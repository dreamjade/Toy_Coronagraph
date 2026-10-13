"""
Target: the original one-image interface (one target image through one coronagraph).

    target = Target()                         # toycoronagraph_para.toml in the working directory
    target.plot_origin()                      # the input image -> origin.png
    target.add_planet([0.5, 0.0], 5e-5, mode="cartesian")
    final = target.plot_final(charge=6)       # the final image -> charge6_with_planets_final.png
    target.contrast(charge=6)                 # planet and dust light in the planet's core (Jy)
    target.planet_video(charge=6)             # a moving planet over one period -> planet_video.webm

The parameters come from the parameter file (toycoronagraph.config, either layout); the
image is in Jy/arcsec^2 on target_pixel x target_pixel pixels of psf_scale arcsec, star at
pixel (N-1)/2. The coronagraph is [coronagraph] type: a vortex (charge: positive, even) or the
Lazuli telescope (toycoronagraph.lazuli; charge 0 = no coronagraph, or an even vortex charge; odd
image size). The PSFs and final images are made and cached by toycoronagraph.observe's
Coronagraph (shared with observe()). Planet positions are (radius arcsec, angle rad from +x
towards +y), with y increasing with the image row as in FITS and in scenes
(tool.arcsec_to_pixel); plots show +y upwards (flip=False turns the y axis over). Planet
brightnesses are stored per pixel area (Jy/arcsec^2 = Jy / psf_scale^2).
"""
import os
import numbers
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from toycoronagraph.psf import cir_psf_contrast, spatial_psf_contrast, radial_source_positions
from toycoronagraph.planet import planet_position, orbit_position, orbit_plot
from toycoronagraph.tool import convert_to_polar, is_positive_even_integer, is_planet_pos_allowed, arcsec_to_pixel
from toycoronagraph.config import load_config, write_used, LAZULI_KEYS
from toycoronagraph.inputs import read_target, target_path
from toycoronagraph import selection
from toycoronagraph.video import write_video

# How to make the image size odd for a lazuli coronagraph with the older layout (see observe.GRID_SIZE_FIX)
LEGACY_SIZE_FIX = ("use an image of %d or %d pixels: this layout does not resample it ([target] target_pixel follows "
                   "the image; [grid] size of the scene layout resamples it)")

class Target(object):
    """A target: a dust ring, a debris disk, ice remains, ... (not necessarily circularly symmetric).

    The parameters come from toycoronagraph_para.toml (see toycoronagraph.config). The
    final image is made with the circular method when the target and the PSFs are both
    circularly symmetric, otherwise with the spatial method (see check_symmetry); the
    parameter file, or the method argument, can choose either.

    Attributes:
        par (SimpleNamespace): The parameters, resolved against the input image (target_pixel N,
            psf_scale in arcsec/pixel, wavelength in m, charge, method, ...).
        file_name (str): The FITS file of the target.
        scene (Scene or None): The scene the image was read through (scene-layout parameter
            files only; see toycoronagraph.scene).
        xpix (np.ndarray): (N,) float64 pixel-centre offsets from the star along the columns, in
            arcsec (x).
        ypix (np.ndarray): (N,) float64 pixel-centre offsets from the star along the rows, in arcsec,
            in row order (equal to xpix; y increases with the row, see tool.arcsec_to_pixel).
        extent (list): [left, right, bottom, top] pixel edges in arcsec, for imshow.
        pre_img (np.ndarray): (N, N) float64 target image in Jy/arcsec^2.
        planets (list): Each planet's position (radius arcsec, angle rad).
        orbits (list): Each planet's orbit: np.ndarray [a (arcsec), e, pa (deg), inc (deg), t (orbital
            periods)] for moving planets, None for static ones.
        planets_brightness (list): Each planet's brightness per pixel area in Jy/arcsec^2 (Jy / psf_scale^2).
        output_prefix (str): Prefix of saved figures and videos: "<target name>_", "" for the example.
    """
    def __init__(self, para_file=None):
        """ Initialization

        Initializes the instance based on the parameter file and the target fits file; writes
        the values used to <parameter file>_used.toml (the parameter file is not modified).

        Args:
            para_file (str): Parameter file; default toycoronagraph_para.toml in the working directory.

        Raises:
            ValueError: The image cannot be read or put in Jy/arcsec^2 (e.g. an unknown unit,
                or a contrast image without the star's flux), or the coronagraph cannot work on
                it (lazuli: an even image size, too few pupil samples for its field, no light at
                its wavelengths, an OPD map that is not valid; see observe.optics_for).
            OSError: A lazuli coronagraph's OPD map cannot be read.
            ImportError: A lazuli coronagraph without prysm or the Lazuli config packages.
        """
        # Read the latest parameters (toycoronagraph_para.toml in the working directory)
        self.par = load_config(para_file)

        # Read the target image (in Jy/arcsec^2); the parameters left "auto" are filled from it,
        # and the values used are saved to <parameter file>_used.toml (the parameter file is not modified).
        # Files in the scene layout are read through toycoronagraph.scene (resampling, units, checks).
        self.file_name = target_path(self.par)
        self.scene = None
        if getattr(self.par, "params", None) is not None:
            F_v = self._image_from_scene()
        else:
            F_v = read_target(self.par)
            write_used(self.par)

        # Pixel centres in arcsec (star at pixel (N-1)/2) and the image extent (pixel edges)
        n, scale = self.par.target_pixel, self.par.psf_scale
        self.xpix = (np.arange(n)-(n-1)/2)*scale
        self.ypix = (np.arange(n)-(n-1)/2)*scale
        self.extent = [-n/2*scale, n/2*scale, -n/2*scale, n/2*scale]

        self.pre_img = F_v.astype(np.float64)
        self.planets = []
        self.orbits = []
        self.planets_brightness = []
        self._image_symmetry = None  # (bool, reason) of pre_img, measured once (see _choose)
        self._last_choice = None  # the method messages printed last (printed again only when they change)
        self._reported = set()    # images whose estimated error was printed already
        # Saved figures and videos start with the target's name (none for the bundled example)
        stem = os.path.basename(self.file_name).split(".")[0]
        self.output_prefix = "" if self.par.fits_filename is None else stem+"_"
        self._optics_of = {}  # charge -> Optics (see _optics)
        self._optics(self.par.charge)  # built now, so the optics are checked against the image at once

    def _image_from_scene(self):
        """The target image through the scene framework.

        Renders the scene of a scene-layout parameter file, fills the flat parameters from the
        scene grid (target_pixel, psf_scale, wavelength, diameter, D, visual_range, psf_range,
        rot_number if <= 0, F_transfer, input_description, input_unit), prints the input summary
        and notes, and writes the values used (<parameter file>_used.toml). Sets self.scene.

        Returns:
            np.ndarray: (N, N) total scene image in Jy/arcsec^2.

        Raises:
            ValueError: The scene is in contrast (the star's flux is unknown).
        """
        from toycoronagraph.scene import Scene
        from toycoronagraph import config
        scene = Scene.from_params(self.par.params)
        scene.render(verbose=False)
        if scene.internal_unit != "Jy/arcsec^2":
            raise ValueError("Target works in Jy/arcsec^2, but the image is in contrast and the star's flux is unknown: "
                             "set [star] temperature (with radius and distance) or flux_jy")
        g, par, p = scene.grid, self.par, scene.params
        par.target_pixel, par.psf_scale, par.wavelength, par.diameter = g.size, g.pixel_scale, g.wavelength, g.diameter
        par.D = g.diameter
        par.visual_range = g.visual_range
        par.psf_range = g.visual_range/2
        if par.rot_number <= 0:
            par.rot_number = g.size
        # nu F_nu (W/m^2/pixel) to F_nu (Jy/arcsec^2), as in config.resolve
        par.F_transfer = config.JY*g.wavelength/config.LIGHT_SPEED/g.pixel_scale**2
        name = next(iter(scene.images))
        par.input_description, par.input_unit = scene.origins[name], str(scene.units[name])
        print("Input: %s, %d x %d, %s; %.5g arcsec/pixel (%s) = %.4g lambda/D across at %.4g m (%s)"
              % (par.input_description, g.size, g.size, par.input_unit, g.pixel_scale, p.origin("grid.pixel_scale"),
                 g.visual_range, g.wavelength, p.origin("observation.wavelength")))
        for note in scene.notes:
            print(note)
        p.write_used(config.used_path(par))
        self.scene = scene
        return scene.total

    # ------------------------------------------------------------ optics, PSFs and method

    def _charge(self, charge):
        """The charge to use.

        Args:
            charge (int or None): The requested charge; None: [coronagraph] charge.

        Returns:
            int: `charge`, or the parameters' charge if it is None or not a positive even integer
                (lazuli: nor 0) (with a message).
        """
        if charge is None:
            return self.par.charge
        if not (is_positive_even_integer(charge) or (self.par.coronagraph_type == 'lazuli' and charge == 0)):
            print("Charge number is not compatible with %s coronagraph type, use %d from the parameters"
                  % (self.par.coronagraph_type, self.par.charge))
            return self.par.charge
        return charge

    def _optics(self, charge):
        """The Optics of a charge on this target's grid ([coronagraph] type and parameters; observe.optics_for),
        built once per charge.

        Args:
            charge (int): Charge (vortex; lazuli: 0 = the telescope alone).

        Returns:
            Optics: N = target_pixel pixels, psf_range (half-width in lambda/D), and the [coronagraph]
                parameters (lazuli: also the target's wavelength and diameter).

        Raises:
            ParamError: A lazuli coronagraph on an image of even size, or with too few pupil samples for it.
            ValueError, OSError, ImportError: A lazuli coronagraph's settings cannot be resolved (see
                observe.optics_for).
        """
        if charge not in self._optics_of:
            from toycoronagraph.observe import optics_for, GRID_SIZE_FIX
            par = self.par
            values = dict(charge=charge, lyot_mask_size=par.lyot_mask_size, symmetric=par.psf_symmetric,
                          **{key: getattr(par, key) for key in LAZULI_KEYS})
            self._optics_of[charge] = optics_for(par.coronagraph_type, values, par.target_pixel, par.psf_range,
                                                 par.wavelength, par.diameter,
                                                 GRID_SIZE_FIX if self.scene is not None else LEGACY_SIZE_FIX)
        return self._optics_of[charge]

    def _coronagraph(self, charge):
        """The coronagraph of a charge on this target's grid (toycoronagraph.observe.Coronagraph).

        The same object for the same charge and settings, so its PSFs (cached on disk in
        psf_cache_dir and in memory) and the final images it has made are reused.

        Args:
            charge (int): Charge (vortex; lazuli: 0 = the telescope alone).

        Returns:
            Coronagraph: Named "<type>_c<charge>", e.g. "vortex_c6" or "lazuli_c0".
        """
        from toycoronagraph.observe import Coronagraph
        par = self.par
        return Coronagraph.get(self._optics(charge), "%s_c%d" % (par.coronagraph_type, charge), num_cores=par.num_cores,
                               error_target=par.error_target, max_kernel_gb=par.max_kernel_gb, use_gpu=par.use_gpu,
                               cache_dir=par.psf_cache_dir, symmetric=par.psf_symmetric, rot_number=par.rot_number,
                               symmetry_tol=par.symmetry_tol)

    def radial_psfs_file(self, charge=None):
        """Radial PSF file

        The file of the radial PSF series (sources along the +x axis) on frames larger than the
        image, calculated and saved in psf_cache_dir if it does not exist. Its name records
        the optics and sampling, so a cache made with other parameters is never reused.

        Args:
            charge (int): Charge number (vortex; lazuli: 0 = the telescope alone); default from the parameters.

        Returns:
            str: The path of the .npy file (psf_cache_dir/radial_<tag>_frame<size>.npy, a
                ((N+1)//2, size, size) float32 array; see Coronagraph.radial_psfs_file).
        """
        return self._coronagraph(self._charge(charge)).radial_psfs_file()

    def radial_psfs(self, charge=None):
        """The radial PSF series (see radial_psfs_file).

        Args:
            charge (int): Charge number (vortex; lazuli: 0 = the telescope alone); default from the parameters.

        Returns:
            np.ndarray: ((N+1)//2, size, size) float32 PSFs of sources along +x, size > N.
        """
        return np.load(self.radial_psfs_file(charge))

    def choose_method(self, charge=None, method=None, spatial_mode=None, verbose=True):
        """Method choice

        The method for the final image: method / spatial_mode if given, else the parameter
        file's; "auto" uses the circular method only when the target and the PSFs are both
        circularly symmetric ([target] / [coronagraph] symmetric, or measured / from the optics model).

        Args:
            charge (int): Charge number (vortex; lazuli: 0 = the telescope alone); default from the parameters.
            method (str): "auto", "circular" or "spatial"; default from the parameters.
            spatial_mode (str): "auto", "rotate" or "direct"; default from the parameters.
            verbose (bool): Print the choice and why.

        Returns:
            tuple: method (str) and spatial_mode (str; None for the circular method).

        Raises:
            ValueError: An unknown method / spatial_mode, or "rotate" for optics that do not allow it.
        """
        charge = self._charge(charge)
        method, spatial_mode, messages = self._choose(charge, method, spatial_mode)
        if verbose:
            print("\n".join(messages))
            self._last_choice = messages
        return method, spatial_mode

    def _method(self, charge, method, spatial_mode):
        """choose_method for the image makers: prints the choice only when it changes.

        Args:
            charge (int): Charge (already checked).
            method (str or None): "auto", "circular" or "spatial"; None: from the parameters.
            spatial_mode (str or None): "auto", "rotate" or "direct"; None: from the parameters.

        Returns:
            tuple: method (str) and spatial_mode (str; None for the circular method).
        """
        method, spatial_mode, messages = self._choose(charge, method, spatial_mode)
        if self._last_choice != messages:
            print("\n".join(messages))
            self._last_choice = messages
        return method, spatial_mode

    def _choose(self, charge, method, spatial_mode):
        """The method, spatial mode and the messages explaining the choice (see selection.choose_method).

        The target's symmetry is measured (or taken from [target] symmetric) once per Target.

        Args:
            charge (int): Charge (already checked).
            method (str or None): "auto", "circular" or "spatial"; None: from the parameters.
            spatial_mode (str or None): "auto", "rotate" or "direct"; None: from the parameters.

        Returns:
            tuple: method (str), spatial_mode (str or None for circular), messages (list of str).
        """
        optics = self._optics(charge)
        if self._image_symmetry is None:
            self._image_symmetry = selection.image_symmetry(self.pre_img, self.par.image_symmetric, self.par.symmetry_tol)
        image_sym, image_why = self._image_symmetry
        psf_sym, psf_why = self._coronagraph(charge).psf_symmetry()
        backend = optics.backend()
        # 'rotate' needs rotationally symmetric optics with a known band limit
        return selection.choose_method(
            method or self.par.method, spatial_mode or self.par.spatial_mode, image_sym, psf_sym,
            bool(backend.rotationally_symmetric and backend.bandlimit), (image_why, psf_why))

    def check_symmetry(self, charge=None):
        """Symmetry check

        Print whether the target and the PSFs are circularly symmetric and which method the
        final image will use (see choose_method).

        Args:
            charge (int): Charge number (vortex; lazuli: 0 = the telescope alone); default from the parameters.

        Returns:
            tuple: method (str) and spatial_mode (str; None for the circular method).
        """
        return self.choose_method(charge)

    def spatial_model(self, charge=None, spatial_mode="rotate"):
        """Spatial PSF model

        The spatially-sampled PSFs for this target (built once per charge and mode): the
        kernel size and the grid are chosen so the final image meets [method] error_target;
        the runs along the x axis are cached in psf_cache_dir. Prints the error analysis.

        Args:
            charge (int): Charge number (vortex; lazuli: 0 = the telescope alone); default from the parameters.
            spatial_mode (str): "rotate" or "direct".

        Returns:
            dict: The model (see psf.spatial_psf_calculation).
        """
        return self._coronagraph(self._charge(charge)).spatial_model(spatial_mode, self.pre_img)

    def _disk_and_planets(self, charge, iwa_ignore, add_planet, method, spatial_mode):
        """The final image (disk, plus planets if add_planet) with the chosen method.

        The disk's estimated error (spatial method) is printed once per (charge, mode, iwa_ignore).

        Args:
            charge (int): Charge (already checked).
            iwa_ignore (bool): Leave out the target inside the IWA.
            add_planet (bool): Add self.planets (see _planets_image).
            method (str): "circular" or "spatial" (already chosen).
            spatial_mode (str or None): "rotate" or "direct" for the spatial method.

        Returns:
            np.ndarray: (N, N) float64 final image in Jy/arcsec^2.
        """
        coro = self._coronagraph(charge)
        final_img, _ = coro.image(self.pre_img, method, spatial_mode, iwa_ignore=iwa_ignore,
                                  report=method == "spatial" and self._first_report((charge, spatial_mode, iwa_ignore)))
        if add_planet and self.planets:
            final_img = final_img+self._planets_image(self.planets, self.planets_brightness, charge, method)
        return final_img

    def _planets_image(self, positions, brightness, charge, method):
        """Planets through the coronagraph: the radial PSFs rotated (circular), else their exact PSFs.

        Args:
            positions (list): Each planet's (radius arcsec, angle rad).
            brightness (list): Each planet's brightness per pixel area in Jy/arcsec^2.
            charge (int): Charge (already checked).
            method (str): "circular" or "spatial".

        Returns:
            np.ndarray: (N, N) image of the planets in Jy/arcsec^2 (planets off the image are left out).
        """
        coro = self._coronagraph(charge)
        if method == "circular":
            return coro.interpolated_point_sources(positions, brightness, self.par.psf_scale)
        return coro.point_sources([arcsec_to_pixel(pos, self.par.psf_scale, self.par.target_pixel) for pos in positions],
                                  brightness)

    def _first_report(self, key):
        """Whether the estimated error of this disk image is still to be printed (once per image).

        Args:
            key (tuple): Identifies the image, e.g. (charge, spatial_mode, iwa_ignore).

        Returns:
            bool: True the first time the key is seen (it is then recorded).
        """
        first = key not in self._reported
        self._reported.add(key)
        return first

    def _save_figure(self, fig, file_name, default_name):
        """Save a figure exactly at file_name, or under the default name.

        Args:
            fig (matplotlib.figure.Figure): The figure.
            file_name (str or os.PathLike or None): "name" or "folder/name" (the folder is created;
                the format follows the extension, PNG if it is not a known one); a message gives
                the path. None: <output_prefix><default_name>.png, without a message.
            default_name (str): The default file name, without the extension.

        Returns:
            str: The path written.
        """
        if file_name is None:
            path = self.output_prefix+default_name+".png"
        else:
            path = os.fspath(file_name)
            folder = os.path.dirname(path)
            if folder:
                os.makedirs(folder, exist_ok=True)
        # The format follows the extension (.png, .pdf, .svg, ...); PNG when there is none
        known = os.path.splitext(path)[1].lstrip(".").lower() in fig.canvas.get_supported_filetypes()
        fig.savefig(path, format=None if known else "png", bbox_inches='tight')
        if file_name is not None:
            print("Saved the image to %s" % path)
        return path

    def _colorbar_label(self):
        """The colour bar label: [target] unit, or "Jy/arcsec^2" when it is "none".

        Returns:
            str: The label.
        """
        return self.par.unit if self.par.unit != "none" else "Jy/arcsec^2"

    def _optics_text(self, charge):
        """The coronagraph in words, after "the" in figure titles, e.g. "charge-6 vortex" or "Lazuli telescope alone".

        Args:
            charge (int): The charge.

        Returns:
            str: The text.
        """
        if self.par.coronagraph_type == "vortex":
            return "charge-%d vortex" % charge
        return "Lazuli telescope with a charge-%d vortex" % charge if charge else "Lazuli telescope alone"

    @staticmethod
    def _colour_scale(image, log, log_range):
        """The colour map and scale of an image plot.

        Args:
            image (np.ndarray): The image shown.
            log (bool): Log scale, from the image's peak down `log_range` decades (zero or negative
                pixels in the darkest colour); else linear from 0.
            log_range (float): Decades shown on a log scale.

        Returns:
            tuple: (Colormap, Normalize) for imshow.
        """
        cmap = matplotlib.colormaps['gnuplot']
        if not log:
            return cmap, matplotlib.colors.Normalize(vmin=0, vmax=np.max(image))
        top = np.max(image)
        return cmap.with_extremes(bad=cmap(0.0), under=cmap(0.0)), matplotlib.colors.LogNorm(top*10.0**-log_range, top)

    # ------------------------------------------------------------ plots

    def plot_origin(self, plot_planets=True, res=1000, plot_dpi=300, boundary=True, flip=True, file_name=None,
                    log=False, log_range=5):
        """Plot original target

        Plot the original image (Jy/arcsec^2) with or without planets: each planet is a
        coloured circle labelled with its brightness in Jy, with its orbit (dashed) if it moves.
        Planets outside the image are reported, greyed, and drawn only without boundary.

        Args:
            plot_planets (bool): Whether to plot planets or not.
            res (int): Number of points used for orbit calculations.
            plot_dpi (int): Dots per inch (DPI) for the plot.
            boundary (bool): Whether to limit the axes to the image.
            flip (bool): True: +y upwards (row 0 at the bottom, as DS9 shows FITS); False: the y axis turned over.
            file_name (str): Save the figure exactly there: "name.png" or "folder/name.png" (the
                folder is created). The format follows the extension (PNG if there is none). None
                (default): origin.png / origin_with_planets.png (with the target's name in front
                for your own targets).
            log (bool): Log colour scale (faint outer parts visible); else linear from 0.
            log_range (float): Log scale: decades shown below the image's peak.

        Returns:
            None: The figure is saved (origin*.png, see file_name) and shown.
        """
        # Create a new figure and axes for plotting
        fig = plt.figure(dpi=plot_dpi)
        axs = plt.gca()

        # Plot the pre-processed image data
        cmap, norm = self._colour_scale(self.pre_img, log, log_range)
        img = axs.imshow(self.pre_img, cmap=cmap, norm=norm, extent=self.extent, origin='lower')
        marked = plot_planets and self.planets != []
        axs.set_title("Input image%s" % (" (planets: brightness in Jy)" if marked else ""), fontsize=10)

        # Calculate boundary limits based on pixel scale (half the image width, arcsec)
        x_range = self.par.target_pixel*self.par.psf_scale/2
        y_range = self.par.target_pixel*self.par.psf_scale/2

        if boundary:
            axs.set_xlim(-x_range, x_range)
            axs.set_ylim(-y_range, y_range)
        if not flip:
            axs.invert_yaxis()

        # Set axis labels
        axs.set_ylabel('y [arcsec]')
        axs.set_xlabel('x [arcsec]')

        # Add color bar
        colorbar=plt.colorbar(img,orientation='vertical')
        colorbar.set_label(self._colorbar_label())

        # Iterate through planets, their brightness, orbits, and orders to add them into the plot
        if plot_planets and self.planets !=[]:
            circle_size = 0.02*min(x_range, y_range)
            colors = [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#bcbd22', u'#17becf']
            for p, b, orbit, order in zip(self.planets, self.planets_brightness, self.orbits, range(len(self.planets))):
                x,y = p[0]*np.cos(p[1]),p[0]*np.sin(p[1])
                color = colors[order%9]

                outside = abs(x) > x_range or abs(y) > y_range
                if outside:
                    print("Planet #"+str(order+1)+" is outside the range of the plot")
                    color = u'#7f7f7f'

                # A planet outside the image is drawn only when the axes are not limited to it
                if not outside or not boundary:
                    circle = plt.Circle((x,y), radius=circle_size, color=color)
                    axs.add_artist(circle)
                    axs.annotate(str(b*self.par.psf_scale**2), (x,y), color=color)  # brightness in Jy

                if orbit is not None:
                    orbit_x, orbit_y = orbit_position(orbit[0], orbit[1], orbit[2], orbit[3], res)
                    plt.plot(orbit_x, orbit_y, color=color, linestyle = '--', alpha = 0.5) # Plot the orbit.

            default_name = "origin_with_planets"
        else:
            default_name = "origin"
        self._save_figure(fig, file_name, default_name)

        # Show the plot
        plt.show()

    def add_planet(self, pos, brightness, mode="moving", res=1000):
        """Add a planet

        Add a planet to the Target object (appended to planets, planets_brightness and orbits).
        Invalid input (a position that tool.is_planet_pos_allowed refuses, e.g. one with an entry
        that is not a number, a brightness that is not a positive number, an unknown mode) is
        reported with a message and nothing is added.

        Args:
            pos (array-like): Planet position information (depends on mode); each entry a finite real
                number (Python or numpy int or float). For "moving" mode, pos = [a, e, pa, inc, t]: semi-major
                axis (arcsec), eccentricity (0 <= e < 1), position angle (deg), inclination (deg)
                and time (in orbital periods); for "polar" mode, pos = [r (arcsec), theta (deg, from
                +x towards +y)]; for "cartesian" mode, pos = [x, y] (arcsec).
            brightness (numbers.Real): Planet brightness in Jy (> 0): a Python or numpy int or float
                (not a bool).
            mode (str): Mode for specifying planet position coordinates ("moving", "polar", "cartesian").
            res (int): Number of points used for orbit calculations.

        Returns:
            None
        """
        brightness_valid = isinstance(brightness, numbers.Real) and not isinstance(brightness, bool) and brightness > 0

        if mode == "moving":
            # Check if planet position and brightness are valid
            if not is_planet_pos_allowed(pos, mode):
                print("Planet position is not allowed")
            elif not brightness_valid:
                print("Brightness is invalid")
            else:
                pos = np.array(pos, dtype=float)  # (integer times would be truncated when the planet moves)
                self.planets.append(planet_position(pos[0],pos[1],pos[2],pos[3],pos[4],mode="polar",res=res))
                self.planets_brightness.append(float(brightness)/self.par.psf_scale**2)
                self.orbits.append(pos) #pos = np.array(a, e, pa, inc, t)

        elif mode == "polar":
            # Check if planet position and brightness are valid
            if not is_planet_pos_allowed(pos, mode):
                print("Planet position is not allowed")
            elif not brightness_valid:
                print("Brightness is invalid")
            else:
                # Convert polar mode angle from degree to radian (in floats: integers would truncate it)
                pos = np.array(pos, dtype=float)
                pos[1] = pos[1]*np.pi/180.0
                self.planets.append(pos)
                self.planets_brightness.append(float(brightness)/self.par.psf_scale**2)
                self.orbits.append(None)

        elif mode == "cartesian":
            # Check if planet position and brightness are valid
            if not is_planet_pos_allowed(pos, mode):
                print("planet position is not allowed")
            elif not brightness_valid:
                print("Brightness is invalid")
            else:
                self.planets.append(convert_to_polar(pos))
                self.planets_brightness.append(float(brightness)/self.par.psf_scale**2)
                self.orbits.append(None)
        else:
            print("no "+str(mode)+" mode")

    def list_planets(self):
        """List planets

        Print the planets present in the Target object with their order numbers (from 1),
        static or moving, position (x, y) in arcsec and brightness in Jy.

        Returns:
            None
        """
        # Check if the planet list is empty
        if len(self.planets)==0:
            print("There is no planet, but you can add planet via Target.add_planet(pos, brightness)")
        else:
            for order, (p, b, o) in enumerate(zip(self.planets, self.planets_brightness, self.orbits), start=1):
                print("{} Planet {:d}: ({}, {}) arcsec, brightness: {:.2e} Jy".format(
                    "Static" if o is None else "Moving", order, p[0]*np.cos(p[1]), p[0]*np.sin(p[1]), b*self.par.psf_scale**2))

    def delete_planet(self, order=1):
        """Delete planet

        Deletes a planet from self.planets, self.orbits, and self.planets_brightness with a given
        order number, then prints the remaining planets.

        Args:
            order (int): The order number of the planet to be deleted (from 1). Defaults to 1.

        Returns:
            None
        """
        # Check if the input order is valid
        if not isinstance(order, int) or order < 1 or order > len(self.planets):
            print("Input order number is not existed, please try again or use Target.list_planets() to list all planets")
        else:
            # Delete planet, brightness, and orbit information at the given order
            del self.planets[order-1]
            del self.planets_brightness[order-1]
            del self.orbits[order-1]
            print("Successfully remove planet #"+str(order)+", here is the latest planet list:")
            self.list_planets()

    def planet_move(self, time, order=1, mode="cumulate", plot_pos=False, plot_dpi=300, res=1000, flip=True, message=True):
        """Planet moving

        Move a moving planet along its orbit (static planets cannot move).

        Args:
            time (float): The time, in orbital periods.
            order (int): The order number of the planet (from 1).
            mode (str): "cumulate": time is added to the planet's current time; "specific": the
                planet's time is set to `time`. Any other mode leaves the time unchanged.
            plot_pos (bool): Whether to plot the new planet's position in its orbit
                (saved as oribit_planet<order>.png and shown, see planet.orbit_plot).
            plot_dpi (int): Dots per inch (DPI) for the plot.
            res (int): Number of points used for orbit calculations.
            flip (bool): True: +y upwards; False: the y axis turned over.
            message (bool): Whether to show the notification once job done (when not plotting).

        Returns:
            None
        """
        # Check if the input order is valid
        if not isinstance(order, int) or order < 1 or order > len(self.planets):
            print("Input order number is not existed, please try again or use Target.list_planets() to list all planets")
        elif self.orbits[order-1] is None or len(self.orbits[order-1]) != 5:
            print("A static planet is unable to move")
        else:
            if mode == "cumulate":
                # Update the time of the planet's orbit by adding to the current time
                self.orbits[order-1][4] += time
            elif mode == "specific":
                # Set a new time for the planet's orbit
                self.orbits[order-1][4] = time

            # Get the updated position of the planet based on the updated time
            pos = self.orbits[order-1]
            self.planets[order-1] = planet_position(pos[0],pos[1],pos[2],pos[3],pos[4],mode="polar",res=res)
            if plot_pos:
                # Plot the new planet position on its orbit
                orbit_plot(pos[0],pos[1],pos[2],pos[3], self.planets[order-1], plot_dpi, res, flip, "polar", '_planet'+str(order))
            elif message:
                print("Planet has moved to new position")

    def planet_video(self, charge=None, order=1, iwa_ignore=False, plot_dpi=300, length=3, fps=20, flip=True,
                     method=None, spatial_mode=None):
        """Planet Video

        Create a movie of the planet moving over one orbital period, starting from its current
        position; afterwards the planet is back at that position. The disk's final image is
        made once; each frame adds the planet's image at its new position.

        The video is lossless (VP9 lossless in WebM, see toycoronagraph.video): every frame decodes
        to exactly the rendered pixels, so colours can be read against the colour bar. It stays
        small (about 100-200 kB for 60 frames) because the frames differ only around the planet.
        It plays in browsers, VS Code and VLC (not QuickTime). Needs ffmpeg with libvpx-vp9.

        Args:
            charge (int): Charge number (vortex; lazuli: 0 = the telescope alone); default from the parameters.
            order (int): The order number of the planet, using self.list_planets() to look up the order.
            iwa_ignore (bool): Whether to ignore the dust inside IWA.
            plot_dpi (int): Dots per inch of the 10 x 8 inch frames (300: 3000 x 2400 pixels).
            length (int): The length of the movie in seconds.
            fps (int): The frames per second of the movie.
            flip (bool): True: +y upwards; False: the y axis turned over.
            method (str): "auto", "circular" or "spatial"; default from the parameters.
            spatial_mode (str): "auto", "rotate" or "direct"; default from the parameters.

        Returns:
            None: The video is saved as <prefix>planet_video[_iwa_ignore][_spatial].webm in the
                working directory (with the target's name in front for your own targets).
        """
        charge = self._charge(charge)

        # Check if the input order is valid
        if not isinstance(order, int) or order < 1 or order > len(self.planets):
            print("Input order number is not existed, please try again or use Target.list_planets() to list all planets")
        elif self.orbits[order-1] is None or len(self.orbits[order-1]) != 5:
            print("A static planet is unable to move")
        else:
            method, spatial_mode = self._method(charge, method, spatial_mode)
            total_frames = length*fps
            target_img = self._disk_and_planets(charge, iwa_ignore, False, method, spatial_mode)
            fig,[ax,cax] = plt.subplots(1,2, width_ratios=[50,1], figsize=(10, 8), dpi=plot_dpi)
            # Set the colormap and norm to correspond to the data for which the colorbar will be used.
            cmap = matplotlib.colormaps['gnuplot']
            planet_brightness = self.planets_brightness[order-1]
            t0 = self.orbits[order-1][-1]
            initial_planet_pos = self.planets[order-1]
            initial_planet_img = self._planets_image([initial_planet_pos], [planet_brightness], charge, method)
            # The colour scale is fixed by the first frame
            norm = matplotlib.colors.Normalize(vmin=0, vmax=np.max(target_img+initial_planet_img))
            colorbar = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap,norm=norm, orientation='vertical')
            colorbar.set_label(self._colorbar_label())
            img = ax.imshow(target_img+initial_planet_img, extent=self.extent, cmap=cmap, norm=norm, origin='lower')
            if not flip:
                ax.invert_yaxis()
            ax.set_ylabel('y [arcsec]')
            ax.set_xlabel('x [arcsec]')
            ax.set_title("Planet #%d over one orbit, through the %s%s" % (order, self._optics_text(charge),
                         ", dust inside the IWA removed" if iwa_ignore else ""))
            def animate(num):
                """Frame num: the planet at time t0 + num/total_frames periods, added to the disk image.

                Args:
                    num (int): Frame number (0 ... total_frames-1).

                Returns:
                    list: The updated artist (the image).
                """
                self.planet_move(t0+num/total_frames, order=order, mode="specific", res=total_frames*10, message=False)
                current_img = target_img+self._planets_image([self.planets[order-1]], [planet_brightness], charge, method)
                img.set_array(current_img)
                return [img]
            anim_name = self.output_prefix + 'planet_video' + ('_iwa_ignore' if iwa_ignore else '') + ('_spatial' if method == "spatial" else '') + '.webm'
            write_video(fig, animate, total_frames, anim_name, fps, plot_dpi)
            plt.close(fig)
            # The video does not move the planet: back to where it was
            self.orbits[order-1][4] = t0
            self.planets[order-1] = initial_planet_pos
            print("Saved the video to %s" % anim_name)

    def plot_orbit(self, order=1, plot_dpi=300, res=1000, flip=True):
        """Orbit plot

        Plot the orbit of a moving planet, with its current position (planet.orbit_plot).

        Args:
            order (int): The order number of the planet, using self.list_planets() to look up the order.
            plot_dpi (int): Dots per inch (DPI) for the plot.
            res (int): Number of points used for orbit calculations.
            flip (bool): True: +y upwards; False: the y axis turned over.

        Returns:
            None: The figure is saved as oribit_planet<order>.png (sic) and shown.
        """
        if not isinstance(order, int) or not 1 <= order <= len(self.planets):
            print("Input order number is not existed, please try again or use Target.list_planets() to list all planets")
            return
        pos = self.orbits[order-1]
        if pos is None:
            print("A static planet has no orbit")
            return
        orbit_plot(pos[0],pos[1],pos[2],pos[3], self.planets[order-1], plot_dpi, res, flip, "polar", '_planet'+str(order))

    def contrast(self, charge=None, order=1, plot_contrast=True, plot_dpi=300, exozodi_limit=10, res=100,
                 method=None, spatial_mode=None):
        """Contrast

        Find the planet brightness and background brightness in the final image: the planet's
        core (the pixels where its final image is above half its peak) is the aperture, and the
        planet and the dust (with and without the dust inside the IWA) are summed inside it.

        Args:
            charge (int): Charge number (vortex; lazuli: 0 = the telescope alone); default from the parameters.
            order (int): The order number of the planet, using self.list_planets() to look up the order.
            plot_contrast (bool): Whether to plot the dust to planet light ratio in the planet's core
                against a factor scaling the dust's brightness (e.g. more exozodiacal dust; shown, not saved).
            plot_dpi (int): Dots per inch (DPI) for the plot.
            exozodi_limit (float): The largest dust brightness factor plotted.
            res (int): The number of points in the exozodi ratio range.
            method (str): "auto", "circular" or "spatial"; default from the parameters.
            spatial_mode (str): "auto", "rotate" or "direct"; default from the parameters.

        Returns:
            tuple or None: planet brightness, background brightness, and background brightness
                (ignored dust inside IWA), floats in Jy (also printed, and plotted if plot_contrast);
                None for an invalid order, or (circular method) a planet at or beyond the last radial
                PSF's source, which the final image leaves out (see psf.cir_psf_planets).
        """
        # Check if the input order is valid
        if not isinstance(order, int) or order < 1 or order > len(self.planets):
            print("Input order number is not existed, please try again or use Target.list_planets() to list all planets")
            return None
        charge = self._charge(charge)
        if self.par.coronagraph_type == 'vortex':
            print("This is a charge-" + str(charge) + " vortex coronagraph.")
        elif charge:
            print("This is the Lazuli telescope with a charge-%d vortex coronagraph." % charge)
        else:
            print("This is the Lazuli telescope without a coronagraph.")
        method, spatial_mode = self._method(charge, method, spatial_mode)
        planet_pos = self.planets[order-1]
        planet_r = planet_pos[0]/self.par.psf_scale  # in pixels
        if method == "circular":
            # Check if the planet is inside the range of the radial PSFs (the limit of psf.cir_psf_planets;
            # the radial PSF file holds all (N+1)//2 of them)
            if planet_r >= radial_source_positions(self.par.target_pixel)[-1]:
                print("Planet is outside the range of the plot")
                return None
            planet_b, dust_b, dust_b_iwa = cir_psf_contrast(self.pre_img, planet_r, planet_pos[1], self.planets_brightness[order-1],
                                                            self.par.psf_scale, self.par.target_pixel, self.par.rot_number,
                                                            self.radial_psfs_file(charge))
        else:
            planet_b, dust_b, dust_b_iwa = spatial_psf_contrast(self.pre_img, planet_pos, self.planets_brightness[order-1],
                                                                self.par.psf_scale, self.spatial_model(charge, spatial_mode),
                                                                self.par.use_gpu, self.par.num_cores,
                                                                report=self._first_report((charge, spatial_mode, False)))
        if plot_contrast:
            # The dust's final image is linear in its brightness: scaled by a factor, so is the dust in the core
            plt.figure(dpi=plot_dpi)
            axs = plt.gca()
            axs.set_ylabel("dust / planet light in the planet's core")
            axs.set_xlabel("dust brightness relative to the target's (1 = as given)")
            axs.set_title("Planet #%d through the %s" % (order, self._optics_text(charge)), fontsize=10)
            exozodi_ratio = np.linspace(0, exozodi_limit, res)
            dp_ratio = (dust_b / planet_b) * exozodi_ratio
            dp_ratio_iwa = (dust_b_iwa / planet_b) * exozodi_ratio
            plt.plot(exozodi_ratio, dp_ratio, label="all the dust")
            plt.plot(exozodi_ratio, dp_ratio_iwa, label="dust inside the IWA left out")
            plt.axhline(1, color="grey", lw=0.8, ls=":")
            plt.legend()
            plt.show()
        print("Planet brightness: %.2e Jy" % planet_b)
        print("Dust brightness: %.2e Jy" % dust_b)
        print("Dust brightness but ignore dust inside IWA: %.2e Jy" % dust_b_iwa)
        return planet_b, dust_b, dust_b_iwa

    def plot_final(self, charge=None, iwa_ignore=False, add_planet=True, plot_dpi=300, flip=True,
                   method=None, spatial_mode=None, file_name=None, log=False, log_range=5):
        """Final image

        Plots the target image after processing with the coronagraph, with the circular or the
        spatial method (see choose_method).

        Args:
            charge (int): Charge number (vortex; lazuli: 0 = the telescope alone); default from the parameters.
            iwa_ignore (bool): Whether to ignore the target inside the IWA.
            add_planet (bool): Whether to include planets in the plot.
            plot_dpi (int): Dots per inch (DPI) for the final plot.
            flip (bool): True: +y upwards; False: the y axis turned over.
            method (str): "auto", "circular" or "spatial"; default from the parameters.
            spatial_mode (str): "auto", "rotate" or "direct"; default from the parameters.
            file_name (str): Save the figure exactly there: "name.png" or "folder/name.png" (the
                folder is created). The format follows the extension (PNG if there is none). None
                (default): charge<c>[_with_planets][_iwa_ignore][_spatial]_final.png ("lazuli_charge<c>..."
                for the Lazuli telescope; with the target's name in front for your own targets).
            log (bool): Log colour scale (faint outer parts visible); else linear from 0.
            log_range (float): Log scale: decades shown below the image's peak.

        Returns:
            np.ndarray: The (N, N) float64 final image in Jy/arcsec^2 (also saved, see file_name, and shown).
        """
        charge = self._charge(charge)
        method, spatial_mode = self._method(charge, method, spatial_mode)

        # Generate the final image of the disk (and the planets)
        final_img = self._disk_and_planets(charge, iwa_ignore, add_planet, method, spatial_mode)

        # Create a new figure and axes for plotting
        fig = plt.figure(dpi=plot_dpi)
        axs = plt.gca()

        # Plot the final image
        cmap, norm = self._colour_scale(final_img, log, log_range)
        img = axs.imshow(final_img, cmap=cmap, norm=norm, extent=self.extent, origin='lower')
        axs.set_title("Through the %s (%s method)%s%s" % (self._optics_text(charge), method,
                                                         ", with planets" if add_planet and self.planets else "",
                                                         ", dust inside the IWA removed" if iwa_ignore else ""), fontsize=10)
        if not flip:
            axs.invert_yaxis()

        axs.set_ylabel('y [arcsec]')
        axs.set_xlabel('x [arcsec]')

        # Add color bar
        colorbar=plt.colorbar(img,orientation='vertical')
        colorbar.set_label(self._colorbar_label())

        # Create a name for the final image
        final_image_name = ("" if self.par.coronagraph_type == "vortex" else "lazuli_")+"charge"+str(charge)
        if add_planet and self.planets != []:
            final_image_name += "_with_planets"

        if iwa_ignore:
            final_image_name += "_iwa_ignore"

        if method == "spatial":
            final_image_name += "_spatial"

        # Save the final image
        self._save_figure(fig, file_name, final_image_name+"_final")

        # Show the plot
        plt.show()
        return final_img
