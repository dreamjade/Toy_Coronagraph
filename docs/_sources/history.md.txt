# Version History

## v1.5.2
published on Aug 28, 2023.

**able to add moving planets**

## v1.5.3
published on Sep 8, 2023.

**normalization of psfs**

**Users can now set the number of cores to use for parallel computing**

**feature allow user to ignore dust inside IWA**

**show the planet brightness unit**

**print planet brightness and dust brightness**

**turn on/off dust inside IWA**

**auto check multithreading module before using it in psf.py**

## v1.5.4
published on Sep 12, 2023.

**flip or not flip y-axis**

**you could choose lyot_mask_size in toycoronagraph_para.py**

**check planet order inputs are positive integers**

## v1.6.0
published on Sep 14, 2023.

**core_throughput plot**

**video of planet moving**

**optimate the parameter number**

## v1.6.1
published on Sep 15, 2023.

**complete python docstrings for all new features**

## v1.6.2
published on Sep 16, 2023.

**add_planet(res=...): number of points used for orbit calculations**

**example.fits is included in the PyPI package**

**docstring fixes (flip, res, return values)**

## v2.0.0
published on Oct 12, 2026.

**Scenes: several inputs on one grid, centred on the star (Scene), each made by a source module: a FITS image, an MCFOST model, a stellar wind (Thomson scattering), a coronal mass ejection (ported from astro_toolbox's cme_imaging), planets, or your own module**

**observe(): a scene through one or more coronagraphs; the Result holds the final image of every input through every coronagraph, with plots (final images, radial profiles, each input's share as 2-D maps, two sets of PSFs compared) and FITS files that load back**

**output units: physical (e.g. Jy/arcsec^2) or relative to the star (contrast per pixel, per arcsec^2, or to the peak of the star's unocculted PSF)**

**the Lazuli telescope's PSFs, computed with prysm: the telescope alone or with a vortex coronagraph**

**optics backends: a new optical system plugs in with register_backend**

**inner working angle of every coronagraph (Result.iwa: where a point source keeps half of its light), marked on profile plots and kept in the run's FITS files**

**plot settings shared by every scene and result plot (unit, log, log_range, axis, dpi), with defaults in [display]; other keywords go to matplotlib (e.g. cmap, lw)**

**spatial method for any target: PSFs that change across the image, sampled where their shape changes, with the kernel size chosen for an error target and an error estimate for each image; "auto" picks the circular or spatial method from the target's and the PSFs' symmetry**

**circular method: the light beyond the image edge, out to the corners, is kept**

**TOML parameter files: one home per quantity, every value's origin recorded, a JSON schema for editors; python -m toycoronagraph new | describe | migrate | schema | reference; old toycoronagraph_para.py files are still read**

**CME time series (observe.video) and planet videos as lossless WebM**

**faster reruns: radial PSFs memory-mapped and the spatial method's models cached on disk, reused only if everything they were made from is the same**

**a warning when the pupil is sampled too coarsely (fewer than 5 pixels per lambda/D)**

**unit tests (1,280) with a coverage badge; documentation rebuilt for every module**

**needs Python 3.11 or newer**
