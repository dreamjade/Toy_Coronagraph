"""
Lossless videos of matplotlib figures (VP9 lossless in WebM).

Every frame decodes to exactly the rendered pixels, so colours can be read against the colour
bar. The files stay small when consecutive frames differ little, since VP9 codes the
differences. They play in browsers, VS Code and VLC (not QuickTime). Needs ffmpeg with libvpx-vp9.
Used by Target.planet_video (toycoronagraph.main) and observe.video.
"""
import os
import matplotlib.animation as animation
import tqdm

# Lossless VP9 on RGB (gbrp: no colour conversion); cpu-used 0 gives the smallest file
LOSSLESS_VP9 = ['-c:v', 'libvpx-vp9', '-lossless', '1', '-pix_fmt', 'gbrp', '-deadline', 'good', '-cpu-used', '0',
                '-row-mt', '1']


def write_video(fig, update, frames, path, fps=20, dpi=100):
    """
    Write a figure animated frame by frame to a lossless WebM file.

    Args:
        fig (matplotlib.figure.Figure): The figure (not closed here).
        update (callable): update(i) changes the figure for frame i (0 ... frames-1); it may be
            called more than once for a frame.
        frames (int): Number of frames.
        path (str): The .webm file to write.
        fps (int): Frames per second.
        dpi (int): Dots per inch of the frames.

    Returns:
        str: path.
    """
    anim = animation.FuncAnimation(fig, update, frames=frames, interval=1000/fps)
    with tqdm.tqdm(total=frames, desc=os.path.basename(path), leave=True, mininterval=2) as bar:
        anim.save(path, fps=fps, dpi=dpi, extra_args=LOSSLESS_VP9, progress_callback=lambda i, n: bar.update(1))
    return path
