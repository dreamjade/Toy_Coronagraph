Toy_Coronagraph
===============

A Python package simulating coronagraphic images: a target (a debris disk image, an MCFOST model, a stellar
wind, a coronal mass ejection, planets, or your own source module) through one or more coronagraphs (vortex
coronagraphs computed with hcipy, or the Lazuli telescope computed with prysm).

Install it with ``pip install toycoronagraph``. The `README <https://github.com/dreamjade/Toy_Coronagraph#readme>`_
is the tutorial; ``tests/test.ipynb`` runs every feature. These pages are the reference, made from the docstrings.

There are two ways in:

* **Scenes** (:doc:`scene`, :doc:`observe`): several inputs on one grid, through several coronagraphs, with one
  result per input and coronagraph, saved as FITS files.
* **Target** (:doc:`main`): the original interface, one target image through one vortex coronagraph, with
  planets and planet videos.

.. toctree::
   :maxdepth: 2

   scenes
   optics_psfs
   parameters
   target
   history

* :ref:`genindex`
* :ref:`modindex`
