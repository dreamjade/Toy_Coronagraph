Scenes
======

A scene is the image before the coronagraph: several inputs, each made by a source module, on one grid centred on
the star (:doc:`scene`). :doc:`observe` sends it through one or more coronagraphs and returns a Result, with the
final image of every input through every coronagraph. Write your own input as a source module (:doc:`sources`).

.. toctree::
   :maxdepth: 1

   scene
   sources
   sources.images
   sources.thomson
   sources.cme
   sources.planets
   observe
   units
   selection
