import os

__version__ = "1.2.3"

# set Python env variable to keep track of example data dir
toycoronagraph_dir = os.path.dirname(__file__)
DATADIR = os.path.join(toycoronagraph_dir, "example_data/")
