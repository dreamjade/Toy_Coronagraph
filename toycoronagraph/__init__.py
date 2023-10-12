import os

# Define the version of the module
__version__ = "1.6.2"

# Set Python environment variable to keep track of example data directory
toycoronagraph_dir = os.path.dirname(__file__)

# Create the full path to the example data directory using the module's location
DATADIR = os.path.join(toycoronagraph_dir, "example_data/")