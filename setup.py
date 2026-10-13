from setuptools import setup, find_packages, Extension
import re

# auto-updating version code stolen from RadVel
def get_property(prop, project):
    result = re.search(
        r'{}\s*=\s*[\'"]([^\'"]*)[\'"]'.format(prop),
        open(project + "/__init__.py").read(),
    )
    return result.group(1)

def get_requires():
    """The requirements in requirements.txt (blank lines and # comments left out)."""
    lines = (line.split("#")[0].strip() for line in open("requirements.txt"))
    return [line for line in lines if line]

setup(
    name="toycoronagraph",
    version=get_property("__version__", "toycoronagraph"),
    description="Toy!!",
    url="https://github.com/dreamjade/Toy_Coronagraph",
    author="dreamjade",
    author_email="",
    license="MIT",
    packages=find_packages(),
    # the bundled example target (Target() with no parameter file); the MCFOST model series in example_data/
    # (186 MB) is too large for PyPI and stays in the repository
    package_data={"toycoronagraph": ["example_data/example.fits"]},
    include_package_data=True,
    #package_data={"": ["kernels/*.cu"]},
    #ext_modules=get_extensions(),
    #include_dirs=[numpy.get_include()],
    #include_package_data=True,
    zip_safe=False,
    classifiers=[
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.11",
    ],
    install_requires=get_requires(),
    # the unit tests and the coverage badge (tests/make_badges.py): pip install -e ".[test]";
    # the documentation (docs/make_docs.py): pip install -e ".[docs]"
    extras_require={"test": ["pytest", "pytest-cov"], "docs": ["sphinx", "pydata-sphinx-theme", "myst-parser"]},
    keywords="Toy Coronagraph"
)
