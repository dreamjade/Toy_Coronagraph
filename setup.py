from setuptools import setup, find_packages, Extension
import re

# auto-updating version code stolen from RadVel
def get_property(prop, project):
    result = re.search(
        r'{}\s*=\s*[\'"]([^\'"]*)[\'"]'.format(prop),
        open(project + "/__init__.py").read(),
    )
    return result.group(1)

REPO = "https://github.com/dreamjade/Toy_Coronagraph"

def get_long_description():
    """
    README.md for the PyPI page. PyPI cannot follow the README's relative links, so they point to GitHub:
    images (./static/...) to their raw files, other files to their GitHub pages, and #section links to the
    README's sections on GitHub.
    """
    def absolute(match):
        before, path = match.groups()
        base = "https://raw.githubusercontent.com/dreamjade/Toy_Coronagraph/main/" if path.startswith("static/") \
            else REPO + "/blob/main/"
        return before + base + path
    text = open("README.md", encoding="utf-8").read()
    text = re.sub(r'(\]\(|src="|href=")\./([^)"]+)', absolute, text)
    return text.replace("](#", "](%s#" % REPO)

def get_requires():
    """The requirements in requirements.txt (blank lines and # comments left out)."""
    lines = (line.split("#")[0].strip() for line in open("requirements.txt"))
    return [line for line in lines if line]

setup(
    name="toycoronagraph",
    version=get_property("__version__", "toycoronagraph"),
    description="Simulated coronagraphic images: scenes of dust disks, stellar winds, CMEs and planets through "
                "vortex coronagraphs and the Lazuli telescope",
    long_description=get_long_description(),
    long_description_content_type="text/markdown",
    url=REPO,
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
