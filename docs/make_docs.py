"""
Build the documentation: the HTML pages in docs/ from the sources in docs/source/.

    pip install -e ".[docs]"        # once: sphinx, pydata-sphinx-theme, myst-parser
    python docs/make_docs.py

The pages are made from the package's docstrings by Sphinx's autodoc (napoleon reads their Google-style
Args:/Returns: sections). docs/source/ holds conf.py, index.rst (the contents), one page per module
(<module>.rst: a title and ".. automodule:: toycoronagraph.<module>") and history.md. GitHub Pages
serves the built docs/ folder (repository Settings > Pages: branch main, folder /docs).

Rules:
  * every module of the package needs a page; the build stops naming any module without one (add
    docs/source/<module>.rst and list it in index.rst);
  * everything in docs/ except source/ and this script is build output: after a successful build it
    is replaced as a whole, so pages of removed modules do not linger (a failed build changes nothing);
  * Sphinx's warnings (mostly docstrings that do not render as written) are printed, and should be fixed.
"""
import os
import sys
import glob
import shutil
import pkgutil
import tempfile
import subprocess

DOCS = os.path.dirname(os.path.abspath(__file__))
SOURCE = os.path.join(DOCS, "source")
ROOT = os.path.dirname(DOCS)
KEEP = {"source", os.path.basename(__file__)}


def undocumented():
    """The package's modules that no page in docs/source/ documents (with ".. automodule::")."""
    text = ""
    for page in glob.glob(os.path.join(SOURCE, "*.rst")):
        with open(page) as f:
            text += f.read()
    package = os.path.join(ROOT, "toycoronagraph")
    sys.path.insert(0, ROOT)                                   # subpackages are imported to list their modules
    modules = [m.name for m in pkgutil.walk_packages([package], "toycoronagraph.")]
    return [m for m in modules if ".. automodule:: %s\n" % m not in text + "\n"]


def build(out):
    """Run Sphinx (HTML) from docs/source/ into out. Returns: the number of warnings."""
    with tempfile.TemporaryDirectory() as tmp:
        warnings = os.path.join(tmp, "warnings.txt")
        run = subprocess.run([sys.executable, "-m", "sphinx", "-b", "html", "-q", "-E", "-d", os.path.join(tmp, "doctrees"),
                              "-w", warnings, SOURCE, out], cwd=ROOT, env=dict(os.environ, MPLBACKEND="Agg"))
        if run.returncode:
            sys.exit("Sphinx failed (above): docs/ is not changed.")
        with open(warnings) as f:
            return sum(1 for line in f if "WARNING" in line or "ERROR" in line)


if __name__ == "__main__":
    missing = undocumented()
    if missing:
        sys.exit("No page documents %s: add docs/source/<module>.rst with \".. automodule:: <module>\" and list it in "
                 "docs/source/index.rst." % ", ".join(missing))
    with tempfile.TemporaryDirectory() as tmp:
        html = os.path.join(tmp, "html")
        n_warnings = build(html)
        for name in os.listdir(DOCS):
            if name not in KEEP:
                path = os.path.join(DOCS, name)
                shutil.rmtree(path) if os.path.isdir(path) and not os.path.islink(path) else os.remove(path)
        for name in os.listdir(html):
            shutil.move(os.path.join(html, name), os.path.join(DOCS, name))
    open(os.path.join(DOCS, ".nojekyll"), "w").close()       # GitHub Pages: serve _static/ etc. as they are
    print("docs/ built: open docs/index.html%s" % ("" if not n_warnings else
          "; %d warnings (above) to fix" % n_warnings))
