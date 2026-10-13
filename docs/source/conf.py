"""
Sphinx configuration of the toycoronagraph documentation.

Build it with  python docs/make_docs.py  (see that script). The pages are made from the package's
docstrings by autodoc; napoleon reads their Google-style sections (Args:, Returns:, ...).
"""
import os
import re
import sys
import textwrap

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)                      # document this checkout, not an installed copy
import toycoronagraph

project = "Toy_Coronagraph"
author = "dreamjade"                        # as in setup.py and LICENSE
copyright = "2023-2026, " + author
version = release = toycoronagraph.__version__

extensions = ["sphinx.ext.autodoc", "sphinx.ext.napoleon", "sphinx.ext.viewcode", "myst_parser"]
source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
exclude_patterns = []

autodoc_default_options = {"members": True, "member-order": "bysource", "show-inheritance": True}
autoclass_content = "both"                    # a class's docstring, then its __init__'s
autodoc_typehints = "none"                    # the docstrings state the types
napoleon_google_docstring = True
napoleon_numpy_docstring = False
napoleon_use_ivar = True                      # Attributes: a compact list, like Args:

html_theme = "pydata_sphinx_theme"
html_title = "Toy_Coronagraph %s" % release
html_theme_options = {
    "github_url": "https://github.com/dreamjade/Toy_Coronagraph",
    "show_toc_level": 2,
    "navigation_with_keys": True,
}
html_show_sourcelink = False


# The docstrings are plain text, written to read well in help() and in the source. Three rules make them render
# as written (reStructuredText would otherwise run indented lines together and read |x| or a*b*c as markup):
#   * a block indented under a paragraph (a code example, an aligned table, a list) is shown as written,
#     highlighted when it is Python;
#   * in a Google-style section (Args:, Returns:, ...), an entry ending with ":" over more deeply indented lines
#     (the parts of a returned tuple or dict, the cases of an argument) becomes that entry and a bullet list,
#     and an entry's other lines are joined to it (napoleon would show a Returns entry's as a definition);
#     an attribute entry naming several things ("x, y (float): ...") keeps them all as its name;
#   * everywhere else, \ | * and a word's trailing _ are shown as typed; `name` is shown as code.
# napoleon then turns the sections into parameter lists (it runs after this hook).
default_role = "code"
SECTION = re.compile(r"\s*(Args|Arguments|Attributes|Examples?|Keyword Args|Keyword Arguments|Methods|Notes?|"
                     r"Other Parameters|Parameters|Raises?|References|Returns?|See Also|Todo|Warnings?|Warns|Yields?):\s*$")
RETURNS = re.compile(r"\s*(Returns?|Yields?):\s*$")
NAMES = re.compile(r"^(\s*)([\w*]+(?:, [\w*]+)+)(?=\s*(\(|:))")     # an entry naming several things: "x, y (float): ..."


def _indent(line):
    return len(line) - len(line.lstrip())


def _plain(line):
    """The line with reStructuredText's markup characters escaped, outside `code` (a leading bullet "* " is kept)."""
    body = line.lstrip()
    head = line[:len(line)-len(body)]
    if body.startswith("* "):
        head, body = head + "* ", body[2:]
    parts = re.split(r"(`[^`]*`)", body)
    for k in range(0, len(parts), 2):
        parts[k] = re.sub(r"_(?![A-Za-z0-9])", r"\\_", parts[k].replace("\\", "\\\\").replace("|", "\\|").replace("*", "\\*"))
    return head + "".join(parts)


def _sections(lines):
    """The lines with the sections' nested parts made bullet lists and the sections' text escaped."""
    out, i, header = [], 0, None                    # header: (indent, name) of the section being read
    while i < len(lines):
        line, indent = lines[i], _indent(lines[i])
        if line.strip() and header and indent <= header[0]:
            header = None
        if line.strip() and SECTION.match(line):
            header = (indent, line.strip()[:-1])
            out.append(line)
            i += 1
            continue
        if header is None or not line.strip():
            out.append(line)
            i += 1
            continue
        if header[1] == "Attributes":                # Sphinx would read "x," as the type of y (napoleon splits Args itself)
            line = NAMES.sub(lambda m: m.group(1) + m.group(2).replace(", ", ",\u200b"), line)
        end = i+1                                    # the lines indented under this one
        while end < len(lines) and (not lines[end].strip() or _indent(lines[end]) > indent):
            end += 1
        nested = [l for l in lines[i+1:end] if l.strip()]
        if not line.rstrip().endswith(":") or not nested:      # an entry and its hanging lines: one line
            out.append(_plain(" ".join([line.rstrip()] + [l.strip() for l in nested])))
            i = end
            continue
        items, n = [], _indent(nested[0])
        for l in nested:                             # a part per line at the parts' indent; deeper lines continue it
            if _indent(l) <= n or not items:
                items.append(l.strip())
            else:
                items[-1] += " " + l.strip()
        at = " "*(indent if RETURNS.match(header[1] + ":") else n)    # napoleon reads a Returns list at the entry's indent
        out += [_plain(line), ""] + [_plain(at + "* " + (t[2:] if t[:2] in ("- ", "* ") else t)) for t in items] + [""]
        i = end
    return out


def _blocks(lines):
    """The lines with each block indented under a paragraph (outside the sections) made a literal block."""
    out, i, prose, section = [], 0, None, None     # prose: indent of the text; section: indent of its header
    while i < len(lines):
        line = lines[i]
        if not line.strip():
            out.append(line)
            i += 1
            continue
        indent = _indent(line)
        if section is not None and indent > section:
            out.append(line)
            i += 1
            continue
        section = indent if SECTION.match(line) else None
        if prose is None or indent <= prose or section is not None:
            prose = indent
            out.append(line if section is not None else _plain(line))
            i += 1
            continue
        end = i                                     # the block: down to the next line back at the text's indent
        while end < len(lines) and (not lines[end].strip() or _indent(lines[end]) > prose):
            end += 1
        while not lines[end-1].strip():
            end -= 1
        if out[-1].strip():
            out.append("")
        try:                                        # Python code is highlighted as Python, the rest shown plain
            compile(textwrap.dedent("\n".join(lines[i:end])), "<docstring>", "exec")
            out += [" "*prose + ".. code-block:: python", ""]
        except SyntaxError:
            out += [" "*prose + ".. code-block:: none", ""]
        out.extend(lines[i:end])
        if end < len(lines) and lines[end].strip():
            out.append("")
        i = end
    return out


def plain_text_docstring(app, what, name, obj, options, lines):
    """autodoc-process-docstring: the three rules above, in place."""
    lines[:] = _blocks(_sections(lines))


def setup(app):
    app.connect("autodoc-process-docstring", plain_text_docstring, priority=400)   # before napoleon (500)
