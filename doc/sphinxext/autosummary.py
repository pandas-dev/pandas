"""
===========
autosummary
===========

Sphinx extension that adds an autosummary:: directive, which can be
used to generate function/method/attribute/etc. summary lists, similar
to those output eg. by Epydoc and other API doc generation tools.

An :autolink: role is also provided.

autosummary directive
---------------------

The autosummary directive has the form::

    .. autosummary::
       :nosignatures:
       :toctree: generated/
       
       module.function_1
       module.function_2
       ...

and it generates an output table (containing signatures, optionally)

    ========================  =============================================
    module.function_1(args)   Summary line from the docstring of function_1
    module.function_2(args)   Summary line from the docstring
    ...
    ========================  =============================================

If the :toctree: option is specified, files matching the function names
are inserted to the toctree with the given prefix:

    generated/module.function_1
    generated/module.function_2
    ...

Note: The file names contain the module:: or currentmodule:: prefixes.

.. seealso:: autosummary_generate.py


autolink role
-------------

The autolink role functions as ``:obj:`` when the name referred can be
resolved to a Python object, and otherwise it becomes simple emphasis.
This can be used as the default role to make links 'smart'.

"""
import sys, os, posixpath, re

from docutils.parsers.rst import directives
from docutils.statemachine import ViewList
from docutils import nodes

import sphinx.addnodes, sphinx.roles
from sphinx.util import patfilter

from docscrape_sphinx import get_doc_object

import warnings
warnings.warn(
    "The numpydoc.autosummary extension can also be found as "
    "sphinx.ext.autosummary in Sphinx >= 0.6, and the version in "
    "Sphinx >= 0.7 is superior to the one in numpydoc. This numpydoc "
    "version of autosummary is no longer maintained.",
    DeprecationWarning, stacklevel=2)

def setup(app):
    app.add_directive('autosummary', autosummary_directive, True, (0, 0, False),
                      toctree=directives.unchanged,
                      nosignatures=directives.flag)
    app.add_role('autolink', autolink_role)
    
    app.add_node(autosummary_toc,
                 html=(autosummary_toc_visit_html, autosummary_toc_depart_noop),
                 latex=(autosummary_toc_visit_latex, autosummary_toc_depart_noop))
    app.connect('doctree-read', process_autosummary_toc)

#------------------------------------------------------------------------------
# autosummary_toc node
#------------------------------------------------------------------------------

class autosummary_toc(nodes.comment):
    pass

def process_autosummary_toc(app, doctree):
    """
    Insert items described in autosummary:: to the TOC tree, but do
    not generate the toctree:: list.

    """
    env = app.builder.env
    crawled = {}
    def crawl_toc(node, depth=1):
        crawled[node] = True
        for j, subnode in enumerate(node):
            try:
                if (isinstance(subnode, autosummary_toc)
                    and isinstance(subnode[0], sphinx.addnodes.toctree)):
                    env.note_toctree(env.docname, subnode[0])
                    continue
            except IndexError:
                continue
            if not isinstance(subnode, nodes.section):
                continue
            if subnode not in crawled:
                crawl_toc(subnode, depth+1)
    crawl_toc(doctree)

def autosummary_toc_visit_html(self, node):
    """Hide autosummary toctree list in HTML output"""
    raise nodes.SkipNode

def autosummary_toc_visit_latex(self, node):
    """Show autosummary toctree (= put the referenced pages here) in Latex"""
    pass

def autosummary_toc_depart_noop(self, node):
    pass

#------------------------------------------------------------------------------
# .. autosummary::
#------------------------------------------------------------------------------

def autosummary_directive(dirname, arguments, options, content, lineno,
                          content_offset, block_text, state, state_machine):
    """
    Pretty table containing short signatures and summaries of functions etc.

    autosummary also generates a (hidden) toctree:: node.

    """

    names = []
    names += [x.strip().split()[0] for x in content
              if x.strip() and re.search(r'^[a-zA-Z_]', x.strip()[0])]

    table, warnings, real_names = get_autosummary(names, state,
                                                  'nosignatures' in options)
    node = table

    env = state.document.settings.env
    suffix = env.config.source_suffix
    all_docnames = env.found_docs.copy()
    dirname = posixpath.dirname(env.docname)

    if 'toctree' in options:
        tree_prefix = options['toctree'].strip()
        docnames = []
        for name in names:
            name = real_names.get(name, name)

            docname = tree_prefix + name
            if docname.endswith(suffix):
                docname = docname[:-len(suffix)]
            docname = posixpath.normpath(posixpath.join(dirname, docname))
            if docname not in env.found_docs:
                warnings.append(state.document.reporter.warning(
                    'toctree references unknown document %r' % docname,
                    line=lineno))
            docnames.append(docname)

        tocnode = sphinx.addnodes.toctree()
        tocnode['includefiles'] = docnames
        tocnode['maxdepth'] = -1
        tocnode['glob'] = None
        tocnode['entries'] = [(None, docname) for docname in docnames]

        tocnode = autosummary_toc('', '', tocnode)
        return warnings + [node] + [tocnode]
    else:
        return warnings + [node]

def get_autosummary(names, state, no_signatures=False):
    """
    Generate a proper table node for autosummary:: directive.

    Parameters
    ----------
    names : list of str
        Names of Python objects to be imported and added to the table.
    document : document
        Docutils document object
    
    """
    document = state.document
    
    real_names = {}
    warnings = []

    prefixes = ['']
    prefixes.insert(0, document.settings.env.currmodule)

    table = nodes.table('')
    group = nodes.tgroup('', cols=2)
    table.append(group)
    group.append(nodes.colspec('', colwidth=10))
    group.append(nodes.colspec('', colwidth=90))
    body = nodes.tbody('')
    group.append(body)

    def append_row(*column_texts):
        row = nodes.row('')
        for text in column_texts:
            node = nodes.paragraph('')
            vl = ViewList()
            vl.append(text, '<autosummary>')
            state.nested_parse(vl, 0, node)
            try:
                if isinstance(node[0], nodes.paragraph):
                    node = node[0]
            except IndexError:
                pass
            row.append(nodes.entry('', node))
        body.append(row)

    for name in names:
        try:
            obj, real_name = import_by_name(name, prefixes=prefixes)
        except ImportError:
            warnings.append(document.reporter.warning(
                'failed to import %s' % name))
            append_row(":obj:`%s`" % name, "")
            continue

        real_names[name] = real_name

        doc = get_doc_object(obj)

        if doc['Summary']:
            title = " ".join(doc['Summary'])
        else:
            title = ""
        
        col1 = u":obj:`%s <%s>`" % (name, real_name)
        if doc['Signature']:
            sig = re.sub('^[^(\[]*', '', doc['Signature'].strip())
            if '=' in sig:
                # abbreviate optional arguments
                sig = re.sub(r', ([a-zA-Z0-9_]+)=', r'[, \1=', sig, count=1)
                sig = re.sub(r'\(([a-zA-Z0-9_]+)=', r'([\1=', sig, count=1)
                sig = re.sub(r'=[^,)]+,', ',', sig)
                sig = re.sub(r'=[^,)]+\)$', '])', sig)
                # shorten long strings
                sig = re.sub(r'(\[.{16,16}[^,]*?),.*?\]\)', r'\1, ...])', sig)
            else:
                sig = re.sub(r'(\(.{16,16}[^,]*?),.*?\)', r'\1, ...)', sig)
            # make signature contain non-breaking spaces
            col1 += u"\\ \u00a0" + unicode(sig).replace(u" ", u"\u00a0")
        col2 = title
        append_row(col1, col2)

    return table, warnings, real_names

def import_by_name(name, prefixes=[None]):
    """
    Import a Python object that has the given name, under one of the prefixes.

    Parameters
    ----------
    name : str
        Name of a Python object, eg. 'numpy.ndarray.view'
    prefixes : list of (str or None), optional
        Prefixes to prepend to the name (None implies no prefix).
        The first prefixed name that results to successful import is used.

    Returns
    -------
    obj
        The imported object
    name
        Name of the imported object (useful if `prefixes` was used)
    
    """
    for prefix in prefixes:
        try:
            if prefix:
                prefixed_name = '.'.join([prefix, name])
            else:
                prefixed_name = name
            return _import_by_name(prefixed_name), prefixed_name
        except ImportError:
            pass
    raise ImportError

def _import_by_name(name):
    """Import a Python object given its full name"""
    try:
        # try first interpret `name` as MODNAME.OBJ
        name_parts = name.split('.')
        try:
            modname = '.'.join(name_parts[:-1])
            __import__(modname)
            return getattr(sys.modules[modname], name_parts[-1])
        except (ImportError, IndexError, AttributeError):
            pass
       
        # ... then as MODNAME, MODNAME.OBJ1, MODNAME.OBJ1.OBJ2, ...
        last_j = 0
        modname = None
        for j in reversed(range(1, len(name_parts)+1)):
            last_j = j
            modname = '.'.join(name_parts[:j])
            try:
                __import__(modname)
            except ImportError:
                continue
            if modname in sys.modules:
                break

        if last_j < len(name_parts):
            obj = sys.modules[modname]
            for obj_name in name_parts[last_j:]:
                obj = getattr(obj, obj_name)
            return obj
        else:
            return sys.modules[modname]
    except (ValueError, ImportError, AttributeError, KeyError), e:
        raise ImportError(e)

#------------------------------------------------------------------------------
# :autolink: (smart default role)
#------------------------------------------------------------------------------

def autolink_role(typ, rawtext, etext, lineno, inliner,
                  options={}, content=[]):
    """
    Smart linking role.

    Expands to ":obj:`text`" if `text` is an object that can be imported;
    otherwise expands to "*text*".
    """
    r = sphinx.roles.xfileref_role('obj', rawtext, etext, lineno, inliner,
                                   options, content)
    pnode = r[0][0]

    prefixes = [None]
    #prefixes.insert(0, inliner.document.settings.env.currmodule)
    try:
        obj, name = import_by_name(pnode['reftarget'], prefixes)
    except ImportError:
        content = pnode[0]
        r[0][0] = nodes.emphasis(rawtext, content[0].astext(),
                                 classes=content['classes'])
    return r
