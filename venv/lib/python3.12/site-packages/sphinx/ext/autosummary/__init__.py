"""Extension that adds an autosummary:: directive.

The directive can be used to generate function/method/attribute/etc. summary
lists, similar to those output eg. by Epydoc and other API doc generation tools.

An :autolink: role is also provided.

autosummary directive
---------------------

The autosummary directive has the form::

    .. autosummary::
       :signatures: none
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

from __future__ import annotations

import functools
import inspect
import operator
import posixpath
import re
import sys
from inspect import Parameter
from types import ModuleType
from typing import TYPE_CHECKING, cast

from docutils import nodes
from docutils.parsers.rst import directives
from docutils.parsers.rst.states import RSTStateMachine, state_classes
from docutils.statemachine import StringList

import sphinx
from sphinx import addnodes
from sphinx.errors import PycodeError
from sphinx.ext.autodoc._directive_options import _AutoDocumenterOptions
from sphinx.ext.autodoc._dynamic._importer import _import_module
from sphinx.ext.autodoc._dynamic._loader import _load_object_by_name
from sphinx.ext.autodoc._dynamic._member_finder import _best_object_type_for_member
from sphinx.ext.autodoc._dynamic._mock import mock
from sphinx.ext.autodoc._sentinels import INSTANCE_ATTR
from sphinx.ext.autodoc._shared import _AutodocAttrGetter, _AutodocConfig
from sphinx.locale import __
from sphinx.pycode import ModuleAnalyzer
from sphinx.util import logging, rst
from sphinx.util.docutils import (
    NullReporter,
    SphinxDirective,
    SphinxRole,
    new_document,
    switch_source_input,
)
from sphinx.util.inspect import getmro, signature_from_str
from sphinx.util.matching import Matcher
from sphinx.util.parsing import nested_parse_to_nodes

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Any, ClassVar

    from docutils.nodes import Node, system_message

    from sphinx.application import Sphinx
    from sphinx.environment import BuildEnvironment
    from sphinx.ext.autodoc._property_types import _AutodocObjType
    from sphinx.util.typing import ExtensionMetadata, OptionSpec
    from sphinx.writers.html5 import HTML5Translator

logger = logging.getLogger(__name__)


periods_re = re.compile(r'\.(?:\s+)')
literal_re = re.compile(r'::\s*$')

WELL_KNOWN_ABBREVIATIONS = ('et al.', 'e.g.', 'i.e.', 'vs.')


# -- autosummary_toc node ------------------------------------------------------


class autosummary_toc(nodes.comment):
    pass


def autosummary_toc_visit_html(self: nodes.NodeVisitor, node: autosummary_toc) -> None:
    """Hide autosummary toctree list in HTML output."""
    raise nodes.SkipNode


def autosummary_noop(self: nodes.NodeVisitor, node: Node) -> None:
    pass


# -- autosummary_table node ----------------------------------------------------


class autosummary_table(nodes.comment):
    pass


def autosummary_table_visit_html(
    self: HTML5Translator, node: autosummary_table
) -> None:
    """Make the first column of the table non-breaking."""
    try:
        table = cast('nodes.table', node[0])
        tgroup = cast('nodes.tgroup', table[0])
        tbody = cast('nodes.tbody', tgroup[-1])
        rows = cast('list[nodes.row]', tbody)
        for row in rows:
            col1_entry = cast('nodes.entry', row[0])
            par = cast('nodes.paragraph', col1_entry[0])
            for j, subnode in enumerate(list(par)):
                if isinstance(subnode, nodes.Text):
                    new_text = subnode.astext().replace(' ', '\u00a0')
                    par[j] = nodes.Text(new_text)
    except IndexError:
        pass


# -- autodoc integration -------------------------------------------------------


def _get_documenter(obj: Any, parent: Any) -> _AutodocObjType:
    """Get the best object type suitable for documenting the given object.

    *obj* is the Python object to be documented, and *parent* is another
    Python object (e.g. a module or a class) to which *obj* belongs.
    """
    if inspect.ismodule(obj):
        return 'module'

    if parent is None or inspect.ismodule(parent):
        parent_obj_type = 'module'
    else:
        parent_opt = _best_object_type_for_member(
            member=parent,
            member_name='',
            is_attr=False,
            parent_obj_type='module',
            parent_props=None,
        )
        parent_obj_type = parent_opt if parent_opt is not None else 'data'

    if obj_type := _best_object_type_for_member(
        member=obj,
        member_name='',
        is_attr=False,
        parent_obj_type=parent_obj_type,
        parent_props=None,
    ):
        return obj_type
    return 'data'


# -- .. autosummary:: ----------------------------------------------------------


class Autosummary(SphinxDirective):
    """Pretty table containing short signatures and summaries of functions etc.

    autosummary can also optionally generate a hidden toctree:: node.
    """

    required_arguments = 0
    optional_arguments = 0
    final_argument_whitespace = False
    has_content = True
    option_spec: ClassVar[OptionSpec] = {
        'caption': directives.unchanged_required,
        'class': directives.class_option,
        'toctree': directives.unchanged,
        'nosignatures': directives.flag,
        'recursive': directives.flag,
        'signatures': directives.unchanged,
        'template': directives.unchanged,
    }

    def run(self) -> list[Node]:
        names = [
            x.strip().split()[0]
            for x in self.content
            if x.strip() and re.search(r'^[~a-zA-Z_]', x.strip()[0])
        ]
        items = self.get_items(names)
        nodes = self.get_table(items)

        if 'toctree' in self.options:
            dirname = posixpath.dirname(self.env.current_document.docname)

            tree_prefix = self.options['toctree'].strip()
            docnames = []
            excluded = Matcher(self.config.exclude_patterns)
            filename_map = self.config.autosummary_filename_map
            for _name, _sig, _summary, real_name in items:
                real_name = filename_map.get(real_name, real_name)
                docname = posixpath.join(tree_prefix, real_name)
                docname = posixpath.normpath(posixpath.join(dirname, docname))
                if docname not in self.env.found_docs:
                    if excluded(str(self.env.doc2path(docname, False))):
                        msg = __(
                            'autosummary references excluded document %r. Ignored.'
                        )
                    else:
                        msg = __(
                            'autosummary: stub file not found %r. '
                            'Check your autosummary_generate setting.'
                        )

                    logger.warning(msg, real_name, location=self.get_location())
                    continue

                docnames.append(docname)

            if docnames:
                tocnode = addnodes.toctree()
                tocnode['includefiles'] = docnames
                tocnode['entries'] = [(None, docn) for docn in docnames]
                tocnode['maxdepth'] = -1
                tocnode['glob'] = None
                tocnode['caption'] = self.options.get('caption')

                nodes.append(autosummary_toc('', '', tocnode))

        if 'toctree' not in self.options and 'caption' in self.options:
            logger.warning(
                __('A captioned autosummary requires :toctree: option. ignored.'),
                location=nodes[-1],
            )

        return nodes

    def import_by_name(
        self, name: str, prefixes: list[str | None]
    ) -> tuple[str, Any, Any, str]:
        with mock(self.config.autosummary_mock_imports):
            try:
                return import_by_name(name, prefixes)
            except ImportExceptionGroup as exc:
                # check existence of instance attribute
                try:
                    return import_ivar_by_name(name, prefixes)
                except ImportError as exc2:
                    if exc2.__cause__:
                        errors: list[BaseException] = [*exc.exceptions, exc2.__cause__]
                    else:
                        errors = [*exc.exceptions, exc2]

                    raise ImportExceptionGroup(exc.args[0], errors) from None

    def get_items(self, names: list[str]) -> list[tuple[str, str | None, str, str]]:
        """Try to import the given names, and return a list of
        ``[(name, signature, summary_string, real_name), ...]``.

        signature is already formatted and is None if :nosignatures: option was given.
        """
        prefixes = get_import_prefixes_from_env(self.env)

        items: list[tuple[str, str | None, str, str]] = []

        signatures_option = self.options.get('signatures')
        if signatures_option is None:
            signatures_option = 'none' if 'nosignatures' in self.options else 'long'
        if signatures_option not in {'none', 'short', 'long'}:
            msg = (
                'Invalid value for autosummary :signatures: option: '
                f"{signatures_option!r}. Valid values are 'none', 'short', 'long'"
            )
            raise ValueError(msg)

        document_settings = self.state.document.settings
        env = self.env
        config = _AutodocConfig.from_config(env.config)
        current_document = env.current_document
        events = env.events
        get_attr = _AutodocAttrGetter(env._registry.autodoc_attrgetters)
        opts = _AutoDocumenterOptions()
        ref_context = env.ref_context
        reread_always = env.reread_always

        max_item_chars = 50

        for name in names:
            display_name = name
            if name.startswith('~'):
                name = name[1:]
                display_name = name.split('.')[-1]

            try:
                real_name, obj, parent, modname = self.import_by_name(
                    name, prefixes=prefixes
                )
            except ImportExceptionGroup as exc:
                errors = list({f'* {type(e).__name__}: {e}' for e in exc.exceptions})
                logger.warning(
                    __('autosummary: failed to import %s.\nPossible hints:\n%s'),
                    name,
                    '\n'.join(errors),
                    location=self.get_location(),
                )
                continue

            obj_type = _get_documenter(obj, parent)
            if isinstance(obj, ModuleType):
                full_name = real_name
            else:
                # give explicitly separated module name, so that members
                # of inner classes can be documented
                full_name = f'{modname}::{real_name[len(modname) + 1 :]}'
            # NB. using full_name here is important, since Documenters
            #     handle module prefixes slightly differently
            props = _load_object_by_name(
                name=full_name,
                objtype=obj_type,
                current_document=current_document,
                config=config,
                events=events,
                get_attr=get_attr,
                options=opts,
                ref_context=ref_context,
                reread_always=reread_always,
            )
            if props is None:
                logger.warning(
                    __('failed to import object %s'),
                    real_name,
                    location=self.get_location(),
                )
                items.append((display_name, '', '', real_name))
                continue

            # -- Grab the signature

            if signatures_option == 'none':
                sig = None
            elif not props.signatures:
                sig = ''
            elif signatures_option == 'short':
                sig = '()' if props.signatures == ('()',) else '(…)'
            else:  # signatures_option == 'long'
                max_chars = max(10, max_item_chars - len(display_name))
                sig = mangle_signature('\n'.join(props.signatures), max_chars=max_chars)

            # -- Grab the summary

            # get content from docstrings or attribute documentation
            summary = extract_summary(props.docstring_lines, document_settings)

            items.append((display_name, sig, summary, real_name))

        return items

    def get_table(self, items: list[tuple[str, str | None, str, str]]) -> list[Node]:
        """Generate a proper list of table nodes for autosummary:: directive.

        *items* is a list produced by :meth:`get_items`.
        """
        table_spec = addnodes.tabular_col_spec()
        table_spec['spec'] = r'\X{1}{2}\X{1}{2}'

        table = autosummary_table('')
        real_table = nodes.table(
            '', classes=['autosummary', 'longtable', *self.options.get('class', ())]
        )
        table.append(real_table)
        group = nodes.tgroup('', cols=2)
        real_table.append(group)
        group.append(nodes.colspec('', colwidth=10))
        group.append(nodes.colspec('', colwidth=90))
        body = nodes.tbody('')
        group.append(body)

        def append_row(*column_texts: str) -> None:
            row = nodes.row('')
            source, line = self.state_machine.get_source_and_line()
            for text in column_texts:
                vl = StringList([text], f'{source}:{line}:<autosummary>')
                with switch_source_input(self.state, vl):
                    col_nodes = nested_parse_to_nodes(
                        self.state, vl, allow_section_headings=False
                    )
                    if col_nodes and isinstance(col_nodes[0], nodes.paragraph):
                        node = col_nodes[0]
                    else:
                        node = nodes.paragraph('')
                    row.append(nodes.entry('', node))
            body.append(row)

        for name, sig, summary, real_name in items:
            qualifier = 'obj'
            if sig is None:
                col1 = f':py:{qualifier}:`{name} <{real_name}>`'
            else:
                col1 = f':py:{qualifier}:`{name} <{real_name}>`\\ {rst.escape(sig)}'

            col2 = summary
            append_row(col1, col2)

        return [table_spec, table]


def strip_arg_typehint(s: str) -> str:
    """Strip a type hint from argument definition."""
    return s.partition(':')[0].strip()


def _cleanup_signature(s: str) -> str:
    """Clean up signature using inspect.signautre() for mangle_signature()"""
    try:
        sig = signature_from_str(s)
        parameters = list(sig.parameters.values())
        for i, param in enumerate(parameters):
            if param.annotation is not Parameter.empty:
                # Remove typehints
                param = param.replace(annotation=Parameter.empty)
            if param.default is not Parameter.empty:
                # Replace default value by "None"
                param = param.replace(default=None)
            parameters[i] = param
        sig = sig.replace(parameters=parameters, return_annotation=Parameter.empty)
        return str(sig)
    except Exception:
        # Return the original signature string if failed to clean (ex. parsing error)
        return s


def mangle_signature(sig: str, max_chars: int = 30) -> str:
    """Reformat a function signature to a more compact form."""
    s = _cleanup_signature(sig)

    # Strip return type annotation
    s = re.sub(r'\)\s*->\s.*$', ')', s)

    # Remove parenthesis
    s = re.sub(r'^\((.*)\)$', r'\1', s).strip()

    # Strip literals (which can contain things that confuse the code below)
    s = re.sub(r'\\\\', '', s)  # escaped backslash (maybe inside string)
    s = re.sub(r"\\'", '', s)  # escaped single quote
    s = re.sub(r'\\"', '', s)  # escaped double quote
    s = re.sub(r"'[^']*'", '', s)  # string literal (w/ single quote)
    s = re.sub(r'"[^"]*"', '', s)  # string literal (w/ double quote)

    # Strip complex objects (maybe default value of arguments)
    while re.search(
        r'\([^)]*\)', s
    ):  # contents of parenthesis (ex. NamedTuple(attr=...))
        s = re.sub(r'\([^)]*\)', '', s)
    while re.search(r'<[^>]*>', s):  # contents of angle brackets (ex. <object>)
        s = re.sub(r'<[^>]*>', '', s)
    while re.search(r'{[^}]*}', s):  # contents of curly brackets (ex. dict)
        s = re.sub(r'{[^}]*}', '', s)

    # Parse the signature to arguments + options
    args: list[str] = []
    opts: list[str] = []

    opt_re = re.compile(r'^(.*, |)([a-zA-Z0-9_*]+)\s*=\s*')
    while s:
        m = opt_re.search(s)
        if not m:
            # The rest are arguments
            args = s.split(', ')
            break

        opts.insert(0, m.group(2))
        s = m.group(1)[:-2]

    # Strip typehints
    for i, arg in enumerate(args):
        args[i] = strip_arg_typehint(arg)

    for i, opt in enumerate(opts):
        opts[i] = strip_arg_typehint(opt)

    # Produce a more compact signature
    sig = limited_join(', ', args, max_chars=max_chars - 2)
    if opts:
        if not sig:
            sig = '[%s]' % limited_join(', ', opts, max_chars=max_chars - 4)
        elif len(sig) < max_chars - 4 - 2 - 3:
            sig += '[, %s]' % limited_join(
                ', ', opts, max_chars=max_chars - len(sig) - 4 - 2
            )

    return '(%s)' % sig


def extract_summary(doc: Sequence[str], settings: Any) -> str:
    """Extract summary from docstring."""
    # Find the first stanza (heading, sentence, paragraph, etc.).
    # If there's a blank line, then we can assume that the stanza has ended,
    # so anything after shouldn't be part of the summary.
    first_stanza = []
    content_started = False
    for line in doc:
        is_blank_line = not line or line.isspace()
        if not content_started:
            # Skip any blank lines at the start
            if is_blank_line:
                continue
            content_started = True
        if content_started:
            if is_blank_line:
                break
            first_stanza.append(line)

    if not first_stanza:
        return ''

    # parse the docstring
    node = _parse_summary(first_stanza, settings)
    if isinstance(node[0], nodes.section):
        # document starts with a section heading, so use that.
        summary = node[0].astext().strip()
    elif not isinstance(node[0], nodes.paragraph):
        # document starts with non-paragraph: pick up the first line
        summary = first_stanza[0].strip()
    else:
        # Try to find the "first sentence", which may span multiple lines
        sentences = periods_re.split(' '.join(first_stanza))
        if len(sentences) == 1:
            summary = sentences[0].strip()
        else:
            summary = ''
            for i in range(len(sentences)):
                summary = '. '.join(sentences[: i + 1]).rstrip('.') + '.'
                node[:] = []
                node = _parse_summary(first_stanza, settings)
                if summary.endswith(WELL_KNOWN_ABBREVIATIONS):
                    pass
                elif not any(node.findall(nodes.system_message)):
                    # considered as that splitting by period does not break inline markups
                    break

    # strip literal notation mark ``::`` from tail of summary
    summary = literal_re.sub('.', summary)

    return summary


def _parse_summary(doc: Sequence[str], settings: Any) -> nodes.document:
    state_machine = RSTStateMachine(state_classes, 'Body')
    node = new_document('', settings)
    node.reporter = NullReporter()
    state_machine.run(doc, node)

    return node


def limited_join(
    sep: str, items: list[str], max_chars: int = 30, overflow_marker: str = '...'
) -> str:
    """Join a number of strings into one, limiting the length to *max_chars*.

    If the string overflows this limit, replace the last fitting item by
    *overflow_marker*.

    Returns: joined_string
    """
    full_str = sep.join(items)
    if len(full_str) < max_chars:
        return full_str

    n_chars = 0
    n_items = 0
    for item in items:
        n_chars += len(item) + len(sep)
        if n_chars < max_chars - len(overflow_marker):
            n_items += 1
        else:
            break

    return sep.join([*list(items[:n_items]), overflow_marker])


# -- Importing items -----------------------------------------------------------


class ImportExceptionGroup(BaseExceptionGroup):
    """Exceptions raised during importing the target objects.

    It contains an error message and a list of exceptions as its arguments.
    """


def get_import_prefixes_from_env(env: BuildEnvironment) -> list[str | None]:
    """Obtain current Python import prefixes (for `import_by_name`)
    from ``document.env``
    """
    prefixes: list[str | None] = [None]

    currmodule = env.ref_context.get('py:module')
    if currmodule:
        prefixes.insert(0, currmodule)

    currclass = env.ref_context.get('py:class')
    if currclass:
        if currmodule:
            prefixes.insert(0, f'{currmodule}.{currclass}')
        else:
            prefixes.insert(0, currclass)

    return prefixes


def import_by_name(
    name: str, prefixes: Sequence[str | None] = (None,)
) -> tuple[str, Any, Any, str]:
    """Import a Python object that has the given *name*, under one of the
    *prefixes*.  The first name that succeeds is used.
    """
    tried = []
    errors: list[ImportExceptionGroup] = []
    for prefix in prefixes:
        if prefix is not None and name.startswith(f'{prefix}.'):
            # Catch and avoid module cycles (e.g., sphinx.ext.sphinx.ext...)
            msg = __(
                'Summarised items should not include the current module. '
                'Replace %r with %r.'
            )
            logger.warning(
                msg,
                name,
                name.removeprefix(f'{prefix}.'),
                type='autosummary',
                subtype='import_cycle',
            )
            continue
        try:
            if prefix:
                prefixed_name = f'{prefix}.{name}'
            else:
                prefixed_name = name
            obj, parent, modname = _import_by_name(
                prefixed_name, grouped_exception=True
            )
            return prefixed_name, obj, parent, modname
        except ImportError:
            tried.append(prefixed_name)
        except ImportExceptionGroup as exc:
            tried.append(prefixed_name)
            errors.append(exc)

    exceptions: list[BaseException] = functools.reduce(
        operator.iadd, (e.exceptions for e in errors), []
    )
    msg = f'could not import {" or ".join(tried)}'
    raise ImportExceptionGroup(msg, exceptions)


def _import_by_name(name: str, grouped_exception: bool = True) -> tuple[Any, Any, str]:
    """Import a Python object given its full name."""
    errors: list[BaseException] = []

    try:
        name_parts = name.split('.')

        # try first interpret `name` as MODNAME.OBJ
        modname = '.'.join(name_parts[:-1])
        if modname:
            try:
                mod = _import_module(modname)
                return getattr(mod, name_parts[-1]), mod, modname
            except (ImportError, IndexError, AttributeError) as exc:
                errors.append(exc.__cause__ or exc)

        # ... then as MODNAME, MODNAME.OBJ1, MODNAME.OBJ1.OBJ2, ...
        last_j = 0
        modname = ''
        for j in reversed(range(1, len(name_parts) + 1)):
            last_j = j
            modname = '.'.join(name_parts[:j])
            try:
                _import_module(modname)
            except ImportError as exc:
                errors.append(exc.__cause__ or exc)

            if modname in sys.modules:
                break

        if last_j < len(name_parts):
            parent = None
            obj = sys.modules[modname]
            for obj_name in name_parts[last_j:]:
                parent = obj
                obj = getattr(obj, obj_name)
            return obj, parent, modname
        else:
            return sys.modules[modname], None, modname
    except (ValueError, ImportError, AttributeError, KeyError) as exc:
        errors.append(exc)
        if grouped_exception:
            raise ImportExceptionGroup('', errors) from None  # NoQA: EM101
        else:
            raise ImportError(*exc.args) from exc


def import_ivar_by_name(
    name: str, prefixes: Sequence[str | None] = (None,), grouped_exception: bool = True
) -> tuple[str, Any, Any, str]:
    """Import an instance variable that has the given *name*, under one of the
    *prefixes*.  The first name that succeeds is used.
    """
    try:
        name, attr = name.rsplit('.', 1)
        real_name, obj, _parent, modname = import_by_name(name, prefixes)

        # Get ancestors of the object (class.__mro__ includes the class itself as
        # the first entry)
        candidate_objects = getmro(obj)
        if len(candidate_objects) == 0:
            candidate_objects = (obj,)

        for candidate_obj in candidate_objects:
            analyzer = ModuleAnalyzer.for_module(
                getattr(candidate_obj, '__module__', modname)
            )
            analyzer.analyze()
            # check for presence in `annotations` to include dataclass attributes
            found_attrs = set()
            found_attrs |= {attr for (qualname, attr) in analyzer.attr_docs}
            found_attrs |= {attr for (qualname, attr) in analyzer.annotations}
            if attr in found_attrs:
                return f'{real_name}.{attr}', INSTANCE_ATTR, obj, modname
    except (ImportError, ValueError, PycodeError) as exc:
        raise ImportError from exc
    except ImportExceptionGroup:
        raise  # pass through it as is

    raise ImportError


# -- :autolink: (smart default role) -------------------------------------------


class AutoLink(SphinxRole):
    """Smart linking role.

    Expands to ':obj:`text`' if `text` is an object that can be imported;
    otherwise expands to '*text*'.
    """

    def run(self) -> tuple[list[Node], list[system_message]]:
        pyobj_role = self.env.domains.python_domain.role('obj')
        assert pyobj_role is not None
        objects, errors = pyobj_role(
            'obj',
            self.rawtext,
            self.text,
            self.lineno,
            self.inliner,
            self.options,
            self.content,
        )
        if errors:
            return objects, errors

        assert len(objects) == 1
        pending_xref = cast('addnodes.pending_xref', objects[0])
        try:
            # try to import object by name
            prefixes = get_import_prefixes_from_env(self.env)
            name = pending_xref['reftarget']
            prefixes = [
                prefix
                for prefix in prefixes
                if prefix is None
                or not (name.startswith(f'{prefix}.') or name == prefix)
            ]
            import_by_name(name, prefixes)
        except ImportExceptionGroup:
            literal = cast('nodes.literal', pending_xref[0])
            objects[0] = nodes.emphasis(
                self.rawtext, literal.astext(), classes=literal['classes']
            )

        return objects, errors


def get_rst_suffix(app: Sphinx) -> str | None:
    def get_supported_format(suffix: str) -> tuple[str, ...]:
        parser_class = app.registry.get_source_parsers().get(suffix.removeprefix('.'))
        if parser_class is None:
            return ('restructuredtext',)
        return parser_class.supported

    suffix = None
    for suffix in app.config.source_suffix:
        if 'restructuredtext' in get_supported_format(suffix):
            return suffix

    return None


def process_generate_options(app: Sphinx) -> None:
    genfiles = app.config.autosummary_generate

    if genfiles is True:
        env = app.env
        genfiles = [
            str(env.doc2path(x, base=False))
            for x in env.found_docs
            if env.doc2path(x).is_file()
        ]
    elif genfiles is False:
        pass
    else:
        ext = list(app.config.source_suffix)
        genfiles = [
            genfile + (ext[0] if not genfile.endswith(tuple(ext)) else '')
            for genfile in genfiles
        ]

        for entry in genfiles[:]:
            if not (app.srcdir / entry).is_file():
                logger.warning(__('autosummary_generate: file not found: %s'), entry)
                genfiles.remove(entry)

    if not genfiles:
        return

    suffix = get_rst_suffix(app)
    if suffix is None:
        logger.warning(
            __(
                'autosummary generates .rst files internally. '
                'But your source_suffix does not contain .rst. Skipped.'
            )
        )
        return

    from sphinx.ext.autosummary.generate import generate_autosummary_docs

    imported_members = app.config.autosummary_imported_members
    with mock(app.config.autosummary_mock_imports):
        generate_autosummary_docs(
            genfiles,
            suffix=suffix,
            base_path=app.srcdir,
            app=app,
            imported_members=imported_members,
            overwrite=app.config.autosummary_generate_overwrite,
            encoding=app.config.source_encoding,
        )


def setup(app: Sphinx) -> ExtensionMetadata:
    # I need autodoc
    app.setup_extension('sphinx.ext.autodoc')
    app.add_node(
        autosummary_toc,
        html=(autosummary_toc_visit_html, autosummary_noop),
        latex=(autosummary_noop, autosummary_noop),
        text=(autosummary_noop, autosummary_noop),
        man=(autosummary_noop, autosummary_noop),
        texinfo=(autosummary_noop, autosummary_noop),
    )
    app.add_node(
        autosummary_table,
        html=(autosummary_table_visit_html, autosummary_noop),
        latex=(autosummary_noop, autosummary_noop),
        text=(autosummary_noop, autosummary_noop),
        man=(autosummary_noop, autosummary_noop),
        texinfo=(autosummary_noop, autosummary_noop),
    )
    app.add_directive('autosummary', Autosummary)
    app.add_role('autolink', AutoLink())
    app.connect('builder-inited', process_generate_options)
    app.add_config_value('autosummary_context', {}, 'env', types=frozenset({dict}))
    app.add_config_value(
        'autosummary_filename_map', {}, 'html', types=frozenset({dict})
    )
    app.add_config_value(
        'autosummary_generate', True, 'env', types=frozenset({bool, list})
    )
    app.add_config_value(
        'autosummary_generate_overwrite', True, '', types=frozenset({bool})
    )
    app.add_config_value(
        'autosummary_mock_imports',
        lambda config: config.autodoc_mock_imports,
        'env',
        types=frozenset({list, tuple}),
    )
    app.add_config_value(
        'autosummary_imported_members', False, '', types=frozenset({bool})
    )
    app.add_config_value(
        'autosummary_ignore_module_all', True, 'env', types=frozenset({bool})
    )

    return {
        'version': sphinx.__display_version__,
        'parallel_read_safe': True,
    }
