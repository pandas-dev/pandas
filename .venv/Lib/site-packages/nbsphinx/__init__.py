"""Jupyter Notebook Tools for Sphinx.

https://nbsphinx.readthedocs.io/

"""
__version__ = '0.9.7'

import collections.abc
import copy
import html
from itertools import chain
import json
import os
import re
import subprocess
import sys
from urllib.parse import unquote
import uuid

import docutils
from docutils.parsers import rst
import jinja2
import nbconvert
import nbformat
import sphinx
import sphinx.directives
import sphinx.directives.other
import sphinx.environment
import sphinx.errors
import sphinx.transforms.post_transforms.images
from sphinx.util.matching import patmatch
try:
    from sphinx.util.display import status_iterator
except ImportError:
    # This will be removed in Sphinx 8:
    from sphinx.util import status_iterator
import traitlets


if sys.version_info >= (3, 8) and sys.platform == 'win32':
    # See: https://github.com/jupyter/jupyter_client/issues/583
    import asyncio

    asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())

_ipynbversion = 4

logger = sphinx.util.logging.getLogger(__name__)

_BROKEN_THUMBNAIL = None

# See nbconvert/exporters/html.py:
DISPLAY_DATA_PRIORITY_HTML = (
    'application/vnd.jupyter.widget-state+json',
    'application/vnd.jupyter.widget-view+json',
    'application/javascript',
    'text/html',
    'text/markdown',
    'image/svg+xml',
    'text/latex',
    'image/png',
    'image/jpeg',
    'text/plain',
)
# See nbconvert/exporters/latex.py:
DISPLAY_DATA_PRIORITY_LATEX = (
    'text/latex',
    'application/pdf',
    'image/png',
    'image/jpeg',
    'image/svg+xml',
    'text/markdown',
    'text/plain',
)

THUMBNAIL_MIME_TYPES = {
    'image/svg+xml': '.svg',
    'image/png': '.png',
    'image/jpeg': '.jpg',
}

# The default rst template name is changing in nbconvert 6, so we substitute
# it in to the *extends* directive.
RST_TEMPLATE = """
{% extends '__RST_DEFAULT_TEMPLATE__' %}

{% macro insert_empty_lines(text) %}
{%- set before, after = text | get_empty_lines %}
{%- if before %}
    :empty-lines-before: {{ before }}
{%- endif %}
{%- if after %}
    :empty-lines-after: {{ after }}
{%- endif %}
{%- endmacro %}


{% block any_cell %}
{%- if cell.metadata.nbsphinx != 'hidden' %}
{{ super() }}
..
{# Empty comment to make sure the preceding directive (if any) is closed #}
{% endif %}
{%- endblock any_cell %}


{% block input -%}
.. nbinput:: {% if cell.metadata.magics_language -%}
{{ cell.metadata.magics_language }}
{%- elif nb.metadata.language_info -%}
{{ nb.metadata.language_info.pygments_lexer or nb.metadata.language_info.name }}
{%- else -%}
{{ resources.codecell_lexer }}
{%- endif -%}
{{ insert_empty_lines(cell.source) }}
{%- if cell.execution_count %}
    :execution-count: {{ cell.execution_count }}
{%- endif %}
{%- if not cell.outputs %}
    :no-output:
{%- endif %}

{{ cell.source.strip('\n') | indent }}
{% endblock input %}


{% macro insert_nboutput(datatype, output, cell) -%}
.. nboutput::
{%- if output.output_type == 'execute_result' and cell.execution_count %}
    :execution-count: {{ cell.execution_count }}
{%- endif %}
{%- if output != cell.outputs[-1] %}
    :more-to-come:
{%- endif %}
{%- if output.name == 'stderr' %}
    :class: stderr
{%- endif %}
{%- if datatype != 'text/plain' %}
    :fancy:
{%- endif %}
{%- if datatype == 'text/plain' %}

    .. rst-class:: highlight

    .. raw:: html

        <pre>
{{ output.data[datatype] | ansi2html | indent | indent }}
        </pre>

    .. raw:: latex

        \\begin{sphinxVerbatim}[commandchars=\\\\\\{\\}]
{{ output.data[datatype] | escape_latex | ansi2latex | indent | indent }}
        \\end{sphinxVerbatim}

{# NB: The "raw" directive doesn't work with empty content #}
{%- if output.data[datatype].strip() %}
    .. raw:: text

{{ output.data[datatype] | indent | indent }}
{%- endif %}

{%- elif datatype in ['image/svg+xml', 'image/png', 'image/jpeg', 'application/pdf'] %}

    .. image:: {{ output.metadata.filenames[datatype] | posix_path }}
{%- if datatype in output.metadata %}
        :class: no-scaled-link
{%- set width = output.metadata[datatype].width %}
{%- if width %}
        :width: {{ width }}
{%- endif %}
{%- set height = output.metadata[datatype].height %}
{%- if height %}
        :height: {{ height }}
{% endif %}
{% endif %}
{%- elif datatype in ['text/markdown'] %}

{{ output.data['text/markdown'] | markdown2rst | indent }}
{%- elif datatype in ['text/latex'] %}

    .. math::
        :nowrap:

{{ output.data['text/latex'] | indent | indent }}
{%- elif datatype == 'text/html' %}
    :class: rendered_html

    .. raw:: html

{{ (output.data['text/html'] or '<!-- empty output -->') | indent | indent }}
{%- elif datatype == 'application/javascript' %}

    .. raw:: html

        <div class="output_javascript"></div>
        <script type="text/javascript">
        var element = document.currentScript.previousSibling.previousSibling;
{{ output.data['application/javascript'] | indent | indent }}
        </script>
{%- elif datatype.startswith('application/vnd.jupyter') and datatype.endswith('+json') %}

    .. raw:: html

        <script type="{{ datatype }}">{{ output.data[datatype] | json_dumps }}</script>
{%- elif datatype == '' %}
{# Empty output data #}

    ..
{% else %}

    .. nbwarning:: Data type cannot be displayed: {{ datatype }}
{%- endif %}
{% endmacro %}


{% block nboutput -%}
..
{# Empty comment to make sure the preceding directive (if any) is closed #}
{%- set html_datatype, latex_datatype = output | get_output_type %}
{%- if html_datatype == latex_datatype %}
{{ insert_nboutput(html_datatype, output, cell) }}
{% else %}
.. only:: html

{{ insert_nboutput(html_datatype, output, cell) | indent }}
.. only:: latex

{{ insert_nboutput(latex_datatype, output, cell) | indent }}
{% endif %}
{% endblock nboutput %}


{% block execute_result %}{{ self.nboutput() }}{% endblock execute_result %}
{% block display_data %}{{ self.nboutput() }}{% endblock display_data %}
{% block stream %}{{ self.nboutput() }}{% endblock stream %}
{% block error %}{{ self.nboutput() }}{% endblock error %}


{% block markdowncell %}
{%- if 'nbsphinx-gallery' in cell.metadata
    or 'nbsphinx-gallery' in cell.metadata.get('tags', [])
    or 'nbsphinx-toctree' in cell.metadata
    or 'nbsphinx-toctree' in cell.metadata.get('tags', [])
    or 'nbsphinx-link-gallery' in cell.metadata
    or 'nbsphinx-link-gallery' in cell.metadata.get('tags', []) %}
{{ cell | extract_gallery_or_toctree }}
{%- else %}
{{ cell | save_attachments or super() | replace_attachments }}
{% endif %}
{% endblock markdowncell %}


{% block rawcell %}
{%- set raw_mimetype = cell.metadata.get('raw_mimetype', '').lower() %}
{%- if raw_mimetype == '' %}
.. raw:: html

{{ (cell.source or '<!-- empty raw cell -->') | indent }}

.. raw:: latex

{{ (cell.source or '% empty raw cell') | indent }}
{%- elif raw_mimetype == 'text/html' %}
.. raw:: html

{{ (cell.source or '<!-- empty raw cell -->') | indent }}
{%- elif raw_mimetype == 'text/latex' %}
.. raw:: latex

{{ (cell.source or '% empty raw cell') | indent }}
{%- elif raw_mimetype == 'text/markdown' %}
{{ cell.source | markdown2rst }}
{%- elif raw_mimetype == 'text/restructuredtext' %}
{{ cell.source }}
{% endif %}
{% endblock rawcell %}


{% block footer %}

{% if 'application/vnd.jupyter.widget-state+json' in nb.metadata.get('widgets', {})%}

.. raw:: html

    <script type="application/vnd.jupyter.widget-state+json">
    {{ nb.metadata.widgets['application/vnd.jupyter.widget-state+json'] | json_dumps }}
    </script>
{% endif %}
{{ super() }}
{% endblock footer %}
""".replace('__RST_DEFAULT_TEMPLATE__', nbconvert.RSTExporter().template_file)


class Exporter(nbconvert.RSTExporter):
    """Convert Jupyter notebooks to reStructuredText.

    Uses nbconvert to convert Jupyter notebooks to a reStructuredText
    string with custom reST directives for input and output cells.

    Notebooks without output cells are automatically executed before
    conversion.

    """

    def __init__(self, execute='auto', kernel_name='', execute_arguments=[],
                 allow_errors=False, timeout=None, codecell_lexer='none'):
        """Initialize the Exporter."""

        # NB: The following stateful Jinja filters are a hack until
        # template-based processing is dropped
        # (https://github.com/spatialaudio/nbsphinx/issues/36) or someone
        # comes up with a better idea.

        # NB: This instance-local state makes the methods non-reentrant!
        attachment_storage = []

        def save_attachments(cell):
            for filename, bundle in cell.get('attachments', {}).items():
                attachment_storage.append((filename, bundle))

        def replace_attachments(text):
            for filename, bundle in attachment_storage:
                # For now, this works only if there is a single MIME bundle
                (mime_type, data), = bundle.items()
                text = re.sub(
                    r'^(\s*\.\. ((\|[^|]*\| )?image|figure)::) attachment:{0}$'
                        .format(filename),
                    r'\1 data:{0};base64,{1}'.format(mime_type, data),
                    text, flags=re.MULTILINE)
            del attachment_storage[:]
            return text

        self._execute = execute
        self._kernel_name = kernel_name
        self._execute_arguments = execute_arguments
        self._allow_errors = allow_errors
        self._timeout = timeout
        self._codecell_lexer = codecell_lexer
        loader = jinja2.DictLoader({'nbsphinx-rst.tpl': RST_TEMPLATE})
        super(Exporter, self).__init__(
            template_file='nbsphinx-rst.tpl', extra_loaders=[loader],
            config=traitlets.config.Config({
                'HighlightMagicsPreprocessor': {'enabled': True},
                # Work around https://github.com/jupyter/nbconvert/issues/720:
                'RegexRemovePreprocessor': {'enabled': False},
            }),
            filters={
                'convert_pandoc': convert_pandoc,
                'markdown2rst': markdown2rst,
                'get_empty_lines': _get_empty_lines,
                'extract_gallery_or_toctree': _extract_gallery_or_toctree,
                'save_attachments': save_attachments,
                'replace_attachments': replace_attachments,
                'get_output_type': _get_output_type,
                'json_dumps': lambda s: re.sub(
                    r'<(/script)',
                    r'<\\\1',
                    json.dumps(s),
                    flags=re.IGNORECASE,
                ),
                'basename': os.path.basename,
                'dirname': os.path.dirname,
            })

    def from_notebook_node(self, nb, resources=None, **kw):
        nb = copy.deepcopy(nb)
        if resources is None:
            resources = {}
        else:
            resources = copy.deepcopy(resources)
        # Set default codecell lexer
        resources['codecell_lexer'] = self._codecell_lexer

        nbsphinx_metadata = nb.metadata.get('nbsphinx', {})

        execute = nbsphinx_metadata.get('execute', self._execute)
        if execute not in ('always', 'never', 'auto'):
            raise ValueError('invalid execute option: {!r}'.format(execute))
        auto_execute = (
            execute == 'auto' and
            # At least one code cell actually containing source code:
            any(c.source for c in nb.cells if c.cell_type == 'code') and
            # No outputs, not even a prompt number:
            not any(c.get('outputs') or c.get('execution_count')
                    for c in nb.cells if c.cell_type == 'code')
        )
        if auto_execute or execute == 'always':
            allow_errors = nbsphinx_metadata.get(
                'allow_errors', self._allow_errors)
            timeout = nbsphinx_metadata.get('timeout', self._timeout)
            pp = nbconvert.preprocessors.ExecutePreprocessor(
                kernel_name=self._kernel_name,
                extra_arguments=self._execute_arguments,
                allow_errors=allow_errors, timeout=timeout)
            nb, resources = pp.preprocess(nb, resources)

        if 'nbsphinx_save_notebook' in resources:
            # Save *executed* notebook *before* the Exporter can change it:
            nbformat.write(nb, resources['nbsphinx_save_notebook'])

        # Call into RSTExporter
        rststr, resources = super(Exporter, self).from_notebook_node(
            nb, resources, **kw)

        orphan = nbsphinx_metadata.get('orphan', False)
        if orphan is True:
            resources['nbsphinx_orphan'] = True
        elif orphan is not False:
            raise ValueError('invalid orphan option: {!r}'.format(orphan))

        if 'application/vnd.jupyter.widget-state+json' in nb.metadata.get(
                'widgets', {}):
            resources['nbsphinx_widgets'] = True

        thumbnail = {}

        def warning(msg, *args):
            logger.warning(
                '"nbsphinx-thumbnail" in cell %s: ' + msg, cell_index, *args,
                location=resources.get('nbsphinx_docname'),
                type='nbsphinx', subtype='thumbnail')
            thumbnail['filename'] = _BROKEN_THUMBNAIL

        for cell_index, cell in enumerate(nb.cells):
            if 'nbsphinx-thumbnail' in cell.metadata:
                data = cell.metadata['nbsphinx-thumbnail'].copy()
                output_index = data.pop('output-index', None)
                tooltip = data.pop('tooltip', '')
                if data:
                    warning('Invalid key(s): %s', set(data))
                    break
            elif 'nbsphinx-thumbnail' in cell.metadata.get('tags', []):
                output_index = None
                tooltip = ''
            else:
                continue

            if cell.cell_type != 'code':
                warning(
                    'Only allowed in code cells; wrong cell type: "%s"',
                    cell.cell_type)
                break

            if thumbnail:
                warning('Only allowed once per notebook')
                break

            if output_index is None:
                output_index = len(cell.outputs) - 1
            try:
                suffix = _extract_thumbnail(cell, output_index)
            except _ExtractThumbnailException as e:
                warning(*e.args)
                break

            thumbnail['filename'] = '{}_{}_{}{}'.format(
                resources['unique_key'],
                cell_index,
                output_index,
                suffix,
            )
            if tooltip:
                thumbnail['tooltip'] = tooltip

        if not thumbnail:
            # No explicit thumbnails were specified in the notebook.
            # Now we are looking for the last output image in the notebook.
            for cell_index, cell in reversed(list(enumerate(nb.cells))):
                if cell.cell_type == 'code':
                    # The following can be quickly skipped if there is
                    # less than 1 image output and 2 stream outputs:
                    if len(cell.outputs) >= 3 and self.config.get(
                            'CoalesceStreamsPreprocessor', {}).get(
                                'enabled', False):
                        # CoalesceStreamsPreprocessor was introduced in
                        # nbconvert version 7.14 and enabled in 7.16.4,
                        # see https://github.com/jupyter/nbconvert/pull/2142.
                        from nbconvert.preprocessors \
                            import CoalesceStreamsPreprocessor
                        pp = CoalesceStreamsPreprocessor()
                        # The CoalesceStreamsPreprocessor will be executed as
                        # part of the RSTExporter, but we already have to call
                        # it here to get the correct output indices:
                        cell, _ = pp.preprocess_cell(cell, {}, cell_index)
                        # NB: This correction doesn't happen for nbconvert<7.14
                    for output_index in reversed(range(len(cell.outputs))):
                        try:
                            suffix = _extract_thumbnail(cell, output_index)
                        except _ExtractThumbnailException:
                            continue
                        thumbnail['filename'] = '{}_{}_{}{}'.format(
                            resources['unique_key'],
                            cell_index,
                            output_index,
                            suffix,
                        )
                        # NB: we use this as marker for implicit thumbnail:
                        thumbnail['tooltip'] = None
                        break
                    else:
                        continue
                    break

        if thumbnail:
            resources['nbsphinx_thumbnail'] = thumbnail
        return rststr, resources


class _ExtractThumbnailException(Exception):
    """Internal exception thrown by _extract_thumbnail()."""


def _extract_thumbnail(cell, output_index):
    if not cell.outputs:
        raise _ExtractThumbnailException('No outputs')
    if output_index not in range(len(cell.outputs)):
        raise _ExtractThumbnailException(
            'Invalid "output-index": %s', output_index)
    out = cell.outputs[output_index]
    if out.output_type not in {'display_data', 'execute_result'}:
        raise _ExtractThumbnailException(
            'Unsupported output type in output %s: "%s"',
            output_index, out.output_type)
    for mime_type in DISPLAY_DATA_PRIORITY_HTML:
        if mime_type not in out.data or mime_type not in THUMBNAIL_MIME_TYPES:
            continue
        return THUMBNAIL_MIME_TYPES[mime_type]
    else:
        raise _ExtractThumbnailException(
            'Unsupported MIME type(s) in output %s: %s',
            output_index, set(out.data))


class NotebookParser(rst.Parser):
    """Sphinx source parser for Jupyter notebooks.

    Uses nbsphinx.Exporter to convert notebook content to a
    reStructuredText string, which is then parsed by Sphinx's built-in
    reST parser.

    """

    supported = 'jupyter_notebook',

    def get_transforms(self):
        """List of transforms for documents parsed by this parser."""
        return rst.Parser.get_transforms(self) + [
            CreateNotebookSectionAnchors,
            ReplaceAlertDivs,
            CopyLinkedFiles,
            ForceEquations,
        ]

    def parse(self, inputstring, document):
        """Parse *inputstring*, write results to *document*.

        *inputstring* is either the JSON representation of a notebook,
        or a paragraph of text coming from the Sphinx translation
        machinery.

        Note: For now, the translation strings use reST formatting,
        because the NotebookParser uses reST as intermediate
        representation.
        However, there are plans to remove this intermediate step
        (https://github.com/spatialaudio/nbsphinx/issues/36), and after
        that, the translated strings will most likely be parsed as
        CommonMark.

        If the configuration value "nbsphinx_custom_formats" is
        specified, the input string is converted to the Jupyter notebook
        format with the given conversion function.

        """
        env = document.settings.env
        formats = {
            '.ipynb': lambda s: nbformat.reads(s, as_version=_ipynbversion)}
        formats.update(env.config.nbsphinx_custom_formats)
        srcfile = str(env.doc2path(env.docname, base=None))
        for format, converter in formats.items():
            if srcfile.endswith(format):
                break
        else:
            raise NotebookError(
                'No converter was found for {!r}'.format(srcfile))
        if (isinstance(converter, collections.abc.Sequence) and
                not isinstance(converter, str)):
            if len(converter) != 2:
                raise NotebookError(
                    'The values of nbsphinx_custom_formats must be '
                    'either strings or 2-element sequences')
            converter, kwargs = converter
        else:
            kwargs = {}
        if isinstance(converter, str):
            converter = sphinx.util.import_object(converter)
        try:
            nb = converter(inputstring, **kwargs)
        except Exception:
            # NB: The use of the RST parser is temporary!
            rst.Parser.parse(self, inputstring, document)
            return

        srcdir = os.path.dirname(env.doc2path(env.docname))
        auxdir = env.nbsphinx_auxdir

        resources = {}
        # Working directory for ExecutePreprocessor
        resources['metadata'] = {'path': srcdir}
        # Sphinx doesn't accept absolute paths in images etc.
        resources['output_files_dir'] = os.path.relpath(auxdir, srcdir)
        resources['unique_key'] = re.sub('[/ ]', '_', env.docname)
        resources['nbsphinx_docname'] = env.docname

        # NB: The source file could have a different suffix
        #     if nbsphinx_custom_formats is used.
        notebookfile = env.docname + '.ipynb'
        env.nbsphinx_notebooks[env.docname] = notebookfile
        auxfile = os.path.join(auxdir, notebookfile)
        sphinx.util.ensuredir(os.path.dirname(auxfile))
        resources['nbsphinx_save_notebook'] = auxfile

        exporter = Exporter(
            execute=env.config.nbsphinx_execute,
            kernel_name=env.config.nbsphinx_kernel_name,
            execute_arguments=env.config.nbsphinx_execute_arguments,
            allow_errors=env.config.nbsphinx_allow_errors,
            timeout=env.config.nbsphinx_timeout,
            codecell_lexer=env.config.nbsphinx_codecell_lexer,
        )

        try:
            rststring, resources = exporter.from_notebook_node(nb, resources)
        except nbconvert.preprocessors.CellExecutionError as e:
            lines = str(e).split('\n')
            lines[0] = 'CellExecutionError in {}:'.format(
                env.doc2path(env.docname, base=None))
            lines.append("You can ignore this error by setting the following "
                         "in conf.py:\n\n    nbsphinx_allow_errors = True\n")
            raise NotebookError('\n'.join(lines))
        except Exception as e:
            raise NotebookError(
                type(e).__name__ + ' in ' +
                str(env.doc2path(env.docname, base=None)) + ':\n' + str(e))

        rststring = """
.. role:: nbsphinx-math(raw)
    :format: latex + html
    :class: math

..

""" + rststring

        # Create additional output files (figures etc.),
        # see nbconvert.writers.FilesWriter.write()
        for filename, data in resources.get('outputs', {}).items():
            dest = os.path.normpath(os.path.join(srcdir, filename))
            with open(dest, 'wb') as f:
                f.write(data)

        if resources.get('nbsphinx_orphan', False):
            rst.Parser.parse(self, ':orphan:', document)
        if env.config.nbsphinx_prolog:
            prolog = exporter.environment.from_string(
                env.config.nbsphinx_prolog).render(env=env)
            rst.Parser.parse(self, prolog, document)
        rst.Parser.parse(self, '.. highlight:: none', document)
        if 'sphinx_codeautolink' in env.config.extensions:
            rst.Parser.parse(self, '.. autolink-concat:: on', document)
        rst.Parser.parse(self, rststring, document)
        if env.config.nbsphinx_epilog:
            epilog = exporter.environment.from_string(
                env.config.nbsphinx_epilog).render(env=env)
            rst.Parser.parse(self, epilog, document)

        if resources.get('nbsphinx_widgets', False):
            env.nbsphinx_widgets.add(env.docname)

        env.nbsphinx_thumbnails[env.docname] = resources.get(
            'nbsphinx_thumbnail', {})


class NotebookError(sphinx.errors.SphinxError):
    """Error during notebook parsing."""

    category = 'Notebook error'


class CodeAreaNode(docutils.nodes.Element):
    """Input area or plain-text output area of a Jupyter notebook code cell."""


class FancyOutputNode(docutils.nodes.Element):
    """A custom node for non-plain-text output of code cells."""


def _create_code_nodes(directive):
    """Create nodes for an input or output code cell."""
    directive.state.document['nbsphinx_code_css'] = True
    execution_count = directive.options.get('execution-count')
    config = directive.state.document.settings.env.config
    if isinstance(directive, NbInput):
        outer_classes = ['nbinput']
        if 'no-output' in directive.options:
            outer_classes.append('nblast')
        inner_classes = ['input_area']
        prompt_template = config.nbsphinx_input_prompt
        if not execution_count:
            execution_count = ' '
    elif isinstance(directive, NbOutput):
        outer_classes = ['nboutput']
        if 'more-to-come' not in directive.options:
            outer_classes.append('nblast')
        inner_classes = ['output_area']
        inner_classes.append(directive.options.get('class', ''))
        prompt_template = config.nbsphinx_output_prompt
    else:
        assert False

    outer_node = docutils.nodes.container(classes=outer_classes)
    if execution_count:
        prompt = prompt_template % (execution_count,)
        prompt_node = docutils.nodes.literal_block(
            prompt, prompt, language='none', classes=['prompt'])
    else:
        prompt = ''
        prompt_node = docutils.nodes.container(classes=['prompt', 'empty'])
    # NB: Prompts are added manually in LaTeX output
    outer_node += sphinx.addnodes.only('', prompt_node, expr='html')

    if isinstance(directive, NbInput):
        text = '\n'.join(directive.content.data)
        if directive.arguments:
            language = directive.arguments[0]
        else:
            language = 'none'
        inner_node = docutils.nodes.literal_block(
            text, text, language=language, classes=inner_classes)
    else:
        inner_node = docutils.nodes.container(classes=inner_classes)
        sphinx.util.nodes.nested_parse_with_titles(
            directive.state, directive.content, inner_node)

    if 'fancy' in directive.options:
        outer_node += FancyOutputNode('', inner_node, prompt=prompt)
    else:
        codearea_node = CodeAreaNode(
            '', inner_node, prompt=prompt, stderr='stderr' in inner_classes)
        # See http://stackoverflow.com/q/34050044/.
        for attr in 'empty-lines-before', 'empty-lines-after':
            value = directive.options.get(attr, 0)
            if value:
                codearea_node[attr] = value
        outer_node += codearea_node
    return [outer_node]


class AdmonitionNode(docutils.nodes.Element):
    """A custom node for info and warning boxes."""


class GalleryToc(docutils.nodes.Element):
    """A wrapper node used for creating galleries."""


class GalleryNode(docutils.nodes.Element):
    """A custom node for thumbnail galleries."""


class DummyTocTree(docutils.nodes.Element):
    """A dummy node to replace and disable sphinx.addnodes.toctree."""


# See http://docutils.sourceforge.net/docs/howto/rst-directives.html

class NbInput(rst.Directive):
    """A notebook input cell with prompt and code area."""

    required_arguments = 0
    optional_arguments = 1  # lexer name
    final_argument_whitespace = False
    option_spec = {
        'execution-count': rst.directives.positive_int,
        'empty-lines-before': rst.directives.nonnegative_int,
        'empty-lines-after': rst.directives.nonnegative_int,
        'no-output': rst.directives.flag,
    }
    has_content = True

    def run(self):
        """This is called by the reST parser."""
        return _create_code_nodes(self)


class NbOutput(rst.Directive):
    """A notebook output cell with optional prompt."""

    required_arguments = 0
    final_argument_whitespace = False
    option_spec = {
        'execution-count': rst.directives.positive_int,
        'more-to-come': rst.directives.flag,
        'fancy': rst.directives.flag,
        'class': rst.directives.unchanged,
    }
    has_content = True

    def run(self):
        """This is called by the reST parser."""
        return _create_code_nodes(self)


class _NbAdmonition(rst.Directive):
    """Base class for NbInfo and NbWarning."""

    required_arguments = 0
    optional_arguments = 0
    option_spec = {}
    has_content = True

    def run(self):
        """This is called by the reST parser."""
        node = AdmonitionNode(classes=['admonition', self._class])
        self.state.nested_parse(self.content, self.content_offset, node)
        return [node]


class NbInfo(_NbAdmonition):
    """An information box."""

    _class = 'note'


class NbWarning(_NbAdmonition):
    """A warning box."""

    _class = 'warning'


class NbGallery(sphinx.directives.other.TocTree):
    """A thumbnail gallery for notebooks."""

    def run(self):
        """Wrap GalleryToc around toctree."""
        ret = super().run()
        try:
            toctree_wrapper = ret[-1]
            toctree, = toctree_wrapper
        except (IndexError, TypeError, ValueError):
            return ret
        if not isinstance(toctree, sphinx.addnodes.toctree):
            return ret
        gallerytoc = GalleryToc()
        gallerytoc.extend(ret)
        if isinstance(self, NbLinkGallery):
            # Disable toctree processing:
            toctree.__class__ = DummyTocTree
            # Make invisible to LaTeX processing
            gallerytoc = sphinx.addnodes.only('', gallerytoc, expr='html')
        return [gallerytoc]


class NbLinkGallery(NbGallery):
    """A thumbnail gallery for notebooks as links."""

    # Not all options of TocTree are allowed:
    option_spec = {
        'name': rst.directives.unchanged,
        'caption': rst.directives.unchanged_required,
        'glob': rst.directives.flag,
        'reversed': rst.directives.flag,
    }


def convert_pandoc(text, from_format, to_format):
    """Simple wrapper for markdown2rst.

    In nbconvert version 5.0, the use of markdown2rst in the RST
    template was replaced by the new filter function convert_pandoc.

    """
    if from_format != 'markdown' and to_format != 'rst':
        raise ValueError('Unsupported conversion')
    return markdown2rst(text)


_CITE_ROLES = {
    'data-cite': 'cite',
    'data-footcite': 'footcite',
}
_CITE_ROLES.update({
    f'{k}-{suffix}': f'{v}:{suffix}'
    for k, v in _CITE_ROLES.items()
    for suffix in ('p', 'ps', 't', 'ts', 'ct', 'cts')
})


class CitationParser(html.parser.HTMLParser):

    def handle_starttag(self, tag, attrs):
        if self._check_cite(attrs):
            self.starttag = tag

    def handle_endtag(self, tag):
        self.endtag = tag

    def handle_startendtag(self, tag, attrs):
        self._check_cite(attrs)

    def _check_cite(self, attrs):
        for name, value in attrs:
            try:
                self.cite = f':{_CITE_ROLES[name]}:`{value}`'
            except KeyError:
                continue
            return True
        return False

    def reset(self):
        super().reset()
        self.starttag = ''
        self.endtag = ''
        self.cite = ''


class ImgParser(html.parser.HTMLParser):
    """Turn HTML <img> tags into raw RST blocks."""

    def handle_starttag(self, tag, attrs):
        self._check_img(tag, attrs)

    def handle_startendtag(self, tag, attrs):
        self._check_img(tag, attrs)

    def _check_img(self, tag, attrs):
        if tag != 'img':
            return
        # NB: attrs is a list of pairs
        attrs = dict(attrs)
        if 'src' not in attrs:
            return
        img_path = nbconvert.filters.posix_path(attrs['src'])
        if img_path.startswith('data'):
            # Allow multi-line data, see issue #474
            img_path = img_path.replace('\n', '')
        lines = ['image:: ' + img_path]
        indent = ' ' * 4
        classes = []
        if 'class' in attrs:
            classes.append(attrs['class'])
        if 'alt' in attrs:
            lines.append(indent + ':alt: ' + attrs['alt'])
        if 'width' in attrs:
            lines.append(indent + ':width: ' + attrs['width'])
        if 'height' in attrs:
            lines.append(indent + ':height: ' + attrs['height'])
        if {'width', 'height'}.intersection(attrs):
            classes.append('no-scaled-link')
        if classes:
            lines.append(indent + ':class: ' + ' '.join(classes))

        definition = '\n'.join(lines)
        hex_id = uuid.uuid4().hex
        definition = '.. |' + hex_id + '| ' + definition
        self.obj = {'t': 'RawInline', 'c': ['rst', '|' + hex_id + '|']}
        self.definition = definition

    def reset(self):
        super().reset()
        self.obj = {}


def markdown2rst(text):
    """Convert a Markdown string to reST via pandoc.

    This is very similar to nbconvert.filters.markdown.markdown2rst(),
    except that it uses a pandoc filter to convert raw LaTeX blocks to
    "math" directives (instead of "raw:: latex" directives).

    NB: At some point, pandoc changed its behavior!  In former times,
    it converted LaTeX math environments to RawBlock ("latex"), at some
    later point this was changed to RawInline ("tex").
    Either way, we convert it to Math/DisplayMath.

    """

    def parse_citation(obj):
        p = CitationParser()
        p.feed(obj['c'][1])
        p.close()
        return p

    def parse_img(obj):
        p = ImgParser()
        p.feed(obj['c'][1])
        p.close()
        return p

    def object_hook(obj):
        if object_hook.open_cite_tag:
            if obj.get('t') == 'RawInline' and obj['c'][0] == 'html':
                p = parse_citation(obj)
                if p.endtag == object_hook.open_cite_tag:
                    object_hook.open_cite_tag = ''
            return {'t': 'Str', 'c': ''}  # Object is replaced by empty string

        if obj.get('t') == 'RawBlock' and obj['c'][0] == 'latex':
            obj['t'] = 'Para'
            obj['c'] = [{
                't': 'Math',
                'c': [
                    {'t': 'DisplayMath', 'c': []},
                    # Special marker characters are removed below:
                    '\x0e:nowrap:\x0f\n\n' + obj['c'][1],
                ]
            }]
        elif obj.get('t') == 'RawInline' and obj['c'][0] == 'tex':
            obj = {'t': 'RawInline',
                   'c': ['rst', ':nbsphinx-math:`{}`'.format(obj['c'][1])]}
        elif obj.get('t') == 'RawInline' and obj['c'][0] == 'html':
            p = parse_citation(obj)
            if p.starttag:
                object_hook.open_cite_tag = p.starttag
            if p.cite:
                obj = {'t': 'RawInline', 'c': ['rst', p.cite]}
            if not p.starttag and not p.cite:
                p = parse_img(obj)
                if p.obj:
                    obj = p.obj
                    object_hook.image_definitions.append(p.definition)
        return obj

    object_hook.open_cite_tag = ''
    object_hook.image_definitions = []

    def filter_func(text):
        json_data = json.loads(text, object_hook=object_hook)
        return json.dumps(json_data)

    input_format = 'markdown'
    input_format += '-implicit_figures'
    input_format += '+lists_without_preceding_blankline'
    v = nbconvert.utils.pandoc.get_pandoc_version()
    if nbconvert.utils.version.check_version(v, '1.13'):
        input_format += '-native_divs+raw_html'
    if nbconvert.utils.version.check_version(v, '2.0'):
        input_format += '-smart'  # Smart quotes etc. are handled by Sphinx

    rststring = pandoc(text, input_format, 'rst', filter_func=filter_func)
    rststring = re.sub(
        r'^\n( *)\x0e:nowrap:\x0f$',
        r'\1:nowrap:',
        rststring,
        flags=re.MULTILINE)
    rststring += '\n\n'
    rststring += '\n'.join(object_hook.image_definitions)
    return rststring


def pandoc(source, fmt, to, filter_func=None):
    """Convert a string in format `from` to format `to` via pandoc.

    This is based on nbconvert.utils.pandoc.pandoc() and extended to
    allow passing a filter function.

    """
    def encode(text):
        return text if isinstance(text, bytes) else text.encode('utf-8')

    def decode(data):
        return data.decode('utf-8') if isinstance(data, bytes) else data

    nbconvert.utils.pandoc.check_pandoc_version()
    v = nbconvert.utils.pandoc.get_pandoc_version()
    cmd = ['pandoc']
    if nbconvert.utils.version.check_version(v, '2.0'):
        # see issue #155
        cmd += ['--eol', 'lf']
    cmd1 = cmd + ['--from', fmt, '--to', 'json']
    cmd2 = cmd + ['--from', 'json', '--to', to]
    cmd2 += ['--columns=500']  # Avoid breaks in tables, see issue #240

    p = subprocess.Popen(cmd1, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    json_data, _ = p.communicate(encode(source))

    if filter_func:
        json_data = encode(filter_func(decode(json_data)))

    p = subprocess.Popen(cmd2, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    out, _ = p.communicate(json_data)
    return decode(out).rstrip('\n')


def _extract_gallery_or_toctree(cell):
    """Extract links from Markdown cell and create gallery/toctree."""
    # If more than one is available, "gallery" takes precedence, then "toctree"
    if 'nbsphinx-gallery' in cell.metadata:
        lines = ['.. nbgallery::']
        options = cell.metadata['nbsphinx-gallery']
    elif 'nbsphinx-gallery' in cell.metadata.get('tags', []):
        lines = ['.. nbgallery::']
        options = {}
    elif 'nbsphinx-toctree' in cell.metadata:
        lines = ['.. toctree::']
        options = cell.metadata['nbsphinx-toctree']
    elif 'nbsphinx-toctree' in cell.metadata.get('tags', []):
        lines = ['.. toctree::']
        options = {}
    elif 'nbsphinx-link-gallery' in cell.metadata:
        lines = ['.. nblinkgallery::']
        options = cell.metadata['nbsphinx-link-gallery']
    elif 'nbsphinx-link-gallery' in cell.metadata.get('tags', []):
        lines = ['.. nblinkgallery::']
        options = {}
    else:
        assert False
    try:
        for option, value in options.items():
            if value is True:
                lines.append(':{}:'.format(option))
            elif value is False:
                pass
            else:
                lines.append(':{}: {}'.format(option, value))
    except AttributeError:
        raise ValueError(
            'invalid nbsphinx-[link-]gallery/nbsphinx-toctree option: {!r}'
            .format(options))

    text = markdown2rst(cell.source)
    settings = docutils.frontend.OptionParser(
        components=(rst.Parser,)).get_default_values()
    node = docutils.utils.new_document('gallery_or_toctree', settings)
    parser = rst.Parser()
    parser.parse(text, node)

    if 'caption' not in options:
        for sec in node.findall(docutils.nodes.section):
            assert sec.children
            assert isinstance(sec.children[0], docutils.nodes.title)
            title = sec.children[0].astext()
            lines.append(':caption: ' + title)
            break
    lines.append('')  # empty line
    for ref in node.findall(docutils.nodes.reference):
        lines.append(ref['name'] + ' <' + unquote(ref.get('refuri')) + '>')
    return '\n    '.join(lines)


def _get_empty_lines(text):
    """Get number of empty lines before and after code."""
    before = len(text) - len(text.lstrip('\n'))
    after = len(text) - len(text.strip('\n')) - before
    return before, after


def _get_output_type(output):
    """Choose appropriate output data types for HTML and LaTeX."""
    if output.output_type == 'stream':
        html_datatype = latex_datatype = 'text/plain'
        text = output.text
        output.data = {'text/plain': text[:-1] if text.endswith('\n') else text}
    elif output.output_type == 'error':
        html_datatype = latex_datatype = 'text/plain'
        output.data = {'text/plain': '\n'.join(output.traceback)}
    else:
        for datatype in DISPLAY_DATA_PRIORITY_HTML:
            if datatype in output.data:
                html_datatype = datatype
                break
        else:
            html_datatype = ', '.join(output.data.keys())
        for datatype in DISPLAY_DATA_PRIORITY_LATEX:
            if datatype in output.data:
                latex_datatype = datatype
                break
        else:
            latex_datatype = ', '.join(output.data.keys())
    return html_datatype, latex_datatype


def _local_file_from_reference(node, document):
    """Get local file path from reference and detect fragment identifier."""
    # NB: Anonymous hyperlinks must be already resolved at this point!
    refuri = node.get('refuri')
    if not refuri:
        refname = node.get('refname')
        if refname:
            refid = document.nameids.get(refname)
        else:
            # NB: This can happen for anonymous hyperlinks
            refid = node.get('refid')
        target = document.ids.get(refid)
        if not target:
            # No corresponding target, Sphinx may warn later
            return '', ''
        refuri = target.get('refuri')
        if not refuri:
            # Target doesn't have URI
            return '', ''
    if '://' in refuri or refuri.startswith('mailto:'):
        # Not a local link
        return '', ''
    elif refuri.startswith('#'):
        # Kinda not a local link
        return '', refuri

    # NB: We look for "fragment identifier" before unquoting
    match = re.match(r'^([^#]*)(#.*)$', refuri)
    if match:
        filename = unquote(match.group(1))
        # NB: The "fragment identifier" is not unquoted
        fragment = match.group(2)
    else:
        filename = unquote(refuri)
        fragment = ''
    return filename, fragment


class RewriteLocalLinks(docutils.transforms.Transform):
    """Turn links to source files into ``:doc:``/``:ref:`` links.

    Links to subsections are possible with ``...#Subsection-Title``.
    These links use the labels created by CreateSectionLabels.

    Links to subsections use ``:ref:``, links to whole source files use
    ``:doc:``.  Latter can be useful if you have an ``index.rst`` but
    also want a distinct ``index.ipynb`` for use with Jupyter.
    In this case you can use such a link in a notebook::

        [Back to main page](index.ipynb)

    In Jupyter, this will create a "normal" link to ``index.ipynb``, but
    in the files generated by Sphinx, this will become a link to the
    main page created from ``index.rst``.

    """

    default_priority = 500  # After AnonymousHyperlinks (440)

    def apply(self):
        env = self.document.settings.env
        for node in self.document.findall(docutils.nodes.reference):
            filename, fragment = _local_file_from_reference(
                node, self.document)
            if not filename:
                if fragment:
                    # This should not be needed if it weren't for
                    # https://github.com/sphinx-doc/sphinx/issues/11336 and
                    # https://github.com/sphinx-doc/sphinx/issues/11335
                    filename = os.path.basename(env.doc2path(env.docname))
                    if fragment == '#':
                        fragment = ''
                else:
                    continue

            for s in env.config.source_suffix:
                if filename.lower().endswith(s.lower()):
                    assert len(s) > 0
                    target = filename[:-len(s)]
                    suffix = filename[-len(s):]
                    if fragment:
                        target_ext = suffix + fragment
                        reftype = 'ref'
                    else:
                        target_ext = ''
                        reftype = 'doc'
                    break
            else:
                continue  # Not a link to a potential Sphinx source file

            target_docname = nbconvert.filters.posix_path(os.path.normpath(
                os.path.join(os.path.dirname(env.docname), target)))
            if target_docname in env.found_docs:
                reftarget = '/' + target_docname + target_ext
                if reftype == 'ref':
                    reftarget = reftarget.lower()
                linktext = node.astext()
                xref = sphinx.addnodes.pending_xref(
                    reftype=reftype, reftarget=reftarget, refdomain='std',
                    refwarn=True, refexplicit=True, refdoc=env.docname)
                xref += docutils.nodes.Text(linktext)
                node.replace_self(xref)
            else:
                # NB: This is a link to an ignored (via exclude_patterns)
                #     source file.
                pass


class CreateNotebookSectionAnchors(docutils.transforms.Transform):
    """Create section anchors for Jupyter notebooks.

    Note: Sphinx lower-cases the HTML section IDs, Jupyter doesn't.
    This transform creates anchors in the Jupyter style.

    """

    default_priority = 200  # Before CreateSectionLabels (250)

    def apply(self):
        all_ids = set()
        for section in self.document.findall(docutils.nodes.section):
            title = section[0].astext()
            link_id = title.replace(' ', '-').replace('"', '%22')
            if link_id in all_ids:
                # Avoid duplicated anchors on the same page
                continue
            all_ids.add(link_id)
            section['ids'] = [link_id]
        if not all_ids:
            logger.warning(
                'Each notebook should have at least one section title',
                location=self.document[0],
                type='nbsphinx', subtype='notebooktitle')


class CreateSectionLabels(docutils.transforms.Transform):
    """Make labels for each document and each section thereof.

    These labels are referenced in RewriteLocalLinks but can also be
    used manually with ``:ref:``.

    """

    default_priority = 250  # Before references.PropagateTargets (260)

    def apply(self):
        env = self.document.settings.env
        file_ext = str(env.doc2path(env.docname, base=None))[len(env.docname):]
        i_still_have_to_create_the_document_label = True
        for section in self.document.findall(docutils.nodes.section):
            assert section.children
            assert isinstance(section.children[0], docutils.nodes.title)
            title = section.children[0].astext()
            link_id = section['ids'][0]
            label = '/' + env.docname + file_ext + '#' + link_id
            label = label.lower()
            env.domaindata['std']['labels'][label] = (
                env.docname, link_id, title)
            env.domaindata['std']['anonlabels'][label] = (
                env.docname, link_id)

            # Create a label for the whole document using the first section:
            if i_still_have_to_create_the_document_label:
                label = '/' + env.docname.lower() + file_ext
                env.domaindata['std']['labels'][label] = (
                    env.docname, '', title)
                env.domaindata['std']['anonlabels'][label] = (
                    env.docname, '')
                i_still_have_to_create_the_document_label = False


class CreateDomainObjectLabels(docutils.transforms.Transform):
    """Create labels for domain-specific object signatures."""

    default_priority = 250  # About the same as CreateSectionLabels

    def apply(self):
        env = self.document.settings.env
        file_ext = str(env.doc2path(env.docname, base=None))[len(env.docname):]
        for sig in self.document.findall(sphinx.addnodes.desc_signature):
            try:
                title = sig['ids'][0]
            except IndexError:
                # Object has same name as another, so skip it
                continue
            link_id = title.replace(' ', '-')
            sig['ids'] = [link_id]
            label = '/' + env.docname + file_ext + '#' + link_id
            label = label.lower()
            env.domaindata['std']['labels'][label] = (
                env.docname, link_id, title)
            env.domaindata['std']['anonlabels'][label] = (
                env.docname, link_id)


class ReplaceAlertDivs(docutils.transforms.Transform):
    """Replace certain <div> elements with AdmonitionNode containers.

    This is a quick-and-dirty work-around until a proper
    Mardown/CommonMark extension for note/warning boxes is available.

    """

    default_priority = 500  # Doesn't really matter

    _start_re = re.compile(
        r'\s*<div\s*class\s*=\s*(?P<q>"|\')([a-z\s-]*)(?P=q)\s*>\s*\Z',
        flags=re.IGNORECASE)
    _class_re = re.compile(r'\s*alert\s*alert-(info|warning)\s*\Z')
    _end_re = re.compile(r'\s*</div\s*>\s*\Z', flags=re.IGNORECASE)

    def apply(self):
        start_tags = []
        for node in self.document.findall(docutils.nodes.raw):
            if node['format'] != 'html':
                continue
            start_match = self._start_re.match(node.astext())
            if not start_match:
                continue
            class_match = self._class_re.match(start_match.group(2))
            if not class_match:
                continue
            admonition_class = class_match.group(1)
            if admonition_class == 'info':
                admonition_class = 'note'
            start_tags.append((node, admonition_class))

        # Reversed order to allow nested <div> elements:
        for node, admonition_class in reversed(start_tags):
            content = []
            for sibling in node.findall(include_self=False, descend=False,
                                         siblings=True, ascend=False):
                end_tag = (isinstance(sibling, docutils.nodes.raw) and
                           sibling['format'] == 'html' and
                           self._end_re.match(sibling.astext()))
                if end_tag:
                    admonition_node = AdmonitionNode(
                        classes=['admonition', admonition_class])
                    admonition_node.extend(content)
                    parent = node.parent
                    parent.replace(node, admonition_node)
                    for n in content:
                        parent.remove(n)
                    parent.remove(sibling)
                    break
                else:
                    content.append(sibling)


class CopyLinkedFiles(docutils.transforms.Transform):
    """Mark linked (local) files to be copied to the HTML output."""

    default_priority = 600  # After RewriteLocalLinks

    def apply(self):
        env = self.document.settings.env
        for node in self.document.findall(docutils.nodes.reference):
            filename, fragment = _local_file_from_reference(
                node, self.document)
            if not filename:
                continue  # Not a local link
            relpath = filename + fragment
            file = os.path.normpath(
                os.path.join(os.path.dirname(env.docname), relpath))
            if not os.path.isfile(os.path.join(env.srcdir, file)):
                logger.warning(
                    'File not found: %r', file, location=node,
                    type='nbsphinx', subtype='localfile')
                continue  # Link is ignored
            elif file.startswith('..'):
                logger.warning(
                    'Link outside source directory: %r', file, location=node,
                    type='nbsphinx', subtype='localfile')
                continue  # Link is ignored
            env.nbsphinx_files.setdefault(env.docname, []).append(file)


class ForceEquations(docutils.transforms.Transform):
    """Unconditionally enable equations on notebooks.

    Except if ``nbsphinx_assume_equations`` is set to ``False``.

    """

    default_priority = 900  # after checking for equations in MathDomain

    def apply(self):
        env = self.document.settings.env
        if env.config.nbsphinx_assume_equations:
            env.get_domain('math').data['has_equations'][env.docname] = True


class GetSizeFromImages(
        sphinx.transforms.post_transforms.images.BaseImageConverter):
    """Get size from images and store it as node attributes.

    This is only done for LaTeX output.

    """

    # After ImageDownloader (100) and DataURIExtractor (150):
    default_priority = 200

    def match(self, node):
        return self.app.builder.format == 'latex'

    def handle(self, node):
        if 'width' not in node and 'height' not in node:
            srcdir = os.path.dirname(self.env.doc2path(self.env.docname))
            image_path = os.path.normpath(os.path.join(srcdir, node['uri']))
            size = sphinx.util.images.get_image_size(image_path)
            if size is not None:
                node['width'], node['height'] = map(str, size)


if hasattr(sphinx.environment.adapters.toctree, '_resolve_toctree'):
    # Since Sphinx 7.2.0
    original_toctree_resolve = \
        sphinx.environment.adapters.toctree._resolve_toctree
else:
    original_toctree_resolve = \
        sphinx.environment.adapters.toctree.TocTree.resolve


def patched_toctree_resolve(self, docname, builder, toctree, *args, **kwargs):
    """Method for monkey-patching Sphinx's TocTree adapter.

    The list of section links is never shown, regardless of the
    ``:hidden:`` option.
    However, this option can still be used to control whether the
    section links are shown in higher-level tables of contents.

    """
    gallery = toctree.get('nbsphinx_gallery', False)
    if gallery:
        toctree = toctree.copy()
        toctree['hidden'] = False
    node = original_toctree_resolve(
        self, docname, builder, toctree, *args, **kwargs)
    if not gallery or node is None:
        return node
    if isinstance(node[0], (docutils.nodes.caption, docutils.nodes.title)):
        del node[1:]
    else:
        del node[:]
    return node


def config_inited(app, config):
    for suffix in config.nbsphinx_custom_formats:
        app.add_source_suffix(suffix, 'jupyter_notebook')
    if '.ipynb' in chain(config.source_suffix, app.registry.source_suffix):
        # We don't interfere if '.ipynb' has been added by the user in conf.py
        # or by another Sphinx extension that has been loaded earlier.
        pass
    else:
        app.add_source_suffix('.ipynb', 'jupyter_notebook')

    if '**.ipynb_checkpoints' not in config.exclude_patterns:
        config.exclude_patterns.append('**.ipynb_checkpoints')

    # Make sure require.js is loaded after all other extensions,
    # see https://github.com/spatialaudio/nbsphinx/issues/409
    app.connect('builder-inited', load_requirejs)

    # http://docs.mathjax.org/en/v3.1-latest/options/document.html
    # http://docs.mathjax.org/en/v2.7-latest/options/preprocessors/tex2jax.html
    mathjax_inline_math = [['$', '$'], ['\\(', '\\)']]
    mathjax_ignore_class = (
        'tex2jax_ignore'  # MathJax 2 default
        '|'
        'mathjax_ignore'  # Mathjax 3 default
        '|'
        'document'  # Main page content
    )
    mathjax_process_class = (
        'tex2jax_process'  # MathJax 2 default
        '|'
        'mathjax_process'  # Mathjax 3 default
        '|'
        'math'  # Used by Sphinx
        '|'
        'output_area'  # Jupyter code cells
    )

    # See also https://github.com/sphinx-doc/sphinx/pull/5504
    if hasattr(config, 'mathjax3_config') and config.mathjax2_config is None:
        # NB: If mathjax_path is used in Sphinx >= 4 to load MathJax v2,
        # this only works if mathjax_config or mathjax2_config is specified.
        if config.mathjax3_config is None:
            config.mathjax3_config = {}
        mathjax3_config = config.mathjax3_config
        tex = {
            'inlineMath': mathjax_inline_math,
            'processEscapes': True,
        }
        tex.update(mathjax3_config.get('tex', {}))
        mathjax3_config['tex'] = tex
        options = {
            'ignoreHtmlClass': mathjax_ignore_class,
            'processHtmlClass': mathjax_process_class,
        }
        options.update(mathjax3_config.get('options', {}))
        mathjax3_config['options'] = options
    else:
        if hasattr(config, 'mathjax2_config'):
            # Sphinx >= 4.0
            if config.mathjax2_config is None:
                config.mathjax2_config = {}
            mathjax2_config = config.mathjax2_config
        else:
            # Sphinx < 4.0
            if config.mathjax_config is None:
                config.mathjax_config = {}
            mathjax2_config = config.mathjax_config
        tex2jax = {
            'inlineMath': mathjax_inline_math,
            'processEscapes': True,
            'ignoreClass': mathjax_ignore_class,
            'processClass': mathjax_process_class,
        }
        tex2jax.update(mathjax2_config.get('tex2jax', {}))
        mathjax2_config['tex2jax'] = tex2jax


def load_requirejs(app):
    config = app.config
    if config.nbsphinx_requirejs_path is None:
        config.nbsphinx_requirejs_path = 'https://cdnjs.cloudflare.com/ajax/libs/require.js/2.3.4/require.min.js'
    if config.nbsphinx_requirejs_options is None:
        config.nbsphinx_requirejs_options = {
            'integrity': 'sha256-Ae2Vz/4ePdIu6ZyI/5ZGsYnb+m0JlOmKPjt6XZ9JJkA=',
            'crossorigin': 'anonymous',
        }
    if config.nbsphinx_requirejs_path:
        # TODO: Add only on pages created from notebooks?
        app.add_js_file(
            config.nbsphinx_requirejs_path,
            **config.nbsphinx_requirejs_options)


def builder_inited(app):
    env = app.env
    env.settings['line_length_limit'] = 100_000_000
    env.nbsphinx_notebooks = {}
    env.nbsphinx_files = {}
    if not hasattr(env, 'nbsphinx_thumbnails'):
        env.nbsphinx_thumbnails = {}
    env.nbsphinx_widgets = set()
    env.nbsphinx_auxdir = os.path.join(env.doctreedir, 'nbsphinx')
    sphinx.util.ensuredir(env.nbsphinx_auxdir)

    if app.builder.format == 'latex':
        sphinx.util.fileutil.copy_asset(
            os.path.join(os.path.dirname(__file__), '_texinputs'),
            os.path.join(app.builder.outdir))

    if app.builder.format == 'html':
        context = {
            'nbsphinx_responsive_width': app.config.nbsphinx_responsive_width,
            'nbsphinx_prompt_width': app.config.nbsphinx_prompt_width,
        }
        assets = (
            'nbsphinx-code-cells.css_t',
            'nbsphinx-gallery.css',
            'nbsphinx-no-thumbnail.svg',
            'nbsphinx-broken-thumbnail.svg',
        )
        for a in assets:
            sphinx.util.fileutil.copy_asset(
                os.path.join(os.path.dirname(__file__), '_static', a),
                os.path.join(app.builder.outdir, '_static'),
                context=context)


def env_merge_info(app, env, docnames, other):
    env.nbsphinx_notebooks.update(other.nbsphinx_notebooks)
    env.nbsphinx_files.update(other.nbsphinx_files)
    env.nbsphinx_thumbnails.update(other.nbsphinx_thumbnails)
    env.nbsphinx_widgets.update(other.nbsphinx_widgets)


def env_purge_doc(app, env, docname):
    """Remove list of local files for a given document."""
    try:
        del env.nbsphinx_notebooks[docname]
    except KeyError:
        pass
    try:
        del env.nbsphinx_files[docname]
    except KeyError:
        pass
    try:
        del env.nbsphinx_thumbnails[docname]
    except KeyError:
        pass
    env.nbsphinx_widgets.discard(docname)


def html_page_context(app, pagename, templatename, context, doctree):
    """Add CSS files for code cells and galleries."""
    # NB: the CSS files are copied in html_collect_pages().
    if doctree and doctree.get('nbsphinx_code_css'):
        app.add_css_file('nbsphinx-code-cells.css')
    if doctree and doctree.get('nbsphinx_gallery_css'):
        app.add_css_file('nbsphinx-gallery.css')


def backwards_compat_overwrite(copyfile=sphinx.util.copyfile):
    """Return kwargs dictionary to pass to copyfile() for consistent behavior

    Before version 8 of Sphinx, the default behavior of the copyfile function
    was to overwrite the file at the target path if it already existed.
    Version 8 requires passing force=True.

    Ref: https://github.com/sphinx-doc/sphinx/pull/12647
    """
    from inspect import signature
    if "force" in signature(copyfile).parameters:
        return {"force": True}
    else:
        return {}


def html_collect_pages(app):
    """This event handler is abused to copy local files around."""
    files = set()
    for file_list in app.env.nbsphinx_files.values():
        files.update(file_list)
    for file in status_iterator(files, 'copying linked files... ',
                                sphinx.util.console.brown, len(files)):
        target = os.path.join(app.builder.outdir, file)
        sphinx.util.ensuredir(os.path.dirname(target))
        try:
            sphinx.util.copyfile(
                os.path.join(app.env.srcdir, file),
                target,
                **backwards_compat_overwrite())
        except OSError as err:
            logger.warning(
                'Cannot copy local file %r: %s', file, err,
                type='nbsphinx', subtype='localfile')
    notebooks = app.env.nbsphinx_notebooks.values()
    for notebook in status_iterator(
            notebooks, 'copying notebooks ... ',
            'brown', len(notebooks)):
        sphinx.util.copyfile(
            os.path.join(app.env.nbsphinx_auxdir, notebook),
            os.path.join(app.builder.outdir, notebook),
            **backwards_compat_overwrite())
    return []  # No new HTML pages are created

def env_updated(app, env):
    widgets_path = app.config.nbsphinx_widgets_path
    if widgets_path is None:
        if env.nbsphinx_widgets:
            try:
                from ipywidgets.embed import DEFAULT_EMBED_REQUIREJS_URL
            except ImportError:
                logger.warning(
                    'nbsphinx_widgets_path not given '
                    'and ipywidgets module unavailable',
                    type='nbsphinx', subtype='ipywidgets')
            else:
                widgets_path = DEFAULT_EMBED_REQUIREJS_URL
        else:
            widgets_path = ''

    if widgets_path:
        app.add_js_file(widgets_path, **app.config.nbsphinx_widgets_options)


def doctree_resolved(app, doctree, fromdocname):
    base = sphinx.util.osutil.relative_uri(
        app.builder.get_target_uri(fromdocname), '')

    # Replace GalleryToc with toctree + GalleryNode
    for node in doctree.findall(GalleryToc):
        toctree_wrapper, = node
        if (len(toctree_wrapper) != 1 or not isinstance(
                toctree_wrapper[0],
                (sphinx.addnodes.toctree, DummyTocTree))):
            # This happens for LaTeX output
            node.replace_self(node.children)
            continue
        toctree, = toctree_wrapper
        entries = []
        for title, doc in toctree['entries']:
            if doc in toctree['includefiles']:
                if title is None:
                    title = app.env.titles[doc].astext()
                uri = app.builder.get_relative_uri(fromdocname, doc)

                # NB: This is how Sphinx implements the "html_sidebars"
                #     config value in StandaloneHTMLBuilder.add_sidebars()

                def has_wildcard(pattern):
                    return any(char in pattern for char in '*?[')

                matched = None
                conf_py_thumbnail = None
                conf_py_thumbnails = app.env.config.nbsphinx_thumbnails.items()
                for pattern, candidate in conf_py_thumbnails:
                    if patmatch(doc, pattern):
                        if matched:
                            if has_wildcard(pattern):
                                # warn if both patterns contain wildcards
                                if has_wildcard(matched):
                                    logger.warning(
                                        'page %s matches two patterns in '
                                        'nbsphinx_thumbnails: %r and %r',
                                        doc, matched, pattern,
                                        type='nbsphinx', subtype='thumbnail')
                                # else the already matched pattern is more
                                # specific than the present one, because it
                                # contains no wildcard
                                continue
                        matched = pattern
                        conf_py_thumbnail = candidate

                thumbnail = app.env.nbsphinx_thumbnails.get(doc, {})
                # NB: "None" is used as marker for implicit thumbnail:
                tooltip = thumbnail.get('tooltip', '')
                filename = thumbnail.get('filename', '')

                # thumbnail priority: broken, explicit in notebook,
                #     from conf.py, implicit in notebook, default
                if filename is _BROKEN_THUMBNAIL:
                    filename = os.path.join(
                        base, '_static', 'nbsphinx-broken-thumbnail.svg')
                elif filename and tooltip is not None:
                    # thumbnail from tagged cell or metadata
                    filename = os.path.join(
                        base, app.builder.imagedir, filename)
                elif conf_py_thumbnail:
                    # NB: Settings from conf.py can be overwritten in notebook
                    filename = os.path.join(base, conf_py_thumbnail)
                elif filename:
                    # implicit thumbnail from the last image in the notebook
                    assert tooltip is None
                    tooltip = ''
                    filename = os.path.join(
                        base, app.builder.imagedir, filename)
                else:
                    filename = os.path.join(
                        base, '_static', 'nbsphinx-no-thumbnail.svg')
                entries.append((title, uri, filename, tooltip))
            else:
                logger.warning(
                    'External links are not supported in gallery: %s', doc,
                    location=fromdocname, type='nbsphinx', subtype='gallery')
        gallery = GalleryNode()
        gallery['entries'] = entries
        if isinstance(toctree, DummyTocTree):
            # NbLinkGallery, toctree is only used to get optional caption.

            # Copied from sphinx/environment/adapters/toctree.py:
            newnode = sphinx.addnodes.compact_paragraph('', '')
            caption = toctree.attributes.get('caption')
            if caption:
                caption_node = docutils.nodes.title(
                    caption, '', *[docutils.nodes.Text(caption)])
                caption_node.line = toctree.line
                caption_node.source = toctree.source
                #caption_node.rawsource = toctree['rawcaption']
                if hasattr(toctree, 'uid'):
                    # move uid to caption_node to translate it
                    caption_node.uid = toctree.uid  # type: ignore
                    del toctree.uid
                newnode += caption_node
            newnode['toctree'] = True

            toctree_wrapper[0] = newnode
        else:
            # NbGallery
            toctree['nbsphinx_gallery'] = True
            # NB: Further processing happens in patched_toctree_resolve()
        toctree_wrapper += gallery
        node.replace_self(toctree_wrapper)
        doctree['nbsphinx_gallery_css'] = True


def depart_codearea_html(self, node):
    """Add empty lines before and after the code."""
    text = self.body[-1]
    text = text.replace('<pre>',
                        '<pre>' + '<br/>' * node.get('empty-lines-before', 0))
    text = text.replace('</pre>',
                        '<br/>' * node.get('empty-lines-after', 0) + '</pre>')
    self.body[-1] = text


def visit_codearea_latex(self, node):
    self.pushbody([])  # See popbody() below


def depart_codearea_latex(self, node):
    """Some changes to code blocks.

    * Change frame color and background color
    * Add empty lines before and after the code
    * Add prompt

    """
    lines = ''.join(self.popbody()).strip('\n').split('\n')
    out = []
    out.append('')
    out.append('{')  # Start a scope for colors
    if 'nbinput' in node.parent['classes']:
        promptcolor = 'nbsphinxin'
    else:
        out.append(r"""
\kern-\sphinxverbatimsmallskipamount\kern-\baselineskip
\kern+\FrameHeightAdjust\kern-\fboxrule
\vspace{\nbsphinxcodecellspacing}
""")
        promptcolor = 'nbsphinxout'
        if node['stderr']:
            out.append(r'\sphinxsetup{VerbatimColor={named}{nbsphinx-stderr}}')
        else:
            out.append(r'\sphinxsetup{VerbatimColor={named}{white}}')
    if lines[0].startswith(r'\fvset{'):  # Sphinx >= 1.6.6 and < 1.8.3
        out.append(lines[0])
        del lines[0]
    # Sphinx 4.1.0 added "sphinxuseclass" environments around "sphinxVerbatim"
    for begin_verbatim, line in enumerate(lines):
        if line.startswith(r'\begin{sphinxVerbatim}'):
            break
    else:
        assert False
    for end_verbatim, line in enumerate(reversed(lines)):
        if line == r'\end{sphinxVerbatim}':
            break
    else:
        assert False
    out.extend(lines[:begin_verbatim + 1])
    code_lines = (
        [''] * node.get('empty-lines-before', 0) +
        lines[begin_verbatim + 1:-end_verbatim - 1] +
        [''] * node.get('empty-lines-after', 0)
    )
    prompt = node['prompt']
    if prompt:
        prompt = nbconvert.filters.latex.escape_latex(prompt)
        prefix = r'\llap{\color{' + promptcolor + '}' + prompt + \
            r'\,\hspace{\fboxrule}\hspace{\fboxsep}}'
        assert code_lines
        code_lines[0] = prefix + code_lines[0]
    out.extend(code_lines)
    out.extend(lines[-end_verbatim - 1:])
    out.append('}')  # End of scope for colors
    out.append('')
    self.body.append('\n'.join(out))


def visit_fancyoutput_latex(self, node):
    out = r"""
\hrule height -\fboxrule\relax
\vspace{\nbsphinxcodecellspacing}
"""
    prompt = node['prompt']
    if prompt:
        prompt = nbconvert.filters.latex.escape_latex(prompt)
        out += r"""
\savebox\nbsphinxpromptbox[0pt][r]{\color{nbsphinxout}\Verb|\strut{%s}\,|}
""" % (prompt,)
    else:
        out += r"""
\makeatletter\setbox\nbsphinxpromptbox\box\voidb@x\makeatother
"""
    out += r"""
\begin{nbsphinxfancyoutput}
"""
    self.body.append(out)


def depart_fancyoutput_latex(self, node):
    self.body.append('\n\\end{nbsphinxfancyoutput}\n')


def visit_admonition_html(self, node):
    self.body.append(self.starttag(node, 'div'))
    if len(node.children) >= 2:
        node[0]['classes'].append('admonition-title')
        html_theme = self.settings.env.config.html_theme
        if html_theme in ('sphinx_rtd_theme', 'julia', 'dask_sphinx_theme'):
            node.children[0]['classes'].extend(['fa', 'fa-exclamation-circle'])


def depart_admonition_html(self, node):
    self.body.append('</div>\n')


def visit_admonition_latex(self, node):
    kind = node['classes'][1]
    if len(node.children) >= 2 and isinstance(
            node.children[0], docutils.nodes.paragraph):
        title = node.children[0].astext()
        del node.children[0]
        self.body.append(
            '\n\\begin{sphinxadmonition}{' + kind + '}{' + title + '}\\par')
    else:
        # See http://tex.stackexchange.com/q/305898/:
        self.body.append(
            '\n\\begin{sphinxadmonition}{' + kind + '}{}\\unskip')


def depart_admonition_latex(self, node):
    self.body.append('\\end{sphinxadmonition}\n')


def visit_admonition_text(self, node):
    self.new_state(0)


def depart_admonition_text(self, node):
    self.end_state()


def depart_gallery_html(self, node):
    self.body.append('<div class="nbsphinx-gallery">\n')
    for title, uri, filename, tooltip in node['entries']:
        if tooltip:
            tooltip = ' title="{}"'.format(html.escape(tooltip))
        self.body.append("""\
<a class="reference internal" href="{uri}"{tooltip}>
  <div><img alt="" src="{filename}"></div>
  <div>{title}</div>
</a>
""".format(
            uri=html.escape(uri),
            title=html.escape(title),
            tooltip=tooltip,
            filename=html.escape(filename),
        ))
    self.body.append('</div>\n')


def do_nothing(self, node):
    pass


def setup(app):
    """Initialize Sphinx extension."""
    app.add_source_parser(NotebookParser)

    app.add_config_value('nbsphinx_execute', 'auto', rebuild='env')
    app.add_config_value('nbsphinx_kernel_name', '', rebuild='env')
    app.add_config_value('nbsphinx_execute_arguments', [], rebuild='env')
    app.add_config_value('nbsphinx_allow_errors', False, rebuild='')
    app.add_config_value('nbsphinx_timeout', None, rebuild='')
    app.add_config_value('nbsphinx_codecell_lexer', 'none', rebuild='env')
    app.add_config_value('nbsphinx_prompt_width', '4.5ex', rebuild='html')
    app.add_config_value('nbsphinx_responsive_width', '540px', rebuild='html')
    app.add_config_value('nbsphinx_prolog', None, rebuild='env')
    app.add_config_value('nbsphinx_epilog', None, rebuild='env')
    app.add_config_value('nbsphinx_input_prompt', '[%s]:', rebuild='env')
    app.add_config_value('nbsphinx_output_prompt', '[%s]:', rebuild='env')
    app.add_config_value('nbsphinx_custom_formats', {}, rebuild='env')
    # Default value is set in config_inited():
    app.add_config_value('nbsphinx_requirejs_path', None, rebuild='html')
    # Default value is set in config_inited():
    app.add_config_value('nbsphinx_requirejs_options', None, rebuild='html')
    # This will be updated in env_updated():
    app.add_config_value('nbsphinx_widgets_path', None, rebuild='html')
    app.add_config_value('nbsphinx_widgets_options', {}, rebuild='html')
    app.add_config_value('nbsphinx_thumbnails', {}, rebuild='html')
    app.add_config_value('nbsphinx_assume_equations', True, rebuild='env')

    app.add_directive('nbinput', NbInput)
    app.add_directive('nboutput', NbOutput)
    app.add_directive('nbinfo', NbInfo)
    app.add_directive('nbwarning', NbWarning)
    app.add_directive('nbgallery', NbGallery)
    app.add_directive('nblinkgallery', NbLinkGallery)
    app.add_node(CodeAreaNode,
                 html=(do_nothing, depart_codearea_html),
                 latex=(visit_codearea_latex, depart_codearea_latex),
                 text=(do_nothing, do_nothing))
    app.add_node(FancyOutputNode,
                 html=(do_nothing, do_nothing),
                 latex=(visit_fancyoutput_latex, depart_fancyoutput_latex),
                 text=(do_nothing, do_nothing))
    app.add_node(AdmonitionNode,
                 html=(visit_admonition_html, depart_admonition_html),
                 latex=(visit_admonition_latex, depart_admonition_latex),
                 text=(visit_admonition_text, depart_admonition_text))
    app.add_node(GalleryNode,
                 html=(do_nothing, depart_gallery_html),
                 latex=(do_nothing, do_nothing),
                 text=(do_nothing, do_nothing))
    app.connect('builder-inited', builder_inited)
    app.connect('config-inited', config_inited)
    app.connect('html-page-context', html_page_context)
    app.connect('html-collect-pages', html_collect_pages)
    app.connect('env-purge-doc', env_purge_doc)
    app.connect('env-updated', env_updated)
    app.connect('doctree-resolved', doctree_resolved)
    app.connect('env-merge-info', env_merge_info)
    app.add_transform(CreateSectionLabels)
    app.add_transform(CreateDomainObjectLabels)
    app.add_transform(RewriteLocalLinks)
    app.add_post_transform(GetSizeFromImages)
    app.add_latex_package('nbsphinx')

    # Make docutils' "code" directive (generated by markdown2rst/pandoc)
    # behave like Sphinx's "code-block",
    # see https://github.com/sphinx-doc/sphinx/issues/2155:
    rst.directives.register_directive('code', sphinx.directives.code.CodeBlock)

    # Monkey-patch Sphinx TocTree adapter
    if hasattr(sphinx.environment.adapters.toctree, '_resolve_toctree'):
        # Since Sphinx 7.2.0
        sphinx.environment.adapters.toctree._resolve_toctree = \
            patched_toctree_resolve
    else:
        sphinx.environment.adapters.toctree.TocTree.resolve = \
            patched_toctree_resolve

    return {
        'version': __version__,
        'parallel_read_safe': True,
        'parallel_write_safe': True,
        'env_version': 4,
    }
