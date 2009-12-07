"""
A special directive for generating a matplotlib plot.

.. warning::

   This is a hacked version of plot_directive.py from Matplotlib.
   It's very much subject to change!


Usage
-----

Can be used like this::

    .. plot:: examples/example.py

    .. plot::

       import matplotlib.pyplot as plt
       plt.plot([1,2,3], [4,5,6])

    .. plot::

       A plotting example:

       >>> import matplotlib.pyplot as plt
       >>> plt.plot([1,2,3], [4,5,6])

The content is interpreted as doctest formatted if it has a line starting
with ``>>>``.

The ``plot`` directive supports the options

    format : {'python', 'doctest'}
        Specify the format of the input

    include-source : bool
        Whether to display the source code. Default can be changed in conf.py
    
and the ``image`` directive options ``alt``, ``height``, ``width``,
``scale``, ``align``, ``class``.

Configuration options
---------------------

The plot directive has the following configuration options:

    plot_include_source
        Default value for the include-source option

    plot_pre_code
        Code that should be executed before each plot.

    plot_basedir
        Base directory, to which plot:: file names are relative to.
        (If None or empty, file names are relative to the directoly where
        the file containing the directive is.)

    plot_formats
        File formats to generate. List of tuples or strings::

            [(suffix, dpi), suffix, ...]

        that determine the file format and the DPI. For entries whose
        DPI was omitted, sensible defaults are chosen.

TODO
----

* Refactor Latex output; now it's plain images, but it would be nice
  to make them appear side-by-side, or in floats.

"""

import sys, os, glob, shutil, imp, warnings, cStringIO, re, textwrap, traceback
import sphinx

import warnings
warnings.warn("A plot_directive module is also available under "
              "matplotlib.sphinxext; expect this numpydoc.plot_directive "
              "module to be deprecated after relevant features have been "
              "integrated there.",
              FutureWarning, stacklevel=2)


#------------------------------------------------------------------------------
# Registration hook
#------------------------------------------------------------------------------

def setup(app):
    setup.app = app
    setup.config = app.config
    setup.confdir = app.confdir
    
    app.add_config_value('plot_pre_code', '', True)
    app.add_config_value('plot_include_source', False, True)
    app.add_config_value('plot_formats', ['png', 'hires.png', 'pdf'], True)
    app.add_config_value('plot_basedir', None, True)

    app.add_directive('plot', plot_directive, True, (0, 1, False),
                      **plot_directive_options)

#------------------------------------------------------------------------------
# plot:: directive
#------------------------------------------------------------------------------
from docutils.parsers.rst import directives
from docutils import nodes

def plot_directive(name, arguments, options, content, lineno,
                   content_offset, block_text, state, state_machine):
    return run(arguments, content, options, state_machine, state, lineno)
plot_directive.__doc__ = __doc__

def _option_boolean(arg):
    if not arg or not arg.strip():
        # no argument given, assume used as a flag
        return True
    elif arg.strip().lower() in ('no', '0', 'false'):
        return False
    elif arg.strip().lower() in ('yes', '1', 'true'):
        return True
    else:
        raise ValueError('"%s" unknown boolean' % arg)

def _option_format(arg):
    return directives.choice(arg, ('python', 'lisp'))

def _option_align(arg):
    return directives.choice(arg, ("top", "middle", "bottom", "left", "center",
                                   "right"))

plot_directive_options = {'alt': directives.unchanged,
                          'height': directives.length_or_unitless,
                          'width': directives.length_or_percentage_or_unitless,
                          'scale': directives.nonnegative_int,
                          'align': _option_align,
                          'class': directives.class_option,
                          'include-source': _option_boolean,
                          'format': _option_format,
                          }

#------------------------------------------------------------------------------
# Generating output
#------------------------------------------------------------------------------

from docutils import nodes, utils

try:
    # Sphinx depends on either Jinja or Jinja2
    import jinja2
    def format_template(template, **kw):
        return jinja2.Template(template).render(**kw)
except ImportError:
    import jinja
    def format_template(template, **kw):
        return jinja.from_string(template, **kw)

TEMPLATE = """
{{ source_code }}

{{ only_html }}

   {% if source_code %}
   (`Source code <{{ source_link }}>`__)

   .. admonition:: Output
      :class: plot-output

   {% endif %}

       {% for img in images %}
       .. figure:: {{ build_dir }}/{{ img.basename }}.png
          {%- for option in options %}
          {{ option }}
          {% endfor %}
    
          (
          {%- if not source_code -%}
            `Source code <{{source_link}}>`__
            {%- for fmt in img.formats -%} 
            , `{{ fmt }} <{{ dest_dir }}/{{ img.basename }}.{{ fmt }}>`__
            {%- endfor -%}
          {%- else -%}
            {%- for fmt in img.formats -%} 
            {%- if not loop.first -%}, {% endif -%}
            `{{ fmt }} <{{ dest_dir }}/{{ img.basename }}.{{ fmt }}>`__
            {%- endfor -%}
          {%- endif -%}
          )
       {% endfor %}

{{ only_latex }}

   {% for img in images %}
   .. image:: {{ build_dir }}/{{ img.basename }}.pdf
   {% endfor %}

"""

class ImageFile(object):
    def __init__(self, basename, dirname):
        self.basename = basename
        self.dirname = dirname
        self.formats = []

    def filename(self, format):
        return os.path.join(self.dirname, "%s.%s" % (self.basename, format))

    def filenames(self):
        return [self.filename(fmt) for fmt in self.formats]

def run(arguments, content, options, state_machine, state, lineno):
    if arguments and content:
        raise RuntimeError("plot:: directive can't have both args and content")

    document = state_machine.document
    config = document.settings.env.config

    options.setdefault('include-source', config.plot_include_source)

    # determine input
    rst_file = document.attributes['source']
    rst_dir = os.path.dirname(rst_file)

    if arguments:
        if not config.plot_basedir:
            source_file_name = os.path.join(rst_dir,
                                            directives.uri(arguments[0]))
        else:
            source_file_name = os.path.join(setup.confdir, config.plot_basedir,
                                            directives.uri(arguments[0]))
        code = open(source_file_name, 'r').read()
        output_base = os.path.basename(source_file_name)
    else:
        source_file_name = rst_file
        code = textwrap.dedent("\n".join(map(str, content)))
        counter = document.attributes.get('_plot_counter', 0) + 1
        document.attributes['_plot_counter'] = counter
        base, ext = os.path.splitext(os.path.basename(source_file_name))
        output_base = '%s-%d.py' % (base, counter)

    base, source_ext = os.path.splitext(output_base)
    if source_ext in ('.py', '.rst', '.txt'):
        output_base = base
    else:
        source_ext = ''

    # ensure that LaTeX includegraphics doesn't choke in foo.bar.pdf filenames
    output_base = output_base.replace('.', '-')

    # is it in doctest format?
    is_doctest = contains_doctest(code)
    if options.has_key('format'):
        if options['format'] == 'python':
            is_doctest = False
        else:
            is_doctest = True

    # determine output directory name fragment
    source_rel_name = relpath(source_file_name, setup.confdir)
    source_rel_dir = os.path.dirname(source_rel_name)
    while source_rel_dir.startswith(os.path.sep):
        source_rel_dir = source_rel_dir[1:]

    # build_dir: where to place output files (temporarily)
    build_dir = os.path.join(os.path.dirname(setup.app.doctreedir),
                             'plot_directive',
                             source_rel_dir)
    if not os.path.exists(build_dir):
        os.makedirs(build_dir)

    # output_dir: final location in the builder's directory
    dest_dir = os.path.abspath(os.path.join(setup.app.builder.outdir,
                                            source_rel_dir))

    # how to link to files from the RST file
    dest_dir_link = os.path.join(relpath(setup.confdir, rst_dir),
                                 source_rel_dir).replace(os.path.sep, '/')
    build_dir_link = relpath(build_dir, rst_dir).replace(os.path.sep, '/')
    source_link = dest_dir_link + '/' + output_base + source_ext

    # make figures
    try:
        images = makefig(code, source_file_name, build_dir, output_base,
                         config)
    except PlotError, err:
        reporter = state.memo.reporter
        sm = reporter.system_message(
            3, "Exception occurred in plotting %s: %s" % (output_base, err),
            line=lineno)
        return [sm]

    # generate output restructuredtext
    if options['include-source']:
        if is_doctest:
            lines = ['']
            lines += [row.rstrip() for row in code.split('\n')]
        else:
            lines = ['.. code-block:: python', '']
            lines += ['    %s' % row.rstrip() for row in code.split('\n')]
        source_code = "\n".join(lines)
    else:
        source_code = ""

    opts = [':%s: %s' % (key, val) for key, val in options.items()
            if key in ('alt', 'height', 'width', 'scale', 'align', 'class')]

    if sphinx.__version__ >= "0.6":
        only_html = ".. only:: html"
        only_latex = ".. only:: latex"
    else:
        only_html = ".. htmlonly::"
        only_latex = ".. latexonly::"

    result = format_template(
        TEMPLATE,
        dest_dir=dest_dir_link,
        build_dir=build_dir_link,
        source_link=source_link,
        only_html=only_html,
        only_latex=only_latex,
        options=opts,
        images=images,
        source_code=source_code)

    lines = result.split("\n")
    if len(lines):
        state_machine.insert_input(
            lines, state_machine.input_lines.source(0))

    # copy image files to builder's output directory
    if not os.path.exists(dest_dir):
        os.makedirs(dest_dir)

    for img in images:
        for fn in img.filenames():
            shutil.copyfile(fn, os.path.join(dest_dir, os.path.basename(fn)))

    # copy script (if necessary)
    if source_file_name == rst_file:
        target_name = os.path.join(dest_dir, output_base + source_ext)
        f = open(target_name, 'w')
        f.write(unescape_doctest(code))
        f.close()

    return []


#------------------------------------------------------------------------------
# Run code and capture figures
#------------------------------------------------------------------------------

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as image
from matplotlib import _pylab_helpers

import exceptions

def contains_doctest(text):
    try:
        # check if it's valid Python as-is
        compile(text, '<string>', 'exec')
        return False
    except SyntaxError:
        pass
    r = re.compile(r'^\s*>>>', re.M)
    m = r.search(text)
    return bool(m)

def unescape_doctest(text):
    """
    Extract code from a piece of text, which contains either Python code
    or doctests.

    """
    if not contains_doctest(text):
        return text

    code = ""
    for line in text.split("\n"):
        m = re.match(r'^\s*(>>>|\.\.\.) (.*)$', line)
        if m:
            code += m.group(2) + "\n"
        elif line.strip():
            code += "# " + line.strip() + "\n"
        else:
            code += "\n"
    return code

class PlotError(RuntimeError):
    pass

def run_code(code, code_path):
    # Change the working directory to the directory of the example, so
    # it can get at its data files, if any.
    pwd = os.getcwd()
    old_sys_path = list(sys.path)
    if code_path is not None:
        dirname = os.path.abspath(os.path.dirname(code_path))
        os.chdir(dirname)
        sys.path.insert(0, dirname)

    # Redirect stdout
    stdout = sys.stdout
    sys.stdout = cStringIO.StringIO()

    # Reset sys.argv
    old_sys_argv = sys.argv
    sys.argv = [code_path]
    
    try:
        try:
            code = unescape_doctest(code)
            ns = {}
            exec setup.config.plot_pre_code in ns
            exec code in ns
        except (Exception, SystemExit), err:
            raise PlotError(traceback.format_exc())
    finally:
        os.chdir(pwd)
        sys.argv = old_sys_argv
        sys.path[:] = old_sys_path
        sys.stdout = stdout
    return ns


#------------------------------------------------------------------------------
# Generating figures
#------------------------------------------------------------------------------

def out_of_date(original, derived):
    """
    Returns True if derivative is out-of-date wrt original,
    both of which are full file paths.
    """
    return (not os.path.exists(derived)
            or os.stat(derived).st_mtime < os.stat(original).st_mtime)


def makefig(code, code_path, output_dir, output_base, config):
    """
    Run a pyplot script *code* and save the images under *output_dir*
    with file names derived from *output_base*

    """

    # -- Parse format list
    default_dpi = {'png': 80, 'hires.png': 200, 'pdf': 50}
    formats = []
    for fmt in config.plot_formats:
        if isinstance(fmt, str):
            formats.append((fmt, default_dpi.get(fmt, 80)))
        elif type(fmt) in (tuple, list) and len(fmt)==2:
            formats.append((str(fmt[0]), int(fmt[1])))
        else:
            raise PlotError('invalid image format "%r" in plot_formats' % fmt)

    # -- Try to determine if all images already exist

    # Look for single-figure output files first
    all_exists = True
    img = ImageFile(output_base, output_dir)
    for format, dpi in formats:
        if out_of_date(code_path, img.filename(format)):
            all_exists = False
            break
        img.formats.append(format)

    if all_exists:
        return [img]

    # Then look for multi-figure output files
    images = []
    all_exists = True
    for i in xrange(1000):
        img = ImageFile('%s_%02d' % (output_base, i), output_dir)
        for format, dpi in formats:
            if out_of_date(code_path, img.filename(format)):
                all_exists = False
                break
            img.formats.append(format)
        
        # assume that if we have one, we have them all
        if not all_exists:
            all_exists = (i > 0)
            break
        images.append(img)

    if all_exists:
        return images

    # -- We didn't find the files, so build them

    # Clear between runs
    plt.close('all')

    # Run code
    run_code(code, code_path)

    # Collect images
    images = []

    fig_managers = _pylab_helpers.Gcf.get_all_fig_managers()
    for i, figman in enumerate(fig_managers):
        if len(fig_managers) == 1:
            img = ImageFile(output_base, output_dir)
        else:
            img = ImageFile("%s_%02d" % (output_base, i), output_dir)
        images.append(img)
        for format, dpi in formats:
            try:
                figman.canvas.figure.savefig(img.filename(format), dpi=dpi)
            except exceptions.BaseException, err:
                raise PlotError(traceback.format_exc())
            img.formats.append(format)

    return images


#------------------------------------------------------------------------------
# Relative pathnames
#------------------------------------------------------------------------------

try:
    from os.path import relpath
except ImportError:
    def relpath(target, base=os.curdir):
        """
        Return a relative path to the target from either the current
        dir or an optional base dir.  Base can be a directory
        specified either as absolute or relative to current dir.
        """

        if not os.path.exists(target):
            raise OSError, 'Target does not exist: '+target

        if not os.path.isdir(base):
            raise OSError, 'Base is not a directory or does not exist: '+base

        base_list = (os.path.abspath(base)).split(os.sep)
        target_list = (os.path.abspath(target)).split(os.sep)

        # On the windows platform the target may be on a completely
        # different drive from the base.
        if os.name in ['nt','dos','os2'] and base_list[0] <> target_list[0]:
            raise OSError, 'Target is on a different drive to base. Target: '+target_list[0].upper()+', base: '+base_list[0].upper()

        # Starting from the filepath root, work out how much of the
        # filepath is shared by base and target.
        for i in range(min(len(base_list), len(target_list))):
            if base_list[i] <> target_list[i]: break
        else:
            # If we broke out of the loop, i is pointing to the first
            # differing path elements.  If we didn't break out of the
            # loop, i is pointing to identical path elements.
            # Increment i so that in all cases it points to the first
            # differing path elements.
            i+=1

        rel_list = [os.pardir] * (len(base_list)-i) + target_list[i:]
        return os.path.join(*rel_list)
