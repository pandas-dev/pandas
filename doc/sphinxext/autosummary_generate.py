#!/usr/bin/env python
r"""
autosummary_generate.py OPTIONS FILES

Generate automatic RST source files for items referred to in
autosummary:: directives.

Each generated RST file contains a single auto*:: directive which
extracts the docstring of the referred item.

Example Makefile rule::

    generate:
            ./ext/autosummary_generate.py -o source/generated source/*.rst

"""
import glob, re, inspect, os, optparse, pydoc
from autosummary import import_by_name

try:
    from phantom_import import import_phantom_module
except ImportError:
    import_phantom_module = lambda x: x

def main():
    p = optparse.OptionParser(__doc__.strip())
    p.add_option("-p", "--phantom", action="store", type="string",
                 dest="phantom", default=None,
                 help="Phantom import modules from a file")
    p.add_option("-o", "--output-dir", action="store", type="string",
                 dest="output_dir", default=None,
                 help=("Write all output files to the given directory (instead "
                       "of writing them as specified in the autosummary:: "
                       "directives)"))
    options, args = p.parse_args()

    if len(args) == 0:
        p.error("wrong number of arguments")

    if options.phantom and os.path.isfile(options.phantom):
        import_phantom_module(options.phantom)

    # read
    names = {}
    for name, loc in get_documented(args).items():
        for (filename, sec_title, keyword, toctree) in loc:
            if toctree is not None:
                path = os.path.join(os.path.dirname(filename), toctree)
                names[name] = os.path.abspath(path)

    # write
    for name, path in sorted(names.items()):
        if options.output_dir is not None:
            path = options.output_dir
        
        if not os.path.isdir(path):
            os.makedirs(path)

        try:
            obj, name = import_by_name(name)
        except ImportError, e:
            print "Failed to import '%s': %s" % (name, e)
            continue

        fn = os.path.join(path, '%s.rst' % name)

        if os.path.exists(fn):
            # skip
            continue

        f = open(fn, 'w')

        try:
            f.write('%s\n%s\n\n' % (name, '='*len(name)))

            if inspect.isclass(obj):
                if issubclass(obj, Exception):
                    f.write(format_modulemember(name, 'autoexception'))
                else:
                    f.write(format_modulemember(name, 'autoclass'))
            elif inspect.ismodule(obj):
                f.write(format_modulemember(name, 'automodule'))
            elif inspect.ismethod(obj) or inspect.ismethoddescriptor(obj):
                f.write(format_classmember(name, 'automethod'))
            elif callable(obj):
                f.write(format_modulemember(name, 'autofunction'))
            elif hasattr(obj, '__get__'):
                f.write(format_classmember(name, 'autoattribute'))
            else:
                f.write(format_modulemember(name, 'autofunction'))
        finally:
            f.close()

def format_modulemember(name, directive):
    parts = name.split('.')
    mod, name = '.'.join(parts[:-1]), parts[-1]
    return ".. currentmodule:: %s\n\n.. %s:: %s\n" % (mod, directive, name)

def format_classmember(name, directive):
    parts = name.split('.')
    mod, name = '.'.join(parts[:-2]), '.'.join(parts[-2:])
    return ".. currentmodule:: %s\n\n.. %s:: %s\n" % (mod, directive, name)

def get_documented(filenames):
    """
    Find out what items are documented in source/*.rst
    See `get_documented_in_lines`.

    """
    documented = {}
    for filename in filenames:
        f = open(filename, 'r')
        lines = f.read().splitlines()
        documented.update(get_documented_in_lines(lines, filename=filename))
        f.close()
    return documented

def get_documented_in_docstring(name, module=None, filename=None):
    """
    Find out what items are documented in the given object's docstring.
    See `get_documented_in_lines`.
    
    """
    try:
        obj, real_name = import_by_name(name)
        lines = pydoc.getdoc(obj).splitlines()
        return get_documented_in_lines(lines, module=name, filename=filename)
    except AttributeError:
        pass
    except ImportError, e:
        print "Failed to import '%s': %s" % (name, e)
    return {}

def get_documented_in_lines(lines, module=None, filename=None):
    """
    Find out what items are documented in the given lines
    
    Returns
    -------
    documented : dict of list of (filename, title, keyword, toctree)
        Dictionary whose keys are documented names of objects.
        The value is a list of locations where the object was documented.
        Each location is a tuple of filename, the current section title,
        the name of the directive, and the value of the :toctree: argument
        (if present) of the directive.

    """
    title_underline_re = re.compile("^[-=*_^#]{3,}\s*$")
    autodoc_re = re.compile(".. auto(function|method|attribute|class|exception|module)::\s*([A-Za-z0-9_.]+)\s*$")
    autosummary_re = re.compile(r'^\.\.\s+autosummary::\s*')
    module_re = re.compile(r'^\.\.\s+(current)?module::\s*([a-zA-Z0-9_.]+)\s*$')
    autosummary_item_re = re.compile(r'^\s+([_a-zA-Z][a-zA-Z0-9_.]*)\s*.*?')
    toctree_arg_re = re.compile(r'^\s+:toctree:\s*(.*?)\s*$')
    
    documented = {}
   
    current_title = []
    last_line = None
    toctree = None
    current_module = module
    in_autosummary = False
    
    for line in lines:
        try:
            if in_autosummary:
                m = toctree_arg_re.match(line)
                if m:
                    toctree = m.group(1)
                    continue

                if line.strip().startswith(':'):
                    continue # skip options

                m = autosummary_item_re.match(line)
                if m:
                    name = m.group(1).strip()
                    if current_module and not name.startswith(current_module + '.'):
                        name = "%s.%s" % (current_module, name)
                    documented.setdefault(name, []).append(
                        (filename, current_title, 'autosummary', toctree))
                    continue
                if line.strip() == '':
                    continue
                in_autosummary = False

            m = autosummary_re.match(line)
            if m:
                in_autosummary = True
                continue

            m = autodoc_re.search(line)
            if m:
                name = m.group(2).strip()
                if m.group(1) == "module":
                    current_module = name
                    documented.update(get_documented_in_docstring(
                        name, filename=filename))
                elif current_module and not name.startswith(current_module+'.'):
                    name = "%s.%s" % (current_module, name)
                documented.setdefault(name, []).append(
                    (filename, current_title, "auto" + m.group(1), None))
                continue

            m = title_underline_re.match(line)
            if m and last_line:
                current_title = last_line.strip()
                continue

            m = module_re.match(line)
            if m:
                current_module = m.group(2)
                continue
        finally:
            last_line = line

    return documented

if __name__ == "__main__":
    main()
