"""
This file contains `__main__` so that it can be run as a commandline tool.

This file contains functions to inspect Numba's support for a given Python
module or a Python package.
"""

import argparse
import pkgutil
import warnings
import types as pytypes

from numba.core import errors
from numba._version import get_versions
from numba.core.registry import cpu_target
from numba.tests.support import captured_stdout


def _get_commit():
    full = get_versions()['full-revisionid']
    if not full:
        warnings.warn(
            "Cannot find git commit hash. Source links could be inaccurate.",
            category=errors.NumbaWarning,
        )
        return 'main'
    return full


commit = _get_commit()
github_url = 'https://github.com/numba/numba/blob/{commit}/{path}#L{firstline}-L{lastline}'    # noqa: E501


def inspect_function(function, target=None):
    """Return information about the support of a function.

    Returns
    -------
    info : dict
        Defined keys:
        - "numba_type": str or None
            The numba type object of the function if supported.
        - "explained": str
            A textual description of the support.
        - "source_infos": dict
            A dictionary containing the source location of each definition.
    """
    target = target or cpu_target
    tyct = target.typing_context
    # Make sure we have loaded all extensions
    tyct.refresh()
    target.target_context.refresh()

    info = {}
    # Try getting the function type
    source_infos = {}
    try:
        nbty = tyct.resolve_value_type(function)
    except ValueError:
        nbty = None
        explained = 'not supported'
    else:
        # Make a longer explanation of the type
        explained = tyct.explain_function_type(nbty)
        for temp in nbty.templates:
            try:
                source_infos[temp] = temp.get_source_info()
            except AttributeError:
                source_infos[temp] = None

    info['numba_type'] = nbty
    info['explained'] = explained
    info['source_infos'] = source_infos
    return info


def inspect_module(module, target=None, alias=None):
    """Inspect a module object and yielding results from `inspect_function()`
    for each function object in the module.
    """
    alias = {} if alias is None else alias
    # Walk the module
    for name in dir(module):
        if name.startswith('_'):
            # Skip
            continue
        obj = getattr(module, name)
        supported_types = (pytypes.FunctionType, pytypes.BuiltinFunctionType)

        if not isinstance(obj, supported_types):
            # Skip if it's not a function
            continue

        info = dict(module=module, name=name, obj=obj)
        if obj in alias:
            info['alias'] = alias[obj]
        else:
            alias[obj] = "{module}.{name}".format(module=module.__name__,
                                                  name=name)
        info.update(inspect_function(obj, target=target))
        yield info


class _Stat(object):
    """For gathering simple statistic of (un)supported functions"""
    def __init__(self):
        self.supported = 0
        self.unsupported = 0

    @property
    def total(self):
        total = self.supported + self.unsupported
        return total

    @property
    def ratio(self):
        ratio = self.supported / self.total * 100
        return ratio

    def describe(self):
        if self.total == 0:
            return "empty"
        return "supported = {supported} / {total} = {ratio:.2f}%".format(
            supported=self.supported,
            total=self.total,
            ratio=self.ratio,
        )

    def __repr__(self):
        return "{clsname}({describe})".format(
            clsname=self.__class__.__name__,
            describe=self.describe(),
        )


def filter_private_module(module_components):
    return not any(x.startswith('_') for x in module_components)


def filter_tests_module(module_components):
    return not any(x == 'tests' for x in module_components)


_default_module_filters = (
    filter_private_module,
    filter_tests_module,
)


def list_modules_in_package(package, module_filters=_default_module_filters):
    """Yield all modules in a given package.

    Recursively walks the package tree.
    """
    onerror_ignore = lambda _: None

    prefix = package.__name__ + "."
    package_walker = pkgutil.walk_packages(
        package.__path__,
        prefix,
        onerror=onerror_ignore,
    )

    def check_filter(modname):
        module_components = modname.split('.')
        return any(not filter_fn(module_components)
                   for filter_fn in module_filters)

    modname = package.__name__
    if not check_filter(modname):
        yield package

    for pkginfo in package_walker:
        modname = pkginfo[1]
        if check_filter(modname):
            continue
        # In case importing of the module print to stdout
        with captured_stdout():
            try:
                # Import the module
                mod = __import__(modname)
            except Exception:
                continue

            # Extract the module
            for part in modname.split('.')[1:]:
                try:
                    mod = getattr(mod, part)
                except AttributeError:
                    # Suppress error in getting the attribute
                    mod = None
                    break

        # Ignore if mod is not a module
        if not isinstance(mod, pytypes.ModuleType):
            # Skip non-module
            continue

        yield mod


class Formatter(object):
    """Base class for formatters.
    """
    def __init__(self, fileobj):
        self._fileobj = fileobj

    def print(self, *args, **kwargs):
        kwargs.setdefault('file', self._fileobj)
        print(*args, **kwargs)


class HTMLFormatter(Formatter):
    """Formatter that outputs HTML
    """

    def escape(self, text):
        import html
        return html.escape(text)

    def title(self, text):
        self.print('<h1>', text, '</h2>')

    def begin_module_section(self, modname):
        self.print('<h2>', modname, '</h2>')
        self.print('<ul>')

    def end_module_section(self):
        self.print('</ul>')

    def write_supported_item(self, modname, itemname, typename, explained,
                             sources, alias):
        self.print('<li>')
        self.print('{}.<b>{}</b>'.format(
            modname,
            itemname,
        ))
        self.print(': <b>{}</b>'.format(typename))
        self.print('<div><pre>', explained, '</pre></div>')

        self.print("<ul>")
        for tcls, source in sources.items():
            if source:
                self.print("<li>")
                impl = source['name']
                sig = source['sig']
                filename = source['filename']
                lines = source['lines']
                self.print(
                    "<p>defined by <b>{}</b>{} at {}:{}-{}</p>".format(
                        self.escape(impl), self.escape(sig),
                        self.escape(filename), lines[0], lines[1],
                    ),
                )
                self.print('<p>{}</p>'.format(
                    self.escape(source['docstring'] or '')
                ))
            else:
                self.print("<li>{}".format(self.escape(str(tcls))))
            self.print("</li>")
        self.print("</ul>")
        self.print('</li>')

    def write_unsupported_item(self, modname, itemname):
        self.print('<li>')
        self.print('{}.<b>{}</b>: UNSUPPORTED'.format(
            modname,
            itemname,
        ))
        self.print('</li>')

    def write_statistic(self, stats):
        self.print('<p>{}</p>'.format(stats.describe()))


class ReSTFormatter(Formatter):
    """Formatter that output ReSTructured text format for Sphinx docs.
    """
    def escape(self, text):
        return text

    def title(self, text):
        self.print(text)
        self.print('=' * len(text))
        self.print()

    def begin_module_section(self, modname):
        self.print(modname)
        self.print('-' * len(modname))
        self.print()

    def end_module_section(self):
        self.print()

    def write_supported_item(self, modname, itemname, typename, explained,
                             sources, alias):
        self.print('.. function:: {}.{}'.format(modname, itemname))
        self.print('   :noindex:')
        self.print()

        if alias:
            self.print("   Alias to: ``{}``".format(alias))
        self.print()

        for tcls, source in sources.items():
            if source:
                impl = source['name']
                sig = source['sig']
                filename = source['filename']
                lines = source['lines']
                source_link = github_url.format(
                    commit=commit,
                    path=filename,
                    firstline=lines[0],
                    lastline=lines[1],
                )
                self.print(
                    "   - defined by ``{}{}`` at `{}:{}-{} <{}>`_".format(
                        impl, sig, filename, lines[0], lines[1], source_link,
                    ),
                )

            else:
                self.print("   - defined by ``{}``".format(str(tcls)))
        self.print()

    def write_unsupported_item(self, modname, itemname):
        pass

    def write_statistic(self, stat):
        if stat.supported == 0:
            self.print("This module is not supported.")
        else:
            msg = "Not showing {} unsupported functions."
            self.print(msg.format(stat.unsupported))
            self.print()
            self.print(stat.describe())
        self.print()


def _format_module_infos(formatter, package_name, mod_sequence, target=None):
    """Format modules.
    """
    formatter.title('Listings for {}'.format(package_name))
    alias_map = {}  # remember object seen to track alias
    for mod in mod_sequence:
        stat = _Stat()
        modname = mod.__name__
        formatter.begin_module_section(formatter.escape(modname))
        for info in inspect_module(mod, target=target, alias=alias_map):
            nbtype = info['numba_type']
            if nbtype is not None:
                stat.supported += 1
                formatter.write_supported_item(
                    modname=formatter.escape(info['module'].__name__),
                    itemname=formatter.escape(info['name']),
                    typename=formatter.escape(str(nbtype)),
                    explained=formatter.escape(info['explained']),
                    sources=info['source_infos'],
                    alias=info.get('alias'),
                )

            else:
                stat.unsupported += 1
                formatter.write_unsupported_item(
                    modname=formatter.escape(info['module'].__name__),
                    itemname=formatter.escape(info['name']),
                )

        formatter.write_statistic(stat)
        formatter.end_module_section()


def write_listings(package_name, filename, output_format):
    """Write listing information into a file.

    Parameters
    ----------
    package_name : str
        Name of the package to inspect.
    filename : str
        Output filename. Always overwrite.
    output_format : str
        Support formats are "html" and "rst".
    """
    package = __import__(package_name)
    if hasattr(package, '__path__'):
        mods = list_modules_in_package(package)
    else:
        mods = [package]

    if output_format == 'html':
        with open(filename + '.html', 'w') as fout:
            fmtr = HTMLFormatter(fileobj=fout)
            _format_module_infos(fmtr, package_name, mods)
    elif output_format == 'rst':
        with open(filename + '.rst', 'w') as fout:
            fmtr = ReSTFormatter(fileobj=fout)
            _format_module_infos(fmtr, package_name, mods)
    else:
        raise ValueError(
            "Output format '{}' is not supported".format(output_format))


program_description = """
Inspect Numba support for a given top-level package.
""".strip()


def main():
    parser = argparse.ArgumentParser(description=program_description)
    parser.add_argument(
        'package', metavar='package', type=str,
        help='Package to inspect',
    )
    parser.add_argument(
        '--format', dest='format', default='html',
        help='Output format; i.e. "html", "rst"',
    )
    parser.add_argument(
        '--file', dest='file', default='inspector_output',
        help='Output filename. Defaults to "inspector_output.<format>"',
    )

    args = parser.parse_args()
    package_name = args.package
    output_format = args.format
    filename = args.file
    write_listings(package_name, filename, output_format)


if __name__ == '__main__':
    main()
