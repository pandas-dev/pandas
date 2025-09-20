# SPDX-FileCopyrightText: 2021 Filipe Laíns <lains@riseup.net>
# SPDX-FileCopyrightText: 2021 Quansight, LLC
# SPDX-FileCopyrightText: 2022 The meson-python developers
#
# SPDX-License-Identifier: MIT

"""Meson Python build backend

Implements PEP 517 hooks.
"""

from __future__ import annotations

import argparse
import collections
import contextlib
import difflib
import functools
import importlib.machinery
import io
import itertools
import json
import os
import pathlib
import platform
import re
import shutil
import subprocess
import sys
import sysconfig
import tarfile
import tempfile
import textwrap
import typing
import warnings

from typing import Dict


if sys.version_info < (3, 11):
    import tomli as tomllib
else:
    import tomllib

import packaging.version
import pyproject_metadata

import mesonpy._compat
import mesonpy._dylib
import mesonpy._elf
import mesonpy._tags
import mesonpy._util
import mesonpy._wheelfile

from mesonpy._compat import Collection, Mapping, cached_property, read_binary


if typing.TYPE_CHECKING:  # pragma: no cover
    from typing import (
        Any, Callable, ClassVar, DefaultDict, List, Literal, Optional, Sequence, TextIO, Tuple, Type, TypeVar, Union
    )

    from mesonpy._compat import Iterator, ParamSpec, Path

    P = ParamSpec('P')
    T = TypeVar('T')


__version__ = '0.13.1'


# XXX: Once Python 3.8 is our minimum supported version, get rid of
#      meson_args_keys and use typing.get_args(MesonArgsKeys) instead.

# Keep both definitions in sync!
_MESON_ARGS_KEYS = ['dist', 'setup', 'compile', 'install']
if typing.TYPE_CHECKING:
    MesonArgsKeys = Literal['dist', 'setup', 'compile', 'install']
    MesonArgs = Mapping[MesonArgsKeys, List[str]]
else:
    MesonArgs = dict


_COLORS = {
    'red': '\33[31m',
    'cyan': '\33[36m',
    'yellow': '\33[93m',
    'light_blue': '\33[94m',
    'bold': '\33[1m',
    'dim': '\33[2m',
    'underline': '\33[4m',
    'reset': '\33[0m',
}
_NO_COLORS = {color: '' for color in _COLORS}

_NINJA_REQUIRED_VERSION = '1.8.2'
_MESON_REQUIRED_VERSION = '0.63.3' # keep in sync with the version requirement in pyproject.toml


class _depstr:
    """Namespace that holds the requirement strings for dependencies we *might*
    need at runtime. Having them in one place makes it easier to update.
    """
    patchelf = 'patchelf >= 0.11.0'
    ninja = f'ninja >= {_NINJA_REQUIRED_VERSION}'


def _init_colors() -> Dict[str, str]:
    """Detect if we should be using colors in the output. We will enable colors
    if running in a TTY, and no environment variable overrides it. Setting the
    NO_COLOR (https://no-color.org/) environment variable force-disables colors,
    and FORCE_COLOR forces color to be used, which is useful for thing like
    Github actions.
    """
    if 'NO_COLOR' in os.environ:
        if 'FORCE_COLOR' in os.environ:
            warnings.warn(
                'Both NO_COLOR and FORCE_COLOR environment variables are set, disabling color',
                stacklevel=1,
            )
        return _NO_COLORS
    elif 'FORCE_COLOR' in os.environ or sys.stdout.isatty():
        return _COLORS
    return _NO_COLORS


_STYLES = _init_colors()  # holds the color values, should be _COLORS or _NO_COLORS


_SUFFIXES = importlib.machinery.all_suffixes()
_EXTENSION_SUFFIXES = importlib.machinery.EXTENSION_SUFFIXES
_EXTENSION_SUFFIX_REGEX = re.compile(r'^\.(?:(?P<abi>[^.]+)\.)?(?:so|pyd|dll)$')
assert all(re.match(_EXTENSION_SUFFIX_REGEX, x) for x in _EXTENSION_SUFFIXES)


# Map Meson installation path placeholders to wheel installation paths.
# See https://docs.python.org/3/library/sysconfig.html#installation-paths
_INSTALLATION_PATH_MAP = {
    '{bindir}': 'scripts',
    '{py_purelib}': 'purelib',
    '{py_platlib}': 'platlib',
    '{moduledir_shared}': 'platlib',
    '{includedir}': 'headers',
    '{datadir}': 'data',
    # custom location
    '{libdir}': 'mesonpy-libs',
    '{libdir_shared}': 'mesonpy-libs',
}


def _map_to_wheel(sources: Dict[str, Dict[str, Any]]) -> DefaultDict[str, List[Tuple[pathlib.Path, str]]]:
    """Map files to the wheel, organized by wheel installation directrory."""
    wheel_files: DefaultDict[str, List[Tuple[pathlib.Path, str]]] = collections.defaultdict(list)
    packages: Dict[str, str] = {}

    for group in sources.values():
        for src, target in group.items():
            destination = pathlib.Path(target['destination'])
            anchor = destination.parts[0]

            path = _INSTALLATION_PATH_MAP.get(anchor)
            if path is None:
                raise BuildError(f'Could not map installation path to an equivalent wheel directory: {str(destination)!r}')

            if path == 'purelib' or path == 'platlib':
                package = destination.parts[1]
                other = packages.setdefault(package, path)
                if other != path:
                    this = os.fspath(pathlib.Path(path, *destination.parts[1:]))
                    that = os.fspath(other / next(d for d, s in wheel_files[other] if d.parts[0] == destination.parts[1]))
                    raise BuildError(
                        f'The {package} package is split between {path} and {other}: '
                        f'{this!r} and {that!r}, a "pure: false" argument may be missing in meson.build. '
                        f'It is recommended to set it in "import(\'python\').find_installation()"')

            wheel_files[path].append((pathlib.Path(*destination.parts[1:]), src))
    return wheel_files


def _showwarning(
    message: Union[Warning, str],
    category: Type[Warning],
    filename: str,
    lineno: int,
    file: Optional[TextIO] = None,
    line: Optional[str] = None,
) -> None:  # pragma: no cover
    """Callable to override the default warning handler, to have colored output."""
    print('{yellow}meson-python: warning:{reset} {}'.format(message, **_STYLES))


def _setup_cli() -> None:
    """Setup CLI stuff (eg. handlers, hooks, etc.). Should only be called when
    actually we are in control of the CLI, not on a normal import.
    """
    warnings.showwarning = _showwarning

    try:  # pragma: no cover
        import colorama
    except ModuleNotFoundError:  # pragma: no cover
        pass
    else:  # pragma: no cover
        colorama.init()  # fix colors on windows


class Error(RuntimeError):
    def __str__(self) -> str:
        return str(self.args[0])


class ConfigError(Error):
    """Error in the backend configuration."""


class BuildError(Error):
    """Error when building the wheel."""


class MesonBuilderError(Error):
    """Error when building the Meson package."""


class _WheelBuilder():
    """Helper class to build wheels from projects."""

    def __init__(
        self,
        project: Project,
        metadata: Optional[pyproject_metadata.StandardMetadata],
        source_dir: pathlib.Path,
        build_dir: pathlib.Path,
        sources: Dict[str, Dict[str, Any]],
    ) -> None:
        self._project = project
        self._metadata = metadata
        self._source_dir = source_dir
        self._build_dir = build_dir
        self._sources = sources

        self._libs_build_dir = self._build_dir / 'mesonpy-wheel-libs'

    @cached_property
    def _wheel_files(self) -> DefaultDict[str, List[Tuple[pathlib.Path, str]]]:
        return _map_to_wheel(self._sources)

    @property
    def _has_internal_libs(self) -> bool:
        return bool(self._wheel_files['mesonpy-libs'])

    @property
    def _has_extension_modules(self) -> bool:
        # Assume that all code installed in {platlib} is Python ABI dependent.
        return bool(self._wheel_files['platlib'])

    @property
    def normalized_name(self) -> str:
        return self._project.name.replace('-', '_')

    @property
    def basename(self) -> str:
        """Normalized wheel name and version (eg. meson_python-1.0.0)."""
        return '{distribution}-{version}'.format(
            distribution=self.normalized_name,
            version=self._project.version,
        )

    @property
    def tag(self) -> mesonpy._tags.Tag:
        """Wheel tags."""
        if self.is_pure:
            return mesonpy._tags.Tag('py3', 'none', 'any')
        if not self._has_extension_modules:
            # The wheel has platform dependent code (is not pure) but
            # does not contain any extension module (does not
            # distribute any file in {platlib}) thus use generic
            # implementation and ABI tags.
            return mesonpy._tags.Tag('py3', 'none', None)
        return mesonpy._tags.Tag(None, self._stable_abi, None)

    @property
    def name(self) -> str:
        """Wheel name, this includes the basename and tag."""
        return '{basename}-{tag}'.format(
            basename=self.basename,
            tag=self.tag,
        )

    @property
    def distinfo_dir(self) -> str:
        return f'{self.basename}.dist-info'

    @property
    def data_dir(self) -> str:
        return f'{self.basename}.data'

    @cached_property
    def is_pure(self) -> bool:
        """Is the wheel "pure" (architecture independent)?"""
        # XXX: I imagine some users might want to force the package to be
        # non-pure, but I think it's better that we evaluate use-cases as they
        # arise and make sure allowing the user to override this is indeed the
        # best option for the use-case.
        if self._wheel_files['platlib']:
            return False
        for _, file in self._wheel_files['scripts']:
            if self._is_native(file):
                return False
        return True

    @property
    def wheel(self) -> bytes:
        """Return WHEEL file for dist-info."""
        return textwrap.dedent('''
            Wheel-Version: 1.0
            Generator: meson
            Root-Is-Purelib: {is_purelib}
            Tag: {tag}
        ''').strip().format(
            is_purelib='true' if self.is_pure else 'false',
            tag=self.tag,
        ).encode()

    @property
    def entrypoints_txt(self) -> bytes:
        """dist-info entry_points.txt."""
        if not self._metadata:
            return b''

        data = self._metadata.entrypoints.copy()
        data.update({
            'console_scripts': self._metadata.scripts,
            'gui_scripts': self._metadata.gui_scripts,
        })

        text = ''
        for entrypoint in data:
            if data[entrypoint]:
                text += f'[{entrypoint}]\n'
                for name, target in data[entrypoint].items():
                    text += f'{name} = {target}\n'
                text += '\n'

        return text.encode()

    @cached_property
    def _stable_abi(self) -> Optional[str]:
        """Determine stabe ABI compatibility.

        Examine all files installed in {platlib} that look like
        extension modules (extension .pyd on Windows, .dll on Cygwin,
        and .so on other platforms) and, if they all share the same
        PEP 3149 filename stable ABI tag, return it.

        Other files are ignored.

        """
        soext = sorted(_EXTENSION_SUFFIXES, key=len)[0]
        abis = []

        for path, _ in self._wheel_files['platlib']:
            # NOTE: When searching for shared objects files, we assume the host
            # and build machines have the same soext, even though that we might
            # be cross compiling.
            if path.suffix == soext:
                match = re.match(r'^[^.]+(.*)$', path.name)
                assert match is not None
                suffix = match.group(1)
                match = _EXTENSION_SUFFIX_REGEX.match(suffix)
                if match:
                    abis.append(match.group('abi'))

        stable = [x for x in abis if x and re.match(r'abi\d+', x)]
        if len(stable) > 0 and len(stable) == len(abis) and all(x == stable[0] for x in stable[1:]):
            return stable[0]
        return None

    @property
    def top_level_modules(self) -> Collection[str]:
        modules = set()
        for type_ in self._wheel_files:
            for path, _ in self._wheel_files[type_]:
                name, dot, ext = path.parts[0].partition('.')
                if dot:
                    # module
                    suffix = dot + ext
                    if suffix in _SUFFIXES:
                        modules.add(name)
                else:
                    # package
                    modules.add(name)
        return modules

    def _is_native(self, file: Union[str, pathlib.Path]) -> bool:
        """Check if file is a native file."""
        self._project.build()  # the project needs to be built for this :/

        with open(file, 'rb') as f:
            if platform.system() == 'Linux':
                return f.read(4) == b'\x7fELF'  # ELF
            elif platform.system() == 'Darwin':
                return f.read(4) in (
                    b'\xfe\xed\xfa\xce',  # 32-bit
                    b'\xfe\xed\xfa\xcf',  # 64-bit
                    b'\xcf\xfa\xed\xfe',  # arm64
                    b'\xca\xfe\xba\xbe',  # universal / fat (same as java class so beware!)
                )
            elif platform.system() == 'Windows':
                return f.read(2) == b'MZ'

        # For unknown platforms, check for file extensions.
        _, ext = os.path.splitext(file)
        if ext in ('.so', '.a', '.out', '.exe', '.dll', '.dylib', '.pyd'):
            return True
        return False

    def _install_path(  # noqa: C901
        self,
        wheel_file: mesonpy._wheelfile.WheelFile,
        counter: mesonpy._util.CLICounter,
        origin: Path,
        destination: pathlib.Path,
    ) -> None:
        """"Install" file or directory into the wheel
        and do the necessary processing before doing so.

        Some files might need to be fixed up to set the RPATH to the internal
        library directory on Linux wheels for eg.
        """
        location = destination.as_posix()
        counter.update(location)

        # fix file
        if os.path.isdir(origin):
            for root, dirnames, filenames in os.walk(str(origin)):
                # Sort the directory names so that `os.walk` will walk them in a
                # defined order on the next iteration.
                dirnames.sort()
                for name in sorted(filenames):
                    path = os.path.normpath(os.path.join(root, name))
                    if os.path.isfile(path):
                        arcname = os.path.join(destination, os.path.relpath(path, origin).replace(os.path.sep, '/'))
                        wheel_file.write(path, arcname)
        else:
            if self._has_internal_libs:
                if platform.system() == 'Linux' or platform.system() == 'Darwin':
                    # add .mesonpy.libs to the RPATH of ELF files
                    if self._is_native(os.fspath(origin)):
                        # copy ELF to our working directory to avoid Meson having to regenerate the file
                        new_origin = self._libs_build_dir / pathlib.Path(origin).relative_to(self._build_dir)
                        os.makedirs(new_origin.parent, exist_ok=True)
                        shutil.copy2(origin, new_origin)
                        origin = new_origin
                        # add our in-wheel libs folder to the RPATH
                        if platform.system() == 'Linux':
                            elf = mesonpy._elf.ELF(origin)
                            libdir_path = \
                                f'$ORIGIN/{os.path.relpath(f".{self._project.name}.mesonpy.libs", destination.parent)}'
                            if libdir_path not in elf.rpath:
                                elf.rpath = [*elf.rpath, libdir_path]
                        elif platform.system() == 'Darwin':
                            dylib = mesonpy._dylib.Dylib(origin)
                            libdir_path = \
                                f'@loader_path/{os.path.relpath(f".{self._project.name}.mesonpy.libs", destination.parent)}'
                            if libdir_path not in dylib.rpath:
                                dylib.rpath = [*dylib.rpath, libdir_path]
                        else:
                            # Internal libraries are currently unsupported on this platform
                            raise NotImplementedError("Bundling libraries in wheel is not supported on platform '{}'"
                                                      .format(platform.system()))

            try:
                wheel_file.write(origin, location)
            except FileNotFoundError:
                # work around for Meson bug, see https://github.com/mesonbuild/meson/pull/11655
                if not os.fspath(origin).endswith('.pdb'):
                    raise

    def _wheel_write_metadata(self, whl: mesonpy._wheelfile.WheelFile) -> None:
        # add metadata
        whl.writestr(f'{self.distinfo_dir}/METADATA', self._project.metadata)
        whl.writestr(f'{self.distinfo_dir}/WHEEL', self.wheel)
        if self.entrypoints_txt:
            whl.writestr(f'{self.distinfo_dir}/entry_points.txt', self.entrypoints_txt)

        # add license (see https://github.com/mesonbuild/meson-python/issues/88)
        if self._project.license_file:
            whl.write(
                self._source_dir / self._project.license_file,
                f'{self.distinfo_dir}/{os.path.basename(self._project.license_file)}',
            )

    def build(self, directory: Path) -> pathlib.Path:
        # ensure project is built
        self._project.build()
        # install the project
        self._project.install()

        wheel_file = pathlib.Path(directory, f'{self.name}.whl')

        with mesonpy._wheelfile.WheelFile(wheel_file, 'w') as whl:
            self._wheel_write_metadata(whl)

            with mesonpy._util.cli_counter(
                len(list(itertools.chain.from_iterable(self._wheel_files.values()))),
            ) as counter:
                # install root scheme files
                root_scheme = 'purelib' if self.is_pure else 'platlib'
                for destination, origin in self._wheel_files[root_scheme]:
                    self._install_path(whl, counter, origin, destination)

                # install bundled libraries
                for destination, origin in self._wheel_files['mesonpy-libs']:
                    destination = pathlib.Path(f'.{self._project.name}.mesonpy.libs', destination)
                    self._install_path(whl, counter, origin, destination)

                # install the other schemes
                for scheme in self._wheel_files.keys():
                    if scheme in (root_scheme, 'mesonpy-libs'):
                        continue
                    for destination, origin in self._wheel_files[scheme]:
                        destination = pathlib.Path(self.data_dir, scheme, destination)
                        self._install_path(whl, counter, origin, destination)

        return wheel_file

    def build_editable(self, directory: Path, verbose: bool = False) -> pathlib.Path:
        # ensure project is built
        self._project.build()

        wheel_file = pathlib.Path(directory, f'{self.name}.whl')

        with mesonpy._wheelfile.WheelFile(wheel_file, 'w') as whl:
            self._wheel_write_metadata(whl)
            whl.writestr(
                f'{self.distinfo_dir}/direct_url.json',
                self._source_dir.as_uri().encode('utf-8'),
            )

            # install loader module
            loader_module_name = f'_{self.normalized_name.replace(".", "_")}_editable_loader'
            whl.writestr(
                f'{loader_module_name}.py',
                read_binary('mesonpy', '_editable.py') + textwrap.dedent(f'''
                   install(
                       {self.top_level_modules!r},
                       {os.fspath(self._build_dir)!r},
                       {self._project._build_command!r},
                       {verbose!r},
                   )''').encode('utf-8'))

            # install .pth file
            whl.writestr(
                f'{self.normalized_name}-editable.pth',
                f'import {loader_module_name}'.encode('utf-8'))

        return wheel_file


def _validate_pyproject_config(pyproject: Dict[str, Any]) -> Dict[str, Any]:

    def _table(scheme: Dict[str, Callable[[Any, str], Any]]) -> Callable[[Any, str], Dict[str, Any]]:
        def func(value: Any, name: str) -> Dict[str, Any]:
            if not isinstance(value, dict):
                raise ConfigError(f'Configuration entry "{name}" must be a table')
            table = {}
            for key, val in value.items():
                check = scheme.get(key)
                if check is None:
                    raise ConfigError(f'Unknown configuration entry "{name}.{key}"')
                table[key] = check(val, f'{name}.{key}')
            return table
        return func

    def _strings(value: Any, name: str) -> List[str]:
        if not isinstance(value, list) or not all(isinstance(x, str) for x in value):
            raise ConfigError(f'Configuration entry "{name}" must be a list of strings')
        return value

    scheme = _table({
        'args': _table({
            name: _strings for name in _MESON_ARGS_KEYS
        })
    })

    table = pyproject.get('tool', {}).get('meson-python', {})
    return scheme(table, 'tool.meson-python')


def _validate_config_settings(config_settings: Dict[str, Any]) -> Dict[str, Any]:
    """Validate options received from build frontend."""

    def _string(value: Any, name: str) -> str:
        if not isinstance(value, str):
            raise ConfigError(f'Only one value for "{name}" can be specified')
        return value

    def _bool(value: Any, name: str) -> bool:
        return True

    def _string_or_strings(value: Any, name: str) -> List[str]:
        return list([value,] if isinstance(value, str) else value)

    options = {
        'builddir': _string,
        'editable-verbose': _bool,
        'dist-args': _string_or_strings,
        'setup-args': _string_or_strings,
        'compile-args': _string_or_strings,
        'install-args': _string_or_strings,
    }
    assert all(f'{name}-args' in options for name in _MESON_ARGS_KEYS)

    config = {}
    for key, value in config_settings.items():
        parser = options.get(key)
        if parser is None:
            matches = difflib.get_close_matches(key, options.keys(), n=2)
            if matches:
                alternatives = ' or '.join(f'"{match}"' for match in matches)
                raise ConfigError(f'Unknown option "{key}". Did you mean {alternatives}?')
            else:
                raise ConfigError(f'Unknown option "{key}"')
        config[key] = parser(value, key)
    return config


class Project():
    """Meson project wrapper to generate Python artifacts."""

    _ALLOWED_DYNAMIC_FIELDS: ClassVar[List[str]] = [
        'version',
    ]
    _metadata: pyproject_metadata.StandardMetadata

    def __init__(
        self,
        source_dir: Path,
        working_dir: Path,
        build_dir: Optional[Path] = None,
        meson_args: Optional[MesonArgs] = None,
        editable_verbose: bool = False,
    ) -> None:
        self._source_dir = pathlib.Path(source_dir).absolute()
        self._working_dir = pathlib.Path(working_dir).absolute()
        self._build_dir = pathlib.Path(build_dir).absolute() if build_dir else (self._working_dir / 'build')
        self._editable_verbose = editable_verbose
        self._install_dir = self._working_dir / 'install'
        self._meson_native_file = self._build_dir / 'meson-python-native-file.ini'
        self._meson_cross_file = self._build_dir / 'meson-python-cross-file.ini'
        self._meson_args: MesonArgs = collections.defaultdict(list)

        _check_meson_version()

        self._ninja = _env_ninja_command()
        if self._ninja is None:
            raise ConfigError(f'Could not find ninja version {_NINJA_REQUIRED_VERSION} or newer.')
        os.environ.setdefault('NINJA', self._ninja)

        # make sure the build dir exists
        self._build_dir.mkdir(exist_ok=True, parents=True)
        self._install_dir.mkdir(exist_ok=True, parents=True)

        # setuptools-like ARCHFLAGS environment variable support
        if sysconfig.get_platform().startswith('macosx-'):
            archflags = os.environ.get('ARCHFLAGS', '').strip()
            if archflags:
                arch, *other = filter(None, (x.strip() for x in archflags.split('-arch')))
                if other:
                    raise ConfigError(f'Multi-architecture builds are not supported but $ARCHFLAGS={archflags!r}')
                macver, _, nativearch = platform.mac_ver()
                if arch != nativearch:
                    x = os.environ.setdefault('_PYTHON_HOST_PLATFORM', f'macosx-{macver}-{arch}')
                    if not x.endswith(arch):
                        raise ConfigError(f'$ARCHFLAGS={archflags!r} and $_PYTHON_HOST_PLATFORM={x!r} do not agree')
                    family = 'aarch64' if arch == 'arm64' else arch
                    cross_file_data = textwrap.dedent(f'''
                        [binaries]
                        c = ['cc', '-arch', {arch!r}]
                        cpp = ['c++', '-arch', {arch!r}]
                        [host_machine]
                        system = 'Darwin'
                        cpu = {arch!r}
                        cpu_family = {family!r}
                        endian = 'little'
                    ''')
                    self._meson_cross_file.write_text(cross_file_data)
                    self._meson_args['setup'].extend(('--cross-file', os.fspath(self._meson_cross_file)))

        # load pyproject.toml
        pyproject = tomllib.loads(self._source_dir.joinpath('pyproject.toml').read_text())

        # load meson args from pyproject.toml
        pyproject_config = _validate_pyproject_config(pyproject)
        for key, value in pyproject_config.get('args', {}).items():
            self._meson_args[key].extend(value)

        # meson arguments from the command line take precedence over
        # arguments from the configuration file thus are added later
        if meson_args:
            for key, value in meson_args.items():
                self._meson_args[key].extend(value)

        # write the native file
        native_file_data = textwrap.dedent(f'''
            [binaries]
            python = '{sys.executable}'
        ''')
        self._meson_native_file.write_text(native_file_data)

        # reconfigure if we have a valid Meson build directory. Meson
        # uses the presence of the 'meson-private/coredata.dat' file
        # in the build directory as indication that the build
        # directory has already been configured and arranges this file
        # to be created as late as possible or deleted if something
        # goes wrong during setup.
        reconfigure = self._build_dir.joinpath('meson-private/coredata.dat').is_file()

        # run meson setup
        self._configure(reconfigure=reconfigure)

        # package metadata
        if 'project' in pyproject:
            self._metadata = pyproject_metadata.StandardMetadata.from_pyproject(pyproject, self._source_dir)
        else:
            self._metadata = pyproject_metadata.StandardMetadata(
                name=self._meson_name, version=packaging.version.Version(self._meson_version))
        self._validate_metadata()

        # set version from meson.build if dynamic
        if 'version' in self._metadata.dynamic:
            self._metadata.version = packaging.version.Version(self._meson_version)

    def _run(self, cmd: Sequence[str]) -> None:
        """Invoke a subprocess."""
        # Flush the line to ensure that the log line with the executed
        # command line appears before the command output. Without it,
        # the lines appear in the wrong order in pip output.
        print('{cyan}{bold}+ {}{reset}'.format(' '.join(cmd), **_STYLES), flush=True)
        r = subprocess.run(cmd, cwd=self._build_dir)
        if r.returncode != 0:
            raise SystemExit(r.returncode)

    def _configure(self, reconfigure: bool = False) -> None:
        """Configure Meson project."""
        setup_args = [
            os.fspath(self._source_dir),
            os.fspath(self._build_dir),
            # default build options
            '-Dbuildtype=release',
            '-Db_ndebug=if-release',
            '-Db_vscrt=md',
            # user build options
            *self._meson_args['setup'],
            # pass native file last to have it override the python
            # interpreter path that may have been specified in user
            # provided native files
            f'--native-file={os.fspath(self._meson_native_file)}',
        ]
        if reconfigure:
            setup_args.insert(0, '--reconfigure')

        self._run(['meson', 'setup', *setup_args])

    def _validate_metadata(self) -> None:
        """Check the pyproject.toml metadata and see if there are any issues."""

        # check for unsupported dynamic fields
        unsupported_dynamic = {
            key for key in self._metadata.dynamic
            if key not in self._ALLOWED_DYNAMIC_FIELDS
        }
        if unsupported_dynamic:
            s = ', '.join(f'"{x}"' for x in unsupported_dynamic)
            raise MesonBuilderError(f'Unsupported dynamic fields: {s}')

        # check if we are running on an unsupported interpreter
        if self._metadata.requires_python:
            self._metadata.requires_python.prereleases = True
            if platform.python_version().rstrip('+') not in self._metadata.requires_python:
                raise MesonBuilderError(
                    f'Unsupported Python version {platform.python_version()}, '
                    f'expected {self._metadata.requires_python}'
                )

    @cached_property
    def _wheel_builder(self) -> _WheelBuilder:
        return _WheelBuilder(
            self,
            self._metadata,
            self._source_dir,
            self._build_dir,
            self._install_plan,
        )

    @property
    def _build_command(self) -> List[str]:
        assert self._ninja is not None  # help mypy out
        if platform.system() == 'Windows':
            # On Windows use 'meson compile' to setup the MSVC compiler
            # environment. Using the --ninja-args option allows to
            # provide the exact same semantics for the compile arguments
            # provided by the users.
            cmd = ['meson', 'compile']
            args = list(self._meson_args['compile'])
            if args:
                cmd.append(f'--ninja-args={args!r}')
            return cmd
        return [self._ninja, *self._meson_args['compile']]

    @functools.lru_cache(maxsize=None)
    def build(self) -> None:
        """Build the Meson project."""
        self._run(self._build_command)

    def install(self) -> None:
        """Install the Meson project."""
        destdir = os.fspath(self._install_dir)
        self._run(['meson', 'install', '--quiet', '--no-rebuild', '--destdir', destdir, *self._meson_args['install']])

    @classmethod
    @contextlib.contextmanager
    def with_temp_working_dir(
        cls,
        source_dir: Path = os.path.curdir,
        build_dir: Optional[Path] = None,
        meson_args: Optional[MesonArgs] = None,
        editable_verbose: bool = False,
    ) -> Iterator[Project]:
        """Creates a project instance pointing to a temporary working directory."""
        with tempfile.TemporaryDirectory(prefix='.mesonpy-', dir=os.fspath(source_dir)) as tmpdir:
            yield cls(source_dir, tmpdir, build_dir, meson_args, editable_verbose)

    @functools.lru_cache()
    def _info(self, name: str) -> Dict[str, Any]:
        """Read info from meson-info directory."""
        file = self._build_dir.joinpath('meson-info', f'{name}.json')
        return typing.cast(
            Dict[str, str],
            json.loads(file.read_text())
        )

    @property
    def _install_plan(self) -> Dict[str, Dict[str, Dict[str, str]]]:
        """Meson install_plan metadata."""

        # copy the install plan so we can modify it
        install_plan = self._info('intro-install_plan').copy()

        # parse install args for install tags (--tags)
        parser = argparse.ArgumentParser()
        parser.add_argument('--tags')
        args, _ = parser.parse_known_args(self._meson_args['install'])

        # filter the install_plan for files that do not fit the install tags
        if args.tags:
            install_tags = args.tags.split(',')

            for files in install_plan.values():
                for file, details in list(files.items()):
                    if details['tag'].strip() not in install_tags:
                        del files[file]

        return install_plan

    @property
    def _meson_name(self) -> str:
        """Name in meson.build."""
        name = self._info('intro-projectinfo')['descriptive_name']
        assert isinstance(name, str)
        return name

    @property
    def _meson_version(self) -> str:
        """Version in meson.build."""
        name = self._info('intro-projectinfo')['version']
        assert isinstance(name, str)
        return name

    @property
    def name(self) -> str:
        """Project name."""
        return str(self._metadata.name).replace('-', '_')

    @property
    def version(self) -> str:
        """Project version."""
        return str(self._metadata.version)

    @cached_property
    def metadata(self) -> bytes:
        """Project metadata as an RFC822 message."""
        return bytes(self._metadata.as_rfc822())

    @property
    def license_file(self) -> Optional[pathlib.Path]:
        if self._metadata:
            license_ = self._metadata.license
            if license_ and license_.file:
                return pathlib.Path(license_.file)
        return None

    @property
    def is_pure(self) -> bool:
        """Is the wheel "pure" (architecture independent)?"""
        return bool(self._wheel_builder.is_pure)

    def sdist(self, directory: Path) -> pathlib.Path:
        """Generates a sdist (source distribution) in the specified directory."""
        # generate meson dist file
        self._run(['meson', 'dist', '--allow-dirty', '--no-tests', '--formats', 'gztar', *self._meson_args['dist']])

        # move meson dist file to output path
        dist_name = f'{self.name}-{self.version}'
        meson_dist_name = f'{self._meson_name}-{self._meson_version}'
        meson_dist_path = pathlib.Path(self._build_dir, 'meson-dist', f'{meson_dist_name}.tar.gz')
        sdist = pathlib.Path(directory, f'{dist_name}.tar.gz')

        with tarfile.open(meson_dist_path, 'r:gz') as meson_dist, mesonpy._util.create_targz(sdist) as (tar, mtime):
            for member in meson_dist.getmembers():
                # calculate the file path in the source directory
                assert member.name, member.name
                member_parts = member.name.split('/')
                if len(member_parts) <= 1:
                    continue
                path = self._source_dir.joinpath(*member_parts[1:])

                if not path.exists() and member.isfile():
                    # File doesn't exists on the source directory but exists on
                    # the Meson dist, so it is generated file, which we need to
                    # include.
                    # See https://mesonbuild.com/Reference-manual_builtin_meson.html#mesonadd_dist_script

                    # MESON_DIST_ROOT could have a different base name
                    # than the actual sdist basename, so we need to rename here
                    file = meson_dist.extractfile(member.name)
                    member.name = str(pathlib.Path(dist_name, *member_parts[1:]).as_posix())
                    tar.addfile(member, file)
                    continue

                if not path.is_file():
                    continue

                info = tarfile.TarInfo(member.name)
                file_stat = os.stat(path)
                info.size = file_stat.st_size
                info.mode = int(oct(file_stat.st_mode)[-3:], 8)

                # rewrite the path if necessary, to match the sdist distribution name
                if dist_name != meson_dist_name:
                    info.name = pathlib.Path(
                        dist_name,
                        path.relative_to(self._source_dir)
                    ).as_posix()

                with path.open('rb') as f:
                    tar.addfile(info, fileobj=f)

            # add PKG-INFO to dist file to make it a sdist
            pkginfo_info = tarfile.TarInfo(f'{dist_name}/PKG-INFO')
            if mtime:
                pkginfo_info.mtime = mtime
            pkginfo_info.size = len(self.metadata)
            tar.addfile(pkginfo_info, fileobj=io.BytesIO(self.metadata))

        return sdist

    def wheel(self, directory: Path) -> pathlib.Path:
        """Generates a wheel (binary distribution) in the specified directory."""
        file = self._wheel_builder.build(directory)
        assert isinstance(file, pathlib.Path)
        return file

    def editable(self, directory: Path) -> pathlib.Path:
        file = self._wheel_builder.build_editable(directory, self._editable_verbose)
        assert isinstance(file, pathlib.Path)
        return file


@contextlib.contextmanager
def _project(config_settings: Optional[Dict[Any, Any]]) -> Iterator[Project]:
    """Create the project given the given config settings."""

    settings = _validate_config_settings(config_settings or {})
    meson_args = {name: settings.get(f'{name}-args', []) for name in _MESON_ARGS_KEYS}

    with Project.with_temp_working_dir(
            build_dir=settings.get('builddir'),
            meson_args=typing.cast(MesonArgs, meson_args),
            editable_verbose=bool(settings.get('editable-verbose'))
    ) as project:
        yield project


def _parse_version_string(string: str) -> Tuple[int, ...]:
    """Parse version string."""
    try:
        return tuple(map(int, string.split('.')[:3]))
    except ValueError:
        return (0, )


def _env_ninja_command(*, version: str = _NINJA_REQUIRED_VERSION) -> Optional[str]:
    """Returns the path to ninja, or None if no ninja found."""
    required_version = _parse_version_string(version)
    env_ninja = os.environ.get('NINJA')
    ninja_candidates = [env_ninja] if env_ninja else ['ninja', 'ninja-build', 'samu']
    for ninja in ninja_candidates:
        ninja_path = shutil.which(ninja)
        if ninja_path is not None:
            version = subprocess.run([ninja_path, '--version'], check=False, text=True, capture_output=True).stdout
            if _parse_version_string(version) >= required_version:
                return ninja_path
    return None


def _check_meson_version(*, version: str = _MESON_REQUIRED_VERSION) -> None:
    """Check that the meson executable in the path has an appropriate version.

    The meson Python package is a dependency of the meson-python
    Python package, however, it may occur that the meson Python
    package is installed but the corresponding meson command is not
    available in $PATH. Implement a runtime check to verify that the
    build environment is setup correcly.

    """
    required_version = _parse_version_string(version)
    meson_version = subprocess.run(['meson', '--version'], check=False, text=True, capture_output=True).stdout
    if _parse_version_string(meson_version) < required_version:
        raise ConfigError(f'Could not find meson version {version} or newer, found {meson_version}.')


def _add_ignore_files(directory: pathlib.Path) -> None:
    directory.joinpath('.gitignore').write_text(textwrap.dedent('''
        # This file is generated by meson-python. It will not be recreated if deleted or modified.
        *
    '''), encoding='utf-8')
    directory.joinpath('.hgignore').write_text(textwrap.dedent('''
        # This file is generated by meson-python. It will not be recreated if deleted or modified.
        syntax: glob
        **/*
    '''), encoding='utf-8')


def _pyproject_hook(func: Callable[P, T]) -> Callable[P, T]:
    @functools.wraps(func)
    def wrapper(*args: P.args, **kwargs: P.kwargs) -> T:
        try:
            return func(*args, **kwargs)
        except (Error, pyproject_metadata.ConfigurationError) as exc:
            prefix = '{red}meson-python: error:{reset} '.format(**_STYLES)
            print('\n' + textwrap.indent(str(exc), prefix))
            raise SystemExit(1) from exc
    return wrapper


@_pyproject_hook
def get_requires_for_build_sdist(
    config_settings: Optional[Dict[str, str]] = None,
) -> List[str]:
    if os.environ.get('NINJA') is None and _env_ninja_command() is None:
        return [_depstr.ninja]
    return []


@_pyproject_hook
def build_sdist(
    sdist_directory: str,
    config_settings: Optional[Dict[Any, Any]] = None,
) -> str:
    _setup_cli()

    out = pathlib.Path(sdist_directory)
    with _project(config_settings) as project:
        return project.sdist(out).name


@_pyproject_hook
def get_requires_for_build_wheel(
    config_settings: Optional[Dict[str, str]] = None,
) -> List[str]:
    dependencies = []

    if os.environ.get('NINJA') is None and _env_ninja_command() is None:
        dependencies.append(_depstr.ninja)

    if sys.platform.startswith('linux'):
        # we may need patchelf
        if not shutil.which('patchelf'):
            # patchelf not already accessible on the system
            if _env_ninja_command() is not None:
                # we have ninja available, so we can run Meson and check if the project needs patchelf
                with _project(config_settings) as project:
                    if not project.is_pure:
                        dependencies.append(_depstr.patchelf)
            else:
                # we can't check if the project needs patchelf, so always add it
                # XXX: wait for https://github.com/mesonbuild/meson/pull/10779
                dependencies.append(_depstr.patchelf)

    return dependencies


@_pyproject_hook
def build_wheel(
    wheel_directory: str,
    config_settings: Optional[Dict[Any, Any]] = None,
    metadata_directory: Optional[str] = None,
) -> str:
    _setup_cli()

    out = pathlib.Path(wheel_directory)
    with _project(config_settings) as project:
        return project.wheel(out).name


@_pyproject_hook
def build_editable(
    wheel_directory: str,
    config_settings: Optional[Dict[Any, Any]] = None,
    metadata_directory: Optional[str] = None,
) -> str:
    _setup_cli()

    # force set a permanent builddir
    if not config_settings:
        config_settings = {}
    if 'builddir' not in config_settings:
        builddir = pathlib.Path('build')
        builddir.mkdir(exist_ok=True)
        if not next(builddir.iterdir(), None):
            _add_ignore_files(builddir)
        config_settings['builddir'] = os.fspath(builddir / str(mesonpy._tags.get_abi_tag()))

    out = pathlib.Path(wheel_directory)
    with _project(config_settings) as project:
        return project.editable(out).name


@_pyproject_hook
def get_requires_for_build_editable(
    config_settings: Optional[Dict[str, str]] = None,
) -> List[str]:
    return get_requires_for_build_wheel()
