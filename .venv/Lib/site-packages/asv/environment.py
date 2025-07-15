# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Manages an environment -- a combination of a version of Python and set
of dependencies.
"""


import hashlib
import os
import re
import sys
import itertools
import subprocess
import importlib
from pathlib import Path

from .console import log
from . import util, build_cache

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

WIN = (os.name == "nt")


def iter_matrix(environment_type, pythons, conf, explicit_selection=False):
    """
    Iterate through all combinations of the given requirement
    matrix and python versions.

    Yields
    ------
    combination : dict of {(key_type, key_name): value, ...}
        Combination of environment settings.
        Possible key types are ('req', 'env', 'env_nobuild', 'python').
    """

    env_classes = {}

    def get_env_type(python):
        env_type = env_classes.get(python)
        if env_type is None:
            cls = get_environment_class(conf, python)
            env_type = cls.tool_name
            env_classes[python] = env_type
        return env_type

    platform_keys = {
        ('environment_type', None): environment_type,
        ('sys_platform', None): sys.platform
    }

    # Parse requirement matrix
    parsed_matrix = _parse_matrix(conf.matrix)
    keys = list(parsed_matrix.keys())
    values = list(parsed_matrix.values())

    # Convert values to lists in the expected format
    values = [value if isinstance(value, list) else [value] for value in values]
    values = [[''] if value == [] else value for value in values]

    # Process excludes
    for python in pythons:
        empty_matrix = True

        # Cartesian product of everything
        all_keys = [('python', None)] + keys
        all_combinations = itertools.product([python], *values)

        for combination in all_combinations:
            target = dict(zip(all_keys, combination))
            target.update(platform_keys)

            if not environment_type:
                try:
                    target[('environment_type', None)] = get_env_type(target[('python', None)])
                except EnvironmentUnavailable as err:
                    log.warning(str(err))
                    continue

            for rule in conf.exclude:
                # check if all fields in the rule match
                rule = _parse_exclude_include_rule(rule)
                if match_rule(target, rule):
                    # rule matched
                    break
            else:
                # not excluded
                empty_matrix = False
                yield dict(item for item in zip(all_keys, combination)
                           if item[1] is not None)

        # If the user explicitly selected environment/python, yield it
        # even if matrix contains no packages to be installed
        if empty_matrix and explicit_selection:
            yield {('python', None): python}

    # Process includes, unless explicit selection
    if explicit_selection:
        return

    for include in conf.include:
        include = _parse_exclude_include_rule(include, is_include=True)

        # Platform keys in include statement act as matching rules
        target = dict(platform_keys)

        if not environment_type:
            try:
                target[('environment_type', None)] = get_env_type(include[('python', None)])
            except EnvironmentUnavailable as err:
                log.warning(str(err))
                continue

        rule = {}

        for key in platform_keys.keys():
            if key in include:
                rule[key] = include.pop(key)

        if match_rule(target, rule):
            # Prune empty keys
            for key in list(include.keys()):
                if include[key] is None:
                    include.pop(key)

            yield include


def _parse_matrix(matrix, bare_keys=()):
    """
    Parse 'matrix' and include/exclude rule configuration entries.

    It is in format::

         {"key_type1": {"key1": value1, "key2, value2, ...},
          ...,
          "nondict_key1": nondict_value1,
          ...}

    or in legacy format::

         {"key1": value1, ..., "nondict_key1": nondict_value1, ...}

    in which the key type is assumed to be "req".

    Parameters
    ----------
    matrix
        Configuration matrix or rule entry
    bare_keys : iterable
        Non-dictionary key values to store as is

    Returns
    -------
    parsed_matrix
        Dictionary {(key_type, key): value, ...}

    """
    matrix = dict(matrix)
    result = {}

    # Insert non-dict ("bare") keys first
    for key in bare_keys:
        if key in matrix:
            result[key, None] = matrix.pop(key)

    # Insert remaining matrix entries
    matrix_types = ('req', 'env', 'env_nobuild')
    if any(t in matrix for t in matrix_types):
        # New-style config
        matrices = []
        for t in matrix_types:
            submatrix = matrix.pop(t, {})
            matrices.append((t, submatrix))

        # Check if spurious keys left
        remaining_keys = tuple(matrix.keys())
        if remaining_keys:
            raise util.UserError('Unknown keys in "matrix" configuration: {}, expected: {}'.format(
                remaining_keys, matrix_types + tuple(bare_keys)))
    else:
        # Backward-compatibility for old-style config
        matrices = [('req', matrix)]

    # Convert values
    for t, m in matrices:
        for key, value in m.items():
            result[t, key] = value

    return result


def _parse_exclude_include_rule(rule, is_include=False):
    """
    Parse exclude/include rule by adding key types.

    Parameters
    ----------
    rule : dict
        Keys must be str, values must be str or None.
        The keys 'python', 'environment_type', 'sys_platform',
        are parsed specially and result to the corresponding key types.

    Returns
    -------
    rule : dict
        Dictionary of {(key_type, key): value, ...}
    """
    if is_include and 'python' not in rule:
        raise util.UserError(f"include rule '{rule}' does not specify Python version")

    bare_keys = ('python', 'environment_type', 'sys_platform')
    return _parse_matrix(rule, bare_keys)


def match_rule(target, rule):
    """
    Match rule to a target.

    Parameters
    ----------
    target : dict
        Dictionary containing [((key_type, key), value), ...].
    rule : dict
        Dictionary containing [((key_type, key), match), ...], to be matched
        to *target*. Match can be str specifying a regexp that must
        match target[key], or None. None matches either None
        or a missing key in *target*. If match is not None,
        and the key is missing in *target*, the rule does not match.
        The key_type must match exactly.

    Returns
    -------
    matched : bool
        Whether the rule matched. The rule matches if
        all keys match.

    """
    for key, value in rule.items():
        if value is None:
            if key in target and target[key] is not None:
                return False
        elif key not in target or target[key] is None:
            return False
        else:
            w = str(target[key])
            m = re.match(str(value), w)
            if m is None or m.end() != len(w):
                return False

    # rule matched
    return True


def get_env_name(tool_name, python, requirements, tagged_env_vars, build=False):
    """
    Get a name to uniquely identify an environment.

    Parameters
    ----------
    build : bool
        Whether to omit non-build environment variables.
        The canonical name of the environment is the name with build=False.
    """
    if tool_name:
        name = [tool_name]
    else:
        # Backward compatibility vs. result file names
        name = []

    name.append(f"py{python}")
    reqs = list(requirements.items())
    reqs.sort()
    for key, val in reqs:
        if val:
            name.append(''.join([key, val]))
        else:
            name.append(key)

    env_vars = _untag_env_vars(tagged_env_vars, build=build)

    for env_var, value in sorted(env_vars.items()):
        name.append(''.join([env_var, value]))

    return util.sanitize_filename('-'.join(name))


def _untag_env_vars(tagged_env_vars, build=False):
    vars = {}

    for (tag, key), value in tagged_env_vars.items():
        if not build or tag == 'build':
            vars[key] = value

    return vars


def get_environments(conf, env_specifiers, verbose=True):
    """
    Iterator returning `Environment` objects for all of the
    permutations of the given versions of Python and a matrix of
    requirements.

    Parameters
    ----------
    conf : dict
        asv configuration object
    env_specifiers : list of str
        List of environment specifiers, in the format
        'env_name:python_spec'. If *env_name* is empty, autodetect
        it. If *python_spec* is missing, use those listed in the
        configuration file. Alternatively, can be the name given
        by *Environment.name* if the environment is in the matrix.
    verbose : bool, optional
        Whether to display warnings about unknown environment types etc.

    """

    if not env_specifiers:
        all_environments = ()
        env_specifiers = [conf.environment_type]
        if not conf.environment_type and verbose:
            log.warning(
                "No `environment_type` specified in asv.conf.json. "
                "This will be required in the future.")
    else:
        all_environments = list(get_environments(conf, None, verbose=verbose))

    for env_spec in env_specifiers:
        env_name_found = False
        for env in all_environments:
            if env_spec == env.name:
                env_name_found = True
                yield env
                break
        if env_name_found:
            continue

        explicit_selection = False

        if env_spec and ':' in env_spec:
            env_type, python_spec = env_spec.split(':', 1)
            pythons = [python_spec]
            explicit_selection = True
        else:
            env_type = env_spec
            if env_type == "existing":
                explicit_selection = True
                pythons = ["same"]
            else:
                pythons = conf.pythons

        if env_type != "existing":
            requirements_iter = iter_matrix(env_type, pythons, conf,
                                            explicit_selection)
        else:
            # Ignore requirement matrix
            requirements_iter = [{('python', None): python} for python in pythons]

        for entries in requirements_iter:
            python, requirements, tagged_env_vars = _parse_matrix_entries(entries)

            try:
                if env_type:
                    cls = get_environment_class_by_name(env_type)
                else:
                    cls = get_environment_class(conf, python)

                yield cls(conf, python, requirements, tagged_env_vars)
            except EnvironmentUnavailable as err:
                if verbose:
                    log.warning(str(err))


def _parse_matrix_entries(entries):
    """
    Parse mixed requirement / environment variable matrix entries
    to requirements and tagged environment variables.
    """
    python = None
    requirements = {}
    tagged_env_vars = {}
    for (key_type, key), value in entries.items():
        if key_type == 'python':
            python = value
        elif key_type == 'env':
            tagged_env_vars[("build", key)] = value
        elif key_type == 'env_nobuild':
            tagged_env_vars[("nobuild", key)] = value
        elif key_type == 'req':
            requirements[key] = value
        else:
            # Shouldn't happen
            raise ValueError(f"Invalid matrix key type {key}")
    return python, requirements, tagged_env_vars


def get_environment_class(conf, python):
    """
    Get a matching environment type class.

    Parameters
    ----------
    conf : dict
        asv configuration object

    python : str
        Python version specifier.  Acceptable values depend on the
        Environment plugins installed but generally are:

        - 'X.Y': A Python version, in which case conda or virtualenv
          will be used to create a new environment.

        - 'python' or '/usr/bin/python': Search for the given
          executable on the search PATH, and use that.  It is assumed
          that all dependencies and the benchmarked project itself are
          already installed.
    """
    if python == 'same':
        return ExistingEnvironment

    # Try the subclasses in reverse order so custom plugins come first
    classes = list(util.iter_subclasses(Environment))[::-1]

    if conf.environment_type:
        cls = get_environment_class_by_name(conf.environment_type)
        classes.remove(cls)
        classes.insert(0, cls)

    for cls in classes:
        if cls.matches_python_fallback and cls.matches(python):
            return cls
    raise EnvironmentUnavailable(
        f"No way to create environment for python='{python}'")


def get_environment_class_by_name(environment_type):
    """
    Find the environment class with the given name.
    """
    for cls in util.iter_subclasses(Environment):
        if cls.tool_name == environment_type:
            return cls
    tool_names = [cls.tool_name for cls in util.iter_subclasses(Environment)]
    raise EnvironmentUnavailable(
        f"Unknown environment type '{environment_type}'. "
        f"Allowed values based on existing plugins are {tool_names}. "
        f"If you are trying to use `mamba`, you may need to install `conda-build`."
    )


def is_existing_only(environments):
    """
    Check if the list of environments only contains ExistingEnvironment
    """
    return all(isinstance(env, ExistingEnvironment) for env in environments)


class EnvironmentUnavailable(BaseException):
    pass


class Environment:
    """
    Manage a single environment -- a combination of a particular
    version of Python and a set of dependencies for the benchmarked
    project.
    """
    tool_name = None
    matches_python_fallback = True

    def __init__(self, conf, python, requirements, tagged_env_vars):
        """
        Get an environment for a given requirement matrix and
        Python version specifier.

        Parameters
        ----------
        conf : dict
            asv configuration object

        python : str
            A Python version specifier.  This is the same as passed to
            the `matches` method, and its exact meaning depends on the
            environment.

        requirements : dict (str -> str)
            Mapping from package names to versions

        tagged_env_vars : dict (tag, key) -> value
            Environment variables, tagged for build vs. non-build

        Raises
        ------
        EnvironmentUnavailable
            The environment for the given combination is not available.

        """
        self._env_dir = conf.env_dir
        self._repo_subdir = conf.repo_subdir
        self._install_timeout = conf.install_timeout  # gh-391
        self._default_benchmark_timeout = conf.default_benchmark_timeout # gh-973
        self._tagged_env_vars = tagged_env_vars
        self._path = os.path.abspath(os.path.join(
            self._env_dir, self.dir_name))
        self._project = conf.project

        self._is_setup = False

        self._cache = build_cache.BuildCache(conf, self._path)
        self._build_root = os.path.abspath(os.path.join(self._path, 'project'))

        self._requirements = requirements
        # These are needed for asv to build and run the project, not part of
        # benchmark name mangling
        self._base_requirements = {}
        # gh-1385
        self._base_requirements["pip+build"] = ""
        # gh-1314
        asv_runner_path = os.getenv("ASV_RUNNER_PATH", "")
        module_path = Path(asv_runner_path) / "asv_runner"

        # Check if the path points to a directory containing the "asv_runner" module
        if module_path.is_dir() and (module_path / "__init__.py").is_file():
            spec = importlib.util.spec_from_file_location("asv_runner",
                                                          module_path / "__init__.py")
            # Attempt to load the module
            asv_runner_module = importlib.util.module_from_spec(spec)
            try:
                spec.loader.exec_module(asv_runner_module)
                self._base_requirements["pip+asv_runner"] = asv_runner_path
            except Exception as e:
                self._base_requirements["pip+asv_runner"] = ""
                log.warning(f"Failed to load module from ASV_RUNNER_PATH: {e}")
        else:
            self._base_requirements["pip+asv_runner"] = ""
            if asv_runner_path:
                log.warning("ASV_RUNNER_PATH does not point"
                            "to a directory containing the 'asv_runner' module")
        if not util.ON_PYPY:
            # XXX: What if pypy installed asv tries to benchmark a cpython
            # python?
            self._base_requirements["pip+pympler"] = ""

        pyproject_path = Path.cwd() / "pyproject.toml"
        if pyproject_path.exists():
            with open(pyproject_path, "rb") as pyproject_file:
                pyproject_data = tomllib.load(pyproject_file)
                requires = pyproject_data.get("build-system", {}).get("requires", [])
                for requirement in requires:
                    self._base_requirements[f"pip+{requirement}"] = ""

        # Update the _base_requirements if needed
        for key in list(self._requirements.keys()):
            if key in self._base_requirements:
                self._base_requirements[key] = self._requirements[key]
                del self._requirements[key]

        self._build_command = conf.build_command
        self._install_command = conf.install_command
        self._uninstall_command = conf.uninstall_command

        self._global_env_vars = {}
        self._global_env_vars['ASV'] = 'true'
        self._global_env_vars['ASV_PROJECT'] = conf.project
        self._global_env_vars['ASV_CONF_DIR'] = os.path.abspath(os.getcwd())
        self._global_env_vars['ASV_ENV_NAME'] = self.name
        self._global_env_vars['ASV_ENV_DIR'] = self._path
        self._global_env_vars['ASV_ENV_TYPE'] = self.tool_name

        installed_commit_hash = self._get_installed_commit_hash()
        self._set_commit_hash(installed_commit_hash)

    def _set_commit_hash(self, commit_hash):
        if commit_hash is None:
            self._global_env_vars.pop('ASV_COMMIT', None)
        else:
            self._global_env_vars['ASV_COMMIT'] = commit_hash

    def _set_build_dirs(self, build_dir, cache_dir):
        if build_dir is None:
            self._global_env_vars.pop('ASV_BUILD_DIR', None)
        else:
            self._global_env_vars['ASV_BUILD_DIR'] = build_dir

        if cache_dir is None:
            self._global_env_vars.pop('ASV_BUILD_CACHE_DIR', None)
        else:
            self._global_env_vars['ASV_BUILD_CACHE_DIR'] = cache_dir

    def _set_installed_commit_hash(self, commit_hash):
        # Save status
        install_checksum = self._get_install_checksum()
        hash_file = os.path.join(self._path, 'asv-install-status.json')
        data = {'commit_hash': commit_hash, 'install_checksum': install_checksum}
        util.write_json(hash_file, data, api_version=1)

    def _get_installed_commit_hash(self):
        hash_file = os.path.join(self._path, 'asv-install-status.json')

        data = {}
        if os.path.isfile(hash_file):
            try:
                data = util.load_json(hash_file, api_version=1)
            except util.UserError:
                pass

        # If configuration changed, force reinstall
        install_checksum = self._get_install_checksum()
        if data.get('install_checksum', None) != install_checksum:
            return None

        return data.get('commit_hash', None)

    def _get_install_checksum(self):
        return [self._repo_subdir,
                self._install_timeout,
                self._project,
                self._build_command,
                self._install_command,
                self._uninstall_command]

    @property
    def installed_commit_hash(self):
        return self._get_installed_commit_hash()

    @classmethod
    def matches(self, python):
        """
        Returns `True` if this environment subclass can handle the
        given Python specifier.
        """
        return False

    @property
    def name(self):
        """
        Get a name to uniquely identify this environment.
        """
        return get_env_name(self.tool_name,
                            self._python,
                            self._requirements,
                            self._tagged_env_vars)

    @property
    def hashname(self):
        """
        Get a hash to uniquely identify this environment.
        """
        return hashlib.md5(self.name.encode('utf-8')).hexdigest()

    @property
    def dir_name(self):
        """
        Get the name of the directory where the environment resides.
        This is not necessarily unique, and may be shared across
        different environments.
        """
        name = get_env_name(self.tool_name,
                            self._python,
                            self._requirements,
                            self._tagged_env_vars,
                            build=True)
        return hashlib.md5(name.encode('utf-8')).hexdigest()

    @property
    def requirements(self):
        """Return the requirements"""
        return self._requirements

    @property
    def env_vars(self):
        """
        All environment variables configured in the matrix.
        """
        return _untag_env_vars(self._tagged_env_vars, build=False)

    @property
    def build_env_vars(self):
        """
        Build-time environment variables configured in the matrix.
        """
        return _untag_env_vars(self._tagged_env_vars, build=True)

    @property
    def python(self):
        return self._python

    def check_presence(self):
        """
        Check whether the environment already exists.
        """
        if not os.path.isdir(self._env_dir):
            return False

        try:
            info = self.load_info_file(self._path)
        except (util.UserError, OSError):
            return False

        expected_info = {
            'tool_name': self.tool_name,
            'python': self._python,
            'requirements': self._requirements,
            'build_env_vars': self.build_env_vars
        }

        if info != expected_info:
            return False

        for executable in ['pip', 'python']:
            try:
                self.find_executable(executable)
            except OSError:
                return False

        try:
            self.run_executable('python', ['-c', 'pass'])
        except (subprocess.CalledProcessError, OSError):
            return False

        return True

    def create(self):
        """
        Create the environment on disk.  If it doesn't exist, it is
        created.  Then, all of the requirements are installed into it.
        """
        if self._is_setup:
            return

        if not self.check_presence():
            if os.path.exists(self._path):
                util.long_path_rmtree(self._path)

            if not os.path.exists(self._env_dir):
                try:
                    os.makedirs(self._env_dir)
                except OSError:
                    # Environment.create may be called in parallel for
                    # environments with different self._path, but same
                    # self._env_dir. This causes a race condition for
                    # the above makedirs() call --- but not for the
                    # rest of the processing. Therefore, we will just
                    # ignore the error here, as things will fail at a
                    # later stage if there is really was a problem.
                    pass

            try:
                self._setup()
            except Exception:
                log.error(f"Failure creating environment for {self.name}")
                if os.path.exists(self._path):
                    util.long_path_rmtree(self._path)
                raise

        self.save_info_file(self._path)

        self._is_setup = True

    def _setup(self):
        """
        Implementation for setting up the environment.
        """
        raise NotImplementedError()

    def run(self, args, **kwargs):
        """
        Start up the environment's python executable with the given
        args.
        """
        raise NotImplementedError()

    def _interpolate_commands(self, commands):
        """
        Parse a command list with interpolated variables to a sequence of commands.

        Parameters
        ----------
        commands : {list of str}
            Commands to execute

        Returns
        -------
        run_commands : list of (cmd, env, return_codes, cwd)
            Parsed commands to run.

        """
        if not commands:
            return []

        if not isinstance(commands, list):
            commands = [commands]

        # All environment variables are available as interpolation variables,
        # lowercased without the prefix.
        kwargs = dict()
        for key, value in self._global_env_vars.items():
            if key == 'ASV':
                continue
            assert key.startswith('ASV_')
            interp_key = key[4:].lower()
            kwargs[interp_key] = value

        # There is an additional {wheel_file} interpolation variable
        if 'build_cache_dir' in kwargs:
            cache_dir = kwargs['build_cache_dir']

            if os.path.isdir(cache_dir):
                files = os.listdir(cache_dir)
                wheels = [fn for fn in files if fn.lower().endswith('.whl')]
                if len(wheels) == 1:
                    kwargs['wheel_file'] = os.path.join(cache_dir, wheels[0])

        # Interpolate, and raise useful error message if it fails
        return [util.interpolate_command(c, kwargs) for c in commands]

    def _interpolate_and_run_commands(self, commands, default_cwd, extra_env=None):
        interpolated = self._interpolate_commands(commands)

        for cmd, env, return_codes, cwd in interpolated:
            environ = dict(os.environ)
            if extra_env is not None:
                environ.update(extra_env)
            environ.update(env)
            if cwd is None:
                cwd = default_cwd
            self.run_executable(cmd[0], cmd[1:], timeout=self._install_timeout, cwd=cwd,
                                env=environ, valid_return_codes=return_codes)

    def checkout_project(self, repo, commit_hash):
        """
        Check out the working tree of the project at given commit hash
        """
        self._set_commit_hash(commit_hash)
        repo.checkout(self._build_root, commit_hash)

    def install_project(self, conf, repo, commit_hash):
        """
        Build and install the benchmarked project into the environment.
        Uninstalls any installed copy of the project first.
        """
        if self._repo_subdir:
            build_dir = os.path.join(self._build_root, self._repo_subdir)
        else:
            build_dir = self._build_root

        # Check first if anything needs to be done
        installed_commit_hash = self._get_installed_commit_hash()
        if installed_commit_hash == commit_hash:
            self._set_commit_hash(installed_commit_hash)
            self._set_build_dirs(None, None)
            return

        # Checkout first, so that uninstall can access build_dir
        # (for e.g. Makefiles)
        self.checkout_project(repo, commit_hash)
        self._set_build_dirs(build_dir, None)

        # Uninstall
        self._uninstall_project()

        # Build if not in cache
        cache_dir = self._cache.get_cache_dir(commit_hash)
        if cache_dir is not None:
            self._set_build_dirs(build_dir, cache_dir)
        else:
            cache_dir = self._cache.create_cache_dir(commit_hash)
            self._set_build_dirs(build_dir, cache_dir)
            self._build_project(repo, commit_hash, build_dir)

        # Install
        self._install_project(repo, commit_hash, build_dir)

        # Mark cached build as valid
        self._cache.finalize_cache_dir(commit_hash)

        # Mark installation as updated
        self._set_installed_commit_hash(commit_hash)

    def _install_project(self, repo, commit_hash, build_dir):
        """
        Run install commands
        """
        cmd = self._install_command
        if cmd is None:
            # Run pip via python -m pip, avoids shebang length limit on Linux.
            # Don't run it in build directory, because it may contain Python packages
            # that pip believes to be already installed.
            # --force-reinstall is needed since versions may not change between
            # asv runs (esp. for compare), e.g. gh-1421
            cmd = ["in-dir={env_dir} python -mpip install {wheel_file} --force-reinstall"]

        if cmd:
            commit_name = repo.get_decorated_hash(commit_hash, 8)
            log.info(f"Installing {commit_name} into {self.name}")
            self._interpolate_and_run_commands(cmd, default_cwd=build_dir,
                                               extra_env=self.build_env_vars)

    def _uninstall_project(self):
        """
        Run uninstall commands
        """
        # Mark installation invalid first
        self._set_installed_commit_hash(None)

        cmd = self._uninstall_command
        if cmd is None:
            # Run pip via python -m pip, avoids shebang length limit on Linux
            # pip uninstall may fail if not installed, so allow any exit code
            cmd = ['return-code=any python -mpip uninstall -y {project}']

        if cmd:
            log.info(f"Uninstalling from {self.name}")
            self._interpolate_and_run_commands(cmd, default_cwd=self._env_dir,
                                               extra_env=self.build_env_vars)

    def _build_project(self, repo, commit_hash, build_dir):
        """
        Run build commands
        """

        cmd = self._build_command

        if cmd is None:
            cmd = [
                "PIP_NO_BUILD_ISOLATION=0 python -m build",
                "python -mpip wheel -w {build_cache_dir} {build_dir}"
            ]

        if cmd:
            commit_name = repo.get_decorated_hash(commit_hash, 8)
            log.info(f"Building {commit_name} for {self.name}")
            self._interpolate_and_run_commands(cmd, default_cwd=build_dir,
                                               extra_env=self.build_env_vars)

    def can_install_project(self):
        """
        Return `True` if specific revisions of the benchmarked project
        can be installed into this environment.
        """
        return True

    def find_executable(self, executable):
        """
        Find an executable (eg. python, pip) in the environment.

        If not found, raises OSError
        """

        # Assume standard virtualenv/Conda layout
        if WIN:
            paths = [self._path,
                     os.path.join(self._path, 'Scripts'),
                     os.path.join(self._path, 'bin')]
        else:
            paths = [os.path.join(self._path, 'bin')]

        return util.which(executable, paths)

    def run_executable(self, executable, args, **kwargs):
        """
        Run a given executable (eg. python, pip) in the environment.
        """
        env = kwargs.pop("env", os.environ).copy()
        env.update(self._global_env_vars)

        # Insert bin dirs to PATH
        if "PATH" in env:
            paths = env["PATH"].split(os.pathsep)
        else:
            paths = []

        if WIN:
            subpaths = ['Library\\mingw-w64\\bin',
                        'Library\\bin',
                        'Library\\usr\\bin',
                        'Scripts']
            for sub in subpaths[::-1]:
                paths.insert(0, os.path.join(self._path, sub))
            paths.insert(0, self._path)
        else:
            paths.insert(0, os.path.join(self._path, "bin"))

        # Discard PYTHONPATH, which can easily break the environment
        # isolation
        if 'ASV_PYTHONPATH' in env:
            env['PYTHONPATH'] = env['ASV_PYTHONPATH']
            env.pop('ASV_PYTHONPATH', None)
        else:
            env.pop('PYTHONPATH', None)

        # When running pip, we need to set PIP_USER to false, as --user (which
        # may have been set from a pip config file) is incompatible with
        # virtualenvs.
        kwargs["env"] = dict(env,
                             PIP_USER=str("false"),
                             PATH=str(os.pathsep.join(paths)))
        exe = self.find_executable(executable)
        if kwargs.get("timeout", None) is None:
            kwargs["timeout"] = self._install_timeout
        return util.check_output([exe] + args, **kwargs)

    def load_info_file(self, path):
        path = os.path.join(path, 'asv-env-info.json')
        return util.load_json(path)

    def save_info_file(self, path):
        """
        Save a file with information about the environment into
        directory `path`.
        """
        path = os.path.join(path, 'asv-env-info.json')
        content = {
            'tool_name': self.tool_name,
            'python': self._python,
            'requirements': self._requirements,
            'build_env_vars': self.build_env_vars
        }
        util.write_json(path, content)


class ExistingEnvironment(Environment):
    tool_name = "existing"

    def __init__(self, conf, executable, requirements, tagged_env_vars):
        if executable == 'same':
            executable = sys.executable

        try:
            executable = os.path.abspath(util.which(executable))

            self._python = util.check_output(
                [executable,
                 '-c',
                 'import sys; '
                 'print(str(sys.version_info[0]) + "." + str(sys.version_info[1]))'
                 ]).strip()
        except (util.ProcessError, OSError):
            raise EnvironmentUnavailable()

        self._executable = executable
        self._requirements = {}

        super(ExistingEnvironment, self).__init__(conf,
                                                  executable,
                                                  requirements,
                                                  tagged_env_vars)
        self._global_env_vars.pop('ASV_ENV_DIR')

    @property
    def installed_commit_hash(self):
        return None

    @classmethod
    def matches(cls, python):
        if python == 'same':
            python = sys.executable

        try:
            util.which(python)
        except OSError:
            return False
        else:
            return True

    @property
    def name(self):
        return get_env_name(self.tool_name,
                            self._executable.replace(os.path.sep, '_'),
                            {},
                            self._tagged_env_vars)

    def check_presence(self):
        return True

    def create(self):
        pass

    def _setup(self):
        pass

    def install_project(self, conf, repo, commit_hash=None):
        pass

    def can_install_project(self):
        return False

    def run(self, args, **kwargs):
        log.debug(f"Running '{' '.join(args)}' in {self.name}")
        return util.check_output([
            self._executable] + args, **kwargs)
