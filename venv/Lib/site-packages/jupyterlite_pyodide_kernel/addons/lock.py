"""a JupyterLite addon for customizing Pyodide lockfiles."""

from __future__ import annotations

import importlib.metadata
import json
import shutil
import urllib.parse
from fnmatch import fnmatch

from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path
from textwrap import indent
from typing import Any, TYPE_CHECKING

from traitlets import Unicode, Bool, Dict
import doit.tools

from jupyterlite_core.trait_types import TypedTuple
from jupyterlite_core.constants import JUPYTERLITE_JSON, UTF8

from ..utils import (
    iter_pep508_specs,
    list_wheels,
    normalize_names,
    get_wheel_name,
    patch_json_path,
    is_pyodide_wheel,
    wheel_to_pep508,
)
from ..constants import (
    PYODIDE_LOCK,
    PYODIDE_LOCK_STEM,
    PYODIDE_LOCK_DEFAULT_URL,
    PYODIDE_UV_WHEELS,
    LOAD_PYODIDE_OPTIONS,
    OPTION_LOCK_FILE_URL,
    OPTION_PACKAGES,
    NOARCH_WHL,
)

from ._base import _BaseAddon


if TYPE_CHECKING:
    from collections.abc import Iterator
    from packaging.utils import NormalizedName
    from pyodide_lock.spec import PackageSpec, PyodideLockSpec
    from jupyterlite_core.manager import LiteManager
    from .pyodide import PyodideAddon

    TTaskGenerator = Iterator[dict[str, Any]]
    TWheels = dict[NormalizedName, Path]


PYODIDE_LOCK_VERSION: str | None

try:
    PYODIDE_LOCK_VERSION = importlib.metadata.version(PYODIDE_LOCK_STEM)
except ImportError:  # pragma: no cover
    PYODIDE_LOCK_VERSION = None


class PyodideLockAddon(_BaseAddon):
    __all__ = ["status", "post_init", "post_build", "check"]

    # CLI
    aliases = {
        "pyodide-lock-url": "PyodideLockAddon.pyodide_lock_url",
        "pyodide-lock-wheels": "PyodideLockAddon.wheels",
        "pyodide-lock-constraints": "PyodideLockAddon.constraints",
        "pyodide-lock-lite-constraints": "PyodideLockAddon.lite_constraints_file",
        "pyodide-lock-specs": "PyodideLockAddon.specs",
        "pyodide-lock-excludes": "PyodideLockAddon.excludes_extra",
        "pyodide-lock-prefetch": "PyodideLockAddon.prefetch_extra",
    }

    flags = {
        "pyodide-lock": (
            {"PyodideLockAddon": {"enabled": True}},
            f"Use pyodide-lock and uv to customize {PYODIDE_LOCK}",
        ),
        "pyodide-lock-no-constrain-extensions": (
            {"PyodideLockAddon": {"constrain_extensions": False}},
            "Add ``LiteBuildConfig.federated_extensions`` to constraints",
        ),
    }

    # traits
    enabled: bool = Bool(
        default_value=False,
        help="whether Pyodide lockfile customization is enabled",
    ).tag(config=True)  # type: ignore[assignment]

    pyodide_lock_url: str | None = Unicode(
        help=f"URL of a remote {PYODIDE_LOCK}: {PYODIDE_LOCK_DEFAULT_URL}",
        allow_none=True,
    ).tag(config=True)  # type: ignore[assignment]

    wheels: tuple[str, ...] = TypedTuple(
        Unicode(), help=f"paths to local wheels or folders to include in {PYODIDE_LOCK}"
    ).tag(config=True)  # type: ignore[assignment]

    specs: tuple[str, ...] = TypedTuple(
        Unicode(),
        help=(
            f"PEP-508 specs for Python packages to include in {PYODIDE_LOCK};"
            " may include ``-r/--requirements`` and ``-g/--group``"
            " relative to ``lite_dir``"
        ),
    ).tag(config=True)  # type: ignore[assignment]

    constraints: tuple[str, ...] = TypedTuple(
        Unicode(),
        help=(
            "PEP-508 specs for Python packages to lock only if required in"
            f" {PYODIDE_LOCK}; may include ``-r/--requirements`` and "
            " ``-g/--group`` relative to ``lite_dir``"
        ),
    ).tag(config=True)  # type: ignore[assignment]

    constrain_extensions: bool = Bool(
        True,
        help=f"constrain wheels from ``LiteBuildConfig.federated_extensions`` in {PYODIDE_LOCK}",
    ).tag(config=True)  # type: ignore[assignment]

    excludes: tuple[str, ...] = TypedTuple(
        Unicode(),
        help=f"Python package names to exclude from {PYODIDE_LOCK}",
        default_value=[
            "jupyter-server",
            "jupyterlab",
            "notebook",
        ],
    ).tag(config=True)  # type: ignore[assignment]

    omit_output_wheels: tuple[str, ...] = TypedTuple(
        help="names of packages to omit from discovery in ``output_dir``",
        default_value=[
            "jupyterlab-widgets",
            "widgetsnbextension",
        ],
    ).tag(config=True)

    excludes_extra: tuple[str, ...] = TypedTuple(
        Unicode(),
        help=f"extra Python package names to exclude from {PYODIDE_LOCK}",
    ).tag(config=True)  # type: ignore[assignment]

    pyodide_lock_uv_options: dict[str, Any] = Dict(
        help="extra options to pass to ``pyodide_lock.uv_pip_compile.UvPipCompile``",
    ).tag(config=True)  # type: ignore[assignment]

    patches: dict[str, Any] = Dict(
        help="partial Pyodide lockfile to merge after the ``uv`` solve and URL rewrites",
    ).tag(config=True)  # type: ignore[assignment]

    prefetch: tuple[str, ...] = TypedTuple(
        Unicode(),
        default_value=[
            "ipykernel",
            "comm",
            "pyodide-kernel",
            "ipython",
        ],
        help=f"Python package names from {PYODIDE_LOCK} to prefetch while initializing Pyodide",
    ).tag(config=True)  # type: ignore[assignment]

    prefetch_extra: tuple[str] = TypedTuple(
        Unicode(),
        help=f"extra Python package names from {PYODIDE_LOCK} to prefetch while initializing Pyodide",
    ).tag(config=True)  # type: ignore[assignment]

    lite_constraints_file: str = Unicode(
        allow_none=True,
        help=(
            "path relative to ``lite_dir`` to a ``requirements.txt``-style"
            f" file with versions of all `none-any.whl` packages from {PYODIDE_LOCK}"
            " written if missing, otherwise added to ``constraints`` and left unchanged"
        ),
    ).tag(config=True)  # type: ignore[assignment]

    # properties
    @property
    def well_known_lock(self) -> Path:
        """a well-known path where pyodide-lock might be stored"""
        return self.manager.lite_dir / "static" / PYODIDE_LOCK_STEM

    @property
    def pyodide_addon(self) -> PyodideAddon:
        addons: dict[str, _BaseAddon] = self.manager._addons  # noqa: SLF001
        return addons["jupyterlite-pyodide-kernel-pyodide"]

    @property
    def output_lock(self) -> Path:
        return self.manager.output_dir / "static" / PYODIDE_LOCK_STEM / PYODIDE_LOCK

    @property
    def cache_dir(self) -> Path:
        """where ``pyodide-lock`` and ``uv`` stuff will go in the cache folder"""
        return self.manager.cache_dir / PYODIDE_LOCK_STEM

    @property
    def cached_remote_lock(self) -> Path:
        """where ``pyodide-lock`` and ``uv`` stuff will go in the cache folder"""
        return self.cache_dir / f"remote-{PYODIDE_LOCK}"

    @property
    def all_prefetch(self) -> list[NormalizedName]:
        """All packages to fetch while ``pyodide`` is initializing."""
        return normalize_names(*self.prefetch, *self.prefetch_extra)

    @property
    def all_excludes(self) -> list[NormalizedName]:
        """All packages to be excluded from the ``uv`` solve, and removed from the lock."""
        return normalize_names(*self.excludes, *self.excludes_extra)

    @property
    def all_specs(self) -> list[str]:
        """All PEP-508 specs, including ``requirements.txt``-style files."""
        return sorted(set(iter_pep508_specs([*self.specs], self.manager.lite_dir)))

    @property
    def all_constraints(self) -> list[str]:
        """All PEP-508 constraints, including ``requirements.txt``-style files."""
        constraints = [*self.constraints]
        lcp = self.lite_constraints_path
        if lcp and lcp.exists():
            constraints = [f"-r {self.lite_constraints_file}", *constraints]
        return sorted(set(iter_pep508_specs(constraints, self.manager.lite_dir)))

    @property
    def lite_constraints_path(self) -> Path | None:
        """Get the effective path of a ``constraints.txt`` to use in a lock, or write."""
        if not self.lite_constraints_file:
            return None
        return self.manager.lite_dir / self.lite_constraints_file

    @property
    def all_extra_uv_args(self) -> list[str]:
        """All arguments to inject for ``uv pip compile``."""
        args: list[str] = []
        if self.manager.source_date_epoch:
            iso = datetime.fromtimestamp(
                self.manager.source_date_epoch, tz=timezone.utc
            ).isoformat()
            args += ["--exclude-newer", iso]
        return args

    @property
    def status_info(self) -> str:
        """The status string, also used for task up-to-date checks."""
        lines = [
            f"pyodide-lock version:  {PYODIDE_LOCK_VERSION or 'not installed'}",
            f"pyodide-lock URL:      {self.effective_lock_url}",
            f"pyodide-lock options:  {self.pyodide_lock_uv_options}",
            "lock:",
            f" - wheels:       {self.wheels}",
            f" - specs:        {self.all_specs}",
            f" - constraints:  {self.all_constraints}",
            f" - excludes:     {self.all_excludes}",
            f" - uv args:      {self.all_extra_uv_args}",
            f" - patches:      {self.patches}",
            "runtime:",
            f" - prefetch packages:   {self.all_prefetch}",
            "output:",
            f" - lite constraints:    {self.lite_constraints_file}",
        ]
        return "\n".join(lines)

    @property
    def effective_lock_url(self) -> str | None:
        """Get an effective remote lock; empty for a local pyodide with no remote lock."""
        if not self.enabled:
            return None

        if self.pyodide_lock_url:
            return self.pyodide_lock_url

        pyodide = self.pyodide_addon

        if pyodide.pyodide_url or pyodide.well_known_pyodide.exists():
            return None

        return PYODIDE_LOCK_DEFAULT_URL

    # JupyterLite API task generators
    def status(self, manager: LiteManager) -> TTaskGenerator:
        """report on the status of pyodide and pyodide-lock"""
        if self.enabled:
            yield self.task(
                name="pyodide-lock",
                actions=[lambda: print(indent(self.status_info, "    "), flush=True)],
            )

    def post_init(self, manager: LiteManager) -> TTaskGenerator:
        """handle downloading of pyodide"""
        lock_url = self.effective_lock_url

        if lock_url:
            yield from self.post_init_cache_pyodide_lock(lock_url)

    def post_build(self, manager: LiteManager) -> TTaskGenerator:
        """configure jupyter-lite.json for Pyodide, potentially after updating a lockfile."""
        if not self.enabled:
            return

        pyodide_addon = self.pyodide_addon

        out = manager.output_dir
        jupyterlite_json = out / JUPYTERLITE_JSON

        patch_kwargs: dict[str, Path] = {}
        patch_uptodate = ""

        wheels_by_name = self.find_wheels_by_name()
        in_lock: Path | None = None
        lcp = self.lite_constraints_path

        candidates = [pyodide_addon.output_pyodide / PYODIDE_LOCK]
        if self.effective_lock_url:
            candidates += [self.cached_remote_lock]
        for candidate in candidates:
            if candidate.is_file():
                in_lock = candidate
                break

        if not (in_lock and in_lock.is_file()):  # pragma: no cover
            self.log.error(
                "A custom %s was requested, but no input lock found in: %s",
                PYODIDE_LOCK,
                candidates,
            )
            return

        yield self.task(
            name="lock:build",
            doc=f"build {PYODIDE_LOCK} with kernel, extension, and user-requested wheels",
            actions=[(self.post_build_lock, [in_lock, wheels_by_name])],
            file_dep=[
                in_lock,
                *wheels_by_name.values(),
                *([lcp] if lcp and lcp.is_file() else []),
            ],
            targets=[self.output_lock],
            uptodate=[doit.tools.config_changed(self.status_info)],
        )

        patch_uptodate += self.status_info

        yield self.task(
            name=f"patch:{JUPYTERLITE_JSON}",
            doc=f"ensure {JUPYTERLITE_JSON} includes pyodide-lock customizations",
            file_dep=[jupyterlite_json, *patch_kwargs.values()],
            actions=[(self.post_build_patch_jupyterlite_json, [jupyterlite_json])],
            uptodate=[doit.tools.config_changed(patch_uptodate)],
        )

        if lcp and not lcp.is_file():
            yield self.task(
                name=f"lite-constraints:{lcp}",
                doc=f"update lite constraints from {PYODIDE_LOCK}",
                actions=[(self.post_build_lite_constraints, [in_lock, lcp])],
                file_dep=[self.output_lock],
                targets=[lcp],
            )

    def check(self, manager: LiteManager) -> TTaskGenerator:
        """ensure the pyodide-lock configuration is sound"""

        if self.enabled:
            yield self.task(
                name="lock",
                doc=f"ensure {PYODIDE_LOCK} and local wheels are consistent",
                actions=[self.check_lock],
            )

    # task actions
    def post_init_cache_pyodide_lock(self, url: str) -> TTaskGenerator:
        """Cache a remote ``pyodide-lock.json`` if requested."""
        yield self.task(
            name=f"fetch:{PYODIDE_LOCK}",
            doc="fetch the pyodide lockfile",
            actions=[(self.fetch_one, [url, self.cached_remote_lock])],
            targets=[self.cached_remote_lock],
        )

    def post_build_lock(self, input_lock: Path, wheels_by_name: TWheels) -> bool:
        """Build a Pyodide lockfile with all kernel and user-requested wheels."""
        from pyodide_lock.uv_pip_compile import UvPipCompile

        tmp_lock = self.cache_dir / PYODIDE_LOCK
        self.copy_one(input_lock, tmp_lock)

        kwargs = deepcopy(self.pyodide_lock_uv_options)
        url_base: str | None = None
        lock_url = self.effective_lock_url
        if not lock_url and self.is_partial_pyodide(input_lock):
            lock_url = PYODIDE_LOCK_DEFAULT_URL
            self.log.warning(
                "Partial Pyodide distribution described in %s, using base URL %s",
                input_lock,
                lock_url,
            )
        if lock_url:
            url_base = lock_url.rsplit("/", 1)[0]
        extra_uv_args = [
            *kwargs.pop("extra_uv_args", []),
            *self.all_extra_uv_args,
        ]

        if self.constrain_extensions:
            extra_uv_args += [*self.get_labextension_uv_args()]

        upc = UvPipCompile(
            input_path=tmp_lock,
            output_path=tmp_lock.parent / f"patched-{tmp_lock.name}",
            input_base_url=url_base,
            wheels=[*wheels_by_name.values()],
            specs=[*self.all_specs],
            constraints=[*self.all_constraints],
            work_dir=self.cache_dir / "_work",
            wheel_dir=self.cache_dir / PYODIDE_UV_WHEELS,
            base_url_for_missing=url_base,
            excludes=[*map(str, self.all_excludes)],
            extra_uv_args=extra_uv_args,
            **kwargs,
        )
        spec = upc.update()

        # start with a clean folder
        shutil.rmtree(self.output_lock.parent, ignore_errors=True)
        self.output_lock.parent.mkdir(parents=True)

        # ensure wheels not already included in the output
        for pkg in spec.packages.values():
            self.ensure_local_spec(pkg, wheels_by_name)

        spec.to_json(path=self.output_lock, indent=2)
        patch_json_path(self.output_lock, self.patches)
        self.maybe_timestamp(self.output_lock)
        return True

    def post_build_lite_constraints(
        self, in_lock_path: Path, lite_constraints: Path
    ) -> bool:
        """Write out a requested constraint file of portable wheels not in base lock."""
        from packaging.utils import canonicalize_name

        in_packages, out_packages = [
            {
                canonicalize_name(pkg["name"]): pkg
                for pkg in json.loads(p.read_text(**UTF8))["packages"].values()
            }
            for p in [in_lock_path, self.output_lock]
        ]
        noarch_fn = f"*{NOARCH_WHL}"
        lines = ["# generated by jupyterlite-pyodide-kernel:PyodideLockAddon"]
        portable_prefixes = (f"../../static/{PYODIDE_LOCK_STEM}/", "http")

        for c_name, pkg in sorted(out_packages.items()):
            fn = str(pkg["file_name"])
            old_sha = in_packages.get(c_name, {}).get("sha256")
            is_portable = fnmatch(fn, noarch_fn) and fn.startswith(portable_prefixes)
            did_change = pkg["sha256"] != old_sha
            if is_portable and did_change:
                lines += [f"""{c_name} =={pkg["version"]}"""]
        lite_constraints.parent.mkdir(parents=True, exist_ok=True)
        lite_constraints.write_text("\n".join([*sorted(lines), ""]), **UTF8)
        return True

    def post_build_patch_jupyterlite_json(
        self,
        config_path: Path,
    ) -> None:
        """update jupyter-lite.json to use the custom Pyodide files"""
        out = self.manager.output_dir
        settings = self.get_pyodide_settings(config_path)

        lpo = settings.setdefault(LOAD_PYODIDE_OPTIONS, {})
        packages = normalize_names(*lpo.get(OPTION_PACKAGES, []), *self.all_prefetch)
        lpo.update(
            {
                OPTION_LOCK_FILE_URL: f"./{self.output_lock.relative_to(out).as_posix()}",
                OPTION_PACKAGES: packages,
            }
        )

        self.set_pyodide_settings(config_path, settings)

    def check_lock(self) -> bool:
        """Check the lock."""
        ok: dict[str, bool] = {}
        spec: PyodideLockSpec | None = None
        try:
            from pyodide_lock.spec import PyodideLockSpec

            spec = PyodideLockSpec.from_json(self.output_lock)
            ok["lock"] = True
        except Exception as err:
            self.log.error(
                "Failed to parse lock with pyodide-lock v%s: %s\n%s\n",
                PYODIDE_LOCK_VERSION,
                self.output_lock,
                err,
            )
            ok["lock"] = False
        if spec:
            c_names = {*normalize_names(*spec.packages)}
            for pkg in spec.packages.values():
                ok.update(self.check_package_spec(pkg, c_names))
        self.log.debug("Lock OK: %s", ok)
        return all(ok.values())

    # helpers
    def get_labextension_uv_args(self) -> list[str]:
        """Get packages from federated extensions to constrain the ``uv`` solve."""
        raw_wheels = []
        for ext in self.manager.federated_extensions:
            url = urllib.parse.urlparse(ext)
            if url.scheme and url.path.endswith(".whl") and is_pyodide_wheel(url.path):
                raw_wheels += [ext]
            elif not url.scheme:
                path = (self.manager.lite_dir / ext).absolute()
                if path.is_dir():
                    raw_wheels += [*map(str, list_wheels(path))]
                elif path.is_file() and is_pyodide_wheel(ext):
                    raw_wheels += [f"{path}"]

        ext_constraints: list[str] = []

        for raw_wheel in raw_wheels:
            pep508 = wheel_to_pep508(raw_wheel)
            if pep508:
                ext_constraints += [pep508]

        if not ext_constraints:
            return []

        tmp_constraints = self.cache_dir / "extension-constraints.txt"
        tmp_constraints.write_text("\n".join(sorted(ext_constraints)))

        return [f"--constraints={tmp_constraints.absolute()}"]

    def find_wheels_by_name(self) -> TWheels:
        """Gather a wheel per canonical name, first-in wins."""

        from packaging.utils import canonicalize_name

        wheels_by_name: TWheels = {}

        # add directly-requested wheels
        for wheel_str in self.wheels:
            wheel_or_dir = self.manager.lite_dir / wheel_str
            if wheel_or_dir.is_dir():
                for wheel_in_dir in list_wheels(wheel_or_dir):
                    self.add_wheel_by_name(wheel_in_dir, wheels_by_name)
            elif wheel_or_dir.is_file():
                self.add_wheel_by_name(wheel_or_dir, wheels_by_name)
            else:  # pragma: no cover
                self.log.warning("Wheel requested, but not found: %s", wheel_or_dir)

        # add wheels from well-known location
        for wheel in list_wheels(self.well_known_lock):
            self.add_wheel_by_name(wheel, wheels_by_name)

        omit = {*map(canonicalize_name, self.omit_output_wheels)}

        # add all wheels already in output, unless omitted
        for wheel in list_wheels(self.output_extensions, recursive=True):
            if get_wheel_name(wheel) in omit:
                continue
            self.add_wheel_by_name(wheel, wheels_by_name)

        return wheels_by_name

    def add_wheel_by_name(self, wheel: Path, wheels_by_name: TWheels) -> None:
        """Add a single wheel."""
        c_name = get_wheel_name(wheel)
        if c_name is None:  # pragma: no cover
            self.log.warning("[???] name cannot be found in wheel: %s", wheel)
            return
        if c_name in self.all_excludes:
            self.log.warning("[%s] local wheel excluded by name: %s", c_name, wheel)
            return
        if c_name in wheels_by_name:  # pragma: no cover
            self.log.warning(
                "[%s] local wheel already collected\n\tfrom %s\n\tdiscarding %s",
                c_name,
                wheels_by_name[c_name],
                wheel,
            )
            return
        wheels_by_name[c_name] = wheel

    def ensure_local_spec(self, pkg: PackageSpec, wheels_by_name: TWheels) -> None:
        """Ensure a ``PackageSpec`` points at a wheel in ``output_dir``."""
        url = urllib.parse.urlparse(pkg.file_name)

        if url.scheme or not url.path.startswith(PYODIDE_UV_WHEELS):
            return

        out = self.manager.output_dir
        just_name = Path(url.path).name
        rel_url: str | None = None
        in_wheel = wheels_by_name.get(normalize_names(pkg.name)[0])
        out_wheel: Path | None = None

        if in_wheel and out in in_wheel.parents:
            out_wheel = in_wheel
            rel_url = out_wheel.relative_to(out).as_posix()
        else:
            uv_wheel = self.cache_dir / PYODIDE_UV_WHEELS / just_name
            out_wheel = self.output_lock.parent / just_name
            if not uv_wheel.is_file():  # pragma: no cover
                self.log.error(
                    "[%s] wheel could not be found from spec: %s", pkg.name, pkg
                )
                raise FileNotFoundError(uv_wheel)
            self.copy_one(uv_wheel, out_wheel)
            rel_url = f"static/{PYODIDE_LOCK_STEM}/{just_name}"

        if not (rel_url and out_wheel and out_wheel.exists()):  # pragma: no cover
            msg = f"Don't know what to do with {just_name} from {pkg}: {out_wheel}"
            raise NotImplementedError(msg)

        # build a relative path from the location the pyodide runtime loads from
        pkg.file_name = f"../../{rel_url}"

    def check_package_spec(
        self, pkg: PackageSpec, c_names: set[NormalizedName]
    ) -> dict[str, bool]:
        """Verify a single package."""
        from packaging.utils import canonicalize_name

        name = canonicalize_name(pkg.name)
        url = urllib.parse.urlparse(pkg.file_name)

        is_ok: dict[str, bool] = {}

        for dep_name in pkg.depends:
            c_dep = canonicalize_name(dep_name)
            dep_ok = is_ok[f"{pkg.name}:depends:{c_dep}"] = c_dep in c_names
            if not dep_ok:
                self.log.error("[%s] missing dependency: %s", name, dep_name)

        if not url.scheme:
            path = self.output_lock.parent / url.path
            file_ok = is_ok[name] = path.is_file()

            if not file_ok:
                self.log.error("[%s] missing wheel: %s", name, path)

        return is_ok

    def is_partial_pyodide(self, lockfile: Path) -> bool:
        """Get whether a lockfile is missing a local file."""
        from pyodide_lock.spec import PyodideLockSpec

        dist_dir = lockfile.parent
        lock = PyodideLockSpec.from_json(lockfile)
        for pkg in lock.packages.values():
            url = urllib.parse.urlparse(pkg.file_name)
            if url.scheme:
                continue
            if not (dist_dir / url.path).is_file():
                return True
        return False
