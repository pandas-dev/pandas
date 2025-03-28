from __future__ import annotations

import sys
from importlib.metadata import EntryPoint, EntryPoints
from unittest import mock

import pytest

from xarray.backends import common, plugins
from xarray.tests import (
    has_h5netcdf,
    has_netCDF4,
    has_pydap,
    has_scipy,
    has_zarr,
)

# Do not import list_engines here, this will break the lazy tests

importlib_metadata_mock = "importlib.metadata"


class DummyBackendEntrypointArgs(common.BackendEntrypoint):
    def open_dataset(filename_or_obj, *args):  # type: ignore[override]
        pass


class DummyBackendEntrypointKwargs(common.BackendEntrypoint):
    def open_dataset(filename_or_obj, **kwargs):  # type: ignore[override]
        pass


class DummyBackendEntrypoint1(common.BackendEntrypoint):
    def open_dataset(self, filename_or_obj, *, decoder):  # type: ignore[override]
        pass


class DummyBackendEntrypoint2(common.BackendEntrypoint):
    def open_dataset(self, filename_or_obj, *, decoder):  # type: ignore[override]
        pass


@pytest.fixture
def dummy_duplicated_entrypoints():
    specs = [
        ["engine1", "xarray.tests.test_plugins:backend_1", "xarray.backends"],
        ["engine1", "xarray.tests.test_plugins:backend_2", "xarray.backends"],
        ["engine2", "xarray.tests.test_plugins:backend_1", "xarray.backends"],
        ["engine2", "xarray.tests.test_plugins:backend_2", "xarray.backends"],
    ]
    eps = [EntryPoint(name, value, group) for name, value, group in specs]
    return eps


@pytest.mark.filterwarnings("ignore:Found")
def test_remove_duplicates(dummy_duplicated_entrypoints) -> None:
    with pytest.warns(RuntimeWarning):
        entrypoints = plugins.remove_duplicates(dummy_duplicated_entrypoints)
    assert len(entrypoints) == 2


def test_broken_plugin() -> None:
    broken_backend = EntryPoint(
        "broken_backend",
        "xarray.tests.test_plugins:backend_1",
        "xarray.backends",
    )
    with pytest.warns(RuntimeWarning) as record:
        _ = plugins.build_engines(EntryPoints([broken_backend]))
    assert len(record) == 1
    message = str(record[0].message)
    assert "Engine 'broken_backend'" in message


def test_remove_duplicates_warnings(dummy_duplicated_entrypoints) -> None:
    with pytest.warns(RuntimeWarning) as record:
        _ = plugins.remove_duplicates(dummy_duplicated_entrypoints)

    assert len(record) == 2
    message0 = str(record[0].message)
    message1 = str(record[1].message)
    assert "entrypoints" in message0
    assert "entrypoints" in message1


@mock.patch(
    f"{importlib_metadata_mock}.EntryPoint.load", mock.MagicMock(return_value=None)
)
def test_backends_dict_from_pkg() -> None:
    specs = [
        ["engine1", "xarray.tests.test_plugins:backend_1", "xarray.backends"],
        ["engine2", "xarray.tests.test_plugins:backend_2", "xarray.backends"],
    ]
    entrypoints = [EntryPoint(name, value, group) for name, value, group in specs]
    engines = plugins.backends_dict_from_pkg(entrypoints)
    assert len(engines) == 2
    assert engines.keys() == {"engine1", "engine2"}


def test_set_missing_parameters() -> None:
    backend_1 = DummyBackendEntrypoint1
    backend_2 = DummyBackendEntrypoint2
    backend_2.open_dataset_parameters = ("filename_or_obj",)
    engines = {"engine_1": backend_1, "engine_2": backend_2}
    plugins.set_missing_parameters(engines)

    assert len(engines) == 2
    assert backend_1.open_dataset_parameters == ("filename_or_obj", "decoder")
    assert backend_2.open_dataset_parameters == ("filename_or_obj",)

    backend_kwargs = DummyBackendEntrypointKwargs
    backend_kwargs.open_dataset_parameters = ("filename_or_obj", "decoder")
    plugins.set_missing_parameters({"engine": backend_kwargs})
    assert backend_kwargs.open_dataset_parameters == ("filename_or_obj", "decoder")

    backend_args = DummyBackendEntrypointArgs
    backend_args.open_dataset_parameters = ("filename_or_obj", "decoder")
    plugins.set_missing_parameters({"engine": backend_args})
    assert backend_args.open_dataset_parameters == ("filename_or_obj", "decoder")

    # reset
    backend_1.open_dataset_parameters = None
    backend_1.open_dataset_parameters = None
    backend_kwargs.open_dataset_parameters = None
    backend_args.open_dataset_parameters = None


def test_set_missing_parameters_raise_error() -> None:
    backend = DummyBackendEntrypointKwargs
    with pytest.raises(TypeError):
        plugins.set_missing_parameters({"engine": backend})

    backend_args = DummyBackendEntrypointArgs
    with pytest.raises(TypeError):
        plugins.set_missing_parameters({"engine": backend_args})


@mock.patch(
    f"{importlib_metadata_mock}.EntryPoint.load",
    mock.MagicMock(return_value=DummyBackendEntrypoint1),
)
def test_build_engines() -> None:
    dummy_pkg_entrypoint = EntryPoint(
        "dummy", "xarray.tests.test_plugins:backend_1", "xarray_backends"
    )
    backend_entrypoints = plugins.build_engines(EntryPoints([dummy_pkg_entrypoint]))

    assert isinstance(backend_entrypoints["dummy"], DummyBackendEntrypoint1)
    assert backend_entrypoints["dummy"].open_dataset_parameters == (
        "filename_or_obj",
        "decoder",
    )


@mock.patch(
    f"{importlib_metadata_mock}.EntryPoint.load",
    mock.MagicMock(return_value=DummyBackendEntrypoint1),
)
def test_build_engines_sorted() -> None:
    dummy_pkg_entrypoints = EntryPoints(
        [
            EntryPoint(
                "dummy2", "xarray.tests.test_plugins:backend_1", "xarray.backends"
            ),
            EntryPoint(
                "dummy1", "xarray.tests.test_plugins:backend_1", "xarray.backends"
            ),
        ]
    )
    backend_entrypoints = list(plugins.build_engines(dummy_pkg_entrypoints))

    indices = []
    for be in plugins.STANDARD_BACKENDS_ORDER:
        try:
            index = backend_entrypoints.index(be)
            backend_entrypoints.pop(index)
            indices.append(index)
        except ValueError:
            pass

    assert set(indices) < {0, -1}
    assert list(backend_entrypoints) == sorted(backend_entrypoints)


@mock.patch(
    "xarray.backends.plugins.list_engines",
    mock.MagicMock(return_value={"dummy": DummyBackendEntrypointArgs()}),
)
def test_no_matching_engine_found() -> None:
    with pytest.raises(ValueError, match=r"did not find a match in any"):
        plugins.guess_engine("not-valid")

    with pytest.raises(ValueError, match=r"found the following matches with the input"):
        plugins.guess_engine("foo.nc")


@mock.patch(
    "xarray.backends.plugins.list_engines",
    mock.MagicMock(return_value={}),
)
def test_engines_not_installed() -> None:
    with pytest.raises(ValueError, match=r"xarray is unable to open"):
        plugins.guess_engine("not-valid")

    with pytest.raises(ValueError, match=r"found the following matches with the input"):
        plugins.guess_engine("foo.nc")


def test_lazy_import() -> None:
    """Test that some modules are imported in a lazy manner.

    When importing xarray these should not be imported as well.
    Only when running code for the first time that requires them.
    """
    deny_list = [
        "cubed",
        "cupy",
        # "dask",  # TODO: backends.locks is not lazy yet :(
        "dask.array",
        "dask.distributed",
        "flox",
        "h5netcdf",
        "matplotlib",
        "nc_time_axis",
        "netCDF4",
        "numbagg",
        "pint",
        "pydap",
        "scipy",
        "sparse",
        "zarr",
    ]
    # ensure that none of the above modules has been imported before
    modules_backup = {}
    for pkg in list(sys.modules.keys()):
        for mod in deny_list + ["xarray"]:
            if pkg.startswith(mod):
                modules_backup[pkg] = sys.modules[pkg]
                del sys.modules[pkg]
                break

    try:
        import xarray  # noqa: F401
        from xarray.backends import list_engines

        list_engines()

        # ensure that none of the modules that are supposed to be
        # lazy loaded are loaded when importing xarray
        is_imported = set()
        for pkg in sys.modules:
            for mod in deny_list:
                if pkg.startswith(mod):
                    is_imported.add(mod)
                    break
        assert (
            len(is_imported) == 0
        ), f"{is_imported} have been imported but should be lazy"

    finally:
        # restore original
        sys.modules.update(modules_backup)


def test_list_engines() -> None:
    from xarray.backends import list_engines

    engines = list_engines()
    assert list_engines.cache_info().currsize == 1

    assert ("scipy" in engines) == has_scipy
    assert ("h5netcdf" in engines) == has_h5netcdf
    assert ("netcdf4" in engines) == has_netCDF4
    assert ("pydap" in engines) == has_pydap
    assert ("zarr" in engines) == has_zarr
    assert "store" in engines


def test_refresh_engines() -> None:
    from xarray.backends import list_engines, refresh_engines

    EntryPointMock1 = mock.MagicMock()
    EntryPointMock1.name = "test1"
    EntryPointMock1.load.return_value = DummyBackendEntrypoint1

    return_value = EntryPoints([EntryPointMock1])

    with mock.patch("xarray.backends.plugins.entry_points", return_value=return_value):
        list_engines.cache_clear()
        engines = list_engines()
    assert "test1" in engines
    assert isinstance(engines["test1"], DummyBackendEntrypoint1)

    EntryPointMock2 = mock.MagicMock()
    EntryPointMock2.name = "test2"
    EntryPointMock2.load.return_value = DummyBackendEntrypoint2

    return_value2 = EntryPoints([EntryPointMock2])

    with mock.patch("xarray.backends.plugins.entry_points", return_value=return_value2):
        refresh_engines()
        engines = list_engines()
    assert "test1" not in engines
    assert "test2" in engines
    assert isinstance(engines["test2"], DummyBackendEntrypoint2)

    # reset to original
    refresh_engines()
