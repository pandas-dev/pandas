from __future__ import annotations

from xarray import DataArray, DataTree, tutorial
from xarray.testing import assert_identical
from xarray.tests import network


@network
class TestLoadDataset:
    def test_download_from_github(self, tmp_path) -> None:
        cache_dir = tmp_path / tutorial._default_cache_dir_name
        ds = tutorial.open_dataset("tiny", cache_dir=cache_dir).load()
        tiny = DataArray(range(5), name="tiny").to_dataset()
        assert_identical(ds, tiny)

    def test_download_from_github_load_without_cache(
        self, tmp_path, monkeypatch
    ) -> None:
        cache_dir = tmp_path / tutorial._default_cache_dir_name

        ds_nocache = tutorial.open_dataset(
            "tiny", cache=False, cache_dir=cache_dir
        ).load()
        ds_cache = tutorial.open_dataset("tiny", cache_dir=cache_dir).load()
        assert_identical(ds_cache, ds_nocache)


@network
class TestLoadDataTree:
    def test_download_from_github(self, tmp_path) -> None:
        cache_dir = tmp_path / tutorial._default_cache_dir_name
        ds = tutorial.open_datatree("tiny", cache_dir=cache_dir).load()
        tiny = DataTree.from_dict({"/": DataArray(range(5), name="tiny").to_dataset()})
        assert_identical(ds, tiny)

    def test_download_from_github_load_without_cache(
        self, tmp_path, monkeypatch
    ) -> None:
        cache_dir = tmp_path / tutorial._default_cache_dir_name

        ds_nocache = tutorial.open_datatree(
            "tiny", cache=False, cache_dir=cache_dir
        ).load()
        ds_cache = tutorial.open_datatree("tiny", cache_dir=cache_dir).load()
        assert_identical(ds_cache, ds_nocache)
