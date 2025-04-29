from __future__ import annotations

import pytest

from xarray import DataArray, tutorial
from xarray.tests import assert_identical, network


@network
class TestLoadDataset:
    @pytest.fixture(autouse=True)
    def setUp(self):
        self.testfile = "tiny"

    def test_download_from_github(self, tmp_path) -> None:
        cache_dir = tmp_path / tutorial._default_cache_dir_name
        ds = tutorial.open_dataset(self.testfile, cache_dir=cache_dir).load()
        tiny = DataArray(range(5), name="tiny").to_dataset()
        assert_identical(ds, tiny)

    def test_download_from_github_load_without_cache(
        self, tmp_path, monkeypatch
    ) -> None:
        cache_dir = tmp_path / tutorial._default_cache_dir_name

        ds_nocache = tutorial.open_dataset(
            self.testfile, cache=False, cache_dir=cache_dir
        ).load()
        ds_cache = tutorial.open_dataset(self.testfile, cache_dir=cache_dir).load()
        assert_identical(ds_cache, ds_nocache)
