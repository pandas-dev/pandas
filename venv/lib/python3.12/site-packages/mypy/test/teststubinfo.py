from __future__ import annotations

import unittest

from mypy.stubinfo import (
    is_module_from_legacy_bundled_package,
    legacy_bundled_packages,
    non_bundled_packages_flat,
    stub_distribution_name,
)


class TestStubInfo(unittest.TestCase):
    def test_is_legacy_bundled_packages(self) -> None:
        assert not is_module_from_legacy_bundled_package("foobar_asdf")
        assert not is_module_from_legacy_bundled_package("PIL")
        assert is_module_from_legacy_bundled_package("pycurl")
        assert is_module_from_legacy_bundled_package("dateparser")

    def test_stub_distribution_name(self) -> None:
        assert stub_distribution_name("foobar_asdf") is None
        assert stub_distribution_name("pycurl") == "types-pycurl"
        assert stub_distribution_name("bs4") == "types-beautifulsoup4"
        assert stub_distribution_name("google.cloud.ndb") == "types-google-cloud-ndb"
        assert stub_distribution_name("google.cloud.ndb.submodule") == "types-google-cloud-ndb"
        assert stub_distribution_name("google.cloud.unknown") is None
        assert stub_distribution_name("google.protobuf") == "types-protobuf"
        assert stub_distribution_name("google.protobuf.submodule") == "types-protobuf"
        assert stub_distribution_name("google") is None

    def test_period_in_top_level(self) -> None:
        for packages in (non_bundled_packages_flat, legacy_bundled_packages):
            for top_level_module in packages:
                assert "." not in top_level_module
