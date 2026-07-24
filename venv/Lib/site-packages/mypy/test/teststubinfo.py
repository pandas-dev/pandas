from __future__ import annotations

import unittest

from mypy.stubinfo import non_bundled_packages_flat, stub_distribution_name


class TestStubInfo(unittest.TestCase):
    def test_stub_distribution_name(self) -> None:
        assert stub_distribution_name("foobar_asdf") is None
        assert stub_distribution_name("pycurl") == "types-pycurl"
        assert stub_distribution_name("psutil") == "types-psutil"
        assert stub_distribution_name("sassutils") == "types-libsass"
        assert stub_distribution_name("google.cloud.ndb") == "types-google-cloud-ndb"
        assert stub_distribution_name("google.cloud.ndb.submodule") == "types-google-cloud-ndb"
        assert stub_distribution_name("google.cloud.unknown") is None
        assert stub_distribution_name("google.protobuf") == "types-protobuf"
        assert stub_distribution_name("google.protobuf.submodule") == "types-protobuf"
        assert stub_distribution_name("google") is None

    def test_period_in_top_level(self) -> None:
        for packages in non_bundled_packages_flat:
            for top_level_module in packages:
                assert "." not in top_level_module
