import os
import unittest

import certifi


class TestCertifi(unittest.TestCase):
    def test_cabundle_exists(self) -> None:
        assert os.path.exists(certifi.where())

    def test_read_contents(self) -> None:
        content = certifi.contents()
        assert "-----BEGIN CERTIFICATE-----" in content

    def test_py_typed_exists(self) -> None:
        assert os.path.exists(
            os.path.join(os.path.dirname(certifi.__file__), 'py.typed')
        )
