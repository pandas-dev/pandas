import os
import unittest

from vortexa_utils.database import DatabaseFactory


class TestEngineFactory(unittest.TestCase):
    def test_create_factory(self):
        db_factory = DatabaseFactory()
        return db_factory

    def test_get_cert(self):
        db_factory = self.test_create_factory()
        cert_file = db_factory.fetch_cert()
        self.assertEqual(cert_file, db_factory.cert_file)
        assert os.path.isfile(cert_file)
