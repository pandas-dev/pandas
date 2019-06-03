from pandas import read_csv
import os
import pytest


class TestFirstRowBom:

    @pytest.fixture(autouse=True)
    def setup_method(self, datapath):
        self.dirpath = datapath('io', 'parser', 'data')
        self.file_path = os.path.join(self.dirpath, 'bom_first_line.txt')

    def test_first_row_bom_python(self):
        assert read_csv(self.file_path,
                        delimiter='\t',
                        engine='python').shape == (8, 22)

    def test_first_row_bom_c(self):
        assert read_csv(self.file_path,
                        delimiter='\t').shape == (8, 22)
