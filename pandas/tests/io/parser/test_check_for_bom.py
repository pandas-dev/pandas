import pytest

from pandas import read_csv


def test_first_row_bom():
    assert read_csv('data/bom_first_line.txt',
                    delimiter='\t',
                    engine='python').shape == (8, 22)
    assert read_csv('data/bom_first_line.txt',
                    delimiter='\t').shape == (8, 22)
