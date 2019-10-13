import pandas as pd


def test_to_string_encoding(float_frame,):
    # GH 28766
    path = "test_to_string_file"
    float_frame.to_string(path, encoding="gbk")
    with open(str(path), "r", encoding="gbk") as f:
        assert float_frame.to_string() == f.read()
