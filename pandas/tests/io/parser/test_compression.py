# -*- coding: utf-8 -*-

"""
Tests compressed data parsing functionality for all
of the parsers defined in parsers.py
"""

import os
import zipfile

import pytest

import pandas as pd
import pandas.util.testing as tm


@pytest.fixture(params=[True, False])
def buffer(request):
    return request.param


@pytest.fixture
def parser_and_data(all_parsers, csv1):
    parser = all_parsers

    with open(csv1, "rb") as f:
        data = f.read()
        expected = parser.read_csv(csv1)

    return parser, data, expected


@pytest.mark.parametrize("compression", ["zip", "infer", "zip2"])
def test_zip(parser_and_data, compression):
    parser, data, expected = parser_and_data

    with tm.ensure_clean("test_file.zip") as path:
        with zipfile.ZipFile(path, mode="w") as tmp:
            tmp.writestr("test_file", data)

        if compression == "zip2":
            with open(path, "rb") as f:
                result = parser.read_csv(f, compression="zip")
        else:
            result = parser.read_csv(path, compression=compression)

        tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("compression", ["zip", "infer"])
def test_zip_error_multiple_files(parser_and_data, compression):
    parser, data, expected = parser_and_data

    with tm.ensure_clean("combined_zip.zip") as path:
        inner_file_names = ["test_file", "second_file"]

        with zipfile.ZipFile(path, mode="w") as tmp:
            for file_name in inner_file_names:
                tmp.writestr(file_name, data)

        with pytest.raises(ValueError, match="Multiple files"):
            parser.read_csv(path, compression=compression)


def test_zip_error_no_files(parser_and_data):
    parser, _, _ = parser_and_data

    with tm.ensure_clean() as path:
        with zipfile.ZipFile(path, mode="w"):
            pass

        with pytest.raises(ValueError, match="Zero files"):
            parser.read_csv(path, compression="zip")


def test_zip_error_invalid_zip(parser_and_data):
    parser, _, _ = parser_and_data

    with tm.ensure_clean() as path:
        with open(path, "wb") as f:
            with pytest.raises(zipfile.BadZipfile,
                               match="File is not a zip file"):
                parser.read_csv(f, compression="zip")


@pytest.mark.parametrize("filename", [None, "test.{ext}"])
def test_compression(parser_and_data, compression_only, buffer, filename):
    parser, data, expected = parser_and_data
    compress_type = compression_only

    ext = "gz" if compress_type == "gzip" else compress_type
    filename = filename if filename is None else filename.format(ext=ext)

    if filename and buffer:
        pytest.skip("Cannot deduce compression from "
                    "buffer of compressed data.")

    with tm.ensure_clean(filename=filename) as path:
        tm.write_to_compressed(compress_type, path, data)
        compression = "infer" if filename else compress_type

        if buffer:
            with open(path, "rb") as f:
                result = parser.read_csv(f, compression=compression)
        else:
            result = parser.read_csv(path, compression=compression)

        tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("ext", [None, "gz", "bz2"])
def test_infer_compression(all_parsers, csv1, buffer, ext):
    # see gh-9770
    parser = all_parsers
    kwargs = dict(index_col=0, parse_dates=True)

    expected = parser.read_csv(csv1, **kwargs)
    kwargs["compression"] = "infer"

    if buffer:
        with open(csv1) as f:
            result = parser.read_csv(f, **kwargs)
    else:
        ext = "." + ext if ext else ""
        result = parser.read_csv(csv1 + ext, **kwargs)

    tm.assert_frame_equal(result, expected)


def test_compression_utf16_encoding(all_parsers, csv_dir_path):
    # see gh-18071
    parser = all_parsers
    path = os.path.join(csv_dir_path, "utf16_ex_small.zip")

    result = parser.read_csv(path, encoding="utf-16",
                             compression="zip", sep="\t")
    expected = pd.DataFrame({
        u"Country": [u"Venezuela", u"Venezuela"],
        u"Twitter": [u"Hugo Chávez Frías", u"Henrique Capriles R."]
    })

    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("invalid_compression", ["sfark", "bz3", "zipper"])
def test_invalid_compression(all_parsers, invalid_compression):
    parser = all_parsers
    compress_kwargs = dict(compression=invalid_compression)

    msg = ("Unrecognized compression "
           "type: {compression}".format(**compress_kwargs))

    with pytest.raises(ValueError, match=msg):
        parser.read_csv("test_file.zip", **compress_kwargs)
