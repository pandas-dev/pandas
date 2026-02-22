# Tests for the `str` accessor are derived from the original
# pandas string accessor tests.

# For reference, here is a copy of the pandas copyright notice:

# (c) 2011-2012, Lambda Foundry, Inc. and PyData Development Team
# All rights reserved.

# Copyright (c) 2008-2011 AQR Capital Management, LLC
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:

#     * Redistributions of source code must retain the above copyright
#        notice, this list of conditions and the following disclaimer.

#     * Redistributions in binary form must reproduce the above
#        copyright notice, this list of conditions and the following
#        disclaimer in the documentation and/or other materials provided
#        with the distribution.

#     * Neither the name of the copyright holder nor the names of any
#        contributors may be used to endorse or promote products derived
#        from this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
from __future__ import annotations

import re
from collections.abc import Callable

import numpy as np
import pytest

import xarray as xr
from xarray.tests import assert_equal, assert_identical, requires_dask


@pytest.fixture(
    params=[pytest.param(np.str_, id="str"), pytest.param(np.bytes_, id="bytes")]
)
def dtype(request):
    return request.param


@requires_dask
def test_dask() -> None:
    import dask.array as da

    arr = da.from_array(["a", "b", "c"], chunks=-1)
    xarr = xr.DataArray(arr)

    result = xarr.str.len().compute()
    expected = xr.DataArray([1, 1, 1])
    assert result.dtype == expected.dtype
    assert_equal(result, expected)


def test_count(dtype) -> None:
    values = xr.DataArray(["foo", "foofoo", "foooofooofommmfoo"]).astype(dtype)
    pat_str = dtype(r"f[o]+")
    pat_re = re.compile(pat_str)

    result_str = values.str.count(pat_str)
    result_re = values.str.count(pat_re)

    expected = xr.DataArray([1, 2, 4])

    assert result_str.dtype == expected.dtype
    assert result_re.dtype == expected.dtype
    assert_equal(result_str, expected)
    assert_equal(result_re, expected)


def test_count_broadcast(dtype) -> None:
    values = xr.DataArray(["foo", "foofoo", "foooofooofommmfoo"]).astype(dtype)
    pat_str = np.array([r"f[o]+", r"o", r"m"]).astype(dtype)
    pat_re = np.array([re.compile(x) for x in pat_str])

    result_str = values.str.count(pat_str)
    result_re = values.str.count(pat_re)

    expected = xr.DataArray([1, 4, 3])

    assert result_str.dtype == expected.dtype
    assert result_re.dtype == expected.dtype
    assert_equal(result_str, expected)
    assert_equal(result_re, expected)


def test_contains(dtype) -> None:
    values = xr.DataArray(["Foo", "xYz", "fOOomMm__fOo", "MMM_"]).astype(dtype)

    # case insensitive using regex
    pat = values.dtype.type("FOO|mmm")
    result = values.str.contains(pat, case=False)
    expected = xr.DataArray([True, False, True, True])
    assert result.dtype == expected.dtype
    assert_equal(result, expected)
    result = values.str.contains(re.compile(pat, flags=re.IGNORECASE))
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    # case sensitive using regex
    pat = values.dtype.type("Foo|mMm")
    result = values.str.contains(pat)
    expected = xr.DataArray([True, False, True, False])
    assert result.dtype == expected.dtype
    assert_equal(result, expected)
    result = values.str.contains(re.compile(pat))
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    # case insensitive without regex
    result = values.str.contains("foo", regex=False, case=False)
    expected = xr.DataArray([True, False, True, False])
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    # case sensitive without regex
    result = values.str.contains("fO", regex=False, case=True)
    expected = xr.DataArray([False, False, True, False])
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    # regex regex=False
    pat_re = re.compile("(/w+)")
    with pytest.raises(
        ValueError,
        match=r"Must use regular expression matching for regular expression object.",
    ):
        values.str.contains(pat_re, regex=False)


def test_contains_broadcast(dtype) -> None:
    values = xr.DataArray(["Foo", "xYz", "fOOomMm__fOo", "MMM_"], dims="X").astype(
        dtype
    )
    pat_str = xr.DataArray(["FOO|mmm", "Foo", "MMM"], dims="Y").astype(dtype)
    pat_re = xr.DataArray([re.compile(x) for x in pat_str.data], dims="Y")

    # case insensitive using regex
    result = values.str.contains(pat_str, case=False)
    expected = xr.DataArray(
        [
            [True, True, False],
            [False, False, False],
            [True, True, True],
            [True, False, True],
        ],
        dims=["X", "Y"],
    )
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    # case sensitive using regex
    result = values.str.contains(pat_str)
    expected = xr.DataArray(
        [
            [False, True, False],
            [False, False, False],
            [False, False, False],
            [False, False, True],
        ],
        dims=["X", "Y"],
    )
    assert result.dtype == expected.dtype
    assert_equal(result, expected)
    result = values.str.contains(pat_re)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    # case insensitive without regex
    result = values.str.contains(pat_str, regex=False, case=False)
    expected = xr.DataArray(
        [
            [False, True, False],
            [False, False, False],
            [False, True, True],
            [False, False, True],
        ],
        dims=["X", "Y"],
    )
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    # case insensitive with regex
    result = values.str.contains(pat_str, regex=False, case=True)
    expected = xr.DataArray(
        [
            [False, True, False],
            [False, False, False],
            [False, False, False],
            [False, False, True],
        ],
        dims=["X", "Y"],
    )
    assert result.dtype == expected.dtype
    assert_equal(result, expected)


def test_starts_ends_with(dtype) -> None:
    values = xr.DataArray(["om", "foo_nom", "nom", "bar_foo", "foo"]).astype(dtype)

    result = values.str.startswith("foo")
    expected = xr.DataArray([False, True, False, False, True])
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    result = values.str.endswith("foo")
    expected = xr.DataArray([False, False, False, True, True])
    assert result.dtype == expected.dtype
    assert_equal(result, expected)


def test_starts_ends_with_broadcast(dtype) -> None:
    values = xr.DataArray(
        ["om", "foo_nom", "nom", "bar_foo", "foo_bar"], dims="X"
    ).astype(dtype)
    pat = xr.DataArray(["foo", "bar"], dims="Y").astype(dtype)

    result = values.str.startswith(pat)
    expected = xr.DataArray(
        [[False, False], [True, False], [False, False], [False, True], [True, False]],
        dims=["X", "Y"],
    )
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    result = values.str.endswith(pat)
    expected = xr.DataArray(
        [[False, False], [False, False], [False, False], [True, False], [False, True]],
        dims=["X", "Y"],
    )
    assert result.dtype == expected.dtype
    assert_equal(result, expected)


def test_case_bytes() -> None:
    value = xr.DataArray(["SOme wOrd"]).astype(np.bytes_)

    exp_capitalized = xr.DataArray(["Some word"]).astype(np.bytes_)
    exp_lowered = xr.DataArray(["some word"]).astype(np.bytes_)
    exp_swapped = xr.DataArray(["soME WoRD"]).astype(np.bytes_)
    exp_titled = xr.DataArray(["Some Word"]).astype(np.bytes_)
    exp_uppered = xr.DataArray(["SOME WORD"]).astype(np.bytes_)

    res_capitalized = value.str.capitalize()
    res_lowered = value.str.lower()
    res_swapped = value.str.swapcase()
    res_titled = value.str.title()
    res_uppered = value.str.upper()

    assert res_capitalized.dtype == exp_capitalized.dtype
    assert res_lowered.dtype == exp_lowered.dtype
    assert res_swapped.dtype == exp_swapped.dtype
    assert res_titled.dtype == exp_titled.dtype
    assert res_uppered.dtype == exp_uppered.dtype

    assert_equal(res_capitalized, exp_capitalized)
    assert_equal(res_lowered, exp_lowered)
    assert_equal(res_swapped, exp_swapped)
    assert_equal(res_titled, exp_titled)
    assert_equal(res_uppered, exp_uppered)


def test_case_str() -> None:
    # This string includes some unicode characters
    # that are common case management corner cases
    value = xr.DataArray(["SOme wOrd Ǆ ß ᾛ ΣΣ ﬃ⁵Å Ç Ⅰ"]).astype(np.str_)

    exp_capitalized = xr.DataArray(["Some word ǆ ß ᾓ σς ﬃ⁵å ç ⅰ"]).astype(np.str_)
    exp_lowered = xr.DataArray(["some word ǆ ß ᾓ σς ﬃ⁵å ç ⅰ"]).astype(np.str_)
    exp_swapped = xr.DataArray(["soME WoRD ǆ SS ᾛ σς FFI⁵å ç ⅰ"]).astype(np.str_)
    exp_titled = xr.DataArray(["Some Word ǅ Ss ᾛ Σς Ffi⁵Å Ç Ⅰ"]).astype(np.str_)
    exp_uppered = xr.DataArray(["SOME WORD Ǆ SS ἫΙ ΣΣ FFI⁵Å Ç Ⅰ"]).astype(np.str_)
    exp_casefolded = xr.DataArray(["some word ǆ ss ἣι σσ ffi⁵å ç ⅰ"]).astype(np.str_)

    exp_norm_nfc = xr.DataArray(["SOme wOrd Ǆ ß ᾛ ΣΣ ﬃ⁵Å Ç Ⅰ"]).astype(np.str_)
    exp_norm_nfkc = xr.DataArray(["SOme wOrd DŽ ß ᾛ ΣΣ ffi5Å Ç I"]).astype(np.str_)
    exp_norm_nfd = xr.DataArray(["SOme wOrd Ǆ ß ᾛ ΣΣ ﬃ⁵Å Ç Ⅰ"]).astype(np.str_)
    exp_norm_nfkd = xr.DataArray(["SOme wOrd DŽ ß ᾛ ΣΣ ffi5Å Ç I"]).astype(np.str_)

    res_capitalized = value.str.capitalize()
    res_casefolded = value.str.casefold()
    res_lowered = value.str.lower()
    res_swapped = value.str.swapcase()
    res_titled = value.str.title()
    res_uppered = value.str.upper()

    res_norm_nfc = value.str.normalize("NFC")
    res_norm_nfd = value.str.normalize("NFD")
    res_norm_nfkc = value.str.normalize("NFKC")
    res_norm_nfkd = value.str.normalize("NFKD")

    assert res_capitalized.dtype == exp_capitalized.dtype
    assert res_casefolded.dtype == exp_casefolded.dtype
    assert res_lowered.dtype == exp_lowered.dtype
    assert res_swapped.dtype == exp_swapped.dtype
    assert res_titled.dtype == exp_titled.dtype
    assert res_uppered.dtype == exp_uppered.dtype

    assert res_norm_nfc.dtype == exp_norm_nfc.dtype
    assert res_norm_nfd.dtype == exp_norm_nfd.dtype
    assert res_norm_nfkc.dtype == exp_norm_nfkc.dtype
    assert res_norm_nfkd.dtype == exp_norm_nfkd.dtype

    assert_equal(res_capitalized, exp_capitalized)
    assert_equal(res_casefolded, exp_casefolded)
    assert_equal(res_lowered, exp_lowered)
    assert_equal(res_swapped, exp_swapped)
    assert_equal(res_titled, exp_titled)
    assert_equal(res_uppered, exp_uppered)

    assert_equal(res_norm_nfc, exp_norm_nfc)
    assert_equal(res_norm_nfd, exp_norm_nfd)
    assert_equal(res_norm_nfkc, exp_norm_nfkc)
    assert_equal(res_norm_nfkd, exp_norm_nfkd)


def test_replace(dtype) -> None:
    values = xr.DataArray(["fooBAD__barBAD"], dims=["x"]).astype(dtype)
    result = values.str.replace("BAD[_]*", "")
    expected = xr.DataArray(["foobar"], dims=["x"]).astype(dtype)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    result = values.str.replace("BAD[_]*", "", n=1)
    expected = xr.DataArray(["foobarBAD"], dims=["x"]).astype(dtype)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    pat = xr.DataArray(["BAD[_]*", "AD[_]*"], dims=["y"]).astype(dtype)
    result = values.str.replace(pat, "")
    expected = xr.DataArray([["foobar", "fooBbarB"]], dims=["x", "y"]).astype(dtype)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    repl = xr.DataArray(["", "spam"], dims=["y"]).astype(dtype)
    result = values.str.replace(pat, repl, n=1)
    expected = xr.DataArray([["foobarBAD", "fooBspambarBAD"]], dims=["x", "y"]).astype(
        dtype
    )
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    values = xr.DataArray(
        ["A", "B", "C", "Aaba", "Baca", "", "CABA", "dog", "cat"]
    ).astype(dtype)
    expected = xr.DataArray(
        ["YYY", "B", "C", "YYYaba", "Baca", "", "CYYYBYYY", "dog", "cat"]
    ).astype(dtype)
    result = values.str.replace("A", "YYY")
    assert result.dtype == expected.dtype
    assert_equal(result, expected)
    result = values.str.replace("A", "YYY", regex=False)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    result = values.str.replace("A", "YYY", case=False)
    expected = xr.DataArray(
        ["YYY", "B", "C", "YYYYYYbYYY", "BYYYcYYY", "", "CYYYBYYY", "dog", "cYYYt"]
    ).astype(dtype)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    result = values.str.replace("^.a|dog", "XX-XX ", case=False)
    expected = xr.DataArray(
        ["A", "B", "C", "XX-XX ba", "XX-XX ca", "", "XX-XX BA", "XX-XX ", "XX-XX t"]
    ).astype(dtype)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)


def test_replace_callable() -> None:
    values = xr.DataArray(["fooBAD__barBAD"])

    # test with callable
    repl = lambda m: m.group(0).swapcase()
    result = values.str.replace("[a-z][A-Z]{2}", repl, n=2)
    exp = xr.DataArray(["foObaD__baRbaD"])
    assert result.dtype == exp.dtype
    assert_equal(result, exp)

    # test regex named groups
    values = xr.DataArray(["Foo Bar Baz"])
    pat = r"(?P<first>\w+) (?P<middle>\w+) (?P<last>\w+)"
    repl = lambda m: m.group("middle").swapcase()
    result = values.str.replace(pat, repl)
    exp = xr.DataArray(["bAR"])
    assert result.dtype == exp.dtype
    assert_equal(result, exp)

    # test broadcast
    values = xr.DataArray(["Foo Bar Baz"], dims=["x"])
    pat = r"(?P<first>\w+) (?P<middle>\w+) (?P<last>\w+)"
    repl2 = xr.DataArray(
        [
            lambda m: m.group("first").swapcase(),
            lambda m: m.group("middle").swapcase(),
            lambda m: m.group("last").swapcase(),
        ],
        dims=["Y"],
    )
    result = values.str.replace(pat, repl2)
    exp = xr.DataArray([["fOO", "bAR", "bAZ"]], dims=["x", "Y"])
    assert result.dtype == exp.dtype
    assert_equal(result, exp)


def test_replace_unicode() -> None:
    # flags + unicode
    values = xr.DataArray([b"abcd,\xc3\xa0".decode("utf-8")])
    expected = xr.DataArray([b"abcd, \xc3\xa0".decode("utf-8")])
    pat = re.compile(r"(?<=\w),(?=\w)", flags=re.UNICODE)
    result = values.str.replace(pat, ", ")
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    # broadcast version
    values = xr.DataArray([b"abcd,\xc3\xa0".decode("utf-8")], dims=["X"])
    expected = xr.DataArray(
        [[b"abcd, \xc3\xa0".decode("utf-8"), b"BAcd,\xc3\xa0".decode("utf-8")]],
        dims=["X", "Y"],
    )
    pat2 = xr.DataArray(
        [re.compile(r"(?<=\w),(?=\w)", flags=re.UNICODE), r"ab"], dims=["Y"]
    )
    repl = xr.DataArray([", ", "BA"], dims=["Y"])
    result = values.str.replace(pat2, repl)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)


def test_replace_compiled_regex(dtype) -> None:
    values = xr.DataArray(["fooBAD__barBAD"], dims=["x"]).astype(dtype)

    # test with compiled regex
    pat = re.compile(dtype("BAD[_]*"))
    result = values.str.replace(pat, "")
    expected = xr.DataArray(["foobar"], dims=["x"]).astype(dtype)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    result = values.str.replace(pat, "", n=1)
    expected = xr.DataArray(["foobarBAD"], dims=["x"]).astype(dtype)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    # broadcast
    pat2 = xr.DataArray(
        [re.compile(dtype("BAD[_]*")), re.compile(dtype("AD[_]*"))], dims=["y"]
    )
    result = values.str.replace(pat2, "")
    expected = xr.DataArray([["foobar", "fooBbarB"]], dims=["x", "y"]).astype(dtype)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    repl = xr.DataArray(["", "spam"], dims=["y"]).astype(dtype)
    result = values.str.replace(pat2, repl, n=1)
    expected = xr.DataArray([["foobarBAD", "fooBspambarBAD"]], dims=["x", "y"]).astype(
        dtype
    )
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    # case and flags provided to str.replace will have no effect
    # and will produce warnings
    values = xr.DataArray(["fooBAD__barBAD__bad"]).astype(dtype)
    pat3 = re.compile(dtype("BAD[_]*"))

    with pytest.raises(
        ValueError, match=r"Flags cannot be set when pat is a compiled regex."
    ):
        result = values.str.replace(pat3, "", flags=re.IGNORECASE)

    with pytest.raises(
        ValueError, match=r"Case cannot be set when pat is a compiled regex."
    ):
        result = values.str.replace(pat3, "", case=False)

    with pytest.raises(
        ValueError, match=r"Case cannot be set when pat is a compiled regex."
    ):
        result = values.str.replace(pat3, "", case=True)

    # test with callable
    values = xr.DataArray(["fooBAD__barBAD"]).astype(dtype)
    repl2 = lambda m: m.group(0).swapcase()
    pat4 = re.compile(dtype("[a-z][A-Z]{2}"))
    result = values.str.replace(pat4, repl2, n=2)
    expected = xr.DataArray(["foObaD__baRbaD"]).astype(dtype)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)


def test_replace_literal(dtype) -> None:
    # GH16808 literal replace (regex=False vs regex=True)
    values = xr.DataArray(["f.o", "foo"], dims=["X"]).astype(dtype)
    expected = xr.DataArray(["bao", "bao"], dims=["X"]).astype(dtype)
    result = values.str.replace("f.", "ba")
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    expected = xr.DataArray(["bao", "foo"], dims=["X"]).astype(dtype)
    result = values.str.replace("f.", "ba", regex=False)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    # Broadcast
    pat = xr.DataArray(["f.", ".o"], dims=["yy"]).astype(dtype)
    expected = xr.DataArray([["bao", "fba"], ["bao", "bao"]], dims=["X", "yy"]).astype(
        dtype
    )
    result = values.str.replace(pat, "ba")
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    expected = xr.DataArray([["bao", "fba"], ["foo", "foo"]], dims=["X", "yy"]).astype(
        dtype
    )
    result = values.str.replace(pat, "ba", regex=False)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    # Cannot do a literal replace if given a callable repl or compiled
    # pattern
    callable_repl = lambda m: m.group(0).swapcase()
    compiled_pat = re.compile("[a-z][A-Z]{2}")

    msg = "Cannot use a callable replacement when regex=False"
    with pytest.raises(ValueError, match=msg):
        values.str.replace("abc", callable_repl, regex=False)

    msg = "Cannot use a compiled regex as replacement pattern with regex=False"
    with pytest.raises(ValueError, match=msg):
        values.str.replace(compiled_pat, "", regex=False)


def test_extract_extractall_findall_empty_raises(dtype) -> None:
    pat_str = dtype(r".*")
    pat_re = re.compile(pat_str)

    value = xr.DataArray([["a"]], dims=["X", "Y"]).astype(dtype)

    with pytest.raises(ValueError, match=r"No capture groups found in pattern."):
        value.str.extract(pat=pat_str, dim="ZZ")

    with pytest.raises(ValueError, match=r"No capture groups found in pattern."):
        value.str.extract(pat=pat_re, dim="ZZ")

    with pytest.raises(ValueError, match=r"No capture groups found in pattern."):
        value.str.extractall(pat=pat_str, group_dim="XX", match_dim="YY")

    with pytest.raises(ValueError, match=r"No capture groups found in pattern."):
        value.str.extractall(pat=pat_re, group_dim="XX", match_dim="YY")

    with pytest.raises(ValueError, match=r"No capture groups found in pattern."):
        value.str.findall(pat=pat_str)

    with pytest.raises(ValueError, match=r"No capture groups found in pattern."):
        value.str.findall(pat=pat_re)


def test_extract_multi_None_raises(dtype) -> None:
    pat_str = r"(\w+)_(\d+)"
    pat_re = re.compile(pat_str)

    value = xr.DataArray([["a_b"]], dims=["X", "Y"]).astype(dtype)

    with pytest.raises(
        ValueError,
        match=r"Dimension must be specified if more than one capture group is given.",
    ):
        value.str.extract(pat=pat_str, dim=None)

    with pytest.raises(
        ValueError,
        match=r"Dimension must be specified if more than one capture group is given.",
    ):
        value.str.extract(pat=pat_re, dim=None)


def test_extract_extractall_findall_case_re_raises(dtype) -> None:
    pat_str = r".*"
    pat_re = re.compile(pat_str)

    value = xr.DataArray([["a"]], dims=["X", "Y"]).astype(dtype)

    with pytest.raises(
        ValueError, match=r"Case cannot be set when pat is a compiled regex."
    ):
        value.str.extract(pat=pat_re, case=True, dim="ZZ")

    with pytest.raises(
        ValueError, match=r"Case cannot be set when pat is a compiled regex."
    ):
        value.str.extract(pat=pat_re, case=False, dim="ZZ")

    with pytest.raises(
        ValueError, match=r"Case cannot be set when pat is a compiled regex."
    ):
        value.str.extractall(pat=pat_re, case=True, group_dim="XX", match_dim="YY")

    with pytest.raises(
        ValueError, match=r"Case cannot be set when pat is a compiled regex."
    ):
        value.str.extractall(pat=pat_re, case=False, group_dim="XX", match_dim="YY")

    with pytest.raises(
        ValueError, match=r"Case cannot be set when pat is a compiled regex."
    ):
        value.str.findall(pat=pat_re, case=True)

    with pytest.raises(
        ValueError, match=r"Case cannot be set when pat is a compiled regex."
    ):
        value.str.findall(pat=pat_re, case=False)


def test_extract_extractall_name_collision_raises(dtype) -> None:
    pat_str = r"(\w+)"
    pat_re = re.compile(pat_str)

    value = xr.DataArray([["a"]], dims=["X", "Y"]).astype(dtype)

    with pytest.raises(KeyError, match=r"Dimension 'X' already present in DataArray."):
        value.str.extract(pat=pat_str, dim="X")

    with pytest.raises(KeyError, match=r"Dimension 'X' already present in DataArray."):
        value.str.extract(pat=pat_re, dim="X")

    with pytest.raises(
        KeyError, match=r"Group dimension 'X' already present in DataArray."
    ):
        value.str.extractall(pat=pat_str, group_dim="X", match_dim="ZZ")

    with pytest.raises(
        KeyError, match=r"Group dimension 'X' already present in DataArray."
    ):
        value.str.extractall(pat=pat_re, group_dim="X", match_dim="YY")

    with pytest.raises(
        KeyError, match=r"Match dimension 'Y' already present in DataArray."
    ):
        value.str.extractall(pat=pat_str, group_dim="XX", match_dim="Y")

    with pytest.raises(
        KeyError, match=r"Match dimension 'Y' already present in DataArray."
    ):
        value.str.extractall(pat=pat_re, group_dim="XX", match_dim="Y")

    with pytest.raises(
        KeyError, match=r"Group dimension 'ZZ' is the same as match dimension 'ZZ'."
    ):
        value.str.extractall(pat=pat_str, group_dim="ZZ", match_dim="ZZ")

    with pytest.raises(
        KeyError, match=r"Group dimension 'ZZ' is the same as match dimension 'ZZ'."
    ):
        value.str.extractall(pat=pat_re, group_dim="ZZ", match_dim="ZZ")


def test_extract_single_case(dtype) -> None:
    pat_str = r"(\w+)_Xy_\d*"
    pat_re: str | bytes = (
        pat_str if dtype == np.str_ else bytes(pat_str, encoding="UTF-8")
    )
    pat_compiled = re.compile(pat_re)

    value = xr.DataArray(
        [
            ["a_Xy_0", "ab_xY_10-bab_Xy_110-baab_Xy_1100", "abc_Xy_01-cbc_Xy_2210"],
            [
                "abcd_Xy_-dcd_Xy_33210-dccd_Xy_332210",
                "",
                "abcdef_Xy_101-fef_Xy_5543210",
            ],
        ],
        dims=["X", "Y"],
    ).astype(dtype)

    targ_none = xr.DataArray(
        [["a", "bab", "abc"], ["abcd", "", "abcdef"]], dims=["X", "Y"]
    ).astype(dtype)
    targ_dim = xr.DataArray(
        [[["a"], ["bab"], ["abc"]], [["abcd"], [""], ["abcdef"]]], dims=["X", "Y", "XX"]
    ).astype(dtype)

    res_str_none = value.str.extract(pat=pat_str, dim=None)
    res_str_dim = value.str.extract(pat=pat_str, dim="XX")
    res_str_none_case = value.str.extract(pat=pat_str, dim=None, case=True)
    res_str_dim_case = value.str.extract(pat=pat_str, dim="XX", case=True)
    res_re_none = value.str.extract(pat=pat_compiled, dim=None)
    res_re_dim = value.str.extract(pat=pat_compiled, dim="XX")

    assert res_str_none.dtype == targ_none.dtype
    assert res_str_dim.dtype == targ_dim.dtype
    assert res_str_none_case.dtype == targ_none.dtype
    assert res_str_dim_case.dtype == targ_dim.dtype
    assert res_re_none.dtype == targ_none.dtype
    assert res_re_dim.dtype == targ_dim.dtype

    assert_equal(res_str_none, targ_none)
    assert_equal(res_str_dim, targ_dim)
    assert_equal(res_str_none_case, targ_none)
    assert_equal(res_str_dim_case, targ_dim)
    assert_equal(res_re_none, targ_none)
    assert_equal(res_re_dim, targ_dim)


def test_extract_single_nocase(dtype) -> None:
    pat_str = r"(\w+)?_Xy_\d*"
    pat_re: str | bytes = (
        pat_str if dtype == np.str_ else bytes(pat_str, encoding="UTF-8")
    )
    pat_compiled = re.compile(pat_re, flags=re.IGNORECASE)

    value = xr.DataArray(
        [
            ["a_Xy_0", "ab_xY_10-bab_Xy_110-baab_Xy_1100", "abc_Xy_01-cbc_Xy_2210"],
            [
                "abcd_Xy_-dcd_Xy_33210-dccd_Xy_332210",
                "_Xy_1",
                "abcdef_Xy_101-fef_Xy_5543210",
            ],
        ],
        dims=["X", "Y"],
    ).astype(dtype)

    targ_none = xr.DataArray(
        [["a", "ab", "abc"], ["abcd", "", "abcdef"]], dims=["X", "Y"]
    ).astype(dtype)
    targ_dim = xr.DataArray(
        [[["a"], ["ab"], ["abc"]], [["abcd"], [""], ["abcdef"]]], dims=["X", "Y", "XX"]
    ).astype(dtype)

    res_str_none = value.str.extract(pat=pat_str, dim=None, case=False)
    res_str_dim = value.str.extract(pat=pat_str, dim="XX", case=False)
    res_re_none = value.str.extract(pat=pat_compiled, dim=None)
    res_re_dim = value.str.extract(pat=pat_compiled, dim="XX")

    assert res_re_dim.dtype == targ_none.dtype
    assert res_str_dim.dtype == targ_dim.dtype
    assert res_re_none.dtype == targ_none.dtype
    assert res_re_dim.dtype == targ_dim.dtype

    assert_equal(res_str_none, targ_none)
    assert_equal(res_str_dim, targ_dim)
    assert_equal(res_re_none, targ_none)
    assert_equal(res_re_dim, targ_dim)


def test_extract_multi_case(dtype) -> None:
    pat_str = r"(\w+)_Xy_(\d*)"
    pat_re: str | bytes = (
        pat_str if dtype == np.str_ else bytes(pat_str, encoding="UTF-8")
    )
    pat_compiled = re.compile(pat_re)

    value = xr.DataArray(
        [
            ["a_Xy_0", "ab_xY_10-bab_Xy_110-baab_Xy_1100", "abc_Xy_01-cbc_Xy_2210"],
            [
                "abcd_Xy_-dcd_Xy_33210-dccd_Xy_332210",
                "",
                "abcdef_Xy_101-fef_Xy_5543210",
            ],
        ],
        dims=["X", "Y"],
    ).astype(dtype)

    expected = xr.DataArray(
        [
            [["a", "0"], ["bab", "110"], ["abc", "01"]],
            [["abcd", ""], ["", ""], ["abcdef", "101"]],
        ],
        dims=["X", "Y", "XX"],
    ).astype(dtype)

    res_str = value.str.extract(pat=pat_str, dim="XX")
    res_re = value.str.extract(pat=pat_compiled, dim="XX")
    res_str_case = value.str.extract(pat=pat_str, dim="XX", case=True)

    assert res_str.dtype == expected.dtype
    assert res_re.dtype == expected.dtype
    assert res_str_case.dtype == expected.dtype

    assert_equal(res_str, expected)
    assert_equal(res_re, expected)
    assert_equal(res_str_case, expected)


def test_extract_multi_nocase(dtype) -> None:
    pat_str = r"(\w+)_Xy_(\d*)"
    pat_re: str | bytes = (
        pat_str if dtype == np.str_ else bytes(pat_str, encoding="UTF-8")
    )
    pat_compiled = re.compile(pat_re, flags=re.IGNORECASE)

    value = xr.DataArray(
        [
            ["a_Xy_0", "ab_xY_10-bab_Xy_110-baab_Xy_1100", "abc_Xy_01-cbc_Xy_2210"],
            [
                "abcd_Xy_-dcd_Xy_33210-dccd_Xy_332210",
                "",
                "abcdef_Xy_101-fef_Xy_5543210",
            ],
        ],
        dims=["X", "Y"],
    ).astype(dtype)

    expected = xr.DataArray(
        [
            [["a", "0"], ["ab", "10"], ["abc", "01"]],
            [["abcd", ""], ["", ""], ["abcdef", "101"]],
        ],
        dims=["X", "Y", "XX"],
    ).astype(dtype)

    res_str = value.str.extract(pat=pat_str, dim="XX", case=False)
    res_re = value.str.extract(pat=pat_compiled, dim="XX")

    assert res_str.dtype == expected.dtype
    assert res_re.dtype == expected.dtype

    assert_equal(res_str, expected)
    assert_equal(res_re, expected)


def test_extract_broadcast(dtype) -> None:
    value = xr.DataArray(
        ["a_Xy_0", "ab_xY_10", "abc_Xy_01"],
        dims=["X"],
    ).astype(dtype)

    pat_str = xr.DataArray(
        [r"(\w+)_Xy_(\d*)", r"(\w+)_xY_(\d*)"],
        dims=["Y"],
    ).astype(dtype)
    pat_compiled = value.str._re_compile(pat=pat_str)

    expected_list = [
        [["a", "0"], ["", ""]],
        [["", ""], ["ab", "10"]],
        [["abc", "01"], ["", ""]],
    ]
    expected = xr.DataArray(expected_list, dims=["X", "Y", "Zz"]).astype(dtype)

    res_str = value.str.extract(pat=pat_str, dim="Zz")
    res_re = value.str.extract(pat=pat_compiled, dim="Zz")

    assert res_str.dtype == expected.dtype
    assert res_re.dtype == expected.dtype

    assert_equal(res_str, expected)
    assert_equal(res_re, expected)


def test_extractall_single_single_case(dtype) -> None:
    pat_str = r"(\w+)_Xy_\d*"
    pat_re: str | bytes = (
        pat_str if dtype == np.str_ else bytes(pat_str, encoding="UTF-8")
    )
    pat_compiled = re.compile(pat_re)

    value = xr.DataArray(
        [["a_Xy_0", "ab_xY_10", "abc_Xy_01"], ["abcd_Xy_", "", "abcdef_Xy_101"]],
        dims=["X", "Y"],
    ).astype(dtype)

    expected = xr.DataArray(
        [[[["a"]], [[""]], [["abc"]]], [[["abcd"]], [[""]], [["abcdef"]]]],
        dims=["X", "Y", "XX", "YY"],
    ).astype(dtype)

    res_str = value.str.extractall(pat=pat_str, group_dim="XX", match_dim="YY")
    res_re = value.str.extractall(pat=pat_compiled, group_dim="XX", match_dim="YY")
    res_str_case = value.str.extractall(
        pat=pat_str, group_dim="XX", match_dim="YY", case=True
    )

    assert res_str.dtype == expected.dtype
    assert res_re.dtype == expected.dtype
    assert res_str_case.dtype == expected.dtype

    assert_equal(res_str, expected)
    assert_equal(res_re, expected)
    assert_equal(res_str_case, expected)


def test_extractall_single_single_nocase(dtype) -> None:
    pat_str = r"(\w+)_Xy_\d*"
    pat_re: str | bytes = (
        pat_str if dtype == np.str_ else bytes(pat_str, encoding="UTF-8")
    )
    pat_compiled = re.compile(pat_re, flags=re.IGNORECASE)

    value = xr.DataArray(
        [["a_Xy_0", "ab_xY_10", "abc_Xy_01"], ["abcd_Xy_", "", "abcdef_Xy_101"]],
        dims=["X", "Y"],
    ).astype(dtype)

    expected = xr.DataArray(
        [[[["a"]], [["ab"]], [["abc"]]], [[["abcd"]], [[""]], [["abcdef"]]]],
        dims=["X", "Y", "XX", "YY"],
    ).astype(dtype)

    res_str = value.str.extractall(
        pat=pat_str, group_dim="XX", match_dim="YY", case=False
    )
    res_re = value.str.extractall(pat=pat_compiled, group_dim="XX", match_dim="YY")

    assert res_str.dtype == expected.dtype
    assert res_re.dtype == expected.dtype

    assert_equal(res_str, expected)
    assert_equal(res_re, expected)


def test_extractall_single_multi_case(dtype) -> None:
    pat_str = r"(\w+)_Xy_\d*"
    pat_re: str | bytes = (
        pat_str if dtype == np.str_ else bytes(pat_str, encoding="UTF-8")
    )
    pat_compiled = re.compile(pat_re)

    value = xr.DataArray(
        [
            ["a_Xy_0", "ab_xY_10-bab_Xy_110-baab_Xy_1100", "abc_Xy_01-cbc_Xy_2210"],
            [
                "abcd_Xy_-dcd_Xy_33210-dccd_Xy_332210",
                "",
                "abcdef_Xy_101-fef_Xy_5543210",
            ],
        ],
        dims=["X", "Y"],
    ).astype(dtype)

    expected = xr.DataArray(
        [
            [[["a"], [""], [""]], [["bab"], ["baab"], [""]], [["abc"], ["cbc"], [""]]],
            [
                [["abcd"], ["dcd"], ["dccd"]],
                [[""], [""], [""]],
                [["abcdef"], ["fef"], [""]],
            ],
        ],
        dims=["X", "Y", "XX", "YY"],
    ).astype(dtype)

    res_str = value.str.extractall(pat=pat_str, group_dim="XX", match_dim="YY")
    res_re = value.str.extractall(pat=pat_compiled, group_dim="XX", match_dim="YY")
    res_str_case = value.str.extractall(
        pat=pat_str, group_dim="XX", match_dim="YY", case=True
    )

    assert res_str.dtype == expected.dtype
    assert res_re.dtype == expected.dtype
    assert res_str_case.dtype == expected.dtype

    assert_equal(res_str, expected)
    assert_equal(res_re, expected)
    assert_equal(res_str_case, expected)


def test_extractall_single_multi_nocase(dtype) -> None:
    pat_str = r"(\w+)_Xy_\d*"
    pat_re: str | bytes = (
        pat_str if dtype == np.str_ else bytes(pat_str, encoding="UTF-8")
    )
    pat_compiled = re.compile(pat_re, flags=re.IGNORECASE)

    value = xr.DataArray(
        [
            ["a_Xy_0", "ab_xY_10-bab_Xy_110-baab_Xy_1100", "abc_Xy_01-cbc_Xy_2210"],
            [
                "abcd_Xy_-dcd_Xy_33210-dccd_Xy_332210",
                "",
                "abcdef_Xy_101-fef_Xy_5543210",
            ],
        ],
        dims=["X", "Y"],
    ).astype(dtype)

    expected = xr.DataArray(
        [
            [
                [["a"], [""], [""]],
                [["ab"], ["bab"], ["baab"]],
                [["abc"], ["cbc"], [""]],
            ],
            [
                [["abcd"], ["dcd"], ["dccd"]],
                [[""], [""], [""]],
                [["abcdef"], ["fef"], [""]],
            ],
        ],
        dims=["X", "Y", "XX", "YY"],
    ).astype(dtype)

    res_str = value.str.extractall(
        pat=pat_str, group_dim="XX", match_dim="YY", case=False
    )
    res_re = value.str.extractall(pat=pat_compiled, group_dim="XX", match_dim="YY")

    assert res_str.dtype == expected.dtype
    assert res_re.dtype == expected.dtype

    assert_equal(res_str, expected)
    assert_equal(res_re, expected)


def test_extractall_multi_single_case(dtype) -> None:
    pat_str = r"(\w+)_Xy_(\d*)"
    pat_re: str | bytes = (
        pat_str if dtype == np.str_ else bytes(pat_str, encoding="UTF-8")
    )
    pat_compiled = re.compile(pat_re)

    value = xr.DataArray(
        [["a_Xy_0", "ab_xY_10", "abc_Xy_01"], ["abcd_Xy_", "", "abcdef_Xy_101"]],
        dims=["X", "Y"],
    ).astype(dtype)

    expected = xr.DataArray(
        [
            [[["a", "0"]], [["", ""]], [["abc", "01"]]],
            [[["abcd", ""]], [["", ""]], [["abcdef", "101"]]],
        ],
        dims=["X", "Y", "XX", "YY"],
    ).astype(dtype)

    res_str = value.str.extractall(pat=pat_str, group_dim="XX", match_dim="YY")
    res_re = value.str.extractall(pat=pat_compiled, group_dim="XX", match_dim="YY")
    res_str_case = value.str.extractall(
        pat=pat_str, group_dim="XX", match_dim="YY", case=True
    )

    assert res_str.dtype == expected.dtype
    assert res_re.dtype == expected.dtype
    assert res_str_case.dtype == expected.dtype

    assert_equal(res_str, expected)
    assert_equal(res_re, expected)
    assert_equal(res_str_case, expected)


def test_extractall_multi_single_nocase(dtype) -> None:
    pat_str = r"(\w+)_Xy_(\d*)"
    pat_re: str | bytes = (
        pat_str if dtype == np.str_ else bytes(pat_str, encoding="UTF-8")
    )
    pat_compiled = re.compile(pat_re, flags=re.IGNORECASE)

    value = xr.DataArray(
        [["a_Xy_0", "ab_xY_10", "abc_Xy_01"], ["abcd_Xy_", "", "abcdef_Xy_101"]],
        dims=["X", "Y"],
    ).astype(dtype)

    expected = xr.DataArray(
        [
            [[["a", "0"]], [["ab", "10"]], [["abc", "01"]]],
            [[["abcd", ""]], [["", ""]], [["abcdef", "101"]]],
        ],
        dims=["X", "Y", "XX", "YY"],
    ).astype(dtype)

    res_str = value.str.extractall(
        pat=pat_str, group_dim="XX", match_dim="YY", case=False
    )
    res_re = value.str.extractall(pat=pat_compiled, group_dim="XX", match_dim="YY")

    assert res_str.dtype == expected.dtype
    assert res_re.dtype == expected.dtype

    assert_equal(res_str, expected)
    assert_equal(res_re, expected)


def test_extractall_multi_multi_case(dtype) -> None:
    pat_str = r"(\w+)_Xy_(\d*)"
    pat_re: str | bytes = (
        pat_str if dtype == np.str_ else bytes(pat_str, encoding="UTF-8")
    )
    pat_compiled = re.compile(pat_re)

    value = xr.DataArray(
        [
            ["a_Xy_0", "ab_xY_10-bab_Xy_110-baab_Xy_1100", "abc_Xy_01-cbc_Xy_2210"],
            [
                "abcd_Xy_-dcd_Xy_33210-dccd_Xy_332210",
                "",
                "abcdef_Xy_101-fef_Xy_5543210",
            ],
        ],
        dims=["X", "Y"],
    ).astype(dtype)

    expected = xr.DataArray(
        [
            [
                [["a", "0"], ["", ""], ["", ""]],
                [["bab", "110"], ["baab", "1100"], ["", ""]],
                [["abc", "01"], ["cbc", "2210"], ["", ""]],
            ],
            [
                [["abcd", ""], ["dcd", "33210"], ["dccd", "332210"]],
                [["", ""], ["", ""], ["", ""]],
                [["abcdef", "101"], ["fef", "5543210"], ["", ""]],
            ],
        ],
        dims=["X", "Y", "XX", "YY"],
    ).astype(dtype)

    res_str = value.str.extractall(pat=pat_str, group_dim="XX", match_dim="YY")
    res_re = value.str.extractall(pat=pat_compiled, group_dim="XX", match_dim="YY")
    res_str_case = value.str.extractall(
        pat=pat_str, group_dim="XX", match_dim="YY", case=True
    )

    assert res_str.dtype == expected.dtype
    assert res_re.dtype == expected.dtype
    assert res_str_case.dtype == expected.dtype

    assert_equal(res_str, expected)
    assert_equal(res_re, expected)
    assert_equal(res_str_case, expected)


def test_extractall_multi_multi_nocase(dtype) -> None:
    pat_str = r"(\w+)_Xy_(\d*)"
    pat_re: str | bytes = (
        pat_str if dtype == np.str_ else bytes(pat_str, encoding="UTF-8")
    )
    pat_compiled = re.compile(pat_re, flags=re.IGNORECASE)

    value = xr.DataArray(
        [
            ["a_Xy_0", "ab_xY_10-bab_Xy_110-baab_Xy_1100", "abc_Xy_01-cbc_Xy_2210"],
            [
                "abcd_Xy_-dcd_Xy_33210-dccd_Xy_332210",
                "",
                "abcdef_Xy_101-fef_Xy_5543210",
            ],
        ],
        dims=["X", "Y"],
    ).astype(dtype)

    expected = xr.DataArray(
        [
            [
                [["a", "0"], ["", ""], ["", ""]],
                [["ab", "10"], ["bab", "110"], ["baab", "1100"]],
                [["abc", "01"], ["cbc", "2210"], ["", ""]],
            ],
            [
                [["abcd", ""], ["dcd", "33210"], ["dccd", "332210"]],
                [["", ""], ["", ""], ["", ""]],
                [["abcdef", "101"], ["fef", "5543210"], ["", ""]],
            ],
        ],
        dims=["X", "Y", "XX", "YY"],
    ).astype(dtype)

    res_str = value.str.extractall(
        pat=pat_str, group_dim="XX", match_dim="YY", case=False
    )
    res_re = value.str.extractall(pat=pat_compiled, group_dim="XX", match_dim="YY")

    assert res_str.dtype == expected.dtype
    assert res_re.dtype == expected.dtype

    assert_equal(res_str, expected)
    assert_equal(res_re, expected)


def test_extractall_broadcast(dtype) -> None:
    value = xr.DataArray(
        ["a_Xy_0", "ab_xY_10", "abc_Xy_01"],
        dims=["X"],
    ).astype(dtype)

    pat_str = xr.DataArray(
        [r"(\w+)_Xy_(\d*)", r"(\w+)_xY_(\d*)"],
        dims=["Y"],
    ).astype(dtype)
    pat_re = value.str._re_compile(pat=pat_str)

    expected_list = [
        [[["a", "0"]], [["", ""]]],
        [[["", ""]], [["ab", "10"]]],
        [[["abc", "01"]], [["", ""]]],
    ]
    expected = xr.DataArray(expected_list, dims=["X", "Y", "ZX", "ZY"]).astype(dtype)

    res_str = value.str.extractall(pat=pat_str, group_dim="ZX", match_dim="ZY")
    res_re = value.str.extractall(pat=pat_re, group_dim="ZX", match_dim="ZY")

    assert res_str.dtype == expected.dtype
    assert res_re.dtype == expected.dtype

    assert_equal(res_str, expected)
    assert_equal(res_re, expected)


def test_findall_single_single_case(dtype) -> None:
    pat_str = r"(\w+)_Xy_\d*"
    pat_re = re.compile(dtype(pat_str))

    value = xr.DataArray(
        [["a_Xy_0", "ab_xY_10", "abc_Xy_01"], ["abcd_Xy_", "", "abcdef_Xy_101"]],
        dims=["X", "Y"],
    ).astype(dtype)

    expected_list: list[list[list]] = [[["a"], [], ["abc"]], [["abcd"], [], ["abcdef"]]]
    expected_dtype = [[[dtype(x) for x in y] for y in z] for z in expected_list]
    expected_np = np.array(expected_dtype, dtype=np.object_)
    expected = xr.DataArray(expected_np, dims=["X", "Y"])

    res_str = value.str.findall(pat=pat_str)
    res_re = value.str.findall(pat=pat_re)
    res_str_case = value.str.findall(pat=pat_str, case=True)

    assert res_str.dtype == expected.dtype
    assert res_re.dtype == expected.dtype
    assert res_str_case.dtype == expected.dtype

    assert_equal(res_str, expected)
    assert_equal(res_re, expected)
    assert_equal(res_str_case, expected)


def test_findall_single_single_nocase(dtype) -> None:
    pat_str = r"(\w+)_Xy_\d*"
    pat_re = re.compile(dtype(pat_str), flags=re.IGNORECASE)

    value = xr.DataArray(
        [["a_Xy_0", "ab_xY_10", "abc_Xy_01"], ["abcd_Xy_", "", "abcdef_Xy_101"]],
        dims=["X", "Y"],
    ).astype(dtype)

    expected_list: list[list[list]] = [
        [["a"], ["ab"], ["abc"]],
        [["abcd"], [], ["abcdef"]],
    ]
    expected_dtype = [[[dtype(x) for x in y] for y in z] for z in expected_list]
    expected_np = np.array(expected_dtype, dtype=np.object_)
    expected = xr.DataArray(expected_np, dims=["X", "Y"])

    res_str = value.str.findall(pat=pat_str, case=False)
    res_re = value.str.findall(pat=pat_re)

    assert res_str.dtype == expected.dtype
    assert res_re.dtype == expected.dtype

    assert_equal(res_str, expected)
    assert_equal(res_re, expected)


def test_findall_single_multi_case(dtype) -> None:
    pat_str = r"(\w+)_Xy_\d*"
    pat_re = re.compile(dtype(pat_str))

    value = xr.DataArray(
        [
            ["a_Xy_0", "ab_xY_10-bab_Xy_110-baab_Xy_1100", "abc_Xy_01-cbc_Xy_2210"],
            [
                "abcd_Xy_-dcd_Xy_33210-dccd_Xy_332210",
                "",
                "abcdef_Xy_101-fef_Xy_5543210",
            ],
        ],
        dims=["X", "Y"],
    ).astype(dtype)

    expected_list: list[list[list]] = [
        [["a"], ["bab", "baab"], ["abc", "cbc"]],
        [
            ["abcd", "dcd", "dccd"],
            [],
            ["abcdef", "fef"],
        ],
    ]
    expected_dtype = [[[dtype(x) for x in y] for y in z] for z in expected_list]
    expected_np = np.array(expected_dtype, dtype=np.object_)
    expected = xr.DataArray(expected_np, dims=["X", "Y"])

    res_str = value.str.findall(pat=pat_str)
    res_re = value.str.findall(pat=pat_re)
    res_str_case = value.str.findall(pat=pat_str, case=True)

    assert res_str.dtype == expected.dtype
    assert res_re.dtype == expected.dtype
    assert res_str_case.dtype == expected.dtype

    assert_equal(res_str, expected)
    assert_equal(res_re, expected)
    assert_equal(res_str_case, expected)


def test_findall_single_multi_nocase(dtype) -> None:
    pat_str = r"(\w+)_Xy_\d*"
    pat_re = re.compile(dtype(pat_str), flags=re.IGNORECASE)

    value = xr.DataArray(
        [
            ["a_Xy_0", "ab_xY_10-bab_Xy_110-baab_Xy_1100", "abc_Xy_01-cbc_Xy_2210"],
            [
                "abcd_Xy_-dcd_Xy_33210-dccd_Xy_332210",
                "",
                "abcdef_Xy_101-fef_Xy_5543210",
            ],
        ],
        dims=["X", "Y"],
    ).astype(dtype)

    expected_list: list[list[list]] = [
        [
            ["a"],
            ["ab", "bab", "baab"],
            ["abc", "cbc"],
        ],
        [
            ["abcd", "dcd", "dccd"],
            [],
            ["abcdef", "fef"],
        ],
    ]
    expected_dtype = [[[dtype(x) for x in y] for y in z] for z in expected_list]
    expected_np = np.array(expected_dtype, dtype=np.object_)
    expected = xr.DataArray(expected_np, dims=["X", "Y"])

    res_str = value.str.findall(pat=pat_str, case=False)
    res_re = value.str.findall(pat=pat_re)

    assert res_str.dtype == expected.dtype
    assert res_re.dtype == expected.dtype

    assert_equal(res_str, expected)
    assert_equal(res_re, expected)


def test_findall_multi_single_case(dtype) -> None:
    pat_str = r"(\w+)_Xy_(\d*)"
    pat_re = re.compile(dtype(pat_str))

    value = xr.DataArray(
        [["a_Xy_0", "ab_xY_10", "abc_Xy_01"], ["abcd_Xy_", "", "abcdef_Xy_101"]],
        dims=["X", "Y"],
    ).astype(dtype)

    expected_list: list[list[list[list]]] = [
        [[["a", "0"]], [], [["abc", "01"]]],
        [[["abcd", ""]], [], [["abcdef", "101"]]],
    ]
    expected_dtype = [
        [[tuple(dtype(x) for x in y) for y in z] for z in w] for w in expected_list
    ]
    expected_np = np.array(expected_dtype, dtype=np.object_)
    expected = xr.DataArray(expected_np, dims=["X", "Y"])

    res_str = value.str.findall(pat=pat_str)
    res_re = value.str.findall(pat=pat_re)
    res_str_case = value.str.findall(pat=pat_str, case=True)

    assert res_str.dtype == expected.dtype
    assert res_re.dtype == expected.dtype
    assert res_str_case.dtype == expected.dtype

    assert_equal(res_str, expected)
    assert_equal(res_re, expected)
    assert_equal(res_str_case, expected)


def test_findall_multi_single_nocase(dtype) -> None:
    pat_str = r"(\w+)_Xy_(\d*)"
    pat_re = re.compile(dtype(pat_str), flags=re.IGNORECASE)

    value = xr.DataArray(
        [["a_Xy_0", "ab_xY_10", "abc_Xy_01"], ["abcd_Xy_", "", "abcdef_Xy_101"]],
        dims=["X", "Y"],
    ).astype(dtype)

    expected_list: list[list[list[list]]] = [
        [[["a", "0"]], [["ab", "10"]], [["abc", "01"]]],
        [[["abcd", ""]], [], [["abcdef", "101"]]],
    ]
    expected_dtype = [
        [[tuple(dtype(x) for x in y) for y in z] for z in w] for w in expected_list
    ]
    expected_np = np.array(expected_dtype, dtype=np.object_)
    expected = xr.DataArray(expected_np, dims=["X", "Y"])

    res_str = value.str.findall(pat=pat_str, case=False)
    res_re = value.str.findall(pat=pat_re)

    assert res_str.dtype == expected.dtype
    assert res_re.dtype == expected.dtype

    assert_equal(res_str, expected)
    assert_equal(res_re, expected)


def test_findall_multi_multi_case(dtype) -> None:
    pat_str = r"(\w+)_Xy_(\d*)"
    pat_re = re.compile(dtype(pat_str))

    value = xr.DataArray(
        [
            ["a_Xy_0", "ab_xY_10-bab_Xy_110-baab_Xy_1100", "abc_Xy_01-cbc_Xy_2210"],
            [
                "abcd_Xy_-dcd_Xy_33210-dccd_Xy_332210",
                "",
                "abcdef_Xy_101-fef_Xy_5543210",
            ],
        ],
        dims=["X", "Y"],
    ).astype(dtype)

    expected_list: list[list[list[list]]] = [
        [
            [["a", "0"]],
            [["bab", "110"], ["baab", "1100"]],
            [["abc", "01"], ["cbc", "2210"]],
        ],
        [
            [["abcd", ""], ["dcd", "33210"], ["dccd", "332210"]],
            [],
            [["abcdef", "101"], ["fef", "5543210"]],
        ],
    ]
    expected_dtype = [
        [[tuple(dtype(x) for x in y) for y in z] for z in w] for w in expected_list
    ]
    expected_np = np.array(expected_dtype, dtype=np.object_)
    expected = xr.DataArray(expected_np, dims=["X", "Y"])

    res_str = value.str.findall(pat=pat_str)
    res_re = value.str.findall(pat=pat_re)
    res_str_case = value.str.findall(pat=pat_str, case=True)

    assert res_str.dtype == expected.dtype
    assert res_re.dtype == expected.dtype
    assert res_str_case.dtype == expected.dtype

    assert_equal(res_str, expected)
    assert_equal(res_re, expected)
    assert_equal(res_str_case, expected)


def test_findall_multi_multi_nocase(dtype) -> None:
    pat_str = r"(\w+)_Xy_(\d*)"
    pat_re = re.compile(dtype(pat_str), flags=re.IGNORECASE)

    value = xr.DataArray(
        [
            ["a_Xy_0", "ab_xY_10-bab_Xy_110-baab_Xy_1100", "abc_Xy_01-cbc_Xy_2210"],
            [
                "abcd_Xy_-dcd_Xy_33210-dccd_Xy_332210",
                "",
                "abcdef_Xy_101-fef_Xy_5543210",
            ],
        ],
        dims=["X", "Y"],
    ).astype(dtype)

    expected_list: list[list[list[list]]] = [
        [
            [["a", "0"]],
            [["ab", "10"], ["bab", "110"], ["baab", "1100"]],
            [["abc", "01"], ["cbc", "2210"]],
        ],
        [
            [["abcd", ""], ["dcd", "33210"], ["dccd", "332210"]],
            [],
            [["abcdef", "101"], ["fef", "5543210"]],
        ],
    ]
    expected_dtype = [
        [[tuple(dtype(x) for x in y) for y in z] for z in w] for w in expected_list
    ]
    expected_np = np.array(expected_dtype, dtype=np.object_)
    expected = xr.DataArray(expected_np, dims=["X", "Y"])

    res_str = value.str.findall(pat=pat_str, case=False)
    res_re = value.str.findall(pat=pat_re)

    assert res_str.dtype == expected.dtype
    assert res_re.dtype == expected.dtype

    assert_equal(res_str, expected)
    assert_equal(res_re, expected)


def test_findall_broadcast(dtype) -> None:
    value = xr.DataArray(
        ["a_Xy_0", "ab_xY_10", "abc_Xy_01"],
        dims=["X"],
    ).astype(dtype)

    pat_str = xr.DataArray(
        [r"(\w+)_Xy_\d*", r"\w+_Xy_(\d*)"],
        dims=["Y"],
    ).astype(dtype)
    pat_re = value.str._re_compile(pat=pat_str)

    expected_list: list[list[list]] = [[["a"], ["0"]], [[], []], [["abc"], ["01"]]]
    expected_dtype = [[[dtype(x) for x in y] for y in z] for z in expected_list]
    expected_np = np.array(expected_dtype, dtype=np.object_)
    expected = xr.DataArray(expected_np, dims=["X", "Y"])

    res_str = value.str.findall(pat=pat_str)
    res_re = value.str.findall(pat=pat_re)

    assert res_str.dtype == expected.dtype
    assert res_re.dtype == expected.dtype

    assert_equal(res_str, expected)
    assert_equal(res_re, expected)


def test_repeat(dtype) -> None:
    values = xr.DataArray(["a", "b", "c", "d"]).astype(dtype)

    result = values.str.repeat(3)
    result_mul = values.str * 3

    expected = xr.DataArray(["aaa", "bbb", "ccc", "ddd"]).astype(dtype)

    assert result.dtype == expected.dtype
    assert result_mul.dtype == expected.dtype

    assert_equal(result_mul, expected)
    assert_equal(result, expected)


def test_repeat_broadcast(dtype) -> None:
    values = xr.DataArray(["a", "b", "c", "d"], dims=["X"]).astype(dtype)
    reps = xr.DataArray([3, 4], dims=["Y"])

    result = values.str.repeat(reps)
    result_mul = values.str * reps

    expected = xr.DataArray(
        [["aaa", "aaaa"], ["bbb", "bbbb"], ["ccc", "cccc"], ["ddd", "dddd"]],
        dims=["X", "Y"],
    ).astype(dtype)

    assert result.dtype == expected.dtype
    assert result_mul.dtype == expected.dtype

    assert_equal(result_mul, expected)
    assert_equal(result, expected)


def test_match(dtype) -> None:
    values = xr.DataArray(["fooBAD__barBAD", "foo"]).astype(dtype)

    # New match behavior introduced in 0.13
    pat = values.dtype.type(".*(BAD[_]+).*(BAD)")
    result = values.str.match(pat)
    expected = xr.DataArray([True, False])
    assert result.dtype == expected.dtype
    assert_equal(result, expected)
    result = values.str.match(re.compile(pat))
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    # Case-sensitive
    pat = values.dtype.type(".*BAD[_]+.*BAD")
    result = values.str.match(pat)
    expected = xr.DataArray([True, False])
    assert result.dtype == expected.dtype
    assert_equal(result, expected)
    result = values.str.match(re.compile(pat))
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    # Case-insensitive
    pat = values.dtype.type(".*bAd[_]+.*bad")
    result = values.str.match(pat, case=False)
    expected = xr.DataArray([True, False])
    assert result.dtype == expected.dtype
    assert_equal(result, expected)
    result = values.str.match(re.compile(pat, flags=re.IGNORECASE))
    assert result.dtype == expected.dtype
    assert_equal(result, expected)


def test_empty_str_methods() -> None:
    empty = xr.DataArray(np.empty(shape=(0,), dtype="U"))
    empty_str = empty
    empty_int = xr.DataArray(np.empty(shape=(0,), dtype=int))
    empty_bool = xr.DataArray(np.empty(shape=(0,), dtype=bool))
    empty_bytes = xr.DataArray(np.empty(shape=(0,), dtype="S"))

    # TODO: Determine why U and S dtype sizes don't match and figure
    # out a reliable way to predict what they should be

    assert empty_bool.dtype == empty.str.contains("a").dtype
    assert empty_bool.dtype == empty.str.endswith("a").dtype
    assert empty_bool.dtype == empty.str.match("^a").dtype
    assert empty_bool.dtype == empty.str.startswith("a").dtype
    assert empty_bool.dtype == empty.str.isalnum().dtype
    assert empty_bool.dtype == empty.str.isalpha().dtype
    assert empty_bool.dtype == empty.str.isdecimal().dtype
    assert empty_bool.dtype == empty.str.isdigit().dtype
    assert empty_bool.dtype == empty.str.islower().dtype
    assert empty_bool.dtype == empty.str.isnumeric().dtype
    assert empty_bool.dtype == empty.str.isspace().dtype
    assert empty_bool.dtype == empty.str.istitle().dtype
    assert empty_bool.dtype == empty.str.isupper().dtype
    assert empty_bytes.dtype.kind == empty.str.encode("ascii").dtype.kind
    assert empty_int.dtype.kind == empty.str.count("a").dtype.kind
    assert empty_int.dtype.kind == empty.str.find("a").dtype.kind
    assert empty_int.dtype.kind == empty.str.len().dtype.kind
    assert empty_int.dtype.kind == empty.str.rfind("a").dtype.kind
    assert empty_str.dtype.kind == empty.str.capitalize().dtype.kind
    assert empty_str.dtype.kind == empty.str.center(42).dtype.kind
    assert empty_str.dtype.kind == empty.str.get(0).dtype.kind
    assert empty_str.dtype.kind == empty.str.lower().dtype.kind
    assert empty_str.dtype.kind == empty.str.lstrip().dtype.kind
    assert empty_str.dtype.kind == empty.str.pad(42).dtype.kind
    assert empty_str.dtype.kind == empty.str.repeat(3).dtype.kind
    assert empty_str.dtype.kind == empty.str.rstrip().dtype.kind
    assert empty_str.dtype.kind == empty.str.slice(step=1).dtype.kind
    assert empty_str.dtype.kind == empty.str.slice(stop=1).dtype.kind
    assert empty_str.dtype.kind == empty.str.strip().dtype.kind
    assert empty_str.dtype.kind == empty.str.swapcase().dtype.kind
    assert empty_str.dtype.kind == empty.str.title().dtype.kind
    assert empty_str.dtype.kind == empty.str.upper().dtype.kind
    assert empty_str.dtype.kind == empty.str.wrap(42).dtype.kind
    assert empty_str.dtype.kind == empty_bytes.str.decode("ascii").dtype.kind

    assert_equal(empty_bool, empty.str.contains("a"))
    assert_equal(empty_bool, empty.str.endswith("a"))
    assert_equal(empty_bool, empty.str.match("^a"))
    assert_equal(empty_bool, empty.str.startswith("a"))
    assert_equal(empty_bool, empty.str.isalnum())
    assert_equal(empty_bool, empty.str.isalpha())
    assert_equal(empty_bool, empty.str.isdecimal())
    assert_equal(empty_bool, empty.str.isdigit())
    assert_equal(empty_bool, empty.str.islower())
    assert_equal(empty_bool, empty.str.isnumeric())
    assert_equal(empty_bool, empty.str.isspace())
    assert_equal(empty_bool, empty.str.istitle())
    assert_equal(empty_bool, empty.str.isupper())
    assert_equal(empty_bytes, empty.str.encode("ascii"))
    assert_equal(empty_int, empty.str.count("a"))
    assert_equal(empty_int, empty.str.find("a"))
    assert_equal(empty_int, empty.str.len())
    assert_equal(empty_int, empty.str.rfind("a"))
    assert_equal(empty_str, empty.str.capitalize())
    assert_equal(empty_str, empty.str.center(42))
    assert_equal(empty_str, empty.str.get(0))
    assert_equal(empty_str, empty.str.lower())
    assert_equal(empty_str, empty.str.lstrip())
    assert_equal(empty_str, empty.str.pad(42))
    assert_equal(empty_str, empty.str.repeat(3))
    assert_equal(empty_str, empty.str.replace("a", "b"))
    assert_equal(empty_str, empty.str.rstrip())
    assert_equal(empty_str, empty.str.slice(step=1))
    assert_equal(empty_str, empty.str.slice(stop=1))
    assert_equal(empty_str, empty.str.strip())
    assert_equal(empty_str, empty.str.swapcase())
    assert_equal(empty_str, empty.str.title())
    assert_equal(empty_str, empty.str.upper())
    assert_equal(empty_str, empty.str.wrap(42))
    assert_equal(empty_str, empty_bytes.str.decode("ascii"))

    table = str.maketrans("a", "b")
    assert empty_str.dtype.kind == empty.str.translate(table).dtype.kind
    assert_equal(empty_str, empty.str.translate(table))


@pytest.mark.parametrize(
    ["func", "expected"],
    [
        pytest.param(
            lambda x: x.str.isalnum(),
            [True, True, True, True, True, False, True, True, False, False],
            id="isalnum",
        ),
        pytest.param(
            lambda x: x.str.isalpha(),
            [True, True, True, False, False, False, True, False, False, False],
            id="isalpha",
        ),
        pytest.param(
            lambda x: x.str.isdigit(),
            [False, False, False, True, False, False, False, True, False, False],
            id="isdigit",
        ),
        pytest.param(
            lambda x: x.str.islower(),
            [False, True, False, False, False, False, False, False, False, False],
            id="islower",
        ),
        pytest.param(
            lambda x: x.str.isspace(),
            [False, False, False, False, False, False, False, False, False, True],
            id="isspace",
        ),
        pytest.param(
            lambda x: x.str.istitle(),
            [True, False, True, False, True, False, False, False, False, False],
            id="istitle",
        ),
        pytest.param(
            lambda x: x.str.isupper(),
            [True, False, False, False, True, False, True, False, False, False],
            id="isupper",
        ),
    ],
)
def test_ismethods(
    dtype, func: Callable[[xr.DataArray], xr.DataArray], expected: list[bool]
) -> None:
    values = xr.DataArray(
        ["A", "b", "Xy", "4", "3A", "", "TT", "55", "-", "  "]
    ).astype(dtype)

    expected_da = xr.DataArray(expected)
    actual = func(values)

    assert actual.dtype == expected_da.dtype
    assert_equal(actual, expected_da)


def test_isnumeric() -> None:
    # 0x00bc: ¼ VULGAR FRACTION ONE QUARTER
    # 0x2605: ★ not number
    # 0x1378: ፸ ETHIOPIC NUMBER SEVENTY
    # 0xFF13: ３ Em 3
    values = xr.DataArray(["A", "3", "¼", "★", "፸", "３", "four"])
    exp_numeric = xr.DataArray([False, True, True, False, True, True, False])
    exp_decimal = xr.DataArray([False, True, False, False, False, True, False])

    res_numeric = values.str.isnumeric()
    res_decimal = values.str.isdecimal()

    assert res_numeric.dtype == exp_numeric.dtype
    assert res_decimal.dtype == exp_decimal.dtype

    assert_equal(res_numeric, exp_numeric)
    assert_equal(res_decimal, exp_decimal)


def test_len(dtype) -> None:
    values = ["foo", "fooo", "fooooo", "fooooooo"]
    result = xr.DataArray(values).astype(dtype).str.len()
    expected = xr.DataArray([len(x) for x in values])
    assert result.dtype == expected.dtype
    assert_equal(result, expected)


def test_find(dtype) -> None:
    values = xr.DataArray(["ABCDEFG", "BCDEFEF", "DEFGHIJEF", "EFGHEF", "XXX"])
    values = values.astype(dtype)

    result_0 = values.str.find("EF")
    result_1 = values.str.find("EF", side="left")
    expected_0 = xr.DataArray([4, 3, 1, 0, -1])
    expected_1 = xr.DataArray([v.find(dtype("EF")) for v in values.values])
    assert result_0.dtype == expected_0.dtype
    assert result_0.dtype == expected_1.dtype
    assert result_1.dtype == expected_0.dtype
    assert result_1.dtype == expected_1.dtype
    assert_equal(result_0, expected_0)
    assert_equal(result_0, expected_1)
    assert_equal(result_1, expected_0)
    assert_equal(result_1, expected_1)

    result_0 = values.str.rfind("EF")
    result_1 = values.str.find("EF", side="right")
    expected_0 = xr.DataArray([4, 5, 7, 4, -1])
    expected_1 = xr.DataArray([v.rfind(dtype("EF")) for v in values.values])
    assert result_0.dtype == expected_0.dtype
    assert result_0.dtype == expected_1.dtype
    assert result_1.dtype == expected_0.dtype
    assert result_1.dtype == expected_1.dtype
    assert_equal(result_0, expected_0)
    assert_equal(result_0, expected_1)
    assert_equal(result_1, expected_0)
    assert_equal(result_1, expected_1)

    result_0 = values.str.find("EF", 3)
    result_1 = values.str.find("EF", 3, side="left")
    expected_0 = xr.DataArray([4, 3, 7, 4, -1])
    expected_1 = xr.DataArray([v.find(dtype("EF"), 3) for v in values.values])
    assert result_0.dtype == expected_0.dtype
    assert result_0.dtype == expected_1.dtype
    assert result_1.dtype == expected_0.dtype
    assert result_1.dtype == expected_1.dtype
    assert_equal(result_0, expected_0)
    assert_equal(result_0, expected_1)
    assert_equal(result_1, expected_0)
    assert_equal(result_1, expected_1)

    result_0 = values.str.rfind("EF", 3)
    result_1 = values.str.find("EF", 3, side="right")
    expected_0 = xr.DataArray([4, 5, 7, 4, -1])
    expected_1 = xr.DataArray([v.rfind(dtype("EF"), 3) for v in values.values])
    assert result_0.dtype == expected_0.dtype
    assert result_0.dtype == expected_1.dtype
    assert result_1.dtype == expected_0.dtype
    assert result_1.dtype == expected_1.dtype
    assert_equal(result_0, expected_0)
    assert_equal(result_0, expected_1)
    assert_equal(result_1, expected_0)
    assert_equal(result_1, expected_1)

    result_0 = values.str.find("EF", 3, 6)
    result_1 = values.str.find("EF", 3, 6, side="left")
    expected_0 = xr.DataArray([4, 3, -1, 4, -1])
    expected_1 = xr.DataArray([v.find(dtype("EF"), 3, 6) for v in values.values])
    assert result_0.dtype == expected_0.dtype
    assert result_0.dtype == expected_1.dtype
    assert result_1.dtype == expected_0.dtype
    assert result_1.dtype == expected_1.dtype
    assert_equal(result_0, expected_0)
    assert_equal(result_0, expected_1)
    assert_equal(result_1, expected_0)
    assert_equal(result_1, expected_1)

    result_0 = values.str.rfind("EF", 3, 6)
    result_1 = values.str.find("EF", 3, 6, side="right")
    expected_0 = xr.DataArray([4, 3, -1, 4, -1])
    expected_1 = xr.DataArray([v.rfind(dtype("EF"), 3, 6) for v in values.values])
    assert result_0.dtype == expected_0.dtype
    assert result_0.dtype == expected_1.dtype
    assert result_1.dtype == expected_0.dtype
    assert result_1.dtype == expected_1.dtype
    assert_equal(result_0, expected_0)
    assert_equal(result_0, expected_1)
    assert_equal(result_1, expected_0)
    assert_equal(result_1, expected_1)


def test_find_broadcast(dtype) -> None:
    values = xr.DataArray(
        ["ABCDEFG", "BCDEFEF", "DEFGHIJEF", "EFGHEF", "XXX"], dims=["X"]
    )
    values = values.astype(dtype)
    sub = xr.DataArray(["EF", "BC", "XX"], dims=["Y"]).astype(dtype)
    start = xr.DataArray([0, 7], dims=["Z"])
    end = xr.DataArray([6, 9], dims=["Z"])

    result_0 = values.str.find(sub, start, end)
    result_1 = values.str.find(sub, start, end, side="left")
    expected = xr.DataArray(
        [
            [[4, -1], [1, -1], [-1, -1]],
            [[3, -1], [0, -1], [-1, -1]],
            [[1, 7], [-1, -1], [-1, -1]],
            [[0, -1], [-1, -1], [-1, -1]],
            [[-1, -1], [-1, -1], [0, -1]],
        ],
        dims=["X", "Y", "Z"],
    )

    assert result_0.dtype == expected.dtype
    assert result_1.dtype == expected.dtype
    assert_equal(result_0, expected)
    assert_equal(result_1, expected)

    result_0 = values.str.rfind(sub, start, end)
    result_1 = values.str.find(sub, start, end, side="right")
    expected = xr.DataArray(
        [
            [[4, -1], [1, -1], [-1, -1]],
            [[3, -1], [0, -1], [-1, -1]],
            [[1, 7], [-1, -1], [-1, -1]],
            [[4, -1], [-1, -1], [-1, -1]],
            [[-1, -1], [-1, -1], [1, -1]],
        ],
        dims=["X", "Y", "Z"],
    )

    assert result_0.dtype == expected.dtype
    assert result_1.dtype == expected.dtype
    assert_equal(result_0, expected)
    assert_equal(result_1, expected)


def test_index(dtype) -> None:
    s = xr.DataArray(["ABCDEFG", "BCDEFEF", "DEFGHIJEF", "EFGHEF"]).astype(dtype)

    result_0 = s.str.index("EF")
    result_1 = s.str.index("EF", side="left")
    expected = xr.DataArray([4, 3, 1, 0])
    assert result_0.dtype == expected.dtype
    assert result_1.dtype == expected.dtype
    assert_equal(result_0, expected)
    assert_equal(result_1, expected)

    result_0 = s.str.rindex("EF")
    result_1 = s.str.index("EF", side="right")
    expected = xr.DataArray([4, 5, 7, 4])
    assert result_0.dtype == expected.dtype
    assert result_1.dtype == expected.dtype
    assert_equal(result_0, expected)
    assert_equal(result_1, expected)

    result_0 = s.str.index("EF", 3)
    result_1 = s.str.index("EF", 3, side="left")
    expected = xr.DataArray([4, 3, 7, 4])
    assert result_0.dtype == expected.dtype
    assert result_1.dtype == expected.dtype
    assert_equal(result_0, expected)
    assert_equal(result_1, expected)

    result_0 = s.str.rindex("EF", 3)
    result_1 = s.str.index("EF", 3, side="right")
    expected = xr.DataArray([4, 5, 7, 4])
    assert result_0.dtype == expected.dtype
    assert result_1.dtype == expected.dtype
    assert_equal(result_0, expected)
    assert_equal(result_1, expected)

    result_0 = s.str.index("E", 4, 8)
    result_1 = s.str.index("E", 4, 8, side="left")
    expected = xr.DataArray([4, 5, 7, 4])
    assert result_0.dtype == expected.dtype
    assert result_1.dtype == expected.dtype
    assert_equal(result_0, expected)
    assert_equal(result_1, expected)

    result_0 = s.str.rindex("E", 0, 5)
    result_1 = s.str.index("E", 0, 5, side="right")
    expected = xr.DataArray([4, 3, 1, 4])
    assert result_0.dtype == expected.dtype
    assert result_1.dtype == expected.dtype
    assert_equal(result_0, expected)
    assert_equal(result_1, expected)

    matchtype = "subsection" if dtype == np.bytes_ else "substring"
    with pytest.raises(ValueError, match=f"{matchtype} not found"):
        s.str.index("DE")


def test_index_broadcast(dtype) -> None:
    values = xr.DataArray(
        ["ABCDEFGEFDBCA", "BCDEFEFEFDBC", "DEFBCGHIEFBC", "EFGHBCEFBCBCBCEF"],
        dims=["X"],
    )
    values = values.astype(dtype)
    sub = xr.DataArray(["EF", "BC"], dims=["Y"]).astype(dtype)
    start = xr.DataArray([0, 6], dims=["Z"])
    end = xr.DataArray([6, 12], dims=["Z"])

    result_0 = values.str.index(sub, start, end)
    result_1 = values.str.index(sub, start, end, side="left")
    expected = xr.DataArray(
        [[[4, 7], [1, 10]], [[3, 7], [0, 10]], [[1, 8], [3, 10]], [[0, 6], [4, 8]]],
        dims=["X", "Y", "Z"],
    )

    assert result_0.dtype == expected.dtype
    assert result_1.dtype == expected.dtype
    assert_equal(result_0, expected)
    assert_equal(result_1, expected)

    result_0 = values.str.rindex(sub, start, end)
    result_1 = values.str.index(sub, start, end, side="right")
    expected = xr.DataArray(
        [[[4, 7], [1, 10]], [[3, 7], [0, 10]], [[1, 8], [3, 10]], [[0, 6], [4, 10]]],
        dims=["X", "Y", "Z"],
    )

    assert result_0.dtype == expected.dtype
    assert result_1.dtype == expected.dtype
    assert_equal(result_0, expected)
    assert_equal(result_1, expected)


def test_translate() -> None:
    values = xr.DataArray(["abcdefg", "abcc", "cdddfg", "cdefggg"])
    table = str.maketrans("abc", "cde")
    result = values.str.translate(table)
    expected = xr.DataArray(["cdedefg", "cdee", "edddfg", "edefggg"])
    assert result.dtype == expected.dtype
    assert_equal(result, expected)


def test_pad_center_ljust_rjust(dtype) -> None:
    values = xr.DataArray(["a", "b", "c", "eeeee"]).astype(dtype)

    result = values.str.center(5)
    expected = xr.DataArray(["  a  ", "  b  ", "  c  ", "eeeee"]).astype(dtype)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)
    result = values.str.pad(5, side="both")
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    result = values.str.ljust(5)
    expected = xr.DataArray(["a    ", "b    ", "c    ", "eeeee"]).astype(dtype)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)
    result = values.str.pad(5, side="right")
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    result = values.str.rjust(5)
    expected = xr.DataArray(["    a", "    b", "    c", "eeeee"]).astype(dtype)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)
    result = values.str.pad(5, side="left")
    assert result.dtype == expected.dtype
    assert_equal(result, expected)


def test_pad_center_ljust_rjust_fillchar(dtype) -> None:
    values = xr.DataArray(["a", "bb", "cccc", "ddddd", "eeeeee"]).astype(dtype)

    result = values.str.center(5, fillchar="X")
    expected = xr.DataArray(["XXaXX", "XXbbX", "Xcccc", "ddddd", "eeeeee"]).astype(
        dtype
    )
    assert result.dtype == expected.dtype
    assert_equal(result, expected)
    result = values.str.pad(5, side="both", fillchar="X")
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    result = values.str.ljust(5, fillchar="X")
    expected = xr.DataArray(["aXXXX", "bbXXX", "ccccX", "ddddd", "eeeeee"]).astype(
        dtype
    )
    assert result.dtype == expected.dtype
    assert_equal(result, expected.astype(dtype))
    result = values.str.pad(5, side="right", fillchar="X")
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    result = values.str.rjust(5, fillchar="X")
    expected = xr.DataArray(["XXXXa", "XXXbb", "Xcccc", "ddddd", "eeeeee"]).astype(
        dtype
    )
    assert result.dtype == expected.dtype
    assert_equal(result, expected.astype(dtype))
    result = values.str.pad(5, side="left", fillchar="X")
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    # If fillchar is not a charatter, normal str raises TypeError
    # 'aaa'.ljust(5, 'XY')
    # TypeError: must be char, not str
    template = "fillchar must be a character, not {dtype}"

    with pytest.raises(TypeError, match=template.format(dtype="str")):
        values.str.center(5, fillchar="XY")

    with pytest.raises(TypeError, match=template.format(dtype="str")):
        values.str.ljust(5, fillchar="XY")

    with pytest.raises(TypeError, match=template.format(dtype="str")):
        values.str.rjust(5, fillchar="XY")

    with pytest.raises(TypeError, match=template.format(dtype="str")):
        values.str.pad(5, fillchar="XY")


def test_pad_center_ljust_rjust_broadcast(dtype) -> None:
    values = xr.DataArray(["a", "bb", "cccc", "ddddd", "eeeeee"], dims="X").astype(
        dtype
    )
    width = xr.DataArray([5, 4], dims="Y")
    fillchar = xr.DataArray(["X", "#"], dims="Y").astype(dtype)

    result = values.str.center(width, fillchar=fillchar)
    expected = xr.DataArray(
        [
            ["XXaXX", "#a##"],
            ["XXbbX", "#bb#"],
            ["Xcccc", "cccc"],
            ["ddddd", "ddddd"],
            ["eeeeee", "eeeeee"],
        ],
        dims=["X", "Y"],
    ).astype(dtype)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)
    result = values.str.pad(width, side="both", fillchar=fillchar)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    result = values.str.ljust(width, fillchar=fillchar)
    expected = xr.DataArray(
        [
            ["aXXXX", "a###"],
            ["bbXXX", "bb##"],
            ["ccccX", "cccc"],
            ["ddddd", "ddddd"],
            ["eeeeee", "eeeeee"],
        ],
        dims=["X", "Y"],
    ).astype(dtype)
    assert result.dtype == expected.dtype
    assert_equal(result, expected.astype(dtype))
    result = values.str.pad(width, side="right", fillchar=fillchar)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    result = values.str.rjust(width, fillchar=fillchar)
    expected = xr.DataArray(
        [
            ["XXXXa", "###a"],
            ["XXXbb", "##bb"],
            ["Xcccc", "cccc"],
            ["ddddd", "ddddd"],
            ["eeeeee", "eeeeee"],
        ],
        dims=["X", "Y"],
    ).astype(dtype)
    assert result.dtype == expected.dtype
    assert_equal(result, expected.astype(dtype))
    result = values.str.pad(width, side="left", fillchar=fillchar)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)


def test_zfill(dtype) -> None:
    values = xr.DataArray(["1", "22", "aaa", "333", "45678"]).astype(dtype)

    result = values.str.zfill(5)
    expected = xr.DataArray(["00001", "00022", "00aaa", "00333", "45678"]).astype(dtype)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    result = values.str.zfill(3)
    expected = xr.DataArray(["001", "022", "aaa", "333", "45678"]).astype(dtype)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)


def test_zfill_broadcast(dtype) -> None:
    values = xr.DataArray(["1", "22", "aaa", "333", "45678"]).astype(dtype)
    width = np.array([4, 5, 0, 3, 8])

    result = values.str.zfill(width)
    expected = xr.DataArray(["0001", "00022", "aaa", "333", "00045678"]).astype(dtype)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)


def test_slice(dtype) -> None:
    arr = xr.DataArray(["aafootwo", "aabartwo", "aabazqux"]).astype(dtype)

    result = arr.str.slice(2, 5)
    exp = xr.DataArray(["foo", "bar", "baz"]).astype(dtype)
    assert result.dtype == exp.dtype
    assert_equal(result, exp)

    for start, stop, step in [(0, 3, -1), (None, None, -1), (3, 10, 2), (3, 0, -1)]:
        try:
            result = arr.str[start:stop:step]
            expected = xr.DataArray([s[start:stop:step] for s in arr.values])
            assert_equal(result, expected.astype(dtype))
        except IndexError:
            print(f"failed on {start}:{stop}:{step}")
            raise


def test_slice_broadcast(dtype) -> None:
    arr = xr.DataArray(["aafootwo", "aabartwo", "aabazqux"]).astype(dtype)
    start = xr.DataArray([1, 2, 3])
    stop = 5

    result = arr.str.slice(start=start, stop=stop)
    exp = xr.DataArray(["afoo", "bar", "az"]).astype(dtype)
    assert result.dtype == exp.dtype
    assert_equal(result, exp)


def test_slice_replace(dtype) -> None:
    da = lambda x: xr.DataArray(x).astype(dtype)
    values = da(["short", "a bit longer", "evenlongerthanthat", ""])

    expected = da(["shrt", "a it longer", "evnlongerthanthat", ""])
    result = values.str.slice_replace(2, 3)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    expected = da(["shzrt", "a zit longer", "evznlongerthanthat", "z"])
    result = values.str.slice_replace(2, 3, "z")
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    expected = da(["shzort", "a zbit longer", "evzenlongerthanthat", "z"])
    result = values.str.slice_replace(2, 2, "z")
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    expected = da(["shzort", "a zbit longer", "evzenlongerthanthat", "z"])
    result = values.str.slice_replace(2, 1, "z")
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    expected = da(["shorz", "a bit longez", "evenlongerthanthaz", "z"])
    result = values.str.slice_replace(-1, None, "z")
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    expected = da(["zrt", "zer", "zat", "z"])
    result = values.str.slice_replace(None, -2, "z")
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    expected = da(["shortz", "a bit znger", "evenlozerthanthat", "z"])
    result = values.str.slice_replace(6, 8, "z")
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    expected = da(["zrt", "a zit longer", "evenlongzerthanthat", "z"])
    result = values.str.slice_replace(-10, 3, "z")
    assert result.dtype == expected.dtype
    assert_equal(result, expected)


def test_slice_replace_broadcast(dtype) -> None:
    values = xr.DataArray(["short", "a bit longer", "evenlongerthanthat", ""]).astype(
        dtype
    )
    start = 2
    stop = np.array([4, 5, None, 7])
    repl = "test"

    expected = xr.DataArray(["shtestt", "a test longer", "evtest", "test"]).astype(
        dtype
    )
    result = values.str.slice_replace(start, stop, repl)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)


def test_strip_lstrip_rstrip(dtype) -> None:
    values = xr.DataArray(["  aa   ", " bb \n", "cc  "]).astype(dtype)

    result = values.str.strip()
    expected = xr.DataArray(["aa", "bb", "cc"]).astype(dtype)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    result = values.str.lstrip()
    expected = xr.DataArray(["aa   ", "bb \n", "cc  "]).astype(dtype)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    result = values.str.rstrip()
    expected = xr.DataArray(["  aa", " bb", "cc"]).astype(dtype)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)


def test_strip_lstrip_rstrip_args(dtype) -> None:
    values = xr.DataArray(["xxABCxx", "xx BNSD", "LDFJH xx"]).astype(dtype)

    result = values.str.strip("x")
    expected = xr.DataArray(["ABC", " BNSD", "LDFJH "]).astype(dtype)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    result = values.str.lstrip("x")
    expected = xr.DataArray(["ABCxx", " BNSD", "LDFJH xx"]).astype(dtype)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    result = values.str.rstrip("x")
    expected = xr.DataArray(["xxABC", "xx BNSD", "LDFJH "]).astype(dtype)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)


def test_strip_lstrip_rstrip_broadcast(dtype) -> None:
    values = xr.DataArray(["xxABCxx", "yy BNSD", "LDFJH zz"]).astype(dtype)
    to_strip = xr.DataArray(["x", "y", "z"]).astype(dtype)

    result = values.str.strip(to_strip)
    expected = xr.DataArray(["ABC", " BNSD", "LDFJH "]).astype(dtype)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    result = values.str.lstrip(to_strip)
    expected = xr.DataArray(["ABCxx", " BNSD", "LDFJH zz"]).astype(dtype)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    result = values.str.rstrip(to_strip)
    expected = xr.DataArray(["xxABC", "yy BNSD", "LDFJH "]).astype(dtype)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)


def test_wrap() -> None:
    # test values are: two words less than width, two words equal to width,
    # two words greater than width, one word less than width, one word
    # equal to width, one word greater than width, multiple tokens with
    # trailing whitespace equal to width
    values = xr.DataArray(
        [
            "hello world",
            "hello world!",
            "hello world!!",
            "abcdefabcde",
            "abcdefabcdef",
            "abcdefabcdefa",
            "ab ab ab ab ",
            "ab ab ab ab a",
            "\t",
        ]
    )

    # expected values
    expected = xr.DataArray(
        [
            "hello world",
            "hello world!",
            "hello\nworld!!",
            "abcdefabcde",
            "abcdefabcdef",
            "abcdefabcdef\na",
            "ab ab ab ab",
            "ab ab ab ab\na",
            "",
        ]
    )

    result = values.str.wrap(12, break_long_words=True)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    # test with pre and post whitespace (non-unicode), NaN, and non-ascii
    # Unicode
    values = xr.DataArray(["  pre  ", "\xac\u20ac\U00008000 abadcafe"])
    expected = xr.DataArray(["  pre", "\xac\u20ac\U00008000 ab\nadcafe"])
    result = values.str.wrap(6)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)


def test_wrap_kwargs_passed() -> None:
    # GH4334

    values = xr.DataArray("  hello world  ")

    result = values.str.wrap(7)
    expected = xr.DataArray("  hello\nworld")
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    result = values.str.wrap(7, drop_whitespace=False)
    expected = xr.DataArray("  hello\n world\n  ")
    assert result.dtype == expected.dtype
    assert_equal(result, expected)


def test_get(dtype) -> None:
    values = xr.DataArray(["a_b_c", "c_d_e", "f_g_h"]).astype(dtype)

    result = values.str[2]
    expected = xr.DataArray(["b", "d", "g"]).astype(dtype)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    # bounds testing
    values = xr.DataArray(["1_2_3_4_5", "6_7_8_9_10", "11_12"]).astype(dtype)

    # positive index
    result = values.str[5]
    expected = xr.DataArray(["_", "_", ""]).astype(dtype)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    # negative index
    result = values.str[-6]
    expected = xr.DataArray(["_", "8", ""]).astype(dtype)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)


def test_get_default(dtype) -> None:
    # GH4334
    values = xr.DataArray(["a_b", "c", ""]).astype(dtype)

    result = values.str.get(2, "default")
    expected = xr.DataArray(["b", "default", "default"]).astype(dtype)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)


def test_get_broadcast(dtype) -> None:
    values = xr.DataArray(["a_b_c", "c_d_e", "f_g_h"], dims=["X"]).astype(dtype)
    inds = xr.DataArray([0, 2], dims=["Y"])

    result = values.str.get(inds)
    expected = xr.DataArray(
        [["a", "b"], ["c", "d"], ["f", "g"]], dims=["X", "Y"]
    ).astype(dtype)
    assert result.dtype == expected.dtype
    assert_equal(result, expected)


def test_encode_decode() -> None:
    data = xr.DataArray(["a", "b", "a\xe4"])
    encoded = data.str.encode("utf-8")
    decoded = encoded.str.decode("utf-8")
    assert data.dtype == decoded.dtype
    assert_equal(data, decoded)


def test_encode_decode_errors() -> None:
    encodeBase = xr.DataArray(["a", "b", "a\x9d"])

    msg = (
        r"'charmap' codec can't encode character '\\x9d' in position 1:"
        " character maps to <undefined>"
    )
    with pytest.raises(UnicodeEncodeError, match=msg):
        encodeBase.str.encode("cp1252")

    f = lambda x: x.encode("cp1252", "ignore")
    result = encodeBase.str.encode("cp1252", "ignore")
    expected = xr.DataArray([f(x) for x in encodeBase.values.tolist()])

    assert result.dtype == expected.dtype
    assert_equal(result, expected)

    decodeBase = xr.DataArray([b"a", b"b", b"a\x9d"])

    msg = (
        "'charmap' codec can't decode byte 0x9d in position 1:"
        " character maps to <undefined>"
    )
    with pytest.raises(UnicodeDecodeError, match=msg):
        decodeBase.str.decode("cp1252")

    f = lambda x: x.decode("cp1252", "ignore")
    result = decodeBase.str.decode("cp1252", "ignore")
    expected = xr.DataArray([f(x) for x in decodeBase.values.tolist()])

    assert result.dtype == expected.dtype
    assert_equal(result, expected)


def test_partition_whitespace(dtype) -> None:
    values = xr.DataArray(
        [
            ["abc def", "spam eggs swallow", "red_blue"],
            ["test0 test1 test2 test3", "", "abra ka da bra"],
        ],
        dims=["X", "Y"],
    ).astype(dtype)

    exp_part_dim_list = [
        [
            ["abc", " ", "def"],
            ["spam", " ", "eggs swallow"],
            ["red_blue", "", ""],
        ],
        [
            ["test0", " ", "test1 test2 test3"],
            ["", "", ""],
            ["abra", " ", "ka da bra"],
        ],
    ]

    exp_rpart_dim_list = [
        [
            ["abc", " ", "def"],
            ["spam eggs", " ", "swallow"],
            ["", "", "red_blue"],
        ],
        [
            ["test0 test1 test2", " ", "test3"],
            ["", "", ""],
            ["abra ka da", " ", "bra"],
        ],
    ]

    exp_part_dim = xr.DataArray(exp_part_dim_list, dims=["X", "Y", "ZZ"]).astype(dtype)
    exp_rpart_dim = xr.DataArray(exp_rpart_dim_list, dims=["X", "Y", "ZZ"]).astype(
        dtype
    )

    res_part_dim = values.str.partition(dim="ZZ")
    res_rpart_dim = values.str.rpartition(dim="ZZ")

    assert res_part_dim.dtype == exp_part_dim.dtype
    assert res_rpart_dim.dtype == exp_rpart_dim.dtype

    assert_equal(res_part_dim, exp_part_dim)
    assert_equal(res_rpart_dim, exp_rpart_dim)


def test_partition_comma(dtype) -> None:
    values = xr.DataArray(
        [
            ["abc, def", "spam, eggs, swallow", "red_blue"],
            ["test0, test1, test2, test3", "", "abra, ka, da, bra"],
        ],
        dims=["X", "Y"],
    ).astype(dtype)

    exp_part_dim_list = [
        [
            ["abc", ", ", "def"],
            ["spam", ", ", "eggs, swallow"],
            ["red_blue", "", ""],
        ],
        [
            ["test0", ", ", "test1, test2, test3"],
            ["", "", ""],
            ["abra", ", ", "ka, da, bra"],
        ],
    ]

    exp_rpart_dim_list = [
        [
            ["abc", ", ", "def"],
            ["spam, eggs", ", ", "swallow"],
            ["", "", "red_blue"],
        ],
        [
            ["test0, test1, test2", ", ", "test3"],
            ["", "", ""],
            ["abra, ka, da", ", ", "bra"],
        ],
    ]

    exp_part_dim = xr.DataArray(exp_part_dim_list, dims=["X", "Y", "ZZ"]).astype(dtype)
    exp_rpart_dim = xr.DataArray(exp_rpart_dim_list, dims=["X", "Y", "ZZ"]).astype(
        dtype
    )

    res_part_dim = values.str.partition(sep=", ", dim="ZZ")
    res_rpart_dim = values.str.rpartition(sep=", ", dim="ZZ")

    assert res_part_dim.dtype == exp_part_dim.dtype
    assert res_rpart_dim.dtype == exp_rpart_dim.dtype

    assert_equal(res_part_dim, exp_part_dim)
    assert_equal(res_rpart_dim, exp_rpart_dim)


def test_partition_empty(dtype) -> None:
    values = xr.DataArray([], dims=["X"]).astype(dtype)
    expected = xr.DataArray(np.zeros((0, 0)), dims=["X", "ZZ"]).astype(dtype)

    res = values.str.partition(sep=", ", dim="ZZ")

    assert res.dtype == expected.dtype
    assert_equal(res, expected)


@pytest.mark.parametrize(
    ["func", "expected"],
    [
        pytest.param(
            lambda x: x.str.split(dim=None),
            [
                [["abc", "def"], ["spam", "eggs", "swallow"], ["red_blue"]],
                [["test0", "test1", "test2", "test3"], [], ["abra", "ka", "da", "bra"]],
            ],
            id="split_full",
        ),
        pytest.param(
            lambda x: x.str.rsplit(dim=None),
            [
                [["abc", "def"], ["spam", "eggs", "swallow"], ["red_blue"]],
                [["test0", "test1", "test2", "test3"], [], ["abra", "ka", "da", "bra"]],
            ],
            id="rsplit_full",
        ),
        pytest.param(
            lambda x: x.str.split(dim=None, maxsplit=1),
            [
                [["abc", "def"], ["spam", "eggs\tswallow"], ["red_blue"]],
                [["test0", "test1\ntest2\n\ntest3"], [], ["abra", "ka\nda\tbra"]],
            ],
            id="split_1",
        ),
        pytest.param(
            lambda x: x.str.rsplit(dim=None, maxsplit=1),
            [
                [["abc", "def"], ["spam\t\teggs", "swallow"], ["red_blue"]],
                [["test0\ntest1\ntest2", "test3"], [], ["abra  ka\nda", "bra"]],
            ],
            id="rsplit_1",
        ),
    ],
)
def test_split_whitespace_nodim(
    dtype, func: Callable[[xr.DataArray], xr.DataArray], expected: xr.DataArray
) -> None:
    values = xr.DataArray(
        [
            ["abc def", "spam\t\teggs\tswallow", "red_blue"],
            ["test0\ntest1\ntest2\n\ntest3", "", "abra  ka\nda\tbra"],
        ],
        dims=["X", "Y"],
    ).astype(dtype)

    expected_dtype = [[[dtype(x) for x in y] for y in z] for z in expected]
    expected_np = np.array(expected_dtype, dtype=np.object_)
    expected_da = xr.DataArray(expected_np, dims=["X", "Y"])

    actual = func(values)

    assert actual.dtype == expected_da.dtype
    assert_equal(actual, expected_da)


@pytest.mark.parametrize(
    ["func", "expected"],
    [
        pytest.param(
            lambda x: x.str.split(dim="ZZ"),
            [
                [
                    ["abc", "def", "", ""],
                    ["spam", "eggs", "swallow", ""],
                    ["red_blue", "", "", ""],
                ],
                [
                    ["test0", "test1", "test2", "test3"],
                    ["", "", "", ""],
                    ["abra", "ka", "da", "bra"],
                ],
            ],
            id="split_full",
        ),
        pytest.param(
            lambda x: x.str.rsplit(dim="ZZ"),
            [
                [
                    ["", "", "abc", "def"],
                    ["", "spam", "eggs", "swallow"],
                    ["", "", "", "red_blue"],
                ],
                [
                    ["test0", "test1", "test2", "test3"],
                    ["", "", "", ""],
                    ["abra", "ka", "da", "bra"],
                ],
            ],
            id="rsplit_full",
        ),
        pytest.param(
            lambda x: x.str.split(dim="ZZ", maxsplit=1),
            [
                [["abc", "def"], ["spam", "eggs\tswallow"], ["red_blue", ""]],
                [["test0", "test1\ntest2\n\ntest3"], ["", ""], ["abra", "ka\nda\tbra"]],
            ],
            id="split_1",
        ),
        pytest.param(
            lambda x: x.str.rsplit(dim="ZZ", maxsplit=1),
            [
                [["abc", "def"], ["spam\t\teggs", "swallow"], ["", "red_blue"]],
                [["test0\ntest1\ntest2", "test3"], ["", ""], ["abra  ka\nda", "bra"]],
            ],
            id="rsplit_1",
        ),
    ],
)
def test_split_whitespace_dim(
    dtype, func: Callable[[xr.DataArray], xr.DataArray], expected: xr.DataArray
) -> None:
    values = xr.DataArray(
        [
            ["abc def", "spam\t\teggs\tswallow", "red_blue"],
            ["test0\ntest1\ntest2\n\ntest3", "", "abra  ka\nda\tbra"],
        ],
        dims=["X", "Y"],
    ).astype(dtype)

    expected_dtype = [[[dtype(x) for x in y] for y in z] for z in expected]
    expected_np = np.array(expected_dtype, dtype=np.object_)
    expected_da = xr.DataArray(expected_np, dims=["X", "Y", "ZZ"]).astype(dtype)

    actual = func(values)

    assert actual.dtype == expected_da.dtype
    assert_equal(actual, expected_da)


@pytest.mark.parametrize(
    ["func", "expected"],
    [
        pytest.param(
            lambda x: x.str.split(sep=",", dim=None),
            [
                [["abc", "def"], ["spam", "", "eggs", "swallow"], ["red_blue"]],
                [
                    ["test0", "test1", "test2", "test3"],
                    [""],
                    ["abra", "ka", "da", "bra"],
                ],
            ],
            id="split_full",
        ),
        pytest.param(
            lambda x: x.str.rsplit(sep=",", dim=None),
            [
                [["abc", "def"], ["spam", "", "eggs", "swallow"], ["red_blue"]],
                [
                    ["test0", "test1", "test2", "test3"],
                    [""],
                    ["abra", "ka", "da", "bra"],
                ],
            ],
            id="rsplit_full",
        ),
        pytest.param(
            lambda x: x.str.split(sep=",", dim=None, maxsplit=1),
            [
                [["abc", "def"], ["spam", ",eggs,swallow"], ["red_blue"]],
                [["test0", "test1,test2,test3"], [""], ["abra", "ka,da,bra"]],
            ],
            id="split_1",
        ),
        pytest.param(
            lambda x: x.str.rsplit(sep=",", dim=None, maxsplit=1),
            [
                [["abc", "def"], ["spam,,eggs", "swallow"], ["red_blue"]],
                [["test0,test1,test2", "test3"], [""], ["abra,ka,da", "bra"]],
            ],
            id="rsplit_1",
        ),
        pytest.param(
            lambda x: x.str.split(sep=",", dim=None, maxsplit=10),
            [
                [["abc", "def"], ["spam", "", "eggs", "swallow"], ["red_blue"]],
                [
                    ["test0", "test1", "test2", "test3"],
                    [""],
                    ["abra", "ka", "da", "bra"],
                ],
            ],
            id="split_10",
        ),
        pytest.param(
            lambda x: x.str.rsplit(sep=",", dim=None, maxsplit=10),
            [
                [["abc", "def"], ["spam", "", "eggs", "swallow"], ["red_blue"]],
                [
                    ["test0", "test1", "test2", "test3"],
                    [""],
                    ["abra", "ka", "da", "bra"],
                ],
            ],
            id="rsplit_10",
        ),
    ],
)
def test_split_comma_nodim(
    dtype, func: Callable[[xr.DataArray], xr.DataArray], expected: xr.DataArray
) -> None:
    values = xr.DataArray(
        [
            ["abc,def", "spam,,eggs,swallow", "red_blue"],
            ["test0,test1,test2,test3", "", "abra,ka,da,bra"],
        ],
        dims=["X", "Y"],
    ).astype(dtype)

    expected_dtype = [[[dtype(x) for x in y] for y in z] for z in expected]
    expected_np = np.array(expected_dtype, dtype=np.object_)
    expected_da = xr.DataArray(expected_np, dims=["X", "Y"])

    actual = func(values)

    assert actual.dtype == expected_da.dtype
    assert_equal(actual, expected_da)


@pytest.mark.parametrize(
    ["func", "expected"],
    [
        pytest.param(
            lambda x: x.str.split(sep=",", dim="ZZ"),
            [
                [
                    ["abc", "def", "", ""],
                    ["spam", "", "eggs", "swallow"],
                    ["red_blue", "", "", ""],
                ],
                [
                    ["test0", "test1", "test2", "test3"],
                    ["", "", "", ""],
                    ["abra", "ka", "da", "bra"],
                ],
            ],
            id="split_full",
        ),
        pytest.param(
            lambda x: x.str.rsplit(sep=",", dim="ZZ"),
            [
                [
                    ["", "", "abc", "def"],
                    ["spam", "", "eggs", "swallow"],
                    ["", "", "", "red_blue"],
                ],
                [
                    ["test0", "test1", "test2", "test3"],
                    ["", "", "", ""],
                    ["abra", "ka", "da", "bra"],
                ],
            ],
            id="rsplit_full",
        ),
        pytest.param(
            lambda x: x.str.split(sep=",", dim="ZZ", maxsplit=1),
            [
                [["abc", "def"], ["spam", ",eggs,swallow"], ["red_blue", ""]],
                [["test0", "test1,test2,test3"], ["", ""], ["abra", "ka,da,bra"]],
            ],
            id="split_1",
        ),
        pytest.param(
            lambda x: x.str.rsplit(sep=",", dim="ZZ", maxsplit=1),
            [
                [["abc", "def"], ["spam,,eggs", "swallow"], ["", "red_blue"]],
                [["test0,test1,test2", "test3"], ["", ""], ["abra,ka,da", "bra"]],
            ],
            id="rsplit_1",
        ),
        pytest.param(
            lambda x: x.str.split(sep=",", dim="ZZ", maxsplit=10),
            [
                [
                    ["abc", "def", "", ""],
                    ["spam", "", "eggs", "swallow"],
                    ["red_blue", "", "", ""],
                ],
                [
                    ["test0", "test1", "test2", "test3"],
                    ["", "", "", ""],
                    ["abra", "ka", "da", "bra"],
                ],
            ],
            id="split_10",
        ),
        pytest.param(
            lambda x: x.str.rsplit(sep=",", dim="ZZ", maxsplit=10),
            [
                [
                    ["", "", "abc", "def"],
                    ["spam", "", "eggs", "swallow"],
                    ["", "", "", "red_blue"],
                ],
                [
                    ["test0", "test1", "test2", "test3"],
                    ["", "", "", ""],
                    ["abra", "ka", "da", "bra"],
                ],
            ],
            id="rsplit_10",
        ),
    ],
)
def test_split_comma_dim(
    dtype, func: Callable[[xr.DataArray], xr.DataArray], expected: xr.DataArray
) -> None:
    values = xr.DataArray(
        [
            ["abc,def", "spam,,eggs,swallow", "red_blue"],
            ["test0,test1,test2,test3", "", "abra,ka,da,bra"],
        ],
        dims=["X", "Y"],
    ).astype(dtype)

    expected_dtype = [[[dtype(x) for x in y] for y in z] for z in expected]
    expected_np = np.array(expected_dtype, dtype=np.object_)
    expected_da = xr.DataArray(expected_np, dims=["X", "Y", "ZZ"]).astype(dtype)

    actual = func(values)

    assert actual.dtype == expected_da.dtype
    assert_equal(actual, expected_da)


def test_splitters_broadcast(dtype) -> None:
    values = xr.DataArray(
        ["ab cd,de fg", "spam, ,eggs swallow", "red_blue"],
        dims=["X"],
    ).astype(dtype)

    sep = xr.DataArray(
        [" ", ","],
        dims=["Y"],
    ).astype(dtype)

    expected_left = xr.DataArray(
        [
            [["ab", "cd,de fg"], ["ab cd", "de fg"]],
            [["spam,", ",eggs swallow"], ["spam", " ,eggs swallow"]],
            [["red_blue", ""], ["red_blue", ""]],
        ],
        dims=["X", "Y", "ZZ"],
    ).astype(dtype)
    expected_right = xr.DataArray(
        [
            [["ab cd,de", "fg"], ["ab cd", "de fg"]],
            [["spam, ,eggs", "swallow"], ["spam, ", "eggs swallow"]],
            [["", "red_blue"], ["", "red_blue"]],
        ],
        dims=["X", "Y", "ZZ"],
    ).astype(dtype)

    res_left = values.str.split(dim="ZZ", sep=sep, maxsplit=1)
    res_right = values.str.rsplit(dim="ZZ", sep=sep, maxsplit=1)

    # assert res_left.dtype == expected_left.dtype
    # assert res_right.dtype == expected_right.dtype

    assert_equal(res_left, expected_left)
    assert_equal(res_right, expected_right)

    expected_left = xr.DataArray(
        [
            [["ab", " ", "cd,de fg"], ["ab cd", ",", "de fg"]],
            [["spam,", " ", ",eggs swallow"], ["spam", ",", " ,eggs swallow"]],
            [["red_blue", "", ""], ["red_blue", "", ""]],
        ],
        dims=["X", "Y", "ZZ"],
    ).astype(dtype)
    expected_right = xr.DataArray(
        [
            [["ab", " ", "cd,de fg"], ["ab cd", ",", "de fg"]],
            [["spam,", " ", ",eggs swallow"], ["spam", ",", " ,eggs swallow"]],
            [["red_blue", "", ""], ["red_blue", "", ""]],
        ],
        dims=["X", "Y", "ZZ"],
    ).astype(dtype)

    res_left = values.str.partition(dim="ZZ", sep=sep)
    res_right = values.str.partition(dim="ZZ", sep=sep)

    # assert res_left.dtype == expected_left.dtype
    # assert res_right.dtype == expected_right.dtype

    assert_equal(res_left, expected_left)
    assert_equal(res_right, expected_right)


def test_split_empty(dtype) -> None:
    values = xr.DataArray([], dims=["X"]).astype(dtype)
    expected = xr.DataArray(np.zeros((0, 0)), dims=["X", "ZZ"]).astype(dtype)

    res = values.str.split(sep=", ", dim="ZZ")

    assert res.dtype == expected.dtype
    assert_equal(res, expected)


def test_get_dummies(dtype) -> None:
    values_line = xr.DataArray(
        [["a|ab~abc|abc", "ab", "a||abc|abcd"], ["abcd|ab|a", "abc|ab~abc", "|a"]],
        dims=["X", "Y"],
    ).astype(dtype)
    values_comma = xr.DataArray(
        [["a~ab|abc~~abc", "ab", "a~abc~abcd"], ["abcd~ab~a", "abc~ab|abc", "~a"]],
        dims=["X", "Y"],
    ).astype(dtype)

    vals_line = np.array(["a", "ab", "abc", "abcd", "ab~abc"]).astype(dtype)
    vals_comma = np.array(["a", "ab", "abc", "abcd", "ab|abc"]).astype(dtype)
    expected_list = [
        [
            [True, False, True, False, True],
            [False, True, False, False, False],
            [True, False, True, True, False],
        ],
        [
            [True, True, False, True, False],
            [False, False, True, False, True],
            [True, False, False, False, False],
        ],
    ]
    expected_np = np.array(expected_list)
    expected = xr.DataArray(expected_np, dims=["X", "Y", "ZZ"])
    targ_line = expected.copy()
    targ_comma = expected.copy()
    targ_line.coords["ZZ"] = vals_line
    targ_comma.coords["ZZ"] = vals_comma

    res_default = values_line.str.get_dummies(dim="ZZ")
    res_line = values_line.str.get_dummies(dim="ZZ", sep="|")
    res_comma = values_comma.str.get_dummies(dim="ZZ", sep="~")

    assert res_default.dtype == targ_line.dtype
    assert res_line.dtype == targ_line.dtype
    assert res_comma.dtype == targ_comma.dtype

    assert_equal(res_default, targ_line)
    assert_equal(res_line, targ_line)
    assert_equal(res_comma, targ_comma)


def test_get_dummies_broadcast(dtype) -> None:
    values = xr.DataArray(
        ["x~x|x~x", "x", "x|x~x", "x~x"],
        dims=["X"],
    ).astype(dtype)

    sep = xr.DataArray(
        ["|", "~"],
        dims=["Y"],
    ).astype(dtype)

    expected_list = [
        [[False, False, True], [True, True, False]],
        [[True, False, False], [True, False, False]],
        [[True, False, True], [True, True, False]],
        [[False, False, True], [True, False, False]],
    ]
    expected_np = np.array(expected_list)
    expected = xr.DataArray(expected_np, dims=["X", "Y", "ZZ"])
    expected.coords["ZZ"] = np.array(["x", "x|x", "x~x"]).astype(dtype)

    res = values.str.get_dummies(dim="ZZ", sep=sep)

    assert res.dtype == expected.dtype

    assert_equal(res, expected)


def test_get_dummies_empty(dtype) -> None:
    values = xr.DataArray([], dims=["X"]).astype(dtype)
    expected = xr.DataArray(np.zeros((0, 0)), dims=["X", "ZZ"]).astype(dtype)

    res = values.str.get_dummies(dim="ZZ")

    assert res.dtype == expected.dtype
    assert_equal(res, expected)


def test_splitters_empty_str(dtype) -> None:
    values = xr.DataArray(
        [["", "", ""], ["", "", ""]],
        dims=["X", "Y"],
    ).astype(dtype)

    targ_partition_dim = xr.DataArray(
        [
            [["", "", ""], ["", "", ""], ["", "", ""]],
            [["", "", ""], ["", "", ""], ["", "", ""]],
        ],
        dims=["X", "Y", "ZZ"],
    ).astype(dtype)

    targ_partition_none_list = [
        [["", "", ""], ["", "", ""], ["", "", ""]],
        [["", "", ""], ["", "", ""], ["", "", "", ""]],
    ]
    targ_partition_none_list = [
        [[dtype(x) for x in y] for y in z] for z in targ_partition_none_list
    ]
    targ_partition_none_np = np.array(targ_partition_none_list, dtype=np.object_)
    del targ_partition_none_np[-1, -1][-1]
    targ_partition_none = xr.DataArray(
        targ_partition_none_np,
        dims=["X", "Y"],
    )

    targ_split_dim = xr.DataArray(
        [[[""], [""], [""]], [[""], [""], [""]]],
        dims=["X", "Y", "ZZ"],
    ).astype(dtype)
    targ_split_none = xr.DataArray(
        np.array([[[], [], []], [[], [], [""]]], dtype=np.object_),
        dims=["X", "Y"],
    )
    del targ_split_none.data[-1, -1][-1]

    res_partition_dim = values.str.partition(dim="ZZ")
    res_rpartition_dim = values.str.rpartition(dim="ZZ")
    res_partition_none = values.str.partition(dim=None)
    res_rpartition_none = values.str.rpartition(dim=None)

    res_split_dim = values.str.split(dim="ZZ")
    res_rsplit_dim = values.str.rsplit(dim="ZZ")
    res_split_none = values.str.split(dim=None)
    res_rsplit_none = values.str.rsplit(dim=None)

    res_dummies = values.str.rsplit(dim="ZZ")

    assert res_partition_dim.dtype == targ_partition_dim.dtype
    assert res_rpartition_dim.dtype == targ_partition_dim.dtype
    assert res_partition_none.dtype == targ_partition_none.dtype
    assert res_rpartition_none.dtype == targ_partition_none.dtype

    assert res_split_dim.dtype == targ_split_dim.dtype
    assert res_rsplit_dim.dtype == targ_split_dim.dtype
    assert res_split_none.dtype == targ_split_none.dtype
    assert res_rsplit_none.dtype == targ_split_none.dtype

    assert res_dummies.dtype == targ_split_dim.dtype

    assert_equal(res_partition_dim, targ_partition_dim)
    assert_equal(res_rpartition_dim, targ_partition_dim)
    assert_equal(res_partition_none, targ_partition_none)
    assert_equal(res_rpartition_none, targ_partition_none)

    assert_equal(res_split_dim, targ_split_dim)
    assert_equal(res_rsplit_dim, targ_split_dim)
    assert_equal(res_split_none, targ_split_none)
    assert_equal(res_rsplit_none, targ_split_none)

    assert_equal(res_dummies, targ_split_dim)


def test_cat_str(dtype) -> None:
    values_1 = xr.DataArray(
        [["a", "bb", "cccc"], ["ddddd", "eeee", "fff"]],
        dims=["X", "Y"],
    ).astype(dtype)
    values_2 = "111"

    targ_blank = xr.DataArray(
        [["a111", "bb111", "cccc111"], ["ddddd111", "eeee111", "fff111"]],
        dims=["X", "Y"],
    ).astype(dtype)

    targ_space = xr.DataArray(
        [["a 111", "bb 111", "cccc 111"], ["ddddd 111", "eeee 111", "fff 111"]],
        dims=["X", "Y"],
    ).astype(dtype)

    targ_bars = xr.DataArray(
        [["a||111", "bb||111", "cccc||111"], ["ddddd||111", "eeee||111", "fff||111"]],
        dims=["X", "Y"],
    ).astype(dtype)

    targ_comma = xr.DataArray(
        [["a, 111", "bb, 111", "cccc, 111"], ["ddddd, 111", "eeee, 111", "fff, 111"]],
        dims=["X", "Y"],
    ).astype(dtype)

    res_blank = values_1.str.cat(values_2)
    res_add = values_1.str + values_2
    res_space = values_1.str.cat(values_2, sep=" ")
    res_bars = values_1.str.cat(values_2, sep="||")
    res_comma = values_1.str.cat(values_2, sep=", ")

    assert res_blank.dtype == targ_blank.dtype
    assert res_add.dtype == targ_blank.dtype
    assert res_space.dtype == targ_space.dtype
    assert res_bars.dtype == targ_bars.dtype
    assert res_comma.dtype == targ_comma.dtype

    assert_equal(res_blank, targ_blank)
    assert_equal(res_add, targ_blank)
    assert_equal(res_space, targ_space)
    assert_equal(res_bars, targ_bars)
    assert_equal(res_comma, targ_comma)


def test_cat_uniform(dtype) -> None:
    values_1 = xr.DataArray(
        [["a", "bb", "cccc"], ["ddddd", "eeee", "fff"]],
        dims=["X", "Y"],
    ).astype(dtype)
    values_2 = xr.DataArray(
        [["11111", "222", "33"], ["4", "5555", "66"]],
        dims=["X", "Y"],
    )

    targ_blank = xr.DataArray(
        [["a11111", "bb222", "cccc33"], ["ddddd4", "eeee5555", "fff66"]],
        dims=["X", "Y"],
    ).astype(dtype)

    targ_space = xr.DataArray(
        [["a 11111", "bb 222", "cccc 33"], ["ddddd 4", "eeee 5555", "fff 66"]],
        dims=["X", "Y"],
    ).astype(dtype)

    targ_bars = xr.DataArray(
        [["a||11111", "bb||222", "cccc||33"], ["ddddd||4", "eeee||5555", "fff||66"]],
        dims=["X", "Y"],
    ).astype(dtype)

    targ_comma = xr.DataArray(
        [["a, 11111", "bb, 222", "cccc, 33"], ["ddddd, 4", "eeee, 5555", "fff, 66"]],
        dims=["X", "Y"],
    ).astype(dtype)

    res_blank = values_1.str.cat(values_2)
    res_add = values_1.str + values_2
    res_space = values_1.str.cat(values_2, sep=" ")
    res_bars = values_1.str.cat(values_2, sep="||")
    res_comma = values_1.str.cat(values_2, sep=", ")

    assert res_blank.dtype == targ_blank.dtype
    assert res_add.dtype == targ_blank.dtype
    assert res_space.dtype == targ_space.dtype
    assert res_bars.dtype == targ_bars.dtype
    assert res_comma.dtype == targ_comma.dtype

    assert_equal(res_blank, targ_blank)
    assert_equal(res_add, targ_blank)
    assert_equal(res_space, targ_space)
    assert_equal(res_bars, targ_bars)
    assert_equal(res_comma, targ_comma)


def test_cat_broadcast_right(dtype) -> None:
    values_1 = xr.DataArray(
        [["a", "bb", "cccc"], ["ddddd", "eeee", "fff"]],
        dims=["X", "Y"],
    ).astype(dtype)
    values_2 = xr.DataArray(
        ["11111", "222", "33"],
        dims=["Y"],
    )

    targ_blank = xr.DataArray(
        [["a11111", "bb222", "cccc33"], ["ddddd11111", "eeee222", "fff33"]],
        dims=["X", "Y"],
    ).astype(dtype)

    targ_space = xr.DataArray(
        [["a 11111", "bb 222", "cccc 33"], ["ddddd 11111", "eeee 222", "fff 33"]],
        dims=["X", "Y"],
    ).astype(dtype)

    targ_bars = xr.DataArray(
        [["a||11111", "bb||222", "cccc||33"], ["ddddd||11111", "eeee||222", "fff||33"]],
        dims=["X", "Y"],
    ).astype(dtype)

    targ_comma = xr.DataArray(
        [["a, 11111", "bb, 222", "cccc, 33"], ["ddddd, 11111", "eeee, 222", "fff, 33"]],
        dims=["X", "Y"],
    ).astype(dtype)

    res_blank = values_1.str.cat(values_2)
    res_add = values_1.str + values_2
    res_space = values_1.str.cat(values_2, sep=" ")
    res_bars = values_1.str.cat(values_2, sep="||")
    res_comma = values_1.str.cat(values_2, sep=", ")

    assert res_blank.dtype == targ_blank.dtype
    assert res_add.dtype == targ_blank.dtype
    assert res_space.dtype == targ_space.dtype
    assert res_bars.dtype == targ_bars.dtype
    assert res_comma.dtype == targ_comma.dtype

    assert_equal(res_blank, targ_blank)
    assert_equal(res_add, targ_blank)
    assert_equal(res_space, targ_space)
    assert_equal(res_bars, targ_bars)
    assert_equal(res_comma, targ_comma)


def test_cat_broadcast_left(dtype) -> None:
    values_1 = xr.DataArray(
        ["a", "bb", "cccc"],
        dims=["Y"],
    ).astype(dtype)
    values_2 = xr.DataArray(
        [["11111", "222", "33"], ["4", "5555", "66"]],
        dims=["X", "Y"],
    )

    targ_blank = (
        xr.DataArray(
            [["a11111", "bb222", "cccc33"], ["a4", "bb5555", "cccc66"]],
            dims=["X", "Y"],
        )
        .astype(dtype)
        .T
    )

    targ_space = (
        xr.DataArray(
            [["a 11111", "bb 222", "cccc 33"], ["a 4", "bb 5555", "cccc 66"]],
            dims=["X", "Y"],
        )
        .astype(dtype)
        .T
    )

    targ_bars = (
        xr.DataArray(
            [["a||11111", "bb||222", "cccc||33"], ["a||4", "bb||5555", "cccc||66"]],
            dims=["X", "Y"],
        )
        .astype(dtype)
        .T
    )

    targ_comma = (
        xr.DataArray(
            [["a, 11111", "bb, 222", "cccc, 33"], ["a, 4", "bb, 5555", "cccc, 66"]],
            dims=["X", "Y"],
        )
        .astype(dtype)
        .T
    )

    res_blank = values_1.str.cat(values_2)
    res_add = values_1.str + values_2
    res_space = values_1.str.cat(values_2, sep=" ")
    res_bars = values_1.str.cat(values_2, sep="||")
    res_comma = values_1.str.cat(values_2, sep=", ")

    assert res_blank.dtype == targ_blank.dtype
    assert res_add.dtype == targ_blank.dtype
    assert res_space.dtype == targ_space.dtype
    assert res_bars.dtype == targ_bars.dtype
    assert res_comma.dtype == targ_comma.dtype

    assert_equal(res_blank, targ_blank)
    assert_equal(res_add, targ_blank)
    assert_equal(res_space, targ_space)
    assert_equal(res_bars, targ_bars)
    assert_equal(res_comma, targ_comma)


def test_cat_broadcast_both(dtype) -> None:
    values_1 = xr.DataArray(
        ["a", "bb", "cccc"],
        dims=["Y"],
    ).astype(dtype)
    values_2 = xr.DataArray(
        ["11111", "4"],
        dims=["X"],
    )

    targ_blank = (
        xr.DataArray(
            [["a11111", "bb11111", "cccc11111"], ["a4", "bb4", "cccc4"]],
            dims=["X", "Y"],
        )
        .astype(dtype)
        .T
    )

    targ_space = (
        xr.DataArray(
            [["a 11111", "bb 11111", "cccc 11111"], ["a 4", "bb 4", "cccc 4"]],
            dims=["X", "Y"],
        )
        .astype(dtype)
        .T
    )

    targ_bars = (
        xr.DataArray(
            [["a||11111", "bb||11111", "cccc||11111"], ["a||4", "bb||4", "cccc||4"]],
            dims=["X", "Y"],
        )
        .astype(dtype)
        .T
    )

    targ_comma = (
        xr.DataArray(
            [["a, 11111", "bb, 11111", "cccc, 11111"], ["a, 4", "bb, 4", "cccc, 4"]],
            dims=["X", "Y"],
        )
        .astype(dtype)
        .T
    )

    res_blank = values_1.str.cat(values_2)
    res_add = values_1.str + values_2
    res_space = values_1.str.cat(values_2, sep=" ")
    res_bars = values_1.str.cat(values_2, sep="||")
    res_comma = values_1.str.cat(values_2, sep=", ")

    assert res_blank.dtype == targ_blank.dtype
    assert res_add.dtype == targ_blank.dtype
    assert res_space.dtype == targ_space.dtype
    assert res_bars.dtype == targ_bars.dtype
    assert res_comma.dtype == targ_comma.dtype

    assert_equal(res_blank, targ_blank)
    assert_equal(res_add, targ_blank)
    assert_equal(res_space, targ_space)
    assert_equal(res_bars, targ_bars)
    assert_equal(res_comma, targ_comma)


def test_cat_multi() -> None:
    values_1 = xr.DataArray(
        ["11111", "4"],
        dims=["X"],
    )

    values_2 = xr.DataArray(
        ["a", "bb", "cccc"],
        dims=["Y"],
    ).astype(np.bytes_)

    values_3 = np.array(3.4)

    values_4 = ""

    values_5 = np.array("", dtype=np.str_)

    sep = xr.DataArray(
        [" ", ", "],
        dims=["ZZ"],
    ).astype(np.str_)

    expected = xr.DataArray(
        [
            [
                ["11111 a 3.4  ", "11111, a, 3.4, , "],
                ["11111 bb 3.4  ", "11111, bb, 3.4, , "],
                ["11111 cccc 3.4  ", "11111, cccc, 3.4, , "],
            ],
            [
                ["4 a 3.4  ", "4, a, 3.4, , "],
                ["4 bb 3.4  ", "4, bb, 3.4, , "],
                ["4 cccc 3.4  ", "4, cccc, 3.4, , "],
            ],
        ],
        dims=["X", "Y", "ZZ"],
    ).astype(np.str_)

    res = values_1.str.cat(values_2, values_3, values_4, values_5, sep=sep)

    assert res.dtype == expected.dtype
    assert_equal(res, expected)


def test_join_scalar(dtype) -> None:
    values = xr.DataArray("aaa").astype(dtype)

    targ = xr.DataArray("aaa").astype(dtype)

    res_blank = values.str.join()
    res_space = values.str.join(sep=" ")

    assert res_blank.dtype == targ.dtype
    assert res_space.dtype == targ.dtype

    assert_identical(res_blank, targ)
    assert_identical(res_space, targ)


def test_join_vector(dtype) -> None:
    values = xr.DataArray(
        ["a", "bb", "cccc"],
        dims=["Y"],
    ).astype(dtype)

    targ_blank = xr.DataArray("abbcccc").astype(dtype)
    targ_space = xr.DataArray("a bb cccc").astype(dtype)

    res_blank_none = values.str.join()
    res_blank_y = values.str.join(dim="Y")

    res_space_none = values.str.join(sep=" ")
    res_space_y = values.str.join(dim="Y", sep=" ")

    assert res_blank_none.dtype == targ_blank.dtype
    assert res_blank_y.dtype == targ_blank.dtype
    assert res_space_none.dtype == targ_space.dtype
    assert res_space_y.dtype == targ_space.dtype

    assert_identical(res_blank_none, targ_blank)
    assert_identical(res_blank_y, targ_blank)
    assert_identical(res_space_none, targ_space)
    assert_identical(res_space_y, targ_space)


def test_join_2d(dtype) -> None:
    values = xr.DataArray(
        [["a", "bb", "cccc"], ["ddddd", "eeee", "fff"]],
        dims=["X", "Y"],
    ).astype(dtype)

    targ_blank_x = xr.DataArray(
        ["addddd", "bbeeee", "ccccfff"],
        dims=["Y"],
    ).astype(dtype)
    targ_space_x = xr.DataArray(
        ["a ddddd", "bb eeee", "cccc fff"],
        dims=["Y"],
    ).astype(dtype)

    targ_blank_y = xr.DataArray(
        ["abbcccc", "dddddeeeefff"],
        dims=["X"],
    ).astype(dtype)
    targ_space_y = xr.DataArray(
        ["a bb cccc", "ddddd eeee fff"],
        dims=["X"],
    ).astype(dtype)

    res_blank_x = values.str.join(dim="X")
    res_blank_y = values.str.join(dim="Y")

    res_space_x = values.str.join(dim="X", sep=" ")
    res_space_y = values.str.join(dim="Y", sep=" ")

    assert res_blank_x.dtype == targ_blank_x.dtype
    assert res_blank_y.dtype == targ_blank_y.dtype
    assert res_space_x.dtype == targ_space_x.dtype
    assert res_space_y.dtype == targ_space_y.dtype

    assert_identical(res_blank_x, targ_blank_x)
    assert_identical(res_blank_y, targ_blank_y)
    assert_identical(res_space_x, targ_space_x)
    assert_identical(res_space_y, targ_space_y)

    with pytest.raises(
        ValueError, match=r"Dimension must be specified for multidimensional arrays."
    ):
        values.str.join()


def test_join_broadcast(dtype) -> None:
    values = xr.DataArray(
        ["a", "bb", "cccc"],
        dims=["X"],
    ).astype(dtype)

    sep = xr.DataArray(
        [" ", ", "],
        dims=["ZZ"],
    ).astype(dtype)

    expected = xr.DataArray(
        ["a bb cccc", "a, bb, cccc"],
        dims=["ZZ"],
    ).astype(dtype)

    res = values.str.join(sep=sep)

    assert res.dtype == expected.dtype
    assert_identical(res, expected)


def test_format_scalar() -> None:
    values = xr.DataArray(
        ["{}.{Y}.{ZZ}", "{},{},{X},{X}", "{X}-{Y}-{ZZ}"],
        dims=["X"],
    ).astype(np.str_)

    pos0 = 1
    pos1 = 1.2
    pos2 = "2.3"
    X = "'test'"
    Y = "X"
    ZZ = None
    W = "NO!"

    expected = xr.DataArray(
        ["1.X.None", "1,1.2,'test','test'", "'test'-X-None"],
        dims=["X"],
    ).astype(np.str_)

    res = values.str.format(pos0, pos1, pos2, X=X, Y=Y, ZZ=ZZ, W=W)

    assert res.dtype == expected.dtype
    assert_equal(res, expected)


def test_format_broadcast() -> None:
    values = xr.DataArray(
        ["{}.{Y}.{ZZ}", "{},{},{X},{X}", "{X}-{Y}-{ZZ}"],
        dims=["X"],
    ).astype(np.str_)

    pos0 = 1
    pos1 = 1.2

    pos2 = xr.DataArray(
        ["2.3", "3.44444"],
        dims=["YY"],
    )

    X = "'test'"
    Y = "X"
    ZZ = None
    W = "NO!"

    expected = xr.DataArray(
        [
            ["1.X.None", "1.X.None"],
            ["1,1.2,'test','test'", "1,1.2,'test','test'"],
            ["'test'-X-None", "'test'-X-None"],
        ],
        dims=["X", "YY"],
    ).astype(np.str_)

    res = values.str.format(pos0, pos1, pos2, X=X, Y=Y, ZZ=ZZ, W=W)

    assert res.dtype == expected.dtype
    assert_equal(res, expected)


def test_mod_scalar() -> None:
    values = xr.DataArray(
        ["%s.%s.%s", "%s,%s,%s", "%s-%s-%s"],
        dims=["X"],
    ).astype(np.str_)

    pos0 = 1
    pos1 = 1.2
    pos2 = "2.3"

    expected = xr.DataArray(
        ["1.1.2.2.3", "1,1.2,2.3", "1-1.2-2.3"],
        dims=["X"],
    ).astype(np.str_)

    res = values.str % (pos0, pos1, pos2)

    assert res.dtype == expected.dtype
    assert_equal(res, expected)


def test_mod_dict() -> None:
    values = xr.DataArray(
        ["%(a)s.%(a)s.%(b)s", "%(b)s,%(c)s,%(b)s", "%(c)s-%(b)s-%(a)s"],
        dims=["X"],
    ).astype(np.str_)

    a = 1
    b = 1.2
    c = "2.3"

    expected = xr.DataArray(
        ["1.1.1.2", "1.2,2.3,1.2", "2.3-1.2-1"],
        dims=["X"],
    ).astype(np.str_)

    res = values.str % {"a": a, "b": b, "c": c}

    assert res.dtype == expected.dtype
    assert_equal(res, expected)


def test_mod_broadcast_single() -> None:
    values = xr.DataArray(
        ["%s_1", "%s_2", "%s_3"],
        dims=["X"],
    ).astype(np.str_)

    pos = xr.DataArray(
        ["2.3", "3.44444"],
        dims=["YY"],
    )

    expected = xr.DataArray(
        [["2.3_1", "3.44444_1"], ["2.3_2", "3.44444_2"], ["2.3_3", "3.44444_3"]],
        dims=["X", "YY"],
    ).astype(np.str_)

    res = values.str % pos

    assert res.dtype == expected.dtype
    assert_equal(res, expected)


def test_mod_broadcast_multi() -> None:
    values = xr.DataArray(
        ["%s.%s.%s", "%s,%s,%s", "%s-%s-%s"],
        dims=["X"],
    ).astype(np.str_)

    pos0 = 1
    pos1 = 1.2

    pos2 = xr.DataArray(
        ["2.3", "3.44444"],
        dims=["YY"],
    )

    expected = xr.DataArray(
        [
            ["1.1.2.2.3", "1.1.2.3.44444"],
            ["1,1.2,2.3", "1,1.2,3.44444"],
            ["1-1.2-2.3", "1-1.2-3.44444"],
        ],
        dims=["X", "YY"],
    ).astype(np.str_)

    res = values.str % (pos0, pos1, pos2)

    assert res.dtype == expected.dtype
    assert_equal(res, expected)
