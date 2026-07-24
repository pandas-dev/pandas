from numpy.f2py import capi_maps


def test_complex_long_double_capi_map():
    assert capi_maps.c2capi_map["complex_long_double"] == "NPY_CLONGDOUBLE"


def test_complex_long_double_is_distinct():
    assert capi_maps.c2pycode_map["complex_long_double"] != capi_maps.c2pycode_map["complex_double"]
    assert capi_maps.c2capi_map["complex_long_double"] != capi_maps.c2capi_map["complex_double"]
