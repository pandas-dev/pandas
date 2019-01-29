import pytest

# A set of Base EA tests that are know to not work for
# the object-dtype PandasArray holding nested data.
skips = {
    'BaseCastingTests.test_astype_str',
    'BaseConstructorsTests.test_array_from_scalars',
    # tuple isn't instance of np.object
    'BaseGetitemTests.test_getitem_scalar',
    # Can't pass tuples to _from_sequence
    'BaseGetitemTests.test_take_series',
    # np.array shape inference
    'BaseInterfaceTests.test_array_interface',
    # Can't construct expected.
    'BaseMethodsTests.test_unique',
    'BaseMethodsTests.test_combine_add',
    'BaseMethodsTests.test_shift_fill_value',
    'BaseMethodsTests.test_where_series',
    'BaseMethodsTests.test_repeat',
    # Can't hash ndarray[tuple]
    'BaseMethodsTests.test_hash_pandas_object_works',
    # Can't construct expected.
    'BaseReshapingTests.test_merge',
    'BaseReshapingTests.test_merge_on_extension_array',
    'BaseReshapingTests.test_merge_on_extension_array_duplicates',

    # ndarray setting
    'BaseSetitemTests.test_setitem_scalar_series',
    'BaseSetitemTests.test_setitem_sequence',
    'BaseSetitemTests.test_setitem_sequence_mismatched_length_raises',
    'BaseSetitemTests.test_setitem_sequence_broadcasts',
    'BaseSetitemTests.test_setitem_sequence_broadcasts',
    'BaseSetitemTests.test_setitem_loc_scalar_mixed',
    'BaseSetitemTests.test_setitem_iloc_scalar_mixed',
    'BaseSetitemTests.test_setitem_loc_scalar_multiple_homogoneous',
    'BaseSetitemTests.test_setitem_iloc_scalar_multiple_homogoneous',
    'BaseSetitemTests.test_setitem_mask_broadcast',
    'BaseSetitemTests.test_setitem_scalar_key_sequence_raise',

    'BaseParsingTests.test_EA_types',
}


def pytest_collection_modifyitems(config, items):
    skip = pytest.mark.skip(reason="Skipping this because ...")
    for item in items:
        # TODO: See if pytest has a better way to resolve the *value*
        # supplied to a fixture. Right now .keywords gets things
        # like 'object' or 'data-object'.
        parts = item.name.split("[")
        if (len(parts) > 1 and 'object' in item.name.split('[')[1]
                and item.obj.__qualname__ in skips):
            item.add_marker(skip)
