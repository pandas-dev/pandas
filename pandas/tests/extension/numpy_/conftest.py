import pytest

# A set of Base EA tests that are know to not work for
# the object-dtype PandasArray holding nested data.
skips = {
    'TestCasting.test_astype_str',
    'TestConstructors.test_array_from_scalars',
    # tuple isn't instance of np.object
    'TestGetitem.test_getitem_scalar',
    # Can't pass tuples to _from_sequence
    'TestGetitem.test_take_series',
    # np.array shape inference
    'TestInterface.test_array_interface',
    # Can't construct expected.
    'TestMethods.test_unique',
    'TestMethods.test_combine_add',
    'TestMethods.test_shift_fill_value',
    'TestMethods.test_where_series',
    'TestMethods.test_repeat',
    # Can't hash ndarray[tuple]
    'TestMethods.test_hash_pandas_object_works',
    # Can't construct expected.
    'TestReshaping.test_merge',
    'TestReshaping.test_merge_on_extension_array',
    'TestReshaping.test_merge_on_extension_array_duplicates',

    # ndarray setting
    'TestSetitem.test_setitem_scalar_series',
    'TestSetitem.test_setitem_sequence',
    'TestSetitem.test_setitem_sequence_mismatched_length_raises',
    'TestSetitem.test_setitem_sequence_broadcasts',
    'TestSetitem.test_setitem_sequence_broadcasts',
    'TestSetitem.test_setitem_loc_scalar_mixed',
    'TestSetitem.test_setitem_iloc_scalar_mixed',
    'TestSetitem.test_setitem_loc_scalar_multiple_homogoneous',
    'TestSetitem.test_setitem_iloc_scalar_multiple_homogoneous',
    'TestSetitem.test_setitem_mask_broadcast',
    'TestSetitem.test_setitem_scalar_key_sequence_raise',

    # parsing differs.
    'TestParsing.test_EA_types',
}


def pytest_collection_modifyitems(config, items):
    skip = pytest.mark.skip(reason="Skipping for nested data.")
    for item in items:
        # TODO: See if pytest has a better way to resolve the *value*
        # supplied to a fixture. Right now .keywords gets things
        # like 'object' or 'data-object'.
        parts = item.name.split("[")
        qualname = item.parent.obj.__class__.__name__ + '.' + item.obj.__name__
        if (len(parts) > 1 and 'object' in item.name.split('[')[1]
                and qualname in skips):
            item.add_marker(skip)
