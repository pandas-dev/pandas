import pytest

from pandas.util._validators import validate_fillna_kwargs

def test_no_value_or_method():

    msg = "Must specify a fill 'value' or 'method'."

    with pytest.raises(ValueError, match=msg):
         validate_fillna_kwargs(None, None, None)


@pytest.mark.parametrize("value", [0, {'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4}])
@pytest.mark.parametrize("method", ['backfill', 'bfill', 'pad', 'ffill'])
def test_value_and_method(value, method):

    msg = "Cannot specify both 'value' and 'method'."
    
    with pytest.raises(ValueError, match=msg):
        validate_fillna_kwargs(value, method, None)

@pytest.mark.parametrize("value", [(1, 2, 3), [1, 2, 3]])
def test_valid_value(value):

        msg = ('"value" parameter must be a scalar or dict, but '
               'you passed a "{0}"'.format(type(value).__name__))

        with pytest.raises(TypeError, match=msg):
                validate_fillna_kwargs(value, None, None)

def test_integer_limit():

        msg = "Limit must be an integer"
        with pytest.raises(ValueError, match=msg):
                validate_fillna_kwargs(0,None, 0.5)
def test_positive_limit():

        msg = "Limit must be greater than 0"
        with pytest.raises(ValueError, match=msg):
                validate_fillna_kwargs(5, None, -5)


