import pytest


@pytest.fixture(params=[None, False])
def sort(request):
    """
    Valid values for 'sort' parameter used in a variety Index methods.

    Caution:
        There are also "sort" fixtures with params [True, False] or [None, True, False].
        For instance these might be used for DataFrame.groupby.

        We can't extend the params here either as sort=True is not permitted in
        many Index methods.
    """
    return request.param
