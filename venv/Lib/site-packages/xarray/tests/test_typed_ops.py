import numpy as np

from xarray import DataArray, Dataset, Variable


def test_variable_typed_ops() -> None:
    """Tests for type checking of typed_ops on Variable"""

    var = Variable(dims=["t"], data=[1, 2, 3])

    def _test(var: Variable) -> None:
        # mypy checks the input type
        assert isinstance(var, Variable)

    _int: int = 1
    _list = [1, 2, 3]
    _ndarray = np.array([1, 2, 3])

    # __add__ as an example of binary ops
    _test(var + _int)
    _test(var + _list)
    _test(var + _ndarray)
    _test(var + var)

    # __radd__ as an example of reflexive binary ops
    _test(_int + var)
    _test(_list + var)
    _test(_ndarray + var)  # type: ignore[arg-type]  # numpy problem

    # __eq__ as an example of cmp ops
    _test(var == _int)
    _test(var == _list)
    _test(var == _ndarray)
    _test(_int == var)  # type: ignore[arg-type]  # typeshed problem
    _test(_list == var)  # type: ignore[arg-type]  # typeshed problem
    _test(_ndarray == var)

    # __lt__ as another example of cmp ops
    _test(var < _int)
    _test(var < _list)
    _test(var < _ndarray)
    _test(_int > var)
    _test(_list > var)
    _test(_ndarray > var)  # type: ignore[arg-type]  # numpy problem

    # __iadd__ as an example of inplace binary ops
    var += _int
    var += _list
    var += _ndarray

    # __neg__ as an example of unary ops
    _test(-var)


def test_dataarray_typed_ops() -> None:
    """Tests for type checking of typed_ops on DataArray"""

    da = DataArray([1, 2, 3], dims=["t"])

    def _test(da: DataArray) -> None:
        # mypy checks the input type
        assert isinstance(da, DataArray)

    _int: int = 1
    _list = [1, 2, 3]
    _ndarray = np.array([1, 2, 3])
    _var = Variable(dims=["t"], data=[1, 2, 3])

    # __add__ as an example of binary ops
    _test(da + _int)
    _test(da + _list)
    _test(da + _ndarray)
    _test(da + _var)
    _test(da + da)

    # __radd__ as an example of reflexive binary ops
    _test(_int + da)
    _test(_list + da)
    _test(_ndarray + da)  # type: ignore[arg-type]  # numpy problem
    _test(_var + da)

    # __eq__ as an example of cmp ops
    _test(da == _int)
    _test(da == _list)
    _test(da == _ndarray)
    _test(da == _var)
    _test(_int == da)  # type: ignore[arg-type]  # typeshed problem
    _test(_list == da)  # type: ignore[arg-type]  # typeshed problem
    _test(_ndarray == da)
    _test(_var == da)

    # __lt__ as another example of cmp ops
    _test(da < _int)
    _test(da < _list)
    _test(da < _ndarray)
    _test(da < _var)
    _test(_int > da)
    _test(_list > da)
    _test(_ndarray > da)  # type: ignore[arg-type]  # numpy problem
    _test(_var > da)

    # __iadd__ as an example of inplace binary ops
    da += _int
    da += _list
    da += _ndarray
    da += _var

    # __neg__ as an example of unary ops
    _test(-da)


def test_dataset_typed_ops() -> None:
    """Tests for type checking of typed_ops on Dataset"""

    ds = Dataset({"a": ("t", [1, 2, 3])})

    def _test(ds: Dataset) -> None:
        # mypy checks the input type
        assert isinstance(ds, Dataset)

    _int: int = 1
    _list = [1, 2, 3]
    _ndarray = np.array([1, 2, 3])
    _var = Variable(dims=["t"], data=[1, 2, 3])
    _da = DataArray([1, 2, 3], dims=["t"])

    # __add__ as an example of binary ops
    _test(ds + _int)
    _test(ds + _list)
    _test(ds + _ndarray)
    _test(ds + _var)
    _test(ds + _da)
    _test(ds + ds)

    # __radd__ as an example of reflexive binary ops
    _test(_int + ds)
    _test(_list + ds)
    _test(_ndarray + ds)
    _test(_var + ds)
    _test(_da + ds)

    # __eq__ as an example of cmp ops
    _test(ds == _int)
    _test(ds == _list)
    _test(ds == _ndarray)
    _test(ds == _var)
    _test(ds == _da)
    _test(_int == ds)  # type: ignore[arg-type]  # typeshed problem
    _test(_list == ds)  # type: ignore[arg-type]  # typeshed problem
    _test(_ndarray == ds)
    _test(_var == ds)
    _test(_da == ds)

    # __lt__ as another example of cmp ops
    _test(ds < _int)
    _test(ds < _list)
    _test(ds < _ndarray)
    _test(ds < _var)
    _test(ds < _da)
    _test(_int > ds)
    _test(_list > ds)
    _test(_ndarray > ds)  # type: ignore[arg-type]  # numpy problem
    _test(_var > ds)
    _test(_da > ds)

    # __iadd__ as an example of inplace binary ops
    ds += _int
    ds += _list
    ds += _ndarray
    ds += _var
    ds += _da

    # __neg__ as an example of unary ops
    _test(-ds)


def test_dataarray_groupy_typed_ops() -> None:
    """Tests for type checking of typed_ops on DataArrayGroupBy"""

    da = DataArray([1, 2, 3], coords={"x": ("t", [1, 2, 2])}, dims=["t"])
    grp = da.groupby("x")

    def _testda(da: DataArray) -> None:
        # mypy checks the input type
        assert isinstance(da, DataArray)

    def _testds(ds: Dataset) -> None:
        # mypy checks the input type
        assert isinstance(ds, Dataset)

    _da = DataArray([5, 6], coords={"x": [1, 2]}, dims="x")
    _ds = _da.to_dataset(name="a")

    # __add__ as an example of binary ops
    _testda(grp + _da)
    _testds(grp + _ds)

    # __radd__ as an example of reflexive binary ops
    _testda(_da + grp)
    _testds(_ds + grp)

    # __eq__ as an example of cmp ops
    _testda(grp == _da)
    _testda(_da == grp)
    _testds(grp == _ds)
    _testds(_ds == grp)

    # __lt__ as another example of cmp ops
    _testda(grp < _da)
    _testda(_da > grp)
    _testds(grp < _ds)
    _testds(_ds > grp)


def test_dataset_groupy_typed_ops() -> None:
    """Tests for type checking of typed_ops on DatasetGroupBy"""

    ds = Dataset({"a": ("t", [1, 2, 3])}, coords={"x": ("t", [1, 2, 2])})
    grp = ds.groupby("x")

    def _test(ds: Dataset) -> None:
        # mypy checks the input type
        assert isinstance(ds, Dataset)

    _da = DataArray([5, 6], coords={"x": [1, 2]}, dims="x")
    _ds = _da.to_dataset(name="a")

    # __add__ as an example of binary ops
    _test(grp + _da)
    _test(grp + _ds)

    # __radd__ as an example of reflexive binary ops
    _test(_da + grp)
    _test(_ds + grp)

    # __eq__ as an example of cmp ops
    _test(grp == _da)
    _test(_da == grp)
    _test(grp == _ds)
    _test(_ds == grp)

    # __lt__ as another example of cmp ops
    _test(grp < _da)
    _test(_da > grp)
    _test(grp < _ds)
    _test(_ds > grp)
