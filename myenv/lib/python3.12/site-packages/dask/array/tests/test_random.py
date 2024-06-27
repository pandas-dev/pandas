from __future__ import annotations

import pytest

pytest.importorskip("numpy")

import numpy as np

import dask
import dask.array as da
from dask.array.core import Array
from dask.array.utils import assert_eq
from dask.multiprocessing import _dumps, _loads
from dask.utils import key_split


@pytest.fixture(params=[da.random.RandomState, da.random.default_rng])
def generator_class(request):
    return request.param


def test_generators(generator_class):
    state = generator_class(5)
    x = state.normal(10, 1, size=10, chunks=5)
    assert_eq(x, x)

    state = generator_class(5)
    y = state.normal(10, 1, size=10, chunks=5)
    assert_eq(x, y)


@pytest.mark.parametrize(
    "sd", [None, 42, np.random.PCG64, da.random.Generator(np.random.PCG64)], ids=type
)
def test_default_rng(sd):
    rng = da.random.default_rng(seed=sd)
    assert isinstance(rng, da.random.Generator)


def test_concurrency(generator_class):
    state = generator_class(5)
    x = state.normal(10, 1, size=10, chunks=2)

    state = generator_class(5)
    y = state.normal(10, 1, size=10, chunks=2)
    assert (x.compute(scheduler="processes") == y.compute(scheduler="processes")).all()


def test_doc_randomstate(generator_class):
    assert "mean" in generator_class(5).normal.__doc__


def test_doc_generator():
    x = da.random.default_rng()
    assert "Generator" in str(x)
    assert "Generator" in repr(x)


def test_serializability(generator_class):
    state = generator_class(5)
    x = state.normal(10, 1, size=10, chunks=5)

    y = _loads(_dumps(x))

    assert_eq(x, y)


def test_determinisim_through_dask_values(generator_class):
    samples_1 = generator_class(42).normal(size=1000, chunks=10)
    samples_2 = generator_class(42).normal(size=1000, chunks=10)

    assert set(samples_1.dask) == set(samples_2.dask)
    assert_eq(samples_1, samples_2)


def test_generator_consistent_names(generator_class):
    state1 = generator_class(42)
    state2 = generator_class(42)
    assert sorted(state1.normal(size=(100, 100), chunks=(10, 10)).dask) == sorted(
        state2.normal(size=(100, 100), chunks=(10, 10)).dask
    )
    assert sorted(
        state1.normal(size=100, loc=4.5, scale=5.0, chunks=10).dask
    ) == sorted(state2.normal(size=100, loc=4.5, scale=5.0, chunks=10).dask)


def test_random(generator_class):
    a = generator_class().random((10, 10), chunks=(5, 5))
    assert isinstance(a, Array)
    assert isinstance(a.name, str) and a.name
    assert a.shape == (10, 10)
    assert a.chunks == ((5, 5), (5, 5))

    x = set(np.array(a).flat)

    assert len(x) > 90


def test_parametrized_random_function(generator_class):
    a = generator_class().exponential(1000, (10, 10), chunks=(5, 5))
    assert isinstance(a, Array)
    assert isinstance(a.name, str) and a.name
    assert a.shape == (10, 10)
    assert a.chunks == ((5, 5), (5, 5))

    x = np.array(a)
    assert 10 < x.mean() < 100000

    y = set(x.flat)
    assert len(y) > 90


def test_kwargs(generator_class):
    a = generator_class().normal(loc=10.0, scale=0.1, size=(10, 10), chunks=(5, 5))
    assert isinstance(a, Array)
    x = np.array(a)
    assert 8 < x.mean() < 12


def test_unique_names(generator_class):
    a = generator_class().random((10, 10), chunks=(5, 5))
    b = generator_class().random((10, 10), chunks=(5, 5))

    assert a.name != b.name


def test_docs(generator_class):
    assert "exponential" in generator_class().exponential.__doc__
    assert "exponential" in generator_class().exponential.__name__
    assert "# doctest: +SKIP" in generator_class().normal.__doc__


def test_can_make_really_big_random_array(generator_class):
    generator_class().normal(10, 1, (1000000, 1000000), chunks=(100000, 100000))


def test_random_seed():
    da.random.seed(123)
    x = da.random.normal(size=10, chunks=5)
    y = da.random.normal(size=10, chunks=5)

    da.random.seed(123)
    a = da.random.normal(size=10, chunks=5)
    b = da.random.normal(size=10, chunks=5)

    assert_eq(x, a)
    assert_eq(y, b)


def test_consistent_across_sizes(generator_class):
    x1 = generator_class(123).random(20, chunks=20)
    x2 = generator_class(123).random(100, chunks=20)[:20]
    x3 = generator_class(123).random(200, chunks=20)[:20]
    assert_eq(x1, x2)
    assert_eq(x1, x3)


@pytest.mark.parametrize("sz", [None, 5, (2, 2)], ids=type)
def test_random_all(sz):
    da.random.beta(1, 2, size=sz, chunks=3).compute()
    da.random.binomial(10, 0.5, size=sz, chunks=3).compute()
    da.random.chisquare(1, size=sz, chunks=3).compute()
    da.random.exponential(1, size=sz, chunks=3).compute()
    da.random.f(1, 2, size=sz, chunks=3).compute()
    da.random.gamma(5, 1, size=sz, chunks=3).compute()
    da.random.geometric(1, size=sz, chunks=3).compute()
    da.random.gumbel(1, size=sz, chunks=3).compute()
    da.random.hypergeometric(1, 2, 3, size=sz, chunks=3).compute()
    da.random.laplace(size=sz, chunks=3).compute()
    da.random.logistic(size=sz, chunks=3).compute()
    da.random.lognormal(size=sz, chunks=3).compute()
    da.random.logseries(0.5, size=sz, chunks=3).compute()
    da.random.multinomial(20, [1 / 6.0] * 6, size=sz, chunks=3).compute()
    da.random.negative_binomial(5, 0.5, size=sz, chunks=3).compute()
    da.random.noncentral_chisquare(2, 2, size=sz, chunks=3).compute()

    da.random.noncentral_f(2, 2, 3, size=sz, chunks=3).compute()
    da.random.normal(2, 2, size=sz, chunks=3).compute()
    da.random.pareto(1, size=sz, chunks=3).compute()
    da.random.poisson(size=sz, chunks=3).compute()

    da.random.power(1, size=sz, chunks=3).compute()
    da.random.rayleigh(size=sz, chunks=3).compute()

    da.random.triangular(1, 2, 3, size=sz, chunks=3).compute()
    da.random.uniform(size=sz, chunks=3).compute()
    da.random.vonmises(2, 3, size=sz, chunks=3).compute()
    da.random.wald(1, 2, size=sz, chunks=3).compute()

    da.random.weibull(2, size=sz, chunks=3).compute()
    da.random.zipf(2, size=sz, chunks=3).compute()

    da.random.standard_cauchy(size=sz, chunks=3).compute()
    da.random.standard_exponential(size=sz, chunks=3).compute()
    da.random.standard_gamma(2, size=sz, chunks=3).compute()
    da.random.standard_normal(size=sz, chunks=3).compute()
    da.random.standard_t(2, size=sz, chunks=3).compute()


def test_RandomState_only_funcs():
    da.random.randint(10, size=5, chunks=3).compute()
    with pytest.warns(DeprecationWarning):
        da.random.random_integers(10, size=5, chunks=3).compute()
    da.random.random_sample(10, chunks=3).compute()


@pytest.mark.parametrize("sz", [None, 5, (2, 2)], ids=type)
def test_Generator_only_funcs(sz):
    da.random.default_rng().integers(5, high=15, size=sz, chunks=3).compute()
    da.random.default_rng().multivariate_hypergeometric(
        [16, 8, 4], 6, size=sz, chunks=6
    ).compute()


@pytest.mark.parametrize("sz", [None, 5, (2, 2)], ids=type)
def test_random_all_with_class_methods(generator_class, sz):
    def rnd_test(func, *args, **kwargs):
        a = func(*args, **kwargs)
        assert type(a._meta) == np.ndarray
        assert_eq(a, a)  # Check that _meta and computed arrays match types

    rnd_test(generator_class().beta, 1, 2, size=sz, chunks=3)
    rnd_test(generator_class().binomial, 10, 0.5, size=sz, chunks=3)
    rnd_test(generator_class().chisquare, 1, size=sz, chunks=3)
    rnd_test(generator_class().exponential, 1, size=sz, chunks=3)
    rnd_test(generator_class().f, 1, 2, size=sz, chunks=3)
    rnd_test(generator_class().gamma, 5, 1, size=sz, chunks=3)
    rnd_test(generator_class().geometric, 1, size=sz, chunks=3)
    rnd_test(generator_class().gumbel, 1, size=sz, chunks=3)
    rnd_test(generator_class().hypergeometric, 1, 2, 3, size=sz, chunks=3)
    rnd_test(generator_class().laplace, size=sz, chunks=3)
    rnd_test(generator_class().logistic, size=sz, chunks=3)
    rnd_test(generator_class().lognormal, size=sz, chunks=3)
    rnd_test(generator_class().logseries, 0.5, size=sz, chunks=3)
    rnd_test(generator_class().multinomial, 20, [1 / 6.0] * 6, size=sz, chunks=3)
    rnd_test(generator_class().negative_binomial, 5, 0.5, size=sz, chunks=3)
    rnd_test(generator_class().noncentral_chisquare, 2, 2, size=sz, chunks=3)
    rnd_test(generator_class().noncentral_f, 2, 2, 3, size=sz, chunks=3)
    rnd_test(generator_class().normal, 2, 2, size=sz, chunks=3)
    rnd_test(generator_class().pareto, 1, size=sz, chunks=3)
    rnd_test(generator_class().poisson, size=sz, chunks=3)
    rnd_test(generator_class().power, 1, size=sz, chunks=3)
    rnd_test(generator_class().rayleigh, size=sz, chunks=3)
    rnd_test(generator_class().triangular, 1, 2, 3, size=sz, chunks=3)
    rnd_test(generator_class().uniform, size=sz, chunks=3)
    rnd_test(generator_class().vonmises, 2, 3, size=sz, chunks=3)
    rnd_test(generator_class().wald, 1, 2, size=sz, chunks=3)
    rnd_test(generator_class().weibull, 2, size=sz, chunks=3)
    rnd_test(generator_class().zipf, 2, size=sz, chunks=3)
    rnd_test(generator_class().standard_cauchy, size=sz, chunks=3)
    rnd_test(generator_class().standard_exponential, size=sz, chunks=3)
    rnd_test(generator_class().standard_gamma, 2, size=sz, chunks=3)
    rnd_test(generator_class().standard_normal, size=sz, chunks=3)
    rnd_test(generator_class().standard_t, 2, size=sz, chunks=3)


@pytest.mark.skipif(
    not hasattr(np, "broadcast_to"), reason='requires numpy 1.10 method "broadcast_to"'
)
def test_array_broadcasting(generator_class):
    arr = np.arange(6).reshape((2, 3))
    daones = da.ones((2, 3, 4), chunks=3)
    assert generator_class().poisson(arr, chunks=3).compute().shape == (2, 3)

    for x in (arr, daones):
        y = generator_class().normal(x, 2, chunks=3)
        assert y.shape == x.shape
        assert y.compute().shape == x.shape

    y = generator_class().normal(daones, 2, chunks=3)
    assert set(daones.dask).issubset(set(y.dask))

    assert generator_class().normal(
        np.ones((1, 4)), da.ones((2, 3, 4), chunks=(2, 3, 4)), chunks=(2, 3, 4)
    ).compute().shape == (2, 3, 4)
    assert generator_class().normal(
        scale=np.ones((1, 4)),
        loc=da.ones((2, 3, 4), chunks=(2, 3, 4)),
        size=(2, 2, 3, 4),
        chunks=(2, 2, 3, 4),
    ).compute().shape == (2, 2, 3, 4)

    with pytest.raises(ValueError):
        generator_class().normal(arr, np.ones((3, 1)), size=(2, 3, 4), chunks=3)

    for o in (np.ones(100), da.ones(100, chunks=(50,)), 1):
        a = generator_class().normal(1000 * o, 0.01, chunks=(50,))
        assert 800 < a.mean().compute() < 1200

    # ensure that mis-matched chunks align well
    x = np.arange(10) ** 3
    y = da.from_array(x, chunks=(1,))
    z = generator_class().normal(y, 0.01, chunks=(10,))

    assert 0.8 < z.mean().compute() / x.mean() < 1.2


def test_multinomial(generator_class):
    for size, chunks in [(5, 3), ((5, 4), (2, 3))]:
        x = generator_class().multinomial(20, [1 / 6.0] * 6, size=size, chunks=chunks)
        y = np.random.multinomial(20, [1 / 6.0] * 6, size=size)

        assert x.shape == y.shape == x.compute().shape


def test_choice(generator_class):
    np_generator = {
        da.random.RandomState: np.random.RandomState,
        da.random.default_rng: np.random.default_rng,
    }
    np_dtype = np_generator[generator_class]().choice(1, size=()).dtype
    size = (10, 3)
    chunks = 4
    x = generator_class().choice(3, size=size, chunks=chunks)
    assert x.dtype == np_dtype
    assert x.shape == size
    res = x.compute()
    assert res.dtype == np_dtype
    assert res.shape == size

    py_a = [1, 3, 5, 7, 9]
    np_a = np.array(py_a, dtype="f8")
    da_a = da.from_array(np_a, chunks=2)

    for a in [py_a, np_a, da_a]:
        x = generator_class().choice(a, size=size, chunks=chunks)
        res = x.compute()
        expected_dtype = np.asarray(a).dtype
        assert x.dtype == expected_dtype
        assert res.dtype == expected_dtype
        assert set(np.unique(res)).issubset(np_a)

    np_p = np.array([0, 0.2, 0.2, 0.3, 0.3])
    da_p = da.from_array(np_p, chunks=2)

    for a, p in [(da_a, np_p), (np_a, da_p)]:
        x = generator_class().choice(a, size=size, chunks=chunks, p=p)
        res = x.compute()
        assert x.dtype == np_a.dtype
        assert res.dtype == np_a.dtype
        assert set(np.unique(res)).issubset(np_a[1:])

    np_dtype = np_generator[generator_class]().choice(1, size=(), p=np.array([1])).dtype
    x = generator_class().choice(5, size=size, chunks=chunks, p=np_p)
    res = x.compute()
    assert x.dtype == np_dtype
    assert res.dtype == np_dtype

    errs = [
        (-1, None),  # negative a
        (np_a[:, None], None),  # a must be 1D
        (np_a, np_p[:, None]),  # p must be 1D
        (np_a, np_p[:-2]),  # a and p must match
        (3, np_p),  # a and p must match
        (4, [0.2, 0.2, 0.3]),
    ]  # p must sum to 1

    for a, p in errs:
        with pytest.raises(ValueError):
            generator_class().choice(a, size=size, chunks=chunks, p=p)

    with pytest.raises(NotImplementedError):
        generator_class().choice(da_a, size=size, chunks=chunks, replace=False)

    # axis needs to be row based since a is 1D
    with pytest.raises(ValueError):
        da.random.default_rng().choice(da_a, size=size, chunks=chunks, axis=1)

    # Want to make sure replace=False works for a single-partition output array
    x = generator_class().choice(da_a, size=da_a.shape[0], chunks=-1, replace=False)
    res = x.compute()
    assert len(res) == len(np.unique(res))


def test_create_with_auto_dimensions():
    with dask.config.set({"array.chunk-size": "128MiB"}):
        x = da.random.random((10000, 10000), chunks=(-1, "auto"))
        assert x.chunks == ((10000,), (1677,) * 5 + (1615,))

        y = da.random.random((10000, 10000), chunks="auto")
        assert y.chunks == ((4096, 4096, 1808),) * 2


def test_names():
    name = da.random.normal(0, 1, size=(1000,), chunks=(500,)).name

    assert name.startswith("normal")
    assert len(key_split(name)) < 10


def test_permutation(generator_class):
    x = da.arange(12, chunks=3)
    y = generator_class().permutation(x)

    assert y.shape == x.shape
    assert y.dtype == x.dtype

    y.compute()  # smoke test

    a = generator_class(0)
    b = generator_class(0)
    r1 = a.permutation(x)
    r2 = b.permutation(x)
    assert_eq(r1, r2)

    x = generator_class().permutation(100)
    assert x.shape == (100,)


def test_auto_chunks(generator_class):
    with dask.config.set({"array.chunk-size": "50 MiB"}):
        x = generator_class().random((10000, 10000))
        assert 4 < x.npartitions < 32


def test_randint_dtype():
    x = da.random.randint(0, 255, size=10, dtype="uint8")
    assert_eq(x, x)
    assert x.dtype == "uint8"
    assert x.compute().dtype == "uint8"


def test_raises_bad_kwarg(generator_class):
    with pytest.raises(Exception) as info:
        generator_class().standard_cauchy(size=(10,), dtype="float64")

    assert "dtype" in str(info.value)


def test_randomstate_kwargs():
    cupy = pytest.importorskip("cupy")

    rs = da.random.RandomState(RandomState=cupy.random.RandomState)
    x = rs.standard_normal((10, 5), dtype=np.float32)
    assert x.dtype == np.float32

    rs = da.random.default_rng(cupy.random.default_rng())
    x = rs.standard_normal((10, 5), dtype=np.float32)
    assert x.dtype == np.float32
