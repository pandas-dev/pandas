from __future__ import annotations

import pytest

pytestmark = pytest.mark.gpu

import numpy as np

import dask.array as da
from dask import config
from dask.array.utils import assert_eq

cupy = pytest.importorskip("cupy")


@pytest.mark.parametrize("backend", ["cupy", "numpy"])
@pytest.mark.parametrize("rs", [None, cupy.random.RandomState, np.random.RandomState])
def test_random_all_RandomState(backend, rs):
    # RandomState argument takes priority over backend
    if rs == cupy.random.RandomState:
        expect = cupy.ndarray
    elif rs == np.random.RandomState:
        expect = np.ndarray
    elif backend == "cupy":
        expect = cupy.ndarray
    else:
        expect = np.ndarray

    def rnd_test(func, *args, **kwargs):
        a = func(*args, **kwargs)
        assert type(a._meta) == expect
        assert_eq(a, a)  # Check that _meta and computed arrays match types

    with config.set({"array.backend": backend}):
        rs = da.random.RandomState(RandomState=rs)

        rnd_test(rs.beta, 1, 2, size=5, chunks=3)
        rnd_test(rs.binomial, 10, 0.5, size=5, chunks=3)
        rnd_test(rs.chisquare, 1, size=5, chunks=3)
        rnd_test(rs.exponential, 1, size=5, chunks=3)
        rnd_test(rs.f, 1, 2, size=5, chunks=3)
        rnd_test(rs.gamma, 5, 1, size=5, chunks=3)
        rnd_test(rs.geometric, 1, size=5, chunks=3)
        rnd_test(rs.gumbel, 1, size=5, chunks=3)
        rnd_test(rs.hypergeometric, 1, 2, 3, size=5, chunks=3)
        rnd_test(rs.laplace, size=5, chunks=3)
        rnd_test(rs.logistic, size=5, chunks=3)
        rnd_test(rs.lognormal, size=5, chunks=3)
        rnd_test(rs.logseries, 0.5, size=5, chunks=3)
        # No RandomState for multinomial in CuPy
        # rnd_test(rs.multinomial, 20, [1 / 6.] * 6, size=5, chunks=3)
        rnd_test(rs.negative_binomial, 5, 0.5, size=5, chunks=3)
        rnd_test(rs.noncentral_chisquare, 2, 2, size=5, chunks=3)

        rnd_test(rs.noncentral_f, 2, 2, 3, size=5, chunks=3)
        rnd_test(rs.normal, 2, 2, size=5, chunks=3)
        rnd_test(rs.pareto, 1, size=5, chunks=3)
        rnd_test(rs.poisson, size=5, chunks=3)

        rnd_test(rs.power, 1, size=5, chunks=3)
        rnd_test(rs.rayleigh, size=5, chunks=3)
        rnd_test(rs.random_sample, size=5, chunks=3)

        rnd_test(rs.triangular, 1, 2, 3, size=5, chunks=3)
        rnd_test(rs.uniform, size=5, chunks=3)
        rnd_test(rs.vonmises, 2, 3, size=5, chunks=3)
        rnd_test(rs.wald, 1, 2, size=5, chunks=3)

        rnd_test(rs.weibull, 2, size=5, chunks=3)
        rnd_test(rs.zipf, 2, size=5, chunks=3)

        rnd_test(rs.standard_cauchy, size=5, chunks=3)
        rnd_test(rs.standard_exponential, size=5, chunks=3)
        rnd_test(rs.standard_gamma, 2, size=5, chunks=3)
        rnd_test(rs.standard_normal, size=5, chunks=3)
        rnd_test(rs.standard_t, 2, size=5, chunks=3)


@pytest.mark.parametrize("backend", ["cupy", "numpy"])
def test_random_all_direct_calls(backend):
    if backend == "cupy":
        expect = cupy.ndarray
    else:
        expect = np.ndarray

    def rnd_test(func, *args, **kwargs):
        a = func(*args, **kwargs)
        assert type(a._meta) == expect
        assert_eq(a, a)  # Check that _meta and computed arrays match types

    with config.set({"array.backend": backend}):
        # check calling random functions directly
        rnd_test(da.random.beta, 1, 2, size=5, chunks=3)
        rnd_test(da.random.binomial, 10, 0.5, size=5, chunks=3)
        rnd_test(da.random.chisquare, 1, size=5, chunks=3)
        rnd_test(da.random.exponential, 1, size=5, chunks=3)
        rnd_test(da.random.f, 1, 2, size=5, chunks=3)
        rnd_test(da.random.gamma, 5, 1, size=5, chunks=3)
        rnd_test(da.random.geometric, 1, size=5, chunks=3)
        rnd_test(da.random.gumbel, 1, size=5, chunks=3)
        rnd_test(da.random.hypergeometric, 1, 2, 3, size=5, chunks=3)
        rnd_test(da.random.laplace, size=5, chunks=3)
        rnd_test(da.random.logistic, size=5, chunks=3)
        rnd_test(da.random.lognormal, size=5, chunks=3)
        rnd_test(da.random.logseries, 0.5, size=5, chunks=3)
        # No RandomState for multinomial in CuPy
        if backend != "cupy":
            rnd_test(da.random.multinomial, 20, [1 / 6.0] * 6, size=5, chunks=3)
        rnd_test(da.random.negative_binomial, 5, 0.5, size=5, chunks=3)
        rnd_test(da.random.noncentral_chisquare, 2, 2, size=5, chunks=3)
        rnd_test(da.random.noncentral_f, 2, 2, 3, size=5, chunks=3)
        rnd_test(da.random.normal, 2, 2, size=5, chunks=3)
        rnd_test(da.random.pareto, 1, size=5, chunks=3)
        rnd_test(da.random.poisson, size=5, chunks=3)
        rnd_test(da.random.power, 1, size=5, chunks=3)
        rnd_test(da.random.rayleigh, size=5, chunks=3)
        rnd_test(da.random.randint, low=10, size=5, chunks=3)
        rnd_test(da.random.random, size=5, chunks=3)
        rnd_test(da.random.random_sample, size=5, chunks=3)
        rnd_test(da.random.triangular, 1, 2, 3, size=5, chunks=3)
        rnd_test(da.random.uniform, size=5, chunks=3)
        rnd_test(da.random.vonmises, 2, 3, size=5, chunks=3)
        rnd_test(da.random.wald, 1, 2, size=5, chunks=3)
        rnd_test(da.random.weibull, 2, size=5, chunks=3)
        rnd_test(da.random.zipf, 2, size=5, chunks=3)
        rnd_test(da.random.standard_cauchy, size=5, chunks=3)
        rnd_test(da.random.standard_exponential, size=5, chunks=3)
        rnd_test(da.random.standard_gamma, 2, size=5, chunks=3)
        rnd_test(da.random.standard_normal, size=5, chunks=3)
        rnd_test(da.random.standard_t, 2, size=5, chunks=3)


@pytest.mark.parametrize("backend", ["cupy", "numpy"])
@pytest.mark.parametrize("gen", [None, cupy.random.default_rng, np.random.default_rng])
@pytest.mark.parametrize("shape", [(2), (2, 3), (2, 3, 4), (2, 3, 4, 2)], ids=type)
def test_random_all_Generator(backend, gen, shape):
    # Generator argument takes priority over backend
    if gen == cupy.random.default_rng:
        expect = cupy.ndarray
    elif gen == np.random.default_rng:
        expect = np.ndarray
    elif backend == "cupy":
        expect = cupy.ndarray
    else:
        expect = np.ndarray

    def rnd_test(func, *args, **kwargs):
        a = func(*args, **kwargs)
        assert type(a._meta) == expect
        assert_eq(a, a)  # Check that _meta and computed arrays match types

    with config.set({"array.backend": backend}):
        generator = gen(5) if gen else None
        rng = da.random.default_rng(generator)

        rnd_test(rng.beta, 1, 2, size=shape, chunks=3)
        rnd_test(rng.binomial, 10, 0.5, size=shape, chunks=3)
        rnd_test(rng.chisquare, 1, size=shape, chunks=3)
        rnd_test(rng.exponential, 1, size=shape, chunks=3)
        rnd_test(rng.f, 1, 2, size=shape, chunks=3)
        rnd_test(rng.gamma, 5, 1, size=shape, chunks=3)
        rnd_test(rng.geometric, 1, size=shape, chunks=3)
        rnd_test(rng.hypergeometric, 1, 2, 3, size=shape, chunks=3)
        rnd_test(rng.integers, 1, high=10, size=shape, chunks=3)
        rnd_test(rng.logseries, 0.5, size=shape, chunks=3)
        rnd_test(rng.poisson, 1, size=shape, chunks=3)
        rnd_test(rng.power, 1, size=shape, chunks=3)
        rnd_test(rng.random, size=shape, chunks=3)
        rnd_test(rng.standard_exponential, size=shape, chunks=3)
        rnd_test(rng.standard_gamma, 2, size=shape, chunks=3)
        rnd_test(rng.standard_normal, size=shape, chunks=3)


@pytest.mark.parametrize("backend", ["cupy", "numpy"])
def test_random_Generator_processes(backend):
    with config.set({"array.backend": backend}):
        # Check that matching seeds produce consistent results
        # with scheduler="processes"
        state = da.random.default_rng(5)
        x = state.standard_normal(size=(2, 3), chunks=3)
        state = da.random.default_rng(5)
        y = state.standard_normal(size=(2, 3), chunks=3)
        assert_eq(x, y, scheduler="processes")


def test_cupy_unsupported():
    with config.set({"array.backend": "cupy"}):
        # permutation supported for np-backed BitGenerator
        x = da.arange(12, chunks=3)
        da.random.default_rng(np.random.PCG64()).permutation(x).compute()

        # permutation not supported for default cupy BitGenerator
        with pytest.raises(NotImplementedError):
            da.random.default_rng().permutation(x).compute()

        # choice not supported for cupy-backed Generator
        with pytest.raises(NotImplementedError):
            da.random.default_rng().choice(10).compute()


@pytest.mark.parametrize("shape", [(2, 3), (2, 3, 4), (2, 3, 4, 2)])
def test_random_shapes(shape):
    rs = da.random.RandomState(RandomState=cupy.random.RandomState)

    x = rs.poisson(size=shape, chunks=3)
    assert type(x._meta) == cupy.ndarray
    assert_eq(x, x)  # Check that _meta and computed arrays match types
    assert x._meta.shape == (0,) * len(shape)
    assert x.shape == shape

    rng = da.random.default_rng(cupy.random.default_rng())

    x = rng.poisson(1.0, size=shape, chunks=3)
    assert type(x._meta) == cupy.ndarray
    assert_eq(x, x)  # Check that _meta and computed arrays match types
    assert x._meta.shape == (0,) * len(shape)
    assert x.shape == shape
