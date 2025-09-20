from __future__ import annotations

import os
from contextlib import contextmanager

import pytest

pytest.importorskip("skimage")
import numpy as np
from skimage.io import imsave

from dask.array.image import imread as da_imread
from dask.utils import tmpdir


@contextmanager
def random_images(n, shape):
    with tmpdir() as dirname:
        for i in range(n):
            fn = os.path.join(dirname, "image.%d.png" % i)
            x = np.random.randint(0, 255, size=shape).astype("u1")
            imsave(fn, x, check_contrast=False)

        yield os.path.join(dirname, "*.png")


def test_imread():
    with random_images(4, (5, 6, 3)) as globstring:
        im = da_imread(globstring)
        assert im.shape == (4, 5, 6, 3)
        assert im.chunks == ((1, 1, 1, 1), (5,), (6,), (3,))
        assert im.dtype == "uint8"

        assert im.compute().shape == (4, 5, 6, 3)
        assert im.compute().dtype == "uint8"


def test_imread_with_custom_function():
    def imread2(fn):
        return np.ones((2, 3, 4), dtype="i1")

    with random_images(4, (5, 6, 3)) as globstring:
        im = da_imread(globstring, imread=imread2)
        assert (im.compute() == np.ones((4, 2, 3, 4), dtype="u1")).all()


def test_preprocess():
    def preprocess(x):
        x[:] = 1
        return x[:, :, 0]

    with random_images(4, (2, 3, 4)) as globstring:
        im = da_imread(globstring, preprocess=preprocess)
        assert (im.compute() == np.ones((4, 2, 3), dtype="u1")).all()
