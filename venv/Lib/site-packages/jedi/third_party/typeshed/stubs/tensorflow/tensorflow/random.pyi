from collections.abc import Sequence
from enum import Enum
from typing import Literal
from typing_extensions import TypeAlias

import numpy as np
import numpy.typing as npt
import tensorflow as tf
from tensorflow._aliases import DTypeLike, ScalarTensorCompatible, ShapeLike
from tensorflow.python.trackable import autotrackable

class Algorithm(Enum):
    PHILOX = 1
    THREEFRY = 2
    AUTO_SELECT = 3

_Alg: TypeAlias = Literal[Algorithm.PHILOX, Algorithm.THREEFRY, Algorithm.AUTO_SELECT, "philox", "threefry", "auto_select"]

class Generator(autotrackable.AutoTrackable):
    @classmethod
    def from_state(cls, state: tf.Variable, alg: _Alg | None) -> Generator: ...
    @classmethod
    def from_seed(cls, seed: int, alg: _Alg | None = None) -> Generator: ...
    @classmethod
    def from_non_deterministic_state(cls, alg: _Alg | None = None) -> Generator: ...
    @classmethod
    def from_key_counter(
        cls, key: ScalarTensorCompatible, counter: Sequence[ScalarTensorCompatible], alg: _Alg | None
    ) -> Generator: ...
    def __init__(self, copy_from: Generator | None = None, state: tf.Variable | None = None, alg: _Alg | None = None) -> None: ...
    def reset(self, state: tf.Variable) -> None: ...
    def reset_from_seed(self, seed: int) -> None: ...
    def reset_from_key_counter(self, key: ScalarTensorCompatible, counter: tf.Variable) -> None: ...
    @property
    def state(self) -> tf.Variable: ...
    @property
    def algorithm(self) -> int: ...
    @property
    def key(self) -> ScalarTensorCompatible: ...
    def skip(self, delta: int) -> tf.Tensor: ...
    def normal(
        self,
        shape: tf.Tensor | Sequence[int],
        mean: ScalarTensorCompatible = 0.0,
        stddev: ScalarTensorCompatible = 1.0,
        dtype: DTypeLike = ...,
        name: str | None = None,
    ) -> tf.Tensor: ...
    def truncated_normal(
        self,
        shape: ShapeLike,
        mean: ScalarTensorCompatible = 0.0,
        stddev: ScalarTensorCompatible = 1.0,
        dtype: DTypeLike = ...,
        name: str | None = None,
    ) -> tf.Tensor: ...
    def uniform(
        self,
        shape: ShapeLike,
        minval: ScalarTensorCompatible = 0,
        maxval: ScalarTensorCompatible | None = None,
        dtype: DTypeLike = ...,
        name: str | None = None,
    ) -> tf.Tensor: ...
    def uniform_full_int(self, shape: ShapeLike, dtype: DTypeLike = ..., name: str | None = None) -> tf.Tensor: ...
    def binomial(
        self, shape: ShapeLike, counts: tf.Tensor, probs: tf.Tensor, dtype: DTypeLike = ..., name: str | None = None
    ) -> tf.Tensor: ...
    def make_seeds(self, count: int = 1) -> tf.Tensor: ...
    def split(self, count: int = 1) -> list[Generator]: ...

def all_candidate_sampler(
    true_classes: tf.Tensor, num_true: int, num_sampled: int, unique: bool, seed: int | None = None, name: str | None = None
) -> tuple[tf.Tensor, tf.Tensor, tf.Tensor]: ...
def categorical(
    logits: tf.Tensor,
    num_samples: int | tf.Tensor,
    dtype: DTypeLike | None = None,
    seed: int | None = None,
    name: str | None = None,
) -> tf.Tensor: ...
def create_rng_state(seed: int, alg: _Alg) -> npt.NDArray[np.int64]: ...
def fixed_unigram_candidate_sampler(
    true_classes: tf.Tensor,
    num_true: int,
    num_sampled: int,
    unique: bool,
    range_max: int,
    vocab_file: str = "",
    distortion: float = 1.0,
    num_reserved_ids: int = 0,
    num_shards: int = 1,
    shard: int = 0,
    unigrams: Sequence[float] = (),
    seed: int | None = None,
    name: str | None = None,
) -> tuple[tf.Tensor, tf.Tensor, tf.Tensor]: ...
def fold_in(seed: tf.Tensor | Sequence[int], data: int, alg: _Alg = "auto_select") -> int: ...
def gamma(
    shape: tf.Tensor | Sequence[int],
    alpha: tf.Tensor | float | Sequence[float],
    beta: tf.Tensor | float | Sequence[float] | None = None,
    dtype: DTypeLike = ...,
    seed: int | None = None,
    name: str | None = None,
) -> tf.Tensor: ...
def get_global_generator() -> Generator: ...
def learned_unigram_candidate_sampler(
    true_classes: tf.Tensor,
    num_true: int,
    num_sampled: int,
    unique: bool,
    range_max: int,
    seed: int | None = None,
    name: str | None = None,
) -> tuple[tf.Tensor, tf.Tensor, tf.Tensor]: ...
def log_uniform_candidate_sampler(
    true_classes: tf.Tensor,
    num_true: int,
    num_sampled: int,
    unique: bool,
    range_max: int,
    seed: int | None = None,
    name: str | None = None,
) -> tuple[tf.Tensor, tf.Tensor, tf.Tensor]: ...
def normal(
    shape: ShapeLike,
    mean: ScalarTensorCompatible = 0.0,
    stddev: ScalarTensorCompatible = 1.0,
    dtype: DTypeLike = ...,
    seed: int | None = None,
    name: str | None = None,
) -> tf.Tensor: ...
def poisson(
    shape: ShapeLike, lam: ScalarTensorCompatible, dtype: DTypeLike = ..., seed: int | None = None, name: str | None = None
) -> tf.Tensor: ...
def set_global_generator(generator: Generator) -> None: ...
def set_seed(seed: int) -> None: ...
def shuffle(value: tf.Tensor, seed: int | None = None, name: str | None = None) -> tf.Tensor: ...
def split(seed: tf.Tensor | Sequence[int], num: int = 2, alg: _Alg = "auto_select") -> tf.Tensor: ...
def stateless_binomial(
    shape: ShapeLike,
    seed: tuple[int, int] | tf.Tensor,
    counts: tf.Tensor,
    probs: tf.Tensor,
    output_dtype: DTypeLike = ...,
    name: str | None = None,
) -> tf.Tensor: ...
def stateless_categorical(
    logits: tf.Tensor,
    num_samples: int | tf.Tensor,
    seed: tuple[int, int] | tf.Tensor,
    dtype: DTypeLike = ...,
    name: str | None = None,
) -> tf.Tensor: ...
def stateless_gamma(
    shape: ShapeLike,
    seed: tuple[int, int] | tf.Tensor,
    alpha: tf.Tensor,
    beta: tf.Tensor | None = None,
    dtype: DTypeLike = ...,
    name: str | None = None,
) -> tf.Tensor: ...
def stateless_normal(
    shape: tf.Tensor | Sequence[int],
    seed: tuple[int, int] | tf.Tensor,
    mean: float | tf.Tensor = 0.0,
    stddev: float | tf.Tensor = 1.0,
    dtype: DTypeLike = ...,
    name: str | None = None,
    alg: _Alg = "auto_select",
) -> tf.Tensor: ...
def stateless_parameterized_truncated_normal(
    shape: tf.Tensor | Sequence[int],
    seed: tuple[int, int] | tf.Tensor,
    means: float | tf.Tensor = 0.0,
    stddevs: float | tf.Tensor = 1.0,
    minvals: tf.Tensor | float = -2.0,
    maxvals: tf.Tensor | float = 2.0,
    name: str | None = None,
) -> tf.Tensor: ...
def stateless_poisson(
    shape: tf.Tensor | Sequence[int],
    seed: tuple[int, int] | tf.Tensor,
    lam: tf.Tensor,
    dtype: DTypeLike = ...,
    name: str | None = None,
) -> tf.Tensor: ...
def stateless_truncated_normal(
    shape: tf.Tensor | Sequence[int],
    seed: tuple[int, int] | tf.Tensor,
    mean: float | tf.Tensor = 0.0,
    stddev: float | tf.Tensor = 1.0,
    dtype: DTypeLike = ...,
    name: str | None = None,
    alg: _Alg = "auto_select",
) -> tf.Tensor: ...
def stateless_uniform(
    shape: tf.Tensor | Sequence[int],
    seed: tuple[int, int] | tf.Tensor,
    minval: float | tf.Tensor = 0,
    maxval: float | tf.Tensor | None = None,
    dtype: DTypeLike = ...,
    name: str | None = None,
    alg: _Alg = "auto_select",
) -> tf.Tensor: ...
def truncated_normal(
    shape: tf.Tensor | Sequence[int],
    mean: float | tf.Tensor = 0.0,
    stddev: float | tf.Tensor = 1.0,
    dtype: DTypeLike = ...,
    seed: int | None = None,
    name: str | None = None,
) -> tf.Tensor: ...
def uniform(
    shape: tf.Tensor | Sequence[int],
    minval: float | tf.Tensor = 0,
    maxval: float | tf.Tensor | None = None,
    dtype: DTypeLike = ...,
    seed: int | None = None,
    name: str | None = None,
) -> tf.Tensor: ...
def uniform_candidate_sampler(
    true_classes: tf.Tensor,
    num_true: int,
    num_sampled: int,
    unique: bool,
    range_max: int,
    seed: int | None = None,
    name: str | None = None,
) -> tuple[tf.Tensor, tf.Tensor, tf.Tensor]: ...
