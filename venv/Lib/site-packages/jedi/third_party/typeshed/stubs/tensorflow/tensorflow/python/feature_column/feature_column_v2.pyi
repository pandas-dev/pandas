# The types here are all undocumented, but all feature columns are return types of the
# public functions in tf.feature_column. As they are undocumented internals while some
# common methods are included, they are incomplete and do not have getattr Incomplete fallback.

from _typeshed import Incomplete
from abc import ABC, ABCMeta, abstractmethod
from collections.abc import Callable, Sequence
from typing import Literal
from typing_extensions import Self, TypeAlias

import tensorflow as tf
from tensorflow._aliases import ShapeLike

_Combiners: TypeAlias = Literal["mean", "sqrtn", "sum"]
_ExampleSpec: TypeAlias = dict[str, tf.io.FixedLenFeature | tf.io.VarLenFeature]

class _FeatureColumn(ABC):
    @property
    @abstractmethod
    def name(self) -> str: ...
    @property
    @abstractmethod
    def parse_example_spec(self) -> _ExampleSpec: ...
    def __lt__(self, other: _FeatureColumn) -> bool: ...
    def __gt__(self, other: _FeatureColumn) -> bool: ...
    @property
    @abstractmethod
    def parents(self) -> list[_FeatureColumn | str]: ...

class DenseColumn(_FeatureColumn, metaclass=ABCMeta): ...
class SequenceDenseColumn(_FeatureColumn, metaclass=ABCMeta): ...

# These classes are mostly subclasses of collections.namedtuple but we can't use
# typing.NamedTuple because they use multiple inheritance with other non namedtuple classes.
# _cls instead of cls is because collections.namedtuple uses _cls for __new__.
class NumericColumn(DenseColumn):
    key: str
    shape: ShapeLike
    default_value: float
    dtype: tf.DType
    normalizer_fn: Callable[[tf.Tensor], tf.Tensor] | None

    def __new__(
        _cls,
        key: str,
        shape: ShapeLike,
        default_value: float,
        dtype: tf.DType,
        normalizer_fn: Callable[[tf.Tensor], tf.Tensor] | None,
    ) -> Self: ...
    @property
    def name(self) -> str: ...
    @property
    def parse_example_spec(self) -> _ExampleSpec: ...
    @property
    def parents(self) -> list[_FeatureColumn | str]: ...

class CategoricalColumn(_FeatureColumn):
    @property
    @abstractmethod
    def num_buckets(self) -> int: ...

class BucketizedColumn(DenseColumn, CategoricalColumn):
    source_column: NumericColumn
    boundaries: list[float] | tuple[float, ...]

    def __new__(_cls, source_column: NumericColumn, boundaries: list[float] | tuple[float, ...]) -> Self: ...
    @property
    def name(self) -> str: ...
    @property
    def num_buckets(self) -> int: ...
    @property
    def parse_example_spec(self) -> _ExampleSpec: ...
    @property
    def parents(self) -> list[_FeatureColumn | str]: ...

class EmbeddingColumn(DenseColumn, SequenceDenseColumn):
    categorical_column: CategoricalColumn
    dimension: int
    combiner: _Combiners
    initializer: Callable[[ShapeLike], tf.Tensor] | None
    ckpt_to_load_from: str | None
    tensor_name_in_ckpt: str | None
    max_norm: float | None
    trainable: bool
    use_safe_embedding_lookup: bool

    # This one subclasses collections.namedtuple and overrides __new__.
    def __new__(
        cls,
        categorical_column: CategoricalColumn,
        dimension: int,
        combiner: _Combiners,
        initializer: Callable[[ShapeLike], tf.Tensor] | None,
        ckpt_to_load_from: str | None,
        tensor_name_in_ckpt: str | None,
        max_norm: float | None,
        trainable: bool,
        use_safe_embedding_lookup: bool = True,
    ) -> Self: ...
    @property
    def name(self) -> str: ...
    @property
    def parse_example_spec(self) -> _ExampleSpec: ...
    @property
    def parents(self) -> list[_FeatureColumn | str]: ...

class SharedEmbeddingColumnCreator:
    def __init__(
        self,
        dimension: int,
        initializer: Callable[[ShapeLike], tf.Tensor] | None,
        ckpt_to_load_from: str | None,
        tensor_name_in_ckpt: str | None,
        num_buckets: int,
        trainable: bool,
        name: str = "shared_embedding_column_creator",
        use_safe_embedding_lookup: bool = True,
    ) -> None: ...
    def __getattr__(self, name: str) -> Incomplete: ...

class SharedEmbeddingColumn(DenseColumn, SequenceDenseColumn):
    categorical_column: CategoricalColumn
    shared_embedding_column_creator: SharedEmbeddingColumnCreator
    combiner: _Combiners
    max_norm: float | None
    use_safe_embedding_lookup: bool

    def __new__(
        cls,
        categorical_column: CategoricalColumn,
        shared_embedding_column_creator: SharedEmbeddingColumnCreator,
        combiner: _Combiners,
        max_norm: float | None,
        use_safe_embedding_lookup: bool = True,
    ) -> Self: ...
    @property
    def name(self) -> str: ...
    @property
    def parse_example_spec(self) -> _ExampleSpec: ...
    @property
    def parents(self) -> list[_FeatureColumn | str]: ...

class CrossedColumn(CategoricalColumn):
    keys: tuple[str, ...]
    hash_bucket_size: int
    hash_key: int | None

    def __new__(_cls, keys: tuple[str, ...], hash_bucket_size: int, hash_key: int | None) -> Self: ...
    @property
    def name(self) -> str: ...
    @property
    def num_buckets(self) -> int: ...
    @property
    def parse_example_spec(self) -> _ExampleSpec: ...
    @property
    def parents(self) -> list[_FeatureColumn | str]: ...

class IdentityCategoricalColumn(CategoricalColumn):
    key: str
    number_buckets: int
    default_value: int | None

    def __new__(_cls, key: str, number_buckets: int, default_value: int | None) -> Self: ...
    @property
    def name(self) -> str: ...
    @property
    def num_buckets(self) -> int: ...
    @property
    def parse_example_spec(self) -> _ExampleSpec: ...
    @property
    def parents(self) -> list[_FeatureColumn | str]: ...

class HashedCategoricalColumn(CategoricalColumn):
    key: str
    hash_bucket_size: int
    dtype: tf.DType

    def __new__(_cls, key: str, hash_bucket_size: int, dtype: tf.DType) -> Self: ...
    @property
    def name(self) -> str: ...
    @property
    def num_buckets(self) -> int: ...
    @property
    def parse_example_spec(self) -> _ExampleSpec: ...
    @property
    def parents(self) -> list[_FeatureColumn | str]: ...

class VocabularyFileCategoricalColumn(CategoricalColumn):
    key: str
    vocabulary_file: str
    vocabulary_size: int | None
    num_oov_buckets: int
    dtype: tf.DType
    default_value: str | int | None
    file_format: str | None

    def __new__(
        cls,
        key: str,
        vocabulary_file: str,
        vocabulary_size: int | None,
        num_oov_buckets: int,
        dtype: tf.DType,
        default_value: str | int | None,
        file_format: str | None = None,
    ) -> Self: ...
    @property
    def name(self) -> str: ...
    @property
    def num_buckets(self) -> int: ...
    @property
    def parse_example_spec(self) -> _ExampleSpec: ...
    @property
    def parents(self) -> list[_FeatureColumn | str]: ...

class VocabularyListCategoricalColumn(CategoricalColumn):
    key: str
    vocabulary_list: Sequence[str] | Sequence[int]
    dtype: tf.DType
    default_value: str | int | None
    num_oov_buckets: int

    def __new__(
        _cls, key: str, vocabulary_list: Sequence[str], dtype: tf.DType, default_value: str | int | None, num_oov_buckets: int
    ) -> Self: ...
    @property
    def name(self) -> str: ...
    @property
    def num_buckets(self) -> int: ...
    @property
    def parse_example_spec(self) -> _ExampleSpec: ...
    @property
    def parents(self) -> list[_FeatureColumn | str]: ...

class WeightedCategoricalColumn(CategoricalColumn):
    categorical_column: CategoricalColumn
    weight_feature_key: str
    dtype: tf.DType

    def __new__(_cls, categorical_column: CategoricalColumn, weight_feature_key: str, dtype: tf.DType) -> Self: ...
    @property
    def name(self) -> str: ...
    @property
    def num_buckets(self) -> int: ...
    @property
    def parse_example_spec(self) -> _ExampleSpec: ...
    @property
    def parents(self) -> list[_FeatureColumn | str]: ...

class IndicatorColumn(DenseColumn, SequenceDenseColumn):
    categorical_column: CategoricalColumn

    def __new__(_cls, categorical_column: CategoricalColumn) -> Self: ...
    @property
    def name(self) -> str: ...
    @property
    def parse_example_spec(self) -> _ExampleSpec: ...
    @property
    def parents(self) -> list[_FeatureColumn | str]: ...

class SequenceCategoricalColumn(CategoricalColumn):
    categorical_column: CategoricalColumn

    def __new__(_cls, categorical_column: CategoricalColumn) -> Self: ...
    @property
    def name(self) -> str: ...
    @property
    def num_buckets(self) -> int: ...
    @property
    def parse_example_spec(self) -> _ExampleSpec: ...
    @property
    def parents(self) -> list[_FeatureColumn | str]: ...
