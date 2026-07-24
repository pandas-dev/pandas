from collections.abc import Callable, Iterable, Sequence

import tensorflow as tf
from tensorflow._aliases import ShapeLike
from tensorflow.python.feature_column import feature_column_v2 as fc, sequence_feature_column as seq_fc

def numeric_column(
    key: str,
    shape: ShapeLike = (1,),
    default_value: float | None = None,
    dtype: tf.DType = ...,
    normalizer_fn: Callable[[tf.Tensor], tf.Tensor] | None = None,
) -> fc.NumericColumn: ...
def bucketized_column(source_column: fc.NumericColumn, boundaries: list[float] | tuple[float, ...]) -> fc.BucketizedColumn: ...
def embedding_column(
    categorical_column: fc.CategoricalColumn,
    dimension: int,
    combiner: fc._Combiners = "mean",
    initializer: Callable[[ShapeLike], tf.Tensor] | None = None,
    ckpt_to_load_from: str | None = None,
    tensor_name_in_ckpt: str | None = None,
    max_norm: float | None = None,
    trainable: bool = True,
    use_safe_embedding_lookup: bool = True,
) -> fc.EmbeddingColumn: ...
def shared_embeddings(
    categorical_columns: Iterable[fc.CategoricalColumn],
    dimension: int,
    combiner: fc._Combiners = "mean",
    initializer: Callable[[ShapeLike], tf.Tensor] | None = None,
    shared_embedding_collection_name: str | None = None,
    ckpt_to_load_from: str | None = None,
    tensor_name_in_ckpt: str | None = None,
    max_norm: float | None = None,
    trainable: bool = True,
    use_safe_embedding_lookup: bool = True,
) -> list[fc.SharedEmbeddingColumn]: ...
def categorical_column_with_identity(
    key: str, num_buckets: int, default_value: int | None = None
) -> fc.IdentityCategoricalColumn: ...
def categorical_column_with_hash_bucket(key: str, hash_bucket_size: int, dtype: tf.DType = ...) -> fc.HashedCategoricalColumn: ...
def categorical_column_with_vocabulary_file(
    key: str,
    vocabulary_file: str,
    vocabulary_size: int | None = None,
    dtype: tf.DType = ...,
    default_value: str | int | None = None,
    num_oov_buckets: int = 0,
    file_format: str | None = None,
) -> fc.VocabularyFileCategoricalColumn: ...
def categorical_column_with_vocabulary_list(
    key: str,
    vocabulary_list: Sequence[str] | Sequence[int],
    dtype: tf.DType | None = None,
    default_value: str | int | None = -1,
    num_oov_buckets: int = 0,
) -> fc.VocabularyListCategoricalColumn: ...
def indicator_column(categorical_column: fc.CategoricalColumn) -> fc.IndicatorColumn: ...
def weighted_categorical_column(
    categorical_column: fc.CategoricalColumn, weight_feature_key: str, dtype: tf.DType = ...
) -> fc.WeightedCategoricalColumn: ...
def crossed_column(
    keys: Iterable[str | fc.CategoricalColumn], hash_bucket_size: int, hash_key: int | None = None
) -> fc.CrossedColumn: ...
def sequence_numeric_column(
    key: str,
    shape: ShapeLike = (1,),
    default_value: float = 0.0,
    dtype: tf.DType = ...,
    normalizer_fn: Callable[[tf.Tensor], tf.Tensor] | None = None,
) -> seq_fc.SequenceNumericColumn: ...
def sequence_categorical_column_with_identity(
    key: str, num_buckets: int, default_value: int | None = None
) -> fc.SequenceCategoricalColumn: ...
def sequence_categorical_column_with_hash_bucket(
    key: str, hash_bucket_size: int, dtype: tf.DType = ...
) -> fc.SequenceCategoricalColumn: ...
def sequence_categorical_column_with_vocabulary_file(
    key: str,
    vocabulary_file: str,
    vocabulary_size: int | None = None,
    num_oov_buckets: int = 0,
    default_value: str | int | None = None,
    dtype: tf.DType = ...,
) -> fc.SequenceCategoricalColumn: ...
def sequence_categorical_column_with_vocabulary_list(
    key: str,
    vocabulary_list: Sequence[str] | Sequence[int],
    dtype: tf.DType | None = None,
    default_value: str | int | None = -1,
    num_oov_buckets: int = 0,
) -> fc.SequenceCategoricalColumn: ...
def make_parse_example_spec(
    feature_columns: Iterable[fc._FeatureColumn],
) -> dict[str, tf.io.FixedLenFeature | tf.io.VarLenFeature]: ...
