from collections.abc import Sequence
from typing import Literal, TypeVar, overload

from tensorflow import RaggedTensor, Tensor
from tensorflow._aliases import StringTensorCompatible, TensorCompatible, UIntTensorCompatible
from tensorflow.dtypes import DType

_TensorOrRaggedTensor = TypeVar("_TensorOrRaggedTensor", Tensor, RaggedTensor)

@overload
def as_string(
    input: TensorCompatible,
    precision: int = -1,
    scientific: bool = False,
    shortest: bool = False,
    width: int = -1,
    fill: str = "",
    name: str | None = None,
) -> Tensor: ...
@overload
def as_string(
    input: RaggedTensor,
    precision: int = -1,
    scientific: bool = False,
    shortest: bool = False,
    width: int = -1,
    fill: str = "",
    name: str | None = None,
) -> RaggedTensor: ...
def bytes_split(input: TensorCompatible | RaggedTensor, name: str | None = None) -> RaggedTensor: ...
def format(
    template: str, inputs: TensorCompatible, placeholder: str = "{}", summarize: int = 3, name: str | None = None
) -> Tensor: ...
def join(inputs: Sequence[TensorCompatible | RaggedTensor], separator: str = "", name: str | None = None) -> Tensor: ...
@overload
def length(input: TensorCompatible, unit: Literal["BYTE", "UTF8_CHAR"] = "BYTE", name: str | None = None) -> Tensor: ...
@overload
def length(input: RaggedTensor, unit: Literal["BYTE", "UTF8_CHAR"] = "BYTE", name: str | None = None) -> RaggedTensor: ...
@overload
def lower(input: TensorCompatible, encoding: Literal["utf-8", ""] = "", name: str | None = None) -> Tensor: ...
@overload
def lower(input: RaggedTensor, encoding: Literal["utf-8", ""] = "", name: str | None = None) -> RaggedTensor: ...
def ngrams(
    data: StringTensorCompatible | RaggedTensor,
    ngram_width: int | Sequence[int],
    separator: str = " ",
    pad_values: tuple[int, int] | str | None = None,
    padding_width: int | None = None,
    preserve_short_sequences: bool = False,
    name: str | None = None,
) -> RaggedTensor: ...
def reduce_join(
    inputs: StringTensorCompatible | RaggedTensor,
    axis: int | None = None,
    keepdims: bool = False,
    separator: str = "",
    name: str | None = None,
) -> Tensor: ...
@overload
def regex_full_match(input: StringTensorCompatible, pattern: StringTensorCompatible, name: str | None = None) -> Tensor: ...
@overload
def regex_full_match(input: RaggedTensor, pattern: StringTensorCompatible, name: str | None = None) -> RaggedTensor: ...
@overload
def regex_replace(
    input: StringTensorCompatible,
    pattern: StringTensorCompatible,
    rewrite: StringTensorCompatible,
    replace_global: bool = True,
    name: str | None = None,
) -> Tensor: ...
@overload
def regex_replace(
    input: RaggedTensor,
    pattern: StringTensorCompatible,
    rewrite: StringTensorCompatible,
    replace_global: bool = True,
    name: str | None = None,
) -> RaggedTensor: ...
def split(
    input: StringTensorCompatible | RaggedTensor,
    sep: StringTensorCompatible | None = None,
    maxsplit: int = -1,
    name: str | None = None,
) -> RaggedTensor: ...
@overload
def strip(input: StringTensorCompatible, name: str | None = None) -> Tensor: ...
@overload
def strip(input: RaggedTensor, name: str | None = None) -> RaggedTensor: ...
@overload
def substr(
    input: StringTensorCompatible,
    pos: TensorCompatible,
    len: TensorCompatible,
    unit: Literal["BYTE", "UTF8_CHAR"] = "BYTE",
    name: str | None = None,
) -> Tensor: ...
@overload
def substr(
    input: RaggedTensor,
    pos: TensorCompatible,
    len: TensorCompatible,
    unit: Literal["BYTE", "UTF8_CHAR"] = "BYTE",
    name: str | None = None,
) -> RaggedTensor: ...
@overload
def to_hash_bucket(input: StringTensorCompatible, num_buckets: int, name: str | None = None) -> Tensor: ...
@overload
def to_hash_bucket(input: RaggedTensor, num_buckets: int, name: str | None = None) -> RaggedTensor: ...
@overload
def to_hash_bucket_fast(input: StringTensorCompatible, num_buckets: int, name: str | None = None) -> Tensor: ...
@overload
def to_hash_bucket_fast(input: RaggedTensor, num_buckets: int, name: str | None = None) -> RaggedTensor: ...
@overload
def to_hash_bucket_strong(
    input: StringTensorCompatible, num_buckets: int, key: Sequence[int], name: str | None = None
) -> Tensor: ...
@overload
def to_hash_bucket_strong(input: RaggedTensor, num_buckets: int, key: Sequence[int], name: str | None = None) -> RaggedTensor: ...
@overload
def to_number(input: StringTensorCompatible, out_type: DType = ..., name: str | None = None) -> Tensor: ...
@overload
def to_number(input: RaggedTensor, out_type: DType = ..., name: str | None = None) -> RaggedTensor: ...
@overload
def unicode_decode(
    input: StringTensorCompatible,
    input_encoding: str,
    errors: Literal["replace", "strict", "ignore"] = "replace",
    replacement_char: int = 65533,
    replace_control_characters: bool = False,
    name: str | None = None,
) -> Tensor | RaggedTensor: ...
@overload
def unicode_decode(
    input: RaggedTensor,
    input_encoding: str,
    errors: Literal["replace", "strict", "ignore"] = "replace",
    replacement_char: int = 65533,
    replace_control_characters: bool = False,
    name: str | None = None,
) -> RaggedTensor: ...
@overload
def unicode_decode_with_offsets(
    input: StringTensorCompatible,
    input_encoding: str,
    errors: Literal["replace", "strict", "ignore"] = "replace",
    replacement_char: int = 65533,
    replace_control_characters: bool = False,
    name: str | None = None,
) -> tuple[_TensorOrRaggedTensor, _TensorOrRaggedTensor]: ...
@overload
def unicode_decode_with_offsets(
    input: RaggedTensor,
    input_encoding: str,
    errors: Literal["replace", "strict", "ignore"] = "replace",
    replacement_char: int = 65533,
    replace_control_characters: bool = False,
    name: str | None = None,
) -> tuple[RaggedTensor, RaggedTensor]: ...
@overload
def unicode_encode(
    input: TensorCompatible,
    output_encoding: Literal["UTF-8", "UTF-16-BE", "UTF-32-BE"],
    errors: Literal["replace", "strict", "ignore"] = "replace",
    replacement_char: int = 65533,
    name: str | None = None,
) -> Tensor: ...
@overload
def unicode_encode(
    input: RaggedTensor,
    output_encoding: Literal["UTF-8", "UTF-16-BE", "UTF-32-BE"],
    errors: Literal["replace", "strict", "ignore"] = "replace",
    replacement_char: int = 65533,
    name: str | None = None,
) -> RaggedTensor: ...
@overload
def unicode_script(input: TensorCompatible, name: str | None = None) -> Tensor: ...
@overload
def unicode_script(input: RaggedTensor, name: str | None = None) -> RaggedTensor: ...
@overload
def unicode_split(
    input: StringTensorCompatible,
    input_encoding: str,
    errors: Literal["replace", "strict", "ignore"] = "replace",
    replacement_char: int = 65533,
    name: str | None = None,
) -> Tensor | RaggedTensor: ...
@overload
def unicode_split(
    input: RaggedTensor,
    input_encoding: str,
    errors: Literal["replace", "strict", "ignore"] = "replace",
    replacement_char: int = 65533,
    name: str | None = None,
) -> RaggedTensor: ...
@overload
def unicode_split_with_offsets(
    input: StringTensorCompatible,
    input_encoding: str,
    errors: Literal["replace", "strict", "ignore"] = "replace",
    replacement_char: int = 65533,
    name: str | None = None,
) -> tuple[_TensorOrRaggedTensor, _TensorOrRaggedTensor]: ...
@overload
def unicode_split_with_offsets(
    input: RaggedTensor,
    input_encoding: str,
    errors: Literal["replace", "strict", "ignore"] = "replace",
    replacement_char: int = 65533,
    name: str | None = None,
) -> tuple[RaggedTensor, RaggedTensor]: ...
@overload
def unicode_transcode(
    input: StringTensorCompatible,
    input_encoding: str,
    output_encoding: Literal["UTF-8", "UTF-16-BE", "UTF-32-BE"],
    errors: Literal["replace", "strict", "ignore"] = "replace",
    replacement_char: int = 65533,
    replace_control_characters: bool = False,
    name: str | None = None,
) -> Tensor: ...
@overload
def unicode_transcode(
    input: RaggedTensor,
    input_encoding: str,
    output_encoding: Literal["UTF-8", "UTF-16-BE", "UTF-32-BE"],
    errors: Literal["replace", "strict", "ignore"] = "replace",
    replacement_char: int = 65533,
    replace_control_characters: bool = False,
    name: str | None = None,
) -> RaggedTensor: ...
def unsorted_segment_join(
    inputs: StringTensorCompatible,
    segment_ids: UIntTensorCompatible,
    num_segments: UIntTensorCompatible,
    separator: str = "",
    name: str | None = None,
) -> Tensor: ...
@overload
def upper(input: TensorCompatible, encoding: Literal["utf-8", ""] = "", name: str | None = None) -> Tensor: ...
@overload
def upper(input: RaggedTensor, encoding: Literal["utf-8", ""] = "", name: str | None = None) -> RaggedTensor: ...
