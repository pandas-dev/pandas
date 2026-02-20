# dialects/oracle/vector.py
# Copyright (C) 2005-2026 the SQLAlchemy authors and contributors
# <see AUTHORS file>
#
# This module is part of SQLAlchemy and is released under
# the MIT License: https://www.opensource.org/licenses/mit-license.php
# mypy: ignore-errors


from __future__ import annotations

import array
from dataclasses import dataclass
from enum import Enum
from typing import Optional
from typing import Union

from ... import types
from ...types import Float


class VectorIndexType(Enum):
    """Enum representing different types of VECTOR index structures.

    See :ref:`oracle_vector_datatype` for background.

    .. versionadded:: 2.0.41

    """

    HNSW = "HNSW"
    """
    The HNSW (Hierarchical Navigable Small World) index type.
    """
    IVF = "IVF"
    """
    The IVF (Inverted File Index) index type
    """


class VectorDistanceType(Enum):
    """Enum representing different types of vector distance metrics.

    See :ref:`oracle_vector_datatype` for background.

    .. versionadded:: 2.0.41

    """

    EUCLIDEAN = "EUCLIDEAN"
    """Euclidean distance (L2 norm).

    Measures the straight-line distance between two vectors in space.
    """
    DOT = "DOT"
    """Dot product similarity.

    Measures the algebraic similarity between two vectors.
    """
    COSINE = "COSINE"
    """Cosine similarity.

    Measures the cosine of the angle between two vectors.
    """
    MANHATTAN = "MANHATTAN"
    """Manhattan distance (L1 norm).

    Calculates the sum of absolute differences across dimensions.
    """


class VectorStorageFormat(Enum):
    """Enum representing the data format used to store vector components.

    See :ref:`oracle_vector_datatype` for background.

    .. versionadded:: 2.0.41

    """

    INT8 = "INT8"
    """
    8-bit integer format.
    """
    BINARY = "BINARY"
    """
    Binary format.
    """
    FLOAT32 = "FLOAT32"
    """
    32-bit floating-point format.
    """
    FLOAT64 = "FLOAT64"
    """
    64-bit floating-point format.
    """


class VectorStorageType(Enum):
    """Enum representing the vector type,

    See :ref:`oracle_vector_datatype` for background.

    .. versionadded:: 2.0.43

    """

    SPARSE = "SPARSE"
    """
    A Sparse vector is a vector which has zero value for
    most of its dimensions.
    """
    DENSE = "DENSE"
    """
    A Dense vector is a vector where most, if not all, elements
    hold meaningful values.
    """


@dataclass
class VectorIndexConfig:
    """Define the configuration for Oracle VECTOR Index.

    See :ref:`oracle_vector_datatype` for background.

    .. versionadded:: 2.0.41

    :param index_type: Enum value from :class:`.VectorIndexType`
     Specifies the indexing method. For HNSW, this must be
     :attr:`.VectorIndexType.HNSW`.

    :param distance: Enum value from :class:`.VectorDistanceType`
     specifies the metric for calculating distance between VECTORS.

    :param accuracy: integer. Should be in the range 0 to 100
     Specifies the accuracy of the nearest neighbor search during
     query execution.

    :param parallel: integer. Specifies degree of parallelism.

    :param hnsw_neighbors: integer. Should be in the range 0 to
     2048. Specifies the number of nearest neighbors considered
     during the search. The attribute :attr:`.VectorIndexConfig.hnsw_neighbors`
     is HNSW index specific.

    :param hnsw_efconstruction: integer. Should be in the range 0
     to 65535. Controls the trade-off between indexing speed and
     recall quality during index construction. The attribute
     :attr:`.VectorIndexConfig.hnsw_efconstruction` is HNSW index
     specific.

    :param ivf_neighbor_partitions: integer. Should be in the range
     0 to 10,000,000. Specifies the number of partitions used to
     divide the dataset. The attribute
     :attr:`.VectorIndexConfig.ivf_neighbor_partitions` is IVF index
     specific.

    :param ivf_sample_per_partition: integer. Should be between 1
     and ``num_vectors / neighbor partitions``. Specifies the
     number of samples used per partition. The attribute
     :attr:`.VectorIndexConfig.ivf_sample_per_partition` is IVF index
     specific.

    :param ivf_min_vectors_per_partition: integer. From 0 (no trimming)
     to the total number of vectors (results in 1 partition). Specifies
     the minimum number of vectors per partition. The attribute
     :attr:`.VectorIndexConfig.ivf_min_vectors_per_partition`
     is IVF index specific.

    """

    index_type: VectorIndexType = VectorIndexType.HNSW
    distance: Optional[VectorDistanceType] = None
    accuracy: Optional[int] = None
    hnsw_neighbors: Optional[int] = None
    hnsw_efconstruction: Optional[int] = None
    ivf_neighbor_partitions: Optional[int] = None
    ivf_sample_per_partition: Optional[int] = None
    ivf_min_vectors_per_partition: Optional[int] = None
    parallel: Optional[int] = None

    def __post_init__(self):
        self.index_type = VectorIndexType(self.index_type)
        for field in [
            "hnsw_neighbors",
            "hnsw_efconstruction",
            "ivf_neighbor_partitions",
            "ivf_sample_per_partition",
            "ivf_min_vectors_per_partition",
            "parallel",
            "accuracy",
        ]:
            value = getattr(self, field)
            if value is not None and not isinstance(value, int):
                raise TypeError(
                    f"{field} must be an integer if"
                    f"provided, got {type(value).__name__}"
                )


class SparseVector:
    """
    Lightweight SQLAlchemy-side version of SparseVector.
    This mimics oracledb.SparseVector.

    .. versionadded:: 2.0.43

    """

    def __init__(
        self,
        num_dimensions: int,
        indices: Union[list, array.array],
        values: Union[list, array.array],
    ):
        if not isinstance(indices, array.array) or indices.typecode != "I":
            indices = array.array("I", indices)
        if not isinstance(values, array.array):
            values = array.array("d", values)
        if len(indices) != len(values):
            raise TypeError("indices and values must be of the same length!")

        self.num_dimensions = num_dimensions
        self.indices = indices
        self.values = values

    def __str__(self):
        return (
            f"SparseVector(num_dimensions={self.num_dimensions}, "
            f"size={len(self.indices)}, typecode={self.values.typecode})"
        )


class VECTOR(types.TypeEngine):
    """Oracle VECTOR datatype.

    For complete background on using this type, see
    :ref:`oracle_vector_datatype`.

    .. versionadded:: 2.0.41

    """

    cache_ok = True

    __visit_name__ = "VECTOR"

    _typecode_map = {
        VectorStorageFormat.INT8: "b",  # Signed int
        VectorStorageFormat.BINARY: "B",  # Unsigned int
        VectorStorageFormat.FLOAT32: "f",  # Float
        VectorStorageFormat.FLOAT64: "d",  # Double
    }

    def __init__(self, dim=None, storage_format=None, storage_type=None):
        """Construct a VECTOR.

        :param dim: integer. The dimension of the VECTOR datatype. This
         should be an integer value.

        :param storage_format: VectorStorageFormat. The VECTOR storage
         type format. This should be Enum values form
         :class:`.VectorStorageFormat` INT8, BINARY, FLOAT32, or FLOAT64.

        :param storage_type: VectorStorageType. The Vector storage type. This
         should be Enum values from :class:`.VectorStorageType` SPARSE or
         DENSE.

        """

        if dim is not None and not isinstance(dim, int):
            raise TypeError("dim must be an integer")
        if storage_format is not None and not isinstance(
            storage_format, VectorStorageFormat
        ):
            raise TypeError(
                "storage_format must be an enum of type VectorStorageFormat"
            )
        if storage_type is not None and not isinstance(
            storage_type, VectorStorageType
        ):
            raise TypeError(
                "storage_type must be an enum of type VectorStorageType"
            )

        self.dim = dim
        self.storage_format = storage_format
        self.storage_type = storage_type

    def _cached_bind_processor(self, dialect):
        """
        Converts a Python-side SparseVector instance into an
        oracledb.SparseVectormor a compatible array format before
        binding it to the database.
        """

        def process(value):
            if value is None or isinstance(value, array.array):
                return value

            # Convert list to a array.array
            elif isinstance(value, list):
                typecode = self._array_typecode(self.storage_format)
                value = array.array(typecode, value)
                return value

            # Convert SqlAlchemy SparseVector to oracledb SparseVector object
            elif isinstance(value, SparseVector):
                return dialect.dbapi.SparseVector(
                    value.num_dimensions,
                    value.indices,
                    value.values,
                )

            else:
                raise TypeError(
                    """
                    Invalid input for VECTOR: expected a list, an array.array,
                    or a SparseVector object.
                    """
                )

        return process

    def _cached_result_processor(self, dialect, coltype):
        """
        Converts database-returned values into Python-native representations.
        If the value is an oracledb.SparseVector, it is converted into the
        SQLAlchemy-side SparseVector class.
        If the value is a array.array, it is converted to a plain Python list.

        """

        def process(value):
            if value is None:
                return None

            elif isinstance(value, array.array):
                return list(value)

            # Convert Oracledb SparseVector to SqlAlchemy SparseVector object
            elif isinstance(value, dialect.dbapi.SparseVector):
                return SparseVector(
                    num_dimensions=value.num_dimensions,
                    indices=value.indices,
                    values=value.values,
                )

        return process

    def _array_typecode(self, typecode):
        """
        Map storage format to array typecode.
        """
        return self._typecode_map.get(typecode, "d")

    class comparator_factory(types.TypeEngine.Comparator):
        def l2_distance(self, other):
            return self.op("<->", return_type=Float)(other)

        def inner_product(self, other):
            return self.op("<#>", return_type=Float)(other)

        def cosine_distance(self, other):
            return self.op("<=>", return_type=Float)(other)
