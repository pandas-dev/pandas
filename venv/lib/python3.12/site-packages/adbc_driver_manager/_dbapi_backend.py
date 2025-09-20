# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the
# specific language governing permissions and limitations
# under the License.

"""
Backend-specific operations for the DB-API layer.

These are mostly functions that convert Python types to/from Arrow types.
They are abstracted so that we can support multiple backends like PyArrow,
polars, and nanoarrow.
"""

import abc
import typing

from . import _lib


class DbapiBackend(abc.ABC):
    """
    Python/Arrow type conversions that the DB-API layer needs.

    The return types can and should vary based on the backend.
    """

    @abc.abstractmethod
    def convert_bind_parameters(self, parameters: typing.Any) -> typing.Any:
        """Convert an arbitrary Python object into bind parameters.

        Parameters
        ----------
        parameters
            A sequence of bind parameters.  For instance: a tuple, where each
            item is a bind parameter in sequence.

        Returns
        -------
        parameters : CapsuleType
            This should be an Arrow stream capsule or an object implementing
            the Arrow PyCapsule interface.

        See Also
        --------
        https://arrow.apache.org/docs/format/CDataInterface/PyCapsuleInterface.html

        """
        ...

    @abc.abstractmethod
    def convert_executemany_parameters(
        self, parameters: typing.Any
    ) -> typing.Tuple[typing.Any, bool]:
        """Convert an arbitrary Python sequence into bind parameters.

        Parameters
        ----------
        parameters
            A sequence of bind parameters.  For instance: an iterable of
            tuples, where each tuple is a row of bind parameters.

        Returns
        -------
        parameters : CapsuleType
            This should be an Arrow stream capsule or an object implementing
            the Arrow PyCapsule interface.
        bind_by_name : bool
            Whether the parameters should be bound by name (e.g. because they
            contain a dictionary).

        See Also
        --------
        https://arrow.apache.org/docs/format/CDataInterface/PyCapsuleInterface.html

        """
        ...

    @abc.abstractmethod
    def import_array_stream(self, handle: _lib.ArrowArrayStreamHandle) -> typing.Any:
        """Import an Arrow stream."""
        ...

    @abc.abstractmethod
    def import_schema(self, handle: _lib.ArrowSchemaHandle) -> typing.Any:
        """Import an Arrow schema."""
        ...


_ALL_BACKENDS: list[DbapiBackend] = []


def default_backend() -> DbapiBackend:
    return _ALL_BACKENDS[-1]


class _NoOpBackend(DbapiBackend):
    def convert_bind_parameters(self, parameters: typing.Any) -> typing.Any:
        raise _lib.ProgrammingError(
            "This API requires PyArrow or another suitable backend to be installed",
            status_code=_lib.AdbcStatusCode.INVALID_STATE,
        )

    def convert_executemany_parameters(
        self, parameters: typing.Any
    ) -> typing.Tuple[typing.Any, bool]:
        raise _lib.ProgrammingError(
            "This API requires PyArrow or another suitable backend to be installed",
            status_code=_lib.AdbcStatusCode.INVALID_STATE,
        )

    def import_array_stream(
        self, handle: _lib.ArrowArrayStreamHandle
    ) -> _lib.ArrowArrayStreamHandle:
        return handle

    def import_schema(self, handle: _lib.ArrowSchemaHandle) -> _lib.ArrowSchemaHandle:
        return handle


def param_iterable_to_dict(parameters: typing.Any) -> typing.Tuple[dict, bool]:
    bind_by_name = False
    cols = {}
    for param in parameters:
        if not cols:
            # First iteration
            if isinstance(param, dict):
                bind_by_name = True
                for k, v in param.items():
                    cols[str(k)] = []
            else:
                for col_idx, v in enumerate(param):
                    cols[str(col_idx)] = []

        if isinstance(param, dict):
            for k, v in param.items():
                cols[str(k)].append(v)
        else:
            for col_idx, v in enumerate(param):
                cols[str(col_idx)].append(v)
    return cols, bind_by_name


_ALL_BACKENDS.append(_NoOpBackend())

try:
    import polars

    class _PolarsBackend(DbapiBackend):
        def convert_bind_parameters(self, parameters: typing.Any) -> polars.DataFrame:
            if isinstance(parameters, dict):
                return polars.DataFrame(
                    {str(k): v for k, v in parameters.items()},
                )

            return polars.DataFrame(
                {str(col_idx): v for col_idx, v in enumerate(parameters)},
            )

        def convert_executemany_parameters(
            self, parameters: typing.Any
        ) -> typing.Tuple[typing.Any, bool]:
            cols, bind_by_name = param_iterable_to_dict(parameters)
            return polars.DataFrame(cols), bind_by_name

        def import_array_stream(
            self, handle: _lib.ArrowArrayStreamHandle
        ) -> typing.Any:
            return polars.from_arrow(handle)

        def import_schema(self, handle: _lib.ArrowSchemaHandle) -> typing.Any:
            raise _lib.NotSupportedError("Polars does not support __arrow_c_schema__")

    _ALL_BACKENDS.append(_PolarsBackend())
except ImportError:
    pass

# Keep PyArrow at the end so it stays default
try:
    import pyarrow

    class _PyArrowBackend(DbapiBackend):
        def convert_bind_parameters(self, parameters: typing.Any) -> typing.Any:
            if isinstance(parameters, dict):
                return pyarrow.record_batch(
                    {str(k): [v] for k, v in parameters.items()},
                )
            return pyarrow.record_batch(
                [[param_value] for param_value in parameters],
                names=[str(i) for i in range(len(parameters))],
            )

        def convert_executemany_parameters(
            self, parameters: typing.Any
        ) -> typing.Tuple[typing.Any, bool]:
            cols, bind_by_name = param_iterable_to_dict(parameters)
            return pyarrow.RecordBatch.from_pydict(cols), bind_by_name

        def import_array_stream(
            self, handle: _lib.ArrowArrayStreamHandle
        ) -> pyarrow.RecordBatchReader:
            return pyarrow.RecordBatchReader._import_from_c(handle.address)

        def import_schema(self, handle: _lib.ArrowSchemaHandle) -> pyarrow.Schema:
            return pyarrow.schema(handle)

    _ALL_BACKENDS.append(_PyArrowBackend())

except ImportError:
    pass
