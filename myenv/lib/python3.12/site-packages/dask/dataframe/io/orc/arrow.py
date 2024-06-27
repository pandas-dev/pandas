from __future__ import annotations

import pyarrow as pa
import pyarrow.orc as orc

from dask.dataframe.io.utils import _get_pyarrow_dtypes, _meta_from_dtypes


class ArrowORCEngine:
    @classmethod
    def read_metadata(
        cls,
        fs,
        paths,
        columns,
        index,
        split_stripes,
        aggregate_files,
        **kwargs,
    ):
        # Convert root directory to file list.
        # TODO: Handle hive-partitioned data
        if len(paths) == 1 and not fs.isfile(paths[0]):
            paths = fs.find(paths[0])

        schema = None
        parts = []

        def _get_schema(_o, schema):
            if schema is None:
                schema = _o.schema
            elif schema != _o.schema:
                raise ValueError("Incompatible schemas while parsing ORC files")
            return schema

        if split_stripes:
            offset = 0
            for path in paths:
                with fs.open(path, "rb") as f:
                    o = orc.ORCFile(f)
                    if schema is None:
                        schema = o.schema
                    elif schema != o.schema:
                        raise ValueError("Incompatible schemas while parsing ORC files")
                    _stripes = list(range(o.nstripes))
                    if offset:
                        parts.append([(path, _stripes[0:offset])])
                    while offset < o.nstripes:
                        parts.append(
                            [(path, _stripes[offset : offset + int(split_stripes)])]
                        )
                        offset += int(split_stripes)
                    if aggregate_files and int(split_stripes) > 1:
                        offset -= o.nstripes
                    else:
                        offset = 0
        else:
            for path in paths:
                if schema is None:
                    with fs.open(paths[0], "rb") as f:
                        o = orc.ORCFile(f)
                        schema = o.schema
                parts.append([(path, None)])

        schema = _get_pyarrow_dtypes(schema, categories=None)
        if columns is not None:
            ex = set(columns) - set(schema)
            if ex:
                raise ValueError(
                    f"Requested columns ({ex}) not in schema ({set(schema)})"
                )

        # Check if we can aggregate adjacent parts together
        parts = cls._aggregate_files(aggregate_files, split_stripes, parts)

        columns = list(schema) if columns is None else columns
        index = [index] if isinstance(index, str) else index
        meta = _meta_from_dtypes(columns, schema, index, [])
        return parts, schema, meta

    @classmethod
    def _aggregate_files(cls, aggregate_files, split_stripes, parts):
        if aggregate_files is True and int(split_stripes) > 1 and len(parts) > 1:
            new_parts = []
            new_part = parts[0]
            nstripes = len(new_part[0][1])
            for part in parts[1:]:
                next_nstripes = len(part[0][1])
                if next_nstripes + nstripes <= split_stripes:
                    new_part.append(part[0])
                    nstripes += next_nstripes
                else:
                    new_parts.append(new_part)
                    new_part = part
                    nstripes = next_nstripes
            new_parts.append(new_part)
            return new_parts
        else:
            return parts

    @classmethod
    def read_partition(cls, fs, parts, schema, columns, **kwargs):
        batches = []
        for path, stripes in parts:
            batches += _read_orc_stripes(fs, path, stripes, schema, columns)
        return pa.Table.from_batches(batches).to_pandas(date_as_object=False)

    @classmethod
    def write_partition(cls, df, path, fs, filename, **kwargs):
        table = pa.Table.from_pandas(df)
        with fs.open(fs.sep.join([path, filename]), "wb") as f:
            orc.write_table(table, f)


def _read_orc_stripes(fs, path, stripes, schema, columns):
    # Construct a list of RecordBatch objects.
    # Each ORC stripe will corresonpond to a single RecordBatch.
    if columns is None:
        columns = list(schema)

    batches = []
    with fs.open(path, "rb") as f:
        o = orc.ORCFile(f)
        _stripes = range(o.nstripes) if stripes is None else stripes
        for stripe in _stripes:
            batches.append(o.read_stripe(stripe, columns))
    return batches
