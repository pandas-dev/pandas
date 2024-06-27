from __future__ import annotations


class ORCEngine:
    """The API necessary to provide a new ORC reader/writer"""

    @classmethod
    def read_metadata(
        cls, fs, paths, columns, index, split_stripes, aggregate_files, **kwargs
    ):
        raise NotImplementedError()

    @classmethod
    def read_partition(cls, fs, part, columns, **kwargs):
        raise NotImplementedError()

    @classmethod
    def write_partition(cls, df, path, fs, filename, **kwargs):
        raise NotImplementedError
