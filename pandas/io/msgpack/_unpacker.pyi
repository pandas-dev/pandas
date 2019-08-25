# flake8: noqa

def unpackb(
    packed,
    object_hook=...,
    list_hook=...,
    use_list=...,
    encoding=...,
    unicode_errors=...,
    object_pairs_hook=...,
    ext_hook=...,
    max_str_len=...,
    max_bin_len=...,
    max_array_len=...,
    max_map_len=...,
    max_ext_len=...,
): ...
def unpack(
    stream,
    object_hook=...,
    list_hook=...,
    use_list=...,
    encoding=...,
    unicode_errors=...,
    object_pairs_hook=...,
): ...

class Unpacker:
    def __cinit__(self): ...
    def __dealloc__(self): ...
    def __init__(
        self,
        file_like=...,
        read_size=...,
        use_list=...,
        object_hook=...,
        object_pairs_hook=...,
        list_hook=...,
        encoding=...,
        unicode_errors=...,
        max_buffer_size: int = ...,
        ext_hook=...,
        max_str_len=...,
        max_bin_len=...,
        max_array_len=...,
        max_map_len=...,
        max_ext_len=...,
    ): ...
    def feed(self, next_bytes): ...
    def append_buffer(self, _buf, _buf_len): ...
    def read_from_file(self): ...
    def _unpack(self, execute, write_bytes, iter=...): ...
    def read_bytes(self, nbytes): ...
    def unpack(self, write_bytes=...): ...
    def skip(self, write_bytes=...): ...
    def read_array_header(self, write_bytes=...): ...
    def read_map_header(self, write_bytes=...): ...
    def __iter__(self): ...
    def __next__(self): ...

