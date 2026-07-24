from __future__ import annotations


def get_short_path_name(long_name: str) -> str:
    """Gets the short path name of a given long path - http://stackoverflow.com/a/23598461/200291."""
    import ctypes  # ruff:ignore[import-outside-top-level]
    from ctypes import wintypes  # ruff:ignore[import-outside-top-level]

    GetShortPathNameW = ctypes.windll.kernel32.GetShortPathNameW  # ruff:ignore[non-lowercase-variable-in-function]  # ty: ignore[unresolved-attribute]
    GetShortPathNameW.argtypes = [wintypes.LPCWSTR, wintypes.LPWSTR, wintypes.DWORD]
    GetShortPathNameW.restype = wintypes.DWORD
    output_buf_size = 0
    while True:
        output_buf = ctypes.create_unicode_buffer(output_buf_size)
        needed = GetShortPathNameW(long_name, output_buf, output_buf_size)
        if output_buf_size >= needed:
            return output_buf.value
        output_buf_size = needed


__all__ = [
    "get_short_path_name",
]
