from importlib.util import find_spec

HAS_CRT = find_spec("awscrt") is not None
HAS_CRC32C = find_spec("crc32c") is not None
