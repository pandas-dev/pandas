
import cramjam
import numpy as np
from fastparquet import parquet_thrift

# TODO: use stream/direct-to-buffer conversions instead of memcopy

compressions = {
    'UNCOMPRESSED': lambda x: x
}
decompressions = {
    'UNCOMPRESSED': lambda x, y: x
}

# Gzip is present regardless
COMPRESSION_LEVEL = 6


def gzip_compress_v3(data, compresslevel=COMPRESSION_LEVEL):
    return cramjam.gzip.compress(data, level=compresslevel)


def gzip_decompress(data, uncompressed_size):
    return cramjam.gzip.decompress(data, output_len=uncompressed_size)


compressions['GZIP'] = gzip_compress_v3
decompressions['GZIP'] = gzip_decompress
compressions['SNAPPY'] = cramjam.snappy.compress_raw
decompressions['SNAPPY'] = cramjam.snappy.decompress_raw

try:
    import lzo
    def lzo_decompress(data, uncompressed_size):
        return lzo.decompress(data)
    compressions['LZO'] = lzo.compress
    decompressions['LZO'] = lzo_decompress
except ImportError:
    pass
compressions['BROTLI'] = cramjam.brotli.compress
decompressions['BROTLI'] = cramjam.brotli.decompress


def lz4_compress(data, **kwargs):
    kwargs['store_size'] = False
    return cramjam.lz4.compress_block(data, **kwargs)


def lz4_decomp(data, size):
    return cramjam.lz4.decompress_block(np.frombuffer(data, 'uint8'), size)


compressions['LZ4'] = lz4_compress
decompressions['LZ4'] = lz4_decomp

# LZ4 is actually LZ4 block, aka "raw", see
# https://github.com/apache/parquet-format/commit/7f06e838cbd1b7dbd722ff2580b9c2525e37fc46
compressions['LZ4_RAW'] = lz4_compress
decompressions['LZ4_RAW'] = lz4_decomp
compressions['ZSTD'] = cramjam.zstd.compress
decompressions['ZSTD'] = cramjam.zstd.decompress
decom_into = {
    "GZIP": cramjam.gzip.decompress_into,
    "SNAPPY": cramjam.snappy.decompress_raw_into,
    "ZSTD": cramjam.zstd.decompress_into,
    "BROTLI": cramjam.brotli.decompress_into
}

compressions = {k.upper(): v for k, v in compressions.items()}
decompressions = {k.upper(): v for k, v in decompressions.items()}

rev_map = {getattr(parquet_thrift.CompressionCodec, key): key for key in
           dir(parquet_thrift.CompressionCodec) if key in
           ['UNCOMPRESSED', 'SNAPPY', 'GZIP', 'LZO', 'BROTLI', 'LZ4', 'ZSTD', 'LZ4_RAW']}


def compress_data(data, compression='gzip'):
    if isinstance(compression, dict):
        algorithm = compression.get('type', 'gzip')
        if isinstance(algorithm, int):
            algorithm = rev_map[compression]
        args = compression.get('args', None)
    else:
        algorithm = compression
        args = None

    if isinstance(algorithm, int):
        algorithm = rev_map[compression]

    if algorithm.upper() not in compressions:
        raise RuntimeError("Compression '%s' not available.  Options: %s" %
                (algorithm, sorted(compressions)))
    if args is None:
        return compressions[algorithm.upper()](data)
    else:
        if not isinstance(args, dict):
            raise ValueError("args dict entry is not a dict")
        return compressions[algorithm.upper()](data, **args)


def decompress_data(data, uncompressed_size, algorithm='gzip'):
    if isinstance(algorithm, int):
        algorithm = rev_map[algorithm]
    if algorithm.upper() not in decompressions:
        raise RuntimeError(
            "Decompression '%s' not available.  Options: %s" %
            (algorithm.upper(), sorted(decompressions))
        )
    if algorithm.upper() in decom_into:
        # ensures writable buffer from cramjam
        x = np.empty(uncompressed_size, dtype='uint8')
        decom_into[algorithm.upper()](np.frombuffer(data, dtype=np.uint8), x)
        return x
    return decompressions[algorithm.upper()](data, uncompressed_size)
