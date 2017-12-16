""" s3 support for remote file interactivity """
from pandas import compat
try:
    import hdfs3
except:
    raise ImportError("The hdfs3 library is required to handle hdfs files")


def get_filepath_or_buffer(filepath_or_buffer, encoding=None,
                           compression=None):
    parsed_url = compat.parse_url(filepath_or_buffer)
    if ":" in parsed_url.netloc:
        host, port = parsed_url.netloc.rsplit(":", 1)
        try:
            port = int(port)
            fs = hdfs3.HDFileSystem(host=host, port=port)
        except ValueError:
            pass
    else:
        fs = hdfs3.HDFileSystem()
    filepath_or_buffer = fs.open(parsed_url.path)
    return filepath_or_buffer, None, compression
