import os
from datetime import datetime, timezone
from typing import List, Optional, Text, Union


def file_exists(path: Union[bytes, str, "os.PathLike[Text]"]) -> bool:
    """Validates file path as existing local file"""
    return os.path.isfile(path)


def get_file_modtime(path: Union[bytes, str, "os.PathLike[Text]"]) -> Text:
    """Returns ISO formatted file modification time in local system timezone"""
    return (
        datetime.fromtimestamp(os.stat(path).st_mtime, timezone.utc)
        .astimezone()
        .isoformat()
    )


def get_tables_argument_list(table_list: Optional[List[Text]]) -> Optional[List[Text]]:
    """Converts a list of OpenType table string into a Python list or
    return None if the table_list was not defined (i.e., it was not included
    in an option on the command line). Tables that are composed of three
    characters must be right padded with a space."""
    if table_list is None:
        return None
    else:
        return [table.ljust(4) for table in table_list]
