import datetime
import json
from json import JSONEncoder
from typing import Any, Optional, Tuple


class _DateTimeEncoder(JSONEncoder):
    def default(self, o):
        if isinstance(o, (datetime.date, datetime.datetime)):
            return o.isoformat()
        else:
            return str(o)


def to_json_str(obj: Any, separators: Optional[Tuple[str, str]] = None) -> str:
    return json.dumps(obj, cls=_DateTimeEncoder, separators=separators)
