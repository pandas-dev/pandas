import io
from typing import Dict, Optional, Tuple

import multipart


def get_body_from_form_data(
    body: bytes, boundary: str
) -> Tuple[Optional[bytes], Dict[str, str]]:
    body_stream = io.BytesIO(body)
    parser = multipart.MultipartParser(body_stream, boundary=boundary)

    data = None
    headers: Dict[str, str] = {}
    for prt in parser.parts():
        if prt.name == "upload_file":
            headers["key"] = prt.name
            data = prt.file.read()
        else:
            val = prt.file.read()
            if prt.name == "file":
                data = val
            else:
                headers[prt.name] = val.decode("utf-8")
    return data, headers
