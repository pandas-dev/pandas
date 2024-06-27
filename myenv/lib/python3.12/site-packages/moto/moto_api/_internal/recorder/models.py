import base64
import io
import json
import os
from typing import Any, Optional, Tuple, Union
from urllib.parse import urlparse

import requests
from botocore.awsrequest import AWSPreparedRequest


class Recorder:
    def __init__(self) -> None:
        self._location = str(os.environ.get("MOTO_RECORDER_FILEPATH", "moto_recording"))
        self._os_enabled = bool(os.environ.get("MOTO_ENABLE_RECORDING", False))
        self._user_enabled = self._os_enabled

    def _record_request(self, request: Any, body: Optional[bytes] = None) -> None:
        """
        Record the current request
        """
        if not self._user_enabled:
            return

        if urlparse(request.url).path.startswith("/moto-api/recorder/"):
            return

        entry = {
            "headers": dict(request.headers),
            "method": request.method,
            "url": request.url,
        }

        if body is None:
            if isinstance(request, AWSPreparedRequest):
                body_str, body_encoded = self._encode_body(body=request.body)
            else:
                try:
                    request_body = None
                    request_body_size = int(request.headers["Content-Length"])
                    request_body = request.environ["wsgi.input"].read(request_body_size)
                    body_str, body_encoded = self._encode_body(body=request_body)
                except (AttributeError, KeyError):
                    body_str = ""
                    body_encoded = False
                finally:
                    if request_body is not None:
                        if isinstance(request_body, str):
                            request_body = request_body.encode("utf-8")
                        request.environ["wsgi.input"] = io.BytesIO(request_body)
        else:
            body_str, body_encoded = self._encode_body(body)
        entry.update({"body": body_str, "body_encoded": body_encoded})

        filepath = self._location
        with open(filepath, "a+") as file:
            file.write(json.dumps(entry))
            file.write("\n")

    def _encode_body(self, body: Any) -> Tuple[str, bool]:
        body_encoded = False
        try:
            if isinstance(body, io.BytesIO):
                body = body.getvalue()
            if isinstance(body, bytes):
                body = base64.b64encode(body).decode("ascii")
                body_encoded = True
        except AttributeError:
            body = None
        return body, body_encoded

    def reset_recording(self) -> None:
        """
        Resets the recording. This will erase any requests made previously.
        """
        filepath = self._location
        with open(filepath, "w"):
            pass

    def start_recording(self) -> None:
        """
        Start the recording, and append incoming requests to the log.
        """
        self._user_enabled = True

    def stop_recording(self) -> None:
        self._user_enabled = False

    def upload_recording(self, data: Union[str, bytes]) -> None:
        """
        Replaces the current log. Remember to replay the recording afterwards.
        """
        filepath = self._location
        if isinstance(data, str):
            data = data.encode("utf-8")
        with open(filepath, "bw") as file:
            file.write(data)

    def download_recording(self) -> str:
        """
        Download the current recording. The result can be uploaded afterwards.
        """
        filepath = self._location
        with open(filepath, "r") as file:
            return file.read()

    def replay_recording(self, target_host: Optional[str] = None) -> None:
        """
        Replays the current log, i.e. replay all requests that were made after the recorder was started.
        Download the recording if you want to manually verify the correct requests will be replayed.
        """
        filepath = self._location

        # do not record the replay itself
        old_setting = self._user_enabled
        self._user_enabled = False

        with open(filepath, "r") as file:
            entries = file.readlines()

        for row in entries:
            row_loaded = json.loads(row)
            body = row_loaded.get("body", "{}")
            if row_loaded.get("body_encoded"):
                body = base64.b64decode(body)
            method = row_loaded.get("method")
            url = row_loaded.get("url")
            if target_host is not None:
                parsed_host = urlparse(target_host)
                parsed_url = urlparse(url)
                url = f"{parsed_host.scheme}://{parsed_host.netloc}{parsed_url.path}"
                if parsed_url.query:
                    url = f"{url}?{parsed_url.query}"
            headers = row_loaded.get("headers")
            requests.request(method=method, url=url, headers=headers, data=body)

        # restore the recording setting
        self._user_enabled = old_setting
