from typing import Final

URL: Final = "url"
METHOD: Final = "method"
USER_AGENT: Final = "user_agent"
CLIENT_IP: Final = "client_ip"
X_FORWARDED_FOR: Final = "x_forwarded_for"
STATUS: Final = "status"
CONTENT_LENGTH: Final = "content_length"
XRAY_HEADER: Final = "X-Amzn-Trace-Id"
ALT_XRAY_HEADER: Final = "HTTP_X_AMZN_TRACE_ID"
request_keys: tuple[str, ...]
response_keys: tuple[str, ...]
