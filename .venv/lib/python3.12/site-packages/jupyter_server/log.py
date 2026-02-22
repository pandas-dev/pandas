"""Log utilities."""

# -----------------------------------------------------------------------------
#  Copyright (c) Jupyter Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file LICENSE, distributed as part of this software.
# -----------------------------------------------------------------------------
import json
from urllib.parse import urlparse, urlunparse

from tornado.log import access_log

from .auth import User
from .prometheus.log_functions import prometheus_log_method

# url params to be scrubbed if seen
# any url param that *contains* one of these
# will be scrubbed from logs
_DEFAULT_SCRUB_PARAM_KEYS = {"token", "auth", "key", "code", "state", "xsrf"}


def _scrub_uri(uri: str, extra_param_keys=None) -> str:
    """scrub auth info from uri"""

    scrub_param_keys = _DEFAULT_SCRUB_PARAM_KEYS.union(set(extra_param_keys or []))

    parsed = urlparse(uri)
    if parsed.query:
        # check for potentially sensitive url params
        # use manual list + split rather than parsing
        # to minimally perturb original
        parts = parsed.query.split("&")
        changed = False
        for i, s in enumerate(parts):
            key, sep, value = s.partition("=")
            for substring in scrub_param_keys:
                if substring in key:
                    parts[i] = f"{key}{sep}[secret]"
                    changed = True
        if changed:
            parsed = parsed._replace(query="&".join(parts))
            return urlunparse(parsed)
    return uri


def log_request(handler, record_prometheus_metrics=True):
    """log a bit more information about each request than tornado's default

    - move static file get success to debug-level (reduces noise)
    - get proxied IP instead of proxy IP
    - log referer for redirect and failed requests
    - log user-agent for failed requests

    if record_prometheus_metrics is true, will record a histogram prometheus
    metric (http_request_duration_seconds) for each request handler
    """
    status = handler.get_status()
    request = handler.request
    try:
        logger = handler.log
    except AttributeError:
        logger = access_log

    extra_param_keys = handler.settings.get("extra_log_scrub_param_keys", [])

    if status < 300 or status == 304:
        # Successes (or 304 FOUND) are debug-level
        log_method = logger.debug
    elif status < 400:
        log_method = logger.info
    elif status < 500:
        log_method = logger.warning
    else:
        log_method = logger.error

    request_time = 1000.0 * handler.request.request_time()
    ns = {
        "status": status,
        "method": request.method,
        "ip": request.remote_ip,
        "uri": _scrub_uri(request.uri, extra_param_keys),
        "request_time": request_time,
    }
    # log username
    # make sure we don't break anything
    # in case mixins cause current_user to not be a User somehow
    try:
        user = handler.current_user
    except Exception:
        user = None
    username = (user.username if isinstance(user, User) else "unknown") if user else ""
    ns["username"] = username

    msg = "{status} {method} {uri} ({username}@{ip}) {request_time:.2f}ms"
    if status >= 400:
        # log bad referrers
        ns["referer"] = _scrub_uri(request.headers.get("Referer", "None"), extra_param_keys)
        msg = msg + " referer={referer}"
    if status >= 500 and status != 502:
        # Log a subset of the headers if it caused an error.
        headers = {}
        for header in ["Host", "Accept", "Referer", "User-Agent"]:
            if header in request.headers:
                headers[header] = request.headers[header]
        log_method(json.dumps(headers, indent=2))
    log_method(msg.format(**ns))
    if record_prometheus_metrics:
        prometheus_log_method(handler)
