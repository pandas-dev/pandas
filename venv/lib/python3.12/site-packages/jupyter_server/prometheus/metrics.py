"""
Prometheus metrics exported by Jupyter Server

Read https://prometheus.io/docs/practices/naming/ for naming
conventions for metrics & labels.
"""

from prometheus_client import Gauge, Histogram, Info

from jupyter_server._version import version_info as server_version_info

try:
    from notebook._version import version_info as notebook_version_info
except ImportError:
    notebook_version_info = None


if (
    notebook_version_info is not None  # No notebook package found
    and notebook_version_info < (7,)  # Notebook package found, is version 6
    # Notebook package found, but its version is the same as jupyter_server
    # version. This means some package (looking at you, nbclassic) has shimmed
    # the notebook package to instead be imports from the jupyter_server package.
    # In such cases, notebook.prometheus.metrics is actually *this file*, so
    # trying to import it will cause a circular import. So we don't.
    and notebook_version_info != server_version_info
):
    # Jupyter Notebook v6 also defined these metrics.  Re-defining them results in a ValueError,
    # so we simply re-export them if we are co-existing with the notebook v6 package.
    # See https://github.com/jupyter/jupyter_server/issues/209
    from notebook.prometheus.metrics import (
        HTTP_REQUEST_DURATION_SECONDS,
        KERNEL_CURRENTLY_RUNNING_TOTAL,
        TERMINAL_CURRENTLY_RUNNING_TOTAL,
    )
else:
    HTTP_REQUEST_DURATION_SECONDS = Histogram(
        "http_request_duration_seconds",
        "duration in seconds for all HTTP requests",
        ["method", "handler", "status_code"],
    )

    TERMINAL_CURRENTLY_RUNNING_TOTAL = Gauge(
        "terminal_currently_running_total",
        "counter for how many terminals are running",
    )

    KERNEL_CURRENTLY_RUNNING_TOTAL = Gauge(
        "kernel_currently_running_total",
        "counter for how many kernels are running labeled by type",
        ["type"],
    )

# New prometheus metrics that do not exist in notebook v6 go here
SERVER_INFO = Info("jupyter_server", "Jupyter Server Version information")
SERVER_EXTENSION_INFO = Info(
    "jupyter_server_extension",
    "Jupyter Server Extension Version Information",
    ["name", "version", "enabled"],
)
LAST_ACTIVITY = Gauge(
    "jupyter_server_last_activity_timestamp_seconds",
    "Timestamp of last seen activity on this Jupyter Server",
)
SERVER_STARTED = Gauge(
    "jupyter_server_started_timestamp_seconds", "Timestamp of when this Jupyter Server was started"
)
ACTIVE_DURATION = Gauge(
    "jupyter_server_active_duration_seconds",
    "Number of seconds this Jupyter Server has been active",
)

__all__ = [
    "HTTP_REQUEST_DURATION_SECONDS",
    "TERMINAL_CURRENTLY_RUNNING_TOTAL",
    "KERNEL_CURRENTLY_RUNNING_TOTAL",
    "SERVER_INFO",
]
