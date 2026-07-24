import sys

from gevent._config import config as config
from gevent._hub_local import get_hub as get_hub
from gevent._hub_primitives import iwait_on_objects as iwait, wait_on_objects as wait
from gevent.greenlet import Greenlet as Greenlet, joinall as joinall, killall as killall
from gevent.hub import (
    GreenletExit as GreenletExit,
    getcurrent as getcurrent,
    idle as idle,
    kill as kill,
    reinit as reinit,
    signal as signal_handler,
    sleep as sleep,
    spawn_raw as spawn_raw,
)
from gevent.timeout import Timeout as Timeout, with_timeout as with_timeout

if sys.platform != "win32":
    from gevent.os import fork

    __all__ = [
        "Greenlet",
        "GreenletExit",
        "Timeout",
        "config",
        "fork",
        "get_hub",
        "getcurrent",
        "getswitchinterval",
        "idle",
        "iwait",
        "joinall",
        "kill",
        "killall",
        "reinit",
        "setswitchinterval",
        "signal_handler",
        "sleep",
        "spawn",
        "spawn_later",
        "spawn_raw",
        "wait",
        "with_timeout",
    ]
else:
    __all__ = [
        "Greenlet",
        "GreenletExit",
        "Timeout",
        "config",
        "get_hub",
        "getcurrent",
        "getswitchinterval",
        "idle",
        "iwait",
        "joinall",
        "kill",
        "killall",
        "reinit",
        "setswitchinterval",
        "signal_handler",
        "sleep",
        "spawn",
        "spawn_later",
        "spawn_raw",
        "wait",
        "with_timeout",
    ]

__version__: str

getswitchinterval = sys.getswitchinterval
setswitchinterval = sys.setswitchinterval
spawn = Greenlet.spawn
spawn_later = Greenlet.spawn_later
