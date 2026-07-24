from _typeshed import Incomplete

class Flake8Exception(Exception): ...
class EarlyQuit(Flake8Exception): ...
class ExecutionError(Flake8Exception): ...

class FailedToLoadPlugin(Flake8Exception):
    FORMAT: str
    plugin_name: Incomplete
    original_exception: Incomplete
    def __init__(self, plugin_name: str, exception: Exception) -> None: ...

class PluginRequestedUnknownParameters(Flake8Exception):
    FORMAT: str
    plugin_name: Incomplete
    original_exception: Incomplete
    def __init__(self, plugin_name: str, exception: Exception) -> None: ...

class PluginExecutionFailed(Flake8Exception):
    FORMAT: str
    filename: Incomplete
    plugin_name: Incomplete
    original_exception: Incomplete
    def __init__(self, filename: str, plugin_name: str, exception: Exception) -> None: ...
