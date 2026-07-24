from argparse import ArgumentParser

from ..options.manager import OptionManager

def stage1_arg_parser() -> ArgumentParser: ...

class JobsArgument:
    is_auto: bool
    n_jobs: int
    def __init__(self, arg: str) -> None: ...

def register_default_options(option_manager: OptionManager) -> None: ...
