from .setopt import option_base

class saveopts(option_base):
    description: str
    def run(self) -> None: ...
