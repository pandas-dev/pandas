from markdown.extensions import Extension

extensions: list[str]

class ExtraExtension(Extension):
    def __init__(self, **kwargs) -> None: ...

def makeExtension(**kwargs) -> ExtraExtension: ...
