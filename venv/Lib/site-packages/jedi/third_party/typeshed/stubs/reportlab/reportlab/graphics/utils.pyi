class RenderPMError(Exception): ...

def setFont(gs, fontName, fontSize) -> None: ...
def pathNumTrunc(n): ...
def text2Path(
    text,
    x: int = 0,
    y: int = 0,
    fontName="Times-Roman",
    fontSize: int = 1000,
    anchor: str = "start",
    truncate: int = 1,
    pathReverse: int = 0,
    gs=None,
    **kwds,
): ...

# NOTE: This only exists on some render backends
def processGlyph(G, truncate=1, pathReverse=0): ...
def text2PathDescription(text, x=0, y=0, fontName=..., fontSize=1000, anchor="start", truncate=1, pathReverse=0, gs=None): ...

__all__ = ("setFont", "pathNumTrunc", "processGlyph", "text2PathDescription", "text2Path", "RenderPMError")
