from tkinter import Frame, Grid, Misc, Pack, Place, Scrollbar, Text
from typing import Any, Optional

# The methods from Pack, Place, and Grid are dynamically added over the parent's impls
class ScrolledText(Text):
    frame: Frame
    vbar: Scrollbar
    def __init__(self, master: Optional[Misc] = ..., **kwargs: Any) -> None: ...
