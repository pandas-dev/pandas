# Implements https://sw.kovidgoyal.net/kitty/graphics-protocol/

from base64 import b64encode, b64decode
import sys
from typing import Union

def _supports_kitty_graphics() -> bool:
    import platform

    if platform.system() not in ("Darwin", "Linux"):
        return False

    isatty = getattr(sys.stdout, "isatty", None)
    if not callable(isatty) or not isatty():
        return False
    # Hardcoding process names instead of using
    # https://sw.kovidgoyal.net/kitty/graphics-protocol/#querying-support-and-available-transmission-mediums
    # to avoid startup slowdown
    supported_terminals = {
        "ghostty",
        "iTerm2",
        "kitty",
        "konsole",
        "warp",
        "wayst",
        "wezterm-gui",
    }
    import psutil

    try:
        process = psutil.Process()
        while process := process.parent():
            if process.name() in supported_terminals:
                return True
    except (psutil.Error, OSError):
        # Walking the process tree can fail when /proc is mounted with
        # ``hidepid`` on shared multi-user systems (common on HPC clusters):
        # ancestor processes owned by other users are inaccessible and psutil
        # raises AccessDenied. Treat as "unsupported" rather than letting it
        # abort the import of IPython.
        return False
    return False


supports_kitty_graphics = _supports_kitty_graphics()


def png_to_kitty_ansi(png: bytes) -> str:
    if not png.startswith(b"\x89PNG\r\n\x1a\n"):
        raise ValueError
    # This simplicity resembles
    # https://sw.kovidgoyal.net/kitty/graphics-protocol/#a-minimal-example
    # but if we need tmux support, we can switch to Unicode like
    # https://github.com/hzeller/timg/blob/main/src/kitty-canvas.cc
    result = ["\033_Ga=T,f=100,", "m=1;"]
    encoded = b64encode(png)
    for i in range(0, len(encoded), 4096):
        result.append(encoded[i : i + 4096].decode("ascii"))
        result.append("\033\\")
        result.append("\033_G")
        result.append("m=1;")
    del result[-2:]
    result[-3] = "m=0;"
    return "".join(result)


def kitty_png_render(png: Union[bytes, str], _md_dict: object) -> None:
    if isinstance(png, str):
        png = png_to_kitty_ansi(b64decode(png))
    else:
        png = png_to_kitty_ansi(png)
    print(png)


display_formatter_default_active_types = [
    "text/plain",
    *(["image/png"] if supports_kitty_graphics else []),
]

terminal_default_mime_renderers = {
    "image/png": kitty_png_render,
}
