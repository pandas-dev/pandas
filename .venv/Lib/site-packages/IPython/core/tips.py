from datetime import datetime
import os
from random import choice
from typing import Any

_tips: Any = {
    # (month, day)
    "every_year": {
        (1, 1): "Happy new year!",
        # European time:
        # [2/8/25, 23:28:24] Fernando Perez: Hi! Yes, this was my first public
        # announcement:
        # [2/8/25, 23:28:25] Fernando Perez:
        # https://mail.python.org/pipermail/python-list/2001-December/093408.html
        # [2/8/25, 23:28:55] Fernando Perez: All that started two months earlier
        # - in October 2001 I read this article:
        # [2/8/25, 23:28:55] Fernando Perez:
        # https://web.archive.org/web/20011202000624/http://www.onlamp.com/pub/a/python/2001/10/11/pythonnews.html
        # [2/8/25, 23:29:05] Fernando Perez: Which is also archived here:
        # [2/8/25, 23:29:05] Fernando Perez:
        # https://docstore.mik.ua/orelly/weblinux2/orn/interactive_python.html
        # [2/8/25, 23:29:48] Fernando Perez: That article put me on the trail of
        # LazyPython, and I also found out (can’t remember where) about IPP
        # (Interactive Python Prompt), another similar projecd by Janko Hauser.
        # [2/8/25, 23:30:20] Fernando Perez: A also read an article around that
        # time, about new features in Python 2.0, which spoke of how
        # sys.displayhook could be programmed to call any object you wanted.
        # [2/8/25, 23:31:01] Fernando Perez: Those things together gave me the
        # idea of implementing a stateful object, the “IPython Prompt”, that
        # could store a  cache of old results, Mathematica-style, and that could
        # have all kinds of other useful tricks up its sleeve.
        # [2/8/25, 23:31:17] Fernando Perez: I can’t remember if I did the
        # prompt stuff first and then read about LazyPython/IPP.
        # [2/8/25, 23:31:41] Fernando Perez: I do know that, implementation
        # wise, at first I just did the tiny IPython 0.0.1 that I posted much
        # later on github as a gist:
        # [2/8/25, 23:31:55] Fernando Perez:
        # https://gist.github.com/fperez/1579699
        # [2/8/25, 23:32:03] Fernando Perez: But I only shared that publicly
        # much later.
        # [2/8/25, 23:33:19] Fernando Perez: I did find out about IPP/LazyPython
        # sometime in October, contacted Janko and Nathan who told me to use
        # their code at will but said they were busy with other things, and got
        # cranking for a few mad weeks on what became the IPython 0.2.0 that I
        # posted about in that comp.lang.python thread of December 2001.
        # [2/8/25, 23:33:52] Fernando Perez: That period from Oct 11 to Dec 9
        # 2001 was maniacal coding, with very little sleep 🙂
        (
            10,
            11,
        ): "IPython first line of code was written {} years ago by Fernando Pérez".format(
            datetime.now().year - 2001
        ),
        (
            12,
            9,
        ): "IPython 0.0.2 was announced {} years ago: https://mail.python.org/pipermail/python-list/2001-December/093408.html".format(
            datetime.now().year - 2001
        ),
        (
            3,
            8,
        ): "Today is international Women's Day: https://www.internationalwomensday.com/",
        (
            3,
            31,
        ): "Happy trans day of visibility! You are valid. You matter. https://en.wikipedia.org/wiki/International_Transgender_Day_of_Visibility",
    },
    "random": [
        "Use `F2` or %edit with no arguments to open an empty editor with a temporary file.",
        "Run your doctests from within IPython for development and debugging. The special %doctest_mode command toggles a mode where the prompt, output and exceptions display matches as closely as possible that of the default Python interpreter.",
        "You can use `files = !ls *.png`",
        "Use the IPython.lib.demo.Demo class to load any Python script as an interactive demo.",
        "Put a ';' at the end of a line to suppress the printing of output.",
        "You can use Ctrl-O to force a new line in terminal IPython",
        "Use `object?` to see the help on `object`, `object??` to view it's source",
        "`?` alone on a line will brings up IPython's help",
        "You can use `%hist` to view history, see the options with `%history?`",
        "You can change the editing mode of IPython to behave more like vi, or emacs.",
        "IPython 9.0+ have hooks to integrate AI/LLM completions.",
        "Use `%timeit` or `%%timeit`, and the  `-r`, `-n`, and `-o` options to easily profile your code.",
        "Use `ipython --help-all | less` to view all the IPython configurations options.",
        "Use `--theme`, or the `%colors` magic to change ipython themes and colors.",
        "The `%timeit` magic has an `-o` flag, which returns the results, making it easy to plot. See `%timeit?`.",
    ],
}

if os.name == "nt":
    _tips["random"].extend(
        [
            "We can't show you all tips on windows as sometime Unicode characters crash windows console, please help us debug it."
        ]
    )
    # unicode may crash windows console, so we filter out tips with non-ASCII characters.
    _tips["every_year"] = {
        k: v
        for k, v in _tips["every_year"].items()
        if all(ord(char) < 128 for char in v)
    }
else:
    _tips["random"].extend(
        [
            "You can use latex or unicode completion, `\\alpha<tab>` will insert the α symbol.",
            "You can find how to type a latex symbol by back completing it `\\θ<tab>` will expand to `\\theta`.",
            "You can find how to type a unicode symbol by back completing it `\\Ⅷ<tab>` will expand to `\\ROMAN NUMERAL EIGHT`.",
            "IPython support combining unicode identifiers, F\\vec<tab> will become F⃗, usefull for physics equations. Play with \\dot \\ddot and others.",
        ]
    )


def pick_tip() -> str:
    current_date = datetime.now()
    month, day = current_date.month, current_date.day

    if (month, day) in _tips["every_year"]:
        return _tips["every_year"][(month, day)]

    return choice(_tips["random"])
