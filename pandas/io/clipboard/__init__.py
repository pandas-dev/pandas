from functools import partial

from pyclip import copy, paste

# pandas aliases
clipboard_get = partial(paste, text=True)
clipboard_set = copy
