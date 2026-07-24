if False:  # MYPY
    from typing import Any, Dict, Optional, List, Union, Optional, Iterator  # NOQA

anchor_attrib = "_yaml_anchor"


class Anchor(object):
    __slots__ = "value", "always_dump"
    attrib = anchor_attrib

    def __init__(self):
        # type: () -> None
        self.value = None
        self.always_dump = False

    def __repr__(self):
        # type: () -> Any
        ad = ", (always dump)" if self.always_dump else ""
        return "Anchor({!r}{})".format(self.value, ad)
