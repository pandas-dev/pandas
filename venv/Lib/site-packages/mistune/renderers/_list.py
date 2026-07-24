from typing import TYPE_CHECKING, Any, Dict, Iterable, cast

from ..util import strip_end

if TYPE_CHECKING:
    from ..core import BaseRenderer, BlockState


def render_list(renderer: "BaseRenderer", token: Dict[str, Any], state: "BlockState") -> str:
    attrs = token["attrs"]
    if attrs["ordered"]:
        children = _render_ordered_list(renderer, token, state)
    else:
        children = _render_unordered_list(renderer, token, state)

    text = "".join(children)
    parent = token.get("parent")
    if parent:
        if parent["tight"]:
            return text
        return text + "\n"
    return strip_end(text) + "\n"


def render_list_item(
    renderer: "BaseRenderer",
    item: Dict[str, Any],
    state: "BlockState",
    marker: str = "",
) -> str:
    parent = item.get("parent")
    if not parent:
        parent = {"leading": "- ", "tight": False}

    leading = cast(str, parent["leading"]) + marker
    text = ""
    prev = None
    for tok in item["children"]:
        if tok["type"] == "list":
            tok["parent"] = parent
        elif tok["type"] == "blank_line":
            continue
        tok["prev"] = prev
        prev = tok
        text += renderer.render_token(tok, state)

    lines = text.splitlines()
    text = (lines[0] if lines else "") + "\n"
    prefix = " " * len(leading)
    for line in lines[1:]:
        if line:
            text += prefix + line + "\n"
        else:
            text += "\n"
    return leading + text


def _render_ordered_list(renderer: "BaseRenderer", token: Dict[str, Any], state: "BlockState") -> Iterable[str]:
    attrs = token["attrs"]
    start = attrs.get("start", 1)
    for item in token["children"]:
        leading = str(start) + token["bullet"] + " "
        item["parent"] = {
            "leading": leading,
            "tight": token["tight"],
        }
        try:
            yield renderer.render_token(item, state)
        finally:
            item.pop("parent", None)
            start += 1


def _render_unordered_list(renderer: "BaseRenderer", token: Dict[str, Any], state: "BlockState") -> Iterable[str]:
    parent = {
        "leading": token["bullet"] + " ",
        "tight": token["tight"],
    }
    for item in token["children"]:
        item["parent"] = parent
        try:
            yield renderer.render_token(item, state)
        finally:
            item.pop("parent", None)
