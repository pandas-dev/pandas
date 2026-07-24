from _typeshed import Incomplete, ReadableBuffer
from array import array
from typing import NamedTuple, type_check_only

is64bit: bool

@type_check_only
class _WMNOTIFY(NamedTuple):
    hwndFrom: Incomplete
    idFrom: Incomplete
    code: Incomplete

def UnpackWMNOTIFY(lparam: int) -> _WMNOTIFY: ...
@type_check_only
class _NMITEMACTIVATE(NamedTuple):
    hwndFrom: Incomplete
    idFrom: Incomplete
    code: Incomplete
    iItem: Incomplete
    iSubItem: Incomplete
    uNewState: Incomplete
    uOldState: Incomplete
    uChanged: Incomplete
    actionx: Incomplete
    actiony: Incomplete
    lParam: Incomplete

def UnpackNMITEMACTIVATE(lparam) -> _NMITEMACTIVATE: ...
def PackMENUITEMINFO(
    fType: Incomplete | None = ...,
    fState: Incomplete | None = ...,
    wID: Incomplete | None = ...,
    hSubMenu: Incomplete | None = ...,
    hbmpChecked: Incomplete | None = ...,
    hbmpUnchecked: Incomplete | None = ...,
    dwItemData: Incomplete | None = ...,
    text: Incomplete | None = ...,
    hbmpItem: Incomplete | None = ...,
    dwTypeData: Incomplete | None = ...,
) -> tuple[array[int], list[Incomplete]]: ...
@type_check_only
class _MENUITEMINFO(NamedTuple):
    fType: int | None
    fState: int | None
    wID: int | None
    hSubMenu: int | None
    hbmpChecked: int | None
    hbmpUnchecked: int | None
    dwItemData: int | None
    text: str | None
    hbmpItem: int | None

def UnpackMENUITEMINFO(s: ReadableBuffer) -> _MENUITEMINFO: ...
def EmptyMENUITEMINFO(mask: Incomplete | None = ..., text_buf_size: int = ...) -> tuple[array[int], list[array[int]]]: ...
def PackMENUINFO(
    dwStyle: Incomplete | None = ...,
    cyMax: Incomplete | None = ...,
    hbrBack: Incomplete | None = ...,
    dwContextHelpID: Incomplete | None = ...,
    dwMenuData: Incomplete | None = ...,
    fMask: int = ...,
) -> array[int]: ...
@type_check_only
class _MENUINFO(NamedTuple):
    dwStyle: Incomplete | None
    cyMax: Incomplete | None
    hbrBack: Incomplete | None
    dwContextHelpID: Incomplete | None
    dwMenuData: Incomplete | None

def UnpackMENUINFO(s: ReadableBuffer) -> _MENUINFO: ...
def EmptyMENUINFO(mask: Incomplete | None = ...) -> array[int]: ...
def PackTVINSERTSTRUCT(parent, insertAfter, tvitem) -> tuple[bytes, list[Incomplete]]: ...
def PackTVITEM(hitem, state, stateMask, text, image, selimage, citems, param) -> tuple[array[int], list[Incomplete]]: ...
def EmptyTVITEM(hitem, mask: Incomplete | None = ..., text_buf_size: int = ...) -> tuple[array[int], list[Incomplete]]: ...
@type_check_only
class _TVITEM(NamedTuple):
    item_hItem: Incomplete
    item_state: Incomplete | None
    item_stateMask: Incomplete | None
    text: Incomplete | None
    item_image: Incomplete | None
    item_selimage: Incomplete | None
    item_cChildren: Incomplete | None
    item_param: Incomplete | None

def UnpackTVITEM(buffer: ReadableBuffer) -> _TVITEM: ...
@type_check_only
class _TVNOTIFY(NamedTuple):
    hwndFrom: Incomplete
    id: Incomplete
    code: Incomplete
    action: Incomplete
    item_old: _TVITEM
    item_new: _TVITEM

def UnpackTVNOTIFY(lparam: int) -> _TVNOTIFY: ...
@type_check_only
class _TVDISPINFO(NamedTuple):
    hwndFrom: Incomplete
    id: Incomplete
    code: Incomplete
    item: _TVITEM

def UnpackTVDISPINFO(lparam: int) -> _TVDISPINFO: ...
def PackLVITEM(
    item: Incomplete | None = ...,
    subItem: Incomplete | None = ...,
    state: Incomplete | None = ...,
    stateMask: Incomplete | None = ...,
    text: Incomplete | None = ...,
    image: Incomplete | None = ...,
    param: Incomplete | None = ...,
    indent: Incomplete | None = ...,
) -> tuple[array[int], list[Incomplete]]: ...
@type_check_only
class _LVITEM(NamedTuple):
    item_item: Incomplete
    item_subItem: Incomplete
    item_state: Incomplete | None
    item_stateMask: Incomplete | None
    text: Incomplete | None
    item_image: Incomplete | None
    item_param: Incomplete | None
    item_indent: Incomplete | None

def UnpackLVITEM(buffer: ReadableBuffer) -> _LVITEM: ...
@type_check_only
class _LVDISPINFO(NamedTuple):
    hwndFrom: Incomplete
    id: Incomplete
    code: Incomplete
    item: _LVITEM

def UnpackLVDISPINFO(lparam: int) -> _LVDISPINFO: ...
@type_check_only
class _UnpackLVNOTIFY(NamedTuple):
    hwndFrom: Incomplete
    id: Incomplete
    code: Incomplete
    item: Incomplete
    subitem: Incomplete
    newstate: Incomplete
    oldstate: Incomplete
    changed: Incomplete
    pt: tuple[Incomplete, Incomplete]
    lparam: Incomplete

def UnpackLVNOTIFY(lparam: int) -> _UnpackLVNOTIFY: ...
def EmptyLVITEM(
    item, subitem, mask: Incomplete | None = ..., text_buf_size: int = ...
) -> tuple[array[int], list[Incomplete]]: ...
def PackLVCOLUMN(
    fmt: Incomplete | None = ...,
    cx: Incomplete | None = ...,
    text: Incomplete | None = ...,
    subItem: Incomplete | None = ...,
    image: Incomplete | None = ...,
    order: Incomplete | None = ...,
) -> tuple[array[int], list[Incomplete]]: ...
@type_check_only
class _LVCOLUMN(NamedTuple):
    fmt: Incomplete | None
    cx: Incomplete | None
    text: Incomplete | None
    subItem: Incomplete | None
    image: Incomplete | None
    order: Incomplete | None

def UnpackLVCOLUMN(lparam: ReadableBuffer) -> _LVCOLUMN: ...
def EmptyLVCOLUMN(mask: Incomplete | None = ..., text_buf_size: int = ...) -> tuple[array[int], list[Incomplete]]: ...
def PackLVHITTEST(pt) -> tuple[array[int], None]: ...
@type_check_only
class _LVHITTEST(NamedTuple):
    pt: tuple[Incomplete, Incomplete]
    flags: Incomplete
    item: Incomplete
    subitem: Incomplete

def UnpackLVHITTEST(buf: ReadableBuffer) -> tuple[tuple[Incomplete, Incomplete], Incomplete, Incomplete, Incomplete]: ...
def PackHDITEM(
    cxy: Incomplete | None = ...,
    text: Incomplete | None = ...,
    hbm: Incomplete | None = ...,
    fmt: Incomplete | None = ...,
    param: Incomplete | None = ...,
    image: Incomplete | None = ...,
    order: Incomplete | None = ...,
) -> tuple[array[int], list[Incomplete]]: ...
def PackDEV_BROADCAST(devicetype, rest_fmt, rest_data, extra_data=...) -> bytes: ...
def PackDEV_BROADCAST_HANDLE(handle, hdevnotify: int = ..., guid=..., name_offset: int = ..., data=...) -> bytes: ...
def PackDEV_BROADCAST_VOLUME(unitmask, flags) -> bytes: ...
def PackDEV_BROADCAST_DEVICEINTERFACE(classguid, name: str = ...) -> bytes: ...

class DEV_BROADCAST_INFO:
    devicetype: Incomplete
    def __init__(self, devicetype, **kw) -> None: ...

def UnpackDEV_BROADCAST(lparam: int) -> DEV_BROADCAST_INFO | None: ...
