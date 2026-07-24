from _typeshed import ConvertibleToFloat, SliceableBuffer, Unused
from collections.abc import Iterable, Sequence
from typing import Final, TypeVar

from Xlib.display import Display
from Xlib.protocol import display, request, rq
from Xlib.xobject import drawable, resource

_T = TypeVar("_T")

extname: Final = "XInputExtension"
PropertyDeleted: Final = 0
PropertyCreated: Final = 1
PropertyModified: Final = 2
NotifyNormal: Final = 0
NotifyGrab: Final = 1
NotifyUngrab: Final = 2
NotifyWhileGrabbed: Final = 3
NotifyPassiveGrab: Final = 4
NotifyPassiveUngrab: Final = 5
NotifyAncestor: Final = 0
NotifyVirtual: Final = 1
NotifyInferior: Final = 2
NotifyNonlinear: Final = 3
NotifyNonlinearVirtual: Final = 4
NotifyPointer: Final = 5
NotifyPointerRoot: Final = 6
NotifyDetailNone: Final = 7
GrabtypeButton: Final = 0
GrabtypeKeycode: Final = 1
GrabtypeEnter: Final = 2
GrabtypeFocusIn: Final = 3
GrabtypeTouchBegin: Final = 4
AnyModifier: Final = 0x80000000
AnyButton: Final = 0
AnyKeycode: Final = 0
AsyncDevice: Final = 0
SyncDevice: Final = 1
ReplayDevice: Final = 2
AsyncPairedDevice: Final = 3
AsyncPair: Final = 4
SyncPair: Final = 5
SlaveSwitch: Final = 1
DeviceChange: Final = 2
MasterAdded: Final = 0x01
MasterRemoved: Final = 0x02
SlaveAdded: Final = 0x04
SlaveRemoved: Final = 0x08
SlaveAttached: Final = 0x10
SlaveDetached: Final = 0x20
DeviceEnabled: Final = 0x40
DeviceDisabled: Final = 0x80
AddMaster: Final = 1
RemoveMaster: Final = 2
AttachSlave: Final = 3
DetachSlave: Final = 4
AttachToMaster: Final = 1
Floating: Final = 2
ModeRelative: Final = 0
ModeAbsolute: Final = 1
MasterPointer: Final = 1
MasterKeyboard: Final = 2
SlavePointer: Final = 3
SlaveKeyboard: Final = 4
FloatingSlave: Final = 5
KeyClass: Final = 0
ButtonClass: Final = 1
ValuatorClass: Final = 2
ScrollClass: Final = 3
TouchClass: Final = 8
KeyRepeat: Final = 0x10000
AllDevices: Final = 0
AllMasterDevices: Final = 1
DeviceChanged: Final = 1
KeyPress: Final = 2
KeyRelease: Final = 3
ButtonPress: Final = 4
ButtonRelease: Final = 5
Motion: Final = 6
Enter: Final = 7
Leave: Final = 8
FocusIn: Final = 9
FocusOut: Final = 10
HierarchyChanged: Final = 11
PropertyEvent: Final = 12
RawKeyPress: Final = 13
RawKeyRelease: Final = 14
RawButtonPress: Final = 15
RawButtonRelease: Final = 16
RawMotion: Final = 17
DeviceChangedMask: Final = 0x00002
KeyPressMask: Final = 0x00004
KeyReleaseMask: Final = 0x00008
ButtonPressMask: Final = 0x00010
ButtonReleaseMask: Final = 0x00020
MotionMask: Final = 0x00040
EnterMask: Final = 0x00080
LeaveMask: Final = 0x00100
FocusInMask: Final = 0x00200
FocusOutMask: Final = 0x00400
HierarchyChangedMask: Final = 0x00800
PropertyEventMask: Final = 0x01000
RawKeyPressMask: Final = 0x02000
RawKeyReleaseMask: Final = 0x04000
RawButtonPressMask: Final = 0x08000
RawButtonReleaseMask: Final = 0x10000
RawMotionMask: Final = 0x20000
GrabModeSync: Final = 0
GrabModeAsync: Final = 1
GrabModeTouch: Final = 2
DEVICEID = rq.Card16
DEVICE = rq.Card16
DEVICEUSE = rq.Card8
PROPERTY_TYPE_FLOAT: Final = "FLOAT"

# ignore[override] because of Liskov substitution principle violations
class FP1616(rq.Int32):
    def check_value(self, value: float) -> int: ...  # type: ignore[override]
    def parse_value(self, value: ConvertibleToFloat, display: Unused) -> float: ...  # type: ignore[override]

class FP3232(rq.ValueField):
    structcode: str
    def check_value(self, value: _T) -> _T: ...  # type: ignore[override]
    def parse_value(self, value: tuple[ConvertibleToFloat, ConvertibleToFloat], display: Unused) -> float: ...  # type: ignore[override]

class XIQueryVersion(rq.ReplyRequest): ...

def query_version(self: Display | resource.Resource) -> XIQueryVersion: ...

class Mask(rq.List):
    def __init__(self, name: str) -> None: ...
    def pack_value(self, val: int | Iterable[int]) -> tuple[bytes, int, None]: ...  # type: ignore[override]

EventMask: rq.Struct

class XISelectEvents(rq.Request): ...

def select_events(self: drawable.Window, event_masks: Sequence[tuple[int, Sequence[int]]]) -> XISelectEvents: ...

AnyInfo: rq.Struct

class ButtonMask:
    def __init__(self, value: int, length: int) -> None: ...
    def __getitem__(self, key: int) -> int: ...
    def __len__(self) -> int: ...

class ButtonState(rq.ValueField):
    structcode: None
    def __init__(self, name: str) -> None: ...
    def parse_binary_value(  # type: ignore[override]  # length: None will error. See: https://github.com/python-xlib/python-xlib/pull/248
        self, data: SliceableBuffer, display: Unused, length: int, fmt: Unused
    ) -> tuple[ButtonMask, SliceableBuffer]: ...

ButtonInfo: rq.Struct
KeyInfo: rq.Struct
ValuatorInfo: rq.Struct
ScrollInfo: rq.Struct
TouchInfo: rq.Struct
INFO_CLASSES: Final[dict[int, rq.Struct]]

class ClassInfoClass:
    structcode: None
    def parse_binary(self, data: SliceableBuffer, display: display.Display | None) -> tuple[rq.DictWrapper, SliceableBuffer]: ...

ClassInfo: ClassInfoClass
DeviceInfo: rq.Struct

class XIQueryDevice(rq.ReplyRequest): ...

def query_device(self: Display | resource.Resource, deviceid: int) -> XIQueryDevice: ...

class XIListProperties(rq.ReplyRequest): ...

def list_device_properties(self: Display | resource.Resource, deviceid: int) -> XIListProperties: ...

class XIGetProperty(rq.ReplyRequest): ...

def get_device_property(
    self: Display | resource.Resource, deviceid: int, property: int, type: int, offset: int, length: int, delete: bool = False
) -> XIGetProperty: ...

class XIChangeProperty(rq.Request): ...

def change_device_property(
    self: Display | resource.Resource, deviceid: int, property: int, type: int, mode: int, value: Sequence[float] | Sequence[str]
) -> XIChangeProperty: ...

class XIDeleteProperty(rq.Request): ...

def delete_device_property(self: Display | resource.Resource, deviceid: int, property: int) -> XIDeleteProperty: ...

class XIGrabDevice(rq.ReplyRequest): ...

def grab_device(
    self: drawable.Window,
    deviceid: int,
    time: int,
    grab_mode: int,
    paired_device_mode: int,
    owner_events: bool,
    event_mask: Sequence[int],
) -> XIGrabDevice: ...

class XIUngrabDevice(rq.Request): ...

def ungrab_device(self: Display | resource.Resource, deviceid: int, time: int) -> XIUngrabDevice: ...

class XIPassiveGrabDevice(rq.ReplyRequest): ...

def passive_grab_device(
    self: drawable.Window,
    deviceid: int,
    time: int,
    detail: int,
    grab_type: int,
    grab_mode: int,
    paired_device_mode: int,
    owner_events: bool,
    event_mask: Sequence[int],
    modifiers: Sequence[int],
) -> XIPassiveGrabDevice: ...
def grab_keycode(
    self: drawable.Window,
    deviceid: int,
    time: int,
    keycode: int,
    grab_mode: int,
    paired_device_mode: int,
    owner_events: bool,
    event_mask: Sequence[int],
    modifiers: Sequence[int],
) -> XIPassiveGrabDevice: ...

class XIPassiveUngrabDevice(rq.Request): ...

def passive_ungrab_device(
    self: drawable.Window, deviceid: int, detail: int, grab_type: int, modifiers: Sequence[int]
) -> XIPassiveUngrabDevice: ...
def ungrab_keycode(self: drawable.Window, deviceid: int, keycode: int, modifiers: Sequence[int]) -> XIPassiveUngrabDevice: ...

HierarchyInfo: rq.Struct
HierarchyEventData: rq.Struct
ModifierInfo: rq.Struct
GroupInfo: rq.Struct
DeviceEventData: rq.Struct
DeviceChangedEventData: rq.Struct
PropertyEventData: rq.Struct

def init(disp: Display, info: request.QueryExtension) -> None: ...
