from collections.abc import Sequence
from typing import Final
from typing_extensions import TypeAlias

from Xlib.display import Display
from Xlib.protocol import request, rq
from Xlib.xobject import drawable, resource

_RandRModeInfo13IntSequence: TypeAlias = Sequence[int]

extname: Final = "RANDR"
RRScreenChangeNotify: Final = 0
RRNotify: Final = 1
RRNotify_CrtcChange: Final = 0
RRNotify_OutputChange: Final = 1
RRNotify_OutputProperty: Final = 2
RRScreenChangeNotifyMask: Final = 0x1
RRCrtcChangeNotifyMask: Final = 0x2
RROutputChangeNotifyMask: Final = 0x4
RROutputPropertyNotifyMask: Final = 0x8
SetConfigSuccess: Final = 0
SetConfigInvalidConfigTime: Final = 1
SetConfigInvalidTime: Final = 2
SetConfigFailed: Final = 3
Rotate_0: Final = 1
Rotate_90: Final = 2
Rotate_180: Final = 4
Rotate_270: Final = 8
Reflect_X: Final = 16
Reflect_Y: Final = 32
HSyncPositive: Final = 0x00000001
HSyncNegative: Final = 0x00000002
VSyncPositive: Final = 0x00000004
VSyncNegative: Final = 0x00000008
Interlace: Final = 0x00000010
DoubleScan: Final = 0x00000020
CSync: Final = 0x00000040
CSyncPositive: Final = 0x00000080
CSyncNegative: Final = 0x00000100
HSkewPresent: Final = 0x00000200
BCast: Final = 0x00000400
PixelMultiplex: Final = 0x00000800
DoubleClock: Final = 0x00001000
ClockDivideBy2: Final = 0x00002000
Connected: Final = 0
Disconnected: Final = 1
UnknownConnection: Final = 2
PROPERTY_RANDR_EDID: Final = "EDID"
PROPERTY_SIGNAL_FORMAT: Final = "SignalFormat"
PROPERTY_SIGNAL_PROPERTIES: Final = "SignalProperties"
PROPERTY_CONNECTOR_TYPE: Final = "ConnectorType"
PROPERTY_CONNECTOR_NUMBER: Final = "ConnectorNumber"
PROPERTY_COMPATIBILITY_LIST: Final = "CompatibilityList"
PROPERTY_CLONE_LIST: Final = "CloneList"
SubPixelUnknown: Final = 0
SubPixelHorizontalRGB: Final = 1
SubPixelHorizontalBGR: Final = 2
SubPixelVerticalRGB: Final = 3
SubPixelVerticalBGR: Final = 4
SubPixelNone: Final = 5
BadRROutput: Final = 0
BadRRCrtc: Final = 1
BadRRMode: Final = 2

class BadRROutputError(Exception): ...
class BadRRCrtcError(Exception): ...
class BadRRModeError(Exception): ...

RandR_ScreenSizes: rq.Struct
RandR_ModeInfo: rq.Struct
RandR_Rates: rq.Struct
Render_Transform: rq.Struct
MonitorInfo: rq.Struct

class QueryVersion(rq.ReplyRequest): ...

def query_version(self: Display | resource.Resource) -> QueryVersion: ...

class _1_0SetScreenConfig(rq.ReplyRequest): ...
class SetScreenConfig(rq.ReplyRequest): ...

def set_screen_config(
    self: drawable.Drawable, size_id: int, rotation: int, config_timestamp: int, rate: int = 0, timestamp: int = 0
) -> SetScreenConfig: ...

class SelectInput(rq.Request): ...

def select_input(self: drawable.Window, mask: int) -> SelectInput: ...

class GetScreenInfo(rq.ReplyRequest): ...

def get_screen_info(self: drawable.Window) -> GetScreenInfo: ...

class GetScreenSizeRange(rq.ReplyRequest): ...

def get_screen_size_range(self: drawable.Window) -> GetScreenSizeRange: ...

class SetScreenSize(rq.Request): ...

def set_screen_size(
    self: drawable.Window,
    width: int,
    height: int,
    width_in_millimeters: int | None = None,
    height_in_millimeters: int | None = None,
) -> SetScreenSize: ...

class GetScreenResources(rq.ReplyRequest): ...

def get_screen_resources(self: drawable.Window) -> GetScreenResources: ...

class GetOutputInfo(rq.ReplyRequest): ...

def get_output_info(self: Display | resource.Resource, output: int, config_timestamp: int) -> GetOutputInfo: ...

class ListOutputProperties(rq.ReplyRequest): ...

def list_output_properties(self: Display | resource.Resource, output: int) -> ListOutputProperties: ...

class QueryOutputProperty(rq.ReplyRequest): ...

def query_output_property(self: Display | resource.Resource, output: int, property: int) -> QueryOutputProperty: ...

class ConfigureOutputProperty(rq.Request): ...

def configure_output_property(self: Display | resource.Resource, output: int, property: int) -> ConfigureOutputProperty: ...

class ChangeOutputProperty(rq.Request): ...

def change_output_property(
    self: Display | resource.Resource, output: int, property: int, type: int, mode: int, value: Sequence[float] | Sequence[str]
) -> ChangeOutputProperty: ...

class DeleteOutputProperty(rq.Request): ...

def delete_output_property(self: Display | resource.Resource, output: int, property: int) -> DeleteOutputProperty: ...

class GetOutputProperty(rq.ReplyRequest): ...

def get_output_property(
    self: Display | resource.Resource,
    output: int,
    property: int,
    type: int,
    long_offset: int,
    long_length: int,
    delete: bool = False,
    pending: bool = False,
) -> GetOutputProperty: ...

class CreateMode(rq.ReplyRequest): ...

def create_mode(self: drawable.Window, mode: _RandRModeInfo13IntSequence, name: str) -> CreateMode: ...

class DestroyMode(rq.Request): ...

def destroy_mode(self: Display | resource.Resource, mode: int) -> DestroyMode: ...

class AddOutputMode(rq.Request): ...

def add_output_mode(self: Display | resource.Resource, output: int, mode: int) -> AddOutputMode: ...

class DeleteOutputMode(rq.Request): ...

def delete_output_mode(self: Display | resource.Resource, output: int, mode: int) -> DeleteOutputMode: ...

class GetCrtcInfo(rq.ReplyRequest): ...

def get_crtc_info(self: Display | resource.Resource, crtc: int, config_timestamp: int) -> GetCrtcInfo: ...

class SetCrtcConfig(rq.ReplyRequest): ...

def set_crtc_config(
    self: Display | resource.Resource,
    crtc: int,
    config_timestamp: int,
    x: int,
    y: int,
    mode: int,
    rotation: int,
    outputs: Sequence[int],
    timestamp: int = 0,
) -> SetCrtcConfig: ...

class GetCrtcGammaSize(rq.ReplyRequest): ...

def get_crtc_gamma_size(self: Display | resource.Resource, crtc: int) -> GetCrtcGammaSize: ...

class GetCrtcGamma(rq.ReplyRequest): ...

def get_crtc_gamma(self: Display | resource.Resource, crtc: int) -> GetCrtcGamma: ...

class SetCrtcGamma(rq.Request): ...

def set_crtc_gamma(
    self: Display | resource.Resource, crtc: int, size: int, red: Sequence[int], green: Sequence[int], blue: Sequence[int]
) -> SetCrtcGamma: ...

class GetScreenResourcesCurrent(rq.ReplyRequest): ...

def get_screen_resources_current(self: drawable.Window) -> GetScreenResourcesCurrent: ...

class SetCrtcTransform(rq.Request): ...

def set_crtc_transform(self: Display | resource.Resource, crtc: int, n_bytes_filter: Sequence[int]) -> SetCrtcTransform: ...

class GetCrtcTransform(rq.ReplyRequest): ...

def get_crtc_transform(self: Display | resource.Resource, crtc: int) -> GetCrtcTransform: ...

class GetPanning(rq.ReplyRequest): ...

def get_panning(self: Display | resource.Resource, crtc: int) -> GetPanning: ...

class SetPanning(rq.ReplyRequest): ...

def set_panning(
    self: Display | resource.Resource,
    crtc: int,
    left: int,
    top: int,
    width: int,
    height: int,
    track_left: int,
    track_top: int,
    track_width: int,
    track_height: int,
    border_left: int,
    border_top: int,
    border_width: int,
    border_height: int,
    timestamp: int = 0,
) -> SetPanning: ...

class SetOutputPrimary(rq.Request): ...

def set_output_primary(self: drawable.Window, output: int) -> SetOutputPrimary: ...

class GetOutputPrimary(rq.ReplyRequest): ...

def get_output_primary(self: drawable.Window) -> GetOutputPrimary: ...

class GetMonitors(rq.ReplyRequest): ...

def get_monitors(self: drawable.Window, is_active: bool = True) -> GetMonitors: ...

class SetMonitor(rq.Request): ...

def set_monitor(
    self: drawable.Window, monitor_info: tuple[int, bool, bool, Sequence[int], int, int, int, int, int]
) -> SetMonitor: ...

class DeleteMonitor(rq.Request): ...

def delete_monitor(self: Display | resource.Resource, name: str) -> DeleteMonitor: ...

class ScreenChangeNotify(rq.Event): ...
class CrtcChangeNotify(rq.Event): ...
class OutputChangeNotify(rq.Event): ...
class OutputPropertyNotify(rq.Event): ...

def init(disp: Display, info: request.QueryExtension) -> None: ...
