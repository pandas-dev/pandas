import _win32typing

def OleCreate(
    clsid,
    clsid1,
    obCLSID: _win32typing.PyIID,
    obIID: _win32typing.PyIID,
    renderopt,
    obFormatEtc,
    obOleClientSite: _win32typing.PyIOleClientSite,
    obStorage: _win32typing.PyIStorage,
    /,
) -> _win32typing.PyIOleObject: ...
def OleLoadPicture(
    stream: _win32typing.PyIStream, size, runMode, arg: _win32typing.PyIID, arg1: _win32typing.PyIID, /
) -> _win32typing.PyIUnknown: ...
def OleLoadPicturePath(
    url_or_path: str, unk, reserved, clr, arg: _win32typing.PyIID, arg1: _win32typing.PyIID, /
) -> _win32typing.PyIUnknown: ...
def OleSetContainedObject(unk: _win32typing.PyIUnknown, fContained, /) -> None: ...
def OleTranslateAccelerator(frame: _win32typing.PyIOleInPlaceFrame, frame_info, msg: _win32typing.PyMSG, /) -> None: ...

EMBDHLP_CREATENOW: int
EMBDHLP_DELAYCREATE: int
EMBDHLP_INPROC_HANDLER: int
EMBDHLP_INPROC_SERVER: int
OLECLOSE_NOSAVE: int
OLECLOSE_PROMPTSAVE: int
OLECLOSE_SAVEIFDIRTY: int
OLECMDF_ENABLED: int
OLECMDF_LATCHED: int
OLECMDF_NINCHED: int
OLECMDF_SUPPORTED: int
OLECMDTEXTF_NAME: int
OLECMDTEXTF_NONE: int
OLECMDTEXTF_STATUS: int
OLECREATE_LEAVERUNNING: int
OLEIVERB_DISCARDUNDOSTATE: int
OLEIVERB_HIDE: int
OLEIVERB_INPLACEACTIVATE: int
OLEIVERB_OPEN: int
OLEIVERB_PRIMARY: int
OLEIVERB_SHOW: int
OLEIVERB_UIACTIVATE: int
IID_IObjectWithSite: _win32typing.PyIID
IID_IOleClientSite: _win32typing.PyIID
IID_IOleCommandTarget: _win32typing.PyIID
IID_IOleControl: _win32typing.PyIID
IID_IOleControlSite: _win32typing.PyIID
IID_IOleInPlaceActiveObject: _win32typing.PyIID
IID_IOleInPlaceFrame: _win32typing.PyIID
IID_IOleInPlaceObject: _win32typing.PyIID
IID_IOleInPlaceSite: _win32typing.PyIID
IID_IOleInPlaceSiteEx: _win32typing.PyIID
IID_IOleInPlaceSiteWindowless: _win32typing.PyIID
IID_IOleInPlaceUIWindow: _win32typing.PyIID
IID_IOleLink: _win32typing.PyIID
IID_IOleObject: _win32typing.PyIID
IID_ISpecifyPropertyPages: _win32typing.PyIID
IID_IViewObject: _win32typing.PyIID
IID_IViewObject2: _win32typing.PyIID
