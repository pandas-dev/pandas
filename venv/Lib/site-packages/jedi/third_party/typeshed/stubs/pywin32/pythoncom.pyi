from _typeshed import Incomplete, Unused
from abc import abstractmethod
from collections.abc import Sequence
from typing import ClassVar, SupportsInt, overload
from typing_extensions import TypeAlias, deprecated, disjoint_base

import _win32typing
from win32.lib.pywintypes import TimeType, com_error as com_error

error: TypeAlias = com_error  # noqa: Y042

class internal_error(Exception): ...

@disjoint_base
class com_record:
    @abstractmethod
    def __init__(self, /, *args, **kwargs) -> None: ...
    TLBID: ClassVar[str]
    MJVER: ClassVar[int]
    MNVER: ClassVar[int]
    LCID: ClassVar[int]
    GUID: ClassVar[str]

def CoCreateFreeThreadedMarshaler(unk: _win32typing.PyIUnknown, /) -> _win32typing.PyIUnknown: ...
def CoCreateInstanceEx(
    clsid: _win32typing.PyIID,
    unkOuter: _win32typing.PyIUnknown,
    context,
    serverInfo: tuple[Incomplete, Incomplete, Incomplete, Incomplete],
    iids: list[_win32typing.PyIID],
    /,
) -> _win32typing.PyIUnknown: ...
def CoCreateInstance(
    clsid: _win32typing.PyIID, unkOuter: _win32typing.PyIUnknown | None, context: int, iid: _win32typing.PyIID, /
) -> _win32typing.PyIUnknown: ...
def CoFreeUnusedLibraries() -> None: ...
def CoInitialize() -> None: ...
def CoInitializeEx(flags, /) -> None: ...
def CoInitializeSecurity(
    sd: _win32typing.PySECURITY_DESCRIPTOR, authSvc, reserved1, authnLevel, impLevel, authInfo, capabilities, reserved2, /
) -> None: ...
def CoGetInterfaceAndReleaseStream(stream: _win32typing.PyIStream, iid: _win32typing.PyIID, /) -> _win32typing.PyIUnknown: ...
def CoMarshalInterThreadInterfaceInStream(iid: _win32typing.PyIID, unk: _win32typing.PyIUnknown, /) -> _win32typing.PyIStream: ...
def CoMarshalInterface(
    Stm: _win32typing.PyIStream, riid: _win32typing.PyIID, Unk: _win32typing.PyIUnknown, DestContext, flags, /
) -> None: ...
def CoUnmarshalInterface(Stm: _win32typing.PyIStream, riid: _win32typing.PyIID, /): ...
def CoReleaseMarshalData(Stm: _win32typing.PyIStream, /) -> None: ...
def CoGetObject(name: str, iid: _win32typing.PyIID, bindOpts: Incomplete | None = ..., /) -> _win32typing.PyIUnknown: ...
def CoUninitialize() -> None: ...
def CoRegisterClassObject(iid: _win32typing.PyIID, factory: _win32typing.PyIUnknown, context, flags, /) -> int: ...
def CoResumeClassObjects() -> None: ...
def CoRevokeClassObject(reg: int, /) -> None: ...
def CoTreatAsClass(clsidold: _win32typing.PyIID, clsidnew: _win32typing.PyIID, /) -> None: ...
def CoWaitForMultipleHandles(Flags, Timeout, Handles: list[int], /): ...
def Connect(cls, /) -> _win32typing.PyIDispatch: ...
def connect(cls, /) -> _win32typing.PyIDispatch: ...
def CreateGuid() -> _win32typing.PyIID: ...
def CreateBindCtx() -> _win32typing.PyIBindCtx: ...
def CreateFileMoniker(filename: str, /) -> _win32typing.PyIMoniker: ...
def CreateItemMoniker(delim: str, item: str, /) -> _win32typing.PyIMoniker: ...
def CreatePointerMoniker(IUnknown: _win32typing.PyIUnknown, /) -> _win32typing.PyIMoniker: ...
def CreateURLMonikerEx(Context, URL, Flags: int = ..., /): ...
def CreateTypeLib(): ...
def CreateTypeLib2(): ...
def CreateStreamOnHGlobal(hGlobal: int | None = ..., DeleteOnRelease: bool = ..., /) -> _win32typing.PyIStream: ...
def CreateILockBytesOnHGlobal(hGlobal: int | None = ..., DeleteOnRelease: bool = ..., /) -> _win32typing.PyILockBytes: ...
def EnableQuitMessage(threadId, /) -> None: ...
def FUNCDESC() -> _win32typing.FUNCDESC: ...
def GetActiveObject(cls, /) -> _win32typing.PyIUnknown: ...
def GetClassFile(fileName, /) -> _win32typing.PyIID: ...
def GetFacilityString(scode, /) -> str: ...
def GetRecordFromGuids(
    iid: str | _win32typing.PyIID,
    verMajor: int,
    verMinor: int,
    lcid: int,
    infoIID: str | _win32typing.PyIID,
    data: Incomplete | None = ...,
    /,
): ...
def GetRecordFromTypeInfo(TypeInfo: _win32typing.PyITypeInfo, /): ...
def GetRunningObjectTable(reserved: int = ..., /) -> _win32typing.PyIRunningObjectTable: ...
def GetScodeString(scode, /) -> str: ...
def GetScodeRangeString(scode, /) -> str: ...
def GetSeverityString(scode, /) -> str: ...
def IsGatewayRegistered(iid: _win32typing.PyIID | None, /) -> int: ...
def LoadRegTypeLib(iid: _win32typing.PyIID, versionMajor, versionMinor, lcid, /) -> _win32typing.PyITypeLib: ...
def LoadTypeLib(libFileName: str, /) -> _win32typing.PyITypeLib: ...
def MakePyFactory(iid: _win32typing.PyIID, /) -> _win32typing.PyIClassFactory: ...
@deprecated("Use pywintypes.IID() instead.")
def MakeIID(iidString: str, is_bytes: bool = ..., /) -> _win32typing.PyIID: ...
@deprecated("Use pywintypes.Time() instead.")
def MakeTime(timeRepr: SupportsInt | Sequence[SupportsInt] | TimeType, /) -> TimeType: ...
def MkParseDisplayName(
    displayName: str, bindCtx: _win32typing.PyIBindCtx | None = ..., /
) -> tuple[_win32typing.PyIMoniker, Incomplete, _win32typing.PyIBindCtx]: ...
def new(iid: _win32typing.PyIID | str, /): ...
def New(cls, /) -> _win32typing.PyIDispatch: ...
def ObjectFromAddress(address, iid: _win32typing.PyIID, /) -> _win32typing.PyIUnknown: ...
def ObjectFromLresult(lresult, iid: _win32typing.PyIID, wparm, /) -> _win32typing.PyIUnknown: ...
def OleInitialize() -> None: ...
def OleGetClipboard() -> _win32typing.PyIDataObject: ...
def OleFlushClipboard() -> None: ...
def OleIsCurrentClipboard(dataObj: _win32typing.PyIDataObject, /): ...
def OleSetClipboard(dataObj: _win32typing.PyIDataObject, /) -> None: ...
def OleLoadFromStream(stream: _win32typing.PyIStream, iid: _win32typing.PyIID, /) -> None: ...
def OleSaveToStream(persist: _win32typing.PyIPersistStream, stream: _win32typing.PyIStream, /) -> None: ...
def OleLoad(storage: _win32typing.PyIStorage, iid: _win32typing.PyIID, site: _win32typing.PyIOleClientSite, /) -> None: ...
def ProgIDFromCLSID(clsid, /) -> str: ...
def PumpWaitingMessages(firstMessage: int = ..., lastMessage: int = ..., /) -> int: ...
def PumpMessages() -> None: ...
def QueryPathOfRegTypeLib(iid: _win32typing.PyIID, versionMajor, versionMinor, lcid, /) -> str: ...
def ReadClassStg(storage: _win32typing.PyIStorage, /) -> _win32typing.PyIID: ...
def ReadClassStm(Stm: _win32typing.PyIStream, /) -> _win32typing.PyIID: ...
def RegisterTypeLib(typelib: _win32typing.PyITypeLib, fullPath: str, lcid, helpDir: str | None = ..., /) -> None: ...
def UnRegisterTypeLib(iid: _win32typing.PyIID, versionMajor, versionMinor, lcid, syskind, /) -> str: ...
def RegisterActiveObject(obUnknown: _win32typing.PyIUnknown, clsid: _win32typing.PyIID, flags, /): ...
def RevokeActiveObject(handle, /) -> None: ...
def RegisterDragDrop(hwnd: int, dropTarget: _win32typing.PyIDropTarget, /) -> None: ...
def RevokeDragDrop(hwnd: int, /) -> None: ...
def DoDragDrop() -> None: ...
def StgCreateDocfile(name: str | None, mode: int, reserved: int = ..., /) -> _win32typing.PyIStorage: ...
def StgCreateDocfileOnILockBytes(lockBytes: _win32typing.PyILockBytes, mode, reserved=..., /) -> _win32typing.PyIStorage: ...
def StgOpenStorageOnILockBytes(
    lockBytes: _win32typing.PyILockBytes,
    stgPriority: _win32typing.PyIStorage,
    mode,
    snbExclude: Incomplete | None = ...,
    reserved: int = ...,
    /,
) -> _win32typing.PyIStorage: ...
def StgIsStorageFile(name: str, /): ...
def STGMEDIUM() -> _win32typing.PySTGMEDIUM: ...
@overload
def StgOpenStorage(
    name: str | None, other: _win32typing.PyIStorage, mode: int, snbExclude: Unused = ..., reserved: int = ..., /
) -> _win32typing.PyIStorage: ...
@overload
def StgOpenStorage(
    name: str, other: _win32typing.PyIStorage | None, mode: int, snbExclude: Unused = ..., reserved: int = ..., /
) -> _win32typing.PyIStorage: ...
def StgOpenStorageEx(
    Name: str, Mode: int, stgfmt: int, Attrs: int, riid: _win32typing.PyIID, StgOptions: Incomplete | None = ...
) -> _win32typing.PyIStorage: ...
def StgCreateStorageEx(
    Name: str,
    Mode: int,
    stgfmt: int,
    Attrs: int,
    riid: _win32typing.PyIID,
    StgOptions: Incomplete | None = ...,
    SecurityDescriptor: _win32typing.PySECURITY_DESCRIPTOR | None = ...,
) -> _win32typing.PyIStorage: ...
def TYPEATTR() -> _win32typing.TYPEATTR: ...
def VARDESC() -> _win32typing.VARDESC: ...
def WrapObject(ob, gatewayIID: _win32typing.PyIID, interfaceIID: _win32typing.PyIID, /) -> _win32typing.PyIUnknown: ...
def WriteClassStg(storage: _win32typing.PyIStorage, iid: _win32typing.PyIID, /) -> None: ...
def WriteClassStm(Stm: _win32typing.PyIStream, clsid: _win32typing.PyIID, /) -> None: ...
def UnwrapObject(ob: _win32typing.PyIUnknown, /) -> _win32typing.PyIDispatch: ...
def FmtIdToPropStgName(fmtid: _win32typing.PyIID, /): ...
def PropStgNameToFmtId(Name: str, /) -> _win32typing.PyIID: ...
def CoGetCallContext(riid: _win32typing.PyIID, /) -> _win32typing.PyIServerSecurity: ...
def CoGetObjectContext(riid: _win32typing.PyIID, /) -> _win32typing.PyIContext: ...
def CoGetCancelObject(riid: _win32typing.PyIID, ThreadID: int = ..., /) -> _win32typing.PyICancelMethodCalls: ...
def CoSetCancelObject(Unk: _win32typing.PyIUnknown, /) -> None: ...
def CoEnableCallCancellation() -> None: ...
def CoDisableCallCancellation() -> None: ...

ACTIVEOBJECT_STRONG: int
ACTIVEOBJECT_WEAK: int
ArgNotFound: _win32typing.ArgNotFound
CLSCTX_ALL: int
CLSCTX_INPROC: int
CLSCTX_INPROC_HANDLER: int
CLSCTX_INPROC_SERVER: int
CLSCTX_LOCAL_SERVER: int
CLSCTX_REMOTE_SERVER: int
CLSCTX_SERVER: int
CLSID_DCOMAccessControl: _win32typing.PyIID
CLSID_StdComponentCategoriesMgr: _win32typing.PyIID
CLSID_StdGlobalInterfaceTable: _win32typing.PyIID
COINIT_APARTMENTTHREADED: int
COINIT_DISABLE_OLE1DDE: int
COINIT_MULTITHREADED: int
COINIT_SPEED_OVER_MEMORY: int
COWAIT_ALERTABLE: int
COWAIT_WAITALL: int
DATADIR_GET: int
DATADIR_SET: int
DESCKIND_FUNCDESC: int
DESCKIND_VARDESC: int
DISPATCH_METHOD: int
DISPATCH_PROPERTYGET: int
DISPATCH_PROPERTYPUT: int
DISPATCH_PROPERTYPUTREF: int
DISPID_COLLECT: int
DISPID_CONSTRUCTOR: int
DISPID_DESTRUCTOR: int
DISPID_EVALUATE: int
DISPID_NEWENUM: int
DISPID_PROPERTYPUT: int
DISPID_STARTENUM: int
DISPID_THIS: int
DISPID_UNKNOWN: int
DISPID_VALUE: int
DVASPECT_CONTENT: int
DVASPECT_DOCPRINT: int
DVASPECT_ICON: int
DVASPECT_THUMBNAIL: int
EOAC_ACCESS_CONTROL: int
EOAC_ANY_AUTHORITY: int
EOAC_APPID: int
EOAC_AUTO_IMPERSONATE: int
EOAC_DEFAULT: int
EOAC_DISABLE_AAA: int
EOAC_DYNAMIC: int
EOAC_DYNAMIC_CLOAKING: int
EOAC_MAKE_FULLSIC: int
EOAC_MUTUAL_AUTH: int
EOAC_NONE: int
EOAC_NO_CUSTOM_MARSHAL: int
EOAC_REQUIRE_FULLSIC: int
EOAC_SECURE_REFS: int
EOAC_STATIC_CLOAKING: int
EXTCONN_CALLABLE: int
EXTCONN_STRONG: int
EXTCONN_WEAK: int
Empty: _win32typing.PyOleEmpty
FMTID_DocSummaryInformation: _win32typing.PyIID
FMTID_SummaryInformation: _win32typing.PyIID
FMTID_UserDefinedProperties: _win32typing.PyIID
FUNCFLAG_FBINDABLE: int
FUNCFLAG_FDEFAULTBIND: int
FUNCFLAG_FDISPLAYBIND: int
FUNCFLAG_FHIDDEN: int
FUNCFLAG_FREQUESTEDIT: int
FUNCFLAG_FRESTRICTED: int
FUNCFLAG_FSOURCE: int
FUNCFLAG_FUSESGETLASTERROR: int
FUNC_DISPATCH: int
FUNC_NONVIRTUAL: int
FUNC_PUREVIRTUAL: int
FUNC_STATIC: int
FUNC_VIRTUAL: int
IDLFLAG_FIN: int
IDLFLAG_FLCID: int
IDLFLAG_FOUT: int
IDLFLAG_FRETVAL: int
IDLFLAG_NONE: int
IID_IBindCtx: _win32typing.PyIID
IID_ICancelMethodCalls: _win32typing.PyIID
IID_ICatInformation: _win32typing.PyIID
IID_ICatRegister: _win32typing.PyIID
IID_IClassFactory: _win32typing.PyIID
IID_IClientSecurity: _win32typing.PyIID
IID_IConnectionPoint: _win32typing.PyIID
IID_IConnectionPointContainer: _win32typing.PyIID
IID_IContext: _win32typing.PyIID
IID_ICreateTypeInfo: _win32typing.PyIID
IID_ICreateTypeLib: _win32typing.PyIID
IID_ICreateTypeLib2: _win32typing.PyIID
IID_IDataObject: _win32typing.PyIID
IID_IDispatch: _win32typing.PyIID
IID_IDispatchEx: _win32typing.PyIID
IID_IDropSource: _win32typing.PyIID
IID_IDropTarget: _win32typing.PyIID
IID_IEnumCATEGORYINFO: _win32typing.PyIID
IID_IEnumConnectionPoints: _win32typing.PyIID
IID_IEnumConnections: _win32typing.PyIID
IID_IEnumContextProps: _win32typing.PyIID
IID_IEnumFORMATETC: _win32typing.PyIID
IID_IEnumGUID: _win32typing.PyIID
IID_IEnumMoniker: _win32typing.PyIID
IID_IEnumSTATPROPSETSTG: _win32typing.PyIID
IID_IEnumSTATPROPSTG: _win32typing.PyIID
IID_IEnumSTATSTG: _win32typing.PyIID
IID_IEnumString: _win32typing.PyIID
IID_IEnumVARIANT: _win32typing.PyIID
IID_IErrorLog: _win32typing.PyIID
IID_IExternalConnection: _win32typing.PyIID
IID_IGlobalInterfaceTable: _win32typing.PyIID
IID_ILockBytes: _win32typing.PyIID
IID_IMarshal: _win32typing.PyIID
IID_IMoniker: _win32typing.PyIID
IID_IOleWindow: _win32typing.PyIID
IID_IPersist: _win32typing.PyIID
IID_IPersistFile: _win32typing.PyIID
IID_IPersistPropertyBag: _win32typing.PyIID
IID_IPersistStorage: _win32typing.PyIID
IID_IPersistStream: _win32typing.PyIID
IID_IPersistStreamInit: _win32typing.PyIID
IID_IPropertyBag: _win32typing.PyIID
IID_IPropertySetStorage: _win32typing.PyIID
IID_IPropertyStorage: _win32typing.PyIID
IID_IProvideClassInfo: _win32typing.PyIID
IID_IProvideClassInfo2: _win32typing.PyIID
IID_IRunningObjectTable: _win32typing.PyIID
IID_IServerSecurity: _win32typing.PyIID
IID_IServiceProvider: _win32typing.PyIID
IID_IStdMarshalInfo: _win32typing.PyIID
IID_IStorage: _win32typing.PyIID
IID_IStream: _win32typing.PyIID
IID_ITypeComp: _win32typing.PyIID
IID_ITypeInfo: _win32typing.PyIID
IID_ITypeLib: _win32typing.PyIID
IID_IUnknown: _win32typing.PyIID
IID_NULL: _win32typing.PyIID
IID_StdOle: _win32typing.PyIID
IMPLTYPEFLAG_FDEFAULT: int
IMPLTYPEFLAG_FRESTRICTED: int
IMPLTYPEFLAG_FSOURCE: int
INVOKE_FUNC: int
INVOKE_PROPERTYGET: int
INVOKE_PROPERTYPUT: int
INVOKE_PROPERTYPUTREF: int
InterfaceNames: dict[str, _win32typing.PyIID]
MKSYS_ANTIMONIKER: int
MKSYS_CLASSMONIKER: int
MKSYS_FILEMONIKER: int
MKSYS_GENERICCOMPOSITE: int
MKSYS_ITEMMONIKER: int
MKSYS_NONE: int
MKSYS_POINTERMONIKER: int
MSHCTX_DIFFERENTMACHINE: int
MSHCTX_INPROC: int
MSHCTX_LOCAL: int
MSHCTX_NOSHAREDMEM: int
MSHLFLAGS_NOPING: int
MSHLFLAGS_NORMAL: int
MSHLFLAGS_TABLESTRONG: int
MSHLFLAGS_TABLEWEAK: int
Missing: _win32typing.PyOleMissing
Nothing: _win32typing.PyOleNothing
PARAMFLAG_FHASDEFAULT: int
PARAMFLAG_FIN: int
PARAMFLAG_FLCID: int
PARAMFLAG_FOPT: int
PARAMFLAG_FOUT: int
PARAMFLAG_FRETVAL: int
PARAMFLAG_NONE: int
REGCLS_MULTIPLEUSE: int
REGCLS_MULTI_SEPARATE: int
REGCLS_SINGLEUSE: int
REGCLS_SUSPENDED: int
ROTFLAGS_ALLOWANYCLIENT: int
ROTFLAGS_REGISTRATIONKEEPSALIVE: int
RPC_C_AUTHN_DCE_PRIVATE: int
RPC_C_AUTHN_DCE_PUBLIC: int
RPC_C_AUTHN_DEC_PUBLIC: int
RPC_C_AUTHN_DEFAULT: int
RPC_C_AUTHN_DPA: int
RPC_C_AUTHN_GSS_KERBEROS: int
RPC_C_AUTHN_GSS_NEGOTIATE: int
RPC_C_AUTHN_GSS_SCHANNEL: int
RPC_C_AUTHN_LEVEL_CALL: int
RPC_C_AUTHN_LEVEL_CONNECT: int
RPC_C_AUTHN_LEVEL_DEFAULT: int
RPC_C_AUTHN_LEVEL_NONE: int
RPC_C_AUTHN_LEVEL_PKT: int
RPC_C_AUTHN_LEVEL_PKT_INTEGRITY: int
RPC_C_AUTHN_LEVEL_PKT_PRIVACY: int
RPC_C_AUTHN_MQ: int
RPC_C_AUTHN_MSN: int
RPC_C_AUTHN_NONE: int
RPC_C_AUTHN_WINNT: int
RPC_C_AUTHZ_DCE: int
RPC_C_AUTHZ_DEFAULT: int
RPC_C_AUTHZ_NAME: int
RPC_C_AUTHZ_NONE: int
RPC_C_IMP_LEVEL_ANONYMOUS: int
RPC_C_IMP_LEVEL_DEFAULT: int
RPC_C_IMP_LEVEL_DELEGATE: int
RPC_C_IMP_LEVEL_IDENTIFY: int
RPC_C_IMP_LEVEL_IMPERSONATE: int
STDOLE2_LCID: int
STDOLE2_MAJORVERNUM: int
STDOLE2_MINORVERNUM: int
STDOLE_LCID: int
STDOLE_MAJORVERNUM: int
STDOLE_MINORVERNUM: int
STREAM_SEEK_CUR: int
STREAM_SEEK_END: int
STREAM_SEEK_SET: int
SYS_MAC: int
SYS_WIN16: int
SYS_WIN32: int
ServerInterfaces: dict[_win32typing.PyIID, bytes]
TKIND_ALIAS: int
TKIND_COCLASS: int
TKIND_DISPATCH: int
TKIND_ENUM: int
TKIND_INTERFACE: int
TKIND_MODULE: int
TKIND_RECORD: int
TKIND_UNION: int
TYMED_ENHMF: int
TYMED_FILE: int
TYMED_GDI: int
TYMED_HGLOBAL: int
TYMED_ISTORAGE: int
TYMED_ISTREAM: int
TYMED_MFPICT: int
TYMED_NULL: int
TYPEFLAG_FAGGREGATABLE: int
TYPEFLAG_FAPPOBJECT: int
TYPEFLAG_FCANCREATE: int
TYPEFLAG_FCONTROL: int
TYPEFLAG_FDISPATCHABLE: int
TYPEFLAG_FDUAL: int
TYPEFLAG_FHIDDEN: int
TYPEFLAG_FLICENSED: int
TYPEFLAG_FNONEXTENSIBLE: int
TYPEFLAG_FOLEAUTOMATION: int
TYPEFLAG_FPREDECLID: int
TYPEFLAG_FREPLACEABLE: int
TYPEFLAG_FRESTRICTED: int
TYPEFLAG_FREVERSEBIND: int
RecordClasses: dict[str, com_record]
TypeIIDs: dict[_win32typing.PyIID, type]
URL_MK_LEGACY: int
URL_MK_UNIFORM: int
VARFLAG_FREADONLY: int
VAR_CONST: int
VAR_DISPATCH: int
VAR_PERINSTANCE: int
VAR_STATIC: int
VT_ARRAY: int
VT_BLOB: int
VT_BLOB_OBJECT: int
VT_BOOL: int
VT_BSTR: int
VT_BSTR_BLOB: int
VT_BYREF: int
VT_CARRAY: int
VT_CF: int
VT_CLSID: int
VT_CY: int
VT_DATE: int
VT_DECIMAL: int
VT_DISPATCH: int
VT_EMPTY: int
VT_ERROR: int
VT_FILETIME: int
VT_HRESULT: int
VT_I1: int
VT_I2: int
VT_I4: int
VT_I8: int
VT_ILLEGAL: int
VT_ILLEGALMASKED: int
VT_INT: int
VT_LPSTR: int
VT_LPWSTR: int
VT_NULL: int
VT_PTR: int
VT_R4: int
VT_R8: int
VT_RECORD: int
VT_RESERVED: int
VT_SAFEARRAY: int
VT_STORAGE: int
VT_STORED_OBJECT: int
VT_STREAM: int
VT_STREAMED_OBJECT: int
VT_TYPEMASK: int
VT_UI1: int
VT_UI2: int
VT_UI4: int
VT_UI8: int
VT_UINT: int
VT_UNKNOWN: int
VT_USERDEFINED: int
VT_VARIANT: int
VT_VECTOR: int
VT_VOID: int

dcom: int
fdexNameCaseInsensitive: int
fdexNameCaseSensitive: int
fdexNameEnsure: int
fdexNameImplicit: int
fdexPropCanCall: int
fdexPropCanConstruct: int
fdexPropCanGet: int
fdexPropCanPut: int
fdexPropCanPutRef: int
fdexPropCanSourceEvents: int
fdexPropCannotCall: int
fdexPropCannotConstruct: int
fdexPropCannotGet: int
fdexPropCannotPut: int
fdexPropCannotPutRef: int
fdexPropCannotSourceEvents: int
fdexPropDynamicType: int
fdexPropNoSideEffects: int
frozen: int
