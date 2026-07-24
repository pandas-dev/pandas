from _typeshed import Incomplete
from typing_extensions import TypeAlias

import _win32typing
from win32.lib.pywintypes import com_error

error: TypeAlias = com_error  # noqa: Y042

def ADsOpenObject(path, username, password, iid: _win32typing.PyIID, reserved: int = ..., /): ...
def ADsGetObject(path, iid: _win32typing.PyIID, /): ...
def ADsBuildEnumerator(container: _win32typing.PyIADsContainer, /): ...
def ADsEnumerateNext(enum, num: int = ..., /): ...
def ADsGetLastError() -> tuple[Incomplete, Incomplete, Incomplete]: ...
def StringAsDS_SELECTION_LIST(buf, /): ...

DSOP_SCOPE_INIT_INFOs = _win32typing.PyDSOP_SCOPE_INIT_INFOs
CLSID_ADsDSOObject: _win32typing.PyIID
CLSID_AccessControlEntry: _win32typing.PyIID
CLSID_AccessControlList: _win32typing.PyIID
CLSID_DsObjectPicker: _win32typing.PyIID
CLSID_SecurityDescriptor: _win32typing.PyIID
DBGUID_LDAPDialect: _win32typing.PyIID
DBPROPSET_ADSISEARCH: _win32typing.PyIID
IID_IADs: _win32typing.PyIID
IID_IADsClass: _win32typing.PyIID
IID_IADsCollection: _win32typing.PyIID
IID_IADsComputer: _win32typing.PyIID
IID_IADsComputerOperations: _win32typing.PyIID
IID_IADsContainer: _win32typing.PyIID
IID_IADsDeleteOps: _win32typing.PyIID
IID_IADsDomain: _win32typing.PyIID
IID_IADsFileService: _win32typing.PyIID
IID_IADsFileServiceOperations: _win32typing.PyIID
IID_IADsFileShare: _win32typing.PyIID
IID_IADsGroup: _win32typing.PyIID
IID_IADsLocality: _win32typing.PyIID
IID_IADsMembers: _win32typing.PyIID
IID_IADsNamespaces: _win32typing.PyIID
IID_IADsO: _win32typing.PyIID
IID_IADsOU: _win32typing.PyIID
IID_IADsOpenDSObject: _win32typing.PyIID
IID_IADsPrintJob: _win32typing.PyIID
IID_IADsPrintJobOperations: _win32typing.PyIID
IID_IADsPrintQueue: _win32typing.PyIID
IID_IADsPrintQueueOperations: _win32typing.PyIID
IID_IADsProperty: _win32typing.PyIID
IID_IADsPropertyList: _win32typing.PyIID
IID_IADsResource: _win32typing.PyIID
IID_IADsSearch: _win32typing.PyIID
IID_IADsService: _win32typing.PyIID
IID_IADsServiceOperations: _win32typing.PyIID
IID_IADsSession: _win32typing.PyIID
IID_IADsSyntax: _win32typing.PyIID
IID_IADsUser: _win32typing.PyIID
IID_IDirectoryObject: _win32typing.PyIID
IID_IDirectorySearch: _win32typing.PyIID
IID_IDsObjectPicker: _win32typing.PyIID
LIBID_ADs: _win32typing.PyIID
