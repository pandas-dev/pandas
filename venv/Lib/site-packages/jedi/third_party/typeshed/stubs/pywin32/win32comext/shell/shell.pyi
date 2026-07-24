from _typeshed import Incomplete
from typing_extensions import TypeAlias

import _win32typing
from win32.lib.pywintypes import com_error

error: TypeAlias = com_error  # noqa: Y042

def AssocCreate() -> _win32typing.PyIQueryAssociations: ...
def AssocCreateForClasses() -> _win32typing.PyIUnknown: ...
def DragQueryFile(hglobal: int, index, /) -> str: ...
def DragQueryFileW(hglobal: int, index, /) -> str: ...
def DragQueryPoint(hglobal: int, /) -> tuple[Incomplete, Incomplete, Incomplete]: ...
def IsUserAnAdmin() -> bool: ...
def SHCreateDataObject(
    parent, children: list[Incomplete], do_inner: _win32typing.PyIDataObject, iid: _win32typing.PyIID, /
) -> _win32typing.PyIUnknown: ...
def SHCreateDefaultContextMenu(dcm, iid: _win32typing.PyIID, /) -> _win32typing.PyIUnknown: ...
def SHCreateDefaultExtractIcon() -> _win32typing.PyIDefaultExtractIconInit: ...
def SHCreateShellFolderView(
    sf: _win32typing.PyIShellFolder, viewOuter: _win32typing.PyIShellView | None = ..., callbacks: Incomplete | None = ..., /
) -> _win32typing.PyIShellView: ...
def SHCreateShellItemArray(
    parent: _win32typing.PyIDL, sf: _win32typing.PyIShellFolder, children: list[_win32typing.PyIDL], /
) -> _win32typing.PyIShellItemArray: ...
def SHCreateShellItemArrayFromDataObject(
    do: _win32typing.PyIDataObject, iid: _win32typing.PyIID, /
) -> _win32typing.PyIShellItemArray: ...
def SHCreateShellItemArrayFromShellItem(
    si: _win32typing.PyIShellItem, riid: _win32typing.PyIID, /
) -> _win32typing.PyIShellItemArray: ...
def SHBrowseForFolder(
    hwndOwner: int | None = ...,
    pidlRoot: _win32typing.PyIDL | None = ...,
    title: str | None = ...,
    flags: int = ...,
    callback: Incomplete | None = ...,
    callback_data: Incomplete | None = ...,
    /,
) -> tuple[_win32typing.PyIDL, Incomplete, Incomplete]: ...
def SHGetFileInfo(
    name: _win32typing.PyIDL | str, dwFileAttributes, uFlags, infoAttrs: int = ..., /
) -> tuple[Incomplete, _win32typing.SHFILEINFO]: ...
def SHGetFolderPath(hwndOwner: int, nFolder, handle: int, flags, /) -> str: ...
def SHSetFolderPath(csidl, Path, hToken: int | None = ..., /) -> None: ...
def SHGetFolderLocation(hwndOwner: int, nFolder, hToken: int | None = ..., reserved=..., /) -> _win32typing.PyIDL: ...
def SHGetSpecialFolderPath(hwndOwner: int, nFolder, bCreate: int = ..., /) -> str: ...
def SHGetSpecialFolderLocation(hwndOwner: int, nFolder, /) -> _win32typing.PyIDL: ...
def SHAddToRecentDocs(Flags, data, /) -> None: ...
def SHEmptyRecycleBin(hwnd: int, path: str, flags, /) -> None: ...
def SHQueryRecycleBin(RootPath: str | None = ..., /) -> tuple[Incomplete, Incomplete]: ...
def SHGetDesktopFolder() -> _win32typing.PyIShellFolder: ...
def SHUpdateImage(HashItem: str, Index, Flags, ImageIndex, /) -> None: ...
def SHChangeNotify(EventId, Flags, Item1, Item2, /) -> None: ...
def SHChangeNotifyRegister(hwnd: int, sources, events, msg, /): ...
def SHChangeNotifyDeregister(_id, /) -> None: ...
def SHCreateItemFromParsingName(name, ctx: _win32typing.PyIBindCtx, riid: _win32typing.PyIID, /) -> _win32typing.PyIShellItem: ...
def SHCreateItemFromRelativeName(
    Parent: _win32typing.PyIShellItem, Name, ctx: _win32typing.PyIBindCtx, riid: _win32typing.PyIID, /
) -> _win32typing.PyIShellItem: ...
def SHCreateItemInKnownFolder(
    FolderId: _win32typing.PyIID, Flags, Name, riid: _win32typing.PyIID, /
) -> _win32typing.PyIShellItem: ...
def SHCreateItemWithParent(
    Parent: _win32typing.PyIDL, sfParent: _win32typing.PyIShellFolder, child: _win32typing.PyIDL, riid: _win32typing.PyIID, /
) -> _win32typing.PyIShellItem: ...
def SHGetInstanceExplorer() -> _win32typing.PyIUnknown: ...
def SHFileOperation(operation: _win32typing.SHFILEOPSTRUCT, /) -> tuple[Incomplete, Incomplete]: ...
def StringAsCIDA(pidl: str, /) -> tuple[_win32typing.PyIDL, Incomplete]: ...
def CIDAAsString(pidl: str, /) -> str: ...
def StringAsPIDL(pidl: str, /) -> _win32typing.PyIDL: ...
def AddressAsPIDL(address, /) -> _win32typing.PyIDL: ...
def PIDLAsString(pidl: _win32typing.PyIDL, /) -> str: ...
def SHGetSettings(mask: int = ..., /): ...
def FILEGROUPDESCRIPTORAsString(descriptors: list[Incomplete], arg, /) -> str: ...
def StringAsFILEGROUPDESCRIPTOR(buf, make_unicode: int = ..., /) -> list[Incomplete]: ...
def ShellExecuteEx(
    fMask: int = ...,
    hwnd: int = ...,
    lpVerb: str = ...,
    lpFile: str = ...,
    lpParameters: str = ...,
    lpDirectory: str = ...,
    nShow: int = ...,
    lpIDlist: _win32typing.PyIDL = ...,
    lpClass: str = ...,
    hkeyClass=...,
    dwHotKey=...,
    hIcon: int = ...,
    hMonitor: int = ...,
): ...
def SHGetViewStatePropertyBag(
    pidl: _win32typing.PyIDL, BagName: str, Flags, riid: _win32typing.PyIID, /
) -> _win32typing.PyIPropertyBag: ...
def SHILCreateFromPath(Path: str, Flags, /) -> tuple[_win32typing.PyIDL, Incomplete]: ...
def SHCreateShellItem(
    pidlParent: _win32typing.PyIDL, sfParent: _win32typing.PyIShellFolder, Child: _win32typing.PyIDL, /
) -> _win32typing.PyIShellItem: ...
def SHOpenFolderAndSelectItems(Folder: _win32typing.PyIDL, Items: tuple[_win32typing.PyIDL, ...], Flags=...) -> None: ...
def SHCreateStreamOnFileEx(
    File: str, Mode: int, Attributes: int, Create: bool, Template: None = None
) -> _win32typing.PyIStream: ...
def SetCurrentProcessExplicitAppUserModelID(AppID: str, /) -> None: ...
def GetCurrentProcessExplicitAppUserModelID() -> str: ...
def SHParseDisplayName(Name, Attributes, BindCtx: _win32typing.PyIBindCtx | None = ...) -> tuple[list[bytes], int]: ...
def SHCreateItemFromIDList(pidl, riid=..., /): ...
def SHCreateShellItemArrayFromIDLists(pidls, /): ...
def SHGetIDListFromObject(unk, /): ...
def SHGetNameFromIDList(pidl, flags: int, /): ...
def SHGetPathFromIDList(pidl, /): ...
def SHGetPathFromIDListW(Pidl, /): ...

BHID_AssociationArray: _win32typing.PyIID
BHID_DataObject: _win32typing.PyIID
BHID_EnumItems: _win32typing.PyIID
BHID_Filter: _win32typing.PyIID
BHID_LinkTargetItem: _win32typing.PyIID
BHID_PropertyStore: _win32typing.PyIID
BHID_SFObject: _win32typing.PyIID
BHID_SFUIObject: _win32typing.PyIID
BHID_SFViewObject: _win32typing.PyIID
BHID_Storage: _win32typing.PyIID
BHID_StorageEnum: _win32typing.PyIID
BHID_Stream: _win32typing.PyIID
BHID_ThumbnailHandler: _win32typing.PyIID
BHID_Transfer: _win32typing.PyIID
CGID_DefView: _win32typing.PyIID
CGID_Explorer: _win32typing.PyIID
CGID_ExplorerBarDoc: _win32typing.PyIID
CGID_ShellDocView: _win32typing.PyIID
CGID_ShellServiceObject: _win32typing.PyIID
CLSID_ActiveDesktop: _win32typing.PyIID
CLSID_ApplicationDestinations: _win32typing.PyIID
CLSID_ApplicationDocumentLists: _win32typing.PyIID
CLSID_ControlPanel: _win32typing.PyIID
CLSID_DestinationList: _win32typing.PyIID
CLSID_DragDropHelper: _win32typing.PyIID
CLSID_EnumerableObjectCollection: _win32typing.PyIID
CLSID_FileOperation: _win32typing.PyIID
CLSID_Internet: _win32typing.PyIID
CLSID_InternetShortcut: _win32typing.PyIID
CLSID_KnownFolderManager: _win32typing.PyIID
CLSID_MyComputer: _win32typing.PyIID
CLSID_MyDocuments: _win32typing.PyIID
CLSID_NetworkDomain: _win32typing.PyIID
CLSID_NetworkPlaces: _win32typing.PyIID
CLSID_NetworkServer: _win32typing.PyIID
CLSID_NetworkShare: _win32typing.PyIID
CLSID_Printers: _win32typing.PyIID
CLSID_RecycleBin: _win32typing.PyIID
CLSID_ShellDesktop: _win32typing.PyIID
CLSID_ShellFSFolder: _win32typing.PyIID
CLSID_ShellItem: _win32typing.PyIID
CLSID_ShellLibrary: _win32typing.PyIID
CLSID_ShellLink: _win32typing.PyIID
CLSID_TaskbarList: _win32typing.PyIID
EP_AdvQueryPane: _win32typing.PyIID
EP_Commands: _win32typing.PyIID
EP_Commands_Organize: _win32typing.PyIID
EP_Commands_View: _win32typing.PyIID
EP_DetailsPane: _win32typing.PyIID
EP_NavPane: _win32typing.PyIID
EP_PreviewPane: _win32typing.PyIID
EP_QueryPane: _win32typing.PyIID
FMTID_AudioSummaryInformation: _win32typing.PyIID
FMTID_Briefcase: _win32typing.PyIID
FMTID_Displaced: _win32typing.PyIID
FMTID_ImageProperties: _win32typing.PyIID
FMTID_ImageSummaryInformation: _win32typing.PyIID
FMTID_InternetSite: _win32typing.PyIID
FMTID_Intshcut: _win32typing.PyIID
FMTID_MediaFileSummaryInformation: _win32typing.PyIID
FMTID_Misc: _win32typing.PyIID
FMTID_Query: _win32typing.PyIID
FMTID_ShellDetails: _win32typing.PyIID
FMTID_Storage: _win32typing.PyIID
FMTID_SummaryInformation: _win32typing.PyIID
FMTID_Volume: _win32typing.PyIID
FMTID_WebView: _win32typing.PyIID
FOLDERID_AddNewPrograms: _win32typing.PyIID
FOLDERID_AdminTools: _win32typing.PyIID
FOLDERID_AppUpdates: _win32typing.PyIID
FOLDERID_CDBurning: _win32typing.PyIID
FOLDERID_ChangeRemovePrograms: _win32typing.PyIID
FOLDERID_CommonAdminTools: _win32typing.PyIID
FOLDERID_CommonOEMLinks: _win32typing.PyIID
FOLDERID_CommonPrograms: _win32typing.PyIID
FOLDERID_CommonStartMenu: _win32typing.PyIID
FOLDERID_CommonStartup: _win32typing.PyIID
FOLDERID_CommonTemplates: _win32typing.PyIID
FOLDERID_ComputerFolder: _win32typing.PyIID
FOLDERID_ConflictFolder: _win32typing.PyIID
FOLDERID_ConnectionsFolder: _win32typing.PyIID
FOLDERID_Contacts: _win32typing.PyIID
FOLDERID_ControlPanelFolder: _win32typing.PyIID
FOLDERID_Cookies: _win32typing.PyIID
FOLDERID_Desktop: _win32typing.PyIID
FOLDERID_DeviceMetadataStore: _win32typing.PyIID
FOLDERID_Documents: _win32typing.PyIID
FOLDERID_DocumentsLibrary: _win32typing.PyIID
FOLDERID_Downloads: _win32typing.PyIID
FOLDERID_Favorites: _win32typing.PyIID
FOLDERID_Fonts: _win32typing.PyIID
FOLDERID_GameTasks: _win32typing.PyIID
FOLDERID_Games: _win32typing.PyIID
FOLDERID_History: _win32typing.PyIID
FOLDERID_HomeGroup: _win32typing.PyIID
FOLDERID_ImplicitAppShortcuts: _win32typing.PyIID
FOLDERID_InternetCache: _win32typing.PyIID
FOLDERID_InternetFolder: _win32typing.PyIID
FOLDERID_Libraries: _win32typing.PyIID
FOLDERID_Links: _win32typing.PyIID
FOLDERID_LocalAppData: _win32typing.PyIID
FOLDERID_LocalAppDataLow: _win32typing.PyIID
FOLDERID_LocalizedResourcesDir: _win32typing.PyIID
FOLDERID_Music: _win32typing.PyIID
FOLDERID_MusicLibrary: _win32typing.PyIID
FOLDERID_NetHood: _win32typing.PyIID
FOLDERID_NetworkFolder: _win32typing.PyIID
FOLDERID_OriginalImages: _win32typing.PyIID
FOLDERID_PhotoAlbums: _win32typing.PyIID
FOLDERID_Pictures: _win32typing.PyIID
FOLDERID_PicturesLibrary: _win32typing.PyIID
FOLDERID_Playlists: _win32typing.PyIID
FOLDERID_PrintHood: _win32typing.PyIID
FOLDERID_PrintersFolder: _win32typing.PyIID
FOLDERID_Profile: _win32typing.PyIID
FOLDERID_ProgramData: _win32typing.PyIID
FOLDERID_ProgramFiles: _win32typing.PyIID
FOLDERID_ProgramFilesCommon: _win32typing.PyIID
FOLDERID_ProgramFilesCommonX64: _win32typing.PyIID
FOLDERID_ProgramFilesCommonX86: _win32typing.PyIID
FOLDERID_ProgramFilesX64: _win32typing.PyIID
FOLDERID_ProgramFilesX86: _win32typing.PyIID
FOLDERID_Programs: _win32typing.PyIID
FOLDERID_Public: _win32typing.PyIID
FOLDERID_PublicDesktop: _win32typing.PyIID
FOLDERID_PublicDocuments: _win32typing.PyIID
FOLDERID_PublicDownloads: _win32typing.PyIID
FOLDERID_PublicGameTasks: _win32typing.PyIID
FOLDERID_PublicLibraries: _win32typing.PyIID
FOLDERID_PublicMusic: _win32typing.PyIID
FOLDERID_PublicPictures: _win32typing.PyIID
FOLDERID_PublicRingtones: _win32typing.PyIID
FOLDERID_PublicVideos: _win32typing.PyIID
FOLDERID_QuickLaunch: _win32typing.PyIID
FOLDERID_Recent: _win32typing.PyIID
FOLDERID_RecordedTVLibrary: _win32typing.PyIID
FOLDERID_RecycleBinFolder: _win32typing.PyIID
FOLDERID_ResourceDir: _win32typing.PyIID
FOLDERID_Ringtones: _win32typing.PyIID
FOLDERID_RoamingAppData: _win32typing.PyIID
FOLDERID_SEARCH_CSC: _win32typing.PyIID
FOLDERID_SEARCH_MAPI: _win32typing.PyIID
FOLDERID_SampleMusic: _win32typing.PyIID
FOLDERID_SamplePictures: _win32typing.PyIID
FOLDERID_SamplePlaylists: _win32typing.PyIID
FOLDERID_SampleVideos: _win32typing.PyIID
FOLDERID_SavedGames: _win32typing.PyIID
FOLDERID_SavedSearches: _win32typing.PyIID
FOLDERID_SearchHome: _win32typing.PyIID
FOLDERID_SendTo: _win32typing.PyIID
FOLDERID_SidebarDefaultParts: _win32typing.PyIID
FOLDERID_SidebarParts: _win32typing.PyIID
FOLDERID_StartMenu: _win32typing.PyIID
FOLDERID_Startup: _win32typing.PyIID
FOLDERID_SyncManagerFolder: _win32typing.PyIID
FOLDERID_SyncResultsFolder: _win32typing.PyIID
FOLDERID_SyncSetupFolder: _win32typing.PyIID
FOLDERID_System: _win32typing.PyIID
FOLDERID_SystemX86: _win32typing.PyIID
FOLDERID_Templates: _win32typing.PyIID
FOLDERID_UserPinned: _win32typing.PyIID
FOLDERID_UserProfiles: _win32typing.PyIID
FOLDERID_UserProgramFiles: _win32typing.PyIID
FOLDERID_UserProgramFilesCommon: _win32typing.PyIID
FOLDERID_UsersFiles: _win32typing.PyIID
FOLDERID_UsersLibraries: _win32typing.PyIID
FOLDERID_Videos: _win32typing.PyIID
FOLDERID_VideosLibrary: _win32typing.PyIID
FOLDERID_Windows: _win32typing.PyIID
FOLDERTYPEID_Communications: _win32typing.PyIID
FOLDERTYPEID_CompressedFolder: _win32typing.PyIID
FOLDERTYPEID_Contacts: _win32typing.PyIID
FOLDERTYPEID_ControlPanelCategory: _win32typing.PyIID
FOLDERTYPEID_ControlPanelClassic: _win32typing.PyIID
FOLDERTYPEID_Documents: _win32typing.PyIID
FOLDERTYPEID_Games: _win32typing.PyIID
FOLDERTYPEID_Generic: _win32typing.PyIID
FOLDERTYPEID_GenericLibrary: _win32typing.PyIID
FOLDERTYPEID_GenericSearchResults: _win32typing.PyIID
FOLDERTYPEID_Invalid: _win32typing.PyIID
FOLDERTYPEID_Music: _win32typing.PyIID
FOLDERTYPEID_NetworkExplorer: _win32typing.PyIID
FOLDERTYPEID_OpenSearch: _win32typing.PyIID
FOLDERTYPEID_OtherUsers: _win32typing.PyIID
FOLDERTYPEID_Pictures: _win32typing.PyIID
FOLDERTYPEID_Printers: _win32typing.PyIID
FOLDERTYPEID_PublishedItems: _win32typing.PyIID
FOLDERTYPEID_RecordedTV: _win32typing.PyIID
FOLDERTYPEID_RecycleBin: _win32typing.PyIID
FOLDERTYPEID_SavedGames: _win32typing.PyIID
FOLDERTYPEID_SearchConnector: _win32typing.PyIID
FOLDERTYPEID_SearchHome: _win32typing.PyIID
FOLDERTYPEID_Searches: _win32typing.PyIID
FOLDERTYPEID_SoftwareExplorer: _win32typing.PyIID
FOLDERTYPEID_StartMenu: _win32typing.PyIID
FOLDERTYPEID_UserFiles: _win32typing.PyIID
FOLDERTYPEID_UsersLibraries: _win32typing.PyIID
FOLDERTYPEID_Videos: _win32typing.PyIID
HOTKEYF_ALT: int
HOTKEYF_CONTROL: int
HOTKEYF_EXT: int
HOTKEYF_SHIFT: int
IID_CDefView: _win32typing.PyIID
IID_IADesktopP2: _win32typing.PyIID
IID_IActiveDesktop: _win32typing.PyIID
IID_IActiveDesktopP: _win32typing.PyIID
IID_IApplicationDestinations: _win32typing.PyIID
IID_IApplicationDocumentLists: _win32typing.PyIID
IID_IAsyncOperation: _win32typing.PyIID
IID_IBrowserFrameOptions: _win32typing.PyIID
IID_ICategorizer: _win32typing.PyIID
IID_ICategoryProvider: _win32typing.PyIID
IID_IColumnProvider: _win32typing.PyIID
IID_IContextMenu: _win32typing.PyIID
IID_IContextMenu2: _win32typing.PyIID
IID_IContextMenu3: _win32typing.PyIID
IID_ICopyHook: _win32typing.PyIID
IID_ICopyHookA: _win32typing.PyIID
IID_ICopyHookW: _win32typing.PyIID
IID_ICurrentItem: _win32typing.PyIID
IID_ICustomDestinationList: _win32typing.PyIID
IID_IDefaultExtractIconInit: _win32typing.PyIID
IID_IDeskBand: _win32typing.PyIID
IID_IDisplayItem: _win32typing.PyIID
IID_IDockingWindow: _win32typing.PyIID
IID_IDropTargetHelper: _win32typing.PyIID
IID_IEmptyVolumeCache: _win32typing.PyIID
IID_IEmptyVolumeCache2: _win32typing.PyIID
IID_IEmptyVolumeCacheCallBack: _win32typing.PyIID
IID_IEnumExplorerCommand: _win32typing.PyIID
IID_IEnumIDList: _win32typing.PyIID
IID_IEnumObjects: _win32typing.PyIID
IID_IEnumResources: _win32typing.PyIID
IID_IEnumShellItems: _win32typing.PyIID
IID_IExplorerBrowser: _win32typing.PyIID
IID_IExplorerBrowserEvents: _win32typing.PyIID
IID_IExplorerCommand: _win32typing.PyIID
IID_IExplorerCommandProvider: _win32typing.PyIID
IID_IExplorerPaneVisibility: _win32typing.PyIID
IID_IExtractIcon: _win32typing.PyIID
IID_IExtractIconW: _win32typing.PyIID
IID_IExtractImage: _win32typing.PyIID
IID_IFileOperation: _win32typing.PyIID
IID_IFileOperationProgressSink: _win32typing.PyIID
IID_IFolderView: _win32typing.PyIID
IID_IIdentityName: _win32typing.PyIID
IID_IKnownFolder: _win32typing.PyIID
IID_IKnownFolderManager: _win32typing.PyIID
IID_INameSpaceTreeControl: _win32typing.PyIID
IID_IObjectArray: _win32typing.PyIID
IID_IObjectCollection: _win32typing.PyIID
IID_IPersistFolder: _win32typing.PyIID
IID_IPersistFolder2: _win32typing.PyIID
IID_IQueryAssociations: _win32typing.PyIID
IID_IRelatedItem: _win32typing.PyIID
IID_IShellBrowser: _win32typing.PyIID
IID_IShellCopyHook: _win32typing.PyIID
IID_IShellCopyHookA: _win32typing.PyIID
IID_IShellCopyHookW: _win32typing.PyIID
IID_IShellExtInit: _win32typing.PyIID
IID_IShellFolder: _win32typing.PyIID
IID_IShellFolder2: _win32typing.PyIID
IID_IShellIcon: _win32typing.PyIID
IID_IShellIconOverlay: _win32typing.PyIID
IID_IShellIconOverlayIdentifier: _win32typing.PyIID
IID_IShellIconOverlayManager: _win32typing.PyIID
IID_IShellItem: _win32typing.PyIID
IID_IShellItem2: _win32typing.PyIID
IID_IShellItemArray: _win32typing.PyIID
IID_IShellItemResources: _win32typing.PyIID
IID_IShellLibrary: _win32typing.PyIID
IID_IShellLink: _win32typing.PyIID
IID_IShellLinkA: _win32typing.PyIID
IID_IShellLinkDataList: _win32typing.PyIID
IID_IShellLinkW: _win32typing.PyIID
IID_IShellView: _win32typing.PyIID
IID_ITaskbarList: _win32typing.PyIID
IID_ITransferAdviseSink: _win32typing.PyIID
IID_ITransferDestination: _win32typing.PyIID
IID_ITransferMediumItem: _win32typing.PyIID
IID_ITransferSource: _win32typing.PyIID
IID_IUniformResourceLocator: _win32typing.PyIID
ResourceTypeStream: _win32typing.PyIID
SID_CtxQueryAssociations: _win32typing.PyIID
SID_DefView: _win32typing.PyIID
SID_LinkSite: _win32typing.PyIID
SID_MenuShellFolder: _win32typing.PyIID
SID_SCommDlgBrowser: _win32typing.PyIID
SID_SGetViewFromViewDual: _win32typing.PyIID
SID_SInternetExplorer: _win32typing.PyIID
SID_SMenuBandBKContextMenu: _win32typing.PyIID
SID_SMenuBandBottom: _win32typing.PyIID
SID_SMenuBandBottomSelected: _win32typing.PyIID
SID_SMenuBandChild: _win32typing.PyIID
SID_SMenuBandContextMenuModifier: _win32typing.PyIID
SID_SMenuBandParent: _win32typing.PyIID
SID_SMenuBandTop: _win32typing.PyIID
SID_SMenuPopup: _win32typing.PyIID
SID_SProgressUI: _win32typing.PyIID
SID_SShellBrowser: _win32typing.PyIID
SID_SShellDesktop: _win32typing.PyIID
SID_STopLevelBrowser: _win32typing.PyIID
SID_STopWindow: _win32typing.PyIID
SID_SUrlHistory: _win32typing.PyIID
SID_SWebBrowserApp: _win32typing.PyIID
SID_ShellFolderViewCB: _win32typing.PyIID
SLGP_RAWPATH: int
SLGP_SHORTPATH: int
SLGP_UNCPRIORITY: int
SLR_ANY_MATCH: int
SLR_INVOKE_MSI: int
SLR_NOLINKINFO: int
SLR_NOSEARCH: int
SLR_NOTRACK: int
SLR_NOUPDATE: int
SLR_NO_UI: int
SLR_UPDATE: int
VID_Details: _win32typing.PyIID
VID_LargeIcons: _win32typing.PyIID
VID_List: _win32typing.PyIID
VID_SmallIcons: _win32typing.PyIID
VID_ThumbStrip: _win32typing.PyIID
VID_Thumbnails: _win32typing.PyIID
VID_Tile: _win32typing.PyIID

def SHGetKnownFolderPath(fid, flags: int = 0, token=None, /): ...
