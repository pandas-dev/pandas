from _typeshed import Incomplete

def WTSOpenServer(ServerName: str) -> int: ...
def WTSCloseServer(Server: int) -> None: ...
def WTSQueryUserConfig(ServerName: str, UserName: str, WTSConfigClass): ...
def WTSSetUserConfig(ServerName: str, UserName: str, WTSConfigClass, Buffer) -> None: ...
def WTSEnumerateServers(DomainName: str | None = ..., Version: int = ..., Reserved=...) -> tuple[str, ...]: ...
def WTSEnumerateSessions(Server: int = ..., Version: int = ..., Reserved=...) -> tuple[dict[str, str | int], ...]: ...
def WTSLogoffSession(Server: int, SessionId: int, Wait: bool) -> None: ...
def WTSDisconnectSession(Server: int, SessionId: int, Wait: bool) -> None: ...
def WTSQuerySessionInformation(Server: int, SessionId: int, WTSInfoClass: int) -> str: ...
def WTSEnumerateProcesses(Server: int = ..., Version: int = ..., Reserved: int = ...) -> tuple[str, ...]: ...
def WTSQueryUserToken(SessionId) -> int: ...
def WTSShutdownSystem(Server: int, ShutdownFlag) -> None: ...
def WTSTerminateProcess(Server: int, ProcessId, ExitCode) -> None: ...
def ProcessIdToSessionId(ProcessId): ...
def WTSGetActiveConsoleSessionId(): ...
def WTSRegisterSessionNotification(Wnd: int, Flags) -> None: ...
def WTSUnRegisterSessionNotification(Wnd: int) -> None: ...
def WTSWaitSystemEvent(Server: int = ..., EventMask=...): ...
def WTSSendMessage(Server: int, SessionId, Title: str, Message: str, Style, Timeout, Wait): ...

NOTIFY_FOR_ALL_SESSIONS: int
NOTIFY_FOR_THIS_SESSION: int
WTSActive: int
WTSApplicationName: int
WTSClientAddress: int
WTSClientBuildNumber: int
WTSClientDirectory: int
WTSClientDisplay: int
WTSClientHardwareId: int
WTSClientName: int
WTSClientProductId: int
WTSClientProtocolType: int
WTSIsRemoteSession: int
WTSConnectQuery: int
WTSConnectState: int
WTSConnected: int
WTSDisconnected: int
WTSDomainName: int
WTSDown: int
WTSIdle: int
WTSInit: int
WTSInitialProgram: int
WTSListen: int
WTSOEMId: int
WTSReset: int
WTSSessionId: int
WTSShadow: int
WTSUserConfigBrokenTimeoutSettings: int
WTSUserConfigInitialProgram: int
WTSUserConfigModemCallbackPhoneNumber: int
WTSUserConfigModemCallbackSettings: int
WTSUserConfigReconnectSettings: int
WTSUserConfigShadowingSettings: int
WTSUserConfigTerminalServerHomeDir: int
WTSUserConfigTerminalServerHomeDirDrive: int
WTSUserConfigTerminalServerProfilePath: int
WTSUserConfigTimeoutSettingsConnections: int
WTSUserConfigTimeoutSettingsDisconnections: int
WTSUserConfigTimeoutSettingsIdle: int
WTSUserConfigWorkingDirectory: int
WTSUserConfigfAllowLogonTerminalServer: int
WTSUserConfigfDeviceClientDefaultPrinter: int
WTSUserConfigfDeviceClientDrives: int
WTSUserConfigfDeviceClientPrinters: int
WTSUserConfigfInheritInitialProgram: int
WTSUserConfigfTerminalServerRemoteHomeDir: int
WTSUserName: int
WTSVirtualClientData: int
WTSVirtualFileHandle: int
WTSWinStationName: int
WTSWorkingDirectory: int
WTS_CURRENT_SERVER: int
WTS_CURRENT_SERVER_HANDLE: int
WTS_CURRENT_SERVER_NAME: Incomplete
WTS_CURRENT_SESSION: int
WTS_EVENT_ALL: int
WTS_EVENT_CONNECT: int
WTS_EVENT_CREATE: int
WTS_EVENT_DELETE: int
WTS_EVENT_DISCONNECT: int
WTS_EVENT_FLUSH: int
WTS_EVENT_LICENSE: int
WTS_EVENT_LOGOFF: int
WTS_EVENT_LOGON: int
WTS_EVENT_NONE: int
WTS_EVENT_RENAME: int
WTS_EVENT_STATECHANGE: int
WTS_PROTOCOL_TYPE_CONSOLE: int
WTS_PROTOCOL_TYPE_ICA: int
WTS_PROTOCOL_TYPE_RDP: int
WTS_WSD_FASTREBOOT: int
WTS_WSD_LOGOFF: int
WTS_WSD_POWEROFF: int
WTS_WSD_REBOOT: int
WTS_WSD_SHUTDOWN: int
