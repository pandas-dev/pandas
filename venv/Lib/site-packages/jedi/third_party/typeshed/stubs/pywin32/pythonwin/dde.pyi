# Can't generate with stubgen because:
# "ImportError: This must be an MFC application - try 'import win32ui' first"
APPCLASS_MONITOR: int
APPCLASS_STANDARD: int
APPCMD_CLIENTONLY: int
APPCMD_FILTERINITS: int
CBF_FAIL_ADVISES: int
CBF_FAIL_ALLSVRXACTIONS: int
CBF_FAIL_CONNECTIONS: int
CBF_FAIL_EXECUTES: int
CBF_FAIL_POKES: int
CBF_FAIL_REQUESTS: int
CBF_FAIL_SELFCONNECTIONS: int
CBF_SKIP_ALLNOTIFICATIONS: int
CBF_SKIP_CONNECT_CONFIRMS: int
CBF_SKIP_DISCONNECTS: int
CBF_SKIP_REGISTRATIONS: int

def CreateConversation(Server, /): ...
def CreateServer(): ...
def CreateServerSystemTopic(): ...
def CreateStringItem(name, /): ...
def CreateTopic(name, /): ...

MF_CALLBACKS: int
MF_CONV: int
MF_ERRORS: int
MF_HSZ_INFO: int
MF_LINKS: int
MF_POSTMSGS: int
MF_SENDMSGS: int

class error(Exception): ...
