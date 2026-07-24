from _typeshed import Incomplete

def CredMarshalCredential(CredType, Credential: str) -> str: ...
def CredUnmarshalCredential(MarshaledCredential: str) -> tuple[Incomplete, str]: ...
def CredIsMarshaledCredential(MarshaledCredential: str) -> bool: ...
def CredEnumerate(Filter: str | None = ..., Flags: int = ...) -> tuple[dict[str, Incomplete], ...]: ...
def CredGetTargetInfo(TargetName: str, Flags: int = ...): ...
def CredGetSessionTypes(MaximumPersistCount: int = 7) -> tuple[int, ...]: ...
def CredWriteDomainCredentials(TargetInfo, Credential, Flags: int = ...) -> None: ...
def CredReadDomainCredentials(TargetInfo, Flags: int = ...) -> tuple[Incomplete, ...]: ...
def CredDelete(TargetName: str, Type, Flags: int = ...) -> None: ...
def CredWrite(Credential, Flags: int = ...) -> None: ...
def CredRead(TargetName: str, Type, Flags: int = ...): ...
def CredRename(OldTargetName: str, NewTargetName: str, Type, Flags: int = ...): ...
def CredUICmdLinePromptForCredentials(
    TargetName: str, AuthError: int = ..., UserName: str | None = ..., Password: str | None = ..., Save: int = ..., Flags=...
) -> tuple[str, str, Incomplete]: ...
def CredUIPromptForCredentials(
    TargetName: str,
    AuthError: int = ...,
    UserName: str | None = ...,
    Password: str | None = ...,
    Save: bool = ...,
    Flags: int = ...,
    UiInfo: Incomplete | None = ...,
) -> tuple[str, str, Incomplete]: ...
def CredUIConfirmCredentials(TargetName: str, Confirm) -> None: ...
def CredUIReadSSOCredW(Realm: str | None = ...) -> str: ...
def CredUIStoreSSOCredW(Realm: str, Username: str, Password: str, Persist) -> None: ...
def CredUIParseUserName(UserName: str) -> tuple[str, str]: ...

CREDUI_FLAGS_ALWAYS_SHOW_UI: int
CREDUI_FLAGS_COMPLETE_USERNAME: int
CREDUI_FLAGS_DO_NOT_PERSIST: int
CREDUI_FLAGS_EXCLUDE_CERTIFICATES: int
CREDUI_FLAGS_EXPECT_CONFIRMATION: int
CREDUI_FLAGS_GENERIC_CREDENTIALS: int
CREDUI_FLAGS_INCORRECT_PASSWORD: int
CREDUI_FLAGS_KEEP_USERNAME: int
CREDUI_FLAGS_PASSWORD_ONLY_OK: int
CREDUI_FLAGS_PERSIST: int
CREDUI_FLAGS_PROMPT_VALID: int
CREDUI_FLAGS_REQUEST_ADMINISTRATOR: int
CREDUI_FLAGS_REQUIRE_CERTIFICATE: int
CREDUI_FLAGS_REQUIRE_SMARTCARD: int
CREDUI_FLAGS_SERVER_CREDENTIAL: int
CREDUI_FLAGS_SHOW_SAVE_CHECK_BOX: int
CREDUI_FLAGS_USERNAME_TARGET_CREDENTIALS: int
CREDUI_FLAGS_VALIDATE_USERNAME: int
CREDUI_MAX_CAPTION_LENGTH: int
CREDUI_MAX_DOMAIN_TARGET_LENGTH: int
CREDUI_MAX_GENERIC_TARGET_LENGTH: int
CREDUI_MAX_MESSAGE_LENGTH: int
CREDUI_MAX_PASSWORD_LENGTH: int
CREDUI_MAX_USERNAME_LENGTH: int
CRED_ALLOW_NAME_RESOLUTION: int
CRED_CACHE_TARGET_INFORMATION: int
CRED_ENUMERATE_ALL_CREDENTIALS: int
CRED_FLAGS_OWF_CRED_BLOB: int
CRED_FLAGS_PASSWORD_FOR_CERT: int
CRED_FLAGS_PROMPT_NOW: int
CRED_FLAGS_USERNAME_TARGET: int
CRED_FLAGS_VALID_FLAGS: int
CRED_MAX_ATTRIBUTES: int
CRED_MAX_DOMAIN_TARGET_NAME_LENGTH: int
CRED_MAX_GENERIC_TARGET_NAME_LENGTH: int
CRED_MAX_STRING_LENGTH: int
CRED_MAX_USERNAME_LENGTH: int
CRED_MAX_VALUE_SIZE: int
CRED_PERSIST_ENTERPRISE: int
CRED_PERSIST_LOCAL_MACHINE: int
CRED_PERSIST_NONE: int
CRED_PERSIST_SESSION: int
CRED_PRESERVE_CREDENTIAL_BLOB: int
CRED_TI_CREATE_EXPLICIT_CRED: int
CRED_TI_DOMAIN_FORMAT_UNKNOWN: int
CRED_TI_ONLY_PASSWORD_REQUIRED: int
CRED_TI_SERVER_FORMAT_UNKNOWN: int
CRED_TI_USERNAME_TARGET: int
CRED_TI_VALID_FLAGS: int
CRED_TI_WORKGROUP_MEMBER: int
CRED_TYPE_DOMAIN_CERTIFICATE: int
CRED_TYPE_DOMAIN_EXTENDED: int
CRED_TYPE_DOMAIN_PASSWORD: int
CRED_TYPE_DOMAIN_VISIBLE_PASSWORD: int
CRED_TYPE_GENERIC: int
CRED_TYPE_GENERIC_CERTIFICATE: int
CRED_TYPE_MAXIMUM: int
CRED_TYPE_MAXIMUM_EX: int
CertCredential: int
UsernameTargetCredential: int
