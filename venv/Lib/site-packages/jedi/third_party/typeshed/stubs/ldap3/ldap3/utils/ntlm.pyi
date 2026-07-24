from _typeshed import Incomplete

oem_encoding: Incomplete
NTLM_SIGNATURE: bytes
NTLM_MESSAGE_TYPE_NTLM_NEGOTIATE: int
NTLM_MESSAGE_TYPE_NTLM_CHALLENGE: int
NTLM_MESSAGE_TYPE_NTLM_AUTHENTICATE: int
FLAG_NEGOTIATE_56: int
FLAG_NEGOTIATE_KEY_EXCH: int
FLAG_NEGOTIATE_128: int
FLAG_NEGOTIATE_VERSION: int
FLAG_NEGOTIATE_TARGET_INFO: int
FLAG_REQUEST_NOT_NT_SESSION_KEY: int
FLAG_NEGOTIATE_IDENTIFY: int
FLAG_NEGOTIATE_EXTENDED_SESSIONSECURITY: int
FLAG_TARGET_TYPE_SERVER: int
FLAG_TARGET_TYPE_DOMAIN: int
FLAG_NEGOTIATE_ALWAYS_SIGN: int
FLAG_NEGOTIATE_OEM_WORKSTATION_SUPPLIED: int
FLAG_NEGOTIATE_OEM_DOMAIN_SUPPLIED: int
FLAG_NEGOTIATE_ANONYMOUS: int
FLAG_NEGOTIATE_NTLM: int
FLAG_NEGOTIATE_LM_KEY: int
FLAG_NEGOTIATE_DATAGRAM: int
FLAG_NEGOTIATE_SEAL: int
FLAG_NEGOTIATE_SIGN: int
FLAG_REQUEST_TARGET: int
FLAG_NEGOTIATE_OEM: int
FLAG_NEGOTIATE_UNICODE: int
FLAG_TYPES: Incomplete
AV_END_OF_LIST: int
AV_NETBIOS_COMPUTER_NAME: int
AV_NETBIOS_DOMAIN_NAME: int
AV_DNS_COMPUTER_NAME: int
AV_DNS_DOMAIN_NAME: int
AV_DNS_TREE_NAME: int
AV_FLAGS: int
AV_TIMESTAMP: int
AV_SINGLE_HOST_DATA: int
AV_TARGET_NAME: int
AV_CHANNEL_BINDINGS: int
AV_TYPES: Incomplete
AV_FLAG_CONSTRAINED: int
AV_FLAG_INTEGRITY: int
AV_FLAG_TARGET_SPN_UNTRUSTED: int
AV_FLAG_TYPES: Incomplete

def pack_windows_version(debug: bool = False): ...
def unpack_windows_version(version_message): ...

class NtlmClient:
    client_config_flags: int
    exported_session_key: Incomplete
    negotiated_flags: Incomplete
    user_name: Incomplete
    user_domain: Incomplete
    no_lm_response_ntlm_v1: Incomplete
    client_blocked: bool
    client_block_exceptions: Incomplete
    client_require_128_bit_encryption: Incomplete
    max_life_time: Incomplete
    client_signing_key: Incomplete
    client_sealing_key: Incomplete
    sequence_number: Incomplete
    server_sealing_key: Incomplete
    server_signing_key: Incomplete
    integrity: bool
    replay_detect: bool
    sequence_detect: bool
    confidentiality: bool
    datagram: bool
    identity: bool
    client_supplied_target_name: Incomplete
    client_channel_binding_unhashed: Incomplete
    unverified_target_name: Incomplete
    server_challenge: Incomplete
    server_target_name: Incomplete
    server_target_info: Incomplete
    server_version: Incomplete
    server_av_netbios_computer_name: Incomplete
    server_av_netbios_domain_name: Incomplete
    server_av_dns_computer_name: Incomplete
    server_av_dns_domain_name: Incomplete
    server_av_dns_forest_name: Incomplete
    server_av_target_name: Incomplete
    server_av_flags: Incomplete
    server_av_timestamp: Incomplete
    server_av_single_host_data: Incomplete
    server_av_channel_bindings: Incomplete
    server_av_flag_constrained: Incomplete
    server_av_flag_integrity: Incomplete
    server_av_flag_target_spn_untrusted: Incomplete
    current_encoding: Incomplete
    client_challenge: Incomplete
    server_target_info_raw: Incomplete
    def __init__(self, domain, user_name, password) -> None: ...
    def get_client_flag(self, flag): ...
    def get_negotiated_flag(self, flag): ...
    def get_server_av_flag(self, flag): ...
    def set_client_flag(self, flags) -> None: ...
    def reset_client_flags(self) -> None: ...
    def unset_client_flag(self, flags) -> None: ...
    def create_negotiate_message(self): ...
    def parse_challenge_message(self, message): ...
    def create_authenticate_message(self): ...
    @staticmethod
    def pack_field(value, offset): ...
    @staticmethod
    def unpack_field(field_message): ...
    @staticmethod
    def unpack_av_info(info): ...
    @staticmethod
    def pack_av_info(avs): ...
    @staticmethod
    def pack_windows_timestamp(): ...
    def compute_nt_response(self): ...
    def ntowf_v2(self): ...
