from _typeshed import Incomplete

class ReverseDnsSetting:
    OFF: Incomplete
    REQUIRE_RESOLVE_ALL_ADDRESSES: Incomplete
    REQUIRE_RESOLVE_IP_ADDRESSES_ONLY: Incomplete
    OPTIONAL_RESOLVE_ALL_ADDRESSES: Incomplete
    OPTIONAL_RESOLVE_IP_ADDRESSES_ONLY: Incomplete
    SUPPORTED_VALUES: Incomplete

def get_hostname_by_addr(addr, success_required: bool = True): ...
def is_ip_addr(addr): ...
