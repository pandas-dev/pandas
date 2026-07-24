from hvac.api.system_backend.audit import Audit as Audit
from hvac.api.system_backend.auth import Auth as Auth
from hvac.api.system_backend.capabilities import Capabilities as Capabilities
from hvac.api.system_backend.health import Health as Health
from hvac.api.system_backend.init import Init as Init
from hvac.api.system_backend.key import Key as Key
from hvac.api.system_backend.leader import Leader as Leader
from hvac.api.system_backend.lease import Lease as Lease
from hvac.api.system_backend.mount import Mount as Mount
from hvac.api.system_backend.namespace import Namespace as Namespace
from hvac.api.system_backend.policies import Policies as Policies
from hvac.api.system_backend.policy import Policy as Policy
from hvac.api.system_backend.quota import Quota as Quota
from hvac.api.system_backend.raft import Raft as Raft
from hvac.api.system_backend.seal import Seal as Seal
from hvac.api.system_backend.system_backend_mixin import SystemBackendMixin as SystemBackendMixin
from hvac.api.system_backend.wrapping import Wrapping as Wrapping
from hvac.api.vault_api_base import VaultApiBase
from hvac.api.vault_api_category import VaultApiCategory

__all__ = (
    "Audit",
    "Auth",
    "Capabilities",
    "Health",
    "Init",
    "Key",
    "Leader",
    "Lease",
    "Mount",
    "Namespace",
    "Policies",
    "Policy",
    "Quota",
    "Raft",
    "Seal",
    "SystemBackend",
    "SystemBackendMixin",
    "Wrapping",
)

class SystemBackend(
    VaultApiCategory,
    Audit,
    Auth,
    Capabilities,
    Health,
    Init,
    Key,
    Leader,
    Lease,
    Mount,
    Namespace,
    Policies,
    Policy,
    Quota,
    Raft,
    Seal,
    Wrapping,
):
    implemented_classes: list[type[VaultApiBase]]
    unimplemented_classes: list[str]
