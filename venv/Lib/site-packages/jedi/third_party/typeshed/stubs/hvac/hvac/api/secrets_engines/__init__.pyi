from hvac.api.secrets_engines.active_directory import ActiveDirectory as ActiveDirectory
from hvac.api.secrets_engines.aws import Aws as Aws
from hvac.api.secrets_engines.azure import Azure as Azure
from hvac.api.secrets_engines.database import Database as Database
from hvac.api.secrets_engines.gcp import Gcp as Gcp
from hvac.api.secrets_engines.identity import Identity as Identity
from hvac.api.secrets_engines.kv import Kv as Kv
from hvac.api.secrets_engines.kv_v1 import KvV1 as KvV1
from hvac.api.secrets_engines.kv_v2 import KvV2 as KvV2
from hvac.api.secrets_engines.ldap import Ldap as Ldap
from hvac.api.secrets_engines.pki import Pki as Pki
from hvac.api.secrets_engines.rabbitmq import RabbitMQ as RabbitMQ
from hvac.api.secrets_engines.ssh import Ssh as Ssh
from hvac.api.secrets_engines.transform import Transform as Transform
from hvac.api.secrets_engines.transit import Transit as Transit
from hvac.api.vault_api_base import VaultApiBase
from hvac.api.vault_api_category import VaultApiCategory

__all__ = (
    "Aws",
    "Azure",
    "Gcp",
    "ActiveDirectory",
    "Identity",
    "Kv",
    "KvV1",
    "KvV2",
    "Ldap",
    "Pki",
    "Transform",
    "Transit",
    "SecretsEngines",
    "Database",
    "RabbitMQ",
    "Ssh",
)

class SecretsEngines(VaultApiCategory):
    implemented_classes: list[type[VaultApiBase]]
    unimplemented_classes: list[str]
