from types import TracebackType
from typing_extensions import Self

from auth0.rest import RestClientOptions

from .actions import Actions
from .attack_protection import AttackProtection
from .blacklists import Blacklists
from .branding import Branding
from .client_credentials import ClientCredentials
from .client_grants import ClientGrants
from .clients import Clients
from .connections import Connections
from .custom_domains import CustomDomains
from .device_credentials import DeviceCredentials
from .email_templates import EmailTemplates
from .emails import Emails
from .grants import Grants
from .guardian import Guardian
from .hooks import Hooks
from .jobs import Jobs
from .log_streams import LogStreams
from .logs import Logs
from .organizations import Organizations
from .prompts import Prompts
from .resource_servers import ResourceServers
from .roles import Roles
from .rules import Rules
from .rules_configs import RulesConfigs
from .stats import Stats
from .tenants import Tenants
from .tickets import Tickets
from .user_blocks import UserBlocks
from .users import Users
from .users_by_email import UsersByEmail

class AsyncAuth0:
    def __init__(self, domain: str, token: str, rest_options: RestClientOptions | None = None) -> None: ...
    def set_session(self, session) -> None: ...
    async def __aenter__(self) -> Self: ...
    async def __aexit__(
        self, exc_type: type[BaseException] | None, exc_val: BaseException | None, exc_tb: TracebackType | None
    ) -> None: ...

    # Same attributes as Auth0
    # See note in stubs/auth0-python/@tests/stubtest_allowlist.txt about _async methods
    actions: Actions
    attack_protection: AttackProtection
    blacklists: Blacklists
    branding: Branding
    client_credentials: ClientCredentials
    client_grants: ClientGrants
    clients: Clients
    connections: Connections
    custom_domains: CustomDomains
    device_credentials: DeviceCredentials
    email_templates: EmailTemplates
    emails: Emails
    grants: Grants
    guardian: Guardian
    hooks: Hooks
    jobs: Jobs
    log_streams: LogStreams
    logs: Logs
    organizations: Organizations
    prompts: Prompts
    resource_servers: ResourceServers
    roles: Roles
    rules_configs: RulesConfigs
    rules: Rules
    stats: Stats
    tenants: Tenants
    tickets: Tickets
    user_blocks: UserBlocks
    users_by_email: UsersByEmail
    users: Users
