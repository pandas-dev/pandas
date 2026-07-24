from hvac.api.auth_methods.approle import AppRole as AppRole
from hvac.api.auth_methods.aws import Aws as Aws
from hvac.api.auth_methods.azure import Azure as Azure
from hvac.api.auth_methods.cert import Cert as Cert
from hvac.api.auth_methods.gcp import Gcp as Gcp
from hvac.api.auth_methods.github import Github as Github
from hvac.api.auth_methods.jwt import JWT as JWT
from hvac.api.auth_methods.kubernetes import Kubernetes as Kubernetes
from hvac.api.auth_methods.ldap import Ldap as Ldap
from hvac.api.auth_methods.legacy_mfa import LegacyMfa as LegacyMfa
from hvac.api.auth_methods.oidc import OIDC as OIDC
from hvac.api.auth_methods.okta import Okta as Okta
from hvac.api.auth_methods.radius import Radius as Radius
from hvac.api.auth_methods.token import Token as Token
from hvac.api.auth_methods.userpass import Userpass as Userpass
from hvac.api.vault_api_base import VaultApiBase
from hvac.api.vault_api_category import VaultApiCategory

__all__ = (
    "AuthMethods",
    "AppRole",
    "Azure",
    "Gcp",
    "Github",
    "JWT",
    "Kubernetes",
    "Ldap",
    "Userpass",
    "LegacyMfa",
    "OIDC",
    "Okta",
    "Radius",
    "Token",
    "Aws",
    "Cert",
)

class AuthMethods(VaultApiCategory):
    implemented_classes: list[type[VaultApiBase]]
    unimplemented_classes: list[str]
