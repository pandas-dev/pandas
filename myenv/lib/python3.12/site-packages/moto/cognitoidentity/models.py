import datetime
import json
import re
from collections import OrderedDict
from typing import Any, Dict, List, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import utcnow

from .exceptions import InvalidNameException, ResourceNotFoundError
from .utils import get_random_identity_id


class CognitoIdentityPool(BaseModel):
    def __init__(self, region: str, identity_pool_name: str, **kwargs: Any):
        self.identity_pool_name = identity_pool_name

        if not re.fullmatch(r"[\w\s+=,.@-]+", identity_pool_name):
            raise InvalidNameException(identity_pool_name)

        self.allow_unauthenticated_identities = kwargs.get(
            "allow_unauthenticated_identities", ""
        )
        self.supported_login_providers = kwargs.get("supported_login_providers", {})
        self.developer_provider_name = kwargs.get("developer_provider_name", "")
        self.open_id_connect_provider_arns = kwargs.get(
            "open_id_connect_provider_arns", []
        )
        self.cognito_identity_providers = kwargs.get("cognito_identity_providers", [])
        self.saml_provider_arns = kwargs.get("saml_provider_arns", [])

        self.identity_pool_id = get_random_identity_id(region)
        self.creation_time = utcnow()

        self.tags = kwargs.get("tags") or {}

    def to_json(self) -> str:
        return json.dumps(
            {
                "IdentityPoolId": self.identity_pool_id,
                "IdentityPoolName": self.identity_pool_name,
                "AllowUnauthenticatedIdentities": self.allow_unauthenticated_identities,
                "SupportedLoginProviders": self.supported_login_providers,
                "DeveloperProviderName": self.developer_provider_name,
                "OpenIdConnectProviderARNs": self.open_id_connect_provider_arns,
                "CognitoIdentityProviders": self.cognito_identity_providers,
                "SamlProviderARNs": self.saml_provider_arns,
                "IdentityPoolTags": self.tags,
            }
        )


class CognitoIdentityBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.identity_pools: Dict[str, CognitoIdentityPool] = OrderedDict()
        self.pools_identities: Dict[str, Dict[str, Any]] = {}

    def describe_identity_pool(self, identity_pool_id: str) -> str:
        identity_pool = self.identity_pools.get(identity_pool_id, None)

        if not identity_pool:
            raise ResourceNotFoundError(identity_pool_id)

        return identity_pool.to_json()

    def create_identity_pool(
        self,
        identity_pool_name: str,
        allow_unauthenticated_identities: bool,
        supported_login_providers: Dict[str, str],
        developer_provider_name: str,
        open_id_connect_provider_arns: List[str],
        cognito_identity_providers: List[Dict[str, Any]],
        saml_provider_arns: List[str],
        tags: Dict[str, str],
    ) -> str:
        new_identity = CognitoIdentityPool(
            self.region_name,
            identity_pool_name,
            allow_unauthenticated_identities=allow_unauthenticated_identities,
            supported_login_providers=supported_login_providers,
            developer_provider_name=developer_provider_name,
            open_id_connect_provider_arns=open_id_connect_provider_arns,
            cognito_identity_providers=cognito_identity_providers,
            saml_provider_arns=saml_provider_arns,
            tags=tags,
        )
        self.identity_pools[new_identity.identity_pool_id] = new_identity
        self.pools_identities.update(
            {
                new_identity.identity_pool_id: {
                    "IdentityPoolId": new_identity.identity_pool_id,
                    "Identities": [],
                }
            }
        )
        return new_identity.to_json()

    def update_identity_pool(
        self,
        identity_pool_id: str,
        identity_pool_name: str,
        allow_unauthenticated: Optional[bool],
        login_providers: Optional[Dict[str, str]],
        provider_name: Optional[str],
        provider_arns: Optional[List[str]],
        identity_providers: Optional[List[Dict[str, Any]]],
        saml_providers: Optional[List[str]],
        tags: Optional[Dict[str, str]],
    ) -> str:
        """
        The AllowClassic-parameter has not yet been implemented
        """
        pool = self.identity_pools[identity_pool_id]
        pool.identity_pool_name = pool.identity_pool_name or identity_pool_name
        if allow_unauthenticated is not None:
            pool.allow_unauthenticated_identities = allow_unauthenticated
        if login_providers is not None:
            pool.supported_login_providers = login_providers
        if provider_name:
            pool.developer_provider_name = provider_name
        if provider_arns is not None:
            pool.open_id_connect_provider_arns = provider_arns
        if identity_providers is not None:
            pool.cognito_identity_providers = identity_providers
        if saml_providers is not None:
            pool.saml_provider_arns = saml_providers
        if tags:
            pool.tags = tags

        return pool.to_json()

    def get_id(self, identity_pool_id: str) -> str:
        identity_id = {"IdentityId": get_random_identity_id(self.region_name)}
        self.pools_identities[identity_pool_id]["Identities"].append(identity_id)
        return json.dumps(identity_id)

    def get_credentials_for_identity(self, identity_id: str) -> str:
        duration = 90
        now = utcnow()
        expiration = now + datetime.timedelta(seconds=duration)
        return json.dumps(
            {
                "Credentials": {
                    "AccessKeyId": "TESTACCESSKEY12345",
                    "Expiration": expiration.timestamp(),
                    "SecretKey": "ABCSECRETKEY",
                    "SessionToken": "ABC12345",
                },
                "IdentityId": identity_id,
            }
        )

    def get_open_id_token_for_developer_identity(self, identity_id: str) -> str:
        return json.dumps(
            {
                "IdentityId": identity_id,
                "Token": get_random_identity_id(self.region_name),
            }
        )

    def get_open_id_token(self, identity_id: str) -> str:
        return json.dumps(
            {
                "IdentityId": identity_id,
                "Token": get_random_identity_id(self.region_name),
            }
        )

    def list_identities(self, identity_pool_id: str) -> str:
        """
        The MaxResults-parameter has not yet been implemented
        """
        return json.dumps(self.pools_identities[identity_pool_id])

    def list_identity_pools(self) -> str:
        """
        The MaxResults-parameter has not yet been implemented
        """
        return json.dumps(
            {
                "IdentityPools": [
                    json.loads(pool.to_json()) for pool in self.identity_pools.values()
                ]
            }
        )


cognitoidentity_backends = BackendDict(CognitoIdentityBackend, "cognito-identity")
