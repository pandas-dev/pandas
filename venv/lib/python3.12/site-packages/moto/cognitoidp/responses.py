import json
from typing import Any

from moto.core.responses import TYPE_RESPONSE, ActionResult, BaseResponse, EmptyResult
from moto.utilities.utils import load_resource

from .exceptions import InvalidParameterException
from .models import (
    CognitoIdpBackend,
    RegionAgnosticBackend,
    UserStatus,
    cognitoidp_backends,
    find_account_region_by_value,
)


class CognitoIdpResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="cognito-idp")

    def _get_region_agnostic_backend(self) -> RegionAgnosticBackend:
        return RegionAgnosticBackend(self.current_account, self.region)

    @property
    def parameters(self) -> dict[str, Any]:  # type: ignore[misc]
        return json.loads(self.body)

    @property
    def backend(self) -> CognitoIdpBackend:
        return cognitoidp_backends[self.current_account][self.region]

    # User pool
    def create_user_pool(self) -> ActionResult:
        name = self.parameters.pop("PoolName")
        user_pool = self.backend.create_user_pool(name, self.parameters)
        return ActionResult({"UserPool": user_pool.to_json(extended=True)})

    def set_user_pool_mfa_config(self) -> ActionResult:
        user_pool_id = self._get_param("UserPoolId")
        sms_config = self._get_param("SmsMfaConfiguration", None)
        token_config = self._get_param("SoftwareTokenMfaConfiguration", None)
        mfa_config = self._get_param("MfaConfiguration")

        if mfa_config not in ["ON", "OFF", "OPTIONAL"]:
            raise InvalidParameterException(
                "[MfaConfiguration] must be one of 'ON', 'OFF', or 'OPTIONAL'."
            )

        if mfa_config in ["ON", "OPTIONAL"]:
            if sms_config is None and token_config is None:
                raise InvalidParameterException(
                    "At least one of [SmsMfaConfiguration] or [SoftwareTokenMfaConfiguration] must be provided."
                )
            if sms_config is not None:
                if "SmsConfiguration" not in sms_config:
                    raise InvalidParameterException(
                        "[SmsConfiguration] is a required member of [SoftwareTokenMfaConfiguration]."
                    )

        response = self.backend.set_user_pool_mfa_config(
            user_pool_id, sms_config, token_config, mfa_config
        )
        return ActionResult(response)

    def get_user_pool_mfa_config(self) -> ActionResult:
        user_pool_id = self._get_param("UserPoolId")
        response = self.backend.get_user_pool_mfa_config(user_pool_id)
        return ActionResult(response)

    def list_user_pools(self) -> ActionResult:
        max_results = self._get_param("MaxResults")
        next_token = self._get_param("NextToken")
        user_pools, next_token = self.backend.list_user_pools(
            max_results=max_results, next_token=next_token
        )
        response: dict[str, Any] = {
            "UserPools": [user_pool.to_json() for user_pool in user_pools]
        }
        if next_token:
            response["NextToken"] = str(next_token)
        return ActionResult(response)

    def describe_user_pool(self) -> ActionResult:
        user_pool_id = self._get_param("UserPoolId")
        user_pool = self.backend.describe_user_pool(user_pool_id)
        return ActionResult({"UserPool": user_pool.to_json(extended=True)})

    def update_user_pool(self) -> None:
        user_pool_id = self._get_param("UserPoolId")
        self.backend.update_user_pool(user_pool_id, self.parameters)

    def delete_user_pool(self) -> ActionResult:
        user_pool_id = self._get_param("UserPoolId")
        self.backend.delete_user_pool(user_pool_id)
        return EmptyResult()

    # User pool domain
    def create_user_pool_domain(self) -> ActionResult:
        domain = self._get_param("Domain")
        user_pool_id = self._get_param("UserPoolId")
        custom_domain_config = self._get_param("CustomDomainConfig")
        user_pool_domain = self.backend.create_user_pool_domain(
            user_pool_id, domain, custom_domain_config
        )
        domain_description = user_pool_domain.to_json(extended=False)
        return ActionResult(domain_description)

    def describe_user_pool_domain(self) -> ActionResult:
        domain = self._get_param("Domain")
        user_pool_domain = self.backend.describe_user_pool_domain(domain)
        domain_description: dict[str, Any] = {}
        if user_pool_domain:
            domain_description = user_pool_domain.to_json()

        return ActionResult({"DomainDescription": domain_description})

    def delete_user_pool_domain(self) -> ActionResult:
        domain = self._get_param("Domain")
        self.backend.delete_user_pool_domain(domain)
        return EmptyResult()

    def update_user_pool_domain(self) -> ActionResult:
        domain = self._get_param("Domain")
        custom_domain_config = self._get_param("CustomDomainConfig")
        user_pool_domain = self.backend.update_user_pool_domain(
            domain, custom_domain_config
        )
        domain_description = user_pool_domain.to_json(extended=False)
        return ActionResult(domain_description)

    # User pool client
    def create_user_pool_client(self) -> ActionResult:
        user_pool_id = self.parameters.pop("UserPoolId")
        generate_secret = self.parameters.pop("GenerateSecret", False)
        user_pool_client = self.backend.create_user_pool_client(
            user_pool_id, generate_secret, self.parameters
        )
        return ActionResult({"UserPoolClient": user_pool_client.to_json(extended=True)})

    def list_user_pool_clients(self) -> ActionResult:
        user_pool_id = self._get_param("UserPoolId")
        max_results = self._get_param("MaxResults")
        next_token = self._get_param("NextToken")
        user_pool_clients, next_token = self.backend.list_user_pool_clients(
            user_pool_id, max_results=max_results, next_token=next_token
        )
        response: dict[str, Any] = {
            "UserPoolClients": [
                user_pool_client.to_json() for user_pool_client in user_pool_clients
            ]
        }
        if next_token:
            response["NextToken"] = str(next_token)
        return ActionResult(response)

    def describe_user_pool_client(self) -> ActionResult:
        user_pool_id = self._get_param("UserPoolId")
        client_id = self._get_param("ClientId")
        user_pool_client = self.backend.describe_user_pool_client(
            user_pool_id, client_id
        )
        return ActionResult({"UserPoolClient": user_pool_client.to_json(extended=True)})

    def update_user_pool_client(self) -> ActionResult:
        user_pool_id = self.parameters.pop("UserPoolId")
        client_id = self.parameters.pop("ClientId")
        user_pool_client = self.backend.update_user_pool_client(
            user_pool_id, client_id, self.parameters
        )
        return ActionResult({"UserPoolClient": user_pool_client.to_json(extended=True)})

    def delete_user_pool_client(self) -> ActionResult:
        user_pool_id = self._get_param("UserPoolId")
        client_id = self._get_param("ClientId")
        self.backend.delete_user_pool_client(user_pool_id, client_id)
        return EmptyResult()

    # Identity provider
    def create_identity_provider(self) -> ActionResult:
        user_pool_id = self._get_param("UserPoolId")
        name = self.parameters.pop("ProviderName")
        identity_provider = self.backend.create_identity_provider(
            user_pool_id, name, self.parameters
        )
        return ActionResult(
            {"IdentityProvider": identity_provider.to_json(extended=True)}
        )

    def list_identity_providers(self) -> ActionResult:
        user_pool_id = self._get_param("UserPoolId")
        max_results = self._get_param("MaxResults")
        next_token = self._get_param("NextToken")
        identity_providers, next_token = self.backend.list_identity_providers(
            user_pool_id, max_results=max_results, next_token=next_token
        )
        response: dict[str, Any] = {
            "Providers": [
                identity_provider.to_json() for identity_provider in identity_providers
            ]
        }
        if next_token:
            response["NextToken"] = str(next_token)
        return ActionResult(response)

    def describe_identity_provider(self) -> ActionResult:
        user_pool_id = self._get_param("UserPoolId")
        name = self._get_param("ProviderName")
        identity_provider = self.backend.describe_identity_provider(user_pool_id, name)
        return ActionResult(
            {"IdentityProvider": identity_provider.to_json(extended=True)}
        )

    def update_identity_provider(self) -> ActionResult:
        user_pool_id = self._get_param("UserPoolId")
        name = self._get_param("ProviderName")
        identity_provider = self.backend.update_identity_provider(
            user_pool_id, name, self.parameters
        )
        return ActionResult(
            {"IdentityProvider": identity_provider.to_json(extended=True)}
        )

    def delete_identity_provider(self) -> ActionResult:
        user_pool_id = self._get_param("UserPoolId")
        name = self._get_param("ProviderName")
        self.backend.delete_identity_provider(user_pool_id, name)
        return EmptyResult()

    # Group
    def create_group(self) -> ActionResult:
        group_name = self._get_param("GroupName")
        user_pool_id = self._get_param("UserPoolId")
        description = self._get_param("Description")
        role_arn = self._get_param("RoleArn")
        precedence = self._get_param("Precedence")

        group = self.backend.create_group(
            user_pool_id, group_name, description, role_arn, precedence
        )

        return ActionResult({"Group": group.to_json()})

    def get_group(self) -> ActionResult:
        group_name = self._get_param("GroupName")
        user_pool_id = self._get_param("UserPoolId")
        group = self.backend.get_group(user_pool_id, group_name)
        return ActionResult({"Group": group.to_json()})

    def list_groups(self) -> ActionResult:
        user_pool_id = self._get_param("UserPoolId")
        limit = self._get_param("Limit")
        token = self._get_param("NextToken")
        groups, token = self.backend.list_groups(
            user_pool_id, limit=limit, next_token=token
        )
        response = {"Groups": [group.to_json() for group in groups], "NextToken": token}
        return ActionResult(response)

    def delete_group(self) -> ActionResult:
        group_name = self._get_param("GroupName")
        user_pool_id = self._get_param("UserPoolId")
        self.backend.delete_group(user_pool_id, group_name)
        return EmptyResult()

    def update_group(self) -> ActionResult:
        group_name = self._get_param("GroupName")
        user_pool_id = self._get_param("UserPoolId")
        description = self._get_param("Description")
        role_arn = self._get_param("RoleArn")
        precedence = self._get_param("Precedence")

        group = self.backend.update_group(
            user_pool_id, group_name, description, role_arn, precedence
        )

        return ActionResult({"Group": group.to_json()})

    def admin_add_user_to_group(self) -> ActionResult:
        user_pool_id = self._get_param("UserPoolId")
        username = self._get_param("Username")
        group_name = self._get_param("GroupName")

        self.backend.admin_add_user_to_group(user_pool_id, group_name, username)

        return EmptyResult()

    def list_users_in_group(self) -> ActionResult:
        user_pool_id = self._get_param("UserPoolId")
        group_name = self._get_param("GroupName")
        limit = self._get_param("Limit")
        token = self._get_param("NextToken")
        users, token = self.backend.list_users_in_group(
            user_pool_id, group_name, limit=limit, next_token=token
        )
        response = {
            "Users": [user.to_json(extended=True) for user in users],
            "NextToken": token,
        }
        return ActionResult(response)

    def admin_list_groups_for_user(self) -> ActionResult:
        username = self._get_param("Username")
        user_pool_id = self._get_param("UserPoolId")
        limit = self._get_param("Limit")
        token = self._get_param("NextToken")
        groups, token = self.backend.admin_list_groups_for_user(
            user_pool_id, username, limit=limit, next_token=token
        )
        response = {"Groups": [group.to_json() for group in groups], "NextToken": token}
        return ActionResult(response)

    def admin_remove_user_from_group(self) -> ActionResult:
        user_pool_id = self._get_param("UserPoolId")
        username = self._get_param("Username")
        group_name = self._get_param("GroupName")

        self.backend.admin_remove_user_from_group(user_pool_id, group_name, username)

        return EmptyResult()

    def admin_reset_user_password(self) -> ActionResult:
        user_pool_id = self._get_param("UserPoolId")
        username = self._get_param("Username")
        self.backend.admin_reset_user_password(user_pool_id, username)
        return EmptyResult()

    # User
    def admin_create_user(self) -> ActionResult:
        user_pool_id = self._get_param("UserPoolId")
        username = self._get_param("Username")
        message_action = self._get_param("MessageAction")
        temporary_password = self._get_param("TemporaryPassword")
        user = self.backend.admin_create_user(
            user_pool_id,
            username,
            message_action,
            temporary_password,
            self._get_param("UserAttributes", []),
        )

        return ActionResult({"User": user.to_json(extended=True)})

    def admin_confirm_sign_up(self) -> ActionResult:
        user_pool_id = self._get_param("UserPoolId")
        username = self._get_param("Username")
        self.backend.admin_confirm_sign_up(user_pool_id, username)
        return EmptyResult()

    def admin_get_user(self) -> ActionResult:
        user_pool_id = self._get_param("UserPoolId")
        username = self._get_param("Username")
        user = self.backend.admin_get_user(user_pool_id, username)
        return ActionResult(
            user.to_json(extended=True, attributes_key="UserAttributes")
        )

    def get_user(self) -> ActionResult:
        access_token = self._get_param("AccessToken")
        user = self._get_region_agnostic_backend().get_user(access_token=access_token)
        return ActionResult(
            user.to_json(extended=True, attributes_key="UserAttributes")
        )

    def list_users(self) -> ActionResult:
        user_pool_id = self._get_param("UserPoolId")
        limit = self._get_param("Limit")
        token = self._get_param("PaginationToken")
        filt = self._get_param("Filter")
        attributes_to_get = self._get_param("AttributesToGet")
        users, token = self.backend.list_users(
            user_pool_id, filt, limit=limit, pagination_token=token
        )
        response: dict[str, Any] = {
            "Users": [
                user.to_json(extended=True, attributes_to_get=attributes_to_get)
                for user in users
            ]
        }
        if token:
            response["PaginationToken"] = str(token)
        return ActionResult(response)

    def admin_disable_user(self) -> ActionResult:
        user_pool_id = self._get_param("UserPoolId")
        username = self._get_param("Username")
        self.backend.admin_disable_user(user_pool_id, username)
        return EmptyResult()

    def admin_enable_user(self) -> ActionResult:
        user_pool_id = self._get_param("UserPoolId")
        username = self._get_param("Username")
        self.backend.admin_enable_user(user_pool_id, username)
        return EmptyResult()

    def admin_delete_user(self) -> ActionResult:
        user_pool_id = self._get_param("UserPoolId")
        username = self._get_param("Username")
        self.backend.admin_delete_user(user_pool_id, username)
        return EmptyResult()

    def admin_initiate_auth(self) -> ActionResult:
        user_pool_id = self._get_param("UserPoolId")
        client_id = self._get_param("ClientId")
        auth_flow = self._get_param("AuthFlow")
        auth_parameters = self._get_param("AuthParameters")

        auth_result = self.backend.admin_initiate_auth(
            user_pool_id, client_id, auth_flow, auth_parameters
        )

        return ActionResult(auth_result)

    def admin_respond_to_auth_challenge(self) -> ActionResult:
        session = self._get_param("Session")
        client_id = self._get_param("ClientId")
        challenge_name = self._get_param("ChallengeName")
        challenge_responses = self._get_param("ChallengeResponses")
        backend = self._get_region_agnostic_backend()
        auth_result = backend.admin_respond_to_auth_challenge(
            session, client_id, challenge_name, challenge_responses
        )

        return ActionResult(auth_result)

    def respond_to_auth_challenge(self) -> ActionResult:
        session = self._get_param("Session")
        client_id = self._get_param("ClientId")
        challenge_name = self._get_param("ChallengeName")
        challenge_responses = self._get_param("ChallengeResponses")
        auth_result = self._get_region_agnostic_backend().respond_to_auth_challenge(
            session, client_id, challenge_name, challenge_responses
        )

        return ActionResult(auth_result)

    def forgot_password(self) -> ActionResult:
        client_id = self._get_param("ClientId")
        username = self._get_param("Username")
        account, region = find_account_region_by_value(
            "client_id", client_id, fallback=(self.current_account, self.region)
        )
        confirmation_code, response = cognitoidp_backends[account][
            region
        ].forgot_password(client_id, username)
        self.response_headers["x-moto-forgot-password-confirmation-code"] = (
            confirmation_code  # type: ignore[assignment]
        )
        return ActionResult(response)

    # This endpoint receives no authorization header, so if moto-server is listening
    # on localhost (doesn't get a region in the host header), it doesn't know what
    # region's backend should handle the traffic, and we use `find_region_by_value` to
    # solve that problem.
    def confirm_forgot_password(self) -> ActionResult:
        client_id = self._get_param("ClientId")
        username = self._get_param("Username")
        password = self._get_param("Password")
        confirmation_code = self._get_param("ConfirmationCode")
        account, region = find_account_region_by_value(
            "client_id", client_id, fallback=(self.current_account, self.region)
        )
        cognitoidp_backends[account][region].confirm_forgot_password(
            client_id, username, password, confirmation_code
        )
        return EmptyResult()

    # Ditto the comment on confirm_forgot_password.
    def change_password(self) -> ActionResult:
        access_token = self._get_param("AccessToken")
        previous_password = self._get_param("PreviousPassword")
        proposed_password = self._get_param("ProposedPassword")
        account, region = find_account_region_by_value(
            "access_token", access_token, fallback=(self.current_account, self.region)
        )
        cognitoidp_backends[account][region].change_password(
            access_token, previous_password, proposed_password
        )
        return EmptyResult()

    def admin_update_user_attributes(self) -> ActionResult:
        user_pool_id = self._get_param("UserPoolId")
        username = self._get_param("Username")
        attributes = self._get_param("UserAttributes")
        self.backend.admin_update_user_attributes(user_pool_id, username, attributes)
        return EmptyResult()

    def admin_delete_user_attributes(self) -> ActionResult:
        user_pool_id = self._get_param("UserPoolId")
        username = self._get_param("Username")
        attributes = self._get_param("UserAttributeNames")
        self.backend.admin_delete_user_attributes(user_pool_id, username, attributes)
        return EmptyResult()

    def admin_user_global_sign_out(self) -> ActionResult:
        user_pool_id = self._get_param("UserPoolId")
        username = self._get_param("Username")
        self.backend.admin_user_global_sign_out(user_pool_id, username)
        return EmptyResult()

    def global_sign_out(self) -> ActionResult:
        access_token = self._get_param("AccessToken")
        self.backend.global_sign_out(access_token)
        return EmptyResult()

    # Resource Server
    def create_resource_server(self) -> ActionResult:
        user_pool_id = self._get_param("UserPoolId")
        identifier = self._get_param("Identifier")
        name = self._get_param("Name")
        scopes = self._get_param("Scopes")
        resource_server = self.backend.create_resource_server(
            user_pool_id, identifier, name, scopes
        )
        return ActionResult({"ResourceServer": resource_server.to_json()})

    def describe_resource_server(self) -> ActionResult:
        user_pool_id = self._get_param("UserPoolId")
        identifier = self._get_param("Identifier")
        resource_server = self.backend.describe_resource_server(
            user_pool_id, identifier
        )
        return ActionResult({"ResourceServer": resource_server.to_json()})

    def list_resource_servers(self) -> ActionResult:
        max_results = self._get_param("MaxResults")
        next_token = self._get_param("NextToken")
        user_pool_id = self._get_param("UserPoolId")
        resource_servers, next_token = self.backend.list_resource_servers(
            user_pool_id, max_results=max_results, next_token=next_token
        )
        response: dict[str, Any] = {
            "ResourceServers": [
                resource_server.to_json() for resource_server in resource_servers
            ]
        }
        if next_token:
            response["NextToken"] = str(next_token)
        return ActionResult(response)

    def sign_up(self) -> ActionResult:
        client_id = self._get_param("ClientId")
        username = self._get_param("Username")
        password = self._get_param("Password")
        user, code_delivery_details = self._get_region_agnostic_backend().sign_up(
            client_id=client_id,
            username=username,
            password=password,
            attributes=self._get_param("UserAttributes", []),
        )
        response = {
            "UserConfirmed": user.status == UserStatus["CONFIRMED"],
            "UserSub": user.id,
        }
        if code_delivery_details:
            response["CodeDeliveryDetails"] = code_delivery_details
        return ActionResult(response)

    def confirm_sign_up(self) -> ActionResult:
        client_id = self._get_param("ClientId")
        username = self._get_param("Username")
        self._get_region_agnostic_backend().confirm_sign_up(
            client_id=client_id, username=username
        )
        return EmptyResult()

    def initiate_auth(self) -> ActionResult:
        client_id = self._get_param("ClientId")
        auth_flow = self._get_param("AuthFlow")
        auth_parameters = self._get_param("AuthParameters")

        auth_result = self._get_region_agnostic_backend().initiate_auth(
            client_id, auth_flow, auth_parameters
        )

        return ActionResult(auth_result)

    def associate_software_token(self) -> ActionResult:
        access_token = self._get_param("AccessToken")
        session = self._get_param("Session")
        result = self._get_region_agnostic_backend().associate_software_token(
            access_token, session
        )
        return ActionResult(result)

    def verify_software_token(self) -> ActionResult:
        access_token = self._get_param("AccessToken")
        session = self._get_param("Session")
        result = self._get_region_agnostic_backend().verify_software_token(
            access_token, session
        )
        return ActionResult(result)

    def set_user_mfa_preference(self) -> ActionResult:
        access_token = self._get_param("AccessToken")
        software_token_mfa_settings = self._get_param("SoftwareTokenMfaSettings")
        sms_mfa_settings = self._get_param("SMSMfaSettings")
        self._get_region_agnostic_backend().set_user_mfa_preference(
            access_token, software_token_mfa_settings, sms_mfa_settings
        )
        return EmptyResult()

    def admin_set_user_mfa_preference(self) -> ActionResult:
        user_pool_id = self._get_param("UserPoolId")
        username = self._get_param("Username")
        software_token_mfa_settings = self._get_param("SoftwareTokenMfaSettings")
        sms_mfa_settings = self._get_param("SMSMfaSettings")
        self.backend.admin_set_user_mfa_preference(
            user_pool_id, username, software_token_mfa_settings, sms_mfa_settings
        )
        return EmptyResult()

    def admin_set_user_password(self) -> ActionResult:
        user_pool_id = self._get_param("UserPoolId")
        username = self._get_param("Username")
        password = self._get_param("Password")
        permanent = self._get_param("Permanent")
        self.backend.admin_set_user_password(
            user_pool_id, username, password, permanent
        )
        return EmptyResult()

    def add_custom_attributes(self) -> ActionResult:
        user_pool_id = self._get_param("UserPoolId")
        custom_attributes = self._get_param("CustomAttributes")
        self.backend.add_custom_attributes(user_pool_id, custom_attributes)
        return EmptyResult()

    def update_user_attributes(self) -> ActionResult:
        access_token = self._get_param("AccessToken")
        attributes = self._get_param("UserAttributes")
        self._get_region_agnostic_backend().update_user_attributes(
            access_token, attributes
        )
        return EmptyResult()


class CognitoIdpJsonWebKeyResponse(BaseResponse):
    json_web_key = json.dumps(
        load_resource(__name__, "resources/jwks-public.json")
    ).encode("utf-8")

    def __init__(self) -> None:
        super().__init__(service_name="cognito-idp")

    @staticmethod
    def serve_json_web_key(*args) -> TYPE_RESPONSE:  # type: ignore
        return (
            200,
            {"Content-Type": "application/json"},
            CognitoIdpJsonWebKeyResponse.json_web_key,
        )
