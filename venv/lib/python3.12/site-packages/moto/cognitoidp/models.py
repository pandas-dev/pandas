import enum
import re
import time
import typing
from collections import OrderedDict
from typing import Any, Dict, List, Optional, Set, Tuple

from joserfc import jwk, jwt

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import utcnow
from moto.moto_api._internal import mock_random as random
from moto.utilities.paginator import paginate
from moto.utilities.utils import get_partition, load_resource, md5_hash

from ..settings import (
    get_cognito_idp_user_pool_client_id_strategy,
    get_cognito_idp_user_pool_id_strategy,
)
from .exceptions import (
    AliasExistsException,
    ExpiredCodeException,
    GroupExistsException,
    InvalidParameterException,
    InvalidPasswordException,
    NotAuthorizedError,
    ResourceNotFoundError,
    UsernameExistsException,
    UserNotConfirmedException,
    UserNotFoundError,
)
from .utils import (
    PAGINATION_MODEL,
    check_secret_hash,
    expand_attrs,
    flatten_attrs,
    generate_id,
    validate_username_format,
)


class UserStatus(str, enum.Enum):
    FORCE_CHANGE_PASSWORD = "FORCE_CHANGE_PASSWORD"
    CONFIRMED = "CONFIRMED"
    UNCONFIRMED = "UNCONFIRMED"
    RESET_REQUIRED = "RESET_REQUIRED"


class AuthFlow(str, enum.Enum):
    # Order follows AWS' order
    ADMIN_NO_SRP_AUTH = "ADMIN_NO_SRP_AUTH"
    ADMIN_USER_PASSWORD_AUTH = "ADMIN_USER_PASSWORD_AUTH"
    USER_SRP_AUTH = "USER_SRP_AUTH"
    REFRESH_TOKEN_AUTH = "REFRESH_TOKEN_AUTH"
    REFRESH_TOKEN = "REFRESH_TOKEN"
    CUSTOM_AUTH = "CUSTOM_AUTH"
    USER_PASSWORD_AUTH = "USER_PASSWORD_AUTH"

    @classmethod
    def list(cls) -> List[str]:
        return [e.value for e in cls]


class CognitoIdpUserPoolAttribute(BaseModel):
    STANDARD_SCHEMA = {
        "sub": {
            "AttributeDataType": "String",
            "Mutable": False,
            "Required": True,
            "StringAttributeConstraints": {"MinLength": "1", "MaxLength": "2048"},
        },
        "name": {
            "AttributeDataType": "String",
            "Mutable": True,
            "Required": False,
            "StringAttributeConstraints": {"MinLength": "0", "MaxLength": "2048"},
        },
        "given_name": {
            "AttributeDataType": "String",
            "Mutable": True,
            "Required": False,
            "StringAttributeConstraints": {"MinLength": "0", "MaxLength": "2048"},
        },
        "family_name": {
            "AttributeDataType": "String",
            "Mutable": True,
            "Required": False,
            "StringAttributeConstraints": {"MinLength": "0", "MaxLength": "2048"},
        },
        "middle_name": {
            "AttributeDataType": "String",
            "Mutable": True,
            "Required": False,
            "StringAttributeConstraints": {"MinLength": "0", "MaxLength": "2048"},
        },
        "nickname": {
            "AttributeDataType": "String",
            "Mutable": True,
            "Required": False,
            "StringAttributeConstraints": {"MinLength": "0", "MaxLength": "2048"},
        },
        "preferred_username": {
            "AttributeDataType": "String",
            "Mutable": True,
            "Required": False,
            "StringAttributeConstraints": {"MinLength": "0", "MaxLength": "2048"},
        },
        "profile": {
            "AttributeDataType": "String",
            "Mutable": True,
            "Required": False,
            "StringAttributeConstraints": {"MinLength": "0", "MaxLength": "2048"},
        },
        "picture": {
            "AttributeDataType": "String",
            "Mutable": True,
            "Required": False,
            "StringAttributeConstraints": {"MinLength": "0", "MaxLength": "2048"},
        },
        "website": {
            "AttributeDataType": "String",
            "Mutable": True,
            "Required": False,
            "StringAttributeConstraints": {"MinLength": "0", "MaxLength": "2048"},
        },
        "email": {
            "AttributeDataType": "String",
            "Mutable": True,
            "Required": False,
            "StringAttributeConstraints": {"MinLength": "0", "MaxLength": "2048"},
        },
        "email_verified": {
            "AttributeDataType": "Boolean",
            "Mutable": True,
            "Required": False,
        },
        "gender": {
            "AttributeDataType": "String",
            "Mutable": True,
            "Required": False,
            "StringAttributeConstraints": {"MinLength": "0", "MaxLength": "2048"},
        },
        "birthdate": {
            "AttributeDataType": "String",
            "Mutable": True,
            "Required": False,
            "StringAttributeConstraints": {"MinLength": "10", "MaxLength": "10"},
        },
        "zoneinfo": {
            "AttributeDataType": "String",
            "Mutable": True,
            "Required": False,
            "StringAttributeConstraints": {"MinLength": "0", "MaxLength": "2048"},
        },
        "locale": {
            "AttributeDataType": "String",
            "Mutable": True,
            "Required": False,
            "StringAttributeConstraints": {"MinLength": "0", "MaxLength": "2048"},
        },
        "phone_number": {
            "AttributeDataType": "String",
            "Mutable": True,
            "Required": False,
            "StringAttributeConstraints": {"MinLength": "0", "MaxLength": "2048"},
        },
        "phone_number_verified": {
            "AttributeDataType": "Boolean",
            "Mutable": True,
            "Required": False,
        },
        "address": {
            "AttributeDataType": "String",
            "Mutable": True,
            "Required": False,
            "StringAttributeConstraints": {"MinLength": "0", "MaxLength": "2048"},
        },
        "updated_at": {
            "AttributeDataType": "Number",
            "Mutable": True,
            "Required": False,
            "NumberAttributeConstraints": {"MinValue": "0"},
        },
    }

    ATTRIBUTE_DATA_TYPES = {"Boolean", "DateTime", "String", "Number"}

    def __init__(self, name: str, custom: bool, schema: Dict[str, Any]):
        self.name = name
        self.custom = custom
        attribute_data_type = schema.get("AttributeDataType", None)
        if (
            attribute_data_type
            and attribute_data_type
            not in CognitoIdpUserPoolAttribute.ATTRIBUTE_DATA_TYPES
        ):
            raise InvalidParameterException(
                f"Validation error detected: Value '{attribute_data_type}' failed to satisfy constraint: Member must satisfy enum value set: [Boolean, Number, String, DateTime]"
            )

        if self.custom:
            self._init_custom(schema)
        else:
            self._init_standard(schema)

    def _init_custom(self, schema: Dict[str, Any]) -> None:
        self.name = "custom:" + self.name
        attribute_data_type = schema.get("AttributeDataType", None)
        if not attribute_data_type:
            raise InvalidParameterException(
                "Invalid AttributeDataType input, consider using the provided AttributeDataType enum."
            )
        self.data_type = attribute_data_type
        self.developer_only = schema.get("DeveloperOnlyAttribute", False)
        if self.developer_only:
            self.name = "dev:" + self.name
        self.mutable = schema.get("Mutable", True)
        if schema.get("Required", False):
            raise InvalidParameterException(
                "Required custom attributes are not supported currently."
            )
        self.required = False
        self._init_constraints(schema, None, show_empty_constraints=True)

    def _init_standard(self, schema: Dict[str, Any]) -> None:
        attribute_data_type = schema.get("AttributeDataType", None)
        default_attribute_data_type = CognitoIdpUserPoolAttribute.STANDARD_SCHEMA[
            self.name
        ]["AttributeDataType"]
        if attribute_data_type and attribute_data_type != default_attribute_data_type:
            raise InvalidParameterException(
                f"You can not change AttributeDataType or set developerOnlyAttribute for standard schema attribute {self.name}"
            )
        self.data_type = default_attribute_data_type
        if schema.get("DeveloperOnlyAttribute", False):
            raise InvalidParameterException(
                f"You can not change AttributeDataType or set developerOnlyAttribute for standard schema attribute {self.name}"
            )
        else:
            self.developer_only = False
        self.mutable = schema.get(
            "Mutable", CognitoIdpUserPoolAttribute.STANDARD_SCHEMA[self.name]["Mutable"]
        )
        self.required = schema.get(
            "Required",
            CognitoIdpUserPoolAttribute.STANDARD_SCHEMA[self.name]["Required"],
        )
        constraints_key = None
        if self.data_type == "Number":
            constraints_key = "NumberAttributeConstraints"
        elif self.data_type == "String":
            constraints_key = "StringAttributeConstraints"
        default_constraints = (
            None
            if not constraints_key
            else CognitoIdpUserPoolAttribute.STANDARD_SCHEMA[self.name][constraints_key]
        )
        self._init_constraints(schema, default_constraints)

    def _init_constraints(
        self,
        schema: Dict[str, Any],
        default_constraints: Any,
        show_empty_constraints: bool = False,
    ) -> None:
        def numeric_limit(num: Optional[str], constraint_type: str) -> Optional[int]:
            if not num:
                return  # type: ignore[return-value]
            parsed = None
            try:
                parsed = int(num)
            except ValueError:
                pass
            if parsed is None or parsed < 0:
                raise InvalidParameterException(
                    f"Invalid {constraint_type} for schema attribute {self.name}"
                )
            return parsed

        self.string_constraints: Optional[Dict[str, Any]] = (
            {} if show_empty_constraints else None
        )
        self.number_constraints = None

        if "AttributeDataType" in schema:
            # Quirk - schema is set/validated only if AttributeDataType is specified
            if self.data_type == "String":
                string_constraints = schema.get(
                    "StringAttributeConstraints", default_constraints
                )
                if not string_constraints:
                    return
                min_len = numeric_limit(
                    string_constraints.get("MinLength", None),
                    "StringAttributeConstraints",
                )
                max_len = numeric_limit(
                    string_constraints.get("MaxLength", None),
                    "StringAttributeConstraints",
                )
                if (min_len and min_len > 2048) or (max_len and max_len > 2048):
                    raise InvalidParameterException(
                        f"user.{self.name}: String attributes cannot have a length of more than 2048"
                    )
                if min_len and max_len and min_len > max_len:
                    raise InvalidParameterException(
                        f"user.{self.name}: Max length cannot be less than min length."
                    )
                self.string_constraints = string_constraints
                self.number_constraints = None
            elif self.data_type == "Number":
                number_constraints = schema.get(
                    "NumberAttributeConstraints", default_constraints
                )
                if not number_constraints:
                    return
                # No limits on either min or max value
                min_val = numeric_limit(
                    number_constraints.get("MinValue", None),
                    "NumberAttributeConstraints",
                )
                max_val = numeric_limit(
                    number_constraints.get("MaxValue", None),
                    "NumberAttributeConstraints",
                )
                if min_val and max_val and min_val > max_val:
                    raise InvalidParameterException(
                        f"user.{self.name}: Max value cannot be less than min value."
                    )
                self.number_constraints = number_constraints
                self.string_constraints = None
            else:
                self.number_constraints = None
                self.string_constraints = None

    def to_json(self) -> Dict[str, Any]:
        return {
            "Name": self.name,
            "AttributeDataType": self.data_type,
            "DeveloperOnlyAttribute": self.developer_only,
            "Mutable": self.mutable,
            "Required": self.required,
            "NumberAttributeConstraints": self.number_constraints,
            "StringAttributeConstraints": self.string_constraints,
        }


DEFAULT_USER_POOL_CONFIG: Dict[str, Any] = {
    "Policies": {
        "PasswordPolicy": {
            "MinimumLength": 8,
            "RequireUppercase": True,
            "RequireLowercase": True,
            "RequireNumbers": True,
            "RequireSymbols": True,
            "TemporaryPasswordValidityDays": 7,
        }
    },
    "AdminCreateUserConfig": {
        "AllowAdminCreateUserOnly": False,
        "UnusedAccountValidityDays": 7,
        "InviteMessageTemplate": {
            "SMSMessage": "Your username is {username} and temporary password is {####}. ",
            "EmailMessage": "Your username is {username} and temporary password is {####}. ",
            "EmailSubject": "Your temporary password",
        },
    },
    "EmailConfiguration": {"EmailSendingAccount": "COGNITO_DEFAULT"},
    "VerificationMessageTemplate": {
        "SmsMessage": "Your verification code is {####}. ",
        "EmailMessage": "Your verification code is {####}. ",
        "EmailSubject": "Your verification code",
        "DefaultEmailOption": "CONFIRM_WITH_CODE",
    },
    "AccountRecoverySetting": {
        "RecoveryMechanisms": [
            {"Priority": 1, "Name": "verified_email"},
            {"Priority": 2, "Name": "verified_phone_number"},
        ]
    },
}


class CognitoIdpUserPool(BaseModel):
    MAX_ID_LENGTH = 55

    def __init__(
        self, account_id: str, region: str, name: str, extended_config: Dict[str, Any]
    ):
        self.account_id = account_id
        self.region = region

        user_pool_id = generate_id(
            get_cognito_idp_user_pool_id_strategy(), region, name, extended_config
        )
        self.id = f"{self.region}_{user_pool_id}"[: self.MAX_ID_LENGTH]
        self.arn = f"arn:{get_partition(region)}:cognito-idp:{region}:{account_id}:userpool/{self.id}"

        self.name = name
        self.status = None

        self.update_extended_config(extended_config)

        self.creation_date = utcnow()
        self.last_modified_date = utcnow()

        self.mfa_config = extended_config.get("MfaConfiguration") or "OFF"
        self.sms_mfa_config: Optional[Dict[str, Any]] = None
        self.token_mfa_config: Optional[Dict[str, bool]] = None

        self.schema_attributes = {}
        for schema in self.extended_config.pop("Schema", {}):
            attribute = CognitoIdpUserPoolAttribute(
                schema["Name"],
                schema["Name"] not in CognitoIdpUserPoolAttribute.STANDARD_SCHEMA,
                schema,
            )
            self.schema_attributes[attribute.name] = attribute
        # If we do not have custom attributes, use the standard schema
        if not self.schema_attributes:
            for (
                standard_attribute_name,
                standard_attribute_schema,
            ) in CognitoIdpUserPoolAttribute.STANDARD_SCHEMA.items():
                self.schema_attributes[standard_attribute_name] = (
                    CognitoIdpUserPoolAttribute(
                        standard_attribute_name, False, standard_attribute_schema
                    )
                )

        self.clients: Dict[str, CognitoIdpUserPoolClient] = OrderedDict()
        self.identity_providers: Dict[str, CognitoIdpIdentityProvider] = OrderedDict()
        self.groups: Dict[str, CognitoIdpGroup] = OrderedDict()
        self.users: Dict[str, CognitoIdpUser] = OrderedDict()
        self.resource_servers: Dict[str, CognitoResourceServer] = OrderedDict()
        self.refresh_tokens: Dict[str, Optional[Tuple[str, str, str]]] = {}
        self.access_tokens: Dict[str, Tuple[str, str]] = {}
        self.id_tokens: Dict[str, Tuple[str, str]] = {}

        jwks_file = load_resource(__name__, "resources/jwks-private.json")
        self.json_web_key = jwk.RSAKey.import_key(jwks_file)

    @property
    def backend(self) -> "CognitoIdpBackend":
        return cognitoidp_backends[self.account_id][self.region]

    @property
    def domain(self) -> Optional["CognitoIdpUserPoolDomain"]:
        return next(
            (
                upd
                for upd in self.backend.user_pool_domains.values()
                if upd.user_pool_id == self.id
            ),
            None,
        )

    def update_extended_config(self, extended_config: Dict[str, Any]) -> None:
        self.extended_config = DEFAULT_USER_POOL_CONFIG.copy()
        self.extended_config.update(extended_config or {})

        message_template = self.extended_config.get("VerificationMessageTemplate")
        if message_template and "SmsVerificationMessage" not in extended_config:
            self.extended_config["SmsVerificationMessage"] = message_template.get(
                "SmsMessage"
            )
        if message_template and "EmailVerificationSubject" not in extended_config:
            self.extended_config["EmailVerificationSubject"] = message_template.get(
                "EmailSubject"
            )
        if message_template and "EmailVerificationMessage" not in extended_config:
            self.extended_config["EmailVerificationMessage"] = message_template.get(
                "EmailMessage"
            )

    def _base_json(self) -> Dict[str, Any]:
        return {
            "Id": self.id,
            "Arn": self.arn,
            "Name": self.name,
            "Status": self.status,
            "CreationDate": self.creation_date,
            "LastModifiedDate": self.last_modified_date,
            "MfaConfiguration": self.mfa_config,
            "EstimatedNumberOfUsers": len(self.users),
        }

    def to_json(self, extended: bool = False) -> Dict[str, Any]:
        user_pool_json = self._base_json()
        if extended:
            user_pool_json.update(self.extended_config)
            user_pool_json.update(
                {
                    "SchemaAttributes": [
                        att.to_json() for att in self.schema_attributes.values()
                    ]
                }
            )
        else:
            user_pool_json["LambdaConfig"] = (
                self.extended_config.get("LambdaConfig") or {}
            )
        if self.domain:
            user_pool_json["Domain"] = self.domain.domain
        return user_pool_json

    def _get_user(self, username: str) -> "CognitoIdpUser":
        """Find a user within a user pool by Username or any UsernameAttributes
        (`email` or `phone_number` or both)"""
        if self.extended_config.get("UsernameAttributes"):
            attribute_types = self.extended_config["UsernameAttributes"]
            for user in self.users.values():
                if username in [
                    flatten_attrs(user.attributes).get(attribute_type)
                    for attribute_type in attribute_types
                ]:
                    return user

        return self.users.get(username)  # type: ignore[return-value]

    def create_jwt(
        self,
        client_id: str,
        username: str,
        token_use: str,
        expires_in: int = 60 * 60,
        extra_data: Optional[Dict[str, Any]] = None,
    ) -> Tuple[str, int]:
        now = int(time.time())
        payload = {
            "iss": f"https://cognito-idp.{self.region}.amazonaws.com/{self.id}",
            "sub": self._get_user(username).id,
            "client_id" if token_use == "access" else "aud": client_id,
            "token_use": token_use,
            "auth_time": now,
            "exp": now + expires_in,
            "jti": str(random.uuid4()),
        }
        username_is_email = "email" in self.extended_config.get(
            "UsernameAttributes", []
        )
        if token_use == "access":
            if username_is_email:
                payload["username"] = payload["sub"]
            else:
                payload["username"] = username
        if token_use == "id":
            if username_is_email:
                payload["cognito:username"] = payload["sub"]
                payload["email"] = username
            else:
                payload["cognito:username"] = username

        payload.update(extra_data or {})
        headers = {"kid": "dummy", "alg": "RS256"}  # KID as present in jwks-public.json

        return (
            jwt.encode(headers, payload, self.json_web_key),
            expires_in,
        )

    def add_custom_attributes(self, custom_attributes: List[Dict[str, str]]) -> None:
        attributes = []
        for attribute_schema in custom_attributes:
            base_name = attribute_schema["Name"]
            target_name = "custom:" + base_name
            if attribute_schema.get("DeveloperOnlyAttribute", False):
                target_name = "dev:" + target_name
            if target_name in self.schema_attributes:
                raise InvalidParameterException(
                    f"custom:{base_name}: Existing attribute already has name {target_name}."
                )
            attribute = CognitoIdpUserPoolAttribute(base_name, True, attribute_schema)
            attributes.append(attribute)
        for attribute in attributes:
            self.schema_attributes[attribute.name] = attribute

    def create_id_token(
        self, client_id: str, username: str, origin_jti: str
    ) -> Tuple[str, int]:
        """
        :returns: (id_token, expires_in)
        """
        extra_data = self.get_user_extra_data_by_client_id(client_id, username)
        extra_data["origin_jti"] = origin_jti
        user = self._get_user(username)
        for attr in user.attributes:
            if attr["Name"].startswith("custom:"):
                extra_data[attr["Name"]] = attr["Value"]
        if len(user.groups) > 0:
            extra_data["cognito:groups"] = [group.group_name for group in user.groups]
        id_token, expires_in = self.create_jwt(
            client_id, username, "id", extra_data=extra_data
        )
        self.id_tokens[id_token] = (client_id, username)
        return id_token, expires_in

    def create_refresh_token(self, client_id: str, username: str) -> Tuple[str, str]:
        """
        :returns: (refresh_token, origin_jti)
        """
        refresh_token = str(random.uuid4())
        origin_jti = str(random.uuid4())
        self.refresh_tokens[refresh_token] = (client_id, username, origin_jti)
        return refresh_token, origin_jti

    def create_access_token(
        self, client_id: str, username: str, origin_jti: str
    ) -> Tuple[str, int]:
        """
        :returns: (access_token, expires_in)
        """
        extra_data: Dict[str, Any] = {
            "origin_jti": origin_jti,
        }
        user = self._get_user(username)
        if len(user.groups) > 0:
            extra_data["cognito:groups"] = [group.group_name for group in user.groups]

        access_token, expires_in = self.create_jwt(
            client_id, username, "access", extra_data=extra_data
        )
        self.access_tokens[access_token] = (client_id, username)
        return access_token, expires_in

    def create_tokens_from_refresh_token(
        self, refresh_token: str
    ) -> Tuple[str, str, int]:
        res = self.refresh_tokens[refresh_token]
        if res is None:
            raise NotAuthorizedError(refresh_token)
        client_id, username, origin_jti = res
        if not username:
            raise NotAuthorizedError(refresh_token)

        access_token, expires_in = self.create_access_token(
            client_id, username, origin_jti=origin_jti
        )
        id_token, _ = self.create_id_token(client_id, username, origin_jti=origin_jti)
        return access_token, id_token, expires_in

    def get_user_extra_data_by_client_id(
        self, client_id: str, username: str
    ) -> Dict[str, Any]:
        extra_data = {}
        current_client = self.clients.get(client_id, None)
        if current_client:
            for readable_field in current_client.get_readable_fields():
                attribute = list(
                    filter(
                        lambda f: f["Name"] == readable_field,
                        self._get_user(username).attributes,
                    )
                )
                if len(attribute) > 0:
                    extra_data.update({attribute[0]["Name"]: attribute[0]["Value"]})
        return extra_data

    def sign_out(self, username: str) -> None:
        for token, token_tuple in list(self.refresh_tokens.items()):
            if token_tuple is None:
                continue
            _, logged_in_user, _ = token_tuple
            if username == logged_in_user:
                self.refresh_tokens[token] = None
        for access_token, (_, logged_in_user) in list(self.access_tokens.items()):
            if username == logged_in_user:
                self.access_tokens.pop(access_token)


class CognitoIdpUserPoolDomain(BaseModel):
    def __init__(
        self,
        user_pool_id: str,
        domain: str,
        custom_domain_config: Optional[Dict[str, Any]] = None,
    ):
        self.user_pool_id = user_pool_id
        self.domain = domain
        self.custom_domain_config = custom_domain_config or {}

    def _distribution_name(self) -> str:
        if self.custom_domain_config and "CertificateArn" in self.custom_domain_config:
            unique_hash = md5_hash(
                self.custom_domain_config["CertificateArn"].encode("utf-8")
            ).hexdigest()
            return f"{unique_hash[:16]}.cloudfront.net"
        unique_hash = md5_hash(self.user_pool_id.encode("utf-8")).hexdigest()
        return f"{unique_hash[:16]}.amazoncognito.com"

    def to_json(self, extended: bool = True) -> Dict[str, Any]:
        distribution = self._distribution_name()
        if extended:
            return {
                "UserPoolId": self.user_pool_id,
                "AWSAccountId": str(random.uuid4()),
                "CloudFrontDistribution": distribution,
                "Domain": self.domain,
                "S3Bucket": None,
                "Status": "ACTIVE",
                "Version": None,
            }
        else:
            return {"CloudFrontDomain": distribution}


class CognitoIdpUserPoolClient(BaseModel):
    MAX_ID_LENGTH = 26

    def __init__(
        self,
        user_pool_id: str,
        generate_secret: bool,
        extended_config: Optional[Dict[str, Any]],
    ):
        self.user_pool_id = user_pool_id
        self.id = generate_id(
            get_cognito_idp_user_pool_client_id_strategy(),
            user_pool_id,
            generate_secret,
            extended_config,
        )[: self.MAX_ID_LENGTH]
        self.secret = str(random.uuid4())
        self.generate_secret = generate_secret or False
        # Some default values - may be overridden by the user
        self.extended_config: Dict[str, Any] = {
            "AllowedOAuthFlowsUserPoolClient": False,
            "AuthSessionValidity": 3,
            "EnablePropagateAdditionalUserContextData": False,
            "EnableTokenRevocation": True,
            "RefreshTokenValidity": 30,
        }
        self.extended_config.update(extended_config or {})

    def _base_json(self) -> Dict[str, Any]:
        return {
            "ClientId": self.id,
            "ClientName": self.extended_config.get("ClientName"),
            "UserPoolId": self.user_pool_id,
        }

    def to_json(self, extended: bool = False) -> Dict[str, Any]:
        user_pool_client_json = self._base_json()
        if self.generate_secret:
            user_pool_client_json.update({"ClientSecret": self.secret})
        if extended:
            user_pool_client_json.update(self.extended_config)

        return user_pool_client_json

    def get_readable_fields(self) -> List[str]:
        return self.extended_config.get("ReadAttributes", [])


class CognitoIdpIdentityProvider(BaseModel):
    def __init__(self, name: str, extended_config: Optional[Dict[str, Any]]):
        self.name = name
        self.extended_config = extended_config or {}
        self.creation_date = utcnow()
        self.last_modified_date = utcnow()

        if "AttributeMapping" not in self.extended_config:
            self.extended_config["AttributeMapping"] = {"username": "sub"}

    def _base_json(self) -> Dict[str, Any]:
        return {
            "ProviderName": self.name,
            "ProviderType": self.extended_config.get("ProviderType"),
            "CreationDate": self.creation_date,
            "LastModifiedDate": self.last_modified_date,
        }

    def to_json(self, extended: bool = False) -> Dict[str, Any]:
        identity_provider_json = self._base_json()
        if extended:
            identity_provider_json.update(self.extended_config)

        return identity_provider_json


class CognitoIdpGroup(BaseModel):
    def __init__(
        self,
        user_pool_id: str,
        group_name: str,
        description: str,
        role_arn: str,
        precedence: int,
    ):
        self.user_pool_id = user_pool_id
        self.group_name = group_name
        self.description = description or ""
        self.role_arn = role_arn
        self.precedence = precedence
        self.last_modified_date = utcnow()
        self.creation_date = self.last_modified_date

        # Users who are members of this group.
        # Note that these links are bidirectional.
        self.users: Set[CognitoIdpUser] = set()

    def update(
        self,
        description: Optional[str],
        role_arn: Optional[str],
        precedence: Optional[int],
    ) -> None:
        if description is not None:
            self.description = description
        if role_arn is not None:
            self.role_arn = role_arn
        if precedence is not None:
            self.precedence = precedence
        self.last_modified_date = utcnow()

    def to_json(self) -> Dict[str, Any]:
        return {
            "GroupName": self.group_name,
            "UserPoolId": self.user_pool_id,
            "Description": self.description,
            "RoleArn": self.role_arn,
            "Precedence": self.precedence,
            "LastModifiedDate": self.last_modified_date,
            "CreationDate": self.creation_date,
        }


class CognitoIdpUser(BaseModel):
    def __init__(
        self,
        user_pool_id: str,
        username: Optional[str],
        password: Optional[str],
        status: str,
        attributes: List[Dict[str, str]],
    ):
        self.id = str(random.uuid4())
        self.user_pool_id = user_pool_id
        # Username is None when users sign up with an email or phone_number,
        # and should be given the value of the internal id generate (sub)
        self.username = username if username else self.id
        self.password = password
        self.status = status
        self.enabled = True
        self.attributes = attributes
        self.attribute_lookup = flatten_attrs(attributes)
        self.create_date = utcnow()
        self.last_modified_date = utcnow()
        self.sms_mfa_enabled = False
        self.software_token_mfa_enabled = False
        self.token_verified = False
        self.confirmation_code: Optional[str] = None
        self.preferred_mfa_setting: Optional[str] = None

        # Groups this user is a member of.
        # Note that these links are bidirectional.
        self.groups: Set[CognitoIdpGroup] = set()

        self.update_attributes([{"Name": "sub", "Value": self.id}])

    def _base_json(self) -> Dict[str, Any]:
        return {
            "UserPoolId": self.user_pool_id,
            "Username": self.username,
            "UserStatus": self.status,
            "UserCreateDate": self.create_date,
            "UserLastModifiedDate": self.last_modified_date,
        }

    # list_users brings back "Attributes" while admin_get_user brings back "UserAttributes".
    def to_json(
        self,
        extended: bool = False,
        attributes_key: str = "Attributes",
        attributes_to_get: Optional[List[str]] = None,
    ) -> Dict[str, Any]:
        user_mfa_setting_list = []
        if self.software_token_mfa_enabled:
            user_mfa_setting_list.append("SOFTWARE_TOKEN_MFA")
        if self.sms_mfa_enabled:
            user_mfa_setting_list.append("SMS_MFA")
        user_json = self._base_json()
        if extended:
            attrs = [
                attr
                for attr in self.attributes
                if not attributes_to_get or attr["Name"] in attributes_to_get
            ]
            user_json.update(
                {
                    "Enabled": self.enabled,
                    attributes_key: attrs,
                    "MFAOptions": [],
                    "UserMFASettingList": user_mfa_setting_list,
                    "PreferredMfaSetting": self.preferred_mfa_setting or "",
                }
            )

        return user_json

    def update_attributes(self, new_attributes: List[Dict[str, Any]]) -> None:
        flat_attributes = flatten_attrs(self.attributes)
        flat_attributes.update(flatten_attrs(new_attributes))
        self.attribute_lookup = flat_attributes
        self.attributes = expand_attrs(flat_attributes)
        self.last_modified_date = utcnow()

    def delete_attributes(self, attrs_to_delete: List[str]) -> None:
        flat_attributes = flatten_attrs(self.attributes)
        wrong_attrs = []
        for attr in attrs_to_delete:
            try:
                flat_attributes.pop(attr)
            except KeyError:
                wrong_attrs.append(attr)
        if wrong_attrs:
            raise InvalidParameterException(
                "Invalid user attributes: "
                + "\n".join(
                    [
                        f"user.{w}: Attribute does not exist in the schema."
                        for w in wrong_attrs
                    ]
                )
                + "\n"
            )
        self.attribute_lookup = flat_attributes
        self.attributes = expand_attrs(flat_attributes)
        self.last_modified_date = utcnow()


class CognitoResourceServer(BaseModel):
    def __init__(
        self,
        user_pool_id: str,
        identifier: str,
        name: str,
        scopes: List[Dict[str, str]],
    ):
        self.user_pool_id = user_pool_id
        self.identifier = identifier
        self.name = name
        self.scopes = scopes

    def to_json(self) -> Dict[str, Any]:
        res: Dict[str, Any] = {
            "UserPoolId": self.user_pool_id,
            "Identifier": self.identifier,
            "Name": self.name,
        }

        if self.scopes:
            res.update({"Scopes": self.scopes})

        return res


class CognitoIdpBackend(BaseBackend):
    """
    Moto mocks the JWK uris.
    If you're using decorators, you can retrieve this information by making a call to `https://cognito-idp.us-west-2.amazonaws.com/someuserpoolid/.well-known/jwks.json`.

    Call `http://localhost:5000/userpoolid/.well-known/jwks.json` instead of you're running Moto in ServerMode or Docker.
    Because Moto cannot determine this is a CognitoIDP-request based on the URL alone, you have to add an Authorization-header instead:
    `Authorization: AWS4-HMAC-SHA256 Credential=mock_access_key/20220524/us-east-1/cognito-idp/aws4_request, SignedHeaders=content-length;content-type;host;x-amz-date, Signature=asdf`

    In some cases, you need to have reproducible IDs for the user pool.
    For example, a single initialization before the start of integration tests.

    This behavior can be enabled by passing the environment variable: MOTO_COGNITO_IDP_USER_POOL_ID_STRATEGY=HASH.
    Passing MOTO_COGNITO_IDP_USER_POOL_CLIENT_ID_STRATEGY=HASH enables the same logic for user pool clients.
    """

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.user_pools: Dict[str, CognitoIdpUserPool] = OrderedDict()
        self.user_pool_domains: Dict[str, CognitoIdpUserPoolDomain] = OrderedDict()
        self.sessions: Dict[str, Tuple[str, CognitoIdpUserPool]] = {}

    # User pool
    def create_user_pool(
        self, name: str, extended_config: Dict[str, Any]
    ) -> CognitoIdpUserPool:
        user_pool = CognitoIdpUserPool(
            self.account_id, self.region_name, name, extended_config
        )
        self.user_pools[user_pool.id] = user_pool
        return user_pool

    def set_user_pool_mfa_config(
        self,
        user_pool_id: str,
        sms_config: Dict[str, Any],
        token_config: Dict[str, bool],
        mfa_config: str,
    ) -> Dict[str, Any]:
        user_pool = self.describe_user_pool(user_pool_id)
        user_pool.mfa_config = mfa_config
        user_pool.sms_mfa_config = sms_config
        user_pool.token_mfa_config = token_config

        return self.get_user_pool_mfa_config(user_pool_id)

    def get_user_pool_mfa_config(self, user_pool_id: str) -> Dict[str, Any]:
        user_pool = self.describe_user_pool(user_pool_id)

        return {
            "SmsMfaConfiguration": user_pool.sms_mfa_config,
            "SoftwareTokenMfaConfiguration": user_pool.token_mfa_config,
            "MfaConfiguration": user_pool.mfa_config,
        }

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_user_pools(self) -> List[CognitoIdpUserPool]:
        return list(self.user_pools.values())

    def describe_user_pool(self, user_pool_id: str) -> CognitoIdpUserPool:
        user_pool = self.user_pools.get(user_pool_id)
        if not user_pool:
            raise ResourceNotFoundError(f"User pool {user_pool_id} does not exist.")

        return user_pool

    def update_user_pool(
        self, user_pool_id: str, extended_config: Dict[str, Any]
    ) -> None:
        user_pool = self.describe_user_pool(user_pool_id)
        user_pool.update_extended_config(extended_config)

    def delete_user_pool(self, user_pool_id: str) -> None:
        self.describe_user_pool(user_pool_id)

        del self.user_pools[user_pool_id]

    # User pool domain
    def create_user_pool_domain(
        self,
        user_pool_id: str,
        domain: str,
        custom_domain_config: Optional[Dict[str, str]] = None,
    ) -> CognitoIdpUserPoolDomain:
        self.describe_user_pool(user_pool_id)

        user_pool_domain = CognitoIdpUserPoolDomain(
            user_pool_id, domain, custom_domain_config=custom_domain_config
        )
        self.user_pool_domains[domain] = user_pool_domain
        return user_pool_domain

    def describe_user_pool_domain(
        self, domain: str
    ) -> Optional[CognitoIdpUserPoolDomain]:
        if domain not in self.user_pool_domains:
            return None

        return self.user_pool_domains[domain]

    def delete_user_pool_domain(self, domain: str) -> None:
        if domain not in self.user_pool_domains:
            raise ResourceNotFoundError(domain)

        del self.user_pool_domains[domain]

    def update_user_pool_domain(
        self, domain: str, custom_domain_config: Dict[str, str]
    ) -> CognitoIdpUserPoolDomain:
        if domain not in self.user_pool_domains:
            raise ResourceNotFoundError(domain)

        user_pool_domain = self.user_pool_domains[domain]
        user_pool_domain.custom_domain_config = custom_domain_config
        return user_pool_domain

    # User pool client
    def create_user_pool_client(
        self, user_pool_id: str, generate_secret: bool, extended_config: Dict[str, str]
    ) -> CognitoIdpUserPoolClient:
        user_pool = self.describe_user_pool(user_pool_id)

        user_pool_client = CognitoIdpUserPoolClient(
            user_pool_id, generate_secret, extended_config
        )
        user_pool.clients[user_pool_client.id] = user_pool_client
        return user_pool_client

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_user_pool_clients(
        self, user_pool_id: str
    ) -> List[CognitoIdpUserPoolClient]:
        user_pool = self.describe_user_pool(user_pool_id)

        return list(user_pool.clients.values())

    def describe_user_pool_client(
        self, user_pool_id: str, client_id: str
    ) -> CognitoIdpUserPoolClient:
        user_pool = self.describe_user_pool(user_pool_id)

        client = user_pool.clients.get(client_id)
        if not client:
            raise ResourceNotFoundError(client_id)

        return client

    def update_user_pool_client(
        self, user_pool_id: str, client_id: str, extended_config: Dict[str, str]
    ) -> CognitoIdpUserPoolClient:
        user_pool = self.describe_user_pool(user_pool_id)

        client = user_pool.clients.get(client_id)
        if not client:
            raise ResourceNotFoundError(client_id)

        client.extended_config.update(extended_config)
        return client

    def delete_user_pool_client(self, user_pool_id: str, client_id: str) -> None:
        user_pool = self.describe_user_pool(user_pool_id)

        if client_id not in user_pool.clients:
            raise ResourceNotFoundError(client_id)

        del user_pool.clients[client_id]

    # Identity provider
    def create_identity_provider(
        self, user_pool_id: str, name: str, extended_config: Dict[str, str]
    ) -> CognitoIdpIdentityProvider:
        user_pool = self.describe_user_pool(user_pool_id)

        identity_provider = CognitoIdpIdentityProvider(name, extended_config)
        user_pool.identity_providers[name] = identity_provider
        return identity_provider

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_identity_providers(
        self, user_pool_id: str
    ) -> List[CognitoIdpIdentityProvider]:
        user_pool = self.describe_user_pool(user_pool_id)

        return list(user_pool.identity_providers.values())

    def describe_identity_provider(
        self, user_pool_id: str, name: str
    ) -> CognitoIdpIdentityProvider:
        user_pool = self.describe_user_pool(user_pool_id)

        identity_provider = user_pool.identity_providers.get(name)
        if not identity_provider:
            raise ResourceNotFoundError(name)

        return identity_provider

    def update_identity_provider(
        self, user_pool_id: str, name: str, extended_config: Dict[str, str]
    ) -> CognitoIdpIdentityProvider:
        user_pool = self.describe_user_pool(user_pool_id)

        identity_provider = user_pool.identity_providers.get(name)
        if not identity_provider:
            raise ResourceNotFoundError(name)

        identity_provider.extended_config.update(extended_config)

        return identity_provider

    def delete_identity_provider(self, user_pool_id: str, name: str) -> None:
        user_pool = self.describe_user_pool(user_pool_id)

        if name not in user_pool.identity_providers:
            raise ResourceNotFoundError(name)

        del user_pool.identity_providers[name]

    # Group
    def create_group(
        self,
        user_pool_id: str,
        group_name: str,
        description: str,
        role_arn: str,
        precedence: int,
    ) -> CognitoIdpGroup:
        user_pool = self.describe_user_pool(user_pool_id)

        group = CognitoIdpGroup(
            user_pool_id, group_name, description, role_arn, precedence
        )
        if group.group_name in user_pool.groups:
            raise GroupExistsException("A group with the name already exists")
        user_pool.groups[group.group_name] = group

        return group

    def get_group(self, user_pool_id: str, group_name: str) -> CognitoIdpGroup:
        user_pool = self.describe_user_pool(user_pool_id)

        if group_name not in user_pool.groups:
            raise ResourceNotFoundError(group_name)

        return user_pool.groups[group_name]

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_groups(self, user_pool_id: str) -> List[CognitoIdpGroup]:
        user_pool = self.describe_user_pool(user_pool_id)

        return list(user_pool.groups.values())

    def delete_group(self, user_pool_id: str, group_name: str) -> None:
        user_pool = self.describe_user_pool(user_pool_id)

        if group_name not in user_pool.groups:
            raise ResourceNotFoundError(group_name)

        group = user_pool.groups[group_name]
        for user in group.users:
            user.groups.remove(group)

        del user_pool.groups[group_name]

    def update_group(
        self,
        user_pool_id: str,
        group_name: str,
        description: str,
        role_arn: str,
        precedence: int,
    ) -> CognitoIdpGroup:
        group = self.get_group(user_pool_id, group_name)

        group.update(description, role_arn, precedence)

        return group

    def admin_add_user_to_group(
        self, user_pool_id: str, group_name: str, username: str
    ) -> None:
        group = self.get_group(user_pool_id, group_name)
        user = self.admin_get_user(user_pool_id, username)

        group.users.add(user)
        user.groups.add(group)

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_users_in_group(
        self, user_pool_id: str, group_name: str
    ) -> List[CognitoIdpUser]:
        user_pool = self.describe_user_pool(user_pool_id)
        group = self.get_group(user_pool_id, group_name)
        return list(filter(lambda user: user in group.users, user_pool.users.values()))

    @paginate(pagination_model=PAGINATION_MODEL)
    def admin_list_groups_for_user(
        self, user_pool_id: str, username: str
    ) -> List[CognitoIdpGroup]:
        user = self.admin_get_user(user_pool_id, username)
        return list(user.groups)

    def admin_remove_user_from_group(
        self, user_pool_id: str, group_name: str, username: str
    ) -> None:
        group = self.get_group(user_pool_id, group_name)
        user = self.admin_get_user(user_pool_id, username)

        group.users.discard(user)
        user.groups.discard(group)

    def admin_reset_user_password(self, user_pool_id: str, username: str) -> None:
        user = self.admin_get_user(user_pool_id, username)
        if not user.enabled:
            raise NotAuthorizedError("User is disabled")
        if user.status is UserStatus.RESET_REQUIRED:
            return
        if user.status is not UserStatus.CONFIRMED:
            raise NotAuthorizedError(
                "User password cannot be reset in the current state."
            )
        if (
            user.attribute_lookup.get("email_verified", "false") == "false"
            and user.attribute_lookup.get("phone_number_verified", "false") == "false"
        ):
            raise InvalidParameterException(
                "Cannot reset password for the user as there is no registered/verified email or phone_number"
            )
        user.status = UserStatus.RESET_REQUIRED

    # User
    def admin_create_user(
        self,
        user_pool_id: str,
        username: str,
        message_action: str,
        temporary_password: str,
        attributes: List[Dict[str, str]],
    ) -> CognitoIdpUser:
        user_pool = self.describe_user_pool(user_pool_id)

        if message_action and message_action == "RESEND":
            self.admin_get_user(user_pool_id, username)
        elif user_pool._get_user(username):
            raise UsernameExistsException(username)

        # UsernameAttributes are attributes (either `email` or `phone_number`
        # or both) than can be used in the place of a unique username. If the
        # user provides an email or phone number when signing up, the user pool
        # performs the following steps:
        # 1. populates the correct field (email, phone_number) with the value
        #    supplied for Username
        # 2. generates a persistent GUID for the user that will be returned as
        #    the value of `Username` in the `get-user` and `list-users`
        #    operations, as well as the value of `sub` in `IdToken` and
        #    `AccessToken`
        #
        # ref: https://docs.aws.amazon.com/cognito/latest/developerguide/user-pool-settings-attributes.html#user-pool-settings-aliases-settings
        has_username_attrs = user_pool.extended_config.get("UsernameAttributes")
        if has_username_attrs:
            username_attributes = user_pool.extended_config["UsernameAttributes"]
            # attribute_type should be one of `email`, `phone_number` or both
            for attribute_type in username_attributes:
                # check if provided username matches one of the attribute types in
                # `UsernameAttributes`
                if attribute_type in username_attributes and validate_username_format(
                    username, _format=attribute_type
                ):
                    # insert provided username into new user's attributes under the
                    # correct key
                    flattened_attrs = flatten_attrs(attributes or [])
                    flattened_attrs.update({attribute_type: username})
                    attributes = expand_attrs(flattened_attrs)

                    # once the username has been validated against a username attribute
                    # type, there is no need to attempt validation against the other
                    # type(s)
                    break

            # The provided username has not matched the required format for any
            # of the possible attributes
            else:
                raise InvalidParameterException(
                    "Username should be either an email or a phone number."
                )

        user = CognitoIdpUser(
            user_pool_id,
            # set username to None so that it will be default to the internal GUID
            # when them user gets created
            None if has_username_attrs else username,
            temporary_password,
            UserStatus.FORCE_CHANGE_PASSWORD,
            attributes,
        )

        user_pool.users[user.username] = user
        return user

    def admin_confirm_sign_up(self, user_pool_id: str, username: str) -> str:
        user = self.admin_get_user(user_pool_id, username)
        user.status = UserStatus["CONFIRMED"]
        return ""

    def admin_get_user(self, user_pool_id: str, username: str) -> CognitoIdpUser:
        user_pool = self.describe_user_pool(user_pool_id)

        user = user_pool._get_user(username)
        if not user:
            raise UserNotFoundError("User does not exist.")
        return user

    def get_user(self, access_token: str) -> CognitoIdpUser:
        for user_pool in self.user_pools.values():
            if access_token in user_pool.access_tokens:
                _, username = user_pool.access_tokens[access_token]
                user = self.admin_get_user(user_pool.id, username)
                if (
                    not user
                    or not user.enabled
                    or user.status is not UserStatus.CONFIRMED
                ):
                    raise NotAuthorizedError("username")
                return user
        raise NotAuthorizedError("Invalid token")

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_users(self, user_pool_id: str, filt: str) -> List[CognitoIdpUser]:
        user_pool = self.describe_user_pool(user_pool_id)
        users = list(user_pool.users.values())
        if filt:
            inherent_attributes: Dict[str, Any] = {
                "cognito:user_status": lambda u: u.status,
                "status": lambda u: "Enabled" if u.enabled else "Disabled",
                "username": lambda u: u.username,
            }
            comparisons: Dict[str, Any] = {
                "=": lambda x, y: x == y,
                "^=": lambda x, y: x.startswith(y),
            }
            allowed_attributes = [
                "username",
                "email",
                "phone_number",
                "name",
                "given_name",
                "family_name",
                "preferred_username",
                "cognito:user_status",
                "status",
                "sub",
            ]

            match = re.match(r"([\w:]+)\s*(=|\^=)\s*\"(.*)\"", filt)
            if match:
                name, op, value = match.groups()
            else:
                raise InvalidParameterException("Error while parsing filter")
            if name not in allowed_attributes:
                raise InvalidParameterException(f"Invalid search attribute: {name}")
            compare = comparisons[op]
            users = [
                user
                for user in users
                if [
                    attr
                    for attr in user.attributes
                    if attr["Name"] == name and compare(attr["Value"], value)
                ]
                or (
                    name in inherent_attributes
                    and compare(inherent_attributes[name](user), value)
                )
            ]
        return users

    def admin_disable_user(self, user_pool_id: str, username: str) -> None:
        user = self.admin_get_user(user_pool_id, username)
        user.enabled = False

    def admin_enable_user(self, user_pool_id: str, username: str) -> None:
        user = self.admin_get_user(user_pool_id, username)
        user.enabled = True

    def admin_delete_user(self, user_pool_id: str, username: str) -> None:
        user_pool = self.describe_user_pool(user_pool_id)
        user = self.admin_get_user(user_pool_id, username)

        for group in user.groups:
            group.users.remove(user)

        # use internal username
        del user_pool.users[user.username]

    def _log_user_in(
        self,
        user_pool: CognitoIdpUserPool,
        client: CognitoIdpUserPoolClient,
        username: str,
    ) -> Dict[str, Dict[str, Any]]:
        refresh_token, _ = user_pool.create_refresh_token(client.id, username)
        access_token, id_token, expires_in = user_pool.create_tokens_from_refresh_token(
            refresh_token
        )

        return {
            "ChallengeParameters": {},
            "AuthenticationResult": {
                "IdToken": id_token,
                "AccessToken": access_token,
                "RefreshToken": refresh_token,
                "ExpiresIn": expires_in,
                "TokenType": "Bearer",
            },
        }

    def _validate_auth_flow(
        self, auth_flow: str, valid_flows: typing.List[AuthFlow]
    ) -> AuthFlow:
        """validate auth_flow value and convert auth_flow to enum"""

        try:
            auth_flow = AuthFlow[auth_flow]
        except KeyError:
            raise InvalidParameterException(
                f"1 validation error detected: Value '{auth_flow}' at 'authFlow' failed to satisfy constraint: "
                f"Member must satisfy enum value set: "
                f"{AuthFlow.list()}"
            )

        if auth_flow not in valid_flows:
            raise InvalidParameterException("Initiate Auth method not supported")

        return auth_flow

    def admin_initiate_auth(
        self,
        user_pool_id: str,
        client_id: str,
        auth_flow: str,
        auth_parameters: Dict[str, str],
    ) -> Dict[str, Any]:
        admin_auth_flows = [
            AuthFlow.ADMIN_NO_SRP_AUTH,
            AuthFlow.ADMIN_USER_PASSWORD_AUTH,
            AuthFlow.REFRESH_TOKEN_AUTH,
            AuthFlow.REFRESH_TOKEN,
        ]
        auth_flow = self._validate_auth_flow(
            auth_flow=auth_flow, valid_flows=admin_auth_flows
        )

        user_pool = self.describe_user_pool(user_pool_id)

        client = user_pool.clients.get(client_id)
        if not client:
            raise ResourceNotFoundError(client_id)

        if auth_flow in (AuthFlow.ADMIN_USER_PASSWORD_AUTH, AuthFlow.ADMIN_NO_SRP_AUTH):
            username: str = auth_parameters.get("USERNAME")  # type: ignore[assignment]
            password: str = auth_parameters.get("PASSWORD")  # type: ignore[assignment]
            user = self.admin_get_user(user_pool_id, username)

            if not user.enabled:
                raise NotAuthorizedError("User is disabled.")

            if user.password != password:
                raise NotAuthorizedError(username)

            if user.status in [
                UserStatus.FORCE_CHANGE_PASSWORD,
                UserStatus.RESET_REQUIRED,
            ]:
                session = str(random.uuid4())
                self.sessions[session] = (username, user_pool)

                return {
                    "ChallengeName": "NEW_PASSWORD_REQUIRED",
                    "ChallengeParameters": {},
                    "Session": session,
                }

            if (
                user.software_token_mfa_enabled
                and user.preferred_mfa_setting == "SOFTWARE_TOKEN_MFA"
            ):
                session = str(random.uuid4())
                self.sessions[session] = (username, user_pool)

                return {
                    "ChallengeName": "SOFTWARE_TOKEN_MFA",
                    "ChallengeParameters": {},
                    "Session": session,
                }

            if user.sms_mfa_enabled and user.preferred_mfa_setting == "SMS_MFA":
                session = str(random.uuid4())
                self.sessions[session] = (username, user_pool)

                return {
                    "ChallengeName": "SMS_MFA",
                    "ChallengeParameters": {},
                    "Session": session,
                }

            return self._log_user_in(user_pool, client, username)
        elif auth_flow in (AuthFlow.REFRESH_TOKEN, AuthFlow.REFRESH_TOKEN_AUTH):
            refresh_token: str = auth_parameters.get("REFRESH_TOKEN")  # type: ignore[assignment]
            (
                access_token,
                id_token,
                expires_in,
            ) = user_pool.create_tokens_from_refresh_token(refresh_token)

            return {
                "AuthenticationResult": {
                    "IdToken": id_token,
                    "AccessToken": access_token,
                    "ExpiresIn": expires_in,
                    "TokenType": "Bearer",
                }
            }
        else:
            # We shouldn't get here due to enum validation of auth_flow
            return None  # type: ignore[return-value]

    def admin_respond_to_auth_challenge(
        self,
        session: str,
        client_id: str,
        challenge_name: str,
        challenge_responses: Dict[str, str],
    ) -> Dict[str, Any]:
        # Responds to an authentication challenge, as an administrator.
        # The only differences between this admin endpoint and public endpoint are not relevant and so we can safely call
        # the public endpoint to do the work:
        #  - The admin endpoint requires a user pool id along with a session; the public endpoint searches across all pools
        #  - ContextData is passed in; we don't use it

        return self.respond_to_auth_challenge(
            session, client_id, challenge_name, challenge_responses
        )

    def respond_to_auth_challenge(
        self,
        session: str,
        client_id: str,
        challenge_name: str,
        challenge_responses: Dict[str, str],
    ) -> Dict[str, Any]:
        if challenge_name == "PASSWORD_VERIFIER":
            session = challenge_responses.get("PASSWORD_CLAIM_SECRET_BLOCK")  # type: ignore[assignment]

        if session not in self.sessions:
            raise ResourceNotFoundError(session)
        _, user_pool = self.sessions[session]

        client = user_pool.clients.get(client_id)
        if not client:
            raise ResourceNotFoundError(client_id)

        if challenge_name == "NEW_PASSWORD_REQUIRED":
            username: str = challenge_responses.get("USERNAME")  # type: ignore[assignment]
            new_password = challenge_responses.get("NEW_PASSWORD")
            if not new_password:
                raise InvalidPasswordException()
            self._validate_password(user_pool.id, new_password)
            user = self.admin_get_user(user_pool.id, username)

            user.password = new_password
            user.status = UserStatus.CONFIRMED

            if user_pool.mfa_config == "ON":
                mfas_can_setup = []
                if user_pool.token_mfa_config == {"Enabled": True}:
                    mfas_can_setup.append("SOFTWARE_TOKEN_MFA")
                if user_pool.sms_mfa_config:
                    mfas_can_setup.append("SMS_MFA")

                if (
                    mfas_can_setup
                    and not user.software_token_mfa_enabled
                    and not user.sms_mfa_enabled
                ):
                    return {
                        "ChallengeName": "MFA_SETUP",
                        "ChallengeParameters": {"MFAS_CAN_SETUP": mfas_can_setup},
                        "Session": session,
                    }

            del self.sessions[session]
            return self._log_user_in(user_pool, client, username)
        elif challenge_name == "PASSWORD_VERIFIER":
            username: str = challenge_responses.get("USERNAME")  # type: ignore[no-redef]
            user = self.admin_get_user(user_pool.id, username)

            password_claim_signature = challenge_responses.get(
                "PASSWORD_CLAIM_SIGNATURE"
            )
            if not password_claim_signature:
                raise ResourceNotFoundError(password_claim_signature)
            password_claim_secret_block = challenge_responses.get(
                "PASSWORD_CLAIM_SECRET_BLOCK"
            )
            if not password_claim_secret_block:
                raise ResourceNotFoundError(password_claim_secret_block)
            timestamp = challenge_responses.get("TIMESTAMP")
            if not timestamp:
                raise ResourceNotFoundError(timestamp)

            if user.status == UserStatus.FORCE_CHANGE_PASSWORD:
                return {
                    "ChallengeName": "NEW_PASSWORD_REQUIRED",
                    "ChallengeParameters": {
                        "USERNAME": username,
                    },
                    "Session": session,
                }
            if user_pool.mfa_config == "ON" and not user.token_verified:
                mfas_can_setup = []
                if user_pool.token_mfa_config == {"Enabled": True}:
                    mfas_can_setup.append("SOFTWARE_TOKEN_MFA")
                if user_pool.sms_mfa_config:
                    mfas_can_setup.append("SMS_MFA")

                return {
                    "ChallengeName": "MFA_SETUP",
                    "ChallengeParameters": {"MFAS_CAN_SETUP": mfas_can_setup},
                    "Session": session,
                }
            if user.software_token_mfa_enabled or (
                user_pool.token_mfa_config == {"Enabled": True} and user.token_verified
            ):
                return {
                    "ChallengeName": "SOFTWARE_TOKEN_MFA",
                    "Session": session,
                    "ChallengeParameters": {},
                }

            if user.sms_mfa_enabled:
                return {
                    "ChallengeName": "SMS_MFA",
                    "Session": session,
                    "ChallengeParameters": {},
                }

            del self.sessions[session]
            return self._log_user_in(user_pool, client, username)
        elif challenge_name == "SOFTWARE_TOKEN_MFA" or challenge_name == "SMS_MFA":
            username: str = challenge_responses.get("USERNAME")  # type: ignore[no-redef]
            self.admin_get_user(user_pool.id, username)

            mfa_code = challenge_responses.get(f"{challenge_name}_CODE")
            if not mfa_code:
                raise ResourceNotFoundError(mfa_code)

            if client.generate_secret:
                secret_hash = challenge_responses.get("SECRET_HASH")
                if not check_secret_hash(
                    client.secret, client.id, username, secret_hash
                ):
                    raise NotAuthorizedError(secret_hash)

            del self.sessions[session]
            return self._log_user_in(user_pool, client, username)

        elif challenge_name == "MFA_SETUP":
            username, user_pool = self.sessions[session]
            return self._log_user_in(user_pool, client, username)

        else:
            return {}

    def confirm_forgot_password(
        self, client_id: str, username: str, password: str, confirmation_code: str
    ) -> None:
        for user_pool in self.user_pools.values():
            if client_id in user_pool.clients and user_pool._get_user(username):
                user = user_pool._get_user(username)
                if (
                    confirmation_code.startswith("moto-confirmation-code:")
                    and user.confirmation_code != confirmation_code
                ):
                    raise ExpiredCodeException(
                        "Invalid code provided, please request a code again."
                    )
                user.password = password
                user.confirmation_code = None
                break
        else:
            raise ResourceNotFoundError(client_id)

    def forgot_password(
        self, client_id: str, username: str
    ) -> Tuple[Optional[str], Dict[str, Any]]:
        """
        The ForgotPassword operation is partially broken in AWS. If the input is 100% correct it works fine.

        Otherwise you get semi-random garbage and HTTP 200 OK, for example:
        - recovery for username which is not registered in any cognito pool
        - recovery for username belonging to a different user pool than the client id is registered to
        - phone-based recovery for a user without phone_number / phone_number_verified attributes
        - same as above, but email / email_verified
        """
        for user_pool in self.user_pools.values():
            if client_id in user_pool.clients:
                recovery_settings = user_pool.extended_config["AccountRecoverySetting"]
                user = user_pool._get_user(username)
                break
        else:
            raise ResourceNotFoundError("Username/client id combination not found.")

        confirmation_code: Optional[str] = None
        if user:
            # An unfortunate bit of magic - confirmation_code is opt-in, as it's returned
            # via a "x-moto-forgot-password-confirmation-code" http header, which is not the AWS way (should be SES, SNS, Cognito built-in email)
            # Verification of user.confirmation_code vs received code will be performed only for codes
            # beginning with 'moto-confirmation-code' prefix. All other codes are considered VALID.
            confirmation_code = (
                f"moto-confirmation-code:{random.randint(100_000, 999_999)}"
            )
            user.confirmation_code = confirmation_code

        code_delivery_details = self._get_code_delivery_details(
            recovery_settings, user, username
        )
        return confirmation_code, {"CodeDeliveryDetails": code_delivery_details}

    def _get_code_delivery_details(
        self, recovery_settings: Any, user: Optional[CognitoIdpUser], username: str
    ) -> Dict[str, str]:
        selected_recovery = min(
            recovery_settings["RecoveryMechanisms"],
            key=lambda recovery_mechanism: recovery_mechanism["Priority"],
        )
        if selected_recovery["Name"] == "admin_only":
            raise NotAuthorizedError("Contact administrator to reset password.")
        if selected_recovery["Name"] == "verified_phone_number":
            number = "+*******9934"
            if user and "phone_number" in user.attribute_lookup:
                number = user.attribute_lookup["phone_number"]
            return {
                "Destination": number,
                "DeliveryMedium": "SMS",
                "AttributeName": "phone_number",
            }
        else:
            email = username + "@h***.com"
            if user and "email" in user.attribute_lookup:
                first, second = user.attribute_lookup["email"].split("@")
                email = f"{first[0]}***@{second[0]}***"
            return {
                "Destination": email,
                "DeliveryMedium": "EMAIL",
                "AttributeName": "email",
            }

    def change_password(
        self, access_token: str, previous_password: str, proposed_password: str
    ) -> None:
        for user_pool in self.user_pools.values():
            if access_token in user_pool.access_tokens:
                self._validate_password(
                    user_pool_id=user_pool.id, password=proposed_password
                )

                _, username = user_pool.access_tokens[access_token]
                user = self.admin_get_user(user_pool.id, username)

                if user.password != previous_password:
                    raise NotAuthorizedError(username)

                user.password = proposed_password
                if user.status in [
                    UserStatus.FORCE_CHANGE_PASSWORD,
                    UserStatus.RESET_REQUIRED,
                ]:
                    user.status = UserStatus.CONFIRMED

                break
        else:
            raise NotAuthorizedError(access_token)

    def admin_update_user_attributes(
        self, user_pool_id: str, username: str, attributes: List[Dict[str, str]]
    ) -> None:
        user = self.admin_get_user(user_pool_id, username)

        email = self._find_attr("email", attributes)
        self._verify_email_is_not_used(user_pool_id, email)

        user.update_attributes(attributes)

    def admin_delete_user_attributes(
        self, user_pool_id: str, username: str, attributes: List[str]
    ) -> None:
        self.admin_get_user(user_pool_id, username).delete_attributes(attributes)

    def admin_user_global_sign_out(self, user_pool_id: str, username: str) -> None:
        user_pool = self.describe_user_pool(user_pool_id)
        self.admin_get_user(user_pool_id, username)

        user_pool.sign_out(username)

    def global_sign_out(self, access_token: str) -> None:
        for user_pool in self.user_pools.values():
            if access_token in user_pool.access_tokens:
                _, username = user_pool.access_tokens[access_token]
                user_pool.sign_out(username)
                return

        raise NotAuthorizedError(access_token)

    def create_resource_server(
        self,
        user_pool_id: str,
        identifier: str,
        name: str,
        scopes: List[Dict[str, str]],
    ) -> CognitoResourceServer:
        user_pool = self.describe_user_pool(user_pool_id)

        if identifier in user_pool.resource_servers:
            raise InvalidParameterException(
                f"{identifier} already exists in user pool {user_pool_id}."
            )

        resource_server = CognitoResourceServer(user_pool_id, identifier, name, scopes)
        user_pool.resource_servers[identifier] = resource_server
        return resource_server

    def describe_resource_server(
        self, user_pool_id: str, identifier: str
    ) -> CognitoResourceServer:
        user_pool = self.user_pools.get(user_pool_id)
        if not user_pool:
            raise ResourceNotFoundError(f"User pool {user_pool_id} does not exist.")

        resource_server = user_pool.resource_servers.get(identifier)
        if not resource_server:
            raise ResourceNotFoundError(f"Resource server {identifier} does not exist.")

        return resource_server

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_resource_servers(self, user_pool_id: str) -> List[CognitoResourceServer]:
        user_pool = self.user_pools[user_pool_id]
        resource_servers = list(user_pool.resource_servers.values())
        return resource_servers

    def sign_up(
        self,
        client_id: str,
        username: str,
        password: str,
        attributes: List[Dict[str, str]],
    ) -> Tuple[CognitoIdpUser, Any]:
        user_pool = None
        for p in self.user_pools.values():
            if client_id in p.clients:
                user_pool = p
        if user_pool is None:
            raise ResourceNotFoundError(client_id)
        elif user_pool._get_user(username):
            raise UsernameExistsException(username)

        # UsernameAttributes are attributes (either `email` or `phone_number`
        # or both) than can be used in the place of a unique username. If the
        # user provides an email or phone number when signing up, the user pool
        # performs the following steps:
        # 1. populates the correct field (email, phone_number) with the value
        #    supplied for Username
        # 2. generates a persistent GUID for the user that will be returned as
        #    the value of `Username` in the `get-user` and `list-users`
        #    operations, as well as the value of `sub` in `IdToken` and
        #    `AccessToken`
        #
        # ref: https://docs.aws.amazon.com/cognito/latest/developerguide/user-pool-settings-attributes.html#user-pool-settings-aliases-settings
        has_username_attrs = user_pool.extended_config.get("UsernameAttributes")
        if has_username_attrs:
            username_attributes = user_pool.extended_config["UsernameAttributes"]
            # attribute_type should be one of `email`, `phone_number` or both
            for attribute_type in username_attributes:
                # check if provided username matches one of the attribute types in
                # `UsernameAttributes`
                if attribute_type in username_attributes and validate_username_format(
                    username, _format=attribute_type
                ):
                    # insert provided username into new user's attributes under the
                    # correct key
                    flattened_attrs = flatten_attrs(attributes or [])
                    flattened_attrs.update({attribute_type: username})
                    attributes = expand_attrs(flattened_attrs)

                    # once the username has been validated against a username attribute
                    # type, there is no need to attempt validation against the other
                    # type(s)
                    break

            else:
                # The provided username has not matched the required format for any
                # of the possible attributes
                raise InvalidParameterException(
                    "Username should be either an email or a phone number."
                )

        self._validate_password(user_pool.id, password)

        user = CognitoIdpUser(
            user_pool_id=user_pool.id,
            # set username to None so that it will be default to the internal GUID
            # when them user gets created
            username=None if has_username_attrs else username,
            password=password,
            attributes=attributes,
            status=UserStatus.UNCONFIRMED,
        )
        user_pool.users[user.username] = user

        has_email = "email" in user.attribute_lookup
        has_phone = "phone_number" in user.attribute_lookup
        verified_attributes = user_pool.extended_config.get(
            "AutoVerifiedAttributes", []
        )
        email_verified = (
            "email_verified" in user.attribute_lookup or "email" in verified_attributes
        )
        phone_verified = (
            "phone_number_verified" in user.attribute_lookup
            or "phone_number" in verified_attributes
        )
        if (has_email and email_verified) or (has_phone and phone_verified):
            recovery_settings = user_pool.extended_config["AccountRecoverySetting"]
            details = self._get_code_delivery_details(
                recovery_settings=recovery_settings, user=user, username=user.username
            )
            return user, details
        else:
            return user, None

    def confirm_sign_up(self, client_id: str, username: str) -> str:
        user_pool = None
        for p in self.user_pools.values():
            if client_id in p.clients:
                user_pool = p
        if user_pool is None:
            raise ResourceNotFoundError(client_id)

        user = self.admin_get_user(user_pool.id, username)

        user.status = UserStatus.CONFIRMED
        return ""

    def initiate_auth(
        self, client_id: str, auth_flow: str, auth_parameters: Dict[str, str]
    ) -> Dict[str, Any]:
        user_auth_flows = [
            AuthFlow.USER_SRP_AUTH,
            AuthFlow.REFRESH_TOKEN_AUTH,
            AuthFlow.REFRESH_TOKEN,
            AuthFlow.CUSTOM_AUTH,
            AuthFlow.USER_PASSWORD_AUTH,
        ]

        auth_flow = self._validate_auth_flow(
            auth_flow=auth_flow, valid_flows=user_auth_flows
        )

        user_pool: Optional[CognitoIdpUserPool] = None
        client: CognitoIdpUserPoolClient = None  # type: ignore[assignment]
        for p in self.user_pools.values():
            if client_id in p.clients:
                user_pool = p
                client = p.clients[client_id]
        if user_pool is None:
            raise ResourceNotFoundError(client_id)

        if auth_flow is AuthFlow.USER_SRP_AUTH:
            username: str = auth_parameters.get("USERNAME")  # type: ignore[assignment]
            srp_a = auth_parameters.get("SRP_A")
            if not srp_a:
                raise ResourceNotFoundError(srp_a)
            if client.generate_secret:
                secret_hash: str = auth_parameters.get("SECRET_HASH")  # type: ignore[assignment]
                if not check_secret_hash(
                    client.secret, client.id, username, secret_hash
                ):
                    raise NotAuthorizedError(secret_hash)

            user = self.admin_get_user(user_pool.id, username)

            if not user.enabled:
                raise NotAuthorizedError("User is disabled.")

            if user.status is UserStatus.UNCONFIRMED:
                raise UserNotConfirmedException("User is not confirmed.")

            session = str(random.uuid4())
            self.sessions[session] = (username, user_pool)

            return {
                "ChallengeName": "PASSWORD_VERIFIER",
                "Session": session,
                "ChallengeParameters": {
                    "SALT": random.uuid4().hex,
                    "SRP_B": random.uuid4().hex,
                    "USERNAME": user.username,
                    "USER_ID_FOR_SRP": user.id,
                    "SECRET_BLOCK": session,
                },
            }
        elif auth_flow is AuthFlow.USER_PASSWORD_AUTH:
            username: str = auth_parameters.get("USERNAME")  # type: ignore[no-redef]
            password: str = auth_parameters.get("PASSWORD")  # type: ignore[assignment]

            user = self.admin_get_user(user_pool.id, username)

            if not user:
                raise UserNotFoundError(username)

            if not user.enabled:
                raise NotAuthorizedError("User is disabled.")

            if user.password != password:
                raise NotAuthorizedError("Incorrect username or password.")

            if user.status is UserStatus.UNCONFIRMED:
                raise UserNotConfirmedException("User is not confirmed.")

            session = str(random.uuid4())
            self.sessions[session] = (username, user_pool)

            if user.status is UserStatus.FORCE_CHANGE_PASSWORD:
                return {
                    "ChallengeName": "NEW_PASSWORD_REQUIRED",
                    "ChallengeParameters": {"USERNAME": user.username},
                    "Session": session,
                }

            if (
                user.software_token_mfa_enabled
                and user.preferred_mfa_setting == "SOFTWARE_TOKEN_MFA"
            ) or (user_pool.mfa_config == "ON" and user.token_verified):
                return {
                    "ChallengeName": "SOFTWARE_TOKEN_MFA",
                    "ChallengeParameters": {},
                    "Session": session,
                }

            if user.sms_mfa_enabled and user.preferred_mfa_setting == "SMS_MFA":
                return {
                    "ChallengeName": "SMS_MFA",
                    "ChallengeParameters": {},
                    "Session": session,
                }

            new_refresh_token, origin_jti = user_pool.create_refresh_token(
                client_id, username
            )
            access_token, expires_in = user_pool.create_access_token(
                client_id, username, origin_jti=origin_jti
            )
            id_token, _ = user_pool.create_id_token(
                client_id, username, origin_jti=origin_jti
            )

            return {
                "AuthenticationResult": {
                    "IdToken": id_token,
                    "AccessToken": access_token,
                    "ExpiresIn": expires_in,
                    "RefreshToken": new_refresh_token,
                    "TokenType": "Bearer",
                }
            }
        elif auth_flow in (AuthFlow.REFRESH_TOKEN, AuthFlow.REFRESH_TOKEN_AUTH):
            refresh_token = auth_parameters.get("REFRESH_TOKEN")
            if not refresh_token:
                raise ResourceNotFoundError(refresh_token)

            res = user_pool.refresh_tokens.get(refresh_token)
            if res is None:
                raise NotAuthorizedError("Refresh Token has been revoked")

            client_id, username, _ = res
            if not username:
                raise ResourceNotFoundError(username)

            if client.generate_secret:
                secret_hash: str = auth_parameters.get("SECRET_HASH")  # type: ignore[no-redef]
                if not check_secret_hash(
                    client.secret, client.id, username, secret_hash
                ):
                    raise NotAuthorizedError(secret_hash)

            (
                access_token,
                id_token,
                expires_in,
            ) = user_pool.create_tokens_from_refresh_token(refresh_token)

            return {
                "AuthenticationResult": {
                    "IdToken": id_token,
                    "AccessToken": access_token,
                    "ExpiresIn": expires_in,
                    "TokenType": "Bearer",
                }
            }
        else:
            # We shouldn't get here due to enum validation of auth_flow
            return None  # type: ignore[return-value]

    def associate_software_token(
        self, access_token: str, session: str
    ) -> Dict[str, str]:
        secret_code = "asdfasdfasdf"
        if session:
            if session in self.sessions:
                return {"SecretCode": secret_code, "Session": session}
            raise NotAuthorizedError(session)

        for user_pool in self.user_pools.values():
            if access_token in user_pool.access_tokens:
                _, username = user_pool.access_tokens[access_token]
                self.admin_get_user(user_pool.id, username)

                return {"SecretCode": secret_code}

        raise NotAuthorizedError(access_token)

    def verify_software_token(self, access_token: str, session: str) -> Dict[str, str]:
        """
        The parameter UserCode has not yet been implemented
        """
        if session:
            if session not in self.sessions:
                raise ResourceNotFoundError(session)

            username, user_pool = self.sessions[session]
            user = self.admin_get_user(user_pool.id, username)
            user.token_verified = True

            session = str(random.uuid4())
            self.sessions[session] = (username, user_pool)
            return {"Status": "SUCCESS", "Session": session}

        for user_pool in self.user_pools.values():
            if access_token in user_pool.access_tokens:
                _, username = user_pool.access_tokens[access_token]
                user = self.admin_get_user(user_pool.id, username)

                user.token_verified = True

                session = str(random.uuid4())
                self.sessions[session] = (username, user_pool)
                return {"Status": "SUCCESS", "Session": session}

        raise NotAuthorizedError(access_token)

    def set_user_mfa_preference(
        self,
        access_token: str,
        software_token_mfa_settings: Dict[str, bool],
        sms_mfa_settings: Dict[str, bool],
    ) -> None:
        for user_pool in self.user_pools.values():
            if access_token in user_pool.access_tokens:
                _, username = user_pool.access_tokens[access_token]

                return self.admin_set_user_mfa_preference(
                    user_pool.id,
                    username,
                    software_token_mfa_settings,
                    sms_mfa_settings,
                )

        raise NotAuthorizedError(access_token)

    def admin_set_user_mfa_preference(
        self,
        user_pool_id: str,
        username: str,
        software_token_mfa_settings: Dict[str, bool],
        sms_mfa_settings: Dict[str, bool],
    ) -> None:
        user = self.admin_get_user(user_pool_id, username)

        if software_token_mfa_settings:
            if software_token_mfa_settings.get("Enabled"):
                if user.token_verified:
                    user.software_token_mfa_enabled = True
                else:
                    raise InvalidParameterException(
                        "User has not verified software token mfa"
                    )
            else:
                user.software_token_mfa_enabled = False

            if software_token_mfa_settings.get("PreferredMfa"):
                user.preferred_mfa_setting = "SOFTWARE_TOKEN_MFA"
            elif user.preferred_mfa_setting != "SMS_MFA":
                user.preferred_mfa_setting = ""

        if sms_mfa_settings:
            if sms_mfa_settings.get("Enabled"):
                user.sms_mfa_enabled = True
            else:
                user.sms_mfa_enabled = False

            if sms_mfa_settings.get("PreferredMfa"):
                user.preferred_mfa_setting = "SMS_MFA"
            elif user.preferred_mfa_setting != "SOFTWARE_TOKEN_MFA":
                user.preferred_mfa_setting = ""

        return None

    def _validate_password(self, user_pool_id: str, password: str) -> None:
        user_pool = self.describe_user_pool(user_pool_id)
        password_policy = user_pool.extended_config.get("Policies", {}).get(
            "PasswordPolicy", {}
        )
        minimum = password_policy.get("MinimumLength", 5)
        maximum = password_policy.get("MaximumLength", 99)
        require_uppercase = password_policy.get("RequireUppercase", True)
        require_lowercase = password_policy.get("RequireLowercase", True)
        require_numbers = password_policy.get("RequireNumbers", True)
        require_symbols = password_policy.get("RequireSymbols", True)

        flagl = minimum <= len(password) < maximum
        flagn = not require_numbers or bool(re.search(r"\d", password))
        # If we require symbols, we assume False - and check a symbol is present
        # If we don't require symbols, we assume True - and we could technically skip the for-loop
        flag_sc = not require_symbols
        sc = "^ $ * . [ ] { } ( ) ? \" ! @ # % & / \\ , > < ' : ; | _ ~ ` = + -"
        for i in password:
            if i in sc:
                flag_sc = True

        flag_u = not require_uppercase or bool(re.search(r"[A-Z]+", password))
        flag_lo = not require_lowercase or bool(re.search(r"[a-z]+", password))
        if not (flagl and flagn and flag_sc and flag_u and flag_lo):
            raise InvalidPasswordException()

    def admin_set_user_password(
        self, user_pool_id: str, username: str, password: str, permanent: bool
    ) -> None:
        user = self.admin_get_user(user_pool_id, username)
        # user.password = password
        self._validate_password(user_pool_id, password)
        user.password = password
        if permanent:
            user.status = UserStatus.CONFIRMED
        else:
            user.status = UserStatus.FORCE_CHANGE_PASSWORD

    def add_custom_attributes(
        self, user_pool_id: str, custom_attributes: List[Dict[str, Any]]
    ) -> None:
        user_pool = self.describe_user_pool(user_pool_id)
        user_pool.add_custom_attributes(custom_attributes)

    def update_user_attributes(
        self, access_token: str, attributes: List[Dict[str, str]]
    ) -> None:
        """
        The parameter ClientMetadata has not yet been implemented. No CodeDeliveryDetails are returned.
        """
        for user_pool in self.user_pools.values():
            if access_token in user_pool.access_tokens:
                _, username = user_pool.access_tokens[access_token]
                user = self.admin_get_user(user_pool.id, username)

                email = self._find_attr("email", attributes)
                self._verify_email_is_not_used(user_pool.id, email)

                user.update_attributes(attributes)
                return

        raise NotAuthorizedError(access_token)

    def _find_attr(self, name: str, attrs: List[Dict[str, str]]) -> Optional[str]:
        return next((a["Value"] for a in attrs if a["Name"] == name), None)

    def _verify_email_is_not_used(
        self, user_pool_id: str, email: Optional[str]
    ) -> None:
        if not email:
            # We're not updating emails
            return
        user_pool = self.describe_user_pool(user_pool_id)
        if "email" not in user_pool.extended_config.get("UsernameAttributes", []):
            # email is not used as a username - duplicate emails are allowed
            return

        for user in user_pool.users.values():
            if user.attribute_lookup.get("email", "") == email:
                raise AliasExistsException


class RegionAgnosticBackend:
    # Some operations are unauthenticated
    # Without authentication-header, we lose the context of which region the request was send to
    # This backend will cycle through all backends as a workaround

    def __init__(self, account_id: str, region_name: str):
        self.account_id = account_id
        self.region_name = region_name

    def _find_backend_by_access_token(self, access_token: str) -> CognitoIdpBackend:
        for account_specific_backends in cognitoidp_backends.values():
            for region, backend in account_specific_backends.items():
                if region == "global":
                    continue
                for p in backend.user_pools.values():
                    if access_token in p.access_tokens:
                        return backend
        return cognitoidp_backends[self.account_id][self.region_name]

    def _find_backend_by_access_token_or_session(
        self, access_token: str, session: str
    ) -> CognitoIdpBackend:
        for account_specific_backends in cognitoidp_backends.values():
            for region, backend in account_specific_backends.items():
                if region == "global":
                    continue
                if session and session in backend.sessions:
                    return backend
                for p in backend.user_pools.values():
                    if access_token and access_token in p.access_tokens:
                        return backend
        return cognitoidp_backends[self.account_id][self.region_name]

    def _find_backend_for_clientid(self, client_id: str) -> CognitoIdpBackend:
        for account_specific_backends in cognitoidp_backends.values():
            for region, backend in account_specific_backends.items():
                if region == "global":
                    continue
                for p in backend.user_pools.values():
                    if client_id in p.clients:
                        return backend
        return cognitoidp_backends[self.account_id][self.region_name]

    def sign_up(
        self,
        client_id: str,
        username: str,
        password: str,
        attributes: List[Dict[str, str]],
    ) -> Tuple[CognitoIdpUser, Any]:
        backend = self._find_backend_for_clientid(client_id)
        return backend.sign_up(client_id, username, password, attributes)

    def initiate_auth(
        self, client_id: str, auth_flow: str, auth_parameters: Dict[str, str]
    ) -> Dict[str, Any]:
        backend = self._find_backend_for_clientid(client_id)
        return backend.initiate_auth(client_id, auth_flow, auth_parameters)

    def confirm_sign_up(self, client_id: str, username: str) -> str:
        backend = self._find_backend_for_clientid(client_id)
        return backend.confirm_sign_up(client_id, username)

    def get_user(self, access_token: str) -> CognitoIdpUser:
        backend = self._find_backend_by_access_token(access_token)
        return backend.get_user(access_token)

    def admin_respond_to_auth_challenge(
        self,
        session: str,
        client_id: str,
        challenge_name: str,
        challenge_responses: Dict[str, str],
    ) -> Dict[str, Any]:
        backend = self._find_backend_for_clientid(client_id)
        return backend.admin_respond_to_auth_challenge(
            session, client_id, challenge_name, challenge_responses
        )

    def respond_to_auth_challenge(
        self,
        session: str,
        client_id: str,
        challenge_name: str,
        challenge_responses: Dict[str, str],
    ) -> Dict[str, Any]:
        backend = self._find_backend_for_clientid(client_id)
        return backend.respond_to_auth_challenge(
            session, client_id, challenge_name, challenge_responses
        )

    def associate_software_token(
        self, access_token: str, session: str
    ) -> Dict[str, str]:
        backend = self._find_backend_by_access_token_or_session(access_token, session)
        return backend.associate_software_token(access_token, session)

    def verify_software_token(self, access_token: str, session: str) -> Dict[str, str]:
        backend = self._find_backend_by_access_token_or_session(access_token, session)
        return backend.verify_software_token(access_token, session)

    def set_user_mfa_preference(
        self,
        access_token: str,
        software_token_mfa_settings: Dict[str, bool],
        sms_mfa_settings: Dict[str, bool],
    ) -> None:
        backend = self._find_backend_by_access_token(access_token)
        return backend.set_user_mfa_preference(
            access_token, software_token_mfa_settings, sms_mfa_settings
        )

    def update_user_attributes(
        self, access_token: str, attributes: List[Dict[str, str]]
    ) -> None:
        backend = self._find_backend_by_access_token(access_token)
        return backend.update_user_attributes(access_token, attributes)


cognitoidp_backends = BackendDict(CognitoIdpBackend, "cognito-idp")


# Hack to help moto-server process requests on localhost, where the region isn't
# specified in the host header. Some endpoints (change password, confirm forgot
# password) have no authorization header from which to extract the region.
def find_account_region_by_value(
    key: str, value: str, fallback: Tuple[str, str]
) -> Tuple[str, str]:
    for account_id, account_specific_backend in cognitoidp_backends.items():
        for region, backend in account_specific_backend.items():
            for user_pool in backend.user_pools.values():
                if key == "client_id" and value in user_pool.clients:
                    return account_id, region

                if key == "access_token" and value in user_pool.access_tokens:
                    return account_id, region
    # If we can't find the `client_id` or `access_token`, we just pass
    # back a default backend region, which will raise the appropriate
    # error message (e.g. NotAuthorized or NotFound).
    return fallback
