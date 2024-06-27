from typing import TYPE_CHECKING, Any, Dict, List, NamedTuple, Optional, Tuple

from botocore.exceptions import ParamValidationError

from moto.core.base_backend import BaseBackend
from moto.core.models import BackendDict
from moto.moto_api._internal import mock_random
from moto.utilities.paginator import paginate

from .exceptions import (
    ConflictException,
    ResourceNotFoundException,
    ValidationException,
)

if TYPE_CHECKING:
    from typing_extensions import Self


class Group(NamedTuple):
    GroupId: str
    DisplayName: str
    ExternalIds: List[Optional[Dict[str, str]]]
    Description: str
    IdentityStoreId: str


class Name(NamedTuple):
    Formatted: Optional[str]
    FamilyName: Optional[str]
    GivenName: Optional[str]
    MiddleName: Optional[str]
    HonorificPrefix: Optional[str]
    HonorificSuffix: Optional[str]

    @classmethod
    def from_dict(cls, name_dict: Dict[str, str]) -> "Optional[Self]":
        if not name_dict:
            return None
        return cls(
            name_dict.get("Formatted"),
            name_dict.get("FamilyName"),
            name_dict.get("GivenName"),
            name_dict.get("MiddleName"),
            name_dict.get("HonorificPrefix"),
            name_dict.get("HonorificSuffix"),
        )


class User(NamedTuple):
    UserId: str
    IdentityStoreId: str
    UserName: str
    Name: Optional[Name]
    DisplayName: str
    NickName: str
    ProfileUrl: str
    Emails: List[Dict[str, str]]
    Addresses: List[Dict[str, str]]
    PhoneNumbers: List[Dict[str, str]]
    UserType: str
    Title: str
    PreferredLanguage: str
    Locale: str
    Timezone: str


class IdentityStoreData:
    def __init__(self) -> None:
        self.groups: Dict[str, Group] = {}
        self.users: Dict[str, User] = {}
        self.group_memberships: Dict[str, Any] = {}


class IdentityStoreBackend(BaseBackend):
    """Implementation of IdentityStore APIs."""

    PAGINATION_MODEL = {
        "list_group_memberships": {
            "input_token": "next_token",
            "limit_key": "max_results",
            "limit_default": 100,
            "unique_attribute": "MembershipId",
        },
        "list_group_memberships_for_member": {
            "input_token": "next_token",
            "limit_key": "max_results",
            "limit_default": 100,
            "unique_attribute": "MembershipId",
        },
        "list_groups": {
            "input_token": "next_token",
            "limit_key": "max_results",
            "limit_default": 100,
            "unique_attribute": "GroupId",
        },
        "list_users": {
            "input_token": "next_token",
            "limit_key": "max_results",
            "limit_default": 100,
            "unique_attribute": "UserId",
        },
    }

    def __init__(self, region_name: str, account_id: str) -> None:
        super().__init__(region_name, account_id)
        self.identity_stores: Dict[str, IdentityStoreData] = {}

    def create_group(
        self, identity_store_id: str, display_name: str, description: str
    ) -> Tuple[str, str]:
        identity_store = self.__get_identity_store(identity_store_id)

        matching = [
            g for g in identity_store.groups.values() if g.DisplayName == display_name
        ]
        if len(matching) > 0:
            raise ConflictException(
                message="Duplicate GroupDisplayName",
                reason="UNIQUENESS_CONSTRAINT_VIOLATION",
            )

        group_id = str(mock_random.uuid4())
        group = Group(
            group_id,
            display_name,
            [],
            description,
            identity_store_id,
        )
        identity_store.groups[group_id] = group
        return group_id, identity_store_id

    def get_group_id(
        self, identity_store_id: str, alternate_identifier: Dict[str, Any]
    ) -> Tuple[str, str]:
        """
        The ExternalId alternate identifier is not yet implemented
        """
        identity_store = self.__get_identity_store(identity_store_id)
        if "UniqueAttribute" in alternate_identifier:
            if (
                "AttributeValue" in alternate_identifier["UniqueAttribute"]
                and alternate_identifier["UniqueAttribute"]["AttributePath"].lower()
                == "displayname"
            ):
                for g in identity_store.groups.values():
                    if (
                        g.DisplayName
                        == alternate_identifier["UniqueAttribute"]["AttributeValue"]
                    ):
                        return g.GroupId, identity_store_id

        raise ResourceNotFoundException(
            message="GROUP not found.", resource_type="GROUP"
        )

    def describe_group(self, identity_store_id: str, group_id: str) -> Group:
        identity_store = self.__get_identity_store(identity_store_id)
        if group_id in identity_store.groups:
            g = identity_store.groups[group_id]
            # External Ids are not implemented
            external_ids: List[Any] = []
            return Group(
                g.GroupId,
                g.DisplayName,
                external_ids,
                g.Description,
                identity_store_id,
            )
        raise ResourceNotFoundException(
            message="GROUP not found.", resource_type="GROUP"
        )

    def delete_group(self, identity_store_id: str, group_id: str) -> None:
        identity_store = self.__get_identity_store(identity_store_id)
        if group_id in identity_store.groups:
            del identity_store.groups[group_id]

    def create_user(
        self,
        identity_store_id: str,
        user_name: str,
        name: Dict[str, str],
        display_name: str,
        nick_name: str,
        profile_url: str,
        emails: List[Dict[str, Any]],
        addresses: List[Dict[str, Any]],
        phone_numbers: List[Dict[str, Any]],
        user_type: str,
        title: str,
        preferred_language: str,
        locale: str,
        timezone: str,
    ) -> Tuple[str, str]:
        identity_store = self.__get_identity_store(identity_store_id)
        user_id = str(mock_random.uuid4())

        new_user = User(
            user_id,
            identity_store_id,
            user_name,
            Name.from_dict(name),
            display_name,
            nick_name,
            profile_url,
            emails,
            addresses,
            phone_numbers,
            user_type,
            title,
            preferred_language,
            locale,
            timezone,
        )
        self.__validate_create_user(new_user, identity_store)

        identity_store.users[user_id] = new_user

        return user_id, identity_store_id

    def describe_user(self, identity_store_id: str, user_id: str) -> User:
        identity_store = self.__get_identity_store(identity_store_id)

        if user_id in identity_store.users:
            return identity_store.users[user_id]

        raise ResourceNotFoundException(message="USER not found.", resource_type="USER")

    def delete_user(self, identity_store_id: str, user_id: str) -> None:
        identity_store = self.__get_identity_store(identity_store_id)

        if user_id in identity_store.users:
            del identity_store.users[user_id]

    def create_group_membership(
        self, identity_store_id: str, group_id: str, member_id: Dict[str, str]
    ) -> Tuple[str, str]:
        identity_store = self.__get_identity_store(identity_store_id)
        user_id = member_id["UserId"]
        if user_id not in identity_store.users:
            raise ResourceNotFoundException(
                message="Member does not exist", resource_type="USER"
            )

        if group_id not in identity_store.groups:
            raise ResourceNotFoundException(
                message="Group does not exist", resource_type="GROUP"
            )

        membership_id = str(mock_random.uuid4())
        identity_store.group_memberships[membership_id] = {
            "IdentityStoreId": identity_store_id,
            "MembershipId": membership_id,
            "GroupId": group_id,
            "MemberId": {"UserId": user_id},
        }

        return membership_id, identity_store_id

    @paginate(pagination_model=PAGINATION_MODEL)  # type: ignore
    def list_group_memberships(
        self, identity_store_id: str, group_id: str
    ) -> List[Any]:
        identity_store = self.__get_identity_store(identity_store_id)

        return [
            m
            for m in identity_store.group_memberships.values()
            if m["GroupId"] == group_id
        ]

    @paginate(pagination_model=PAGINATION_MODEL)  # type: ignore
    def list_group_memberships_for_member(
        self, identity_store_id: str, member_id: Dict[str, str]
    ) -> List[Any]:
        identity_store = self.__get_identity_store(identity_store_id)
        user_id = member_id["UserId"]

        return [
            m
            for m in identity_store.group_memberships.values()
            if m["MemberId"]["UserId"] == user_id
        ]

    @paginate(pagination_model=PAGINATION_MODEL)  # type: ignore
    def list_groups(
        self, identity_store_id: str, filters: List[Dict[str, str]]
    ) -> List[Dict[str, Any]]:
        identity_store = self.__get_identity_store(identity_store_id)

        if filters:
            if filters[0].get("AttributePath") == "DisplayName":
                displayname = filters[0].get("AttributeValue")
                return [
                    m._asdict()
                    for m in identity_store.groups.values()
                    if m.DisplayName == displayname
                ]

        return [m._asdict() for m in identity_store.groups.values()]

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_users(
        self, identity_store_id: str, filters: List[Dict[str, str]]
    ) -> List[Dict[str, str]]:
        identity_store = self.__get_identity_store(identity_store_id)

        users = []
        if filters:
            if filters[0].get("AttributePath") == "UserName":
                username = filters[0].get("AttributeValue")

                for m in identity_store.users.values():
                    if m.UserName == username:
                        user = m._asdict()
                        if user.get("Name"):
                            user["Name"] = m.Name._asdict()  # type: ignore
                        users.append(user)
                return users

        for m in identity_store.users.values():
            user = m._asdict()
            if user.get("Name"):
                user["Name"] = m.Name._asdict()  # type: ignore
            users.append(user)
        return users

    def delete_group_membership(
        self, identity_store_id: str, membership_id: str
    ) -> None:
        identity_store = self.__get_identity_store(identity_store_id)
        if membership_id in identity_store.group_memberships:
            del identity_store.group_memberships[membership_id]

    def __get_identity_store(self, store_id: str) -> IdentityStoreData:
        if len(store_id) < 1:
            raise ParamValidationError(
                msg="Invalid length for parameter IdentityStoreId, value: 0, valid min length: 1"
            )
        if store_id not in self.identity_stores:
            self.identity_stores[store_id] = IdentityStoreData()
        return self.identity_stores[store_id]

    def __validate_create_user(
        self, new_user: User, identity_store: IdentityStoreData
    ) -> None:
        if not new_user.UserName:
            raise ValidationException(message="userName is a required attribute")

        missing = []

        if not new_user.DisplayName:
            missing.append("displayname")
        if not new_user.Name:
            missing.append("name")
        else:
            if not new_user.Name.GivenName:
                missing.append("givenname")
            if not new_user.Name.FamilyName:
                missing.append("familyname")

        if len(missing) > 0:
            message = ", ".join(
                [f"{att}: The attribute {att} is required" for att in missing]
            )
            raise ValidationException(message=message)

        matching = [
            u for u in identity_store.users.values() if u.UserName == new_user.UserName
        ]
        if len(matching) > 0:
            raise ConflictException(
                message="Duplicate UserName", reason="UNIQUENESS_CONSTRAINT_VIOLATION"
            )


identitystore_backends = BackendDict(IdentityStoreBackend, "identitystore")
