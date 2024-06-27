"""Handles incoming identitystore requests, invokes methods, returns responses."""

import json
from typing import Any, Dict, NamedTuple, Optional

from moto.core.responses import BaseResponse

from .models import IdentityStoreBackend, identitystore_backends


class IdentityStoreResponse(BaseResponse):
    """Handler for IdentityStore requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="identitystore")

    @property
    def identitystore_backend(self) -> IdentityStoreBackend:
        """Return backend instance specific for this region."""
        return identitystore_backends[self.current_account][self.region]

    def create_group(self) -> str:
        identity_store_id = self._get_param("IdentityStoreId")
        display_name = self._get_param("DisplayName")
        description = self._get_param("Description")
        group_id, identity_store_id = self.identitystore_backend.create_group(
            identity_store_id=identity_store_id,
            display_name=display_name,
            description=description,
        )
        return json.dumps(dict(GroupId=group_id, IdentityStoreId=identity_store_id))

    def create_group_membership(self) -> str:
        identity_store_id = self._get_param("IdentityStoreId")
        group_id = self._get_param("GroupId")
        member_id = self._get_param("MemberId")
        (
            membership_id,
            identity_store_id,
        ) = self.identitystore_backend.create_group_membership(
            identity_store_id=identity_store_id,
            group_id=group_id,
            member_id=member_id,
        )

        return json.dumps(
            dict(MembershipId=membership_id, IdentityStoreId=identity_store_id)
        )

    def create_user(self) -> str:
        user_id, identity_store_id = self.identitystore_backend.create_user(
            self._get_param("IdentityStoreId"),
            self._get_param("UserName"),
            self._get_param("Name"),
            self._get_param("DisplayName"),
            self._get_param("NickName"),
            self._get_param("ProfileUrl"),
            self._get_param("Emails"),
            self._get_param("Addresses"),
            self._get_param("PhoneNumbers"),
            self._get_param("UserType"),
            self._get_param("Title"),
            self._get_param("PreferredLanguage"),
            self._get_param("Locale"),
            self._get_param("Timezone"),
        )
        return json.dumps(dict(UserId=user_id, IdentityStoreId=identity_store_id))

    def get_group_id(self) -> str:
        identity_store_id = self._get_param("IdentityStoreId")
        alternate_identifier = self._get_param("AlternateIdentifier")
        group_id, identity_store_id = self.identitystore_backend.get_group_id(
            identity_store_id=identity_store_id,
            alternate_identifier=alternate_identifier,
        )
        return json.dumps(dict(GroupId=group_id, IdentityStoreId=identity_store_id))

    def describe_user(self) -> str:
        identity_store_id = self._get_param("IdentityStoreId")
        user_id = self._get_param("UserId")
        (
            user_id,
            identity_store_id,
            user_name,
            name,
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
        ) = self.identitystore_backend.describe_user(
            identity_store_id=identity_store_id,
            user_id=user_id,
        )
        return json.dumps(
            dict(
                UserName=user_name,
                UserId=user_id,
                ExternalIds=None,
                Name=self.named_tuple_to_dict(name),
                DisplayName=display_name,
                NickName=nick_name,
                ProfileUrl=profile_url,
                Emails=emails,
                Addresses=addresses,
                PhoneNumbers=phone_numbers,
                UserType=user_type,
                Title=title,
                PreferredLanguage=preferred_language,
                Locale=locale,
                Timezone=timezone,
                IdentityStoreId=identity_store_id,
            )
        )

    def list_group_memberships(self) -> str:
        identity_store_id = self._get_param("IdentityStoreId")
        group_id = self._get_param("GroupId")
        max_results = self._get_param("MaxResults")
        next_token = self._get_param("NextToken")
        (
            group_memberships,
            next_token,
        ) = self.identitystore_backend.list_group_memberships(
            identity_store_id=identity_store_id,
            group_id=group_id,
            max_results=max_results,
            next_token=next_token,
        )

        return json.dumps(
            dict(GroupMemberships=group_memberships, NextToken=next_token)
        )

    def list_group_memberships_for_member(self) -> str:
        identity_store_id = self._get_param("IdentityStoreId")
        member_id = self._get_param("MemberId")
        max_results = self._get_param("MaxResults")
        next_token = self._get_param("NextToken")
        (
            group_memberships,
            next_token,
        ) = self.identitystore_backend.list_group_memberships_for_member(
            identity_store_id=identity_store_id,
            member_id=member_id,
            max_results=max_results,
            next_token=next_token,
        )

        return json.dumps(
            dict(GroupMemberships=group_memberships, NextToken=next_token)
        )

    def list_groups(self) -> str:
        identity_store_id = self._get_param("IdentityStoreId")
        max_results = self._get_param("MaxResults")
        next_token = self._get_param("NextToken")
        filters = self._get_param("Filters")
        (
            groups,
            next_token,
        ) = self.identitystore_backend.list_groups(
            identity_store_id=identity_store_id,
            max_results=max_results,
            next_token=next_token,
            filters=filters,
        )

        return json.dumps(dict(Groups=groups, NextToken=next_token))

    def describe_group(self) -> str:
        identity_store_id = self._get_param("IdentityStoreId")
        group_id = self._get_param("GroupId")
        (
            group_id,
            display_name,
            external_ids,
            description,
            identity_store_id,
        ) = self.identitystore_backend.describe_group(
            identity_store_id=identity_store_id,
            group_id=group_id,
        )
        return json.dumps(
            dict(
                GroupId=group_id,
                DisplayName=display_name,
                ExternalIds=external_ids,
                Description=description,
                IdentityStoreId=identity_store_id,
            )
        )

    def list_users(self) -> str:
        identity_store_id = self._get_param("IdentityStoreId")
        max_results = self._get_param("MaxResults")
        next_token = self._get_param("NextToken")
        filters = self._get_param("Filters")
        (
            users,
            next_token,
        ) = self.identitystore_backend.list_users(
            identity_store_id=identity_store_id,
            max_results=max_results,
            next_token=next_token,
            filters=filters,
        )

        return json.dumps(dict(Users=users, NextToken=next_token))

    def delete_group(self) -> str:
        identity_store_id = self._get_param("IdentityStoreId")
        group_id = self._get_param("GroupId")
        self.identitystore_backend.delete_group(
            identity_store_id=identity_store_id,
            group_id=group_id,
        )
        return json.dumps(dict())

    def delete_group_membership(self) -> str:
        identity_store_id = self._get_param("IdentityStoreId")
        membership_id = self._get_param("MembershipId")
        self.identitystore_backend.delete_group_membership(
            identity_store_id=identity_store_id,
            membership_id=membership_id,
        )
        return json.dumps(dict())

    def delete_user(self) -> str:
        identity_store_id = self._get_param("IdentityStoreId")
        user_id = self._get_param("UserId")
        self.identitystore_backend.delete_user(
            identity_store_id=identity_store_id,
            user_id=user_id,
        )
        return json.dumps(dict())

    def named_tuple_to_dict(
        self, value: Optional[NamedTuple]
    ) -> Optional[Dict[str, Any]]:
        if value:
            return value._asdict()
        return None
