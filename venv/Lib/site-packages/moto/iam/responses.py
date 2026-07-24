from __future__ import annotations

from moto.core.responses import ActionResult, BaseResponse, EmptyResult
from moto.core.serialize import never_return, return_if_not_empty, url_encode

from .exceptions import NotFoundException
from .models import IAMBackend, User, iam_backends
from .utils import is_role_resource


class IamResponse(BaseResponse):
    # IAM is inconsistent in which attributes it returns for different operations,
    # and some attributes require specialized processing, so we need to explicitly
    # transform some individual response elements using specialized functions.
    #
    # IAM resource-listing operations, for example, return a subset of the available
    # attributes for the resource.
    # See: https://docs.aws.amazon.com/IAM/latest/APIReference/API_ListRoles.html
    RESPONSE_KEY_PATH_TO_TRANSFORMER = {
        "CreatePolicyResponse.Policy.Description": never_return,
        "CreatePolicyResponse.Policy.Tags": return_if_not_empty,
        "CreateRoleResponse.Role.AssumeRolePolicyDocument": url_encode,
        "CreateRoleResponse.Role.PermissionsBoundary": return_if_not_empty,
        "CreateRoleResponse.Role.Tags": return_if_not_empty,
        "GetGroupPolicyResponse.PolicyDocument": url_encode,
        "GetPolicyResponse.Policy.Description": return_if_not_empty,
        "GetRolePolicyResponse.PolicyDocument": url_encode,
        "GetRoleResponse.Role.AssumeRolePolicyDocument": url_encode,
        "GetRoleResponse.Role.PermissionsBoundary": return_if_not_empty,
        "GetRoleResponse.Role.Tags": return_if_not_empty,
        "GetUserPolicyResponse.PolicyDocument": url_encode,
        "ListPoliciesResponse.Policies.Policy.Description": never_return,
        "ListPoliciesResponse.Policies.Policy.Tags": never_return,
        "ListRolesResponse.Roles.Role.PermissionsBoundary": never_return,
        "ListRolesResponse.Roles.Role.RoleLastUsed": never_return,
        "ListRolesResponse.Roles.Role.Tags": never_return,
        "ListUsersResponse.Users.User.Tags": never_return,
    }

    def __init__(self) -> None:
        super().__init__(service_name="iam")
        self.automated_parameter_parsing = True

    @property
    def backend(self) -> IAMBackend:
        return iam_backends[self.current_account][self.partition]

    def _determine_resource(self) -> str:
        if is_role_resource(self.data):
            role_name = self.data["RoleName"][0]
            return self._resolve_resource_arn(role_name)

        return "*"

    def _resolve_resource_arn(self, role_name: str) -> str:
        try:
            self.backend.get_role(role_name)
            role_object = self.backend.get_role(role_name)

            return role_object.arn
        except NotFoundException:
            return "*"

    def attach_role_policy(self) -> ActionResult:
        policy_arn = self._get_param("PolicyArn")
        role_name = self._get_param("RoleName")
        self.backend.attach_role_policy(policy_arn, role_name)
        return EmptyResult()

    def detach_role_policy(self) -> ActionResult:
        role_name = self._get_param("RoleName")
        policy_arn = self._get_param("PolicyArn")
        self.backend.detach_role_policy(policy_arn, role_name)
        return EmptyResult()

    def attach_group_policy(self) -> ActionResult:
        policy_arn = self._get_param("PolicyArn")
        group_name = self._get_param("GroupName")
        self.backend.attach_group_policy(policy_arn, group_name)
        return EmptyResult()

    def detach_group_policy(self) -> ActionResult:
        policy_arn = self._get_param("PolicyArn")
        group_name = self._get_param("GroupName")
        self.backend.detach_group_policy(policy_arn, group_name)
        return EmptyResult()

    def attach_user_policy(self) -> ActionResult:
        policy_arn = self._get_param("PolicyArn")
        user_name = self._get_param("UserName")
        self.backend.attach_user_policy(policy_arn, user_name)
        return EmptyResult()

    def detach_user_policy(self) -> ActionResult:
        policy_arn = self._get_param("PolicyArn")
        user_name = self._get_param("UserName")
        self.backend.detach_user_policy(policy_arn, user_name)
        return EmptyResult()

    def create_policy(self) -> ActionResult:
        description = self._get_param("Description")
        path = self._get_param("Path")
        policy_document = self._get_param("PolicyDocument")
        policy_name = self._get_param("PolicyName")
        tags = self._get_param("Tags", [])
        policy = self.backend.create_policy(
            description, path, policy_document, policy_name, tags
        )
        result = {"Policy": policy}
        return ActionResult(result)

    def get_policy(self) -> ActionResult:
        policy_arn = self._get_param("PolicyArn")
        policy = self.backend.get_policy(policy_arn)
        result = {"Policy": policy}
        return ActionResult(result)

    def list_attached_role_policies(self) -> ActionResult:
        marker = self._get_param("Marker")
        max_items = self._get_int_param("MaxItems", 100)
        path_prefix = self._get_param("PathPrefix", "/")
        role_name = self._get_param("RoleName")
        policies, marker = self.backend.list_attached_role_policies(
            role_name, marker=marker, max_items=max_items, path_prefix=path_prefix
        )
        result = {
            "AttachedPolicies": [
                {"PolicyName": p.name, "PolicyArn": p.arn} for p in policies
            ],
            "Marker": marker,
            "IsTruncated": bool(marker),
        }
        return ActionResult(result)

    def list_attached_group_policies(self) -> ActionResult:
        marker = self._get_param("Marker")
        max_items = self._get_int_param("MaxItems", 100)
        path_prefix = self._get_param("PathPrefix", "/")
        group_name = self._get_param("GroupName")
        policies, marker = self.backend.list_attached_group_policies(
            group_name, marker=marker, max_items=max_items, path_prefix=path_prefix
        )
        result = {
            "AttachedPolicies": [
                {"PolicyName": p.name, "PolicyArn": p.arn} for p in policies
            ],
            "Marker": marker,
            "IsTruncated": bool(marker),
        }
        return ActionResult(result)

    def list_attached_user_policies(self) -> ActionResult:
        marker = self._get_param("Marker")
        max_items = self._get_int_param("MaxItems", 100)
        path_prefix = self._get_param("PathPrefix", "/")
        user_name = self._get_param("UserName")
        policies, marker = self.backend.list_attached_user_policies(
            user_name, marker=marker, max_items=max_items, path_prefix=path_prefix
        )
        result = {
            "AttachedPolicies": [
                {"PolicyName": p.name, "PolicyArn": p.arn} for p in policies
            ],
            "Marker": marker,
            "IsTruncated": bool(marker),
        }
        return ActionResult(result)

    def list_policies(self) -> ActionResult:
        marker = self._get_param("Marker")
        max_items = self._get_int_param("MaxItems", 100)
        only_attached = self._get_bool_param("OnlyAttached", False)
        path_prefix = self._get_param("PathPrefix", "/")
        scope = self._get_param("Scope", "All")
        policies, marker = self.backend.list_policies(
            marker, max_items, only_attached, path_prefix, scope
        )
        result = {"Policies": policies, "Marker": marker, "IsTruncated": bool(marker)}
        return ActionResult(result)

    def list_entities_for_policy(self) -> ActionResult:
        policy_arn = self._get_param("PolicyArn")

        # Options 'User'|'Role'|'Group'|'LocalManagedPolicy'|'AWSManagedPolicy
        entity = self._get_param("EntityFilter")
        path_prefix = self._get_param("PathPrefix")
        # policy_usage_filter = self._get_param('PolicyUsageFilter')
        marker = self._get_param("Marker")
        max_items = self._get_param("MaxItems")

        entity_roles = []
        entity_groups = []
        entity_users = []

        if not entity or entity == "User":
            users = self.backend.list_users(path_prefix, marker, max_items)
            if users:
                for user in users:
                    for p in user.managed_policies:
                        if p == policy_arn:
                            entity_users.append(user)

        if not entity or entity == "Role":
            roles, _ = self.backend.list_roles(path_prefix, marker, max_items)
            if roles:
                for role in roles:
                    for p in role.managed_policies:
                        if p == policy_arn:
                            entity_roles.append(role)

        if not entity or entity == "Group":
            groups = self.backend.list_groups()
            if groups:
                for group in groups:
                    for p in group.managed_policies:
                        if p == policy_arn:
                            entity_groups.append(group)

        if entity == "LocalManagedPolicy" or entity == "AWSManagedPolicy":
            users = self.backend.list_users(path_prefix, marker, max_items)
            if users:
                for user in users:
                    for p in user.managed_policies:
                        if p == policy_arn:
                            entity_users.append(user)

            roles, _ = self.backend.list_roles(path_prefix, marker, max_items)
            if roles:
                for role in roles:
                    for p in role.managed_policies:
                        if p == policy_arn:
                            entity_roles.append(role)

            groups = self.backend.list_groups()
            if groups:
                for group in groups:
                    for p in group.managed_policies:
                        if p == policy_arn:
                            entity_groups.append(group)

        result = {
            "PolicyUsers": entity_users,
            "PolicyRoles": entity_roles,
            "PolicyGroups": entity_groups,
            "IsTruncated": False,
        }
        return ActionResult(result)

    def set_default_policy_version(self) -> ActionResult:
        policy_arn = self._get_param("PolicyArn")
        version_id = self._get_param("VersionId")
        self.backend.set_default_policy_version(policy_arn, version_id)
        return EmptyResult()

    def create_role(self) -> ActionResult:
        role_name = self._get_param("RoleName")
        path = self._get_param("Path")
        assume_role_policy_document = self._get_param("AssumeRolePolicyDocument")
        permissions_boundary = self._get_param("PermissionsBoundary")
        description = self._get_param("Description")
        tags = self._get_param("Tags", [])
        max_session_duration = self._get_param("MaxSessionDuration", 3600)

        role = self.backend.create_role(
            role_name,
            assume_role_policy_document,
            path,
            permissions_boundary,
            description,
            tags,
            max_session_duration,
        )
        result = {"Role": role}
        return ActionResult(result)

    def get_role(self) -> ActionResult:
        role_name = self._get_param("RoleName")
        role = self.backend.get_role(role_name)
        result = {"Role": role}
        return ActionResult(result)

    def delete_role(self) -> ActionResult:
        role_name = self._get_param("RoleName")
        self.backend.delete_role(role_name)
        return EmptyResult()

    def list_role_policies(self) -> ActionResult:
        role_name = self._get_param("RoleName")
        role_policies_names = self.backend.list_role_policies(role_name)
        result = {"PolicyNames": role_policies_names, "IsTruncated": False}
        return ActionResult(result)

    def put_role_policy(self) -> ActionResult:
        role_name = self._get_param("RoleName")
        policy_name = self._get_param("PolicyName")
        policy_document = self._get_param("PolicyDocument")
        self.backend.put_role_policy(role_name, policy_name, policy_document)
        return EmptyResult()

    def delete_role_policy(self) -> ActionResult:
        role_name = self._get_param("RoleName")
        policy_name = self._get_param("PolicyName")
        self.backend.delete_role_policy(role_name, policy_name)
        return EmptyResult()

    def get_role_policy(self) -> ActionResult:
        role_name = self._get_param("RoleName")
        policy_name = self._get_param("PolicyName")
        policy_name, policy_document = self.backend.get_role_policy(
            role_name, policy_name
        )
        result = {
            "RoleName": role_name,
            "PolicyName": policy_name,
            "PolicyDocument": policy_document,
        }
        return ActionResult(result)

    def update_assume_role_policy(self) -> ActionResult:
        role_name = self._get_param("RoleName")
        policy_document = self._get_param("PolicyDocument")
        self.backend.update_assume_role_policy(role_name, policy_document)
        return EmptyResult()

    def update_role_description(self) -> ActionResult:
        role_name = self._get_param("RoleName")
        description = self._get_param("Description")
        role = self.backend.update_role_description(role_name, description)
        result = {"Role": role}
        return ActionResult(result)

    def update_role(self) -> ActionResult:
        role_name = self._get_param("RoleName")
        description = self._get_param("Description")
        max_session_duration = self._get_param("MaxSessionDuration", 3600)
        self.backend.update_role(role_name, description, max_session_duration)
        return EmptyResult()

    def put_role_permissions_boundary(self) -> ActionResult:
        permissions_boundary = self._get_param("PermissionsBoundary")
        role_name = self._get_param("RoleName")
        self.backend.put_role_permissions_boundary(role_name, permissions_boundary)
        return EmptyResult()

    def delete_role_permissions_boundary(self) -> ActionResult:
        role_name = self._get_param("RoleName")
        self.backend.delete_role_permissions_boundary(role_name)
        return EmptyResult()

    def create_policy_version(self) -> ActionResult:
        policy_arn = self._get_param("PolicyArn")
        policy_document = self._get_param("PolicyDocument")
        set_as_default = self._get_param("SetAsDefault", False)
        policy_version = self.backend.create_policy_version(
            policy_arn, policy_document, set_as_default
        )
        result = {"PolicyVersion": policy_version}
        return ActionResult(result)

    def get_policy_version(self) -> ActionResult:
        policy_arn = self._get_param("PolicyArn")
        version_id = self._get_param("VersionId")
        policy_version = self.backend.get_policy_version(policy_arn, version_id)
        result = {"PolicyVersion": policy_version}
        return ActionResult(result)

    def list_policy_versions(self) -> ActionResult:
        policy_arn = self._get_param("PolicyArn")
        policy_versions = self.backend.list_policy_versions(policy_arn)
        result = {"Versions": policy_versions}
        return ActionResult(result)

    def list_policy_tags(self) -> ActionResult:
        policy_arn = self._get_param("PolicyArn")
        marker = self._get_param("Marker")
        max_items = self._get_param("MaxItems", 100)

        tags, marker = self.backend.list_policy_tags(policy_arn, marker, max_items)
        result = {"Tags": tags, "Marker": marker, "IsTruncated": bool(marker)}
        return ActionResult(result)

    def tag_policy(self) -> ActionResult:
        policy_arn = self._get_param("PolicyArn")
        tags = self._get_param("Tags", [])

        self.backend.tag_policy(policy_arn, tags)

        return EmptyResult()

    def untag_policy(self) -> ActionResult:
        policy_arn = self._get_param("PolicyArn")
        tag_keys = self._get_param("TagKeys", [])

        self.backend.untag_policy(policy_arn, tag_keys)

        return EmptyResult()

    def delete_policy_version(self) -> ActionResult:
        policy_arn = self._get_param("PolicyArn")
        version_id = self._get_param("VersionId")

        self.backend.delete_policy_version(policy_arn, version_id)
        return EmptyResult()

    def create_instance_profile(self) -> ActionResult:
        profile_name = self._get_param("InstanceProfileName")
        path = self._get_param("Path", "/")
        tags = self._get_param("Tags", [])

        profile = self.backend.create_instance_profile(
            profile_name, path, role_names=[], tags=tags
        )
        result = {"InstanceProfile": profile}
        return ActionResult(result)

    def delete_instance_profile(self) -> ActionResult:
        profile_name = self._get_param("InstanceProfileName")

        self.backend.delete_instance_profile(profile_name)
        return EmptyResult()

    def get_instance_profile(self) -> ActionResult:
        profile_name = self._get_param("InstanceProfileName")
        profile = self.backend.get_instance_profile(profile_name)

        result = {"InstanceProfile": profile}
        return ActionResult(result)

    def add_role_to_instance_profile(self) -> ActionResult:
        profile_name = self._get_param("InstanceProfileName")
        role_name = self._get_param("RoleName")

        self.backend.add_role_to_instance_profile(profile_name, role_name)
        return EmptyResult()

    def remove_role_from_instance_profile(self) -> ActionResult:
        profile_name = self._get_param("InstanceProfileName")
        role_name = self._get_param("RoleName")

        self.backend.remove_role_from_instance_profile(profile_name, role_name)
        return EmptyResult()

    def list_roles(self) -> ActionResult:
        path_prefix = self._get_param("PathPrefix", "/")
        marker = self._get_param("Marker", "0")
        max_items = self._get_param("MaxItems", 100)

        roles, marker = self.backend.list_roles(path_prefix, marker, max_items)
        result = {"Roles": roles, "IsTruncated": bool(marker), "Marker": marker}
        return ActionResult(result)

    def list_instance_profiles(self) -> ActionResult:
        profiles = self.backend.get_instance_profiles()
        result = {"InstanceProfiles": profiles}
        return ActionResult(result)

    def list_instance_profiles_for_role(self) -> ActionResult:
        role_name = self._get_param("RoleName")
        profiles = self.backend.get_instance_profiles_for_role(role_name=role_name)
        result = {"InstanceProfiles": profiles}
        return ActionResult(result)

    def upload_server_certificate(self) -> ActionResult:
        cert_name = self._get_param("ServerCertificateName")
        cert_body = self._get_param("CertificateBody")
        path = self._get_param("Path")
        private_key = self._get_param("PrivateKey")
        cert_chain = self._get_param("CertificateChain")

        cert = self.backend.upload_server_certificate(
            cert_name, cert_body, private_key, cert_chain=cert_chain, path=path
        )
        result = {
            "ServerCertificateMetadata": {
                "ServerCertificateName": cert.cert_name,
                "Path": cert.path,
                "Arn": cert.arn,
                "UploadDate": "2010-05-08T01:02:03.004Z",
                "ServerCertificateId": "ASCACKCEVSQ6C2EXAMPLE",
                "Expiration": "2012-05-08T01:02:03.004Z",
            }
        }
        return ActionResult(result)

    def list_server_certificates(self) -> ActionResult:
        certs = self.backend.list_server_certificates()
        result = {
            "ServerCertificateMetadataList": [
                {
                    "ServerCertificateName": cert.cert_name,
                    "Path": cert.path,
                    "Arn": cert.arn,
                    "UploadDate": "2010-05-08T01:02:03.004Z",
                    "ServerCertificateId": "ASCACKCEVSQ6C2EXAMPLE",
                    "Expiration": "2012-05-08T01:02:03.004Z",
                }
                for cert in certs
            ],
            "IsTruncated": False,
        }
        return ActionResult(result)

    def get_server_certificate(self) -> ActionResult:
        cert_name = self._get_param("ServerCertificateName")
        cert = self.backend.get_server_certificate(cert_name)
        server_cert = {
            "ServerCertificateMetadata": {
                "ServerCertificateName": cert.cert_name,
                "Path": cert.path,
                "Arn": cert.arn,
                "UploadDate": "2010-05-08T01:02:03.004Z",
                "ServerCertificateId": "ASCACKCEVSQ6C2EXAMPLE",
                "Expiration": "2012-05-08T01:02:03.004Z",
            },
            "CertificateBody": cert.cert_body,
        }
        if cert.cert_chain:
            server_cert["CertificateChain"] = cert.cert_chain
        result = {"ServerCertificate": server_cert}
        return ActionResult(result)

    def delete_server_certificate(self) -> ActionResult:
        cert_name = self._get_param("ServerCertificateName")
        self.backend.delete_server_certificate(cert_name)
        return EmptyResult()

    def create_group(self) -> ActionResult:
        group_name = self._get_param("GroupName")
        path = self._get_param("Path", "/")

        group = self.backend.create_group(group_name, path)
        result = {"Group": group}
        return ActionResult(result)

    def get_group(self) -> ActionResult:
        group_name = self._get_param("GroupName")

        group = self.backend.get_group(group_name)
        result = {"Group": group, "Users": group.users}
        return ActionResult(result)

    def list_groups(self) -> ActionResult:
        groups = self.backend.list_groups()
        result = {"Groups": groups}
        return ActionResult(result)

    def list_groups_for_user(self) -> ActionResult:
        user_name = self._get_param("UserName")

        groups = self.backend.get_groups_for_user(user_name)
        result = {"Groups": groups, "IsTruncated": False}
        return ActionResult(result)

    def put_group_policy(self) -> ActionResult:
        group_name = self._get_param("GroupName")
        policy_name = self._get_param("PolicyName")
        policy_document = self._get_param("PolicyDocument")
        self.backend.put_group_policy(group_name, policy_name, policy_document)
        return EmptyResult()

    def list_group_policies(self) -> ActionResult:
        group_name = self._get_param("GroupName")
        policies = self.backend.list_group_policies(group_name)
        result = {"PolicyNames": policies}
        return ActionResult(result)

    def get_group_policy(self) -> ActionResult:
        group_name = self._get_param("GroupName")
        policy_name = self._get_param("PolicyName")
        policy_result = self.backend.get_group_policy(group_name, policy_name)
        result = {
            "PolicyName": policy_result["policy_name"],
            "PolicyDocument": policy_result["policy_document"],
            "GroupName": group_name,
        }
        return ActionResult(result)

    def delete_group_policy(self) -> ActionResult:
        group_name = self._get_param("GroupName")
        policy_name = self._get_param("PolicyName")
        self.backend.delete_group_policy(group_name, policy_name)
        return EmptyResult()

    def delete_group(self) -> ActionResult:
        group_name = self._get_param("GroupName")
        self.backend.delete_group(group_name)
        return EmptyResult()

    def update_group(self) -> ActionResult:
        group_name = self._get_param("GroupName")
        new_group_name = self._get_param("NewGroupName")
        new_path = self._get_param("NewPath")
        self.backend.update_group(group_name, new_group_name, new_path)
        return EmptyResult()

    def create_user(self) -> ActionResult:
        user_name = self._get_param("UserName")
        path = self._get_param("Path")
        tags = self._get_param("Tags", [])
        user = self.backend.create_user(
            self.region, user_name=user_name, path=path, tags=tags
        )
        result = {"User": user}
        return ActionResult(result)

    def get_user(self) -> ActionResult:
        user_name = self._get_param("UserName")
        if not user_name:
            access_key_id = self.get_access_key()
            user = self.backend.get_user_from_access_key_id(access_key_id)
            if user is None:
                user = User(
                    self.current_account, region_name=self.region, name="default_user"
                )
        else:
            user = self.backend.get_user(user_name)
        result = {"User": user}
        return ActionResult(result)

    def list_users(self) -> ActionResult:
        path_prefix = self._get_param("PathPrefix")
        marker = self._get_param("Marker")
        max_items = self._get_param("MaxItems")
        users = self.backend.list_users(path_prefix, marker, max_items)
        result = {"Users": users, "IsTruncated": False}
        return ActionResult(result)

    def update_user(self) -> ActionResult:
        user_name = self._get_param("UserName")
        new_path = self._get_param("NewPath")
        new_user_name = self._get_param("NewUserName")
        self.backend.update_user(user_name, new_path, new_user_name)
        return EmptyResult()

    def create_login_profile(self) -> ActionResult:
        user_name = self._get_param("UserName")
        password = self._get_param("Password")
        user = self.backend.create_login_profile(user_name, password)
        result = {
            "LoginProfile": {
                "UserName": user.name,
                "CreateDate": user.create_date,
                "PasswordResetRequired": user.password_reset_required,
            }
        }
        return ActionResult(result)

    def get_login_profile(self) -> ActionResult:
        user_name = self._get_param("UserName")
        user = self.backend.get_login_profile(user_name)
        result = {
            "LoginProfile": {
                "UserName": user.name,
                "CreateDate": user.create_date,
                "PasswordResetRequired": user.password_reset_required,
            }
        }
        return ActionResult(result)

    def update_login_profile(self) -> ActionResult:
        user_name = self._get_param("UserName")
        password = self._get_param("Password")
        password_reset_required = self._get_param("PasswordResetRequired")
        self.backend.update_login_profile(user_name, password, password_reset_required)
        return EmptyResult()

    def add_user_to_group(self) -> ActionResult:
        group_name = self._get_param("GroupName")
        user_name = self._get_param("UserName")

        self.backend.add_user_to_group(group_name, user_name)
        return EmptyResult()

    def remove_user_from_group(self) -> ActionResult:
        group_name = self._get_param("GroupName")
        user_name = self._get_param("UserName")

        self.backend.remove_user_from_group(group_name, user_name)
        return EmptyResult()

    def get_user_policy(self) -> ActionResult:
        user_name = self._get_param("UserName")
        policy_name = self._get_param("PolicyName")

        policy = self.backend.get_user_policy(user_name, policy_name)
        result = {
            "PolicyName": policy["policy_name"],
            "PolicyDocument": policy["policy_document"],
            "UserName": policy["user_name"],
        }
        return ActionResult(result)

    def list_user_policies(self) -> ActionResult:
        user_name = self._get_param("UserName")
        policies = self.backend.list_user_policies(user_name)
        result = {"PolicyNames": policies, "IsTruncated": False}
        return ActionResult(result)

    def list_user_tags(self) -> ActionResult:
        user_name = self._get_param("UserName")
        tags = self.backend.list_user_tags(user_name)
        result = {"Tags": tags, "IsTruncated": False}
        return ActionResult(result)

    def put_user_policy(self) -> ActionResult:
        user_name = self._get_param("UserName")
        policy_name = self._get_param("PolicyName")
        policy_document = self._get_param("PolicyDocument")

        self.backend.put_user_policy(user_name, policy_name, policy_document)
        return EmptyResult()

    def delete_user_policy(self) -> ActionResult:
        user_name = self._get_param("UserName")
        policy_name = self._get_param("PolicyName")

        self.backend.delete_user_policy(user_name, policy_name)
        return EmptyResult()

    def create_access_key(self) -> ActionResult:
        user_name = self._get_param("UserName")
        if not user_name:
            access_key_id = self.get_access_key()
            access_key = self.backend.get_access_key_last_used(access_key_id)
            user_name = access_key["user_name"]

        key = self.backend.create_access_key(user_name)
        result = {"AccessKey": key}
        return ActionResult(result)

    def update_access_key(self) -> ActionResult:
        user_name = self._get_param("UserName")
        access_key_id = self._get_param("AccessKeyId")
        status = self._get_param("Status")
        if not user_name:
            access_key = self.backend.get_access_key_last_used(access_key_id)
            user_name = access_key["user_name"]

        self.backend.update_access_key(user_name, access_key_id, status)
        return EmptyResult()

    def get_access_key_last_used(self) -> ActionResult:
        access_key_id = self._get_param("AccessKeyId")
        last_used_response = self.backend.get_access_key_last_used(access_key_id)
        user_name = last_used_response["user_name"]
        last_used = last_used_response["last_used"]
        if last_used:
            last_used = {
                "LastUsedDate": last_used.timestamp,
                "Region": last_used.region,
                "ServiceName": last_used.service,
            }
        else:
            last_used = {"LastUsedDate": None, "Region": "N/A", "ServiceName": "N/A"}
        result = {
            "UserName": user_name,
            "AccessKeyLastUsed": last_used,
        }
        return ActionResult(result)

    def list_access_keys(self) -> ActionResult:
        user_name = self._get_param("UserName")
        if not user_name:
            access_key_id = self.get_access_key()
            access_key = self.backend.get_access_key_last_used(access_key_id)
            user_name = access_key["user_name"]
        keys = self.backend.list_access_keys(user_name)
        result = {"UserName": user_name, "AccessKeyMetadata": keys}
        return ActionResult(result)

    def delete_access_key(self) -> ActionResult:
        user_name = self._get_param("UserName")
        access_key_id = self._get_param("AccessKeyId")
        if not user_name:
            access_key = self.backend.get_access_key_last_used(access_key_id)
            user_name = access_key["user_name"]

        self.backend.delete_access_key(access_key_id, user_name)
        return EmptyResult()

    def upload_ssh_public_key(self) -> ActionResult:
        user_name = self._get_param("UserName")
        ssh_public_key_body = self._get_param("SSHPublicKeyBody")

        key = self.backend.upload_ssh_public_key(user_name, ssh_public_key_body)
        result = {"SSHPublicKey": key}
        return ActionResult(result)

    def get_ssh_public_key(self) -> ActionResult:
        user_name = self._get_param("UserName")
        ssh_public_key_id = self._get_param("SSHPublicKeyId")

        key = self.backend.get_ssh_public_key(user_name, ssh_public_key_id)
        result = {"SSHPublicKey": key}
        return ActionResult(result)

    def list_ssh_public_keys(self) -> ActionResult:
        user_name = self._get_param("UserName")

        keys = self.backend.get_all_ssh_public_keys(user_name)
        result = {"SSHPublicKeys": keys}
        return ActionResult(result)

    def update_ssh_public_key(self) -> ActionResult:
        user_name = self._get_param("UserName")
        ssh_public_key_id = self._get_param("SSHPublicKeyId")
        status = self._get_param("Status")

        self.backend.update_ssh_public_key(user_name, ssh_public_key_id, status)
        return EmptyResult()

    def delete_ssh_public_key(self) -> ActionResult:
        user_name = self._get_param("UserName")
        ssh_public_key_id = self._get_param("SSHPublicKeyId")

        self.backend.delete_ssh_public_key(user_name, ssh_public_key_id)
        return EmptyResult()

    def deactivate_mfa_device(self) -> ActionResult:
        user_name = self._get_param("UserName")
        serial_number = self._get_param("SerialNumber")

        self.backend.deactivate_mfa_device(user_name, serial_number)
        return EmptyResult()

    def enable_mfa_device(self) -> ActionResult:
        user_name = self._get_param("UserName")
        serial_number = self._get_param("SerialNumber")
        authentication_code_1 = self._get_param("AuthenticationCode1")
        authentication_code_2 = self._get_param("AuthenticationCode2")

        self.backend.enable_mfa_device(
            user_name, serial_number, authentication_code_1, authentication_code_2
        )
        return EmptyResult()

    def list_mfa_devices(self) -> ActionResult:
        user_name = self._get_param("UserName")
        devices = self.backend.list_mfa_devices(user_name)
        result = {
            "MFADevices": [
                {
                    "UserName": user_name,
                    "SerialNumber": device.serial_number,
                    "EnableDate": device.enable_date,
                }
                for device in devices
            ],
            "IsTruncated": False,
        }
        return ActionResult(result)

    def create_virtual_mfa_device(self) -> ActionResult:
        path = self._get_param("Path")
        virtual_mfa_device_name = self._get_param("VirtualMFADeviceName")

        virtual_mfa_device = self.backend.create_virtual_mfa_device(
            virtual_mfa_device_name, path
        )

        result = {"VirtualMFADevice": virtual_mfa_device}
        return ActionResult(result)

    def delete_virtual_mfa_device(self) -> ActionResult:
        serial_number = self._get_param("SerialNumber")

        self.backend.delete_virtual_mfa_device(serial_number)

        return EmptyResult()

    def list_virtual_mfa_devices(self) -> ActionResult:
        assignment_status = self._get_param("AssignmentStatus", "Any")
        marker = self._get_param("Marker")
        max_items = self._get_param("MaxItems", 100)

        devices, marker = self.backend.list_virtual_mfa_devices(
            assignment_status, marker, max_items
        )

        result = {
            "VirtualMFADevices": devices,
            "IsTruncated": bool(marker),
            "Marker": marker,
        }
        return ActionResult(result)

    def delete_user(self) -> ActionResult:
        user_name = self._get_param("UserName")
        self.backend.delete_user(user_name)
        return EmptyResult()

    def delete_policy(self) -> ActionResult:
        policy_arn = self._get_param("PolicyArn")
        self.backend.delete_policy(policy_arn)
        return EmptyResult()

    def delete_login_profile(self) -> ActionResult:
        user_name = self._get_param("UserName")
        self.backend.delete_login_profile(user_name)
        return EmptyResult()

    def generate_credential_report(self) -> ActionResult:
        if self.backend.report_generated():
            result = {"State": "COMPLETE"}
        else:
            result = {
                "State": "STARTED",
                "Description": "No report exists. Starting a new report generation task",
            }
        self.backend.generate_report()
        return ActionResult(result)

    def get_credential_report(self) -> ActionResult:
        report = self.backend.get_credential_report()
        result = {
            "Content": report,
            "ReportFormat": "text/csv",
            "GeneratedTime": "2015-02-02T20:02:02Z",
        }
        return ActionResult(result)

    def list_account_aliases(self) -> ActionResult:
        aliases = self.backend.list_account_aliases()
        result = {"AccountAliases": aliases}
        return ActionResult(result)

    def create_account_alias(self) -> ActionResult:
        alias = self._get_param("AccountAlias")
        self.backend.create_account_alias(alias)
        return EmptyResult()

    def delete_account_alias(self) -> ActionResult:
        self.backend.delete_account_alias()
        return EmptyResult()

    def get_account_authorization_details(self) -> ActionResult:
        filter_param = self._get_param("Filter", [])
        account_details = self.backend.get_account_authorization_details(filter_param)
        result = {
            "UserDetailList": account_details["users"],
            "GroupDetailList": account_details["groups"],
            "RoleDetailList": account_details["roles"],
            "Policies": account_details["managed_policies"],
            "IsTruncated": False,
        }
        return ActionResult(result)

    def create_saml_provider(self) -> ActionResult:
        saml_provider_name = self._get_param("Name")
        saml_metadata_document = self._get_param("SAMLMetadataDocument")
        saml_provider = self.backend.create_saml_provider(
            saml_provider_name, saml_metadata_document
        )
        result = {"SAMLProviderArn": saml_provider.arn}
        return ActionResult(result)

    def update_saml_provider(self) -> ActionResult:
        saml_provider_arn = self._get_param("SAMLProviderArn")
        saml_metadata_document = self._get_param("SAMLMetadataDocument")
        saml_provider = self.backend.update_saml_provider(
            saml_provider_arn, saml_metadata_document
        )
        result = {"SAMLProviderArn": saml_provider.arn}
        return ActionResult(result)

    def delete_saml_provider(self) -> ActionResult:
        saml_provider_arn = self._get_param("SAMLProviderArn")
        self.backend.delete_saml_provider(saml_provider_arn)

        return EmptyResult()

    def list_saml_providers(self) -> ActionResult:
        saml_providers = self.backend.list_saml_providers()

        result = {"SAMLProviderList": saml_providers}
        return ActionResult(result)

    def get_saml_provider(self) -> ActionResult:
        saml_provider_arn = self._get_param("SAMLProviderArn")
        saml_provider = self.backend.get_saml_provider(saml_provider_arn)

        return ActionResult(saml_provider)

    def upload_signing_certificate(self) -> ActionResult:
        user_name = self._get_param("UserName")
        cert_body = self._get_param("CertificateBody")

        cert = self.backend.upload_signing_certificate(user_name, cert_body)
        result = {"Certificate": cert}
        return ActionResult(result)

    def update_signing_certificate(self) -> ActionResult:
        user_name = self._get_param("UserName")
        cert_id = self._get_param("CertificateId")
        status = self._get_param("Status")

        self.backend.update_signing_certificate(user_name, cert_id, status)
        return EmptyResult()

    def delete_signing_certificate(self) -> ActionResult:
        user_name = self._get_param("UserName")
        cert_id = self._get_param("CertificateId")

        self.backend.delete_signing_certificate(user_name, cert_id)
        return EmptyResult()

    def list_signing_certificates(self) -> ActionResult:
        user_name = self._get_param("UserName")

        certs = self.backend.list_signing_certificates(user_name)
        result = {"UserName": user_name, "Certificates": certs}
        return ActionResult(result)

    def list_role_tags(self) -> ActionResult:
        role_name = self._get_param("RoleName")
        marker = self._get_param("Marker")
        max_items = self._get_param("MaxItems", 100)

        tags, marker = self.backend.list_role_tags(role_name, marker, max_items)

        result = {
            "Tags": tags,
            "Marker": marker,
            "IsTruncated": bool(marker),
        }
        return ActionResult(result)

    def tag_role(self) -> ActionResult:
        role_name = self._get_param("RoleName")
        tags = self._get_param("Tags", [])

        self.backend.tag_role(role_name, tags)

        return EmptyResult()

    def untag_role(self) -> ActionResult:
        role_name = self._get_param("RoleName")
        tag_keys = self._get_param("TagKeys", [])

        self.backend.untag_role(role_name, tag_keys)

        return EmptyResult()

    def create_open_id_connect_provider(self) -> ActionResult:
        open_id_provider_url = self._get_param("Url")
        thumbprint_list = self._get_param("ThumbprintList", [])
        client_id_list = self._get_param("ClientIDList", [])
        tags = self._get_param("Tags", [])

        open_id_provider = self.backend.create_open_id_connect_provider(
            open_id_provider_url, thumbprint_list, client_id_list, tags
        )

        result = {"OpenIDConnectProviderArn": open_id_provider.arn}
        return ActionResult(result)

    def update_open_id_connect_provider_thumbprint(self) -> ActionResult:
        open_id_provider_arn = self._get_param("OpenIDConnectProviderArn")
        thumbprint_list = self._get_param("ThumbprintList", [])

        self.backend.update_open_id_connect_provider_thumbprint(
            open_id_provider_arn, thumbprint_list
        )

        return EmptyResult()

    def tag_open_id_connect_provider(self) -> ActionResult:
        open_id_provider_arn = self._get_param("OpenIDConnectProviderArn")
        tags = self._get_param("Tags", [])

        self.backend.tag_open_id_connect_provider(open_id_provider_arn, tags)

        return EmptyResult()

    def untag_open_id_connect_provider(self) -> ActionResult:
        open_id_provider_arn = self._get_param("OpenIDConnectProviderArn")
        tag_keys = self._get_param("TagKeys", [])

        self.backend.untag_open_id_connect_provider(open_id_provider_arn, tag_keys)

        return EmptyResult()

    def list_open_id_connect_provider_tags(self) -> ActionResult:
        open_id_provider_arn = self._get_param("OpenIDConnectProviderArn")
        marker = self._get_param("Marker")
        max_items = self._get_param("MaxItems", 100)
        tags, marker = self.backend.list_open_id_connect_provider_tags(
            open_id_provider_arn, marker, max_items
        )
        result = {"Tags": tags, "Marker": marker, "IsTruncated": bool(marker)}
        return ActionResult(result)

    def delete_open_id_connect_provider(self) -> ActionResult:
        open_id_provider_arn = self._get_param("OpenIDConnectProviderArn")

        self.backend.delete_open_id_connect_provider(open_id_provider_arn)

        return EmptyResult()

    def get_open_id_connect_provider(self) -> ActionResult:
        open_id_provider_arn = self._get_param("OpenIDConnectProviderArn")

        open_id_provider = self.backend.get_open_id_connect_provider(
            open_id_provider_arn
        )

        return ActionResult(open_id_provider)

    def list_open_id_connect_providers(self) -> ActionResult:
        open_id_provider_arns = self.backend.list_open_id_connect_providers()

        result = {
            "OpenIDConnectProviderList": [{"Arn": arn} for arn in open_id_provider_arns]
        }
        return ActionResult(result)

    def update_account_password_policy(self) -> ActionResult:
        allow_change_password = self._get_bool_param(
            "AllowUsersToChangePassword", False
        )
        hard_expiry = self._get_bool_param("HardExpiry", False)
        max_password_age = self._get_int_param("MaxPasswordAge", 0)
        minimum_password_length = self._get_int_param("MinimumPasswordLength", 6)
        password_reuse_prevention = self._get_int_param("PasswordReusePrevention")
        require_lowercase_characters = self._get_bool_param(
            "RequireLowercaseCharacters", False
        )
        require_numbers = self._get_bool_param("RequireNumbers", False)
        require_symbols = self._get_bool_param("RequireSymbols", False)
        require_uppercase_characters = self._get_bool_param(
            "RequireUppercaseCharacters", False
        )

        self.backend.update_account_password_policy(
            allow_change_password,
            hard_expiry,
            max_password_age,
            minimum_password_length,
            password_reuse_prevention,
            require_lowercase_characters,
            require_numbers,
            require_symbols,
            require_uppercase_characters,
        )

        return EmptyResult()

    def get_account_password_policy(self) -> ActionResult:
        account_password_policy = self.backend.get_account_password_policy()

        result = {"PasswordPolicy": account_password_policy}
        return ActionResult(result)

    def delete_account_password_policy(self) -> ActionResult:
        self.backend.delete_account_password_policy()

        return EmptyResult()

    def get_account_summary(self) -> ActionResult:
        account_summary = self.backend.get_account_summary()

        result = {"SummaryMap": account_summary.summary_map}
        return ActionResult(result)

    def tag_user(self) -> ActionResult:
        name = self._get_param("UserName")
        tags = self._get_param("Tags", [])

        self.backend.tag_user(name, tags)

        return EmptyResult()

    def untag_user(self) -> ActionResult:
        name = self._get_param("UserName")
        tag_keys = self._get_param("TagKeys", [])

        self.backend.untag_user(name, tag_keys)

        return EmptyResult()

    def create_service_linked_role(self) -> ActionResult:
        service_name = self._get_param("AWSServiceName")
        description = self._get_param("Description")
        suffix = self._get_param("CustomSuffix")

        role = self.backend.create_service_linked_role(
            service_name, description, suffix
        )

        result = {"Role": role}
        return ActionResult(result)

    def delete_service_linked_role(self) -> ActionResult:
        role_name = self._get_param("RoleName")

        deletion_task_id = self.backend.delete_service_linked_role(role_name)

        result = {"DeletionTaskId": deletion_task_id}
        return ActionResult(result)

    def get_service_linked_role_deletion_status(self) -> ActionResult:
        self.backend.get_service_linked_role_deletion_status()

        result = {"Status": "SUCCEEDED"}
        return ActionResult(result)

    def tag_instance_profile(self) -> ActionResult:
        instance_profile_name = self._get_param("InstanceProfileName")
        tags = self._get_param("Tags", [])

        self.backend.tag_instance_profile(
            instance_profile_name=instance_profile_name,
            tags=tags,
        )
        return EmptyResult()

    def untag_instance_profile(self) -> ActionResult:
        instance_profile_name = self._get_param("InstanceProfileName")
        tags = self._get_param("TagKeys", [])

        self.backend.untag_instance_profile(
            instance_profile_name=instance_profile_name,
            tagKeys=tags,
        )
        return EmptyResult()
