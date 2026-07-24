from typing import Any

from moto.core.responses import (
    ActionResult,
    BaseResponse,
    EmptyResult,
)

from .models import S3ControlBackend, s3control_backends


class S3ControlResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="s3control")
        self.automated_parameter_parsing = True

    @property
    def backend(self) -> S3ControlBackend:
        return s3control_backends[self.current_account][self.partition]

    def get_public_access_block(self) -> ActionResult:
        account_id = self._get_param("AccountId")
        public_block_config = self.backend.get_public_access_block(
            account_id=account_id
        )
        return ActionResult({"PublicAccessBlockConfiguration": public_block_config})

    def put_public_access_block(self) -> EmptyResult:
        account_id = self._get_param("AccountId")
        pab_config = self._get_param("PublicAccessBlockConfiguration", {})
        self.backend.put_public_access_block(account_id, pab_config)
        return EmptyResult()

    def delete_public_access_block(self) -> EmptyResult:
        account_id = self._get_param("AccountId")
        self.backend.delete_public_access_block(account_id=account_id)
        return EmptyResult()

    def create_access_point(self) -> ActionResult:
        account_id = self._get_param("AccountId")
        name = self._get_param("Name")
        bucket = self._get_param("Bucket")
        vpc_configuration = self._get_param("VpcConfiguration")
        public_access_block_configuration = self._get_param(
            "PublicAccessBlockConfiguration"
        )
        access_point = self.backend.create_access_point(
            account_id=account_id,
            name=name,
            bucket=bucket,
            vpc_configuration=vpc_configuration,
            public_access_block_configuration=public_access_block_configuration,
        )
        return ActionResult(
            {
                "AccessPointArn": access_point.arn,
                "Alias": access_point.alias,
            }
        )

    def get_access_point(self) -> ActionResult:
        account_id = self._get_param("AccountId")
        name = self._get_param("Name")
        access_point = self.backend.get_access_point(account_id=account_id, name=name)
        return ActionResult(
            {
                "Name": access_point.name,
                "Bucket": access_point.bucket,
                "NetworkOrigin": access_point.network_origin,
                "VpcConfiguration": {"VpcId": access_point.vpc_id}
                if access_point.vpc_id
                else None,
                "PublicAccessBlockConfiguration": access_point.pubc,
                "CreationDate": access_point.created,
                "Alias": access_point.alias,
                "AccessPointArn": access_point.arn,
                "Endpoints": {
                    "ipv4": "s3-accesspoint.us-east-1.amazonaws.com",
                    "fips": "s3-accesspoint-fips.us-east-1.amazonaws.com",
                    "fips_dualstack": "s3-accesspoint-fips.dualstack.us-east-1.amazonaws.com",
                    "dualstack": "s3-accesspoint.dualstack.us-east-1.amazonaws.com",
                },
            }
        )

    def delete_access_point(self) -> EmptyResult:
        account_id = self._get_param("AccountId")
        name = self._get_param("Name")
        self.backend.delete_access_point(account_id=account_id, name=name)
        return EmptyResult()

    def put_access_point_policy(self) -> EmptyResult:
        account_id = self._get_param("AccountId")
        name = self._get_param("Name")
        policy = self._get_param("Policy")
        self.backend.put_access_point_policy(account_id, name, policy)
        return EmptyResult()

    def get_access_point_policy(self) -> ActionResult:
        account_id = self._get_param("AccountId")
        name = self._get_param("Name")
        policy = self.backend.get_access_point_policy(account_id, name)
        return ActionResult({"Policy": policy})

    def delete_access_point_policy(self) -> EmptyResult:
        account_id = self._get_param("AccountId")
        name = self._get_param("Name")
        self.backend.delete_access_point_policy(account_id=account_id, name=name)
        return EmptyResult()

    def get_access_point_policy_status(self) -> ActionResult:
        account_id = self._get_param("AccountId")
        name = self._get_param("Name")
        self.backend.get_access_point_policy_status(account_id, name)
        return ActionResult({"PolicyStatus": {"IsPublic": True}})

    def put_storage_lens_configuration(self) -> EmptyResult:
        account_id = self._get_param("AccountId")
        config_id = self._get_param("ConfigId")
        storage_lens_configuration = self._get_param("StorageLensConfiguration")
        tags = self._get_param("Tags")
        self.backend.put_storage_lens_configuration(
            config_id=config_id,
            account_id=account_id,
            storage_lens_configuration=storage_lens_configuration,
            tags=tags,
        )
        return EmptyResult()

    def get_storage_lens_configuration(self) -> ActionResult:
        account_id = self._get_param("AccountId")
        config_id = self._get_param("ConfigId")
        storage_lens_configuration = self.backend.get_storage_lens_configuration(
            config_id=config_id,
            account_id=account_id,
        )
        return ActionResult(
            {"StorageLensConfiguration": storage_lens_configuration.config}
        )

    def delete_storage_lens_configuration(self) -> EmptyResult:
        account_id = self._get_param("AccountId")
        config_id = self._get_param("ConfigId")
        self.backend.delete_storage_lens_configuration(
            config_id=config_id,
            account_id=account_id,
        )
        return EmptyResult()

    def list_storage_lens_configurations(self) -> ActionResult:
        account_id = self._get_param("AccountId")
        next_token = self._get_param("NextToken")
        storage_lens_configuration_list, next_token = (
            self.backend.list_storage_lens_configurations(
                account_id=account_id,
                next_token=next_token,
            )
        )
        configs = [
            {
                "Id": config.config.get("Id"),
                "IsEnabled": config.config.get("IsEnabled"),
                "StorageLensArn": config.arn,
            }
            for config in storage_lens_configuration_list
        ]
        result = {"StorageLensConfigurationList": configs, "NextToken": next_token}
        return ActionResult(result)

    def put_storage_lens_configuration_tagging(self) -> EmptyResult:
        account_id = self._get_param("AccountId")
        config_id = self._get_param("ConfigId")
        tags = self._get_param("Tags")
        self.backend.put_storage_lens_configuration_tagging(
            config_id=config_id,
            account_id=account_id,
            tags=tags,
        )
        return EmptyResult()

    def get_storage_lens_configuration_tagging(self) -> ActionResult:
        account_id = self._get_param("AccountId")
        config_id = self._get_param("ConfigId")
        storage_lens_tags = self.backend.get_storage_lens_configuration_tagging(
            config_id=config_id,
            account_id=account_id,
        )
        return ActionResult({"Tags": storage_lens_tags})

    def list_access_points(self) -> ActionResult:
        account_id = self._get_param("AccountId")
        bucket = self._get_param("Bucket")
        max_results = self._get_int_param("MaxResults")
        next_token = self._get_param("NextToken")

        access_points, next_token = self.backend.list_access_points(
            account_id=account_id,
            bucket=bucket,
            max_results=max_results,
            next_token=next_token,
        )

        ap_list = [
            {
                "Name": ap.name,
                "NetworkOrigin": ap.network_origin,
                "VpcConfiguration": {"VpcId": ap.vpc_id} if ap.vpc_id else None,
                "Bucket": ap.bucket,
                "AccessPointArn": ap.arn,
                "Alias": ap.alias,
            }
            for ap in access_points
        ]
        result: dict[str, Any] = {"AccessPointList": ap_list}
        if next_token:
            result["NextToken"] = next_token
        return ActionResult(result)

    def create_multi_region_access_point(self) -> ActionResult:
        account_id = self._get_param("AccountId")
        name = self._get_param("Details.Name")
        regions = self._get_param("Details.Regions", [])
        public_access_block = self._get_param("Details.PublicAccessBlock", {})
        operation = self.backend.create_multi_region_access_point(
            account_id=account_id,
            name=name,
            public_access_block=public_access_block,
            regions=regions,
            region_name=self.region,
        )
        return ActionResult({"RequestTokenARN": operation.request_token_arn})

    def delete_multi_region_access_point(self) -> ActionResult:
        account_id = self._get_param("AccountId")
        details = self._get_param("Details", {})
        name = details.get("Name")

        operation = self.backend.delete_multi_region_access_point(
            account_id=account_id,
            name=name,
            region_name=self.region,
        )
        return ActionResult({"RequestTokenARN": operation.request_token_arn})

    def describe_multi_region_access_point_operation(self) -> ActionResult:
        account_id = self._get_param("AccountId")
        request_token = self._get_param("RequestTokenARN")

        operation = self.backend.describe_multi_region_access_point_operation(
            account_id=account_id,
            request_token_arn=request_token,
        )
        return ActionResult({"AsyncOperation": operation.to_dict()})

    def get_multi_region_access_point(self) -> ActionResult:
        account_id = self._get_param("AccountId")
        name = self._get_param("Name")

        mrap = self.backend.get_multi_region_access_point(
            account_id=account_id,
            name=name,
        )
        return ActionResult({"AccessPoint": mrap.to_dict()})

    def get_multi_region_access_point_policy(self) -> ActionResult:
        account_id = self._get_param("AccountId")
        name = self._get_param("Name")

        policy = self.backend.get_multi_region_access_point_policy(
            account_id=account_id,
            name=name,
        )
        return ActionResult({"Policy": {"Established": {"Policy": policy}}})

    def get_multi_region_access_point_policy_status(self) -> ActionResult:
        account_id = self._get_param("AccountId")
        name = self._get_param("Name")

        policy_status = self.backend.get_multi_region_access_point_policy_status(
            account_id=account_id,
            name=name,
        )
        return ActionResult({"Established": {"IsPublic": policy_status["IsPublic"]}})

    def list_multi_region_access_points(self) -> ActionResult:
        account_id = self._get_param("AccountId")
        max_results = self._get_int_param("MaxResults")
        next_token = self._get_param("NextToken")

        mraps, next_token = self.backend.list_multi_region_access_points(
            account_id=account_id,
            max_results=max_results,
            next_token=next_token,
        )

        result: dict[str, Any] = {
            "AccessPoints": [mrap.to_dict() for mrap in mraps],
        }
        if next_token:
            result["NextToken"] = next_token
        return ActionResult(result)

    def put_multi_region_access_point_policy(self) -> ActionResult:
        account_id = self._get_param("AccountId")
        details = self._get_param("Details", {})
        name = details.get("Name")
        policy = details.get("Policy")

        operation = self.backend.put_multi_region_access_point_policy(
            account_id=account_id,
            name=name,
            policy=policy,
            region_name=self.region,
        )
        return ActionResult({"RequestTokenARN": operation.request_token_arn})

    def list_tags_for_resource(self) -> ActionResult:
        resource_arn = self._get_param("ResourceArn")
        tags = self.backend.list_tags_for_resource(resource_arn)
        return ActionResult(result={"Tags": tags})

    def tag_resource(self) -> EmptyResult:
        resource_arn = self._get_param("ResourceArn")
        tags = self._get_param("Tags", [])
        self.backend.tag_resource(resource_arn, tags=tags)
        return EmptyResult()

    def untag_resource(self) -> EmptyResult:
        resource_arn = self._get_param("ResourceArn")
        tag_keys = self._get_param("TagKeys", [])
        self.backend.untag_resource(resource_arn, tag_keys=tag_keys)
        return EmptyResult()
