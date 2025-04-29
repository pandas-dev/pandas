from typing import Any, List, Optional, Tuple

from botocore.awsrequest import AWSPreparedRequest
from werkzeug.wrappers import Request

from moto.core.common_types import TYPE_RESPONSE
from moto.core.responses import ActionResult, BaseResponse
from moto.core.utils import get_service_model
from moto.ec2.models import ec2_backends

from .exceptions import DBParameterGroupNotFoundError
from .models import RDSBackend, rds_backends
from .parser import QueryParser, XFormedDict


def normalize_request(request: AWSPreparedRequest) -> Request:
    from urllib.parse import urlparse

    parsed_url = urlparse(request.url)
    normalized_request = Request.from_values(
        method=request.method,
        base_url=f"{parsed_url.scheme}://{parsed_url.netloc}",
        path=parsed_url.path,
        query_string=parsed_url.query,
        data=request.body,
        headers=[(k, v) for k, v in request.headers.items()],
    )
    return normalized_request


class RDSResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="rds")

    @property
    def backend(self) -> RDSBackend:
        return rds_backends[self.current_account][self.region]

    @property
    def global_backend(self) -> RDSBackend:
        """Return backend instance of the region that stores Global Clusters"""
        return rds_backends[self.current_account]["us-east-1"]

    def _dispatch(self, request: Any, full_url: str, headers: Any) -> TYPE_RESPONSE:
        self.setup_class(request, full_url, headers)

        if isinstance(request, AWSPreparedRequest):
            request = normalize_request(request)

        self.action = request.values["Action"]

        service_model = get_service_model(self.service_name)
        self.operation_model = service_model.operation_model(self.action)

        parser = QueryParser(map_type=XFormedDict)  # type: ignore[no-untyped-call]
        self.parameters = parser.parse(
            {"query_params": request.values},
            self.operation_model,  # type: ignore[no-untyped-call]
        )
        return self.call_action()

    def create_db_instance(self) -> ActionResult:
        db_kwargs = self.parameters
        database = self.backend.create_db_instance(db_kwargs)
        result = {"DBInstance": database}
        return ActionResult(result)

    def create_db_instance_read_replica(self) -> ActionResult:
        db_kwargs = self.parameters
        database = self.backend.create_db_instance_read_replica(db_kwargs)
        result = {"DBInstance": database}
        return ActionResult(result)

    def describe_db_instances(self) -> ActionResult:
        db_instance_identifier = self.parameters.get("DBInstanceIdentifier")
        filters = self.parameters.get("Filters", [])
        filter_dict = {f["Name"]: f["Values"] for f in filters}
        all_instances = list(
            self.backend.describe_db_instances(
                db_instance_identifier, filters=filter_dict
            )
        )
        instances_resp, next_marker = self._paginate(all_instances)
        result = {
            "DBInstances": instances_resp,
            "Marker": next_marker,
        }
        return ActionResult(result)

    def modify_db_instance(self) -> ActionResult:
        db_instance_identifier = self.parameters.get("DBInstanceIdentifier")
        db_kwargs = self.parameters
        # This is a hack because the backend code expects the parameter to be
        # lowercased before validation is performed.  We need to move the
        # lowercasing to a backend setter (in one place) and then do the validation.
        if "PreferredMaintenanceWindow" in db_kwargs:
            db_kwargs["PreferredMaintenanceWindow"] = db_kwargs[
                "PreferredMaintenanceWindow"
            ].lower()
        new_db_instance_identifier = self.parameters.get("NewDBInstanceIdentifier")
        if new_db_instance_identifier:
            db_kwargs["new_db_instance_identifier"] = new_db_instance_identifier
        database = self.backend.modify_db_instance(db_instance_identifier, db_kwargs)
        result = {"DBInstance": database}
        return ActionResult(result)

    def delete_db_instance(self) -> ActionResult:
        db_snapshot_name = self.parameters.get("FinalDBSnapshotIdentifier")
        if db_snapshot_name is not None:
            self.backend.validate_db_snapshot_identifier(
                db_snapshot_name, parameter_name="FinalDBSnapshotIdentifier"
            )
        database = self.backend.delete_db_instance(**self.parameters)
        result = {"DBInstance": database}
        return ActionResult(result)

    def reboot_db_instance(self) -> ActionResult:
        db_instance_identifier = self.parameters.get("DBInstanceIdentifier")
        database = self.backend.reboot_db_instance(db_instance_identifier)
        result = {"DBInstance": database}
        return ActionResult(result)

    def create_db_snapshot(self) -> ActionResult:
        self.backend.validate_db_snapshot_identifier(
            self.parameters["DBSnapshotIdentifier"],
            parameter_name="DBSnapshotIdentifier",
        )
        snapshot = self.backend.create_db_snapshot(**self.parameters)
        result = {"DBSnapshot": snapshot}
        return ActionResult(result)

    def copy_db_snapshot(self) -> ActionResult:
        target_snapshot_identifier = self.parameters.get("TargetDBSnapshotIdentifier")
        self.backend.validate_db_snapshot_identifier(
            target_snapshot_identifier, parameter_name="TargetDBSnapshotIdentifier"
        )
        snapshot = self.backend.copy_db_snapshot(**self.parameters)
        result = {"DBSnapshot": snapshot}
        return ActionResult(result)

    def describe_db_snapshots(self) -> ActionResult:
        db_instance_identifier = self.parameters.get("DBInstanceIdentifier")
        db_snapshot_identifier = self.parameters.get("DBSnapshotIdentifier")
        snapshot_type = self.parameters.get("SnapshotType")
        filters = self.parameters.get("Filters", [])
        filter_dict = {f["Name"]: f["Values"] for f in filters}
        snapshots = self.backend.describe_db_snapshots(
            db_instance_identifier, db_snapshot_identifier, snapshot_type, filter_dict
        )
        result = {"DBSnapshots": snapshots}
        return ActionResult(result)

    def promote_read_replica(self) -> ActionResult:
        db_kwargs = self.parameters
        database = self.backend.promote_read_replica(db_kwargs)
        result = {"DBInstance": database}
        return ActionResult(result)

    def delete_db_snapshot(self) -> ActionResult:
        db_snapshot_identifier = self.parameters.get("DBSnapshotIdentifier")
        snapshot = self.backend.delete_db_snapshot(db_snapshot_identifier)
        result = {"DBSnapshot": snapshot}
        return ActionResult(result)

    def restore_db_instance_from_db_snapshot(self) -> ActionResult:
        db_snapshot_identifier = self.parameters.get("DBSnapshotIdentifier")
        db_kwargs = self.parameters
        new_instance = self.backend.restore_db_instance_from_db_snapshot(
            db_snapshot_identifier, db_kwargs
        )
        result = {"DBInstance": new_instance}
        return ActionResult(result)

    def restore_db_instance_to_point_in_time(self) -> ActionResult:
        source_db_identifier = self.parameters.get("SourceDBInstanceIdentifier")
        target_db_identifier = self.parameters.get("TargetDBInstanceIdentifier")

        db_kwargs = self.parameters
        new_instance = self.backend.restore_db_instance_to_point_in_time(
            source_db_identifier, target_db_identifier, db_kwargs
        )
        result = {"DBInstance": new_instance}
        return ActionResult(result)

    def restore_db_cluster_to_point_in_time(self) -> ActionResult:
        cluster = self.backend.restore_db_cluster_to_point_in_time(**self.parameters)
        result = {"DBCluster": cluster}
        return ActionResult(result)

    def failover_db_cluster(self) -> ActionResult:
        cluster = self.backend.failover_db_cluster(**self.parameters)
        result = {"DBCluster": cluster}
        return ActionResult(result)

    def list_tags_for_resource(self) -> ActionResult:
        arn = self.parameters.get("ResourceName")
        tags = self.backend.list_tags_for_resource(arn)
        result = {"TagList": tags}
        return ActionResult(result)

    def add_tags_to_resource(self) -> ActionResult:
        arn = self.parameters.get("ResourceName")
        tags = self.parameters.get("Tags", [])
        self.backend.add_tags_to_resource(arn, tags)
        return ActionResult({})

    def remove_tags_from_resource(self) -> ActionResult:
        arn = self.parameters.get("ResourceName")
        tag_keys = self.parameters.get("TagKeys")
        self.backend.remove_tags_from_resource(arn, tag_keys)
        return ActionResult({})

    def stop_db_instance(self) -> ActionResult:
        db_instance_identifier = self.parameters.get("DBInstanceIdentifier")
        db_snapshot_identifier = self.parameters.get("DBSnapshotIdentifier")
        if db_snapshot_identifier is not None:
            self.backend.validate_db_snapshot_identifier(
                db_snapshot_identifier, parameter_name="DBSnapshotIdentifier"
            )

        database = self.backend.stop_db_instance(
            db_instance_identifier, db_snapshot_identifier
        )
        result = {"DBInstance": database}
        return ActionResult(result)

    def start_db_instance(self) -> ActionResult:
        db_instance_identifier = self.parameters.get("DBInstanceIdentifier")
        database = self.backend.start_db_instance(db_instance_identifier)
        result = {"DBInstance": database}
        return ActionResult(result)

    def create_db_security_group(self) -> ActionResult:
        group_name = self.parameters.get("DBSecurityGroupName")
        description = self.parameters.get("DBSecurityGroupDescription")
        tags = self.parameters.get("Tags", [])
        security_group = self.backend.create_db_security_group(
            group_name, description, tags
        )
        result = {"DBSecurityGroup": security_group}
        return ActionResult(result)

    def describe_db_security_groups(self) -> ActionResult:
        security_group_name = self.parameters.get("DBSecurityGroupName")
        security_groups = self.backend.describe_security_groups(security_group_name)
        result = {"DBSecurityGroups": security_groups}
        return ActionResult(result)

    def delete_db_security_group(self) -> ActionResult:
        security_group_name = self.parameters.get("DBSecurityGroupName")
        security_group = self.backend.delete_security_group(security_group_name)
        result = {"DBSecurityGroup": security_group}
        return ActionResult(result)

    def authorize_db_security_group_ingress(self) -> ActionResult:
        security_group_name = self.parameters.get("DBSecurityGroupName")
        cidr_ip = self.parameters.get("CIDRIP")
        security_group = self.backend.authorize_security_group(
            security_group_name, cidr_ip
        )
        result = {"DBSecurityGroup": security_group}
        return ActionResult(result)

    def create_db_subnet_group(self) -> ActionResult:
        subnet_name = self.parameters.get("DBSubnetGroupName")
        description = self.parameters.get("DBSubnetGroupDescription")
        subnet_ids = self.parameters.get("SubnetIds", [])
        tags = self.parameters.get("Tags", [])
        subnets = [
            ec2_backends[self.current_account][self.region].get_subnet(subnet_id)
            for subnet_id in subnet_ids
        ]
        subnet_group = self.backend.create_subnet_group(
            subnet_name, description, subnets, tags
        )
        result = {"DBSubnetGroup": subnet_group}
        return ActionResult(result)

    def describe_db_subnet_groups(self) -> ActionResult:
        subnet_name = self.parameters.get("DBSubnetGroupName")
        subnet_groups = self.backend.describe_db_subnet_groups(subnet_name)
        result = {"DBSubnetGroups": subnet_groups}
        return ActionResult(result)

    def modify_db_subnet_group(self) -> ActionResult:
        subnet_name = self.parameters.get("DBSubnetGroupName")
        description = self.parameters.get("DBSubnetGroupDescription")
        subnet_ids = self.parameters.get("SubnetIds", [])
        subnets = [
            ec2_backends[self.current_account][self.region].get_subnet(subnet_id)
            for subnet_id in subnet_ids
        ]
        subnet_group = self.backend.modify_db_subnet_group(
            subnet_name, description, subnets
        )
        result = {"DBSubnetGroup": subnet_group}
        return ActionResult(result)

    def delete_db_subnet_group(self) -> ActionResult:
        subnet_name = self.parameters.get("DBSubnetGroupName")
        subnet_group = self.backend.delete_subnet_group(subnet_name)
        result = {"DBSubnetGroup": subnet_group}
        return ActionResult(result)

    def create_option_group(self) -> ActionResult:
        kwargs = self.parameters
        option_group = self.backend.create_option_group(kwargs)
        result = {"OptionGroup": option_group}
        return ActionResult(result)

    def delete_option_group(self) -> ActionResult:
        name = self.parameters["OptionGroupName"]
        option_group = self.backend.delete_option_group(name)
        result = {"OptionGroup": option_group}
        return ActionResult(result)

    def describe_option_groups(self) -> ActionResult:
        kwargs = self.parameters
        option_groups = self.backend.describe_option_groups(kwargs)
        option_groups, marker = self._paginate(option_groups)
        result = {
            "OptionGroupsList": option_groups,
            "Marker": marker,
        }
        return ActionResult(result)

    def describe_option_group_options(self) -> str:
        engine_name = self.parameters.get("EngineName")
        major_engine_version = self.parameters.get("MajorEngineVersion")
        return self.backend.describe_option_group_options(
            engine_name, major_engine_version
        )

    def modify_option_group(self) -> ActionResult:
        option_group_name = self.parameters.get("OptionGroupName")
        options_to_include = self.parameters.get("OptionsToInclude", [])
        options_to_remove = self.parameters.get("OptionsToRemove", [])
        option_group = self.backend.modify_option_group(
            option_group_name, options_to_include, options_to_remove
        )
        result = {"OptionGroup": option_group}
        return ActionResult(result)

    def create_db_parameter_group(self) -> ActionResult:
        kwargs = self.parameters
        db_parameter_group = self.backend.create_db_parameter_group(kwargs)
        result = {"DBParameterGroup": db_parameter_group}
        return ActionResult(result)

    def copy_db_parameter_group(self) -> ActionResult:
        target_db_parameter_group = self.backend.copy_db_parameter_group(
            **self.parameters
        )
        result = {"DBParameterGroup": target_db_parameter_group}
        return ActionResult(result)

    def describe_db_parameter_groups(self) -> ActionResult:
        kwargs = self.parameters
        db_parameter_groups = self.backend.describe_db_parameter_groups(kwargs)
        db_parameter_groups, _ = self._paginate(db_parameter_groups)
        result = {"DBParameterGroups": db_parameter_groups}
        return ActionResult(result)

    def modify_db_parameter_group(self) -> ActionResult:
        db_parameter_group_name = self.parameters.get("DBParameterGroupName")
        param_list = self.parameters.get("Parameters", [])
        # Raw dict is stored on the backend, so we need the original PascalCase items.
        db_parameter_group_parameters = [
            dict(param.original_items()) for param in param_list
        ]
        db_parameter_group = self.backend.modify_db_parameter_group(
            db_parameter_group_name, db_parameter_group_parameters
        )
        result = {"DBParameterGroupName": db_parameter_group.name}
        return ActionResult(result)

    def describe_db_parameters(self) -> ActionResult:
        db_parameter_group_name = self.parameters.get("DBParameterGroupName")
        db_parameter_groups = self.backend.describe_db_parameter_groups(
            {"db_parameter_group_name": db_parameter_group_name}
        )
        if not db_parameter_groups:
            raise DBParameterGroupNotFoundError(db_parameter_group_name)
        parameters = db_parameter_groups[0].parameters.values()
        result = {"Parameters": parameters}
        return ActionResult(result)

    def delete_db_parameter_group(self) -> ActionResult:
        name = self.parameters["DBParameterGroupName"]
        db_parameter_group = self.backend.delete_db_parameter_group(name)
        return ActionResult(db_parameter_group)

    def describe_db_cluster_parameters(self) -> ActionResult:
        # TODO: This never worked at all...
        db_parameter_group_name = self.parameters.get("DBParameterGroupName")
        db_parameter_groups = self.backend.describe_db_cluster_parameters()
        if db_parameter_groups is None:
            raise DBParameterGroupNotFoundError(db_parameter_group_name)
        result = {"Parameters": db_parameter_groups}
        return ActionResult(result)

    def create_db_cluster(self) -> ActionResult:
        kwargs = self.parameters
        cluster = self.backend.create_db_cluster(kwargs)
        result = {"DBCluster": cluster}
        return ActionResult(result)

    def modify_db_cluster(self) -> ActionResult:
        kwargs = self.parameters
        cluster = self.backend.modify_db_cluster(kwargs)
        result = {"DBCluster": cluster}
        return ActionResult(result)

    def describe_db_clusters(self) -> ActionResult:
        _id = self.parameters.get("DBClusterIdentifier")
        filters = self.parameters.get("Filters", [])
        filter_dict = {f["Name"]: f["Values"] for f in filters}
        clusters = self.backend.describe_db_clusters(
            db_cluster_identifier=_id, filters=filter_dict
        )
        result = {"DBClusters": clusters}
        return ActionResult(result)

    def delete_db_cluster(self) -> ActionResult:
        _id = self.parameters.get("DBClusterIdentifier")
        snapshot_name = self.parameters.get("FinalDBSnapshotIdentifier")
        cluster = self.backend.delete_db_cluster(
            cluster_identifier=_id, snapshot_name=snapshot_name
        )
        result = {"DBCluster": cluster}
        return ActionResult(result)

    def start_db_cluster(self) -> ActionResult:
        _id = self.parameters.get("DBClusterIdentifier")
        cluster = self.backend.start_db_cluster(cluster_identifier=_id)
        result = {"DBCluster": cluster}
        return ActionResult(result)

    def stop_db_cluster(self) -> ActionResult:
        _id = self.parameters.get("DBClusterIdentifier")
        cluster = self.backend.stop_db_cluster(cluster_identifier=_id)
        result = {"DBCluster": cluster}
        return ActionResult(result)

    def create_db_cluster_snapshot(self) -> ActionResult:
        snapshot = self.backend.create_db_cluster_snapshot(**self.parameters)
        result = {"DBClusterSnapshot": snapshot}
        return ActionResult(result)

    def copy_db_cluster_snapshot(self) -> ActionResult:
        snapshot = self.backend.copy_db_cluster_snapshot(**self.parameters)
        result = {"DBClusterSnapshot": snapshot}
        return ActionResult(result)

    def describe_db_cluster_snapshots(self) -> ActionResult:
        db_cluster_identifier = self.parameters.get("DBClusterIdentifier")
        db_snapshot_identifier = self.parameters.get("DBClusterSnapshotIdentifier")
        snapshot_type = self.parameters.get("SnapshotType")
        filters = self.parameters.get("Filters", [])
        filter_dict = {f["Name"]: f["Values"] for f in filters}
        snapshots = self.backend.describe_db_cluster_snapshots(
            db_cluster_identifier, db_snapshot_identifier, snapshot_type, filter_dict
        )
        results = {"DBClusterSnapshots": snapshots}
        return ActionResult(results)

    def delete_db_cluster_snapshot(self) -> ActionResult:
        db_snapshot_identifier = self.parameters["DBClusterSnapshotIdentifier"]
        snapshot = self.backend.delete_db_cluster_snapshot(db_snapshot_identifier)
        result = {"DBClusterSnapshot": snapshot}
        return ActionResult(result)

    def restore_db_cluster_from_snapshot(self) -> ActionResult:
        db_snapshot_identifier = self.parameters.get("SnapshotIdentifier")
        db_kwargs = self.parameters
        new_cluster = self.backend.restore_db_cluster_from_snapshot(
            db_snapshot_identifier, db_kwargs
        )
        result = {"DBCluster": new_cluster}
        return ActionResult(result)

    def start_export_task(self) -> ActionResult:
        kwargs = self.parameters
        export_task = self.backend.start_export_task(kwargs)
        return ActionResult(export_task)

    def cancel_export_task(self) -> ActionResult:
        export_task_identifier = self.parameters.get("ExportTaskIdentifier")
        export_task = self.backend.cancel_export_task(export_task_identifier)
        return ActionResult(export_task)

    def describe_export_tasks(self) -> ActionResult:
        export_task_identifier = self.parameters.get("ExportTaskIdentifier")
        tasks = self.backend.describe_export_tasks(export_task_identifier)
        result = {"ExportTasks": tasks}
        return ActionResult(result)

    def create_event_subscription(self) -> ActionResult:
        kwargs = self.parameters
        subscription = self.backend.create_event_subscription(kwargs)
        result = {"EventSubscription": subscription}
        return ActionResult(result)

    def delete_event_subscription(self) -> ActionResult:
        subscription_name = self.parameters.get("SubscriptionName")
        subscription = self.backend.delete_event_subscription(subscription_name)
        result = {"EventSubscription": subscription}
        return ActionResult(result)

    def describe_event_subscriptions(self) -> ActionResult:
        subscription_name = self.parameters.get("SubscriptionName")
        subscriptions = self.backend.describe_event_subscriptions(subscription_name)
        result = {"EventSubscriptionsList": subscriptions}
        return ActionResult(result)

    def describe_orderable_db_instance_options(self) -> ActionResult:
        engine = self.parameters.get("Engine")
        engine_version = self.parameters.get("EngineVersion")
        options = self.backend.describe_orderable_db_instance_options(
            engine, engine_version
        )
        result = {"OrderableDBInstanceOptions": options}
        return ActionResult(result)

    def describe_global_clusters(self) -> ActionResult:
        clusters = self.global_backend.describe_global_clusters()
        result = {"GlobalClusters": clusters}
        return ActionResult(result)

    def create_global_cluster(self) -> ActionResult:
        params = self.parameters
        cluster = self.global_backend.create_global_cluster(
            global_cluster_identifier=params["GlobalClusterIdentifier"],
            source_db_cluster_identifier=params.get("SourceDBClusterIdentifier"),
            engine=params.get("Engine"),
            engine_version=params.get("EngineVersion"),
            storage_encrypted=params.get("StorageEncrypted"),
            deletion_protection=params.get("DeletionProtection"),
        )
        result = {"GlobalCluster": cluster}
        return ActionResult(result)

    def delete_global_cluster(self) -> ActionResult:
        params = self.parameters
        cluster = self.global_backend.delete_global_cluster(
            global_cluster_identifier=params["GlobalClusterIdentifier"],
        )
        result = {"GlobalCluster": cluster}
        return ActionResult(result)

    def remove_from_global_cluster(self) -> ActionResult:
        params = self.parameters
        global_cluster = self.backend.remove_from_global_cluster(
            global_cluster_identifier=params["GlobalClusterIdentifier"],
            db_cluster_identifier=params["DbClusterIdentifier"],
        )
        result = {"GlobalCluster": global_cluster}
        return ActionResult(result)

    def create_db_cluster_parameter_group(self) -> ActionResult:
        group_name = self.parameters.get("DBClusterParameterGroupName")
        family = self.parameters.get("DBParameterGroupFamily")
        desc = self.parameters.get("Description")
        db_cluster_parameter_group = self.backend.create_db_cluster_parameter_group(
            group_name=group_name,
            family=family,
            description=desc,
        )
        result = {"DBClusterParameterGroup": db_cluster_parameter_group}
        return ActionResult(result)

    def describe_db_cluster_parameter_groups(self) -> ActionResult:
        group_name = self.parameters.get("DBClusterParameterGroupName")
        db_parameter_groups = self.backend.describe_db_cluster_parameter_groups(
            group_name=group_name,
        )
        result = {"DBClusterParameterGroups": db_parameter_groups}
        return ActionResult(result)

    def delete_db_cluster_parameter_group(self) -> ActionResult:
        group_name = self.parameters.get("DBClusterParameterGroupName")
        self.backend.delete_db_cluster_parameter_group(
            group_name=group_name,
        )
        return ActionResult({})

    def promote_read_replica_db_cluster(self) -> ActionResult:
        db_cluster_identifier = self.parameters.get("DBClusterIdentifier")
        cluster = self.backend.promote_read_replica_db_cluster(db_cluster_identifier)
        result = {"DBCluster": cluster}
        return ActionResult(result)

    def describe_db_snapshot_attributes(self) -> ActionResult:
        params = self.parameters
        db_snapshot_identifier = params["DBSnapshotIdentifier"]
        db_snapshot_attributes_result = self.backend.describe_db_snapshot_attributes(
            db_snapshot_identifier=db_snapshot_identifier,
        )
        result = {
            "DBSnapshotAttributesResult": {
                "DBSnapshotIdentifier": db_snapshot_identifier,
                "DBSnapshotAttributes": [
                    {"AttributeName": name, "AttributeValues": values}
                    for name, values in db_snapshot_attributes_result.items()
                ],
            }
        }
        return ActionResult(result)

    def modify_db_snapshot_attribute(self) -> ActionResult:
        params = self.parameters
        db_snapshot_identifier = params["DBSnapshotIdentifier"]
        db_snapshot_attributes_result = self.backend.modify_db_snapshot_attribute(
            db_snapshot_identifier=db_snapshot_identifier,
            attribute_name=params["AttributeName"],
            values_to_add=params.get("ValuesToAdd"),
            values_to_remove=params.get("ValuesToRemove"),
        )
        result = {
            "DBSnapshotAttributesResult": {
                "DBSnapshotIdentifier": db_snapshot_identifier,
                "DBSnapshotAttributes": [
                    {"AttributeName": name, "AttributeValues": values}
                    for name, values in db_snapshot_attributes_result.items()
                ],
            }
        }
        return ActionResult(result)

    def describe_db_cluster_snapshot_attributes(self) -> ActionResult:
        params = self.parameters
        db_cluster_snapshot_identifier = params["DBClusterSnapshotIdentifier"]
        db_cluster_snapshot_attributes_result = (
            self.backend.describe_db_cluster_snapshot_attributes(
                db_cluster_snapshot_identifier=db_cluster_snapshot_identifier,
            )
        )
        result = {
            "DBClusterSnapshotAttributesResult": {
                "DBClusterSnapshotIdentifier": db_cluster_snapshot_identifier,
                "DBClusterSnapshotAttributes": [
                    {"AttributeName": name, "AttributeValues": values}
                    for name, values in db_cluster_snapshot_attributes_result.items()
                ],
            }
        }
        return ActionResult(result)

    def modify_db_cluster_snapshot_attribute(self) -> ActionResult:
        params = self.parameters
        db_cluster_snapshot_identifier = params["DBClusterSnapshotIdentifier"]
        db_cluster_snapshot_attributes_result = (
            self.backend.modify_db_cluster_snapshot_attribute(
                db_cluster_snapshot_identifier=db_cluster_snapshot_identifier,
                attribute_name=params["AttributeName"],
                values_to_add=params.get("ValuesToAdd"),
                values_to_remove=params.get("ValuesToRemove"),
            )
        )
        result = {
            "DBClusterSnapshotAttributesResult": {
                "DBClusterSnapshotIdentifier": db_cluster_snapshot_identifier,
                "DBClusterSnapshotAttributes": [
                    {"AttributeName": name, "AttributeValues": values}
                    for name, values in db_cluster_snapshot_attributes_result.items()
                ],
            }
        }
        return ActionResult(result)

    def describe_db_proxies(self) -> ActionResult:
        params = self.parameters
        db_proxy_name = params.get("DBProxyName")
        # filters = params.get("Filters")
        marker = params.get("Marker")
        db_proxies = self.backend.describe_db_proxies(
            db_proxy_name=db_proxy_name,
            # filters=filters,
        )
        result = {
            "DBProxies": db_proxies,
            "Marker": marker,
        }
        return ActionResult(result)

    def create_db_proxy(self) -> ActionResult:
        params = self.parameters
        db_proxy_name = params["DBProxyName"]
        engine_family = params["EngineFamily"]
        auth = params["Auth"]
        role_arn = params["RoleArn"]
        vpc_subnet_ids = params["VpcSubnetIds"]
        vpc_security_group_ids = params.get("VpcSecurityGroupIds")
        require_tls = params.get("RequireTLS")
        idle_client_timeout = params.get("IdleClientTimeout")
        debug_logging = params.get("DebugLogging")
        tags = self.parameters.get("Tags", [])
        db_proxy = self.backend.create_db_proxy(
            db_proxy_name=db_proxy_name,
            engine_family=engine_family,
            auth=auth,
            role_arn=role_arn,
            vpc_subnet_ids=vpc_subnet_ids,
            vpc_security_group_ids=vpc_security_group_ids,
            require_tls=require_tls,
            idle_client_timeout=idle_client_timeout,
            debug_logging=debug_logging,
            tags=tags,
        )
        result = {"DBProxy": db_proxy}
        return ActionResult(result)

    def register_db_proxy_targets(self) -> ActionResult:
        db_proxy_name = self.parameters.get("DBProxyName")
        target_group_name = self.parameters.get("TargetGroupName")
        db_cluster_identifiers = self.parameters.get("DBClusterIdentifiers", [])
        db_instance_identifiers = self.parameters.get("DBInstanceIdentifiers", [])
        targets = self.backend.register_db_proxy_targets(
            db_proxy_name=db_proxy_name,
            target_group_name=target_group_name,
            db_cluster_identifiers=db_cluster_identifiers,
            db_instance_identifiers=db_instance_identifiers,
        )
        result = {"DBProxyTargets": targets}
        return ActionResult(result)

    def deregister_db_proxy_targets(self) -> ActionResult:
        db_proxy_name = self.parameters.get("DBProxyName")
        target_group_name = self.parameters.get("TargetGroupName")
        db_cluster_identifiers = self.parameters.get("DBClusterIdentifiers", [])
        db_instance_identifiers = self.parameters.get("DBInstanceIdentifiers", [])
        self.backend.deregister_db_proxy_targets(
            db_proxy_name=db_proxy_name,
            target_group_name=target_group_name,
            db_cluster_identifiers=db_cluster_identifiers,
            db_instance_identifiers=db_instance_identifiers,
        )
        return ActionResult({})

    def describe_db_proxy_targets(self) -> ActionResult:
        proxy_name = self.parameters.get("DBProxyName")
        targets = self.backend.describe_db_proxy_targets(proxy_name=proxy_name)
        result = {"Targets": targets}
        return ActionResult(result)

    def delete_db_proxy(self) -> ActionResult:
        proxy_name = self.parameters.get("DBProxyName")
        proxy = self.backend.delete_db_proxy(proxy_name=proxy_name)
        result = {"DBProxy": proxy}
        return ActionResult(result)

    def describe_db_proxy_target_groups(self) -> ActionResult:
        proxy_name = self.parameters.get("DBProxyName")
        groups = self.backend.describe_db_proxy_target_groups(proxy_name=proxy_name)
        result = {"TargetGroups": groups}
        return ActionResult(result)

    def modify_db_proxy_target_group(self) -> ActionResult:
        proxy_name = self.parameters.get("DBProxyName")
        config = self.parameters.get("ConnectionPoolConfig", {})
        group = self.backend.modify_db_proxy_target_group(
            proxy_name=proxy_name, config=config
        )
        result = {"DBProxyTargetGroup": group}
        return ActionResult(result)

    def describe_db_instance_automated_backups(self) -> ActionResult:
        automated_backups = self.backend.describe_db_instance_automated_backups(
            **self.parameters
        )
        result = {"DBInstanceAutomatedBackups": automated_backups}
        return ActionResult(result)

    def describe_events(self) -> ActionResult:
        events = self.backend.describe_events(**self.parameters)
        result = {"Events": events}
        return ActionResult(result)

    def describe_db_log_files(self) -> ActionResult:
        log_files = self.backend.describe_db_log_files(**self.parameters)
        result = {"DescribeDBLogFiles": log_files}
        return ActionResult(result)

    def _paginate(self, resources: List[Any]) -> Tuple[List[Any], Optional[str]]:
        from moto.rds.exceptions import InvalidParameterValue

        marker = self.parameters.get("Marker")
        # Default was originally set to 50 instead of 100 for ease of testing.  Should fix.
        page_size = self.parameters.get("MaxRecords", 50)
        if page_size < 20 or page_size > 100:
            msg = (
                f"Invalid value {page_size} for MaxRecords. Must be between 20 and 100"
            )
            raise InvalidParameterValue(msg)
        all_resources = list(resources)
        all_ids = [resource.resource_id for resource in all_resources]
        if marker:
            start = all_ids.index(marker) + 1
        else:
            start = 0
        paginated_resources = all_resources[start : start + page_size]
        next_marker = None
        if len(all_resources) > start + page_size:
            next_marker = paginated_resources[-1].resource_id
        return paginated_resources, next_marker
