from typing import Any, Dict, Iterator, List, Optional, Tuple

from moto.acm.models import AWSCertificateManagerBackend, acm_backends
from moto.appsync.models import AppSyncBackend, appsync_backends
from moto.athena.models import athena_backends
from moto.awslambda.models import LambdaBackend, lambda_backends
from moto.backup.models import BackupBackend, backup_backends
from moto.clouddirectory import CloudDirectoryBackend, clouddirectory_backends
from moto.cloudfront.models import CloudFrontBackend, cloudfront_backends
from moto.cloudwatch.models import CloudWatchBackend, cloudwatch_backends
from moto.comprehend.models import ComprehendBackend, comprehend_backends
from moto.connectcampaigns.models import (
    ConnectCampaignServiceBackend,
    connectcampaigns_backends,
)
from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.exceptions import RESTError
from moto.directconnect.models import DirectConnectBackend, directconnect_backends
from moto.dms.models import DatabaseMigrationServiceBackend, dms_backends
from moto.dynamodb.models import DynamoDBBackend, dynamodb_backends
from moto.ec2 import ec2_backends
from moto.ecs.models import EC2ContainerServiceBackend, ecs_backends
from moto.efs.models import EFSBackend, efs_backends
from moto.elb.models import ELBBackend, elb_backends
from moto.elbv2.models import ELBv2Backend, elbv2_backends
from moto.emr.models import ElasticMapReduceBackend, emr_backends
from moto.events.models import EventsBackend, events_backends
from moto.glacier.models import GlacierBackend, glacier_backends
from moto.glue.models import GlueBackend, glue_backends
from moto.kafka.models import KafkaBackend, kafka_backends
from moto.kinesis.models import KinesisBackend, kinesis_backends
from moto.kinesisanalyticsv2.models import (
    KinesisAnalyticsV2Backend,
    kinesisanalyticsv2_backends,
)
from moto.kms.models import KmsBackend, kms_backends
from moto.lexv2models.models import LexModelsV2Backend, lexv2models_backends
from moto.logs.models import LogsBackend, logs_backends
from moto.moto_api._internal import mock_random
from moto.quicksight.models import QuickSightBackend, quicksight_backends
from moto.rds.models import RDSBackend, rds_backends
from moto.redshift.models import RedshiftBackend, redshift_backends
from moto.s3.models import S3Backend, s3_backends
from moto.sagemaker.models import SageMakerModelBackend, sagemaker_backends
from moto.secretsmanager import secretsmanager_backends
from moto.secretsmanager.models import ReplicaSecret, SecretsManagerBackend
from moto.sns.models import SNSBackend, sns_backends
from moto.sqs.models import SQSBackend, sqs_backends
from moto.ssm.models import SimpleSystemManagerBackend, ssm_backends
from moto.stepfunctions.models import StepFunctionBackend, stepfunctions_backends
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import get_partition
from moto.workspaces.models import WorkSpacesBackend, workspaces_backends
from moto.workspacesweb.models import WorkSpacesWebBackend, workspacesweb_backends

# Left: EC2 ElastiCache RDS ELB Lambda EMR Glacier Kinesis Redshift Route53
# StorageGateway DynamoDB MachineLearning ACM DirectConnect DirectoryService CloudHSM
# Inspector Elasticsearch


class ResourceGroupsTaggingAPIBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)

        self._pages: Dict[str, Any] = {}
        # Like 'someuuid': {'gen': <generator>, 'misc': None}
        # Misc is there for peeking from a generator and it cant
        # fit in the current request. As we only store generators
        # there is really no point cleaning up

    @property
    def appsync_backend(self) -> AppSyncBackend:
        return appsync_backends[self.account_id][self.region_name]

    @property
    def s3_backend(self) -> S3Backend:
        return s3_backends[self.account_id][self.partition]

    @property
    def directconnect_backend(self) -> DirectConnectBackend:
        return directconnect_backends[self.account_id][self.region_name]

    @property
    def dms_backend(self) -> DatabaseMigrationServiceBackend:
        return dms_backends[self.account_id][self.region_name]

    @property
    def ec2_backend(self) -> Any:  # type: ignore[misc]
        return ec2_backends[self.account_id][self.region_name]

    @property
    def efs_backend(self) -> EFSBackend:
        return efs_backends[self.account_id][self.region_name]

    @property
    def elb_backend(self) -> ELBBackend:
        return elb_backends[self.account_id][self.region_name]

    @property
    def elbv2_backend(self) -> ELBv2Backend:
        return elbv2_backends[self.account_id][self.region_name]

    @property
    def events_backend(self) -> EventsBackend:
        return events_backends[self.account_id][self.region_name]

    @property
    def glue_backend(self) -> GlueBackend:
        return glue_backends[self.account_id][self.region_name]

    @property
    def kinesis_backend(self) -> KinesisBackend:
        return kinesis_backends[self.account_id][self.region_name]

    @property
    def kinesisanalyticsv2_backend(self) -> KinesisAnalyticsV2Backend:
        return kinesisanalyticsv2_backends[self.account_id][self.region_name]

    @property
    def kms_backend(self) -> KmsBackend:
        return kms_backends[self.account_id][self.region_name]

    @property
    def logs_backend(self) -> LogsBackend:
        return logs_backends[self.account_id][self.region_name]

    @property
    def rds_backend(self) -> RDSBackend:
        return rds_backends[self.account_id][self.region_name]

    @property
    def glacier_backend(self) -> GlacierBackend:
        return glacier_backends[self.account_id][self.region_name]

    @property
    def emr_backend(self) -> ElasticMapReduceBackend:
        return emr_backends[self.account_id][self.region_name]

    @property
    def redshift_backend(self) -> RedshiftBackend:
        return redshift_backends[self.account_id][self.region_name]

    @property
    def lambda_backend(self) -> LambdaBackend:
        return lambda_backends[self.account_id][self.region_name]

    @property
    def ecs_backend(self) -> EC2ContainerServiceBackend:
        return ecs_backends[self.account_id][self.region_name]

    @property
    def acm_backend(self) -> AWSCertificateManagerBackend:
        return acm_backends[self.account_id][self.region_name]

    @property
    def secretsmanager_backend(self) -> SecretsManagerBackend:
        return secretsmanager_backends[self.account_id][self.region_name]

    @property
    def sns_backend(self) -> SNSBackend:
        return sns_backends[self.account_id][self.region_name]

    @property
    def ssm_backend(self) -> SimpleSystemManagerBackend:
        return ssm_backends[self.account_id][self.region_name]

    @property
    def sqs_backend(self) -> SQSBackend:
        return sqs_backends[self.account_id][self.region_name]

    @property
    def stepfunctions_backend(self) -> StepFunctionBackend:
        return stepfunctions_backends[self.account_id][self.region_name]

    @property
    def backup_backend(self) -> BackupBackend:
        return backup_backends[self.account_id][self.region_name]

    @property
    def dynamodb_backend(self) -> DynamoDBBackend:
        return dynamodb_backends[self.account_id][self.region_name]

    @property
    def workspaces_backend(self) -> Optional[WorkSpacesBackend]:
        # Workspaces service has limited region availability
        if self.region_name in workspaces_backends[self.account_id].regions:
            return workspaces_backends[self.account_id][self.region_name]
        return None

    @property
    def workspacesweb_backends(self) -> Optional[WorkSpacesWebBackend]:
        # Workspaces service has limited region availability
        if self.region_name in workspaces_backends[self.account_id].regions:
            return workspacesweb_backends[self.account_id][self.region_name]
        return None

    @property
    def comprehend_backend(self) -> Optional[ComprehendBackend]:
        # aws Comprehend has limited region availability
        if self.region_name in comprehend_backends[self.account_id].regions:
            return comprehend_backends[self.account_id][self.region_name]
        return None

    @property
    def kafka_backend(self) -> KafkaBackend:
        return kafka_backends[self.account_id][self.region_name]

    @property
    def sagemaker_backend(self) -> SageMakerModelBackend:
        return sagemaker_backends[self.account_id][self.region_name]

    @property
    def lexv2_backend(self) -> Optional[LexModelsV2Backend]:
        if self.region_name in lexv2models_backends[self.account_id].regions:
            return lexv2models_backends[self.account_id][self.region_name]
        return None

    @property
    def clouddirectory_backend(self) -> Optional[CloudDirectoryBackend]:
        if self.region_name in clouddirectory_backends[self.account_id].regions:
            return clouddirectory_backends[self.account_id][self.region_name]
        return None

    @property
    def cloudfront_backend(self) -> CloudFrontBackend:
        return cloudfront_backends[self.account_id][self.partition]

    @property
    def cloudwatch_backend(self) -> CloudWatchBackend:
        return cloudwatch_backends[self.account_id][self.region_name]

    @property
    def connectcampaigns_backend(self) -> Optional[ConnectCampaignServiceBackend]:
        # Connect Campaigns service has limited region availability
        if self.region_name in connectcampaigns_backends[self.account_id].regions:
            return connectcampaigns_backends[self.account_id][self.region_name]
        return None

    @property
    def quicksight_backend(self) -> Optional[QuickSightBackend]:
        if self.region_name in quicksight_backends[self.account_id].regions:
            return quicksight_backends[self.account_id][self.region_name]
        return None

    def _get_resources_generator(
        self,
        tag_filters: Optional[List[Dict[str, Any]]] = None,
        resource_type_filters: Optional[List[str]] = None,
    ) -> Iterator[Dict[str, Any]]:
        # Look at
        # https://docs.aws.amazon.com/general/latest/gr/aws-arns-and-namespaces.html

        # TODO move these to their respective backends
        filters = []
        for tag_filter_dict in tag_filters:  # type: ignore
            values = tag_filter_dict.get("Values", [])
            if len(values) == 0:
                # Check key matches
                filters.append(lambda t, v, key=tag_filter_dict["Key"]: t == key)
            elif len(values) == 1:
                # Check it's exactly the same as key, value
                filters.append(
                    lambda t, v, key=tag_filter_dict["Key"], value=values[0]: t == key  # type: ignore
                    and v == value
                )
            else:
                # Check key matches and value is one of the provided values
                filters.append(
                    lambda t, v, key=tag_filter_dict["Key"], vl=values: t == key  # type: ignore
                    and v in vl
                )

        def tag_filter(tag_list: List[Dict[str, Any]]) -> bool:
            result = []
            if tag_filters:
                for f in filters:
                    temp_result = []
                    for tag in tag_list:
                        f_result = f(tag["Key"], tag["Value"])  # type: ignore
                        temp_result.append(f_result)
                    result.append(any(temp_result))
                return all(result)
            else:
                return True

        def format_tags(tags: Dict[str, Any]) -> List[Dict[str, Any]]:
            result = []
            for key, value in tags.items():
                result.append({"Key": key, "Value": value})
            return result

        def format_tag_keys(
            tags: List[Dict[str, Any]], keys: List[str]
        ) -> List[Dict[str, Any]]:
            result = []
            for tag in tags:
                result.append({"Key": tag[keys[0]], "Value": tag[keys[1]]})
            return result

        # ACM
        if not resource_type_filters or "acm" in resource_type_filters:
            for certificate in self.acm_backend._certificates.values():
                tags = format_tags(certificate.tags)
                if not tags or not tag_filter(tags):
                    continue
                yield {"ResourceARN": f"{certificate.arn}", "Tags": tags}

        # AppSync
        if not resource_type_filters or "appsync" in resource_type_filters:
            for api in self.appsync_backend.graphql_apis.values():
                tags = self.appsync_backend.tagger.list_tags_for_resource(api.arn)[
                    "Tags"
                ]
                if not tags or not tag_filter(
                    tags
                ):  # Skip if no tags, or invalid filter
                    continue

                yield {"ResourceARN": f"{api.arn}", "Tags": tags}

        # Athena
        if not resource_type_filters or "athena" in resource_type_filters:
            athena_backend = athena_backends[self.account_id][self.region_name]

            # Capacity Reservations
            for capacity_reservation in athena_backend.capacity_reservations.values():
                tags = athena_backend.tagger.list_tags_for_resource(
                    capacity_reservation.name
                )["Tags"]
                if not tags or not tag_filter(tags):
                    continue
                yield {"ResourceARN": f"{capacity_reservation.name}", "Tags": tags}

            # Workgroups
            for work_group in athena_backend.work_groups.values():
                tags = athena_backend.tagger.list_tags_for_resource(work_group.name)[
                    "Tags"
                ]
                if not tags or not tag_filter(tags):
                    continue
                yield {"ResourceARN": f"{work_group.name}", "Tags": tags}

            # Data Catalogs
            for data_catalog in athena_backend.data_catalogs.values():
                tags = athena_backend.tagger.list_tags_for_resource(data_catalog.name)[
                    "Tags"
                ]
                if not tags or not tag_filter(tags):
                    continue
                yield {"ResourceARN": f"{data_catalog.name}", "Tags": tags}

        # Backup
        if not resource_type_filters or "backup" in resource_type_filters:
            for vault in self.backup_backend.vaults.values():
                tags = self.backup_backend.tagger.list_tags_for_resource(
                    vault.backup_vault_arn
                )["Tags"]
                if not tags or not tag_filter(
                    tags
                ):  # Skip if no tags, or invalid filter
                    continue

                yield {"ResourceARN": f"{vault.backup_vault_arn}", "Tags": tags}

        # Comprehend
        if self.comprehend_backend:
            if not resource_type_filters or "comprehend" in resource_type_filters:
                for document_classifier in self.comprehend_backend.classifiers.values():
                    tags = self.comprehend_backend.tagger.list_tags_for_resource(
                        document_classifier.arn
                    )["Tags"]
                    if not tags or not tag_filter(tags):
                        continue
                    yield {
                        "ResourceARN": f"{document_classifier.arn}",
                        "Tags": tags,
                    }

                for entity_recognizer in self.comprehend_backend.recognizers.values():
                    tags = self.comprehend_backend.tagger.list_tags_for_resource(
                        entity_recognizer.arn
                    )["Tags"]
                    if not tags or not tag_filter(tags):
                        continue
                    yield {
                        "ResourceARN": f"{entity_recognizer.arn}",
                        "Tags": tags,
                    }

        # S3
        if (
            not resource_type_filters
            or "s3" in resource_type_filters
            or "s3:bucket" in resource_type_filters
        ):
            for bucket in self.s3_backend.buckets.values():
                tags = self.s3_backend.tagger.list_tags_for_resource(bucket.arn)["Tags"]
                if not tags or not tag_filter(
                    tags
                ):  # Skip if no tags, or invalid filter
                    continue
                yield {"ResourceARN": bucket.arn, "Tags": tags}

        # Cloud Directory
        if self.clouddirectory_backend:
            if not resource_type_filters or "clouddirectory" in resource_type_filters:
                clouddirectory_backend = clouddirectory_backends[self.account_id][
                    self.region_name
                ]
                for directory in clouddirectory_backend.directories.values():
                    tags = clouddirectory_backend.tagger.list_tags_for_resource(
                        directory.directory_arn
                    )["Tags"]
                    if not tags or not tag_filter(tags):
                        continue
                    yield {"ResourceARN": f"{directory.directory_arn}", "Tags": tags}

        # CloudFormation
        if not resource_type_filters or "cloudformation:stack" in resource_type_filters:
            try:
                from moto.cloudformation import cloudformation_backends

                backend = cloudformation_backends[self.account_id][self.region_name]

                for stack in backend.stacks.values():
                    tags = format_tags(stack.tags)
                    if not tag_filter(tags):
                        continue
                    yield {"ResourceARN": f"{stack.stack_id}", "Tags": tags}

            except Exception:
                pass

        # Cloudfront
        if (
            not resource_type_filters
            or "cloudfront" in resource_type_filters
            or "cloudfront:distributions" in resource_type_filters
        ):
            for dist in self.cloudfront_backend.distributions.values():
                tags = self.cloudfront_backend.tagger.list_tags_for_resource(dist.arn)[
                    "Tags"
                ]
                if (
                    not tag_filter(tags) or len(tags) == 0
                ):  # Skip if no tags, or invalid filter
                    continue
                yield {"ResourceARN": f"{dist.arn}", "Tags": tags}

        if self.cloudwatch_backend:
            cloudwatch_resource_map: Dict[str, Dict[str, Any]] = {
                "cloudwatch:alarm": self.cloudwatch_backend.alarms,
                "cloudwatch:insight-rule": self.cloudwatch_backend.insight_rules,
            }
            for resource_type, resource_source in cloudwatch_resource_map.items():
                if (
                    not resource_type_filters
                    or "cloudwatch" in resource_type_filters
                    or resource_type in resource_type_filters
                ):
                    for resource in resource_source.values():
                        if resource_type == "cloudwatch:alarm":
                            end_of_arn = f"alarm:{resource.name}"
                        elif resource_type == "cloudwatch:insight-rule":
                            end_of_arn = f"insight-rule/{resource.name}"

                        arn = f"arn:{get_partition(self.region_name)}:cloudwatch:{self.region_name}:{self.account_id}:{end_of_arn}"

                        cw_tags = self.cloudwatch_backend.list_tags_for_resource(arn)

                        tags = format_tags(cw_tags)
                        if not tags or not tag_filter(tags):
                            continue
                        yield {
                            "ResourceARN": arn,
                            "Tags": tags,
                        }

        # Connect Campaigns v1
        if self.connectcampaigns_backend:
            if not resource_type_filters or "connectcampaigns" in resource_type_filters:
                connectcampaigns_backend: ConnectCampaignServiceBackend = (
                    connectcampaigns_backends[self.account_id][self.region_name]
                )
                for campaign in connectcampaigns_backend.campaigns.values():
                    tags = connectcampaigns_backend.tagger.list_tags_for_resource(
                        campaign.arn
                    )["Tags"]
                    if not tags or not tag_filter(tags):
                        continue
                    yield {"ResourceARN": f"{campaign.arn}", "Tags": tags}

        # Direct Connect
        if self.directconnect_backend:
            if not resource_type_filters or "directconnect" in resource_type_filters:
                directconnect_backend = directconnect_backends[self.account_id][
                    self.region_name
                ]

                # Connections
                for connection in directconnect_backend.connections.values():
                    tags = directconnect_backend.tagger.list_tags_for_resource(
                        connection.connection_id
                    )["Tags"]
                    tags = format_tag_keys(tags, ["key", "value"])
                    if not tags or not tag_filter(tags):
                        continue
                    yield {"ResourceARN": f"{connection.connection_id}", "Tags": tags}

                # LAGs
                for lag in directconnect_backend.lags.values():
                    tags = directconnect_backend.tagger.list_tags_for_resource(
                        lag.lag_id
                    )["Tags"]
                    tags = format_tag_keys(tags, ["key", "value"])
                    if not tags or not tag_filter(tags):
                        continue
                    yield {"ResourceARN": f"{lag.lag_id}", "Tags": tags}

        # DMS
        if not resource_type_filters or "dms:endpoint" in resource_type_filters:
            for endpoint in self.dms_backend.endpoints.values():
                tags = self.dms_backend.tagger.list_tags_for_resource(
                    endpoint.endpoint_arn
                )["Tags"]
                if not tag_filter(tags):
                    continue
                yield {"ResourceARN": f"{endpoint.endpoint_arn}", "Tags": tags}

        if (
            not resource_type_filters
            or "dms:replication-instance" in resource_type_filters
        ):
            for replication_instance in self.dms_backend.replication_instances.values():
                tags = self.dms_backend.tagger.list_tags_for_resource(
                    replication_instance.arn
                )["Tags"]
                if not tag_filter(tags):
                    continue
                yield {"ResourceARN": f"{replication_instance.arn}", "Tags": tags}

        # ECS
        if not resource_type_filters or "ecs:service" in resource_type_filters:
            for service in self.ecs_backend.services.values():
                tags = format_tag_keys(service.tags, ["key", "value"])
                if not tag_filter(tags):
                    continue
                yield {"ResourceARN": f"{service.physical_resource_id}", "Tags": tags}

        if not resource_type_filters or "ecs:cluster" in resource_type_filters:
            for cluster in self.ecs_backend.clusters.values():
                tags = format_tag_keys(cluster.tags, ["key", "value"])  # type: ignore[arg-type]
                if not tag_filter(tags):
                    continue
                yield {"ResourceARN": f"{cluster.arn}", "Tags": tags}

        if not resource_type_filters or "ecs:task" in resource_type_filters:
            for task_dict in self.ecs_backend.tasks.values():
                for task in task_dict.values():
                    tags = format_tag_keys(task.tags, ["key", "value"])
                    if not tag_filter(tags):
                        continue
                    yield {"ResourceARN": f"{task.task_arn}", "Tags": tags}

        if not resource_type_filters or "ecs:task-definition" in resource_type_filters:
            for task_definition_dict in self.ecs_backend.task_definitions.values():
                for task_definition in task_definition_dict.values():
                    tags = format_tag_keys(task_definition.tags, ["key", "value"])
                    if not tag_filter(tags):
                        continue
                    yield {"ResourceARN": f"{task_definition.arn}", "Tags": tags}

        # EC2 Resources
        ec2_resource_types = {
            "ec2:image": self.ec2_backend.amis.values(),
            "ec2:instance": (
                instance
                for reservation in self.ec2_backend.reservations.values()
                for instance in reservation.instances
            ),
            "ec2:network-interface": self.ec2_backend.enis.values(),
            "ec2:security-group": (
                sg for vpc in self.ec2_backend.groups.values() for sg in vpc.values()
            ),
            "ec2:snapshot": self.ec2_backend.snapshots.values(),
            "ec2:volume": self.ec2_backend.volumes.values(),
            "ec2:vpc": self.ec2_backend.vpcs.values(),
            "ec2:subnet": (
                subnet
                for subnet in self.ec2_backend.subnets.values()
                for subnet in subnet.values()
            ),
            "ec2:vpc-peering-connection": self.ec2_backend.vpc_pcxs.values(),
            "ec2:transit-gateway": self.ec2_backend.transit_gateways.values(),
            "ec2:transit-gateway-attachment": self.ec2_backend.transit_gateway_attachments.values(),
            "ec2:route-table": self.ec2_backend.route_tables.values(),
            "ec2:customer-gateway": self.ec2_backend.customer_gateways.values(),
            "ec2:vpn-connection": self.ec2_backend.vpn_connections.values(),
            "ec2:natgateway": self.ec2_backend.nat_gateways.values(),
            "ec2:internet-gateway": self.ec2_backend.internet_gateways.values(),
            "ec2:managed-prefix-lists": self.ec2_backend.managed_prefix_lists.values(),
            "ec2:flow-logs": self.ec2_backend.flow_logs.values(),
            "ec2:spot-instance-request": self.ec2_backend.spot_instance_requests.values(),
            # TODO: "ec2:reserved-instance": ...,
        }

        for resource_type, resources in ec2_resource_types.items():
            if (
                not resource_type_filters
                or "ec2" in resource_type_filters
                or resource_type in resource_type_filters
            ):
                for resource in resources:
                    tags = format_tags(self.ec2_backend.tags.get(resource.id, {}))
                    if not tags or not tag_filter(tags):
                        continue
                    yield {
                        "ResourceARN": f"arn:{self.partition}:ec2:{self.region_name}:{self.account_id}:{resource_type}/{resource.id}",
                        "Tags": tags,
                    }

        # EFS, resource type elasticfilesystem:access-point
        if (
            not resource_type_filters
            or "elasticfilesystem" in resource_type_filters
            or "elasticfilesystem:access-point" in resource_type_filters
        ):
            for ap in self.efs_backend.access_points.values():
                tags = self.efs_backend.list_tags_for_resource(ap.access_point_id)
                if not tag_filter(tags):
                    continue
                yield {"ResourceARN": f"{ap.access_point_arn}", "Tags": tags}

        # EFS, resource type elasticfilesystem:file-system
        if (
            not resource_type_filters
            or "elasticfilesystem" in resource_type_filters
            or "elasticfilesystem:file-system" in resource_type_filters
        ):
            for fs in self.efs_backend.file_systems_by_id.values():
                tags = self.efs_backend.list_tags_for_resource(fs.file_system_id)
                if not tag_filter(tags):
                    continue
                yield {"ResourceARN": f"{fs.file_system_arn}", "Tags": tags}

        # ELB (Classic Load Balancers)
        if (
            not resource_type_filters
            or "elb" in resource_type_filters
            or "elb:loadbalancer" in resource_type_filters
        ):
            for elb in self.elb_backend.load_balancers.values():
                tags = format_tags(elb.tags)
                if not tags or not tag_filter(
                    tags
                ):  # Skip if no tags, or invalid filter
                    continue

                yield {
                    "ResourceARN": f"arn:{get_partition(self.region_name)}:elasticloadbalancing:{self.region_name}:{self.account_id}:loadbalancer/{elb.name}",
                    "Tags": tags,
                }

        # TODO add these to the keys and values functions / combine functions
        # ELB, resource type elasticloadbalancing:loadbalancer
        if (
            not resource_type_filters
            or "elasticloadbalancing" in resource_type_filters
            or "elasticloadbalancing:loadbalancer" in resource_type_filters
        ):
            for elbv2 in self.elbv2_backend.load_balancers.values():
                tags = self.elbv2_backend.tagging_service.list_tags_for_resource(
                    elbv2.arn
                )["Tags"]
                if not tag_filter(tags):  # Skip if no tags, or invalid filter
                    continue

                yield {"ResourceARN": f"{elbv2.arn}", "Tags": tags}

        # ELB Target Group, resource type elasticloadbalancing:targetgroup
        if (
            not resource_type_filters
            or "elasticloadbalancing" in resource_type_filters
            or "elasticloadbalancing:targetgroup" in resource_type_filters
        ):
            for target_group in self.elbv2_backend.target_groups.values():
                tags = self.elbv2_backend.tagging_service.list_tags_for_resource(
                    target_group.arn
                )["Tags"]
                if not tag_filter(tags):  # Skip if no tags, or invalid filter
                    continue

                yield {"ResourceARN": f"{target_group.arn}", "Tags": tags}

        # EMR Cluster

        # Events
        if (
            not resource_type_filters
            or "events" in resource_type_filters
            or "events:event-bus" in resource_type_filters
        ):
            for bus in self.events_backend.event_buses.values():
                tags = self.events_backend.tagger.list_tags_for_resource(bus.arn)[
                    "Tags"
                ]
                if (
                    not tag_filter(tags) or len(tags) == 0
                ):  # Skip if no tags, or invalid filter
                    continue
                yield {"ResourceARN": f"{bus.arn}", "Tags": tags}

        # Glacier Vault

        # Glue
        if not resource_type_filters or any(
            ("glue" in _type) for _type in resource_type_filters
        ):
            if not resource_type_filters or "glue" in resource_type_filters:
                arns_starting_with = [
                    f"arn:{get_partition(self.region_name)}:glue:{self.region_name}:{self.account_id}:"
                ]
            else:
                arns_starting_with = []
                for resource_type in resource_type_filters:
                    if resource_type.startswith("glue:"):
                        glue_type = resource_type.split(":")[-1]
                        arns_starting_with.append(
                            f"arn:{get_partition(self.region_name)}:glue:{self.region_name}:{self.account_id}:{glue_type}"
                        )
            for glue_arn in self.glue_backend.tagger.tags.keys():
                if any(glue_arn.startswith(arn) for arn in arns_starting_with):
                    tags = self.glue_backend.tagger.list_tags_for_resource(glue_arn)[
                        "Tags"
                    ]
                    yield {"ResourceARN": glue_arn, "Tags": tags}

        # Kinesis

        # KinesisAnalyticsV2
        if self.kinesisanalyticsv2_backend and (
            not resource_type_filters or "kinesisanalyticsv2" in resource_type_filters
        ):
            for application in self.kinesisanalyticsv2_backend.applications.values():
                tags = self.kinesisanalyticsv2_backend.tagger.list_tags_for_resource(
                    application.application_arn
                )["Tags"]
                if not tags or not tag_filter(tags):
                    continue
                yield {
                    "ResourceARN": application.application_arn,
                    "Tags": tags,
                }

        # KMS
        if not resource_type_filters or "kms" in resource_type_filters:
            for kms_key in self.kms_backend.list_keys():
                tags = format_tag_keys(
                    self.kms_backend.list_resource_tags(kms_key.id).get("Tags", []),
                    ["TagKey", "TagValue"],
                )
                if not tag_filter(tags):  # Skip if no tags, or invalid filter
                    continue

                yield {"ResourceARN": f"{kms_key.arn}", "Tags": tags}

        # LexV2
        if self.lexv2_backend:
            lex_v2_resource_map: Dict[str, Dict[str, Any]] = {
                "lexv2:bot": self.lexv2_backend.bots,
                "lexv2:bot-alias": self.lexv2_backend.bot_aliases,
            }
            for resource_type, resource_source in lex_v2_resource_map.items():
                if (
                    not resource_type_filters
                    or "lexv2" in resource_type_filters
                    or resource_type in resource_type_filters
                ):
                    for resource in resource_source.values():
                        bot_tags = self.lexv2_backend.list_tags_for_resource(
                            resource.arn
                        )

                        tags = format_tags(bot_tags)
                        if not tags or not tag_filter(tags):
                            continue
                        yield {
                            "ResourceARN": resource.arn,
                            "Tags": tags,
                        }

        # LOGS
        if (
            not resource_type_filters
            or "logs" in resource_type_filters
            or "logs:loggroup" in resource_type_filters
        ):
            for group in self.logs_backend.groups.values():
                log_tags = self.logs_backend.list_tags_for_resource(group.arn)
                tags = format_tags(log_tags)

                if not log_tags or not tag_filter(tags):
                    # Skip if no tags, or invalid filter
                    continue
                yield {"ResourceARN": group.arn, "Tags": tags}

        # Quicksight
        if self.quicksight_backend:
            quicksight_resource_map: dict[str, dict[str, Any]] = {
                "quicksight:dashboards": dict(self.quicksight_backend.dashboards),
                "quicksight:data_sources": dict(self.quicksight_backend.data_sources),
                "quicksight:data_sets": dict(self.quicksight_backend.data_sets),
                "quicksight:users": dict(self.quicksight_backend.users),
            }

            for resource_type, resource_source in quicksight_resource_map.items():
                if (
                    not resource_type_filters
                    or "quicksight" in resource_type_filters
                    or resource_type in resource_type_filters
                ):
                    for resource in resource_source.values():
                        tags = self.quicksight_backend.tagger.list_tags_for_resource(
                            resource.arn
                        )["Tags"]

                        if not tags or not tag_filter(tags):
                            continue

                        yield {
                            "ResourceARN": resource.arn,
                            "Tags": tags,
                        }

        # RDS resources
        resource_map: dict[str, dict[str, Any]] = {
            "rds:cluster": dict(self.rds_backend.clusters),
            "rds:db": dict(self.rds_backend.databases),
            "rds:snapshot": dict(self.rds_backend.database_snapshots),
            "rds:cluster-snapshot": dict(self.rds_backend.cluster_snapshots),
            "rds:db-proxy": self.rds_backend.db_proxies,
        }
        for resource_type, resource_source in resource_map.items():
            if (
                not resource_type_filters
                or "rds" in resource_type_filters
                or resource_type in resource_type_filters
            ):
                for resource in resource_source.values():
                    tags = resource.get_tags()
                    if not tags or not tag_filter(tags):
                        continue
                    yield {
                        "ResourceARN": resource.arn,
                        "Tags": tags,
                    }

        # RDS Reserved Database Instance
        # RDS Option Group
        # RDS Parameter Group
        # RDS Security Group
        # RDS Subnet Group
        # RDS Event Subscription

        # RedShift Cluster
        # RedShift Hardware security module (HSM) client certificate
        # RedShift HSM connection
        # RedShift Parameter group
        # RedShift Snapshot
        # RedShift Subnet group

        # Secrets Manager
        if (
            not resource_type_filters
            or "secretsmanager" in resource_type_filters
            or "secretsmanager:secret" in resource_type_filters
        ):
            for secret in self.secretsmanager_backend.secrets.values():
                if isinstance(secret, ReplicaSecret):
                    secret_tags = secret.source.tags
                else:
                    secret_tags = secret.tags

                if secret_tags:
                    formated_tags = format_tag_keys(secret_tags, ["Key", "Value"])
                    if not formated_tags or not tag_filter(formated_tags):
                        continue
                    yield {"ResourceARN": f"{secret.arn}", "Tags": formated_tags}

        # SQS
        if not resource_type_filters or "sqs" in resource_type_filters:
            for queue in self.sqs_backend.queues.values():
                tags = format_tags(queue.tags)
                if not tags or not tag_filter(
                    tags
                ):  # Skip if no tags, or invalid filter
                    continue

                yield {"ResourceARN": f"{queue.queue_arn}", "Tags": tags}

        # SNS
        if not resource_type_filters or "sns" in resource_type_filters:
            for topic in self.sns_backend.topics.values():
                tags = format_tags(topic._tags)
                if not tags or not tag_filter(
                    tags
                ):  # Skip if no tags, or invalid filter
                    continue
                yield {"ResourceARN": f"{topic.arn}", "Tags": tags}

        # SSM
        if not resource_type_filters or "ssm" in resource_type_filters:
            for document in self.ssm_backend._documents.values():
                doc_name = document.describe()["Name"]
                tags = self.ssm_backend._get_documents_tags(doc_name)
                if not tags or not tag_filter(
                    tags
                ):  # Skip if no tags, or invalid filter
                    continue
                yield {
                    "ResourceARN": f"arn:{get_partition(self.region_name)}:ssm:{self.region_name}:{self.account_id}:document/{doc_name}",
                    "Tags": tags,
                }
        # Step Functions
        if not resource_type_filters or "states:stateMachine" in resource_type_filters:
            for state_machine in self.stepfunctions_backend.state_machines:
                tags = format_tag_keys(
                    state_machine.backend.get_tags_list_for_state_machine(
                        state_machine.arn
                    ),
                    [
                        state_machine.backend.tagger.key_name,
                        state_machine.backend.tagger.value_name,
                    ],
                )
                if not tags or not tag_filter(tags):
                    continue
                yield {"ResourceARN": state_machine.arn, "Tags": tags}

        # Workspaces
        if self.workspaces_backend and (
            not resource_type_filters or "workspaces" in resource_type_filters
        ):
            for ws in self.workspaces_backend.workspaces.values():
                tags = format_tag_keys(ws.tags, ["Key", "Value"])
                if not tags or not tag_filter(
                    tags
                ):  # Skip if no tags, or invalid filter
                    continue

                yield {
                    "ResourceARN": f"arn:{get_partition(self.region_name)}:workspaces:{self.region_name}:{self.account_id}:workspace/{ws.workspace_id}",
                    "Tags": tags,
                }

        # Workspace Directories
        if self.workspaces_backend and (
            not resource_type_filters or "workspaces-directory" in resource_type_filters
        ):
            for wd in self.workspaces_backend.workspace_directories.values():
                tags = format_tag_keys(wd.tags, ["Key", "Value"])
                if not tags or not tag_filter(
                    tags
                ):  # Skip if no tags, or invalid filter
                    continue

                yield {
                    "ResourceARN": f"arn:{get_partition(self.region_name)}:workspaces:{self.region_name}:{self.account_id}:directory/{wd.directory_id}",
                    "Tags": tags,
                }

        # Workspace Images
        if self.workspaces_backend and (
            not resource_type_filters or "workspaces-image" in resource_type_filters
        ):
            for wi in self.workspaces_backend.workspace_images.values():
                tags = format_tag_keys(wi.tags, ["Key", "Value"])
                if not tags or not tag_filter(
                    tags
                ):  # Skip if no tags, or invalid filter
                    continue

                yield {
                    "ResourceARN": f"arn:{get_partition(self.region_name)}:workspaces:{self.region_name}:{self.account_id}:workspaceimage/{wi.image_id}",
                    "Tags": tags,
                }

        # Kafka (MSK)
        if self.kafka_backend and (
            not resource_type_filters or "kafka" in resource_type_filters
        ):
            for msk_cluster in self.kafka_backend.clusters.values():
                tag_dict = self.kafka_backend.list_tags_for_resource(msk_cluster.arn)
                tags = [{"Key": key, "Value": value} for key, value in tag_dict.items()]

                if not tags or not tag_filter(tags):
                    continue

                yield {
                    "ResourceARN": msk_cluster.arn,
                    "Tags": tags,
                }

        # Workspaces Web
        if self.workspacesweb_backends and (
            not resource_type_filters or "workspaces-web" in resource_type_filters
        ):
            for portal in self.workspacesweb_backends.portals.values():
                tags = self.workspacesweb_backends.tagger.list_tags_for_resource(
                    portal.arn
                )["Tags"]
                if not tags or not tag_filter(tags):
                    continue
                yield {
                    "ResourceARN": f"arn:{get_partition(self.region_name)}:workspaces-web:{self.region_name}:{self.account_id}:portal/{portal.portal_id}",
                    "Tags": tags,
                }

        # Lambda Instance
        if (
            not resource_type_filters
            or "lambda" in resource_type_filters
            or "lambda:function" in resource_type_filters
        ):
            for f in self.lambda_backend.list_functions():
                tags = format_tags(f.tags)
                if not tags or not tag_filter(tags):
                    continue
                yield {
                    "ResourceARN": f.function_arn,
                    "Tags": tags,
                }

        if (
            not resource_type_filters
            or "dynamodb" in resource_type_filters
            or "dynamodb:table" in resource_type_filters
        ):
            for table in self.dynamodb_backend.tables.values():
                tags = table.tags

                if not tags or not tag_filter(tags):
                    continue
                yield {
                    "ResourceARN": table.table_arn,
                    "Tags": tags,
                }

        # sagemaker cluster, automljob, compilation-job, domain, model-explainability-job-definition, model-quality-job-definition, and hyper-parameter-tuning-job currently supported
        sagemaker_resource_map: Dict[str, Dict[str, Any]] = {
            "sagemaker:cluster": self.sagemaker_backend.clusters,
            "sagemaker:automl-job": self.sagemaker_backend.auto_ml_jobs,
            "sagemaker:compilation-job": self.sagemaker_backend.compilation_jobs,
            "sagemaker:domain": self.sagemaker_backend.domains,
            "sagemaker:model-explainability-job-definition": self.sagemaker_backend.model_explainability_job_definitions,
            "sagemaker:model-quality-job-definition": self.sagemaker_backend.model_quality_job_definitions,
            "sagemaker:hyper-parameter-tuning-job": self.sagemaker_backend.hyper_parameter_tuning_jobs,
            "sagemaker:model-bias-job-definition": self.sagemaker_backend.model_bias_job_definitions,
            "sagemaker:data-quality-job-definition": self.sagemaker_backend.data_quality_job_definitions,
            "sagemaker:model": self.sagemaker_backend._models,
            "sagemaker:notebook-instance": self.sagemaker_backend.notebook_instances,
            "sagemaker:endpoint-config": self.sagemaker_backend.endpoint_configs,
            "sagemaker:endpoint": self.sagemaker_backend.endpoints,
            "sagemaker:experiment": self.sagemaker_backend.experiments,
            "sagemaker:pipeline": self.sagemaker_backend.pipelines,
            "sagemaker:pipeline-execution": self.sagemaker_backend.pipeline_executions,
            "sagemaker:processing-job": self.sagemaker_backend.processing_jobs,
            "sagemaker:trial": self.sagemaker_backend.trials,
            "sagemaker:trial-component": self.sagemaker_backend.trial_components,
            "sagemaker:training-job": self.sagemaker_backend.training_jobs,
            "sagemaker:transform-job": self.sagemaker_backend.transform_jobs,
            "sagemaker:notebook-instance-lifecycle-config": self.sagemaker_backend.notebook_instance_lifecycle_configurations,
            "sagemaker:model-card": self.sagemaker_backend.model_cards,
            "sagemaker:model-package-group": self.sagemaker_backend.model_package_groups,
            "sagemaker:model-package": self.sagemaker_backend.model_packages,
            "sagemaker:feature-group": self.sagemaker_backend.feature_groups,
        }
        for resource_type, resource_source in sagemaker_resource_map.items():
            if (
                not resource_type_filters
                or "sagemaker" in resource_type_filters
                or resource_type in resource_type_filters
            ):
                for resource in resource_source.values():
                    tags = self.sagemaker_backend.list_tags(resource.arn)[0]
                    if not tags or not tag_filter(tags):
                        continue
                    yield {
                        "ResourceARN": resource.arn,
                        "Tags": tags,
                    }

    def _get_tag_keys_generator(self) -> Iterator[str]:
        # Look at
        # https://docs.aws.amazon.com/general/latest/gr/aws-arns-and-namespaces.html

        # S3
        for bucket in self.s3_backend.buckets.values():
            tags = self.s3_backend.tagger.get_tag_dict_for_resource(bucket.arn)
            for key, _ in tags.items():
                yield key

        # EC2 tags
        def get_ec2_keys(res_id: str) -> List[Dict[str, str]]:
            result = []
            for key in self.ec2_backend.tags.get(res_id, {}):
                result.append(key)
            return result

        # EC2 AMI, resource type ec2:image
        for ami in self.ec2_backend.amis.values():
            for key in get_ec2_keys(ami.id):  # type: ignore[assignment]
                yield key

        # EC2 Instance, resource type ec2:instance
        for reservation in self.ec2_backend.reservations.values():
            for instance in reservation.instances:
                for key in get_ec2_keys(instance.id):  # type: ignore[assignment]
                    yield key

        # EC2 NetworkInterface, resource type ec2:network-interface
        for eni in self.ec2_backend.enis.values():
            for key in get_ec2_keys(eni.id):  # type: ignore[assignment]
                yield key

        # TODO EC2 ReservedInstance

        # EC2 SecurityGroup, resource type ec2:security-group
        for vpc in self.ec2_backend.groups.values():
            for sg in vpc.values():
                for key in get_ec2_keys(sg.id):  # type: ignore[assignment]
                    yield key

        # EC2 Snapshot, resource type ec2:snapshot
        for snapshot in self.ec2_backend.snapshots.values():
            for key in get_ec2_keys(snapshot.id):  # type: ignore[assignment]
                yield key

        # TODO EC2 SpotInstanceRequest

        # EC2 Volume, resource type ec2:volume
        for volume in self.ec2_backend.volumes.values():
            for key in get_ec2_keys(volume.id):  # type: ignore[assignment]
                yield key

        # Glue
        for tag_dict in self.glue_backend.tagger.tags.values():
            for tag_key in tag_dict.keys():
                yield tag_key

    def _get_tag_values_generator(self, tag_key: str) -> Iterator[str]:
        # Look at
        # https://docs.aws.amazon.com/general/latest/gr/aws-arns-and-namespaces.html

        # Do S3, resource type s3
        for bucket in self.s3_backend.buckets.values():
            tags = self.s3_backend.tagger.get_tag_dict_for_resource(bucket.arn)
            for key, value in tags.items():
                if key == tag_key:
                    yield value

        # EC2 tags
        def get_ec2_values(res_id: str) -> List[Dict[str, str]]:
            result = []
            for key, value in self.ec2_backend.tags.get(res_id, {}).items():
                if key == tag_key:
                    result.append(value)
            return result

        # EC2 AMI, resource type ec2:image
        for ami in self.ec2_backend.amis.values():
            for value in get_ec2_values(ami.id):  # type: ignore[assignment]
                yield value

        # EC2 Instance, resource type ec2:instance
        for reservation in self.ec2_backend.reservations.values():
            for instance in reservation.instances:
                for value in get_ec2_values(instance.id):  # type: ignore[assignment]
                    yield value

        # EC2 NetworkInterface, resource type ec2:network-interface
        for eni in self.ec2_backend.enis.values():
            for value in get_ec2_values(eni.id):  # type: ignore[assignment]
                yield value

        # TODO EC2 ReservedInstance

        # EC2 SecurityGroup, resource type ec2:security-group
        for vpc in self.ec2_backend.groups.values():
            for sg in vpc.values():
                for value in get_ec2_values(sg.id):  # type: ignore[assignment]
                    yield value

        # EC2 Snapshot, resource type ec2:snapshot
        for snapshot in self.ec2_backend.snapshots.values():
            for value in get_ec2_values(snapshot.id):  # type: ignore[assignment]
                yield value

        # TODO EC2 SpotInstanceRequest

        # EC2 Volume, resource type ec2:volume
        for volume in self.ec2_backend.volumes.values():
            for value in get_ec2_values(volume.id):  # type: ignore[assignment]
                yield value

        # Glue
        for tag_dict in self.glue_backend.tagger.tags.values():
            for key, tag_value in tag_dict.items():
                if key == tag_key and tag_value is not None:
                    yield tag_value

    def get_resources(
        self,
        pagination_token: Optional[str] = None,
        resources_per_page: int = 50,
        tags_per_page: int = 100,
        tag_filters: Optional[List[Dict[str, Any]]] = None,
        resource_type_filters: Optional[List[str]] = None,
    ) -> Tuple[Optional[str], List[Dict[str, Any]]]:
        # Simple range checking
        if 100 >= tags_per_page >= 500:
            raise RESTError(
                "InvalidParameterException", "TagsPerPage must be between 100 and 500"
            )
        if 1 >= resources_per_page >= 50:
            raise RESTError(
                "InvalidParameterException", "ResourcesPerPage must be between 1 and 50"
            )

        # If we have a token, go and find the respective generator, or error
        if pagination_token:
            if pagination_token not in self._pages:
                raise RESTError(
                    "PaginationTokenExpiredException", "Token does not exist"
                )

            generator = self._pages[pagination_token]["gen"]
            left_over = self._pages[pagination_token]["misc"]
        else:
            generator = self._get_resources_generator(
                tag_filters=tag_filters, resource_type_filters=resource_type_filters
            )
            left_over = None

        result = []
        current_tags = 0
        current_resources = 0
        if left_over:
            result.append(left_over)
            current_resources += 1
            current_tags += len(left_over["Tags"])

        try:
            while True:
                # Generator format: [{'ResourceARN': str, 'Tags': [{'Key': str, 'Value': str]}, ...]
                next_item = next(generator)
                resource_tags = len(next_item["Tags"])

                if current_resources >= resources_per_page:
                    break
                if current_tags + resource_tags >= tags_per_page:
                    break

                current_resources += 1
                current_tags += resource_tags

                result.append(next_item)

        except StopIteration:
            # Finished generator before invalidating page limiting constraints
            return None, result

        # Didn't hit StopIteration so there's stuff left in generator
        new_token = str(mock_random.uuid4())
        self._pages[new_token] = {"gen": generator, "misc": next_item}

        # Token used up, might as well bin now, if you call it again you're an idiot
        if pagination_token:
            del self._pages[pagination_token]
        return new_token, result

    def get_tag_keys(
        self, pagination_token: Optional[str] = None
    ) -> Tuple[Optional[str], List[str]]:
        if pagination_token:
            if pagination_token not in self._pages:
                raise RESTError(
                    "PaginationTokenExpiredException", "Token does not exist"
                )

            generator = self._pages[pagination_token]["gen"]
            left_over = self._pages[pagination_token]["misc"]
        else:
            generator = self._get_tag_keys_generator()
            left_over = None

        result = []
        current_tags = 0
        if left_over:
            result.append(left_over)
            current_tags += 1

        try:
            while True:
                # Generator format: ['tag', 'tag', 'tag', ...]
                next_item = next(generator)

                if current_tags + 1 >= 128:
                    break

                current_tags += 1

                result.append(next_item)

        except StopIteration:
            # Finished generator before invalidating page limiting constraints
            return None, result

        # Didn't hit StopIteration so there's stuff left in generator
        new_token = str(mock_random.uuid4())
        self._pages[new_token] = {"gen": generator, "misc": next_item}

        # Token used up, might as well bin now, if you call it again your an idiot
        if pagination_token:
            del self._pages[pagination_token]

        return new_token, result

    def get_tag_values(
        self, pagination_token: Optional[str], key: str
    ) -> Tuple[Optional[str], List[str]]:
        if pagination_token:
            if pagination_token not in self._pages:
                raise RESTError(
                    "PaginationTokenExpiredException", "Token does not exist"
                )

            generator = self._pages[pagination_token]["gen"]
            left_over = self._pages[pagination_token]["misc"]
        else:
            generator = self._get_tag_values_generator(key)
            left_over = None

        result = []
        current_tags = 0
        if left_over:
            result.append(left_over)
            current_tags += 1

        try:
            while True:
                # Generator format: ['value', 'value', 'value', ...]
                next_item = next(generator)

                if current_tags + 1 >= 128:
                    break

                current_tags += 1

                result.append(next_item)

        except StopIteration:
            # Finished generator before invalidating page limiting constraints
            return None, result

        # Didn't hit StopIteration so there's stuff left in generator
        new_token = str(mock_random.uuid4())
        self._pages[new_token] = {"gen": generator, "misc": next_item}

        # Token used up, might as well bin now, if you call it again your an idiot
        if pagination_token:
            del self._pages[pagination_token]

        return new_token, result

    def tag_resources(
        self, resource_arns: List[str], tags: Dict[str, str]
    ) -> Dict[str, Dict[str, Any]]:
        """
        Only DynamoDB, EFS, Lambda Logs, Quicksight RDS, and SageMaker resources are currently supported
        """
        missing_resources = []
        missing_error: Dict[str, Any] = {
            "StatusCode": 404,
            "ErrorCode": "InternalServiceException",
            "ErrorMessage": "Service not yet supported",
        }
        for arn in resource_arns:
            if arn.startswith(
                f"arn:{get_partition(self.region_name)}:rds:"
            ) or arn.startswith(f"arn:{get_partition(self.region_name)}:snapshot:"):
                self.rds_backend.add_tags_to_resource(
                    arn, TaggingService.convert_dict_to_tags_input(tags)
                )
            elif arn.startswith(
                f"arn:{get_partition(self.region_name)}:workspaces-web:"
            ):
                resource_id = arn.split("/")[-1]
                self.workspacesweb_backends.create_tags(  # type: ignore[union-attr]
                    resource_id, TaggingService.convert_dict_to_tags_input(tags)
                )
            elif arn.startswith(f"arn:{get_partition(self.region_name)}:workspaces:"):
                resource_id = arn.split("/")[-1]
                self.workspaces_backend.create_tags(  # type: ignore[union-attr]
                    resource_id, TaggingService.convert_dict_to_tags_input(tags)
                )
            elif arn.startswith(f"arn:{get_partition(self.region_name)}:logs:"):
                self.logs_backend.tag_resource(arn, tags)
            elif arn.startswith(f"arn:{get_partition(self.region_name)}:dynamodb"):
                self.dynamodb_backend.tag_resource(
                    arn, TaggingService.convert_dict_to_tags_input(tags)
                )
            elif arn.startswith(f"arn:{get_partition(self.region_name)}:sagemaker:"):
                self.sagemaker_backend.add_tags(
                    arn, TaggingService.convert_dict_to_tags_input(tags)
                )
            elif arn.startswith(f"arn:{get_partition(self.region_name)}:lambda:"):
                self.lambda_backend.tag_resource(arn, tags)
            elif arn.startswith(
                f"arn:{get_partition(self.region_name)}:elasticfilesystem:"
            ):
                resource_id = arn.split("/")[-1]
                self.efs_backend.tag_resource(
                    resource_id, TaggingService.convert_dict_to_tags_input(tags)
                )
            elif arn.startswith(f"arn:{get_partition(self.region_name)}:quicksight:"):
                assert self.quicksight_backend is not None
                self.quicksight_backend.tag_resource(
                    arn, TaggingService.convert_dict_to_tags_input(tags)
                )
            else:
                missing_resources.append(arn)
        return {arn: missing_error for arn in missing_resources}

    def untag_resources(
        self, resource_arn_list: List[str], tag_keys: List[str]
    ) -> Dict[str, Dict[str, Any]]:
        """
        Only EFS, Lambda, and Quicksight resources are currently supported
        """
        missing_resources = []
        missing_error: Dict[str, Any] = {
            "StatusCode": 404,
            "ErrorCode": "InternalServiceException",
            "ErrorMessage": "Service not yet supported",
        }

        for arn in resource_arn_list:
            if arn.startswith(f"arn:{get_partition(self.region_name)}:lambda:"):
                self.lambda_backend.untag_resource(arn, tag_keys)
            elif arn.startswith(
                f"arn:{get_partition(self.region_name)}:elasticfilesystem:"
            ):
                resource_id = arn.split("/")[-1]
                self.efs_backend.untag_resource(resource_id, tag_keys)
            elif arn.startswith(f"arn:{get_partition(self.region_name)}:quicksight:"):
                assert self.quicksight_backend is not None
                self.quicksight_backend.untag_resource(arn, tag_keys)
            else:
                missing_resources.append(arn)

        return {arn: missing_error for arn in missing_resources}


resourcegroupstaggingapi_backends = BackendDict(
    ResourceGroupsTaggingAPIBackend, "resourcegroupstaggingapi"
)
