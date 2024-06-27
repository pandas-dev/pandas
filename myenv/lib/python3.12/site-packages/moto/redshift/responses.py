import json
from typing import Any, Dict, List

import xmltodict
from jinja2 import Template

from moto.core.common_types import TYPE_RESPONSE
from moto.core.responses import BaseResponse

from .models import RedshiftBackend, redshift_backends


def convert_json_error_to_xml(json_error: Any) -> str:
    error = json.loads(json_error)
    code = error["Error"]["Code"]
    message = error["Error"]["Message"]
    template = Template(
        """
        <RedshiftClientError>
            <Error>
              <Code>{{ code }}</Code>
              <Message>{{ message }}</Message>
              <Type>Sender</Type>
            </Error>
            <RequestId>6876f774-7273-11e4-85dc-39e55ca848d1</RequestId>
        </RedshiftClientError>"""
    )
    return template.render(code=code, message=message)


def itemize(data: Any) -> Dict[str, Any]:
    """
    The xmltodict.unparse requires we modify the shape of the input dictionary slightly. Instead of a dict of the form:
        {'key': ['value1', 'value2']}
    We must provide:
        {'key': {'item': ['value1', 'value2']}}
    """
    if isinstance(data, dict):
        ret = {}
        for key in data:
            ret[key] = itemize(data[key])
        return ret
    elif isinstance(data, list):
        return {"item": [itemize(value) for value in data]}
    else:
        return data


class RedshiftResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="redshift")

    @property
    def redshift_backend(self) -> RedshiftBackend:
        return redshift_backends[self.current_account][self.region]

    def get_response(self, response: Any) -> str:
        if self.request_json:
            return json.dumps(response)
        else:
            xml = xmltodict.unparse(itemize(response), full_document=False)
            if hasattr(xml, "decode"):
                xml = xml.decode("utf-8")
            return xml

    def call_action(self) -> TYPE_RESPONSE:
        status, headers, body = super().call_action()
        if status >= 400 and not self.request_json:
            body = convert_json_error_to_xml(body)
        return status, headers, body

    def unpack_list_params(self, label: str, child_label: str) -> Any:
        root = self._get_multi_param_dict(label) or {}
        return root.get(child_label, [])

    def _get_cluster_security_groups(self) -> List[str]:
        cluster_security_groups = self._get_multi_param("ClusterSecurityGroups.member")
        if not cluster_security_groups:
            cluster_security_groups = self._get_multi_param(
                "ClusterSecurityGroups.ClusterSecurityGroupName"
            )
        return cluster_security_groups

    def _get_vpc_security_group_ids(self) -> List[str]:
        vpc_security_group_ids = self._get_multi_param("VpcSecurityGroupIds.member")
        if not vpc_security_group_ids:
            vpc_security_group_ids = self._get_multi_param(
                "VpcSecurityGroupIds.VpcSecurityGroupId"
            )
        return vpc_security_group_ids

    def _get_iam_roles(self) -> List[str]:
        iam_roles = self._get_multi_param("IamRoles.member")
        if not iam_roles:
            iam_roles = self._get_multi_param("IamRoles.IamRoleArn")
        return iam_roles

    def _get_subnet_ids(self) -> List[str]:
        subnet_ids = self._get_multi_param("SubnetIds.member")
        if not subnet_ids:
            subnet_ids = self._get_multi_param("SubnetIds.SubnetIdentifier")
        return subnet_ids

    def create_cluster(self) -> str:
        cluster_kwargs = {
            "cluster_identifier": self._get_param("ClusterIdentifier"),
            "node_type": self._get_param("NodeType"),
            "master_username": self._get_param("MasterUsername"),
            "master_user_password": self._get_param("MasterUserPassword"),
            "db_name": self._get_param("DBName"),
            "cluster_type": self._get_param("ClusterType"),
            "cluster_security_groups": self._get_cluster_security_groups(),
            "vpc_security_group_ids": self._get_vpc_security_group_ids(),
            "cluster_subnet_group_name": self._get_param("ClusterSubnetGroupName"),
            "availability_zone": self._get_param("AvailabilityZone"),
            "preferred_maintenance_window": self._get_param(
                "PreferredMaintenanceWindow"
            ),
            "cluster_parameter_group_name": self._get_param(
                "ClusterParameterGroupName"
            ),
            "automated_snapshot_retention_period": self._get_int_param(
                "AutomatedSnapshotRetentionPeriod"
            ),
            "port": self._get_int_param("Port"),
            "cluster_version": self._get_param("ClusterVersion"),
            "allow_version_upgrade": self._get_bool_param("AllowVersionUpgrade"),
            "number_of_nodes": self._get_int_param("NumberOfNodes"),
            "publicly_accessible": self._get_param("PubliclyAccessible"),
            "encrypted": self._get_param("Encrypted"),
            "region_name": self.region,
            "tags": self.unpack_list_params("Tags", "Tag"),
            "iam_roles_arn": self._get_iam_roles(),
            "enhanced_vpc_routing": self._get_param("EnhancedVpcRouting"),
            "kms_key_id": self._get_param("KmsKeyId"),
        }
        cluster = self.redshift_backend.create_cluster(**cluster_kwargs).to_json()
        cluster["ClusterStatus"] = "creating"
        return self.get_response(
            {
                "CreateClusterResponse": {
                    "CreateClusterResult": {"Cluster": cluster},
                    "ResponseMetadata": {
                        "RequestId": "384ac68d-3775-11df-8963-01868b7c937a"
                    },
                }
            }
        )

    def pause_cluster(self) -> str:
        cluster_id = self._get_param("ClusterIdentifier")
        cluster = self.redshift_backend.pause_cluster(cluster_id).to_json()
        return self.get_response(
            {
                "PauseClusterResponse": {
                    "PauseClusterResult": {"Cluster": cluster},
                    "ResponseMetadata": {
                        "RequestId": "384ac68d-3775-11df-8963-01868b7c937a"
                    },
                }
            }
        )

    def resume_cluster(self) -> str:
        cluster_id = self._get_param("ClusterIdentifier")
        cluster = self.redshift_backend.resume_cluster(cluster_id).to_json()
        return self.get_response(
            {
                "ResumeClusterResponse": {
                    "ResumeClusterResult": {"Cluster": cluster},
                    "ResponseMetadata": {
                        "RequestId": "384ac68d-3775-11df-8963-01868b7c937a"
                    },
                }
            }
        )

    def restore_from_cluster_snapshot(self) -> str:
        enhanced_vpc_routing = self._get_bool_param("EnhancedVpcRouting")
        node_type = self._get_param("NodeType")
        number_of_nodes = self._get_int_param("NumberOfNodes")
        restore_kwargs = {
            "snapshot_identifier": self._get_param("SnapshotIdentifier"),
            "cluster_identifier": self._get_param("ClusterIdentifier"),
            "port": self._get_int_param("Port"),
            "availability_zone": self._get_param("AvailabilityZone"),
            "allow_version_upgrade": self._get_bool_param("AllowVersionUpgrade"),
            "cluster_subnet_group_name": self._get_param("ClusterSubnetGroupName"),
            "publicly_accessible": self._get_param("PubliclyAccessible"),
            "cluster_parameter_group_name": self._get_param(
                "ClusterParameterGroupName"
            ),
            "cluster_security_groups": self._get_cluster_security_groups(),
            "vpc_security_group_ids": self._get_vpc_security_group_ids(),
            "preferred_maintenance_window": self._get_param(
                "PreferredMaintenanceWindow"
            ),
            "automated_snapshot_retention_period": self._get_int_param(
                "AutomatedSnapshotRetentionPeriod"
            ),
            "region_name": self.region,
            "iam_roles_arn": self._get_iam_roles(),
        }
        if enhanced_vpc_routing is not None:
            restore_kwargs["enhanced_vpc_routing"] = enhanced_vpc_routing
        if node_type is not None:
            restore_kwargs["node_type"] = node_type
        if number_of_nodes is not None:
            restore_kwargs["number_of_nodes"] = number_of_nodes
        cluster = self.redshift_backend.restore_from_cluster_snapshot(
            **restore_kwargs
        ).to_json()
        cluster["ClusterStatus"] = "creating"
        return self.get_response(
            {
                "RestoreFromClusterSnapshotResponse": {
                    "RestoreFromClusterSnapshotResult": {"Cluster": cluster},
                    "ResponseMetadata": {
                        "RequestId": "384ac68d-3775-11df-8963-01868b7c937a"
                    },
                }
            }
        )

    def describe_clusters(self) -> str:
        cluster_identifier = self._get_param("ClusterIdentifier")
        clusters = self.redshift_backend.describe_clusters(cluster_identifier)

        return self.get_response(
            {
                "DescribeClustersResponse": {
                    "DescribeClustersResult": {
                        "Clusters": [cluster.to_json() for cluster in clusters]
                    },
                    "ResponseMetadata": {
                        "RequestId": "384ac68d-3775-11df-8963-01868b7c937a"
                    },
                }
            }
        )

    def modify_cluster(self) -> str:
        request_kwargs = {
            "cluster_identifier": self._get_param("ClusterIdentifier"),
            "new_cluster_identifier": self._get_param("NewClusterIdentifier"),
            "node_type": self._get_param("NodeType"),
            "master_user_password": self._get_param("MasterUserPassword"),
            "cluster_type": self._get_param("ClusterType"),
            "cluster_security_groups": self._get_cluster_security_groups(),
            "vpc_security_group_ids": self._get_vpc_security_group_ids(),
            "cluster_subnet_group_name": self._get_param("ClusterSubnetGroupName"),
            "preferred_maintenance_window": self._get_param(
                "PreferredMaintenanceWindow"
            ),
            "cluster_parameter_group_name": self._get_param(
                "ClusterParameterGroupName"
            ),
            "automated_snapshot_retention_period": self._get_int_param(
                "AutomatedSnapshotRetentionPeriod"
            ),
            "cluster_version": self._get_param("ClusterVersion"),
            "allow_version_upgrade": self._get_bool_param("AllowVersionUpgrade"),
            "number_of_nodes": self._get_int_param("NumberOfNodes"),
            "publicly_accessible": self._get_param("PubliclyAccessible"),
            "encrypted": self._get_param("Encrypted"),
            "iam_roles_arn": self._get_iam_roles(),
            "enhanced_vpc_routing": self._get_param("EnhancedVpcRouting"),
        }
        cluster_kwargs = {}
        # We only want parameters that were actually passed in, otherwise
        # we'll stomp all over our cluster metadata with None values.
        for key, value in request_kwargs.items():
            if value is not None and value != []:
                cluster_kwargs[key] = value

        cluster = self.redshift_backend.modify_cluster(**cluster_kwargs)

        return self.get_response(
            {
                "ModifyClusterResponse": {
                    "ModifyClusterResult": {"Cluster": cluster.to_json()},
                    "ResponseMetadata": {
                        "RequestId": "384ac68d-3775-11df-8963-01868b7c937a"
                    },
                }
            }
        )

    def delete_cluster(self) -> str:
        request_kwargs = {
            "cluster_identifier": self._get_param("ClusterIdentifier"),
            "final_cluster_snapshot_identifier": self._get_param(
                "FinalClusterSnapshotIdentifier"
            ),
            "skip_final_snapshot": self._get_bool_param("SkipFinalClusterSnapshot"),
        }

        cluster = self.redshift_backend.delete_cluster(**request_kwargs)

        return self.get_response(
            {
                "DeleteClusterResponse": {
                    "DeleteClusterResult": {"Cluster": cluster.to_json()},
                    "ResponseMetadata": {
                        "RequestId": "384ac68d-3775-11df-8963-01868b7c937a"
                    },
                }
            }
        )

    def create_cluster_subnet_group(self) -> str:
        cluster_subnet_group_name = self._get_param("ClusterSubnetGroupName")
        description = self._get_param("Description")
        subnet_ids = self._get_subnet_ids()
        tags = self.unpack_list_params("Tags", "Tag")

        subnet_group = self.redshift_backend.create_cluster_subnet_group(
            cluster_subnet_group_name=cluster_subnet_group_name,
            description=description,
            subnet_ids=subnet_ids,
            region_name=self.region,
            tags=tags,
        )

        return self.get_response(
            {
                "CreateClusterSubnetGroupResponse": {
                    "CreateClusterSubnetGroupResult": {
                        "ClusterSubnetGroup": subnet_group.to_json()
                    },
                    "ResponseMetadata": {
                        "RequestId": "384ac68d-3775-11df-8963-01868b7c937a"
                    },
                }
            }
        )

    def describe_cluster_subnet_groups(self) -> str:
        subnet_identifier = self._get_param("ClusterSubnetGroupName")
        subnet_groups = self.redshift_backend.describe_cluster_subnet_groups(
            subnet_identifier
        )

        return self.get_response(
            {
                "DescribeClusterSubnetGroupsResponse": {
                    "DescribeClusterSubnetGroupsResult": {
                        "ClusterSubnetGroups": [
                            subnet_group.to_json() for subnet_group in subnet_groups
                        ]
                    },
                    "ResponseMetadata": {
                        "RequestId": "384ac68d-3775-11df-8963-01868b7c937a"
                    },
                }
            }
        )

    def delete_cluster_subnet_group(self) -> str:
        subnet_identifier = self._get_param("ClusterSubnetGroupName")
        self.redshift_backend.delete_cluster_subnet_group(subnet_identifier)

        return self.get_response(
            {
                "DeleteClusterSubnetGroupResponse": {
                    "ResponseMetadata": {
                        "RequestId": "384ac68d-3775-11df-8963-01868b7c937a"
                    }
                }
            }
        )

    def create_cluster_security_group(self) -> str:
        cluster_security_group_name = self._get_param("ClusterSecurityGroupName")
        description = self._get_param("Description")
        tags = self.unpack_list_params("Tags", "Tag")

        security_group = self.redshift_backend.create_cluster_security_group(
            cluster_security_group_name=cluster_security_group_name,
            description=description,
            tags=tags,
        )

        return self.get_response(
            {
                "CreateClusterSecurityGroupResponse": {
                    "CreateClusterSecurityGroupResult": {
                        "ClusterSecurityGroup": security_group.to_json()
                    },
                    "ResponseMetadata": {
                        "RequestId": "384ac68d-3775-11df-8963-01868b7c937a"
                    },
                }
            }
        )

    def describe_cluster_security_groups(self) -> str:
        cluster_security_group_name = self._get_param("ClusterSecurityGroupName")
        security_groups = self.redshift_backend.describe_cluster_security_groups(
            cluster_security_group_name
        )

        return self.get_response(
            {
                "DescribeClusterSecurityGroupsResponse": {
                    "DescribeClusterSecurityGroupsResult": {
                        "ClusterSecurityGroups": [
                            security_group.to_json()
                            for security_group in security_groups
                        ]
                    },
                    "ResponseMetadata": {
                        "RequestId": "384ac68d-3775-11df-8963-01868b7c937a"
                    },
                }
            }
        )

    def delete_cluster_security_group(self) -> str:
        security_group_identifier = self._get_param("ClusterSecurityGroupName")
        self.redshift_backend.delete_cluster_security_group(security_group_identifier)

        return self.get_response(
            {
                "DeleteClusterSecurityGroupResponse": {
                    "ResponseMetadata": {
                        "RequestId": "384ac68d-3775-11df-8963-01868b7c937a"
                    }
                }
            }
        )

    def authorize_cluster_security_group_ingress(self) -> str:
        cluster_security_group_name = self._get_param("ClusterSecurityGroupName")
        cidr_ip = self._get_param("CIDRIP")

        security_group = self.redshift_backend.authorize_cluster_security_group_ingress(
            cluster_security_group_name, cidr_ip
        )

        return self.get_response(
            {
                "AuthorizeClusterSecurityGroupIngressResponse": {
                    "AuthorizeClusterSecurityGroupIngressResult": {
                        "ClusterSecurityGroup": {
                            "ClusterSecurityGroupName": cluster_security_group_name,
                            "Description": security_group.description,
                            "IPRanges": [
                                {
                                    "Status": "authorized",
                                    "CIDRIP": cidr_ip,
                                    "Tags": security_group.tags,
                                },
                            ],
                        }
                    }
                }
            }
        )

    def create_cluster_parameter_group(self) -> str:
        cluster_parameter_group_name = self._get_param("ParameterGroupName")
        group_family = self._get_param("ParameterGroupFamily")
        description = self._get_param("Description")
        tags = self.unpack_list_params("Tags", "Tag")

        parameter_group = self.redshift_backend.create_cluster_parameter_group(
            cluster_parameter_group_name, group_family, description, self.region, tags
        )

        return self.get_response(
            {
                "CreateClusterParameterGroupResponse": {
                    "CreateClusterParameterGroupResult": {
                        "ClusterParameterGroup": parameter_group.to_json()
                    },
                    "ResponseMetadata": {
                        "RequestId": "384ac68d-3775-11df-8963-01868b7c937a"
                    },
                }
            }
        )

    def describe_cluster_parameter_groups(self) -> str:
        cluster_parameter_group_name = self._get_param("ParameterGroupName")
        parameter_groups = self.redshift_backend.describe_cluster_parameter_groups(
            cluster_parameter_group_name
        )

        return self.get_response(
            {
                "DescribeClusterParameterGroupsResponse": {
                    "DescribeClusterParameterGroupsResult": {
                        "ParameterGroups": [
                            parameter_group.to_json()
                            for parameter_group in parameter_groups
                        ]
                    },
                    "ResponseMetadata": {
                        "RequestId": "384ac68d-3775-11df-8963-01868b7c937a"
                    },
                }
            }
        )

    def delete_cluster_parameter_group(self) -> str:
        cluster_parameter_group_name = self._get_param("ParameterGroupName")
        self.redshift_backend.delete_cluster_parameter_group(
            cluster_parameter_group_name
        )

        return self.get_response(
            {
                "DeleteClusterParameterGroupResponse": {
                    "ResponseMetadata": {
                        "RequestId": "384ac68d-3775-11df-8963-01868b7c937a"
                    }
                }
            }
        )

    def create_cluster_snapshot(self) -> str:
        cluster_identifier = self._get_param("ClusterIdentifier")
        snapshot_identifier = self._get_param("SnapshotIdentifier")
        tags = self.unpack_list_params("Tags", "Tag")

        snapshot = self.redshift_backend.create_cluster_snapshot(
            cluster_identifier, snapshot_identifier, self.region, tags
        )
        return self.get_response(
            {
                "CreateClusterSnapshotResponse": {
                    "CreateClusterSnapshotResult": {"Snapshot": snapshot.to_json()},
                    "ResponseMetadata": {
                        "RequestId": "384ac68d-3775-11df-8963-01868b7c937a"
                    },
                }
            }
        )

    def describe_cluster_snapshots(self) -> str:
        cluster_identifier = self._get_param("ClusterIdentifier")
        snapshot_identifier = self._get_param("SnapshotIdentifier")
        snapshot_type = self._get_param("SnapshotType")
        snapshots = self.redshift_backend.describe_cluster_snapshots(
            cluster_identifier, snapshot_identifier, snapshot_type
        )
        return self.get_response(
            {
                "DescribeClusterSnapshotsResponse": {
                    "DescribeClusterSnapshotsResult": {
                        "Snapshots": [snapshot.to_json() for snapshot in snapshots]
                    },
                    "ResponseMetadata": {
                        "RequestId": "384ac68d-3775-11df-8963-01868b7c937a"
                    },
                }
            }
        )

    def delete_cluster_snapshot(self) -> str:
        snapshot_identifier = self._get_param("SnapshotIdentifier")
        snapshot = self.redshift_backend.delete_cluster_snapshot(snapshot_identifier)

        return self.get_response(
            {
                "DeleteClusterSnapshotResponse": {
                    "DeleteClusterSnapshotResult": {"Snapshot": snapshot.to_json()},
                    "ResponseMetadata": {
                        "RequestId": "384ac68d-3775-11df-8963-01868b7c937a"
                    },
                }
            }
        )

    def create_snapshot_copy_grant(self) -> str:
        copy_grant_kwargs = {
            "snapshot_copy_grant_name": self._get_param("SnapshotCopyGrantName"),
            "kms_key_id": self._get_param("KmsKeyId"),
            "region_name": self._get_param("Region"),
        }

        copy_grant = self.redshift_backend.create_snapshot_copy_grant(
            **copy_grant_kwargs
        )
        return self.get_response(
            {
                "CreateSnapshotCopyGrantResponse": {
                    "CreateSnapshotCopyGrantResult": {
                        "SnapshotCopyGrant": copy_grant.to_json()
                    },
                    "ResponseMetadata": {
                        "RequestId": "384ac68d-3775-11df-8963-01868b7c937a"
                    },
                }
            }
        )

    def delete_snapshot_copy_grant(self) -> str:
        copy_grant_kwargs = {
            "snapshot_copy_grant_name": self._get_param("SnapshotCopyGrantName")
        }
        self.redshift_backend.delete_snapshot_copy_grant(**copy_grant_kwargs)
        return self.get_response(
            {
                "DeleteSnapshotCopyGrantResponse": {
                    "ResponseMetadata": {
                        "RequestId": "384ac68d-3775-11df-8963-01868b7c937a"
                    }
                }
            }
        )

    def describe_snapshot_copy_grants(self) -> str:
        copy_grant_kwargs = {
            "snapshot_copy_grant_name": self._get_param("SnapshotCopyGrantName")
        }

        copy_grants = self.redshift_backend.describe_snapshot_copy_grants(
            **copy_grant_kwargs
        )
        return self.get_response(
            {
                "DescribeSnapshotCopyGrantsResponse": {
                    "DescribeSnapshotCopyGrantsResult": {
                        "SnapshotCopyGrants": [
                            copy_grant.to_json() for copy_grant in copy_grants
                        ]
                    },
                    "ResponseMetadata": {
                        "RequestId": "384ac68d-3775-11df-8963-01868b7c937a"
                    },
                }
            }
        )

    def create_tags(self) -> str:
        resource_name = self._get_param("ResourceName")
        tags = self.unpack_list_params("Tags", "Tag")

        self.redshift_backend.create_tags(resource_name, tags)

        return self.get_response(
            {
                "CreateTagsResponse": {
                    "ResponseMetadata": {
                        "RequestId": "384ac68d-3775-11df-8963-01868b7c937a"
                    }
                }
            }
        )

    def describe_tags(self) -> str:
        resource_name = self._get_param("ResourceName")
        resource_type = self._get_param("ResourceType")

        tagged_resources = self.redshift_backend.describe_tags(
            resource_name, resource_type
        )
        return self.get_response(
            {
                "DescribeTagsResponse": {
                    "DescribeTagsResult": {"TaggedResources": tagged_resources},
                    "ResponseMetadata": {
                        "RequestId": "384ac68d-3775-11df-8963-01868b7c937a"
                    },
                }
            }
        )

    def delete_tags(self) -> str:
        resource_name = self._get_param("ResourceName")
        tag_keys = self.unpack_list_params("TagKeys", "TagKey")

        self.redshift_backend.delete_tags(resource_name, tag_keys)

        return self.get_response(
            {
                "DeleteTagsResponse": {
                    "ResponseMetadata": {
                        "RequestId": "384ac68d-3775-11df-8963-01868b7c937a"
                    }
                }
            }
        )

    def enable_snapshot_copy(self) -> str:
        snapshot_copy_kwargs = {
            "cluster_identifier": self._get_param("ClusterIdentifier"),
            "destination_region": self._get_param("DestinationRegion"),
            "retention_period": self._get_param("RetentionPeriod", 7),
            "snapshot_copy_grant_name": self._get_param("SnapshotCopyGrantName"),
        }
        cluster = self.redshift_backend.enable_snapshot_copy(**snapshot_copy_kwargs)

        return self.get_response(
            {
                "EnableSnapshotCopyResponse": {
                    "EnableSnapshotCopyResult": {"Cluster": cluster.to_json()},
                    "ResponseMetadata": {
                        "RequestId": "384ac68d-3775-11df-8963-01868b7c937a"
                    },
                }
            }
        )

    def disable_snapshot_copy(self) -> str:
        snapshot_copy_kwargs = {
            "cluster_identifier": self._get_param("ClusterIdentifier")
        }
        cluster = self.redshift_backend.disable_snapshot_copy(**snapshot_copy_kwargs)

        return self.get_response(
            {
                "DisableSnapshotCopyResponse": {
                    "DisableSnapshotCopyResult": {"Cluster": cluster.to_json()},
                    "ResponseMetadata": {
                        "RequestId": "384ac68d-3775-11df-8963-01868b7c937a"
                    },
                }
            }
        )

    def modify_snapshot_copy_retention_period(self) -> str:
        snapshot_copy_kwargs = {
            "cluster_identifier": self._get_param("ClusterIdentifier"),
            "retention_period": self._get_param("RetentionPeriod"),
        }
        cluster = self.redshift_backend.modify_snapshot_copy_retention_period(
            **snapshot_copy_kwargs
        )

        return self.get_response(
            {
                "ModifySnapshotCopyRetentionPeriodResponse": {
                    "ModifySnapshotCopyRetentionPeriodResult": {
                        "Clusters": [cluster.to_json()]
                    },
                    "ResponseMetadata": {
                        "RequestId": "384ac68d-3775-11df-8963-01868b7c937a"
                    },
                }
            }
        )

    def get_cluster_credentials(self) -> str:
        cluster_identifier = self._get_param("ClusterIdentifier")
        db_user = self._get_param("DbUser")
        auto_create = self._get_bool_param("AutoCreate", False)
        duration_seconds = self._get_int_param("DurationSeconds", 900)

        cluster_credentials = self.redshift_backend.get_cluster_credentials(
            cluster_identifier, db_user, auto_create, duration_seconds
        )

        return self.get_response(
            {
                "GetClusterCredentialsResponse": {
                    "GetClusterCredentialsResult": cluster_credentials,
                    "ResponseMetadata": {
                        "RequestId": "384ac68d-3775-11df-8963-01868b7c937a"
                    },
                }
            }
        )
