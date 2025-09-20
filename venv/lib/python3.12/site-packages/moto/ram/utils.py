from datetime import datetime, timezone
from typing import Optional, TypedDict


class ManagedPermissionDict(TypedDict):
    arn: str
    version: str
    defaultVersion: bool
    name: str
    resourceType: str
    status: str
    creationTime: str
    lastUpdatedTime: str
    isResourceTypeDefault: bool
    permissionType: str


def format_ram_permission(
    name: str,
    resource_type: str,
    version: str = "1",
    arn_prefix: str = "arn:aws:ram::aws:permission/",
    status: str = "ATTACHABLE",
    creation_time: Optional[str] = None,
    last_updated_time: Optional[str] = None,
    is_resource_type_default: bool = True,
    permission_type: str = "AWS_MANAGED",
    default_version: bool = True,
) -> ManagedPermissionDict:
    """
    Format a RAM (Resource Access Manager) permission dictionary with the
    specified attributes.

    Args:
        name (str): The name of the permission.
        resource_type (str): The type of resource the permission applies to.
        version (str): The version of the permission. Defaults to "1".
        arn_prefix (str): The prefix for the permission ARN. Defaults
            to "arn:aws:ram::aws:permission/".
        status (str): The status of the permission. Defaults to
            "ATTACHABLE".
        creation_time (str, optional): The creation time in UTC. If not
            provided, uses the current time.
        last_updated_time (str, optional): The last updated time in UTC. If
            not provided, uses the creation time.
        is_resource_type_default (bool): Whether this is the default
            permission for the resource type. Defaults to True.
        permission_type (str): The type of permission (e.g.,
            "AWS_MANAGED"). Defaults to "AWS_MANAGED".
        default_version (bool): Whether this is the default version.
            Defaults to True.

    Returns:
        ManagedPermissionDict: A dictionary representing the formatted RAM permission.
    """
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S.%f")[:-3]
    creation_time = creation_time or now
    last_updated_time = last_updated_time or creation_time
    arn = f"{arn_prefix}{name}"
    return {
        "arn": arn,
        "version": version,
        "defaultVersion": default_version,
        "name": name,
        "resourceType": resource_type,
        "status": status,
        "creationTime": creation_time,
        "lastUpdatedTime": last_updated_time,
        "isResourceTypeDefault": is_resource_type_default,
        "permissionType": permission_type,
    }


# List of tuples with resourceType, serviceName, and resourceRegionScope based
# on the SHARED_RESOURCE_TYPES in models.py
_RESOURCE_TYPE_DEFINITIONS = [
    ("rds:Cluster", "rds", "REGIONAL"),
    ("imagebuilder:Component", "imagebuilder", "REGIONAL"),
    ("networkmanager:CoreNetwork", "networkmanager", "GLOBAL"),
    ("resource-groups:Group", "resource-groups", "REGIONAL"),
    ("imagebuilder:Image", "imagebuilder", "REGIONAL"),
    ("imagebuilder:ImageRecipe", "imagebuilder", "REGIONAL"),
    ("license-manager:LicenseConfiguration", "license-manager", "REGIONAL"),
    ("appmesh:Mesh", "appmesh", "REGIONAL"),
    ("ec2:PrefixList", "ec2", "REGIONAL"),
    ("codebuild:Project", "codebuild", "REGIONAL"),
    ("codebuild:ReportGroup", "codebuild", "REGIONAL"),
    ("route53resolver:ResolverRule", "route53resolver", "REGIONAL"),
    ("ec2:Subnet", "ec2", "REGIONAL"),
    ("ec2:TransitGatewayMulticastDomain", "ec2", "REGIONAL"),
    ("glue:Database", "glue", "REGIONAL"),
    ("glue:Table", "glue", "REGIONAL"),
    ("glue:Catalog", "glue", "REGIONAL"),
]

# List of dictionaries representing RAM resource types based on the definitions
RAM_RESOURCE_TYPES = [
    {
        "resourceType": resource_type,
        "serviceName": service_name,
        "resourceRegionScope": region_scope,
    }
    for resource_type, service_name, region_scope in _RESOURCE_TYPE_DEFINITIONS
]

# List of AWS managed permissions for RAM, formatted using the
# format_ram_permission function.
AWS_MANAGED_PERMISSIONS = [
    format_ram_permission(
        name="AWSRAMDefaultPermissionRDSCluster",
        resource_type="rds:Cluster",
        creation_time="2022-06-30 17:04:16.671000",
    ),
    format_ram_permission(
        name="AWSRAMDefaultPermissionImageBuilderComponent",
        resource_type="imagebuilder:Component",
        creation_time="2022-06-30 17:04:21.977000",
    ),
    format_ram_permission(
        name="AWSRAMDefaultPermissionsNetworkManagerCoreNetwork",
        resource_type="networkmanager:CoreNetwork",
        version="3",
        creation_time="2024-10-23 16:34:51.604000",
    ),
    format_ram_permission(
        name="AWSRAMTransitGatewayPermissionsNetworkManagerCoreNetwork",
        resource_type="networkmanager:CoreNetwork",
        is_resource_type_default=False,
        creation_time="2022-06-30 17:03:48.071000",
    ),
    format_ram_permission(
        name="AWSRAMVPCPermissionsNetworkManagerCoreNetwork",
        resource_type="networkmanager:CoreNetwork",
        is_resource_type_default=False,
        creation_time="2022-06-30 17:03:46.477000",
    ),
    format_ram_permission(
        name="AWSRAMDefaultPermissionResourceGroup",
        resource_type="resource-groups:Group",
        creation_time="2022-06-30 17:04:10.914000",
    ),
    format_ram_permission(
        name="AWSRAMDefaultPermissionImageBuilderImage",
        resource_type="imagebuilder:Image",
        creation_time="2022-06-30 17:04:21.201000",
    ),
    format_ram_permission(
        name="AWSRAMDefaultPermissionImageBuilderImageRecipe",
        resource_type="imagebuilder:ImageRecipe",
        creation_time="2022-06-30 17:04:22.832000",
    ),
    format_ram_permission(
        name="AWSRAMDefaultPermissionLicenseConfiguration",
        resource_type="license-manager:LicenseConfiguration",
        creation_time="2024-01-31 16:36:13.664000",
    ),
    format_ram_permission(
        name="AWSRAMDefaultPermissionAppMesh",
        resource_type="appmesh:Mesh",
        version="3",
        creation_time="2023-03-28 12:27:55.910000",
    ),
    format_ram_permission(
        name="AWSRAMDefaultPermissionPrefixList",
        resource_type="ec2:PrefixList",
        creation_time="2022-06-30 17:04:20.594000",
    ),
    format_ram_permission(
        name="AWSRAMDefaultPermissionCodeBuildProject",
        resource_type="codebuild:Project",
        creation_time="2022-06-30 17:04:17.595000",
    ),
    format_ram_permission(
        name="AWSRAMDefaultPermissionCodeBuildReportGroup",
        resource_type="codebuild:ReportGroup",
        creation_time="2022-06-30 17:04:17.761000",
    ),
    format_ram_permission(
        name="AWSRAMDefaultPermissionResolverRule",
        resource_type="route53resolver:ResolverRule",
        creation_time="2022-06-30 17:04:14.752000",
    ),
    format_ram_permission(
        name="AWSRAMDefaultPermissionSubnet",
        resource_type="ec2:Subnet",
        creation_time="2024-10-30 15:11:41.358000",
    ),
    format_ram_permission(
        name="AWSRAMDefaultPermissionTransitGatewayMulticastDomain",
        resource_type="ec2:TransitGatewayMulticastDomain",
        creation_time="2022-06-30 17:04:13.461000",
    ),
    format_ram_permission(
        name="AWSRAMDefaultPermissionGlueDatabase",
        resource_type="glue:Database",
        creation_time="2022-06-30 17:04:25.912000",
    ),
    format_ram_permission(
        name="AWSRAMLFEnabledGlueAllTablesReadWriteForDatabase",
        resource_type="glue:Database",
        is_resource_type_default=False,
        creation_time="2023-06-27 17:05:24.471000",
    ),
    format_ram_permission(
        name="AWSRAMLFEnabledGlueDatabaseReadWrite",
        resource_type="glue:Database",
        is_resource_type_default=False,
        creation_time="2023-06-27 17:05:27.153000",
    ),
    format_ram_permission(
        name="AWSRAMPermissionGlueAllTablesReadWriteForDatabase",
        resource_type="glue:Database",
        version="3",
        is_resource_type_default=False,
        creation_time="2022-10-27 14:19:43.345000",
    ),
    format_ram_permission(
        name="AWSRAMPermissionGlueDatabaseReadWrite",
        resource_type="glue:Database",
        version="3",
        is_resource_type_default=False,
        creation_time="2022-10-27 14:20:08.415000",
    ),
    format_ram_permission(
        name="AWSRAMPermissionGlueTableReadWriteForDatabase",
        resource_type="glue:Database",
        version="3",
        is_resource_type_default=False,
        creation_time="2022-10-27 14:19:17.062000",
    ),
    format_ram_permission(
        name="AWSRAMPermissionLFTagGlueDatabaseReadWrite",
        resource_type="glue:Database",
        is_resource_type_default=False,
        creation_time="2022-10-27 14:22:14.138000",
    ),
    format_ram_permission(
        name="AWSRAMPermissionLFTagGlueTableReadWriteForDatabase",
        resource_type="glue:Database",
        is_resource_type_default=False,
        creation_time="2022-10-27 14:21:49.920000",
    ),
    format_ram_permission(
        name="AWSRAMDefaultPermissionGlueTable",
        resource_type="glue:Table",
        creation_time="2022-06-30 17:04:24.372000",
    ),
    format_ram_permission(
        name="AWSRAMLFEnabledGlueDatabaseReadWriteForTable",
        resource_type="glue:Table",
        is_resource_type_default=False,
        creation_time="2023-06-27 17:05:33.534000",
    ),
    format_ram_permission(
        name="AWSRAMLFEnabledGlueTableReadWrite",
        resource_type="glue:Table",
        is_resource_type_default=False,
        creation_time="2023-06-27 17:05:21.095000",
    ),
    format_ram_permission(
        name="AWSRAMPermissionGlueDatabaseReadWriteForTable",
        resource_type="glue:Table",
        version="2",
        is_resource_type_default=False,
        creation_time="2022-06-30 17:04:09.034000",
    ),
    format_ram_permission(
        name="AWSRAMPermissionGlueTableReadWrite",
        resource_type="glue:Table",
        version="2",
        is_resource_type_default=False,
        creation_time="2022-06-30 17:04:03.947000",
    ),
    format_ram_permission(
        name="AWSRAMPermissionLFTagGlueDatabaseReadWriteForTable",
        resource_type="glue:Table",
        is_resource_type_default=False,
        creation_time="2022-10-27 14:22:26.191000",
    ),
    format_ram_permission(
        name="AWSRAMPermissionLFTagGlueTableReadWrite",
        resource_type="glue:Table",
        is_resource_type_default=False,
        creation_time="2022-10-27 14:22:01.908000",
    ),
    format_ram_permission(
        name="AWSRAMDefaultPermissionGlueCatalog",
        resource_type="glue:Catalog",
        creation_time="2022-06-30 17:04:26.612000",
    ),
    format_ram_permission(
        name="AWSRAMPermissionGlueAllTablesReadWriteForCatalog",
        resource_type="glue:Catalog",
        version="3",
        is_resource_type_default=False,
        creation_time="2022-10-27 14:19:56.250000",
    ),
    format_ram_permission(
        name="AWSRAMPermissionGlueDatabaseReadWriteForCatalog",
        resource_type="glue:Catalog",
        version="3",
        is_resource_type_default=False,
        creation_time="2022-10-27 14:20:20.322000",
    ),
    format_ram_permission(
        name="AWSRAMPermissionGlueTableReadWriteForCatalog",
        resource_type="glue:Catalog",
        version="3",
        is_resource_type_default=False,
        creation_time="2022-10-27 14:19:30.644000",
    ),
    format_ram_permission(
        name="AWSRAMPermissionLFTagGlueDatabaseReadWriteForCatalog",
        resource_type="glue:Catalog",
        is_resource_type_default=False,
        creation_time="2022-10-27 14:22:37.435000",
    ),
    format_ram_permission(
        name="AWSRAMPermissionLFTagGlueTableReadWriteForCatalog",
        resource_type="glue:Catalog",
        is_resource_type_default=False,
        creation_time="2022-10-27 14:21:37.593000",
    ),
]
