import base64
import calendar
import copy
import hashlib
import io
import json
import logging
import os
import re
import tarfile
import threading
import time
import warnings
import weakref
import zipfile
from collections import defaultdict
from datetime import datetime
from gzip import GzipFile
from sys import platform
from typing import Any, Dict, Iterable, List, Optional, Tuple, TypedDict, Union

import requests.exceptions

from moto import settings
from moto.awslambda.policy import Policy
from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel, CloudFormationModel
from moto.core.exceptions import RESTError
from moto.core.utils import iso_8601_datetime_with_nanoseconds, unix_time_millis, utcnow
from moto.dynamodb import dynamodb_backends
from moto.dynamodbstreams import dynamodbstreams_backends
from moto.ecr.exceptions import ImageNotFoundException
from moto.ecr.models import ecr_backends
from moto.iam.exceptions import IAMNotFoundException
from moto.iam.models import iam_backends
from moto.kinesis.models import KinesisBackend, kinesis_backends
from moto.logs.models import logs_backends
from moto.moto_api._internal import mock_random as random
from moto.s3.exceptions import MissingBucket, MissingKey
from moto.s3.models import FakeKey, s3_backends
from moto.sqs.models import sqs_backends
from moto.utilities.docker_utilities import DockerModel
from moto.utilities.utils import (
    ARN_PARTITION_REGEX,
    get_partition,
    load_resource_as_bytes,
)

from .exceptions import (
    ConflictException,
    CrossAccountNotAllowed,
    GenericResourcNotFound,
    InvalidParameterValueException,
    InvalidRoleFormat,
    UnknownAliasException,
    UnknownEventConfig,
    UnknownFunctionException,
    UnknownLayerException,
    UnknownLayerVersionException,
    ValidationException,
)
from .utils import (
    make_function_arn,
    make_function_ver_arn,
    make_layer_arn,
    make_layer_ver_arn,
    split_layer_arn,
)

logger = logging.getLogger(__name__)


class LayerDataType(TypedDict):
    Arn: str
    CodeSize: int


def zip2tar(zip_bytes: bytes) -> io.BytesIO:
    tarstream = io.BytesIO()
    timeshift = int((datetime.now() - utcnow()).total_seconds())
    tarf = tarfile.TarFile(fileobj=tarstream, mode="w")
    with zipfile.ZipFile(io.BytesIO(zip_bytes), "r") as zipf:
        for zipinfo in zipf.infolist():
            if zipinfo.is_dir():
                continue

            tarinfo = tarfile.TarInfo(name=zipinfo.filename)
            tarinfo.size = zipinfo.file_size
            tarinfo.mtime = calendar.timegm(zipinfo.date_time) - timeshift
            infile = zipf.open(zipinfo.filename)
            tarf.addfile(tarinfo, infile)

    tarstream.seek(0)
    return tarstream


def file2tar(file_content: bytes, file_name: str) -> io.BytesIO:
    tarstream = io.BytesIO()
    tarf = tarfile.TarFile(fileobj=tarstream, mode="w")
    tarinfo = tarfile.TarInfo(name=file_name)
    tarinfo.size = len(file_content)
    tarf.addfile(tarinfo, io.BytesIO(file_content))

    tarstream.seek(0)
    return tarstream


class _VolumeRefCount:
    __slots__ = "refcount", "volume"

    def __init__(self, refcount: int, volume: Any):
        self.refcount = refcount
        self.volume = volume


class _DockerDataVolumeContext:
    # {sha256: _VolumeRefCount}
    _data_vol_map: Dict[str, _VolumeRefCount] = defaultdict(
        lambda: _VolumeRefCount(0, None)
    )
    _lock = threading.Lock()

    def __init__(self, lambda_func: "LambdaFunction"):
        self._lambda_func = lambda_func
        self._vol_ref: Optional[_VolumeRefCount] = None

    @property
    def name(self) -> str:
        return self._vol_ref.volume.name  # type: ignore[union-attr]

    def __enter__(self) -> "_DockerDataVolumeContext":
        # See if volume is already known
        with self.__class__._lock:
            self._vol_ref = self.__class__._data_vol_map[self._lambda_func.code_digest]
            self._vol_ref.refcount += 1
            if self._vol_ref.refcount > 1:
                return self

            # See if the volume already exists
            for vol in self._lambda_func.docker_client.volumes.list():
                if vol.name == self._lambda_func.code_digest:
                    self._vol_ref.volume = vol
                    return self

            # It doesn't exist so we need to create it
            self._vol_ref.volume = self._lambda_func.docker_client.volumes.create(
                self._lambda_func.code_digest
            )
            volumes = {self.name: {"bind": settings.LAMBDA_DATA_DIR, "mode": "rw"}}

            self._lambda_func.ensure_image_exists("busybox")
            container = self._lambda_func.docker_client.containers.run(
                "busybox", "sleep 100", volumes=volumes, detach=True
            )
            try:
                with zip2tar(self._lambda_func.code_bytes) as stream:
                    container.put_archive(settings.LAMBDA_DATA_DIR, stream)
                if settings.is_test_proxy_mode():
                    ca_cert = load_resource_as_bytes(__name__, "../moto_proxy/ca.crt")
                    with file2tar(ca_cert, "ca.crt") as cert_stream:
                        container.put_archive(settings.LAMBDA_DATA_DIR, cert_stream)
            finally:
                container.remove(force=True)

        return self

    def __exit__(self, exc_type: Any, exc_val: Any, exc_tb: Any) -> None:
        with self.__class__._lock:
            self._vol_ref.refcount -= 1  # type: ignore[union-attr]
            if self._vol_ref.refcount == 0:  # type: ignore[union-attr]
                import docker.errors

                try:
                    self._vol_ref.volume.remove()  # type: ignore[union-attr]
                except docker.errors.APIError as e:
                    if e.status_code != 409:
                        raise

                    raise  # multiple processes trying to use same volume?


class _DockerDataVolumeLayerContext:
    _data_vol_map: Dict[str, _VolumeRefCount] = defaultdict(
        lambda: _VolumeRefCount(0, None)
    )
    _lock = threading.Lock()

    def __init__(self, lambda_func: "LambdaFunction"):
        self._lambda_func = lambda_func
        self._layers: List[LayerDataType] = self._lambda_func.layers
        self._vol_ref: Optional[_VolumeRefCount] = None

    @property
    def name(self) -> str:
        return self._vol_ref.volume.name  # type: ignore[union-attr]

    @property
    def hash(self) -> str:
        return "-".join(
            [
                layer["Arn"].split("layer:")[-1].replace(":", "_")
                for layer in self._layers
            ]
        )

    def __enter__(self) -> "_DockerDataVolumeLayerContext":
        # See if volume is already known
        with self.__class__._lock:
            self._vol_ref = self.__class__._data_vol_map[self.hash]
            self._vol_ref.refcount += 1
            if self._vol_ref.refcount > 1:
                return self

            # See if the volume already exists
            for vol in self._lambda_func.docker_client.volumes.list():
                if vol.name == self.hash:
                    self._vol_ref.volume = vol
                    return self

            # It doesn't exist so we need to create it
            self._vol_ref.volume = self._lambda_func.docker_client.volumes.create(
                self.hash
            )
            # If we don't have any layers to apply, just return at this point
            # When invoking the function, we will bind this empty volume
            if len(self._layers) == 0:
                return self
            volumes = {self.name: {"bind": "/opt", "mode": "rw"}}

            self._lambda_func.ensure_image_exists("busybox")
            container = self._lambda_func.docker_client.containers.run(
                "busybox", "sleep 100", volumes=volumes, detach=True
            )
            backend: "LambdaBackend" = lambda_backends[self._lambda_func.account_id][
                self._lambda_func.region
            ]
            try:
                for layer in self._layers:
                    try:
                        layer_zip = backend.layers_versions_by_arn(  # type: ignore[union-attr]
                            layer["Arn"]
                        ).code_bytes
                        layer_tar = zip2tar(layer_zip)
                        container.put_archive("/opt", layer_tar)
                    except zipfile.BadZipfile as e:
                        warnings.warn(f"Error extracting layer to Lambda: {e}")
            finally:
                container.remove(force=True)

        return self

    def __exit__(self, exc_type: Any, exc_val: Any, exc_tb: Any) -> None:
        with self.__class__._lock:
            self._vol_ref.refcount -= 1  # type: ignore[union-attr]
            if self._vol_ref.refcount == 0:  # type: ignore[union-attr]
                import docker.errors

                try:
                    self._vol_ref.volume.remove()  # type: ignore[union-attr]
                except docker.errors.APIError as e:
                    if e.status_code != 409:
                        raise

                    raise  # multiple processes trying to use same volume?


def _zipfile_content(zipfile_content: Union[str, bytes]) -> Tuple[bytes, int, str, str]:
    try:
        to_unzip_code = base64.b64decode(bytes(zipfile_content, "utf-8"))  # type: ignore[arg-type]
    except Exception:
        to_unzip_code = base64.b64decode(zipfile_content)

    sha_code = hashlib.sha256(to_unzip_code)
    base64ed_sha = base64.b64encode(sha_code.digest()).decode("utf-8")
    sha_hex_digest = sha_code.hexdigest()
    return to_unzip_code, len(to_unzip_code), base64ed_sha, sha_hex_digest


def _s3_content(key: Any) -> Tuple[bytes, int, str, str]:
    sha_code = hashlib.sha256(key.value)
    base64ed_sha = base64.b64encode(sha_code.digest()).decode("utf-8")
    sha_hex_digest = sha_code.hexdigest()
    return key.value, key.size, base64ed_sha, sha_hex_digest


def _validate_s3_bucket_and_key(
    account_id: str, partition: str, data: Dict[str, Any]
) -> Optional[FakeKey]:
    key = None
    try:
        # FIXME: does not validate bucket region
        key = s3_backends[account_id][partition].get_object(
            data["S3Bucket"], data["S3Key"]
        )
    except MissingBucket:
        if do_validate_s3():
            raise InvalidParameterValueException(
                "Error occurred while GetObject. S3 Error Code: NoSuchBucket. S3 Error Message: The specified bucket does not exist"
            )
    except MissingKey:
        if do_validate_s3():
            raise ValueError(
                "InvalidParameterValueException",
                "Error occurred while GetObject. S3 Error Code: NoSuchKey. S3 Error Message: The specified key does not exist.",
            )
    return key


class EventInvokeConfig:
    def __init__(self, arn: str, last_modified: str, config: Dict[str, Any]) -> None:
        self.config = config
        self.validate_max()
        self.validate()
        self.arn = arn
        self.last_modified = last_modified

    def validate_max(self) -> None:
        if "MaximumRetryAttempts" in self.config:
            mra = self.config["MaximumRetryAttempts"]
            if mra > 2:
                raise ValidationException(
                    mra,
                    "maximumRetryAttempts",
                    "Member must have value less than or equal to 2",
                )

            # < 0 validation done by botocore
        if "MaximumEventAgeInSeconds" in self.config:
            mra = self.config["MaximumEventAgeInSeconds"]
            if mra > 21600:
                raise ValidationException(
                    mra,
                    "maximumEventAgeInSeconds",
                    "Member must have value less than or equal to 21600",
                )

            # < 60 validation done by botocore

    def validate(self) -> None:
        # https://docs.aws.amazon.com/lambda/latest/dg/API_OnSuccess.html
        regex = r"^$|arn:(aws[a-zA-Z0-9-]*):([a-zA-Z0-9\-])+:([a-z]{2}(-gov)?-[a-z]+-\d{1})?:(\d{12})?:(.*)"
        pattern = re.compile(regex)

        if self.config["DestinationConfig"]:
            destination_config = self.config["DestinationConfig"]
            if (
                "OnSuccess" in destination_config
                and "Destination" in destination_config["OnSuccess"]
            ):
                contents = destination_config["OnSuccess"]["Destination"]
                if not pattern.match(contents):
                    raise ValidationException(
                        contents,
                        "destinationConfig.onSuccess.destination",
                        f"Member must satisfy regular expression pattern: {regex}",
                    )
            if (
                "OnFailure" in destination_config
                and "Destination" in destination_config["OnFailure"]
            ):
                contents = destination_config["OnFailure"]["Destination"]
                if not pattern.match(contents):
                    raise ValidationException(
                        contents,
                        "destinationConfig.onFailure.destination",
                        f"Member must satisfy regular expression pattern: {regex}",
                    )

    def response(self) -> Dict[str, Any]:
        response = {"FunctionArn": self.arn, "LastModified": self.last_modified}
        response.update(self.config)
        return response


class ImageConfig:
    def __init__(self, config: Dict[str, Any]) -> None:
        self.cmd = config.get("Command", [])
        self.entry_point = config.get("EntryPoint", [])
        self.working_directory = config.get("WorkingDirectory", None)

    def response(self) -> Dict[str, Any]:
        content = {
            "Command": self.cmd,
            "EntryPoint": self.entry_point,
        }
        if self.working_directory is not None:
            content["WorkingDirectory"] = self.working_directory
        return dict(content)


class Permission(CloudFormationModel):
    def __init__(self, region: str):
        self.region = region

    @staticmethod
    def cloudformation_name_type() -> str:
        return "Permission"

    @staticmethod
    def cloudformation_type() -> str:
        return "AWS::Lambda::Permission"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Dict[str, Any],
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "Permission":
        properties = cloudformation_json["Properties"]
        backend = lambda_backends[account_id][region_name]
        fn = backend.get_function(properties["FunctionName"])
        fn.policy.add_statement(raw=json.dumps(properties))
        return Permission(region=region_name)


class LayerVersion(CloudFormationModel):
    def __init__(self, spec: Dict[str, Any], account_id: str, region: str):
        # required
        self.account_id = account_id
        self.region = region
        self.name = spec["LayerName"]
        self.content = spec["Content"]

        # optional
        self.description = spec.get("Description", "")
        self.compatible_architectures = spec.get("CompatibleArchitectures", [])
        self.compatible_runtimes = spec.get("CompatibleRuntimes", [])
        self.license_info = spec.get("LicenseInfo", "")

        # auto-generated
        self.created_date = utcnow().strftime("%Y-%m-%d %H:%M:%S")
        self.version: Optional[int] = None
        self._attached = False
        self._layer: Optional["Layer"] = None

        if "ZipFile" in self.content:
            (
                self.code_bytes,
                self.code_size,
                self.code_sha_256,
                self.code_digest,
            ) = _zipfile_content(self.content["ZipFile"])
        else:
            partition = get_partition(region)
            key = _validate_s3_bucket_and_key(
                account_id, partition=partition, data=self.content
            )
            if key:
                (
                    self.code_bytes,
                    self.code_size,
                    self.code_sha_256,
                    self.code_digest,
                ) = _s3_content(key)
            else:
                self.code_bytes = b""
                self.code_size = 0
                self.code_sha_256 = ""
                self.code_digest = ""

    @property
    def arn(self) -> str:
        if self.version:
            return make_layer_ver_arn(
                self.region, self.account_id, self.name, self.version
            )
        raise ValueError("Layer version is not set")

    def attach(self, layer: "Layer", version: int) -> None:
        self._attached = True
        self._layer = layer
        self.version = version

    def get_layer_version(self) -> Dict[str, Any]:
        return {
            "Content": {
                "Location": "s3://",
                "CodeSha256": self.code_sha_256,
                "CodeSize": self.code_size,
            },
            "Version": self.version,
            "LayerArn": self._layer.layer_arn,  # type: ignore[union-attr]
            "LayerVersionArn": self.arn,
            "CreatedDate": self.created_date,
            "CompatibleArchitectures": self.compatible_architectures,
            "CompatibleRuntimes": self.compatible_runtimes,
            "Description": self.description,
            "LicenseInfo": self.license_info,
        }

    @staticmethod
    def cloudformation_name_type() -> str:
        return "LayerVersion"

    @staticmethod
    def cloudformation_type() -> str:
        return "AWS::Lambda::LayerVersion"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Dict[str, Any],
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "LayerVersion":
        properties = cloudformation_json["Properties"]
        optional_properties = ("Description", "CompatibleRuntimes", "LicenseInfo")

        # required
        spec = {
            "Content": properties["Content"],
            "LayerName": resource_name,
        }
        for prop in optional_properties:
            if prop in properties:
                spec[prop] = properties[prop]

        backend = lambda_backends[account_id][region_name]
        layer_version = backend.publish_layer_version(spec)
        return layer_version


class LambdaAlias(BaseModel):
    def __init__(
        self,
        account_id: str,
        region: str,
        name: str,
        function_name: str,
        function_version: str,
        description: str,
        routing_config: str,
    ):
        self.arn = f"arn:{get_partition(region)}:lambda:{region}:{account_id}:function:{function_name}:{name}"
        self.name = name
        self.function_version = function_version
        self.description = description
        self.routing_config = routing_config
        self.revision_id = str(random.uuid4())

    def update(
        self,
        description: Optional[str],
        function_version: Optional[str],
        routing_config: Optional[str],
    ) -> None:
        if description is not None:
            self.description = description
        if function_version is not None:
            self.function_version = function_version
        if routing_config is not None:
            self.routing_config = routing_config

    def to_json(self) -> Dict[str, Any]:
        return {
            "AliasArn": self.arn,
            "Description": self.description,
            "FunctionVersion": self.function_version,
            "Name": self.name,
            "RevisionId": self.revision_id,
            "RoutingConfig": self.routing_config or None,
        }


class Layer(object):
    def __init__(self, layer_version: LayerVersion):
        self.region = layer_version.region
        self.name = layer_version.name

        self.layer_arn = make_layer_arn(
            self.region, layer_version.account_id, self.name
        )
        self._latest_version = 0
        self.layer_versions: Dict[str, LayerVersion] = {}

    def attach_version(self, layer_version: LayerVersion) -> None:
        self._latest_version += 1
        layer_version.attach(self, self._latest_version)
        self.layer_versions[str(self._latest_version)] = layer_version

    def delete_version(self, layer_version: str) -> None:
        self.layer_versions.pop(str(layer_version), None)

    def to_dict(self) -> Dict[str, Any]:
        if not self.layer_versions:
            return {}

        last_key = sorted(self.layer_versions.keys(), key=lambda version: int(version))[
            -1
        ]
        return {
            "LayerName": self.name,
            "LayerArn": self.layer_arn,
            "LatestMatchingVersion": self.layer_versions[last_key].get_layer_version(),
        }


class LambdaFunction(CloudFormationModel, DockerModel):
    def __init__(
        self,
        account_id: str,
        spec: Dict[str, Any],
        region: str,
        version: Union[str, int] = 1,
    ):
        DockerModel.__init__(self)
        # required
        self.account_id = account_id
        self.region = region
        self.partition = get_partition(region)
        self.code = spec["Code"]
        self.function_name = spec["FunctionName"]
        self.handler = spec.get("Handler")
        self.role = spec["Role"]
        self.run_time = spec.get("Runtime")
        self.logs_backend = logs_backends[account_id][self.region]
        self.environment_vars = spec.get("Environment", {}).get("Variables", {})
        self.policy = Policy(self)
        self.url_config: Optional[FunctionUrlConfig] = None
        self.state = "Active"
        self.reserved_concurrency = spec.get("ReservedConcurrentExecutions", None)

        # optional
        self.ephemeral_storage: str
        self.code_digest: str
        self.code_bytes: bytes
        self.event_invoke_config: List[EventInvokeConfig] = []

        self.description = spec.get("Description", "")
        self.memory_size = spec.get("MemorySize", 128)
        self.package_type = spec.get("PackageType", "Zip")
        self.publish = spec.get("Publish", False)  # this is ignored currently
        self.timeout = spec.get("Timeout", 3)
        self.layers: List[LayerDataType] = self._get_layers_data(spec.get("Layers", []))
        self.signing_profile_version_arn = spec.get("SigningProfileVersionArn")
        self.signing_job_arn = spec.get("SigningJobArn")
        self.code_signing_config_arn = spec.get("CodeSigningConfigArn")
        self.tracing_config = spec.get("TracingConfig") or {"Mode": "PassThrough"}
        self.architectures: List[str] = spec.get("Architectures", ["x86_64"])
        self.image_config: ImageConfig = ImageConfig(spec.get("ImageConfig", {}))
        _es = spec.get("EphemeralStorage")
        if _es:
            self.ephemeral_storage = _es["Size"]
        else:
            self.ephemeral_storage = 512

        self.logs_group_name = f"/aws/lambda/{self.function_name}"

        # this isn't finished yet. it needs to find out the VpcId value
        self._vpc_config = spec.get(
            "VpcConfig", {"SubnetIds": [], "SecurityGroupIds": []}
        )

        # auto-generated
        self.version = version
        self.last_modified = iso_8601_datetime_with_nanoseconds()

        self._set_function_code(self.code)

        self.function_arn: str = make_function_arn(
            self.region, self.account_id, self.function_name
        )

        self.tags = spec.get("Tags") or dict()

    def __getstate__(self) -> Dict[str, Any]:
        return {
            k: v
            for (k, v) in self.__dict__.items()
            if k != "_DockerModel__docker_client"
        }

    def set_version(self, version: int) -> None:
        self.function_arn = make_function_ver_arn(
            self.region, self.account_id, self.function_name, version
        )
        self.version = version
        self.last_modified = iso_8601_datetime_with_nanoseconds()

    @property
    def architectures(self) -> List[str]:
        return self._architectures

    @architectures.setter
    def architectures(self, architectures: List[str]) -> None:
        if (
            len(architectures) > 1
            or not architectures
            or architectures[0] not in ("x86_64", "arm64")
        ):
            raise ValidationException(
                str(architectures),
                "architectures",
                "Member must satisfy constraint: "
                "[Member must satisfy enum value set: [x86_64, arm64], Member must not be null]",
            )
        self._architectures = architectures

    @property
    def ephemeral_storage(self) -> int:
        return self._ephemeral_storage

    @ephemeral_storage.setter
    def ephemeral_storage(self, ephemeral_storage: int) -> None:
        if ephemeral_storage > 10240:
            raise ValidationException(
                str(ephemeral_storage),
                "ephemeralStorage.size",
                "Member must have value less than or equal to 10240",
            )

        # ephemeral_storage < 512 is handled by botocore 1.30.0
        self._ephemeral_storage = ephemeral_storage

    @property
    def vpc_config(self) -> Dict[str, Any]:  # type: ignore[misc]
        config = self._vpc_config.copy()
        if config["SecurityGroupIds"]:
            config.update({"VpcId": "vpc-123abc"})
        return config

    @property
    def physical_resource_id(self) -> str:
        return self.function_name

    def __repr__(self) -> str:
        return json.dumps(self.get_configuration())

    def _get_layers_data(self, layers_versions_arns: List[str]) -> List[LayerDataType]:
        backend = lambda_backends[self.account_id][self.region]
        layer_versions = [
            backend.layers_versions_by_arn(layer_version)
            for layer_version in layers_versions_arns
        ]
        if not all(layer_versions):
            raise UnknownLayerVersionException(layers_versions_arns)
        # The `if lv` part is not necessary - we know there are no None's, because of the `all()`-check earlier
        # But MyPy does not seem to understand this
        return [
            {"Arn": lv.arn, "CodeSize": lv.code_size} for lv in layer_versions if lv
        ]

    def get_function_code_signing_config(self) -> Dict[str, Any]:
        return {
            "CodeSigningConfigArn": self.code_signing_config_arn,
            "FunctionName": self.function_name,
        }

    def get_configuration(self, on_create: bool = False) -> Dict[str, Any]:
        config = {
            "CodeSha256": self.code_sha_256,
            "CodeSize": self.code_size,
            "Description": self.description,
            "FunctionArn": self.function_arn,
            "FunctionName": self.function_name,
            "Handler": self.handler,
            "LastModified": self.last_modified,
            "MemorySize": self.memory_size,
            "Role": self.role,
            "Runtime": self.run_time,
            "State": self.state,
            "PackageType": self.package_type,
            "Timeout": self.timeout,
            "Version": str(self.version),
            "VpcConfig": self.vpc_config,
            "Layers": self.layers,
            "SigningProfileVersionArn": self.signing_profile_version_arn,
            "SigningJobArn": self.signing_job_arn,
            "TracingConfig": self.tracing_config,
            "Architectures": self.architectures,
            "EphemeralStorage": {
                "Size": self.ephemeral_storage,
            },
            "SnapStart": {"ApplyOn": "None", "OptimizationStatus": "Off"},
        }
        if self.package_type == "Image":
            config["ImageConfigResponse"] = {
                "ImageConfig": self.image_config.response(),
            }
        if not on_create:
            # Only return this variable after the first creation
            config["LastUpdateStatus"] = "Successful"
        if self.environment_vars:
            config["Environment"] = {"Variables": self.environment_vars}

        return config

    def get_code(self) -> Dict[str, Any]:
        resp = {"Configuration": self.get_configuration()}
        if "S3Key" in self.code:
            resp["Code"] = {
                "Location": f"s3://awslambda-{self.region}-tasks.s3-{self.region}.amazonaws.com/{self.code['S3Key']}",
                "RepositoryType": "S3",
            }
        elif "ImageUri" in self.code:
            resp["Code"] = {
                "RepositoryType": "ECR",
                "ImageUri": self.code.get("ImageUri"),
                "ResolvedImageUri": self.code.get("ImageUri").split(":")[0]
                + "@sha256:"
                + self.code_sha_256,
            }
        if self.tags:
            resp["Tags"] = self.tags
        if self.reserved_concurrency:
            resp.update(
                {
                    "Concurrency": {
                        "ReservedConcurrentExecutions": self.reserved_concurrency
                    }
                }
            )
        return resp

    def update_configuration(self, config_updates: Dict[str, Any]) -> Dict[str, Any]:
        for key, value in config_updates.items():
            if key == "Description":
                self.description = value
            elif key == "Handler":
                self.handler = value
            elif key == "MemorySize":
                self.memory_size = value
            elif key == "Role":
                self.role = value
            elif key == "Runtime":
                self.run_time = value
            elif key == "Timeout":
                self.timeout = value
            elif key == "VpcConfig":
                self._vpc_config = value
            elif key == "Environment":
                self.environment_vars = value["Variables"]
            elif key == "Layers":
                self.layers = self._get_layers_data(value)

        return self.get_configuration()

    def _set_function_code(self, updated_spec: Dict[str, Any]) -> None:
        from_update = updated_spec is not self.code

        # "DryRun" is only used for UpdateFunctionCode
        if from_update and "DryRun" in updated_spec and updated_spec["DryRun"]:
            return

        if "ZipFile" in updated_spec:
            if from_update:
                self.code["ZipFile"] = updated_spec["ZipFile"]

            (
                self.code_bytes,
                self.code_size,
                self.code_sha_256,
                self.code_digest,
            ) = _zipfile_content(updated_spec["ZipFile"])

            # TODO: we should be putting this in a lambda bucket
            self.code["UUID"] = str(random.uuid4())
            self.code["S3Key"] = f"{self.function_name}-{self.code['UUID']}"
        elif "S3Bucket" in updated_spec and "S3Key" in updated_spec:
            key = None
            try:
                if from_update:
                    # FIXME: does not validate bucket region
                    key = s3_backends[self.account_id][self.partition].get_object(
                        updated_spec["S3Bucket"], updated_spec["S3Key"]
                    )
                else:
                    key = _validate_s3_bucket_and_key(
                        self.account_id, partition=self.partition, data=self.code
                    )
            except MissingBucket:
                if do_validate_s3():
                    raise ValueError(
                        "InvalidParameterValueException",
                        "Error occurred while GetObject. S3 Error Code: NoSuchBucket. S3 Error Message: The specified bucket does not exist",
                    )
            except MissingKey:
                if do_validate_s3():
                    raise ValueError(
                        "InvalidParameterValueException",
                        "Error occurred while GetObject. S3 Error Code: NoSuchKey. S3 Error Message: The specified key does not exist.",
                    )
            if key:
                (
                    self.code_bytes,
                    self.code_size,
                    self.code_sha_256,
                    self.code_digest,
                ) = _s3_content(key)
            else:
                self.code_bytes = b""
                self.code_size = 0
                self.code_sha_256 = ""
            if from_update:
                self.code["S3Bucket"] = updated_spec["S3Bucket"]
                self.code["S3Key"] = updated_spec["S3Key"]
        elif "ImageUri" in updated_spec:
            if settings.lambda_stub_ecr():
                self.code_sha_256 = hashlib.sha256(
                    updated_spec["ImageUri"].encode("utf-8")
                ).hexdigest()
                self.code_size = 0
            else:
                if "@" in updated_spec["ImageUri"]:
                    # deploying via digest
                    uri, digest = updated_spec["ImageUri"].split("@")
                    image_id = {"imageDigest": digest}
                else:
                    # deploying via tag
                    uri, tag = updated_spec["ImageUri"].split(":")
                    image_id = {"imageTag": tag}

                repo_name = uri.split("/")[-1]
                ecr_backend = ecr_backends[self.account_id][self.region]
                registry_id = ecr_backend.describe_registry()["registryId"]
                images = ecr_backend.batch_get_image(
                    repository_name=repo_name, image_ids=[image_id]
                )["images"]

                if len(images) == 0:
                    raise ImageNotFoundException(image_id, repo_name, registry_id)  # type: ignore
                else:
                    manifest = json.loads(images[0]["imageManifest"])
                    self.code_sha_256 = images[0]["imageId"]["imageDigest"].replace(
                        "sha256:", ""
                    )
                    self.code_size = manifest["config"]["size"]
            if from_update:
                self.code["ImageUri"] = updated_spec["ImageUri"]

    def update_function_code(self, updated_spec: Dict[str, Any]) -> Dict[str, Any]:
        self._set_function_code(updated_spec)
        return self.get_configuration()

    @staticmethod
    def convert(s: Any) -> str:  # type: ignore[misc]
        try:
            return str(s, encoding="utf-8")
        except Exception:
            return s

    def _invoke_lambda(self, event: Optional[str] = None) -> Tuple[str, bool, str]:
        import docker
        import docker.errors

        # Create the LogGroup if necessary, to write the result to
        self.logs_backend.ensure_log_group(self.logs_group_name)
        # TODO: context not yet implemented
        if event is None:
            event = dict()  # type: ignore[assignment]
        output = None

        try:
            # TODO: I believe we can keep the container running and feed events as needed
            #       also need to hook it up to the other services so it can make kws/s3 etc calls
            #  Should get invoke_id /RequestId from invocation
            env_vars = {
                "_HANDLER": self.handler,
                "AWS_EXECUTION_ENV": f"AWS_Lambda_{self.run_time}",
                "AWS_LAMBDA_FUNCTION_TIMEOUT": self.timeout,
                "AWS_LAMBDA_FUNCTION_NAME": self.function_name,
                "AWS_LAMBDA_FUNCTION_MEMORY_SIZE": self.memory_size,
                "AWS_LAMBDA_FUNCTION_VERSION": self.version,
                "AWS_REGION": self.region,
                "AWS_ACCESS_KEY_ID": "role-account-id",
                "AWS_SECRET_ACCESS_KEY": "role-secret-key",
                "AWS_SESSION_TOKEN": "session-token",
            }

            env_vars.update(self.environment_vars)
            env_vars["MOTO_HOST"] = settings.moto_server_host()
            moto_port = settings.moto_server_port()
            env_vars["MOTO_PORT"] = moto_port
            env_vars["MOTO_HTTP_ENDPOINT"] = f'{env_vars["MOTO_HOST"]}:{moto_port}'

            if settings.is_test_proxy_mode():
                env_vars["HTTPS_PROXY"] = env_vars["MOTO_HTTP_ENDPOINT"]
                env_vars["AWS_CA_BUNDLE"] = "/var/task/ca.crt"

            container = exit_code = None
            log_config = docker.types.LogConfig(type=docker.types.LogConfig.types.JSON)

            with _DockerDataVolumeContext(
                self
            ) as data_vol, _DockerDataVolumeLayerContext(self) as layer_context:
                try:
                    run_kwargs: Dict[str, Any] = dict()
                    network_name = settings.moto_network_name()
                    network_mode = settings.moto_network_mode()
                    if network_name:
                        run_kwargs["network"] = network_name
                    elif network_mode:
                        run_kwargs["network_mode"] = network_mode
                    elif settings.TEST_SERVER_MODE:
                        # AWSLambda can make HTTP requests to a Docker container called 'motoserver'
                        # Only works if our Docker-container is named 'motoserver'
                        # TODO: should remove this and rely on 'network_mode' instead, as this is too tightly coupled with our own test setup
                        run_kwargs["links"] = {"motoserver": "motoserver"}

                    # add host.docker.internal host on linux to emulate Mac + Windows behavior
                    #   for communication with other mock AWS services running on localhost
                    if platform == "linux" or platform == "linux2":
                        run_kwargs["extra_hosts"] = {
                            "host.docker.internal": "host-gateway"
                        }

                    # Change entry point if requested
                    if self.image_config.entry_point:
                        run_kwargs["entrypoint"] = self.image_config.entry_point

                    # The requested image can be found in one of a few repos:
                    # - User-provided repo
                    # - mlupin/docker-lambda (the repo with up-to-date AWSLambda images
                    # - lambci/lambda (the repo with older/outdated AWSLambda images
                    #
                    # We'll cycle through all of them - when we find the repo that contains our image, we use it
                    image_repos = set(
                        [
                            settings.moto_lambda_image(),
                            "mlupin/docker-lambda",
                            "lambci/lambda",
                        ]
                    )
                    for image_repo in image_repos:
                        image_ref = f"{image_repo}:{self.run_time}"
                        try:
                            self.ensure_image_exists(image_ref)
                            break
                        except docker.errors.NotFound:
                            pass
                    volumes = {
                        data_vol.name: {"bind": "/var/task", "mode": "rw"},
                        layer_context.name: {"bind": "/opt", "mode": "rw"},
                    }
                    container = self.docker_client.containers.run(
                        image_ref,
                        [self.handler, json.dumps(event)],
                        remove=False,
                        mem_limit=f"{self.memory_size}m",
                        volumes=volumes,
                        environment=env_vars,
                        detach=True,
                        log_config=log_config,
                        **run_kwargs,
                    )
                finally:
                    if container:
                        try:
                            exit_code = container.wait(timeout=300)["StatusCode"]
                        except requests.exceptions.ReadTimeout:
                            exit_code = -1
                            container.stop()
                            container.kill()

                        output = container.logs(stdout=False, stderr=True)
                        output += container.logs(stdout=True, stderr=False)
                        container.remove()

            output = output.decode("utf-8")  # type: ignore[union-attr]

            self.save_logs(output)

            # We only care about the response from the lambda
            # Which is the last line of the output, according to https://github.com/lambci/docker-lambda/issues/25
            resp = output.splitlines()[-1]
            logs = os.linesep.join(
                [line for line in self.convert(output).splitlines()[:-1]]
            )
            invocation_error = exit_code != 0
            return resp, invocation_error, logs
        except docker.errors.DockerException as e:
            # Docker itself is probably not running - there will be no Lambda-logs to handle
            msg = f"error running docker: {e}"
            logger.error(msg)
            self.save_logs(msg)
            return msg, True, ""

    def save_logs(self, output: str) -> None:
        # Send output to "logs" backend
        invoke_id = random.uuid4().hex
        date = utcnow()
        log_stream_name = (
            f"{date.year}/{date.month:02d}/{date.day:02d}/[{self.version}]{invoke_id}"
        )
        self.logs_backend.create_log_stream(self.logs_group_name, log_stream_name)
        log_events = [
            {"timestamp": unix_time_millis(), "message": line}
            for line in output.splitlines()
        ]
        self.logs_backend.put_log_events(
            self.logs_group_name, log_stream_name, log_events
        )

    def invoke(
        self, body: str, request_headers: Any, response_headers: Any
    ) -> Union[str, bytes]:
        if body:
            body = json.loads(body)
        else:
            body = "{}"

        # Get the invocation type:
        res, errored, logs = self._invoke_lambda(event=body)
        if errored:
            response_headers["x-amz-function-error"] = "Handled"

        inv_type = request_headers.get("x-amz-invocation-type", "RequestResponse")
        if inv_type == "RequestResponse":
            encoded = base64.b64encode(logs.encode("utf-8"))
            response_headers["x-amz-log-result"] = encoded.decode("utf-8")
            return res.encode("utf-8")
        else:
            return res

    @staticmethod
    def cloudformation_name_type() -> str:
        return "FunctionName"

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-lambda-function.html
        return "AWS::Lambda::Function"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Dict[str, Any],
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "LambdaFunction":
        properties = cloudformation_json["Properties"]
        optional_properties = (
            "Description",
            "MemorySize",
            "Publish",
            "Timeout",
            "VpcConfig",
            "Environment",
            "ReservedConcurrentExecutions",
        )

        # required
        spec = {
            "Code": properties["Code"],
            "FunctionName": resource_name,
            "Handler": properties["Handler"],
            "Role": properties["Role"],
            "Runtime": properties["Runtime"],
        }

        # NOTE: Not doing `properties.get(k, DEFAULT)` to avoid duplicating the
        # default logic
        for prop in optional_properties:
            if prop in properties:
                spec[prop] = properties[prop]

        # when ZipFile is present in CloudFormation, per the official docs,
        # the code it's a plaintext code snippet up to 4096 bytes.
        # this snippet converts this plaintext code to a proper base64-encoded ZIP file.
        if "ZipFile" in properties["Code"]:
            spec["Code"]["ZipFile"] = base64.b64encode(
                cls._create_zipfile_from_plaintext_code(spec["Code"]["ZipFile"])
            )

        backend = lambda_backends[account_id][region_name]
        fn = backend.create_function(spec)
        return fn

    @classmethod
    def has_cfn_attr(cls, attr: str) -> bool:
        return attr in ["Arn"]

    def get_cfn_attribute(self, attribute_name: str) -> str:
        from moto.cloudformation.exceptions import UnformattedGetAttTemplateException

        if attribute_name == "Arn":
            return make_function_arn(self.region, self.account_id, self.function_name)
        raise UnformattedGetAttTemplateException()

    @classmethod
    def update_from_cloudformation_json(  # type: ignore[misc]
        cls,
        original_resource: "LambdaFunction",
        new_resource_name: str,
        cloudformation_json: Dict[str, Any],
        account_id: str,
        region_name: str,
    ) -> "LambdaFunction":
        updated_props = cloudformation_json["Properties"]
        original_resource.update_configuration(updated_props)
        original_resource.update_function_code(updated_props["Code"])
        return original_resource

    @staticmethod
    def _create_zipfile_from_plaintext_code(code: str) -> bytes:
        zip_output = io.BytesIO()
        zip_file = zipfile.ZipFile(zip_output, "w", zipfile.ZIP_DEFLATED)
        zip_file.writestr("index.py", code)
        # This should really be part of the 'lambci' docker image
        from moto.packages.cfnresponse import cfnresponse

        with open(cfnresponse.__file__) as cfn:
            zip_file.writestr("cfnresponse.py", cfn.read())
        zip_file.close()
        zip_output.seek(0)
        return zip_output.read()

    def delete(self, account_id: str, region: str) -> None:
        lambda_backends[account_id][region].delete_function(self.function_name)

    def create_url_config(self, config: Dict[str, Any]) -> "FunctionUrlConfig":
        self.url_config = FunctionUrlConfig(function=self, config=config)
        return self.url_config

    def delete_url_config(self) -> None:
        self.url_config = None

    def get_url_config(self) -> "FunctionUrlConfig":
        if not self.url_config:
            raise GenericResourcNotFound()
        return self.url_config

    def update_url_config(self, config: Dict[str, Any]) -> "FunctionUrlConfig":
        self.url_config.update(config)  # type: ignore[union-attr]
        return self.url_config  # type: ignore[return-value]


class FunctionUrlConfig:
    def __init__(self, function: LambdaFunction, config: Dict[str, Any]):
        self.function = function
        self.config = config
        self.url = f"https://{random.uuid4().hex}.lambda-url.{function.region}.on.aws"
        self.created = utcnow().strftime("%Y-%m-%dT%H:%M:%S.000+0000")
        self.last_modified = self.created

    def to_dict(self) -> Dict[str, Any]:
        return {
            "FunctionUrl": self.url,
            "FunctionArn": self.function.function_arn,
            "AuthType": self.config.get("AuthType"),
            "Cors": self.config.get("Cors"),
            "CreationTime": self.created,
            "LastModifiedTime": self.last_modified,
            "InvokeMode": self.config.get("InvokeMode") or "Buffered",
        }

    def update(self, new_config: Dict[str, Any]) -> None:
        if new_config.get("Cors"):
            self.config["Cors"] = new_config["Cors"]
        if new_config.get("AuthType"):
            self.config["AuthType"] = new_config["AuthType"]
        self.last_modified = utcnow().strftime("%Y-%m-%dT%H:%M:%S")


class EventSourceMapping(CloudFormationModel):
    def __init__(self, spec: Dict[str, Any]):
        # required
        self.function_name = spec["FunctionName"]
        self.event_source_arn = spec["EventSourceArn"]

        # optional
        self.batch_size = spec.get("BatchSize")  # type: ignore[assignment]
        self.starting_position = spec.get("StartingPosition", "TRIM_HORIZON")
        self.enabled = spec.get("Enabled", True)
        self.starting_position_timestamp = spec.get("StartingPositionTimestamp", None)

        self.function_arn: str = spec["FunctionArn"]
        self.uuid = str(random.uuid4())
        self.last_modified = time.mktime(utcnow().timetuple())

    def _get_service_source_from_arn(self, event_source_arn: str) -> str:
        return event_source_arn.split(":")[2].lower()

    def _validate_event_source(self, event_source_arn: str) -> bool:
        valid_services = ("dynamodb", "kinesis", "sqs")
        service = self._get_service_source_from_arn(event_source_arn)
        return service in valid_services

    @property
    def event_source_arn(self) -> str:
        return self._event_source_arn

    @event_source_arn.setter
    def event_source_arn(self, event_source_arn: str) -> None:
        if not self._validate_event_source(event_source_arn):
            raise ValueError(
                "InvalidParameterValueException", "Unsupported event source type"
            )
        self._event_source_arn = event_source_arn

    @property
    def batch_size(self) -> int:
        return self._batch_size

    @batch_size.setter
    def batch_size(self, batch_size: Optional[int]) -> None:
        batch_size_service_map = {
            "kinesis": (100, 10000),
            "dynamodb": (100, 1000),
            "sqs": (10, 10),
        }

        source_type = self._get_service_source_from_arn(self.event_source_arn)
        batch_size_for_source = batch_size_service_map[source_type]

        if batch_size is None:
            self._batch_size = batch_size_for_source[0]
        elif batch_size > batch_size_for_source[1]:
            error_message = (
                f"BatchSize {batch_size} exceeds the max of {batch_size_for_source[1]}"
            )
            raise ValueError("InvalidParameterValueException", error_message)
        else:
            self._batch_size = int(batch_size)

    def get_configuration(self) -> Dict[str, Any]:
        return {
            "UUID": self.uuid,
            "BatchSize": self.batch_size,
            "EventSourceArn": self.event_source_arn,
            "FunctionArn": self.function_arn,
            "LastModified": self.last_modified,
            "LastProcessingResult": None,
            "State": "Enabled" if self.enabled else "Disabled",
            "StateTransitionReason": "User initiated",
            "StartingPosition": self.starting_position,
        }

    def delete(self, account_id: str, region_name: str) -> None:
        lambda_backend = lambda_backends[account_id][region_name]
        lambda_backend.delete_event_source_mapping(self.uuid)

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-lambda-eventsourcemapping.html
        return "AWS::Lambda::EventSourceMapping"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Dict[str, Any],
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "EventSourceMapping":
        properties = cloudformation_json["Properties"]
        lambda_backend = lambda_backends[account_id][region_name]
        return lambda_backend.create_event_source_mapping(properties)

    @classmethod
    def update_from_cloudformation_json(  # type: ignore[misc]
        cls,
        original_resource: Any,
        new_resource_name: str,
        cloudformation_json: Dict[str, Any],
        account_id: str,
        region_name: str,
    ) -> "EventSourceMapping":
        properties = cloudformation_json["Properties"]
        event_source_uuid = original_resource.uuid
        lambda_backend = lambda_backends[account_id][region_name]
        return lambda_backend.update_event_source_mapping(event_source_uuid, properties)  # type: ignore[return-value]

    @classmethod
    def delete_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Dict[str, Any],
        account_id: str,
        region_name: str,
    ) -> None:
        properties = cloudformation_json["Properties"]
        lambda_backend = lambda_backends[account_id][region_name]
        esms = lambda_backend.list_event_source_mappings(
            event_source_arn=properties["EventSourceArn"],
            function_name=properties["FunctionName"],
        )

        for esm in esms:
            if esm.uuid == resource_name:
                esm.delete(account_id, region_name)

    @property
    def physical_resource_id(self) -> str:
        return self.uuid


class LambdaVersion(CloudFormationModel):
    def __init__(self, spec: Dict[str, Any]):
        self.version = spec["Version"]

    def __repr__(self) -> str:
        return str(self.logical_resource_id)  # type: ignore[attr-defined]

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-lambda-version.html
        return "AWS::Lambda::Version"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Dict[str, Any],
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "LambdaVersion":
        properties = cloudformation_json["Properties"]
        function_name = properties["FunctionName"]
        func = lambda_backends[account_id][region_name].publish_version(function_name)
        spec = {"Version": func.version}  # type: ignore[union-attr]
        return LambdaVersion(spec)


class LambdaStorage(object):
    def __init__(self, region_name: str, account_id: str):
        # Format 'func_name' {'versions': []}
        self._functions: Dict[str, Any] = {}
        self._arns: weakref.WeakValueDictionary[str, LambdaFunction] = (
            weakref.WeakValueDictionary()
        )
        self.region_name = region_name
        self.account_id = account_id
        self.partition = get_partition(region_name)

        # function-arn -> alias -> LambdaAlias
        self._aliases: Dict[str, Dict[str, LambdaAlias]] = defaultdict(lambda: {})

    def _get_latest(self, name: str) -> LambdaFunction:
        return self._functions[name]["latest"]

    def _get_version(self, name: str, version: str) -> Optional[LambdaFunction]:
        for config in self._functions[name]["versions"]:
            if str(config.version) == version:
                return config

        return None

    def _get_function_aliases(self, function_name: str) -> Dict[str, LambdaAlias]:
        fn = self.get_function_by_name_or_arn_with_qualifier(function_name)
        # Split ARN to retrieve an ARN without a qualifier present
        [arn, _, _] = self.split_function_arn(fn.function_arn)
        return self._aliases[arn]

    def delete_alias(self, name: str, function_name: str) -> None:
        aliases = self._get_function_aliases(function_name)
        aliases.pop(name, None)

    def get_alias(self, name: str, function_name: str) -> LambdaAlias:
        aliases = self._get_function_aliases(function_name)
        if name in aliases:
            return aliases[name]

        arn = f"arn:{get_partition(self.region_name)}:lambda:{self.region_name}:{self.account_id}:function:{function_name}:{name}"
        raise UnknownAliasException(arn)

    def put_alias(
        self,
        name: str,
        function_name: str,
        function_version: str,
        description: str,
        routing_config: str,
    ) -> LambdaAlias:
        fn = self.get_function_by_name_or_arn_with_qualifier(
            function_name, function_version
        )
        aliases = self._get_function_aliases(function_name)
        if name in aliases:
            arn = f"arn:{get_partition(self.region_name)}:lambda:{self.region_name}:{self.account_id}:function:{function_name}:{name}"
            raise ConflictException(f"Alias already exists: {arn}")

        alias = LambdaAlias(
            account_id=self.account_id,
            region=self.region_name,
            name=name,
            function_name=fn.function_name,
            function_version=function_version,
            description=description,
            routing_config=routing_config,
        )
        aliases[name] = alias
        return alias

    def update_alias(
        self,
        name: str,
        function_name: str,
        function_version: str,
        description: str,
        routing_config: str,
    ) -> LambdaAlias:
        alias = self.get_alias(name, function_name)

        # errors if new function version doesn't exist
        self.get_function_by_name_or_arn_with_qualifier(function_name, function_version)

        alias.update(description, function_version, routing_config)
        return alias

    def get_function_by_name_forbid_qualifier(self, name: str) -> LambdaFunction:
        """
        Get function by name forbidding a qualifier
        :raises: UnknownFunctionException if function not found
        :raises: InvalidParameterValue if qualifier is provided
        """

        if name.count(":") == 1:
            raise InvalidParameterValueException("Cannot provide qualifier")

        if name not in self._functions:
            raise self.construct_unknown_function_exception(name)

        return self._get_latest(name)

    def get_function_by_name_with_qualifier(
        self, name: str, qualifier: Optional[str] = None
    ) -> LambdaFunction:
        """
        Get function by name with an optional qualifier
        :raises: UnknownFunctionException if function not found
        """

        # Function name may contain an alias
        # <fn_name>:<alias_name>
        if ":" in name:
            # Prefer qualifier in name over qualifier arg
            [name, qualifier] = name.split(":")

        # Find without qualifier
        if qualifier is None:
            return self.get_function_by_name_forbid_qualifier(name)

        if name not in self._functions:
            raise self.construct_unknown_function_exception(name, qualifier)

        # Find by latest
        if qualifier.lower() == "$latest":
            return self._functions[name]["latest"]

        # Find by version
        found_version = self._get_version(name, qualifier)
        if found_version:
            return found_version

        # Find by alias
        aliases = self._get_function_aliases(name)
        if qualifier in aliases:
            alias_version = aliases[qualifier].function_version

            # Find by alias pointing to latest
            if alias_version.lower() == "$latest":
                return self._functions[name]["latest"]

            # Find by alias pointing to version
            found_alias = self._get_version(name, alias_version)
            if found_alias:
                return found_alias

        raise self.construct_unknown_function_exception(name, qualifier)

    def list_versions_by_function(self, name: str) -> Iterable[LambdaFunction]:
        if name not in self._functions:
            return []

        latest = copy.copy(self._functions[name]["latest"])
        latest.function_arn += ":$LATEST"
        return [latest] + self._functions[name]["versions"]

    def list_aliases(self, function_name: str) -> Iterable[LambdaAlias]:
        aliases = self._get_function_aliases(function_name)
        return sorted(aliases.values(), key=lambda alias: alias.name)

    def get_arn(self, arn: str) -> Optional[LambdaFunction]:
        [arn_without_qualifier, _, _] = self.split_function_arn(arn)
        return self._arns.get(arn_without_qualifier, None)

    def split_function_arn(self, arn: str) -> Tuple[str, str, Optional[str]]:
        """
        Handy utility to parse an ARN into:
        - ARN without qualifier
        - Function name
        - Optional qualifier
        """
        qualifier = None
        # Function ARN may contain an alias
        # arn:aws:lambda:region:account_id:function:<fn_name>:<alias_name>
        if ":" in arn.split(":function:")[-1]:
            qualifier = arn.split(":")[-1]
            # arn = arn:aws:lambda:region:account_id:function:<fn_name>
            arn = ":".join(arn.split(":")[0:-1])
        name = arn.split(":")[-1]
        return arn, name, qualifier

    def get_function_by_name_or_arn_forbid_qualifier(
        self, name_or_arn: str
    ) -> LambdaFunction:
        """
        Get function by name or arn forbidding a qualifier
        :raises: UnknownFunctionException if function not found
        :raises: InvalidParameterValue if qualifier is provided
        """

        if re.match(ARN_PARTITION_REGEX, name_or_arn):
            [_, name, qualifier] = self.split_function_arn(name_or_arn)

            if qualifier is not None:
                raise InvalidParameterValueException("Cannot provide qualifier")

            return self.get_function_by_name_forbid_qualifier(name)
        else:
            # name_or_arn is not an arn
            return self.get_function_by_name_forbid_qualifier(name_or_arn)

    def get_function_by_name_or_arn_with_qualifier(
        self, name_or_arn: str, qualifier: Optional[str] = None
    ) -> LambdaFunction:
        """
        Get function by name or arn with an optional qualifier
        :raises: UnknownFunctionException if function not found
        """

        if re.match(ARN_PARTITION_REGEX, name_or_arn):
            [_, name, qualifier_in_arn] = self.split_function_arn(name_or_arn)
            return self.get_function_by_name_with_qualifier(
                name, qualifier_in_arn or qualifier
            )
        else:
            return self.get_function_by_name_with_qualifier(name_or_arn, qualifier)

    def construct_unknown_function_exception(
        self, name_or_arn: str, qualifier: Optional[str] = None
    ) -> UnknownFunctionException:
        if re.match(ARN_PARTITION_REGEX, name_or_arn):
            arn = name_or_arn
        else:
            # name_or_arn is a function name with optional qualifier <func_name>[:<qualifier>]
            arn = make_function_arn(self.region_name, self.account_id, name_or_arn)
            # Append explicit qualifier to arn only if the name doesn't already have it
            if qualifier and ":" not in name_or_arn:
                arn = f"{arn}:{qualifier}"
        return UnknownFunctionException(arn)

    def put_function(self, fn: LambdaFunction) -> None:
        valid_role = re.match(InvalidRoleFormat.pattern, fn.role)
        if valid_role:
            account = valid_role.group(2)
            if account != self.account_id:
                raise CrossAccountNotAllowed()
            try:
                iam_backend = iam_backends[self.account_id][self.partition]
                iam_backend.get_role_by_arn(fn.role)
            except IAMNotFoundException:
                raise InvalidParameterValueException(
                    "The role defined for the function cannot be assumed by Lambda."
                )
        else:
            raise InvalidRoleFormat(fn.role)
        if fn.function_name in self._functions:
            self._functions[fn.function_name]["latest"] = fn
        else:
            self._functions[fn.function_name] = {"latest": fn, "versions": []}
        self._arns[fn.function_arn] = fn

    def publish_version(
        self, name_or_arn: str, description: str = ""
    ) -> Optional[LambdaFunction]:
        function = self.get_function_by_name_or_arn_forbid_qualifier(name_or_arn)
        name = function.function_name
        if name not in self._functions:
            return None
        if not self._functions[name]["latest"]:
            return None

        all_versions = self._functions[name]["versions"]
        if all_versions:
            latest_published = all_versions[-1]
            if latest_published.code_sha_256 == function.code_sha_256:
                # Nothing has changed, don't publish
                return latest_published

        new_version = len(all_versions) + 1
        fn = copy.copy(self._functions[name]["latest"])
        fn.set_version(new_version)
        if description:
            fn.description = description

        self._functions[name]["versions"].append(fn)
        self._arns[fn.function_arn] = fn
        return fn

    def del_function(self, name_or_arn: str, qualifier: Optional[str] = None) -> None:
        # Qualifier may be explicitly passed or part of function name or ARN, extract it here
        if re.match(ARN_PARTITION_REGEX, name_or_arn):
            # Extract from ARN
            if ":" in name_or_arn.split(":function:")[-1]:
                qualifier = name_or_arn.split(":")[-1]
        else:
            # Extract from function name
            if ":" in name_or_arn:
                qualifier = name_or_arn.split(":")[1]

        function = self.get_function_by_name_or_arn_with_qualifier(
            name_or_arn, qualifier
        )
        name = function.function_name
        if not qualifier:
            # Something is still reffing this so delete all arns
            latest = self._functions[name]["latest"].function_arn
            del self._arns[latest]

            for fn in self._functions[name]["versions"]:
                del self._arns[fn.function_arn]

            del self._functions[name]

        else:
            if qualifier == "$LATEST":
                self._functions[name]["latest"] = None
            else:
                self._functions[name]["versions"].remove(function)

            # If theres no functions left
            if (
                not self._functions[name]["versions"]
                and not self._functions[name]["latest"]
            ):
                del self._functions[name]

        self._aliases[function.function_arn] = {}

    def all(self) -> Iterable[LambdaFunction]:
        result = []

        for function_group in list(self._functions.values()):
            latest = copy.deepcopy(function_group["latest"])
            latest.function_arn = f"{latest.function_arn}:$LATEST"
            result.append(latest)

            result.extend(function_group["versions"])

        return result

    def latest(self) -> Iterable[LambdaFunction]:
        """
        Return the list of functions with version @LATEST
        :return:
        """
        result = []
        for function_group in self._functions.values():
            if function_group["latest"] is not None:
                result.append(function_group["latest"])

        return result


class LayerStorage(object):
    def __init__(self) -> None:
        self._layers: Dict[str, Layer] = {}
        self._arns: weakref.WeakValueDictionary[str, LambdaFunction] = (
            weakref.WeakValueDictionary()
        )

    def _find_layer_by_name_or_arn(self, name_or_arn: str) -> Layer:
        if name_or_arn in self._layers:
            return self._layers[name_or_arn]
        for layer in self._layers.values():
            if layer.layer_arn == name_or_arn:
                return layer
        raise UnknownLayerException()

    def put_layer_version(self, layer_version: LayerVersion) -> None:
        """
        :param layer_version: LayerVersion
        """
        if layer_version.name not in self._layers:
            self._layers[layer_version.name] = Layer(layer_version)
        self._layers[layer_version.name].attach_version(layer_version)

    def list_layers(self) -> Iterable[Dict[str, Any]]:
        return [
            layer.to_dict() for layer in self._layers.values() if layer.layer_versions
        ]

    def delete_layer_version(self, layer_name: str, layer_version: str) -> None:
        layer = self._find_layer_by_name_or_arn(layer_name)
        layer.delete_version(layer_version)

    def get_layer_version(self, layer_name: str, layer_version: str) -> LayerVersion:
        layer = self._find_layer_by_name_or_arn(layer_name)
        for lv in layer.layer_versions.values():
            if lv.version == int(layer_version):
                return lv
        raise UnknownLayerException()

    def get_layer_versions(self, layer_name: str) -> List[LayerVersion]:
        if layer_name in self._layers:
            return list(iter(self._layers[layer_name].layer_versions.values()))
        return []

    def get_layer_version_by_arn(
        self, layer_version_arn: str
    ) -> Optional[LayerVersion]:
        split_arn = split_layer_arn(layer_version_arn)
        if split_arn.layer_name in self._layers:
            return self._layers[split_arn.layer_name].layer_versions.get(
                split_arn.version, None
            )
        return None


class LambdaBackend(BaseBackend):
    """
    Implementation of the AWS Lambda endpoint.
    Invoking functions is supported - they will run inside a Docker container, emulating the real AWS behaviour as closely as possible.

    .. warning:: When invoking a function using the decorators, the created Docker container cannot reach Moto (or it's in-memory state). Any AWS SDK-requests within your Lambda will try to connect to AWS instead.

    It is possible to connect from AWS Lambdas to other services, as long as you are running MotoProxy or the MotoServer in a Docker container.

    When running the MotoProxy, calls to other AWS services are automatically proxied.

    When running MotoServer, the Lambda has access to environment variables `MOTO_HOST` and `MOTO_PORT`, which can be used to build the url that MotoServer runs on:

    .. sourcecode:: python

        def lambda_handler(event, context):
            host = os.environ.get("MOTO_HOST")
            port = os.environ.get("MOTO_PORT")
            url = host + ":" + port
            ec2 = boto3.client('ec2', region_name='us-west-2', endpoint_url=url)

            # Or even simpler:
            full_url = os.environ.get("MOTO_HTTP_ENDPOINT")
            ec2 = boto3.client("ec2", region_name="eu-west-1", endpoint_url=full_url)

            ec2.do_whatever_inside_the_existing_moto_server()

    Moto will run on port 5000 by default. This can be overwritten by setting an environment variable when starting Moto:

    .. sourcecode:: bash

        # This env var will be propagated to the Docker container running the Lambda functions
        MOTO_PORT=5000 moto_server

    The Docker container uses the default network mode, `bridge`.
    The following environment variables are available for fine-grained control over the Docker connection options:

    .. sourcecode:: bash

        # Provide the name of a custom network to connect to
        MOTO_DOCKER_NETWORK_NAME=mycustomnetwork moto_server

        # Override the network mode
        # For example, network_mode=host would use the network of the host machine
        # Note that this option will be ignored if MOTO_DOCKER_NETWORK_NAME is also set
        MOTO_DOCKER_NETWORK_MODE=host moto_server

    _-_-_-_

    The Docker images used by Moto are taken from the following repositories:

    - `mlupin/docker-lambda` (for recent versions)
    - `lambci/lambda` (for older/outdated versions)

    Use the following environment variable to configure Moto to look for images in an additional repository:

    .. sourcecode:: bash

        MOTO_DOCKER_LAMBDA_IMAGE=mlupin/docker-lambda

    _-_-_-_

    Use the following environment variable if you want to configure the data directory used by the Docker containers:

    .. sourcecode:: bash

        MOTO_LAMBDA_DATA_DIR=/tmp/data

    _-_-_-_

    If you want to mock the Lambda-containers invocation without starting a Docker-container, use the simple decorator:

    .. sourcecode:: python

        @mock_aws(config={"lambda": {"use_docker": False}})

    """

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self._lambdas = LambdaStorage(region_name=region_name, account_id=account_id)
        self._event_source_mappings: Dict[str, EventSourceMapping] = {}
        self._layers = LayerStorage()

    def create_alias(
        self,
        name: str,
        function_name: str,
        function_version: str,
        description: str,
        routing_config: str,
    ) -> LambdaAlias:
        return self._lambdas.put_alias(
            name, function_name, function_version, description, routing_config
        )

    def delete_alias(self, name: str, function_name: str) -> None:
        return self._lambdas.delete_alias(name, function_name)

    def get_alias(self, name: str, function_name: str) -> LambdaAlias:
        return self._lambdas.get_alias(name, function_name)

    def update_alias(
        self,
        name: str,
        function_name: str,
        function_version: str,
        description: str,
        routing_config: str,
    ) -> LambdaAlias:
        """
        The RevisionId parameter is not yet implemented
        """
        return self._lambdas.update_alias(
            name, function_name, function_version, description, routing_config
        )

    def create_function(self, spec: Dict[str, Any]) -> LambdaFunction:
        """
        The Code.ImageUri is not validated by default. Set environment variable MOTO_LAMBDA_STUB_ECR=false if you want to validate the image exists in our mocked ECR.
        """
        function_name = spec.get("FunctionName", None)
        if function_name is None:
            raise RESTError("InvalidParameterValueException", "Missing FunctionName")

        fn = LambdaFunction(
            account_id=self.account_id,
            spec=spec,
            region=self.region_name,
            version="$LATEST",
        )

        self._lambdas.put_function(fn)

        if spec.get("Publish"):
            ver = self.publish_version(function_name)
            fn = copy.deepcopy(
                fn
            )  # We don't want to change the actual version - just the return value
            fn.version = ver.version  # type: ignore[union-attr]
        return fn

    def create_function_url_config(
        self, name_or_arn: str, config: Dict[str, Any]
    ) -> FunctionUrlConfig:
        """
        The Qualifier-parameter is not yet implemented.
        Function URLs are not yet mocked, so invoking them will fail
        """
        function = self._lambdas.get_function_by_name_or_arn_forbid_qualifier(
            name_or_arn
        )
        return function.create_url_config(config)

    def delete_function_url_config(self, name_or_arn: str) -> None:
        """
        The Qualifier-parameter is not yet implemented
        """
        function = self._lambdas.get_function_by_name_or_arn_forbid_qualifier(
            name_or_arn
        )
        function.delete_url_config()

    def get_function_url_config(self, name_or_arn: str) -> FunctionUrlConfig:
        """
        The Qualifier-parameter is not yet implemented
        """
        function = self._lambdas.get_function_by_name_or_arn_forbid_qualifier(
            name_or_arn
        )
        if not function:
            raise UnknownFunctionException(arn=name_or_arn)
        return function.get_url_config()

    def update_function_url_config(
        self, name_or_arn: str, config: Dict[str, Any]
    ) -> FunctionUrlConfig:
        """
        The Qualifier-parameter is not yet implemented
        """
        function = self._lambdas.get_function_by_name_or_arn_forbid_qualifier(
            name_or_arn
        )
        return function.update_url_config(config)

    def create_event_source_mapping(self, spec: Dict[str, Any]) -> EventSourceMapping:
        required = ["EventSourceArn", "FunctionName"]
        for param in required:
            if not spec.get(param):
                raise RESTError("InvalidParameterValueException", f"Missing {param}")

        # Validate function name
        func = self._lambdas.get_function_by_name_or_arn_with_qualifier(
            spec.get("FunctionName", "")
        )

        # Validate queue
        sqs_backend = sqs_backends[self.account_id][self.region_name]
        for queue in sqs_backend.queues.values():
            if queue.queue_arn == spec["EventSourceArn"]:
                if queue.lambda_event_source_mappings.get("func.function_arn"):
                    # TODO: Correct exception?
                    raise RESTError(
                        "ResourceConflictException", "The resource already exists."
                    )
                spec.update({"FunctionArn": func.function_arn})
                esm = EventSourceMapping(spec)
                self._event_source_mappings[esm.uuid] = esm

                # Set backend function on queue
                queue.lambda_event_source_mappings[esm.function_arn] = esm

                return esm

        ddbstream_backend = dynamodbstreams_backends[self.account_id][self.region_name]
        ddb_backend = dynamodb_backends[self.account_id][self.region_name]
        for stream in json.loads(ddbstream_backend.list_streams())["Streams"]:
            if stream["StreamArn"] == spec["EventSourceArn"]:
                spec.update({"FunctionArn": func.function_arn})
                esm = EventSourceMapping(spec)
                self._event_source_mappings[esm.uuid] = esm
                table_name = stream["TableName"]
                table = ddb_backend.get_table(table_name)
                table.lambda_event_source_mappings[esm.function_arn] = esm
                return esm

        kinesis_backend: KinesisBackend = kinesis_backends[self.account_id][
            self.region_name
        ]
        for stream in kinesis_backend.streams.values():
            if stream.arn == spec["EventSourceArn"]:
                spec.update({"FunctionArn": func.function_arn})
                esm = EventSourceMapping(spec)
                self._event_source_mappings[esm.uuid] = esm
                stream.lambda_event_source_mappings[esm.event_source_arn] = esm
                return esm

        raise RESTError("ResourceNotFoundException", "Invalid EventSourceArn")

    def publish_layer_version(self, spec: Dict[str, Any]) -> LayerVersion:
        required = ["LayerName", "Content"]
        for param in required:
            if not spec.get(param):
                raise InvalidParameterValueException(f"Missing {param}")
        layer_version = LayerVersion(
            spec, account_id=self.account_id, region=self.region_name
        )
        self._layers.put_layer_version(layer_version)
        return layer_version

    def list_layers(self) -> Iterable[Dict[str, Any]]:
        return self._layers.list_layers()

    def delete_layer_version(self, layer_name: str, layer_version: str) -> None:
        return self._layers.delete_layer_version(layer_name, layer_version)

    def get_layer_version(self, layer_name: str, layer_version: str) -> LayerVersion:
        return self._layers.get_layer_version(layer_name, layer_version)

    def list_layer_versions(self, layer_name: str) -> Iterable[LayerVersion]:
        return self._layers.get_layer_versions(layer_name)

    def layers_versions_by_arn(self, layer_version_arn: str) -> Optional[LayerVersion]:
        return self._layers.get_layer_version_by_arn(layer_version_arn)

    def publish_version(
        self, function_name: str, description: str = ""
    ) -> Optional[LambdaFunction]:
        return self._lambdas.publish_version(function_name, description)

    def get_function(
        self, function_name_or_arn: str, qualifier: Optional[str] = None
    ) -> LambdaFunction:
        return self._lambdas.get_function_by_name_or_arn_with_qualifier(
            function_name_or_arn, qualifier
        )

    def list_versions_by_function(self, function_name: str) -> Iterable[LambdaFunction]:
        return self._lambdas.list_versions_by_function(function_name)

    def list_aliases(self, function_name: str) -> Iterable[LambdaAlias]:
        return self._lambdas.list_aliases(function_name)

    def get_event_source_mapping(self, uuid: str) -> Optional[EventSourceMapping]:
        return self._event_source_mappings.get(uuid)

    def delete_event_source_mapping(self, uuid: str) -> Optional[EventSourceMapping]:
        return self._event_source_mappings.pop(uuid, None)

    def update_event_source_mapping(
        self, uuid: str, spec: Dict[str, Any]
    ) -> Optional[EventSourceMapping]:
        esm = self.get_event_source_mapping(uuid)
        if not esm:
            return None

        for key in spec.keys():
            if key == "FunctionName":
                func = self._lambdas.get_function_by_name_or_arn_with_qualifier(
                    spec[key]
                )
                esm.function_arn = func.function_arn
            elif key == "BatchSize":
                esm.batch_size = spec[key]
            elif key == "Enabled":
                esm.enabled = spec[key]

        esm.last_modified = time.mktime(utcnow().timetuple())
        return esm

    def list_event_source_mappings(
        self, event_source_arn: str, function_name: str
    ) -> Iterable[EventSourceMapping]:
        esms = list(self._event_source_mappings.values())
        if event_source_arn:
            esms = list(filter(lambda x: x.event_source_arn == event_source_arn, esms))
        if function_name:
            esms = list(filter(lambda x: x.function_name == function_name, esms))
        return esms

    def get_function_by_arn(self, function_arn: str) -> Optional[LambdaFunction]:
        return self._lambdas.get_arn(function_arn)

    def delete_function(
        self, function_name: str, qualifier: Optional[str] = None
    ) -> None:
        self._lambdas.del_function(function_name, qualifier)

    def list_functions(
        self, func_version: Optional[str] = None
    ) -> Iterable[LambdaFunction]:
        if func_version == "ALL":
            return self._lambdas.all()
        return self._lambdas.latest()

    def send_sqs_batch(self, function_arn: str, messages: Any, queue_arn: str) -> bool:
        success = True
        for message in messages:
            result = self._send_sqs_message(function_arn, message, queue_arn)
            if not result:
                success = False
        return success

    def _send_sqs_message(
        self, function_arn: str, message: Any, queue_arn: str
    ) -> bool:
        event = {
            "Records": [
                {
                    "messageId": message.id,
                    "receiptHandle": message.receipt_handle,
                    "body": message.body,
                    "attributes": {
                        "ApproximateReceiveCount": "1",
                        "SentTimestamp": "1545082649183",
                        "SenderId": "AIDAIENQZJOLO23YVJ4VO",
                        "ApproximateFirstReceiveTimestamp": "1545082649185",
                    },
                    "messageAttributes": {},
                    "md5OfBody": "098f6bcd4621d373cade4e832627b4f6",
                    "eventSource": "aws:sqs",
                    "eventSourceARN": queue_arn,
                    "awsRegion": self.region_name,
                }
            ]
        }
        if queue_arn.endswith(".fifo"):
            # Messages from FIFO queue have additional attributes
            event["Records"][0]["attributes"].update(
                {
                    "MessageGroupId": message.group_id,
                    "MessageDeduplicationId": message.deduplication_id,
                }
            )

        request_headers: Dict[str, Any] = {}
        response_headers: Dict[str, Any] = {}
        self.invoke(
            function_name=function_arn,
            qualifier=None,
            body=json.dumps(event),
            headers=request_headers,
            response_headers=response_headers,
        )
        return "x-amz-function-error" not in response_headers

    def send_kinesis_message(
        self,
        function_name: str,
        kinesis_stream: str,
        kinesis_partition_key: str,
        kinesis_sequence_number: str,
        kinesis_data: str,
        kinesis_shard_id: str,
    ) -> None:
        func = self._lambdas.get_function_by_name_or_arn_with_qualifier(
            function_name, qualifier=None
        )
        event = {
            "Records": [
                {
                    "kinesis": {
                        "kinesisSchemaVersion": "1.0",
                        "partitionKey": kinesis_partition_key,
                        "sequenceNumber": kinesis_sequence_number,
                        "data": kinesis_data,
                        "approximateArrivalTimestamp": round(time.time(), 3),
                    },
                    "eventSource": "aws:kinesis",
                    "eventVersion": "1.0",
                    "eventID": f"{kinesis_shard_id}:{kinesis_sequence_number}",
                    "eventName": "aws:kinesis:record",
                    "invokeIdentityArn": func.role,
                    "awsRegion": self.region_name,
                    "eventSourceARN": kinesis_stream,
                }
            ]
        }
        func.invoke(json.dumps(event), {}, {})

    def send_sns_message(
        self,
        function_name: str,
        message: str,
        subject: Optional[str] = None,
        qualifier: Optional[str] = None,
    ) -> None:
        event = {
            "Records": [
                {
                    "EventVersion": "1.0",
                    "EventSubscriptionArn": "arn:aws:sns:EXAMPLE",
                    "EventSource": "aws:sns",
                    "Sns": {
                        "SignatureVersion": "1",
                        "Timestamp": "1970-01-01T00:00:00.000Z",
                        "Signature": "EXAMPLE",
                        "SigningCertUrl": "EXAMPLE",
                        "MessageId": "95df01b4-ee98-5cb9-9903-4c221d41eb5e",
                        "Message": message,
                        "MessageAttributes": {
                            "Test": {"Type": "String", "Value": "TestString"},
                            "TestBinary": {"Type": "Binary", "Value": "TestBinary"},
                        },
                        "Type": "Notification",
                        "UnsubscribeUrl": "EXAMPLE",
                        "TopicArn": "arn:aws:sns:EXAMPLE",
                        "Subject": subject or "TestInvoke",
                    },
                }
            ]
        }
        func = self._lambdas.get_function_by_name_or_arn_with_qualifier(
            function_name, qualifier
        )
        func.invoke(json.dumps(event), {}, {})

    def send_dynamodb_items(
        self, function_arn: str, items: List[Any], source: str
    ) -> Union[str, bytes]:
        event = {
            "Records": [
                {
                    "eventID": item.to_json()["eventID"],
                    "eventName": "INSERT",
                    "eventVersion": item.to_json()["eventVersion"],
                    "eventSource": item.to_json()["eventSource"],
                    "awsRegion": self.region_name,
                    "dynamodb": item.to_json()["dynamodb"],
                    "eventSourceARN": source,
                }
                for item in items
            ]
        }
        func = self._lambdas.get_arn(function_arn)
        return func.invoke(json.dumps(event), {}, {})  # type: ignore[union-attr]

    def send_log_event(
        self,
        function_arn: str,
        filter_name: str,
        log_group_name: str,
        log_stream_name: str,
        log_events: Any,
    ) -> None:
        data = {
            "messageType": "DATA_MESSAGE",
            "owner": self.account_id,
            "logGroup": log_group_name,
            "logStream": log_stream_name,
            "subscriptionFilters": [filter_name],
            "logEvents": log_events,
        }

        output = io.BytesIO()
        with GzipFile(fileobj=output, mode="w") as f:
            f.write(json.dumps(data, separators=(",", ":")).encode("utf-8"))
        payload_gz_encoded = base64.b64encode(output.getvalue()).decode("utf-8")

        event = {"awslogs": {"data": payload_gz_encoded}}

        func = self._lambdas.get_arn(function_arn)
        func.invoke(json.dumps(event), {}, {})  # type: ignore[union-attr]

    def list_tags(self, resource: str) -> Dict[str, str]:
        return self._lambdas.get_function_by_name_or_arn_with_qualifier(resource).tags

    def tag_resource(self, resource: str, tags: Dict[str, str]) -> None:
        fn = self._lambdas.get_function_by_name_or_arn_with_qualifier(resource)
        fn.tags.update(tags)

    def untag_resource(self, resource: str, tagKeys: List[str]) -> None:
        fn = self._lambdas.get_function_by_name_or_arn_with_qualifier(resource)
        for key in tagKeys:
            fn.tags.pop(key, None)

    def add_permission(
        self, function_name: str, qualifier: str, raw: str
    ) -> Dict[str, Any]:
        fn = self.get_function(function_name, qualifier)
        return fn.policy.add_statement(raw, qualifier)

    def remove_permission(
        self, function_name: str, sid: str, revision: str = ""
    ) -> None:
        fn = self.get_function(function_name)
        fn.policy.del_statement(sid, revision)

    def get_function_code_signing_config(self, function_name: str) -> Dict[str, Any]:
        fn = self.get_function(function_name)
        return fn.get_function_code_signing_config()

    def get_policy(self, function_name: str, qualifier: Optional[str] = None) -> str:
        fn = self._lambdas.get_function_by_name_or_arn_with_qualifier(
            function_name, qualifier
        )
        return fn.policy.wire_format()

    def update_function_code(
        self, function_name: str, qualifier: str, body: Dict[str, Any]
    ) -> Optional[Dict[str, Any]]:
        fn: LambdaFunction = self.get_function(function_name, qualifier)
        fn.update_function_code(body)

        if body.get("Publish", False):
            fn = self.publish_version(function_name)  # type: ignore[assignment]

        return fn.update_function_code(body)

    def update_function_configuration(
        self, function_name: str, qualifier: str, body: Dict[str, Any]
    ) -> Optional[Dict[str, Any]]:
        fn = self.get_function(function_name, qualifier)

        return fn.update_configuration(body)

    def invoke(
        self,
        function_name: str,
        qualifier: Optional[str],
        body: Any,
        headers: Any,
        response_headers: Any,
    ) -> Optional[Union[str, bytes]]:
        """
        Invoking a Function with PackageType=Image is not yet supported.

        Invoking a Funcation against Lambda without docker now supports customised responses, the default being `Simple Lambda happy path OK`.
        You can use a dedicated API to override this, by configuring a queue of expected results.

        A request to `invoke` will take the first result from that queue.

        Configure this queue by making an HTTP request to `/moto-api/static/lambda-simple/response`. An example invocation looks like this:

        .. sourcecode:: python

            expected_results = {"results": ["test", "test 2"], "region": "us-east-1"}
            resp = requests.post(
                "http://motoapi.amazonaws.com/moto-api/static/lambda-simple/response",
                json=expected_results
            )
            assert resp.status_code == 201

            client = boto3.client("lambda", region_name="us-east-1")
            resp = client.invoke(...) # resp["Payload"].read().decode() == "test"
            resp = client.invoke(...) # resp["Payload"].read().decode() == "test2"
        """
        fn = self.get_function(function_name, qualifier)
        payload = fn.invoke(body, headers, response_headers)
        response_headers["Content-Length"] = str(len(payload))
        return payload

    def put_function_concurrency(
        self, function_name: str, reserved_concurrency: str
    ) -> str:
        """Establish concurrency limit/reservations for a function

        Actual lambda restricts concurrency to 1000 (default) per region/account
        across all functions; we approximate that behavior by summing across all
        functions (hopefully all in the same account and region) and allowing the
        caller to simulate an increased quota.

        By default, no quota is enforced in order to preserve compatibility with
        existing code that assumes it can do as many things as it likes. To model
        actual AWS behavior, define the MOTO_LAMBDA_CONCURRENCY_QUOTA environment
        variable prior to testing.
        """

        quota: Optional[str] = os.environ.get("MOTO_LAMBDA_CONCURRENCY_QUOTA")
        if quota is not None:
            # Enforce concurrency limits as described above
            available = int(quota) - int(reserved_concurrency)
            for fnx in self.list_functions():
                if fnx.reserved_concurrency and fnx.function_name != function_name:
                    available -= int(fnx.reserved_concurrency)
            if available < 100:
                raise InvalidParameterValueException(
                    "Specified ReservedConcurrentExecutions for function decreases account's UnreservedConcurrentExecution below its minimum value of [100]."
                )

        fn = self.get_function(function_name)
        fn.reserved_concurrency = reserved_concurrency
        return fn.reserved_concurrency

    def delete_function_concurrency(self, function_name: str) -> Optional[str]:
        fn = self.get_function(function_name)
        fn.reserved_concurrency = None
        return fn.reserved_concurrency

    def get_function_concurrency(self, function_name: str) -> str:
        fn = self.get_function(function_name)
        return fn.reserved_concurrency

    def put_function_event_invoke_config(
        self, function_name: str, config: Dict[str, Any]
    ) -> Dict[str, Any]:
        fn = self.get_function(function_name)
        event_config = EventInvokeConfig(fn.function_arn, fn.last_modified, config)
        fn.event_invoke_config.append(event_config)
        return event_config.response()

    def update_function_event_invoke_config(
        self, function_name: str, config: Dict[str, Any]
    ) -> Dict[str, Any]:
        # partial functionality, the update function just does a put
        # instead of partial update
        return self.put_function_event_invoke_config(function_name, config)

    def get_function_event_invoke_config(self, function_name: str) -> Dict[str, Any]:
        fn = self.get_function(function_name)
        if fn.event_invoke_config:
            response = fn.event_invoke_config[0]
            return response.response()
        else:
            raise UnknownEventConfig(fn.function_arn)

    def delete_function_event_invoke_config(self, function_name: str) -> None:
        if self.get_function_event_invoke_config(function_name):
            fn = self.get_function(function_name)
            fn.event_invoke_config = []

    def list_function_event_invoke_configs(self, function_name: str) -> Dict[str, Any]:
        response: Dict[str, List[Dict[str, Any]]] = {"FunctionEventInvokeConfigs": []}
        try:
            response["FunctionEventInvokeConfigs"] = [
                self.get_function_event_invoke_config(function_name)
            ]
            return response
        except UnknownEventConfig:
            return response


def do_validate_s3() -> bool:
    return os.environ.get("VALIDATE_LAMBDA_S3", "") in ["", "1", "true"]


lambda_backends = BackendDict(LambdaBackend, "lambda")
