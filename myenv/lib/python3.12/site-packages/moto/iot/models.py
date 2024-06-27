import hashlib
import re
import time
from collections import OrderedDict
from datetime import datetime, timedelta
from typing import TYPE_CHECKING, Any, Dict, Iterable, List, Optional, Pattern, Tuple

from cryptography import x509
from cryptography.hazmat._oid import NameOID
from cryptography.hazmat.backends import default_backend
from cryptography.hazmat.primitives import hashes, serialization
from cryptography.hazmat.primitives.asymmetric import rsa

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import utcnow
from moto.moto_api._internal import mock_random as random
from moto.utilities.paginator import paginate
from moto.utilities.utils import get_partition

from .exceptions import (
    CertificateStateException,
    DeleteConflictException,
    InvalidRequestException,
    InvalidStateTransitionException,
    ResourceAlreadyExistsException,
    ResourceNotFoundException,
    ThingStillAttached,
    VersionConflictException,
    VersionsLimitExceededException,
)
from .utils import PAGINATION_MODEL

if TYPE_CHECKING:
    from moto.iotdata.models import FakeShadow


class FakeThing(BaseModel):
    def __init__(
        self,
        thing_name: str,
        thing_type: Optional["FakeThingType"],
        attributes: Dict[str, Any],
        account_id: str,
        region_name: str,
    ):
        self.region_name = region_name
        self.thing_id = str(random.uuid4())
        self.thing_name = thing_name
        self.thing_type = thing_type
        self.attributes = attributes
        self.arn = f"arn:{get_partition(region_name)}:iot:{region_name}:{account_id}:thing/{thing_name}"
        self.version = 1
        # TODO: we need to handle "version"?

        # for iot-data
        self.thing_shadows: Dict[Optional[str], FakeShadow] = {}

    def matches(self, query_string: str) -> bool:
        if query_string == "*":
            return True
        if query_string.startswith("thingTypeName:"):
            if not self.thing_type:
                return False
            qs = query_string[14:].replace("*", ".*").replace("?", ".")
            return re.search(f"^{qs}$", self.thing_type.thing_type_name) is not None
        if query_string.startswith("thingName:"):
            qs = query_string[10:].replace("*", ".*").replace("?", ".")
            return re.search(f"^{qs}$", self.thing_name) is not None
        if query_string.startswith("attributes."):
            k, v = query_string[11:].split(":")
            return self.attributes.get(k) == v
        return query_string in self.thing_name

    def to_dict(
        self,
        include_default_client_id: bool = False,
        include_connectivity: bool = False,
        include_thing_id: bool = False,
    ) -> Dict[str, Any]:
        obj = {
            "thingName": self.thing_name,
            "thingArn": self.arn,
            "attributes": self.attributes,
            "version": self.version,
        }
        if self.thing_type:
            obj["thingTypeName"] = self.thing_type.thing_type_name
        if include_default_client_id:
            obj["defaultClientId"] = self.thing_name
        if include_connectivity:
            obj["connectivity"] = {
                "connected": True,
                "timestamp": time.mktime(utcnow().timetuple()),
            }
        if include_thing_id:
            obj["thingId"] = self.thing_id
        return obj


class FakeThingType(BaseModel):
    def __init__(
        self,
        thing_type_name: str,
        thing_type_properties: Optional[Dict[str, Any]],
        account_id: str,
        region_name: str,
    ):
        self.account_id = account_id
        self.region_name = region_name
        self.thing_type_name = thing_type_name
        self.thing_type_properties = thing_type_properties
        self.thing_type_id = str(random.uuid4())  # I don't know the rule of id
        t = time.time()
        self.metadata = {"deprecated": False, "creationDate": int(t * 1000) / 1000.0}
        self.arn = f"arn:{get_partition(self.region_name)}:iot:{self.region_name}:{self.account_id}:thingtype/{thing_type_name}"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "thingTypeName": self.thing_type_name,
            "thingTypeId": self.thing_type_id,
            "thingTypeProperties": self.thing_type_properties,
            "thingTypeMetadata": self.metadata,
            "thingTypeArn": self.arn,
        }


class FakeThingGroup(BaseModel):
    def __init__(
        self,
        thing_group_name: str,
        parent_group_name: str,
        thing_group_properties: Dict[str, str],
        account_id: str,
        region_name: str,
        thing_groups: Dict[str, "FakeThingGroup"],
    ):
        self.account_id = account_id
        self.region_name = region_name
        self.thing_group_name = thing_group_name
        self.thing_group_id = str(random.uuid4())  # I don't know the rule of id
        self.version = 1  # TODO: tmp
        self.parent_group_name = parent_group_name
        self.thing_group_properties = thing_group_properties or {}
        t = time.time()
        self.metadata: Dict[str, Any] = {"creationDate": int(t * 1000) / 1000.0}
        if parent_group_name:
            self.metadata["parentGroupName"] = parent_group_name
            # initilize rootToParentThingGroups
            if "rootToParentThingGroups" not in self.metadata:
                self.metadata["rootToParentThingGroups"] = []
            # search for parent arn
            for thing_group in thing_groups.values():
                if thing_group.thing_group_name == parent_group_name:
                    parent_thing_group_structure = thing_group
                    break
            # if parent arn found (should always be found)
            if parent_thing_group_structure:
                # copy parent's rootToParentThingGroups
                if "rootToParentThingGroups" in parent_thing_group_structure.metadata:
                    self.metadata["rootToParentThingGroups"].extend(
                        parent_thing_group_structure.metadata["rootToParentThingGroups"]
                    )
                self.metadata["rootToParentThingGroups"].extend(
                    [
                        {
                            "groupName": parent_group_name,
                            "groupArn": parent_thing_group_structure.arn,  # type: ignore
                        }
                    ]
                )
        self.arn = f"arn:{get_partition(self.region_name)}:iot:{self.region_name}:{self.account_id}:thinggroup/{thing_group_name}"
        self.things: Dict[str, FakeThing] = OrderedDict()

    def to_dict(self) -> Dict[str, Any]:
        return {
            "thingGroupName": self.thing_group_name,
            "thingGroupId": self.thing_group_id,
            "version": self.version,
            "thingGroupProperties": self.thing_group_properties,
            "thingGroupMetadata": self.metadata,
            "thingGroupArn": self.arn,
        }


class FakeCertificate(BaseModel):
    def __init__(
        self,
        certificate_pem: str,
        status: str,
        account_id: str,
        region_name: str,
        ca_certificate_id: Optional[str] = None,
    ):
        m = hashlib.sha256()
        m.update(certificate_pem.encode("utf-8"))
        self.certificate_id = m.hexdigest()
        self.arn = f"arn:{get_partition(region_name)}:iot:{region_name}:{account_id}:cert/{self.certificate_id}"
        self.certificate_pem = certificate_pem
        self.status = status

        self.owner = account_id
        self.transfer_data: Dict[str, str] = {}
        self.creation_date = time.time()
        self.last_modified_date = self.creation_date
        self.validity_not_before = time.time() - 86400
        self.validity_not_after = time.time() + 86400
        self.ca_certificate_id = ca_certificate_id

    def to_dict(self) -> Dict[str, Any]:
        return {
            "certificateArn": self.arn,
            "certificateId": self.certificate_id,
            "caCertificateId": self.ca_certificate_id,
            "status": self.status,
            "creationDate": self.creation_date,
        }

    def to_description_dict(self) -> Dict[str, Any]:
        """
        You might need keys below in some situation
          - caCertificateId
          - previousOwnedBy
        """
        return {
            "certificateArn": self.arn,
            "certificateId": self.certificate_id,
            "status": self.status,
            "certificatePem": self.certificate_pem,
            "ownedBy": self.owner,
            "creationDate": self.creation_date,
            "lastModifiedDate": self.last_modified_date,
            "validity": {
                "notBefore": self.validity_not_before,
                "notAfter": self.validity_not_after,
            },
            "transferData": self.transfer_data,
        }


class FakeCaCertificate(FakeCertificate):
    def __init__(
        self,
        ca_certificate: str,
        status: str,
        account_id: str,
        region_name: str,
        registration_config: Dict[str, str],
    ):
        super().__init__(
            certificate_pem=ca_certificate,
            status=status,
            account_id=account_id,
            region_name=region_name,
            ca_certificate_id=None,
        )
        self.registration_config = registration_config


class FakePolicy(BaseModel):
    def __init__(
        self,
        name: str,
        document: Dict[str, Any],
        account_id: str,
        region_name: str,
        default_version_id: str = "1",
    ):
        self.name = name
        self.document = document
        self.arn = f"arn:{get_partition(region_name)}:iot:{region_name}:{account_id}:policy/{name}"
        self.default_version_id = default_version_id
        self.versions = [
            FakePolicyVersion(self.name, document, True, account_id, region_name)
        ]
        self._max_version_id = self.versions[0]._version_id

    def to_get_dict(self) -> Dict[str, Any]:
        return {
            "policyName": self.name,
            "policyArn": self.arn,
            "policyDocument": self.document,
            "defaultVersionId": self.default_version_id,
        }

    def to_dict_at_creation(self) -> Dict[str, Any]:
        return {
            "policyName": self.name,
            "policyArn": self.arn,
            "policyDocument": self.document,
            "policyVersionId": self.default_version_id,
        }

    def to_dict(self) -> Dict[str, str]:
        return {"policyName": self.name, "policyArn": self.arn}


class FakePolicyVersion:
    def __init__(
        self,
        policy_name: str,
        document: Dict[str, Any],
        is_default: bool,
        account_id: str,
        region_name: str,
        version_id: int = 1,
    ):
        self.name = policy_name
        self.arn = f"arn:{get_partition(region_name)}:iot:{region_name}:{account_id}:policy/{policy_name}"
        self.document = document or {}
        self.is_default = is_default
        self._version_id = version_id

        self.create_datetime = time.mktime(datetime(2015, 1, 1).timetuple())
        self.last_modified_datetime = time.mktime(datetime(2015, 1, 2).timetuple())

    @property
    def version_id(self) -> str:
        return str(self._version_id)

    def to_get_dict(self) -> Dict[str, Any]:
        return {
            "policyName": self.name,
            "policyArn": self.arn,
            "policyDocument": self.document,
            "policyVersionId": self.version_id,
            "isDefaultVersion": self.is_default,
            "creationDate": self.create_datetime,
            "lastModifiedDate": self.last_modified_datetime,
            "generationId": self.version_id,
        }

    def to_dict_at_creation(self) -> Dict[str, Any]:
        return {
            "policyArn": self.arn,
            "policyDocument": self.document,
            "policyVersionId": self.version_id,
            "isDefaultVersion": self.is_default,
        }

    def to_dict(self) -> Dict[str, Any]:
        return {
            "versionId": self.version_id,
            "isDefaultVersion": self.is_default,
            "createDate": self.create_datetime,
        }


class FakeJob(BaseModel):
    JOB_ID_REGEX_PATTERN = "[a-zA-Z0-9_-]"
    JOB_ID_REGEX = re.compile(JOB_ID_REGEX_PATTERN)

    def __init__(
        self,
        job_id: str,
        targets: List[str],
        document_source: str,
        document: str,
        description: str,
        presigned_url_config: Dict[str, Any],
        target_selection: str,
        job_executions_rollout_config: Dict[str, Any],
        document_parameters: Dict[str, str],
        account_id: str,
        region_name: str,
    ):
        if not self._job_id_matcher(self.JOB_ID_REGEX, job_id):
            raise InvalidRequestException()

        self.account_id = account_id
        self.region_name = region_name
        self.job_id = job_id
        self.job_arn = f"arn:{get_partition(self.region_name)}:iot:{self.region_name}:{self.account_id}:job/{job_id}"
        self.targets = targets
        self.document_source = document_source
        self.document = document
        self.force = False
        self.description = description
        self.presigned_url_config = presigned_url_config
        self.target_selection = target_selection
        self.job_executions_rollout_config = job_executions_rollout_config
        self.status = "QUEUED"  # IN_PROGRESS | CANCELED | COMPLETED
        self.comment: Optional[str] = None
        self.reason_code: Optional[str] = None
        self.created_at = time.mktime(datetime(2015, 1, 1).timetuple())
        self.last_updated_at = time.mktime(datetime(2015, 1, 1).timetuple())
        self.completed_at = None
        self.job_process_details = {
            "processingTargets": targets,
            "numberOfQueuedThings": 1,
            "numberOfCanceledThings": 0,
            "numberOfSucceededThings": 0,
            "numberOfFailedThings": 0,
            "numberOfRejectedThings": 0,
            "numberOfInProgressThings": 0,
            "numberOfRemovedThings": 0,
        }
        self.document_parameters = document_parameters

    def to_dict(self) -> Dict[str, Any]:
        obj = {
            "jobArn": self.job_arn,
            "jobId": self.job_id,
            "targets": self.targets,
            "description": self.description,
            "presignedUrlConfig": self.presigned_url_config,
            "targetSelection": self.target_selection,
            "jobExecutionsRolloutConfig": self.job_executions_rollout_config,
            "status": self.status,
            "comment": self.comment,
            "forceCanceled": self.force,
            "reasonCode": self.reason_code,
            "createdAt": self.created_at,
            "lastUpdatedAt": self.last_updated_at,
            "completedAt": self.completed_at,
            "jobProcessDetails": self.job_process_details,
            "documentParameters": self.document_parameters,
            "document": self.document,
            "documentSource": self.document_source,
        }

        return obj

    def _job_id_matcher(self, regex: Pattern[str], argument: str) -> bool:
        regex_match = regex.match(argument) is not None
        length_match = len(argument) <= 64
        return regex_match and length_match


class FakeJobExecution(BaseModel):
    def __init__(
        self,
        job_id: str,
        thing_arn: str,
        status: str = "QUEUED",
        force_canceled: bool = False,
        status_details_map: Optional[Dict[str, Any]] = None,
    ):
        self.job_id = job_id
        self.status = status  # IN_PROGRESS | CANCELED | COMPLETED
        self.force_canceled = force_canceled
        self.status_details_map = status_details_map or {}
        self.thing_arn = thing_arn
        self.queued_at = time.mktime(datetime(2015, 1, 1).timetuple())
        self.started_at = time.mktime(datetime(2015, 1, 1).timetuple())
        self.last_updated_at = time.mktime(datetime(2015, 1, 1).timetuple())
        self.execution_number = 123
        self.version_number = 123
        self.approximate_seconds_before_time_out = 123

    def to_get_dict(self) -> Dict[str, Any]:
        return {
            "jobId": self.job_id,
            "status": self.status,
            "forceCanceled": self.force_canceled,
            "statusDetails": {"detailsMap": self.status_details_map},
            "thingArn": self.thing_arn,
            "queuedAt": self.queued_at,
            "startedAt": self.started_at,
            "lastUpdatedAt": self.last_updated_at,
            "executionNumber": self.execution_number,
            "versionNumber": self.version_number,
            "approximateSecondsBeforeTimedOut": self.approximate_seconds_before_time_out,
        }

    def to_dict(self) -> Dict[str, Any]:
        return {
            "jobId": self.job_id,
            "thingArn": self.thing_arn,
            "jobExecutionSummary": {
                "status": self.status,
                "queuedAt": self.queued_at,
                "startedAt": self.started_at,
                "lastUpdatedAt": self.last_updated_at,
                "executionNumber": self.execution_number,
            },
        }


class FakeEndpoint(BaseModel):
    def __init__(self, endpoint_type: str, region_name: str):
        if endpoint_type not in [
            "iot:Data",
            "iot:Data-ATS",
            "iot:CredentialProvider",
            "iot:Jobs",
        ]:
            raise InvalidRequestException(
                " An error occurred (InvalidRequestException) when calling the DescribeEndpoint "
                f"operation: Endpoint type {endpoint_type} not recognized."
            )
        self.region_name = region_name
        identifier = random.get_random_string(length=14, lower_case=True)
        if endpoint_type == "iot:Data":
            self.endpoint = f"{identifier}.iot.{self.region_name}.amazonaws.com"
        elif "iot:Data-ATS" in endpoint_type:
            self.endpoint = f"{identifier}-ats.iot.{self.region_name}.amazonaws.com"
        elif "iot:CredentialProvider" in endpoint_type:
            self.endpoint = (
                f"{identifier}.credentials.iot.{self.region_name}.amazonaws.com"
            )
        elif "iot:Jobs" in endpoint_type:
            self.endpoint = f"{identifier}.jobs.iot.{self.region_name}.amazonaws.com"
        self.endpoint_type = endpoint_type

    def to_get_dict(self) -> Dict[str, str]:
        return {
            "endpointAddress": self.endpoint,
        }

    def to_dict(self) -> Dict[str, str]:
        return {
            "endpointAddress": self.endpoint,
        }


class FakeRule(BaseModel):
    def __init__(
        self,
        rule_name: str,
        description: str,
        created_at: int,
        rule_disabled: bool,
        topic_pattern: Optional[str],
        actions: List[Dict[str, Any]],
        error_action: Dict[str, Any],
        sql: str,
        aws_iot_sql_version: str,
        account_id: str,
        region_name: str,
    ):
        self.account_id = account_id
        self.region_name = region_name
        self.rule_name = rule_name
        self.description = description or ""
        self.created_at = created_at
        self.rule_disabled = bool(rule_disabled)
        self.topic_pattern = topic_pattern
        self.actions = actions or []
        self.error_action = error_action or {}
        self.sql = sql
        self.aws_iot_sql_version = aws_iot_sql_version or "2016-03-23"
        self.arn = f"arn:{get_partition(self.region_name)}:iot:{self.region_name}:{self.account_id}:rule/{rule_name}"

    def to_get_dict(self) -> Dict[str, Any]:
        return {
            "rule": {
                "actions": self.actions,
                "awsIotSqlVersion": self.aws_iot_sql_version,
                "createdAt": self.created_at,
                "description": self.description,
                "errorAction": self.error_action,
                "ruleDisabled": self.rule_disabled,
                "ruleName": self.rule_name,
                "sql": self.sql,
            },
            "ruleArn": self.arn,
        }

    def to_dict(self) -> Dict[str, Any]:
        return {
            "ruleName": self.rule_name,
            "createdAt": self.created_at,
            "ruleArn": self.arn,
            "ruleDisabled": self.rule_disabled,
            "topicPattern": self.topic_pattern,
        }


class FakeDomainConfiguration(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        domain_configuration_name: str,
        domain_name: str,
        server_certificate_arns: List[str],
        domain_configuration_status: str,
        service_type: str,
        authorizer_config: Optional[Dict[str, Any]],
        domain_type: str,
    ):
        if service_type and service_type not in ["DATA", "CREDENTIAL_PROVIDER", "JOBS"]:
            raise InvalidRequestException(
                "An error occurred (InvalidRequestException) when calling the DescribeDomainConfiguration "
                f"operation: Service type {service_type} not recognized."
            )
        self.domain_configuration_name = domain_configuration_name
        self.domain_configuration_arn = f"arn:{get_partition(region_name)}:iot:{region_name}:{account_id}:domainconfiguration/{domain_configuration_name}/{random.get_random_string(length=5)}"
        self.domain_name = domain_name
        self.server_certificates = []
        if server_certificate_arns:
            for sc in server_certificate_arns:
                self.server_certificates.append(
                    {"serverCertificateArn": sc, "serverCertificateStatus": "VALID"}
                )
        self.domain_configuration_status = domain_configuration_status
        self.service_type = service_type
        self.authorizer_config = authorizer_config
        self.domain_type = domain_type
        self.last_status_change_date = time.time()

    def to_description_dict(self) -> Dict[str, Any]:
        return {
            "domainConfigurationName": self.domain_configuration_name,
            "domainConfigurationArn": self.domain_configuration_arn,
            "domainName": self.domain_name,
            "serverCertificates": self.server_certificates,
            "authorizerConfig": self.authorizer_config,
            "domainConfigurationStatus": self.domain_configuration_status,
            "serviceType": self.service_type,
            "domainType": self.domain_type,
            "lastStatusChangeDate": self.last_status_change_date,
        }

    def to_dict(self) -> Dict[str, Any]:
        return {
            "domainConfigurationName": self.domain_configuration_name,
            "domainConfigurationArn": self.domain_configuration_arn,
        }


class IoTBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.things: Dict[str, FakeThing] = OrderedDict()
        self.jobs: Dict[str, FakeJob] = OrderedDict()
        self.job_executions: Dict[Tuple[str, str], FakeJobExecution] = OrderedDict()
        self.thing_types: Dict[str, FakeThingType] = OrderedDict()
        self.thing_groups: Dict[str, FakeThingGroup] = OrderedDict()
        self.ca_certificates: Dict[str, FakeCaCertificate] = OrderedDict()
        self.certificates: Dict[str, FakeCertificate] = OrderedDict()
        self.policies: Dict[str, FakePolicy] = OrderedDict()
        self.principal_policies: Dict[Tuple[str, str], Tuple[str, FakePolicy]] = (
            OrderedDict()
        )
        self.principal_things: Dict[Tuple[str, str], Tuple[str, FakeThing]] = (
            OrderedDict()
        )
        self.rules: Dict[str, FakeRule] = OrderedDict()
        self.endpoint: Optional[FakeEndpoint] = None
        self.domain_configurations: Dict[str, FakeDomainConfiguration] = OrderedDict()

    @staticmethod
    def default_vpc_endpoint_service(
        service_region: str, zones: List[str]
    ) -> List[Dict[str, str]]:
        """Default VPC endpoint service."""
        return BaseBackend.default_vpc_endpoint_service_factory(
            service_region, zones, "iot"
        ) + BaseBackend.default_vpc_endpoint_service_factory(
            service_region,
            zones,
            "data.iot",
            private_dns_names=False,
            special_service_name="iot.data",
            policy_supported=False,
        )

    def create_certificate_from_csr(
        self, csr: str, set_as_active: bool
    ) -> FakeCertificate:
        cert = x509.load_pem_x509_csr(csr.encode("utf-8"), default_backend())
        pem = self._generate_certificate_pem(
            domain_name="example.com", subject=cert.subject
        )
        return self.register_certificate(
            pem, ca_certificate_pem=None, set_as_active=set_as_active, status="INACTIVE"
        )

    def _generate_certificate_pem(
        self,
        domain_name: str,
        subject: x509.Name,
        key: Optional[rsa.RSAPrivateKey] = None,
    ) -> str:
        sans = [x509.DNSName(domain_name)]

        key = key or rsa.generate_private_key(
            public_exponent=65537, key_size=2048, backend=default_backend()
        )
        issuer = x509.Name(
            [  # C = US, O = Moto, OU = Server CA 1B, CN = Moto
                x509.NameAttribute(x509.NameOID.COUNTRY_NAME, "US"),
                x509.NameAttribute(x509.NameOID.ORGANIZATION_NAME, "Moto"),
                x509.NameAttribute(
                    x509.NameOID.ORGANIZATIONAL_UNIT_NAME, "Server CA 1B"
                ),
                x509.NameAttribute(x509.NameOID.COMMON_NAME, "Moto"),
            ]
        )
        cert = (
            x509.CertificateBuilder()
            .subject_name(subject)
            .issuer_name(issuer)
            .public_key(key.public_key())
            .serial_number(x509.random_serial_number())
            .not_valid_before(utcnow())
            .not_valid_after(utcnow() + timedelta(days=365))
            .add_extension(x509.SubjectAlternativeName(sans), critical=False)
            .sign(key, hashes.SHA512(), default_backend())
        )

        return cert.public_bytes(serialization.Encoding.PEM).decode("utf-8")

    def create_thing(
        self,
        thing_name: str,
        thing_type_name: str,
        attribute_payload: Optional[Dict[str, Any]],
    ) -> FakeThing:
        thing_types = self.list_thing_types()
        thing_type = None
        if thing_type_name:
            filtered_thing_types = [
                _ for _ in thing_types if _.thing_type_name == thing_type_name
            ]
            if len(filtered_thing_types) == 0:
                raise ResourceNotFoundException()
            thing_type = filtered_thing_types[0]

            if thing_type.metadata["deprecated"]:
                # note - typo (depreated) exists also in the original exception.
                raise InvalidRequestException(
                    msg=f"Can not create new thing with depreated thing type:{thing_type_name}"
                )
        if attribute_payload is None:
            attributes: Dict[str, Any] = {}
        elif "attributes" not in attribute_payload:
            attributes = {}
        else:
            attributes = attribute_payload["attributes"]
        thing = FakeThing(
            thing_name, thing_type, attributes, self.account_id, self.region_name
        )
        self.things[thing.arn] = thing
        return thing

    def create_thing_type(
        self, thing_type_name: str, thing_type_properties: Dict[str, Any]
    ) -> Tuple[str, str]:
        if thing_type_properties is None:
            thing_type_properties = {}
        thing_type = FakeThingType(
            thing_type_name, thing_type_properties, self.account_id, self.region_name
        )
        self.thing_types[thing_type.arn] = thing_type
        return thing_type.thing_type_name, thing_type.arn

    def list_thing_types(
        self, thing_type_name: Optional[str] = None
    ) -> Iterable[FakeThingType]:
        if thing_type_name:
            # It's weird but thing_type_name is filtered by forward match, not complete match
            return [
                _
                for _ in self.thing_types.values()
                if _.thing_type_name.startswith(thing_type_name)
            ]
        return self.thing_types.values()

    def list_things(
        self,
        attribute_name: str,
        attribute_value: str,
        thing_type_name: str,
        max_results: int,
        token: Optional[str],
    ) -> Tuple[Iterable[FakeThing], Optional[str]]:
        all_things = [_.to_dict() for _ in self.things.values()]
        if attribute_name is not None and thing_type_name is not None:
            filtered_things = list(
                filter(
                    lambda elem: attribute_name in elem["attributes"]
                    and elem["attributes"][attribute_name] == attribute_value
                    and "thingTypeName" in elem
                    and elem["thingTypeName"] == thing_type_name,
                    all_things,
                )
            )
        elif attribute_name is not None and thing_type_name is None:
            filtered_things = list(
                filter(
                    lambda elem: attribute_name in elem["attributes"]
                    and elem["attributes"][attribute_name] == attribute_value,
                    all_things,
                )
            )
        elif attribute_name is None and thing_type_name is not None:
            filtered_things = list(
                filter(
                    lambda elem: "thingTypeName" in elem
                    and elem["thingTypeName"] == thing_type_name,
                    all_things,
                )
            )
        else:
            filtered_things = all_things

        if token is None:
            things = filtered_things[0:max_results]
            next_token = (
                str(max_results) if len(filtered_things) > max_results else None
            )
        else:
            int_token = int(token)
            things = filtered_things[int_token : int_token + max_results]
            next_token = (
                str(int_token + max_results)
                if len(filtered_things) > int_token + max_results
                else None
            )

        return things, next_token

    def describe_thing(self, thing_name: str) -> FakeThing:
        things = [_ for _ in self.things.values() if _.thing_name == thing_name]
        if len(things) == 0:
            raise ResourceNotFoundException()
        return things[0]

    def describe_thing_type(self, thing_type_name: str) -> FakeThingType:
        thing_types = [
            _ for _ in self.thing_types.values() if _.thing_type_name == thing_type_name
        ]
        if len(thing_types) == 0:
            raise ResourceNotFoundException()
        return thing_types[0]

    def describe_endpoint(self, endpoint_type: str) -> FakeEndpoint:
        self.endpoint = FakeEndpoint(endpoint_type, self.region_name)
        return self.endpoint

    def delete_thing(self, thing_name: str) -> None:
        """
        The ExpectedVersion-parameter is not yet implemented
        """

        # can raise ResourceNotFoundError
        thing = self.describe_thing(thing_name)

        for k in list(self.principal_things.keys()):
            if k[1] == thing_name:
                raise ThingStillAttached(thing_name)

        del self.things[thing.arn]

    def delete_thing_type(self, thing_type_name: str) -> None:
        # can raise ResourceNotFoundError
        thing_type = self.describe_thing_type(thing_type_name)
        del self.thing_types[thing_type.arn]

    def deprecate_thing_type(
        self, thing_type_name: str, undo_deprecate: bool
    ) -> FakeThingType:
        thing_types = [
            _ for _ in self.thing_types.values() if _.thing_type_name == thing_type_name
        ]
        if len(thing_types) == 0:
            raise ResourceNotFoundException()
        thing_types[0].metadata["deprecated"] = not undo_deprecate
        return thing_types[0]

    def update_thing(
        self,
        thing_name: str,
        thing_type_name: str,
        attribute_payload: Optional[Dict[str, Any]],
        remove_thing_type: bool,
    ) -> None:
        """
        The ExpectedVersion-parameter is not yet implemented
        """
        # if attributes payload = {}, nothing
        thing = self.describe_thing(thing_name)

        if remove_thing_type and thing_type_name:
            raise InvalidRequestException()

        # thing_type
        if thing_type_name:
            thing_types = self.list_thing_types()
            filtered_thing_types = [
                _ for _ in thing_types if _.thing_type_name == thing_type_name
            ]
            if len(filtered_thing_types) == 0:
                raise ResourceNotFoundException()
            thing_type = filtered_thing_types[0]

            if thing_type.metadata["deprecated"]:
                raise InvalidRequestException(
                    msg=f"Can not update a thing to use deprecated thing type: {thing_type_name}"
                )

            thing.thing_type = thing_type

        if remove_thing_type:
            thing.thing_type = None

        # attribute
        if attribute_payload is not None and "attributes" in attribute_payload:
            do_merge = attribute_payload.get("merge", False)
            attributes = attribute_payload["attributes"]
            if not do_merge:
                thing.attributes = attributes
            else:
                thing.attributes.update(attributes)
            thing.attributes = {k: v for k, v in thing.attributes.items() if v}

    def create_keys_and_certificate(
        self, set_as_active: bool
    ) -> Tuple[FakeCertificate, Dict[str, str]]:
        # implement here
        # caCertificate can be blank
        private_key = rsa.generate_private_key(
            public_exponent=65537, key_size=2048, backend=default_backend()
        )
        key_pair = {
            "PublicKey": private_key.public_key()
            .public_bytes(
                encoding=serialization.Encoding.PEM,
                format=serialization.PublicFormat.SubjectPublicKeyInfo,
            )
            .decode("utf-8"),
            "PrivateKey": private_key.private_bytes(
                encoding=serialization.Encoding.PEM,
                format=serialization.PrivateFormat.TraditionalOpenSSL,
                encryption_algorithm=serialization.NoEncryption(),
            ).decode("utf-8"),
        }
        subject = x509.Name(
            [x509.NameAttribute(NameOID.COMMON_NAME, "AWS IoT Certificate")]
        )
        certificate_pem = self._generate_certificate_pem(
            "getmoto.org", subject, key=private_key
        )
        status = "ACTIVE" if set_as_active else "INACTIVE"
        certificate = FakeCertificate(
            certificate_pem, status, self.account_id, self.region_name
        )
        self.certificates[certificate.certificate_id] = certificate
        return certificate, key_pair

    def delete_ca_certificate(self, certificate_id: str) -> None:
        cert = self.describe_ca_certificate(certificate_id)
        self._validation_delete(cert)
        del self.ca_certificates[certificate_id]

    def delete_certificate(self, certificate_id: str, force_delete: bool) -> None:
        cert = self.describe_certificate(certificate_id)
        self._validation_delete(cert, force_delete)
        del self.certificates[certificate_id]

    def _validation_delete(
        self, cert: FakeCertificate, force_delete: bool = False
    ) -> None:
        if cert.status == "ACTIVE":
            raise CertificateStateException(
                "Certificate must be deactivated (not ACTIVE) before deletion.",
                cert.certificate_id,
            )

        certs = [
            k[0]
            for k, v in self.principal_things.items()
            if self._get_principal(k[0]).certificate_id == cert.certificate_id
        ]
        if len(certs) > 0:
            raise DeleteConflictException(
                f"Things must be detached before deletion (arn: {certs[0]})"
            )

        certs = [
            k[0]
            for k, v in self.principal_policies.items()
            if self._get_principal(k[0]).certificate_id == cert.certificate_id
        ]
        if len(certs) > 0 and not force_delete:
            raise DeleteConflictException(
                "Certificate policies must be detached before deletion (arn: %s)"
                % certs[0]
            )

    def describe_ca_certificate(self, certificate_id: str) -> FakeCaCertificate:
        if certificate_id not in self.ca_certificates:
            raise ResourceNotFoundException()
        return self.ca_certificates[certificate_id]

    def describe_certificate(self, certificate_id: str) -> FakeCertificate:
        certs = [
            _ for _ in self.certificates.values() if _.certificate_id == certificate_id
        ]
        if len(certs) == 0:
            raise ResourceNotFoundException()
        return certs[0]

    def get_registration_code(self) -> str:
        return str(random.uuid4())

    def list_certificates(self) -> Iterable[FakeCertificate]:
        """
        Pagination is not yet implemented
        """
        return self.certificates.values()

    def list_certificates_by_ca(self, ca_certificate_id: str) -> List[FakeCertificate]:
        """
        Pagination is not yet implemented
        """
        return [
            cert
            for cert in self.certificates.values()
            if cert.ca_certificate_id == ca_certificate_id
        ]

    def __raise_if_certificate_already_exists(
        self, certificate_id: str, certificate_arn: str
    ) -> None:
        if certificate_id in self.certificates:
            raise ResourceAlreadyExistsException(
                "The certificate is already provisioned or registered",
                certificate_id,
                certificate_arn,
            )

    def register_ca_certificate(
        self,
        ca_certificate: str,
        set_as_active: bool,
        registration_config: Dict[str, str],
    ) -> FakeCaCertificate:
        """
        The VerificationCertificate-parameter is not yet implemented
        """
        certificate = FakeCaCertificate(
            ca_certificate=ca_certificate,
            status="ACTIVE" if set_as_active else "INACTIVE",
            account_id=self.account_id,
            region_name=self.region_name,
            registration_config=registration_config,
        )

        self.ca_certificates[certificate.certificate_id] = certificate
        return certificate

    def _find_ca_certificate(self, ca_certificate_pem: Optional[str]) -> Optional[str]:
        for ca_cert in self.ca_certificates.values():
            if ca_cert.certificate_pem == ca_certificate_pem:
                return ca_cert.certificate_id
        return None

    def register_certificate(
        self,
        certificate_pem: str,
        ca_certificate_pem: Optional[str],
        set_as_active: bool,
        status: str,
    ) -> FakeCertificate:
        ca_certificate_id = self._find_ca_certificate(ca_certificate_pem)
        certificate = FakeCertificate(
            certificate_pem,
            "ACTIVE" if set_as_active else status,
            self.account_id,
            self.region_name,
            ca_certificate_id,
        )
        self.__raise_if_certificate_already_exists(
            certificate.certificate_id, certificate_arn=certificate.arn
        )

        self.certificates[certificate.certificate_id] = certificate
        return certificate

    def register_certificate_without_ca(
        self, certificate_pem: str, status: str
    ) -> FakeCertificate:
        certificate = FakeCertificate(
            certificate_pem, status, self.account_id, self.region_name
        )
        self.__raise_if_certificate_already_exists(
            certificate.certificate_id, certificate_arn=certificate.arn
        )

        self.certificates[certificate.certificate_id] = certificate
        return certificate

    def update_ca_certificate(
        self,
        certificate_id: str,
        new_status: Optional[str],
        config: Optional[Dict[str, str]],
    ) -> None:
        """
        The newAutoRegistrationStatus and removeAutoRegistration-parameters are not yet implemented
        """
        cert = self.describe_ca_certificate(certificate_id)
        if new_status is not None:
            cert.status = new_status
        if config is not None:
            cert.registration_config = config

    def update_certificate(self, certificate_id: str, new_status: str) -> None:
        cert = self.describe_certificate(certificate_id)
        # TODO: validate new_status
        cert.status = new_status

    def create_policy(
        self, policy_name: str, policy_document: Dict[str, Any]
    ) -> FakePolicy:
        if policy_name in self.policies:
            current_policy = self.policies[policy_name]
            raise ResourceAlreadyExistsException(
                f"Policy cannot be created - name already exists (name={policy_name})",
                current_policy.name,
                current_policy.arn,
            )
        policy = FakePolicy(
            policy_name, policy_document, self.account_id, self.region_name
        )
        self.policies[policy.name] = policy
        return policy

    def attach_policy(self, policy_name: str, target: str) -> None:
        principal = self._get_principal(target)
        policy = self.get_policy(policy_name)
        k = (target, policy_name)
        if k in self.principal_policies:
            return
        self.principal_policies[k] = (principal, policy)

    def detach_policy(self, policy_name: str, target: str) -> None:
        # this may raise ResourceNotFoundException
        self._get_principal(target)
        self.get_policy(policy_name)

        k = (target, policy_name)
        if k not in self.principal_policies:
            raise ResourceNotFoundException()
        del self.principal_policies[k]

    def list_attached_policies(self, target: str) -> List[FakePolicy]:
        """
        Pagination is not yet implemented
        """
        return [v[1] for k, v in self.principal_policies.items() if k[0] == target]

    def list_policies(self) -> Iterable[FakePolicy]:
        """
        Pagination is not yet implemented
        """
        return self.policies.values()

    def get_policy(self, policy_name: str) -> FakePolicy:
        policies = [_ for _ in self.policies.values() if _.name == policy_name]
        if len(policies) == 0:
            raise ResourceNotFoundException()
        return policies[0]

    def delete_policy(self, policy_name: str) -> None:
        policies = [
            k[1] for k, v in self.principal_policies.items() if k[1] == policy_name
        ]
        if len(policies) > 0:
            raise DeleteConflictException(
                "The policy cannot be deleted as the policy is attached to one or more principals (name=%s)"
                % policy_name
            )

        policy = self.get_policy(policy_name)
        if len(policy.versions) > 1:
            raise DeleteConflictException(
                "Cannot delete the policy because it has one or more policy versions attached to it (name=%s)"
                % policy_name
            )
        del self.policies[policy.name]

    def create_policy_version(
        self, policy_name: str, policy_document: Dict[str, Any], set_as_default: bool
    ) -> FakePolicyVersion:
        policy = self.get_policy(policy_name)
        if not policy:
            raise ResourceNotFoundException()
        if len(policy.versions) >= 5:
            raise VersionsLimitExceededException(policy_name)

        policy._max_version_id += 1
        version = FakePolicyVersion(
            policy_name,
            policy_document,
            set_as_default,
            self.account_id,
            self.region_name,
            version_id=policy._max_version_id,
        )
        policy.versions.append(version)
        if set_as_default:
            self.set_default_policy_version(policy_name, version.version_id)
        return version

    def set_default_policy_version(self, policy_name: str, version_id: str) -> None:
        policy = self.get_policy(policy_name)
        if not policy:
            raise ResourceNotFoundException()
        for version in policy.versions:
            if version.version_id == version_id:
                version.is_default = True
                policy.default_version_id = version.version_id
                policy.document = version.document
            else:
                version.is_default = False

    def get_policy_version(
        self, policy_name: str, version_id: str
    ) -> FakePolicyVersion:
        policy = self.get_policy(policy_name)
        if not policy:
            raise ResourceNotFoundException()
        for version in policy.versions:
            if version.version_id == version_id:
                return version
        raise ResourceNotFoundException()

    def list_policy_versions(self, policy_name: str) -> Iterable[FakePolicyVersion]:
        policy = self.get_policy(policy_name)
        if not policy:
            raise ResourceNotFoundException()
        return policy.versions

    def delete_policy_version(self, policy_name: str, version_id: str) -> None:
        policy = self.get_policy(policy_name)
        if not policy:
            raise ResourceNotFoundException()
        if version_id == policy.default_version_id:
            raise InvalidRequestException(
                "Cannot delete the default version of a policy"
            )
        for i, v in enumerate(policy.versions):
            if v.version_id == version_id:
                del policy.versions[i]
                return
        raise ResourceNotFoundException()

    def _get_principal(self, principal_arn: str) -> Any:
        """
        raise ResourceNotFoundException
        """
        if ":cert/" in principal_arn:
            certs = [_ for _ in self.certificates.values() if _.arn == principal_arn]
            if len(certs) == 0:
                raise ResourceNotFoundException()
            principal = certs[0]
            return principal
        if ":thinggroup/" in principal_arn:
            try:
                return self.thing_groups[principal_arn]
            except KeyError:
                raise ResourceNotFoundException(
                    f"No thing group with ARN {principal_arn} exists"
                )
        from moto.cognitoidentity import cognitoidentity_backends

        cognito = cognitoidentity_backends[self.account_id][self.region_name]
        identities = []
        for identity_pool in cognito.identity_pools:
            pool_identities = cognito.pools_identities.get(identity_pool, {})
            identities.extend(
                [pi["IdentityId"] for pi in pool_identities.get("Identities", [])]
            )
            if principal_arn in identities:
                return {"IdentityId": principal_arn}

        raise ResourceNotFoundException()

    def attach_principal_policy(self, policy_name: str, principal_arn: str) -> None:
        principal = self._get_principal(principal_arn)
        policy = self.get_policy(policy_name)
        k = (principal_arn, policy_name)
        if k in self.principal_policies:
            return
        self.principal_policies[k] = (principal, policy)

    def detach_principal_policy(self, policy_name: str, principal_arn: str) -> None:
        # this may raises ResourceNotFoundException
        self._get_principal(principal_arn)
        self.get_policy(policy_name)

        k = (principal_arn, policy_name)
        if k not in self.principal_policies:
            raise ResourceNotFoundException()
        del self.principal_policies[k]

    def list_principal_policies(self, principal_arn: str) -> List[FakePolicy]:
        """
        Pagination is not yet implemented
        """
        policies = [
            v[1] for k, v in self.principal_policies.items() if k[0] == principal_arn
        ]
        return policies

    def list_policy_principals(self, policy_name: str) -> List[str]:
        """
        Pagination is not yet implemented
        """
        # this action is deprecated
        # https://docs.aws.amazon.com/iot/latest/apireference/API_ListTargetsForPolicy.html
        # should use ListTargetsForPolicy instead
        principals = [
            k[0] for k, v in self.principal_policies.items() if k[1] == policy_name
        ]
        return principals

    def list_targets_for_policy(self, policy_name: str) -> List[str]:
        """
        Pagination is not yet implemented
        """
        # This behaviour is different to list_policy_principals which will just return an empty list
        if policy_name not in self.policies:
            raise ResourceNotFoundException("Policy not found")
        return self.list_policy_principals(policy_name=policy_name)

    def attach_thing_principal(self, thing_name: str, principal_arn: str) -> None:
        principal = self._get_principal(principal_arn)
        thing = self.describe_thing(thing_name)
        k = (principal_arn, thing_name)
        if k in self.principal_things:
            return
        self.principal_things[k] = (principal, thing)

    def detach_thing_principal(self, thing_name: str, principal_arn: str) -> None:
        # this may raises ResourceNotFoundException
        self._get_principal(principal_arn)
        self.describe_thing(thing_name)

        k = (principal_arn, thing_name)
        if k not in self.principal_things:
            raise ResourceNotFoundException()
        del self.principal_things[k]

    def list_principal_things(self, principal_arn: str) -> List[str]:
        thing_names = [
            k[1] for k, v in self.principal_things.items() if k[0] == principal_arn
        ]
        return thing_names

    def list_thing_principals(self, thing_name: str) -> List[str]:
        things = [_ for _ in self.things.values() if _.thing_name == thing_name]
        if len(things) == 0:
            raise ResourceNotFoundException(
                "Failed to list principals for thing %s because the thing does not exist in your account"
                % thing_name
            )

        principals = [
            k[0] for k, v in self.principal_things.items() if k[1] == thing_name
        ]
        return principals

    def describe_thing_group(self, thing_group_name: str) -> FakeThingGroup:
        thing_groups = [
            _
            for _ in self.thing_groups.values()
            if _.thing_group_name == thing_group_name
        ]
        if len(thing_groups) == 0:
            raise ResourceNotFoundException()
        return thing_groups[0]

    def create_thing_group(
        self,
        thing_group_name: str,
        parent_group_name: str,
        thing_group_properties: Dict[str, Any],
    ) -> Tuple[str, str, str]:
        thing_group = FakeThingGroup(
            thing_group_name,
            parent_group_name,
            thing_group_properties,
            self.account_id,
            self.region_name,
            self.thing_groups,
        )
        # this behavior is not documented, but AWS does it like that
        # if a thing group with the same name exists, it's properties are compared
        # if they differ, an error is returned.
        # Otherwise, the old thing group is returned
        if thing_group.arn in self.thing_groups:
            current_thing_group = self.thing_groups[thing_group.arn]
            if current_thing_group.thing_group_properties != thing_group_properties:
                raise ResourceAlreadyExistsException(
                    msg=f"Thing Group {thing_group_name} already exists in current account with different properties",
                    resource_arn=thing_group.arn,
                    resource_id=current_thing_group.thing_group_id,
                )
            thing_group = current_thing_group
        else:
            self.thing_groups[thing_group.arn] = thing_group
        return thing_group.thing_group_name, thing_group.arn, thing_group.thing_group_id

    def delete_thing_group(self, thing_group_name: str) -> None:
        """
        The ExpectedVersion-parameter is not yet implemented
        """
        child_groups = [
            thing_group
            for _, thing_group in self.thing_groups.items()
            if thing_group.parent_group_name == thing_group_name
        ]
        if len(child_groups) > 0:
            raise InvalidRequestException(
                " Cannot delete thing group : "
                + thing_group_name
                + " when there are still child groups attached to it"
            )
        try:
            thing_group = self.describe_thing_group(thing_group_name)
            del self.thing_groups[thing_group.arn]
        except ResourceNotFoundException:
            # AWS returns success even if the thing group does not exist.
            pass

    def list_thing_groups(
        self,
        parent_group: Optional[str],
        name_prefix_filter: Optional[str],
        recursive: Optional[bool],
    ) -> List[FakeThingGroup]:
        if recursive is None:
            recursive = True
        if name_prefix_filter is None:
            name_prefix_filter = ""
        if parent_group and parent_group not in [
            _.thing_group_name for _ in self.thing_groups.values()
        ]:
            raise ResourceNotFoundException()
        thing_groups = [
            _ for _ in self.thing_groups.values() if _.parent_group_name == parent_group
        ]
        if recursive:
            for g in thing_groups:
                thing_groups.extend(
                    self.list_thing_groups(
                        parent_group=g.thing_group_name,
                        name_prefix_filter=None,
                        recursive=False,
                    )
                )
        # thing_groups = groups_to_process.values()
        return [
            _ for _ in thing_groups if _.thing_group_name.startswith(name_prefix_filter)
        ]

    def update_thing_group(
        self,
        thing_group_name: str,
        thing_group_properties: Dict[str, Any],
        expected_version: int,
    ) -> int:
        thing_group = self.describe_thing_group(thing_group_name)
        if expected_version and expected_version != thing_group.version:
            raise VersionConflictException(thing_group_name)
        attribute_payload = thing_group_properties.get("attributePayload", None)
        if attribute_payload is not None and "attributes" in attribute_payload:
            do_merge = attribute_payload.get("merge", False)
            attributes = attribute_payload["attributes"]
            if attributes:
                # might not exist yet, for example when the thing group was created without attributes
                current_attribute_payload = (
                    thing_group.thing_group_properties.setdefault(
                        "attributePayload",
                        {"attributes": {}},  # type: ignore
                    )
                )
                if not do_merge:
                    current_attribute_payload["attributes"] = attributes  # type: ignore
                else:
                    current_attribute_payload["attributes"].update(attributes)  # type: ignore
        elif attribute_payload is not None and "attributes" not in attribute_payload:
            thing_group.attributes = {}  # type: ignore
        if "thingGroupDescription" in thing_group_properties:
            thing_group.thing_group_properties["thingGroupDescription"] = (
                thing_group_properties["thingGroupDescription"]
            )
        thing_group.version = thing_group.version + 1
        return thing_group.version

    def _identify_thing_group(
        self, thing_group_name: Optional[str], thing_group_arn: Optional[str]
    ) -> FakeThingGroup:
        # identify thing group
        if thing_group_name is None and thing_group_arn is None:
            raise InvalidRequestException(
                " Both thingGroupArn and thingGroupName are empty. Need to specify at least one of them"
            )
        if thing_group_name is not None:
            thing_group = self.describe_thing_group(thing_group_name)
            if thing_group_arn and thing_group.arn != thing_group_arn:
                raise InvalidRequestException(
                    "ThingGroupName thingGroupArn does not match specified thingGroupName in request"
                )
        elif thing_group_arn is not None:
            if thing_group_arn not in self.thing_groups:
                raise InvalidRequestException()
            thing_group = self.thing_groups[thing_group_arn]
        return thing_group

    def _identify_thing(
        self, thing_name: Optional[str], thing_arn: Optional[str]
    ) -> FakeThing:
        # identify thing
        if thing_name is None and thing_arn is None:
            raise InvalidRequestException(
                "Both thingArn and thingName are empty. Need to specify at least one of them"
            )
        if thing_name is not None:
            thing = self.describe_thing(thing_name)
            if thing_arn and thing.arn != thing_arn:
                raise InvalidRequestException(
                    "ThingName thingArn does not match specified thingName in request"
                )
        elif thing_arn is not None:
            if thing_arn not in self.things:
                raise InvalidRequestException()
            thing = self.things[thing_arn]
        return thing

    def add_thing_to_thing_group(
        self,
        thing_group_name: str,
        thing_group_arn: Optional[str],
        thing_name: str,
        thing_arn: Optional[str],
    ) -> None:
        thing_group = self._identify_thing_group(thing_group_name, thing_group_arn)
        thing = self._identify_thing(thing_name, thing_arn)
        if thing.arn in thing_group.things:
            # aws ignores duplicate registration
            return
        thing_group.things[thing.arn] = thing

    def remove_thing_from_thing_group(
        self,
        thing_group_name: str,
        thing_group_arn: Optional[str],
        thing_name: str,
        thing_arn: Optional[str],
    ) -> None:
        thing_group = self._identify_thing_group(thing_group_name, thing_group_arn)
        thing = self._identify_thing(thing_name, thing_arn)
        if thing.arn not in thing_group.things:
            # aws ignores non-registered thing
            return
        del thing_group.things[thing.arn]

    def list_things_in_thing_group(self, thing_group_name: str) -> Iterable[FakeThing]:
        """
        Pagination and the recursive-parameter is not yet implemented
        """
        thing_group = self.describe_thing_group(thing_group_name)
        return thing_group.things.values()

    def list_thing_groups_for_thing(self, thing_name: str) -> List[Dict[str, str]]:
        """
        Pagination is not yet implemented
        """
        thing = self.describe_thing(thing_name)
        all_thing_groups = self.list_thing_groups(None, None, None)
        ret = []
        for thing_group in all_thing_groups:
            if thing.arn in thing_group.things:
                ret.append(
                    {
                        "groupName": thing_group.thing_group_name,
                        "groupArn": thing_group.arn,
                    }
                )
        return ret

    def update_thing_groups_for_thing(
        self,
        thing_name: str,
        thing_groups_to_add: List[str],
        thing_groups_to_remove: List[str],
    ) -> None:
        thing = self.describe_thing(thing_name)
        for thing_group_name in thing_groups_to_add:
            thing_group = self.describe_thing_group(thing_group_name)
            self.add_thing_to_thing_group(
                thing_group.thing_group_name, None, thing.thing_name, None
            )
        for thing_group_name in thing_groups_to_remove:
            thing_group = self.describe_thing_group(thing_group_name)
            self.remove_thing_from_thing_group(
                thing_group.thing_group_name, None, thing.thing_name, None
            )

    def create_job(
        self,
        job_id: str,
        targets: List[str],
        document_source: str,
        document: str,
        description: str,
        presigned_url_config: Dict[str, Any],
        target_selection: str,
        job_executions_rollout_config: Dict[str, Any],
        document_parameters: Dict[str, str],
    ) -> Tuple[str, str, str]:
        job = FakeJob(
            job_id,
            targets,
            document_source,
            document,
            description,
            presigned_url_config,
            target_selection,
            job_executions_rollout_config,
            document_parameters,
            self.account_id,
            self.region_name,
        )
        self.jobs[job_id] = job

        for thing_arn in targets:
            thing_name = thing_arn.split(":")[-1].split("/")[-1]
            job_execution = FakeJobExecution(job_id, thing_arn)
            self.job_executions[(job_id, thing_name)] = job_execution
        return job.job_arn, job_id, description

    def describe_job(self, job_id: str) -> FakeJob:
        jobs = [_ for _ in self.jobs.values() if _.job_id == job_id]
        if len(jobs) == 0:
            raise ResourceNotFoundException()
        return jobs[0]

    def delete_job(self, job_id: str, force: bool) -> None:
        job = self.jobs[job_id]

        if job.status == "IN_PROGRESS" and force:
            del self.jobs[job_id]
        elif job.status != "IN_PROGRESS":
            del self.jobs[job_id]
        else:
            raise InvalidStateTransitionException()

    def cancel_job(
        self, job_id: str, reason_code: str, comment: str, force: bool
    ) -> FakeJob:
        job = self.jobs[job_id]

        job.reason_code = reason_code if reason_code is not None else job.reason_code
        job.comment = comment if comment is not None else job.comment
        job.force = force if force is not None and force != job.force else job.force
        job.status = "CANCELED"

        if job.status == "IN_PROGRESS" and force:
            self.jobs[job_id] = job
        elif job.status != "IN_PROGRESS":
            self.jobs[job_id] = job
        else:
            raise InvalidStateTransitionException()

        return job

    def get_job_document(self, job_id: str) -> FakeJob:
        return self.jobs[job_id]

    def list_jobs(
        self, max_results: int, token: Optional[str]
    ) -> Tuple[List[Dict[str, Any]], Optional[str]]:
        """
        The following parameter are not yet implemented: Status, TargetSelection, ThingGroupName, ThingGroupId
        """
        all_jobs = [_.to_dict() for _ in self.jobs.values()]
        filtered_jobs = all_jobs

        if token is None:
            jobs = filtered_jobs[0:max_results]
            next_token = str(max_results) if len(filtered_jobs) > max_results else None
        else:
            int_token = int(token)
            jobs = filtered_jobs[int_token : int_token + max_results]
            next_token = (
                str(int_token + max_results)
                if len(filtered_jobs) > int_token + max_results
                else None
            )

        return jobs, next_token

    def describe_job_execution(
        self, job_id: str, thing_name: str, execution_number: int
    ) -> FakeJobExecution:
        try:
            job_execution = self.job_executions[(job_id, thing_name)]
        except KeyError:
            raise ResourceNotFoundException()

        if job_execution is None or (
            execution_number is not None
            and job_execution.execution_number != execution_number
        ):
            raise ResourceNotFoundException()

        return job_execution

    def cancel_job_execution(self, job_id: str, thing_name: str, force: bool) -> None:
        """
        The parameters ExpectedVersion and StatusDetails are not yet implemented
        """
        job_execution = self.job_executions[(job_id, thing_name)]

        if job_execution is None:
            raise ResourceNotFoundException()

        job_execution.force_canceled = (
            force if force is not None else job_execution.force_canceled
        )
        # TODO: implement expected_version and status_details (at most 10 can be specified)

        if job_execution.status == "IN_PROGRESS" and force:
            job_execution.status = "CANCELED"
            self.job_executions[(job_id, thing_name)] = job_execution
        elif job_execution.status != "IN_PROGRESS":
            job_execution.status = "CANCELED"
            self.job_executions[(job_id, thing_name)] = job_execution
        else:
            raise InvalidStateTransitionException()

    def delete_job_execution(
        self, job_id: str, thing_name: str, execution_number: int, force: bool
    ) -> None:
        job_execution = self.job_executions[(job_id, thing_name)]

        if job_execution.execution_number != execution_number:
            raise ResourceNotFoundException()

        if job_execution.status == "IN_PROGRESS" and force:
            del self.job_executions[(job_id, thing_name)]
        elif job_execution.status != "IN_PROGRESS":
            del self.job_executions[(job_id, thing_name)]
        else:
            raise InvalidStateTransitionException()

    def list_job_executions_for_job(
        self, job_id: str, status: str, max_results: int, token: Optional[str]
    ) -> Tuple[List[Dict[str, Any]], Optional[str]]:
        job_executions = [
            self.job_executions[je].to_dict()
            for je in self.job_executions
            if je[0] == job_id
        ]

        if status is not None:
            job_executions = list(
                filter(
                    lambda elem: elem["jobExecutionSummary"].get("status") == status,
                    job_executions,
                )
            )

        if token is None:
            job_executions = job_executions[0:max_results]
            next_token = str(max_results) if len(job_executions) > max_results else None
        else:
            int_token = int(token)
            job_executions = job_executions[int_token : int_token + max_results]
            next_token = (
                str(int_token + max_results)
                if len(job_executions) > int_token + max_results
                else None
            )

        return job_executions, next_token

    @paginate(PAGINATION_MODEL)  # type: ignore[misc]
    def list_job_executions_for_thing(
        self, thing_name: str, status: Optional[str]
    ) -> List[Dict[str, Any]]:
        job_executions = [
            self.job_executions[je].to_dict()
            for je in self.job_executions
            if je[1] == thing_name
        ]

        if status is not None:
            job_executions = list(
                filter(
                    lambda elem: elem["jobExecutionSummary"].get("status") == status,
                    job_executions,
                )
            )

        return job_executions

    def list_topic_rules(self) -> List[Dict[str, Any]]:
        return [r.to_dict() for r in self.rules.values()]

    def get_topic_rule(self, rule_name: str) -> Dict[str, Any]:
        if rule_name not in self.rules:
            raise ResourceNotFoundException()
        return self.rules[rule_name].to_get_dict()

    def create_topic_rule(self, rule_name: str, sql: str, **kwargs: Any) -> None:
        if not re.match("^[a-zA-Z0-9_]+$", rule_name):
            msg = f"1 validation error detected: Value '{rule_name}' at 'ruleName' failed to satisfy constraint: Member must satisfy regular expression pattern: ^[a-zA-Z0-9_]+$"
            raise InvalidRequestException(msg)
        if rule_name in self.rules:
            raise ResourceAlreadyExistsException(
                "Rule with given name already exists", "", self.rules[rule_name].arn
            )
        result = re.search(r"FROM\s+([^\s]*)", sql)
        topic = result.group(1).strip("'") if result else None
        self.rules[rule_name] = FakeRule(
            rule_name=rule_name,
            created_at=int(time.time()),
            topic_pattern=topic,
            sql=sql,
            account_id=self.account_id,
            region_name=self.region_name,
            **kwargs,
        )

    def replace_topic_rule(self, rule_name: str, **kwargs: Any) -> None:
        self.delete_topic_rule(rule_name)
        self.create_topic_rule(rule_name, **kwargs)

    def delete_topic_rule(self, rule_name: str) -> None:
        if rule_name not in self.rules:
            raise ResourceNotFoundException()
        del self.rules[rule_name]

    def enable_topic_rule(self, rule_name: str) -> None:
        if rule_name not in self.rules:
            raise ResourceNotFoundException()
        self.rules[rule_name].rule_disabled = False

    def disable_topic_rule(self, rule_name: str) -> None:
        if rule_name not in self.rules:
            raise ResourceNotFoundException()
        self.rules[rule_name].rule_disabled = True

    def create_domain_configuration(
        self,
        domain_configuration_name: str,
        domain_name: str,
        server_certificate_arns: List[str],
        authorizer_config: Dict[str, Any],
        service_type: str,
    ) -> FakeDomainConfiguration:
        """
        The ValidationCertificateArn-parameter is not yet implemented
        """
        if domain_configuration_name in self.domain_configurations:
            raise ResourceAlreadyExistsException(
                "Domain configuration with given name already exists.",
                self.domain_configurations[
                    domain_configuration_name
                ].domain_configuration_name,
                self.domain_configurations[
                    domain_configuration_name
                ].domain_configuration_arn,
            )
        self.domain_configurations[domain_configuration_name] = FakeDomainConfiguration(
            self.account_id,
            self.region_name,
            domain_configuration_name,
            domain_name,
            server_certificate_arns,
            "ENABLED",
            service_type,
            authorizer_config,
            "CUSTOMER_MANAGED",
        )
        return self.domain_configurations[domain_configuration_name]

    def delete_domain_configuration(self, domain_configuration_name: str) -> None:
        if domain_configuration_name not in self.domain_configurations:
            raise ResourceNotFoundException("The specified resource does not exist.")
        del self.domain_configurations[domain_configuration_name]

    def describe_domain_configuration(
        self, domain_configuration_name: str
    ) -> FakeDomainConfiguration:
        if domain_configuration_name not in self.domain_configurations:
            raise ResourceNotFoundException("The specified resource does not exist.")
        return self.domain_configurations[domain_configuration_name]

    def list_domain_configurations(self) -> List[Dict[str, Any]]:
        return [_.to_dict() for _ in self.domain_configurations.values()]

    def update_domain_configuration(
        self,
        domain_configuration_name: str,
        authorizer_config: Dict[str, Any],
        domain_configuration_status: str,
        remove_authorizer_config: Optional[bool],
    ) -> FakeDomainConfiguration:
        if domain_configuration_name not in self.domain_configurations:
            raise ResourceNotFoundException("The specified resource does not exist.")
        domain_configuration = self.domain_configurations[domain_configuration_name]
        if authorizer_config is not None:
            domain_configuration.authorizer_config = authorizer_config
        if domain_configuration_status is not None:
            domain_configuration.domain_configuration_status = (
                domain_configuration_status
            )
        if remove_authorizer_config is not None and remove_authorizer_config is True:
            domain_configuration.authorizer_config = None
        return domain_configuration

    def search_index(self, query_string: str) -> List[Dict[str, Any]]:
        """
        Pagination is not yet implemented. Only basic search queries are supported for now.
        """
        things = [
            thing for thing in self.things.values() if thing.matches(query_string)
        ]
        return [t.to_dict(include_connectivity=True) for t in things]


iot_backends = BackendDict(IoTBackend, "iot")
