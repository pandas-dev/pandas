import base64
import bz2
import codecs
import copy
import datetime
import gzip
import itertools
import json
import os
import string
import sys
import tempfile
import threading
import urllib.parse
from bisect import insort
from importlib import reload
from io import BytesIO
from typing import Any, Dict, Iterator, List, Optional, Set, Tuple, Union

from moto.cloudwatch.models import MetricDatum
from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import (
    BaseModel,
    CloudFormationModel,
    CloudWatchMetricProvider,
)
from moto.core.utils import (
    iso_8601_datetime_without_milliseconds_s3,
    rfc_1123_datetime,
    unix_time,
    unix_time_millis,
    utcnow,
)
from moto.moto_api._internal import mock_random as random
from moto.moto_api._internal.managed_state_model import ManagedState
from moto.s3.exceptions import (
    AccessDeniedByLock,
    BadRequest,
    BucketAlreadyExists,
    BucketNeedsToBeNew,
    CopyObjectMustChangeSomething,
    CrossLocationLoggingProhibitted,
    DaysMustNotProvidedForSelectRequest,
    DaysMustProvidedExceptForSelectRequest,
    EntityTooSmall,
    HeadOnDeleteMarker,
    InvalidBucketName,
    InvalidNotificationDestination,
    InvalidNotificationEvent,
    InvalidObjectState,
    InvalidPart,
    InvalidPublicAccessBlockConfiguration,
    InvalidRequest,
    InvalidStorageClass,
    InvalidTagError,
    InvalidTargetBucketForLogging,
    MalformedXML,
    MissingBucket,
    MissingKey,
    NoSuchPublicAccessBlockConfiguration,
    NoSuchUpload,
    ObjectLockConfigurationNotFoundError,
)
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import PARTITION_NAMES, LowercaseDict, get_partition, md5_hash

from ..events.notifications import send_notification as events_send_notification
from ..settings import (
    S3_UPLOAD_PART_MIN_SIZE,
    get_s3_default_key_buffer_size,
    s3_allow_crossdomain_access,
)
from . import notifications
from .cloud_formation import cfn_to_api_encryption, is_replacement_update
from .select_object_content import parse_query
from .utils import (
    ARCHIVE_STORAGE_CLASSES,
    LOGGING_SERVICE_PRINCIPAL,
    STORAGE_CLASS,
    CaseInsensitiveDict,
    _VersionedKeyStore,
    compute_checksum,
)

MAX_BUCKET_NAME_LENGTH = 63
MIN_BUCKET_NAME_LENGTH = 3
UPLOAD_ID_BYTES = 43
DEFAULT_TEXT_ENCODING = sys.getdefaultencoding()
OWNER = "75aa57f09aa0c8caeab4f8c24e99d10f8e7faeebf76c078efc7c6caea54ba06a"


class FakeDeleteMarker(BaseModel):
    def __init__(self, key: "FakeKey"):
        self.key = key
        self.name = key.name
        self.last_modified = utcnow()
        self._version_id = str(random.uuid4())

    @property
    def last_modified_ISO8601(self) -> str:
        return iso_8601_datetime_without_milliseconds_s3(self.last_modified)  # type: ignore

    @property
    def version_id(self) -> str:
        return self._version_id


class FakeKey(BaseModel, ManagedState):
    def __init__(
        self,
        name: str,
        value: bytes,
        account_id: str,
        region_name: str,
        storage: Optional[str] = "STANDARD",
        etag: Optional[str] = None,
        is_versioned: bool = False,
        version_id: Optional[str] = None,
        max_buffer_size: Optional[int] = None,
        multipart: Optional["FakeMultipart"] = None,
        bucket_name: Optional[str] = None,
        encryption: Optional[str] = None,
        kms_key_id: Optional[str] = None,
        bucket_key_enabled: Any = None,
        lock_mode: Optional[str] = None,
        lock_legal_status: Optional[str] = None,
        lock_until: Optional[str] = None,
        checksum_value: Optional[str] = None,
    ):
        ManagedState.__init__(
            self,
            "s3::keyrestore",
            transitions=[
                (None, "IN_PROGRESS"),
                ("IN_PROGRESS", "RESTORED"),
            ],
        )
        self.name = name
        self.account_id = account_id
        self.region_name = region_name
        self.partition = get_partition(region_name)
        self.last_modified = utcnow()
        self.acl: Optional[FakeAcl] = get_canned_acl("private")
        self.website_redirect_location: Optional[str] = None
        self.checksum_algorithm = None
        self._storage_class: Optional[str] = storage if storage else "STANDARD"
        self._metadata = LowercaseDict()
        self._expiry: Optional[datetime.datetime] = None
        self._etag = etag
        self._version_id = version_id
        self._is_versioned = is_versioned
        self.multipart = multipart
        self.bucket_name = bucket_name

        self._max_buffer_size = (
            max_buffer_size if max_buffer_size else get_s3_default_key_buffer_size()
        )
        self._value_buffer = tempfile.SpooledTemporaryFile(self._max_buffer_size)
        self.disposed = False
        self.value = value  # type: ignore
        self.lock = threading.Lock()

        self.encryption = encryption
        self.kms_key_id = kms_key_id
        self.bucket_key_enabled = bucket_key_enabled

        self.lock_mode = lock_mode
        self.lock_legal_status = lock_legal_status
        self.lock_until = lock_until
        self.checksum_value = checksum_value

        # Default metadata values
        self._metadata["Content-Type"] = "binary/octet-stream"

    def safe_name(self, encoding_type: Optional[str] = None) -> str:
        if encoding_type == "url":
            return urllib.parse.quote(self.name)
        return self.name

    @property
    def version_id(self) -> str:
        return self._version_id or "null"

    @property
    def value(self) -> bytes:
        with self.lock:
            self._value_buffer.seek(0)
            r = self._value_buffer.read()
            r = copy.copy(r)
            return r

    @property
    def arn(self) -> str:
        # S3 Objects don't have an ARN, but we do need something unique when creating tags against this resource
        return f"arn:{self.partition}:s3:::{self.bucket_name}/{self.name}/{self.version_id}"

    @value.setter  # type: ignore
    def value(self, new_value: bytes) -> None:
        self._value_buffer.seek(0)
        self._value_buffer.truncate()

        # Hack for working around moto's own unit tests; this probably won't
        # actually get hit in normal use.
        if isinstance(new_value, str):
            new_value = new_value.encode(DEFAULT_TEXT_ENCODING)
        self._value_buffer.write(new_value)
        self.contentsize = len(new_value)

    @property
    def status(self) -> Optional[str]:
        previous = self._status
        new_status = super().status
        if previous != "RESTORED" and new_status == "RESTORED":
            s3_backend = s3_backends[self.account_id][self.partition]
            bucket = s3_backend.get_bucket(self.bucket_name)  # type: ignore
            notifications.send_event(
                self.account_id,
                notifications.S3NotificationEvent.OBJECT_RESTORE_COMPLETED_EVENT,
                bucket,
                key=self,
            )
        return new_status

    @status.setter
    def status(self, value: str) -> None:
        self._status = value

    def set_metadata(self, metadata: Any, replace: bool = False) -> None:
        if replace:
            self._metadata = {}  # type: ignore
        self._metadata.update(metadata)

    def set_storage_class(self, storage: Optional[str]) -> None:
        if storage is not None and storage not in STORAGE_CLASS:
            raise InvalidStorageClass(storage=storage)
        self._storage_class = storage

    def set_expiry(self, expiry: Optional[datetime.datetime]) -> None:
        self._expiry = expiry

    def set_acl(self, acl: Optional["FakeAcl"]) -> None:
        self.acl = acl

    def restore(self, days: int) -> None:
        self._expiry = utcnow() + datetime.timedelta(days)
        s3_backend = s3_backends[self.account_id][self.partition]
        bucket = s3_backend.get_bucket(self.bucket_name)  # type: ignore
        notifications.send_event(
            self.account_id,
            notifications.S3NotificationEvent.OBJECT_RESTORE_POST_EVENT,
            bucket,
            key=self,
        )

    @property
    def etag(self) -> str:
        if self._etag is None:
            value_md5 = md5_hash()
            self._value_buffer.seek(0)
            while True:
                block = self._value_buffer.read(16 * 1024 * 1024)  # read in 16MB chunks
                if not block:
                    break
                value_md5.update(block)

            self._etag = value_md5.hexdigest()
        return f'"{self._etag}"'

    @property
    def last_modified_ISO8601(self) -> str:
        return iso_8601_datetime_without_milliseconds_s3(self.last_modified)  # type: ignore

    @property
    def last_modified_RFC1123(self) -> str:
        # Different datetime formats depending on how the key is obtained
        # https://github.com/boto/boto/issues/466
        return rfc_1123_datetime(self.last_modified)

    @property
    def metadata(self) -> LowercaseDict:
        return self._metadata

    @property
    def response_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        res: Dict[str, Any] = {
            "ETag": self.etag,
            "last-modified": self.last_modified_RFC1123,
            "content-length": str(self.size),
        }
        if self.encryption is not None:
            res["x-amz-server-side-encryption"] = self.encryption
            if self.encryption == "aws:kms" and self.kms_key_id is not None:
                res["x-amz-server-side-encryption-aws-kms-key-id"] = self.kms_key_id
        if self.encryption == "aws:kms" and self.bucket_key_enabled is not None:
            res["x-amz-server-side-encryption-bucket-key-enabled"] = (
                self.bucket_key_enabled
            )
        if self._storage_class != "STANDARD":
            res["x-amz-storage-class"] = self._storage_class
        if self._expiry is not None:
            if self.status == "IN_PROGRESS":
                header = 'ongoing-request="true"'
            else:
                header = f'ongoing-request="false", expiry-date="{self.expiry_date}"'
            res["x-amz-restore"] = header

        if self._is_versioned:
            res["x-amz-version-id"] = str(self.version_id)

        if self.checksum_algorithm is not None:
            res["x-amz-sdk-checksum-algorithm"] = self.checksum_algorithm
        if self.website_redirect_location:
            res["x-amz-website-redirect-location"] = self.website_redirect_location
        if self.lock_legal_status:
            res["x-amz-object-lock-legal-hold"] = self.lock_legal_status
        if self.lock_until:
            res["x-amz-object-lock-retain-until-date"] = self.lock_until
        if self.lock_mode:
            res["x-amz-object-lock-mode"] = self.lock_mode

        if self.lock_legal_status:
            res["x-amz-object-lock-legal-hold"] = self.lock_legal_status
        if self.lock_until:
            res["x-amz-object-lock-retain-until-date"] = self.lock_until
        if self.lock_mode:
            res["x-amz-object-lock-mode"] = self.lock_mode

        backend = s3_backends[self.account_id][self.partition]
        tags = backend.tagger.get_tag_dict_for_resource(self.arn)
        if tags:
            res["x-amz-tagging-count"] = str(len(tags.keys()))

        return res

    @property
    def size(self) -> int:
        return self.contentsize

    @property
    def storage_class(self) -> Optional[str]:
        return self._storage_class

    @property
    def expiry_date(self) -> Optional[str]:
        if self._expiry is not None:
            return self._expiry.strftime("%a, %d %b %Y %H:%M:%S GMT")
        return None

    # Keys need to be pickleable due to some implementation details of boto3.
    # Since file objects aren't pickleable, we need to override the default
    # behavior. The following is adapted from the Python docs:
    # https://docs.python.org/3/library/pickle.html#handling-stateful-objects
    def __getstate__(self) -> Dict[str, Any]:
        state = self.__dict__.copy()
        try:
            state["value"] = self.value
        except ValueError:
            # Buffer is already closed, so we can't reach the data
            # Only happens if the key was deleted
            state["value"] = ""
        del state["_value_buffer"]
        del state["lock"]
        return state

    def __setstate__(self, state: Dict[str, Any]) -> None:
        self.__dict__.update({k: v for k, v in state.items() if k != "value"})

        self._value_buffer = tempfile.SpooledTemporaryFile(
            max_size=self._max_buffer_size
        )
        self.value = state["value"]  # type: ignore
        self.lock = threading.Lock()

    @property
    def is_locked(self) -> bool:
        if self.lock_legal_status == "ON":
            return True

        if self.lock_mode == "COMPLIANCE":
            now = utcnow()
            try:
                until = datetime.datetime.strptime(
                    self.lock_until,  # type: ignore
                    "%Y-%m-%dT%H:%M:%SZ",
                )
            except ValueError:
                until = datetime.datetime.strptime(
                    self.lock_until,  # type: ignore
                    "%Y-%m-%dT%H:%M:%S.%fZ",
                )

            if until > now:
                return True

        return False

    def dispose(self, garbage: bool = False) -> None:
        if garbage and not self.disposed:
            import warnings

            warnings.warn("S3 key was not disposed of in time", ResourceWarning)
        try:
            self._value_buffer.close()
            if self.multipart:
                self.multipart.dispose()
        except:  # noqa: E722 Do not use bare except
            pass
        self.disposed = True

    def __del__(self) -> None:
        self.dispose(garbage=True)


class FakeMultipart(BaseModel):
    def __init__(
        self,
        key_name: str,
        metadata: CaseInsensitiveDict,  # type: ignore
        account_id: str,
        region_name: str,
        storage: Optional[str] = None,
        tags: Optional[Dict[str, str]] = None,
        acl: Optional["FakeAcl"] = None,
        sse_encryption: Optional[str] = None,
        kms_key_id: Optional[str] = None,
    ):
        self.key_name = key_name
        self.metadata = metadata
        self.account_id = account_id
        self.region_name = region_name
        self.storage = storage
        self.tags = tags
        self.acl = acl
        self.parts: Dict[int, FakeKey] = {}
        self.partlist: List[int] = []  # ordered list of part ID's
        rand_b64 = base64.b64encode(os.urandom(UPLOAD_ID_BYTES))
        self.id = (
            rand_b64.decode("utf-8").replace("=", "").replace("+", "").replace("/", "")
        )
        self.sse_encryption = sse_encryption
        self.kms_key_id = kms_key_id

    def complete(
        self, body: Iterator[Tuple[int, str]]
    ) -> Tuple[bytes, str, Optional[str]]:
        checksum_algo = self.metadata.get("x-amz-checksum-algorithm")
        decode_hex = codecs.getdecoder("hex_codec")
        total = bytearray()
        md5s = bytearray()
        checksum = bytearray()

        last = None
        count = 0
        for pn, etag in body:
            part = self.parts.get(pn)
            part_etag = None
            if part is not None:
                part_etag = part.etag.replace('"', "")
                etag = etag.replace('"', "")
            if part is None or part_etag != etag:
                raise InvalidPart()
            if last is not None and last.contentsize < S3_UPLOAD_PART_MIN_SIZE:
                raise EntityTooSmall()
            md5s.extend(decode_hex(part_etag)[0])  # type: ignore
            total.extend(part.value)
            if checksum_algo:
                checksum.extend(
                    compute_checksum(part.value, checksum_algo, encode_base64=False)
                )
            last = part
            count += 1

        if count == 0:
            raise MalformedXML

        full_etag = md5_hash()
        full_etag.update(bytes(md5s))
        if checksum_algo:
            encoded_checksum = compute_checksum(checksum, checksum_algo).decode("utf-8")
        else:
            encoded_checksum = None
        return total, f"{full_etag.hexdigest()}-{count}", encoded_checksum

    def set_part(self, part_id: int, value: bytes) -> FakeKey:
        if part_id < 1:
            raise NoSuchUpload(upload_id=part_id)

        key = FakeKey(
            part_id,  # type: ignore
            value,
            account_id=self.account_id,
            region_name=self.region_name,
            encryption=self.sse_encryption,
            kms_key_id=self.kms_key_id,
        )
        if part_id in self.parts:
            # We're overwriting the current part - dispose of it first
            self.parts[part_id].dispose()
        self.parts[part_id] = key
        if part_id not in self.partlist:
            insort(self.partlist, part_id)
        return key

    def list_parts(self, part_number_marker: int, max_parts: int) -> Iterator[FakeKey]:
        max_marker = part_number_marker + max_parts
        for part_id in self.partlist[part_number_marker:max_marker]:
            yield self.parts[part_id]

    def dispose(self) -> None:
        for part in self.parts.values():
            part.dispose()


class FakeGrantee(BaseModel):
    def __init__(self, grantee_id: str = "", uri: str = "", display_name: str = ""):
        self.id = grantee_id
        self.uri = uri
        self.display_name = display_name

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, FakeGrantee):
            return False
        return (
            self.id == other.id
            and self.uri == other.uri
            and self.display_name == other.display_name
        )

    @property
    def type(self) -> str:
        return "Group" if self.uri else "CanonicalUser"

    def __repr__(self) -> str:
        return f"FakeGrantee(display_name: '{self.display_name}', id: '{self.id}', uri: '{self.uri}')"


ALL_USERS_GRANTEE = FakeGrantee(uri="http://acs.amazonaws.com/groups/global/AllUsers")
AUTHENTICATED_USERS_GRANTEE = FakeGrantee(
    uri="http://acs.amazonaws.com/groups/global/AuthenticatedUsers"
)
LOG_DELIVERY_GRANTEE = FakeGrantee(uri="http://acs.amazonaws.com/groups/s3/LogDelivery")

PERMISSION_FULL_CONTROL = "FULL_CONTROL"
PERMISSION_WRITE = "WRITE"
PERMISSION_READ = "READ"
PERMISSION_WRITE_ACP = "WRITE_ACP"
PERMISSION_READ_ACP = "READ_ACP"

CAMEL_CASED_PERMISSIONS = {
    "FULL_CONTROL": "FullControl",
    "WRITE": "Write",
    "READ": "Read",
    "WRITE_ACP": "WriteAcp",
    "READ_ACP": "ReadAcp",
}


class FakeGrant(BaseModel):
    def __init__(self, grantees: List[FakeGrantee], permissions: List[str]):
        self.grantees = grantees
        self.permissions = permissions

    def __repr__(self) -> str:
        return f"FakeGrant(grantees: {self.grantees}, permissions: {self.permissions})"


class FakeAcl(BaseModel):
    def __init__(self, grants: Optional[List[FakeGrant]] = None):
        self.grants = grants or []

    @property
    def public_read(self) -> bool:
        for grant in self.grants:
            if ALL_USERS_GRANTEE in grant.grantees:
                if PERMISSION_READ in grant.permissions:
                    return True
                if PERMISSION_FULL_CONTROL in grant.permissions:
                    return True
        return False

    def __repr__(self) -> str:
        return f"FakeAcl(grants: {self.grants})"

    def to_config_dict(self) -> Dict[str, Any]:
        """Returns the object into the format expected by AWS Config"""
        data: Dict[str, Any] = {
            "grantSet": None,  # Always setting this to None. Feel free to change.
            "owner": {"displayName": None, "id": OWNER},
        }

        # Add details for each Grant:
        grant_list = []
        for grant in self.grants:
            permissions = (
                grant.permissions
                if isinstance(grant.permissions, list)
                else [grant.permissions]  # type: ignore
            )
            for permission in permissions:
                for grantee in grant.grantees:
                    if grantee.uri:
                        grant_list.append(
                            {
                                "grantee": grantee.uri.split(
                                    "http://acs.amazonaws.com/groups/s3/"
                                )[1],
                                "permission": CAMEL_CASED_PERMISSIONS[permission],
                            }
                        )
                    else:
                        grant_list.append(
                            {
                                "grantee": {  # type: ignore
                                    "id": grantee.id,
                                    "displayName": None
                                    if not grantee.display_name
                                    else grantee.display_name,
                                },
                                "permission": CAMEL_CASED_PERMISSIONS[permission],
                            }
                        )

        if grant_list:
            data["grantList"] = grant_list

        return data


def get_canned_acl(acl: str) -> FakeAcl:
    owner_grantee = FakeGrantee(grantee_id=OWNER)
    grants = [FakeGrant([owner_grantee], [PERMISSION_FULL_CONTROL])]
    if acl == "private":
        pass  # no other permissions
    elif acl == "public-read":
        grants.append(FakeGrant([ALL_USERS_GRANTEE], [PERMISSION_READ]))
    elif acl == "public-read-write":
        grants.append(
            FakeGrant([ALL_USERS_GRANTEE], [PERMISSION_READ, PERMISSION_WRITE])
        )
    elif acl == "authenticated-read":
        grants.append(FakeGrant([AUTHENTICATED_USERS_GRANTEE], [PERMISSION_READ]))
    elif acl == "bucket-owner-read":
        pass  # TODO: bucket owner ACL
    elif acl == "bucket-owner-full-control":
        pass  # TODO: bucket owner ACL
    elif acl == "aws-exec-read":
        pass  # TODO: bucket owner, EC2 Read
    elif acl == "log-delivery-write":
        grants.append(
            FakeGrant([LOG_DELIVERY_GRANTEE], [PERMISSION_READ_ACP, PERMISSION_WRITE])
        )
    else:
        assert False, f"Unknown canned acl: {acl}"
    return FakeAcl(grants=grants)


class LifecycleFilter(BaseModel):
    def __init__(
        self,
        prefix: Optional[str] = None,
        tag: Optional[Tuple[str, str]] = None,
        and_filter: Optional["LifecycleAndFilter"] = None,
    ):
        self.prefix = prefix
        (self.tag_key, self.tag_value) = tag if tag else (None, None)
        self.and_filter = and_filter

    def to_config_dict(self) -> Dict[str, Any]:
        if self.prefix is not None:
            return {
                "predicate": {"type": "LifecyclePrefixPredicate", "prefix": self.prefix}
            }

        elif self.tag_key:
            return {
                "predicate": {
                    "type": "LifecycleTagPredicate",
                    "tag": {"key": self.tag_key, "value": self.tag_value},
                }
            }

        else:
            return {
                "predicate": {
                    "type": "LifecycleAndOperator",
                    "operands": self.and_filter.to_config_dict(),  # type: ignore
                }
            }


class LifecycleAndFilter(BaseModel):
    def __init__(
        self, prefix: Optional[str] = None, tags: Optional[Dict[str, str]] = None
    ):
        self.prefix = prefix
        self.tags = tags or {}

    def to_config_dict(self) -> List[Dict[str, Any]]:
        data: List[Dict[str, Any]] = []

        if self.prefix is not None:
            data.append({"type": "LifecyclePrefixPredicate", "prefix": self.prefix})

        for key, value in self.tags.items():
            data.append(
                {"type": "LifecycleTagPredicate", "tag": {"key": key, "value": value}}
            )

        return data


class LifecycleTransition(BaseModel):
    def __init__(
        self,
        date: Optional[str] = None,
        days: Optional[int] = None,
        storage_class: Optional[str] = None,
    ):
        self.date = date
        self.days = days
        self.storage_class = storage_class

    def to_config_dict(self) -> Dict[str, Any]:
        config: Dict[str, Any] = {}
        if self.date is not None:
            config["date"] = self.date
        if self.days is not None:
            config["days"] = self.days
        if self.storage_class is not None:
            config["storageClass"] = self.storage_class
        return config


class LifeCycleNoncurrentVersionTransition(BaseModel):
    def __init__(
        self, days: int, storage_class: str, newer_versions: Optional[int] = None
    ):
        self.newer_versions = newer_versions
        self.days = days
        self.storage_class = storage_class

    def to_config_dict(self) -> Dict[str, Any]:
        config: Dict[str, Any] = {}
        if self.newer_versions is not None:
            config["newerNoncurrentVersions"] = self.newer_versions
        if self.days is not None:
            config["noncurrentDays"] = self.days
        if self.storage_class is not None:
            config["storageClass"] = self.storage_class
        return config


class LifecycleRule(BaseModel):
    def __init__(
        self,
        rule_id: Optional[str] = None,
        prefix: Optional[str] = None,
        lc_filter: Optional[LifecycleFilter] = None,
        status: Optional[str] = None,
        expiration_days: Optional[str] = None,
        expiration_date: Optional[str] = None,
        transitions: Optional[List[LifecycleTransition]] = None,
        expired_object_delete_marker: Optional[str] = None,
        nve_noncurrent_days: Optional[str] = None,
        noncurrent_version_transitions: Optional[
            List[LifeCycleNoncurrentVersionTransition]
        ] = None,
        aimu_days: Optional[str] = None,
    ):
        self.id = rule_id
        self.prefix = prefix
        self.filter = lc_filter
        self.status = status
        self.expiration_days = expiration_days
        self.expiration_date = expiration_date
        self.transitions = transitions
        self.expired_object_delete_marker = expired_object_delete_marker
        self.nve_noncurrent_days = nve_noncurrent_days
        self.noncurrent_version_transitions = noncurrent_version_transitions
        self.aimu_days = aimu_days

    def to_config_dict(self) -> Dict[str, Any]:
        """Converts the object to the AWS Config data dict.

        :param kwargs:
        :return:
        """

        lifecycle_dict: Dict[str, Any] = {
            "id": self.id,
            "prefix": self.prefix,
            "status": self.status,
            "expirationInDays": int(self.expiration_days)
            if self.expiration_days
            else None,
            "expiredObjectDeleteMarker": self.expired_object_delete_marker,
            "noncurrentVersionExpirationInDays": -1 or int(self.nve_noncurrent_days),  # type: ignore
            "expirationDate": self.expiration_date,
        }

        if self.transitions:
            lifecycle_dict["transitions"] = [
                t.to_config_dict() for t in self.transitions
            ]
        else:
            lifecycle_dict["transitions"] = None

        if self.noncurrent_version_transitions:
            lifecycle_dict["noncurrentVersionTransitions"] = [
                t.to_config_dict() for t in self.noncurrent_version_transitions
            ]
        else:
            lifecycle_dict["noncurrentVersionTransitions"] = None

        if self.aimu_days:
            lifecycle_dict["abortIncompleteMultipartUpload"] = {
                "daysAfterInitiation": self.aimu_days
            }
        else:
            lifecycle_dict["abortIncompleteMultipartUpload"] = None

        # Format the filter:
        if self.prefix is None and self.filter is None:
            lifecycle_dict["filter"] = {"predicate": None}

        elif self.prefix:
            lifecycle_dict["filter"] = None
        else:
            lifecycle_dict["filter"] = self.filter.to_config_dict()  # type: ignore

        return lifecycle_dict


class CorsRule(BaseModel):
    def __init__(
        self,
        allowed_methods: Any,
        allowed_origins: Any,
        allowed_headers: Any = None,
        expose_headers: Any = None,
        max_age_seconds: Any = None,
    ):
        self.allowed_methods = (
            [allowed_methods] if isinstance(allowed_methods, str) else allowed_methods
        )
        self.allowed_origins = (
            [allowed_origins] if isinstance(allowed_origins, str) else allowed_origins
        )
        self.allowed_headers = (
            [allowed_headers] if isinstance(allowed_headers, str) else allowed_headers
        )
        self.exposed_headers = (
            [expose_headers] if isinstance(expose_headers, str) else expose_headers
        )
        self.max_age_seconds = max_age_seconds


class Notification(BaseModel):
    def __init__(
        self,
        arn: str,
        events: List[str],
        filters: Optional[Dict[str, Any]] = None,
        notification_id: Optional[str] = None,
    ):
        self.id = notification_id or "".join(
            random.choice(string.ascii_letters + string.digits) for _ in range(50)
        )
        self.arn = arn
        for event_name in events:
            if not notifications.S3NotificationEvent.is_event_valid(event_name):
                raise InvalidNotificationEvent(event_name)
        self.events = events
        self.filters = filters if filters else {}

    def _event_matches(self, event_name: str) -> bool:
        if event_name in self.events:
            return True
        # s3:ObjectCreated:Put --> s3:ObjectCreated:*
        wildcard = ":".join(event_name.rsplit(":")[0:2]) + ":*"
        if wildcard in self.events:
            return True
        return False

    def _key_matches(self, key_name: str) -> bool:
        if "S3Key" not in self.filters:
            return True
        _filters = {f["Name"]: f["Value"] for f in self.filters["S3Key"]["FilterRule"]}
        prefix_matches = "prefix" not in _filters or key_name.startswith(
            _filters["prefix"]
        )
        suffix_matches = "suffix" not in _filters or key_name.endswith(
            _filters["suffix"]
        )
        return prefix_matches and suffix_matches

    def matches(self, event_name: str, key_name: str) -> bool:
        if self._event_matches(event_name):
            if self._key_matches(key_name):
                return True
        return False

    def to_config_dict(self) -> Dict[str, Any]:
        # Type and ARN will be filled in by NotificationConfiguration's to_config_dict:
        data: Dict[str, Any] = {"events": [event for event in self.events]}

        if self.filters:
            data["filter"] = {
                "s3KeyFilter": {
                    "filterRules": [
                        {"name": fr["Name"], "value": fr["Value"]}
                        for fr in self.filters["S3Key"]["FilterRule"]
                    ]
                }
            }
        else:
            data["filter"] = None

        # Not sure why this is a thing since AWS just seems to return this as filters ¯\_(ツ)_/¯
        data["objectPrefixes"] = []

        return data


class NotificationConfiguration(BaseModel):
    def __init__(
        self,
        topic: Optional[List[Dict[str, Any]]] = None,
        queue: Optional[List[Dict[str, Any]]] = None,
        cloud_function: Optional[List[Dict[str, Any]]] = None,
        event_bridge: Optional[Dict[str, Any]] = None,
    ):
        self.topic = (
            [
                Notification(
                    t["Topic"],
                    t["Event"],
                    filters=t.get("Filter"),
                    notification_id=t.get("Id"),
                )
                for t in topic
            ]
            if topic
            else []
        )
        self.queue = (
            [
                Notification(
                    q["Queue"],
                    q["Event"],
                    filters=q.get("Filter"),
                    notification_id=q.get("Id"),
                )
                for q in queue
            ]
            if queue
            else []
        )
        self.cloud_function = (
            [
                Notification(
                    c["CloudFunction"],
                    c["Event"],
                    filters=c.get("Filter"),
                    notification_id=c.get("Id"),
                )
                for c in cloud_function
            ]
            if cloud_function
            else []
        )
        self.event_bridge = event_bridge

    def to_config_dict(self) -> Dict[str, Any]:
        data: Dict[str, Any] = {"configurations": {}}

        for topic in self.topic:
            topic_config = topic.to_config_dict()
            topic_config["topicARN"] = topic.arn
            topic_config["type"] = "TopicConfiguration"
            data["configurations"][topic.id] = topic_config

        for queue in self.queue:
            queue_config = queue.to_config_dict()
            queue_config["queueARN"] = queue.arn
            queue_config["type"] = "QueueConfiguration"
            data["configurations"][queue.id] = queue_config

        for cloud_function in self.cloud_function:
            cf_config = cloud_function.to_config_dict()
            cf_config["queueARN"] = cloud_function.arn
            cf_config["type"] = "LambdaConfiguration"
            data["configurations"][cloud_function.id] = cf_config

        if self.event_bridge is not None:
            data["configurations"]["EventBridgeConfiguration"] = self.event_bridge
        return data


def convert_str_to_bool(item: Any) -> bool:
    """Converts a boolean string to a boolean value"""
    if isinstance(item, str):
        return item.lower() == "true"

    return False


class PublicAccessBlock(BaseModel):
    def __init__(
        self,
        block_public_acls: Optional[str],
        ignore_public_acls: Optional[str],
        block_public_policy: Optional[str],
        restrict_public_buckets: Optional[str],
    ):
        # The boto XML appears to expect these values to exist as lowercase strings...
        self.block_public_acls = block_public_acls or "false"
        self.ignore_public_acls = ignore_public_acls or "false"
        self.block_public_policy = block_public_policy or "false"
        self.restrict_public_buckets = restrict_public_buckets or "false"

    def to_config_dict(self) -> Dict[str, bool]:
        # Need to make the string values booleans for Config:
        return {
            "blockPublicAcls": convert_str_to_bool(self.block_public_acls),
            "ignorePublicAcls": convert_str_to_bool(self.ignore_public_acls),
            "blockPublicPolicy": convert_str_to_bool(self.block_public_policy),
            "restrictPublicBuckets": convert_str_to_bool(self.restrict_public_buckets),
        }


class MultipartDict(Dict[str, FakeMultipart]):
    def __delitem__(self, key: str) -> None:
        if key in self:
            self[key].dispose()
        super().__delitem__(key)


class FakeBucket(CloudFormationModel):
    def __init__(self, name: str, account_id: str, region_name: str):
        self.name = name
        self.account_id = account_id
        self.region_name = region_name
        self.partition = get_partition(region_name)
        self.keys = _VersionedKeyStore()
        self.multiparts = MultipartDict()
        self.versioning_status: Optional[str] = None
        self.rules: List[LifecycleRule] = []
        self.policy: Optional[bytes] = None
        self.website_configuration: Optional[Dict[str, Any]] = None
        self.acl: Optional[FakeAcl] = get_canned_acl("private")
        self.cors: List[CorsRule] = []
        self.logging: Dict[str, Any] = {}
        self.notification_configuration: Optional[NotificationConfiguration] = None
        self.accelerate_configuration: Optional[str] = None
        self.payer = "BucketOwner"
        self.creation_date = datetime.datetime.now(tz=datetime.timezone.utc)
        self.public_access_block: Optional[PublicAccessBlock] = None
        self.encryption: Optional[Dict[str, Any]] = None
        self.object_lock_enabled = False
        self.default_lock_mode: Optional[str] = ""
        self.default_lock_days: Optional[int] = 0
        self.default_lock_years: Optional[int] = 0
        self.ownership_rule: Optional[Dict[str, Any]] = None
        s3_backends.bucket_accounts[name] = (self.partition, account_id)

    @property
    def location(self) -> str:
        return self.region_name

    @property
    def creation_date_ISO8601(self) -> str:
        return iso_8601_datetime_without_milliseconds_s3(self.creation_date)  # type: ignore

    @property
    def is_versioned(self) -> bool:
        return self.versioning_status == "Enabled"

    def get_permission(self, action: str, resource: str) -> Any:
        from moto.iam.access_control import IAMPolicy, PermissionResult

        if self.policy is None:
            return PermissionResult.NEUTRAL

        iam_policy = IAMPolicy(self.policy.decode())
        return iam_policy.is_action_permitted(action, resource)

    def set_lifecycle(self, rules: List[Dict[str, Any]]) -> None:
        self.rules = []
        for rule in rules:
            # Extract and validate actions from Lifecycle rule
            expiration = rule.get("Expiration")

            transitions_input = rule.get("Transition", [])
            if transitions_input and not isinstance(transitions_input, list):
                transitions_input = [rule.get("Transition")]

            transitions = [
                LifecycleTransition(
                    date=transition.get("Date"),
                    days=transition.get("Days"),
                    storage_class=transition.get("StorageClass"),
                )
                for transition in transitions_input
            ]

            try:
                top_level_prefix = (
                    rule["Prefix"] or ""
                )  # If it's `None` the set to the empty string
            except KeyError:
                top_level_prefix = None

            nve_noncurrent_days = None
            if rule.get("NoncurrentVersionExpiration") is not None:
                if rule["NoncurrentVersionExpiration"].get("NoncurrentDays") is None:
                    raise MalformedXML()
                nve_noncurrent_days = rule["NoncurrentVersionExpiration"][
                    "NoncurrentDays"
                ]

            nv_transitions_input = rule.get("NoncurrentVersionTransition", [])
            if nv_transitions_input and not isinstance(nv_transitions_input, list):
                nv_transitions_input = [rule.get("NoncurrentVersionTransition")]

            noncurrent_version_transitions = []
            for nvt in nv_transitions_input:
                if nvt.get("NoncurrentDays") is None or nvt.get("StorageClass") is None:
                    raise MalformedXML()

                transition = LifeCycleNoncurrentVersionTransition(
                    newer_versions=nvt.get("NewerNoncurrentVersions"),
                    days=nvt.get("NoncurrentDays"),
                    storage_class=nvt.get("StorageClass"),
                )
                noncurrent_version_transitions.append(transition)

            aimu_days = None
            if rule.get("AbortIncompleteMultipartUpload") is not None:
                if (
                    rule["AbortIncompleteMultipartUpload"].get("DaysAfterInitiation")
                    is None
                ):
                    raise MalformedXML()
                aimu_days = rule["AbortIncompleteMultipartUpload"][
                    "DaysAfterInitiation"
                ]

            eodm = None
            if expiration and expiration.get("ExpiredObjectDeleteMarker") is not None:
                # This cannot be set if Date or Days is set:
                if expiration.get("Days") or expiration.get("Date"):
                    raise MalformedXML()
                eodm = expiration["ExpiredObjectDeleteMarker"]

            # Pull out the filter:
            lc_filter = None
            if rule.get("Filter"):
                # Can't have both `Filter` and `Prefix` (need to check for the presence of the key):
                try:
                    # 'Prefix' cannot be outside of a Filter:
                    if rule["Prefix"] or not rule["Prefix"]:
                        raise MalformedXML()
                except KeyError:
                    pass

                filters = 0
                try:
                    prefix_filter = (
                        rule["Filter"]["Prefix"] or ""
                    )  # If it's `None` the set to the empty string
                    filters += 1
                except KeyError:
                    prefix_filter = None

                and_filter = None
                if rule["Filter"].get("And"):
                    filters += 1
                    and_tags = {}
                    if rule["Filter"]["And"].get("Tag"):
                        if not isinstance(rule["Filter"]["And"]["Tag"], list):
                            rule["Filter"]["And"]["Tag"] = [
                                rule["Filter"]["And"]["Tag"]
                            ]

                        for t in rule["Filter"]["And"]["Tag"]:
                            and_tags[t["Key"]] = t.get("Value", "")

                    try:
                        and_prefix = (
                            rule["Filter"]["And"]["Prefix"] or ""
                        )  # If it's `None` then set to the empty string
                    except KeyError:
                        and_prefix = None

                    and_filter = LifecycleAndFilter(prefix=and_prefix, tags=and_tags)

                filter_tag = None
                if rule["Filter"].get("Tag"):
                    filters += 1
                    filter_tag = (
                        rule["Filter"]["Tag"]["Key"],
                        rule["Filter"]["Tag"].get("Value", ""),
                    )

                # Can't have more than 1 filter:
                if filters > 1:
                    raise MalformedXML()

                lc_filter = LifecycleFilter(
                    prefix=prefix_filter, tag=filter_tag, and_filter=and_filter
                )

            # If no top level prefix and no filter is present, then this is invalid:
            if top_level_prefix is None:
                try:
                    rule["Filter"]
                except KeyError:
                    raise MalformedXML()

            self.rules.append(
                LifecycleRule(
                    rule_id=rule.get("ID"),
                    prefix=top_level_prefix,
                    lc_filter=lc_filter,
                    status=rule["Status"],
                    expiration_days=expiration.get("Days") if expiration else None,
                    expiration_date=expiration.get("Date") if expiration else None,
                    transitions=transitions,
                    expired_object_delete_marker=eodm,
                    nve_noncurrent_days=nve_noncurrent_days,
                    noncurrent_version_transitions=noncurrent_version_transitions,
                    aimu_days=aimu_days,
                )
            )

    def delete_lifecycle(self) -> None:
        self.rules = []

    def set_cors(self, rules: List[Dict[str, Any]]) -> None:
        self.cors = []

        if len(rules) > 100:
            raise MalformedXML()

        for rule in rules:
            assert isinstance(rule["AllowedMethod"], list) or isinstance(
                rule["AllowedMethod"], str
            )
            assert isinstance(rule["AllowedOrigin"], list) or isinstance(
                rule["AllowedOrigin"], str
            )
            assert isinstance(rule.get("AllowedHeader", []), list) or isinstance(
                rule.get("AllowedHeader", ""), str
            )
            assert isinstance(rule.get("ExposeHeader", []), list) or isinstance(
                rule.get("ExposeHeader", ""), str
            )
            assert isinstance(rule.get("MaxAgeSeconds", "0"), str)

            if isinstance(rule["AllowedMethod"], str):
                methods = [rule["AllowedMethod"]]
            else:
                methods = rule["AllowedMethod"]

            for method in methods:
                if method not in ["GET", "PUT", "HEAD", "POST", "DELETE"]:
                    raise InvalidRequest(method)

            self.cors.append(
                CorsRule(
                    rule["AllowedMethod"],
                    rule["AllowedOrigin"],
                    rule.get("AllowedHeader"),
                    rule.get("ExposeHeader"),
                    rule.get("MaxAgeSeconds"),
                )
            )

    def delete_cors(self) -> None:
        self.cors = []

    @staticmethod
    def _log_permissions_enabled_policy(
        target_bucket: "FakeBucket", target_prefix: Optional[str]
    ) -> bool:
        target_bucket_policy = target_bucket.policy
        if target_bucket_policy:
            target_bucket_policy_json = json.loads(target_bucket_policy.decode())
            for stmt in target_bucket_policy_json["Statement"]:
                if (
                    stmt.get("Principal", {}).get("Service")
                    != LOGGING_SERVICE_PRINCIPAL
                ):
                    continue
                if stmt.get("Effect", "") != "Allow":
                    continue
                if "s3:PutObject" not in stmt.get("Action", []):
                    continue
                if (
                    stmt.get("Resource")
                    != f"arn:{target_bucket.partition}:s3:::{target_bucket.name}/{target_prefix if target_prefix else ''}*"
                    and stmt.get("Resource")
                    != f"arn:{target_bucket.partition}:s3:::{target_bucket.name}/*"
                    and stmt.get("Resource")
                    != f"arn:{target_bucket.partition}:s3:::{target_bucket.name}"
                ):
                    continue
                return True

        return False

    @staticmethod
    def _log_permissions_enabled_acl(target_bucket: "FakeBucket") -> bool:
        write = read_acp = False
        for grant in target_bucket.acl.grants:  # type: ignore
            # Must be granted to: http://acs.amazonaws.com/groups/s3/LogDelivery
            for grantee in grant.grantees:
                if grantee.uri == "http://acs.amazonaws.com/groups/s3/LogDelivery":
                    if (
                        "WRITE" in grant.permissions
                        or "FULL_CONTROL" in grant.permissions
                    ):
                        write = True

                    if (
                        "READ_ACP" in grant.permissions
                        or "FULL_CONTROL" in grant.permissions
                    ):
                        read_acp = True
                    break

        return write and read_acp

    def set_logging(
        self, logging_config: Optional[Dict[str, Any]], bucket_backend: "S3Backend"
    ) -> None:
        if not logging_config:
            self.logging = {}
            return

        # Target bucket must exist in the same account (assuming all moto buckets are in the same account):
        target_bucket = bucket_backend.buckets.get(logging_config["TargetBucket"])
        if not target_bucket:
            raise InvalidTargetBucketForLogging(
                "The target bucket for logging does not exist."
            )

        target_prefix = logging_config.get("TargetPrefix", None)
        has_policy_permissions = self._log_permissions_enabled_policy(
            target_bucket=target_bucket, target_prefix=target_prefix
        )
        has_acl_permissions = self._log_permissions_enabled_acl(
            target_bucket=target_bucket
        )
        if not (has_policy_permissions or has_acl_permissions):
            raise InvalidTargetBucketForLogging(
                "You must either provide the necessary permissions to the logging service using a bucket "
                "policy or give the log-delivery group WRITE and READ_ACP permissions to the target bucket"
            )

        # Buckets must also exist within the same region:
        if target_bucket.region_name != self.region_name:
            raise CrossLocationLoggingProhibitted()

        # Checks pass -- set the logging config:
        self.logging = logging_config

    def set_notification_configuration(
        self, notification_config: Optional[Dict[str, Any]]
    ) -> None:
        if not notification_config:
            self.notification_configuration = None
            return

        self.notification_configuration = NotificationConfiguration(
            topic=notification_config.get("TopicConfiguration"),
            queue=notification_config.get("QueueConfiguration"),
            cloud_function=notification_config.get("CloudFunctionConfiguration"),
            event_bridge=notification_config.get("EventBridgeConfiguration"),
        )

        # Validate that the region is correct:
        for thing in ["topic", "queue", "cloud_function"]:
            for t in getattr(self.notification_configuration, thing):
                region = t.arn.split(":")[3]
                if region != self.region_name:
                    raise InvalidNotificationDestination()

        # Send test events so the user can verify these notifications were set correctly
        notifications.send_test_event(account_id=self.account_id, bucket=self)

    def set_accelerate_configuration(self, accelerate_config: str) -> None:
        if self.accelerate_configuration is None and accelerate_config == "Suspended":
            # Cannot "suspend" a not active acceleration. Leaves it undefined
            return

        self.accelerate_configuration = accelerate_config

    @classmethod
    def has_cfn_attr(cls, attr: str) -> bool:
        return attr in [
            "Arn",
            "DomainName",
            "DualStackDomainName",
            "RegionalDomainName",
            "WebsiteURL",
        ]

    def get_cfn_attribute(self, attribute_name: str) -> Any:
        from moto.cloudformation.exceptions import UnformattedGetAttTemplateException

        if attribute_name == "Arn":
            return self.arn
        elif attribute_name == "DomainName":
            return self.domain_name
        elif attribute_name == "DualStackDomainName":
            return self.dual_stack_domain_name
        elif attribute_name == "RegionalDomainName":
            return self.regional_domain_name
        elif attribute_name == "WebsiteURL":
            return self.website_url
        raise UnformattedGetAttTemplateException()

    def set_acl(self, acl: Optional[FakeAcl]) -> None:
        self.acl = acl

    @property
    def arn(self) -> str:
        return f"arn:{self.partition}:s3:::{self.name}"

    @property
    def domain_name(self) -> str:
        return f"{self.name}.s3.amazonaws.com"

    @property
    def dual_stack_domain_name(self) -> str:
        return f"{self.name}.s3.dualstack.{self.region_name}.amazonaws.com"

    @property
    def regional_domain_name(self) -> str:
        return f"{self.name}.s3.{self.region_name}.amazonaws.com"

    @property
    def website_url(self) -> str:
        return f"http://{self.name}.s3-website.{self.region_name}.amazonaws.com"

    @property
    def physical_resource_id(self) -> str:
        return self.name

    @staticmethod
    def cloudformation_name_type() -> str:
        return "BucketName"

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-s3-bucket.html
        return "AWS::S3::Bucket"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "FakeBucket":
        partition = get_partition(region_name)
        backend = s3_backends[account_id][partition]
        bucket = backend.create_bucket(resource_name, region_name)

        properties = cloudformation_json.get("Properties", {})

        if "BucketEncryption" in properties:
            bucket_encryption = cfn_to_api_encryption(properties["BucketEncryption"])
            backend.put_bucket_encryption(
                bucket_name=resource_name, encryption=bucket_encryption
            )

        return bucket

    @classmethod
    def update_from_cloudformation_json(  # type: ignore[misc]
        cls,
        original_resource: Any,
        new_resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> "FakeBucket":
        properties = cloudformation_json["Properties"]

        if is_replacement_update(properties):
            resource_name_property = cls.cloudformation_name_type()
            if resource_name_property not in properties:
                properties[resource_name_property] = new_resource_name
            new_resource = cls.create_from_cloudformation_json(
                properties[resource_name_property],
                cloudformation_json,
                account_id,
                region_name,
            )
            properties[resource_name_property] = original_resource.name
            cls.delete_from_cloudformation_json(
                original_resource.name, cloudformation_json, account_id, region_name
            )
            return new_resource

        else:  # No Interruption
            if "BucketEncryption" in properties:
                bucket_encryption = cfn_to_api_encryption(
                    properties["BucketEncryption"]
                )
                s3_backends[account_id][
                    get_partition(region_name)
                ].put_bucket_encryption(
                    bucket_name=original_resource.name, encryption=bucket_encryption
                )
            return original_resource

    @classmethod
    def delete_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> None:
        s3_backends[account_id][get_partition(region_name)].delete_bucket(resource_name)

    def to_config_dict(self) -> Dict[str, Any]:
        """Return the AWS Config JSON format of this S3 bucket.

        Note: The following features are not implemented and will need to be if you care about them:
        - Bucket Accelerate Configuration
        """
        backend = s3_backends[self.account_id][get_partition(self.region_name)]
        config_dict: Dict[str, Any] = {
            "version": "1.3",
            "configurationItemCaptureTime": str(self.creation_date),
            "configurationItemStatus": "ResourceDiscovered",
            "configurationStateId": str(int(unix_time())),
            "configurationItemMD5Hash": "",
            "arn": self.arn,
            "resourceType": "AWS::S3::Bucket",
            "resourceId": self.name,
            "resourceName": self.name,
            "awsRegion": self.region_name,
            "availabilityZone": "Regional",
            "resourceCreationTime": str(self.creation_date),
            "relatedEvents": [],
            "relationships": [],
            "tags": backend.tagger.get_tag_dict_for_resource(self.arn),
            "configuration": {
                "name": self.name,
                "owner": {"id": OWNER},
                "creationDate": self.creation_date.isoformat(),
            },
        }

        # Make the supplementary configuration:
        # This is a dobule-wrapped JSON for some reason...
        s_config: Dict[str, Any] = {
            "AccessControlList": json.dumps(json.dumps(self.acl.to_config_dict()))  # type: ignore
        }

        if self.public_access_block:
            s_config["PublicAccessBlockConfiguration"] = json.dumps(
                self.public_access_block.to_config_dict()
            )

        # Tagging is special:
        if config_dict["tags"]:
            s_config["BucketTaggingConfiguration"] = json.dumps(
                {"tagSets": [{"tags": config_dict["tags"]}]}
            )

        # TODO implement Accelerate Configuration:
        s_config["BucketAccelerateConfiguration"] = {"status": None}

        if self.rules:
            s_config["BucketLifecycleConfiguration"] = {
                "rules": [rule.to_config_dict() for rule in self.rules]
            }

        s_config["BucketLoggingConfiguration"] = {
            "destinationBucketName": self.logging.get("TargetBucket", None),
            "logFilePrefix": self.logging.get("TargetPrefix", None),
        }

        s_config["BucketPolicy"] = {
            "policyText": self.policy.decode("utf-8") if self.policy else None
        }

        s_config["IsRequesterPaysEnabled"] = (
            "false" if self.payer == "BucketOwner" else "true"
        )

        if self.notification_configuration:
            s_config["BucketNotificationConfiguration"] = (
                self.notification_configuration.to_config_dict()
            )
        else:
            s_config["BucketNotificationConfiguration"] = {"configurations": {}}

        config_dict["supplementaryConfiguration"] = s_config

        return config_dict

    @property
    def has_default_lock(self) -> bool:
        if not self.object_lock_enabled:
            return False

        if self.default_lock_mode:
            return True

        return False

    def default_retention(self) -> str:
        now = utcnow()
        now += datetime.timedelta(self.default_lock_days)  # type: ignore
        now += datetime.timedelta(self.default_lock_years * 365)  # type: ignore
        return now.strftime("%Y-%m-%dT%H:%M:%SZ")


class S3Backend(BaseBackend, CloudWatchMetricProvider):
    """
    Custom S3 endpoints are supported, if you are using a S3-compatible storage solution like Ceph.
    Example usage:

    .. sourcecode:: python

        os.environ["MOTO_S3_CUSTOM_ENDPOINTS"] = "http://custom.internal.endpoint,http://custom.other.endpoint"

        @mock_aws
        def test_my_custom_endpoint():
            boto3.client("s3", endpoint_url="http://custom.internal.endpoint")
            ...

    Note that this only works if the environment variable is set **before** the mock is initialized.

    _-_-_-_

    When using the MultiPart-API manually, the minimum part size is 5MB, just as with AWS. Use the following environment variable to lower this:

    .. sourcecode:: bash

        S3_UPLOAD_PART_MIN_SIZE=256

    _-_-_-_

    When listing objects, the default max-keys value is 1000, just as with AWS. Use the following environment variable to configure this this:

    .. sourcecode:: bash

        MOTO_S3_DEFAULT_MAX_KEYS=256

    _-_-_-_

    CrossAccount access is allowed by default. If you want Moto to throw an AccessDenied-error when accessing a bucket in another account, use this environment variable:

    .. sourcecode:: bash

        MOTO_S3_ALLOW_CROSSACCOUNT_ACCESS=false

    _-_-_-_

    Install `moto[s3crc32c]` if you use the CRC32C algorithm, and absolutely need the correct value. Alternatively, you can install the `crc32c` dependency manually.

    If this dependency is not installed, Moto will fall-back to the CRC32-computation when computing checksums.

    """

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.buckets: Dict[str, FakeBucket] = {}
        self.tagger = TaggingService()
        self._pagination_tokens: Dict[str, str] = {}

    def reset(self) -> None:
        # For every key and multipart, Moto opens a TemporaryFile to write the value of those keys
        # Ensure that these TemporaryFile-objects are closed, and leave no filehandles open
        #
        # First, check all known buckets/keys
        for bucket in self.buckets.values():
            for key in bucket.keys.values():  # type: ignore
                if isinstance(key, FakeKey):
                    key.dispose()
            for part in bucket.multiparts.values():
                part.dispose()
            s3_backends.bucket_accounts.pop(bucket.name, None)
        #
        # Second, go through the list of instances
        # It may contain FakeKeys created earlier, which are no longer tracked
        for mp in FakeMultipart.instances:  # type: ignore
            mp.dispose()
        for key in FakeKey.instances:  # type: ignore
            key.dispose()
        super().reset()

    def log_incoming_request(self, request: Any, bucket_name: str) -> None:
        """
        Process incoming requests
        If the request is made to a bucket with logging enabled, logs will be persisted in the appropriate bucket
        """
        try:
            bucket = self.get_bucket(bucket_name)
            target_bucket = bucket.logging["TargetBucket"]
            prefix = bucket.logging.get("TargetPrefix", "")

            now = datetime.datetime.now()
            file_name = now.strftime(
                f"%Y-%m-%d-%H-%M-%S-{random.get_random_hex(16).upper()}"
            )
            date = now.strftime("%d/%b/%Y:%H:%M:%S +0000")
            source_ip = "0.0.0.0"
            source_iam = "-"  # Can be the user ARN, or empty
            unknown_hex = random.get_random_hex(16)
            source = f"REST.{request.method}.BUCKET"  # REST/CLI/CONSOLE
            key_name = "-"
            path = urllib.parse.urlparse(request.url).path or "-"
            http_line = f"{request.method} {path} HTTP/1.1"
            response = '200 - - 1 2 "-"'
            user_agent = f"{request.headers.get('User-Agent')} prompt/off command/s3api.put-object"
            content = f"{random.get_random_hex(64)} originbucket [{date}] {source_ip} {source_iam} {unknown_hex} {source} {key_name} {http_line} {response} {user_agent} - c29tZSB1bmtub3duIGRhdGE= SigV4 ECDHE-RSA-AES128-GCM-SHA256 AuthHeader {request.url.split('amazonaws.com')[0]}amazonaws.com TLSv1.2 - -"
            self.put_object(target_bucket, prefix + file_name, value=content)  # type: ignore
        except:  # noqa: E722 Do not use bare except
            # log delivery is not guaranteed in AWS, so if anything goes wrong, it's 'safe' to just ignore it
            # Realistically, we should only get here when the bucket does not exist, or logging is not enabled
            pass

    @property
    def _url_module(self) -> Any:  # type: ignore
        # The urls-property can be different depending on env variables
        # Force a reload, to retrieve the correct set of URLs
        import moto.s3.urls as backend_urls_module

        reload(backend_urls_module)
        return backend_urls_module

    @staticmethod
    def default_vpc_endpoint_service(
        service_region: str, zones: List[str]
    ) -> List[Dict[str, str]]:
        """List of dicts representing default VPC endpoints for this service."""
        accesspoint = {
            "AcceptanceRequired": False,
            "AvailabilityZones": zones,
            "BaseEndpointDnsNames": [
                f"accesspoint.s3-global.{service_region}.vpce.amazonaws.com",
            ],
            "ManagesVpcEndpoints": False,
            "Owner": "amazon",
            "PrivateDnsName": "*.accesspoint.s3-global.amazonaws.com",
            "PrivateDnsNameVerificationState": "verified",
            "PrivateDnsNames": [
                {"PrivateDnsName": "*.accesspoint.s3-global.amazonaws.com"}
            ],
            "ServiceId": f"vpce-svc-{BaseBackend.vpce_random_number()}",
            "ServiceName": "com.amazonaws.s3-global.accesspoint",
            "ServiceType": [{"ServiceType": "Interface"}],
            "Tags": [],
            "VpcEndpointPolicySupported": True,
        }
        return (
            BaseBackend.default_vpc_endpoint_service_factory(
                service_region, zones, "s3", "Interface"
            )
            + BaseBackend.default_vpc_endpoint_service_factory(
                service_region, zones, "s3", "Gateway"
            )
            + [accesspoint]
        )

    @classmethod
    def get_cloudwatch_metrics(cls, account_id: str, region: str) -> List[MetricDatum]:
        metrics = []
        for name, bucket in s3_backends[account_id][
            get_partition(region)
        ].buckets.items():
            metrics.append(
                MetricDatum(
                    namespace="AWS/S3",
                    name="BucketSizeBytes",
                    value=bucket.keys.item_size(),
                    dimensions=[
                        {"Name": "StorageType", "Value": "StandardStorage"},
                        {"Name": "BucketName", "Value": name},
                    ],
                    timestamp=datetime.datetime.now(tz=datetime.timezone.utc).replace(
                        hour=0, minute=0, second=0, microsecond=0
                    ),
                    unit="Bytes",
                )
            )
            metrics.append(
                MetricDatum(
                    namespace="AWS/S3",
                    name="NumberOfObjects",
                    value=len(bucket.keys),
                    dimensions=[
                        {"Name": "StorageType", "Value": "AllStorageTypes"},
                        {"Name": "BucketName", "Value": name},
                    ],
                    timestamp=datetime.datetime.now(tz=datetime.timezone.utc).replace(
                        hour=0, minute=0, second=0, microsecond=0
                    ),
                    unit="Count",
                )
            )
        return metrics

    def create_bucket(self, bucket_name: str, region_name: str) -> FakeBucket:
        if bucket_name in s3_backends.bucket_accounts.keys():
            raise BucketAlreadyExists(bucket=bucket_name)
        if not MIN_BUCKET_NAME_LENGTH <= len(bucket_name) <= MAX_BUCKET_NAME_LENGTH:
            raise InvalidBucketName()
        new_bucket = FakeBucket(
            name=bucket_name, account_id=self.account_id, region_name=region_name
        )

        self.buckets[bucket_name] = new_bucket

        notification_detail = {
            "version": "0",
            "bucket": {"name": bucket_name},
            "request-id": "N4N7GDK58NMKJ12R",
            "requester": self.account_id,
            "source-ip-address": "1.2.3.4",
            "reason": "PutObject",
        }
        events_send_notification(
            source="aws.s3",
            event_name="CreateBucket",
            region=region_name,
            resources=[f"arn:{new_bucket.partition}:s3:::{bucket_name}"],
            detail=notification_detail,
        )

        return new_bucket

    def list_buckets(self) -> List[FakeBucket]:
        return list(self.buckets.values())

    def get_bucket(self, bucket_name: str) -> FakeBucket:
        if bucket_name in self.buckets:
            return self.buckets[bucket_name]

        if bucket_name in s3_backends.bucket_accounts:
            if not s3_allow_crossdomain_access():
                raise AccessDeniedByLock
            (partition, account_id) = s3_backends.bucket_accounts[bucket_name]
            return s3_backends[account_id][partition].get_bucket(bucket_name)

        raise MissingBucket(bucket=bucket_name)

    def head_bucket(self, bucket_name: str) -> FakeBucket:
        return self.get_bucket(bucket_name)

    def delete_bucket(self, bucket_name: str) -> Optional[FakeBucket]:
        bucket = self.get_bucket(bucket_name)
        if bucket.keys:
            # Can't delete a bucket with keys
            return None
        else:
            s3_backends.bucket_accounts.pop(bucket_name, None)
            return self.buckets.pop(bucket_name)

    def put_bucket_versioning(self, bucket_name: str, status: str) -> None:
        self.get_bucket(bucket_name).versioning_status = status

    def get_bucket_versioning(self, bucket_name: str) -> Optional[str]:
        return self.get_bucket(bucket_name).versioning_status

    def get_bucket_encryption(self, bucket_name: str) -> Optional[Dict[str, Any]]:
        return self.get_bucket(bucket_name).encryption

    def list_object_versions(
        self,
        bucket_name: str,
        delimiter: Optional[str] = None,
        key_marker: Optional[str] = None,
        max_keys: Optional[int] = 1000,
        prefix: str = "",
        version_id_marker: Optional[str] = None,
    ) -> Tuple[
        List[FakeKey], List[str], List[FakeDeleteMarker], Optional[str], Optional[str]
    ]:
        """
        The default value for the MaxKeys-argument is 100. This can be configured with an environment variable:

        MOTO_S3_DEFAULT_MAX_KEYS=5
        """
        bucket = self.get_bucket(bucket_name)

        common_prefixes: Set[str] = set()
        requested_versions: List[FakeKey] = []
        delete_markers: List[FakeDeleteMarker] = []
        next_key_marker: Optional[str] = None
        next_version_id_marker: Optional[str] = None
        all_versions = list(
            itertools.chain(
                *(
                    copy.deepcopy(version_key)
                    for key, version_key in bucket.keys.iterlists()
                )
            )
        )
        # sort by name, revert last-modified-date
        all_versions.sort(key=lambda r: (r.name, -unix_time_millis(r.last_modified)))
        last_name = None

        skip_versions = True
        last_item_added: Union[None, FakeKey, FakeDeleteMarker] = None
        for version in all_versions:
            # Pagination
            if skip_versions:
                if key_marker is None:
                    # If KeyMarker is not supplied, we do not skip anything
                    skip_versions = False
                elif not version_id_marker:
                    # Only KeyMarker is supplied, and it will be set to the last item of the previous page
                    # We skip all versions with ``name < key_marker``
                    # Because our list is ordered, we keep everything where ``name >= key_marker`` (i.e.: the next page)
                    skip_versions = version.name < key_marker
                    continue
                elif (
                    version.name == key_marker
                    and version.version_id == version_id_marker
                ):
                    # KeyMarker and VersionIdMarker are set to the last item of the previous page
                    # Which means we should still skip the current version
                    # But continue processing all subsequent versions
                    skip_versions = False
                    continue
                else:
                    continue

            name = version.name
            # guaranteed to be sorted - so the first key with this name will be the latest
            version.is_latest = name != last_name
            if version.is_latest:
                last_name = name

            # Filter for keys that start with prefix
            if not name.startswith(prefix):
                continue
            # separate keys that contain the same string between the prefix and the first occurrence of the delimiter
            is_common_prefix = False
            if delimiter and delimiter in name[len(prefix) :]:
                end_of_delimiter = (
                    len(prefix) + name[len(prefix) :].index(delimiter) + len(delimiter)
                )
                prefix_including_delimiter = name[0:end_of_delimiter]

                # Skip already-processed common prefix.
                if prefix_including_delimiter == key_marker:
                    continue

                common_prefixes.add(prefix_including_delimiter)
                name = prefix_including_delimiter
                is_common_prefix = True
            elif last_item_added:
                name = last_item_added.name

            # Only return max_keys items.
            if (
                max_keys is not None
                and len(requested_versions) + len(delete_markers) + len(common_prefixes)
                >= max_keys
            ):
                next_key_marker = name
                if is_common_prefix:
                    # No NextToken when returning common prefixes
                    next_version_id_marker = None
                elif last_item_added is not None:
                    # NextToken is set to the (version of the) latest item
                    next_version_id_marker = last_item_added.version_id
                else:
                    # Should only happen when max_keys == 0, so when we do not have a last item
                    next_version_id_marker = None
                break

            if not is_common_prefix:
                last_item_added = version
                # Differentiate between FakeKey and FakeDeleteMarkers
                if not isinstance(version, FakeKey):
                    delete_markers.append(version)
                else:
                    requested_versions.append(version)

        return (
            requested_versions,
            sorted(common_prefixes),
            delete_markers,
            next_key_marker,
            next_version_id_marker,
        )

    def get_bucket_policy(self, bucket_name: str) -> Optional[bytes]:
        return self.get_bucket(bucket_name).policy

    def put_bucket_policy(self, bucket_name: str, policy: bytes) -> None:
        """
        Basic policy enforcement is in place.

        Restrictions:
         - Only statements with principal=* are taken into account
         - Conditions are not taken into account
        """
        self.get_bucket(bucket_name).policy = policy

    def delete_bucket_policy(self, bucket_name: str) -> None:
        bucket = self.get_bucket(bucket_name)
        bucket.policy = None

    def put_bucket_encryption(
        self, bucket_name: str, encryption: Dict[str, Any]
    ) -> None:
        self.get_bucket(bucket_name).encryption = encryption

    def delete_bucket_encryption(self, bucket_name: str) -> None:
        self.get_bucket(bucket_name).encryption = None

    def get_bucket_ownership_controls(
        self, bucket_name: str
    ) -> Optional[Dict[str, Any]]:
        return self.get_bucket(bucket_name).ownership_rule

    def put_bucket_ownership_controls(
        self, bucket_name: str, ownership: Dict[str, Any]
    ) -> None:
        self.get_bucket(bucket_name).ownership_rule = ownership

    def delete_bucket_ownership_controls(self, bucket_name: str) -> None:
        self.get_bucket(bucket_name).ownership_rule = None

    def get_bucket_replication(self, bucket_name: str) -> Optional[Dict[str, Any]]:
        bucket = self.get_bucket(bucket_name)
        return getattr(bucket, "replication", None)

    def put_bucket_replication(
        self, bucket_name: str, replication: Dict[str, Any]
    ) -> None:
        if isinstance(replication["Rule"], dict):
            replication["Rule"] = [replication["Rule"]]
        for rule in replication["Rule"]:
            if "Priority" not in rule:
                rule["Priority"] = 1
            if "ID" not in rule:
                rule["ID"] = "".join(
                    random.choice(string.ascii_letters + string.digits)
                    for _ in range(30)
                )
        bucket = self.get_bucket(bucket_name)
        bucket.replication = replication  # type: ignore

    def delete_bucket_replication(self, bucket_name: str) -> None:
        bucket = self.get_bucket(bucket_name)
        bucket.replication = None  # type: ignore

    def put_bucket_lifecycle(
        self, bucket_name: str, rules: List[Dict[str, Any]]
    ) -> None:
        bucket = self.get_bucket(bucket_name)
        bucket.set_lifecycle(rules)

    def delete_bucket_lifecycle(self, bucket_name: str) -> None:
        bucket = self.get_bucket(bucket_name)
        bucket.delete_lifecycle()

    def set_bucket_website_configuration(
        self, bucket_name: str, website_configuration: Dict[str, Any]
    ) -> None:
        bucket = self.get_bucket(bucket_name)
        bucket.website_configuration = website_configuration

    def get_bucket_website_configuration(
        self, bucket_name: str
    ) -> Optional[Dict[str, Any]]:
        bucket = self.get_bucket(bucket_name)
        return bucket.website_configuration

    def delete_bucket_website(self, bucket_name: str) -> None:
        bucket = self.get_bucket(bucket_name)
        bucket.website_configuration = None

    def get_public_access_block(self, bucket_name: str) -> PublicAccessBlock:
        bucket = self.get_bucket(bucket_name)

        if not bucket.public_access_block:
            raise NoSuchPublicAccessBlockConfiguration()

        return bucket.public_access_block

    def put_object(
        self,
        bucket_name: str,
        key_name: str,
        value: bytes,
        storage: Optional[str] = None,
        etag: Optional[str] = None,
        multipart: Optional[FakeMultipart] = None,
        encryption: Optional[str] = None,
        kms_key_id: Optional[str] = None,
        bucket_key_enabled: Any = None,
        lock_mode: Optional[str] = None,
        lock_legal_status: Optional[str] = None,
        lock_until: Optional[str] = None,
        checksum_value: Optional[str] = None,
        # arguments to handle notification
        request_method: Optional[str] = "PUT",
        disable_notification: Optional[bool] = False,
    ) -> FakeKey:
        if storage is not None and storage not in STORAGE_CLASS:
            raise InvalidStorageClass(storage=storage)

        bucket = self.get_bucket(bucket_name)

        # getting default config from bucket if not included in put request
        if bucket.encryption:
            bucket_key_enabled = bucket_key_enabled or bucket.encryption["Rule"].get(
                "BucketKeyEnabled", False
            )
            kms_key_id = kms_key_id or bucket.encryption["Rule"][
                "ApplyServerSideEncryptionByDefault"
            ].get("KMSMasterKeyID")
            encryption = (
                encryption
                or bucket.encryption["Rule"]["ApplyServerSideEncryptionByDefault"][
                    "SSEAlgorithm"
                ]
            )

        new_key = FakeKey(
            name=key_name,
            bucket_name=bucket_name,
            value=value,
            account_id=self.account_id,
            region_name=self.region_name,
            storage=storage,
            etag=etag,
            is_versioned=bucket.is_versioned,
            # AWS uses VersionId=null in both requests and responses
            version_id=str(random.uuid4()) if bucket.is_versioned else "null",
            multipart=multipart,
            encryption=encryption,
            kms_key_id=kms_key_id,
            bucket_key_enabled=bucket_key_enabled,
            lock_mode=lock_mode,
            lock_legal_status=lock_legal_status,
            lock_until=lock_until,
            checksum_value=checksum_value,
        )

        existing_keys = bucket.keys.getlist(key_name, [])
        if bucket.is_versioned:
            keys = existing_keys + [new_key]
        else:
            keys = [new_key]
        bucket.keys.setlist(key_name, keys)

        if not disable_notification:
            # Send event notification
            if request_method == "POST":
                notify_event_name = (
                    notifications.S3NotificationEvent.OBJECT_CREATED_POST_EVENT
                )
            else:  # PUT request
                notify_event_name = (
                    notifications.S3NotificationEvent.OBJECT_CREATED_PUT_EVENT
                )
            notifications.send_event(
                self.account_id,
                notify_event_name,
                bucket,
                new_key,
            )

        return new_key

    def put_object_acl(
        self,
        bucket_name: str,
        key_name: str,
        acl: Optional[FakeAcl],
    ) -> None:
        key = self.get_object(bucket_name, key_name)
        # TODO: Support the XML-based ACL format
        if key is not None:
            key.set_acl(acl)
        else:
            raise MissingKey(key=key_name)

    def put_object_legal_hold(
        self,
        bucket_name: str,
        key_name: str,
        version_id: Optional[str],
        legal_hold_status: Dict[str, Any],
    ) -> None:
        key = self.get_object(bucket_name, key_name, version_id=version_id)
        key.lock_legal_status = legal_hold_status  # type: ignore

    def put_object_retention(
        self,
        bucket_name: str,
        key_name: str,
        version_id: Optional[str],
        retention: Tuple[Optional[str], Optional[str]],
    ) -> None:
        key = self.get_object(bucket_name, key_name, version_id=version_id)
        key.lock_mode = retention[0]  # type: ignore
        key.lock_until = retention[1]  # type: ignore

    def get_object_attributes(
        self,
        key: FakeKey,
        attributes_to_get: List[str],
    ) -> Dict[str, Any]:
        """
        The following attributes are not yet returned: DeleteMarker, RequestCharged, ObjectParts
        """
        response_keys: Dict[str, Any] = {
            "etag": None,
            "checksum": None,
            "size": None,
            "storage_class": None,
        }
        if "ETag" in attributes_to_get:
            response_keys["etag"] = key.etag.replace('"', "")
        if "Checksum" in attributes_to_get and key.checksum_value is not None:
            response_keys["checksum"] = {key.checksum_algorithm: key.checksum_value}
        if "ObjectSize" in attributes_to_get:
            response_keys["size"] = key.size
        if "StorageClass" in attributes_to_get:
            response_keys["storage_class"] = key.storage_class
        return response_keys

    def get_object(
        self,
        bucket_name: str,
        key_name: str,
        version_id: Optional[str] = None,
        part_number: Optional[str] = None,
        return_delete_marker: bool = False,
    ) -> Optional[FakeKey]:
        bucket = self.get_bucket(bucket_name)

        key = None

        if bucket:
            if version_id is None:
                if key_name in bucket.keys:
                    key = bucket.keys[key_name]
            else:
                for key_version in bucket.keys.getlist(key_name, default=[]):
                    if str(key_version.version_id) == str(version_id):
                        key = key_version
                        break

            if part_number and key and key.multipart:
                key = key.multipart.parts[part_number]

        if isinstance(key, FakeKey):
            key.advance()
            return key
        else:
            if return_delete_marker and isinstance(key, FakeDeleteMarker):
                return key  # type: ignore
            return None

    def head_object(
        self,
        bucket_name: str,
        key_name: str,
        version_id: Optional[str] = None,
        part_number: Optional[str] = None,
    ) -> Optional[FakeKey]:
        obj = self.get_object(
            bucket_name, key_name, version_id, part_number, return_delete_marker=True
        )
        if isinstance(obj, FakeDeleteMarker):
            raise HeadOnDeleteMarker(obj)
        return obj

    def get_object_acl(self, key: FakeKey) -> Optional[FakeAcl]:
        return key.acl

    def get_object_legal_hold(self, key: FakeKey) -> Optional[str]:
        return key.lock_legal_status

    def get_object_lock_configuration(
        self, bucket_name: str
    ) -> Tuple[bool, Optional[str], Optional[int], Optional[int]]:
        bucket = self.get_bucket(bucket_name)
        if not bucket.object_lock_enabled:
            raise ObjectLockConfigurationNotFoundError
        return (
            bucket.object_lock_enabled,
            bucket.default_lock_mode,
            bucket.default_lock_days,
            bucket.default_lock_years,
        )

    def get_object_tagging(self, key: FakeKey) -> Dict[str, List[Dict[str, str]]]:
        return self.tagger.list_tags_for_resource(key.arn)

    def put_object_tagging(
        self,
        key: Optional[FakeKey],
        tags: Optional[Dict[str, str]],
        key_name: Optional[str] = None,
    ) -> FakeKey:
        if key is None:
            raise MissingKey(key=key_name)

        # get bucket for eventbridge notification
        # we can assume that the key has its bucket
        bucket = self.get_bucket(key.bucket_name)  # type: ignore

        tags_input = self.tagger.convert_dict_to_tags_input(tags)
        # Validation custom to S3
        if tags:
            if len(tags_input) > 10:
                raise BadRequest("Object tags cannot be greater than 10")
            if any([tagkey.startswith("aws") for tagkey in tags.keys()]):
                raise InvalidTagError("Your TagKey cannot be prefixed with aws:")
        # Validation shared across all services
        errmsg = self.tagger.validate_tags(tags_input)
        if errmsg:
            raise InvalidTagError(errmsg)
        self.tagger.delete_all_tags_for_resource(key.arn)
        self.tagger.tag_resource(key.arn, tags_input)
        notifications.send_event(
            self.account_id,
            notifications.S3NotificationEvent.OBJECT_TAGGING_PUT_EVENT,
            bucket,
            key,
        )
        return key

    def get_bucket_tagging(self, bucket_name: str) -> Dict[str, List[Dict[str, str]]]:
        bucket = self.get_bucket(bucket_name)
        return self.tagger.list_tags_for_resource(bucket.arn)

    def put_bucket_tagging(self, bucket_name: str, tags: Dict[str, str]) -> None:
        bucket = self.get_bucket(bucket_name)
        self.tagger.delete_all_tags_for_resource(bucket.arn)
        self.tagger.tag_resource(
            bucket.arn, [{"Key": key, "Value": value} for key, value in tags.items()]
        )

    def put_object_lock_configuration(
        self,
        bucket_name: str,
        lock_enabled: bool,
        mode: Optional[str] = None,
        days: Optional[int] = None,
        years: Optional[int] = None,
    ) -> None:
        bucket = self.get_bucket(bucket_name)

        if bucket.keys.item_size() > 0:
            raise BucketNeedsToBeNew

        if lock_enabled:
            bucket.object_lock_enabled = True
            bucket.versioning_status = "Enabled"

        bucket.default_lock_mode = mode
        bucket.default_lock_days = days
        bucket.default_lock_years = years

    def delete_bucket_tagging(self, bucket_name: str) -> None:
        bucket = self.get_bucket(bucket_name)
        self.tagger.delete_all_tags_for_resource(bucket.arn)

    def put_bucket_cors(
        self, bucket_name: str, cors_rules: List[Dict[str, Any]]
    ) -> None:
        bucket = self.get_bucket(bucket_name)
        bucket.set_cors(cors_rules)

    def put_bucket_logging(
        self, bucket_name: str, logging_config: Dict[str, Any]
    ) -> None:
        bucket = self.get_bucket(bucket_name)
        bucket.set_logging(logging_config, self)

    def delete_bucket_cors(self, bucket_name: str) -> None:
        bucket = self.get_bucket(bucket_name)
        bucket.delete_cors()

    def delete_public_access_block(self, bucket_name: str) -> None:
        bucket = self.get_bucket(bucket_name)
        bucket.public_access_block = None

    def put_bucket_notification_configuration(
        self, bucket_name: str, notification_config: Dict[str, Any]
    ) -> None:
        """
        The configuration can be persisted, but at the moment we only send notifications to the following targets:

         - AWSLambda
         - SNS
         - SQS
         - EventBridge

        For the following events:
         - 's3:ObjectCreated:CompleteMultipartUpload'
         - 's3:ObjectCreated:Copy'
         - 's3:ObjectCreated:Post'
         - 's3:ObjectCreated:Put'
         - 's3:ObjectDeleted'
         - 's3:ObjectRestore:Post'
        """
        bucket = self.get_bucket(bucket_name)
        bucket.set_notification_configuration(notification_config)

    def put_bucket_accelerate_configuration(
        self, bucket_name: str, accelerate_configuration: str
    ) -> None:
        if accelerate_configuration not in ["Enabled", "Suspended"]:
            raise MalformedXML()

        bucket = self.get_bucket(bucket_name)
        if bucket.name.find(".") != -1:
            raise InvalidRequest("PutBucketAccelerateConfiguration")
        bucket.set_accelerate_configuration(accelerate_configuration)

    def put_public_access_block(
        self, bucket_name: str, pub_block_config: Optional[Dict[str, Any]]
    ) -> None:
        bucket = self.get_bucket(bucket_name)

        if not pub_block_config:
            raise InvalidPublicAccessBlockConfiguration()

        bucket.public_access_block = PublicAccessBlock(
            pub_block_config.get("BlockPublicAcls"),
            pub_block_config.get("IgnorePublicAcls"),
            pub_block_config.get("BlockPublicPolicy"),
            pub_block_config.get("RestrictPublicBuckets"),
        )

    def abort_multipart_upload(self, bucket_name: str, multipart_id: str) -> None:
        bucket = self.get_bucket(bucket_name)
        multipart_data = bucket.multiparts.get(multipart_id, None)
        if not multipart_data:
            raise NoSuchUpload(upload_id=multipart_id)
        del bucket.multiparts[multipart_id]

    def list_parts(
        self,
        bucket_name: str,
        multipart_id: str,
        part_number_marker: int = 0,
        max_parts: int = 1000,
    ) -> List[FakeKey]:
        bucket = self.get_bucket(bucket_name)
        if multipart_id not in bucket.multiparts:
            raise NoSuchUpload(upload_id=multipart_id)
        return list(
            bucket.multiparts[multipart_id].list_parts(part_number_marker, max_parts)
        )

    def is_truncated(
        self, bucket_name: str, multipart_id: str, next_part_number_marker: int
    ) -> bool:
        bucket = self.get_bucket(bucket_name)
        return len(bucket.multiparts[multipart_id].parts) > next_part_number_marker

    def create_multipart_upload(
        self,
        bucket_name: str,
        key_name: str,
        metadata: CaseInsensitiveDict,  # type: ignore
        storage_type: str,
        tags: Dict[str, str],
        acl: Optional[FakeAcl],
        sse_encryption: str,
        kms_key_id: str,
    ) -> str:
        multipart = FakeMultipart(
            key_name,
            metadata,
            account_id=self.account_id,
            region_name=self.region_name,
            storage=storage_type,
            tags=tags,
            acl=acl,
            sse_encryption=sse_encryption,
            kms_key_id=kms_key_id,
        )

        bucket = self.get_bucket(bucket_name)
        bucket.multiparts[multipart.id] = multipart
        return multipart.id

    def complete_multipart_upload(
        self, bucket_name: str, multipart_id: str, body: Iterator[Tuple[int, str]]
    ) -> Optional[FakeKey]:
        bucket = self.get_bucket(bucket_name)
        multipart = bucket.multiparts[multipart_id]
        value, etag, checksum = multipart.complete(body)
        if value is not None:
            del bucket.multiparts[multipart_id]

        if value is None:
            return None

        key = self.put_object(
            bucket_name,
            multipart.key_name,
            value,
            storage=multipart.storage,
            etag=etag,
            multipart=multipart,
            encryption=multipart.sse_encryption,
            kms_key_id=multipart.kms_key_id,
        )
        key.set_metadata(multipart.metadata)

        if checksum:
            key.checksum_algorithm = multipart.metadata.get("x-amz-checksum-algorithm")
            key.checksum_value = checksum

        self.put_object_tagging(key, multipart.tags)
        self.put_object_acl(
            bucket_name=bucket_name,
            key_name=key.name,
            acl=multipart.acl,
        )

        notifications.send_event(
            self.account_id,
            notifications.S3NotificationEvent.OBJECT_CREATED_COMPLETE_MULTIPART_UPLOAD_EVENT,
            bucket,
            bucket.keys.get(multipart.key_name),
        )
        return key

    def get_all_multiparts(self, bucket_name: str) -> Dict[str, FakeMultipart]:
        bucket = self.get_bucket(bucket_name)
        return bucket.multiparts

    def upload_part(
        self, bucket_name: str, multipart_id: str, part_id: int, value: bytes
    ) -> FakeKey:
        bucket = self.get_bucket(bucket_name)
        multipart = bucket.multiparts[multipart_id]
        return multipart.set_part(part_id, value)

    def upload_part_copy(
        self,
        dest_bucket_name: str,
        multipart_id: str,
        part_id: int,
        src_bucket_name: str,
        src_key_name: str,
        src_version_id: Optional[str],
        start_byte: int,
        end_byte: int,
    ) -> FakeKey:
        dest_bucket = self.get_bucket(dest_bucket_name)
        multipart = dest_bucket.multiparts[multipart_id]

        src_value = self.get_object(  # type: ignore
            src_bucket_name, src_key_name, version_id=src_version_id
        ).value
        if start_byte is not None:
            src_value = src_value[start_byte : end_byte + 1]
        return multipart.set_part(part_id, src_value)

    def list_objects(
        self,
        bucket: FakeBucket,
        prefix: Optional[str],
        delimiter: Optional[str],
        marker: Optional[str],
        max_keys: Optional[int],
    ) -> Tuple[Set[FakeKey], Set[str], bool, Optional[str]]:
        """
        The default value for the MaxKeys-argument is 100. This can be configured with an environment variable:

        MOTO_S3_DEFAULT_MAX_KEYS=5
        """
        key_results = set()
        folder_results = set()
        if prefix:
            for key_name, key in bucket.keys.items():  # type: ignore
                if key_name.startswith(prefix):
                    key_without_prefix = key_name.replace(prefix, "", 1)
                    if delimiter and delimiter in key_without_prefix:
                        # If delimiter, we need to split out folder_results
                        key_without_delimiter = key_without_prefix.split(delimiter)[0]
                        folder_results.add(
                            f"{prefix}{key_without_delimiter}{delimiter}"
                        )
                    else:
                        key_results.add(key)
        else:
            for key_name, key in bucket.keys.items():  # type: ignore
                if delimiter and delimiter in key_name:
                    # If delimiter, we need to split out folder_results
                    folder_results.add(key_name.split(delimiter)[0] + delimiter)
                else:
                    key_results.add(key)

        key_results = filter(  # type: ignore
            lambda key: not isinstance(key, FakeDeleteMarker), key_results
        )
        key_results = sorted(key_results, key=lambda key: key.name)  # type: ignore
        folder_results = [  # type: ignore
            folder_name for folder_name in sorted(folder_results, key=lambda key: key)
        ]

        if marker:
            limit = self._pagination_tokens.get(marker) or marker
            key_results = self._get_results_from_token(key_results, limit)

        if max_keys is not None:
            key_results, is_truncated, next_marker = self._truncate_result(
                key_results, max_keys
            )
        else:
            is_truncated = False
            next_marker = None

        return key_results, folder_results, is_truncated, next_marker

    def list_objects_v2(
        self,
        bucket: FakeBucket,
        prefix: Optional[str],
        delimiter: Optional[str],
        continuation_token: Optional[str],
        start_after: Optional[str],
        max_keys: int,
    ) -> Tuple[Set[Union[FakeKey, str]], bool, Optional[str]]:
        """
        The default value for the MaxKeys-argument is 100. This can be configured with an environment variable:

        MOTO_S3_DEFAULT_MAX_KEYS=5
        """
        result_keys, result_folders, _, _ = self.list_objects(
            bucket, prefix, delimiter, marker=None, max_keys=None
        )
        # sort the combination of folders and keys into lexicographical order
        all_keys = result_keys + result_folders  # type: ignore
        all_keys.sort(key=self._get_name)

        if continuation_token or start_after:
            limit = (
                self._pagination_tokens.get(continuation_token)
                if continuation_token
                else start_after
            )
            all_keys = self._get_results_from_token(all_keys, limit)

        truncated_keys, is_truncated, next_continuation_token = self._truncate_result(
            all_keys, max_keys
        )

        return truncated_keys, is_truncated, next_continuation_token

    def _get_results_from_token(self, result_keys: Any, token: Any) -> Any:
        continuation_index = 0
        for key in result_keys:
            if (key.name if isinstance(key, FakeKey) else key) > token:
                break
            continuation_index += 1
        return result_keys[continuation_index:]

    def _truncate_result(self, result_keys: Any, max_keys: int) -> Any:
        if max_keys == 0:
            result_keys = []
            is_truncated = True
            next_continuation_token = None
        elif len(result_keys) > max_keys:
            is_truncated = "true"  # type: ignore
            result_keys = result_keys[:max_keys]
            item = result_keys[-1]
            key_id = item.name if isinstance(item, FakeKey) else item
            next_continuation_token = md5_hash(key_id.encode("utf-8")).hexdigest()
            self._pagination_tokens[next_continuation_token] = key_id
        else:
            is_truncated = "false"  # type: ignore
            next_continuation_token = None
        return result_keys, is_truncated, next_continuation_token

    @staticmethod
    def _get_name(key: Union[str, FakeKey]) -> str:
        if isinstance(key, FakeKey):
            return key.name
        else:
            return key

    def _set_delete_marker(self, bucket_name: str, key_name: str) -> FakeDeleteMarker:
        bucket = self.get_bucket(bucket_name)
        delete_marker = FakeDeleteMarker(key=bucket.keys[key_name])
        bucket.keys[key_name] = delete_marker
        return delete_marker

    def delete_object_tagging(
        self, bucket_name: str, key_name: str, version_id: Optional[str] = None
    ) -> None:
        key = self.get_object(bucket_name, key_name, version_id=version_id)
        bucket = self.get_bucket(bucket_name)

        self.tagger.delete_all_tags_for_resource(key.arn)  # type: ignore

        notifications.send_event(
            self.account_id,
            notifications.S3NotificationEvent.OBJECT_TAGGING_DELETE_EVENT,
            bucket,
            key,
        )

    def delete_object(
        self,
        bucket_name: str,
        key_name: str,
        version_id: Optional[str] = None,
        bypass: bool = False,
    ) -> Tuple[bool, Optional[Dict[str, Any]]]:
        bucket = self.get_bucket(bucket_name)

        response_meta = {}
        delete_key = bucket.keys.get(key_name)

        try:
            if not bucket.is_versioned:
                bucket.keys.pop(key_name)
                notifications.send_event(
                    self.account_id,
                    notifications.S3NotificationEvent.OBJECT_REMOVED_DELETE_EVENT,
                    bucket,
                    delete_key,
                )
            else:
                if version_id is None:
                    delete_marker = self._set_delete_marker(bucket_name, key_name)
                    response_meta["version-id"] = delete_marker.version_id
                    response_meta["delete-marker"] = "true"
                else:
                    if key_name not in bucket.keys:
                        raise KeyError

                    response_meta["version-id"] = version_id

                    for key in bucket.keys.getlist(key_name):
                        if str(key.version_id) == str(version_id):
                            if (
                                hasattr(key, "is_locked")
                                and key.is_locked
                                and not bypass
                            ):
                                raise AccessDeniedByLock

                            if type(key) is FakeDeleteMarker:
                                if type(key.key) is FakeDeleteMarker:  # type: ignore
                                    # Our key is a DeleteMarker, that usually contains a link to the actual FakeKey
                                    # But: If we have deleted the FakeKey multiple times,
                                    # We have a DeleteMarker linking to a DeleteMarker (etc..) linking to a FakeKey
                                    response_meta["delete-marker"] = "true"
                                # The alternative is that we're deleting the DeleteMarker that points directly to a FakeKey
                                # In this scenario, AWS does not return the `delete-marker` header

                            break

                    bucket.keys.setlist(
                        key_name,
                        [
                            key
                            for key in bucket.keys.getlist(key_name)
                            if str(key.version_id) != str(version_id)
                        ],
                    )

                    if not bucket.keys.getlist(key_name):
                        bucket.keys.pop(key_name)
                        notifications.send_event(
                            self.account_id,
                            notifications.S3NotificationEvent.OBJECT_REMOVED_DELETE_EVENT,
                            bucket,
                            delete_key,
                        )

            return True, response_meta
        except KeyError:
            return False, None

    def delete_objects(
        self, bucket_name: str, objects: List[Dict[str, Any]]
    ) -> List[Tuple[str, Optional[str]]]:
        deleted_objects = []
        for object_ in objects:
            key_name = object_["Key"]
            version_id = object_.get("VersionId", None)

            self.delete_object(bucket_name, key_name, version_id=version_id)
            deleted_objects.append((key_name, version_id))
        return deleted_objects

    def copy_object(
        self,
        src_key: FakeKey,
        dest_bucket_name: str,
        dest_key_name: str,
        storage: Optional[str] = None,
        encryption: Optional[str] = None,
        kms_key_id: Optional[str] = None,
        bucket_key_enabled: Any = None,
        mdirective: Optional[str] = None,
        metadata: Optional[Any] = None,
        website_redirect_location: Optional[str] = None,
        lock_mode: Optional[str] = None,
        lock_legal_status: Optional[str] = None,
        lock_until: Optional[str] = None,
    ) -> None:
        bucket = self.get_bucket(dest_bucket_name)
        if src_key.name == dest_key_name and src_key.bucket_name == dest_bucket_name:
            if src_key.encryption and src_key.encryption != "AES256" and not encryption:
                # this a special case, as now S3 default to AES256 when not provided
                # if the source key had encryption, and we did not specify it for the destination, S3 will accept a
                # copy in place even without any required attributes
                encryption = "AES256"

            if not any(
                (
                    storage,
                    encryption,
                    mdirective == "REPLACE",
                    website_redirect_location,
                    bucket.encryption,  # S3 will allow copy in place if the bucket has encryption configured
                    src_key._version_id and bucket.is_versioned,
                )
            ):
                raise CopyObjectMustChangeSomething

        new_key = self.put_object(
            bucket_name=dest_bucket_name,
            key_name=dest_key_name,
            value=src_key.value,
            storage=storage,
            multipart=src_key.multipart,
            encryption=encryption,
            kms_key_id=kms_key_id,  # TODO: use aws managed key if not provided
            bucket_key_enabled=bucket_key_enabled,
            lock_mode=lock_mode,
            lock_legal_status=lock_legal_status,
            lock_until=lock_until,
            disable_notification=True,  # avoid sending PutObject events here
        )
        self.tagger.copy_tags(src_key.arn, new_key.arn)
        if mdirective != "REPLACE":
            new_key.set_metadata(src_key.metadata)
        else:
            new_key.set_metadata(metadata)

        if website_redirect_location:
            new_key.website_redirect_location = website_redirect_location

        if src_key.storage_class in ARCHIVE_STORAGE_CLASSES:
            # Object copied from Glacier object should not have expiry
            new_key.set_expiry(None)

        if src_key.checksum_value:
            new_key.checksum_value = src_key.checksum_value
            new_key.checksum_algorithm = src_key.checksum_algorithm

        # Send notifications that an object was copied
        notifications.send_event(
            self.account_id,
            notifications.S3NotificationEvent.OBJECT_CREATED_COPY_EVENT,
            bucket,
            new_key,
        )

    def put_bucket_acl(self, bucket_name: str, acl: Optional[FakeAcl]) -> None:
        bucket = self.get_bucket(bucket_name)
        bucket.set_acl(acl)

    def get_bucket_acl(self, bucket_name: str) -> Optional[FakeAcl]:
        bucket = self.get_bucket(bucket_name)
        return bucket.acl

    def get_bucket_cors(self, bucket_name: str) -> List[CorsRule]:
        bucket = self.get_bucket(bucket_name)
        return bucket.cors

    def get_bucket_lifecycle(self, bucket_name: str) -> List[LifecycleRule]:
        bucket = self.get_bucket(bucket_name)
        return bucket.rules

    def get_bucket_location(self, bucket_name: str) -> str:
        bucket = self.get_bucket(bucket_name)

        return bucket.location

    def get_bucket_logging(self, bucket_name: str) -> Dict[str, Any]:
        bucket = self.get_bucket(bucket_name)
        return bucket.logging

    def get_bucket_notification_configuration(
        self, bucket_name: str
    ) -> Optional[NotificationConfiguration]:
        bucket = self.get_bucket(bucket_name)
        return bucket.notification_configuration

    def select_object_content(
        self,
        bucket_name: str,
        key_name: str,
        select_query: str,
        input_details: Dict[str, Any],
        output_details: Dict[str, Any],
    ) -> List[bytes]:
        """
        Highly experimental. Please raise an issue if you find any inconsistencies/bugs.

        Known missing features:
         - Function aliases (count(*) as cnt)
         - Most functions (only count() is supported)
         - Result is always in JSON
         - FieldDelimiters are ignored
        """
        self.get_bucket(bucket_name)
        key = self.get_object(bucket_name, key_name)
        if key is None:
            raise MissingKey(key=key_name)
        if input_details.get("CompressionType") == "GZIP":
            with gzip.open(BytesIO(key.value), "rt") as f:
                query_input = f.read()
        elif input_details.get("CompressionType") == "BZIP2":
            query_input = bz2.decompress(key.value).decode("utf-8")
        else:
            query_input = key.value.decode("utf-8")
        if "CSV" in input_details:
            # input is in CSV - we need to convert it to JSON before parsing
            from py_partiql_parser import csv_to_json

            use_headers = (input_details.get("CSV") or {}).get(
                "FileHeaderInfo", ""
            ) == "USE"
            query_input = csv_to_json(query_input, use_headers)
        query_result = parse_query(query_input, select_query)  # type: ignore

        record_delimiter = "\n"
        if "JSON" in output_details:
            record_delimiter = (output_details.get("JSON") or {}).get(
                "RecordDelimiter"
            ) or "\n"
        elif "CSV" in output_details:
            record_delimiter = (output_details.get("CSV") or {}).get(
                "RecordDelimiter"
            ) or "\n"

        if "CSV" in output_details:
            field_delim = (output_details.get("CSV") or {}).get("FieldDelimiter") or ","

            from py_partiql_parser import json_to_csv

            query_result = json_to_csv(query_result, field_delim, record_delimiter)
            return [query_result.encode("utf-8")]  # type: ignore

        else:
            from py_partiql_parser import SelectEncoder

            return [
                (
                    json.dumps(x, indent=None, separators=(",", ":"), cls=SelectEncoder)
                    + record_delimiter
                ).encode("utf-8")
                for x in query_result
            ]

    def restore_object(
        self, bucket_name: str, key_name: str, days: Optional[str], type_: Optional[str]
    ) -> bool:
        key = self.get_object(bucket_name, key_name)
        if not key:
            raise MissingKey

        if days is None and type_ is None:
            raise DaysMustProvidedExceptForSelectRequest()

        if days and type_:
            raise DaysMustNotProvidedForSelectRequest()

        if key.storage_class not in ARCHIVE_STORAGE_CLASSES:
            raise InvalidObjectState(storage_class=key.storage_class)
        had_expiry_date = key.expiry_date is not None
        if days:
            key.restore(int(days))
        return had_expiry_date

    def upload_file(self) -> None:
        # Listed for the implementation coverage
        # Implementation part of responses.py
        pass

    def upload_fileobj(self) -> None:
        # Listed for the implementation coverage
        # Implementation part of responses.py
        pass


class S3BackendDict(BackendDict[S3Backend]):
    """
    Encapsulation class to hold S3 backends.

    This is specialised to include additional attributes to help multi-account support in S3
    but is otherwise identical to the superclass.
    """

    def __init__(
        self,
        backend: Any,
        service_name: str,
        use_boto3_regions: bool = True,
        additional_regions: Optional[List[str]] = None,
    ):
        super().__init__(backend, service_name, use_boto3_regions, additional_regions)

        # Maps bucket names to (partition, account IDs). This is used to locate the exact S3Backend
        # holding the bucket and to maintain the common bucket namespace.
        self.bucket_accounts: Dict[str, Tuple[str, str]] = {}


s3_backends = S3BackendDict(
    S3Backend,
    service_name="s3",
    use_boto3_regions=False,
    additional_regions=PARTITION_NAMES,
)
