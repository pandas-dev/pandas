import io
import re
import urllib.parse
from collections.abc import Iterator
from typing import Any
from urllib.parse import parse_qs, unquote, urlencode, urlparse, urlunparse
from xml.dom import minidom
from xml.parsers.expat import ExpatError

import xmltodict

from moto import settings
from moto.core.common_types import TYPE_RESPONSE
from moto.core.mime_types import APP_XML
from moto.core.responses import ActionResult, BaseResponse, EmptyResult
from moto.core.utils import (
    ALT_DOMAIN_SUFFIXES,
    ensure_boolean,
    extract_region_from_aws_authorization,
    path_url,
    str_to_rfc_1123_datetime,
)
from moto.s3bucket_path.utils import (
    bucket_name_from_url as bucketpath_bucket_name_from_url,
)
from moto.s3bucket_path.utils import (
    parse_key_name as bucketpath_parse_key_name,
)
from moto.utilities.utils import PARTITION_NAMES, get_partition

from ..core.exceptions import SignatureDoesNotMatchError
from .exceptions import (
    AccessForbidden,
    BucketAccessDeniedError,
    BucketAlreadyExists,
    BucketAlreadyOwnedByYou,
    BucketNotEmpty,
    DuplicateTagKeys,
    HeadOnDeleteMarker,
    IllegalLocationConstraintException,
    IncompatibleLocationConstraintException,
    InvalidContinuationToken,
    InvalidLocationConstraintException,
    InvalidMaxPartArgument,
    InvalidMaxPartNumberArgument,
    InvalidNamespaceHeaderException,
    InvalidNotificationARN,
    InvalidObjectState,
    InvalidPartOrder,
    InvalidRange,
    LockNotEnabled,
    MalformedACLError,
    MalformedXML,
    MethodNotAllowed,
    MissingBucket,
    MissingKey,
    MissingRequestBody,
    MissingUploadObjectWithObjectLockHeaders,
    MissingVersion,
    NoSuchBucketPolicy,
    NoSuchCORSConfiguration,
    NoSuchLifecycleConfiguration,
    NoSuchTagSet,
    NoSuchWebsiteConfiguration,
    NoSystemTags,
    NotAnIntegerException,
    ObjectNotInActiveTierError,
    OwnershipControlsNotFoundError,
    PreconditionFailed,
    RangeNotSatisfiable,
    ReplicationConfigurationNotFoundError,
    S3AclAndGrantError,
    S3ClientError,
    ServerSideEncryptionConfigurationNotFoundError,
    VersioningNotEnabledForReplication,
)
from .models import (
    FakeAcl,
    FakeBucket,
    FakeDeleteMarker,
    FakeGrant,
    FakeGrantee,
    FakeKey,
    S3Backend,
    TransitionDefaultMinimumObjectSize,
    get_canned_acl,
    s3_backends,
)
from .utils import (
    ARCHIVE_STORAGE_CLASSES,
    bucket_name_from_url,
    compute_checksum,
    cors_matches_origin,
    metadata_from_headers,
    parse_region_from_url,
)

DEFAULT_REGION_NAME = "us-east-1"

ACTION_MAP = {
    "BUCKET": {
        "HEAD": {"DEFAULT": "HeadBucket"},
        "GET": {
            "uploads": "ListBucketMultipartUploads",
            "location": "GetBucketLocation",
            "lifecycle": "GetLifecycleConfiguration",
            "versioning": "GetBucketVersioning",
            "policy": "GetBucketPolicy",
            "website": "GetBucketWebsite",
            "acl": "GetBucketAcl",
            "tagging": "GetBucketTagging",
            "logging": "GetBucketLogging",
            "cors": "GetBucketCORS",
            "notification": "GetBucketNotification",
            "accelerate": "GetAccelerateConfiguration",
            "versions": "ListBucketVersions",
            "public_access_block": "GetPublicAccessBlock",
            "DEFAULT": "ListBucket",
        },
        "PUT": {
            "lifecycle": "PutLifecycleConfiguration",
            "versioning": "PutBucketVersioning",
            "policy": "PutBucketPolicy",
            "website": "PutBucketWebsite",
            "acl": "PutBucketAcl",
            "tagging": "PutBucketTagging",
            "logging": "PutBucketLogging",
            "cors": "PutBucketCORS",
            "notification": "PutBucketNotification",
            "accelerate": "PutAccelerateConfiguration",
            "public_access_block": "PutPublicAccessBlock",
            "DEFAULT": "CreateBucket",
        },
        "DELETE": {
            "lifecycle": "PutLifecycleConfiguration",
            "policy": "DeleteBucketPolicy",
            "website": "DeleteBucketWebsite",
            "tagging": "PutBucketTagging",
            "cors": "PutBucketCORS",
            "public_access_block": "DeletePublicAccessBlock",
            "DEFAULT": "DeleteBucket",
        },
    },
    "KEY": {
        "HEAD": {"DEFAULT": "GetObject"},
        "GET": {
            "uploadId": "ListMultipartUploadParts",
            "acl": "GetObjectAcl",
            "tagging": "GetObjectTagging",
            "versionId": "GetObjectVersion",
            "DEFAULT": "GetObject",
        },
        "PUT": {
            "acl": "PutObjectAcl",
            "tagging": "PutObjectTagging",
            "DEFAULT": "PutObject",
        },
        "DELETE": {
            "uploadId": "AbortMultipartUpload",
            "versionId": "DeleteObjectVersion",
            "DEFAULT": "DeleteObject",
        },
        "POST": {
            "uploads": "PutObject",
            "restore": "RestoreObject",
            "uploadId": "PutObject",
            "select": "SelectObject",
        },
    },
    "CONTROL": {
        "GET": {"publicAccessBlock": "GetPublicAccessBlock"},
        "PUT": {"publicAccessBlock": "PutPublicAccessBlock"},
        "DELETE": {"publicAccessBlock": "DeletePublicAccessBlock"},
    },
}

ALLOWED_HEADER_OVERRIDES = {
    "response-content-type": "content-type",
    "response-content-language": "Content-Language",
    "response-expires": "Expires",
    "response-cache-control": "Cache-Control",
    "response-content-disposition": "Content-Disposition",
    "response-content-encoding": "Content-Encoding",
}


def parse_key_name(pth: str) -> str:
    # strip the first '/' left by urlparse
    return pth[1:] if pth.startswith("/") else pth


class S3Response(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="s3")
        # Whatever format requests come in, we should never touch them
        # There are some nuances here - this decision should be method-specific, instead of service-specific
        # E.G.: we don't want to touch put_object(), but we might have to decompress put_object_configuration()
        # Taking the naive approach to never decompress anything from S3 for now
        self.allow_request_decompression = False

    def setup_class(self, request: Any, full_url: str, headers: Any) -> None:  # type: ignore[override]
        super().setup_class(request, full_url, headers, use_raw_body=True)
        self.region = parse_region_from_url(full_url, use_default_region=False)
        if self.region is None:
            self.region = (
                extract_region_from_aws_authorization(
                    request.headers.get("Authorization", "")
                )
                or DEFAULT_REGION_NAME
            )
        self.bucket_name = self.parse_bucket_name_from_url(request, full_url)
        self.request = request
        if (
            self.request.headers.get("x-amz-content-sha256")
            == "STREAMING-UNSIGNED-PAYLOAD-TRAILER"
        ):
            self.body = self._handle_encoded_body(self.body)

        if request.headers.get("x-amz-content-sha256") in [
            "STREAMING-AWS4-HMAC-SHA256-PAYLOAD",
            "STREAMING-AWS4-HMAC-SHA256-PAYLOAD-TRAILER",
        ]:
            self.body = self._handle_v4_chunk_signatures(
                self.raw_body, int(request.headers["x-amz-decoded-content-length"])
            )

    def get_safe_path(self) -> str:
        return unquote(self.raw_path)

    @property
    def is_access_point(self) -> bool:
        return ".s3-accesspoint." in self.headers["host"]

    @property
    def backend(self) -> S3Backend:
        return s3_backends[self.current_account][get_partition(self.region)]

    @property
    def should_autoescape(self) -> bool:
        return True

    def abort_multipart_upload(self) -> TYPE_RESPONSE:
        upload_id = self._get_param("uploadId")
        self.backend.abort_multipart_upload(self.bucket_name, upload_id)
        return 204, {}, ""

    def all_buckets(self) -> TYPE_RESPONSE:
        self.data["Action"] = "ListAllMyBuckets"
        self._authenticate_and_authorize_s3_action()

        # No bucket specified. Listing all buckets
        prefix = self._get_param("prefix")
        bucket_region = self._get_param("bucket-region")
        max_buckets = self._get_param("max-buckets")
        all_buckets = self.backend.list_buckets(
            prefix=prefix, bucket_region=bucket_region, max_buckets=max_buckets
        )
        buckets = [
            {"Name": bucket.name, "CreationDate": bucket.creation_date_ISO8601}
            for bucket in all_buckets
        ]
        self.data["Action"] = "ListBuckets"
        return self.serialized(
            ActionResult(
                {
                    "Owner": {
                        "ID": "bcaf1ffd86f41161ca5fb16fd081034f",
                        "DisplayName": "webfile",
                    },
                    "Buckets": buckets,
                }
            )
        )

    def subdomain_based_buckets(self, request: Any) -> bool:
        if settings.S3_IGNORE_SUBDOMAIN_BUCKETNAME:
            return False
        host = request.headers.get("host", request.headers.get("Host"))
        if not host:
            host = urlparse(request.url).netloc

        custom_endpoints = settings.get_s3_custom_endpoints()
        if (
            host
            and custom_endpoints
            and any(host in endpoint for endpoint in custom_endpoints)
        ):
            # Default to path-based buckets for S3-compatible SDKs (Ceph, DigitalOcean Spaces, etc)
            return False

        if (
            not host
            or host.startswith("localhost")
            or host.startswith("localstack")
            or host.startswith("host.docker.internal")
            or re.match(r"^[^.]+$", host)
            or re.match(r"^.*\.svc\.cluster\.local:?\d*$", host)
        ):
            # Default to path-based buckets for (1) localhost, (2) localstack hosts (e.g. localstack.dev),
            # (3) local host names that do not contain a "." (e.g., Docker container host names), or
            # (4) kubernetes host names
            return False

        match = re.match(r"^([^\[\]:]+)(:\d+)?$", host)
        if match:
            match = re.match(
                r"((25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)(\.|$)){4}", match.groups()[0]
            )
            if match:
                return False

        match = re.match(r"^\[(.+)\](:\d+)?$", host)
        if match:
            match = re.match(
                r"^(((?=.*(::))(?!.*\3.+\3))\3?|[\dA-F]{1,4}:)([\dA-F]{1,4}(\3|:\b)|\2){5}(([\dA-F]{1,4}(\3|:\b|$)|\2){2}|(((2[0-4]|1\d|[1-9])?\d|25[0-5])\.?\b){4})\Z",
                match.groups()[0],
                re.IGNORECASE,
            )
            if match:
                return False

        path_based = (
            host == "s3.amazonaws.com"
            or re.match(
                r"s3[\.\-](?:(?:dualstack|fips)[\.\-])?([^.]*)\.amazonaws\.com", host
            )
            or any(
                re.match(
                    r"s3[\.\-](?:(?:dualstack|fips)[\.\-])?([^.]*)\." + suffix, host
                )
                for suffix in ALT_DOMAIN_SUFFIXES
            )
        )
        return not path_based

    def is_delete_keys(self) -> bool:
        qs = parse_qs(urlparse(self.path).query, keep_blank_values=True)
        return "delete" in qs

    def parse_bucket_name_from_url(self, request: Any, url: str) -> str:
        bucket_name = ""
        if self.subdomain_based_buckets(request):
            bucket_name = bucket_name_from_url(url)  # type: ignore
        else:
            bucket_name = bucketpath_bucket_name_from_url(url)  # type: ignore

        if self.is_access_point:
            # import here to avoid circular dependency error
            from moto.s3control import s3control_backends

            ap_name = bucket_name[: -(len(self.current_account) + 1)]
            ap = s3control_backends[self.current_account][
                self.partition
            ].get_access_point(self.current_account, ap_name)
            bucket_name = ap.bucket

        return bucket_name

    def parse_key_name(self) -> str:
        url = self.get_safe_path()
        if self.subdomain_based_buckets(self.request):
            return parse_key_name(url)
        else:
            return bucketpath_parse_key_name(url)

    @classmethod
    def ambiguous_dispatch(cls, *args: Any, **kwargs: Any) -> TYPE_RESPONSE:
        return cls().ambiguous_response(*args, **kwargs)

    @classmethod
    def bucket_dispatch(cls, *args: Any, **kwargs: Any) -> TYPE_RESPONSE:
        return cls().bucket_response(*args, **kwargs)

    @classmethod
    def key_dispatch(cls, *args: Any, **kwargs: Any) -> TYPE_RESPONSE:
        return cls().key_response(*args, **kwargs)

    def ambiguous_response(
        self, request: Any, full_url: str, headers: Any
    ) -> TYPE_RESPONSE:
        # Depending on which calling format the client is using, we don't know
        # if this is a bucket or key request so we have to check
        if self.subdomain_based_buckets(request):
            return self.key_response(request, full_url, headers)
        else:
            # Using path-based buckets
            return self.bucket_response(request, full_url, headers)

    def bucket_response(
        self, request: Any, full_url: str, headers: Any
    ) -> TYPE_RESPONSE:
        self.setup_class(request, full_url, headers)
        bucket_name = self.parse_bucket_name_from_url(request, full_url)
        self.backend.log_incoming_request(request, bucket_name)
        try:
            response = self._bucket_response(request, full_url)
        except S3ClientError as s3error:
            response = self.serialized(ActionResult(s3error))

        return self._send_response(response)

    @classmethod
    def _send_response(cls, response: TYPE_RESPONSE | str | bytes) -> TYPE_RESPONSE:  # type: ignore
        if isinstance(response, (str, bytes)):
            status_code = 200
            headers: dict[str, str] = {}
            body = response
        else:
            status_code, headers, body = response

        headers, body = cls._enrich_response(headers, body)

        if isinstance(body, str):
            body = body.encode("utf-8")

        if body and "content-type" not in (header.lower() for header in headers.keys()):
            headers["content-type"] = APP_XML

        return status_code, headers, body

    def _bucket_response(self, request: Any, full_url: str) -> str | TYPE_RESPONSE:
        querystring = self._get_querystring(request, full_url)
        method = request.method
        bucket_name = self.parse_bucket_name_from_url(request, full_url)
        if not bucket_name:
            # If no bucket specified, list all buckets
            return self.all_buckets()

        self.data["BucketName"] = bucket_name

        if method == "HEAD":
            self._set_action("BUCKET", "HEAD", querystring)
            self._authenticate_and_authorize_s3_action(bucket_name=bucket_name)
            return self.head_bucket()
        elif method == "GET":
            return self._bucket_response_get(bucket_name, querystring)
        elif method == "PUT":
            return self._bucket_response_put(request, bucket_name, querystring)
        elif method == "DELETE":
            return self._bucket_response_delete(bucket_name, querystring)
        elif method == "POST":
            return self._bucket_response_post(request, bucket_name)
        elif method == "OPTIONS":
            return self._response_options(request.headers, bucket_name)
        else:
            raise NotImplementedError(
                f"Method {method} has not been implemented in the S3 backend yet"
            )

    def _get_querystring(self, request: Any, full_url: str) -> dict[str, Any]:
        # Flask's Request has the querystring already parsed
        # In ServerMode, we can use this, instead of manually parsing this
        if hasattr(request, "args"):
            query_dict = {}
            for key, val in dict(request.args).items():
                # The parse_qs-method returns List[str, List[Any]]
                # Ensure that we confirm to the same response-type here
                query_dict[key] = val if isinstance(val, list) else [val]
            return query_dict

        parsed_url = urlparse(full_url)
        # full_url can be one of two formats, depending on the version of werkzeug used:
        # http://foobaz.localhost:5000/?prefix=bar%2Bbaz
        # http://foobaz.localhost:5000/?prefix=bar+baz
        # Werkzeug helpfully encodes the plus-sign for us, from >= 2.1.0
        # However, the `parse_qs` method will (correctly) replace '+' with a space
        #
        # Workaround - manually reverse the encoding.
        # Keep the + encoded, ensuring that parse_qsl doesn't replace it, and parse_qsl will unquote it afterwards
        qs = (parsed_url.query or "").replace("+", "%2B")
        return parse_qs(qs, keep_blank_values=True)

    def head_bucket(self) -> TYPE_RESPONSE:
        try:
            bucket = self.backend.head_bucket(self.bucket_name)
        except MissingBucket:
            # Unless we do this, boto3 does not raise ClientError on
            # HEAD (which the real API responds with), and instead
            # raises NoSuchBucket, leading to inconsistency in
            # error response between real and mocked responses.
            return 404, {}, ""
        headers = {"x-amz-bucket-region": bucket.region_name, "content-type": APP_XML}
        return 200, headers, ""

    def _set_cors_headers_options(
        self, headers: dict[str, str], bucket: FakeBucket
    ) -> None:
        """
        TODO: smarter way of matching the right CORS rule:
        See https://docs.aws.amazon.com/AmazonS3/latest/userguide/cors.html

        "When Amazon S3 receives a preflight request from a browser, it evaluates
        the CORS configuration for the bucket and uses the first CORSRule rule
        that matches the incoming browser request to enable a cross-origin request."
        This here just uses all rules and the last rule will override the previous ones
        if they are re-defining the same headers.
        """

        def _to_string(header: list[str] | str) -> str:
            # We allow list and strs in header values. Transform lists in comma-separated strings
            if isinstance(header, list):
                return ", ".join(header)
            return header

        for cors_rule in bucket.cors:
            if cors_rule.allowed_methods is not None:
                self.response_headers["Access-Control-Allow-Methods"] = _to_string(
                    cors_rule.allowed_methods
                )
            if cors_rule.allowed_origins is not None:
                origin = headers.get("Origin")
                if cors_matches_origin(origin, cors_rule.allowed_origins):  # type: ignore
                    self.response_headers["Access-Control-Allow-Origin"] = origin  # type: ignore
                else:
                    raise AccessForbidden(
                        "CORSResponse: This CORS request is not allowed. This is usually because the evalution of Origin, request method / Access-Control-Request-Method or Access-Control-Request-Headers are not whitelisted by the resource's CORS spec."
                    )
            if cors_rule.allowed_headers is not None:
                self.response_headers["Access-Control-Allow-Headers"] = _to_string(
                    cors_rule.allowed_headers
                )
            if cors_rule.exposed_headers is not None:
                self.response_headers["Access-Control-Expose-Headers"] = _to_string(
                    cors_rule.exposed_headers
                )
            if cors_rule.max_age_seconds is not None:
                self.response_headers["Access-Control-Max-Age"] = _to_string(
                    cors_rule.max_age_seconds
                )

    def _response_options(
        self, headers: dict[str, str], bucket_name: str
    ) -> TYPE_RESPONSE:
        # Return 200 with the headers from the bucket CORS configuration
        self._authenticate_and_authorize_s3_action(bucket_name=bucket_name)
        try:
            bucket = self.backend.head_bucket(bucket_name)
        except MissingBucket:
            # AWS S3 seems to return 403 on OPTIONS and 404 on GET/HEAD
            return 403, {}, ""

        self._set_cors_headers_options(headers, bucket)

        return 200, self.response_headers, ""

    @staticmethod
    def _acl_to_dict(acl: Any) -> dict[str, Any]:
        """Convert an ACL object to a dict matching the botocore shape.

        Mirrors the logic and hardcoded values from the original
        S3_OBJECT_ACL_RESPONSE Jinja template.
        """
        grants = [
            {
                "Grantee": {
                    "Type": grant.grantees[0].type,
                    "URI": grant.grantees[0].uri or None,
                    "ID": grant.grantees[0].id or None,
                    "DisplayName": grant.grantees[0].display_name or None,
                },
                "Permission": grant.permissions[0],
            }
            for grant in acl.grants
        ]
        return {
            "Owner": {
                "ID": "75aa57f09aa0c8caeab4f8c24e99d10f8e7faeebf76c078efc7c6caea54ba06a",
                "DisplayName": "webfile",
            },
            "Grants": grants,
        }

    def _get_cors_headers_other(self) -> dict[str, Any]:
        """
        Returns a dictionary with the appropriate CORS headers
        Should be used for non-OPTIONS requests only
        Applicable if the 'Origin' header matches one of a CORS-rules - returns an empty dictionary otherwise
        """
        response_headers: dict[str, Any] = {}
        try:
            origin = self.headers.get("Origin")
            if not origin:
                return response_headers
            bucket = self.backend.get_bucket(self.bucket_name)

            def _to_string(header: list[str] | str) -> str:
                # We allow list and strs in header values. Transform lists in comma-separated strings
                if isinstance(header, list):
                    return ", ".join(header)
                return header

            for cors_rule in bucket.cors:
                if cors_rule.allowed_origins is not None:
                    if cors_matches_origin(origin, cors_rule.allowed_origins):
                        response_headers["Access-Control-Allow-Origin"] = origin
                        if cors_rule.allowed_methods is not None:
                            response_headers["Access-Control-Allow-Methods"] = (
                                _to_string(cors_rule.allowed_methods)
                            )
                        if cors_rule.allowed_headers is not None:
                            response_headers["Access-Control-Allow-Headers"] = (
                                _to_string(cors_rule.allowed_headers)
                            )
                        if cors_rule.exposed_headers is not None:
                            response_headers["Access-Control-Expose-Headers"] = (
                                _to_string(cors_rule.exposed_headers)
                            )
                        if cors_rule.max_age_seconds is not None:
                            response_headers["Access-Control-Max-Age"] = _to_string(
                                cors_rule.max_age_seconds
                            )

                        return response_headers
        except S3ClientError:
            pass
        return response_headers

    def _bucket_response_get(
        self, bucket_name: str, querystring: dict[str, Any]
    ) -> str | TYPE_RESPONSE:
        self._set_action("BUCKET", "GET", querystring)
        self._authenticate_and_authorize_s3_action(bucket_name=bucket_name)

        if "object-lock" in querystring:
            return self.get_object_lock_configuration()

        if "uploads" in querystring:
            return self.list_multipart_uploads()
        elif "location" in querystring:
            return self.get_bucket_location()
        elif "lifecycle" in querystring:
            return self.get_bucket_lifecycle()
        elif "versioning" in querystring:
            return self.get_bucket_versioning()
        elif "policy" in querystring:
            return self.get_bucket_policy()
        elif "website" in querystring:
            return self.get_bucket_website()
        elif "acl" in querystring:
            return self.get_bucket_acl()
        elif "tagging" in querystring:
            return self.get_bucket_tags()
        elif "logging" in querystring:
            return self.get_bucket_logging()
        elif "cors" in querystring:
            return self.get_bucket_cors()
        elif "notification" in querystring:
            return self.get_bucket_notification()
        elif "accelerate" in querystring:
            return self.get_bucket_accelerate_configuration()
        elif "publicAccessBlock" in querystring:
            return self.get_public_access_block()
        elif "inventory" in querystring:
            # Only GET includes "id" in the querystring, LIST does not
            if "id" in querystring:
                return self.get_bucket_inventory_configuration()
            else:
                return self.list_bucket_inventory_configurations()
        elif "versions" in querystring:
            return self.list_object_versions()
        elif "encryption" in querystring:
            return self.get_bucket_encryption()
        elif querystring.get("list-type", [None])[0] == "2":
            return self.list_objects_v2()
        elif "replication" in querystring:
            return self.get_bucket_replication()
        elif "ownershipControls" in querystring:
            return self.get_bucket_ownership_controls()

        return self.list_objects()

    def _set_action(
        self, action_resource_type: str, method: str, querystring: dict[str, Any]
    ) -> None:
        action_set = False
        for action_in_querystring, action in ACTION_MAP[action_resource_type][
            method
        ].items():
            if action_in_querystring in querystring:
                self.data["Action"] = action
                action_set = True
        if not action_set:
            self.data["Action"] = ACTION_MAP[action_resource_type][method]["DEFAULT"]

    def list_multipart_uploads(self) -> TYPE_RESPONSE:
        multiparts = list(
            self.backend.list_multipart_uploads(self.bucket_name).values()
        )
        if "prefix" in self.querystring:
            prefix = self.querystring.get("prefix", [None])[0]
            multiparts = [
                upload for upload in multiparts if upload.key_name.startswith(prefix)
            ]
        uploads = [
            {
                "Key": upload.key_name,
                "UploadId": upload.id,
                "Initiator": {
                    "ID": f"arn:aws:iam::{self.current_account}:user/user1-11111a31-17b5-4fb7-9df5-b111111f13de",
                    "DisplayName": "user1-11111a31-17b5-4fb7-9df5-b111111f13de",
                },
                "Owner": {
                    "ID": "75aa57f09aa0c8caeab4f8c24e99d10f8e7faeebf76c078efc7c6caea54ba06a",
                    "DisplayName": "webfile",
                },
                "StorageClass": "STANDARD",
                "Initiated": "2010-11-10T20:48:33.000Z",
                "ChecksumAlgorithm": upload.metadata.get("x-amz-checksum-algorithm"),
                "ChecksumType": upload.metadata.get("x-amz-checksum-type"),
            }
            for upload in multiparts
        ]
        result = {
            "Bucket": self.bucket_name,
            "KeyMarker": "",
            "UploadIdMarker": "",
            "MaxUploads": 1000,
            "IsTruncated": False,
            "Uploads": uploads,
        }
        self.data["Action"] = "ListMultipartUploads"
        return self.serialized(ActionResult(result))

    def list_objects(self) -> TYPE_RESPONSE:
        bucket = self.backend.get_bucket(self.bucket_name)
        querystring = self._get_querystring(self.request, self.uri)
        prefix = querystring.get("prefix", [None])[0]
        if prefix and isinstance(prefix, bytes):
            prefix = prefix.decode("utf-8")
        delimiter = self.querystring.get("delimiter", [None])[0]
        max_keys = int(
            self.querystring.get("max-keys", [settings.get_s3_default_max_keys()])[0]
        )
        marker = self.querystring.get("marker", [None])[0]
        encoding_type = self.querystring.get("encoding-type", [None])[0]

        (
            result_keys,
            result_folders,
            is_truncated,
            next_marker,
        ) = self.backend.list_objects(
            bucket=bucket,
            prefix=prefix,
            delimiter=delimiter,
            marker=marker,
            max_keys=max_keys,
        )

        contents = [
            {
                "Key": key.safe_name(encoding_type),
                "LastModified": key.last_modified_ISO8601,
                "ETag": key.etag,
                "Size": key.size,
                "StorageClass": key.storage_class,
                "Owner": {
                    "ID": "75aa57f09aa0c8caeab4f8c24e99d10f8e7faeebf76c078efc7c6caea54ba06a",
                    "DisplayName": "webfile",
                },
            }
            for key in result_keys
        ]

        result: dict[str, Any] = {
            "Name": bucket.name,
            "MaxKeys": max_keys,
            "IsTruncated": is_truncated,
            "Contents": contents,
        }
        if prefix is not None:
            result["Prefix"] = prefix
        if delimiter:
            result["Delimiter"] = delimiter
            result["CommonPrefixes"] = [{"Prefix": folder} for folder in result_folders]
        if encoding_type:
            result["EncodingType"] = encoding_type
        if next_marker:
            result["NextMarker"] = next_marker

        self.data["Action"] = "ListObjects"
        return self.serialized(ActionResult(result))

    def list_objects_v2(self) -> TYPE_RESPONSE:
        bucket = self.backend.get_bucket(self.bucket_name)

        continuation_token = self.querystring.get("continuation-token", [None])[0]
        if continuation_token is not None and continuation_token == "":
            raise InvalidContinuationToken()

        querystring = self._get_querystring(self.request, self.uri)
        prefix = querystring.get("prefix", [None])[0]
        if prefix and isinstance(prefix, bytes):
            prefix = prefix.decode("utf-8")
        delimiter = self.querystring.get("delimiter", [None])[0]

        fetch_owner = self.querystring.get("fetch-owner", [False])[0]
        max_keys = int(
            self.querystring.get("max-keys", [settings.get_s3_default_max_keys()])[0]
        )
        start_after = self.querystring.get("start-after", [None])[0]
        encoding_type = self.querystring.get("encoding-type", [None])[0]

        (
            truncated_keys,
            is_truncated,
            next_continuation_token,
        ) = self.backend.list_objects_v2(
            bucket=bucket,
            prefix=prefix,
            delimiter=delimiter,
            continuation_token=continuation_token,
            start_after=start_after,
            max_keys=max_keys,
        )

        result_keys, result_folders = self._split_truncated_keys(truncated_keys)

        key_count = len(result_keys) + len(result_folders)

        if encoding_type == "url":
            prefix = urllib.parse.quote(prefix) if prefix else ""
            result_folders = [urllib.parse.quote(folder) for folder in result_folders]

        contents = []
        for key in result_keys:
            item: dict[str, Any] = {
                "Key": key.safe_name(encoding_type),
                "LastModified": key.last_modified_ISO8601,
                "ETag": key.etag,
                "Size": key.size,
                "StorageClass": key.storage_class,
            }
            if fetch_owner:
                item["Owner"] = {
                    "ID": "75aa57f09aa0c8caeab4f8c24e99d10f8e7faeebf76c078efc7c6caea54ba06a",
                    "DisplayName": "webfile",
                }
            if key.checksum_algorithm:
                item["ChecksumAlgorithm"] = [key.checksum_algorithm]
            contents.append(item)

        result: dict[str, Any] = {
            "Name": bucket.name,
            "Prefix": prefix or "",
            "MaxKeys": max_keys,
            "KeyCount": key_count,
            "IsTruncated": is_truncated,
            "Contents": contents,
        }
        if delimiter:
            result["Delimiter"] = delimiter
            result["CommonPrefixes"] = [{"Prefix": folder} for folder in result_folders]
        if encoding_type:
            result["EncodingType"] = encoding_type
        if next_continuation_token:
            result["NextContinuationToken"] = next_continuation_token
        if not continuation_token and start_after:
            result["StartAfter"] = start_after

        self.data["Action"] = "ListObjectsV2"
        return self.serialized(ActionResult(result))

    def list_object_versions(self) -> TYPE_RESPONSE:
        delimiter = self.querystring.get("delimiter", [None])[0]
        key_marker = self.querystring.get("key-marker", [None])[0]
        max_keys = int(
            self.querystring.get("max-keys", [settings.get_s3_default_max_keys()])[0]
        )
        querystring = self._get_querystring(self.request, self.uri)
        prefix = querystring.get("prefix", [""])[0]
        version_id_marker = self.querystring.get("version-id-marker", [None])[0]

        bucket = self.backend.get_bucket(self.bucket_name)
        (
            versions,
            common_prefixes,
            delete_markers,
            next_key_marker,
            next_version_id_marker,
        ) = self.backend.list_object_versions(
            self.bucket_name,
            delimiter=delimiter,
            key_marker=key_marker,
            max_keys=max_keys,
            prefix=prefix,
            version_id_marker=version_id_marker,
        )
        key_list = versions

        is_truncated = next_key_marker is not None

        owner = {
            "ID": "75aa57f09aa0c8caeab4f8c24e99d10f8e7faeebf76c078efc7c6caea54ba06a",
            "DisplayName": "webfile",
        }
        version_list = [
            {
                "Key": key.name,
                "VersionId": key.version_id if key.version_id is not None else "null",
                "IsLatest": getattr(key, "is_latest", False),
                "LastModified": key.last_modified_ISO8601,
                "ETag": key.etag,
                "Size": key.size,
                "StorageClass": key.storage_class,
                "Owner": owner,
            }
            for key in key_list
        ]
        delete_marker_list = [
            {
                "Key": marker.name,
                "VersionId": marker.version_id,
                "IsLatest": getattr(marker, "is_latest", False),
                "LastModified": marker.last_modified_ISO8601,
                "Owner": owner,
            }
            for marker in delete_markers
        ]

        result: dict[str, Any] = {
            "Name": bucket.name,
            "Prefix": prefix if prefix is not None else None,
            "Delimiter": delimiter,
            "KeyMarker": key_marker or "",
            "VersionIdMarker": version_id_marker or "",
            "MaxKeys": max_keys,
            "IsTruncated": is_truncated,
            "Versions": version_list,
            "DeleteMarkers": delete_marker_list,
        }
        if common_prefixes:
            result["CommonPrefixes"] = [{"Prefix": p} for p in common_prefixes]
        if is_truncated:
            result["NextKeyMarker"] = next_key_marker
            result["NextVersionIdMarker"] = next_version_id_marker

        self.data["Action"] = "ListObjectVersions"
        return self.serialized(ActionResult(result))

    @staticmethod
    def _split_truncated_keys(truncated_keys: Any) -> Any:  # type: ignore[misc]
        result_keys = []
        result_folders = []
        for key in truncated_keys:
            if isinstance(key, FakeKey):
                result_keys.append(key)
            else:
                result_folders.append(key)
        return result_keys, result_folders

    def _get_location_constraint(self) -> str | None:
        try:
            if self.body:
                return xmltodict.parse(self.body)["CreateBucketConfiguration"][
                    "LocationConstraint"
                ]
        except KeyError:
            pass
        return None

    def _parse_pab_config(self) -> dict[str, Any]:
        parsed_xml = xmltodict.parse(self.body)
        parsed_xml["PublicAccessBlockConfiguration"].pop("@xmlns", None)

        return parsed_xml

    def _bucket_response_put(
        self,
        request: Any,
        bucket_name: str,
        querystring: dict[str, Any],
    ) -> str | TYPE_RESPONSE:
        if querystring and not request.headers.get("Content-Length"):
            return 411, {}, "Content-Length required"

        self._set_action("BUCKET", "PUT", querystring)
        self._authenticate_and_authorize_s3_action(bucket_name=bucket_name)

        if "object-lock" in querystring:
            config = self._process_lock_config_from_body()

            self.backend.put_object_lock_configuration(
                bucket_name,
                config.get("enabled"),  # type: ignore
                config.get("mode"),
                config.get("days"),
                config.get("years"),
            )
            return 200, {}, ""

        if "versioning" in querystring:
            body = self.body.decode("utf-8")
            ver = re.search(r"<Status>([A-Za-z]+)</Status>", body)
            if ver:
                self.backend.put_bucket_versioning(bucket_name, ver.group(1))
                self.data["Action"] = "PutBucketVersioning"
                return self.serialized(EmptyResult())
            else:
                return 404, {}, ""
        elif "lifecycle" in querystring:
            lifecycle_config = xmltodict.parse(self.body)["LifecycleConfiguration"]
            # Handle empty lifecycle configuration (no rules)
            if "Rule" not in lifecycle_config:
                rules = []
            else:
                rules = lifecycle_config["Rule"]
                if not isinstance(rules, list):
                    # If there is only one rule, xmldict returns just the item
                    rules = [rules]
            self.backend.put_bucket_lifecycle(bucket_name, rules)
            return ""
        elif "policy" in querystring:
            self.backend.put_bucket_policy(bucket_name, self.body)
            return "True"
        elif "acl" in querystring:
            # Headers are first. If not set, then look at the body (consistent with the documentation):
            acls = self._acl_from_headers(request.headers)
            if not acls:
                acls = self._acl_from_body()
            self.backend.put_bucket_acl(bucket_name, acls)
            return ""
        elif "tagging" in querystring:
            tagging = self._bucket_tagging_from_body()
            self.backend.put_bucket_tagging(bucket_name, tagging)
            return 204, {}, ""
        elif "website" in querystring:
            self.backend.put_bucket_website(bucket_name, self.body)
            return ""
        elif "cors" in querystring:
            try:
                self.backend.put_bucket_cors(bucket_name, self._cors_from_body())
                return ""
            except KeyError:
                raise MalformedXML()
        elif "logging" in querystring:
            try:
                self.backend.put_bucket_logging(bucket_name, self._logging_from_body())
                return ""
            except KeyError:
                raise MalformedXML()
        elif "notification" in querystring:
            try:
                self.backend.put_bucket_notification_configuration(
                    bucket_name, self._notification_config_from_body()
                )
                return ""
            except KeyError:
                raise MalformedXML()
        elif "accelerate" in querystring:
            try:
                accelerate_status = self._accelerate_config_from_body()
                self.backend.put_bucket_accelerate_configuration(
                    bucket_name, accelerate_status
                )
                return ""
            except KeyError:
                raise MalformedXML()

        elif "publicAccessBlock" in querystring:
            pab_config = self._parse_pab_config()
            self.backend.put_public_access_block(
                bucket_name, pab_config["PublicAccessBlockConfiguration"]
            )
            return ""
        elif "encryption" in querystring:
            try:
                self.backend.put_bucket_encryption(
                    bucket_name, self._encryption_config_from_body()
                )
                return ""
            except KeyError:
                raise MalformedXML()
        elif "replication" in querystring:
            bucket = self.backend.get_bucket(bucket_name)
            if not bucket.is_versioned:
                raise VersioningNotEnabledForReplication(bucket_name=bucket_name)
            replication_config = self._replication_config_from_xml(self.body)
            self.backend.put_bucket_replication(bucket_name, replication_config)
            return ""
        elif "ownershipControls" in querystring:
            ownership_rule = self._ownership_rule_from_body()
            self.backend.put_bucket_ownership_controls(
                bucket_name, ownership=ownership_rule
            )
            return ""
        elif "inventory" in querystring:
            inventory_config = self._inventory_config_from_body()
            self.backend.put_bucket_inventory_configuration(
                bucket_name, inventory_config
            )
            return ""
        else:
            # us-east-1, the default AWS region behaves a bit differently
            # - you should not use any location constraint
            location_constraint = self._get_location_constraint()
            if self.region == DEFAULT_REGION_NAME:
                # REGION = us-east-1 - we should never receive a LocationConstraint for `us-east-1`
                # However, you can create buckets in other regions from us-east-1
                if location_constraint == DEFAULT_REGION_NAME:
                    raise InvalidLocationConstraintException

            else:
                # Non-Standard region - the LocationConstraint must be equal to the actual region the request was send to
                if not location_constraint:
                    raise IllegalLocationConstraintException()
                if location_constraint != self.region:
                    raise IncompatibleLocationConstraintException(location_constraint)

            bucket_region = (
                location_constraint if location_constraint else DEFAULT_REGION_NAME
            )
            bucket_namespace = request.headers.get("x-amz-bucket-namespace")
            if bucket_namespace == "account-regional":
                expected_suffix = f"-{self.get_current_account()}-{bucket_region}-an"
                if not bucket_name.endswith(expected_suffix):
                    raise InvalidNamespaceHeaderException(
                        self.get_current_account(), bucket_region
                    )
            try:
                new_bucket = self.backend.create_bucket(
                    bucket_name, bucket_region, bucket_namespace=bucket_namespace
                )
            except BucketAlreadyExists:
                new_bucket = self.backend.get_bucket(bucket_name)
                if new_bucket.account_id == self.get_current_account():
                    # special cases when the bucket belongs to self
                    if (
                        new_bucket.region_name == DEFAULT_REGION_NAME
                        and self.region == DEFAULT_REGION_NAME
                    ):
                        # us-east-1 has different behavior - creating a bucket there is an idempotent operation
                        pass
                    else:
                        raise BucketAlreadyOwnedByYou(bucket_name=bucket_name)
                else:
                    raise

            if "x-amz-acl" in request.headers:
                # TODO: Support the XML-based ACL format
                self.backend.put_bucket_acl(
                    bucket_name, self._acl_from_headers(request.headers)
                )

            if (
                request.headers.get("x-amz-bucket-object-lock-enabled", "").lower()
                == "true"
            ):
                new_bucket.object_lock_enabled = True
                new_bucket.versioning_status = "Enabled"

            ownership_rule = request.headers.get("x-amz-object-ownership")
            if ownership_rule:
                new_bucket.ownership_rule = ownership_rule

            result: dict[str, Any] = {"Location": f"/{new_bucket.name}"}
            if new_bucket.bucket_namespace == "account-regional":
                result["BucketArn"] = new_bucket.arn
            self.data["Action"] = "CreateBucket"
            return self.serialized(ActionResult(result))

    def get_bucket_accelerate_configuration(self) -> TYPE_RESPONSE:
        accelerate_configuration = self.backend.get_bucket_accelerate_configuration(
            self.bucket_name
        )
        result = {"Status": accelerate_configuration}
        self.data["Action"] = "GetBucketAccelerateConfiguration"
        return self.serialized(ActionResult(result))

    def get_bucket_acl(self) -> TYPE_RESPONSE:
        acl = self.backend.get_bucket_acl(self.bucket_name)
        self.data["Action"] = "GetBucketAcl"
        return self.serialized(ActionResult(self._acl_to_dict(acl)))

    def get_bucket_cors(self) -> str | TYPE_RESPONSE:
        self.data["Action"] = "GetBucketCors"
        cors = self.backend.get_bucket_cors(self.bucket_name)
        if len(cors) == 0:
            raise NoSuchCORSConfiguration(bucket_name=self.bucket_name)
        cors_rules = [
            {
                "ID": rule.id_,
                "AllowedOrigins": rule.allowed_origins,
                "AllowedMethods": rule.allowed_methods,
                "AllowedHeaders": rule.allowed_headers,
                "ExposeHeaders": rule.exposed_headers,
                "MaxAgeSeconds": rule.max_age_seconds,
            }
            for rule in cors
        ]
        return self.serialized(ActionResult({"CORSRules": cors_rules}))

    def get_bucket_encryption(self) -> str | TYPE_RESPONSE:
        self.data["Action"] = "GetBucketEncryption"
        encryption = self.backend.get_bucket_encryption(self.bucket_name)
        if not encryption:
            raise ServerSideEncryptionConfigurationNotFoundError(
                bucket_name=self.bucket_name
            )
        rule = encryption["Rule"]
        result_rule: dict[str, Any] = {
            "ApplyServerSideEncryptionByDefault": rule[
                "ApplyServerSideEncryptionByDefault"
            ],
            "BucketKeyEnabled": ensure_boolean(rule.get("BucketKeyEnabled")),
        }
        self.data["Action"] = "GetBucketEncryption"
        return self.serialized(
            ActionResult(
                {
                    "ServerSideEncryptionConfiguration": {"Rules": [result_rule]},
                }
            )
        )

    def get_bucket_lifecycle(self) -> TYPE_RESPONSE:
        self.data["Action"] = "GetBucketLifecycleConfiguration"
        rules = self.backend.get_bucket_lifecycle(self.bucket_name)
        if not rules:
            raise NoSuchLifecycleConfiguration(bucket_name=self.bucket_name)

        lifecycle_rules = []
        for rule in rules:
            r: dict[str, Any] = {"ID": rule.id, "Status": rule.status}

            if rule.filter:
                f: dict[str, Any] = {}
                if rule.filter.prefix is not None:
                    f["Prefix"] = rule.filter.prefix
                if rule.filter.tag_key:
                    f["Tag"] = {
                        "Key": rule.filter.tag_key,
                        "Value": rule.filter.tag_value,
                    }
                if rule.filter.and_filter:
                    and_f: dict[str, Any] = {}
                    if rule.filter.and_filter.prefix is not None:
                        and_f["Prefix"] = rule.filter.and_filter.prefix
                    if rule.filter.and_filter.tags:
                        and_f["Tags"] = [
                            {"Key": k, "Value": v}
                            for k, v in rule.filter.and_filter.tags.items()
                        ]
                    f["And"] = and_f
                r["Filter"] = f
            else:
                if rule.prefix is not None:
                    r["Prefix"] = rule.prefix
                else:
                    r["Filter"] = {}

            if rule.transitions:
                r["Transitions"] = [
                    {
                        k: v
                        for k, v in [
                            ("Days", t.days),
                            ("Date", t.date),
                            ("StorageClass", t.storage_class),
                        ]
                        if v is not None
                    }
                    for t in rule.transitions
                ]

            if (
                rule.expiration_days
                or rule.expiration_date
                or rule.expired_object_delete_marker
            ):
                exp: dict[str, Any] = {}
                if rule.expiration_days:
                    exp["Days"] = rule.expiration_days
                if rule.expiration_date:
                    exp["Date"] = rule.expiration_date
                if rule.expired_object_delete_marker:
                    exp["ExpiredObjectDeleteMarker"] = rule.expired_object_delete_marker
                r["Expiration"] = exp

            if rule.noncurrent_version_transitions:
                r["NoncurrentVersionTransitions"] = [
                    {
                        k: v
                        for k, v in [
                            ("NewerNoncurrentVersions", nvt.newer_versions),
                            ("NoncurrentDays", nvt.days),
                            ("StorageClass", nvt.storage_class),
                        ]
                        if v is not None
                    }
                    for nvt in rule.noncurrent_version_transitions
                ]

            if rule.nve_noncurrent_days:
                r["NoncurrentVersionExpiration"] = {
                    "NoncurrentDays": rule.nve_noncurrent_days
                }

            if rule.aimu_days:
                r["AbortIncompleteMultipartUpload"] = {
                    "DaysAfterInitiation": rule.aimu_days
                }

            lifecycle_rules.append(r)

        result = {
            "Rules": lifecycle_rules,
            "TransitionDefaultMinimumObjectSize": TransitionDefaultMinimumObjectSize.ALL_STORAGE_CLASSES_128K.value,
        }
        return self.serialized(ActionResult(result))

    def get_bucket_location(self) -> TYPE_RESPONSE:
        location: str | None = self.backend.get_bucket_location(self.bucket_name)

        # us-east-1 is different - returns a None location
        if location == DEFAULT_REGION_NAME:
            location = None

        self.data["Action"] = "GetBucketLocation"
        return self.serialized(ActionResult({"LocationConstraint": location}))

    def get_bucket_logging(self) -> TYPE_RESPONSE:
        logging = self.backend.get_bucket_logging(self.bucket_name)
        result: dict[str, Any] = {}
        if logging:
            logging_enabled: dict[str, Any] = {
                "TargetBucket": logging["TargetBucket"],
                "TargetPrefix": logging["TargetPrefix"],
            }
            if logging.get("TargetGrants"):
                logging_enabled["TargetGrants"] = [
                    {
                        "Grantee": {
                            "Type": grant.grantees[0].type,
                            "URI": grant.grantees[0].uri,
                            "ID": grant.grantees[0].id,
                            "DisplayName": grant.grantees[0].display_name,
                        },
                        "Permission": grant.permissions[0],
                    }
                    for grant in logging["TargetGrants"]
                ]
            result["LoggingEnabled"] = logging_enabled
        self.data["Action"] = "GetBucketLogging"
        return self.serialized(ActionResult(result))

    def get_bucket_notification(self) -> TYPE_RESPONSE:
        notification_configuration = self.backend.get_bucket_notification_configuration(
            self.bucket_name
        )
        if not notification_configuration:
            return self.serialized(EmptyResult())

        def _build_config(item: Any, arn_key: str) -> dict[str, Any]:
            config: dict[str, Any] = {
                "Id": item.id,
                arn_key: item.arn,
                "Events": item.events,
            }
            if item.filters:
                config["Filter"] = {
                    "Key": {"FilterRules": item.filters["S3Key"]["FilterRule"]}
                }
            return config

        result: dict[str, Any] = {
            "TopicConfigurations": [
                _build_config(t, "TopicArn") for t in notification_configuration.topic
            ],
            "QueueConfigurations": [
                _build_config(q, "QueueArn") for q in notification_configuration.queue
            ],
            "LambdaFunctionConfigurations": [
                _build_config(cf, "LambdaFunctionArn")
                for cf in notification_configuration.cloud_function
            ],
        }
        self.data["Action"] = "GetBucketNotificationConfiguration"
        return self.serialized(ActionResult(result))

    def get_bucket_ownership_controls(self) -> str | TYPE_RESPONSE:
        self.data["Action"] = "GetBucketOwnershipControls"
        ownership_rule = self.backend.get_bucket_ownership_controls(self.bucket_name)
        if not ownership_rule:
            raise OwnershipControlsNotFoundError(bucket_name=self.bucket_name)
        return self.serialized(
            ActionResult(
                {
                    "OwnershipControls": {
                        "Rules": [{"ObjectOwnership": ownership_rule}],
                    },
                }
            )
        )

    def get_bucket_policy(self) -> str | TYPE_RESPONSE:
        policy = self.backend.get_bucket_policy(self.bucket_name)
        if not policy:
            raise NoSuchBucketPolicy(bucket_name=self.bucket_name)
        return 200, {}, policy

    def get_bucket_replication(self) -> str | TYPE_RESPONSE:
        self.data["Action"] = "GetBucketReplication"
        replication = self.backend.get_bucket_replication(self.bucket_name)
        if not replication:
            raise ReplicationConfigurationNotFoundError(bucket_name=self.bucket_name)
        rules = [
            {
                "ID": rule["ID"],
                "Priority": rule["Priority"],
                "Status": rule["Status"],
                "DeleteMarkerReplication": {"Status": "Disabled"},
                "Filter": {"Prefix": ""},
                "Destination": rule["Destination"],
            }
            for rule in replication["Rule"]
        ]
        return self.serialized(
            ActionResult(
                {
                    "ReplicationConfiguration": {
                        "Role": replication["Role"],
                        "Rules": rules,
                    },
                }
            )
        )

    def get_bucket_tags(self) -> str | TYPE_RESPONSE:
        self.data["Action"] = "GetBucketTagging"
        tags = self.backend.get_bucket_tagging(self.bucket_name)["Tags"]
        # "Special Error" if no tags:
        if len(tags) == 0:
            raise NoSuchTagSet(bucket_name=self.bucket_name)
        return self.serialized(ActionResult({"TagSet": tags}))

    def get_bucket_versioning(self) -> TYPE_RESPONSE:
        versioning = self.backend.get_bucket_versioning(self.bucket_name)
        result = {"Status": versioning}
        self.data["Action"] = "GetBucketVersioning"
        return self.serialized(ActionResult(result))

    def get_bucket_website(self) -> TYPE_RESPONSE:
        website_configuration = self.backend.get_bucket_website_configuration(
            self.bucket_name
        )
        if not website_configuration:
            raise NoSuchWebsiteConfiguration(bucket_name=self.bucket_name)
        return 200, {}, website_configuration

    def get_object_acl(self) -> TYPE_RESPONSE:
        response_headers = self._get_cors_headers_other()
        key, not_modified = self._get_key()
        if key.version_id != "null":
            response_headers["x-amz-version-id"] = key.version_id
        if not_modified:
            return 304, response_headers, "Not Modified"

        acl = self.backend.get_object_acl(key)
        self.data["Action"] = "GetObjectAcl"
        status, headers, body = self.serialized(ActionResult(self._acl_to_dict(acl)))
        headers.update(response_headers)
        return status, headers, body

    def get_object_lock_configuration(self) -> TYPE_RESPONSE:
        (
            lock_enabled,
            mode,
            days,
            years,
        ) = self.backend.get_object_lock_configuration(self.bucket_name)
        config: dict[str, Any] = {
            "ObjectLockEnabled": "Enabled" if lock_enabled else "Disabled",
        }
        if mode:
            config["Rule"] = {
                "DefaultRetention": {
                    "Mode": mode,
                    "Days": days,
                    "Years": years,
                },
            }
        self.data["Action"] = "GetObjectLockConfiguration"
        return self.serialized(ActionResult({"ObjectLockConfiguration": config}))

    def get_public_access_block(self) -> TYPE_RESPONSE:
        public_block_config = self.backend.get_public_access_block(self.bucket_name)
        self.data["Action"] = "GetPublicAccessBlock"
        return self.serialized(
            ActionResult(
                {
                    "PublicAccessBlockConfiguration": public_block_config,
                }
            )
        )

    def get_bucket_inventory_configuration(self) -> TYPE_RESPONSE:
        config_id = self.querystring["id"][0]
        inventory_configuration = self.backend.get_bucket_inventory_configuration(
            bucket_name=self.bucket_name, id=config_id
        )
        config: dict[str, Any] = {
            "Destination": inventory_configuration.destination,
            "IsEnabled": ensure_boolean(inventory_configuration.is_enabled),
            "Id": inventory_configuration.id,
            "IncludedObjectVersions": "All",
            "Schedule": inventory_configuration.schedule,
        }
        if inventory_configuration.filters:
            config["Filter"] = inventory_configuration.filters
        if inventory_configuration.optional_fields:
            # TODO: fix input processing to make sure this is parsed as list and not a dict.
            optional_fields: list[str] = (
                inventory_configuration.optional_fields.get("Field", [])
                if isinstance(inventory_configuration.optional_fields, dict)
                else inventory_configuration.optional_fields
            )
            config["OptionalFields"] = optional_fields
        self.data["Action"] = "GetBucketInventoryConfiguration"
        return self.serialized(ActionResult({"InventoryConfiguration": config}))

    def list_bucket_inventory_configurations(self) -> TYPE_RESPONSE:
        inventory_configuration_list = (
            self.backend.list_bucket_inventory_configurations(
                bucket_name=self.bucket_name,
            )
        )
        configs = []
        for inv in inventory_configuration_list:
            config: dict[str, Any] = {
                "Destination": inv.destination,
                "IsEnabled": inv.is_enabled,
                "Id": inv.id,
                "IncludedObjectVersions": "All",
                "Schedule": inv.schedule,
            }
            if inv.filters:
                config["Filter"] = inv.filters
            if inv.optional_fields:
                config["OptionalFields"] = inv.optional_fields
            configs.append(config)

        # TODO: Add support for pagination/ continuation tokens
        result = {
            "InventoryConfigurationList": configs,
            "IsTruncated": False,
        }
        self.data["Action"] = "ListBucketInventoryConfigurations"
        return self.serialized(ActionResult(result))

    def _bucket_response_delete(
        self, bucket_name: str, querystring: dict[str, Any]
    ) -> TYPE_RESPONSE:
        self._set_action("BUCKET", "DELETE", querystring)
        self._authenticate_and_authorize_s3_action(bucket_name=bucket_name)

        if "policy" in querystring:
            return self.delete_bucket_policy()
        elif "tagging" in querystring:
            return self.delete_bucket_tagging()
        elif "website" in querystring:
            return self.delete_bucket_website()
        elif "cors" in querystring:
            return self.delete_bucket_cors()
        elif "lifecycle" in querystring:
            return self.delete_bucket_lifecycle()
        elif "publicAccessBlock" in querystring:
            return self.delete_public_access_block()
        elif "encryption" in querystring:
            return self.delete_bucket_encryption()
        elif "replication" in querystring:
            return self.delete_bucket_replication()
        elif "ownershipControls" in querystring:
            return self.delete_bucket_ownership_controls()

        return self.delete_bucket()

    def delete_bucket(self) -> TYPE_RESPONSE:
        self.data["Action"] = "DeleteBucket"
        removed_bucket = self.backend.delete_bucket(self.bucket_name)
        if removed_bucket:
            # Bucket exists
            return self.serialized(EmptyResult())
        else:
            # Tried to delete a bucket that still has keys
            raise BucketNotEmpty(bucket_name=self.bucket_name)

    def delete_bucket_ownership_controls(self) -> TYPE_RESPONSE:
        self.backend.delete_bucket_ownership_controls(self.bucket_name)
        return 204, {}, ""

    def delete_bucket_replication(self) -> TYPE_RESPONSE:
        self.backend.delete_bucket_replication(self.bucket_name)
        return 204, {}, ""

    def delete_bucket_encryption(self) -> TYPE_RESPONSE:
        self.backend.delete_bucket_encryption(self.bucket_name)
        return 204, {}, ""

    def delete_public_access_block(self) -> TYPE_RESPONSE:
        self.backend.delete_public_access_block(self.bucket_name)
        return 204, {}, ""

    def delete_bucket_lifecycle(self) -> TYPE_RESPONSE:
        self.backend.delete_bucket_lifecycle(self.bucket_name)
        return 204, {}, ""

    def delete_bucket_cors(self) -> TYPE_RESPONSE:
        self.backend.delete_bucket_cors(self.bucket_name)
        return 204, {}, ""

    def delete_bucket_website(self) -> TYPE_RESPONSE:
        self.backend.delete_bucket_website(self.bucket_name)
        return 204, {}, ""

    def delete_bucket_tagging(self) -> TYPE_RESPONSE:
        self.backend.delete_bucket_tagging(self.bucket_name)
        return 204, {}, ""

    def delete_bucket_policy(self) -> TYPE_RESPONSE:
        self.backend.delete_bucket_policy(self.bucket_name)
        return 204, {}, ""

    def _bucket_response_post(self, request: Any, bucket_name: str) -> TYPE_RESPONSE:
        response_headers = {}
        if not request.headers.get("Content-Length"):
            return 411, {}, "Content-Length required"

        self.path = self._get_path(request)

        if self.is_delete_keys():
            self.data["Action"] = "DeleteObject"
            try:
                self._authenticate_and_authorize_s3_action(bucket_name=bucket_name)
                return self._bucket_response_delete_keys(bucket_name)
            except BucketAccessDeniedError:
                return self._bucket_response_delete_keys(
                    bucket_name, authenticated=False
                )

        self.data["Action"] = "PutObject"
        self._authenticate_and_authorize_s3_action(bucket_name=bucket_name)

        key = self.querystring["key"][0]
        f = self.body

        if "success_action_redirect" in self.querystring:
            redirect = self.querystring["success_action_redirect"][0]
            parts = urlparse(redirect)
            queryargs: dict[str, Any] = parse_qs(parts.query)
            queryargs["key"] = key
            queryargs["bucket"] = bucket_name
            redirect_queryargs = urlencode(queryargs, doseq=True)
            newparts = (
                parts.scheme,
                parts.netloc,
                parts.path,
                parts.params,
                redirect_queryargs,
                parts.fragment,
            )
            fixed_redirect = urlunparse(newparts)

            response_headers["Location"] = fixed_redirect

        if "success_action_status" in self.querystring:
            status_code = self.querystring["success_action_status"][0]
        elif "success_action_redirect" in self.querystring:
            status_code = 303
        else:
            status_code = 204

        new_key = self.backend.put_object(
            bucket_name, key, f, request_method=request.method
        )

        if self.querystring.get("acl"):
            acl = get_canned_acl(self.querystring["acl"][0])
            new_key.set_acl(acl)

        # Metadata
        metadata = metadata_from_headers(self.form_data)
        new_key.set_metadata(metadata)

        return status_code, response_headers, ""

    @staticmethod
    def _get_path(request: Any) -> str:  # type: ignore[misc]
        return (
            request.full_path
            if hasattr(request, "full_path")
            else path_url(request.url)
        )

    def _bucket_response_delete_keys(
        self, bucket_name: str, authenticated: bool = True
    ) -> TYPE_RESPONSE:
        body_dict = xmltodict.parse(self.body, strip_whitespace=False)
        bypass_retention = (
            self.headers.get("x-amz-bypass-governance-retention", "").lower() == "true"
        )

        objects = body_dict["Delete"].get("Object", [])
        if not isinstance(objects, list):
            # We expect a list of objects, but when there is a single <Object> node xmltodict does not
            # return a list.
            objects = [objects]
        if len(objects) == 0:
            raise MalformedXML()

        if authenticated:
            objects_to_delete = []
            errors = []
            for obj in objects:
                from moto.iam.access_control import PermissionResult

                bucket = self.backend.get_bucket(bucket_name)
                perm = bucket.get_permission(
                    "s3:DeleteObject", f"arn:aws:s3:::{bucket_name}/{obj['Key']}"
                )
                if perm == PermissionResult.DENIED:
                    errors.append((obj["Key"], "AccessDenied", "Access Denied"))
                else:
                    objects_to_delete.append(obj)
            deleted, errored = self.backend.delete_objects(
                bucket_name, objects_to_delete, bypass_retention
            )
            errors.extend(
                [
                    (
                        err,
                        "AccessDenied",
                        "Access Denied because object protected by object lock.",
                    )
                    for err in errored
                ]
            )
        else:
            deleted = []
            # [(key_name, errorcode, 'error message'), ..]
            errors = [(o["Key"], "AccessDenied", "Access Denied") for o in objects]

        deleted_objects = [
            {
                "Key": k,
                "VersionId": v or None,
                "DeleteMarkerVersionId": dv or None,
                "DeleteMarker": True if dv else None,
            }
            for k, v, dv in deleted
        ]
        error_objects = [{"Key": k, "Code": c, "Message": m} for k, c, m in errors]
        self.data["Action"] = "DeleteObjects"
        return self.serialized(
            ActionResult(
                {
                    "Deleted": deleted_objects,
                    "Errors": error_objects,
                }
            )
        )

    def _handle_range_header(
        self, request: Any, response_headers: dict[str, Any], response_content: Any
    ) -> TYPE_RESPONSE:
        if request.method == "HEAD":
            length = int(response_headers["content-length"])
        else:
            length = len(response_content)
        last = length - 1

        _, rspec = request.headers.get("range").split("=")
        if "," in rspec:
            return 200, response_headers, response_content

        try:
            begin, end = [int(i) if i else None for i in rspec.split("-")]
        except ValueError:
            # if we can't parse the Range header, S3 just treat the request as a non-range request
            return 200, response_headers, response_content

        if (begin is None and end == 0) or (begin is not None and begin > last):
            if request.method == "HEAD":
                return 416, {"Content-Type": APP_XML}, b""
            raise InvalidRange(
                actual_size=str(length), range_requested=request.headers.get("range")
            )

        if begin is not None:  # byte range
            end = last if end is None else min(end, last)
        elif end is not None:  # suffix byte range
            begin = length - min(end, length)
            end = last
        else:
            # Treat as non-range request
            return 200, response_headers, response_content

        if begin > min(end, last):
            # Treat as non-range request if after the logic is applied, the start of the range is greater than the end
            return 200, response_headers, response_content

        if begin or end < last:
            # range requests do not return the checksum
            for key in [h for h in response_headers if h.startswith("x-amz-checksum-")]:
                del response_headers[key]

        response_headers["content-range"] = f"bytes {begin}-{end}/{length}"
        content = response_content[begin : end + 1]
        if request.method == "HEAD":
            response_headers["content-length"] = str((end - begin) + 1)
        else:
            response_headers["content-length"] = str(len(content))
        return 206, response_headers, content

    def _handle_v4_chunk_signatures(self, body: bytes, content_length: int) -> bytes:
        body_io = io.BytesIO(body)
        new_body = bytearray(content_length)
        pos = 0
        line = body_io.readline()
        while line:
            # https://docs.aws.amazon.com/AmazonS3/latest/API/sigv4-streaming.html#sigv4-chunked-body-definition
            # str(hex(chunk-size)) + ";chunk-signature=" + signature + \r\n + chunk-data + \r\n
            chunk_size = int(line[: line.find(b";")].decode("utf8"), 16)
            if chunk_size == 0:
                # Signature 'STREAMING-AWS4-HMAC-SHA256-PAYLOAD-TRAILER' will have trailing parts that define the checksum
                # We can safely ignore those
                break
            new_body[pos : pos + chunk_size] = body_io.read(chunk_size)
            pos = pos + chunk_size
            body_io.read(2)  # skip trailing \r\n
            line = body_io.readline()
        return bytes(new_body)

    def _handle_encoded_body(self, body: io.BufferedIOBase | bytes) -> bytes:
        decoded_body = b""
        if not body:
            return decoded_body
        if isinstance(body, bytes):
            body = io.BytesIO(body)
        # first line should equal '{content_length}\r\n' while the content_length is a hex number
        content_length = int(body.readline().strip(), 16)
        while content_length > 0:
            # read the content_length of the actual data
            decoded_body += body.read(content_length)
            # next is line with just '\r\n' so we skip it
            body.readline()
            # read another line with '{content_length}\r\n'
            content_length = int(body.readline().strip(), 16)

        return decoded_body
        # last line should equal
        # amz-checksum-sha256:<..>\r\n

    def key_response(
        self, request: Any, full_url: str, headers: dict[str, Any]
    ) -> TYPE_RESPONSE:
        # Key and Control are lumped in because splitting out the regex is too much of a pain :/
        self.setup_class(request, full_url, headers)
        bucket_name = self.parse_bucket_name_from_url(request, full_url)
        self.backend.log_incoming_request(request, bucket_name)

        try:
            response = self._key_response(request, full_url)
        except S3ClientError as s3error:
            response = self.serialized(ActionResult(s3error))

        status_code, response_headers, response_content = self._send_response(response)

        if (
            status_code == 200
            and "range" in request.headers
            and request.headers["range"] != ""
        ):
            try:
                return self._handle_range_header(
                    request, response_headers, response_content
                )
            except S3ClientError as s3error:
                return self.serialized(ActionResult(s3error))
        return status_code, response_headers, response_content

    def _key_response(self, request: Any, full_url: str) -> TYPE_RESPONSE:
        parsed_url = urlparse(full_url)
        query = parse_qs(parsed_url.query, keep_blank_values=True)
        method = request.method

        key_name = self.parse_key_name()
        bucket_name = self.parse_bucket_name_from_url(request, full_url)

        # SDK requests tend to have Authorization set automatically
        # If users make an HTTP-request, such as `requests.get("https://bucket-name.s3.amazonaws.com/file-name")`,
        # The authorization-header may not be set
        authorized_request = "Authorization" in request.headers
        if hasattr(request, "url"):
            signed_url = "Signature=" in request.url
        elif hasattr(request, "requestline"):
            signed_url = "Signature=" in request.path
        try:
            key = self.backend.get_object(bucket_name, key_name)
            bucket = self.backend.get_bucket(bucket_name)
        except S3ClientError:
            key = bucket = None

        # Enforces policy when creating new objects, `if key` does not https://github.com/getmoto/moto/issues/7837
        if bucket:
            resource = f"arn:{bucket.partition}:s3:::{bucket_name}/{key_name}"

            # Authorization Workflow
            # https://docs.aws.amazon.com/AmazonS3/latest/userguide/access-control-auth-workflow-object-operation.html

            # A bucket can deny all actions, regardless of who makes the request
            from moto.iam.access_control import PermissionResult

            action = f"s3:{method.upper()[0]}{method.lower()[1:]}Object"
            bucket_permissions = bucket.get_permission(action, resource)
            if bucket_permissions == PermissionResult.DENIED:
                return 403, {}, ""

            if key:
                # If the request is not authorized, and not signed,
                # that means that the action should be allowed for anonymous users
                if not authorized_request and not signed_url:
                    # We already know that the bucket permissions do not explicitly deny this
                    # So bucket permissions are either not set, or do not explicitly allow
                    # Next check is to see if the ACL of the individual key allows this action
                    if bucket_permissions != PermissionResult.PERMITTED and (
                        key.acl and not key.acl.public_read
                    ):
                        return 403, {}, ""
        if not key and signed_url and not authorized_request:
            # coming in from requests.get(s3.generate_presigned_url())
            if self._invalid_headers(request.url, dict(request.headers)):
                raise SignatureDoesNotMatchError()

        body = self.body or b""

        if method == "GET":
            return self._key_response_get(query, key_name)
        elif method == "PUT":
            return self._key_response_put(query, key_name)
        elif method == "HEAD":
            return self.head_object(
                bucket_name, query, key_name, headers=request.headers
            )
        elif method == "DELETE":
            return self._key_response_delete(bucket_name, query, key_name)
        elif method == "POST":
            return self._key_response_post(request, body, bucket_name, query, key_name)
        elif method == "OPTIONS":
            # OPTIONS response doesn't depend on the key_name: always return 200 with CORS headers
            return self._response_options(request.headers, bucket_name)
        else:
            raise NotImplementedError(
                f"Method {method} has not been implemented in the S3 backend yet"
            )

    def _key_response_get(self, query: dict[str, Any], key_name: str) -> TYPE_RESPONSE:
        self._set_action("KEY", "GET", query)
        self._authenticate_and_authorize_s3_action(
            bucket_name=self.bucket_name, key_name=key_name
        )

        if query.get("uploadId"):
            return self.list_parts()

        if "acl" in query:
            return self.get_object_acl()
        if "tagging" in query:
            return self.get_object_tagging()
        if "legal-hold" in query:
            return self.get_object_legal_hold()
        if "attributes" in query:
            return self.get_object_attributes()

        return self.get_object()

    def get_object(self) -> TYPE_RESPONSE:
        key, not_modified = self._get_key()
        response_headers = self._get_cors_headers_other()

        if key.version_id != "null":
            response_headers["x-amz-version-id"] = key.version_id

        response_headers.update(key.response_dict)

        if not_modified:
            # Real S3 omits any content-* headers for a 304
            for header in list(response_headers.keys()):
                if header.startswith("content-"):
                    response_headers.pop(header)
            return 304, response_headers, "Not Modified"

        # set the checksum after not_modified has been checked
        if (
            self.headers.get("x-amz-checksum-mode") == "ENABLED"
            and key.checksum_algorithm
        ):
            qualified_checksum = key.checksum_value
            if key.checksum_type == "COMPOSITE":
                qualified_checksum = f"{key.checksum_value}-{key.checksum_parts}"

            response_headers[f"x-amz-checksum-{key.checksum_algorithm.lower()}"] = (
                qualified_checksum
            )

            if key.checksum_type:
                response_headers["x-amz-checksum-type"] = key.checksum_type

        response_headers.update(key.metadata)
        response_headers.update({"Accept-Ranges": "bytes"})

        for param_name, header_name in ALLOWED_HEADER_OVERRIDES.items():
            if param_name in self.querystring:
                response_headers[header_name] = self.querystring[param_name][0]

        part_number = self._get_int_param("partNumber")
        if part_number:
            if key.multipart:
                # TODO
                pass
            else:
                if part_number > 1:
                    raise RangeNotSatisfiable
                response_headers["content-range"] = f"bytes 0-{key.size - 1}/{key.size}"
                return 206, response_headers, key.value
        return 200, response_headers, key.value

    def get_object_attributes(self) -> TYPE_RESPONSE:
        # Get the Key, but do not validate StorageClass - we can retrieve the attributes of Glacier-objects
        key, not_modified = self._get_key(validate_storage_class=False)
        response_headers = self._get_cors_headers_other()
        if not_modified:
            return 304, response_headers, "Not Modified"

        if key.version_id != "null":
            response_headers["x-amz-version-id"] = key.version_id

        attributes_to_get = self.headers.get("x-amz-object-attributes", "").split(",")
        response_keys = self.backend.get_object_attributes(key, attributes_to_get)
        response_headers["Last-Modified"] = key.last_modified_ISO8601

        checksum = response_keys["checksum"]
        result: dict[str, Any] = {
            "ETag": response_keys["etag"],
            "Checksum": {f"Checksum{algo}": value for algo, value in checksum.items()}
            if checksum
            else None,
            "ObjectSize": response_keys["size"],
            "StorageClass": response_keys["storage_class"],
        }
        if checksum:
            result["Checksum"]["ChecksumType"] = checksum.get("type")
        self.data["Action"] = "GetObjectAttributes"
        status, headers, body = self.serialized(ActionResult(result))
        headers.update(response_headers)
        return status, headers, body

    def get_object_legal_hold(self) -> TYPE_RESPONSE:
        key, not_modified = self._get_key()
        response_headers = self._get_cors_headers_other()
        if not_modified:
            return 304, response_headers, "Not Modified"

        if key.version_id != "null":
            response_headers["x-amz-version-id"] = key.version_id

        legal_hold = self.backend.get_object_legal_hold(key)
        self.data["Action"] = "GetObjectLegalHold"
        status, headers, body = self.serialized(
            ActionResult(
                {
                    "LegalHold": {"Status": legal_hold},
                }
            )
        )
        headers.update(response_headers)
        return status, headers, body

    def get_object_tagging(self) -> TYPE_RESPONSE:
        key, not_modified = self._get_key()
        response_headers = self._get_cors_headers_other()
        if not_modified:
            return 304, response_headers, "Not Modified"

        if key.version_id != "null":
            response_headers["x-amz-version-id"] = key.version_id

        tags = self.backend.get_object_tagging(key)["Tags"]

        self.data["Action"] = "GetObjectTagging"
        status, headers, body = self.serialized(ActionResult({"TagSet": tags}))
        headers.update(response_headers)
        return status, headers, body

    def list_parts(self) -> TYPE_RESPONSE:
        response_headers = self._get_cors_headers_other()

        upload_id = self._get_param("uploadId")

        # 0 <= PartNumberMarker <= 2,147,483,647
        part_number_marker = self._get_int_param("part-number-marker", 0)
        if part_number_marker > 2147483647:
            raise NotAnIntegerException(
                name="part-number-marker", value=part_number_marker
            )
        if not (0 <= part_number_marker <= 2147483647):
            raise InvalidMaxPartArgument("part-number-marker", 0, 2147483647)

        # 0 <= MaxParts <= 2,147,483,647 (default is 1,000)
        max_parts = self._get_int_param("max-parts", 1000)
        if max_parts > 2147483647:
            raise NotAnIntegerException(name="max-parts", value=max_parts)
        if not (0 <= max_parts <= 2147483647):
            raise InvalidMaxPartArgument("max-parts", 0, 2147483647)

        parts = self.backend.list_parts(
            self.bucket_name,
            upload_id,
            part_number_marker=part_number_marker,
            max_parts=max_parts,
        )
        next_part_number_marker = parts[-1].name if parts else 0
        is_truncated = len(parts) != 0 and self.backend.is_truncated(
            self.bucket_name,
            upload_id,
            next_part_number_marker,  # type: ignore
        )

        key_name = self.parse_key_name()
        parts_list = [
            {
                "PartNumber": part.name,
                "LastModified": part.last_modified_ISO8601,
                "ETag": part.etag,
                "Size": part.size,
            }
            for part in parts
        ]
        result = {
            "Bucket": self.bucket_name,
            "Key": key_name,
            "UploadId": upload_id,
            "StorageClass": "STANDARD",
            "Initiator": {
                "ID": "75aa57f09aa0c8caeab4f8c24e99d10f8e7faeebf76c078efc7c6caea54ba06a",
                "DisplayName": "webfile",
            },
            "Owner": {
                "ID": "75aa57f09aa0c8caeab4f8c24e99d10f8e7faeebf76c078efc7c6caea54ba06a",
                "DisplayName": "webfile",
            },
            "PartNumberMarker": part_number_marker,
            "NextPartNumberMarker": next_part_number_marker,
            "MaxParts": max_parts,
            "IsTruncated": is_truncated,
            "Parts": parts_list,
        }
        self.data["Action"] = "ListParts"
        status, headers, body = self.serialized(ActionResult(result))
        headers.update(response_headers)
        return status, headers, body

    def _get_key(self, validate_storage_class: bool = True) -> tuple[FakeKey, bool]:
        key_name = self.parse_key_name()
        version_id = self.querystring.get("versionId", [None])[0]
        if_modified_since = self.headers.get("If-Modified-Since")
        if_match = self.headers.get("If-Match")
        if_none_match = self.headers.get("If-None-Match")
        if_unmodified_since = self.headers.get("If-Unmodified-Since")

        key = self.backend.get_object(self.bucket_name, key_name, version_id=version_id)
        if key is None and version_id is None:
            raise MissingKey(key=key_name)
        elif key is None:
            raise MissingVersion()
        if validate_storage_class and key.storage_class in ARCHIVE_STORAGE_CLASSES:
            if 'ongoing-request="false"' not in key.response_dict.get(
                "x-amz-restore", ""
            ):
                raise InvalidObjectState(storage_class=key.storage_class)
        if if_unmodified_since:
            if_unmodified_since = str_to_rfc_1123_datetime(if_unmodified_since)
            if key.last_modified.replace(microsecond=0) > if_unmodified_since:
                raise PreconditionFailed("If-Unmodified-Since")
        if if_match and key.etag not in [if_match, f'"{if_match}"']:
            raise PreconditionFailed("If-Match")
        not_modified = False
        if if_modified_since:
            if_modified_since = str_to_rfc_1123_datetime(if_modified_since)
            if key.last_modified.replace(microsecond=0) <= if_modified_since:
                not_modified = True
        if if_none_match and key.etag in [if_none_match, f'"{if_none_match}"']:
            not_modified = True
        return key, not_modified

    def _key_response_put(self, query: dict[str, Any], key_name: str) -> TYPE_RESPONSE:
        self._set_action("KEY", "PUT", query)
        self._authenticate_and_authorize_s3_action(
            bucket_name=self.bucket_name, key_name=key_name
        )

        if query.get("uploadId") and query.get("partNumber"):
            if "x-amz-copy-source" in self.headers:
                return self.upload_part_copy()
            else:
                return self.upload_part()

        bucket = self.backend.get_bucket(self.bucket_name)
        lock_enabled = bucket.object_lock_enabled

        lock_mode = self.headers.get("x-amz-object-lock-mode")
        lock_until = self.headers.get("x-amz-object-lock-retain-until-date", None)
        legal_hold = self.headers.get("x-amz-object-lock-legal-hold")

        if lock_mode or lock_until or legal_hold == "ON":
            if not self.headers.get("Content-Md5") and not self.headers.get(
                "x-amz-sdk-checksum-algorithm"
            ):
                raise MissingUploadObjectWithObjectLockHeaders
            if not lock_enabled:
                raise LockNotEnabled

        elif lock_enabled and bucket.has_default_lock:
            if not self.headers.get("Content-Md5") and not self.headers.get(
                "x-amz-sdk-checksum-algorithm"
            ):
                raise MissingUploadObjectWithObjectLockHeaders

        if "retention" in query:
            return self.put_object_retention()

        if "legal-hold" in query:
            return self.put_object_legal_hold()

        if "acl" in query:
            return self.put_object_acl()

        if "tagging" in query:
            return self.put_object_tagging()

        if "x-amz-copy-source" in self.headers:
            return self.copy_object()

        return self.put_object()

    def put_object(self) -> TYPE_RESPONSE:
        key_name = self.parse_key_name()
        response_headers = self._get_cors_headers_other()

        storage_class = self.headers.get("x-amz-storage-class", "STANDARD")
        encryption = self.headers.get("x-amz-server-side-encryption")
        kms_key_id = self.headers.get("x-amz-server-side-encryption-aws-kms-key-id")
        if_match = self.headers.get("If-Match")
        if_none_match = self.headers.get("If-None-Match")
        bucket_key_enabled = self.headers.get(
            "x-amz-server-side-encryption-bucket-key-enabled"
        )
        if bucket_key_enabled is not None:
            bucket_key_enabled = str(bucket_key_enabled).lower()

        if if_match:
            if not (obj := self.backend.get_object(self.bucket_name, key_name)):
                raise MissingKey
            # Check if the ETags are the same. S3 doesn't seem to care about quotes, so we shouldn't either
            elif if_match.replace('"', "") != obj.etag.replace('"', ""):
                raise PreconditionFailed("If-Match")
        if (
            if_none_match == "*"
            and self.backend.get_object(self.bucket_name, key_name) is not None
        ):
            raise PreconditionFailed("If-None-Match")

        checksum_algorithm, checksum_value = self._get_checksum(response_headers)

        bucket = self.backend.get_bucket(self.bucket_name)
        lock_enabled = bucket.object_lock_enabled

        legal_hold, lock_mode, lock_until = self._get_lock_details(bucket, lock_enabled)

        acl = self._acl_from_headers(self.headers)
        if acl is None:
            acl = bucket.acl
        tagging = self._tagging_from_headers(self.headers)

        new_key = self.backend.put_object(
            self.bucket_name,
            key_name,
            self.body or b"",
            storage=storage_class,
            encryption=encryption,
            kms_key_id=kms_key_id,
            bucket_key_enabled=bucket_key_enabled,
            lock_mode=lock_mode,
            lock_legal_status=legal_hold,
            lock_until=lock_until,
            checksum_value=checksum_value,
        )
        metadata = metadata_from_headers(self.headers)
        metadata.update(metadata_from_headers(self.querystring))
        new_key.set_metadata(metadata)
        new_key.set_acl(acl)
        new_key.website_redirect_location = self.headers.get(
            "x-amz-website-redirect-location"
        )
        if checksum_algorithm:
            new_key.checksum_algorithm = checksum_algorithm
        self.backend.put_object_tagging(new_key, tagging)
        response_headers.update(new_key.response_dict)
        # Remove content-length - the response body is empty for this request
        response_headers.pop("content-length", None)
        return 200, response_headers, ""

    def _get_checksum(self, response_headers: dict[str, Any]) -> tuple[str, str | None]:
        checksum_algorithm = self.headers.get("x-amz-sdk-checksum-algorithm", "")
        checksum_header = f"x-amz-checksum-{checksum_algorithm.lower()}"
        checksum_value = self.headers.get(checksum_header)
        if not checksum_value and checksum_algorithm:
            # Extract the checksum-value from the body first
            search = re.search(rb"x-amz-checksum-\w+:(.+={1,2})", self.raw_body)
            checksum_value = search.group(1) if search else None
        if checksum_value:
            # TODO: AWS computes the provided value and verifies it's the same
            # Afterwards, it should be returned in every subsequent call
            if isinstance(checksum_value, bytes):
                checksum_value = checksum_value.decode("utf-8")
            response_headers.update({checksum_header: checksum_value})
        elif checksum_algorithm:
            # If the value is not provided, we compute it and only return it as part of this request
            checksum_value = compute_checksum(
                self.raw_body, algorithm=checksum_algorithm
            )
            if isinstance(checksum_value, bytes):
                checksum_value = checksum_value.decode("utf-8")
            response_headers.update({checksum_header: checksum_value})
        return checksum_algorithm, checksum_value

    def copy_object(self) -> TYPE_RESPONSE:
        # you can have a quoted ?version=abc with a version Id, so work on
        # we need to parse the unquoted string first
        copy_source = self.headers.get("x-amz-copy-source")
        if isinstance(copy_source, bytes):
            copy_source = copy_source.decode("utf-8")
        copy_source_parsed = urlparse(copy_source)
        src_bucket, src_key = unquote(copy_source_parsed.path).lstrip("/").split("/", 1)
        src_version_id = parse_qs(copy_source_parsed.query).get("versionId", [None])[0]

        key_to_copy = self.backend.get_object(
            src_bucket, src_key, version_id=src_version_id
        )
        key_name = self.parse_key_name()

        bucket = self.backend.get_bucket(self.bucket_name)
        lock_enabled = bucket.object_lock_enabled

        encryption = self.headers.get("x-amz-server-side-encryption")
        kms_key_id = self.headers.get("x-amz-server-side-encryption-aws-kms-key-id")
        bucket_key_enabled = self.headers.get(
            "x-amz-server-side-encryption-bucket-key-enabled"
        )
        if bucket_key_enabled is not None:
            bucket_key_enabled = str(bucket_key_enabled).lower()

        legal_hold, lock_mode, lock_until = self._get_lock_details(bucket, lock_enabled)

        if key_to_copy is not None:
            if "x-amz-copy-source-if-none-match" in self.headers:
                requested_etag = self.headers["x-amz-copy-source-if-none-match"]
                if requested_etag in [key_to_copy.etag, key_to_copy.etag[1:-1]]:
                    raise PreconditionFailed(
                        failed_condition="x-amz-copy-source-If-None-Match"
                    )

            if key_to_copy.storage_class in ARCHIVE_STORAGE_CLASSES:
                if key_to_copy.response_dict.get(
                    "x-amz-restore"
                ) is None or 'ongoing-request="true"' in key_to_copy.response_dict.get(  # type: ignore
                    "x-amz-restore"
                ):
                    raise ObjectNotInActiveTierError(key_to_copy)

            website_redirect_location = self.headers.get(
                "x-amz-website-redirect-location"
            )

            mdirective = self.headers.get("x-amz-metadata-directive")
            metadata = metadata_from_headers(self.headers)
            self.backend.copy_object(
                key_to_copy,
                self.bucket_name,
                key_name,
                storage=self.headers.get("x-amz-storage-class"),
                kms_key_id=kms_key_id,
                encryption=encryption,
                bucket_key_enabled=bucket_key_enabled,
                mdirective=mdirective,
                metadata=metadata,
                website_redirect_location=website_redirect_location,
                lock_mode=lock_mode,
                lock_legal_status=legal_hold,
                lock_until=lock_until,
                provided_version_id=src_version_id,
            )
        else:
            if src_version_id:
                raise MissingVersion()
            raise MissingKey(key=src_key)

        new_key: FakeKey = self.backend.get_object(self.bucket_name, key_name)  # type: ignore

        acl = self._acl_from_headers(self.headers)
        if acl is None:
            acl = bucket.acl
        if acl is not None:
            new_key.set_acl(acl)

        response_headers = self._get_cors_headers_other()
        tdirective = self.headers.get("x-amz-tagging-directive")
        if tdirective == "REPLACE":
            tagging = self._tagging_from_headers(self.headers)
            self.backend.put_object_tagging(new_key, tagging)
        if key_to_copy.version_id != "null":
            response_headers["x-amz-copy-source-version-id"] = key_to_copy.version_id

        # checksum stuff, do we need to compute hash of the copied object
        checksum_algorithm = self.headers.get("x-amz-checksum-algorithm")
        if checksum_algorithm:
            checksum_value = compute_checksum(
                new_key.value, algorithm=checksum_algorithm
            ).decode("utf-8")
            response_headers.update(
                {"Checksum": {f"Checksum{checksum_algorithm}": checksum_value}}
            )
            # By default, the checksum-details for the copy will be the same as the original
            # But if another algorithm is provided during the copy-operation, we override the values
            new_key.checksum_algorithm = checksum_algorithm
            new_key.checksum_value = checksum_value

        copy_result: dict[str, Any] = {
            "ETag": new_key.etag,
            "LastModified": new_key.last_modified_ISO8601,
        }
        if new_key.checksum_value:
            copy_result[f"Checksum{new_key.checksum_algorithm}"] = (
                new_key.checksum_value
            )

        self.data["Action"] = "CopyObject"
        status, headers, body = self.serialized(
            ActionResult({"CopyObjectResult": copy_result})
        )
        headers.update(response_headers)
        # TODO: Remove this in favor of setting the appropriate keys in the
        # result dict, e.g. "ServerSideEncryption": new_key.encryption
        headers.update(new_key.response_dict)
        headers["content-length"] = str(len(body))
        return status, headers, body

    def _get_lock_details(
        self, bucket: "FakeBucket", lock_enabled: bool
    ) -> tuple[str | None, str | None, str | None]:
        lock_mode = self.headers.get("x-amz-object-lock-mode")
        lock_until = self.headers.get("x-amz-object-lock-retain-until-date")
        legal_hold = self.headers.get("x-amz-object-lock-legal-hold")
        if lock_mode or lock_until or legal_hold == "ON":
            pass
        elif lock_enabled and bucket.has_default_lock:
            lock_until = bucket.default_retention()
            lock_mode = bucket.default_lock_mode
        return legal_hold, lock_mode, lock_until

    def put_object_tagging(self) -> TYPE_RESPONSE:
        key_name = self.parse_key_name()
        version_id = self._get_param("versionId")
        key_to_tag = self.backend.get_object(
            self.bucket_name, key_name, version_id=version_id, return_delete_marker=True
        )

        if isinstance(key_to_tag, FakeDeleteMarker):
            raise MethodNotAllowed(method="PUT", resource_type="DeleteMarker")
        elif key_to_tag is None and version_id is None:
            raise MissingKey(key=key_name)
        elif key_to_tag is None:
            raise MissingVersion(key=key_name, version_id=version_id)

        tagging = self._tagging_from_xml()
        self.backend.put_object_tagging(key=key_to_tag, tags=tagging)

        response_headers = self._get_cors_headers_other()
        self._get_checksum(response_headers)
        return 200, response_headers, ""

    def put_object_acl(self) -> TYPE_RESPONSE:
        acl = self._acl_from_headers(self.headers)
        if acl is None:
            bucket = self.backend.get_bucket(self.bucket_name)
            acl = bucket.acl
        key_name = self.parse_key_name()
        self.backend.put_object_acl(self.bucket_name, key_name, acl)

        response_headers = self._get_cors_headers_other()
        self._get_checksum(response_headers)
        return 200, response_headers, ""

    def put_object_legal_hold(self) -> TYPE_RESPONSE:
        version_id = self._get_param("versionId")
        bucket = self.backend.get_bucket(self.bucket_name)
        lock_enabled = bucket.object_lock_enabled
        key_name = self.parse_key_name()

        if not lock_enabled:
            raise LockNotEnabled
        legal_hold_status = self._legal_hold_status_from_xml(self.body)
        self.backend.put_object_legal_hold(
            self.bucket_name, key_name, version_id, legal_hold_status
        )

        response_headers = self._get_cors_headers_other()
        self._get_checksum(response_headers)
        return 200, response_headers, ""

    def put_object_retention(self) -> TYPE_RESPONSE:
        version_id = self._get_param("versionId")
        bucket = self.backend.get_bucket(self.bucket_name)
        lock_enabled = bucket.object_lock_enabled
        key_name = self.parse_key_name()

        if not lock_enabled:
            raise LockNotEnabled
        retention = self._mode_until_from_body()
        self.backend.put_object_retention(
            self.bucket_name, key_name, version_id=version_id, retention=retention
        )

        response_headers = self._get_cors_headers_other()
        self._get_checksum(response_headers)
        return 200, response_headers, ""

    def upload_part(self) -> TYPE_RESPONSE:
        upload_id = self._get_param("uploadId")
        part_number = self._get_int_param("partNumber")

        if part_number > 10000:
            raise InvalidMaxPartNumberArgument(part_number)
        key = self.backend.upload_part(
            self.bucket_name, upload_id, part_number, self.body
        )

        response_headers = self._get_cors_headers_other()
        self._get_checksum(response_headers)
        response_headers.update(key.response_dict)
        response_headers["content-length"] = "0"
        return 200, response_headers, ""

    def upload_part_copy(self) -> TYPE_RESPONSE:
        response_headers = self._get_cors_headers_other()
        self._get_checksum(response_headers)

        upload_id = self._get_param("uploadId")
        part_number = self._get_int_param("partNumber")

        copy_source = self.headers.get("x-amz-copy-source")
        if isinstance(copy_source, bytes):
            copy_source = copy_source.decode("utf-8")
        copy_source_parsed = urlparse(copy_source)
        src_bucket, src_key = unquote(copy_source_parsed.path).lstrip("/").split("/", 1)
        src_version_id = parse_qs(copy_source_parsed.query).get("versionId", [None])[0]
        src_range = self.headers.get("x-amz-copy-source-range", "").split("bytes=")[-1]
        try:
            start_byte, end_byte = src_range.split("-")
            start_byte, end_byte = int(start_byte), int(end_byte)
        except ValueError:
            start_byte, end_byte = None, None
        if self.backend.get_object(src_bucket, src_key, version_id=src_version_id):
            key = self.backend.upload_part_copy(
                self.bucket_name,
                upload_id,
                part_number,
                src_bucket_name=src_bucket,
                src_key_name=src_key,
                src_version_id=src_version_id,
                start_byte=start_byte,
                end_byte=end_byte,
            )
        else:
            return 404, response_headers, ""
        self.data["Action"] = "UploadPartCopy"
        status, headers, body = self.serialized(
            ActionResult(
                {
                    "CopyPartResult": {
                        "LastModified": key.last_modified_ISO8601,
                        "ETag": key.etag,
                    },
                }
            )
        )
        headers.update(response_headers)
        headers.update(key.response_dict)
        headers["content-length"] = str(len(body))
        return status, headers, body

    def head_object(
        self,
        bucket_name: str,
        query: dict[str, Any],
        key_name: str,
        headers: dict[str, Any],
    ) -> TYPE_RESPONSE:
        self._set_action("KEY", "HEAD", query)
        self._authenticate_and_authorize_s3_action(
            bucket_name=bucket_name, key_name=key_name
        )

        response_headers: dict[str, Any] = {}
        version_id = query.get("versionId", [None])[0]
        if version_id and not self.backend.get_bucket(bucket_name).is_versioned:
            return 400, response_headers, ""

        part_number = query.get("partNumber", [None])[0]
        if part_number:
            part_number = int(part_number)

        checksum_mode = headers.get("x-amz-checksum-mode") == "ENABLED"
        if_modified_since = headers.get("If-Modified-Since", None)
        if_match = headers.get("If-Match", None)
        if_none_match = headers.get("If-None-Match", None)
        if_unmodified_since = headers.get("If-Unmodified-Since", None)

        try:
            key = self.backend.head_object(
                bucket_name, key_name, version_id=version_id, part_number=part_number
            )
        except HeadOnDeleteMarker as exc:
            headers = {
                "x-amz-delete-marker": "true",
                "x-amz-version-id": version_id,
                "content-type": APP_XML,
            }
            if version_id:
                headers["allow"] = "DELETE"
                return 405, headers, "Method Not Allowed"
            else:
                headers["x-amz-version-id"] = exc.marker.version_id
                return 404, headers, "Not Found"
        if key:
            response_headers.update(key.metadata)
            response_headers.update(key.response_dict)
            response_headers.update({"Accept-Ranges": "bytes"})

            if if_unmodified_since:
                if_unmodified_since = str_to_rfc_1123_datetime(if_unmodified_since)
                if key.last_modified.replace(microsecond=0) > if_unmodified_since:
                    return 412, response_headers, ""
            if if_match and key.etag != if_match:
                return 412, response_headers, ""

            if if_modified_since:
                if_modified_since = str_to_rfc_1123_datetime(if_modified_since)
                if key.last_modified.replace(microsecond=0) <= if_modified_since:
                    return 304, response_headers, "Not Modified"
            if if_none_match and key.etag == if_none_match:
                return 304, response_headers, "Not Modified"

            if checksum_mode and key.checksum_algorithm:
                response_headers[f"x-amz-checksum-{key.checksum_algorithm.lower()}"] = (
                    key.checksum_value
                )

            if part_number:
                full_key = self.backend.head_object(bucket_name, key_name, version_id)
                if full_key.multipart:  # type: ignore
                    mp_part_count = str(len(full_key.multipart.partlist))  # type: ignore
                    response_headers["x-amz-mp-parts-count"] = mp_part_count
                else:
                    if part_number > 1:
                        raise RangeNotSatisfiable
                    response_headers["content-range"] = (
                        f"bytes 0-{full_key.size - 1}/{full_key.size}"  # type: ignore
                    )
                    return 206, response_headers, ""

            return 200, response_headers, ""
        else:
            return 404, response_headers, ""

    def _process_lock_config_from_body(self) -> dict[str, Any]:
        try:
            return self._lock_config_from_body()
        except (TypeError, ExpatError):
            raise MissingRequestBody
        except KeyError:
            raise MalformedXML

    def _lock_config_from_body(self) -> dict[str, Any]:
        response_dict: dict[str, Any] = {
            "enabled": False,
            "mode": None,
            "days": None,
            "years": None,
        }
        parsed_xml = xmltodict.parse(self.body)
        enabled = (
            parsed_xml["ObjectLockConfiguration"]["ObjectLockEnabled"] == "Enabled"
        )
        response_dict["enabled"] = enabled

        default_retention = parsed_xml.get("ObjectLockConfiguration").get("Rule")
        if default_retention:
            default_retention = default_retention.get("DefaultRetention")
            mode = default_retention["Mode"]
            days = int(default_retention.get("Days", 0))
            years = int(default_retention.get("Years", 0))

            if days and years:
                raise MalformedXML
            response_dict["mode"] = mode
            response_dict["days"] = days
            response_dict["years"] = years

        return response_dict

    def _acl_from_body(self) -> FakeAcl | None:
        parsed_xml = xmltodict.parse(self.body)
        if not parsed_xml.get("AccessControlPolicy"):
            raise MalformedACLError()

        # The owner is needed for some reason...
        if not parsed_xml["AccessControlPolicy"].get("Owner"):
            # TODO: Validate that the Owner is actually correct.
            raise MalformedACLError()

        # If empty, then no ACLs:
        if parsed_xml["AccessControlPolicy"].get("AccessControlList") is None:
            return None

        if not parsed_xml["AccessControlPolicy"]["AccessControlList"].get("Grant"):
            raise MalformedACLError()

        permissions = ["READ", "WRITE", "READ_ACP", "WRITE_ACP", "FULL_CONTROL"]

        if not isinstance(
            parsed_xml["AccessControlPolicy"]["AccessControlList"]["Grant"], list
        ):
            parsed_xml["AccessControlPolicy"]["AccessControlList"]["Grant"] = [
                parsed_xml["AccessControlPolicy"]["AccessControlList"]["Grant"]
            ]

        grants = self._get_grants_from_xml(
            parsed_xml["AccessControlPolicy"]["AccessControlList"]["Grant"],
            MalformedACLError,
            permissions,
        )
        return FakeAcl(grants)

    def _get_grants_from_xml(
        self,
        grant_list: list[dict[str, Any]],
        exception_type: type[S3ClientError],
        permissions: list[str],
    ) -> list[FakeGrant]:
        grants = []
        for grant in grant_list:
            if grant.get("Permission", "") not in permissions:
                raise exception_type()

            if grant["Grantee"].get("@xsi:type", "") not in [
                "CanonicalUser",
                "AmazonCustomerByEmail",
                "Group",
            ]:
                raise exception_type()

            # TODO: Verify that the proper grantee data is supplied based on the type.

            grants.append(
                FakeGrant(
                    [
                        FakeGrantee(
                            grantee_id=grant["Grantee"].get("ID", ""),
                            display_name=grant["Grantee"].get("DisplayName", ""),
                            uri=grant["Grantee"].get("URI", ""),
                        )
                    ],
                    [grant["Permission"]],
                )
            )

        return grants

    def _acl_from_headers(self, headers: dict[str, str]) -> FakeAcl | None:
        canned_acl = headers.get("x-amz-acl", "")

        grants = []
        for header, value in headers.items():
            header = header.lower()
            if not header.startswith("x-amz-grant-"):
                continue

            permission = {
                "read": "READ",
                "write": "WRITE",
                "read-acp": "READ_ACP",
                "write-acp": "WRITE_ACP",
                "full-control": "FULL_CONTROL",
            }[header[len("x-amz-grant-") :]]

            grantees = []
            for key_and_value in value.split(","):
                key, value = re.match(  # type: ignore
                    '([^=]+)="?([^"]+)"?', key_and_value.strip()
                ).groups()
                if key.lower() == "id":
                    grantees.append(FakeGrantee(grantee_id=value))
                else:
                    grantees.append(FakeGrantee(uri=value))
            grants.append(FakeGrant(grantees, [permission]))

        if canned_acl and grants:
            raise S3AclAndGrantError()
        if canned_acl:
            return get_canned_acl(canned_acl)
        if grants:
            return FakeAcl(grants)
        else:
            return None

    def _tagging_from_headers(self, headers: dict[str, Any]) -> dict[str, str]:
        tags = {}
        if headers.get("x-amz-tagging"):
            parsed_header = parse_qs(headers["x-amz-tagging"], keep_blank_values=True)
            for tag in parsed_header.items():
                tags[tag[0]] = tag[1][0]
        return tags

    def _tagging_from_xml(self) -> dict[str, str]:
        parsed_xml = xmltodict.parse(self.body, force_list={"Tag": True})

        tags = {}
        for tag in parsed_xml["Tagging"]["TagSet"]["Tag"]:
            tags[tag["Key"]] = tag["Value"] or ""

        return tags

    def _bucket_tagging_from_body(self) -> dict[str, str]:
        parsed_xml = xmltodict.parse(self.body)

        tags = {}
        # Optional if no tags are being sent:
        if parsed_xml["Tagging"].get("TagSet"):
            # If there is only 1 tag, then it's not a list:
            if not isinstance(parsed_xml["Tagging"]["TagSet"]["Tag"], list):
                tags[parsed_xml["Tagging"]["TagSet"]["Tag"]["Key"]] = parsed_xml[
                    "Tagging"
                ]["TagSet"]["Tag"]["Value"]
            else:
                for tag in parsed_xml["Tagging"]["TagSet"]["Tag"]:
                    if tag["Key"] in tags:
                        raise DuplicateTagKeys()
                    tags[tag["Key"]] = tag["Value"]

        # Verify that "aws:" is not in the tags. If so, then this is a problem:
        for key, _ in tags.items():
            if key.startswith("aws:"):
                raise NoSystemTags()

        return tags

    def _cors_from_body(self) -> list[dict[str, Any]]:
        parsed_xml = xmltodict.parse(self.body)

        if isinstance(parsed_xml["CORSConfiguration"]["CORSRule"], list):
            return list(parsed_xml["CORSConfiguration"]["CORSRule"])

        return [parsed_xml["CORSConfiguration"]["CORSRule"]]

    def _mode_until_from_body(self) -> tuple[str | None, str | None]:
        parsed_xml = xmltodict.parse(self.body)
        return (
            parsed_xml.get("Retention", None).get("Mode", None),
            parsed_xml.get("Retention", None).get("RetainUntilDate", None),
        )

    def _legal_hold_status_from_xml(self, xml: bytes) -> dict[str, Any]:
        parsed_xml = xmltodict.parse(xml)
        return parsed_xml["LegalHold"]["Status"]

    def _encryption_config_from_body(self) -> dict[str, Any]:
        parsed_xml = xmltodict.parse(self.body)

        if (
            not parsed_xml["ServerSideEncryptionConfiguration"].get("Rule")
            or not parsed_xml["ServerSideEncryptionConfiguration"]["Rule"].get(
                "ApplyServerSideEncryptionByDefault"
            )
            or not parsed_xml["ServerSideEncryptionConfiguration"]["Rule"][
                "ApplyServerSideEncryptionByDefault"
            ].get("SSEAlgorithm")
        ):
            raise MalformedXML()

        return parsed_xml["ServerSideEncryptionConfiguration"]

    def _ownership_rule_from_body(self) -> dict[str, Any]:
        parsed_xml = xmltodict.parse(self.body)

        if not parsed_xml["OwnershipControls"]["Rule"].get("ObjectOwnership"):
            raise MalformedXML()

        return parsed_xml["OwnershipControls"]["Rule"]["ObjectOwnership"]

    def _logging_from_body(self) -> dict[str, Any]:
        parsed_xml = xmltodict.parse(self.body)

        if not parsed_xml["BucketLoggingStatus"].get("LoggingEnabled"):
            return {}

        if not parsed_xml["BucketLoggingStatus"]["LoggingEnabled"].get("TargetBucket"):
            raise MalformedXML()

        if not parsed_xml["BucketLoggingStatus"]["LoggingEnabled"].get("TargetPrefix"):
            parsed_xml["BucketLoggingStatus"]["LoggingEnabled"]["TargetPrefix"] = ""

        # Get the ACLs:
        if parsed_xml["BucketLoggingStatus"]["LoggingEnabled"].get("TargetGrants"):
            permissions = ["READ", "WRITE", "FULL_CONTROL"]
            if not isinstance(
                parsed_xml["BucketLoggingStatus"]["LoggingEnabled"]["TargetGrants"][
                    "Grant"
                ],
                list,
            ):
                target_grants = self._get_grants_from_xml(
                    [
                        parsed_xml["BucketLoggingStatus"]["LoggingEnabled"][
                            "TargetGrants"
                        ]["Grant"]
                    ],
                    MalformedXML,
                    permissions,
                )
            else:
                target_grants = self._get_grants_from_xml(
                    parsed_xml["BucketLoggingStatus"]["LoggingEnabled"]["TargetGrants"][
                        "Grant"
                    ],
                    MalformedXML,
                    permissions,
                )

            parsed_xml["BucketLoggingStatus"]["LoggingEnabled"]["TargetGrants"] = (
                target_grants
            )

        return parsed_xml["BucketLoggingStatus"]["LoggingEnabled"]

    def _notification_config_from_body(self) -> dict[str, Any]:
        parsed_xml = xmltodict.parse(self.body)

        if not len(parsed_xml["NotificationConfiguration"]):
            return {}

        # The types of notifications, and their required fields (apparently lambda is categorized by the API as
        # "CloudFunction"):
        notification_fields = [
            ("Topic", "sns"),
            ("Queue", "sqs"),
            ("CloudFunction", "lambda"),
            ("EventBridge", "events"),
        ]

        found_notifications = (
            0  # Tripwire -- if this is not ever set, then there were no notifications
        )
        for name, arn_string in notification_fields:
            # EventBridgeConfiguration is passed as an empty dict.
            if name == "EventBridge":
                events_field = f"{name}Configuration"
                if events_field in parsed_xml["NotificationConfiguration"]:
                    parsed_xml["NotificationConfiguration"][events_field] = {}
                    found_notifications += 1
            # 1st verify that the proper notification configuration has been passed in (with an ARN that is close
            # to being correct -- nothing too complex in the ARN logic):
            the_notification = parsed_xml["NotificationConfiguration"].get(
                f"{name}Configuration"
            )
            if the_notification:
                found_notifications += 1
                if not isinstance(the_notification, list):
                    the_notification = parsed_xml["NotificationConfiguration"][
                        f"{name}Configuration"
                    ] = [the_notification]

                for n in the_notification:
                    if not any(
                        n[name].startswith(f"arn:{p}:{arn_string}:")
                        for p in PARTITION_NAMES
                    ):
                        raise InvalidNotificationARN()

                    # 2nd, verify that the Events list is correct:
                    assert n["Event"]
                    if not isinstance(n["Event"], list):
                        n["Event"] = [n["Event"]]

                    # Parse out the filters:
                    if n.get("Filter"):
                        # Error if S3Key is blank:
                        if not n["Filter"]["S3Key"]:
                            raise KeyError()

                        if not isinstance(n["Filter"]["S3Key"]["FilterRule"], list):
                            n["Filter"]["S3Key"]["FilterRule"] = [
                                n["Filter"]["S3Key"]["FilterRule"]
                            ]

                        for filter_rule in n["Filter"]["S3Key"]["FilterRule"]:
                            assert filter_rule["Name"] in ["suffix", "prefix"]
                            assert filter_rule["Value"]

        if not found_notifications:
            return {}

        return parsed_xml["NotificationConfiguration"]

    def _accelerate_config_from_body(self) -> str:
        parsed_xml = xmltodict.parse(self.body)
        config = parsed_xml["AccelerateConfiguration"]
        return config["Status"]

    def _replication_config_from_xml(self, xml: str) -> dict[str, Any]:
        parsed_xml = xmltodict.parse(xml, dict_constructor=dict)
        config = parsed_xml["ReplicationConfiguration"]
        return config

    def _inventory_config_from_body(self) -> dict[str, Any]:
        parsed_xml = xmltodict.parse(self.body)
        config = parsed_xml["InventoryConfiguration"]
        return config

    def _key_response_delete(
        self, bucket_name: str, query: dict[str, Any], key_name: str
    ) -> TYPE_RESPONSE:
        self._set_action("KEY", "DELETE", query)
        self._authenticate_and_authorize_s3_action(
            bucket_name=bucket_name, key_name=key_name
        )

        if query.get("uploadId"):
            return self.abort_multipart_upload()
        if "tagging" in query:
            return self.delete_object_tagging()
        return self.delete_object()

    def delete_object(self) -> TYPE_RESPONSE:
        bypass = self.headers.get("X-Amz-Bypass-Governance-Retention")
        key_name = self.parse_key_name()
        version_id = self._get_param("versionId")
        if_match = self.headers.get("If-Match")

        if if_match:
            if not (obj := self.backend.get_object(self.bucket_name, key_name)):
                raise MissingKey
            # Check if the ETags are the same. S3 doesn't seem to care about quotes, so we shouldn't either
            elif if_match.replace('"', "") != obj.etag.replace('"', ""):
                raise PreconditionFailed("If-Match")

        _, response_meta = self.backend.delete_object(
            self.bucket_name, key_name, version_id=version_id, bypass=bypass
        )
        response_headers = {}
        for k in response_meta:
            response_headers[f"x-amz-{k}"] = response_meta[k]
        return 204, response_headers, ""

    def delete_object_tagging(self) -> TYPE_RESPONSE:
        key_name = self.parse_key_name()
        version_id = self._get_param("versionId")
        self.backend.delete_object_tagging(
            self.bucket_name, key_name, version_id=version_id
        )
        self.data["Action"] = "DeleteObjectTagging"
        return self.serialized(ActionResult({"VersionId": version_id}))

    def _complete_multipart_body(self, body: bytes) -> Iterator[tuple[int, str]]:
        ps = minidom.parseString(body).getElementsByTagName("Part")
        prev = 0
        for p in ps:
            pn = int(p.getElementsByTagName("PartNumber")[0].firstChild.wholeText)  # type: ignore[union-attr]
            if pn <= prev:
                raise InvalidPartOrder()
            yield pn, p.getElementsByTagName("ETag")[0].firstChild.wholeText  # type: ignore[union-attr]

    def _key_response_post(
        self,
        request: Any,
        body: bytes,
        bucket_name: str,
        query: dict[str, Any],
        key_name: str,
    ) -> TYPE_RESPONSE:
        self._set_action("KEY", "POST", query)
        self._authenticate_and_authorize_s3_action(
            bucket_name=bucket_name, key_name=key_name
        )

        encryption = request.headers.get("x-amz-server-side-encryption")
        kms_key_id = request.headers.get("x-amz-server-side-encryption-aws-kms-key-id")

        if body == b"" and "uploads" in query:
            response_headers = {}
            metadata = metadata_from_headers(request.headers)
            tagging = self._tagging_from_headers(request.headers)
            storage_type = request.headers.get("x-amz-storage-class", "STANDARD")
            acl = self._acl_from_headers(request.headers)

            multipart = self.backend.create_multipart_upload(
                bucket_name,
                key_name,
                metadata,
                storage_type,
                tagging,
                acl,
                encryption,
                kms_key_id,
            )
            if encryption:
                response_headers["x-amz-server-side-encryption"] = encryption
            if kms_key_id:
                response_headers["x-amz-server-side-encryption-aws-kms-key-id"] = (
                    kms_key_id
                )

            self.data["Action"] = "CreateMultipartUpload"
            status, resp_headers, resp_body = self.serialized(
                ActionResult(
                    {
                        "Bucket": bucket_name,
                        "Key": key_name,
                        "UploadId": multipart.id,
                        "ChecksumAlgorithm": multipart.metadata.get(
                            "x-amz-checksum-algorithm"
                        ),
                        "ChecksumType": multipart.metadata.get("x-amz-checksum-type"),
                    }
                )
            )
            resp_headers.update(response_headers)
            return status, resp_headers, resp_body

        if query.get("uploadId"):
            existing = self.backend.get_object(self.bucket_name, key_name)

            multipart_id = query["uploadId"][0]

            if (
                existing is not None
                and existing.multipart
                and existing.multipart.id == multipart_id
            ):
                # Operation is idempotent
                key: FakeKey | None = existing
            else:
                if self.headers.get("If-None-Match") == "*" and existing is not None:
                    raise PreconditionFailed("If-None-Match")

                key = self.backend.complete_multipart_upload(
                    bucket_name, multipart_id, self._complete_multipart_body(body)
                )
            if key is None:
                return 400, {}, ""

            headers: dict[str, Any] = {}

            if key.version_id:
                headers["x-amz-version-id"] = key.version_id

            if key.encryption:
                headers["x-amz-server-side-encryption"] = key.encryption

            if key.kms_key_id:
                headers["x-amz-server-side-encryption-aws-kms-key-id"] = key.kms_key_id

            self.data["Action"] = "CompleteMultipartUpload"
            status, resp_headers, resp_body = self.serialized(
                ActionResult(
                    {
                        "Location": f"http://{bucket_name}.s3.amazonaws.com/{key.name}",
                        "Bucket": bucket_name,
                        "Key": key.name,
                        "ETag": key.etag,
                    }
                )
            )
            resp_headers.update(headers)
            return status, resp_headers, resp_body

        elif "restore" in query:
            params = xmltodict.parse(body)["RestoreRequest"]
            previously_restored = self.backend.restore_object(
                bucket_name,
                key_name,
                params.get("Days", None),
                params.get("Type", None),
            )
            status_code = 200 if previously_restored else 202
            return status_code, {}, ""
        elif "select" in query:
            request = xmltodict.parse(body)["SelectObjectContentRequest"]
            select_query = request["Expression"]
            input_details = request["InputSerialization"]
            output_details = request["OutputSerialization"]
            results, bytes_scanned = self.backend.select_object_content(
                bucket_name, key_name, select_query, input_details, output_details
            )
            bytes_returned = sum(len(line) for line in results)
            event_stream = [
                {"Records": {"Payload": b"".join(results)}},
                {
                    "Stats": {
                        "Details": {
                            "BytesScanned": bytes_scanned,
                            "BytesProcessed": bytes_scanned,
                            "BytesReturned": bytes_returned,
                        }
                    }
                },
                {"End": {}},
            ]
            self.data["Action"] = "SelectObjectContent"
            return self.serialized(ActionResult({"Payload": event_stream}))

        else:
            raise NotImplementedError(
                "Method POST had only been implemented for multipart uploads and restore operations, so far"
            )

    def _invalid_headers(self, url: str, headers: dict[str, str]) -> bool:
        """
        Verify whether the provided metadata in the URL is also present in the headers
        :param url: .../file.txt&content-type=app%2Fjson&Signature=..
        :param headers: Content-Type=app/json
        :return: True or False
        """
        metadata_to_check = {
            "content-disposition": "Content-Disposition",
            "content-encoding": "Content-Encoding",
            "content-language": "Content-Language",
            "content-length": "Content-Length",
            "content-md5": "Content-MD5",
            "content-type": "Content-Type",
        }
        for url_key, header_key in metadata_to_check.items():
            metadata_in_url = re.search(url_key + "=(.+?)(&.+$|$)", url)
            if metadata_in_url:
                url_value = unquote(metadata_in_url.group(1))
                if header_key not in headers or (url_value != headers[header_key]):
                    return True
        return False


S3ResponseInstance = S3Response()
