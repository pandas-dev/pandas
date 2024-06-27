import io
import re
import urllib.parse
from typing import Any, Dict, Iterator, List, Optional, Tuple, Type, Union
from urllib.parse import parse_qs, unquote, urlencode, urlparse, urlunparse
from xml.dom import minidom

import xmltodict

from moto import settings
from moto.core.common_types import TYPE_RESPONSE
from moto.core.responses import BaseResponse
from moto.core.utils import (
    ALT_DOMAIN_SUFFIXES,
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

from .exceptions import (
    AccessForbidden,
    BucketAccessDeniedError,
    BucketAlreadyExists,
    BucketMustHaveLockeEnabled,
    DuplicateTagKeys,
    HeadOnDeleteMarker,
    IllegalLocationConstraintException,
    InvalidContentMD5,
    InvalidContinuationToken,
    InvalidMaxPartArgument,
    InvalidMaxPartNumberArgument,
    InvalidNotificationARN,
    InvalidObjectState,
    InvalidPartOrder,
    InvalidRange,
    LockNotEnabled,
    MalformedACLError,
    MalformedXML,
    MissingBucket,
    MissingKey,
    MissingVersion,
    NoSystemTags,
    NotAnIntegerException,
    ObjectNotInActiveTierError,
    PreconditionFailed,
    S3AclAndGrantError,
    S3ClientError,
)
from .models import (
    FakeAcl,
    FakeBucket,
    FakeGrant,
    FakeGrantee,
    FakeKey,
    S3Backend,
    get_canned_acl,
    s3_backends,
)
from .select_object_content import serialize_select
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
        "HEAD": {"DEFAULT": "HeadObject"},
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

    def all_buckets(self) -> str:
        self.data["Action"] = "ListAllMyBuckets"
        self._authenticate_and_authorize_s3_action()

        # No bucket specified. Listing all buckets
        all_buckets = self.backend.list_buckets()
        template = self.response_template(S3_ALL_BUCKETS)
        return template.render(buckets=all_buckets)

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
            and any([host in endpoint for endpoint in custom_endpoints])
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
            or re.match(r"s3[\.\-]([^.]*)\.amazonaws\.com", host)
            or any(
                re.match(r"s3[\.\-]([^.]*)\." + suffix, host)
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

    def parse_key_name(self, request: Any, url: str) -> str:
        if self.subdomain_based_buckets(request):
            return parse_key_name(url)
        else:
            return bucketpath_parse_key_name(url)

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
        self.setup_class(request, full_url, headers, use_raw_body=True)
        bucket_name = self.parse_bucket_name_from_url(request, full_url)
        self.backend.log_incoming_request(request, bucket_name)
        try:
            response = self._bucket_response(request, full_url)
        except S3ClientError as s3error:
            response = s3error.code, {}, s3error.description

        return self._send_response(response)

    @staticmethod
    def _send_response(response: Any) -> TYPE_RESPONSE:  # type: ignore
        if isinstance(response, str):
            return 200, {}, response.encode("utf-8")
        else:
            status_code, headers, response_content = response
            if not isinstance(response_content, bytes):
                response_content = response_content.encode("utf-8")

            return status_code, headers, response_content

    def _bucket_response(
        self, request: Any, full_url: str
    ) -> Union[str, TYPE_RESPONSE]:
        querystring = self._get_querystring(request, full_url)
        method = request.method
        region_name = parse_region_from_url(full_url, use_default_region=False)
        if region_name is None:
            region_name = extract_region_from_aws_authorization(
                request.headers.get("Authorization", "")
            )
            region_name = region_name or DEFAULT_REGION_NAME

        bucket_name = self.parse_bucket_name_from_url(request, full_url)
        if not bucket_name:
            # If no bucket specified, list all buckets
            return self.all_buckets()

        self.data["BucketName"] = bucket_name

        if method == "HEAD":
            return self._bucket_response_head(bucket_name, querystring)
        elif method == "GET":
            return self._bucket_response_get(bucket_name, querystring)
        elif method == "PUT":
            return self._bucket_response_put(
                request, region_name, bucket_name, querystring
            )
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

    def _get_querystring(self, request: Any, full_url: str) -> Dict[str, Any]:
        # Flask's Request has the querystring already parsed
        # In ServerMode, we can use this, instead of manually parsing this
        if hasattr(request, "args"):
            query_dict = dict()
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

    def _bucket_response_head(
        self, bucket_name: str, querystring: Dict[str, Any]
    ) -> TYPE_RESPONSE:
        self._set_action("BUCKET", "HEAD", querystring)
        self._authenticate_and_authorize_s3_action(bucket_name=bucket_name)

        try:
            bucket = self.backend.head_bucket(bucket_name)
        except MissingBucket:
            # Unless we do this, boto3 does not raise ClientError on
            # HEAD (which the real API responds with), and instead
            # raises NoSuchBucket, leading to inconsistency in
            # error response between real and mocked responses.
            return 404, {}, ""
        return 200, {"x-amz-bucket-region": bucket.region_name}, ""

    def _set_cors_headers_options(
        self, headers: Dict[str, str], bucket: FakeBucket
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

        def _to_string(header: Union[List[str], str]) -> str:
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
        self, headers: Dict[str, str], bucket_name: str
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

    def _get_cors_headers_other(
        self, headers: Dict[str, str], bucket_name: str
    ) -> Dict[str, Any]:
        """
        Returns a dictionary with the appropriate CORS headers
        Should be used for non-OPTIONS requests only
        Applicable if the 'Origin' header matches one of a CORS-rules - returns an empty dictionary otherwise
        """
        response_headers: Dict[str, Any] = dict()
        try:
            origin = headers.get("Origin")
            if not origin:
                return response_headers
            bucket = self.backend.get_bucket(bucket_name)

            def _to_string(header: Union[List[str], str]) -> str:
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
        self, bucket_name: str, querystring: Dict[str, Any]
    ) -> Union[str, TYPE_RESPONSE]:
        self._set_action("BUCKET", "GET", querystring)
        self._authenticate_and_authorize_s3_action(bucket_name=bucket_name)

        if "object-lock" in querystring:
            (
                lock_enabled,
                mode,
                days,
                years,
            ) = self.backend.get_object_lock_configuration(bucket_name)
            template = self.response_template(S3_BUCKET_LOCK_CONFIGURATION)

            return template.render(
                lock_enabled=lock_enabled, mode=mode, days=days, years=years
            )

        if "uploads" in querystring:
            for unsup in ("delimiter", "max-uploads"):
                if unsup in querystring:
                    raise NotImplementedError(
                        f"Listing multipart uploads with {unsup} has not been implemented yet."
                    )
            multiparts = list(self.backend.get_all_multiparts(bucket_name).values())
            if "prefix" in querystring:
                prefix = querystring.get("prefix", [None])[0]
                multiparts = [
                    upload
                    for upload in multiparts
                    if upload.key_name.startswith(prefix)
                ]
            template = self.response_template(S3_ALL_MULTIPARTS)
            return template.render(
                bucket_name=bucket_name,
                uploads=multiparts,
                account_id=self.current_account,
            )
        elif "location" in querystring:
            location: Optional[str] = self.backend.get_bucket_location(bucket_name)
            template = self.response_template(S3_BUCKET_LOCATION)

            # us-east-1 is different - returns a None location
            if location == DEFAULT_REGION_NAME:
                location = None

            return template.render(location=location)
        elif "lifecycle" in querystring:
            rules = self.backend.get_bucket_lifecycle(bucket_name)
            if not rules:
                template = self.response_template(S3_NO_LIFECYCLE)
                return 404, {}, template.render(bucket_name=bucket_name)
            template = self.response_template(S3_BUCKET_LIFECYCLE_CONFIGURATION)
            return template.render(rules=rules)
        elif "versioning" in querystring:
            versioning = self.backend.get_bucket_versioning(bucket_name)
            template = self.response_template(S3_BUCKET_GET_VERSIONING)
            return template.render(status=versioning)
        elif "policy" in querystring:
            policy = self.backend.get_bucket_policy(bucket_name)
            if not policy:
                template = self.response_template(S3_NO_POLICY)
                return 404, {}, template.render(bucket_name=bucket_name)
            return 200, {}, policy
        elif "website" in querystring:
            website_configuration = self.backend.get_bucket_website_configuration(
                bucket_name
            )
            if not website_configuration:
                template = self.response_template(S3_NO_BUCKET_WEBSITE_CONFIG)
                return 404, {}, template.render(bucket_name=bucket_name)
            return 200, {}, website_configuration  # type: ignore
        elif "acl" in querystring:
            acl = self.backend.get_bucket_acl(bucket_name)
            template = self.response_template(S3_OBJECT_ACL_RESPONSE)
            return template.render(acl=acl)
        elif "tagging" in querystring:
            tags = self.backend.get_bucket_tagging(bucket_name)["Tags"]
            # "Special Error" if no tags:
            if len(tags) == 0:
                template = self.response_template(S3_NO_BUCKET_TAGGING)
                return 404, {}, template.render(bucket_name=bucket_name)
            template = self.response_template(S3_OBJECT_TAGGING_RESPONSE)
            return template.render(tags=tags)
        elif "logging" in querystring:
            logging = self.backend.get_bucket_logging(bucket_name)
            if not logging:
                template = self.response_template(S3_NO_LOGGING_CONFIG)
                return 200, {}, template.render()
            template = self.response_template(S3_LOGGING_CONFIG)
            return 200, {}, template.render(logging=logging)
        elif "cors" in querystring:
            cors = self.backend.get_bucket_cors(bucket_name)
            if len(cors) == 0:
                template = self.response_template(S3_NO_CORS_CONFIG)
                return 404, {}, template.render(bucket_name=bucket_name)
            template = self.response_template(S3_BUCKET_CORS_RESPONSE)
            return template.render(cors=cors)
        elif "notification" in querystring:
            notification_configuration = (
                self.backend.get_bucket_notification_configuration(bucket_name)
            )
            if not notification_configuration:
                return 200, {}, ""
            template = self.response_template(S3_GET_BUCKET_NOTIFICATION_CONFIG)
            return template.render(config=notification_configuration)
        elif "accelerate" in querystring:
            bucket = self.backend.get_bucket(bucket_name)
            if bucket.accelerate_configuration is None:
                template = self.response_template(S3_BUCKET_ACCELERATE_NOT_SET)
                return 200, {}, template.render()
            template = self.response_template(S3_BUCKET_ACCELERATE)
            return template.render(bucket=bucket)
        elif "publicAccessBlock" in querystring:
            public_block_config = self.backend.get_public_access_block(bucket_name)
            template = self.response_template(S3_PUBLIC_ACCESS_BLOCK_CONFIGURATION)
            return template.render(public_block_config=public_block_config)

        elif "versions" in querystring:
            delimiter = querystring.get("delimiter", [None])[0]
            key_marker = querystring.get("key-marker", [None])[0]
            max_keys = int(
                querystring.get("max-keys", [settings.get_s3_default_max_keys()])[0]
            )
            prefix = querystring.get("prefix", [""])[0]
            version_id_marker = querystring.get("version-id-marker", [None])[0]

            bucket = self.backend.get_bucket(bucket_name)
            (
                versions,
                common_prefixes,
                delete_markers,
                next_key_marker,
                next_version_id_marker,
            ) = self.backend.list_object_versions(
                bucket_name,
                delimiter=delimiter,
                key_marker=key_marker,
                max_keys=max_keys,
                prefix=prefix,
                version_id_marker=version_id_marker,
            )
            key_list = versions

            is_truncated = False
            if next_key_marker is not None:
                is_truncated = True

            template = self.response_template(S3_BUCKET_GET_VERSIONS)

            return (
                200,
                {},
                template.render(
                    common_prefixes=common_prefixes,
                    key_list=key_list,
                    delete_marker_list=delete_markers,
                    bucket=bucket,
                    prefix=prefix,
                    max_keys=max_keys,
                    delimiter=delimiter,
                    key_marker=key_marker,
                    version_id_marker=version_id_marker,
                    is_truncated=is_truncated,
                    next_key_marker=next_key_marker,
                    next_version_id_marker=next_version_id_marker,
                ),
            )
        elif "encryption" in querystring:
            encryption = self.backend.get_bucket_encryption(bucket_name)
            if not encryption:
                template = self.response_template(S3_NO_ENCRYPTION)
                return 404, {}, template.render(bucket_name=bucket_name)
            template = self.response_template(S3_ENCRYPTION_CONFIG)
            return 200, {}, template.render(encryption=encryption)
        elif querystring.get("list-type", [None])[0] == "2":
            return 200, {}, self._handle_list_objects_v2(bucket_name, querystring)
        elif "replication" in querystring:
            replication = self.backend.get_bucket_replication(bucket_name)
            if not replication:
                template = self.response_template(S3_NO_REPLICATION)
                return 404, {}, template.render(bucket_name=bucket_name)
            template = self.response_template(S3_REPLICATION_CONFIG)
            return 200, {}, template.render(replication=replication)
        elif "ownershipControls" in querystring:
            ownership_rule = self.backend.get_bucket_ownership_controls(bucket_name)
            if not ownership_rule:
                template = self.response_template(S3_ERROR_BUCKET_ONWERSHIP_NOT_FOUND)
                return 404, {}, template.render(bucket_name=bucket_name)
            template = self.response_template(S3_BUCKET_GET_OWNERSHIP_RULE)
            return 200, {}, template.render(ownership_rule=ownership_rule)

        bucket = self.backend.get_bucket(bucket_name)
        prefix = querystring.get("prefix", [None])[0]
        if prefix and isinstance(prefix, bytes):
            prefix = prefix.decode("utf-8")
        delimiter = querystring.get("delimiter", [None])[0]
        max_keys = int(
            querystring.get("max-keys", [settings.get_s3_default_max_keys()])[0]
        )
        marker = querystring.get("marker", [None])[0]
        encoding_type = querystring.get("encoding-type", [None])[0]

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

        template = self.response_template(S3_BUCKET_GET_RESPONSE)
        return (
            200,
            {},
            template.render(
                bucket=bucket,
                prefix=prefix,
                delimiter=delimiter,
                result_keys=result_keys,
                result_folders=result_folders,
                is_truncated=is_truncated,
                next_marker=next_marker,
                max_keys=max_keys,
                encoding_type=encoding_type,
            ),
        )

    def _set_action(
        self, action_resource_type: str, method: str, querystring: Dict[str, Any]
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

    def _handle_list_objects_v2(
        self, bucket_name: str, querystring: Dict[str, Any]
    ) -> str:
        template = self.response_template(S3_BUCKET_GET_RESPONSE_V2)
        bucket = self.backend.get_bucket(bucket_name)

        continuation_token = querystring.get("continuation-token", [None])[0]
        if continuation_token is not None and continuation_token == "":
            raise InvalidContinuationToken()

        prefix = querystring.get("prefix", [None])[0]
        if prefix and isinstance(prefix, bytes):
            prefix = prefix.decode("utf-8")
        delimiter = querystring.get("delimiter", [None])[0]

        fetch_owner = querystring.get("fetch-owner", [False])[0]
        max_keys = int(
            querystring.get("max-keys", [settings.get_s3_default_max_keys()])[0]
        )
        start_after = querystring.get("start-after", [None])[0]
        encoding_type = querystring.get("encoding-type", [None])[0]

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
            result_folders = list(
                map(lambda folder: urllib.parse.quote(folder), result_folders)
            )

        return template.render(
            bucket=bucket,
            prefix=prefix or "",
            delimiter=delimiter,
            key_count=key_count,
            result_keys=result_keys,
            result_folders=result_folders,
            fetch_owner=fetch_owner,
            max_keys=max_keys,
            is_truncated=is_truncated,
            next_continuation_token=next_continuation_token,
            start_after=None if continuation_token else start_after,
            encoding_type=encoding_type,
        )

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

    def _body_contains_location_constraint(self, body: bytes) -> bool:
        if body:
            try:
                xmltodict.parse(body)["CreateBucketConfiguration"]["LocationConstraint"]
                return True
            except KeyError:
                pass
        return False

    def _create_bucket_configuration_is_empty(self, body: bytes) -> bool:
        if body:
            try:
                create_bucket_configuration = xmltodict.parse(body)[
                    "CreateBucketConfiguration"
                ]
                del create_bucket_configuration["@xmlns"]
                if len(create_bucket_configuration) == 0:
                    return True
            except KeyError:
                pass
        return False

    def _parse_pab_config(self) -> Dict[str, Any]:
        parsed_xml = xmltodict.parse(self.body)
        parsed_xml["PublicAccessBlockConfiguration"].pop("@xmlns", None)

        return parsed_xml

    def _bucket_response_put(
        self,
        request: Any,
        region_name: str,
        bucket_name: str,
        querystring: Dict[str, Any],
    ) -> Union[str, TYPE_RESPONSE]:
        if querystring and not request.headers.get("Content-Length"):
            return 411, {}, "Content-Length required"

        self._set_action("BUCKET", "PUT", querystring)
        self._authenticate_and_authorize_s3_action(bucket_name=bucket_name)

        if "object-lock" in querystring:
            config = self._lock_config_from_body()

            if not self.backend.get_bucket(bucket_name).object_lock_enabled:
                raise BucketMustHaveLockeEnabled

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
                template = self.response_template(S3_BUCKET_VERSIONING)
                return template.render(bucket_versioning_status=ver.group(1))
            else:
                return 404, {}, ""
        elif "lifecycle" in querystring:
            rules = xmltodict.parse(self.body)["LifecycleConfiguration"]["Rule"]
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
            self.backend.set_bucket_website_configuration(bucket_name, self.body)
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
            except Exception as e:
                raise e
        elif "accelerate" in querystring:
            try:
                accelerate_status = self._accelerate_config_from_body()
                self.backend.put_bucket_accelerate_configuration(
                    bucket_name, accelerate_status
                )
                return ""
            except KeyError:
                raise MalformedXML()
            except Exception as e:
                raise e

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
            except Exception as e:
                raise e
        elif "replication" in querystring:
            bucket = self.backend.get_bucket(bucket_name)
            if not bucket.is_versioned:
                template = self.response_template(S3_NO_VERSIONING_ENABLED)
                return 400, {}, template.render(bucket_name=bucket_name)
            replication_config = self._replication_config_from_xml(self.body)
            self.backend.put_bucket_replication(bucket_name, replication_config)
            return ""
        elif "ownershipControls" in querystring:
            ownership_rule = self._ownership_rule_from_body()
            self.backend.put_bucket_ownership_controls(
                bucket_name, ownership=ownership_rule
            )
            return ""

        else:
            # us-east-1, the default AWS region behaves a bit differently
            # - you should not use it as a location constraint --> it fails
            # - querying the location constraint returns None
            # - LocationConstraint has to be specified if outside us-east-1
            if (
                region_name != DEFAULT_REGION_NAME
                and not self._body_contains_location_constraint(self.body)
            ):
                raise IllegalLocationConstraintException()
            if self.body:
                if self._create_bucket_configuration_is_empty(self.body):
                    raise MalformedXML()
                try:
                    forced_region = xmltodict.parse(self.body)[
                        "CreateBucketConfiguration"
                    ]["LocationConstraint"]

                    if forced_region == DEFAULT_REGION_NAME:
                        raise S3ClientError(
                            "InvalidLocationConstraint",
                            "The specified location-constraint is not valid",
                        )
                    else:
                        region_name = forced_region
                except KeyError:
                    pass

            try:
                new_bucket = self.backend.create_bucket(bucket_name, region_name)
            except BucketAlreadyExists:
                new_bucket = self.backend.get_bucket(bucket_name)
                if new_bucket.account_id == self.get_current_account():
                    # special cases when the bucket belongs to self
                    if (
                        new_bucket.region_name == DEFAULT_REGION_NAME
                        and region_name == DEFAULT_REGION_NAME
                    ):
                        # us-east-1 has different behavior - creating a bucket there is an idempotent operation
                        pass
                    else:
                        template = self.response_template(S3_DUPLICATE_BUCKET_ERROR)
                        return 409, {}, template.render(bucket_name=bucket_name)
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

            template = self.response_template(S3_BUCKET_CREATE_RESPONSE)
            return 200, {}, template.render(bucket=new_bucket)

    def _bucket_response_delete(
        self, bucket_name: str, querystring: Dict[str, Any]
    ) -> TYPE_RESPONSE:
        self._set_action("BUCKET", "DELETE", querystring)
        self._authenticate_and_authorize_s3_action(bucket_name=bucket_name)

        if "policy" in querystring:
            self.backend.delete_bucket_policy(bucket_name)
            return 204, {}, ""
        elif "tagging" in querystring:
            self.backend.delete_bucket_tagging(bucket_name)
            return 204, {}, ""
        elif "website" in querystring:
            self.backend.delete_bucket_website(bucket_name)
            return 204, {}, ""
        elif "cors" in querystring:
            self.backend.delete_bucket_cors(bucket_name)
            return 204, {}, ""
        elif "lifecycle" in querystring:
            self.backend.delete_bucket_lifecycle(bucket_name)
            return 204, {}, ""
        elif "publicAccessBlock" in querystring:
            self.backend.delete_public_access_block(bucket_name)
            return 204, {}, ""
        elif "encryption" in querystring:
            self.backend.delete_bucket_encryption(bucket_name)
            return 204, {}, ""
        elif "replication" in querystring:
            self.backend.delete_bucket_replication(bucket_name)
            return 204, {}, ""
        elif "ownershipControls" in querystring:
            self.backend.delete_bucket_ownership_controls(bucket_name)
            return 204, {}, ""

        removed_bucket = self.backend.delete_bucket(bucket_name)

        if removed_bucket:
            # Bucket exists
            template = self.response_template(S3_DELETE_BUCKET_SUCCESS)
            return 204, {}, template.render(bucket=removed_bucket)
        else:
            # Tried to delete a bucket that still has keys
            template = self.response_template(S3_DELETE_BUCKET_WITH_ITEMS_ERROR)
            return 409, {}, template.render(bucket=removed_bucket)

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
            queryargs: Dict[str, Any] = parse_qs(parts.query)
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
        template = self.response_template(S3_DELETE_KEYS_RESPONSE)
        body_dict = xmltodict.parse(self.body, strip_whitespace=False)

        objects = body_dict["Delete"].get("Object", [])
        if not isinstance(objects, list):
            # We expect a list of objects, but when there is a single <Object> node xmltodict does not
            # return a list.
            objects = [objects]
        if len(objects) == 0:
            raise MalformedXML()

        if authenticated:
            deleted_objects = self.backend.delete_objects(bucket_name, objects)
            errors = []
        else:
            deleted_objects = []
            # [(key_name, errorcode, 'error message'), ..]
            errors = [(o["Key"], "AccessDenied", "Access Denied") for o in objects]

        return (
            200,
            {},
            template.render(deleted=deleted_objects, delete_errors=errors),
        )

    def _handle_range_header(
        self, request: Any, response_headers: Dict[str, Any], response_content: Any
    ) -> TYPE_RESPONSE:
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

        response_headers["content-range"] = f"bytes {begin}-{end}/{length}"
        content = response_content[begin : end + 1]
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
            new_body[pos : pos + chunk_size] = body_io.read(chunk_size)
            pos = pos + chunk_size
            body_io.read(2)  # skip trailing \r\n
            line = body_io.readline()
        return bytes(new_body)

    def _handle_encoded_body(self, body: bytes, content_length: int) -> bytes:
        body_io = io.BytesIO(body)
        # first line should equal '{content_length}\r\n
        body_io.readline()
        # Body contains actual data next
        return body_io.read(content_length)
        # last line should equal
        # amz-checksum-sha256:<..>\r\n

    def key_response(
        self, request: Any, full_url: str, headers: Dict[str, Any]
    ) -> TYPE_RESPONSE:
        # Key and Control are lumped in because splitting out the regex is too much of a pain :/
        self.setup_class(request, full_url, headers, use_raw_body=True)
        bucket_name = self.parse_bucket_name_from_url(request, full_url)
        self.backend.log_incoming_request(request, bucket_name)
        response_headers: Dict[str, Any] = {}

        try:
            response = self._key_response(request, full_url, self.headers)
        except S3ClientError as s3error:
            response = s3error.code, {}, s3error.description

        if isinstance(response, str):
            status_code = 200
            response_content = response
        else:
            status_code, response_headers, response_content = response

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
                return s3error.code, {}, s3error.description
        return status_code, response_headers, response_content

    def _key_response(
        self, request: Any, full_url: str, headers: Dict[str, Any]
    ) -> TYPE_RESPONSE:
        parsed_url = urlparse(full_url)
        url_path = self.get_safe_path()
        query = parse_qs(parsed_url.query, keep_blank_values=True)
        method = request.method

        key_name = self.parse_key_name(request, url_path)
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
        if key:
            resource = f"arn:{bucket.partition}:s3:::{bucket_name}/{key_name}"  # type: ignore[union-attr]

            # Authorization Workflow
            # https://docs.aws.amazon.com/AmazonS3/latest/userguide/access-control-auth-workflow-object-operation.html

            # A bucket can deny all actions, regardless of who makes the request
            from moto.iam.access_control import PermissionResult

            action = f"s3:{method.upper()[0]}{method.lower()[1:]}Object"
            bucket_permissions = bucket.get_permission(action, resource)  # type: ignore
            if bucket_permissions == PermissionResult.DENIED:
                return 403, {}, ""

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

        elif signed_url and not authorized_request:
            # coming in from requests.get(s3.generate_presigned_url())
            if self._invalid_headers(request.url, dict(request.headers)):
                return 403, {}, S3_INVALID_PRESIGNED_PARAMETERS

        body = self.body or b""

        if (
            request.headers.get("x-amz-content-sha256", None)
            == "STREAMING-AWS4-HMAC-SHA256-PAYLOAD"
        ):
            body = self._handle_v4_chunk_signatures(
                body, int(request.headers["x-amz-decoded-content-length"])
            )

        if method == "GET":
            return self._key_response_get(
                bucket_name, query, key_name, headers=request.headers
            )
        elif method == "PUT":
            return self._key_response_put(request, body, bucket_name, query, key_name)
        elif method == "HEAD":
            return self._key_response_head(
                bucket_name, query, key_name, headers=request.headers
            )
        elif method == "DELETE":
            return self._key_response_delete(headers, bucket_name, query, key_name)
        elif method == "POST":
            return self._key_response_post(request, body, bucket_name, query, key_name)
        elif method == "OPTIONS":
            # OPTIONS response doesn't depend on the key_name: always return 200 with CORS headers
            return self._response_options(request.headers, bucket_name)
        else:
            raise NotImplementedError(
                f"Method {method} has not been implemented in the S3 backend yet"
            )

    def _key_response_get(
        self,
        bucket_name: str,
        query: Dict[str, Any],
        key_name: str,
        headers: Dict[str, Any],
    ) -> TYPE_RESPONSE:
        self._set_action("KEY", "GET", query)
        self._authenticate_and_authorize_s3_action(
            bucket_name=bucket_name, key_name=key_name
        )

        response_headers = self._get_cors_headers_other(headers, bucket_name)
        if query.get("uploadId"):
            upload_id = query["uploadId"][0]

            # 0 <= PartNumberMarker <= 2,147,483,647
            part_number_marker = int(query.get("part-number-marker", [0])[0])
            if part_number_marker > 2147483647:
                raise NotAnIntegerException(
                    name="part-number-marker", value=part_number_marker
                )
            if not (0 <= part_number_marker <= 2147483647):
                raise InvalidMaxPartArgument("part-number-marker", 0, 2147483647)

            # 0 <= MaxParts <= 2,147,483,647 (default is 1,000)
            max_parts = int(query.get("max-parts", [1000])[0])
            if max_parts > 2147483647:
                raise NotAnIntegerException(name="max-parts", value=max_parts)
            if not (0 <= max_parts <= 2147483647):
                raise InvalidMaxPartArgument("max-parts", 0, 2147483647)

            parts = self.backend.list_parts(
                bucket_name,
                upload_id,
                part_number_marker=part_number_marker,
                max_parts=max_parts,
            )
            next_part_number_marker = parts[-1].name if parts else 0
            is_truncated = len(parts) != 0 and self.backend.is_truncated(
                bucket_name,
                upload_id,
                next_part_number_marker,  # type: ignore
            )

            template = self.response_template(S3_MULTIPART_LIST_RESPONSE)
            return (
                200,
                response_headers,
                template.render(
                    bucket_name=bucket_name,
                    key_name=key_name,
                    upload_id=upload_id,
                    is_truncated=str(is_truncated).lower(),
                    max_parts=max_parts,
                    next_part_number_marker=next_part_number_marker,
                    parts=parts,
                    part_number_marker=part_number_marker,
                ),
            )
        version_id = query.get("versionId", [None])[0]
        if_modified_since = headers.get("If-Modified-Since", None)
        if_match = headers.get("If-Match", None)
        if_none_match = headers.get("If-None-Match", None)
        if_unmodified_since = headers.get("If-Unmodified-Since", None)

        key = self.backend.get_object(bucket_name, key_name, version_id=version_id)
        if key is None and version_id is None:
            raise MissingKey(key=key_name)
        elif key is None:
            raise MissingVersion()

        if key.version_id != "null":
            response_headers["x-amz-version-id"] = key.version_id

        if key.storage_class in ARCHIVE_STORAGE_CLASSES:
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

        if if_modified_since:
            if_modified_since = str_to_rfc_1123_datetime(if_modified_since)
            if key.last_modified.replace(microsecond=0) <= if_modified_since:
                return 304, response_headers, "Not Modified"
        if if_none_match and key.etag in [if_none_match, f'"{if_none_match}"']:
            return 304, response_headers, "Not Modified"

        if "acl" in query:
            acl = self.backend.get_object_acl(key)
            template = self.response_template(S3_OBJECT_ACL_RESPONSE)
            return 200, response_headers, template.render(acl=acl)
        if "tagging" in query:
            tags = self.backend.get_object_tagging(key)["Tags"]
            template = self.response_template(S3_OBJECT_TAGGING_RESPONSE)
            return 200, response_headers, template.render(tags=tags)
        if "legal-hold" in query:
            legal_hold = self.backend.get_object_legal_hold(key)
            template = self.response_template(S3_OBJECT_LEGAL_HOLD)
            return 200, response_headers, template.render(legal_hold=legal_hold)
        if "attributes" in query:
            attributes_to_get = headers.get("x-amz-object-attributes", "").split(",")
            response_keys = self.backend.get_object_attributes(key, attributes_to_get)

            response_headers["Last-Modified"] = key.last_modified_ISO8601

            template = self.response_template(S3_OBJECT_ATTRIBUTES_RESPONSE)
            return 200, response_headers, template.render(**response_keys)

        response_headers.update(key.metadata)
        response_headers.update(key.response_dict)
        response_headers.update({"Accept-Ranges": "bytes"})
        return 200, response_headers, key.value

    def _key_response_put(
        self,
        request: Any,
        body: bytes,
        bucket_name: str,
        query: Dict[str, Any],
        key_name: str,
    ) -> TYPE_RESPONSE:
        self._set_action("KEY", "PUT", query)
        self._authenticate_and_authorize_s3_action(
            bucket_name=bucket_name, key_name=key_name
        )

        response_headers = self._get_cors_headers_other(request.headers, bucket_name)
        if query.get("uploadId") and query.get("partNumber"):
            upload_id = query["uploadId"][0]
            part_number = int(query["partNumber"][0])
            if "x-amz-copy-source" in request.headers:
                copy_source = request.headers.get("x-amz-copy-source")
                if isinstance(copy_source, bytes):
                    copy_source = copy_source.decode("utf-8")
                copy_source_parsed = urlparse(copy_source)
                src_bucket, src_key = (
                    unquote(copy_source_parsed.path).lstrip("/").split("/", 1)
                )
                src_version_id = parse_qs(copy_source_parsed.query).get(
                    "versionId", [None]
                )[0]
                src_range = request.headers.get("x-amz-copy-source-range", "").split(
                    "bytes="
                )[-1]

                try:
                    start_byte, end_byte = src_range.split("-")
                    start_byte, end_byte = int(start_byte), int(end_byte)
                except ValueError:
                    start_byte, end_byte = None, None

                if self.backend.get_object(
                    src_bucket, src_key, version_id=src_version_id
                ):
                    key = self.backend.upload_part_copy(
                        bucket_name,
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

                template = self.response_template(S3_MULTIPART_UPLOAD_RESPONSE)
                response = template.render(part=key)
            else:
                if part_number > 10000:
                    raise InvalidMaxPartNumberArgument(part_number)
                key = self.backend.upload_part(
                    bucket_name, upload_id, part_number, body
                )
                response = ""
            response_headers.update(key.response_dict)
            response_headers["content-length"] = str(len(response))
            return 200, response_headers, response

        storage_class = request.headers.get("x-amz-storage-class", "STANDARD")
        encryption = request.headers.get("x-amz-server-side-encryption", None)
        kms_key_id = request.headers.get(
            "x-amz-server-side-encryption-aws-kms-key-id", None
        )

        checksum_algorithm = request.headers.get("x-amz-sdk-checksum-algorithm", "")
        checksum_header = f"x-amz-checksum-{checksum_algorithm.lower()}"
        checksum_value = request.headers.get(checksum_header)
        if not checksum_value and checksum_algorithm:
            # Extract the checksum-value from the body first
            search = re.search(rb"x-amz-checksum-\w+:(.+={1,2})", body)
            checksum_value = search.group(1) if search else None

        if checksum_value:
            # TODO: AWS computes the provided value and verifies it's the same
            # Afterwards, it should be returned in every subsequent call
            if isinstance(checksum_value, bytes):
                checksum_value = checksum_value.decode("utf-8")
            response_headers.update({checksum_header: checksum_value})
        elif checksum_algorithm:
            # If the value is not provided, we compute it and only return it as part of this request
            checksum_value = compute_checksum(body, algorithm=checksum_algorithm)
            response_headers.update({checksum_header: checksum_value})

        # Extract the actual data from the body second
        if (
            request.headers.get("x-amz-content-sha256", None)
            == "STREAMING-UNSIGNED-PAYLOAD-TRAILER"
        ):
            body = self._handle_encoded_body(
                body, int(request.headers["x-amz-decoded-content-length"])
            )

        bucket_key_enabled = request.headers.get(
            "x-amz-server-side-encryption-bucket-key-enabled", None
        )
        if bucket_key_enabled is not None:
            bucket_key_enabled = str(bucket_key_enabled).lower()

        bucket = self.backend.get_bucket(bucket_name)
        lock_enabled = bucket.object_lock_enabled

        lock_mode = request.headers.get("x-amz-object-lock-mode", None)
        lock_until = request.headers.get("x-amz-object-lock-retain-until-date", None)
        legal_hold = request.headers.get("x-amz-object-lock-legal-hold", None)

        if lock_mode or lock_until or legal_hold == "ON":
            if not request.headers.get("Content-Md5"):
                raise InvalidContentMD5
            if not lock_enabled:
                raise LockNotEnabled

        elif lock_enabled and bucket.has_default_lock:
            if not request.headers.get("Content-Md5"):
                raise InvalidContentMD5
            lock_until = bucket.default_retention()
            lock_mode = bucket.default_lock_mode

        acl = self._acl_from_headers(request.headers)
        if acl is None:
            acl = bucket.acl
        tagging = self._tagging_from_headers(request.headers)

        if "versionId" in query:
            version_id = query["versionId"][0]
        else:
            version_id = None

        if "retention" in query:
            if not lock_enabled:
                raise LockNotEnabled
            retention = self._mode_until_from_body()
            self.backend.put_object_retention(
                bucket_name, key_name, version_id=version_id, retention=retention
            )
            return 200, response_headers, ""

        if "legal-hold" in query:
            if not lock_enabled:
                raise LockNotEnabled
            legal_hold_status = self._legal_hold_status_from_xml(body)
            self.backend.put_object_legal_hold(
                bucket_name, key_name, version_id, legal_hold_status
            )
            return 200, response_headers, ""

        if "acl" in query:
            self.backend.put_object_acl(bucket_name, key_name, acl)
            return 200, response_headers, ""

        if "tagging" in query:
            key_to_tag = self.backend.get_object(
                bucket_name, key_name, version_id=version_id
            )
            tagging = self._tagging_from_xml(body)
            self.backend.put_object_tagging(key_to_tag, tagging, key_name)
            return 200, response_headers, ""

        if "x-amz-copy-source" in request.headers:
            # Copy key
            # you can have a quoted ?version=abc with a version Id, so work on
            # we need to parse the unquoted string first
            copy_source = request.headers.get("x-amz-copy-source")
            if isinstance(copy_source, bytes):
                copy_source = copy_source.decode("utf-8")
            copy_source_parsed = urlparse(copy_source)
            src_bucket, src_key = (
                unquote(copy_source_parsed.path).lstrip("/").split("/", 1)
            )
            src_version_id = parse_qs(copy_source_parsed.query).get(
                "versionId", [None]
            )[0]

            key_to_copy = self.backend.get_object(
                src_bucket, src_key, version_id=src_version_id
            )

            if key_to_copy is not None:
                if key_to_copy.storage_class in ARCHIVE_STORAGE_CLASSES:
                    if (
                        key_to_copy.response_dict.get("x-amz-restore") is None
                        or 'ongoing-request="true"'
                        in key_to_copy.response_dict.get(  # type: ignore
                            "x-amz-restore"
                        )
                    ):
                        raise ObjectNotInActiveTierError(key_to_copy)

                website_redirect_location = request.headers.get(
                    "x-amz-website-redirect-location"
                )

                mdirective = request.headers.get("x-amz-metadata-directive")
                metadata = metadata_from_headers(request.headers)
                self.backend.copy_object(
                    key_to_copy,
                    bucket_name,
                    key_name,
                    storage=request.headers.get("x-amz-storage-class"),
                    kms_key_id=kms_key_id,
                    encryption=encryption,
                    bucket_key_enabled=bucket_key_enabled,
                    mdirective=mdirective,
                    metadata=metadata,
                    website_redirect_location=website_redirect_location,
                    lock_mode=lock_mode,
                    lock_legal_status=legal_hold,
                    lock_until=lock_until,
                )
            else:
                raise MissingKey(key=src_key)

            new_key: FakeKey = self.backend.get_object(bucket_name, key_name)  # type: ignore

            if acl is not None:
                new_key.set_acl(acl)

            tdirective = request.headers.get("x-amz-tagging-directive")
            if tdirective == "REPLACE":
                tagging = self._tagging_from_headers(request.headers)
                self.backend.put_object_tagging(new_key, tagging)
            if key_to_copy.version_id != "null":
                response_headers["x-amz-copy-source-version-id"] = (
                    key_to_copy.version_id
                )

            # checksum stuff, do we need to compute hash of the copied object
            checksum_algorithm = request.headers.get("x-amz-checksum-algorithm")
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

            template = self.response_template(S3_OBJECT_COPY_RESPONSE)
            response_headers.update(new_key.response_dict)
            response = template.render(key=new_key)
            response_headers["content-length"] = str(len(response))
            return 200, response_headers, response

        # Initial data
        new_key = self.backend.put_object(
            bucket_name,
            key_name,
            body,
            storage=storage_class,
            encryption=encryption,
            kms_key_id=kms_key_id,
            bucket_key_enabled=bucket_key_enabled,
            lock_mode=lock_mode,
            lock_legal_status=legal_hold,
            lock_until=lock_until,
            checksum_value=checksum_value,
        )

        metadata = metadata_from_headers(request.headers)
        metadata.update(metadata_from_headers(query))
        new_key.set_metadata(metadata)
        new_key.set_acl(acl)
        new_key.website_redirect_location = request.headers.get(
            "x-amz-website-redirect-location"
        )
        if checksum_algorithm:
            new_key.checksum_algorithm = checksum_algorithm
        self.backend.put_object_tagging(new_key, tagging)

        response_headers.update(new_key.response_dict)
        # Remove content-length - the response body is empty for this request
        response_headers.pop("content-length", None)
        return 200, response_headers, ""

    def _key_response_head(
        self,
        bucket_name: str,
        query: Dict[str, Any],
        key_name: str,
        headers: Dict[str, Any],
    ) -> TYPE_RESPONSE:
        self._set_action("KEY", "HEAD", query)
        self._authenticate_and_authorize_s3_action(
            bucket_name=bucket_name, key_name=key_name
        )

        response_headers: Dict[str, Any] = {}
        version_id = query.get("versionId", [None])[0]
        if version_id and not self.backend.get_bucket(bucket_name).is_versioned:
            return 400, response_headers, ""

        part_number = query.get("partNumber", [None])[0]
        if part_number:
            part_number = int(part_number)

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
                "content-type": "application/xml",
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

            if part_number:
                full_key = self.backend.head_object(bucket_name, key_name, version_id)
                if full_key.multipart:  # type: ignore
                    mp_part_count = str(len(full_key.multipart.partlist))  # type: ignore
                    response_headers["x-amz-mp-parts-count"] = mp_part_count

            return 200, response_headers, ""
        else:
            return 404, response_headers, ""

    def _lock_config_from_body(self) -> Dict[str, Any]:
        response_dict: Dict[str, Any] = {
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

    def _acl_from_body(self) -> Optional[FakeAcl]:
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
        grant_list: List[Dict[str, Any]],
        exception_type: Type[S3ClientError],
        permissions: List[str],
    ) -> List[FakeGrant]:
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

    def _acl_from_headers(self, headers: Dict[str, str]) -> Optional[FakeAcl]:
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

    def _tagging_from_headers(self, headers: Dict[str, Any]) -> Dict[str, str]:
        tags = {}
        if headers.get("x-amz-tagging"):
            parsed_header = parse_qs(headers["x-amz-tagging"], keep_blank_values=True)
            for tag in parsed_header.items():
                tags[tag[0]] = tag[1][0]
        return tags

    def _tagging_from_xml(self, xml: bytes) -> Dict[str, str]:
        parsed_xml = xmltodict.parse(xml, force_list={"Tag": True})

        tags = {}
        for tag in parsed_xml["Tagging"]["TagSet"]["Tag"]:
            tags[tag["Key"]] = tag["Value"] or ""

        return tags

    def _bucket_tagging_from_body(self) -> Dict[str, str]:
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

    def _cors_from_body(self) -> List[Dict[str, Any]]:
        parsed_xml = xmltodict.parse(self.body)

        if isinstance(parsed_xml["CORSConfiguration"]["CORSRule"], list):
            return [cors for cors in parsed_xml["CORSConfiguration"]["CORSRule"]]

        return [parsed_xml["CORSConfiguration"]["CORSRule"]]

    def _mode_until_from_body(self) -> Tuple[Optional[str], Optional[str]]:
        parsed_xml = xmltodict.parse(self.body)
        return (
            parsed_xml.get("Retention", None).get("Mode", None),
            parsed_xml.get("Retention", None).get("RetainUntilDate", None),
        )

    def _legal_hold_status_from_xml(self, xml: bytes) -> Dict[str, Any]:
        parsed_xml = xmltodict.parse(xml)
        return parsed_xml["LegalHold"]["Status"]

    def _encryption_config_from_body(self) -> Dict[str, Any]:
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

    def _ownership_rule_from_body(self) -> Dict[str, Any]:
        parsed_xml = xmltodict.parse(self.body)

        if not parsed_xml["OwnershipControls"]["Rule"].get("ObjectOwnership"):
            raise MalformedXML()

        return parsed_xml["OwnershipControls"]["Rule"]["ObjectOwnership"]

    def _logging_from_body(self) -> Dict[str, Any]:
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

    def _notification_config_from_body(self) -> Dict[str, Any]:
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

    def _replication_config_from_xml(self, xml: str) -> Dict[str, Any]:
        parsed_xml = xmltodict.parse(xml, dict_constructor=dict)
        config = parsed_xml["ReplicationConfiguration"]
        return config

    def _key_response_delete(
        self, headers: Any, bucket_name: str, query: Dict[str, Any], key_name: str
    ) -> TYPE_RESPONSE:
        self._set_action("KEY", "DELETE", query)
        self._authenticate_and_authorize_s3_action(
            bucket_name=bucket_name, key_name=key_name
        )

        if query.get("uploadId"):
            upload_id = query["uploadId"][0]
            self.backend.abort_multipart_upload(bucket_name, upload_id)
            return 204, {}, ""
        version_id = query.get("versionId", [None])[0]
        if "tagging" in query:
            self.backend.delete_object_tagging(
                bucket_name, key_name, version_id=version_id
            )
            template = self.response_template(S3_DELETE_KEY_TAGGING_RESPONSE)
            return 204, {}, template.render(version_id=version_id)
        bypass = headers.get("X-Amz-Bypass-Governance-Retention")
        _, response_meta = self.backend.delete_object(
            bucket_name, key_name, version_id=version_id, bypass=bypass
        )
        response_headers = {}
        if response_meta is not None:
            for k in response_meta:
                response_headers[f"x-amz-{k}"] = response_meta[k]
        return 204, response_headers, ""

    def _complete_multipart_body(self, body: bytes) -> Iterator[Tuple[int, str]]:
        ps = minidom.parseString(body).getElementsByTagName("Part")
        prev = 0
        for p in ps:
            pn = int(p.getElementsByTagName("PartNumber")[0].firstChild.wholeText)  # type: ignore[union-attr]
            if pn <= prev:
                raise InvalidPartOrder()
            yield (pn, p.getElementsByTagName("ETag")[0].firstChild.wholeText)  # type: ignore[union-attr]

    def _key_response_post(
        self,
        request: Any,
        body: bytes,
        bucket_name: str,
        query: Dict[str, Any],
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

            multipart_id = self.backend.create_multipart_upload(
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

            template = self.response_template(S3_MULTIPART_INITIATE_RESPONSE)
            response = template.render(
                bucket_name=bucket_name, key_name=key_name, upload_id=multipart_id
            )
            return 200, response_headers, response

        if query.get("uploadId"):
            multipart_id = query["uploadId"][0]

            key = self.backend.complete_multipart_upload(
                bucket_name, multipart_id, self._complete_multipart_body(body)
            )
            if key is None:
                return 400, {}, ""

            headers: Dict[str, Any] = {}

            template = self.response_template(S3_MULTIPART_COMPLETE_RESPONSE)
            if key.version_id:
                headers["x-amz-version-id"] = key.version_id

            if key.encryption:
                headers["x-amz-server-side-encryption"] = key.encryption

            if key.kms_key_id:
                headers["x-amz-server-side-encryption-aws-kms-key-id"] = key.kms_key_id

            return (
                200,
                headers,
                template.render(
                    bucket_name=bucket_name, key_name=key.name, etag=key.etag
                ),
            )

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
            results = self.backend.select_object_content(
                bucket_name, key_name, select_query, input_details, output_details
            )
            return 200, {}, serialize_select(results)

        else:
            raise NotImplementedError(
                "Method POST had only been implemented for multipart uploads and restore operations, so far"
            )

    def _invalid_headers(self, url: str, headers: Dict[str, str]) -> bool:
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

S3_ALL_BUCKETS = """<ListAllMyBucketsResult xmlns="http://s3.amazonaws.com/doc/2006-03-01">
  <Owner>
    <ID>bcaf1ffd86f41161ca5fb16fd081034f</ID>
    <DisplayName>webfile</DisplayName>
  </Owner>
  <Buckets>
    {% for bucket in buckets %}
      <Bucket>
        <Name>{{ bucket.name }}</Name>
        <CreationDate>{{ bucket.creation_date_ISO8601 }}</CreationDate>
      </Bucket>
    {% endfor %}
 </Buckets>
</ListAllMyBucketsResult>"""

S3_BUCKET_GET_RESPONSE = """<?xml version="1.0" encoding="UTF-8"?>
<ListBucketResult xmlns="http://s3.amazonaws.com/doc/2006-03-01/">
  <Name>{{ bucket.name }}</Name>
  {% if prefix != None %}
  <Prefix>{{ prefix }}</Prefix>
  {% endif %}
  <MaxKeys>{{ max_keys }}</MaxKeys>
  {% if delimiter %}
    <Delimiter>{{ delimiter }}</Delimiter>
  {% endif %}
  {% if encoding_type %}
    <EncodingType>{{ encoding_type }}</EncodingType>
  {% endif %}
  <IsTruncated>{{ is_truncated }}</IsTruncated>
  {% if next_marker %}
    <NextMarker>{{ next_marker }}</NextMarker>
  {% endif %}
  {% for key in result_keys %}
    <Contents>
      <Key>{{ key.safe_name(encoding_type) }}</Key>
      <LastModified>{{ key.last_modified_ISO8601 }}</LastModified>
      <ETag>{{ key.etag }}</ETag>
      <Size>{{ key.size }}</Size>
      <StorageClass>{{ key.storage_class }}</StorageClass>
      <Owner>
        <ID>75aa57f09aa0c8caeab4f8c24e99d10f8e7faeebf76c078efc7c6caea54ba06a</ID>
        <DisplayName>webfile</DisplayName>
      </Owner>
    </Contents>
  {% endfor %}
  {% if delimiter %}
    {% for folder in result_folders %}
      <CommonPrefixes>
        <Prefix>{{ folder }}</Prefix>
      </CommonPrefixes>
    {% endfor %}
  {% endif %}
  </ListBucketResult>"""

S3_BUCKET_GET_RESPONSE_V2 = """<?xml version="1.0" encoding="UTF-8"?>
<ListBucketResult xmlns="http://s3.amazonaws.com/doc/2006-03-01/">
  <Name>{{ bucket.name }}</Name>
{% if prefix != None %}
  <Prefix>{{ prefix }}</Prefix>
{% endif %}
  <MaxKeys>{{ max_keys }}</MaxKeys>
  <KeyCount>{{ key_count }}</KeyCount>
{% if delimiter %}
  <Delimiter>{{ delimiter }}</Delimiter>
{% endif %}
{% if encoding_type %}
  <EncodingType>{{ encoding_type }}</EncodingType>
{% endif %}
  <IsTruncated>{{ is_truncated }}</IsTruncated>
{% if next_continuation_token %}
  <NextContinuationToken>{{ next_continuation_token }}</NextContinuationToken>
{% endif %}
{% if start_after %}
  <StartAfter>{{ start_after }}</StartAfter>
{% endif %}
  {% for key in result_keys %}
    <Contents>
      <Key>{{ key.safe_name(encoding_type) }}</Key>
      <LastModified>{{ key.last_modified_ISO8601 }}</LastModified>
      <ETag>{{ key.etag }}</ETag>
      <Size>{{ key.size }}</Size>
      <StorageClass>{{ key.storage_class }}</StorageClass>
      {% if fetch_owner %}
      <Owner>
        <ID>75aa57f09aa0c8caeab4f8c24e99d10f8e7faeebf76c078efc7c6caea54ba06a</ID>
        <DisplayName>webfile</DisplayName>
      </Owner>
      {% endif %}
      {% if key.checksum_algorithm %}
      <ChecksumAlgorithm>{{ key.checksum_algorithm }}</ChecksumAlgorithm>
      {% endif %}
    </Contents>
  {% endfor %}
  {% if delimiter %}
    {% for folder in result_folders %}
      <CommonPrefixes>
        <Prefix>{{ folder }}</Prefix>
      </CommonPrefixes>
    {% endfor %}
  {% endif %}
  </ListBucketResult>"""

S3_BUCKET_CREATE_RESPONSE = """<CreateBucketResponse xmlns="http://s3.amazonaws.com/doc/2006-03-01">
  <CreateBucketResponse>
    <Bucket>{{ bucket.name }}</Bucket>
  </CreateBucketResponse>
</CreateBucketResponse>"""

S3_DELETE_BUCKET_SUCCESS = """<DeleteBucketResponse xmlns="http://s3.amazonaws.com/doc/2006-03-01">
  <DeleteBucketResponse>
    <Code>204</Code>
    <Description>No Content</Description>
  </DeleteBucketResponse>
</DeleteBucketResponse>"""

S3_DELETE_BUCKET_WITH_ITEMS_ERROR = """<?xml version="1.0" encoding="UTF-8"?>
<Error><Code>BucketNotEmpty</Code>
<Message>The bucket you tried to delete is not empty</Message>
<BucketName>{{ bucket.name }}</BucketName>
<RequestId>asdfasdfsdafds</RequestId>
<HostId>sdfgdsfgdsfgdfsdsfgdfs</HostId>
</Error>"""

S3_BUCKET_LOCATION = """<?xml version="1.0" encoding="UTF-8"?>
<LocationConstraint xmlns="http://s3.amazonaws.com/doc/2006-03-01/">{% if location != None %}{{ location }}{% endif %}</LocationConstraint>"""

S3_BUCKET_LIFECYCLE_CONFIGURATION = """<?xml version="1.0" encoding="UTF-8"?>
<LifecycleConfiguration xmlns="http://s3.amazonaws.com/doc/2006-03-01/">
    {% for rule in rules %}
    <Rule>
        <ID>{{ rule.id }}</ID>
        {% if rule.filter %}
        <Filter>
            {% if rule.filter.prefix != None %}
            <Prefix>{{ rule.filter.prefix }}</Prefix>
            {% endif %}
            {% if rule.filter.tag_key %}
            <Tag>
                <Key>{{ rule.filter.tag_key }}</Key>
                <Value>{{ rule.filter.tag_value }}</Value>
            </Tag>
            {% endif %}
            {% if rule.filter.and_filter %}
            <And>
                {% if rule.filter.and_filter.prefix != None %}
                <Prefix>{{ rule.filter.and_filter.prefix }}</Prefix>
                {% endif %}
                {% for key, value in rule.filter.and_filter.tags.items() %}
                <Tag>
                    <Key>{{ key }}</Key>
                    <Value>{{ value }}</Value>
                </Tag>
                {% endfor %}
            </And>
            {% endif %}
        </Filter>
        {% else %}
            {% if rule.prefix != None %}
            <Prefix>{{ rule.prefix }}</Prefix>
            {% endif %}
        {% endif %}
        <Status>{{ rule.status }}</Status>
        {% for transition in rule.transitions %}
            <Transition>
                {% if transition.days %}
                <Days>{{ transition.days }}</Days>
                {% endif %}
                {% if transition.date %}
                <Date>{{ transition.date }}</Date>
                {% endif %}
                {% if transition.storage_class %}
                <StorageClass>{{ transition.storage_class }}</StorageClass>
                {% endif %}
            </Transition>
        {% endfor %}
        {% if rule.expiration_days or rule.expiration_date or rule.expired_object_delete_marker %}
        <Expiration>
            {% if rule.expiration_days %}
               <Days>{{ rule.expiration_days }}</Days>
            {% endif %}
            {% if rule.expiration_date %}
               <Date>{{ rule.expiration_date }}</Date>
            {% endif %}
            {% if rule.expired_object_delete_marker %}
                <ExpiredObjectDeleteMarker>{{ rule.expired_object_delete_marker }}</ExpiredObjectDeleteMarker>
            {% endif %}
        </Expiration>
        {% endif %}
        {% for nvt in rule.noncurrent_version_transitions %}
            <NoncurrentVersionTransition>
                {% if nvt.newer_versions %}
                <NewerNoncurrentVersions>{{ nvt.newer_versions }}</NewerNoncurrentVersions>
                {% endif %}
                {% if nvt.days %}
                <NoncurrentDays>{{ nvt.days }}</NoncurrentDays>
                {% endif %}
                {% if nvt.storage_class %}
                <StorageClass>{{ nvt.storage_class }}</StorageClass>
                {% endif %}
            </NoncurrentVersionTransition>
        {% endfor %}
        {% if rule.nve_noncurrent_days %}
        <NoncurrentVersionExpiration>
           <NoncurrentDays>{{ rule.nve_noncurrent_days }}</NoncurrentDays>
        </NoncurrentVersionExpiration>
        {% endif %}
        {% if rule.aimu_days %}
        <AbortIncompleteMultipartUpload>
           <DaysAfterInitiation>{{ rule.aimu_days }}</DaysAfterInitiation>
        </AbortIncompleteMultipartUpload>
        {% endif %}
    </Rule>
    {% endfor %}
</LifecycleConfiguration>
"""

S3_BUCKET_VERSIONING = """<?xml version="1.0" encoding="UTF-8"?>
<VersioningConfiguration xmlns="http://s3.amazonaws.com/doc/2006-03-01/">
    <Status>{{ bucket_versioning_status }}</Status>
</VersioningConfiguration>
"""

S3_BUCKET_GET_VERSIONING = """<?xml version="1.0" encoding="UTF-8"?>
{% if status is none %}
    <VersioningConfiguration xmlns="http://s3.amazonaws.com/doc/2006-03-01/"/>
{% else %}
    <VersioningConfiguration xmlns="http://s3.amazonaws.com/doc/2006-03-01/">
    <Status>{{ status }}</Status>
    </VersioningConfiguration>
{% endif %}
"""

S3_BUCKET_GET_VERSIONS = """<?xml version="1.0" encoding="UTF-8"?>
<ListVersionsResult xmlns="http://s3.amazonaws.com/doc/2006-03-01">
    <Name>{{ bucket.name }}</Name>
    {% if prefix != None %}
    <Prefix>{{ prefix }}</Prefix>
    {% endif %}
    {% if common_prefixes %}
    {% for prefix in common_prefixes %}
    <CommonPrefixes>
        <Prefix>{{ prefix }}</Prefix>
    </CommonPrefixes>
    {% endfor %}
    {% endif %}
    <Delimiter>{{ delimiter }}</Delimiter>
    <KeyMarker>{{ key_marker or "" }}</KeyMarker>
    <VersionIdMarker>{{ version_id_marker or "" }}</VersionIdMarker>
    <MaxKeys>{{ max_keys }}</MaxKeys>
    {% if is_truncated %}
    <IsTruncated>true</IsTruncated>
    <NextKeyMarker>{{ next_key_marker }}</NextKeyMarker>
    {% if next_version_id_marker %}
    <NextVersionIdMarker>{{ next_version_id_marker }}</NextVersionIdMarker>
    {% endif %}
    {% else %}
    <IsTruncated>false</IsTruncated>
    {% endif %}
    {% for key in key_list %}
    <Version>
        <Key>{{ key.name }}</Key>
        <VersionId>{% if key.version_id is none %}null{% else %}{{ key.version_id }}{% endif %}</VersionId>
        <IsLatest>{{ 'true' if key.is_latest else 'false' }}</IsLatest>
        <LastModified>{{ key.last_modified_ISO8601 }}</LastModified>
        <ETag>{{ key.etag }}</ETag>
        <Size>{{ key.size }}</Size>
        <StorageClass>{{ key.storage_class }}</StorageClass>
        <Owner>
            <ID>75aa57f09aa0c8caeab4f8c24e99d10f8e7faeebf76c078efc7c6caea54ba06a</ID>
            <DisplayName>webfile</DisplayName>
        </Owner>
    </Version>
    {% endfor %}
    {% for marker in delete_marker_list %}
    <DeleteMarker>
        <Key>{{ marker.name }}</Key>
        <VersionId>{{ marker.version_id }}</VersionId>
        <IsLatest>{{ 'true' if marker.is_latest else 'false' }}</IsLatest>
        <LastModified>{{ marker.last_modified_ISO8601 }}</LastModified>
        <Owner>
            <ID>75aa57f09aa0c8caeab4f8c24e99d10f8e7faeebf76c078efc7c6caea54ba06a</ID>
            <DisplayName>webfile</DisplayName>
        </Owner>
    </DeleteMarker>
    {% endfor %}
</ListVersionsResult>
"""

S3_DELETE_KEYS_RESPONSE = """<?xml version="1.0" encoding="UTF-8"?>
<DeleteResult xmlns="http://s3.amazonaws.com/doc/2006-03-01">
{% for k, v in deleted %}
<Deleted>
<Key>{{k}}</Key>
{% if v %}<VersionId>{{v}}</VersionId>{% endif %}
</Deleted>
{% endfor %}
{% for k,c,m in delete_errors %}
<Error>
<Key>{{k}}</Key>
<Code>{{c}}</Code>
<Message>{{m}}</Message>
</Error>
{% endfor %}
</DeleteResult>"""

S3_DELETE_KEY_TAGGING_RESPONSE = """<?xml version="1.0" encoding="UTF-8"?>
<DeleteObjectTaggingResult xmlns="http://s3.amazonaws.com/doc/2006-03-01">
<VersionId>{{version_id}}</VersionId>
</DeleteObjectTaggingResult>
"""

S3_OBJECT_ACL_RESPONSE = """<?xml version="1.0" encoding="UTF-8"?>
    <AccessControlPolicy xmlns="http://s3.amazonaws.com/doc/2006-03-01/">
      <Owner>
        <ID>75aa57f09aa0c8caeab4f8c24e99d10f8e7faeebf76c078efc7c6caea54ba06a</ID>
        <DisplayName>webfile</DisplayName>
      </Owner>
      <AccessControlList>
        {% for grant in acl.grants %}
        <Grant>
          {% for grantee in grant.grantees %}
          <Grantee xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
                   xsi:type="{{ grantee.type }}">
            {% if grantee.uri %}
            <URI>{{ grantee.uri }}</URI>
            {% endif %}
            {% if grantee.id %}
            <ID>{{ grantee.id }}</ID>
            {% endif %}
            {% if grantee.display_name %}
            <DisplayName>{{ grantee.display_name }}</DisplayName>
            {% endif %}
          </Grantee>
          {% endfor %}
          {% for permission in grant.permissions %}
          <Permission>{{ permission }}</Permission>
          {% endfor %}
        </Grant>
        {% endfor %}
      </AccessControlList>
    </AccessControlPolicy>"""

S3_OBJECT_LEGAL_HOLD = """<?xml version="1.0" encoding="UTF-8"?>
<LegalHold>
   <Status>{{ legal_hold }}</Status>
</LegalHold>
"""

S3_OBJECT_TAGGING_RESPONSE = """\
<?xml version="1.0" encoding="UTF-8"?>
<Tagging xmlns="http://s3.amazonaws.com/doc/2006-03-01/">
  <TagSet>
    {% for tag in tags %}
    <Tag>
      <Key>{{ tag.Key }}</Key>
      <Value>{{ tag.Value }}</Value>
    </Tag>
    {% endfor %}
  </TagSet>
</Tagging>"""

S3_BUCKET_CORS_RESPONSE = """<?xml version="1.0" encoding="UTF-8"?>
<CORSConfiguration>
  {% for cors in cors %}
  <CORSRule>
    {% for origin in cors.allowed_origins %}
    <AllowedOrigin>{{ origin }}</AllowedOrigin>
    {% endfor %}
    {% for method in cors.allowed_methods %}
    <AllowedMethod>{{ method }}</AllowedMethod>
    {% endfor %}
    {% if cors.allowed_headers is not none %}
      {% for header in cors.allowed_headers %}
      <AllowedHeader>{{ header }}</AllowedHeader>
      {% endfor %}
    {% endif %}
    {% if cors.exposed_headers is not none %}
      {% for header in cors.exposed_headers %}
      <ExposedHeader>{{ header }}</ExposedHeader>
      {% endfor %}
    {% endif %}
    {% if cors.max_age_seconds is not none %}
    <MaxAgeSeconds>{{ cors.max_age_seconds }}</MaxAgeSeconds>
    {% endif %}
  </CORSRule>
  {% endfor %}
  </CORSConfiguration>
"""

# https://docs.aws.amazon.com/AmazonS3/latest/API/API_CopyObject.html
S3_OBJECT_COPY_RESPONSE = """\
<CopyObjectResult xmlns="http://doc.s3.amazonaws.com/2006-03-01">
    <ETag>{{ key.etag }}</ETag>
    <LastModified>{{ key.last_modified_ISO8601 }}</LastModified>
    {% if key.checksum_value  %}
      {% if "CRC32" in key.checksum_algorithm %}<ChecksumCRC32>{{ key.checksum_value }}</ChecksumCRC32>{% endif %}
      {% if "CRC32C" in key.checksum_algorithm %}<ChecksumCRC32C>{{ key.checksum_value }}</ChecksumCRC32C>{% endif %}
      {% if "SHA1" in key.checksum_algorithm %}<ChecksumSHA1>{{ key.checksum_value }}</ChecksumSHA1>{% endif %}
      {% if "SHA256" in key.checksum_algorithm %}<ChecksumSHA256>{{ key.checksum_value }}</ChecksumSHA256>{% endif %}
    {% endif %}
</CopyObjectResult>"""

S3_MULTIPART_INITIATE_RESPONSE = """<?xml version="1.0" encoding="UTF-8"?>
<InitiateMultipartUploadResult xmlns="http://s3.amazonaws.com/doc/2006-03-01/">
  <Bucket>{{ bucket_name }}</Bucket>
  <Key>{{ key_name }}</Key>
  <UploadId>{{ upload_id }}</UploadId>
</InitiateMultipartUploadResult>"""

S3_MULTIPART_UPLOAD_RESPONSE = """<?xml version="1.0" encoding="UTF-8"?>
<CopyPartResult xmlns="http://s3.amazonaws.com/doc/2006-03-01/">
  <LastModified>{{ part.last_modified_ISO8601 }}</LastModified>
  <ETag>{{ part.etag }}</ETag>
</CopyPartResult>"""

S3_MULTIPART_LIST_RESPONSE = """<?xml version="1.0" encoding="UTF-8"?>
<ListPartsResult xmlns="http://s3.amazonaws.com/doc/2006-03-01/">
  <Bucket>{{ bucket_name }}</Bucket>
  <Key>{{ key_name }}</Key>
  <UploadId>{{ upload_id }}</UploadId>
  <StorageClass>STANDARD</StorageClass>
  <Initiator>
    <ID>75aa57f09aa0c8caeab4f8c24e99d10f8e7faeebf76c078efc7c6caea54ba06a</ID>
    <DisplayName>webfile</DisplayName>
  </Initiator>
  <Owner>
    <ID>75aa57f09aa0c8caeab4f8c24e99d10f8e7faeebf76c078efc7c6caea54ba06a</ID>
    <DisplayName>webfile</DisplayName>
  </Owner>
  <PartNumberMarker>{{ part_number_marker }}</PartNumberMarker>
  <NextPartNumberMarker>{{ next_part_number_marker }}</NextPartNumberMarker>
  <MaxParts>{{ max_parts }}</MaxParts>
  <IsTruncated>{{ is_truncated }}</IsTruncated>
  {% for part in parts %}
  <Part>
    <PartNumber>{{ part.name }}</PartNumber>
    <LastModified>{{ part.last_modified_ISO8601 }}</LastModified>
    <ETag>{{ part.etag }}</ETag>
    <Size>{{ part.size }}</Size>
  </Part>
  {% endfor %}
</ListPartsResult>"""

S3_MULTIPART_COMPLETE_RESPONSE = """<?xml version="1.0" encoding="UTF-8"?>
<CompleteMultipartUploadResult xmlns="http://s3.amazonaws.com/doc/2006-03-01/">
  <Location>http://{{ bucket_name }}.s3.amazonaws.com/{{ key_name }}</Location>
  <Bucket>{{ bucket_name }}</Bucket>
  <Key>{{ key_name }}</Key>
  <ETag>{{ etag }}</ETag>
</CompleteMultipartUploadResult>
"""

S3_ALL_MULTIPARTS = """<?xml version="1.0" encoding="UTF-8"?>
<ListMultipartUploadsResult xmlns="http://s3.amazonaws.com/doc/2006-03-01/">
  <Bucket>{{ bucket_name }}</Bucket>
  <KeyMarker></KeyMarker>
  <UploadIdMarker></UploadIdMarker>
  <MaxUploads>1000</MaxUploads>
  <IsTruncated>false</IsTruncated>
  {% for upload in uploads %}
  <Upload>
    <Key>{{ upload.key_name }}</Key>
    <UploadId>{{ upload.id }}</UploadId>
    <Initiator>
      <ID>arn:aws:iam::{{ account_id }}:user/user1-11111a31-17b5-4fb7-9df5-b111111f13de</ID>
      <DisplayName>user1-11111a31-17b5-4fb7-9df5-b111111f13de</DisplayName>
    </Initiator>
    <Owner>
      <ID>75aa57f09aa0c8caeab4f8c24e99d10f8e7faeebf76c078efc7c6caea54ba06a</ID>
      <DisplayName>webfile</DisplayName>
    </Owner>
    <StorageClass>STANDARD</StorageClass>
    <Initiated>2010-11-10T20:48:33.000Z</Initiated>
  </Upload>
  {% endfor %}
</ListMultipartUploadsResult>
"""

S3_NO_POLICY = """<?xml version="1.0" encoding="UTF-8"?>
<Error>
  <Code>NoSuchBucketPolicy</Code>
  <Message>The bucket policy does not exist</Message>
  <BucketName>{{ bucket_name }}</BucketName>
  <RequestId>0D68A23BB2E2215B</RequestId>
  <HostId>9Gjjt1m+cjU4OPvX9O9/8RuvnG41MRb/18Oux2o5H5MY7ISNTlXN+Dz9IG62/ILVxhAGI0qyPfg=</HostId>
</Error>
"""

S3_NO_LIFECYCLE = """<?xml version="1.0" encoding="UTF-8"?>
<Error>
  <Code>NoSuchLifecycleConfiguration</Code>
  <Message>The lifecycle configuration does not exist</Message>
  <BucketName>{{ bucket_name }}</BucketName>
  <RequestId>44425877V1D0A2F9</RequestId>
  <HostId>9Gjjt1m+cjU4OPvX9O9/8RuvnG41MRb/18Oux2o5H5MY7ISNTlXN+Dz9IG62/ILVxhAGI0qyPfg=</HostId>
</Error>
"""

S3_NO_BUCKET_TAGGING = """<?xml version="1.0" encoding="UTF-8"?>
<Error>
  <Code>NoSuchTagSet</Code>
  <Message>The TagSet does not exist</Message>
  <BucketName>{{ bucket_name }}</BucketName>
  <RequestId>44425877V1D0A2F9</RequestId>
  <HostId>9Gjjt1m+cjU4OPvX9O9/8RuvnG41MRb/18Oux2o5H5MY7ISNTlXN+Dz9IG62/ILVxhAGI0qyPfg=</HostId>
</Error>
"""

S3_NO_BUCKET_WEBSITE_CONFIG = """<?xml version="1.0" encoding="UTF-8"?>
<Error>
  <Code>NoSuchWebsiteConfiguration</Code>
  <Message>The specified bucket does not have a website configuration</Message>
  <BucketName>{{ bucket_name }}</BucketName>
  <RequestId>44425877V1D0A2F9</RequestId>
  <HostId>9Gjjt1m+cjU4OPvX9O9/8RuvnG41MRb/18Oux2o5H5MY7ISNTlXN+Dz9IG62/ILVxhAGI0qyPfg=</HostId>
</Error>
"""

S3_INVALID_CORS_REQUEST = """<?xml version="1.0" encoding="UTF-8"?>
<Error>
  <Code>NoSuchWebsiteConfiguration</Code>
  <Message>The specified bucket does not have a website configuration</Message>
  <BucketName>{{ bucket_name }}</BucketName>
  <RequestId>44425877V1D0A2F9</RequestId>
  <HostId>9Gjjt1m+cjU4OPvX9O9/8RuvnG41MRb/18Oux2o5H5MY7ISNTlXN+Dz9IG62/ILVxhAGI0qyPfg=</HostId>
</Error>
"""

S3_NO_CORS_CONFIG = """<?xml version="1.0" encoding="UTF-8"?>
<Error>
  <Code>NoSuchCORSConfiguration</Code>
  <Message>The CORS configuration does not exist</Message>
  <BucketName>{{ bucket_name }}</BucketName>
  <RequestId>44425877V1D0A2F9</RequestId>
  <HostId>9Gjjt1m+cjU4OPvX9O9/8RuvnG41MRb/18Oux2o5H5MY7ISNTlXN+Dz9IG62/ILVxhAGI0qyPfg=</HostId>
</Error>
"""

S3_LOGGING_CONFIG = """<?xml version="1.0" encoding="UTF-8"?>
<BucketLoggingStatus xmlns="http://doc.s3.amazonaws.com/2006-03-01">
  <LoggingEnabled>
    <TargetBucket>{{ logging["TargetBucket"] }}</TargetBucket>
    <TargetPrefix>{{ logging["TargetPrefix"] }}</TargetPrefix>
    {% if logging.get("TargetGrants") %}
    <TargetGrants>
      {% for grant in logging["TargetGrants"] %}
      <Grant>
        <Grantee xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
                 xsi:type="{{ grant.grantees[0].type }}">
          {% if grant.grantees[0].uri %}
          <URI>{{ grant.grantees[0].uri }}</URI>
          {% endif %}
          {% if grant.grantees[0].id %}
          <ID>{{ grant.grantees[0].id }}</ID>
          {% endif %}
          {% if grant.grantees[0].display_name %}
          <DisplayName>{{ grant.grantees[0].display_name }}</DisplayName>
          {% endif %}
        </Grantee>
        <Permission>{{ grant.permissions[0] }}</Permission>
      </Grant>
      {% endfor %}
    </TargetGrants>
    {% endif %}
  </LoggingEnabled>
</BucketLoggingStatus>
"""

S3_NO_LOGGING_CONFIG = """<?xml version="1.0" encoding="UTF-8"?>
<BucketLoggingStatus xmlns="http://doc.s3.amazonaws.com/2006-03-01" />
"""

S3_ENCRYPTION_CONFIG = """<?xml version="1.0" encoding="UTF-8"?>
<ServerSideEncryptionConfiguration xmlns="http://doc.s3.amazonaws.com/2006-03-01">
    {% if encryption %}
        <Rule>
            <ApplyServerSideEncryptionByDefault>
                <SSEAlgorithm>{{ encryption["Rule"]["ApplyServerSideEncryptionByDefault"]["SSEAlgorithm"] }}</SSEAlgorithm>
                {% if encryption["Rule"]["ApplyServerSideEncryptionByDefault"].get("KMSMasterKeyID") %}
                <KMSMasterKeyID>{{ encryption["Rule"]["ApplyServerSideEncryptionByDefault"]["KMSMasterKeyID"] }}</KMSMasterKeyID>
                {% endif %}
            </ApplyServerSideEncryptionByDefault>
            <BucketKeyEnabled>{{ 'true' if encryption["Rule"].get("BucketKeyEnabled") == 'true' else 'false' }}</BucketKeyEnabled>
        </Rule>
    {% endif %}
</ServerSideEncryptionConfiguration>
"""

S3_INVALID_PRESIGNED_PARAMETERS = """<?xml version="1.0" encoding="UTF-8"?>
<Error>
  <Code>SignatureDoesNotMatch</Code>
  <Message>The request signature we calculated does not match the signature you provided. Check your key and signing method.</Message>
  <RequestId>0D68A23BB2E2215B</RequestId>
  <HostId>9Gjjt1m+cjU4OPvX9O9/8RuvnG41MRb/18Oux2o5H5MY7ISNTlXN+Dz9IG62/ILVxhAGI0qyPfg=</HostId>
</Error>
"""

S3_NO_ENCRYPTION = """<?xml version="1.0" encoding="UTF-8"?>
<Error>
  <Code>ServerSideEncryptionConfigurationNotFoundError</Code>
  <Message>The server side encryption configuration was not found</Message>
  <BucketName>{{ bucket_name }}</BucketName>
  <RequestId>0D68A23BB2E2215B</RequestId>
  <HostId>9Gjjt1m+cjU4OPvX9O9/8RuvnG41MRb/18Oux2o5H5MY7ISNTlXN+Dz9IG62/ILVxhAGI0qyPfg=</HostId>
</Error>
"""

S3_GET_BUCKET_NOTIFICATION_CONFIG = """<?xml version="1.0" encoding="UTF-8"?>
<NotificationConfiguration xmlns="http://s3.amazonaws.com/doc/2006-03-01/">
  {% for topic in config.topic %}
  <TopicConfiguration>
    <Id>{{ topic.id }}</Id>
    <Topic>{{ topic.arn }}</Topic>
    {% for event in topic.events %}
    <Event>{{ event }}</Event>
    {% endfor %}
    {% if topic.filters %}
      <Filter>
        <S3Key>
          {% for rule in topic.filters["S3Key"]["FilterRule"] %}
          <FilterRule>
            <Name>{{ rule["Name"] }}</Name>
            <Value>{{ rule["Value"] }}</Value>
          </FilterRule>
          {% endfor %}
        </S3Key>
      </Filter>
    {% endif %}
  </TopicConfiguration>
  {% endfor %}
  {% for queue in config.queue %}
  <QueueConfiguration>
    <Id>{{ queue.id }}</Id>
    <Queue>{{ queue.arn }}</Queue>
    {% for event in queue.events %}
    <Event>{{ event }}</Event>
    {% endfor %}
    {% if queue.filters %}
      <Filter>
        <S3Key>
          {% for rule in queue.filters["S3Key"]["FilterRule"] %}
          <FilterRule>
            <Name>{{ rule["Name"] }}</Name>
            <Value>{{ rule["Value"] }}</Value>
          </FilterRule>
          {% endfor %}
        </S3Key>
      </Filter>
    {% endif %}
  </QueueConfiguration>
  {% endfor %}
  {% for cf in config.cloud_function %}
  <CloudFunctionConfiguration>
    <Id>{{ cf.id }}</Id>
    <CloudFunction>{{ cf.arn }}</CloudFunction>
    {% for event in cf.events %}
    <Event>{{ event }}</Event>
    {% endfor %}
    {% if cf.filters %}
      <Filter>
        <S3Key>
          {% for rule in cf.filters["S3Key"]["FilterRule"] %}
          <FilterRule>
            <Name>{{ rule["Name"] }}</Name>
            <Value>{{ rule["Value"] }}</Value>
          </FilterRule>
          {% endfor %}
        </S3Key>
      </Filter>
    {% endif %}
  </CloudFunctionConfiguration>
  {% endfor %}
</NotificationConfiguration>
"""

S3_BUCKET_ACCELERATE = """
<AccelerateConfiguration xmlns="http://s3.amazonaws.com/doc/2006-03-01/">
  <Status>{{ bucket.accelerate_configuration }}</Status>
</AccelerateConfiguration>
"""

S3_BUCKET_ACCELERATE_NOT_SET = """
<AccelerateConfiguration xmlns="http://s3.amazonaws.com/doc/2006-03-01/"/>
"""

S3_PUBLIC_ACCESS_BLOCK_CONFIGURATION = """
<PublicAccessBlockConfiguration>
  <BlockPublicAcls>{{public_block_config.block_public_acls}}</BlockPublicAcls>
  <IgnorePublicAcls>{{public_block_config.ignore_public_acls}}</IgnorePublicAcls>
  <BlockPublicPolicy>{{public_block_config.block_public_policy}}</BlockPublicPolicy>
  <RestrictPublicBuckets>{{public_block_config.restrict_public_buckets}}</RestrictPublicBuckets>
</PublicAccessBlockConfiguration>
"""

S3_BUCKET_LOCK_CONFIGURATION = """
<ObjectLockConfiguration xmlns="http://s3.amazonaws.com/doc/2006-03-01/">
    {%if lock_enabled %}
    <ObjectLockEnabled>Enabled</ObjectLockEnabled>
    {% else %}
    <ObjectLockEnabled>Disabled</ObjectLockEnabled>
    {% endif %}
    {% if mode %}
    <Rule>
        <DefaultRetention>
            <Mode>{{mode}}</Mode>
            <Days>{{days}}</Days>
            <Years>{{years}}</Years>
        </DefaultRetention>
    </Rule>
    {% endif %}
</ObjectLockConfiguration>
"""

S3_DUPLICATE_BUCKET_ERROR = """<?xml version="1.0" encoding="UTF-8"?>
<Error>
  <Code>BucketAlreadyOwnedByYou</Code>
  <Message>Your previous request to create the named bucket succeeded and you already own it.</Message>
  <BucketName>{{ bucket_name }}</BucketName>
  <RequestId>44425877V1D0A2F9</RequestId>
  <HostId>9Gjjt1m+cjU4OPvX9O9/8RuvnG41MRb/18Oux2o5H5MY7ISNTlXN+Dz9IG62/ILVxhAGI0qyPfg=</HostId>
</Error>
"""

S3_NO_REPLICATION = """<?xml version="1.0" encoding="UTF-8"?>
<Error>
  <Code>ReplicationConfigurationNotFoundError</Code>
  <Message>The replication configuration was not found</Message>
  <BucketName>{{ bucket_name }}</BucketName>
  <RequestId>ZM6MA8EGCZ1M9EW9</RequestId>
  <HostId>SMUZFedx1CuwjSaZQnM2bEVpet8UgX9uD/L7e MlldClgtEICTTVFz3C66cz8Bssci2OsWCVlog=</HostId>
</Error>
"""

S3_NO_VERSIONING_ENABLED = """<?xml version="1.0" encoding="UTF-8"?>
<Error>
  <Code>InvalidRequest</Code>
  <Message>Versioning must be 'Enabled' on the bucket to apply a replication configuration</Message>
  <BucketName>{{ bucket_name }}</BucketName>
  <RequestId>ZM6MA8EGCZ1M9EW9</RequestId>
  <HostId>SMUZFedx1CuwjSaZQnM2bEVpet8UgX9uD/L7e MlldClgtEICTTVFz3C66cz8Bssci2OsWCVlog=</HostId>
</Error>
"""

S3_REPLICATION_CONFIG = """<?xml version="1.0" encoding="UTF-8"?>
<ReplicationConfiguration xmlns="http://s3.amazonaws.com/doc/2006-03-01/">
{% for rule in replication["Rule"] %}
<Rule>
  <ID>{{ rule["ID"] }}</ID>
  <Priority>{{ rule["Priority"] }}</Priority>
  <Status>{{ rule["Status"] }}</Status>
  <DeleteMarkerReplication>
    <Status>Disabled</Status>
  </DeleteMarkerReplication>
  <Filter>
    <Prefix></Prefix>
  </Filter>
  <Destination>
    <Bucket>{{ rule["Destination"]["Bucket"] }}</Bucket>
  </Destination>
</Rule>
{% endfor %}
<Role>{{ replication["Role"] }}</Role>
</ReplicationConfiguration>
"""

S3_BUCKET_GET_OWNERSHIP_RULE = """<?xml version="1.0" encoding="UTF-8"?>
<OwnershipControls xmlns="http://s3.amazonaws.com/doc/2006-03-01/">
    <Rule>
        <ObjectOwnership>{{ownership_rule}}</ObjectOwnership>
    </Rule>
</OwnershipControls>
"""

S3_ERROR_BUCKET_ONWERSHIP_NOT_FOUND = """
<Error>
    <Code>OwnershipControlsNotFoundError</Code>
    <Message>The bucket ownership controls were not found</Message>
    <BucketName>{{bucket_name}}</BucketName>
    <RequestId>294PFVCB9GFVXY2S</RequestId>
    <HostId>l/tqqyk7HZbfvFFpdq3+CAzA9JXUiV4ZajKYhwolOIpnmlvZrsI88AKsDLsgQI6EvZ9MuGHhk7M=</HostId>
</Error>
"""


S3_OBJECT_ATTRIBUTES_RESPONSE = """<?xml version="1.0" encoding="UTF-8"?>
<GetObjectAttributesOutput xmlns="http://s3.amazonaws.com/doc/2006-03-01/">
    {% if etag is not none %}<ETag>{{ etag }}</ETag>{% endif %}
    {% if checksum is not none %}
      <Checksum>
      {% if "CRC32" in checksum %}<ChecksumCRC32>{{ checksum["CRC32"] }}</ChecksumCRC32>{% endif %}
      {% if "CRC32C" in checksum %}<ChecksumCRC32C>{{ checksum["CRC32C"] }}</ChecksumCRC32C>{% endif %}
      {% if "SHA1" in checksum %}<ChecksumSHA1>{{ checksum["SHA1"] }}</ChecksumSHA1>{% endif %}
      {% if "SHA256" in checksum %}<ChecksumSHA256>{{ checksum["SHA256"] }}</ChecksumSHA256>{% endif %}
      </Checksum>
    {% endif %}
    {% if size is not none %}<ObjectSize>{{ size }}</ObjectSize>{% endif %}
    {% if storage_class is not none %}<StorageClass>{{ storage_class }}</StorageClass>{% endif %}
</GetObjectAttributesOutput>
"""
