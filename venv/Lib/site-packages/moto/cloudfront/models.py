import string
from typing import Any, Dict, Iterable, List, Optional, Tuple

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import iso_8601_datetime_with_milliseconds
from moto.moto_api._internal import mock_random as random
from moto.moto_api._internal.managed_state_model import ManagedState
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import PARTITION_NAMES, get_partition

from .exceptions import (
    DistributionAlreadyExists,
    DomainNameNotAnS3Bucket,
    InvalidIfMatchVersion,
    NoSuchDistribution,
    NoSuchInvalidation,
    NoSuchOriginAccessControl,
    OriginDoesNotExist,
)


def random_id(uppercase: bool = True, length: int = 13) -> str:
    ascii_set = string.ascii_uppercase if uppercase else string.ascii_lowercase
    chars = list(range(10)) + list(ascii_set)
    resource_id = random.choice(ascii_set) + "".join(
        str(random.choice(chars)) for _ in range(length - 1)
    )
    return resource_id


class ActiveTrustedSigners:
    def __init__(self) -> None:
        self.enabled = False
        self.quantity = 0
        self.signers: List[Any] = []


class ActiveTrustedKeyGroups:
    def __init__(self) -> None:
        self.enabled = False
        self.quantity = 0
        self.kg_key_pair_ids: List[Any] = []


class LambdaFunctionAssociation:
    def __init__(self) -> None:
        self.arn = ""
        self.event_type = ""
        self.include_body = False


class ForwardedValues:
    def __init__(self, config: Dict[str, Any]):
        self.query_string = config.get("QueryString", "false")
        self.cookie_forward = config.get("Cookies", {}).get("Forward") or "none"
        self.whitelisted_names = (
            config.get("Cookies", {}).get("WhitelistedNames", {}).get("Items") or {}
        )
        self.whitelisted_names = self.whitelisted_names.get("Name") or []
        if isinstance(self.whitelisted_names, str):
            self.whitelisted_names = [self.whitelisted_names]
        self.headers: List[Any] = []
        self.query_string_cache_keys: List[Any] = []
        self.cookies: List[Dict[str, Any]] = config.get("Cookies") or []


class TrustedKeyGroups:
    def __init__(self, config: Dict[str, Any]):
        items = config.get("Items") or {}
        self.group_ids = items.get("KeyGroup") or []
        if isinstance(self.group_ids, str):
            self.group_ids = [self.group_ids]


class DefaultCacheBehaviour:
    def __init__(self, config: Dict[str, Any]):
        self.target_origin_id = config["TargetOriginId"]
        self.trusted_signers_enabled = False
        self.trusted_signers: List[Any] = []
        self.trusted_key_groups = TrustedKeyGroups(config.get("TrustedKeyGroups") or {})
        self.viewer_protocol_policy = config["ViewerProtocolPolicy"]
        methods = config.get("AllowedMethods", {})
        self.allowed_methods = methods.get("Items", {}).get("Method", ["HEAD", "GET"])
        self.cached_methods = (
            methods.get("CachedMethods", {})
            .get("Items", {})
            .get("Method", ["GET", "HEAD"])
        )
        self.smooth_streaming = config.get("SmoothStreaming") or True
        self.compress = config.get("Compress", "true").lower() == "true"
        self.lambda_function_associations: List[Any] = []
        self.function_associations: List[Any] = []
        self.field_level_encryption_id = config.get("FieldLevelEncryptionId") or ""
        self.forwarded_values = ForwardedValues(config.get("ForwardedValues", {}))
        self.min_ttl = config.get("MinTTL") or 0
        self.default_ttl = config.get("DefaultTTL") or 0
        self.max_ttl = config.get("MaxTTL") or 0
        self.realtime_log_config_arn = config.get("RealtimeLogConfigArn") or ""


class CacheBehaviour(DefaultCacheBehaviour):
    def __init__(self, config: Dict[str, Any]):
        super().__init__(config)
        self.path_pattern: str = config.get("PathPattern", "")
        methods = config.get("AllowedMethods", {})
        self.cached_methods: List[str] = (
            methods.get("CachedMethods", {}).get("Items", {}).get("Method", [])
        )
        self.allowed_methods: List[str] = methods.get("Items", {}).get("Method", [])
        self.cache_policy_id = config.get("CachePolicyId", "")
        self.origin_request_policy_id = config.get("OriginRequestPolicyId", "")


class Logging:
    def __init__(self, config: Dict[str, Any]) -> None:
        self.enabled = config.get("Enabled") or False
        self.include_cookies = config.get("IncludeCookies") or False
        self.bucket = config.get("Bucket") or ""
        self.prefix = config.get("Prefix") or ""


class ViewerCertificate:
    def __init__(self) -> None:
        self.cloud_front_default_certificate = True
        self.min_protocol_version = "TLSv1"
        self.certificate_source = "cloudfront"


class CustomOriginConfig:
    def __init__(self, config: Dict[str, Any]):
        self.http_port = config.get("HTTPPort")
        self.https_port = config.get("HTTPSPort")
        self.keep_alive = config.get("OriginKeepaliveTimeout") or 5
        self.protocol_policy = config.get("OriginProtocolPolicy")
        self.read_timeout = config.get("OriginReadTimeout") or 30
        protocols = config.get("OriginSslProtocols", {}).get("Items") or {}
        self.ssl_protocols = protocols.get("SslProtocol") or []
        if isinstance(self.ssl_protocols, str):
            self.ssl_protocols = [self.ssl_protocols]


class Origin:
    def __init__(self, origin: Dict[str, Any]):
        self.id = origin["Id"]
        self.domain_name = origin["DomainName"]
        self.origin_path = origin.get("OriginPath") or ""
        self.s3_access_identity = ""
        self.custom_origin = None
        self.origin_shield = origin.get("OriginShield")
        self.connection_attempts = origin.get("ConnectionAttempts") or 3
        self.connection_timeout = origin.get("ConnectionTimeout") or 10

        if "S3OriginConfig" in origin:
            # Very rough validation
            if not self.domain_name.endswith("amazonaws.com"):
                raise DomainNameNotAnS3Bucket
            self.s3_access_identity = origin["S3OriginConfig"]["OriginAccessIdentity"]

        if "CustomOriginConfig" in origin:
            self.custom_origin = CustomOriginConfig(origin["CustomOriginConfig"])

        custom_headers = origin.get("CustomHeaders") or {}
        custom_headers = custom_headers.get("Items") or {}
        custom_headers = custom_headers.get("OriginCustomHeader") or []
        if isinstance(custom_headers, dict):
            # Happens if user only sends a single header
            custom_headers = [custom_headers]
        self.custom_headers = custom_headers


class GeoRestrictions:
    def __init__(self, config: Dict[str, Any]):
        config = config.get("GeoRestriction") or {}
        self._type = config.get("RestrictionType", "none")
        self.restrictions = (config.get("Items") or {}).get("Location") or []


class DistributionConfig:
    def __init__(self, config: Dict[str, Any]):
        self.config = config
        self.aliases = ((config.get("Aliases") or {}).get("Items") or {}).get(
            "CNAME"
        ) or []
        if isinstance(self.aliases, str):
            self.aliases = [self.aliases]
        self.comment = config.get("Comment") or ""
        self.default_cache_behavior = DefaultCacheBehaviour(
            config["DefaultCacheBehavior"]
        )
        self.cache_behaviors: List[Any] = []
        if config.get("CacheBehaviors", {}).get("Items"):
            for _, v in config.get("CacheBehaviors", {}).get("Items").items():
                self.cache_behaviors.append(CacheBehaviour(v))
        self.custom_error_responses: List[Any] = []
        self.logging = Logging(config.get("Logging") or {})
        self.enabled = config.get("Enabled") or False
        self.viewer_certificate = ViewerCertificate()
        self.geo_restriction = GeoRestrictions(config.get("Restrictions") or {})
        self.caller_reference = config.get("CallerReference", str(random.uuid4()))
        self.origins = config["Origins"]["Items"]["Origin"]
        if not isinstance(self.origins, list):
            self.origins = [self.origins]

        # This check happens before any other Origins-validation
        if self.default_cache_behavior.target_origin_id not in [
            o["Id"] for o in self.origins
        ]:
            raise OriginDoesNotExist

        self.origins = [Origin(o) for o in self.origins]
        self.price_class = config.get("PriceClass", "PriceClass_All")
        self.http_version = config.get("HttpVersion", "http2")
        self.is_ipv6_enabled = config.get("IsIPV6Enabled", "true").lower() == "true"
        self.default_root_object = config.get("DefaultRootObject") or ""
        self.web_acl_id = config.get("WebACLId") or ""


class Distribution(BaseModel, ManagedState):
    def __init__(self, account_id: str, region_name: str, config: Dict[str, Any]):
        # Configured ManagedState
        super().__init__(
            "cloudfront::distribution", transitions=[("InProgress", "Deployed")]
        )
        # Configure internal properties
        self.distribution_id = random_id()
        self.arn = f"arn:{get_partition(region_name)}:cloudfront:{account_id}:distribution/{self.distribution_id}"
        self.distribution_config = DistributionConfig(config)
        self.active_trusted_signers = ActiveTrustedSigners()
        self.active_trusted_key_groups = ActiveTrustedKeyGroups()
        self.origin_groups: List[Any] = []
        self.alias_icp_recordals: List[Any] = []
        self.last_modified_time = "2021-11-27T10:34:26.802Z"
        self.in_progress_invalidation_batches = 0
        self.has_active_trusted_key_groups = False
        self.domain_name = f"{random_id(uppercase=False)}.cloudfront.net"
        self.etag = random_id()

    @property
    def location(self) -> str:
        return f"https://cloudfront.amazonaws.com/2020-05-31/distribution/{self.distribution_id}"


class OriginAccessControl(BaseModel):
    def __init__(self, config_dict: Dict[str, str]):
        self.id = random_id()
        self.name = config_dict.get("Name")
        self.description = config_dict.get("Description")
        self.signing_protocol = config_dict.get("SigningProtocol")
        self.signing_behaviour = config_dict.get("SigningBehavior")
        self.origin_type = config_dict.get("OriginAccessControlOriginType")
        self.etag = random_id()

    def update(self, config: Dict[str, str]) -> None:
        if "Name" in config:
            self.name = config["Name"]
        if "Description" in config:
            self.description = config["Description"]
        if "SigningProtocol" in config:
            self.signing_protocol = config["SigningProtocol"]
        if "SigningBehavior" in config:
            self.signing_behaviour = config["SigningBehavior"]
        if "OriginAccessControlOriginType" in config:
            self.origin_type = config["OriginAccessControlOriginType"]


class Invalidation(BaseModel):
    def __init__(
        self, distribution: Distribution, paths: Dict[str, Any], caller_ref: str
    ):
        self.invalidation_id = random_id()
        self.create_time = iso_8601_datetime_with_milliseconds()
        self.distribution = distribution
        self.status = "COMPLETED"

        self.paths = paths
        self.caller_ref = caller_ref

    @property
    def location(self) -> str:
        return self.distribution.location + f"/invalidation/{self.invalidation_id}"


class PublicKey(BaseModel):
    def __init__(self, caller_ref: str, name: str, encoded_key: str):
        self.id = random_id(length=14)
        self.caller_ref = caller_ref
        self.name = name
        self.encoded_key = encoded_key
        self.created = iso_8601_datetime_with_milliseconds()
        self.etag = random_id(length=14)
        self.location = (
            f"https://cloudfront.amazonaws.com/2020-05-31/public-key/{self.id}"
        )

        # Last newline-separator is lost in the XML->Python transformation, but should exist
        if not self.encoded_key.endswith("\n"):
            self.encoded_key += "\n"


class KeyGroup(BaseModel):
    def __init__(self, name: str, items: List[str]):
        self.id = random_id(length=14)
        self.name = name
        self.items = items
        self.etag = random_id(length=14)
        self.location = (
            f"https://cloudfront.amazonaws.com/2020-05-31/key-group/{self.id}"
        )


class CloudFrontBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.distributions: Dict[str, Distribution] = dict()
        self.invalidations: Dict[str, List[Invalidation]] = dict()
        self.origin_access_controls: Dict[str, OriginAccessControl] = dict()
        self.public_keys: Dict[str, PublicKey] = dict()
        self.key_groups: Dict[str, KeyGroup] = dict()
        self.tagger = TaggingService()

    def create_distribution(
        self, distribution_config: Dict[str, Any], tags: List[Dict[str, str]]
    ) -> Tuple[Distribution, str, str]:
        """
        Not all configuration options are supported yet.  Please raise an issue if
        we're not persisting/returning the correct attributes for your
        use-case.
        """
        # We'll always call dist_with_tags, as the incoming request is the same
        return self.create_distribution_with_tags(distribution_config, tags)

    def create_distribution_with_tags(
        self, distribution_config: Dict[str, Any], tags: List[Dict[str, str]]
    ) -> Tuple[Distribution, str, str]:
        dist = Distribution(self.account_id, self.region_name, distribution_config)
        caller_reference = dist.distribution_config.caller_reference
        existing_dist = self._distribution_with_caller_reference(caller_reference)
        if existing_dist is not None:
            raise DistributionAlreadyExists(existing_dist.distribution_id)
        self.distributions[dist.distribution_id] = dist
        self.tagger.tag_resource(dist.arn, tags)
        return dist, dist.location, dist.etag

    def get_distribution(self, distribution_id: str) -> Tuple[Distribution, str]:
        if distribution_id not in self.distributions:
            raise NoSuchDistribution
        dist = self.distributions[distribution_id]
        dist.advance()
        return dist, dist.etag

    def get_distribution_config(self, distribution_id: str) -> Tuple[Distribution, str]:
        if distribution_id not in self.distributions:
            raise NoSuchDistribution
        dist = self.distributions[distribution_id]
        dist.advance()
        return dist, dist.etag

    def delete_distribution(self, distribution_id: str, if_match: bool) -> None:
        """
        The IfMatch-value is ignored - any value is considered valid.
        Calling this function without a value is invalid, per AWS' behaviour
        """
        if not if_match:
            raise InvalidIfMatchVersion
        if distribution_id not in self.distributions:
            raise NoSuchDistribution
        del self.distributions[distribution_id]

    def list_distributions(self) -> Iterable[Distribution]:
        """
        Pagination is not supported yet.
        """
        for dist in self.distributions.values():
            dist.advance()
        return self.distributions.values()

    def _distribution_with_caller_reference(
        self, reference: str
    ) -> Optional[Distribution]:
        for dist in self.distributions.values():
            config = dist.distribution_config
            if config.caller_reference == reference:
                return dist
        return None

    def update_distribution(
        self, dist_config: Dict[str, Any], _id: str, if_match: bool
    ) -> Tuple[Distribution, str, str]:
        """
        The IfMatch-value is ignored - any value is considered valid.
        Calling this function without a value is invalid, per AWS' behaviour
        """
        if _id not in self.distributions or _id is None:
            raise NoSuchDistribution
        if not if_match:
            raise InvalidIfMatchVersion
        if not dist_config:
            raise NoSuchDistribution
        dist = self.distributions[_id]

        dist.distribution_config = DistributionConfig(dist_config)
        self.distributions[_id] = dist
        dist.advance()
        return dist, dist.location, dist.etag

    def create_invalidation(
        self, dist_id: str, paths: Dict[str, Any], caller_ref: str
    ) -> Invalidation:
        dist, _ = self.get_distribution(dist_id)
        invalidation = Invalidation(dist, paths, caller_ref)
        try:
            self.invalidations[dist_id].append(invalidation)
        except KeyError:
            self.invalidations[dist_id] = [invalidation]

        return invalidation

    def list_invalidations(self, dist_id: str) -> Iterable[Invalidation]:
        """
        Pagination is not yet implemented
        """
        return self.invalidations.get(dist_id) or []

    def get_invalidation(self, dist_id: str, id: str) -> Invalidation:
        if dist_id not in self.distributions:
            raise NoSuchDistribution
        try:
            invalidations = self.invalidations[dist_id]
            if invalidations:
                for invalidation in invalidations:
                    if invalidation.invalidation_id == id:
                        return invalidation
        except KeyError:
            pass
        raise NoSuchInvalidation

    def list_tags_for_resource(self, resource: str) -> Dict[str, List[Dict[str, str]]]:
        return self.tagger.list_tags_for_resource(resource)

    def create_origin_access_control(
        self, config_dict: Dict[str, str]
    ) -> OriginAccessControl:
        control = OriginAccessControl(config_dict)
        self.origin_access_controls[control.id] = control
        return control

    def get_origin_access_control(self, control_id: str) -> OriginAccessControl:
        if control_id not in self.origin_access_controls:
            raise NoSuchOriginAccessControl
        return self.origin_access_controls[control_id]

    def update_origin_access_control(
        self, control_id: str, config: Dict[str, str]
    ) -> OriginAccessControl:
        """
        The IfMatch-parameter is not yet implemented
        """
        control = self.get_origin_access_control(control_id)
        control.update(config)
        return control

    def list_origin_access_controls(self) -> Iterable[OriginAccessControl]:
        """
        Pagination is not yet implemented
        """
        return self.origin_access_controls.values()

    def delete_origin_access_control(self, control_id: str) -> None:
        """
        The IfMatch-parameter is not yet implemented
        """
        self.origin_access_controls.pop(control_id)

    def create_public_key(
        self, caller_ref: str, name: str, encoded_key: str
    ) -> PublicKey:
        key = PublicKey(name=name, caller_ref=caller_ref, encoded_key=encoded_key)
        self.public_keys[key.id] = key
        return key

    def get_public_key(self, key_id: str) -> PublicKey:
        return self.public_keys[key_id]

    def delete_public_key(self, key_id: str) -> None:
        """
        IfMatch is not yet implemented - deletion always succeeds
        """
        self.public_keys.pop(key_id, None)

    def list_public_keys(self) -> List[PublicKey]:
        """
        Pagination is not yet implemented
        """
        return list(self.public_keys.values())

    def create_key_group(self, name: str, items: List[str]) -> KeyGroup:
        key_group = KeyGroup(name=name, items=items)
        self.key_groups[key_group.id] = key_group
        return key_group

    def get_key_group(self, group_id: str) -> KeyGroup:
        return self.key_groups[group_id]

    def list_key_groups(self) -> List[KeyGroup]:
        """
        Pagination is not yet implemented
        """
        return list(self.key_groups.values())


cloudfront_backends = BackendDict(
    CloudFrontBackend,
    "cloudfront",
    use_boto3_regions=False,
    additional_regions=PARTITION_NAMES,
)
