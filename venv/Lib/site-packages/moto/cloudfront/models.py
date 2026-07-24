import string
from collections.abc import Iterator
from typing import Any

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.resource_tagging import TaggableResourcesMixin, TaggedResource
from moto.core.utils import utcnow
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
        self.items: list[Any] = []


class ActiveTrustedKeyGroups:
    def __init__(self) -> None:
        self.enabled = False
        self.quantity = 0
        self.items: list[Any] = []


class LambdaFunctionAssociation:
    def __init__(self) -> None:
        self.arn = ""
        self.event_type = ""
        self.include_body = False


class ForwardedValues:
    def __init__(self, config: dict[str, Any]):
        if "QueryString" not in config:
            config["QueryString"] = False
        if "Cookies" not in config:
            config["Cookies"] = {"Forward": "none"}
        self.__dict__.update(config)


class TrustedSigners:
    def __init__(self, config: dict[str, Any]):
        self.items = config.get("Items", [])
        self.quantity = len(self.items)
        self.enabled = True if self.quantity else False


class TrustedKeyGroups:
    def __init__(self, config: dict[str, Any]):
        self.items = config.get("Items", [])
        self.quantity = len(self.items)
        self.enabled = True if self.quantity > 0 else False


class AllowedMethods:
    def __init__(self, config: dict[str, Any]):
        self.items = config.get("Items", ["HEAD", "GET"])
        self.quantity = len(self.items)
        cached_methods = config.get("CachedMethods", {})
        cached_methods_items = cached_methods.get("Items", ["GET", "HEAD"])
        self.cached_methods = {
            "Items": cached_methods_items,
            "Quantity": len(cached_methods_items),
        }


class DefaultCacheBehaviour:
    def __init__(self, config: dict[str, Any]):
        self.target_origin_id = config["TargetOriginId"]
        self.trusted_signers_enabled = False
        self.trusted_signers = TrustedSigners(config.get("TrustedSigners") or {})
        self.trusted_key_groups = TrustedKeyGroups(config.get("TrustedKeyGroups") or {})
        self.viewer_protocol_policy = config["ViewerProtocolPolicy"]
        methods = config.get("AllowedMethods", {})
        self.allowed_methods = AllowedMethods(methods)
        self.smooth_streaming = config.get("SmoothStreaming", False)
        self.compress = config.get("Compress", True)
        self.lambda_function_associations = {"Quantity": 0}
        self.function_associations = {"Quantity": 0}
        self.field_level_encryption_id = config.get("FieldLevelEncryptionId") or ""
        self.forwarded_values = ForwardedValues(config.get("ForwardedValues", {}))
        self.min_ttl = config.get("MinTTL") or 0
        self.default_ttl = config.get("DefaultTTL") or 0
        self.max_ttl = config.get("MaxTTL") or 0
        self.realtime_log_config_arn = config.get("RealtimeLogConfigArn") or ""
        self.cache_policy_id = config.get("CachePolicyId", "")
        self.origin_request_policy_id = config.get("OriginRequestPolicyId")
        self.response_headers_policy_id = config.get("ResponseHeadersPolicyId")


class CacheBehaviour(DefaultCacheBehaviour):
    def __init__(self, config: dict[str, Any]):
        super().__init__(config)
        self.path_pattern: str = config.get("PathPattern", "")
        methods = config.get("AllowedMethods", {})
        self.allowed_methods = AllowedMethods(methods)
        self.cache_policy_id = config.get("CachePolicyId", "")
        self.origin_request_policy_id = config.get("OriginRequestPolicyId", "")


class Logging:
    def __init__(self, config: dict[str, Any]) -> None:
        self.enabled = config.get("Enabled") or False
        self.include_cookies = config.get("IncludeCookies") or False
        self.bucket = config.get("Bucket") or ""
        self.prefix = config.get("Prefix") or ""


class ViewerCertificate:
    def __init__(self, config: dict[str, Any]) -> None:
        self.cloud_front_default_certificate = config.get(
            "CloudFrontDefaultCertificate", True
        )
        self.iam_certificate_id = config.get("IAMCertificateId") or ""
        self.acm_certificate_arn = config.get("ACMCertificateArn") or ""
        self.ssl_support_method = config.get("SSLSupportMethod") or "sni-only"
        self.minimum_protocol_version = config.get("MinimumProtocolVersion") or "TLSv1"
        self.certificate_source = "cloudfront"
        self.certificate = config.get("Certificate", "")


class CustomOriginConfig:
    def __init__(self, config: dict[str, Any]):
        self.http_port = config.get("HTTPPort")
        self.https_port = config.get("HTTPSPort")
        self.origin_keepalive_timeout = config.get("OriginKeepaliveTimeout") or 5
        self.origin_protocol_policy = config.get("OriginProtocolPolicy")
        self.origin_read_timeout = config.get("OriginReadTimeout") or 30
        protocols = config.get("OriginSslProtocols", {}).get("Items", [])
        self.origin_ssl_protocols = {
            "Quantity": len(protocols),
            "Items": protocols,
        }


class Origin:
    def __init__(self, origin: dict[str, Any]):
        self.id = origin["Id"]
        self.domain_name = origin["DomainName"]
        self.origin_path = origin.get("OriginPath") or ""
        self.s3_access_identity = ""
        self.custom_origin = None
        if "OriginShield" not in origin:
            origin["OriginShield"] = {"Enabled": False}
        self.origin_shield = origin["OriginShield"]
        self.connection_attempts = origin.get("ConnectionAttempts") or 3
        self.connection_timeout = origin.get("ConnectionTimeout") or 10

        if "S3OriginConfig" in origin:
            # Very rough validation
            if not self.domain_name.endswith("amazonaws.com"):
                raise DomainNameNotAnS3Bucket
            self.s3_access_identity = origin["S3OriginConfig"]["OriginAccessIdentity"]

        if "CustomOriginConfig" in origin:
            self.custom_origin_config = CustomOriginConfig(origin["CustomOriginConfig"])

        if "CustomHeaders" not in origin:
            origin["CustomHeaders"] = {"Quantity": 0, "Items": []}
        self.custom_headers = origin["CustomHeaders"]


class GeoRestriction:
    def __init__(self, config: dict[str, Any]):
        self.restriction_type = config.get("RestrictionType", "none")
        self.items = config.get("Items", [])
        self.quantity = len(self.items)
        if not self.quantity:
            self.items = None


class DistributionConfig:
    def __init__(self, config: dict[str, Any]):
        if "Aliases" not in config:
            config["Aliases"] = {"Quantity": 0}
        else:
            config["Aliases"]["Quantity"] = len(config["Aliases"].get("Items", []))
        if "OriginGroups" not in config:
            config["OriginGroups"] = {"Quantity": 0}
        if "CustomErrorResponses" not in config:
            config["CustomErrorResponses"] = {"Quantity": 0}
        if "ViewerCertificate" not in config:
            config["ViewerCertificate"] = ViewerCertificate({})
        else:
            config["ViewerCertificate"] = ViewerCertificate(config["ViewerCertificate"])
        config["Origins"]["Items"] = [Origin(o) for o in config["Origins"]["Items"]]
        if "CacheBehaviors" not in config:
            config["CacheBehaviors"] = {"Quantity": 0}
        elif config["CacheBehaviors"].get("Quantity"):
            config["CacheBehaviors"]["Items"] = [
                CacheBehaviour(cb) for cb in config["CacheBehaviors"]["Items"]
            ]
        if "Restrictions" not in config:
            config["Restrictions"] = {"GeoRestriction": GeoRestriction({})}
        elif config.get("Restrictions", {}).get("GeoRestriction"):
            config["Restrictions"]["GeoRestriction"] = GeoRestriction(
                config["Restrictions"]["GeoRestriction"]
            )
        config["Logging"] = Logging(config.get("Logging", {}))
        config["DefaultCacheBehavior"] = DefaultCacheBehaviour(
            config.get("DefaultCacheBehavior", {})
        )
        config["PriceClass"] = config.get("PriceClass", "PriceClass_All")
        config["HttpVersion"] = config.get("HttpVersion", "http2")
        config["IsIPV6Enabled"] = config.get("IsIPV6Enabled", True)
        config["DefaultRootObject"] = config.get("DefaultRootObject", "")
        config["WebACLId"] = config.get("WebACLId", "")
        if config["DefaultCacheBehavior"].target_origin_id not in [
            o.id for o in config["Origins"]["Items"]
        ]:
            raise OriginDoesNotExist
        self.__dict__.update(config)
        # HACK: this attribute is referenced in backend methods.
        self.caller_reference = config["CallerReference"]


class Distribution(BaseModel, ManagedState):
    def __init__(self, account_id: str, region_name: str, config: dict[str, Any]):
        # Configured ManagedState
        super().__init__(
            "cloudfront::distribution", transitions=[("InProgress", "Deployed")]
        )
        # Configure internal properties
        self.distribution_id = random_id()
        self.id = self.distribution_id
        self.arn = f"arn:{get_partition(region_name)}:cloudfront::{account_id}:distribution/{self.distribution_id}"
        self.distribution_config = DistributionConfig(config)
        self.active_trusted_signers = ActiveTrustedSigners()
        self.active_trusted_key_groups = ActiveTrustedKeyGroups()
        self.origin_groups: list[Any] = []
        self.alias_icp_recordals: list[Any] = []
        self.last_modified_time = "2021-11-27T10:34:26.802Z"
        self.in_progress_invalidation_batches = 0
        self.has_active_trusted_key_groups = False
        self.domain_name = f"{random_id(uppercase=False)}.cloudfront.net"
        self.etag = random_id()

    @property
    def location(self) -> str:
        return f"https://cloudfront.amazonaws.com/2020-05-31/distribution/{self.distribution_id}"


class OriginAccessControl(BaseModel):
    def __init__(self, config_dict: dict[str, str]):
        self.id = random_id()
        self.name = config_dict.get("Name")
        self.description = config_dict.get("Description")
        self.signing_protocol = config_dict.get("SigningProtocol")
        self.signing_behavior = config_dict.get("SigningBehavior")
        self.origin_type = config_dict.get("OriginAccessControlOriginType")
        self.etag = random_id()

    def update(self, config: dict[str, str]) -> None:
        if "Name" in config:
            self.name = config["Name"]
        if "Description" in config:
            self.description = config["Description"]
        if "SigningProtocol" in config:
            self.signing_protocol = config["SigningProtocol"]
        if "SigningBehavior" in config:
            self.signing_behavior = config["SigningBehavior"]
        if "OriginAccessControlOriginType" in config:
            self.origin_type = config["OriginAccessControlOriginType"]


class Invalidation(BaseModel):
    def __init__(self, distribution: Distribution, paths: list[str], caller_ref: str):
        self.id = random_id()
        self.create_time = utcnow()
        self.distribution = distribution
        self.status = "COMPLETED"
        self.paths = paths
        self.caller_ref = caller_ref

    @property
    def location(self) -> str:
        return self.distribution.location + f"/invalidation/{self.id}"

    @property
    def invalidation_batch(self) -> dict[str, Any]:
        return {
            "Paths": {"Quantity": len(self.paths), "Items": self.paths},
            "CallerReference": self.caller_ref,
        }


class PublicKey(BaseModel):
    def __init__(self, caller_ref: str, name: str, encoded_key: str):
        self.id = random_id(length=14)
        self.caller_ref = caller_ref
        self.name = name
        self.encoded_key = encoded_key
        self.created_time = utcnow()
        self.comment = ""
        self.etag = random_id(length=14)
        self.location = (
            f"https://cloudfront.amazonaws.com/2020-05-31/public-key/{self.id}"
        )

        # Last newline-separator is lost in the XML->Python transformation, but should exist
        if not self.encoded_key.endswith("\n"):
            self.encoded_key += "\n"

    @property
    def public_key_config(self) -> dict[str, str]:
        return {
            "CallerReference": self.caller_ref,
            "Name": self.name,
            "EncodedKey": self.encoded_key,
            "Comment": self.comment,
        }


class KeyGroup(BaseModel):
    def __init__(self, name: str, items: list[str]):
        self.id = random_id(length=14)
        self.name = name
        self.items = items
        self.etag = random_id(length=14)
        self.location = (
            f"https://cloudfront.amazonaws.com/2020-05-31/key-group/{self.id}"
        )

    @property
    def key_group_config(self) -> dict[str, Any]:
        return {
            "Items": self.items,
            "Name": self.name,
        }


class CloudFrontBackend(BaseBackend, TaggableResourcesMixin):
    SERVICE_NAMESPACE = "cloudfront"

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.distributions: dict[str, Distribution] = {}
        self.invalidations: dict[str, list[Invalidation]] = {}
        self.origin_access_controls: dict[str, OriginAccessControl] = {}
        self.public_keys: dict[str, PublicKey] = {}
        self.key_groups: dict[str, KeyGroup] = {}
        self.tagger = TaggingService()

    def create_distribution(
        self, distribution_config: dict[str, Any], tags: list[dict[str, str]]
    ) -> tuple[Distribution, str, str]:
        """
        Not all configuration options are supported yet.  Please raise an issue if
        we're not persisting/returning the correct attributes for your
        use-case.
        """
        # We'll always call dist_with_tags, as the incoming request is the same
        return self.create_distribution_with_tags(distribution_config, tags)

    def create_distribution_with_tags(
        self, distribution_config: dict[str, Any], tags: list[dict[str, str]]
    ) -> tuple[Distribution, str, str]:
        dist = Distribution(self.account_id, self.region_name, distribution_config)
        caller_reference = dist.distribution_config.caller_reference
        existing_dist = self._distribution_with_caller_reference(caller_reference)
        if existing_dist is not None:
            raise DistributionAlreadyExists(existing_dist.distribution_id)
        self.distributions[dist.distribution_id] = dist
        self.tagger.tag_resource(dist.arn, tags)
        return dist, dist.location, dist.etag

    def get_distribution(self, distribution_id: str) -> tuple[Distribution, str]:
        if distribution_id not in self.distributions:
            raise NoSuchDistribution
        dist = self.distributions[distribution_id]
        dist.advance()
        return dist, dist.etag

    def get_distribution_config(
        self, distribution_id: str
    ) -> tuple[DistributionConfig, str]:
        if distribution_id not in self.distributions:
            raise NoSuchDistribution
        dist = self.distributions[distribution_id]
        dist.advance()
        return dist.distribution_config, dist.etag

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

    def list_distributions(self) -> list[Distribution]:
        """
        Pagination is not supported yet.
        """
        for dist in self.distributions.values():
            dist.advance()
        return list(self.distributions.values())

    def _distribution_with_caller_reference(
        self, reference: str
    ) -> Distribution | None:
        for dist in self.distributions.values():
            config = dist.distribution_config
            if config.caller_reference == reference:
                return dist
        return None

    def update_distribution(
        self, dist_config: dict[str, Any], _id: str, if_match: bool
    ) -> tuple[Distribution, str, str]:
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
        self, dist_id: str, paths: list[str], caller_ref: str
    ) -> Invalidation:
        dist, _ = self.get_distribution(dist_id)
        invalidation = Invalidation(dist, paths, caller_ref)
        try:
            self.invalidations[dist_id].append(invalidation)
        except KeyError:
            self.invalidations[dist_id] = [invalidation]

        return invalidation

    def list_invalidations(self, dist_id: str) -> list[Invalidation]:
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
                    if invalidation.id == id:
                        return invalidation
        except KeyError:
            pass
        raise NoSuchInvalidation

    def list_tags_for_resource(self, resource: str) -> dict[str, list[dict[str, str]]]:
        return self.tagger.list_tags_for_resource(resource)

    def create_origin_access_control(
        self, config_dict: dict[str, str]
    ) -> OriginAccessControl:
        control = OriginAccessControl(config_dict)
        self.origin_access_controls[control.id] = control
        return control

    def get_origin_access_control(self, control_id: str) -> OriginAccessControl:
        if control_id not in self.origin_access_controls:
            raise NoSuchOriginAccessControl
        return self.origin_access_controls[control_id]

    def update_origin_access_control(
        self, control_id: str, config: dict[str, str]
    ) -> OriginAccessControl:
        """
        The IfMatch-parameter is not yet implemented
        """
        control = self.get_origin_access_control(control_id)
        control.update(config)
        return control

    def list_origin_access_controls(self) -> list[OriginAccessControl]:
        """
        Pagination is not yet implemented
        """
        return list(self.origin_access_controls.values())

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

    def list_public_keys(self) -> list[PublicKey]:
        """
        Pagination is not yet implemented
        """
        return list(self.public_keys.values())

    def create_key_group(self, name: str, items: list[str]) -> KeyGroup:
        key_group = KeyGroup(name=name, items=items)
        self.key_groups[key_group.id] = key_group
        return key_group

    def get_key_group(self, group_id: str) -> KeyGroup:
        return self.key_groups[group_id]

    def list_key_groups(self) -> list[KeyGroup]:
        """
        Pagination is not yet implemented
        """
        return list(self.key_groups.values())

    # Resource Groups Tagging API (TaggableResourcesMixin method overrides)
    def iter_tagged_resources(self) -> Iterator[TaggedResource]:
        for dist in self.distributions.values():
            yield TaggedResource(
                arn=dist.arn,
                tags=self.tagger.get_tag_dict_for_resource(dist.arn),
                resource_type="cloudfront:distribution",
            )

    def tag_resource(self, arn: str, tags: dict[str, str]) -> None:
        self.tagger.tag_resource(arn, self.tagger.convert_dict_to_tags_input(tags))

    def untag_resource(self, arn: str, tag_keys: list[str]) -> None:
        self.tagger.untag_resource_using_names(arn, tag_keys)


cloudfront_backends = BackendDict(
    CloudFrontBackend,
    "cloudfront",
    use_boto3_regions=False,
    additional_regions=PARTITION_NAMES,
)
