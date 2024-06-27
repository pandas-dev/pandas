from typing import Any, List, Optional

from moto.core.exceptions import JsonRESTError


class NameTooLongException(JsonRESTError):
    code = 400

    def __init__(self, name: str, location: str, max_limit: int = 256):
        message = (
            f"1 validation error detected: Value '{name}' at '{location}' "
            f"failed to satisfy constraint: Member must have length less "
            f"than or equal to {max_limit}"
        )
        super().__init__("ValidationException", message)


class InvalidConfigurationRecorderNameException(JsonRESTError):
    code = 400

    def __init__(self, name: Optional[str]):
        message = (
            f"The configuration recorder name '{name}' is not valid, blank string."
        )
        super().__init__("InvalidConfigurationRecorderNameException", message)


class MaxNumberOfConfigurationRecordersExceededException(JsonRESTError):
    code = 400

    def __init__(self, name: str):
        message = (
            f"Failed to put configuration recorder '{name}' because the maximum number of "
            "configuration recorders: 1 is reached."
        )
        super().__init__("MaxNumberOfConfigurationRecordersExceededException", message)


class InvalidRecordingGroupException(JsonRESTError):
    code = 400

    def __init__(self) -> None:
        message = "The recording group provided is not valid"
        super().__init__("InvalidRecordingGroupException", message)


class InvalidResourceTypeException(JsonRESTError):
    code = 400

    def __init__(self, bad_list: List[str], good_list: Any):
        message = (
            f"{len(bad_list)} validation error detected: Value '{bad_list}' at "
            "'configurationRecorder.recordingGroup.resourceTypes' failed to satisfy constraint: "
            f"Member must satisfy constraint: [Member must satisfy enum value set: {good_list}]"
        )

        super().__init__("ValidationException", message)


class NoSuchConfigurationAggregatorException(JsonRESTError):
    code = 400

    def __init__(self, number: int = 1):
        if number == 1:
            message = "The configuration aggregator does not exist. Check the configuration aggregator name and try again."
        else:
            message = (
                "At least one of the configuration aggregators does not exist. Check the configuration aggregator"
                " names and try again."
            )
        super().__init__("NoSuchConfigurationAggregatorException", message)


class NoSuchConfigurationRecorderException(JsonRESTError):
    code = 400

    def __init__(self, name: str):
        message = (
            f"Cannot find configuration recorder with the specified name '{name}'."
        )
        super().__init__("NoSuchConfigurationRecorderException", message)


class InvalidDeliveryChannelNameException(JsonRESTError):
    code = 400

    def __init__(self, name: Optional[str]):
        message = f"The delivery channel name '{name}' is not valid, blank string."
        super().__init__("InvalidDeliveryChannelNameException", message)


class NoSuchBucketException(JsonRESTError):
    """We are *only* validating that there is value that is not '' here."""

    code = 400

    def __init__(self) -> None:
        message = "Cannot find a S3 bucket with an empty bucket name."
        super().__init__("NoSuchBucketException", message)


class InvalidNextTokenException(JsonRESTError):
    code = 400

    def __init__(self) -> None:
        message = "The nextToken provided is invalid"
        super().__init__("InvalidNextTokenException", message)


class InvalidS3KeyPrefixException(JsonRESTError):
    code = 400

    def __init__(self) -> None:
        message = "The s3 key prefix '' is not valid, empty s3 key prefix."
        super().__init__("InvalidS3KeyPrefixException", message)


class InvalidS3KmsKeyArnException(JsonRESTError):
    code = 400

    def __init__(self) -> None:
        message = "The arn '' is not a valid kms key or alias arn."
        super().__init__("InvalidS3KmsKeyArnException", message)


class InvalidSNSTopicARNException(JsonRESTError):
    """We are *only* validating that there is value that is not '' here."""

    code = 400

    def __init__(self) -> None:
        message = "The sns topic arn '' is not valid."
        super().__init__("InvalidSNSTopicARNException", message)


class InvalidDeliveryFrequency(JsonRESTError):
    code = 400

    def __init__(self, value: str, good_list: Any):
        message = (
            f"1 validation error detected: Value '{value}' at "
            "'deliveryChannel.configSnapshotDeliveryProperties.deliveryFrequency' failed to satisfy "
            f"constraint: Member must satisfy enum value set: {good_list}"
        )
        super().__init__("InvalidDeliveryFrequency", message)


class MaxNumberOfDeliveryChannelsExceededException(JsonRESTError):
    code = 400

    def __init__(self, name: str):
        message = f"Failed to put delivery channel '{name}' because the maximum number of delivery channels: 1 is reached."
        super().__init__("MaxNumberOfDeliveryChannelsExceededException", message)


class NoSuchDeliveryChannelException(JsonRESTError):
    code = 400

    def __init__(self, name: str):
        message = f"Cannot find delivery channel with specified name '{name}'."
        super().__init__("NoSuchDeliveryChannelException", message)


class NoAvailableConfigurationRecorderException(JsonRESTError):
    code = 400

    def __init__(self) -> None:
        message = "Configuration recorder is not available to put delivery channel."
        super().__init__("NoAvailableConfigurationRecorderException", message)


class NoAvailableDeliveryChannelException(JsonRESTError):
    code = 400

    def __init__(self) -> None:
        message = "Delivery channel is not available to start configuration recorder."
        super().__init__("NoAvailableDeliveryChannelException", message)


class LastDeliveryChannelDeleteFailedException(JsonRESTError):
    code = 400

    def __init__(self, name: str):
        message = (
            f"Failed to delete last specified delivery channel with name '{name}', because there, "
            "because there is a running configuration recorder."
        )
        super().__init__("LastDeliveryChannelDeleteFailedException", message)


class TooManyAccountSources(JsonRESTError):
    code = 400

    def __init__(self, length: int):
        locations = ["com.amazonaws.xyz"] * length

        locs = ", ".join(locations)
        message = (
            f"Value '[{locs}]' at 'accountAggregationSources' failed to satisfy constraint: "
            "Member must have length less than or equal to 1"
        )
        super().__init__("ValidationException", message)


class DuplicateTags(JsonRESTError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "InvalidInput",
            "Duplicate tag keys found. Please note that Tag keys are case insensitive.",
        )


class TagKeyTooBig(JsonRESTError):
    code = 400

    def __init__(self, tag: str, param: str = "tags.X.member.key"):
        super().__init__(
            "ValidationException",
            f"1 validation error detected: Value '{tag}' at '{param}' failed to satisfy "
            "constraint: Member must have length less than or equal to 128",
        )


class TagValueTooBig(JsonRESTError):
    code = 400

    def __init__(self, tag: str, param: str = "tags.X.member.value"):
        super().__init__(
            "ValidationException",
            f"1 validation error detected: Value '{tag}' at '{param}' failed to satisfy "
            "constraint: Member must have length less than or equal to 256",
        )


class InvalidParameterValueException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("InvalidParameterValueException", message)


class InvalidTagCharacters(JsonRESTError):
    code = 400

    def __init__(self, tag: str, param: str = "tags.X.member.key"):
        message = f"1 validation error detected: Value '{tag}' at '{param}' failed to satisfy "
        message += "constraint: Member must satisfy regular expression pattern: [\\\\p{L}\\\\p{Z}\\\\p{N}_.:/=+\\\\-@]+"

        super().__init__("ValidationException", message)


class TooManyTags(JsonRESTError):
    code = 400

    def __init__(self, tags: Any, param: str = "tags"):
        super().__init__(
            "ValidationException",
            f"1 validation error detected: Value '{tags}' at '{param}' failed to satisfy "
            "constraint: Member must have length less than or equal to 50.",
        )


class InvalidResourceParameters(JsonRESTError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "ValidationException",
            "Both Resource ID and Resource Name " "cannot be specified in the request",
        )


class InvalidLimitException(JsonRESTError):
    code = 400

    def __init__(self, value: int):
        super().__init__(
            "InvalidLimitException",
            f"Value '{value}' at 'limit' failed to satisfy constraint: Member"
            " must have value less than or equal to 100",
        )


class TooManyResourceIds(JsonRESTError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "ValidationException",
            "The specified list had more than 20 resource ID's. "
            "It must have '20' or less items",
        )


class ResourceNotDiscoveredException(JsonRESTError):
    code = 400

    def __init__(self, resource_type: str, resource: str):
        super().__init__(
            "ResourceNotDiscoveredException",
            f"Resource {resource} of resourceType:{resource_type} is unknown or has not been discovered",
        )


class ResourceNotFoundException(JsonRESTError):
    code = 400

    def __init__(self, resource_arn: str):
        super().__init__(
            "ResourceNotFoundException", f"ResourceArn '{resource_arn}' does not exist"
        )


class TooManyResourceKeys(JsonRESTError):
    code = 400

    def __init__(self, bad_list: List[str]):
        message = (
            f"1 validation error detected: Value '{bad_list}' at "
            "'resourceKeys' failed to satisfy constraint: "
            "Member must have length less than or equal to 100"
        )
        super().__init__("ValidationException", message)


class InvalidResultTokenException(JsonRESTError):
    code = 400

    def __init__(self) -> None:
        message = "The resultToken provided is invalid"
        super().__init__("InvalidResultTokenException", message)


class ValidationException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("ValidationException", message)


class NoSuchOrganizationConformancePackException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("NoSuchOrganizationConformancePackException", message)


class MaxNumberOfConfigRulesExceededException(JsonRESTError):
    code = 400

    def __init__(self, name: str, max_limit: int):
        message = f"Failed to put config rule '{name}' because the maximum number of config rules: {max_limit} is reached."
        super().__init__("MaxNumberOfConfigRulesExceededException", message)


class ResourceInUseException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("ResourceInUseException", message)


class InsufficientPermissionsException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("InsufficientPermissionsException", message)


class NoSuchConfigRuleException(JsonRESTError):
    code = 400

    def __init__(self, rule_name: str):
        message = f"The ConfigRule '{rule_name}' provided in the request is invalid. Please check the configRule name"
        super().__init__("NoSuchConfigRuleException", message)


class MissingRequiredConfigRuleParameterException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("ParamValidationError", message)


class NoSuchRetentionConfigurationException(JsonRESTError):
    code = 400

    def __init__(self, name: str):
        message = (
            f"Cannot find retention configuration with the specified name '{name}'."
        )
        super().__init__("NoSuchRetentionConfigurationException", message)
