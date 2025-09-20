import json
import re
from typing import Any, Dict, Optional, Pattern, Tuple

from .exceptions import (
    DisabledCompatibilityVersioningException,
    GeneralGSRAlreadyExistsException,
    GeneralResourceNumberLimitExceededException,
    InvalidCompatibilityException,
    InvalidDataFormatException,
    InvalidNumberOfTagsException,
    InvalidRegistryIdBothParamsProvidedException,
    InvalidSchemaDefinitionException,
    InvalidSchemaIdBothParamsProvidedException,
    InvalidSchemaIdNotProvidedException,
    InvalidSchemaVersionIdProvidedWithOtherParamsException,
    InvalidSchemaVersionNumberBothParamsProvidedException,
    InvalidSchemaVersionNumberNotProvidedException,
    ParamValueContainsInvalidCharactersException,
    RegistryNotFoundException,
    ResourceNameTooLongException,
    SchemaNotFoundException,
    SchemaVersionMetadataLimitExceededException,
)
from .glue_schema_registry_constants import (
    ARN_PATTERN,
    DEFAULT_REGISTRY_NAME,
    DESCRIPTION_PATTERN,
    LATEST_VERSION,
    MAX_ARN_LENGTH,
    MAX_DESCRIPTION_LENGTH,
    MAX_REGISTRIES_ALLOWED,
    MAX_REGISTRY_NAME_LENGTH,
    MAX_SCHEMA_DEFINITION_LENGTH,
    MAX_SCHEMA_NAME_LENGTH,
    MAX_SCHEMA_VERSION_METADATA_ALLOWED,
    MAX_SCHEMA_VERSION_METADATA_LENGTH,
    MAX_SCHEMA_VERSIONS_ALLOWED,
    MAX_SCHEMAS_ALLOWED,
    MAX_TAGS_ALLOWED,
    METADATA_KEY,
    METADATA_VALUE,
    REGISTRY_ARN,
    REGISTRY_NAME,
    RESOURCE_NAME_PATTERN,
    SCHEMA_ARN,
    SCHEMA_DEFINITION,
    SCHEMA_NAME,
    SCHEMA_VERSION_ID,
    SCHEMA_VERSION_ID_PATTERN,
    SCHEMA_VERSION_METADATA_PATTERN,
    VERSION_NUMBER,
)


def validate_registry_name_pattern_and_length(param_value: str) -> None:
    validate_param_pattern_and_length(
        param_value,
        param_name="registryName",
        max_name_length=MAX_REGISTRY_NAME_LENGTH,
        pattern=RESOURCE_NAME_PATTERN,
    )


def validate_arn_pattern_and_length(param_value: str) -> None:
    validate_param_pattern_and_length(
        param_value,
        param_name="registryArn",
        max_name_length=MAX_ARN_LENGTH,
        pattern=ARN_PATTERN,
    )


def validate_description_pattern_and_length(param_value: str) -> None:
    validate_param_pattern_and_length(
        param_value,
        param_name="description",
        max_name_length=MAX_DESCRIPTION_LENGTH,
        pattern=DESCRIPTION_PATTERN,
    )


def validate_schema_name_pattern_and_length(param_value: str) -> None:
    validate_param_pattern_and_length(
        param_value,
        param_name="schemaName",
        max_name_length=MAX_SCHEMA_NAME_LENGTH,
        pattern=RESOURCE_NAME_PATTERN,
    )


def validate_schema_version_metadata_key_pattern_and_length(param_value: str) -> None:
    validate_param_pattern_and_length(
        param_value,
        param_name="key",
        max_name_length=MAX_SCHEMA_VERSION_METADATA_LENGTH,
        pattern=SCHEMA_VERSION_METADATA_PATTERN,
    )


def validate_schema_version_metadata_value_pattern_and_length(param_value: str) -> None:
    validate_param_pattern_and_length(
        param_value,
        param_name="value",
        max_name_length=MAX_SCHEMA_VERSION_METADATA_LENGTH,
        pattern=SCHEMA_VERSION_METADATA_PATTERN,
    )


def validate_param_pattern_and_length(
    param_value: str, param_name: str, max_name_length: int, pattern: Pattern[str]
) -> None:
    if len(param_value.encode("utf-8")) > max_name_length:
        raise ResourceNameTooLongException(param_name)

    if re.match(pattern, param_value) is None:
        raise ParamValueContainsInvalidCharactersException(param_name)


def validate_schema_definition(schema_definition: str, data_format: str) -> None:
    validate_schema_definition_length(schema_definition)
    if data_format in ["AVRO", "JSON"]:
        try:
            json.loads(schema_definition)
        except ValueError as err:
            raise InvalidSchemaDefinitionException(data_format, err)


def validate_schema_definition_length(schema_definition: str) -> None:
    if len(schema_definition) > MAX_SCHEMA_DEFINITION_LENGTH:
        param_name = SCHEMA_DEFINITION
        raise ResourceNameTooLongException(param_name)


def validate_schema_version_id_pattern(schema_version_id: str) -> None:
    if re.match(SCHEMA_VERSION_ID_PATTERN, schema_version_id) is None:
        raise ParamValueContainsInvalidCharactersException(SCHEMA_VERSION_ID)


def validate_number_of_tags(tags: Dict[str, str]) -> None:
    if len(tags) > MAX_TAGS_ALLOWED:
        raise InvalidNumberOfTagsException()


def validate_registry_id(
    registry_id: Dict[str, Any], registries: Dict[str, Any]
) -> str:
    if not registry_id:
        return DEFAULT_REGISTRY_NAME

    if registry_id.get(REGISTRY_NAME) and registry_id.get(REGISTRY_ARN):
        raise InvalidRegistryIdBothParamsProvidedException()

    if registry_id.get(REGISTRY_NAME):
        registry_name = registry_id.get(REGISTRY_NAME)
        validate_registry_name_pattern_and_length(registry_name)  # type: ignore

    elif registry_id.get(REGISTRY_ARN):
        registry_arn = registry_id.get(REGISTRY_ARN)
        validate_arn_pattern_and_length(registry_arn)  # type: ignore
        registry_name = registry_arn.split("/")[-1]  # type: ignore

    if registry_name != DEFAULT_REGISTRY_NAME and registry_name not in registries:
        if registry_id.get(REGISTRY_NAME):
            raise RegistryNotFoundException(
                resource="Registry",
                param_name=REGISTRY_NAME,
                param_value=registry_name,
            )
        if registry_id.get(REGISTRY_ARN):
            raise RegistryNotFoundException(
                resource="Registry",
                param_name=REGISTRY_ARN,
                param_value=registry_arn,
            )

    return registry_name  # type: ignore


def validate_registry_params(
    registries: Any,
    registry_name: str,
    description: Optional[str] = None,
    tags: Optional[Dict[str, str]] = None,
) -> None:
    validate_registry_name_pattern_and_length(registry_name)

    if description:
        validate_description_pattern_and_length(description)

    if tags:
        validate_number_of_tags(tags)

    if len(registries) >= MAX_REGISTRIES_ALLOWED:
        raise GeneralResourceNumberLimitExceededException(resource="registries")

    if registry_name in registries:
        raise GeneralGSRAlreadyExistsException(
            resource="Registry",
            param_name=REGISTRY_NAME,
            param_value=registry_name,
        )


def validate_schema_id(
    schema_id: Dict[str, str], registries: Dict[str, Any]
) -> Tuple[str, str, Optional[str]]:
    schema_arn = schema_id.get(SCHEMA_ARN)
    registry_name = schema_id.get(REGISTRY_NAME)
    schema_name = schema_id.get(SCHEMA_NAME)
    if schema_arn:
        if registry_name or schema_name:
            raise InvalidSchemaIdBothParamsProvidedException()
        validate_arn_pattern_and_length(schema_arn)
        arn_components = schema_arn.split("/")
        schema_name = arn_components[-1]
        registry_name = arn_components[-2]

    else:
        if registry_name is None or schema_name is None:
            raise InvalidSchemaIdNotProvidedException()
        validate_registry_name_pattern_and_length(registry_name)
        validate_schema_name_pattern_and_length(schema_name)

    if (
        registry_name not in registries
        or schema_name not in registries[registry_name].schemas
    ):
        raise SchemaNotFoundException(schema_name, registry_name, schema_arn)

    return registry_name, schema_name, schema_arn


def validate_schema_params(
    registry: Any,
    schema_name: str,
    data_format: str,
    compatibility: str,
    schema_definition: str,
    num_schemas: int,
    description: Optional[str] = None,
    tags: Optional[Dict[str, str]] = None,
) -> None:
    validate_schema_name_pattern_and_length(schema_name)

    if data_format not in ["AVRO", "JSON", "PROTOBUF"]:
        raise InvalidDataFormatException()

    if compatibility not in [
        "NONE",
        "DISABLED",
        "BACKWARD",
        "BACKWARD_ALL",
        "FORWARD",
        "FORWARD_ALL",
        "FULL",
        "FULL_ALL",
    ]:
        raise InvalidCompatibilityException()

    if description:
        validate_description_pattern_and_length(description)

    if tags:
        validate_number_of_tags(tags)

    validate_schema_definition(schema_definition, data_format)

    if num_schemas >= MAX_SCHEMAS_ALLOWED:
        raise GeneralResourceNumberLimitExceededException(resource="schemas")

    if schema_name in registry.schemas:
        raise GeneralGSRAlreadyExistsException(
            resource="Schema",
            param_name=SCHEMA_NAME,
            param_value=schema_name,
        )


def validate_register_schema_version_params(
    registry_name: str,
    schema_name: str,
    schema_arn: Optional[str],
    num_schema_versions: int,
    schema_definition: str,
    compatibility: str,
    data_format: str,
) -> None:
    if compatibility == "DISABLED":
        raise DisabledCompatibilityVersioningException(
            schema_name, registry_name, schema_arn
        )

    validate_schema_definition(schema_definition, data_format)

    if num_schema_versions >= MAX_SCHEMA_VERSIONS_ALLOWED:
        raise GeneralResourceNumberLimitExceededException(resource="schema versions")


def validate_schema_version_params(  # type: ignore[return]
    registries: Dict[str, Any],
    schema_id: Optional[Dict[str, Any]],
    schema_version_id: Optional[str],
    schema_version_number: Optional[Dict[str, Any]],
) -> Tuple[
    Optional[str],
    Optional[str],
    Optional[str],
    Optional[str],
    Optional[str],
    Optional[str],
]:
    if not schema_version_id and not schema_id and not schema_version_number:
        raise InvalidSchemaIdNotProvidedException()

    if schema_version_id and (schema_id or schema_version_number):
        raise InvalidSchemaVersionIdProvidedWithOtherParamsException()

    if schema_version_id:
        validate_schema_version_id_pattern(schema_version_id)

        # returns schema_version_id, registry_name, schema_name, schema_arn, version_number, latest_version
        return schema_version_id, None, None, None, None, None

    if schema_id and schema_version_number:
        registry_name, schema_name, schema_arn = validate_schema_id(
            schema_id, registries
        )
        version_number, latest_version = validate_schema_version_number(
            registries, registry_name, schema_name, schema_version_number
        )
        return (
            None,
            registry_name,
            schema_name,
            schema_arn,
            version_number,
            latest_version,
        )

    if not schema_id:
        raise InvalidSchemaIdNotProvidedException()

    if not schema_version_number:
        raise InvalidSchemaVersionNumberNotProvidedException()


def validate_schema_version_number(
    registries: Dict[str, Any],
    registry_name: str,
    schema_name: str,
    schema_version_number: Dict[str, str],
) -> Tuple[str, str]:
    latest_version = schema_version_number.get(LATEST_VERSION)
    version_number = schema_version_number.get(VERSION_NUMBER)
    schema = registries[registry_name].schemas[schema_name]
    if latest_version:
        if version_number:
            raise InvalidSchemaVersionNumberBothParamsProvidedException()
        return schema.latest_schema_version, latest_version

    return version_number, latest_version  # type: ignore


def validate_schema_version_metadata_pattern_and_length(
    metadata_key_value: Dict[str, str],
) -> Tuple[str, str]:
    metadata_key = metadata_key_value.get(METADATA_KEY)
    metadata_value = metadata_key_value.get(METADATA_VALUE)

    validate_schema_version_metadata_key_pattern_and_length(metadata_key)  # type: ignore
    validate_schema_version_metadata_value_pattern_and_length(metadata_value)  # type: ignore

    return metadata_key, metadata_value  # type: ignore[return-value]


def validate_number_of_schema_version_metadata_allowed(
    metadata: Dict[str, Any],
) -> None:
    num_metadata_key_value_pairs = 0
    for m in metadata.values():
        num_metadata_key_value_pairs += len(m)

    if num_metadata_key_value_pairs >= MAX_SCHEMA_VERSION_METADATA_ALLOWED:
        raise SchemaVersionMetadataLimitExceededException()


def get_schema_version_if_definition_exists(
    schema_versions: Any, data_format: str, schema_definition: str
) -> Optional[Dict[str, Any]]:
    if data_format in ["AVRO", "JSON"]:
        for schema_version in schema_versions:
            if json.loads(schema_definition) == json.loads(
                schema_version.schema_definition
            ):
                return schema_version.as_dict()
    else:
        for schema_version in schema_versions:
            if schema_definition == schema_version.schema_definition:
                return schema_version.as_dict()
    return None


def get_put_schema_version_metadata_response(
    schema_id: Dict[str, Any],
    schema_version_number: Optional[Dict[str, str]],
    schema_version_id: str,
    metadata_key_value: Dict[str, str],
) -> Dict[str, Any]:
    put_schema_version_metadata_response_dict: Dict[str, Any] = {}
    if schema_version_id:
        put_schema_version_metadata_response_dict[SCHEMA_VERSION_ID] = schema_version_id
    if schema_id:
        schema_arn = schema_id.get(SCHEMA_ARN)
        registry_name = schema_id.get(REGISTRY_NAME)
        schema_name = schema_id.get(SCHEMA_NAME)
        if schema_arn:
            put_schema_version_metadata_response_dict[SCHEMA_ARN] = schema_arn
        if registry_name:
            put_schema_version_metadata_response_dict[REGISTRY_NAME] = registry_name
        if schema_name:
            put_schema_version_metadata_response_dict[SCHEMA_NAME] = schema_name

    if schema_version_number:
        latest_version = schema_version_number.get(LATEST_VERSION)
        version_number = schema_version_number.get(VERSION_NUMBER)
        if latest_version:
            put_schema_version_metadata_response_dict[LATEST_VERSION] = latest_version
        else:
            put_schema_version_metadata_response_dict[LATEST_VERSION] = False

        if version_number:
            put_schema_version_metadata_response_dict[VERSION_NUMBER] = version_number
    else:
        put_schema_version_metadata_response_dict[LATEST_VERSION] = False
        put_schema_version_metadata_response_dict[VERSION_NUMBER] = 0

    metadata_key = metadata_key_value.get(METADATA_KEY)
    metadata_value = metadata_key_value.get(METADATA_VALUE)
    put_schema_version_metadata_response_dict[METADATA_KEY] = metadata_key
    put_schema_version_metadata_response_dict[METADATA_VALUE] = metadata_value

    return put_schema_version_metadata_response_dict


def delete_schema_response(
    schema_name: str, schema_arn: str, status: str
) -> Dict[str, Any]:
    return {
        "SchemaName": schema_name,
        "SchemaArn": schema_arn,
        "Status": status,
    }
