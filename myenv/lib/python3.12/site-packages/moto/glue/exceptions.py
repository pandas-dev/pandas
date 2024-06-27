from typing import Optional

from moto.core.exceptions import JsonRESTError


class GlueClientError(JsonRESTError):
    code = 400


class AlreadyExistsException(GlueClientError):
    def __init__(self, typ: str):
        super().__init__("AlreadyExistsException", f"{typ} already exists.")


class DatabaseAlreadyExistsException(AlreadyExistsException):
    def __init__(self) -> None:
        super().__init__("Database")


class TableAlreadyExistsException(AlreadyExistsException):
    def __init__(self) -> None:
        super().__init__("Table")


class PartitionAlreadyExistsException(AlreadyExistsException):
    def __init__(self) -> None:
        super().__init__("Partition")


class CrawlerAlreadyExistsException(AlreadyExistsException):
    def __init__(self) -> None:
        super().__init__("Crawler")


class SessionAlreadyExistsException(AlreadyExistsException):
    def __init__(self) -> None:
        super().__init__("Session")


class EntityNotFoundException(GlueClientError):
    def __init__(self, msg: str):
        super().__init__("EntityNotFoundException", msg)


class DatabaseNotFoundException(EntityNotFoundException):
    def __init__(self, db: str):
        super().__init__(f"Database {db} not found.")


class TableNotFoundException(EntityNotFoundException):
    def __init__(self, tbl: str):
        super().__init__(f"Table {tbl} not found.")


class PartitionNotFoundException(EntityNotFoundException):
    def __init__(self) -> None:
        super().__init__("Cannot find partition.")


class CrawlerNotFoundException(EntityNotFoundException):
    def __init__(self, crawler: str):
        super().__init__(f"Crawler {crawler} not found.")


class JobNotFoundException(EntityNotFoundException):
    def __init__(self, job: str):
        super().__init__(f"Job {job} not found.")


class JobRunNotFoundException(EntityNotFoundException):
    def __init__(self, job_run: str):
        super().__init__(f"Job run {job_run} not found.")


class VersionNotFoundException(EntityNotFoundException):
    def __init__(self) -> None:
        super().__init__("Version not found.")


class SchemaNotFoundException(EntityNotFoundException):
    def __init__(
        self,
        schema_name: str,
        registry_name: str,
        schema_arn: Optional[str],
        null: str = "null",
    ):
        super().__init__(
            f"Schema is not found. RegistryName: {registry_name if registry_name else null}, SchemaName: {schema_name if schema_name else null}, SchemaArn: {schema_arn if schema_arn else null}",
        )


class SchemaVersionNotFoundFromSchemaIdException(EntityNotFoundException):
    def __init__(
        self,
        registry_name: Optional[str],
        schema_name: Optional[str],
        schema_arn: Optional[str],
        version_number: Optional[str],
        latest_version: Optional[str],
        null: str = "null",
        false: str = "false",
    ):
        super().__init__(
            f"Schema version is not found. RegistryName: {registry_name if registry_name else null}, SchemaName: {schema_name if schema_name else null}, SchemaArn: {schema_arn if schema_arn else null}, VersionNumber: {version_number if version_number else null}, isLatestVersion: {latest_version if latest_version else false}",
        )


class SchemaVersionNotFoundFromSchemaVersionIdException(EntityNotFoundException):
    def __init__(self, schema_version_id: str):
        super().__init__(
            f"Schema version is not found. SchemaVersionId: {schema_version_id}",
        )


class SessionNotFoundException(EntityNotFoundException):
    def __init__(self, session: str):
        super().__init__(f"Session {session} not found.")


class IllegalSessionStateException(GlueClientError):
    def __init__(self, msg: str):
        super().__init__("IllegalSessionStateException", msg)


class RegistryNotFoundException(EntityNotFoundException):
    def __init__(self, resource: str, param_name: str, param_value: Optional[str]):
        super().__init__(
            resource + " is not found. " + param_name + ": " + param_value,  # type: ignore
        )


class TriggerNotFoundException(EntityNotFoundException):
    def __init__(self, trigger: str):
        super().__init__(f"Trigger {trigger} not found.")


class CrawlerRunningException(GlueClientError):
    def __init__(self, msg: str):
        super().__init__("CrawlerRunningException", msg)


class CrawlerNotRunningException(GlueClientError):
    def __init__(self, msg: str):
        super().__init__("CrawlerNotRunningException", msg)


class ConcurrentRunsExceededException(GlueClientError):
    def __init__(self, msg: str):
        super().__init__("ConcurrentRunsExceededException", msg)


class ResourceNumberLimitExceededException(GlueClientError):
    def __init__(self, msg: str):
        super().__init__(
            "ResourceNumberLimitExceededException",
            msg,
        )


class GeneralResourceNumberLimitExceededException(ResourceNumberLimitExceededException):
    def __init__(self, resource: str):
        super().__init__(
            "More "
            + resource
            + " cannot be created. The maximum limit has been reached.",
        )


class SchemaVersionMetadataLimitExceededException(ResourceNumberLimitExceededException):
    def __init__(self) -> None:
        super().__init__(
            "Your resource limits for Schema Version Metadata have been exceeded.",
        )


class GSRAlreadyExistsException(GlueClientError):
    def __init__(self, msg: str):
        super().__init__(
            "AlreadyExistsException",
            msg,
        )


class SchemaVersionMetadataAlreadyExistsException(GSRAlreadyExistsException):
    def __init__(self, schema_version_id: str, metadata_key: str, metadata_value: str):
        super().__init__(
            f"Resource already exist for schema version id: {schema_version_id}, metadata key: {metadata_key}, metadata value: {metadata_value}",
        )


class GeneralGSRAlreadyExistsException(GSRAlreadyExistsException):
    def __init__(self, resource: str, param_name: str, param_value: str):
        super().__init__(
            resource + " already exists. " + param_name + ": " + param_value,
        )


class _InvalidOperationException(GlueClientError):
    def __init__(self, error_type: str, op: str, msg: str):
        super().__init__(
            error_type,
            "An error occurred (%s) when calling the %s operation: %s"
            % (error_type, op, msg),
        )


class InvalidStateException(_InvalidOperationException):
    def __init__(self, op: str, msg: str):
        super().__init__("InvalidStateException", op, msg)


class InvalidInputException(_InvalidOperationException):
    def __init__(self, op: str, msg: str):
        super().__init__("InvalidInputException", op, msg)


class GSRInvalidInputException(GlueClientError):
    def __init__(self, msg: str):
        super().__init__("InvalidInputException", msg)


class ResourceNameTooLongException(GSRInvalidInputException):
    def __init__(self, param_name: str):
        super().__init__(
            "The resource name contains too many or too few characters. Parameter Name: "
            + param_name,
        )


class ParamValueContainsInvalidCharactersException(GSRInvalidInputException):
    def __init__(self, param_name: str):
        super().__init__(
            "The parameter value contains one or more characters that are not valid. Parameter Name: "
            + param_name,
        )


class InvalidNumberOfTagsException(GSRInvalidInputException):
    def __init__(self) -> None:
        super().__init__(
            "New Tags cannot be empty or more than 50",
        )


class InvalidDataFormatException(GSRInvalidInputException):
    def __init__(self) -> None:
        super().__init__(
            "Data format is not valid.",
        )


class InvalidCompatibilityException(GSRInvalidInputException):
    def __init__(self) -> None:
        super().__init__(
            "Compatibility is not valid.",
        )


class InvalidSchemaDefinitionException(GSRInvalidInputException):
    def __init__(self, data_format_name: str, err: ValueError):
        super().__init__(
            "Schema definition of "
            + data_format_name
            + " data format is invalid: "
            + str(err),
        )


class InvalidRegistryIdBothParamsProvidedException(GSRInvalidInputException):
    def __init__(self) -> None:
        super().__init__(
            "One of registryName or registryArn has to be provided, both cannot be provided.",
        )


class InvalidSchemaIdBothParamsProvidedException(GSRInvalidInputException):
    def __init__(self) -> None:
        super().__init__(
            "One of (registryName and schemaName) or schemaArn has to be provided, both cannot be provided.",
        )


class InvalidSchemaIdNotProvidedException(GSRInvalidInputException):
    def __init__(self) -> None:
        super().__init__(
            "At least one of (registryName and schemaName) or schemaArn has to be provided.",
        )


class InvalidSchemaVersionNumberBothParamsProvidedException(GSRInvalidInputException):
    def __init__(self) -> None:
        super().__init__("Only one of VersionNumber or LatestVersion is required.")


class InvalidSchemaVersionNumberNotProvidedException(GSRInvalidInputException):
    def __init__(self) -> None:
        super().__init__("One of version number (or) latest version is required.")


class InvalidSchemaVersionIdProvidedWithOtherParamsException(GSRInvalidInputException):
    def __init__(self) -> None:
        super().__init__(
            "No other input parameters can be specified when fetching by SchemaVersionId."
        )


class DisabledCompatibilityVersioningException(GSRInvalidInputException):
    def __init__(
        self,
        schema_name: str,
        registry_name: str,
        schema_arn: Optional[str],
        null: str = "null",
    ):
        super().__init__(
            f"Compatibility DISABLED does not allow versioning. SchemaId: SchemaId(schemaArn={schema_arn if schema_arn else null}, schemaName={schema_name if schema_name else null}, registryName={registry_name if registry_name else null})"
        )
