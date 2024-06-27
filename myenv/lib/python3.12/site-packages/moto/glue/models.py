import json
import re
import time
from collections import OrderedDict
from datetime import datetime
from typing import Any, Dict, List, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import unix_time, utcnow
from moto.moto_api._internal import mock_random
from moto.moto_api._internal.managed_state_model import ManagedState
from moto.utilities.utils import get_partition

from ..utilities.paginator import paginate
from ..utilities.tagging_service import TaggingService
from .exceptions import (
    ConcurrentRunsExceededException,
    CrawlerAlreadyExistsException,
    CrawlerNotFoundException,
    CrawlerNotRunningException,
    CrawlerRunningException,
    DatabaseAlreadyExistsException,
    DatabaseNotFoundException,
    IllegalSessionStateException,
    JobNotFoundException,
    JobRunNotFoundException,
    JsonRESTError,
    PartitionAlreadyExistsException,
    PartitionNotFoundException,
    SchemaNotFoundException,
    SchemaVersionMetadataAlreadyExistsException,
    SchemaVersionNotFoundFromSchemaIdException,
    SchemaVersionNotFoundFromSchemaVersionIdException,
    SessionAlreadyExistsException,
    SessionNotFoundException,
    TableAlreadyExistsException,
    TableNotFoundException,
    TriggerNotFoundException,
    VersionNotFoundException,
)
from .glue_schema_registry_constants import (
    AVAILABLE_STATUS,
    DEFAULT_REGISTRY_NAME,
    DELETING_STATUS,
)
from .glue_schema_registry_utils import (
    delete_schema_response,
    get_put_schema_version_metadata_response,
    get_schema_version_if_definition_exists,
    validate_number_of_schema_version_metadata_allowed,
    validate_register_schema_version_params,
    validate_registry_id,
    validate_registry_params,
    validate_schema_definition_length,
    validate_schema_id,
    validate_schema_params,
    validate_schema_version_metadata_pattern_and_length,
    validate_schema_version_params,
)
from .utils import PartitionFilter


class GlueBackend(BaseBackend):
    PAGINATION_MODEL = {
        "get_jobs": {
            "input_token": "next_token",
            "limit_key": "max_results",
            "limit_default": 100,
            "unique_attribute": "name",
        },
        "get_triggers": {
            "input_token": "next_token",
            "limit_key": "max_results",
            "limit_default": 100,
            "unique_attribute": "name",
        },
        "list_crawlers": {
            "input_token": "next_token",
            "limit_key": "max_results",
            "limit_default": 100,
            "unique_attribute": "name",
        },
        "list_jobs": {
            "input_token": "next_token",
            "limit_key": "max_results",
            "limit_default": 100,
            "unique_attribute": "name",
        },
        "list_sessions": {
            "input_token": "next_token",
            "limit_key": "max_results",
            "limit_default": 100,
            "unique_attribute": "session_id",
        },
        "list_triggers": {
            "input_token": "next_token",
            "limit_key": "max_results",
            "limit_default": 100,
            "unique_attribute": "name",
        },
    }

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.databases: Dict[str, FakeDatabase] = OrderedDict()
        self.crawlers: Dict[str, FakeCrawler] = OrderedDict()
        self.jobs: Dict[str, FakeJob] = OrderedDict()
        self.job_runs: Dict[str, FakeJobRun] = OrderedDict()
        self.sessions: Dict[str, FakeSession] = OrderedDict()
        self.tagger = TaggingService()
        self.triggers: Dict[str, FakeTrigger] = OrderedDict()
        self.registries: Dict[str, FakeRegistry] = OrderedDict()
        self.num_schemas = 0
        self.num_schema_versions = 0

    def create_database(
        self,
        database_name: str,
        database_input: Dict[str, Any],
        tags: Optional[Dict[str, str]] = None,
    ) -> "FakeDatabase":
        if database_name in self.databases:
            raise DatabaseAlreadyExistsException()

        database = FakeDatabase(
            database_name, database_input, catalog_id=self.account_id
        )
        self.databases[database_name] = database
        resource_arn = f"arn:{get_partition(self.region_name)}:glue:{self.region_name}:{self.account_id}:database/{database_name}"
        self.tag_resource(resource_arn, tags)
        return database

    def get_database(self, database_name: str) -> "FakeDatabase":
        try:
            return self.databases[database_name]
        except KeyError:
            raise DatabaseNotFoundException(database_name)

    def update_database(
        self, database_name: str, database_input: Dict[str, Any]
    ) -> None:
        if database_name not in self.databases:
            raise DatabaseNotFoundException(database_name)

        self.databases[database_name].input = database_input

    def get_databases(self) -> List["FakeDatabase"]:
        return [self.databases[key] for key in self.databases] if self.databases else []

    def delete_database(self, database_name: str) -> None:
        if database_name not in self.databases:
            raise DatabaseNotFoundException(database_name)
        del self.databases[database_name]

    def create_table(
        self, database_name: str, table_name: str, table_input: Dict[str, Any]
    ) -> "FakeTable":
        database = self.get_database(database_name)

        if table_name in database.tables:
            raise TableAlreadyExistsException()

        table = FakeTable(
            database_name, table_name, table_input, catalog_id=self.account_id
        )
        database.tables[table_name] = table
        return table

    def get_table(self, database_name: str, table_name: str) -> "FakeTable":
        database = self.get_database(database_name)
        try:
            return database.tables[table_name]
        except KeyError:
            raise TableNotFoundException(table_name)

    def get_tables(
        self, database_name: str, expression: Optional[str]
    ) -> List["FakeTable"]:
        database = self.get_database(database_name)
        if expression:
            # sanitise expression, * is treated as a glob-like wildcard
            # so we make it a valid regex
            if "*" in expression:
                if expression.endswith(".*"):
                    expression = (
                        f"{expression[:-2].replace('*', '.*')}{expression[-2:]}"
                    )
                else:
                    expression = expression.replace("*", ".*")
            return [
                table
                for table_name, table in database.tables.items()
                if re.match(expression, table_name)
            ]
        else:
            return [table for table_name, table in database.tables.items()]

    def delete_table(self, database_name: str, table_name: str) -> None:
        database = self.get_database(database_name)
        try:
            del database.tables[table_name]
        except KeyError:
            raise TableNotFoundException(table_name)

    def update_table(
        self, database_name: str, table_name: str, table_input: Dict[str, Any]
    ) -> None:
        table = self.get_table(database_name, table_name)
        table.update(table_input)

    def get_table_version(
        self, database_name: str, table_name: str, ver_id: str
    ) -> str:
        table = self.get_table(database_name, table_name)

        return json.dumps(
            {
                "TableVersion": {
                    "Table": table.as_dict(version=ver_id),
                    "VersionId": ver_id,
                }
            }
        )

    def get_table_versions(
        self, database_name: str, table_name: str
    ) -> Dict[str, Dict[str, Any]]:
        table = self.get_table(database_name, table_name)
        return {version: table.as_dict(version) for version in table.versions.keys()}

    def delete_table_version(
        self, database_name: str, table_name: str, version_id: str
    ) -> None:
        table = self.get_table(database_name, table_name)
        table.delete_version(version_id)

    def create_partition(
        self, database_name: str, table_name: str, part_input: Dict[str, Any]
    ) -> None:
        table = self.get_table(database_name, table_name)
        table.create_partition(part_input)

    def get_partition(
        self, database_name: str, table_name: str, values: str
    ) -> "FakePartition":
        table = self.get_table(database_name, table_name)
        return table.get_partition(values)

    def get_partitions(
        self, database_name: str, table_name: str, expression: str
    ) -> List["FakePartition"]:
        """
        See https://docs.aws.amazon.com/glue/latest/webapi/API_GetPartitions.html
        for supported expressions.

        Expression caveats:

        - Column names must consist of UPPERCASE, lowercase, dots and underscores only.
        - Literal dates and timestamps must be valid, i.e. no support for February 31st.
        - LIKE expressions are converted to Python regexes, escaping special characters.
          Only % and _ wildcards are supported, and SQL escaping using [] does not work.
        """
        table = self.get_table(database_name, table_name)
        return table.get_partitions(expression)

    def update_partition(
        self,
        database_name: str,
        table_name: str,
        part_input: Dict[str, Any],
        part_to_update: str,
    ) -> None:
        table = self.get_table(database_name, table_name)
        table.update_partition(part_to_update, part_input)

    def delete_partition(
        self, database_name: str, table_name: str, part_to_delete: int
    ) -> None:
        table = self.get_table(database_name, table_name)
        table.delete_partition(part_to_delete)

    def create_crawler(
        self,
        name: str,
        role: str,
        database_name: str,
        description: str,
        targets: Dict[str, Any],
        schedule: str,
        classifiers: List[str],
        table_prefix: str,
        schema_change_policy: Dict[str, str],
        recrawl_policy: Dict[str, str],
        lineage_configuration: Dict[str, str],
        configuration: str,
        crawler_security_configuration: str,
        tags: Dict[str, str],
    ) -> None:
        if name in self.crawlers:
            raise CrawlerAlreadyExistsException()

        crawler = FakeCrawler(
            name=name,
            role=role,
            database_name=database_name,
            description=description,
            targets=targets,
            schedule=schedule,
            classifiers=classifiers,
            table_prefix=table_prefix,
            schema_change_policy=schema_change_policy,
            recrawl_policy=recrawl_policy,
            lineage_configuration=lineage_configuration,
            configuration=configuration,
            crawler_security_configuration=crawler_security_configuration,
            tags=tags,
            backend=self,
        )
        self.crawlers[name] = crawler

    def get_crawler(self, name: str) -> "FakeCrawler":
        try:
            return self.crawlers[name]
        except KeyError:
            raise CrawlerNotFoundException(name)

    def get_crawlers(self) -> List["FakeCrawler"]:
        return [self.crawlers[key] for key in self.crawlers] if self.crawlers else []

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_crawlers(self) -> List["FakeCrawler"]:
        return [crawler for _, crawler in self.crawlers.items()]

    def start_crawler(self, name: str) -> None:
        crawler = self.get_crawler(name)
        crawler.start_crawler()

    def stop_crawler(self, name: str) -> None:
        crawler = self.get_crawler(name)
        crawler.stop_crawler()

    def delete_crawler(self, name: str) -> None:
        try:
            del self.crawlers[name]
        except KeyError:
            raise CrawlerNotFoundException(name)

    def create_job(
        self,
        name: str,
        role: str,
        command: str,
        description: str,
        log_uri: str,
        execution_property: Dict[str, int],
        default_arguments: Dict[str, str],
        non_overridable_arguments: Dict[str, str],
        connections: Dict[str, List[str]],
        max_retries: int,
        allocated_capacity: int,
        timeout: int,
        max_capacity: float,
        security_configuration: str,
        tags: Dict[str, str],
        notification_property: Dict[str, int],
        glue_version: str,
        number_of_workers: int,
        worker_type: str,
        code_gen_configuration_nodes: Dict[str, Any],
        execution_class: str,
        source_control_details: Dict[str, str],
    ) -> None:
        self.jobs[name] = FakeJob(
            name,
            role,
            command,
            description,
            log_uri,
            execution_property,
            default_arguments,
            non_overridable_arguments,
            connections,
            max_retries,
            allocated_capacity,
            timeout,
            max_capacity,
            security_configuration,
            tags,
            notification_property,
            glue_version,
            number_of_workers,
            worker_type,
            code_gen_configuration_nodes,
            execution_class,
            source_control_details,
            backend=self,
        )

    def get_job(self, name: str) -> "FakeJob":
        try:
            return self.jobs[name]
        except KeyError:
            raise JobNotFoundException(name)

    @paginate(pagination_model=PAGINATION_MODEL)
    def get_jobs(self) -> List["FakeJob"]:
        return [job for _, job in self.jobs.items()]

    def start_job_run(self, name: str) -> str:
        job = self.get_job(name)
        return job.start_job_run()

    def get_job_run(self, name: str, run_id: str) -> "FakeJobRun":
        job = self.get_job(name)
        return job.get_job_run(run_id)

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_jobs(self) -> List["FakeJob"]:
        return [job for _, job in self.jobs.items()]

    def delete_job(self, name: str) -> None:
        if name in self.jobs:
            del self.jobs[name]

    def get_tags(self, resource_id: str) -> Dict[str, str]:
        return self.tagger.get_tag_dict_for_resource(resource_id)

    def tag_resource(self, resource_arn: str, tags: Optional[Dict[str, str]]) -> None:
        tag_list = TaggingService.convert_dict_to_tags_input(tags or {})
        self.tagger.tag_resource(resource_arn, tag_list)

    def untag_resource(self, resource_arn: str, tag_keys: List[str]) -> None:
        self.tagger.untag_resource_using_names(resource_arn, tag_keys)

    def create_registry(
        self,
        registry_name: str,
        description: Optional[str] = None,
        tags: Optional[Dict[str, str]] = None,
    ) -> Dict[str, Any]:
        # If registry name id default-registry, create default-registry
        if registry_name == DEFAULT_REGISTRY_NAME:
            registry = FakeRegistry(self, registry_name, description, tags)
            self.registries[registry_name] = registry
            return registry  # type: ignore

        # Validate Registry Parameters
        validate_registry_params(self.registries, registry_name, description, tags)

        registry = FakeRegistry(self, registry_name, description, tags)
        self.registries[registry_name] = registry
        return registry.as_dict()

    def delete_registry(self, registry_id: Dict[str, Any]) -> Dict[str, Any]:
        registry_name = validate_registry_id(registry_id, self.registries)
        return self.registries.pop(registry_name).as_dict()

    def get_registry(self, registry_id: Dict[str, Any]) -> Dict[str, Any]:
        registry_name = validate_registry_id(registry_id, self.registries)
        return self.registries[registry_name].as_dict()

    def list_registries(self) -> List[Dict[str, Any]]:
        return [reg.as_dict() for reg in self.registries.values()]

    def create_schema(
        self,
        registry_id: Dict[str, Any],
        schema_name: str,
        data_format: str,
        compatibility: str,
        schema_definition: str,
        description: Optional[str] = None,
        tags: Optional[Dict[str, str]] = None,
    ) -> Dict[str, Any]:
        """
        The following parameters/features are not yet implemented: Glue Schema Registry: compatibility checks NONE | BACKWARD | BACKWARD_ALL | FORWARD | FORWARD_ALL | FULL | FULL_ALL and  Data format parsing and syntax validation.
        """

        # Validate Registry Id
        registry_name = validate_registry_id(registry_id, self.registries)
        if (
            registry_name == DEFAULT_REGISTRY_NAME
            and DEFAULT_REGISTRY_NAME not in self.registries
        ):
            self.create_registry(registry_name)
        registry = self.registries[registry_name]

        # Validate Schema Parameters
        validate_schema_params(
            registry,
            schema_name,
            data_format,
            compatibility,
            schema_definition,
            self.num_schemas,
            description,
            tags,
        )

        # Create Schema
        schema_version = FakeSchemaVersion(
            self,
            registry_name,
            schema_name,
            schema_definition,
            version_number=1,
        )
        schema_version_id = schema_version.get_schema_version_id()
        schema = FakeSchema(
            self,
            registry_name,
            schema_name,
            data_format,
            compatibility,
            schema_version_id,
            description,
        )
        self.tagger.tag_resource(
            schema.schema_arn, self.tagger.convert_dict_to_tags_input(tags)
        )
        registry.schemas[schema_name] = schema
        self.num_schemas += 1

        schema.schema_versions[schema.schema_version_id] = schema_version
        self.num_schema_versions += 1

        resp = schema.as_dict()
        if tags:
            resp.update({"Tags": tags})
        return resp

    def register_schema_version(
        self, schema_id: Dict[str, Any], schema_definition: str
    ) -> Dict[str, Any]:
        # Validate Schema Id
        registry_name, schema_name, schema_arn = validate_schema_id(
            schema_id, self.registries
        )

        compatibility = (
            self.registries[registry_name].schemas[schema_name].compatibility
        )
        data_format = self.registries[registry_name].schemas[schema_name].data_format

        validate_register_schema_version_params(
            registry_name,
            schema_name,
            schema_arn,
            self.num_schema_versions,
            schema_definition,
            compatibility,
            data_format,
        )

        # If the same schema definition is already stored in Schema Registry as a version,
        # the schema ID of the existing schema is returned to the caller.
        schema_versions = (
            self.registries[registry_name].schemas[schema_name].schema_versions.values()
        )
        existing_schema_version = get_schema_version_if_definition_exists(
            schema_versions, data_format, schema_definition
        )
        if existing_schema_version:
            return existing_schema_version

        # Register Schema Version
        version_number = (
            self.registries[registry_name]
            .schemas[schema_name]
            .get_next_schema_version()
        )
        self.registries[registry_name].schemas[schema_name].update_next_schema_version()

        self.registries[registry_name].schemas[
            schema_name
        ].update_latest_schema_version()
        self.num_schema_versions += 1

        schema_version = FakeSchemaVersion(
            self,
            registry_name,
            schema_name,
            schema_definition,
            version_number,
        )
        self.registries[registry_name].schemas[schema_name].schema_versions[
            schema_version.schema_version_id
        ] = schema_version

        return schema_version.as_dict()

    def get_schema_version(
        self,
        schema_id: Optional[Dict[str, str]] = None,
        schema_version_id: Optional[str] = None,
        schema_version_number: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        # Validate Schema Parameters
        (
            schema_version_id,
            registry_name,
            schema_name,
            schema_arn,
            version_number,
            latest_version,
        ) = validate_schema_version_params(
            self.registries, schema_id, schema_version_id, schema_version_number
        )

        # GetSchemaVersion using SchemaVersionId
        if schema_version_id:
            for registry in self.registries.values():
                for schema in registry.schemas.values():
                    if (
                        schema.schema_versions.get(schema_version_id, None)
                        and schema.schema_versions[
                            schema_version_id
                        ].schema_version_status
                        != DELETING_STATUS
                    ):
                        get_schema_version_dict = schema.schema_versions[
                            schema_version_id
                        ].get_schema_version_as_dict()
                        get_schema_version_dict["DataFormat"] = schema.data_format
                        return get_schema_version_dict
            raise SchemaVersionNotFoundFromSchemaVersionIdException(schema_version_id)

        # GetSchemaVersion using VersionNumber
        schema = self.registries[registry_name].schemas[schema_name]  # type: ignore
        for schema_version in schema.schema_versions.values():
            if (
                version_number == schema_version.version_number  # type: ignore
                and schema_version.schema_version_status != DELETING_STATUS
            ):
                get_schema_version_dict = schema_version.get_schema_version_as_dict()
                get_schema_version_dict["DataFormat"] = schema.data_format
                return get_schema_version_dict
        raise SchemaVersionNotFoundFromSchemaIdException(
            registry_name, schema_name, schema_arn, version_number, latest_version
        )

    def get_schema_by_definition(
        self, schema_id: Dict[str, str], schema_definition: str
    ) -> Dict[str, Any]:
        # Validate SchemaId
        validate_schema_definition_length(schema_definition)
        registry_name, schema_name, schema_arn = validate_schema_id(
            schema_id, self.registries
        )

        # Get Schema By Definition
        schema = self.registries[registry_name].schemas[schema_name]
        for schema_version in schema.schema_versions.values():
            if (
                schema_definition == schema_version.schema_definition
                and schema_version.schema_version_status != DELETING_STATUS
            ):
                get_schema_by_definition_dict = (
                    schema_version.get_schema_by_definition_as_dict()
                )
                get_schema_by_definition_dict["DataFormat"] = schema.data_format
                return get_schema_by_definition_dict
        raise SchemaNotFoundException(schema_name, registry_name, schema_arn)

    def put_schema_version_metadata(
        self,
        schema_id: Dict[str, Any],
        schema_version_number: Dict[str, str],
        schema_version_id: str,
        metadata_key_value: Dict[str, str],
    ) -> Dict[str, Any]:
        # Validate metadata_key_value and schema version params
        (
            metadata_key,
            metadata_value,
        ) = validate_schema_version_metadata_pattern_and_length(metadata_key_value)
        (
            schema_version_id,
            registry_name,
            schema_name,
            schema_arn,
            version_number,
            latest_version,
        ) = validate_schema_version_params(  # type: ignore
            self.registries, schema_id, schema_version_id, schema_version_number
        )

        # PutSchemaVersionMetadata using SchemaVersionId
        if schema_version_id:
            for registry in self.registries.values():
                for schema in registry.schemas.values():
                    if schema.schema_versions.get(schema_version_id, None):
                        metadata = schema.schema_versions[schema_version_id].metadata
                        validate_number_of_schema_version_metadata_allowed(metadata)

                        if metadata_key in metadata:
                            if metadata_value in metadata[metadata_key]:
                                raise SchemaVersionMetadataAlreadyExistsException(
                                    schema_version_id, metadata_key, metadata_value
                                )
                            metadata[metadata_key].append(metadata_value)
                        else:
                            metadata[metadata_key] = [metadata_value]
                        return get_put_schema_version_metadata_response(
                            schema_id,
                            schema_version_number,
                            schema_version_id,
                            metadata_key_value,
                        )

            raise SchemaVersionNotFoundFromSchemaVersionIdException(schema_version_id)

        # PutSchemaVersionMetadata using VersionNumber
        schema = self.registries[registry_name].schemas[schema_name]  # type: ignore
        for schema_version in schema.schema_versions.values():
            if version_number == schema_version.version_number:  # type: ignore
                validate_number_of_schema_version_metadata_allowed(
                    schema_version.metadata
                )
                if metadata_key in schema_version.metadata:
                    if metadata_value in schema_version.metadata[metadata_key]:
                        raise SchemaVersionMetadataAlreadyExistsException(
                            schema_version.schema_version_id,
                            metadata_key,
                            metadata_value,
                        )
                    schema_version.metadata[metadata_key].append(metadata_value)
                else:
                    schema_version.metadata[metadata_key] = [metadata_value]
                return get_put_schema_version_metadata_response(
                    schema_id,
                    schema_version_number,
                    schema_version_id,
                    metadata_key_value,
                )

        raise SchemaVersionNotFoundFromSchemaIdException(
            registry_name, schema_name, schema_arn, version_number, latest_version
        )

    def get_schema(self, schema_id: Dict[str, str]) -> Dict[str, Any]:
        registry_name, schema_name, _ = validate_schema_id(schema_id, self.registries)
        schema = self.registries[registry_name].schemas[schema_name]
        return schema.as_dict()

    def delete_schema(self, schema_id: Dict[str, str]) -> Dict[str, Any]:
        # Validate schema_id
        registry_name, schema_name, _ = validate_schema_id(schema_id, self.registries)

        # delete schema pre-processing
        schema = self.registries[registry_name].schemas[schema_name]
        num_schema_version_in_schema = len(schema.schema_versions)
        schema.schema_status = DELETING_STATUS
        response = delete_schema_response(
            schema.schema_name, schema.schema_arn, schema.schema_status
        )

        # delete schema
        del self.registries[registry_name].schemas[schema_name]
        self.num_schemas -= 1
        self.num_schema_versions -= num_schema_version_in_schema

        return response

    def update_schema(
        self, schema_id: Dict[str, str], compatibility: str, description: str
    ) -> Dict[str, Any]:
        """
        The SchemaVersionNumber-argument is not yet implemented
        """
        registry_name, schema_name, _ = validate_schema_id(schema_id, self.registries)
        schema = self.registries[registry_name].schemas[schema_name]

        if compatibility is not None:
            schema.compatibility = compatibility
        if description is not None:
            schema.description = description

        return schema.as_dict()

    def create_session(
        self,
        session_id: str,
        description: str,
        role: str,
        command: Dict[str, str],
        timeout: int,
        idle_timeout: int,
        default_arguments: Dict[str, str],
        connections: Dict[str, List[str]],
        max_capacity: float,
        number_of_workers: int,
        worker_type: str,
        security_configuration: str,
        glue_version: str,
        tags: Dict[str, str],
        request_origin: str,
    ) -> None:
        if session_id in self.sessions:
            raise SessionAlreadyExistsException()

        session = FakeSession(
            session_id,
            description,
            role,
            command,
            timeout,
            idle_timeout,
            default_arguments,
            connections,
            max_capacity,
            number_of_workers,
            worker_type,
            security_configuration,
            glue_version,
            tags,
            request_origin,
            backend=self,
        )
        self.sessions[session_id] = session

    def get_session(self, session_id: str) -> "FakeSession":
        try:
            return self.sessions[session_id]
        except KeyError:
            raise SessionNotFoundException(session_id)

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_sessions(self) -> List["FakeSession"]:
        return [session for _, session in self.sessions.items()]

    def stop_session(self, session_id: str) -> None:
        session = self.get_session(session_id)
        session.stop_session()

    def delete_session(self, session_id: str) -> None:
        try:
            del self.sessions[session_id]
        except KeyError:
            raise SessionNotFoundException(session_id)

    def create_trigger(
        self,
        name: str,
        workflow_name: str,
        trigger_type: str,
        schedule: str,
        predicate: Dict[str, Any],
        actions: List[Dict[str, Any]],
        description: str,
        start_on_creation: bool,
        tags: Dict[str, str],
        event_batching_condition: Dict[str, Any],
    ) -> None:
        self.triggers[name] = FakeTrigger(
            name=name,
            workflow_name=workflow_name,
            trigger_type=trigger_type,
            schedule=schedule,
            predicate=predicate,
            actions=actions,
            description=description,
            start_on_creation=start_on_creation,
            tags=tags,
            event_batching_condition=event_batching_condition,
            backend=self,
        )

    def get_trigger(self, name: str) -> "FakeTrigger":
        try:
            return self.triggers[name]
        except KeyError:
            raise TriggerNotFoundException(name)

    def start_trigger(self, name: str) -> None:
        trigger = self.get_trigger(name)
        trigger.start_trigger()

    def stop_trigger(self, name: str) -> None:
        trigger = self.get_trigger(name)
        trigger.stop_trigger()

    @paginate(pagination_model=PAGINATION_MODEL)
    def get_triggers(self, dependent_job_name: str) -> List["FakeTrigger"]:
        if dependent_job_name:
            triggers = []
            for trigger in self.triggers.values():
                for action in trigger.actions:
                    if ("JobName" in action) and (
                        action["JobName"] == dependent_job_name
                    ):
                        triggers.append(trigger)
            return triggers

        return list(self.triggers.values())

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_triggers(self, dependent_job_name: str) -> List["FakeTrigger"]:
        if dependent_job_name:
            triggers = []
            for trigger in self.triggers.values():
                for action in trigger.actions:
                    if ("JobName" in action) and (
                        action["JobName"] == dependent_job_name
                    ):
                        triggers.append(trigger)
            return triggers

        return list(self.triggers.values())

    def delete_trigger(self, name: str) -> None:
        if name in self.triggers:
            del self.triggers[name]

    def batch_delete_table(
        self, database_name: str, tables: List[str]
    ) -> List[Dict[str, Any]]:
        errors = []
        for table_name in tables:
            try:
                self.delete_table(database_name, table_name)
            except TableNotFoundException:
                errors.append(
                    {
                        "TableName": table_name,
                        "ErrorDetail": {
                            "ErrorCode": "EntityNotFoundException",
                            "ErrorMessage": "Table not found",
                        },
                    }
                )
        return errors

    def batch_get_partition(
        self,
        database_name: str,
        table_name: str,
        partitions_to_get: List[Dict[str, str]],
    ) -> List[Dict[str, Any]]:
        table = self.get_table(database_name, table_name)

        partitions = []
        for values in partitions_to_get:
            try:
                p = table.get_partition(values=values["Values"])
                partitions.append(p.as_dict())
            except PartitionNotFoundException:
                continue
        return partitions

    def batch_create_partition(
        self, database_name: str, table_name: str, partition_input: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        table = self.get_table(database_name, table_name)

        errors_output = []
        for part_input in partition_input:
            try:
                table.create_partition(part_input)
            except PartitionAlreadyExistsException:
                errors_output.append(
                    {
                        "PartitionValues": part_input["Values"],
                        "ErrorDetail": {
                            "ErrorCode": "AlreadyExistsException",
                            "ErrorMessage": "Partition already exists.",
                        },
                    }
                )
        return errors_output

    def batch_update_partition(
        self, database_name: str, table_name: str, entries: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        table = self.get_table(database_name, table_name)

        errors_output = []
        for entry in entries:
            part_to_update = entry["PartitionValueList"]
            part_input = entry["PartitionInput"]

            try:
                table.update_partition(part_to_update, part_input)
            except PartitionNotFoundException:
                errors_output.append(
                    {
                        "PartitionValueList": part_to_update,
                        "ErrorDetail": {
                            "ErrorCode": "EntityNotFoundException",
                            "ErrorMessage": "Partition not found.",
                        },
                    }
                )
        return errors_output

    def batch_delete_partition(
        self, database_name: str, table_name: str, parts: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        table = self.get_table(database_name, table_name)

        errors_output = []
        for part_input in parts:
            values = part_input.get("Values")
            try:
                table.delete_partition(values)  # type: ignore
            except PartitionNotFoundException:
                errors_output.append(
                    {
                        "PartitionValues": values,
                        "ErrorDetail": {
                            "ErrorCode": "EntityNotFoundException",
                            "ErrorMessage": "Partition not found",
                        },
                    }
                )
        return errors_output

    def batch_get_crawlers(self, crawler_names: List[str]) -> List[Dict[str, Any]]:
        crawlers = []
        for crawler in self.get_crawlers():
            if crawler.as_dict()["Name"] in crawler_names:
                crawlers.append(crawler.as_dict())
        return crawlers

    def batch_get_jobs(self, job_names: List[str]) -> List[Dict[str, Any]]:
        jobs = []
        for job_name in job_names:
            if job_name in self.jobs:
                jobs.append(self.jobs[job_name].as_dict())
        return jobs

    def batch_get_triggers(self, trigger_names: List[str]) -> List[Dict[str, Any]]:
        triggers = []
        for trigger_name in trigger_names:
            if trigger_name in self.triggers:
                triggers.append(self.triggers[trigger_name].as_dict())
        return triggers


class FakeDatabase(BaseModel):
    def __init__(
        self, database_name: str, database_input: Dict[str, Any], catalog_id: str
    ):
        self.name = database_name
        self.input = database_input
        self.catalog_id = catalog_id
        self.created_time = utcnow()
        self.tables: Dict[str, FakeTable] = OrderedDict()

    def as_dict(self) -> Dict[str, Any]:
        return {
            "Name": self.name,
            "Description": self.input.get("Description"),
            "LocationUri": self.input.get("LocationUri"),
            "Parameters": self.input.get("Parameters"),
            "CreateTime": unix_time(self.created_time),
            "CreateTableDefaultPermissions": self.input.get(
                "CreateTableDefaultPermissions"
            ),
            "TargetDatabase": self.input.get("TargetDatabase"),
            "CatalogId": self.input.get("CatalogId") or self.catalog_id,
        }


class FakeTable(BaseModel):
    def __init__(
        self,
        database_name: str,
        table_name: str,
        table_input: Dict[str, Any],
        catalog_id: str,
    ):
        self.database_name = database_name
        self.name = table_name
        self.catalog_id = catalog_id
        self.partitions: Dict[str, FakePartition] = OrderedDict()
        self.created_time = utcnow()
        self.updated_time: Optional[datetime] = None
        self._current_version = 1
        self.versions: Dict[str, Dict[str, Any]] = {
            str(self._current_version): table_input
        }

    def update(self, table_input: Dict[str, Any]) -> None:
        self.versions[str(self._current_version + 1)] = table_input
        self._current_version += 1
        self.updated_time = utcnow()

    def get_version(self, ver: str) -> Dict[str, Any]:
        try:
            int(ver)
        except ValueError as e:
            raise JsonRESTError("InvalidInputException", str(e))

        try:
            return self.versions[ver]
        except KeyError:
            raise VersionNotFoundException()

    def delete_version(self, version_id: str) -> None:
        self.versions.pop(version_id)

    def as_dict(self, version: Optional[str] = None) -> Dict[str, Any]:
        version = version or self._current_version  # type: ignore
        obj = {
            "DatabaseName": self.database_name,
            "Name": self.name,
            "CreateTime": unix_time(self.created_time),
            **self.get_version(str(version)),
            # Add VersionId after we get the version-details, just to make sure that it's a valid version (int)
            "VersionId": str(version),
            "CatalogId": self.catalog_id,
        }
        if self.updated_time is not None:
            obj["UpdateTime"] = unix_time(self.updated_time)
        return obj

    def create_partition(self, partiton_input: Dict[str, Any]) -> None:
        partition = FakePartition(self.database_name, self.name, partiton_input)
        key = str(partition.values)
        if key in self.partitions:
            raise PartitionAlreadyExistsException()
        self.partitions[str(partition.values)] = partition

    def get_partitions(self, expression: str) -> List["FakePartition"]:
        return list(filter(PartitionFilter(expression, self), self.partitions.values()))

    def get_partition(self, values: str) -> "FakePartition":
        try:
            return self.partitions[str(values)]
        except KeyError:
            raise PartitionNotFoundException()

    def update_partition(self, old_values: str, partiton_input: Dict[str, Any]) -> None:
        partition = FakePartition(self.database_name, self.name, partiton_input)
        key = str(partition.values)
        if old_values == partiton_input["Values"]:
            # Altering a partition in place. Don't remove it so the order of
            # returned partitions doesn't change
            if key not in self.partitions:
                raise PartitionNotFoundException()
        else:
            removed = self.partitions.pop(str(old_values), None)
            if removed is None:
                raise PartitionNotFoundException()
            if key in self.partitions:
                # Trying to update to overwrite a partition that exists
                raise PartitionAlreadyExistsException()
        self.partitions[key] = partition

    def delete_partition(self, values: int) -> None:
        try:
            del self.partitions[str(values)]
        except KeyError:
            raise PartitionNotFoundException()


class FakePartition(BaseModel):
    def __init__(
        self, database_name: str, table_name: str, partiton_input: Dict[str, Any]
    ):
        self.creation_time = time.time()
        self.database_name = database_name
        self.table_name = table_name
        self.partition_input = partiton_input
        self.values = self.partition_input.get("Values", [])

    def as_dict(self) -> Dict[str, Any]:
        obj = {
            "DatabaseName": self.database_name,
            "TableName": self.table_name,
            "CreationTime": self.creation_time,
        }
        obj.update(self.partition_input)
        return obj


class FakeCrawler(BaseModel):
    def __init__(
        self,
        name: str,
        role: str,
        database_name: str,
        description: str,
        targets: Dict[str, Any],
        schedule: str,
        classifiers: List[str],
        table_prefix: str,
        schema_change_policy: Dict[str, str],
        recrawl_policy: Dict[str, str],
        lineage_configuration: Dict[str, str],
        configuration: str,
        crawler_security_configuration: str,
        tags: Dict[str, str],
        backend: GlueBackend,
    ):
        self.name = name
        self.role = role
        self.database_name = database_name
        self.description = description
        self.targets = targets
        self.schedule = schedule
        self.classifiers = classifiers
        self.table_prefix = table_prefix
        self.schema_change_policy = schema_change_policy
        self.recrawl_policy = recrawl_policy
        self.lineage_configuration = lineage_configuration
        self.configuration = configuration
        self.crawler_security_configuration = crawler_security_configuration
        self.state = "READY"
        self.creation_time = utcnow()
        self.last_updated = self.creation_time
        self.version = 1
        self.crawl_elapsed_time = 0
        self.last_crawl_info = None
        self.arn = f"arn:{get_partition(backend.region_name)}:glue:{backend.region_name}:{backend.account_id}:crawler/{self.name}"
        self.backend = backend
        self.backend.tag_resource(self.arn, tags)

    def get_name(self) -> str:
        return self.name

    def as_dict(self) -> Dict[str, Any]:
        last_crawl = self.last_crawl_info.as_dict() if self.last_crawl_info else None  # type: ignore
        data = {
            "Name": self.name,
            "Role": self.role,
            "Targets": self.targets,
            "DatabaseName": self.database_name,
            "Description": self.description,
            "Classifiers": self.classifiers,
            "RecrawlPolicy": self.recrawl_policy,
            "SchemaChangePolicy": self.schema_change_policy,
            "LineageConfiguration": self.lineage_configuration,
            "State": self.state,
            "TablePrefix": self.table_prefix,
            "CrawlElapsedTime": self.crawl_elapsed_time,
            "CreationTime": self.creation_time.isoformat(),
            "LastUpdated": self.last_updated.isoformat(),
            "LastCrawl": last_crawl,
            "Version": self.version,
            "Configuration": self.configuration,
            "CrawlerSecurityConfiguration": self.crawler_security_configuration,
        }

        if self.schedule:
            data["Schedule"] = {
                "ScheduleExpression": self.schedule,
                "State": "SCHEDULED",
            }

        if self.last_crawl_info:
            data["LastCrawl"] = self.last_crawl_info.as_dict()

        return data

    def start_crawler(self) -> None:
        if self.state == "RUNNING":
            raise CrawlerRunningException(
                f"Crawler with name {self.name} has already started"
            )
        self.state = "RUNNING"

    def stop_crawler(self) -> None:
        if self.state != "RUNNING":
            raise CrawlerNotRunningException(
                f"Crawler with name {self.name} isn't running"
            )
        self.state = "STOPPING"


class LastCrawlInfo(BaseModel):
    def __init__(
        self,
        error_message: str,
        log_group: str,
        log_stream: str,
        message_prefix: str,
        start_time: str,
        status: str,
    ):
        self.error_message = error_message
        self.log_group = log_group
        self.log_stream = log_stream
        self.message_prefix = message_prefix
        self.start_time = start_time
        self.status = status

    def as_dict(self) -> Dict[str, Any]:
        return {
            "ErrorMessage": self.error_message,
            "LogGroup": self.log_group,
            "LogStream": self.log_stream,
            "MessagePrefix": self.message_prefix,
            "StartTime": self.start_time,
            "Status": self.status,
        }


class FakeJob:
    def __init__(
        self,
        name: str,
        role: str,
        command: str,
        description: str,
        log_uri: str,
        execution_property: Dict[str, int],
        default_arguments: Dict[str, str],
        non_overridable_arguments: Dict[str, str],
        connections: Dict[str, List[str]],
        max_retries: int,
        allocated_capacity: int,
        timeout: int,
        max_capacity: float,
        security_configuration: str,
        tags: Dict[str, str],
        notification_property: Dict[str, int],
        glue_version: str,
        number_of_workers: int,
        worker_type: str,
        code_gen_configuration_nodes: Dict[str, Any],
        execution_class: str,
        source_control_details: Dict[str, str],
        backend: GlueBackend,
    ):
        self.name = name
        self.description = description
        self.log_uri = log_uri
        self.role = role
        self.execution_property = execution_property or {}
        self.command = command
        self.default_arguments = default_arguments
        self.non_overridable_arguments = non_overridable_arguments
        self.connections = connections
        self.max_retries = max_retries
        self.allocated_capacity = allocated_capacity
        self.timeout = timeout
        self.max_capacity = max_capacity
        self.security_configuration = security_configuration
        self.notification_property = notification_property
        self.glue_version = glue_version
        self.number_of_workers = number_of_workers
        self.worker_type = worker_type
        self.code_gen_configuration_nodes = code_gen_configuration_nodes
        self.execution_class = execution_class or "STANDARD"
        self.source_control_details = source_control_details
        self.created_on = utcnow()
        self.last_modified_on = utcnow()
        self.arn = f"arn:{get_partition(backend.region_name)}:glue:{backend.region_name}:{backend.account_id}:job/{self.name}"
        self.backend = backend
        self.backend.tag_resource(self.arn, tags)

        self.job_runs: List[FakeJobRun] = []

    def get_name(self) -> str:
        return self.name

    def as_dict(self) -> Dict[str, Any]:
        return {
            "Name": self.name,
            "Description": self.description,
            "LogUri": self.log_uri,
            "Role": self.role,
            "CreatedOn": self.created_on.isoformat(),
            "LastModifiedOn": self.last_modified_on.isoformat(),
            "ExecutionProperty": self.execution_property,
            "Command": self.command,
            "DefaultArguments": self.default_arguments,
            "NonOverridableArguments": self.non_overridable_arguments,
            "Connections": self.connections,
            "MaxRetries": self.max_retries,
            "AllocatedCapacity": self.allocated_capacity,
            "Timeout": self.timeout,
            "MaxCapacity": self.max_capacity,
            "WorkerType": self.worker_type,
            "NumberOfWorkers": self.number_of_workers,
            "SecurityConfiguration": self.security_configuration,
            "NotificationProperty": self.notification_property,
            "GlueVersion": self.glue_version,
            "CodeGenConfigurationNodes": self.code_gen_configuration_nodes,
            "ExecutionClass": self.execution_class,
            "SourceControlDetails": self.source_control_details,
        }

    def start_job_run(self) -> str:
        running_jobs = len(
            [jr for jr in self.job_runs if jr.status in ["STARTING", "RUNNING"]]
        )
        if running_jobs >= self.execution_property.get("MaxConcurrentRuns", 1):
            raise ConcurrentRunsExceededException(
                f"Job with name {self.name} already running"
            )
        fake_job_run = FakeJobRun(job_name=self.name)
        self.job_runs.append(fake_job_run)
        return fake_job_run.job_run_id

    def get_job_run(self, run_id: str) -> "FakeJobRun":
        for job_run in self.job_runs:
            if job_run.job_run_id == run_id:
                job_run.advance()
                return job_run
        raise JobRunNotFoundException(run_id)


class FakeJobRun(ManagedState):
    def __init__(
        self,
        job_name: str,
        job_run_id: str = "01",
        arguments: Optional[Dict[str, Any]] = None,
        allocated_capacity: Optional[int] = None,
        timeout: Optional[int] = None,
        worker_type: str = "Standard",
    ):
        ManagedState.__init__(
            self,
            model_name="glue::job_run",
            transitions=[("STARTING", "RUNNING"), ("RUNNING", "SUCCEEDED")],
        )
        self.job_name = job_name
        self.job_run_id = job_run_id
        self.arguments = arguments
        self.allocated_capacity = allocated_capacity
        self.timeout = timeout
        self.worker_type = worker_type
        self.started_on = utcnow()
        self.modified_on = utcnow()
        self.completed_on = utcnow()

    def get_name(self) -> str:
        return self.job_name

    def as_dict(self) -> Dict[str, Any]:
        return {
            "Id": self.job_run_id,
            "Attempt": 1,
            "PreviousRunId": "01",
            "TriggerName": "test_trigger",
            "JobName": self.job_name,
            "StartedOn": self.started_on.isoformat(),
            "LastModifiedOn": self.modified_on.isoformat(),
            "CompletedOn": self.completed_on.isoformat(),
            "JobRunState": self.status,
            "Arguments": self.arguments or {"runSpark": "spark -f test_file.py"},
            "ErrorMessage": "",
            "PredecessorRuns": [
                {"JobName": "string", "RunId": "string"},
            ],
            "AllocatedCapacity": self.allocated_capacity or 123,
            "ExecutionTime": 123,
            "Timeout": self.timeout or 123,
            "MaxCapacity": 123.0,
            "WorkerType": self.worker_type,
            "NumberOfWorkers": 123,
            "SecurityConfiguration": "string",
            "LogGroupName": "test/log",
            "NotificationProperty": {"NotifyDelayAfter": 123},
            "GlueVersion": "0.9",
        }


class FakeRegistry(BaseModel):
    def __init__(
        self,
        backend: GlueBackend,
        registry_name: str,
        description: Optional[str] = None,
        tags: Optional[Dict[str, str]] = None,
    ):
        self.name = registry_name
        self.description = description
        self.tags = tags
        self.created_time = utcnow()
        self.updated_time = utcnow()
        self.status = "AVAILABLE"
        self.registry_arn = f"arn:{get_partition(backend.region_name)}:glue:{backend.region_name}:{backend.account_id}:registry/{self.name}"
        self.schemas: Dict[str, FakeSchema] = OrderedDict()

    def as_dict(self) -> Dict[str, Any]:
        return {
            "RegistryArn": self.registry_arn,
            "RegistryName": self.name,
            "Description": self.description,
            "Tags": self.tags,
        }


class FakeSchema(BaseModel):
    def __init__(
        self,
        backend: GlueBackend,
        registry_name: str,
        schema_name: str,
        data_format: str,
        compatibility: str,
        schema_version_id: str,
        description: Optional[str] = None,
    ):
        self.registry_name = registry_name
        self.registry_arn = f"arn:{get_partition(backend.region_name)}:glue:{backend.region_name}:{backend.account_id}:registry/{self.registry_name}"
        self.schema_name = schema_name
        self.schema_arn = f"arn:{get_partition(backend.region_name)}:glue:{backend.region_name}:{backend.account_id}:schema/{self.registry_name}/{self.schema_name}"
        self.description = description
        self.data_format = data_format
        self.compatibility = compatibility
        self.schema_checkpoint = 1
        self.latest_schema_version = 1
        self.next_schema_version = 2
        self.schema_status = AVAILABLE_STATUS
        self.schema_version_id = schema_version_id
        self.schema_version_status = AVAILABLE_STATUS
        self.created_time = utcnow()
        self.updated_time = utcnow()
        self.schema_versions: Dict[str, FakeSchemaVersion] = OrderedDict()

    def update_next_schema_version(self) -> None:
        self.next_schema_version += 1

    def update_latest_schema_version(self) -> None:
        self.latest_schema_version += 1

    def get_next_schema_version(self) -> int:
        return self.next_schema_version

    def as_dict(self) -> Dict[str, Any]:
        return {
            "RegistryArn": self.registry_arn,
            "RegistryName": self.registry_name,
            "SchemaName": self.schema_name,
            "SchemaArn": self.schema_arn,
            "DataFormat": self.data_format,
            "Compatibility": self.compatibility,
            "SchemaCheckpoint": self.schema_checkpoint,
            "LatestSchemaVersion": self.latest_schema_version,
            "NextSchemaVersion": self.next_schema_version,
            "SchemaStatus": self.schema_status,
            "SchemaVersionId": self.schema_version_id,
            "SchemaVersionStatus": self.schema_version_status,
            "Description": self.description,
        }


class FakeSchemaVersion(BaseModel):
    def __init__(
        self,
        backend: GlueBackend,
        registry_name: str,
        schema_name: str,
        schema_definition: str,
        version_number: int,
    ):
        self.registry_name = registry_name
        self.schema_name = schema_name
        self.schema_arn = f"arn:{get_partition(backend.region_name)}:glue:{backend.region_name}:{backend.account_id}:schema/{self.registry_name}/{self.schema_name}"
        self.schema_definition = schema_definition
        self.schema_version_status = AVAILABLE_STATUS
        self.version_number = version_number
        self.schema_version_id = str(mock_random.uuid4())
        self.created_time = utcnow()
        self.updated_time = utcnow()
        self.metadata: Dict[str, Any] = {}

    def get_schema_version_id(self) -> str:
        return self.schema_version_id

    def as_dict(self) -> Dict[str, Any]:
        return {
            "SchemaVersionId": self.schema_version_id,
            "VersionNumber": self.version_number,
            "Status": self.schema_version_status,
        }

    def get_schema_version_as_dict(self) -> Dict[str, Any]:
        # add data_format for full return dictionary of get_schema_version
        return {
            "SchemaVersionId": self.schema_version_id,
            "SchemaDefinition": self.schema_definition,
            "SchemaArn": self.schema_arn,
            "VersionNumber": self.version_number,
            "Status": self.schema_version_status,
            "CreatedTime": str(self.created_time),
        }

    def get_schema_by_definition_as_dict(self) -> Dict[str, Any]:
        # add data_format for full return dictionary of get_schema_by_definition
        return {
            "SchemaVersionId": self.schema_version_id,
            "SchemaArn": self.schema_arn,
            "Status": self.schema_version_status,
            "CreatedTime": str(self.created_time),
        }


class FakeSession(BaseModel):
    def __init__(
        self,
        session_id: str,
        description: str,
        role: str,
        command: Dict[str, str],
        timeout: int,
        idle_timeout: int,
        default_arguments: Dict[str, str],
        connections: Dict[str, List[str]],
        max_capacity: float,
        number_of_workers: int,
        worker_type: str,
        security_configuration: str,
        glue_version: str,
        tags: Dict[str, str],
        request_origin: str,
        backend: GlueBackend,
    ):
        self.session_id = session_id
        self.description = description
        self.role = role
        self.command = command
        self.timeout = timeout
        self.idle_timeout = idle_timeout
        self.default_arguments = default_arguments
        self.connections = connections
        self.max_capacity = max_capacity
        self.number_of_workers = number_of_workers
        self.worker_type = worker_type
        self.security_configuration = security_configuration
        self.glue_version = glue_version
        self.tags = tags
        self.request_origin = request_origin
        self.creation_time = utcnow()
        self.last_updated = self.creation_time
        self.arn = f"arn:{get_partition(backend.region_name)}:glue:{backend.region_name}:{backend.account_id}:session/{self.session_id}"
        self.backend = backend
        self.backend.tag_resource(self.arn, tags)
        self.state = "READY"

    def get_id(self) -> str:
        return self.session_id

    def as_dict(self) -> Dict[str, Any]:
        return {
            "Id": self.session_id,
            "CreatedOn": self.creation_time.isoformat(),
            "Status": self.state,
            "ErrorMessage": "string",
            "Description": self.description,
            "Role": self.role,
            "Command": self.command,
            "DefaultArguments": self.default_arguments,
            "Connections": self.connections,
            "Progress": 12.3,
            "MaxCapacity": 123.0,
            "SecurityConfiguration": self.security_configuration,
            "GlueVersion": self.glue_version,
        }

    def stop_session(self) -> None:
        if self.state != "READY":
            raise IllegalSessionStateException(f"Session is in {self.state} status")
        self.state = "STOPPING"


class FakeTrigger(BaseModel):
    def __init__(
        self,
        backend: GlueBackend,
        name: str,
        workflow_name: str,
        trigger_type: str,  # to avoid any issues with built-in function type()
        schedule: str,
        predicate: Dict[str, Any],
        actions: List[Dict[str, Any]],
        description: str,
        start_on_creation: bool,
        tags: Dict[str, str],
        event_batching_condition: Dict[str, Any],
    ):
        self.name = name
        self.workflow_name = workflow_name
        self.trigger_type = trigger_type
        self.schedule = schedule
        self.predicate = predicate
        self.actions = actions
        self.description = description
        if start_on_creation:
            self.state = "ACTIVATED"
        else:
            self.state = "CREATED"
        self.event_batching_condition = event_batching_condition
        self.arn = f"arn:{get_partition(backend.region_name)}:glue:{backend.region_name}:{backend.account_id}:trigger/{self.name}"
        self.backend = backend
        self.backend.tag_resource(self.arn, tags)

    def get_name(self) -> str:
        return self.name

    def start_trigger(self) -> None:
        self.state = "ACTIVATED"

    def stop_trigger(self) -> None:
        self.state = "DEACTIVATED"

    def as_dict(self) -> Dict[str, Any]:
        data: Dict[str, Any] = {
            "Name": self.name,
            "Type": self.trigger_type,
            "Actions": self.actions,
            "State": self.state,
        }

        if self.workflow_name:
            data["WorkflowName"] = self.workflow_name

        if self.trigger_type == "SCHEDULED":
            data["Schedule"] = self.schedule

        if self.predicate:
            data["Predicate"] = self.predicate

        if self.description:
            data["Description"] = self.description

        if self.event_batching_condition:
            data["EventBatchingCondition"] = self.event_batching_condition

        return data


glue_backends = BackendDict(GlueBackend, "glue")
