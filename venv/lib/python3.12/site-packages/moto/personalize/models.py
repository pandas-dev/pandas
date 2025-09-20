from typing import Any, Dict, Iterable

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import unix_time
from moto.utilities.utils import get_partition

from .exceptions import ResourceNotFoundException


class Schema(BaseModel):
    def __init__(
        self,
        account_id: str,
        region: str,
        name: str,
        schema: Dict[str, Any],
        domain: str,
    ):
        self.name = name
        self.schema = schema
        self.domain = domain
        self.arn = f"arn:{get_partition(region)}:personalize:{region}:{account_id}:schema/{name}"
        self.created = unix_time()

    def to_dict(self, full: bool = True) -> Dict[str, Any]:
        d: Dict[str, Any] = {
            "name": self.name,
            "schemaArn": self.arn,
            "domain": self.domain,
            "creationDateTime": self.created,
            "lastUpdatedDateTime": self.created,
        }
        if full:
            d["schema"] = self.schema
        return d


class PersonalizeBackend(BaseBackend):
    """Implementation of Personalize APIs."""

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.schemas: Dict[str, Schema] = dict()

    def create_schema(self, name: str, schema_dict: Dict[str, Any], domain: str) -> str:
        schema = Schema(
            region=self.region_name,
            account_id=self.account_id,
            name=name,
            schema=schema_dict,
            domain=domain,
        )
        self.schemas[schema.arn] = schema
        return schema.arn

    def delete_schema(self, schema_arn: str) -> None:
        if schema_arn not in self.schemas:
            raise ResourceNotFoundException(schema_arn)
        self.schemas.pop(schema_arn, None)

    def describe_schema(self, schema_arn: str) -> Schema:
        if schema_arn not in self.schemas:
            raise ResourceNotFoundException(schema_arn)
        return self.schemas[schema_arn]

    def list_schemas(self) -> Iterable[Schema]:
        """
        Pagination is not yet implemented
        """
        return self.schemas.values()


personalize_backends = BackendDict(PersonalizeBackend, "personalize")
