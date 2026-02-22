import json

from moto.core.responses import BaseResponse

from .models import RDSDataServiceBackend, rdsdata_backends


class RDSDataServiceResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="rds-data")

    @property
    def rdsdata_backend(self) -> RDSDataServiceBackend:
        """Return backend instance specific for this region."""
        return rdsdata_backends[self.current_account][self.region]

    def execute_statement(self) -> str:
        resource_arn = self._get_param("resourceArn")
        sql = self._get_param("sql")
        query_result = self.rdsdata_backend.execute_statement(
            resource_arn=resource_arn, sql=sql
        )
        return json.dumps(query_result.to_json())
