import json

from moto.core.common_types import TYPE_RESPONSE
from moto.core.responses import BaseResponse

from .models import RedshiftDataAPIServiceBackend, redshiftdata_backends


class RedshiftDataAPIServiceResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="redshift-data")

    @property
    def redshiftdata_backend(self) -> RedshiftDataAPIServiceBackend:
        return redshiftdata_backends[self.current_account][self.region]

    def cancel_statement(self) -> TYPE_RESPONSE:
        statement_id = self._get_param("Id")
        self.redshiftdata_backend.cancel_statement(statement_id=statement_id)
        return 200, {}, json.dumps({"Status": True})

    def describe_statement(self) -> TYPE_RESPONSE:
        statement_id = self._get_param("Id")
        statement = self.redshiftdata_backend.describe_statement(
            statement_id=statement_id
        )
        return 200, {}, json.dumps(dict(statement))

    def execute_statement(self) -> TYPE_RESPONSE:
        cluster_identifier = self._get_param("ClusterIdentifier")
        database = self._get_param("Database")
        db_user = self._get_param("DbUser")
        parameters = self._get_param("Parameters")
        secret_arn = self._get_param("SecretArn")
        sql = self._get_param("Sql")
        statement = self.redshiftdata_backend.execute_statement(
            cluster_identifier=cluster_identifier,
            database=database,
            db_user=db_user,
            parameters=parameters,
            secret_arn=secret_arn,
            sql=sql,
        )

        return (
            200,
            {},
            json.dumps(
                {
                    "ClusterIdentifier": statement.cluster_identifier,
                    "CreatedAt": statement.created_at,
                    "Database": statement.database,
                    "DbUser": statement.db_user,
                    "Id": statement.id,
                    "SecretArn": statement.secret_arn,
                }
            ),
        )

    def get_statement_result(self) -> TYPE_RESPONSE:
        statement_id = self._get_param("Id")
        statement_result = self.redshiftdata_backend.get_statement_result(
            statement_id=statement_id
        )

        return 200, {}, json.dumps(dict(statement_result))
