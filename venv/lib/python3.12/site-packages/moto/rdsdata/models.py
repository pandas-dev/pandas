from typing import Any, Dict, List, Optional, Tuple

from moto.core.base_backend import BackendDict, BaseBackend


class QueryResults:
    def __init__(
        self,
        records: Optional[List[List[Dict[str, Any]]]] = None,
        column_metadata: Optional[List[Dict[str, Any]]] = None,
        number_of_records_updated: Optional[int] = None,
        generated_fields: Optional[List[Dict[str, Any]]] = None,
        formatted_records: Optional[str] = None,
    ):
        self.records = records
        self.column_metadata = column_metadata
        self.number_of_records_updated = number_of_records_updated
        self.generated_fields = generated_fields
        self.formatted_records = formatted_records

    def to_json(self) -> Dict[str, Any]:
        return {
            "records": self.records,
            "columnMetadata": self.column_metadata,
            "numberOfRecordsUpdated": self.number_of_records_updated,
            "generatedFields": self.generated_fields,
            "formattedRecords": self.formatted_records,
        }


class RDSDataServiceBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.results_queue: List[QueryResults] = []
        self.sql_results: Dict[Tuple[str, str], QueryResults] = dict()

    def execute_statement(self, resource_arn: str, sql: str) -> QueryResults:
        """
        There is no validation yet on any of the input parameters.

        SQL statements are not executed by Moto, so this call will always return 0 records by default.

        You can use a dedicated API to override this, by configuring a queue of expected results.

        A request to `execute_statement` will take the first result from that queue, and assign it to the provided SQL-query. Subsequent requests using the same SQL-query will return the same result. Other requests using a different SQL-query will take the next result from the queue, or return an empty result if the queue is empty.

        Configure this queue by making an HTTP request to `/moto-api/static/rds-data/statement-results`. An example invocation looks like this:

        .. sourcecode:: python

            expected_results = {
                "account_id": "123456789012",  # This is the default - can be omitted
                "region": "us-east-1",  # This is the default - can be omitted
                "results": [
                    {
                        "records": [...],
                        "columnMetadata": [...],
                        "numberOfRecordsUpdated": 42,
                        "generatedFields": [...],
                        "formattedRecords": "some json"
                    },
                    # other results as required
                ],
            }
            resp = requests.post(
                "http://motoapi.amazonaws.com/moto-api/static/rds-data/statement-results",
                json=expected_results,
            )
            assert resp.status_code == 201

            rdsdata = boto3.client("rds-data", region_name="us-east-1")
            resp = rdsdata.execute_statement(resourceArn="not applicable", secretArn="not applicable", sql="SELECT some FROM thing")

        """

        if (resource_arn, sql) in self.sql_results:
            return self.sql_results[(resource_arn, sql)]
        elif self.results_queue:
            self.sql_results[(resource_arn, sql)] = self.results_queue.pop()
            return self.sql_results[(resource_arn, sql)]
        else:
            return QueryResults(records=[])


rdsdata_backends = BackendDict(RDSDataServiceBackend, "rds-data")
