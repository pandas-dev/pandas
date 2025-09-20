import json

from moto.core.common_types import TYPE_RESPONSE
from moto.core.responses import BaseResponse

from .models import ForecastBackend, forecast_backends


class ForecastResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="forecast")

    @property
    def forecast_backend(self) -> ForecastBackend:
        return forecast_backends[self.current_account][self.region]

    def create_dataset_group(self) -> TYPE_RESPONSE:
        dataset_group_name = self._get_param("DatasetGroupName")
        domain = self._get_param("Domain")
        dataset_arns = self._get_param("DatasetArns")
        tags = self._get_param("Tags")

        dataset_group = self.forecast_backend.create_dataset_group(
            dataset_group_name=dataset_group_name,
            domain=domain,
            dataset_arns=dataset_arns,
            tags=tags,
        )
        response = {"DatasetGroupArn": dataset_group.arn}
        return 200, {}, json.dumps(response)

    def describe_dataset_group(self) -> TYPE_RESPONSE:
        dataset_group_arn = self._get_param("DatasetGroupArn")

        dataset_group = self.forecast_backend.describe_dataset_group(
            dataset_group_arn=dataset_group_arn
        )
        response = {
            "CreationTime": dataset_group.creation_date,
            "DatasetArns": dataset_group.dataset_arns,
            "DatasetGroupArn": dataset_group.arn,
            "DatasetGroupName": dataset_group.dataset_group_name,
            "Domain": dataset_group.domain,
            "LastModificationTime": dataset_group.modified_date,
            "Status": "ACTIVE",
        }
        return 200, {}, json.dumps(response)

    def delete_dataset_group(self) -> TYPE_RESPONSE:
        dataset_group_arn = self._get_param("DatasetGroupArn")
        self.forecast_backend.delete_dataset_group(dataset_group_arn)
        return 200, {}, ""

    def update_dataset_group(self) -> TYPE_RESPONSE:
        dataset_group_arn = self._get_param("DatasetGroupArn")
        dataset_arns = self._get_param("DatasetArns")
        self.forecast_backend.update_dataset_group(dataset_group_arn, dataset_arns)
        return 200, {}, ""

    def list_dataset_groups(self) -> TYPE_RESPONSE:
        list_all = sorted(
            [
                {
                    "DatasetGroupArn": dsg.arn,
                    "DatasetGroupName": dsg.dataset_group_name,
                    "CreationTime": dsg.creation_date,
                    "LastModificationTime": dsg.creation_date,
                }
                for dsg in self.forecast_backend.list_dataset_groups()
            ],
            key=lambda x: x["LastModificationTime"],
            reverse=True,
        )
        response = {"DatasetGroups": list_all}
        return 200, {}, json.dumps(response)
