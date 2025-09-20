import re
from datetime import datetime
from typing import Dict, List, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.utils import iso_8601_datetime_without_milliseconds
from moto.utilities.utils import get_partition

from .exceptions import (
    InvalidInputException,
    ResourceAlreadyExistsException,
    ResourceNotFoundException,
    ValidationException,
)


class DatasetGroup:
    accepted_dataset_group_name_format = re.compile(r"^[a-zA-Z][a-z-A-Z0-9_]*")
    accepted_dataset_group_arn_format = re.compile(r"^[a-zA-Z0-9\-\_\.\/\:]+$")
    accepted_dataset_types = [
        "INVENTORY_PLANNING",
        "METRICS",
        "RETAIL",
        "EC2_CAPACITY",
        "CUSTOM",
        "WEB_TRAFFIC",
        "WORK_FORCE",
    ]

    def __init__(
        self,
        account_id: str,
        region_name: str,
        dataset_arns: List[str],
        dataset_group_name: str,
        domain: str,
        tags: Optional[List[Dict[str, str]]] = None,
    ):
        self.creation_date = iso_8601_datetime_without_milliseconds(datetime.now())
        self.modified_date = self.creation_date

        self.arn = f"arn:{get_partition(region_name)}:forecast:{region_name}:{account_id}:dataset-group/{dataset_group_name}"
        self.dataset_arns = dataset_arns if dataset_arns else []
        self.dataset_group_name = dataset_group_name
        self.domain = domain
        self.tags = tags
        self._validate()

    def update(self, dataset_arns: List[str]) -> None:
        self.dataset_arns = dataset_arns
        self.last_modified_date = iso_8601_datetime_without_milliseconds(datetime.now())

    def _validate(self) -> None:
        errors = []

        errors.extend(self._validate_dataset_group_name())
        errors.extend(self._validate_dataset_group_name_len())
        errors.extend(self._validate_dataset_group_domain())

        if errors:
            err_count = len(errors)
            message = str(err_count) + " validation error"
            message += "s" if err_count > 1 else ""
            message += " detected: "
            message += "; ".join(errors)
            raise ValidationException(message)

    def _validate_dataset_group_name(self) -> List[str]:
        errors = []
        if not re.match(
            self.accepted_dataset_group_name_format, self.dataset_group_name
        ):
            errors.append(
                "Value '"
                + self.dataset_group_name
                + "' at 'datasetGroupName' failed to satisfy constraint: Member must satisfy regular expression pattern "
                + self.accepted_dataset_group_name_format.pattern
            )
        return errors

    def _validate_dataset_group_name_len(self) -> List[str]:
        errors = []
        if len(self.dataset_group_name) >= 64:
            errors.append(
                "Value '"
                + self.dataset_group_name
                + "' at 'datasetGroupName' failed to satisfy constraint: Member must have length less than or equal to 63"
            )
        return errors

    def _validate_dataset_group_domain(self) -> List[str]:
        errors = []
        if self.domain not in self.accepted_dataset_types:
            errors.append(
                "Value '"
                + self.domain
                + "' at 'domain' failed to satisfy constraint: Member must satisfy enum value set "
                + str(self.accepted_dataset_types)
            )
        return errors


class ForecastBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.dataset_groups: Dict[str, DatasetGroup] = {}
        self.datasets: Dict[str, str] = {}

    def create_dataset_group(
        self,
        dataset_group_name: str,
        domain: str,
        dataset_arns: List[str],
        tags: List[Dict[str, str]],
    ) -> DatasetGroup:
        dataset_group = DatasetGroup(
            account_id=self.account_id,
            region_name=self.region_name,
            dataset_group_name=dataset_group_name,
            domain=domain,
            dataset_arns=dataset_arns,
            tags=tags,
        )

        if dataset_arns:
            for dataset_arn in dataset_arns:
                if dataset_arn not in self.datasets:
                    raise InvalidInputException(
                        "Dataset arns: [" + dataset_arn + "] are not found"
                    )

        if self.dataset_groups.get(dataset_group.arn):
            raise ResourceAlreadyExistsException(
                "A dataset group already exists with the arn: " + dataset_group.arn
            )

        self.dataset_groups[dataset_group.arn] = dataset_group
        return dataset_group

    def describe_dataset_group(self, dataset_group_arn: str) -> DatasetGroup:
        try:
            return self.dataset_groups[dataset_group_arn]
        except KeyError:
            raise ResourceNotFoundException("No resource found " + dataset_group_arn)

    def delete_dataset_group(self, dataset_group_arn: str) -> None:
        try:
            del self.dataset_groups[dataset_group_arn]
        except KeyError:
            raise ResourceNotFoundException("No resource found " + dataset_group_arn)

    def update_dataset_group(
        self, dataset_group_arn: str, dataset_arns: List[str]
    ) -> None:
        try:
            dsg = self.dataset_groups[dataset_group_arn]
        except KeyError:
            raise ResourceNotFoundException("No resource found " + dataset_group_arn)

        for dataset_arn in dataset_arns:
            if dataset_arn not in dsg.dataset_arns:
                raise InvalidInputException(
                    "Dataset arns: [" + dataset_arn + "] are not found"
                )

        dsg.update(dataset_arns)

    def list_dataset_groups(self) -> List[DatasetGroup]:
        return [v for (_, v) in self.dataset_groups.items()]


forecast_backends = BackendDict(ForecastBackend, "forecast")
