"""CostExplorerBackend class with methods for supported APIs."""

from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import iso_8601_datetime_without_milliseconds
from moto.moto_api._internal import mock_random
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import PARTITION_NAMES, get_partition

from .exceptions import CostCategoryNotFound


def first_day() -> str:
    as_date = (
        datetime.today()
        .replace(day=1)
        .replace(hour=0)
        .replace(minute=0)
        .replace(second=0)
    )
    return iso_8601_datetime_without_milliseconds(as_date)


class CostCategoryDefinition(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        name: str,
        effective_start: Optional[str],
        rule_version: str,
        rules: List[Dict[str, Any]],
        default_value: str,
        split_charge_rules: List[Dict[str, Any]],
    ):
        self.name = name
        self.rule_version = rule_version
        self.rules = rules
        self.default_value = default_value
        self.split_charge_rules = split_charge_rules
        self.arn = f"arn:{get_partition(region_name)}:ce::{account_id}:costcategory/{str(mock_random.uuid4())}"
        self.effective_start: str = effective_start or first_day()

    def update(
        self,
        rule_version: str,
        effective_start: Optional[str],
        rules: List[Dict[str, Any]],
        default_value: str,
        split_charge_rules: List[Dict[str, Any]],
    ) -> None:
        self.rule_version = rule_version
        self.rules = rules
        self.default_value = default_value
        self.split_charge_rules = split_charge_rules
        self.effective_start = effective_start or first_day()

    def to_json(self) -> Dict[str, Any]:
        return {
            "CostCategoryArn": self.arn,
            "Name": self.name,
            "EffectiveStart": self.effective_start,
            "RuleVersion": self.rule_version,
            "Rules": self.rules,
            "DefaultValue": self.default_value,
            "SplitChargeRules": self.split_charge_rules,
        }


class CostExplorerBackend(BaseBackend):
    """Implementation of CostExplorer APIs."""

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.cost_categories: Dict[str, CostCategoryDefinition] = dict()
        self.cost_usage_results_queue: List[Dict[str, Any]] = []
        self.cost_usage_results: Dict[str, Dict[str, Any]] = {}
        self.tagger = TaggingService()

    def create_cost_category_definition(
        self,
        name: str,
        effective_start: Optional[str],
        rule_version: str,
        rules: List[Dict[str, Any]],
        default_value: str,
        split_charge_rules: List[Dict[str, Any]],
        tags: List[Dict[str, str]],
    ) -> Tuple[str, str]:
        """
        The EffectiveOn and ResourceTags-parameters are not yet implemented
        """
        ccd = CostCategoryDefinition(
            account_id=self.account_id,
            region_name=self.region_name,
            name=name,
            effective_start=effective_start,
            rule_version=rule_version,
            rules=rules,
            default_value=default_value,
            split_charge_rules=split_charge_rules,
        )
        self.cost_categories[ccd.arn] = ccd
        self.tag_resource(ccd.arn, tags)
        return ccd.arn, ccd.effective_start

    def describe_cost_category_definition(
        self, cost_category_arn: str
    ) -> CostCategoryDefinition:
        """
        The EffectiveOn-parameter is not yet implemented
        """
        if cost_category_arn not in self.cost_categories:
            ccd_id = cost_category_arn.split("/")[-1]
            raise CostCategoryNotFound(ccd_id)
        return self.cost_categories[cost_category_arn]

    def delete_cost_category_definition(
        self, cost_category_arn: str
    ) -> Tuple[str, str]:
        """
        The EffectiveOn-parameter is not yet implemented
        """
        self.cost_categories.pop(cost_category_arn, None)
        return cost_category_arn, ""

    def update_cost_category_definition(
        self,
        cost_category_arn: str,
        effective_start: Optional[str],
        rule_version: str,
        rules: List[Dict[str, Any]],
        default_value: str,
        split_charge_rules: List[Dict[str, Any]],
    ) -> Tuple[str, str]:
        """
        The EffectiveOn-parameter is not yet implemented
        """
        cost_category = self.describe_cost_category_definition(cost_category_arn)
        cost_category.update(
            rule_version=rule_version,
            rules=rules,
            default_value=default_value,
            split_charge_rules=split_charge_rules,
            effective_start=effective_start,
        )

        return cost_category_arn, cost_category.effective_start

    def list_tags_for_resource(self, resource_arn: str) -> List[Dict[str, str]]:
        return self.tagger.list_tags_for_resource(arn=resource_arn)["Tags"]

    def tag_resource(self, resource_arn: str, tags: List[Dict[str, str]]) -> None:
        self.tagger.tag_resource(resource_arn, tags)

    def untag_resource(self, resource_arn: str, tag_keys: List[str]) -> None:
        self.tagger.untag_resource_using_names(resource_arn, tag_keys)

    def get_cost_and_usage(self, body: str) -> Dict[str, Any]:
        """
        There is no validation yet on any of the input parameters.

        Cost or usage is not tracked by Moto, so this call will return nothing by default.

        You can use a dedicated API to override this, by configuring a queue of expected results.

        A request to `get_cost_and_usage` will take the first result from that queue, and assign it to the provided parameters. Subsequent requests using the same parameters will return the same result. Other requests using different parameters will take the next result from the queue, or return an empty result if the queue is empty.

        Configure this queue by making an HTTP request to `/moto-api/static/ce/cost-and-usage-results`. An example invocation looks like this:

        .. sourcecode:: python

            result = {
                "results": [
                    {
                        "ResultsByTime": [
                            {
                                "TimePeriod": {"Start": "2024-01-01", "End": "2024-01-02"},
                                "Total": {
                                    "BlendedCost": {"Amount": "0.0101516483", "Unit": "USD"}
                                },
                                "Groups": [],
                                "Estimated": False
                            }
                        ],
                        "DimensionValueAttributes": [{"Value": "v", "Attributes": {"a": "b"}}]
                    },
                    {
                        ...
                    },
                ]
            }
            resp = requests.post(
                "http://motoapi.amazonaws.com/moto-api/static/ce/cost-and-usage-results",
                json=expected_results,
            )
            assert resp.status_code == 201

            ce = boto3.client("ce", region_name="us-east-1")
            resp = ce.get_cost_and_usage(...)
        """
        default_result: Dict[str, Any] = {
            "ResultsByTime": [],
            "DimensionValueAttributes": [],
        }
        if body not in self.cost_usage_results and self.cost_usage_results_queue:
            self.cost_usage_results[body] = self.cost_usage_results_queue.pop(0)
        return self.cost_usage_results.get(body, default_result)


ce_backends = BackendDict(
    CostExplorerBackend,
    "ce",
    use_boto3_regions=False,
    additional_regions=PARTITION_NAMES,
)
