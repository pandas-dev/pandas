"""Handles incoming ce requests, invokes methods, returns responses."""

import json

from moto.core.responses import BaseResponse

from .models import CostExplorerBackend, ce_backends


class CostExplorerResponse(BaseResponse):
    """Handler for CostExplorer requests and responses."""

    @property
    def ce_backend(self) -> CostExplorerBackend:
        """Return backend instance specific for this region."""
        return ce_backends[self.current_account][self.partition]

    def create_cost_category_definition(self) -> str:
        params = json.loads(self.body)
        name = params.get("Name")
        rule_version = params.get("RuleVersion")
        rules = params.get("Rules")
        default_value = params.get("DefaultValue")
        split_charge_rules = params.get("SplitChargeRules")
        effective_start = params.get("EffectiveStart")
        tags = params.get("ResourceTags")
        (
            cost_category_arn,
            effective_start,
        ) = self.ce_backend.create_cost_category_definition(
            name=name,
            effective_start=effective_start,
            rule_version=rule_version,
            rules=rules,
            default_value=default_value,
            split_charge_rules=split_charge_rules,
            tags=tags,
        )
        return json.dumps(
            dict(CostCategoryArn=cost_category_arn, EffectiveStart=effective_start)
        )

    def describe_cost_category_definition(self) -> str:
        params = json.loads(self.body)
        cost_category_arn = params.get("CostCategoryArn")
        cost_category = self.ce_backend.describe_cost_category_definition(
            cost_category_arn=cost_category_arn
        )
        return json.dumps(dict(CostCategory=cost_category.to_json()))

    def delete_cost_category_definition(self) -> str:
        params = json.loads(self.body)
        cost_category_arn = params.get("CostCategoryArn")
        (
            cost_category_arn,
            effective_end,
        ) = self.ce_backend.delete_cost_category_definition(
            cost_category_arn=cost_category_arn,
        )
        return json.dumps(
            dict(CostCategoryArn=cost_category_arn, EffectiveEnd=effective_end)
        )

    def update_cost_category_definition(self) -> str:
        params = json.loads(self.body)
        cost_category_arn = params.get("CostCategoryArn")
        effective_start = params.get("EffectiveStart")
        rule_version = params.get("RuleVersion")
        rules = params.get("Rules")
        default_value = params.get("DefaultValue")
        split_charge_rules = params.get("SplitChargeRules")
        (
            cost_category_arn,
            effective_start,
        ) = self.ce_backend.update_cost_category_definition(
            cost_category_arn=cost_category_arn,
            effective_start=effective_start,
            rule_version=rule_version,
            rules=rules,
            default_value=default_value,
            split_charge_rules=split_charge_rules,
        )
        return json.dumps(
            dict(CostCategoryArn=cost_category_arn, EffectiveStart=effective_start)
        )

    def list_tags_for_resource(self) -> str:
        params = json.loads(self.body)
        resource_arn = params.get("ResourceArn")
        tags = self.ce_backend.list_tags_for_resource(resource_arn)
        return json.dumps({"ResourceTags": tags})

    def tag_resource(self) -> str:
        params = json.loads(self.body)
        resource_arn = params.get("ResourceArn")
        tags = params.get("ResourceTags")
        self.ce_backend.tag_resource(resource_arn, tags)
        return json.dumps({})

    def untag_resource(self) -> str:
        params = json.loads(self.body)
        resource_arn = params.get("ResourceArn")
        tag_names = params.get("ResourceTagKeys")
        self.ce_backend.untag_resource(resource_arn, tag_names)
        return json.dumps({})

    def get_cost_and_usage(self) -> str:
        resp = self.ce_backend.get_cost_and_usage(self.body)
        return json.dumps(resp)
