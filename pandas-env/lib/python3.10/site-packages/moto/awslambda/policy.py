import json
from typing import TYPE_CHECKING, Any, Callable, Dict, List, Optional, TypeVar

from moto.awslambda.exceptions import (
    GenericResourcNotFound,
    PreconditionFailedException,
    UnknownPolicyException,
)
from moto.moto_api._internal import mock_random

if TYPE_CHECKING:
    from .models import LambdaFunction


TYPE_IDENTITY = TypeVar("TYPE_IDENTITY")


class Policy:
    def __init__(self, parent: "LambdaFunction"):
        self.revision = str(mock_random.uuid4())
        self.statements: List[Dict[str, Any]] = []
        self.parent = parent

    def wire_format(self) -> str:
        p = self.get_policy()
        p["Policy"] = json.dumps(p["Policy"])
        return json.dumps(p)

    def get_policy(self) -> Dict[str, Any]:
        if not self.statements:
            raise GenericResourcNotFound()
        return {
            "Policy": {
                "Version": "2012-10-17",
                "Id": "default",
                "Statement": self.statements,
            },
            "RevisionId": self.revision,
        }

    # adds the raw JSON statement to the policy
    def add_statement(
        self, raw: str, qualifier: Optional[str] = None
    ) -> Dict[str, Any]:
        policy = json.loads(raw, object_hook=self.decode_policy)
        if len(policy.revision) > 0 and self.revision != policy.revision:
            raise PreconditionFailedException(
                "The RevisionId provided does not match the latest RevisionId"
                " for the Lambda function or alias. Call the GetFunction or the GetAlias API to retrieve"
                " the latest RevisionId for your resource."
            )
        # Remove #LATEST from the Resource (Lambda ARN)
        if policy.statements[0].get("Resource", "").endswith("$LATEST"):
            policy.statements[0]["Resource"] = policy.statements[0]["Resource"][0:-8]
        if qualifier:
            policy.statements[0]["Resource"] = (
                policy.statements[0]["Resource"] + ":" + qualifier
            )
        self.statements.append(policy.statements[0])
        self.revision = str(mock_random.uuid4())
        return policy.statements[0]

    # removes the statement that matches 'sid' from the policy
    def del_statement(self, sid: str, revision: str = "") -> None:
        if len(revision) > 0 and self.revision != revision:
            raise PreconditionFailedException(
                "The RevisionId provided does not match the latest RevisionId"
                " for the Lambda function or alias. Call the GetFunction or the GetAlias API to retrieve"
                " the latest RevisionId for your resource."
            )
        for statement in self.statements:
            if "Sid" in statement and statement["Sid"] == sid:
                self.statements.remove(statement)
                break
        else:
            raise UnknownPolicyException()

    # converts AddPermission request to PolicyStatement
    # https://docs.aws.amazon.com/lambda/latest/dg/API_AddPermission.html
    def decode_policy(self, obj: Dict[str, Any]) -> "Policy":
        policy = Policy(self.parent)
        policy.revision = obj.get("RevisionId", "")

        # set some default values if these keys are not set
        self.ensure_set(obj, "Effect", "Allow")
        self.ensure_set(obj, "Resource", self.parent.function_arn + ":$LATEST")
        self.ensure_set(obj, "StatementId", str(mock_random.uuid4()))

        # transform field names and values
        self.transform_property(obj, "StatementId", "Sid", self.nop_formatter)
        self.transform_property(obj, "Principal", "Principal", self.principal_formatter)

        self.transform_property(
            obj, "SourceArn", "SourceArn", self.source_arn_formatter
        )
        self.transform_property(
            obj, "SourceAccount", "SourceAccount", self.source_account_formatter
        )
        self.transform_property(
            obj, "PrincipalOrgID", "Condition", self.principal_org_id_formatter
        )

        # remove RevisionId and EventSourceToken if they are set
        self.remove_if_set(obj, ["RevisionId", "EventSourceToken"])

        # merge conditional statements into a single map under the Condition key
        self.condition_merge(obj)

        # append resulting statement to policy.statements
        policy.statements.append(obj)

        return policy

    def nop_formatter(self, obj: TYPE_IDENTITY) -> TYPE_IDENTITY:
        return obj

    def ensure_set(self, obj: Dict[str, Any], key: str, value: Any) -> None:
        if key not in obj:
            obj[key] = value

    def principal_formatter(self, obj: Dict[str, Any]) -> Dict[str, Any]:
        if isinstance(obj, str):
            if obj.endswith(".amazonaws.com"):
                return {"Service": obj}
            if obj.endswith(":root"):
                return {"AWS": obj}
        return obj

    def source_account_formatter(
        self, obj: TYPE_IDENTITY
    ) -> Dict[str, Dict[str, TYPE_IDENTITY]]:
        return {"StringEquals": {"AWS:SourceAccount": obj}}

    def source_arn_formatter(
        self, obj: TYPE_IDENTITY
    ) -> Dict[str, Dict[str, TYPE_IDENTITY]]:
        return {"ArnLike": {"AWS:SourceArn": obj}}

    def principal_org_id_formatter(
        self, obj: TYPE_IDENTITY
    ) -> Dict[str, Dict[str, TYPE_IDENTITY]]:
        return {"StringEquals": {"aws:PrincipalOrgID": obj}}

    def transform_property(
        self,
        obj: Dict[str, Any],
        old_name: str,
        new_name: str,
        formatter: Callable[..., Any],
    ) -> None:
        if old_name in obj:
            obj[new_name] = formatter(obj[old_name])
            if new_name != old_name:
                del obj[old_name]

    def remove_if_set(self, obj: Dict[str, Any], keys: List[str]) -> None:
        for key in keys:
            if key in obj:
                del obj[key]

    def condition_merge(self, obj: Dict[str, Any]) -> None:
        if "SourceArn" in obj:
            if "Condition" not in obj:
                obj["Condition"] = {}
            obj["Condition"].update(obj["SourceArn"])
            del obj["SourceArn"]

        if "SourceAccount" in obj:
            if "Condition" not in obj:
                obj["Condition"] = {}
            obj["Condition"].update(obj["SourceAccount"])
            del obj["SourceAccount"]
