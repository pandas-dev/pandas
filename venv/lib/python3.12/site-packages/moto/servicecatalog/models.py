"""ServiceCatalogBackend class with methods for supported APIs."""

import uuid
from datetime import datetime
from typing import Any, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.utilities.paginator import paginate
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import get_partition

from .exceptions import (
    InvalidParametersException,
    ResourceNotFoundException,
)

PAGINATION_MODEL = {
    "list_portfolio_access": {
        "input_token": "page_token",
        "limit_key": "page_size",
        "limit_default": 20,
        "unique_attribute": "account_id",
    }
}


class Portfolio(BaseModel):
    """Portfolio resource."""

    def __init__(
        self,
        portfolio_id: str,
        display_name: str,
        description: str,
        provider_name: str,
        tags: list[dict[str, str]],
        region_name: str,
        account_id: str,
        backend: "ServiceCatalogBackend",
    ) -> None:
        self.id = portfolio_id
        self.display_name = display_name
        self.description = description
        self.provider_name = provider_name
        self.created_time = datetime.now()
        self.tags = tags
        self.region_name = region_name
        self.account_id = account_id
        self.backend = backend
        if tags:
            self.backend._tag_resource(self.arn, tags)

    @property
    def arn(self) -> str:
        return f"arn:{get_partition(self.region_name)}:catalog:{self.region_name}:{self.account_id}:portfolio/{self.id}"

    def to_dict(self) -> dict[str, str]:
        return {
            "Id": self.id,
            "ARN": self.arn,
            "DisplayName": self.display_name,
            "Description": self.description,
            "CreatedTime": self.created_time.isoformat(),
            "ProviderName": self.provider_name,
        }


class Product(BaseModel):
    """Service Catalog Product resource."""

    def __init__(
        self,
        name: str,
        owner: str,
        description: Optional[str],
        distributor: Optional[str],
        support_description: Optional[str],
        support_email: Optional[str],
        support_url: Optional[str],
        product_type: str,
        tags: Optional[list[dict[str, str]]],
        provisioning_artifact_parameters: Optional[dict[str, Any]],
        source_connection: Optional[dict[str, Any]],
        accept_language: Optional[str],
        region_name: str,
        account_id: str,
        backend: "ServiceCatalogBackend",
    ) -> None:
        self.id = str(uuid.uuid4())
        self.name = name
        self.owner = owner
        self.description = description
        self.distributor = distributor
        self.support_description = support_description
        self.support_email = support_email
        self.support_url = support_url
        self.product_type = product_type
        self.tags = tags or []
        self.provisioning_artifact_parameters = provisioning_artifact_parameters or {}
        self.source_connection = source_connection or {}
        self.accept_language = accept_language
        self.created_time = datetime.now()
        self.region_name = region_name
        self.account_id = account_id
        self.backend = backend
        if tags:
            self.backend._tag_resource(self.arn, tags)

    @property
    def arn(self) -> str:
        return f"arn:{get_partition(self.region_name)}:catalog:{self.region_name}:{self.account_id}:product/{self.id}"

    def to_dict(self) -> dict[str, Any]:
        product_view_summary = {
            "Id": self.id,  # ?
            "ProductId": self.id,  # ?
            "Name": self.name,
            "Owner": self.owner,
            "Type": self.product_type,
            "Distributor": self.distributor,
            "HasDefaultPath": False,
            "ShortDescription": self.description,
            "SupportDescription": self.support_description,
            "SupportEmail": self.support_email,
            "SupportUrl": self.support_url,
        }

        return product_view_summary


class ServiceCatalogBackend(BaseBackend):
    """Implementation of ServiceCatalog APIs."""

    def __init__(self, region_name: str, account_id: str) -> None:
        super().__init__(region_name, account_id)
        self.portfolio_access: dict[str, list[str]] = {}
        self.portfolios: dict[str, Portfolio] = {}
        self.idempotency_tokens: dict[str, str] = {}
        self.portfolio_share_tokens: dict[str, list[str]] = {}
        self.products: dict[str, Product] = {}
        self.tagger = TaggingService()

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_portfolio_access(
        self,
        accept_language: Optional[str],
        portfolio_id: str,
        organization_parent_id: Optional[str],
        page_token: Optional[str],
        page_size: Optional[int] = None,
    ) -> list[dict[str, str]]:
        # TODO: Implement organization_parent_id and accept_language

        account_ids = self.portfolio_access.get(portfolio_id, [])
        return [{"account_id": account_id} for account_id in account_ids]

    def delete_portfolio(self, accept_language: Optional[str], id: str) -> None:
        # TODO: Implement accept_language

        if id in self.portfolio_access:
            del self.portfolio_access[id]

        if id in self.portfolios:
            del self.portfolios[id]

        return None

    def delete_portfolio_share(
        self,
        accept_language: Optional[str],
        portfolio_id: str,
        account_id: Optional[str],
        organization_node: Optional[dict[str, str]],
    ) -> Optional[str]:
        # TODO: Implement accept_language

        if (
            portfolio_id in self.portfolio_access
            and account_id in self.portfolio_access[portfolio_id]
        ):
            self.portfolio_access[portfolio_id].remove(account_id)

        portfolio_share_token = None
        if organization_node:
            org_type = organization_node.get("Type", "")
            org_value = organization_node.get("Value", "")

            # Arbitrary naming for the portfolio share token
            portfolio_share_token = f"share-{portfolio_id}-{org_type}-{org_value}"

            if portfolio_id in self.portfolio_share_tokens:
                tokens_to_remove = []

                for token in self.portfolio_share_tokens[portfolio_id]:
                    if org_type in token and org_value in token:
                        tokens_to_remove.append(token)

                for token in tokens_to_remove:
                    if token in self.portfolio_share_tokens[portfolio_id]:
                        self.portfolio_share_tokens[portfolio_id].remove(token)

        return portfolio_share_token

    def create_portfolio(
        self,
        accept_language: Optional[str],
        display_name: str,
        description: Optional[str],
        provider_name: str,
        tags: Optional[list[dict[str, str]]],
        idempotency_token: Optional[str],
    ) -> tuple[dict[str, str], list[dict[str, str]]]:
        # TODO: Implement accept_language

        if idempotency_token and idempotency_token in self.idempotency_tokens:
            portfolio_id = self.idempotency_tokens[idempotency_token]
            portfolio = self.portfolios[portfolio_id]
            return portfolio.to_dict(), portfolio.tags

        portfolio_id = str(uuid.uuid4())

        portfolio = Portfolio(
            portfolio_id=portfolio_id,
            display_name=display_name,
            description=description or "",
            provider_name=provider_name,
            tags=tags or [],
            region_name=self.region_name,
            account_id=self.account_id,
            backend=self,
        )

        self.portfolios[portfolio_id] = portfolio

        if idempotency_token:
            self.idempotency_tokens[idempotency_token] = portfolio_id

        self.portfolio_access[portfolio_id] = []

        return portfolio.to_dict(), portfolio.tags

    def create_portfolio_share(
        self,
        accept_language: Optional[str],
        portfolio_id: str,
        account_id: Optional[str],
        organization_node: Optional[dict[str, str]],
        share_tag_options: bool,
        share_principals: bool,
    ) -> Optional[str]:
        # TODO: Implement accept_language

        if portfolio_id not in self.portfolios:
            return None

        if account_id:
            if portfolio_id not in self.portfolio_access:
                self.portfolio_access[portfolio_id] = []

            if account_id not in self.portfolio_access[portfolio_id]:
                self.portfolio_access[portfolio_id].append(account_id)

            return None

        portfolio_share_token = None
        if organization_node:
            org_type = organization_node.get("Type", "")
            org_value = organization_node.get("Value", "")

            # Arbitrary org for the portfolio share token
            portfolio_share_token = f"share-{portfolio_id}-{org_type}-{org_value}"

            if share_tag_options:
                portfolio_share_token += "-tags"
            if share_principals:
                portfolio_share_token += "-principals"

            if portfolio_id not in self.portfolio_share_tokens:
                self.portfolio_share_tokens[portfolio_id] = []

            if portfolio_share_token not in self.portfolio_share_tokens[portfolio_id]:
                self.portfolio_share_tokens[portfolio_id].append(portfolio_share_token)

        return portfolio_share_token

    def list_portfolios(
        self,
        accept_language: Optional[str] = None,
        page_token: Optional[str] = None,
        page_size: Optional[int] = None,
    ) -> tuple[list[dict[str, str]], Optional[str]]:
        """TODO: Implement pagination and accept_language"""
        portfolio_details = [
            portfolio.to_dict() for portfolio in self.portfolios.values()
        ]
        return portfolio_details, None

    def describe_portfolio_shares(
        self,
        portfolio_id: str,
        type: str,
        page_token: Optional[str] = None,
        page_size: Optional[int] = None,
    ) -> tuple[Optional[str], list[dict[str, Any]]]:
        """TODO: Implement pagination"""

        portfolio_share_details = []

        if portfolio_id not in self.portfolios:
            return None, []

        if type == "ACCOUNT":
            account_ids = self.portfolio_access.get(portfolio_id, [])
            for account_id in account_ids:
                portfolio_share_details.append(
                    {
                        "PrincipalId": account_id,
                        "Type": "ACCOUNT",
                        "Accepted": True,
                        "ShareTagOptions": False,
                        "SharePrincipals": False,
                    }
                )

        elif type == "ORGANIZATION":
            tokens = self.portfolio_share_tokens.get(portfolio_id, [])

            for token in tokens:
                if "ORGANIZATION" in token and "-o-" in token:
                    org_id_start = token.find("-ORGANIZATION-") + len("-ORGANIZATION-")

                    if "-tags" in token:
                        org_id_end = token.find("-tags", org_id_start)
                    elif "-principals" in token:
                        org_id_end = token.find("-principals", org_id_start)
                    else:
                        org_id_end = len(token)

                    org_id = token[org_id_start:org_id_end]

                    portfolio_share_details.append(
                        {
                            "PrincipalId": org_id,
                            "Type": "ORGANIZATION",
                            "Accepted": True,
                            "ShareTagOptions": "tags" in token,
                            "SharePrincipals": "principals" in token,
                        }
                    )

        return None, portfolio_share_details

    def describe_portfolio(
        self, accept_language: Optional[str], id: str
    ) -> tuple[
        dict[str, Any], list[dict[str, str]], list[dict[str, Any]], list[dict[str, str]]
    ]:
        # TODO: Implement accept_language

        if id not in self.portfolios:
            return {}, [], [], []

        portfolio = self.portfolios[id]
        portfolio_detail = portfolio.to_dict()

        tags = portfolio.tags

        tag_options: list[dict[str, Any]] = []
        budgets: list[dict[str, Any]] = []

        return portfolio_detail, tags, tag_options, budgets

    def create_product(
        self,
        name: str,
        owner: str,
        description: Optional[str],
        distributor: Optional[str],
        support_description: Optional[str],
        support_email: Optional[str],
        support_url: Optional[str],
        product_type: str,
        tags: list[dict[str, str]],
        provisioning_artifact_parameters: Optional[dict[str, Any]],
        idempotency_token: Optional[str] = None,
        source_connection: Optional[dict[str, Any]] = None,
        accept_language: Optional[str] = None,
    ) -> Product:
        token = idempotency_token or str(uuid.uuid4())
        existing_id = self.idempotency_tokens.get(token)
        if existing_id:
            return self.products[existing_id]

        product = Product(
            name=name,
            owner=owner,
            description=description,
            distributor=distributor,
            support_description=support_description,
            support_email=support_email,
            support_url=support_url,
            product_type=product_type,
            tags=tags,
            provisioning_artifact_parameters=provisioning_artifact_parameters,
            source_connection=source_connection,
            accept_language=accept_language,
            region_name=self.region_name,
            account_id=self.account_id,
            backend=self,
        )

        self.products[product.id] = product
        self.idempotency_tokens[token] = product.id

        return product

    def describe_product(
        self, accept_language: Optional[str], id: str, name: str
    ) -> Product:
        if not id and not name:
            raise InvalidParametersException("Either Id or Name must be specified.")

        if id:
            product = self.products.get(id)
            if not product:
                raise ResourceNotFoundException(f"Product with Id '{id}' not found.")
        else:
            product = self.lookup_by_name(name)
            if not product:
                raise ResourceNotFoundException(
                    f"Product with Name '{name}' not found."
                )

        return product

    def lookup_by_name(self, name: str) -> Optional[Product]:
        for product in self.products.values():
            if name == product.name:
                return product

        return None

    def delete_product(self, accept_language: Optional[str], id: str) -> None:
        if id in self.products:
            del self.products[id]

        return None

    def _tag_resource(self, resource_arn: str, tags: list[dict[str, str]]) -> None:
        self.tagger.tag_resource(resource_arn, tags)

    def _list_tags_for_resource(self, resource_arn: str) -> list[dict[str, str]]:
        if self.tagger.has_tags(resource_arn):
            return self.tagger.list_tags_for_resource(resource_arn)["Tags"]
        return []


servicecatalog_backends = BackendDict(ServiceCatalogBackend, "servicecatalog")
