"""Handles incoming servicecatalog requests, invokes methods, returns responses."""

import json

from moto.core.responses import BaseResponse

from .models import ServiceCatalogBackend, servicecatalog_backends


class ServiceCatalogResponse(BaseResponse):
    """Handler for ServiceCatalog requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="servicecatalog")

    @property
    def servicecatalog_backend(self) -> ServiceCatalogBackend:
        """Return backend instance specific for this region."""
        # TODO
        # servicecatalog_backends is not yet typed
        # Please modify moto/backends.py to add the appropriate type annotations for this service
        return servicecatalog_backends[self.current_account][self.region]

    def list_portfolio_access(self) -> str:
        params = json.loads(self.body)
        accept_language = params.get("AcceptLanguage")
        portfolio_id = params.get("PortfolioId")
        organization_parent_id = params.get("OrganizationParentId")
        page_token = params.get("PageToken")
        page_size = params.get("PageSize")

        account_id_objects, next_page_token = (
            self.servicecatalog_backend.list_portfolio_access(
                accept_language=accept_language,
                portfolio_id=portfolio_id,
                organization_parent_id=organization_parent_id,
                page_token=page_token,
                page_size=page_size,
            )
        )

        account_ids = [obj["account_id"] for obj in account_id_objects]

        return json.dumps({"AccountIds": account_ids, "NextPageToken": next_page_token})

    def delete_portfolio(self) -> str:
        params = json.loads(self.body)
        accept_language = params.get("AcceptLanguage")
        id = params.get("Id")

        self.servicecatalog_backend.delete_portfolio(
            accept_language=accept_language,
            id=id,
        )

        return json.dumps({})

    def delete_portfolio_share(self) -> str:
        params = json.loads(self.body)
        accept_language = params.get("AcceptLanguage")
        portfolio_id = params.get("PortfolioId")
        account_id = params.get("AccountId")
        organization_node = params.get("OrganizationNode")

        portfolio_share_token = self.servicecatalog_backend.delete_portfolio_share(
            accept_language=accept_language,
            portfolio_id=portfolio_id,
            account_id=account_id,
            organization_node=organization_node,
        )

        response = {}
        if portfolio_share_token:
            response["PortfolioShareToken"] = portfolio_share_token

        return json.dumps(response)

    def create_portfolio(self) -> str:
        params = json.loads(self.body)

        accept_language = params.get("AcceptLanguage")
        display_name = params.get("DisplayName")
        description = params.get("Description")
        provider_name = params.get("ProviderName")
        tags = params.get("Tags")
        idempotency_token = params.get("IdempotencyToken")

        portfolio_detail, tags = self.servicecatalog_backend.create_portfolio(
            accept_language=accept_language,
            display_name=display_name,
            description=description,
            provider_name=provider_name,
            tags=tags,
            idempotency_token=idempotency_token,
        )

        return json.dumps({"PortfolioDetail": portfolio_detail, "Tags": tags})

    def create_portfolio_share(self) -> str:
        params = json.loads(self.body)
        accept_language = params.get("AcceptLanguage")
        portfolio_id = params.get("PortfolioId")
        account_id = params.get("AccountId")
        organization_node = params.get("OrganizationNode")
        share_tag_options = params.get("ShareTagOptions", False)
        share_principals = params.get("SharePrincipals", False)

        portfolio_share_token = self.servicecatalog_backend.create_portfolio_share(
            accept_language=accept_language,
            portfolio_id=portfolio_id,
            account_id=account_id,
            organization_node=organization_node,
            share_tag_options=share_tag_options,
            share_principals=share_principals,
        )

        response = {}
        if portfolio_share_token:
            response["PortfolioShareToken"] = portfolio_share_token

        return json.dumps(response)

    def list_portfolios(self) -> str:
        params = self._get_params()
        accept_language = params.get("AcceptLanguage")
        page_token = params.get("PageToken")
        page_size = params.get("PageSize")

        portfolio_details, next_page_token = (
            self.servicecatalog_backend.list_portfolios(
                accept_language=accept_language,
                page_token=page_token,
                page_size=page_size,
            )
        )

        response = {
            "PortfolioDetails": portfolio_details,
            "NextPageToken": next_page_token,
        }
        return json.dumps(response)

    def describe_portfolio_shares(self) -> str:
        params = json.loads(self.body)
        portfolio_id = params.get("PortfolioId")
        type = params.get("Type")
        page_token = params.get("PageToken")
        page_size = params.get("PageSize")

        next_page_token, portfolio_share_details = (
            self.servicecatalog_backend.describe_portfolio_shares(
                portfolio_id=portfolio_id,
                type=type,
                page_token=page_token,
                page_size=page_size,
            )
        )

        response = {
            "NextPageToken": next_page_token,
            "PortfolioShareDetails": portfolio_share_details,
        }

        return json.dumps(response)

    def describe_portfolio(self) -> str:
        """Handle describe_portfolio request."""
        params = json.loads(self.body)
        accept_language = params.get("AcceptLanguage")
        id = params.get("Id")

        portfolio_detail, tags, tag_options, budgets = (
            self.servicecatalog_backend.describe_portfolio(
                accept_language=accept_language,
                id=id,
            )
        )

        response = {
            "PortfolioDetail": portfolio_detail,
            "Tags": tags,
            "TagOptions": tag_options,
            "Budgets": budgets,
        }

        return json.dumps(response)

    def create_product(self) -> str:
        params = json.loads(self.body)
        name = params.get("Name")
        owner = params.get("Owner")
        description = params.get("Description")
        distributor = params.get("Distributor")
        support_description = params.get("SupportDescription")
        support_email = params.get("SupportEmail")
        support_url = params.get("SupportUrl")
        product_type = params.get("ProductType")
        tags = params.get("Tags")
        provisioning_artifact_parameters = params.get("ProvisioningArtifactParameters")
        idempotency_token = params.get("IdempotencyToken")
        source_connection = params.get("SourceConnection")
        accept_language = params.get("AcceptLanguage")

        response = self.servicecatalog_backend.create_product(
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
            idempotency_token=idempotency_token,
            source_connection=source_connection,
            accept_language=accept_language,
        )

        product_view = {
            "ProductViewSummary": response.to_dict(),
            "Status": "AVAILABLE",
            "ProductARN": response.arn,
            "CreatedTime": response.created_time.isoformat(),
            "SourceConnection": response.source_connection,
        }

        return json.dumps(
            {
                "ProductViewDetail": product_view,
                "ProvisioningArtifactDetail": {},
                "Tags": response.tags,
            }
        )

    def describe_product(self) -> str:
        params = json.loads(self.body)
        accept_language = params.get("AcceptLanguage")
        id = params.get("Id")
        name = params.get("Name")

        response = self.servicecatalog_backend.describe_product(
            accept_language=accept_language, id=id, name=name
        )

        result = {
            "ProductViewSummary": response.to_dict(),
            "ProvisioningArtifacts": [],
            "Budgets": [],
            "LaunchPaths": [],
        }

        return json.dumps(result)

    def delete_product(self) -> str:
        params = json.loads(self.body)
        accept_language = params.get("AcceptLanguage")
        id = params.get("Id")

        self.servicecatalog_backend.delete_product(
            accept_language=accept_language,
            id=id,
        )

        return json.dumps({})
