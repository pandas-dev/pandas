"""Handles incoming servicequotas requests, invokes methods, returns responses."""

import json

from moto.core.responses import BaseResponse

from .models import ServiceQuotasBackend, servicequotas_backends


class ServiceQuotasResponse(BaseResponse):
    """Handler for ServiceQuotas requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="service-quotas")

    @property
    def backend(self) -> ServiceQuotasBackend:
        """Return backend instance specific for this region."""
        return servicequotas_backends[self.current_account][self.region]

    def list_aws_default_service_quotas(self) -> str:
        params = json.loads(self.body)
        service_code = str(params.get("ServiceCode"))
        quotas = self.backend.list_aws_default_service_quotas(service_code)
        return json.dumps(dict(Quotas=quotas))

    def get_service_quota(self) -> str:
        params = json.loads(self.body)
        service_code = str(params.get("ServiceCode"))
        quota_code = str(params.get("QuotaCode"))
        quota = self.backend.get_service_quota(
            service_code=service_code,
            quota_code=quota_code,
        )
        return json.dumps(dict(Quota=quota))
