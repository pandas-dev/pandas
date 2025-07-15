"""ServiceQuotasBackend class with methods for supported APIs."""

from typing import Any, Dict, List

from moto.core.base_backend import BackendDict, BaseBackend

from .exceptions import NoSuchResource
from .resources.default_quotas.vpc import VPC_DEFAULT_QUOTAS


class ServiceQuotasBackend(BaseBackend):
    """Implementation of ServiceQuotas APIs."""

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)

    def list_aws_default_service_quotas(
        self, service_code: str
    ) -> List[Dict[str, Any]]:
        """
        The ServiceCodes that are currently implemented are: vpc
        Pagination is not yet implemented.
        """
        if service_code == "vpc":
            return VPC_DEFAULT_QUOTAS
        raise NoSuchResource

    def get_service_quota(self, service_code: str, quota_code: str) -> Dict[str, Any]:
        if service_code == "vpc":
            for quota in VPC_DEFAULT_QUOTAS:
                if quota["QuotaCode"] == quota_code:
                    return quota
        raise NoSuchResource


servicequotas_backends = BackendDict(ServiceQuotasBackend, "service-quotas")
