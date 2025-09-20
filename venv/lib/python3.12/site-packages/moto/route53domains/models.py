from datetime import datetime, timedelta, timezone
from typing import Any, Dict, List, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.route53 import route53_backends
from moto.route53.models import Route53Backend
from moto.utilities.paginator import paginate
from moto.utilities.utils import PARTITION_NAMES

from .exceptions import (
    DomainLimitExceededException,
    DuplicateRequestException,
    InvalidInputException,
)
from .validators import (
    DOMAIN_OPERATION_STATUSES,
    DOMAIN_OPERATION_TYPES,
    DomainFilterField,
    DomainsFilter,
    DomainSortOrder,
    DomainsSortCondition,
    NameServer,
    Route53Domain,
    Route53DomainsContactDetail,
    Route53DomainsOperation,
    ValidationException,
)


class Route53DomainsBackend(BaseBackend):
    """Implementation of Route53Domains APIs."""

    DEFAULT_MAX_DOMAINS_COUNT = 20
    PAGINATION_MODEL = {
        "list_domains": {
            "input_token": "marker",
            "limit_key": "max_items",
            "limit_default": 20,
            "unique_attribute": "domain_name",
        },
        "list_operations": {
            "input_token": "marker",
            "limit_key": "max_items",
            "limit_default": 20,
            "unique_attribute": "id",
        },
    }

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.__route53_backend: Route53Backend = route53_backends[account_id][
            self.partition
        ]
        self.__domains: Dict[str, Route53Domain] = {}
        self.__operations: Dict[str, Route53DomainsOperation] = {}

    def register_domain(
        self,
        domain_name: str,
        duration_in_years: int,
        auto_renew: bool,
        admin_contact: Dict[str, Any],
        registrant_contact: Dict[str, Any],
        tech_contact: Dict[str, Any],
        private_protect_admin_contact: bool,
        private_protect_registrant_contact: bool,
        private_protect_tech_contact: bool,
        extra_params: List[Dict[str, Any]],
    ) -> Route53DomainsOperation:
        """Register a domain"""

        if len(self.__domains) == self.DEFAULT_MAX_DOMAINS_COUNT:
            raise DomainLimitExceededException()

        requested_operation = Route53DomainsOperation.validate(
            domain_name=domain_name, status="SUCCESSFUL", type_="REGISTER_DOMAIN"
        )

        self.__validate_duplicate_operations(requested_operation)

        expiration_date = datetime.now(timezone.utc) + timedelta(
            days=365 * duration_in_years
        )

        try:
            domain = Route53Domain.validate(
                domain_name=domain_name,
                auto_renew=auto_renew,
                admin_contact=Route53DomainsContactDetail.validate_dict(admin_contact),
                registrant_contact=Route53DomainsContactDetail.validate_dict(
                    registrant_contact
                ),
                tech_contact=Route53DomainsContactDetail.validate_dict(tech_contact),
                admin_privacy=private_protect_admin_contact,
                registrant_privacy=private_protect_registrant_contact,
                tech_privacy=private_protect_tech_contact,
                expiration_date=expiration_date,
                extra_params=extra_params,
            )

        except ValidationException as e:
            raise InvalidInputException(e.errors)
        self.__operations[requested_operation.id] = requested_operation

        self.__route53_backend.create_hosted_zone(
            name=domain.domain_name, private_zone=False
        )

        self.__domains[domain_name] = domain
        return requested_operation

    def delete_domain(self, domain_name: str) -> Route53DomainsOperation:
        requested_operation = Route53DomainsOperation.validate(
            domain_name=domain_name, status="SUCCESSFUL", type_="DELETE_DOMAIN"
        )
        self.__validate_duplicate_operations(requested_operation)

        input_errors: List[str] = []
        Route53Domain.validate_domain_name(domain_name, input_errors)

        if input_errors:
            raise InvalidInputException(input_errors)

        if domain_name not in self.__domains:
            raise InvalidInputException(
                [f"Domain {domain_name} isn't registered in the current account"]
            )

        self.__operations[requested_operation.id] = requested_operation
        del self.__domains[domain_name]
        return requested_operation

    def __validate_duplicate_operations(
        self, requested_operation: Route53DomainsOperation
    ) -> None:
        for operation in self.__operations.values():
            if (
                operation.domain_name == requested_operation.domain_name
                and operation.type == requested_operation.type
            ):
                raise DuplicateRequestException()

    def get_domain(self, domain_name: str) -> Route53Domain:
        input_errors: List[str] = []
        Route53Domain.validate_domain_name(domain_name, input_errors)
        if input_errors:
            raise InvalidInputException(input_errors)

        if domain_name not in self.__domains:
            raise InvalidInputException(["Domain is not associated with this account"])

        return self.__domains[domain_name]

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_domains(
        self,
        filter_conditions: Optional[List[Dict[str, Any]]] = None,
        sort_condition: Optional[Dict[str, Any]] = None,
    ) -> List[Route53Domain]:
        try:
            filters: List[DomainsFilter] = (
                [DomainsFilter.validate_dict(f) for f in filter_conditions]
                if filter_conditions
                else []
            )
            sort: Optional[DomainsSortCondition] = (
                DomainsSortCondition.validate_dict(sort_condition)
                if sort_condition
                else None
            )
        except ValidationException as e:
            raise InvalidInputException(e.errors)

        filter_fields = [f.name for f in filters]
        if sort and filter_fields and sort.name not in filter_fields:
            raise InvalidInputException(
                ["Sort condition must be the same as the filter condition"]
            )

        domains_to_return: List[Route53Domain] = []

        for domain in self.__domains.values():
            if all([f.filter(domain) for f in filters]):
                domains_to_return.append(domain)

        if sort:
            if sort.name == DomainFilterField.DOMAIN_NAME:
                domains_to_return.sort(
                    key=lambda d: d.domain_name,
                    reverse=(sort.sort_order == DomainSortOrder.DESCENDING),
                )
            else:
                domains_to_return.sort(
                    key=lambda d: d.expiration_date,
                    reverse=(sort.sort_order == DomainSortOrder.DESCENDING),
                )

        return domains_to_return

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_operations(
        self,
        submitted_since_timestamp: Optional[int] = None,
        statuses: Optional[List[str]] = None,
        types: Optional[List[str]] = None,
        sort_by: Optional[str] = None,
        sort_order: Optional[str] = None,
    ) -> List[Route53DomainsOperation]:
        input_errors: List[str] = []
        statuses = statuses or []
        types = types or []

        if any(status not in DOMAIN_OPERATION_STATUSES for status in statuses):
            input_errors.append("Status is invalid")
        if any(type_ not in DOMAIN_OPERATION_TYPES for type_ in types):
            input_errors.append("Type is invalid")

        if input_errors:
            raise InvalidInputException(input_errors)

        submitted_since = (
            datetime.fromtimestamp(submitted_since_timestamp, timezone.utc)
            if submitted_since_timestamp
            else None
        )

        operations_to_return: List[Route53DomainsOperation] = []

        for operation in self.__operations.values():
            if statuses and operation.status not in statuses:
                continue

            if types and operation.type not in types:
                continue

            if submitted_since and operation.submitted_date < submitted_since:
                continue

            operations_to_return.append(operation)

        if sort_by == "SubmittedDate":
            operations_to_return.sort(
                key=lambda op: op.submitted_date,
                reverse=sort_order == DomainSortOrder.ASCENDING,
            )

        return operations_to_return

    def get_operation(self, operation_id: str) -> Route53DomainsOperation:
        if operation_id not in self.__operations:
            raise InvalidInputException(
                [f"Operation with id {operation_id} doesn't exist"]
            )

        return self.__operations[operation_id]

    def update_domain_nameservers(
        self, domain_name: str, nameservers: List[Dict[str, Any]]
    ) -> Route53DomainsOperation:
        input_errors: List[str] = []
        Route53Domain.validate_domain_name(domain_name, input_errors)
        if len(nameservers) < 1:
            input_errors.append("Must supply nameservers")

        servers: List[NameServer] = []
        try:
            servers = [NameServer.validate_dict(obj) for obj in nameservers]
        except ValidationException as e:
            input_errors += e.errors

        for server in servers:
            if domain_name in server.name and not server.glue_ips:
                input_errors.append(
                    f"Must supply glue IPs for name server {server.name} because it is a subdomain of "
                    f"the domain"
                )

        if input_errors:
            raise InvalidInputException(input_errors)

        if domain_name not in self.__domains:
            raise InvalidInputException(
                [f"Domain {domain_name} is not registered to the current AWS account"]
            )

        requested_operation = Route53DomainsOperation.validate(
            domain_name=domain_name, status="SUCCESSFUL", type_="UPDATE_NAMESERVER"
        )
        self.__validate_duplicate_operations(requested_operation)

        domain = self.__domains[domain_name]
        domain.nameservers = servers
        self.__operations[requested_operation.id] = requested_operation

        return requested_operation


route53domains_backends = BackendDict(
    Route53DomainsBackend,
    "route53domains",
    use_boto3_regions=False,
    additional_regions=PARTITION_NAMES,
)
