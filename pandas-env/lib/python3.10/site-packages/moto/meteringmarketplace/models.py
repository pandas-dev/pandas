import collections
from typing import Any, Deque, Dict, List

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.moto_api._internal import mock_random


class UsageRecord(BaseModel, Dict[str, Any]):  # type: ignore[misc]
    def __init__(
        self,
        timestamp: str,
        customer_identifier: str,
        dimension: str,
        quantity: int = 0,
    ):
        super().__init__()
        self.timestamp = timestamp
        self.customer_identifier = customer_identifier
        self.dimension = dimension
        self.quantity = quantity
        self.metering_record_id = mock_random.uuid4().hex

    @property
    def timestamp(self) -> str:
        return self["Timestamp"]

    @timestamp.setter
    def timestamp(self, value: str) -> None:
        self["Timestamp"] = value

    @property
    def customer_identifier(self) -> str:
        return self["CustomerIdentifier"]

    @customer_identifier.setter
    def customer_identifier(self, value: str) -> None:
        self["CustomerIdentifier"] = value

    @property
    def dimension(self) -> str:
        return self["Dimension"]

    @dimension.setter
    def dimension(self, value: str) -> None:
        self["Dimension"] = value

    @property
    def quantity(self) -> int:
        return self["Quantity"]

    @quantity.setter
    def quantity(self, value: int) -> None:
        self["Quantity"] = value


class Result(BaseModel, Dict[str, Any]):  # type: ignore[misc]
    SUCCESS = "Success"
    CUSTOMER_NOT_SUBSCRIBED = "CustomerNotSubscribed"
    DUPLICATE_RECORD = "DuplicateRecord"

    def __init__(self, **kwargs: Any):
        self.usage_record = UsageRecord(
            timestamp=kwargs["Timestamp"],
            customer_identifier=kwargs["CustomerIdentifier"],
            dimension=kwargs["Dimension"],
            quantity=kwargs["Quantity"],
        )
        self.status = Result.SUCCESS
        self["MeteringRecordId"] = self.usage_record.metering_record_id

    @property
    def metering_record_id(self) -> str:
        return self["MeteringRecordId"]

    @property
    def status(self) -> str:
        return self["Status"]

    @status.setter
    def status(self, value: str) -> None:
        self["Status"] = value

    @property
    def usage_record(self) -> UsageRecord:
        return self["UsageRecord"]

    @usage_record.setter
    def usage_record(self, value: UsageRecord) -> None:
        self["UsageRecord"] = value

    def is_duplicate(self, other: Any) -> bool:
        """
        DuplicateRecord - Indicates that the UsageRecord was invalid and not honored.
        A previously metered UsageRecord had the same customer, dimension, and time,
        but a different quantity.
        """
        assert isinstance(other, Result), "Needs to be a Result type"
        usage_record, other = other.usage_record, self.usage_record
        return (
            other.customer_identifier == usage_record.customer_identifier
            and other.dimension == usage_record.dimension
            and other.timestamp == usage_record.timestamp
            and other.quantity != usage_record.quantity
        )


class CustomerDeque(Deque[str]):
    def is_subscribed(self, customer: str) -> bool:
        return customer in self


class ResultDeque(Deque[Result]):
    def is_duplicate(self, result: Result) -> bool:
        return any(record.is_duplicate(result) for record in self)


class MeteringMarketplaceBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.customers_by_product: Dict[str, CustomerDeque] = collections.defaultdict(
            CustomerDeque
        )
        self.records_by_product: Dict[str, ResultDeque] = collections.defaultdict(
            ResultDeque
        )

    def batch_meter_usage(
        self, product_code: str, usage_records: List[Dict[str, Any]]
    ) -> List[Result]:
        results = []
        for usage in usage_records:
            result = Result(**usage)
            if not self.customers_by_product[product_code].is_subscribed(
                result.usage_record.customer_identifier
            ):
                result.status = result.CUSTOMER_NOT_SUBSCRIBED
            elif self.records_by_product[product_code].is_duplicate(result):
                result.status = result.DUPLICATE_RECORD
            results.append(result)
        return results


meteringmarketplace_backends = BackendDict(
    MeteringMarketplaceBackend, "meteringmarketplace"
)
