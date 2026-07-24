from typing import Any

from .base import DictType

class Healthcheck(DictType[Any]):
    def __init__(
        self,
        *,
        test: str | list[str] | None = ...,
        Test: str | list[str] | None = ...,
        interval: int | None = ...,
        Interval: int | None = ...,
        timeout: int | None = ...,
        Timeout: int | None = ...,
        retries: int | None = ...,
        Retries: int | None = ...,
        start_period: int | None = ...,
        StartPeriod: int | None = ...,
    ) -> None: ...
    @property
    def test(self) -> list[str] | None: ...
    @test.setter
    def test(self, value: str | list[str] | None) -> None: ...
    @property
    def interval(self) -> int | None: ...
    @interval.setter
    def interval(self, value: int | None) -> None: ...
    @property
    def timeout(self) -> int | None: ...
    @timeout.setter
    def timeout(self, value: int | None) -> None: ...
    @property
    def retries(self) -> int | None: ...
    @retries.setter
    def retries(self, value: int | None) -> None: ...
    @property
    def start_period(self) -> int | None: ...
    @start_period.setter
    def start_period(self, value: int | None) -> None: ...
