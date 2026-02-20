import re
import string
from collections.abc import Callable, Iterator
from functools import lru_cache
from threading import RLock
from typing import (
    TYPE_CHECKING,
    Any,
    ClassVar,
    Optional,
    TypeVar,
)
from uuid import uuid4

from boto3 import Session

from moto.settings import allow_unknown_region, enable_iso_regions
from moto.utilities.utils import get_partition

from .model_instances import model_data
from .responses import TYPE_RESPONSE
from .utils import ISO_REGION_DOMAINS, convert_regex_to_flask_path

if TYPE_CHECKING:
    from moto.core.common_models import BaseModel


class InstanceTrackerMeta(type):
    def __new__(meta, name: str, bases: Any, dct: dict[str, Any]) -> type:
        cls = super().__new__(meta, name, bases, dct)
        if name == "BaseModel":
            return cls

        service = cls.__module__.split(".")[1]
        if name not in model_data[service]:
            model_data[service][name] = cls
        cls.instances_tracked: ClassVar[list[BaseModel]] = []  # type: ignore
        return cls


class BaseBackend:
    def __init__(self, region_name: str, account_id: str):
        self.region_name = region_name
        self.account_id = account_id
        self.partition = get_partition(region_name)

    def reset(self) -> None:
        region_name = self.region_name
        account_id = self.account_id
        self.__dict__ = {}
        self.__init__(region_name, account_id)  # type: ignore[misc]

    @property
    def _url_module(self) -> Any:  # type: ignore[misc]
        backend_module = self.__class__.__module__
        # moto.service.models --> moto.service.urls
        backend_urls_module_name = backend_module.rstrip("models") + "urls"
        backend_urls_module = __import__(
            backend_urls_module_name, fromlist=["url_bases", "url_paths"]
        )
        return backend_urls_module

    @property
    def urls(self) -> dict[str, Callable[[Any, str, Any], TYPE_RESPONSE]]:  # type: ignore[misc]
        """
        A dictionary of the urls to be mocked with this service and the handlers
        that should be called in their place
        """
        url_bases = self.url_bases
        unformatted_paths = self._url_module.url_paths

        urls = {}
        for url_base in url_bases:
            # The default URL_base will look like: http://service.[..].amazonaws.com/...
            # This extension ensures support for the China & ISO regions
            alt_dns_suffixes = {"cn": "amazonaws.com.cn"}
            if enable_iso_regions():
                alt_dns_suffixes.update(ISO_REGION_DOMAINS)

            for url_path, handler in unformatted_paths.items():
                url = url_path.format(url_base)
                urls[url] = handler
                for dns_suffix in alt_dns_suffixes.values():
                    alt_url_base = re.sub(r"amazonaws\\?.com$", dns_suffix, url_base)
                    alt_url = url_path.format(alt_url_base)
                    urls[alt_url] = handler

        return urls

    @property
    def url_paths(self) -> dict[str, Callable[[Any, str, Any], TYPE_RESPONSE]]:  # type: ignore[misc]
        """
        A dictionary of the paths of the urls to be mocked with this service and
        the handlers that should be called in their place
        """
        unformatted_paths = self._url_module.url_paths

        paths = {}
        for unformatted_path, handler in unformatted_paths.items():
            path = unformatted_path.format("")
            paths[path] = handler

        return paths

    @property
    def url_bases(self) -> list[str]:
        """
        A list containing the url_bases extracted from urls.py
        """
        return self._url_module.url_bases

    @property
    def flask_paths(self) -> dict[str, Callable[[Any, str, Any], TYPE_RESPONSE]]:  # type: ignore[misc]
        """
        The url paths that will be used for the flask server
        """
        paths = {}
        for url_path, handler in self.url_paths.items():
            url_path = convert_regex_to_flask_path(url_path)
            paths[url_path] = handler

        return paths

    @staticmethod
    def default_vpc_endpoint_service(
        service_region: str,
        zones: list[str],
    ) -> list[dict[str, str]]:
        """Invoke the factory method for any VPC endpoint(s) services."""
        return []

    @staticmethod
    def vpce_random_number() -> str:
        from moto.moto_api._internal import mock_random as random

        """Return random number for a VPC endpoint service ID."""
        return "".join([random.choice(string.hexdigits.lower()) for i in range(17)])

    @staticmethod
    def default_vpc_endpoint_service_factory(  # type: ignore[misc]
        service_region: str,
        zones: list[str],
        service: str = "",
        service_type: str = "Interface",
        private_dns_names: bool = True,
        special_service_name: str = "",
        policy_supported: bool = True,
        base_endpoint_dns_names: Optional[list[str]] = None,
    ) -> list[dict[str, Any]]:
        """List of dicts representing default VPC endpoints for this service."""
        if special_service_name:
            service_name = f"com.amazonaws.{service_region}.{special_service_name}"
        else:
            service_name = f"com.amazonaws.{service_region}.{service}"

        if not base_endpoint_dns_names:
            base_endpoint_dns_names = [f"{service}.{service_region}.vpce.amazonaws.com"]

        endpoint_service = {
            "AcceptanceRequired": False,
            "AvailabilityZones": zones,
            "BaseEndpointDnsNames": base_endpoint_dns_names,
            "ManagesVpcEndpoints": False,
            "Owner": "amazon",
            "ServiceId": f"vpce-svc-{BaseBackend.vpce_random_number()}",
            "ServiceName": service_name,
            "ServiceType": [{"ServiceType": service_type}],
            "Tags": [],
            "VpcEndpointPolicySupported": policy_supported,
        }

        # Don't know how private DNS names are different, so for now just
        # one will be added.
        if private_dns_names:
            endpoint_service["PrivateDnsName"] = (
                f"{service}.{service_region}.amazonaws.com"
            )
            endpoint_service["PrivateDnsNameVerificationState"] = "verified"
            endpoint_service["PrivateDnsNames"] = [
                {"PrivateDnsName": f"{service}.{service_region}.amazonaws.com"}
            ]
        return [endpoint_service]

    # def list_config_service_resources(self, resource_ids, resource_name, limit, next_token):
    #     """For AWS Config. This will list all of the resources of the given type and optional resource name and region"""
    #     raise NotImplementedError()


backend_lock = RLock()
SERVICE_BACKEND = TypeVar("SERVICE_BACKEND", bound=BaseBackend)


class AccountSpecificBackend(dict[str, SERVICE_BACKEND]):
    """
    Dictionary storing the data for a service in a specific account.
    Data access pattern:
      account_specific_backend[region: str] = backend: BaseBackend
    """

    def __init__(
        self,
        service_name: str,
        account_id: str,
        backend: type,
        use_boto3_regions: bool,
        additional_regions: Optional[list[str]],
    ):
        self._id = str(uuid4())
        self.service_name = service_name
        self.account_id = account_id
        self.backend = backend
        self.regions = []
        if use_boto3_regions:
            self.regions.extend(self._generate_regions(service_name))
        self.regions.extend(additional_regions or [])

    @lru_cache
    def _generate_regions(self, service_name: str) -> list[str]:
        regions = []
        for (
            partition
        ) in AccountSpecificBackend.get_session().get_available_partitions():
            partition_regions = (
                AccountSpecificBackend.get_session().get_available_regions(
                    service_name, partition_name=partition
                )
            )
            regions.extend(partition_regions)
        return regions

    @classmethod
    @lru_cache
    def get_session(cls) -> Session:  # type: ignore[misc]
        # Only instantiate Session when we absolutely need it
        # This gives the user time to remove any env variables that break botocore, like AWS_PROFILE
        # See https://github.com/getmoto/moto/issues/5469#issuecomment-2474897120
        return Session()

    def __hash__(self) -> int:  # type: ignore[override]
        return hash(self._id)

    def __eq__(self, other: Any) -> bool:
        return (
            other
            and isinstance(other, AccountSpecificBackend)
            and other._id == self._id
        )

    def __ne__(self, other: Any) -> bool:
        return not self.__eq__(other)

    def reset(self) -> None:
        for region_specific_backend in self.values():
            region_specific_backend.reset()

    def __contains__(self, region: str) -> bool:  # type: ignore[override]
        if region == "global":
            region = "aws"
        return region in self.regions or region in self.keys()

    def __delitem__(self, key: str) -> None:
        super().__delitem__(key)

    def __iter__(self) -> Iterator[str]:
        return super().__iter__()

    def __len__(self) -> int:
        return super().__len__()

    def __setitem__(self, key: str, value: SERVICE_BACKEND) -> None:
        super().__setitem__(key, value)

    def __getitem__(self, region_name: str) -> SERVICE_BACKEND:
        # Some services, like S3, used to be truly global - meaning one Backend serving all
        # Now that we support partitions (AWS, AWS-CN, AWS-GOV, etc), there will be one backend per partition
        # Because the concept of 'region' doesn't exist in a global service, we use the partition name to keep the backends separate
        # We used to use the term 'global' in lieu of a region name, and users may still use this
        # It should resolve to 'aws', to ensure consistency
        if region_name == "global":
            region_name = "aws"
        if region_name in self.keys():
            return super().__getitem__(region_name)
        # Create the backend for a specific region
        with backend_lock:
            if region_name in self.regions and region_name not in self.keys():
                super().__setitem__(
                    region_name, self.backend(region_name, account_id=self.account_id)
                )
            if region_name not in self.regions and allow_unknown_region():
                super().__setitem__(
                    region_name, self.backend(region_name, account_id=self.account_id)
                )
        return super().__getitem__(region_name)


class BackendDict(dict[str, AccountSpecificBackend[SERVICE_BACKEND]]):
    """
    Data Structure to store everything related to a specific service.
    Format:
      [account_id: str]: AccountSpecificBackend
      [account_id: str][region: str] = BaseBackend
    """

    # We keep track of the BackendDict's that were:
    # - instantiated
    # - contain at least one AccountSpecificBackend
    #
    # In other words, this is the list of backends which are in use by the user
    #   making it trivial to determine which backends to reset when the mocks end
    _instances: list["BackendDict[SERVICE_BACKEND]"] = []

    @classmethod
    def reset(cls) -> None:
        with backend_lock:
            for backend in BackendDict._instances:  # type: ignore[misc]
                for account_specific_backend in backend.values():
                    account_specific_backend.reset()
                    # account_specific_backend.__getitem__.cache_clear()
                backend.clear()
            # https://github.com/getmoto/moto/issues/6592
            # Could be fixed by removing the cache, forcing all data to be regenerated every reset
            # But this also incurs a significant performance hit
            # backend.__getitem__.cache_clear()
            BackendDict._instances.clear()  # type: ignore[misc]

    def __init__(
        self,
        backend: type[SERVICE_BACKEND],
        service_name: str,
        use_boto3_regions: bool = True,
        additional_regions: Optional[list[str]] = None,
    ):
        self.backend = backend
        self.service_name = service_name
        self._use_boto3_regions = use_boto3_regions
        self._additional_regions = additional_regions
        self._id = str(uuid4())

    def __hash__(self) -> int:  # type: ignore[override]
        # Required for the LRUcache to work.
        # service_name is enough to determine uniqueness - other properties are dependent
        return hash(self._id)

    def __eq__(self, other: Any) -> bool:
        return other and isinstance(other, BackendDict) and other._id == self._id

    def __ne__(self, other: Any) -> bool:
        return not self.__eq__(other)

    def __getitem__(self, account_id: str) -> AccountSpecificBackend[SERVICE_BACKEND]:
        self._create_account_specific_backend(account_id)
        return super().__getitem__(account_id)

    def __delitem__(self, key: str) -> None:
        super().__delitem__(key)

    def __iter__(self) -> Iterator[str]:
        return super().__iter__()

    def __len__(self) -> int:
        return super().__len__()

    def __setitem__(
        self,
        key: str,
        value: AccountSpecificBackend[SERVICE_BACKEND],
    ) -> None:
        super().__setitem__(key, value)

    def _create_account_specific_backend(self, account_id: str) -> None:
        with backend_lock:
            if account_id not in list(self.keys()):
                self[account_id] = AccountSpecificBackend(
                    service_name=self.service_name,
                    account_id=account_id,
                    backend=self.backend,
                    use_boto3_regions=self._use_boto3_regions,
                    additional_regions=self._additional_regions,
                )
                BackendDict._instances.append(self)  # type: ignore[misc]

    def iter_backends(self) -> Iterator[tuple[str, str, BaseBackend]]:
        """
        Iterate over a flattened view of all base backends in a BackendDict.
        Each record is a tuple of account id, region name, and the base backend within that account and region.
        """
        for account_id, account_specific_backend in self.items():
            for region_name, backend in account_specific_backend.items():
                yield account_id, region_name, backend
