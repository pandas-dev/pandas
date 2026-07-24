from __future__ import annotations

import itertools
import random
from typing import Any, cast


class ReservedInstancesOffering(dict[str, Any]):
    """Lightweight container for a Reserved Instance Offering.

    This is a simplified representation that captures commonly-used fields
    for testing without modeling full AWS behavior.
    """

    _filter_attributes = {
        "availability-zone": ["AvailabilityZone"],
        "instance-type": ["InstanceType"],
        "product-description": ["ProductDescription"],
        "offering-type": ["OfferingType"],
    }

    def __init__(self, **kwargs: Any):
        super().__init__(**kwargs)

    def get_filter_value(self, filter_name: str) -> Any:
        path = self._filter_attributes.get(filter_name)
        value: Any = self
        if not path:
            # Unknown filter -> don't match
            return None
        for key in path:
            value = (value or {}).get(key)
        return value


class ReservedInstancesBackend:
    """Backend mixin to support DescribeReservedInstancesOfferings.

    We generate a small, static set of offerings per region/zone to make
    local tests useful while keeping scope limited.
    """

    _ri_offerings_generated: bool = False
    _ri_offerings: dict[str, list[ReservedInstancesOffering]]

    def _ensure_reserved_offerings(self) -> None:
        if getattr(self, "_ri_offerings_generated", False):
            return

        # Build a handful of offerings for the current region based on its AZs
        region = self.region_name  # type: ignore[attr-defined]
        zones = [z.name for z in self.describe_availability_zones()]  # type: ignore[attr-defined]
        instance_types = [
            "t2.micro",
            "t3.small",
            "m5.large",
        ]
        offering_types = ["No Upfront", "Partial Upfront", "All Upfront"]
        product = "Linux/UNIX"

        offerings: list[ReservedInstancesOffering] = []
        rid_counter = 1
        for itype, zone, otype in itertools.product(
            instance_types, zones, offering_types
        ):
            rid = f"{region}-ri-off-{rid_counter:04d}"
            rid_counter += 1
            duration = 31536000  # 1 year in seconds
            fixed = 0.0 if otype == "No Upfront" else (float(random.randint(50, 300)))
            usage = 0.0 if otype == "All Upfront" else 0.05

            offerings.append(
                ReservedInstancesOffering(
                    ReservedInstancesOfferingId=rid,
                    InstanceType=itype,
                    AvailabilityZone=zone,
                    ProductDescription=product,
                    Duration=duration,
                    UsagePrice=usage,
                    FixedPrice=fixed,
                    CurrencyCode="USD",
                    InstanceTenancy="default",
                    OfferingClass="standard",
                    OfferingType=otype,
                    Marketplace=False,
                )
            )

        # Also add region-wide offerings (no specific AZ) for convenience
        for itype, otype in itertools.product(instance_types, offering_types):
            rid = f"{region}-ri-off-{rid_counter:04d}"
            rid_counter += 1
            offerings.append(
                ReservedInstancesOffering(
                    ReservedInstancesOfferingId=rid,
                    InstanceType=itype,
                    AvailabilityZone=None,
                    ProductDescription=product,
                    Duration=31536000,
                    UsagePrice=0.05,
                    FixedPrice=0.0 if otype == "No Upfront" else 100.0,
                    CurrencyCode="USD",
                    InstanceTenancy="default",
                    OfferingClass="standard",
                    OfferingType=otype,
                    Marketplace=False,
                )
            )

        if not hasattr(self, "_ri_offerings"):
            self._ri_offerings = {}
        self._ri_offerings[region] = offerings
        self._ri_offerings_generated = True

    def describe_reserved_instances_offerings(
        self, filters: dict[str, list[str]] | None = None
    ) -> list[dict[str, Any]]:
        self._ensure_reserved_offerings()
        region = self.region_name  # type: ignore[attr-defined]
        offerings = self._ri_offerings.get(region, [])

        if not filters:
            return cast(list[dict[str, Any]], offerings)

        # Simple filter predicate across supported keys
        def predicate(obj: ReservedInstancesOffering, ff: dict[str, list[str]]) -> bool:
            for name, values in ff.items():
                val = obj.get_filter_value(name)
                if val is None:
                    # Unsupported filter -> treat as non-match
                    return False
                if isinstance(val, list):
                    if not any(v in val for v in values):
                        return False
                else:
                    if val not in values:
                        return False
            return True

        return cast(
            list[dict[str, Any]], [o for o in offerings if predicate(o, filters)]
        )
