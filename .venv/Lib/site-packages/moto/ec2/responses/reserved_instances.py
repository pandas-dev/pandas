from moto.core.responses import BaseResponse


class ReservedInstances(BaseResponse):
    def cancel_reserved_instances_listing(self) -> None:
        self.error_on_dryrun()

        raise NotImplementedError(
            "ReservedInstances.cancel_reserved_instances_listing is not yet implemented"
        )

    def create_reserved_instances_listing(self) -> None:
        self.error_on_dryrun()

        raise NotImplementedError(
            "ReservedInstances.create_reserved_instances_listing is not yet implemented"
        )

    def describe_reserved_instances(self) -> None:
        raise NotImplementedError(
            "ReservedInstances.describe_reserved_instances is not yet implemented"
        )

    def describe_reserved_instances_listings(self) -> None:
        raise NotImplementedError(
            "ReservedInstances.describe_reserved_instances_listings is not yet implemented"
        )

    def describe_reserved_instances_offerings(self) -> None:
        raise NotImplementedError(
            "ReservedInstances.describe_reserved_instances_offerings is not yet implemented"
        )

    def purchase_reserved_instances_offering(self) -> None:
        self.error_on_dryrun()

        raise NotImplementedError(
            "ReservedInstances.purchase_reserved_instances_offering is not yet implemented"
        )
