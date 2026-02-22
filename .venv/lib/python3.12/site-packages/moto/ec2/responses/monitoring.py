from ._base_response import EC2BaseResponse


class Monitoring(EC2BaseResponse):
    def monitor_instances(self) -> str:
        self.error_on_dryrun()

        raise NotImplementedError("Monitoring.monitor_instances is not yet implemented")

    def unmonitor_instances(self) -> str:
        self.error_on_dryrun()

        raise NotImplementedError(
            "Monitoring.unmonitor_instances is not yet implemented"
        )
