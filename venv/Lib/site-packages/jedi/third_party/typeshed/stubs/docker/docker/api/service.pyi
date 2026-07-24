class ServiceApiMixin:
    def create_service(
        self,
        task_template,
        name=None,
        labels=None,
        mode=None,
        update_config=None,
        networks=None,
        endpoint_config=None,
        endpoint_spec=None,
        rollback_config=None,
    ): ...
    def inspect_service(self, service, insert_defaults=None): ...
    def inspect_task(self, task): ...
    def remove_service(self, service): ...
    def services(self, filters=None, status=None): ...
    def service_logs(
        self,
        service,
        details: bool = False,
        follow: bool = False,
        stdout: bool = False,
        stderr: bool = False,
        since: int = 0,
        timestamps: bool = False,
        tail: str = "all",
        is_tty=None,
    ): ...
    def tasks(self, filters=None): ...
    def update_service(
        self,
        service,
        version,
        task_template=None,
        name=None,
        labels=None,
        mode=None,
        update_config=None,
        networks=None,
        endpoint_config=None,
        endpoint_spec=None,
        fetch_current_spec: bool = False,
        rollback_config=None,
    ): ...
