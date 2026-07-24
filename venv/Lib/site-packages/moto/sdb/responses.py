from moto.core.responses import ActionResult, BaseResponse, EmptyResult

from .models import SimpleDBBackend, sdb_backends


class SimpleDBResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="sdb")
        self.automated_parameter_parsing = True

    @property
    def sdb_backend(self) -> SimpleDBBackend:
        return sdb_backends[self.current_account][self.region]

    def create_domain(self) -> ActionResult:
        domain_name = self._get_param("DomainName")
        self.sdb_backend.create_domain(domain_name=domain_name)
        return EmptyResult()

    def delete_domain(self) -> ActionResult:
        domain_name = self._get_param("DomainName")
        self.sdb_backend.delete_domain(domain_name=domain_name)
        return EmptyResult()

    def list_domains(self) -> ActionResult:
        domain_names = self.sdb_backend.list_domains()
        result = {"DomainNames": domain_names}
        return ActionResult(result)

    def get_attributes(self) -> ActionResult:
        domain_name = self._get_param("DomainName")
        item_name = self._get_param("ItemName")
        attribute_names = self._get_param("AttributeNames")
        attributes = self.sdb_backend.get_attributes(
            domain_name=domain_name,
            item_name=item_name,
            attribute_names=attribute_names,
        )
        result = {"Attributes": attributes}
        return ActionResult(result)

    def put_attributes(self) -> ActionResult:
        domain_name = self._get_param("DomainName")
        item_name = self._get_param("ItemName")
        attributes = self._get_param("Attributes")
        self.sdb_backend.put_attributes(
            domain_name=domain_name, item_name=item_name, attributes=attributes
        )
        return EmptyResult()
