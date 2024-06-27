from moto.core.responses import BaseResponse

from .models import SimpleDBBackend, sdb_backends


class SimpleDBResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="sdb")

    @property
    def sdb_backend(self) -> SimpleDBBackend:
        return sdb_backends[self.current_account][self.region]

    def create_domain(self) -> str:
        domain_name = self._get_param("DomainName")
        self.sdb_backend.create_domain(domain_name=domain_name)
        template = self.response_template(CREATE_DOMAIN_TEMPLATE)
        return template.render()

    def delete_domain(self) -> str:
        domain_name = self._get_param("DomainName")
        self.sdb_backend.delete_domain(domain_name=domain_name)
        template = self.response_template(DELETE_DOMAIN_TEMPLATE)
        return template.render()

    def list_domains(self) -> str:
        domain_names = self.sdb_backend.list_domains()
        template = self.response_template(LIST_DOMAINS_TEMPLATE)
        return template.render(domain_names=domain_names, next_token=None)

    def get_attributes(self) -> str:
        domain_name = self._get_param("DomainName")
        item_name = self._get_param("ItemName")
        attribute_names = self._get_multi_param("AttributeName.")
        attributes = self.sdb_backend.get_attributes(
            domain_name=domain_name,
            item_name=item_name,
            attribute_names=attribute_names,
        )
        template = self.response_template(GET_ATTRIBUTES_TEMPLATE)
        return template.render(attributes=attributes)

    def put_attributes(self) -> str:
        domain_name = self._get_param("DomainName")
        item_name = self._get_param("ItemName")
        attributes = self._get_list_prefix("Attribute")
        self.sdb_backend.put_attributes(
            domain_name=domain_name, item_name=item_name, attributes=attributes
        )
        template = self.response_template(PUT_ATTRIBUTES_TEMPLATE)
        return template.render()


CREATE_DOMAIN_TEMPLATE = """<?xml version="1.0" encoding="UTF-8"?>
<CreateDomainResult  xmlns="http://sdb.amazonaws.com/doc/2009-04-15/"></CreateDomainResult>
"""


LIST_DOMAINS_TEMPLATE = """<?xml version="1.0" encoding="UTF-8"?>
<ListDomainsResponse  xmlns="http://sdb.amazonaws.com/doc/2009-04-15/">
    <ListDomainsResult>
        {% for name in domain_names %}
            <DomainName>{{ name }}</DomainName>
        {% endfor %}
        <NextToken>{{ next_token }}</NextToken>
    </ListDomainsResult>
</ListDomainsResponse>
"""

DELETE_DOMAIN_TEMPLATE = """<?xml version="1.0"?>
<DeleteDomainResponse xmlns="http://sdb.amazonaws.com/doc/2009-04-15/">
  <ResponseMetadata>
    <RequestId>64d9c3ac-ef19-2e3d-7a03-9ea46205eb71</RequestId>
    <BoxUsage>0.0055590278</BoxUsage>
  </ResponseMetadata>
</DeleteDomainResponse>"""

PUT_ATTRIBUTES_TEMPLATE = """<?xml version="1.0" encoding="UTF-8"?>
<PutAttributesResult xmlns="http://sdb.amazonaws.com/doc/2009-04-15/"></PutAttributesResult>
"""

GET_ATTRIBUTES_TEMPLATE = """<GetAttributesResponse xmlns="http://sdb.amazonaws.com/doc/2009-04-15/">
  <ResponseMetadata>
    <RequestId>1549581b-12b7-11e3-895e-1334aEXAMPLE</RequestId>
  </ResponseMetadata>
  <GetAttributesResult>
{% for attribute in attributes %}
      <Attribute>
        <Name>{{ attribute["name"] }}</Name>
        <Value>{{ attribute["value"] }}</Value>
      </Attribute>
{% endfor %}
  </GetAttributesResult>
</GetAttributesResponse>"""
