from jinja2 import Template

from moto.core.exceptions import RESTError


class NeptuneClientError(RESTError):
    def __init__(self, code: str, message: str):
        super().__init__(error_type=code, message=message)
        template = Template(
            """
        <ErrorResponse>
            <Error>
              <Code>{{ code }}</Code>
              <Message>{{ message }}</Message>
              <Type>Sender</Type>
            </Error>
            <RequestId>6876f774-7273-11e4-85dc-39e55ca848d1</RequestId>
        </ErrorResponse>"""
        )
        self.description = template.render(code=code, message=message)


class DBClusterNotFoundError(NeptuneClientError):
    def __init__(self, cluster_identifier: str):
        super().__init__(
            "DBClusterNotFoundFault", f"DBCluster {cluster_identifier} not found."
        )
