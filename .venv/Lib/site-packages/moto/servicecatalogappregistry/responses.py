"""Handles incoming servicecatalogappregistry requests, invokes methods, returns responses."""

import json
from urllib.parse import unquote

from moto.core.responses import BaseResponse
from moto.servicecatalogappregistry.exceptions import ValidationException

from .models import (
    Application,
    AppRegistryBackend,
    servicecatalogappregistry_backends,
)


class AppRegistryResponse(BaseResponse):
    """Handler for AppRegistry requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="servicecatalog-appregistry")

    @property
    def servicecatalogappregistry_backend(self) -> AppRegistryBackend:
        """Return backend instance specific for this region."""
        return servicecatalogappregistry_backends[self.current_account][self.region]

    def create_application(self) -> str:
        name = self._get_param("name")
        description = self._get_param("description")
        tags = self._get_param("tags")
        client_token = self._get_param("clientToken")
        application = self.servicecatalogappregistry_backend.create_application(
            name=name,
            description=description,
            tags=tags,
            client_token=client_token,
        )
        return json.dumps({"application": application.to_json()})

    def list_applications(self) -> str:
        applications = self.servicecatalogappregistry_backend.list_applications()
        json_list = []
        for app in applications:
            json_list.append(app.to_json())
        return json.dumps({"applications": json_list})

    def associate_resource(self) -> str:
        application = unquote(self._get_param("application"))
        resource_type = self._get_param("resourceType")
        resource = unquote(self._get_param("resource"))
        options = self._get_param("options")
        app = None
        app = self._find_app_by_any_value(application)
        if options is None:
            options = []
        new_resource = self.servicecatalogappregistry_backend.associate_resource(
            app.arn,
            resource_type,
            resource,
            options,
        )
        return json.dumps(new_resource)

    def list_associated_resources(self) -> str:
        application = unquote(self._get_param("application"))
        app = self._find_app_by_any_value(application)
        return_list = []
        for resource in app.associated_resources.values():
            return_list.append(resource.to_json())
        return json.dumps({"resources": return_list})

    def _find_app_by_any_value(self, search: str) -> Application:
        app = None
        if search in self.servicecatalogappregistry_backend.applications:
            app = self.servicecatalogappregistry_backend.applications[search]
        else:
            for a in self.servicecatalogappregistry_backend.applications.values():
                if search == a.id:
                    app = a
                elif search == a.name:
                    app = a
        if app is None:
            raise ValidationException
        return app
