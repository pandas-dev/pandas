"""Handles incoming sesv2 requests, invokes methods, returns responses."""

import base64
import json
from typing import List
from urllib.parse import unquote

from moto.core.responses import BaseResponse

from ..ses.responses import SEND_EMAIL_RESPONSE
from .models import SESV2Backend, sesv2_backends


class SESV2Response(BaseResponse):
    """Handler for SESV2 requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="sesv2")

    @property
    def sesv2_backend(self) -> SESV2Backend:
        """Return backend instance specific for this region."""
        return sesv2_backends[self.current_account][self.region]

    def send_email(self) -> str:
        """Piggy back on functionality from v1 mostly"""

        params = json.loads(self.body)
        from_email_address = params.get("FromEmailAddress")
        destination = params.get("Destination", {})
        content = params.get("Content")
        if "Raw" in content:
            all_destinations: List[str] = []
            if "ToAddresses" in destination:
                all_destinations = all_destinations + destination["ToAddresses"]
            if "CcAddresses" in destination:
                all_destinations = all_destinations + destination["CcAddresses"]
            if "BccAddresses" in destination:
                all_destinations = all_destinations + destination["BccAddresses"]
            message = self.sesv2_backend.send_raw_email(
                source=from_email_address,
                destinations=all_destinations,
                raw_data=base64.b64decode(content["Raw"]["Data"]).decode("utf-8"),
            )
        elif "Simple" in content:
            content_body = content["Simple"]["Body"]
            if "Html" in content_body:
                body = content_body["Html"]["Data"]
            else:
                body = content_body["Text"]["Data"]
            message = self.sesv2_backend.send_email(  # type: ignore
                source=from_email_address,
                destinations=destination,
                subject=content["Simple"]["Subject"]["Data"],
                body=body,
            )
        elif "Template" in content:
            raise NotImplementedError("Template functionality not ready")

        # use v1 templates as response same in v1 and v2
        template = self.response_template(SEND_EMAIL_RESPONSE)
        return template.render(message=message)

    def create_contact_list(self) -> str:
        params = json.loads(self.body)
        self.sesv2_backend.create_contact_list(params)
        return json.dumps({})

    def get_contact_list(self) -> str:
        contact_list_name = self._get_param("ContactListName")
        contact_list = self.sesv2_backend.get_contact_list(contact_list_name)
        return json.dumps(contact_list.response_object)

    def list_contact_lists(self) -> str:
        contact_lists = self.sesv2_backend.list_contact_lists()
        return json.dumps(dict(ContactLists=[c.response_object for c in contact_lists]))

    def delete_contact_list(self) -> str:
        name = self._get_param("ContactListName")
        self.sesv2_backend.delete_contact_list(name)
        return json.dumps({})

    def create_contact(self) -> str:
        contact_list_name = self._get_param("ContactListName")
        params = json.loads(self.body)
        self.sesv2_backend.create_contact(contact_list_name, params)
        return json.dumps({})

    def get_contact(self) -> str:
        email = unquote(self._get_param("EmailAddress"))
        contact_list_name = self._get_param("ContactListName")
        contact = self.sesv2_backend.get_contact(email, contact_list_name)
        return json.dumps(contact.response_object)

    def list_contacts(self) -> str:
        contact_list_name = self._get_param("ContactListName")
        contacts = self.sesv2_backend.list_contacts(contact_list_name)
        return json.dumps(dict(Contacts=[c.response_object for c in contacts]))

    def delete_contact(self) -> str:
        email = self._get_param("EmailAddress")
        contact_list_name = self._get_param("ContactListName")
        self.sesv2_backend.delete_contact(unquote(email), contact_list_name)
        return json.dumps({})
