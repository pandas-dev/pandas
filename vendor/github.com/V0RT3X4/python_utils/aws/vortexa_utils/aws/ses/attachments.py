import email
from .application_mapper import application_mapping


class Attachment(object):

    def __init__(self, attachment: email.message.EmailMessage):
        self.attachment = attachment

    def to_df(self):
        content_type = self.attachment.get_content_type()
        reader = application_mapping.get(content_type)
        if reader is None:
            raise TypeError(f"unknown content_type {content_type}")
        return reader(self.attachment.get_content())
