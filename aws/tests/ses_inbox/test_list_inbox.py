# cd aws/vortexa_utils
# cd ..
from typing import Iterable
from vortexa_utils.aws.ses.inbox import Inbox
from email.message import EmailMessage
from itertools import islice


Path = 'incoming_email/'

inbox = Inbox(default_bucket='ops-data.incoming-emails')


def test_list_inbox():
    inbox = Inbox(default_bucket='ops-data.incoming-emails')
    emails: Iterable[EmailMessage] = islice(
        inbox.list_emails(Path=Path),
        10
    )

    for email in emails:
        # print(email.as_string())
        attachments = list(email.iter_attachments())
        print(list(a.get_filename() for a in attachments))
        print(list(a.get_content_type() for a in attachments))
