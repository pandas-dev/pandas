from vortexa_utils.aws.utils.dataclasses import nested_dataclass
from .action import Action
from .verdicts import (DKIMVerdict,
                       DMARCVerdict,
                       SPFVerdict,
                       SpamVerdict,
                       VirusVerdict)


@nested_dataclass
class Receipt:
    """SNS Recipt object.

    Attributes
    ----------
    action : Action
        Encapsulates information about the action that was executed.

    dkimVerdict : DKIMVerdict
        Indicates whether the DomainKeys Identified Mail (DKIM) check passed.

    dmarcPolicy : str
        Domain-based Message Authentication, Reporting & Conformance (DMARC)
        settings for the sending domain.
        This field only appears if the message fails DMARC authentication.
        Possible values for this field are:
          - none: no specific action be taken on messages that fail DMARC.
          - quarantine: messages that fail DMARC be treated as suspicious.
          - reject: messages that fail DMARC authentication be rejected.

    dmarcVerdict : DMARCVerdict
        Indicates whether the DMARC check passed.

    processingTimeMillis : str
        `str` specifies the period, in milliseconds, from the time Amazon SES
        received the message to the time it triggered the action.

    recipients : list[str]
        list of recipients that were matched by the active receipt rule.
        The addresses may differ from those listed by the destination field
        in the mail Object.

    spamVerdict : SpamVerdict
        Indicates whether the message is spam

    spfVerdict : SPFVerdict
        Whether the Sender Policy Framework (SPF) check passed

    timestamp : str
        ISO 8601 format string representing when the action was triggered.

    virusVerdict : VirusVerdict
        Whether the message contains a virus.
         For a list of possible values, see virusVerdict Object.
    """
    action: Action
    processingTimeMillis: str
    recipients: str
    timestamp: str
    dmarcPolicy: str = None
    dmarcVerdict: DMARCVerdict = None
    dkimVerdict: DKIMVerdict = None
    spamVerdict: SpamVerdict = None
    spfVerdict: SPFVerdict = None
    virusVerdict: VirusVerdict = None
