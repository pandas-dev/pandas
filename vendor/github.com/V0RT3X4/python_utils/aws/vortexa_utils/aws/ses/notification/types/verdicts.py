from dataclasses import dataclass


@dataclass
class Verdict(object):
    """Verdict object.

    Attributes
    ----------
    status : str
        String that contains the verdict. Possible values are:
          - PASS: The message passed the given test.
          - FAIL: The message failed the given test.
          - GRAY: The message failed the given test,
          - PROCESSING_FAILED: There is an issue that prevents Amazon SES
                from providing a verdict to the given test.
    """
    status: str


@dataclass
class DKIMVerdict(Verdict):
    ...


@dataclass
class DMARCVerdict(Verdict):
    ...


@dataclass
class SpamVerdict(Verdict):
    ...


@dataclass
class SPFVerdict(Verdict):
    ...


@dataclass
class VirusVerdict(Verdict):
    ...
