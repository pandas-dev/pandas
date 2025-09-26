from typing import Final

from moto.stepfunctions.parser.asl.component.intrinsic.member import Member


class MemberAccess(Member):
    def __init__(self, subject: Member, target: Member):
        self.subject: Final[Member] = subject
        self.target: Final[Member] = target
