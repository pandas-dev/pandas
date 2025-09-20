from typing import Final, List

from moto.stepfunctions.parser.asl.component.component import Component


class CSVHeaders(Component):
    header_names: Final[List[str]]

    def __init__(self, header_names: List[str]):
        self.header_names = header_names
