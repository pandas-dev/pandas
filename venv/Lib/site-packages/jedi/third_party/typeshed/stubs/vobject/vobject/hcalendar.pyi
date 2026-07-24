from .icalendar import VCalendar2_0

class HCalendar(VCalendar2_0):
    name: str
    @classmethod
    def serialize(cls, obj, buf=None, lineLength=None, validate: bool = True): ...
