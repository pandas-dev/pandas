import serial

class Serial(serial.Serial):
    def from_url(self, url: str) -> str: ...
