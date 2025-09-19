from threading import Timer as ThreadingTimer

from moto.core.common_models import BaseModel


class Timer(BaseModel):
    def __init__(self, background_timer: ThreadingTimer, started_event_id: int):
        self.background_timer = background_timer
        self.started_event_id = started_event_id

    def start(self) -> None:
        return self.background_timer.start()

    def is_alive(self) -> bool:
        return self.background_timer.is_alive()

    def cancel(self) -> None:
        return self.background_timer.cancel()
