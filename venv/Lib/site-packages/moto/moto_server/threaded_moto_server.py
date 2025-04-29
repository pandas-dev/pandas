from threading import Event, Thread
from typing import Optional, Tuple

from werkzeug.serving import BaseWSGIServer, make_server

from .werkzeug_app import DomainDispatcherApplication, create_backend_app


class ThreadedMotoServer:
    def __init__(
        self, ip_address: str = "0.0.0.0", port: int = 5000, verbose: bool = True
    ):
        self._port = port

        self._thread: Optional[Thread] = None
        self._ip_address = ip_address
        self._server: Optional[BaseWSGIServer] = None
        self._server_ready_event = Event()
        self._verbose = verbose

    def _server_entry(self) -> None:
        app = DomainDispatcherApplication(create_backend_app)

        self._server = make_server(self._ip_address, self._port, app, True)
        self._server_ready_event.set()
        self._server.serve_forever()

    def start(self) -> None:
        if self._verbose:
            print(  # noqa
                f"Starting a new Thread with MotoServer running on {self._ip_address}:{self._port}..."
            )
        self._thread = Thread(target=self._server_entry, daemon=True)
        self._thread.start()
        self._server_ready_event.wait()

    def get_host_and_port(self) -> Tuple[str, int]:
        assert self._server is not None, "Make sure to call start() first"
        host, port = self._server.server_address
        return (str(host), port)

    def stop(self) -> None:
        self._server_ready_event.clear()
        if self._server:
            self._server.shutdown()

        self._thread.join()  # type: ignore[union-attr]
