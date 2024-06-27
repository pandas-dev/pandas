# Copyright (C) 2011 Sebastian Rahlf <basti at redtoad dot de>
# with some ideas from http://code.activestate.com/recipes/440690/
# SmtpMailsink Copyright 2005 Aviarc Corporation
# Written by Adam Feuer, Matt Branthwaite, and Troy Frever
# which is Licensed under the PSF License
import email

import aiosmtpd.controller


class MessageDetails:
    def __init__(self, peer, mailfrom, rcpttos, *, mail_options=None, rcpt_options=None):
        self.peer = peer
        self.mailfrom = mailfrom
        self.rcpttos = rcpttos
        if mail_options:
            self.mail_options = mail_options
        if rcpt_options:
            self.rcpt_options = rcpt_options


class Handler:
    def __init__(self):
        self.outbox = []

    async def handle_DATA(self, server, session, envelope):
        message = email.message_from_bytes(envelope.content)
        message.details = MessageDetails(session.peer, envelope.mail_from, envelope.rcpt_tos)
        self.outbox.append(message)
        return "250 OK"


class Server(aiosmtpd.controller.Controller):

    """
    Small SMTP test server.

    This is little more than a wrapper around aiosmtpd.controller.Controller
    which offers a slightly different interface for backward compatibility with
    earlier versions of pytest-localserver. You can just as well use a standard
    Controller and pass it a Handler instance.

    Here is how to use this class for sending an email, if you really need to::

        server = Server(port=8080)
        server.start()
        print 'SMTP server is running on %s:%i' % server.addr

        # any e-mail sent to localhost:8080 will end up in server.outbox
        # ...

        server.stop()

    """

    def __init__(self, host="localhost", port=0):
        try:
            super().__init__(Handler(), hostname=host, port=port, server_hostname=host)
        except TypeError:
            # for aiosmtpd <1.3
            super().__init__(Handler(), hostname=host, port=port)

    @property
    def outbox(self):
        return self.handler.outbox

    def _set_server_socket_attributes(self):
        """
        Set the addr and port attributes on this Server instance, if they're not
        already set.
        """

        # I split this out into its own method to allow running this code in
        # aiosmtpd <1.4, which doesn't have the _trigger_server() method on
        # the Controller class. If I put it directly in _trigger_server(), it
        # would fail when calling super()._trigger_server(). In the future, when
        # we can safely require aiosmtpd >=1.4, this method can be inlined
        # directly into _trigger_server().
        if hasattr(self, "addr"):
            assert hasattr(self, "port")
            return

        self.addr = self.server.sockets[0].getsockname()[:2]

        # Work around a bug/missing feature in aiosmtpd (https://github.com/aio-libs/aiosmtpd/issues/276)
        if self.port == 0:
            self.port = self.addr[1]
            assert self.port != 0

    def _trigger_server(self):
        self._set_server_socket_attributes()
        super()._trigger_server()

    def is_alive(self):
        return self._thread is not None and self._thread.is_alive()

    @property
    def accepting(self):
        try:
            return self.server.is_serving()
        except AttributeError:
            # asyncio.base_events.Server.is_serving() only exists in Python 3.6
            # and up. For Python 3.5, asyncio.base_events.BaseEventLoop.is_running()
            # is a close approximation; it should mostly return the same value
            # except for brief periods when the server is starting up or shutting
            # down. Once we drop support for Python 3.5, this branch becomes
            # unnecessary.
            return self.loop.is_running()

    # for aiosmtpd <1.4
    if not hasattr(aiosmtpd.controller.Controller, "_trigger_server"):

        def start(self):
            super().start()
            self._set_server_socket_attributes()

    def stop(self, timeout=None):
        """
        Stops test server.
        :param timeout: When the timeout argument is present and not None, it
        should be a floating point number specifying a timeout for the
        operation in seconds (or fractions thereof).
        """

        # This mostly copies the implementation from Controller.stop(), with two
        # differences:
        # - It removes the assertion that the thread exists, allowing stop() to
        #   be called more than once safely
        # - It passes the timeout argument to Thread.join()
        if self.loop.is_running():
            try:
                self.loop.call_soon_threadsafe(self.cancel_tasks)
            except AttributeError:
                # for aiosmtpd < 1.4.3
                self.loop.call_soon_threadsafe(self._stop)
        if self._thread is not None:
            self._thread.join(timeout)
            self._thread = None
        self._thread_exception = None
        self._factory_invoked = None
        self.server_coro = None
        self.server = None
        self.smtpd = None

    def __del__(self):
        # This is just for backward compatibility, to preserve the behavior that
        # the server is stopped when this object is finalized. But it seems
        # sketchy to rely on this to stop the server. Typically, the server
        # should be stopped "manually", before it gets deleted.
        if self.is_alive():
            self.stop()

    def __repr__(self):  # pragma: no cover
        return "<smtp.Server %s:%s>" % self.addr


def main():
    import time

    server = Server()
    server.start()

    print("SMTP server is running on %s:%i" % server.addr)
    print("Type <Ctrl-C> to stop")

    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        pass
    finally:
        print("\rstopping...")
        server.stop()


if __name__ == "__main__":  # pragma: no cover
    main()
