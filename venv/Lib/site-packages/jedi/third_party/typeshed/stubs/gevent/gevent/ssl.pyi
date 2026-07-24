import sys
from _typeshed import StrOrBytesPath
from ssl import *

import gevent.socket

# for simplicity we trust that gevent's implementation matches the stdlib version exactly
# for the most part they just copy all the symbols anyways and re-implment the few that
# need to work differently. The only potentially problematic symbol is SSLSocket, since
# it derives from gevent's socket, rather than the stdlib one. SSLContext derives from
# the stdlib SSLContext. Since we already punted on socket, we don't need to change
# anything here either, until we decide that we can't punt on socket.

if sys.version_info >= (3, 12):
    # FIXME: wrap_socket has been removed in 3.12, gevent implements its own, so it
    #        will probably still be there in 3.12, but until we stub out gevent.ssl
    #        properly we will have to just pretend it still exists
    def wrap_socket(
        sock: gevent.socket.socket,
        keyfile: StrOrBytesPath | None = None,
        certfile: StrOrBytesPath | None = None,
        server_side: bool = False,
        cert_reqs: int = ...,
        ssl_version: int = ...,
        ca_certs: str | None = None,
        do_handshake_on_connect: bool = True,
        suppress_ragged_eofs: bool = True,
        ciphers: str | None = None,
    ) -> SSLSocket: ...
