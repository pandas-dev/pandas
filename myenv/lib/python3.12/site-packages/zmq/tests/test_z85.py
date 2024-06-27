"""Test Z85 encoding

confirm values and roundtrip with test values from the reference implementation.
"""

# Copyright (C) PyZMQ Developers
# Distributed under the terms of the Modified BSD License.

from unittest import TestCase

from zmq.utils import z85


class TestZ85(TestCase):
    def test_client_public(self):
        client_public = (
            b"\xbb\x88\x47\x1d\x65\xe2\x65\x9b"
            b"\x30\xc5\x5a\x53\x21\xce\xbb\x5a"
            b"\xab\x2b\x70\xa3\x98\x64\x5c\x26"
            b"\xdc\xa2\xb2\xfc\xb4\x3f\xc5\x18"
        )
        encoded = z85.encode(client_public)

        assert encoded == b"Yne@$w-vo<fVvi]a<NY6T1ed:M$fCG*[IaLV{hID"
        decoded = z85.decode(encoded)
        assert decoded == client_public

    def test_client_secret(self):
        client_secret = (
            b"\x7b\xb8\x64\xb4\x89\xaf\xa3\x67"
            b"\x1f\xbe\x69\x10\x1f\x94\xb3\x89"
            b"\x72\xf2\x48\x16\xdf\xb0\x1b\x51"
            b"\x65\x6b\x3f\xec\x8d\xfd\x08\x88"
        )
        encoded = z85.encode(client_secret)

        assert encoded == b"D:)Q[IlAW!ahhC2ac:9*A}h:p?([4%wOTJ%JR%cs"
        decoded = z85.decode(encoded)
        assert decoded == client_secret

    def test_server_public(self):
        server_public = (
            b"\x54\xfc\xba\x24\xe9\x32\x49\x96"
            b"\x93\x16\xfb\x61\x7c\x87\x2b\xb0"
            b"\xc1\xd1\xff\x14\x80\x04\x27\xc5"
            b"\x94\xcb\xfa\xcf\x1b\xc2\xd6\x52"
        )
        encoded = z85.encode(server_public)

        assert encoded == b"rq:rM>}U?@Lns47E1%kR.o@n%FcmmsL/@{H8]yf7"
        decoded = z85.decode(encoded)
        assert decoded == server_public

    def test_server_secret(self):
        server_secret = (
            b"\x8e\x0b\xdd\x69\x76\x28\xb9\x1d"
            b"\x8f\x24\x55\x87\xee\x95\xc5\xb0"
            b"\x4d\x48\x96\x3f\x79\x25\x98\x77"
            b"\xb4\x9c\xd9\x06\x3a\xea\xd3\xb7"
        )
        encoded = z85.encode(server_secret)

        assert encoded == b"JTKVSB%%)wK0E.X)V>+}o?pNmC{O&4W4b!Ni{Lh6"
        decoded = z85.decode(encoded)
        assert decoded == server_secret
