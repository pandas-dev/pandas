# Copyright (C) 2011 Sebastian Rahlf <basti at redtoad dot de>
#
# This program is release under the MIT license. You can find the full text of
# the license in the LICENSE file.
import pytest


@pytest.fixture
def httpserver(request):
    """The returned ``httpserver`` provides a threaded HTTP server instance
    running on a randomly assigned port on localhost. It can be taught which
    content (i.e. string) to serve with which response code and comes with
    following attributes:

    * ``code`` - HTTP response code (int)
    * ``content`` - content of next response (str)
    * ``headers`` - response headers (dict)

    Once these attribute are set, all subsequent requests will be answered with
    these values until they are changed or the server is stopped. A more
    convenient way to change these is ::

        httpserver.serve_content(
            content='My content', code=200,
            headers={'content-type': 'text/plain'})

    The server address can be found in property

    * ``url``

    which is the string representation of tuple ``server_address`` (host as
    str, port as int).

    Example::

        import requests
        def scrape(url):
            html = requests.get(url).text
            # some parsing happens here
            # ...
            return result

        def test_retrieve_some_content(httpserver):
            httpserver.serve_content(open('cached-content.html').read())
            assert scrape(httpserver.url) == 'Found it!'

    """
    from pytest_localserver import http

    server = http.ContentServer()
    server.start()
    request.addfinalizer(server.stop)
    return server


@pytest.fixture
def httpsserver(request):
    """The returned ``httpsserver`` (note the additional S!) provides a
    threaded HTTP server instance similar to funcarg ``httpserver`` but with
    SSL encryption.
    """
    from pytest_localserver import https

    server = https.SecureContentServer()
    server.start()
    request.addfinalizer(server.stop)
    return server


@pytest.fixture
def smtpserver(request):
    """The returned ``smtpserver`` provides a threaded instance of
    ``smtpd.SMTPServer`` running on localhost.  It has the following
    attributes:

    * ``addr`` - server address as tuple (host as str, port as int)
    """
    from pytest_localserver import smtp

    server = smtp.Server()
    server.start()
    request.addfinalizer(server.stop)
    return server
