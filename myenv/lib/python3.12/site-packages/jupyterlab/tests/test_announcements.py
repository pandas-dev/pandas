# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

import hashlib
import json
from unittest.mock import patch

from . import fake_client_factory

FAKE_ATOM_FEED = b"""<?xml version="1.0" encoding="utf-8"?><feed xmlns="http://www.w3.org/2005/Atom" ><generator uri="https://jekyllrb.com/" version="3.9.2">Jekyll</generator><link href="https://jupyterlab.github.io/assets/feed.xml" rel="self" type="application/atom+xml" /><link href="https://jupyterlab.github.io/assets/" rel="alternate" type="text/html" /><updated>2022-11-02T15:14:50+00:00</updated><id>https://jupyterlab.github.io/assets/feed.xml</id><title type="html">JupyterLab News</title><subtitle>Subscribe to get news about JupyterLab.</subtitle><entry><title type="html">Thanks for using JupyterLab</title><link href="https://jupyterlab.github.io/assets/posts/2022/11/02/demo.html" rel="alternate" type="text/html" title="Thanks for using JupyterLab" /><published>2022-11-02T14:00:00+00:00</published><updated>2022-11-02T14:00:00+00:00</updated><id>https://jupyterlab.github.io/assets/posts/2022/11/02/demo</id><content type="html" xml:base="https://jupyterlab.github.io/assets/posts/2022/11/02/demo.html">&lt;h1 id=&quot;welcome&quot;&gt;Welcome&lt;/h1&gt;

&lt;p&gt;Thanks a lot for your interest in JupyterLab.&lt;/p&gt;</content><author><name></name></author><category term="posts" /><summary type="html">Big thanks to you, beloved JupyterLab user.</summary></entry></feed>"""

FAKE_JUPYTERLAB_PYPI_JSON = b"""{ "info": { "version": "1000.0.0" } }"""


@patch("tornado.httpclient.AsyncHTTPClient", new_callable=fake_client_factory)
async def test_NewsHandler_get_success(mock_client, labserverapp, jp_fetch):
    mock_client.body = FAKE_ATOM_FEED

    response = await jp_fetch("lab", "api", "news", method="GET")

    assert response.code == 200
    payload = json.loads(response.body)
    assert payload["news"] == [
        {
            "createdAt": 1667397600000.0,
            "message": "Thanks for using JupyterLab\nBig thanks to you, beloved JupyterLab user.",
            "modifiedAt": 1667397600000.0,
            "type": "info",
            "link": [
                "Open full post",
                "https://jupyterlab.github.io/assets/posts/2022/11/02/demo.html",
            ],
            "options": {
                "data": {
                    "id": "https://jupyterlab.github.io/assets/posts/2022/11/02/demo",
                    "tags": ["news"],
                }
            },
        }
    ]


@patch("tornado.httpclient.AsyncHTTPClient", new_callable=fake_client_factory)
async def test_NewsHandler_get_failure(mock_client, labserverapp, jp_fetch):
    response = await jp_fetch("lab", "api", "news", method="GET")

    assert response.code == 200
    payload = json.loads(response.body)
    assert payload["news"] == []


@patch("tornado.httpclient.AsyncHTTPClient", new_callable=fake_client_factory)
async def test_CheckForUpdateHandler_get_pypi_success(mock_client, labserverapp, jp_fetch):
    mock_client.body = FAKE_JUPYTERLAB_PYPI_JSON

    response = await jp_fetch("lab", "api", "update", method="GET")

    message = "A newer version (1000.0.0) of JupyterLab is available."
    assert response.code == 200
    payload = json.loads(response.body)
    assert payload["notification"]["message"] == message
    assert payload["notification"]["link"] == [
        "Open changelog",
        "https://github.com/jupyterlab/jupyterlab/releases/tag/v1000.0.0",
    ]
    assert payload["notification"]["options"] == {
        "data": {"id": hashlib.sha1(message.encode()).hexdigest(), "tags": ["update"]}  # noqa: S324
    }


@patch("tornado.httpclient.AsyncHTTPClient", new_callable=fake_client_factory)
async def test_CheckForUpdateHandler_get_failure(mock_client, labserverapp, jp_fetch):
    response = await jp_fetch("lab", "api", "update", method="GET")

    assert response.code == 200
    payload = json.loads(response.body)
    assert payload["notification"] is None


FAKE_NO_SUMMARY_ATOM_FEED = b"""<?xml version='1.0' encoding='UTF-8'?><feed xmlns="http://www.w3.org/2005/Atom" xml:lang="en"><id>https://jupyterlab.github.io/assets/feed.xml</id><title>JupyterLab News</title><updated>2023-05-02T19:01:33.669598+00:00</updated><author><name>John Doe</name><email>john@example.de</email></author><link href="https://jupyterlab.github.io/assets/feed.xml" rel="self" type="application/atom+xml"/><link href="https://jupyterlab.github.io/assets/" rel="alternate" type="text/html"/><generator uri="https://lkiesow.github.io/python-feedgen" version="0.9.0">python-feedgen</generator><logo>http://ex.com/logo.jpg</logo><subtitle>Subscribe to get news about JupyterLab.</subtitle><entry><id>https://jupyterlab.github.io/assets/posts/2022/11/02/demo</id><title>Thanks for using JupyterLab</title><updated>2022-11-02T14:00:00+00:00</updated><link href="https://jupyterlab.github.io/assets/posts/2022/11/02/demo.html" rel="alternate" type="text/html" title="Thanks for using JupyterLab"/><published>2022-11-02T14:00:00+00:00</published></entry></feed>"""


@patch("tornado.httpclient.AsyncHTTPClient", new_callable=fake_client_factory)
async def test_NewsHandler_get_missing_summary(mock_client, labserverapp, jp_fetch):
    mock_client.body = FAKE_NO_SUMMARY_ATOM_FEED

    response = await jp_fetch("lab", "api", "news", method="GET")

    assert response.code == 200
    payload = json.loads(response.body)
    assert payload["news"] == [
        {
            "createdAt": 1667397600000.0,
            "message": "Thanks for using JupyterLab",
            "modifiedAt": 1667397600000.0,
            "type": "info",
            "link": [
                "Open full post",
                "https://jupyterlab.github.io/assets/posts/2022/11/02/demo.html",
            ],
            "options": {
                "data": {
                    "id": "https://jupyterlab.github.io/assets/posts/2022/11/02/demo",
                    "tags": ["news"],
                }
            },
        }
    ]


FAKE_MULTI_ENTRY_LINKS_ATOM_FEED = b"""<?xml version='1.0' encoding='UTF-8'?><feed xmlns="http://www.w3.org/2005/Atom" xml:lang="en"><id>https://jupyterlab.github.io/assets/feed.xml</id><title>JupyterLab News</title><updated>2023-05-02T19:59:44.332080+00:00</updated><author><name>John Doe</name><email>john@example.de</email></author><link href="https://jupyterlab.github.io/assets/feed.xml" rel="self" type="application/atom+xml"/><link href="https://jupyterlab.github.io/assets/" rel="alternate" type="text/html"/><generator uri="https://lkiesow.github.io/python-feedgen" version="0.9.0">python-feedgen</generator><logo>http://ex.com/logo.jpg</logo><subtitle>Subscribe to get news about JupyterLab.</subtitle><entry><id>https://jupyterlab.github.io/assets/posts/2022/11/02/demo</id><title>Thanks for using JupyterLab</title><updated>2022-11-02T14:00:00+00:00</updated><link href="https://jupyterlab.github.io/assets/posts/2022/11/02/demo_self.html" rel="self" type="text/html" title="Thanks for using JupyterLab"/><link href="https://jupyterlab.github.io/assets/posts/2022/11/02/demo.html" rel="alternate" type="text/html" title="Thanks for using JupyterLab"/><summary>Big thanks to you, beloved JupyterLab user.</summary><published>2022-11-02T14:00:00+00:00</published></entry></feed>"""


@patch("tornado.httpclient.AsyncHTTPClient", new_callable=fake_client_factory)
async def test_NewsHandler_multi_entry_links(mock_client, labserverapp, jp_fetch):
    mock_client.body = FAKE_MULTI_ENTRY_LINKS_ATOM_FEED

    response = await jp_fetch("lab", "api", "news", method="GET")

    assert response.code == 200
    payload = json.loads(response.body)
    assert payload["news"] == [
        {
            "createdAt": 1667397600000.0,
            "message": "Thanks for using JupyterLab\nBig thanks to you, beloved JupyterLab user.",
            "modifiedAt": 1667397600000.0,
            "type": "info",
            "link": [
                "Open full post",
                "https://jupyterlab.github.io/assets/posts/2022/11/02/demo.html",
            ],
            "options": {
                "data": {
                    "id": "https://jupyterlab.github.io/assets/posts/2022/11/02/demo",
                    "tags": ["news"],
                }
            },
        }
    ]


FAKE_NO_PUBLISHED_ATOM_FEED = b"""<?xml version='1.0' encoding='UTF-8'?><feed xmlns="http://www.w3.org/2005/Atom" xml:lang="en"><id>https://jupyterlab.github.io/assets/feed.xml</id><title>JupyterLab News</title><updated>2023-05-02T19:32:08.566055+00:00</updated><author><name>John Doe</name><email>john@example.de</email></author><link href="https://jupyterlab.github.io/assets/feed.xml" rel="self" type="application/atom+xml"/><link href="https://jupyterlab.github.io/assets/" rel="alternate" type="text/html"/><generator uri="https://lkiesow.github.io/python-feedgen" version="0.9.0">python-feedgen</generator><logo>http://ex.com/logo.jpg</logo><subtitle>Subscribe to get news about JupyterLab.</subtitle><entry><id>https://jupyterlab.github.io/assets/posts/2022/11/02/demo</id><title>Thanks for using JupyterLab</title><updated>2022-11-02T14:00:00+00:00</updated><link href="https://jupyterlab.github.io/assets/posts/2022/11/02/demo.html" rel="alternate" type="text/html" title="Thanks for using JupyterLab"/><summary>Big thanks to you, beloved JupyterLab user.</summary></entry></feed>"""


@patch("tornado.httpclient.AsyncHTTPClient", new_callable=fake_client_factory)
async def test_NewsHandler_no_published(mock_client, labserverapp, jp_fetch):
    mock_client.body = FAKE_NO_PUBLISHED_ATOM_FEED

    response = await jp_fetch("lab", "api", "news", method="GET")

    assert response.code == 200
    payload = json.loads(response.body)
    assert payload["news"] == [
        {
            "createdAt": 1667397600000.0,
            "message": "Thanks for using JupyterLab\nBig thanks to you, beloved JupyterLab user.",
            "modifiedAt": 1667397600000.0,
            "type": "info",
            "link": [
                "Open full post",
                "https://jupyterlab.github.io/assets/posts/2022/11/02/demo.html",
            ],
            "options": {
                "data": {
                    "id": "https://jupyterlab.github.io/assets/posts/2022/11/02/demo",
                    "tags": ["news"],
                }
            },
        }
    ]


FAKE_LINK_NO_REL_ATOM_FEED = b"""<?xml version='1.0' encoding='UTF-8'?><feed xmlns="http://www.w3.org/2005/Atom" xml:lang="en"><id>https://jupyterlab.github.io/assets/feed.xml</id><title>JupyterLab News</title><updated>2023-05-03T17:06:43.950978+00:00</updated><author><name>John Doe</name><email>john@example.de</email></author><link href="https://jupyterlab.github.io/assets/feed.xml" rel="self" type="application/atom+xml"/><link href="https://jupyterlab.github.io/assets/" rel="alternate" type="text/html"/><generator uri="https://lkiesow.github.io/python-feedgen" version="0.9.0">python-feedgen</generator><logo>http://ex.com/logo.jpg</logo><subtitle>Subscribe to get news about JupyterLab.</subtitle><entry><id>https://jupyterlab.github.io/assets/posts/2022/11/02/demo</id><title>Thanks for using JupyterLab</title><updated>2022-11-02T14:00:00+00:00</updated><link href="https://jupyterlab.github.io/assets/posts/2022/11/02/demo.html" type="text/html" title="Thanks for using JupyterLab"/><summary>Big thanks to you, beloved JupyterLab user.</summary><published>2022-11-02T14:00:00+00:00</published></entry></feed>"""


@patch("tornado.httpclient.AsyncHTTPClient", new_callable=fake_client_factory)
async def test_NewsHandler_link_no_rel(mock_client, labserverapp, jp_fetch):
    mock_client.body = FAKE_LINK_NO_REL_ATOM_FEED

    response = await jp_fetch("lab", "api", "news", method="GET")

    assert response.code == 200
    payload = json.loads(response.body)
    assert payload["news"] == [
        {
            "createdAt": 1667397600000.0,
            "message": "Thanks for using JupyterLab\nBig thanks to you, beloved JupyterLab user.",
            "modifiedAt": 1667397600000.0,
            "type": "info",
            "link": [
                "Open full post",
                "https://jupyterlab.github.io/assets/posts/2022/11/02/demo.html",
            ],
            "options": {
                "data": {
                    "id": "https://jupyterlab.github.io/assets/posts/2022/11/02/demo",
                    "tags": ["news"],
                }
            },
        }
    ]


FAKE_NO_LINK_ATOM_FEED = b"""<?xml version='1.0' encoding='UTF-8'?><feed xmlns="http://www.w3.org/2005/Atom" xml:lang="en"><id>https://jupyterlab.github.io/assets/feed.xml</id><title>JupyterLab News</title><updated>2023-05-03T17:06:43.950978+00:00</updated><author><name>John Doe</name><email>john@example.de</email></author><link href="https://jupyterlab.github.io/assets/feed.xml" rel="self" type="application/atom+xml"/><link href="https://jupyterlab.github.io/assets/" rel="alternate" type="text/html"/><generator uri="https://lkiesow.github.io/python-feedgen" version="0.9.0">python-feedgen</generator><logo>http://ex.com/logo.jpg</logo><subtitle>Subscribe to get news about JupyterLab.</subtitle><entry><id>https://jupyterlab.github.io/assets/posts/2022/11/02/demo</id><title>Thanks for using JupyterLab</title><updated>2022-11-02T14:00:00+00:00</updated><summary>Big thanks to you, beloved JupyterLab user.</summary><published>2022-11-02T14:00:00+00:00</published></entry></feed>"""


@patch("tornado.httpclient.AsyncHTTPClient", new_callable=fake_client_factory)
async def test_NewsHandler_no_links(mock_client, labserverapp, jp_fetch):
    mock_client.body = FAKE_NO_LINK_ATOM_FEED

    response = await jp_fetch("lab", "api", "news", method="GET")

    assert response.code == 200
    payload = json.loads(response.body)
    assert payload["news"] == [
        {
            "createdAt": 1667397600000.0,
            "message": "Thanks for using JupyterLab\nBig thanks to you, beloved JupyterLab user.",
            "modifiedAt": 1667397600000.0,
            "type": "info",
            "link": None,
            "options": {
                "data": {
                    "id": "https://jupyterlab.github.io/assets/posts/2022/11/02/demo",
                    "tags": ["news"],
                }
            },
        }
    ]
