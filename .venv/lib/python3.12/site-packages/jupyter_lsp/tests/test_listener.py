import asyncio
import re

import pytest
import traitlets
from tornado.queues import Queue

from jupyter_lsp import lsp_message_listener


@pytest.mark.parametrize("bad_string", ["not-a-function", "jupyter_lsp.__version__"])
@pytest.mark.asyncio
async def test_listener_bad_traitlets(bad_string, handlers):
    handler, ws_handler = handlers
    manager = handler.manager

    with pytest.raises(traitlets.TraitError):
        manager.all_listeners = [bad_string]


@pytest.mark.asyncio
async def test_listeners(known_server, handlers, jsonrpc_init_msg):
    """will some listeners listen?"""
    handler, ws_handler = handlers
    manager = handler.manager

    manager.all_listeners = ["jupyter_lsp.tests.listener.dummy_listener"]

    manager.initialize()
    manager._listeners["client"] = []  # hide predefined client listeners

    assert len(manager._listeners["all"]) == 1

    dummy_listener = manager._listeners["all"][0]
    assert re.match(
        (
            "<MessageListener listener=<function dummy_listener at .*?>,"
            " method=None, language_server=None>"
        ),
        repr(dummy_listener),
    )

    handler_listened = Queue()
    server_listened = Queue()
    all_listened = Queue()

    # some client listeners
    @lsp_message_listener("client", language_server=known_server, method="initialize")
    async def client_listener(scope, message, language_server, manager):
        await handler_listened.put(message)

    @lsp_message_listener("client", method=r"not-a-method")
    async def other_client_listener(
        scope, message, language_server, manager
    ):  # pragma: no cover
        await handler_listened.put(message)
        raise NotImplementedError("shouldn't get here")

    # some server listeners
    @lsp_message_listener("server", language_server=None, method=None)
    async def server_listener(scope, message, language_server, manager):
        await server_listened.put(message)

    @lsp_message_listener("server", language_server=r"not-a-language-server")
    async def other_server_listener(
        scope, message, language_server, manager
    ):  # pragma: no cover
        await handler_listened.put(message)
        raise NotImplementedError("shouldn't get here")

    # an all listener
    @lsp_message_listener("all")
    async def all_listener(
        scope, message, language_server, manager
    ):  # pragma: no cover
        await all_listened.put(message)

    assert len(manager._listeners["server"]) == 2
    assert len(manager._listeners["client"]) == 2
    assert len(manager._listeners["all"]) == 2

    await ws_handler.open(known_server)

    await ws_handler.on_message(jsonrpc_init_msg)

    results = await asyncio.wait_for(
        asyncio.gather(
            handler_listened.get(),
            server_listened.get(),
            all_listened.get(),
            all_listened.get(),
            return_exceptions=True,
        ),
        240 if known_server == "julia-language-server" else 20,
    )
    assert all([isinstance(res, dict) for res in results])

    ws_handler.on_close()

    handler_listened.task_done()
    server_listened.task_done()
    all_listened.task_done()
    all_listened.task_done()

    [
        manager.unregister_message_listener(listener)
        for listener in [
            client_listener,
            other_client_listener,
            server_listener,
            other_server_listener,
            all_listener,
        ]
    ]

    assert not manager._listeners["server"]
    assert not manager._listeners["client"]
    assert len(manager._listeners["all"]) == 1
