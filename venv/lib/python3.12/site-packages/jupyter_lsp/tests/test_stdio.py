import asyncio
import subprocess
import sys

import pytest
from tornado.queues import Queue

from jupyter_lsp.stdio import LspStdIoReader

WRITER_TEMPLATE = """
from time import sleep

print('Content-Length: {length}')
print()

for repeat in range({repeats}):
    sleep({interval})
    print('{message}', end='')

if {add_excess}:
    print("extra", end='')

print()
"""


class CommunicatorSpawner:
    def __init__(self, tmp_path):
        self.tmp_path = tmp_path

    def spawn_writer(
        self, message: str, repeats: int = 1, interval=None, add_excess=False
    ):
        length = len(message) * repeats
        commands_file = self.tmp_path / "writer.py"
        commands_file.write_text(
            WRITER_TEMPLATE.format(
                length=length,
                repeats=repeats,
                interval=interval or 0,
                message=message,
                add_excess=add_excess,
            )
        )
        return subprocess.Popen(
            [sys.executable, "-u", str(commands_file)],
            stdout=subprocess.PIPE,
            bufsize=0,
        )


@pytest.fixture
def communicator_spawner(tmp_path):
    return CommunicatorSpawner(tmp_path)


async def join_process(process: subprocess.Popen, headstart=1, timeout=1):
    await asyncio.sleep(headstart)
    result = process.wait(timeout=timeout)
    if process.stdout:
        process.stdout.close()
    return result


@pytest.mark.parametrize(
    "message,repeats,interval,add_excess",
    [
        ["short", 1, None, False],
        ["ab" * 10_0000, 1, None, False],
        ["ab", 2, 0.01, False],
        ["ab", 45, 0.01, False],
        ["message", 2, 0.01, True],
    ],
    ids=["short", "long", "intermittent", "intensive-intermittent", "with-excess"],
)
@pytest.mark.asyncio
async def test_reader(message, repeats, interval, add_excess, communicator_spawner):
    queue = Queue()

    process = communicator_spawner.spawn_writer(
        message=message, repeats=repeats, interval=interval, add_excess=add_excess
    )
    reader = LspStdIoReader(stream=process.stdout, queue=queue)

    await asyncio.gather(join_process(process, headstart=3, timeout=1), reader.read())

    result = queue.get_nowait()
    assert result == message * repeats
