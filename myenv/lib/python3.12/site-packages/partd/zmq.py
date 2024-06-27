import zmq
import logging
from itertools import chain
from bisect import bisect
import socket
from operator import add
from time import sleep, time
from toolz import accumulate, topk, pluck, merge, keymap
import uuid
from collections import defaultdict
from contextlib import contextmanager, suppress
from threading import Thread, Lock
from datetime import datetime
from multiprocessing import Process
import traceback
import sys
from .dict import Dict
from .file import File
from .buffer import Buffer
from . import core


tuple_sep = b'-|-'

logger = logging.getLogger(__name__)


@contextmanager
def logerrors():
    try:
        yield
    except Exception as e:
        logger.exception(e)
        raise


class Server:
    def __init__(self, partd=None, bind=None, start=True, block=False,
            hostname=None):
        self.context = zmq.Context()
        if partd is None:
            partd = Buffer(Dict(), File())
        self.partd = partd

        self.socket = self.context.socket(zmq.ROUTER)

        if hostname is None:
            hostname = socket.gethostname()
        if isinstance(bind, str):
            bind = bind.encode()
        if bind is None:
            port = self.socket.bind_to_random_port('tcp://*')
        else:
            self.socket.bind(bind)
            port = int(bind.split(':')[-1].rstrip('/'))
        self.address = ('tcp://%s:%d' % (hostname, port)).encode()

        self.status = 'created'

        self.partd.lock.acquire()
        self._lock = Lock()
        self._socket_lock = Lock()

        if start:
            self.start()

        if block:
            self.block()

    def start(self):
        if self.status != 'run':
            self.status = 'run'
            self._listen_thread = Thread(target=self.listen)
            self._listen_thread.start()
            logger.debug('Start server at %s', self.address)

    def block(self):
        """ Block until all threads close """
        try:
            self._listen_thread.join()
        except AttributeError:
            pass

    def listen(self):
        with logerrors():
            logger.debug('Start listening %s', self.address)
            while self.status != 'closed':
                if not self.socket.poll(100):
                    continue

                with self._socket_lock:
                    payload = self.socket.recv_multipart()

                address, command, payload = payload[0], payload[1], payload[2:]
                logger.debug('Server receives %s %s', address, command)
                if command == b'close':
                    logger.debug('Server closes')
                    self.ack(address)
                    self.status = 'closed'
                    break
                    # self.close()

                elif command == b'append':
                    keys, values = payload[::2], payload[1::2]
                    keys = list(map(deserialize_key, keys))
                    data = dict(zip(keys, values))
                    self.partd.append(data, lock=False)
                    logger.debug('Server appends %d keys', len(data))
                    self.ack(address)

                elif command == b'iset':
                    key, value = payload
                    key = deserialize_key(key)
                    self.partd.iset(key, value, lock=False)
                    self.ack(address)

                elif command == b'get':
                    keys = list(map(deserialize_key, payload))
                    logger.debug('get %s', keys)
                    result = self.get(keys)
                    self.send_to_client(address, result)
                    self.ack(address, flow_control=False)

                elif command == b'delete':
                    keys = list(map(deserialize_key, payload))
                    logger.debug('delete %s', keys)
                    self.partd.delete(keys, lock=False)
                    self.ack(address, flow_control=False)

                elif command == b'syn':
                    self.ack(address)

                elif command == b'drop':
                    self.drop()
                    self.ack(address)

                else:
                    logger.debug("Unknown command: %s", command)
                    raise ValueError("Unknown command: " + command)

    def send_to_client(self, address, result):
        with logerrors():
            if not isinstance(result, list):
                result = [result]
            with self._socket_lock:
                self.socket.send_multipart([address] + result)

    def ack(self, address, flow_control=True):
        with logerrors():
            logger.debug('Server sends ack')
            self.send_to_client(address, b'ack')

    def append(self, data):
        self.partd.append(data, lock=False)
        logger.debug('Server appends %d keys', len(data))

    def drop(self):
        with logerrors():
            self.partd.drop()

    def get(self, keys):
        with logerrors():
            logger.debug('Server gets keys: %s', keys)
            with self._lock:
                result = self.partd.get(keys, lock=False)
            return result

    def close(self):
        logger.debug('Server closes')
        self.status = 'closed'
        self.block()
        with suppress(zmq.error.ZMQError):
            self.socket.close(1)
        with suppress(zmq.error.ZMQError):
            self.context.destroy(3)
        self.partd.lock.release()

    def __enter__(self):
        self.start()
        return self

    def __exit__(self, *args):
        self.close()
        self.partd.__exit__(*args)


def keys_to_flush(lengths, fraction=0.1, maxcount=100000):
    """ Which keys to remove

    >>> lengths = {'a': 20, 'b': 10, 'c': 15, 'd': 15,
    ...            'e': 10, 'f': 25, 'g': 5}
    >>> keys_to_flush(lengths, 0.5)
    ['f', 'a']
    """
    top = topk(max(len(lengths) // 2, 1),
               lengths.items(),
               key=1)
    total = sum(lengths.values())
    cutoff = min(maxcount, max(1,
                   bisect(list(accumulate(add, pluck(1, top))),
                          total * fraction)))
    result = [k for k, v in top[:cutoff]]
    assert result
    return result


def serialize_key(key):
    """

    >>> serialize_key('x')
    b'x'
    >>> serialize_key(('a', 'b', 1))
    b'a-|-b-|-1'
    """
    if isinstance(key, tuple):
        return tuple_sep.join(map(serialize_key, key))
    if isinstance(key, bytes):
        return key
    if isinstance(key, str):
        return key.encode()
    return str(key).encode()


def deserialize_key(text):
    """

    >>> deserialize_key(b'x')
    b'x'
    >>> deserialize_key(b'a-|-b-|-1')
    (b'a', b'b', b'1')
    """
    if tuple_sep in text:
        return tuple(text.split(tuple_sep))
    else:
        return text


from .core import Interface
from .file import File


class Client(Interface):
    def __init__(self, address=None, create_server=False, **kwargs):
        self.address = address
        self.context = zmq.Context()
        self.socket = self.context.socket(zmq.DEALER)
        logger.debug('Client connects to %s', address)
        self.socket.connect(address)
        self.send(b'syn', [], ack_required=False)
        self.lock = NotALock()  # Server sequentializes everything
        Interface.__init__(self)

    def __getstate__(self):
        return {'address': self.address}

    def __setstate__(self, state):
        self.__init__(state['address'])
        logger.debug('Reconstruct client from pickled state')

    def send(self, command, payload, recv=False, ack_required=True):
        if ack_required:
            ack = self.socket.recv_multipart()
            assert ack == [b'ack']
        logger.debug('Client sends command: %s', command)
        self.socket.send_multipart([command] + payload)
        if recv:
            result = self.socket.recv_multipart()
        else:
            result = None
        return result

    def _get(self, keys, lock=None):
        """

        Lock argument is ignored.  Everything is sequential (I think)
        """
        logger.debug('Client gets %s %s', self.address, keys)
        keys = list(map(serialize_key, keys))
        return self.send(b'get', keys, recv=True)

    def append(self, data, lock=None):
        logger.debug('Client appends %s %s', self.address, str(len(data)) + ' keys')
        data = keymap(serialize_key, data)
        payload = list(chain.from_iterable(data.items()))
        self.send(b'append', payload)

    def _delete(self, keys, lock=None):
        logger.debug('Client deletes %s %s', self.address, str(len(keys)) + ' keys')
        keys = list(map(serialize_key, keys))
        self.send(b'delete', keys)

    def _iset(self, key, value):
        self.send(b'iset', [serialize_key(key), value])

    def drop(self):
        self.send(b'drop', [])
        sleep(0.05)

    def close_server(self):
        self.send(b'close', [])

    def close(self):
        if hasattr(self, 'server_process'):
            with suppress(zmq.error.ZMQError):
                self.close_server()
            self.server_process.join()
        with suppress(zmq.error.ZMQError):
            self.socket.close(1)
        with suppress(zmq.error.ZMQError):
            self.context.destroy(1)

    def __exit__(self, type, value, traceback):
        self.drop()
        self.close()

    def __del__(self):
        self.close()


class NotALock:
    def acquire(self): pass
    def release(self): pass

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass
