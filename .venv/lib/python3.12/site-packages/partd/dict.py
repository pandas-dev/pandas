from .core import Interface
from threading import Lock


class Dict(Interface):
    def __init__(self):
        self.lock = Lock()
        self.data = dict()
        Interface.__init__(self)

    def __getstate__(self):
        return {'data': self.data}

    def __setstate__(self, state):
        Interface.__setstate__(self, state)
        Dict.__init__(self)
        self.data = state['data']

    def append(self, data, lock=True, **kwargs):
        if lock: self.lock.acquire()
        try:
            for k, v in data.items():
                if k not in self.data:
                    self.data[k] = []
                self.data[k].append(v)
        finally:
            if lock: self.lock.release()

    def _get(self, keys, lock=True, **kwargs):
        assert isinstance(keys, (list, tuple, set))
        if lock:
            self.lock.acquire()
        try:
            result = [b''.join(self.data.get(key, [])) for key in keys]
        finally:
            if lock:
                self.lock.release()
        return result

    def _iset(self, key, value, lock=True):
        """ Idempotent set """
        if lock:
            self.lock.acquire()
        try:
            self.data[key] = [value]
        finally:
            if lock:
                self.lock.release()

    def _delete(self, keys, lock=True):
        if lock:
            self.lock.acquire()
        try:
            for key in keys:
                if key in self.data:
                    del self.data[key]
        finally:
            if lock:
                self.lock.release()

    def drop(self):
        self._iset_seen.clear()
        self.data.clear()

    def __exit__(self, *args):
        self.drop()
