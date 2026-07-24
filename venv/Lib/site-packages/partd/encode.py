from .core import Interface
from .file import File
from toolz import valmap
from .utils import frame, framesplit


class Encode(Interface):
    def __init__(self, encode, decode, join, partd=None):
        if not partd or isinstance(partd, str):
            partd = File(partd)
        self.partd = partd
        self.encode = encode
        self.decode = decode
        self.join = join
        Interface.__init__(self)

    def __getstate__(self):
        return self.__dict__

    __setstate__ = Interface.__setstate__

    def append(self, data, **kwargs):
        data = valmap(self.encode, data)
        data = valmap(frame, data)
        self.partd.append(data, **kwargs)

    def _get(self, keys, **kwargs):
        raw = self.partd._get(keys, **kwargs)
        return [self.join([self.decode(frame) for frame in framesplit(chunk)])
                for chunk in raw]

    def delete(self, keys, **kwargs):
        return self.partd.delete(keys, **kwargs)

    def _iset(self, key, value, **kwargs):
        return self.partd.iset(key, frame(self.encode(value)), **kwargs)

    def drop(self):
        return self.partd.drop()

    @property
    def lock(self):
        return self.partd.lock

    def __exit__(self, *args):
        self.drop()
        self.partd.__exit__(*args)
