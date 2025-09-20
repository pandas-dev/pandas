from collections import defaultdict


class DuplicatedNameError(NameError):
    pass


class NameScope(object):
    def __init__(self):
        self._useset = set([''])
        self._basenamemap = defaultdict(int)

    def is_used(self, name):
        return name in self._useset

    def register(self, name, deduplicate=False):
        if deduplicate:
            name = self.deduplicate(name)
        elif self.is_used(name):
            raise DuplicatedNameError(name)
        self._useset.add(name)
        return name

    def deduplicate(self, name):
        basename = name
        while self.is_used(name):
            ident = self._basenamemap[basename] + 1
            self._basenamemap[basename] = ident
            name = "{0}.{1}".format(basename, ident)
        return name

    def get_child(self):
        return type(self)(parent=self)


class _StrCaching(object):

    def _clear_string_cache(self):
        try:
            del self.__cached_str
        except AttributeError:
            pass

    def __str__(self):
        try:
            return self.__cached_str
        except AttributeError:
            s = self.__cached_str = self._to_string()
            return s


class _StringReferenceCaching(object):

    def get_reference(self):
        try:
            return self.__cached_refstr
        except AttributeError:
            s = self.__cached_refstr = self._get_reference()
            return s


class _HasMetadata(object):

    def set_metadata(self, name, node):
        """
        Attach unnamed metadata *node* to the metadata slot *name* of this
        value.
        """
        self.metadata[name] = node

    def _stringify_metadata(self, leading_comma=False):
        if self.metadata:
            buf = []
            if leading_comma:
                buf.append("")
            buf += ["!{0} {1}".format(k, v.get_reference())
                    for k, v in self.metadata.items()]
            return ', '.join(buf)
        else:
            return ''
