from .distribution import Distribution

class Index(dict):

    def __setitem__(self, key, value):
        if not isinstance(value, Distribution):
            raise ValueError('Not a distribution: %r.' % value)
        if key != '%s-%s' % (value.name, value.version):
            raise ValueError('Key must match <name>-<version>.')
        super(Index, self).__setitem__(key, value)

    def add(self, distribution):
        key = '%s-%s' % (distribution.name, distribution.version)
        self[key] = distribution

