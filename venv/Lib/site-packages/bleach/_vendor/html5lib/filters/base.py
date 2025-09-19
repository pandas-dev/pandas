class Filter:
    def __init__(self, source):
        self.source = source

    def __iter__(self):
        return iter(self.source)

    def __getattr__(self, name):
        return getattr(self.source, name)
