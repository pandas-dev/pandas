def setup():
    # module-level
    pass


def setup_cache():
    # module-level
    pass


def track_test():
    # module-level 難
    return 0


def track_pretty_source_test():
    return 0


track_pretty_source_test.pretty_source = '''
    int track_pretty_source_test() {
        return 0;
    }'''


class MyClass:
    def setup(self):
        # class-level
        pass

    def setup_cache(self):
        # class-level
        pass

    def track_test(self):
        # class-level 難
        return 0
