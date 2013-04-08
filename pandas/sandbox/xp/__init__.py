
# monkey patched onto something, when the config options is enabled
def foo(self):
    print "Hi, I'm an experimental new method. grrrr!" % type(self)
