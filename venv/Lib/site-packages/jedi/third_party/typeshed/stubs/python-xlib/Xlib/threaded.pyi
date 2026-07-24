# This isn't just a re-export from from Xlib.support import lock
# Importing from this module will cause the lock.allocate_lock function to
# return a basic Python lock, instead of the default dummy lock
from Xlib.support import lock as lock
