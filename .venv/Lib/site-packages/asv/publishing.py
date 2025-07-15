# Licensed under a 3-clause BSD style license - see LICENSE.rst

import math


class OutputPublisher:
    """
    A base class for pages displaying output in the JS application
    """
    name = None
    button_label = None
    description = None
    order = math.inf

    @classmethod
    def publish(cls, conf, repo, benchmarks, graphs, revisions):
        pass
