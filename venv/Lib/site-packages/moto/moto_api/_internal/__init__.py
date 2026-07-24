from .models import moto_api_backend
from .moto_random import MotoRandom
from .recorder.models import Recorder  # noqa
from .state_manager import StateManager  # noqa

moto_api_backends = {"global": moto_api_backend}

"""
Random-object used throughout Moto. Can be seeded to make identifiers deterministic.
"""
mock_random = MotoRandom()
