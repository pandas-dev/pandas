### IMPORTS
### ============================================================================
## Future

## Standard Library
import warnings

## Installed

## Application
from . import json
from . import utils

### CONSTANTS
### ============================================================================
ORJSON_AVAILABLE = utils.package_is_available("orjson")
MSGSPEC_AVAILABLE = utils.package_is_available("msgspec")
