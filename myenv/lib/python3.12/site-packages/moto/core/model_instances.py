from collections import defaultdict
from typing import Dict

"""
Storage of all instances that extend BaseModel
This allows us to easily expose all internal state using the MotoServer dashboard
"""
model_data: Dict[str, Dict[str, object]] = defaultdict(dict)


def reset_model_data() -> None:
    # Remove all references to the models stored
    for models in model_data.values():
        for model in models.values():
            model.instances = []  # type: ignore[attr-defined]
