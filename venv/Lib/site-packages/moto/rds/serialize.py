from __future__ import annotations

from collections import defaultdict
from functools import lru_cache
from typing import Any, Dict, List, Optional

from botocore import xform_name

from moto.core.utils import get_service_model


class XFormedAttributeAccessMixin:
    """Mixin allowing access to "xformed" attributes:

    obj.DBInstanceIdentifier will retrieve the value of obj.db_instance_identifier

    """

    BOTOCORE_MODEL: Optional[str] = None

    _model_attribute_aliases: Dict[str, List[str]] = {}
    _xform_cache: Dict[str, str] = {}

    def __getattr__(self, name: str) -> Any:
        if name in self.model_attributes:
            return self.get_modeled_attribute(name)
        raise AttributeError(f"Attribute '{name}' not found!")

    def get_modeled_attribute(self, attr_name: str) -> Any:
        for attr_alias in self.model_attribute_aliases[attr_name]:
            try:
                return super().__getattribute__(attr_alias)
            except AttributeError:
                pass
        else:
            raise AttributeError

    @property
    def model_attributes(self) -> List[str]:
        return list(self.model_attribute_aliases.keys())

    @property
    def model_attribute_aliases(self) -> Dict[str, List[str]]:
        if not self._model_attribute_aliases:
            self._model_attribute_aliases = self.get_model_attributes_info()
        return self._model_attribute_aliases

    @classmethod
    @lru_cache()
    def get_model_attributes_info(cls) -> Dict[str, List[str]]:
        service_name = cls.__module__.split(".")[1]
        model_name = cls.BOTOCORE_MODEL or cls.__name__
        service_model = get_service_model(service_name)
        model_shape = service_model.shape_for(model_name)
        valid_attributes: Dict[str, List[str]] = defaultdict(list)
        for member_name, member_shape in model_shape.members.items():  # type: ignore[attr-defined]
            aliases = valid_attributes[member_name]
            if member_shape.type_name == "list":
                if member_name != member_shape.name:
                    xformed_name = cls._xform_name(member_shape.name)
                    aliases.append(xformed_name)
            xformed_member_name = cls._xform_name(member_name)
            aliases.append(xformed_member_name)
            if member_name.startswith(model_name):
                short_name = member_name[len(model_name) :]
                xformed_short_name = cls._xform_name(short_name)
                aliases.append(xformed_short_name)
        return valid_attributes

    @classmethod
    def _xform_name(cls, name: str) -> str:
        return xform_name(name, _xform_cache=cls._xform_cache)
