from typing import Any, Dict, List, Optional


class Policy(object):
    def __init__(self, policy_name: str, policy_type_name: str):
        self.policy_name = policy_name
        self.policy_type_name = policy_type_name


class AppCookieStickinessPolicy(Policy):
    def __init__(self, policy_name: str, cookie_name: str):
        super().__init__(policy_name, policy_type_name="AppCookieStickinessPolicy")
        self.cookie_name = cookie_name


class LbCookieStickinessPolicy(Policy):
    def __init__(self, policy_name: str, cookie_expiration_period: Optional[int]):
        super().__init__(policy_name, policy_type_name="LbCookieStickinessPolicy")
        self.cookie_expiration_period = cookie_expiration_period


class OtherPolicy(Policy):
    def __init__(
        self,
        policy_name: str,
        policy_type_name: str,
        policy_attrs: List[Dict[str, Any]],
    ):
        super().__init__(policy_name, policy_type_name=policy_type_name)
        self.attributes = policy_attrs or []
