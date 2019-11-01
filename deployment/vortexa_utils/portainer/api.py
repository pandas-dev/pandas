# @Author: richard
# @Date:   2018-12-04T18:05:38+00:00
# @Last modified by:   richard
# @Last modified time: 2018-12-04T18:14:40+00:00
import os
import requests
from functools import wraps
from urllib.parse import urlparse
from .stacks import Stacks
from . import notNone


class RequestHelper(object):
    def __init__(self, api, base_url='api'):
        self.api = api
        self.base_url = base_url

    def wrapper(self, func):
        @wraps(func)
        def caller(url=None, *args, **kwargs):
            parts = filter(notNone, (self.api.host, self.base_url, url))
            parts = map(str, parts)
            headers = kwargs.get("headers", {})
            headers.update(self.api.get_header())
            kwargs["headers"] = headers
            return func(os.path.join(*parts),
                        *args, **kwargs).json()
        return caller

    def __getattr__(self, name, *args, **kwargs):
        method = getattr(requests, name, *args, **kwargs)
        return self.wrapper(method)


class PortainerAPI(object):
    def __init__(self, host, user=None, pw=None):
        self.host = urlparse(host, scheme='http').geturl()
        self.user = user
        self.pw = pw
        if any(ting is not None for ting in (host, user, pw)):
            self.get_jwt()
        self.requests = RequestHelper(self)
        self.stacks = Stacks(self)

    def get_jwt(self):
        """
        http POST :9000/api/auth Username="admin" Password="adminpassword"
        """
        url = f'{self.host}/api/auth'
        resp = requests.post(url, json=dict(Username=self.user,
                                            Password=self.pw))
        self.token = resp.json().get('jwt')
        return self.token

    def get_header(self):
        return {"Authorization": f"Bearer {self.token}"}
