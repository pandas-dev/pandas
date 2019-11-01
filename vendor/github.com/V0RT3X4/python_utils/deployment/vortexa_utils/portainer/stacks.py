# @Author: richard
# @Date:   2018-12-04T18:04:55+00:00
# @Last modified by:   richard
# @Last modified time: 2018-12-04T18:34:09+00:00
from .api import RequestHelper


class Stacks(object):
    def __init__(self, api):
        self.api = api
        self.requests = RequestHelper(api, 'api/stacks')

    def get(self, stack_id):
        return self.requests.get(stack_id)

    def list(self):
        return self.get(None)

    def filter(self, **kwargs):
        def filter_kwargs(stack):
            return all(str(stack[k]) == str(v) for k, v in kwargs.items())

        return filter(filter_kwargs, self.list())

    def first(self, **kwargs):
        return next(self.filter(**kwargs))

    def get_file(self, stack_id):
        return self.requests.get(f'{stack_id}/file')

    def update(self, stack_id=None, endpointId=None, name=None,
               Env=None, StackFileContent=None, Prune=False):
        # get the stack by filtering on name or stack_id
        if name is not None:
            stack = self.first(Name=name)
            stack_id = stack['Id']
        elif stack_id is not None:
            stack = self.get(stack_id)

        endpointId = stack.get('EndpointId', endpointId)
        if endpointId is None:
            raise Exception("no entrypointID found or set")

        # update the old Env with the new Env
        old_Env = stack.get('Env')
        if old_Env is not None:
            update_keys = set(e['name'] for e in Env)
            old_Env = list(e for e in old_Env if e['name'] not in update_keys)
        Env += old_Env

        if StackFileContent is None:
            StackFileContent = self.get_file(stack_id)['StackFileContent']
        body = dict(StackFileContent=StackFileContent,
                    Env=Env,
                    Prune=Prune)

        return self.requests.put(
            stack_id,
            params=dict(endpointId=endpointId),
            json=body
        )
