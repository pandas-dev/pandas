#!/user/bin/env python3
# @Author: richard
# @Date:   2018-12-04T18:10:07+00:00
# @Last modified by:   richard
# @Last modified time: 2018-12-04T18:31:45+00:00
import argparse
from pprint import pprint
from .import notNone
from .api import PortainerAPI


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Portainer API CLI')
    parser.add_argument('--host', '-H', type=str,
                        help='Host name of Portainer API',
                        default='https://lawgiver.vortexa.com:9000')
    parser.add_argument('--user', '-u', type=str,
                        help='User name',
                        default='kraftwork_updater')
    parser.add_argument('--pass', '-p', type=str, dest='password',
                        help='Password name')

    parser.add_argument('--name', '-n', type=str,
                        help='Stack name to filter')

    parser.add_argument('--env', '-e', nargs=2, action='append',
                        help='key value pairs of confic to update')

    parser.add_argument('--filter', '-f', nargs=2, action='append',
                        help='key value pairs of confic to update')

    def add_cmd(flag):
        def command(func):
            parser.add_argument(
                flag,
                action='store_const',
                const=func,
                dest='cmd'
            )
            return func

    def get_filter():
        Filter = {}
        if args.filter is not None:
            Filter.update(args.filter)
        if args.name is not None:
            Filter.update(Name=args.name)
        return Filter

    @add_cmd('--list')
    def list_stacks():
        if any(map(notNone, ((args.name, args.filter)))):
            Filter = get_filter()
            return list(api.stacks.filter(**Filter))
        else:
            return api.stacks.list()

    @add_cmd('--update')
    def update_stacks():
        env = [dict(name=k, value=v) for k, v in args.env]
        return api.stacks.update(name=args.name, Env=env)

    args = parser.parse_args()

    api = PortainerAPI(host=args.host,
                       user=args.user,
                       pw=args.password)

    pprint(args.cmd())

#    api.stacks.list()
#    api.stacks.update(
#        1, 1,
#        Env=[{
#            "name": "KFAFTWERK_BUILD_NUM",
#            "value": '376'
#        }]
#    )
#
#
#    content = Path('docker/scripts/docker-compose.yml').read_text()
#
#    api.requests.post('stacks?type=1&method=string&endpointId=1',
#        json=dict(
#            Name="myStack",
#            StackFileContent=content,
#            Env=[dict(name="Hello",value="world")],
#            SwarmID='729a4f2h5kj2sd42x34pl3uu1'
#        )
#    )
