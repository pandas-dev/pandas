"""pyzmq log watcher.

Easily view log messages published by the PUBHandler in zmq.log.handlers

Designed to be run as an executable module - try this to see options:
    python -m zmq.log -h

Subscribes to the '' (empty string) topic by default which means it will work
out-of-the-box with a PUBHandler object instantiated with default settings.
If you change the root topic with PUBHandler.setRootTopic() you must pass
the value to this script with the --topic argument.

Note that the default formats for the PUBHandler object selectively include
the log level in the message. This creates redundancy in this script as it
always prints the topic of the message, which includes the log level.
Consider overriding the default formats with PUBHandler.setFormat() to
avoid this issue.

"""

# encoding: utf-8

# Copyright (C) PyZMQ Developers
# Distributed under the terms of the Modified BSD License.

import argparse
from datetime import datetime
from typing import Dict

import zmq

parser = argparse.ArgumentParser('ZMQ Log Watcher')
parser.add_argument('zmq_pub_url', type=str, help='URL to a ZMQ publisher socket.')
parser.add_argument(
    '-t',
    '--topic',
    type=str,
    default='',
    help='Only receive messages that start with this topic.',
)
parser.add_argument(
    '--timestamp', action='store_true', help='Append local time to the log messages.'
)
parser.add_argument(
    '--separator',
    type=str,
    default=' | ',
    help='String to print between topic and message.',
)
parser.add_argument(
    '--dateformat',
    type=str,
    default='%Y-%d-%m %H:%M',
    help='Set alternative date format for use with --timestamp.',
)
parser.add_argument(
    '--align',
    action='store_true',
    default=False,
    help='Try to align messages by the width of their topics.',
)
parser.add_argument(
    '--color',
    action='store_true',
    default=False,
    help='Color the output based on the error level. Requires the colorama module.',
)
args = parser.parse_args()


if args.color:
    import colorama

    colorama.init()
    colors = {
        'DEBUG': colorama.Fore.LIGHTCYAN_EX,
        'INFO': colorama.Fore.LIGHTWHITE_EX,
        'WARNING': colorama.Fore.YELLOW,
        'ERROR': colorama.Fore.LIGHTRED_EX,
        'CRITICAL': colorama.Fore.LIGHTRED_EX,
        '__RESET__': colorama.Fore.RESET,
    }
else:
    colors = {}


ctx = zmq.Context()
sub = ctx.socket(zmq.SUB)
sub.subscribe(args.topic.encode("utf8"))
sub.connect(args.zmq_pub_url)

topic_widths: Dict[int, int] = {}

while True:
    try:
        if sub.poll(10, zmq.POLLIN):
            topic, msg = sub.recv_multipart()
            topics = topic.decode('utf8').strip().split('.')

            if args.align:
                topics.extend(' ' for extra in range(len(topics), len(topic_widths)))
                aligned_parts = []
                for key, part in enumerate(topics):
                    topic_widths[key] = max(len(part), topic_widths.get(key, 0))
                    fmt = ''.join(('{:<', str(topic_widths[key]), '}'))
                    aligned_parts.append(fmt.format(part))

            if len(topics) == 1:
                level = topics[0]
            else:
                level = topics[1]

            fields = {
                'msg': msg.decode('utf8').strip(),
                'ts': (
                    datetime.now().strftime(args.dateformat) + ' '
                    if args.timestamp
                    else ''
                ),
                'aligned': (
                    '.'.join(aligned_parts)
                    if args.align
                    else topic.decode('utf8').strip()
                ),
                'color': colors.get(level, ''),
                'color_rst': colors.get('__RESET__', ''),
                'sep': args.separator,
            }
            print('{ts}{color}{aligned}{sep}{msg}{color_rst}'.format(**fields))
    except KeyboardInterrupt:
        break

sub.disconnect(args.zmq_pub_url)
if args.color:
    print(colorama.Fore.RESET)
