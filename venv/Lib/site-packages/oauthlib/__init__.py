"""
    oauthlib
    ~~~~~~~~

    A generic, spec-compliant, thorough implementation of the OAuth
    request-signing logic.

    :copyright: (c) 2019 by The OAuthlib Community
    :license: BSD, see LICENSE for details.
"""
import logging
from logging import NullHandler

__author__ = 'The OAuthlib Community'
__version__ = '3.2.2'

logging.getLogger('oauthlib').addHandler(NullHandler())

_DEBUG = False

def set_debug(debug_val):
	"""Set value of debug flag
	
    :param debug_val: Value to set. Must be a bool value.
	"""
	global _DEBUG
	_DEBUG = debug_val

def get_debug():
	"""Get debug mode value. 
	
	:return: `True` if debug mode is on, `False` otherwise
	"""
	return _DEBUG	
