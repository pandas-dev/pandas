# For usage of lark with PyInstaller. See https://pyinstaller-sample-hook.readthedocs.io/en/latest/index.html

import os

def get_hook_dirs():
    return [os.path.dirname(__file__)]
