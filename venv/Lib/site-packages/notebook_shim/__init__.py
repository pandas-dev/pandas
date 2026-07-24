def _jupyter_server_extension_points():
    return [
        {
            'module': 'notebook_shim.nbserver',
        }
    ]
