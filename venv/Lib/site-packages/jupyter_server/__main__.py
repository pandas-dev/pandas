"""The main entry point for Jupyter Server."""

if __name__ == "__main__":
    from jupyter_server import serverapp as app

    app.launch_new_instance()
