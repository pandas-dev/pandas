"""Extension config."""

from jupyter_server.services.config.manager import ConfigManager

DEFAULT_SECTION_NAME = "jupyter_server_config"


class ExtensionConfigManager(ConfigManager):
    """A manager class to interface with Jupyter Server Extension config
    found in a `config.d` folder. It is assumed that all configuration
    files in this directory are JSON files.
    """

    def get_jpserver_extensions(self, section_name=DEFAULT_SECTION_NAME):
        """Return the jpserver_extensions field from all
        config files found."""
        data = self.get(section_name)
        return data.get("ServerApp", {}).get("jpserver_extensions", {})

    def enabled(self, name, section_name=DEFAULT_SECTION_NAME, include_root=True):
        """Is the extension enabled?"""
        extensions = self.get_jpserver_extensions(section_name)
        try:
            return extensions[name]
        except KeyError:
            return False

    def enable(self, name):
        """Enable an extension by name."""
        data = {"ServerApp": {"jpserver_extensions": {name: True}}}
        self.update(name, data)

    def disable(self, name):
        """Disable an extension by name."""
        data = {"ServerApp": {"jpserver_extensions": {name: False}}}
        self.update(name, data)
