"""CLI entrypoint for managing piplite wheels"""

from pathlib import Path

from jupyter_core.application import JupyterApp
from jupyterlite_core.app import DescribedMixin
from jupyterlite_core.trait_types import CPath

from ._version import __version__
from .addons.piplite import list_wheels


class PipliteIndex(DescribedMixin, JupyterApp):
    """index a directory of wheels for piplite into an all.json

    this file is suitable for including in a pre-built lab extension and will be
    found by adding to the extension's ``package.json``:

    .. code-block: json

        {
            "name": "my-extension",
            "jupyterlab": {
                "extension": true
            },
            "piplite": {
                "wheelDir": "./pypi"
            }
        }
    """

    version = __version__

    wheel_dir = CPath(Path.cwd(), help="a path of wheels")

    def parse_command_line(self, argv=None):
        super(PipliteIndex, self).parse_command_line(argv)

        if self.extra_args:
            self.wheel_dir = Path(self.extra_args[0])

    def start(self):
        if not self.wheel_dir.exists():
            raise ValueError(f"{self.wheel_dir} does not exist")
        if not list_wheels(self.wheel_dir):
            raise ValueError(f"no supported wheels found in {self.wheel_dir}")
        from .addons.piplite import write_wheel_index

        write_wheel_index(self.wheel_dir)


class PipliteApp(DescribedMixin, JupyterApp):
    """tools for working with piplite"""

    subcommands = {
        k: (v, v.__doc__.splitlines()[0].strip())
        for k, v in dict(
            index=PipliteIndex,
        ).items()
    }


main = launch_new_instance = PipliteApp.launch_instance

if __name__ == "__main__":  # pragma: no cover
    main()
