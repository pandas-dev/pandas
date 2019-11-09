from argparse import ArgumentParser
from dataclasses import dataclass, field
from vortexa_utils.versioning.versioner import Versioner


@dataclass
class VersionCLI(object):
    versioner: Versioner
    parser: ArgumentParser = field(default=None, init=False)

    def __post_init__(self):
        self.parser = ArgumentParser(
            description='Package Version Tool.'
        )
        self.specs = self.versioner.VERSION_SPEC.split(
            self.versioner.VERSION_SEP
        )
        for spec in self.specs:
            self.parser.add_argument(
                f'--bump-{spec.lower()}',
                f'-{spec[0]}',
                action='store_true'
            )

    def parse_args(self):
        args = self.parser.parse_args()
        spec_flags = list(
            getattr(args, f'bump_{spec.lower()}')
            for spec in self.specs
        )
        if any(spec_flags):
            print(f"Current Version: {self.versioner}")
            if sum(spec_flags) > 1:
                print("You can only bump one spec at a time")
                self.parser.print_help()
            else:
                self.versioner.update_version(spec_flags)
                print(f"New Version {self.versioner}")
        else:
            print(f"{self.versioner}")


if __name__ == "__main__":
    version = Versioner()
    cli = VersionCLI(version)
    cli.parse_args()
