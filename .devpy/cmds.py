import click
from devpy.cmds.util import run


@click.command()
@click.option("--clean", is_flag=True, help="Clean build directory")
@click.option(
    "--debug",
    is_flag=True,
    help="Build C extensions with debugging symbols and without optimization",
)
def build(clean, debug):
    extra_args = []

    if debug:
        # Need to clean first, since there is a bug in
        # Cython, where it can't tell extensions need to be recompiled
        # when rebuilding in debug mode
        extra_args.append("--with-debugging-symbols")
        clean = True
    if clean:
        print("Running clean.")
        run(["python", "setup.py", "clean"])

    # TODO: Rely on pip for backend agnostic build
    build_cmd = ["python", "setup.py", "develop"]

    build_cmd += extra_args

    run(build_cmd)


@click.command()
@click.argument("pytest_args", nargs=-1)
def test(pytest_args):
    run(["python", "-m", "pytest"] + list(pytest_args))
