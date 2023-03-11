import argparse
import string


def generate_env(template_file, outfile, python_version, pyarrow_version, bad_deps):
    str_env_file = ""
    with open(template_file) as tmpl_file:
        for line in tmpl_file:
            append = True
            for bad_dep in bad_deps:
                if bad_dep in line:
                    append = False
                    break
            if pyarrow_version not in [None, "latest"] and "pyarrow" in line:
                line = line.split("pyarrow")[
                    0
                ]  # Keep the stuff before we specify the version
                line += f"pyarrow={pyarrow_version}\n"
            if append:
                str_env_file += line

    tmpl = string.Template(str_env_file)
    str_env_file = tmpl.substitute(python_version=python_version)
    with open(outfile, "w") as out_file:
        out_file.write(str_env_file)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "infile", type=str, help="Path to the template environment file"
    )
    parser.add_argument("outfile", type=str, help="Path to write the file to")
    parser.add_argument(
        "--python_version", help="Python version to generate an environment file for"
    )
    parser.add_argument(
        "--pyarrow_version", help="Pyarrow version to use in the environment"
    )
    parser.add_argument(
        "--disable_deps",
        default=[],
        nargs="+",
        help="Space separated deps to remove from the environment file",
    )
    args = parser.parse_args()

    generate_env(
        args.infile,
        args.outfile,
        args.python_version,
        args.pyarrow_version,
        args.disable_deps,
    )


main()
