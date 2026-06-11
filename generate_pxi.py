import argparse
import os

from Cython import Tempita


def process_tempita(pxifile, outfile) -> None:
    with open(pxifile, encoding="utf-8") as f:
        tmpl = f.read()
    pyxcontent = Tempita.sub(tmpl)

    with open(outfile, "w", encoding="utf-8") as f:
        f.write(pyxcontent)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", type=str, help="Path to the input file")
    parser.add_argument("-o", "--outdir", type=str, help="Path to the output directory")
    args = parser.parse_args()

    if not args.infile.endswith(".in"):
        raise ValueError(f"Unexpected extension: {args.infile}")

    outdir_abs = os.path.join(os.getcwd(), args.outdir)
    outfile = os.path.join(
        outdir_abs, os.path.splitext(os.path.split(args.infile)[1])[0]
    )

    process_tempita(args.infile, outfile)


main()
