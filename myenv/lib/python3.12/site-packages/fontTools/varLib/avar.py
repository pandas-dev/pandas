from fontTools.varLib import _add_avar, load_designspace
from fontTools.misc.cliTools import makeOutputFileName
import logging

log = logging.getLogger("fontTools.varLib.avar")


def main(args=None):
    """Add `avar` table from designspace file to variable font."""

    if args is None:
        import sys

        args = sys.argv[1:]

    from fontTools import configLogger
    from fontTools.ttLib import TTFont
    from fontTools.designspaceLib import DesignSpaceDocument
    import argparse

    parser = argparse.ArgumentParser(
        "fonttools varLib.avar",
        description="Add `avar` table from designspace file to variable font.",
    )
    parser.add_argument("font", metavar="varfont.ttf", help="Variable-font file.")
    parser.add_argument(
        "designspace", metavar="family.designspace", help="Designspace file."
    )
    parser.add_argument(
        "-o",
        "--output-file",
        type=str,
        help="Output font file name.",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Run more verbosely."
    )

    options = parser.parse_args(args)

    configLogger(level=("INFO" if options.verbose else "WARNING"))

    font = TTFont(options.font)
    if not "fvar" in font:
        log.error("Not a variable font.")
        return 1

    axisTags = [a.axisTag for a in font["fvar"].axes]

    ds = load_designspace(options.designspace)

    if "avar" in font:
        log.warning("avar table already present, overwriting.")
        del font["avar"]

    _add_avar(font, ds.axes, ds.axisMappings, axisTags)

    if options.output_file is None:
        outfile = makeOutputFileName(options.font, overWrite=True, suffix=".avar")
    else:
        outfile = options.output_file
    if outfile:
        log.info("Saving %s", outfile)
        font.save(outfile)


if __name__ == "__main__":
    import sys

    sys.exit(main())
