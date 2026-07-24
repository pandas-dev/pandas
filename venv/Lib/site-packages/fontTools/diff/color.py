from typing import Dict, Text

ansicolors: Dict[Text, Text] = {
    "BLACK": "\033[30m",
    "RED": "\033[31m",
    "GREEN": "\033[32m",
    "YELLOW": "\033[33m",
    "BLUE": "\033[34m",
    "MAGENTA": "\033[35m",
    "CYAN": "\033[36m",
    "WHITE": "\033[37m",
    "BOLD": "\033[1m",
    "RESET": "\033[0m",
}

green_start: Text = ansicolors["GREEN"]
red_start: Text = ansicolors["RED"]
cyan_start: Text = ansicolors["CYAN"]
reset: Text = ansicolors["RESET"]


def color_unified_diff_line(line: Text) -> Text:
    """Returns an ANSI escape code colored string with color based
    on the unified diff line type."""
    if line[0:2] == "+ ":
        return f"{green_start}{line}{reset}"
    elif line == "+\n":
        # some lines are formatted as hyphen only with no other characters
        # this indicates an added empty line
        return f"{green_start}{line}{reset}"
    elif line[0:2] == "- ":
        return f"{red_start}{line}{reset}"
    elif line == "-\n":
        # some lines are formatted as hyphen only with no other characters
        # this indicates a deleted empty line
        return f"{red_start}{line}{reset}"
    elif line[0:3] == "@@ ":
        return f"{cyan_start}{line}{reset}"
    elif line[0:4] == "--- ":
        return f"{red_start}{line}{reset}"
    elif line[0:4] == "+++ ":
        return f"{green_start}{line}{reset}"
    else:
        return line
