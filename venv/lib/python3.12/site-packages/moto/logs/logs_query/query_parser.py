from typing import List, Optional, Tuple

from moto.utilities.tokenizer import GenericTokenizer


class ParsedQuery:
    def __init__(self) -> None:
        self.limit: Optional[int] = None
        self.fields: List[str] = []
        self.sort: List[Tuple[str, str]] = []

    def sort_reversed(self) -> bool:
        # Descending is the default
        if self.sort:
            # sort_reversed is True if we want to sort in ascending order
            return self.sort[-1][-1] == "asc"
        return False


def parse_query(query: str) -> ParsedQuery:
    tokenizer = GenericTokenizer(query)
    state = "COMMAND"
    characters = ""
    parsed_query = ParsedQuery()

    for char in tokenizer:
        if char.isspace():
            if state == "SORT":
                parsed_query.sort.append((characters, "desc"))
                characters = ""
                state = "SORT_ORDER"
            if state == "COMMAND":
                if characters.lower() in ["fields", "limit", "sort"]:
                    state = characters.upper()
                else:
                    # Unknown/Unsupported command
                    pass
                characters = ""
            tokenizer.skip_white_space()
            continue

        if char == "|":
            if state == "FIELDS":
                parsed_query.fields.append(characters)
                characters = ""
            if state == "LIMIT":
                parsed_query.limit = int(characters)
                characters = ""
            if state == "SORT_ORDER":
                if characters != "":
                    parsed_query.sort[-1] = (parsed_query.sort[-1][0], characters)
                    characters = ""
            state = "COMMAND"
            tokenizer.skip_white_space()
            continue

        if char == ",":
            if state == "FIELDS":
                parsed_query.fields.append(characters)
                characters = ""
                continue

        characters += char

    if state == "FIELDS":
        parsed_query.fields.append(characters)
    if state == "LIMIT":
        parsed_query.limit = int(characters)
    if state == "SORT":
        parsed_query.sort.append((characters, "desc"))
    if state == "SORT_ORDER":
        parsed_query.sort[-1] = (parsed_query.sort[-1][0], characters)

    return parsed_query
