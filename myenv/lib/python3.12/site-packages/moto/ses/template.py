from typing import Any, Dict, Optional, Type

from moto.utilities.tokenizer import GenericTokenizer

from .exceptions import MissingRenderingAttributeException


class BlockProcessor:
    def __init__(
        self, template: str, template_data: Dict[str, Any], tokenizer: GenericTokenizer
    ):
        self.template = template
        self.template_data = template_data
        self.tokenizer = tokenizer

    def parse(self) -> str:
        # Added to make MyPy happy
        # Not all implementations have this method
        # It's up to the caller to know whether to call this method
        raise NotImplementedError


class EachBlockProcessor(BlockProcessor):
    def __init__(
        self, template: str, template_data: Dict[str, Any], tokenizer: GenericTokenizer
    ):
        self.template = template
        self.tokenizer = tokenizer

        self.tokenizer.skip_characters("#each")
        self.tokenizer.skip_white_space()
        var_name = self.tokenizer.read_until("}}").strip()
        self.tokenizer.skip_characters("}}")
        self.template_data = template_data.get(var_name, [])

    def parse(self) -> str:
        parsed = ""
        current_pos = self.tokenizer.token_pos

        for template_data in self.template_data:
            self.tokenizer.token_pos = current_pos
            for char in self.tokenizer:
                if char == "{" and self.tokenizer.peek() == "{":
                    self.tokenizer.skip_characters("{")
                    self.tokenizer.skip_white_space()

                    _processor = get_processor(self.tokenizer)(
                        self.template,
                        template_data,  # type: ignore
                        self.tokenizer,
                    )
                    # If we've reached the end, we should stop processing
                    # Our parent will continue with whatever comes after {{/each}}
                    if type(_processor) == EachEndBlockProcessor:
                        break
                    # If we've encountered another processor, they can continue
                    parsed += _processor.parse()

                    continue

                parsed += char

        return parsed


class EachEndBlockProcessor(BlockProcessor):
    def __init__(
        self, template: str, template_data: Dict[str, Any], tokenizer: GenericTokenizer
    ):
        super().__init__(template, template_data, tokenizer)

        self.tokenizer.skip_characters("/each")
        self.tokenizer.skip_white_space()
        self.tokenizer.skip_characters("}}")


class IfBlockProcessor(BlockProcessor):
    def __init__(
        self, template: str, template_data: Dict[str, Any], tokenizer: GenericTokenizer
    ):
        super().__init__(template, template_data, tokenizer)

        self.tokenizer.skip_characters("#if")
        self.tokenizer.skip_white_space()
        condition = self.tokenizer.read_until("}}").strip()
        self.tokenizer.skip_characters("}}")
        self.parse_contents = template_data.get(condition)

    def parse(self) -> str:
        parsed = ""

        for char in self.tokenizer:
            if char == "{" and self.tokenizer.peek() == "{":
                self.tokenizer.skip_characters("{")
                self.tokenizer.skip_white_space()

                _processor = get_processor(self.tokenizer)(
                    self.template, self.template_data, self.tokenizer
                )
                if type(_processor) == IfEndBlockProcessor:
                    break
                elif type(_processor) == ElseBlockProcessor:
                    self.parse_contents = not self.parse_contents
                    continue
                if self.parse_contents:
                    parsed += _processor.parse()

                continue

            if self.parse_contents:
                parsed += char

        return parsed


class IfEndBlockProcessor(BlockProcessor):
    def __init__(
        self, template: str, template_data: Dict[str, Any], tokenizer: GenericTokenizer
    ):
        super().__init__(template, template_data, tokenizer)

        self.tokenizer.skip_characters("/if")
        self.tokenizer.skip_white_space()
        self.tokenizer.skip_characters("}}")


class ElseBlockProcessor(BlockProcessor):
    def __init__(
        self, template: str, template_data: Dict[str, Any], tokenizer: GenericTokenizer
    ):
        super().__init__(template, template_data, tokenizer)

        self.tokenizer.skip_characters("else")
        self.tokenizer.skip_white_space()
        self.tokenizer.skip_characters("}}")


class VarBlockProcessor(BlockProcessor):
    def parse(self) -> str:
        var_name = self.tokenizer.read_until("}}").strip()
        if self.template_data.get(var_name) is None:
            raise MissingRenderingAttributeException(var_name)
        data: str = self.template_data.get(var_name)  # type: ignore
        self.tokenizer.skip_white_space()
        self.tokenizer.skip_characters("}}")
        return data


def get_processor(tokenizer: GenericTokenizer) -> Type[BlockProcessor]:
    if tokenizer.peek(5) == "#each":
        return EachBlockProcessor
    if tokenizer.peek(5) == "/each":
        return EachEndBlockProcessor
    if tokenizer.peek(3) == "#if":
        return IfBlockProcessor
    if tokenizer.peek(3) == "/if":
        return IfEndBlockProcessor
    if tokenizer.peek(4) == "else":
        return ElseBlockProcessor
    return VarBlockProcessor


def parse_template(
    template: str,
    template_data: Dict[str, Any],
    tokenizer: Optional[GenericTokenizer] = None,
) -> str:
    tokenizer = tokenizer or GenericTokenizer(template)

    parsed = ""
    for char in tokenizer:
        if char == "{" and tokenizer.peek() == "{":
            # Two braces next to each other indicate a variable/language construct such as for-each
            # We have different processors handling different constructs
            tokenizer.skip_characters("{")
            tokenizer.skip_white_space()

            _processor = get_processor(tokenizer)(template, template_data, tokenizer)
            parsed += _processor.parse()
            continue

        parsed += char
    return parsed
