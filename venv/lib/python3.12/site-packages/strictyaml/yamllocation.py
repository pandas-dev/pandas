from strictyaml.ruamel.comments import CommentedSeq, CommentedMap
from strictyaml.exceptions import YAMLValidationError
from strictyaml.yamlpointer import YAMLPointer
from strictyaml import utils
from copy import deepcopy, copy
import sys

if sys.version_info[0] == 3:
    unicode = str


class YAMLChunk(object):
    """
    Represents a section of the document with references to the ruamel
    parsed document and the strictparsed document.

    Most operations done by validators on the document are done using this object.

    Before validation the strictparsed document will be identical to the
    ruamelparsed document. After it will contain CommentedMaps, CommentedSeqs
    and YAML objects.
    """

    def __init__(
        self,
        ruamelparsed,
        pointer=None,
        label=None,
        strictparsed=None,
        key_association=None,
    ):
        self._ruamelparsed = ruamelparsed
        self._strictparsed = (
            deepcopy(ruamelparsed) if strictparsed is None else strictparsed
        )
        self._pointer = pointer if pointer is not None else YAMLPointer()
        self._label = label

        # Associates strictparsed key names with ruamelparsed key names
        # E.g. "my-key-name" -> "My Key name"
        self._key_association = {} if key_association is None else key_association

    def expecting_but_found(self, expecting, found=None):
        raise YAMLValidationError(
            expecting,
            found if found is not None else "found {0}".format(self.found()),
            self,
        )

    def while_parsing_found(self, what, found=None):
        self.expecting_but_found("while parsing {0}".format(what), found=found)

    def process(self, new_item):
        strictparsed = self.pointer.parent().get(self._strictparsed, strictdoc=True)
        current_parsed = (
            strictparsed._value if hasattr(strictparsed, "_value") else strictparsed
        )

        def actual_key_from_string_key(string_key):
            if string_key in current_parsed.keys():
                return string_key
            else:
                for key in current_parsed.keys():
                    if hasattr(key, "_value"):
                        if key.text == string_key:
                            return key

        if self.pointer.is_index():
            current_parsed[self.pointer.last_index] = new_item
        elif self.pointer.is_val():
            current_parsed[
                actual_key_from_string_key(self.pointer.last_regularkey)
            ] = new_item
        elif self.pointer.is_key():
            key = actual_key_from_string_key(self.pointer.last_regularkey)
            existing_val = current_parsed[key]
            del current_parsed[key]
            current_parsed[new_item] = existing_val

    def is_sequence(self):
        return isinstance(self.contents, CommentedSeq)

    def is_mapping(self):
        return isinstance(self.contents, CommentedMap)

    def is_scalar(self):
        return not isinstance(self.contents, (CommentedMap, CommentedSeq))

    def found(self):
        if self.is_sequence():
            return "a sequence"
        elif self.is_mapping():
            return "a mapping"
        elif self.contents == "":
            return "a blank string"
        elif utils.is_integer(self.contents):
            return "an arbitrary integer"
        elif utils.is_decimal(self.contents):
            return "an arbitrary number"
        else:
            return "arbitrary text"

    def expect_sequence(self, expecting="when expecting a sequence"):
        if not self.is_sequence():
            self.expecting_but_found(expecting, "found {0}".format(self.found()))
        return [self.index(i) for i in range(len(self.contents))]

    def expect_mapping(self):
        if not self.is_mapping():
            self.expecting_but_found(
                "when expecting a mapping", "found {0}".format(self.found())
            )
        return [
            (
                self.key(regular_key, unicode(validated_key)),
                self.val(unicode(validated_key)),
            )
            for (regular_key, validated_key) in zip(
                self.contents.keys(), self.strictparsed().keys()
            )
        ]

    def expect_scalar(self, what):
        if not self.is_scalar():
            self.expecting_but_found(
                "when expecting {0}".format(what), "found {0}".format(self.found())
            )

    @property
    def label(self):
        return self._label

    @property
    def whole_document(self):
        return self._ruamelparsed

    @property
    def pointer(self):
        return self._pointer

    def fork(self, strictindex, new_value):
        """
        Return a chunk referring to the same location in a duplicated document.

        Used when modifying a YAML chunk so that the modification can be validated
        before changing it.
        """
        forked_chunk = YAMLChunk(
            deepcopy(self._ruamelparsed),
            pointer=self.pointer,
            label=self.label,
            key_association=copy(self._key_association),
        )
        if self.is_scalar():
            # Necessary for e.g. EmptyDict, which reports as a scalar.
            forked_chunk.pointer.set(forked_chunk, "_ruamelparsed", CommentedMap())
            forked_chunk.pointer.set(
                forked_chunk, "_strictparsed", CommentedMap(), strictdoc=True
            )
        forked_chunk.contents[self.ruamelindex(strictindex)] = new_value.as_marked_up()
        forked_chunk.strictparsed()[strictindex] = deepcopy(new_value.as_marked_up())
        return forked_chunk

    def add_key_association(self, unprocessed_key, processed_key):
        self._key_association[processed_key] = unprocessed_key

    @property
    def key_association(self):
        return self._key_association

    def make_child_of(self, chunk):
        """
        Link one YAML chunk to another.

        Used when inserting a chunk of YAML into another chunk.
        """
        if self.is_mapping():
            for key, value in self.contents.items():
                self.key(key, key).pointer.make_child_of(chunk.pointer)
                self.val(key).make_child_of(chunk)
        elif self.is_sequence():
            for index, item in enumerate(self.contents):
                self.index(index).make_child_of(chunk)
        else:
            self.pointer.make_child_of(chunk.pointer)

    def _select(self, pointer):
        """
        Get a YAMLChunk referenced by a pointer.
        """
        return YAMLChunk(
            self._ruamelparsed,
            pointer=pointer,
            label=self._label,
            strictparsed=self._strictparsed,
            key_association=copy(self._key_association),
        )

    def index(self, strictindex):
        """
        Return a chunk in a sequence referenced by index.
        """
        return self._select(self._pointer.index(self.ruamelindex(strictindex)))

    def ruamelindex(self, strictindex):
        """
        Get the ruamel equivalent of a strict parsed index.

        E.g. 0 -> 0, 1 -> 2, parsed-via-slugify -> Parsed via slugify
        """
        return (
            self.key_association.get(strictindex, strictindex)
            if self.is_mapping()
            else strictindex
        )

    def val(self, strictkey):
        """
        Return a chunk referencing a value in a mapping with the key 'key'.
        """
        ruamelkey = self.ruamelindex(strictkey)
        return self._select(self._pointer.val(ruamelkey, strictkey))

    def key(self, key, strictkey=None):
        """
        Return a chunk referencing a key in a mapping with the name 'key'.
        """
        return self._select(self._pointer.key(key, strictkey))

    def textslice(self, start, end):
        """
        Return a chunk referencing a slice of a scalar text value.
        """
        return self._select(self._pointer.textslice(start, end))

    def start_line(self):
        return self._pointer.start_line(self._ruamelparsed)

    def end_line(self):
        return self._pointer.end_line(self._ruamelparsed)

    def lines(self):
        return self._pointer.lines(self._ruamelparsed)

    def lines_before(self, how_many):
        return self._pointer.lines_before(self._ruamelparsed, how_many)

    def lines_after(self, how_many):
        return self._pointer.lines_after(self._ruamelparsed, how_many)

    @property
    def contents(self):
        return self._pointer.get(self._ruamelparsed)

    def strictparsed(self):
        return self._pointer.get(self._strictparsed, strictdoc=True)
