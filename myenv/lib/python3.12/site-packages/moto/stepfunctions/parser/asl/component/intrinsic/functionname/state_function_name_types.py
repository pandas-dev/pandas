from enum import Enum

from moto.stepfunctions.parser.asl.antlr.runtime.ASLIntrinsicLexer import (
    ASLIntrinsicLexer,
)


class StatesFunctionNameType(Enum):
    Format = ASLIntrinsicLexer.Format
    StringToJson = ASLIntrinsicLexer.StringToJson
    JsonToString = ASLIntrinsicLexer.JsonToString
    Array = ASLIntrinsicLexer.Array
    ArrayPartition = ASLIntrinsicLexer.ArrayPartition
    ArrayContains = ASLIntrinsicLexer.ArrayContains
    ArrayRange = ASLIntrinsicLexer.ArrayRange
    ArrayGetItem = ASLIntrinsicLexer.ArrayGetItem
    ArrayLength = ASLIntrinsicLexer.ArrayLength
    ArrayUnique = ASLIntrinsicLexer.ArrayUnique
    Base64Encode = ASLIntrinsicLexer.Base64Encode
    Base64Decode = ASLIntrinsicLexer.Base64Decode
    Hash = ASLIntrinsicLexer.Hash
    JsonMerge = ASLIntrinsicLexer.JsonMerge
    MathRandom = ASLIntrinsicLexer.MathRandom
    MathAdd = ASLIntrinsicLexer.MathAdd
    StringSplit = ASLIntrinsicLexer.StringSplit
    UUID = ASLIntrinsicLexer.UUID

    def name(self) -> str:  # noqa  # pylint: disable=function-redefined
        return ASLIntrinsicLexer.literalNames[self.value][1:-1]
