"""Message constants for generating error messages during type checking.

Literal messages should be defined as constants in this module so they won't get out of sync
if used in more than one place, and so that they can be easily introspected. These messages are
ultimately consumed by messages.MessageBuilder.fail(). For more non-trivial message generation,
add a method to MessageBuilder and call this instead.
"""

from typing_extensions import Final

# Invalid types
INVALID_TYPE_RAW_ENUM_VALUE = "Invalid type: try using Literal[{}.{}] instead?"  # type: Final

# Type checker error message constants
NO_RETURN_VALUE_EXPECTED = 'No return value expected'  # type: Final
MISSING_RETURN_STATEMENT = 'Missing return statement'  # type: Final
INVALID_IMPLICIT_RETURN = 'Implicit return in function which does not return'  # type: Final
INCOMPATIBLE_RETURN_VALUE_TYPE = 'Incompatible return value type'  # type: Final
RETURN_VALUE_EXPECTED = 'Return value expected'  # type: Final
NO_RETURN_EXPECTED = 'Return statement in function which does not return'  # type: Final
INVALID_EXCEPTION = 'Exception must be derived from BaseException'  # type: Final
INVALID_EXCEPTION_TYPE = 'Exception type must be derived from BaseException'  # type: Final
RETURN_IN_ASYNC_GENERATOR = '"return" with value in async generator is not allowed'  # type: Final
INVALID_RETURN_TYPE_FOR_GENERATOR = \
    'The return type of a generator function should be "Generator"' \
    ' or one of its supertypes'  # type: Final
INVALID_RETURN_TYPE_FOR_ASYNC_GENERATOR = \
    'The return type of an async generator function should be "AsyncGenerator" or one of its ' \
    'supertypes'  # type: Final
INVALID_GENERATOR_RETURN_ITEM_TYPE = \
    'The return type of a generator function must be None in' \
    ' its third type parameter in Python 2'  # type: Final
YIELD_VALUE_EXPECTED = 'Yield value expected'  # type: Final
INCOMPATIBLE_TYPES = 'Incompatible types'  # type: Final
INCOMPATIBLE_TYPES_IN_ASSIGNMENT = 'Incompatible types in assignment'  # type: Final
INCOMPATIBLE_REDEFINITION = 'Incompatible redefinition'  # type: Final
INCOMPATIBLE_TYPES_IN_AWAIT = 'Incompatible types in "await"'  # type: Final
INCOMPATIBLE_TYPES_IN_ASYNC_WITH_AENTER = \
    'Incompatible types in "async with" for "__aenter__"'  # type: Final
INCOMPATIBLE_TYPES_IN_ASYNC_WITH_AEXIT = \
    'Incompatible types in "async with" for "__aexit__"'  # type: Final
INCOMPATIBLE_TYPES_IN_ASYNC_FOR = 'Incompatible types in "async for"'  # type: Final

INCOMPATIBLE_TYPES_IN_YIELD = 'Incompatible types in "yield"'  # type: Final
INCOMPATIBLE_TYPES_IN_YIELD_FROM = 'Incompatible types in "yield from"'  # type: Final
INCOMPATIBLE_TYPES_IN_STR_INTERPOLATION = \
    'Incompatible types in string interpolation'  # type: Final
MUST_HAVE_NONE_RETURN_TYPE = 'The return type of "{}" must be None'  # type: Final
INVALID_TUPLE_INDEX_TYPE = 'Invalid tuple index type'  # type: Final
TUPLE_INDEX_OUT_OF_RANGE = 'Tuple index out of range'  # type: Final
INVALID_SLICE_INDEX = 'Slice index must be an integer or None'  # type: Final
CANNOT_INFER_LAMBDA_TYPE = 'Cannot infer type of lambda'  # type: Final
CANNOT_ACCESS_INIT = 'Cannot access "__init__" directly'  # type: Final
NON_INSTANCE_NEW_TYPE = '"__new__" must return a class instance (got {})'  # type: Final
INVALID_NEW_TYPE = 'Incompatible return type for "__new__"'  # type: Final
BAD_CONSTRUCTOR_TYPE = 'Unsupported decorated constructor type'  # type: Final
CANNOT_ASSIGN_TO_METHOD = 'Cannot assign to a method'  # type: Final
CANNOT_ASSIGN_TO_TYPE = 'Cannot assign to a type'  # type: Final
INCONSISTENT_ABSTRACT_OVERLOAD = \
    'Overloaded method has both abstract and non-abstract variants'  # type: Final
MULTIPLE_OVERLOADS_REQUIRED = 'Single overload definition, multiple required'  # type: Final
READ_ONLY_PROPERTY_OVERRIDES_READ_WRITE = \
    'Read-only property cannot override read-write property'  # type: Final
FORMAT_REQUIRES_MAPPING = 'Format requires a mapping'  # type: Final
RETURN_TYPE_CANNOT_BE_CONTRAVARIANT = \
    "Cannot use a contravariant type variable as return type"  # type: Final
FUNCTION_PARAMETER_CANNOT_BE_COVARIANT = \
    "Cannot use a covariant type variable as a parameter"  # type: Final
INCOMPATIBLE_IMPORT_OF = "Incompatible import of"  # type: Final
FUNCTION_TYPE_EXPECTED = "Function is missing a type annotation"  # type: Final
ONLY_CLASS_APPLICATION = "Type application is only supported for generic classes"  # type: Final
RETURN_TYPE_EXPECTED = "Function is missing a return type annotation"  # type: Final
ARGUMENT_TYPE_EXPECTED = \
    "Function is missing a type annotation for one or more arguments"  # type: Final
KEYWORD_ARGUMENT_REQUIRES_STR_KEY_TYPE = \
    'Keyword argument only valid with "str" key type in call to "dict"'  # type: Final
ALL_MUST_BE_SEQ_STR = 'Type of __all__ must be {}, not {}'  # type: Final
INVALID_TYPEDDICT_ARGS = \
    'Expected keyword arguments, {...}, or dict(...) in TypedDict constructor'  # type: Final
TYPEDDICT_KEY_MUST_BE_STRING_LITERAL = \
    'Expected TypedDict key to be string literal'  # type: Final
MALFORMED_ASSERT = 'Assertion is always true, perhaps remove parentheses?'  # type: Final
DUPLICATE_TYPE_SIGNATURES = 'Function has duplicate type signatures'  # type: Final
DESCRIPTOR_SET_NOT_CALLABLE = "{}.__set__ is not callable"  # type: Final
DESCRIPTOR_GET_NOT_CALLABLE = "{}.__get__ is not callable"  # type: Final
MODULE_LEVEL_GETATTRIBUTE = '__getattribute__ is not valid at the module level'  # type: Final

# Generic
GENERIC_INSTANCE_VAR_CLASS_ACCESS = \
    'Access to generic instance variables via class is ambiguous'  # type: Final
GENERIC_CLASS_VAR_ACCESS = \
    'Access to generic class variables is ambiguous'  # type: Final
BARE_GENERIC = 'Missing type parameters for generic type {}'  # type: Final
IMPLICIT_GENERIC_ANY_BUILTIN = \
    'Implicit generic "Any". Use "{}" and specify generic parameters'  # type: Final

# TypeVar
INCOMPATIBLE_TYPEVAR_VALUE = 'Value of type variable "{}" of {} cannot be {}'  # type: Final
CANNOT_USE_TYPEVAR_AS_EXPRESSION = \
    'Type variable "{}.{}" cannot be used as an expression'  # type: Final

# Super
TOO_MANY_ARGS_FOR_SUPER = 'Too many arguments for "super"'  # type: Final
TOO_FEW_ARGS_FOR_SUPER = 'Too few arguments for "super"'  # type: Final
SUPER_WITH_SINGLE_ARG_NOT_SUPPORTED = '"super" with a single argument not supported'  # type: Final
UNSUPPORTED_ARG_1_FOR_SUPER = 'Unsupported argument 1 for "super"'  # type: Final
UNSUPPORTED_ARG_2_FOR_SUPER = 'Unsupported argument 2 for "super"'  # type: Final
SUPER_VARARGS_NOT_SUPPORTED = 'Varargs not supported with "super"'  # type: Final
SUPER_POSITIONAL_ARGS_REQUIRED = '"super" only accepts positional arguments'  # type: Final
SUPER_ARG_2_NOT_INSTANCE_OF_ARG_1 = \
    'Argument 2 for "super" not an instance of argument 1'  # type: Final
TARGET_CLASS_HAS_NO_BASE_CLASS = 'Target class has no base class'  # type: Final
SUPER_OUTSIDE_OF_METHOD_NOT_SUPPORTED = \
    'super() outside of a method is not supported'  # type: Final
SUPER_ENCLOSING_POSITIONAL_ARGS_REQUIRED = \
    'super() requires one or more positional arguments in enclosing function'  # type: Final

# Self-type
MISSING_OR_INVALID_SELF_TYPE = \
    "Self argument missing for a non-static method (or an invalid type for self)"  # type: Final
ERASED_SELF_TYPE_NOT_SUPERTYPE = \
    'The erased type of self "{}" is not a supertype of its class "{}"'  # type: Final
INVALID_SELF_TYPE_OR_EXTRA_ARG = \
    "Invalid type for self, or extra argument type in function annotation"  # type: Final

# Final
CANNOT_INHERIT_FROM_FINAL = 'Cannot inherit from final class "{}"'  # type: Final
DEPENDENT_FINAL_IN_CLASS_BODY = \
    "Final name declared in class body cannot depend on type variables"  # type: Final
CANNOT_ACCESS_FINAL_INSTANCE_ATTR = \
    'Cannot access final instance attribute "{}" on class object'  # type: Final

# ClassVar
CANNOT_OVERRIDE_INSTANCE_VAR = \
    'Cannot override instance variable (previously declared on base class "{}") with class ' \
    'variable'  # type: Final
CANNOT_OVERRIDE_CLASS_VAR = \
    'Cannot override class variable (previously declared on base class "{}") with instance ' \
    'variable'  # type: Final

# Protocol
RUNTIME_PROTOCOL_EXPECTED = \
    'Only @runtime_checkable protocols can be used with instance and class checks'  # type: Final
CANNOT_INSTANTIATE_PROTOCOL = 'Cannot instantiate protocol class "{}"'  # type: Final
