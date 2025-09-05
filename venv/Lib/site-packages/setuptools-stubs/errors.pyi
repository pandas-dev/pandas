from ._distutils import errors as _distutils_errors

ByteCompileError = _distutils_errors.DistutilsByteCompileError
CCompilerError = _distutils_errors.CCompilerError
ClassError = _distutils_errors.DistutilsClassError
CompileError = _distutils_errors.CompileError
ExecError = _distutils_errors.DistutilsExecError
FileError = _distutils_errors.DistutilsFileError
InternalError = _distutils_errors.DistutilsInternalError
LibError = _distutils_errors.LibError
LinkError = _distutils_errors.LinkError
ModuleError = _distutils_errors.DistutilsModuleError
OptionError = _distutils_errors.DistutilsOptionError
PlatformError = _distutils_errors.DistutilsPlatformError
PreprocessError = _distutils_errors.PreprocessError
SetupError = _distutils_errors.DistutilsSetupError
TemplateError = _distutils_errors.DistutilsTemplateError
UnknownFileError = _distutils_errors.UnknownFileError
BaseError = _distutils_errors.DistutilsError

class InvalidConfigError(OptionError): ...
class RemovedConfigError(OptionError): ...
class RemovedCommandError(BaseError, RuntimeError): ...
class PackageDiscoveryError(BaseError, RuntimeError): ...
