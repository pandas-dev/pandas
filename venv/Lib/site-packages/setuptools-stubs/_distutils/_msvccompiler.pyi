from .compilers.C import msvc

__all__ = ["MSVCCompiler"]

MSVCCompiler = msvc.Compiler
