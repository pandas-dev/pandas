# Check if SSE/AVX instructions are available on the machine where
# the project is compiled.

if(CMAKE_SYSTEM_NAME MATCHES "Linux")
   exec_program(cat ARGS "/proc/cpuinfo" OUTPUT_VARIABLE CPUINFO)

   string(REGEX REPLACE "^.*(sse2).*$" "\\1" SSE_THERE "${CPUINFO}")
   string(COMPARE EQUAL "sse2" "${SSE_THERE}" SSE2_TRUE)
   if(SSE2_TRUE)
      set(SSE2_FOUND true CACHE BOOL "SSE2 available on host")
   else()
      set(SSE2_FOUND false CACHE BOOL "SSE2 available on host")
   endif()

   string(REGEX REPLACE "^.*(avx2).*$" "\\1" SSE_THERE "${CPUINFO}")
   string(COMPARE EQUAL "avx2" "${SSE_THERE}" AVX2_TRUE)
   if(AVX2_TRUE)
      set(AVX2_FOUND true CACHE BOOL "AVX2 available on host")
   else()
      set(AVX2_FOUND false CACHE BOOL "AVX2 available on host")
   endif()

elseif(CMAKE_SYSTEM_NAME MATCHES "Darwin")
   exec_program("/usr/sbin/sysctl -a | grep machdep.cpu.features" OUTPUT_VARIABLE CPUINFO)
   string(REGEX REPLACE "^.*[^S](SSE2).*$" "\\1" SSE_THERE "${CPUINFO}")
   string(COMPARE EQUAL "SSE2" "${SSE_THERE}" SSE2_TRUE)
   if(SSE2_TRUE)
      set(SSE2_FOUND true CACHE BOOL "SSE2 available on host")
   else()
      set(SSE2_FOUND false CACHE BOOL "SSE2 available on host")
   endif()

   exec_program("/usr/sbin/sysctl -a | grep machdep.cpu.leaf7_features" OUTPUT_VARIABLE CPUINFO)
   string(REGEX REPLACE "^.*(AVX2).*$" "\\1" SSE_THERE "${CPUINFO}")
   string(COMPARE EQUAL "AVX2" "${SSE_THERE}" AVX2_TRUE)
   if(AVX2_TRUE)
      set(AVX2_FOUND true CACHE BOOL "AVX2 available on host")
   else()
      set(AVX2_FOUND false CACHE BOOL "AVX2 available on host")
   endif()

elseif(CMAKE_SYSTEM_NAME MATCHES "Windows")
   # TODO.  For now supposing SSE2 is safe enough
   set(SSE2_FOUND true  CACHE BOOL "SSE2 available on host")
   set(AVX2_FOUND false CACHE BOOL "AVX2 available on host")
else()
   set(SSE2_FOUND true  CACHE BOOL "SSE2 available on host")
   set(AVX2_FOUND false CACHE BOOL "AVX2 available on host")
endif()

if(NOT SSE2_FOUND)
   message(STATUS "Could not find hardware support for SSE2 on this machine.")
endif()
if(NOT AVX2_FOUND)
   message(STATUS "Could not find hardware support for AVX2 on this machine.")
endif()

mark_as_advanced(SSE2_FOUND AVX2_FOUND)
