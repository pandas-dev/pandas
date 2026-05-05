# All architectures we might support
# Key is the architecture name used in file suffixes and macros
is_msvc_syntax = cxx.get_argument_syntax() == 'msvc'
simd_arch_flags = {
    'sse2': is_msvc_syntax ? ['/arch:SSE2'] : ['-msse2'],
    'avx2': is_msvc_syntax ? ['/arch:AVX2'] : ['-mavx2'],
    'avx512cd': is_msvc_syntax ? ['/arch:AVX512'] : ['-mavx512cd'],
    'neon': [],
}

simd_config = configuration_data()
supported_simd_archs = {}
if get_option('simd').allowed()
    foreach name, flags : simd_arch_flags
        if host_machine.cpu_family() in ['x86', 'x86_64']
            if name in ['sse2', 'avx2', 'avx512cd']
                if cxx.has_multi_arguments(flags)
                    supported_simd_archs += {name: flags}
                endif
            endif
        elif host_machine.cpu_family() == 'aarch64'
            if name == 'neon'
                supported_simd_archs += {name: flags}
            endif
        endif
    endforeach
endif

foreach arch_name, arch_flags : supported_simd_archs
    simd_config.set('PANDAS_HAVE_@0@'.format(arch_name.to_upper()), 1)
endforeach

# Ensure scalar version on all architectures for now...
simd_config.set('PANDAS_HAVE_SCALAR', 1)

configure_file(
    output: 'pandas_simd_config.h',
    configuration: simd_config,
)

simd_config_inc = include_directories('.')
