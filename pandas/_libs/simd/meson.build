# All architectures we might support
# Key is the architecture name used in file suffixes and macros
is_msvc_syntax = cxx.get_argument_syntax() == 'msvc'
simd_x86_flags = {
    'sse2': is_msvc_syntax ? ['/arch:SSE2'] : ['-msse2'],
    'avx2': is_msvc_syntax ? ['/arch:AVX2'] : ['-mavx2'],
    'avx512cd': is_msvc_syntax ? ['/arch:AVX512'] : ['-mavx512cd'],
}

simd_config = configuration_data()
supported_simd_archs = {}
if get_option('simd').allowed()
    if host_machine.cpu_family() == 'aarch64'
        supported_simd_archs += {'neon': []}
    elif host_machine.cpu_family() in ['x86', 'x86_64']
        foreach name, flags : simd_x86_flags
            if cxx.has_multi_arguments(flags)
                supported_simd_archs += {name: flags}
            endif
        endforeach
    endif
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
