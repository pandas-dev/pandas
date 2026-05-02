host_cpu = host_machine.cpu_family()

simd_config = configuration_data()
if host_cpu == 'aarch64'
    simd_target = 'neon'
    simd_config.set('PANDAS_HAVE_NEON', 1)
elif host_cpu == 'x86_64'
    simd_target = 'sse2'
    simd_config.set('PANDAS_HAVE_SSE2', 1)
else
    simd_target = 'scalar'
    simd_config.set('PANDAS_HAVE_SCALAR', 1)
endif

configure_file(
    output: 'pandas_simd_config.h',
    configuration: simd_config,
)

simd_config_inc = include_directories('.')

moments_libs = []
foreach arch_name, arch_flags : supported_simd_archs
    moments_libs += static_library(
        'moments_simd_@0@'.format(arch_name),
        'moments_inst_@0@.cpp'.format(arch_name),
        include_directories: [inc_pd],
        dependencies: [xsimd_dep],
        cpp_args: arch_flags,
    )
endforeach

moments_libs += static_library(
    'moments_base',
    'moments.cpp',
    include_directories: [inc_pd],
    dependencies: [xsimd_dep],
)

moments_simd_dep = declare_dependency(
    link_with: moments_libs,
    include_directories: [inc_pd],
)
