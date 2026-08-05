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

moments_lib = static_library(
    'moments_@0@'.format(simd_target),
    'moments_@0@.cpp'.format(simd_target),
    include_directories: [inc_pd],
    dependencies: [xsimd_dep],
)

moments_simd_dep = declare_dependency(
    link_with: moments_lib,
    include_directories: [inc_pd],
)
