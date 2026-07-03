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
