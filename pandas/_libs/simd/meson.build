host_cpu = host_machine.cpu_family()

host_baseline = {'aarch64': 'neon', 'x86_64': 'sse2'}

simd_config = configuration_data()
if get_option('simd-mode') == 'scalar' or host_cpu not in host_baseline
    simd_target = 'scalar'
else
    simd_target = host_baseline[host_cpu]
    simd_config.set('PANDAS_HAVE_' + simd_target.to_upper(), 1)
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
    cpp_args: strict_c_args,
)

moments_simd_dep = declare_dependency(
    link_with: moments_lib,
    include_directories: [inc_pd],
)
