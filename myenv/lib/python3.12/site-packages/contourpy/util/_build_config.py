# _build_config.py.in is converted into _build_config.py during the meson build process.

from __future__ import annotations


def build_config() -> dict[str, str]:
    """
    Return a dictionary containing build configuration settings.

    All dictionary keys and values are strings, for example ``False`` is
    returned as ``"False"``.

        .. versionadded:: 1.1.0
    """
    return dict(
        # Python settings
        python_version="3.12",
        python_install_dir=r"/usr/local/lib/python3.12/site-packages/",
        python_path=r"/private/var/folders/24/8k48jl6d249_n_qfxwsl6xvm0000gn/T/build-env-a_k0k2bv/bin/python",

        # Package versions
        contourpy_version="1.2.1",
        meson_version="1.4.0",
        mesonpy_version="0.15.0",
        pybind11_version="2.12.0",

        # Misc meson settings
        meson_backend="ninja",
        build_dir=r"/Users/runner/work/contourpy/contourpy/.mesonpy-cqejv7k5/lib/contourpy/util",
        source_dir=r"/Users/runner/work/contourpy/contourpy/lib/contourpy/util",
        cross_build="False",

        # Build options
        build_options=r"-Dbuildtype=release -Db_ndebug=if-release -Db_vscrt=md -Dvsenv=True --native-file=/Users/runner/work/contourpy/contourpy/.mesonpy-cqejv7k5/meson-python-native-file.ini",
        buildtype="release",
        cpp_std="c++17",
        debug="False",
        optimization="3",
        vsenv="True",
        b_ndebug="if-release",
        b_vscrt="from_buildtype",

        # C++ compiler
        compiler_name="clang",
        compiler_version="13.0.0",
        linker_id="ld64",
        compile_command="c++",

        # Host machine
        host_cpu="x86_64",
        host_cpu_family="x86_64",
        host_cpu_endian="little",
        host_cpu_system="darwin",

        # Build machine, same as host machine if not a cross_build
        build_cpu="x86_64",
        build_cpu_family="x86_64",
        build_cpu_endian="little",
        build_cpu_system="darwin",
    )
