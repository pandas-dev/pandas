import glob
import os
import shutil
import subprocess
import sys

if os.name == "nt":
    py_ver = f"{sys.version_info.major}.{sys.version_info.minor}"
    is_32_bit = os.getenv("IS_32_BIT") == "true"
    try:
        wheel_dir = sys.argv[1]
        wheel_path = glob.glob(f"{wheel_dir}/*.whl")[0]
    except IndexError:
        # Not passed
        wheel_path = None
    print(f"IS_32_BIT is {is_32_bit}")
    print(f"Path to built wheel is {wheel_path}")
    if is_32_bit:
        sys.exit(0)  # No way to test Windows 32-bit(no docker image)
    if wheel_path is None:
        raise ValueError("Wheel path must be passed in if on 64-bit Windows")
    print(f"Pulling docker image to test Windows 64-bit Python {py_ver}")
    subprocess.run(f"docker pull python:{py_ver}-windowsservercore", check=True)
    pandas_base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    print(f"pandas project dir is {pandas_base_dir}")
    dist_dir = os.path.join(pandas_base_dir, "dist")
    print(f"Copying wheel into pandas_base_dir/dist ({dist_dir})")
    os.mkdir(dist_dir)
    shutil.copy(wheel_path, dist_dir)
    print(os.listdir(dist_dir))
    subprocess.run(
        rf"docker run -v %cd%:c:\pandas "
        f"python:{py_ver}-windowsservercore /pandas/ci/test_wheels_windows.bat",
        check=True,
        shell=True,
        cwd=pandas_base_dir,
    )
else:
    import pandas as pd

    pd.test(
        extra_args=[
            "-m not clipboard and not single_cpu",
            "--skip-slow",
            "--skip-network",
            "--skip-db",
            "-n=2",
        ]
    )
    pd.test(
        extra_args=[
            "-m not clipboard and single_cpu",
            "--skip-slow",
            "--skip-network",
            "--skip-db",
        ]
    )
