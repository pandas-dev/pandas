import glob
import os
import shutil
import subprocess
import sys

if os.name == "nt":
    is_32_bit = os.getenv("IS_32_BIT") == "true"
    print(f"IS_32_BIT is {is_32_bit}")

    try:
        wheel_dir = sys.argv[1]
        wheel_path = glob.glob(f"{wheel_dir}/*.whl")[0]
    except IndexError:
        # Not passed
        wheel_path = None
    print(f"Path to built wheel is {wheel_path}")

    print("Verifying file hashes in wheel RECORD file")
    try:
        tmp_dir = "tmp"
        subprocess.run(
            ["wheel", "unpack", wheel_path, "-d", tmp_dir],
            check=True,
            capture_output=True,
        )
    except subprocess.CalledProcessError as err:
        print("wheel RECORD file hash verification failed.")
        print(err.stdout)
        print(err.stderr)
        sys.exit(1)
    finally:
        shutil.rmtree(tmp_dir)

    if is_32_bit:
        print("Cannot test on Windows 32-bit(no docker image)")
        sys.exit(0)
    if wheel_path is None:
        raise ValueError("Wheel path must be passed in if on 64-bit Windows")

    pandas_base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    print(f"pandas project dir is {pandas_base_dir}")
    dist_dir = os.path.join(pandas_base_dir, "dist")
    print(f"Copying wheel into pandas_base_dir/dist ({dist_dir})")
    os.mkdir(dist_dir)
    shutil.copy(wheel_path, dist_dir)
    print(os.listdir(dist_dir))

    py_ver = f"{sys.version_info.major}.{sys.version_info.minor}"
    print(f"Pulling docker image to test Windows 64-bit Python {py_ver}")
    subprocess.run(f"docker pull python:{py_ver}-windowsservercore", check=True)
    subprocess.run(
        rf"docker run -v %cd%:c:\pandas "
        f"python:{py_ver}-windowsservercore /pandas/ci/test_wheels_windows.bat",
        check=True,
        shell=True,
        cwd=pandas_base_dir,
    )
else:
    import pandas as pd

    multi_args = [
        "-m not clipboard and not single_cpu and not slow and not network and not db",
        "-n 2",
    ]
    pd.test(extra_args=multi_args)
    pd.test(
        extra_args=[
            "-m not clipboard and single_cpu and not slow and not network and not db",
        ]
    )
