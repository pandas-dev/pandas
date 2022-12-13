import os
import shutil
import sys
import zipfile

try:
    if len(sys.argv) != 3:
        raise ValueError(
            "User must pass the path to the wheel and the destination directory."
        )
    wheel_path = sys.argv[1]
    dest_dir = sys.argv[2]
    # Figure out whether we are building on 32 or 64 bit python
    is_32 = sys.maxsize <= 2**32
    PYTHON_ARCH = "x86" if is_32 else "x64"
except ValueError:
    # Too many/little values to unpack
    raise ValueError(
        "User must pass the path to the wheel and the destination directory."
    )
# Wheels are zip files
if not os.path.isdir(dest_dir):
    print(f"Created directory {dest_dir}")
    os.mkdir(dest_dir)
shutil.copy(wheel_path, dest_dir)  # Remember to delete if process fails
wheel_name = os.path.basename(wheel_path)
success = True
exception = None
repaired_wheel_path = os.path.join(dest_dir, wheel_name)
with zipfile.ZipFile(repaired_wheel_path, "a") as zipf:
    try:
        # TODO: figure out how licensing works for the redistributables
        base_redist_dir = (
            f"C:/Program Files (x86)/Microsoft Visual Studio/2019/"
            f"Enterprise/VC/Redist/MSVC/14.29.30133/{PYTHON_ARCH}/"
            f"Microsoft.VC142.CRT/"
        )
        zipf.write(
            os.path.join(base_redist_dir, "msvcp140.dll"),
            "pandas/_libs/window/msvcp140.dll",
        )
        zipf.write(
            os.path.join(base_redist_dir, "concrt140.dll"),
            "pandas/_libs/window/concrt140.dll",
        )
        if not is_32:
            zipf.write(
                os.path.join(base_redist_dir, "vcruntime140_1.dll"),
                "pandas/_libs/window/vcruntime140_1.dll",
            )
    except Exception as e:
        success = False
        exception = e

if not success:
    os.remove(repaired_wheel_path)
    raise exception
print(f"Successfully repaired wheel was written to {repaired_wheel_path}")
