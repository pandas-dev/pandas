"""
This file "repairs" our Windows wheels by copying the necessary DLLs for pandas to run
on a barebones Windows installation into the wheel.

NOTE: The paths for the DLLs are hard-coded to the location of the Visual Studio
redistributables
"""
import os
import shutil
import subprocess
import sys
import zipfile

if len(sys.argv) != 3:
    raise ValueError(
        "User must pass the path to the wheel and the destination directory."
    )
wheel_path = sys.argv[1]
dest_dir = sys.argv[2]

# Figure out whether we are building on 32 or 64 bit python
is_32 = sys.maxsize <= 2**32
PYTHON_ARCH = "x86" if is_32 else "x64"

if not os.path.isdir(dest_dir):
    print(f"Creating directory {dest_dir}")
    os.mkdir(dest_dir)

# Use the wheel CLI for zipping up the wheel since the CLI will
# take care of rebuilding the hashes found in the record file
tmp_dir = os.path.join(dest_dir, "tmp")
with zipfile.ZipFile(wheel_path, "r") as f:
    # Extracting all the members of the zip
    # into a specific location.
    f.extractall(path=tmp_dir)
base_redist_dir = (
    f"C:/Program Files (x86)/Microsoft Visual Studio/2019/"
    f"Enterprise/VC/Redist/MSVC/14.29.30133/{PYTHON_ARCH}/"
    f"Microsoft.VC142.CRT/"
)
dest_dll_dir = os.path.join(tmp_dir, "pandas/_libs/window")
required_dlls = ["msvcp140.dll", "concrt140.dll"]
if not is_32:
    required_dlls += ["vcruntime140_1.dll"]
for dll in required_dlls:
    src = os.path.join(base_redist_dir, dll)
    shutil.copy(src, dest_dll_dir)
try:
    subprocess.run(
        ["wheel", "pack", tmp_dir, "-d", dest_dir], check=True, capture_output=True
    )
    print("Successfully repaired wheel")
except subprocess.CalledProcessError as err:
    print("Failed to add DLLS to wheel.")
    print(err.stdout)
    print(err.stderr)
    sys.exit(1)
