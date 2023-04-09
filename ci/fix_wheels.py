"""
This file "repairs" our Windows wheels by copying the necessary DLLs for pandas to run
on a barebones Windows installation() into the wheel.

NOTE: The paths for the DLLs are hard-coded to the location of the Visual Studio
redistributables
"""
import os
import shutil
import subprocess
import sys

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
repaired_wheel_path = os.path.join(dest_dir, wheel_name)

try:
# Use the wheel CLI instead of manipulating zipfiles, since the CLI will
# take care of rebuilding the hashes found in the record file
tmp_dir = os.path.join(dest_dir, "tmp")
subprocess.run(["wheel", "unpack", f"-d {tmp_dir}", wheel_path], check=True)
base_redist_dir = (
    f"C:/Program Files (x86)/Microsoft Visual Studio/2019/"
    f"Enterprise/VC/Redist/MSVC/14.29.30133/{PYTHON_ARCH}/"
    f"Microsoft.VC142.CRT/"
)
required_dlls = ["msvcp140.dll", "concrt140.dll"]
if not is_32:
    required_dlls += "vcruntime140_1.dll"
for dll in required_dlls:
    src = os.path.join(base_redist_dir, dll)
    shutil.copy(src, tmp_dir)
subprocess.run(["wheel", "pack", tmp_dir, f"-d {dest_dir}"], check=True)

if not success:
    os.remove(repaired_wheel_path)
    sys.exit(1)
print(f"Successfully repaired wheel was written to {repaired_wheel_path}")
