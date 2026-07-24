# Copyright 2017 Virgil Dupras

# This software is licensed under the "BSD" License as described in the "LICENSE" file,
# which should be included with this package. The terms are also available at
# http://www.hardcoded.net/licenses/bsd_license

# This is a reimplementation of plat_other.py with reference to the
# freedesktop.org trash specification:
#   [1] http://www.freedesktop.org/wiki/Specifications/trash-spec
#   [2] http://www.ramendik.ru/docs/trashspec.html
# See also:
#   [3] http://standards.freedesktop.org/basedir-spec/basedir-spec-latest.html
#
# For external volumes this implementation will raise an exception if it can't
# find or create the user's trash directory.

from __future__ import unicode_literals

import errno
import shutil
import os
import os.path as op
from datetime import datetime
import stat
from urllib.parse import quote
from send2trash.util import preprocess_paths
from send2trash.exceptions import TrashPermissionError

FILES_DIR = b"files"
INFO_DIR = b"info"
INFO_SUFFIX = b".trashinfo"

# Default of ~/.local/share [3]
XDG_DATA_HOME = op.expanduser(os.environb.get(b"XDG_DATA_HOME", b"~/.local/share"))
HOMETRASH_B = op.join(XDG_DATA_HOME, b"Trash")
HOMETRASH = os.fsdecode(HOMETRASH_B)

uid = os.getuid()
TOPDIR_TRASH = b".Trash"
TOPDIR_FALLBACK = b".Trash-" + str(uid).encode("ascii")


def is_parent(parent, path):
    path = op.realpath(path)  # In case it's a symlink
    if isinstance(path, str):
        path = os.fsencode(path)
    parent = op.realpath(parent)
    if isinstance(parent, str):
        parent = os.fsencode(parent)
    return path.startswith(parent)


def format_date(date):
    return date.strftime("%Y-%m-%dT%H:%M:%S")


def info_for(src, topdir):
    # ...it MUST not include a ".." directory, and for files not "under" that
    # directory, absolute pathnames must be used. [2]
    if topdir is None or not is_parent(topdir, src):
        src = op.abspath(src)
    else:
        src = op.relpath(src, topdir)

    info = "[Trash Info]\n"
    info += "Path=" + quote(src) + "\n"
    info += "DeletionDate=" + format_date(datetime.now()) + "\n"
    return info


def check_create(folder):
    # use 0700 for paths [3]
    if not op.exists(folder):
        os.makedirs(folder, 0o700)


def trash_move(src, dst, topdir=None, cross_dev=False):
    file_name = op.basename(src)
    files_path = op.join(dst, FILES_DIR)
    info_path = op.join(dst, INFO_DIR)
    base_name, ext = op.splitext(file_name)

    counter = 0
    dest_name = file_name
    while op.exists(op.join(files_path, dest_name)) or op.exists(op.join(info_path, dest_name + INFO_SUFFIX)):
        counter += 1
        dest_name = base_name + b" " + str(counter).encode("ascii") + ext

    check_create(files_path)
    check_create(info_path)

    with open(op.join(info_path, dest_name + INFO_SUFFIX), "w") as f:
        f.write(info_for(src, topdir))
    dest_path = op.join(files_path, dest_name)
    if cross_dev:
        shutil.move(os.fsdecode(src), os.fsdecode(dest_path))
    else:
        os.rename(src, dest_path)


def find_mount_point(path):
    # Even if something's wrong, "/" is a mount point, so the loop will exit.
    # Use realpath in case it's a symlink
    path = op.realpath(path)  # Required to avoid infinite loop
    while not op.ismount(path):  # Note ismount() does not always detect mounts
        path = op.split(path)[0]
    return path


def find_ext_volume_global_trash(volume_root):
    # from [2] Trash directories (1) check for a .Trash dir with the right
    # permissions set.
    trash_dir = op.join(volume_root, TOPDIR_TRASH)
    if not op.exists(trash_dir):
        return None

    mode = os.lstat(trash_dir).st_mode
    # vol/.Trash must be a directory, cannot be a symlink, and must have the
    # sticky bit set.
    if not op.isdir(trash_dir) or op.islink(trash_dir) or not mode & stat.S_ISVTX:
        return None

    trash_dir = op.join(trash_dir, str(uid).encode("ascii"))
    try:
        check_create(trash_dir)
    except OSError:
        return None
    return trash_dir


def find_ext_volume_fallback_trash(volume_root):
    # from [2] Trash directories (1) create a .Trash-$uid dir.
    trash_dir = op.join(volume_root, TOPDIR_FALLBACK)
    # Try to make the directory, if we lack permission, raise TrashPermissionError
    try:
        check_create(trash_dir)
    except OSError as e:
        if e.errno == errno.EACCES:
            raise TrashPermissionError(e.filename) from e
        raise
    return trash_dir


def find_ext_volume_trash(volume_root):
    trash_dir = find_ext_volume_global_trash(volume_root)
    if trash_dir is None:
        trash_dir = find_ext_volume_fallback_trash(volume_root)
    return trash_dir


# Pull this out so it's easy to stub (to avoid stubbing lstat itself)
def get_dev(path):
    return os.lstat(path).st_dev


def send2trash(paths):
    paths = preprocess_paths(paths)
    for path in paths:
        if isinstance(path, str):
            path_b = os.fsencode(path)
        elif isinstance(path, bytes):
            path_b = path
        else:
            raise TypeError(f"str, bytes or PathLike expected, not {type(path)}")

        if not op.exists(path_b):
            raise OSError(errno.ENOENT, f"File not found: {path}")
        # ...should check whether the user has the necessary permissions to delete
        # it, before starting the trashing operation itself. [2]
        if not os.access(path_b, os.W_OK):
            raise OSError(errno.EACCES, f"Permission denied: {path}")

        path_dev = get_dev(path_b)
        # If XDG_DATA_HOME or HOMETRASH do not yet exist we need to stat the
        # home directory, and these paths will be created further on if needed.
        trash_dev = get_dev(op.expanduser(b"~"))

        # if the file to be trashed is on the same device as HOMETRASH we
        # want to move it there.
        if path_dev == trash_dev:
            topdir = XDG_DATA_HOME
            dest_trash = HOMETRASH_B
        else:
            topdir = find_mount_point(path_b)
            trash_dev = get_dev(topdir)
            if trash_dev != path_dev:
                raise OSError(f"Couldn't find mount point for {path}")
            dest_trash = find_ext_volume_trash(topdir)
        try:
            trash_move(path_b, dest_trash, topdir)
        except OSError as error:
            # Cross link errors default back to HOMETRASH
            if error.errno == errno.EXDEV:
                trash_move(path_b, HOMETRASH_B, XDG_DATA_HOME, cross_dev=True)
            else:
                raise
