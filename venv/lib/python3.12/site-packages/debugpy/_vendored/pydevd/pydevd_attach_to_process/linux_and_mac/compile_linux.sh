set -e

ARCH="$(uname -m)"
case $ARCH in
    i*86) SUFFIX=x86;;
    x86_64*) SUFFIX=amd64;;
    *) echo >&2 "unsupported: $ARCH, this script may not work";;
esac

SRC="$(dirname "$0")/.."
g++ -std=c++11 -shared -fPIC -O2 -D_FORTIFY_SOURCE=2 -nostartfiles -fstack-protector-strong $SRC/linux_and_mac/attach.cpp -o $SRC/attach_linux_$SUFFIX.so
