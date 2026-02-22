# This file must be used with "source bin/activate" *from bash*
# you cannot run it directly


if [ "${BASH_SOURCE-}" = "$0" ]; then
    echo "You must source this script: \$ source $0" >&2
    exit 33
fi

deactivate () {
    unset -f pydoc >/dev/null 2>&1 || true

    # reset old environment variables
    if [ -n "${_OLD_VIRTUAL_PATH:-}" ] ; then
        PATH="$_OLD_VIRTUAL_PATH"
        export PATH
        unset _OLD_VIRTUAL_PATH
    fi
    if [ -n "${_OLD_VIRTUAL_PYTHONHOME:-}" ] ; then
        PYTHONHOME="$_OLD_VIRTUAL_PYTHONHOME"
        export PYTHONHOME
        unset _OLD_VIRTUAL_PYTHONHOME
    fi

    if [ -n "${_OLD_VIRTUAL_TCL_LIBRARY:-}" ]; then
        TCL_LIBRARY="$_OLD_VIRTUAL_TCL_LIBRARY"
        export TCL_LIBRARY
        unset _OLD_VIRTUAL_TCL_LIBRARY
    fi
    if [ -n "${_OLD_VIRTUAL_TK_LIBRARY:-}" ]; then
        TK_LIBRARY="$_OLD_VIRTUAL_TK_LIBRARY"
        export TK_LIBRARY
        unset _OLD_VIRTUAL_TK_LIBRARY
    fi
    if [ -n "${_OLD_PKG_CONFIG_PATH:-}" ]; then
        PKG_CONFIG_PATH="$_OLD_PKG_CONFIG_PATH"
        export PKG_CONFIG_PATH
        unset _OLD_PKG_CONFIG_PATH
    fi

    # The hash command must be called to get it to forget past
    # commands. Without forgetting past commands the $PATH changes
    # we made may not be respected
    hash -r 2>/dev/null

    if [ -n "${_OLD_VIRTUAL_PS1:-}" ] ; then
        PS1="$_OLD_VIRTUAL_PS1"
        export PS1
        unset _OLD_VIRTUAL_PS1
    fi

    unset VIRTUAL_ENV
    unset VIRTUAL_ENV_PROMPT
    if [ ! "${1-}" = "nondestructive" ] ; then
    # Self destruct!
        unset -f deactivate
    fi
}

# unset irrelevant variables
deactivate nondestructive

if [ ! -d __VIRTUAL_ENV__ ]; then
    echo "Virtual environment directory __VIRTUAL_ENV__ does not exist!" >&2
    CURRENT_PATH=$(realpath "${0}")
    CURRENT_DIR=$(dirname "${CURRENT_PATH}")
    VIRTUAL_ENV="$(realpath "${CURRENT_DIR}/../")"
else
    VIRTUAL_ENV=__VIRTUAL_ENV__
fi

case "$(uname)" in
    CYGWIN*|MSYS*|MINGW*)
        VIRTUAL_ENV=$(cygpath "$VIRTUAL_ENV")
        ;;
esac
export VIRTUAL_ENV

_OLD_VIRTUAL_PATH="$PATH"
PATH="$VIRTUAL_ENV/"__BIN_NAME__":$PATH"
export PATH

_OLD_PKG_CONFIG_PATH="${PKG_CONFIG_PATH:-}"
PKG_CONFIG_PATH="${VIRTUAL_ENV}/lib/pkgconfig${PKG_CONFIG_PATH:+:${PKG_CONFIG_PATH}}"
export PKG_CONFIG_PATH

if [ "x"__VIRTUAL_PROMPT__ != x ] ; then
    VIRTUAL_ENV_PROMPT=__VIRTUAL_PROMPT__
else
    VIRTUAL_ENV_PROMPT=$(basename "$VIRTUAL_ENV")
fi
export VIRTUAL_ENV_PROMPT

# unset PYTHONHOME if set
if [ -n "${PYTHONHOME:-}" ] ; then
    _OLD_VIRTUAL_PYTHONHOME="$PYTHONHOME"
    unset PYTHONHOME
fi

if [ __TCL_LIBRARY__ != "" ]; then
    if [ -n "${TCL_LIBRARY:-}" ] ; then
        _OLD_VIRTUAL_TCL_LIBRARY="$TCL_LIBRARY"
    fi
    TCL_LIBRARY=__TCL_LIBRARY__
    export TCL_LIBRARY
fi

if [ __TK_LIBRARY__ != "" ]; then
    if [ -n "${TK_LIBRARY:-}" ] ; then
        _OLD_VIRTUAL_TK_LIBRARY="$TK_LIBRARY"
    fi
    TK_LIBRARY=__TK_LIBRARY__
    export TK_LIBRARY
fi

if [ -z "${VIRTUAL_ENV_DISABLE_PROMPT-}" ] ; then
    _OLD_VIRTUAL_PS1="${PS1-}"
    PS1="(${VIRTUAL_ENV_PROMPT}) ${PS1-}"
    export PS1
fi

# Make sure to unalias pydoc if it's already there
alias pydoc 2>/dev/null >/dev/null && unalias pydoc || true

pydoc () {
    python -m pydoc "$@"
}

# The hash command must be called to get it to forget past
# commands. Without forgetting past commands the $PATH changes
# we made may not be respected
hash -r 2>/dev/null || true
