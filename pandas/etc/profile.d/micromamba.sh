# Copyright (C) 2012 Anaconda, Inc
# SPDX-License-Identifier: BSD-3-Clause

__mamba_exe() (
    "$MAMBA_EXE" "${@}"
)

__mamba_hashr() {
    if [ -n "${ZSH_VERSION:+x}" ]; then
        \rehash
    elif [ -n "${POSH_VERSION:+x}" ]; then
        :  # pass
    else
        \hash -r
    fi
}

__mamba_xctivate() {
    \local ask_conda
    ask_conda="$(PS1="${PS1:-}" __mamba_exe shell "${@}" --shell bash)" || \return
    \eval "${ask_conda}"
    __mamba_hashr
}

micromamba() {
    \local cmd="${1-__missing__}"
    case "${cmd}" in
        activate|reactivate|deactivate)
            __mamba_xctivate "${@}"
            ;;
        install|update|upgrade|remove|uninstall)
            __mamba_exe "${@}" || \return
            __mamba_xctivate reactivate
            ;;
        self-update)
            __mamba_exe "${@}" || \return

            # remove leftover backup file on Windows
            if [ -f "$MAMBA_EXE.bkup" ]; then
                rm -f "$MAMBA_EXE.bkup"
            fi
            ;;
        *)
            __mamba_exe "${@}"
            ;;
    esac
}

if [ -z "${CONDA_SHLVL+x}" ]; then
    \export CONDA_SHLVL=0
    # In dev-mode MAMBA_EXE is python.exe and on Windows
    # it is in a different relative location to condabin.
    if [ -n "${_CE_CONDA+x}" ] && [ -n "${WINDIR+x}" ]; then
        PATH="${MAMBA_ROOT_PREFIX}/condabin:${PATH}"
    else
        PATH="${MAMBA_ROOT_PREFIX}/condabin:${PATH}"
    fi
    \export PATH

    # We're not allowing PS1 to be unbound. It must at least be set.
    # However, we're not exporting it, which can cause problems when starting a second shell
    # via a first shell (i.e. starting zsh from bash).
    if [ -z "${PS1+x}" ]; then
        PS1=
    fi
fi
