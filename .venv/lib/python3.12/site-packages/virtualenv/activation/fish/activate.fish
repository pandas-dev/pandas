# This file must be used using `source bin/activate.fish` *within a running fish ( http://fishshell.com ) session*.
# Do not run it directly.

function deactivate -d 'Exit virtualenv mode and return to the normal environment.'
    # reset old environment variables
    if test -n "$_OLD_VIRTUAL_PATH"
        set -gx PATH $_OLD_VIRTUAL_PATH
        set -e _OLD_VIRTUAL_PATH
    end

    if test -n __TCL_LIBRARY__
      if test -n "$_OLD_VIRTUAL_TCL_LIBRARY";
        set -gx TCL_LIBRARY "$_OLD_VIRTUAL_TCL_LIBRARY";
        set -e _OLD_VIRTUAL_TCL_LIBRARY;
      else;
        set -e TCL_LIBRARY;
      end
    end
    if test -n __TK_LIBRARY__
      if test -n "$_OLD_VIRTUAL_TK_LIBRARY";
        set -gx TK_LIBRARY "$_OLD_VIRTUAL_TK_LIBRARY";
        set -e _OLD_VIRTUAL_TK_LIBRARY;
      else;
        set -e TK_LIBRARY;
      end
    end

    if test -n "$_OLD_PKG_CONFIG_PATH"
        set -gx PKG_CONFIG_PATH "$_OLD_PKG_CONFIG_PATH"
        set -e _OLD_PKG_CONFIG_PATH
    end

    if test -n "$_OLD_VIRTUAL_PYTHONHOME"
        set -gx PYTHONHOME "$_OLD_VIRTUAL_PYTHONHOME"
        set -e _OLD_VIRTUAL_PYTHONHOME
    end

    if test -n "$_OLD_FISH_PROMPT_OVERRIDE"
       and functions -q _old_fish_prompt
        # Set an empty local `$fish_function_path` to allow the removal of `fish_prompt` using `functions -e`.
        set -l fish_function_path

        # Erase virtualenv's `fish_prompt` and restore the original.
        functions -e fish_prompt
        functions -c _old_fish_prompt fish_prompt
        functions -e _old_fish_prompt
        set -e _OLD_FISH_PROMPT_OVERRIDE
    end

    set -e VIRTUAL_ENV
    set -e VIRTUAL_ENV_PROMPT

    if test "$argv[1]" != 'nondestructive'
        # Self-destruct!
        functions -e pydoc
        functions -e deactivate
    end
end

# Unset irrelevant variables.
deactivate nondestructive

set -gx VIRTUAL_ENV __VIRTUAL_ENV__

set -gx _OLD_PKG_CONFIG_PATH "$PKG_CONFIG_PATH"
set -gx PKG_CONFIG_PATH "$VIRTUAL_ENV/lib/pkgconfig:$PKG_CONFIG_PATH"

set -gx _OLD_VIRTUAL_PATH $PATH
set -gx PATH "$VIRTUAL_ENV"'/'__BIN_NAME__ $PATH

if test -n __TCL_LIBRARY__
  if set -q TCL_LIBRARY;
    set -gx _OLD_VIRTUAL_TCL_LIBRARY $TCL_LIBRARY;
  end
  set -gx TCL_LIBRARY '__TCL_LIBRARY__'
end
if test -n __TK_LIBRARY__
  if set -q TK_LIBRARY;
    set -gx _OLD_VIRTUAL_TK_LIBRARY $TK_LIBRARY;
  end
  set -gx TK_LIBRARY '__TK_LIBRARY__'
end

# Prompt override provided?
# If not, just use the environment name.
if test -n __VIRTUAL_PROMPT__
    set -gx VIRTUAL_ENV_PROMPT __VIRTUAL_PROMPT__
else
    set -gx VIRTUAL_ENV_PROMPT (basename "$VIRTUAL_ENV")
end

# Unset `$PYTHONHOME` if set.
if set -q PYTHONHOME
    set -gx _OLD_VIRTUAL_PYTHONHOME $PYTHONHOME
    set -e PYTHONHOME
end

function pydoc
    python -m pydoc $argv
end

if test -z "$VIRTUAL_ENV_DISABLE_PROMPT"
    # Copy the current `fish_prompt` function as `_old_fish_prompt`.
    functions -c fish_prompt _old_fish_prompt

    function fish_prompt
        set -l old_status $status
        # Run the user's prompt first; it might depend on (pipe)status.
        set -l prompt (_old_fish_prompt)

        printf '(%s) ' $VIRTUAL_ENV_PROMPT

        string join -- \n $prompt # handle multi-line prompts
        echo "exit $old_status" | .
    end

    set -gx _OLD_FISH_PROMPT_OVERRIDE "$VIRTUAL_ENV"
end
