<#
.Synopsis
Activate a Python virtual environment for the current PowerShell session.

.Description
Pushes the python executable for a virtual environment to the front of the
$Env:PATH environment variable and sets the prompt to signify that you are
in a Python virtual environment. Makes use of the command line switches as
well as the `pyvenv.cfg` file values present in the virtual environment.

.Parameter VenvDir
Path to the directory that contains the virtual environment to activate. The
default value for this is the parent of the directory that the activate.ps1
script is located within.

.Parameter Prompt
The prompt prefix to display when this virtual environment is activated. By
default, this prompt is the name of the virtual environment folder (VenvDir)
surrounded by parentheses and followed by a single space (ie. '(.venv) ').

.Example
activate.ps1
Activates the Python virtual environment that contains the activate.ps1 script.

.Example
activate.ps1 -Verbose
Activates the Python virtual environment that contains the activate.ps1 script,
and shows extra information about the activation as it executes.

.Example
activate.ps1 -VenvDir C:\Users\MyUser\Common\.venv
Activates the Python virtual environment located in the specified location.

.Example
activate.ps1 -Prompt "MyPython"
Activates the Python virtual environment that contains the activate.ps1 script,
and prefixes the current prompt with the specified string (surrounded in
parentheses) while the virtual environment is active.

.Notes
On Windows, it may be required to enable this activate.ps1 script by setting the
execution policy for the user. You can do this by issuing the following PowerShell
command:

PS C:\> Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser

For more information on Execution Policies:
https://go.microsoft.com/fwlink/?LinkID=135170

#>
Param(
    [Parameter(Mandatory = $false)]
    [String]
    $VenvDir,
    [Parameter(Mandatory = $false)]
    [String]
    $Prompt
)

function Get-PyVenvConfig(
    [String]
    $ConfigDir
) {
    Write-Verbose "Given ConfigDir=$ConfigDir, obtain values in pyvenv.cfg"

    $pyvenvConfigPath = Join-Path -Resolve -Path $ConfigDir -ChildPath 'pyvenv.cfg' -ErrorAction Continue

    $pyvenvConfig = @{ }

    if ($pyvenvConfigPath) {
        Write-Verbose "File exists, parse ``key = value`` lines"
        $pyvenvConfigContent = Get-Content -Path $pyvenvConfigPath

        $pyvenvConfigContent | ForEach-Object {
            $keyval = $PSItem -split "\s*=\s*", 2
            if ($keyval[0] -and $keyval[1]) {
                $val = $keyval[1]

                if ("'""".Contains($val.Substring(0, 1))) {
                    $val = $val.Substring(1, $val.Length - 2)
                }

                $pyvenvConfig[$keyval[0]] = $val
                Write-Verbose "Adding Key: '$($keyval[0])'='$val'"
            }
        }
    }
    return $pyvenvConfig
}

function global:deactivate([switch] $NonDestructive) {
    if (Test-Path variable:_OLD_VIRTUAL_PATH) {
        $env:PATH = $variable:_OLD_VIRTUAL_PATH
        Remove-Variable "_OLD_VIRTUAL_PATH" -Scope global
    }

    if (Test-Path variable:_OLD_VIRTUAL_TCL_LIBRARY) {
        $env:TCL_LIBRARY = $variable:_OLD_VIRTUAL_TCL_LIBRARY
        Remove-Variable "_OLD_VIRTUAL_TCL_LIBRARY" -Scope global
    } else {
        if (Test-Path env:TCL_LIBRARY) {
            Remove-Item env:TCL_LIBRARY -ErrorAction SilentlyContinue
        }
    }

    if (Test-Path variable:_OLD_VIRTUAL_TK_LIBRARY) {
        $env:TK_LIBRARY = $variable:_OLD_VIRTUAL_TK_LIBRARY
        Remove-Variable "_OLD_VIRTUAL_TK_LIBRARY" -Scope global
    } else {
        if (Test-Path env:TK_LIBRARY) {
            Remove-Item env:TK_LIBRARY -ErrorAction SilentlyContinue
        }
    }

    if (Test-Path variable:_OLD_PKG_CONFIG_PATH) {
        $env:PKG_CONFIG_PATH = $variable:_OLD_PKG_CONFIG_PATH
        Remove-Variable "_OLD_PKG_CONFIG_PATH" -Scope global
    }

    if (Test-Path function:_old_virtual_prompt) {
        $function:prompt = $function:_old_virtual_prompt
        Remove-Item function:\_old_virtual_prompt
    }

    if ($env:VIRTUAL_ENV) {
        Remove-Item env:VIRTUAL_ENV -ErrorAction SilentlyContinue
    }

    if ($env:VIRTUAL_ENV_PROMPT) {
        Remove-Item env:VIRTUAL_ENV_PROMPT -ErrorAction SilentlyContinue
    }

    if (Get-Variable -Name "_PYTHON_VENV_PROMPT_PREFIX" -ErrorAction SilentlyContinue) {
        Remove-Variable -Name _PYTHON_VENV_PROMPT_PREFIX -Scope Global -Force
    }

    if (!$NonDestructive) {
        Remove-Item function:deactivate
        Remove-Item function:pydoc
    }
}

function global:pydoc {
    & python -m pydoc $args
}

$VenvExecPath = Split-Path -Parent $MyInvocation.MyCommand.Definition
$VenvExecDir = Get-Item -Path $VenvExecPath

Write-Verbose "Activation script is located in path: '$VenvExecPath'"

if ($VenvDir) {
    Write-Verbose "VenvDir given as parameter, using '$VenvDir' to determine values"
} else {
    $VenvDir = $VenvExecDir.Parent.FullName.TrimEnd("\\/")
    Write-Verbose "VenvDir=$VenvDir"
}

$pyvenvCfg = Get-PyVenvConfig -ConfigDir $VenvDir

if ($Prompt) {
    Write-Verbose "Prompt specified as argument, using '$Prompt'"
} else {
    if ($pyvenvCfg -and $pyvenvCfg['prompt']) {
        Write-Verbose "Setting based on value in pyvenv.cfg='$($pyvenvCfg['prompt'])'"
        $Prompt = $pyvenvCfg['prompt']
    } elseif (__VIRTUAL_PROMPT__ -ne "") {
        $Prompt = __VIRTUAL_PROMPT__
    } else {
        $Prompt = Split-Path -Path $VenvDir -Leaf
    }
}

deactivate -nondestructive

$env:VIRTUAL_ENV = $VenvDir
$env:VIRTUAL_ENV_PROMPT = $Prompt

if (__TCL_LIBRARY__ -ne "") {
    if (Test-Path env:TCL_LIBRARY) {
        New-Variable -Scope global -Name _OLD_VIRTUAL_TCL_LIBRARY -Value $env:TCL_LIBRARY
    }
    $env:TCL_LIBRARY = __TCL_LIBRARY__
}

if (__TK_LIBRARY__ -ne "") {
    if (Test-Path env:TK_LIBRARY) {
        New-Variable -Scope global -Name _OLD_VIRTUAL_TK_LIBRARY -Value $env:TK_LIBRARY
    }
    $env:TK_LIBRARY = __TK_LIBRARY__
}

New-Variable -Scope global -Name _OLD_VIRTUAL_PATH -Value $env:PATH

if (Test-Path env:PKG_CONFIG_PATH) {
    New-Variable -Scope global -Name _OLD_PKG_CONFIG_PATH -Value $env:PKG_CONFIG_PATH
}
$env:PKG_CONFIG_PATH = "$env:VIRTUAL_ENV\lib\pkgconfig;$env:PKG_CONFIG_PATH"

$env:PATH = "$env:VIRTUAL_ENV/" + __BIN_NAME__ + __PATH_SEP__ + $env:PATH

if (!$env:VIRTUAL_ENV_DISABLE_PROMPT) {
    function global:_old_virtual_prompt {
        ""
    }
    $function:_old_virtual_prompt = $function:prompt
    New-Variable -Name _PYTHON_VENV_PROMPT_PREFIX -Description "Python virtual environment prompt prefix" -Scope Global -Option ReadOnly -Visibility Public -Value $Prompt

    function global:prompt {
        Write-Host -NoNewline -ForegroundColor Green "($($_PYTHON_VENV_PROMPT_PREFIX)) "
        _old_virtual_prompt
    }
}
