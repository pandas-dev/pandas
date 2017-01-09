# Sample script to install Miniconda under Windows
# Authors: Olivier Grisel, Jonathan Helmus and Kyle Kastner, Robert McGibbon
# License: CC0 1.0 Universal: http://creativecommons.org/publicdomain/zero/1.0/

$MINICONDA_URL = "http://repo.continuum.io/miniconda/"


function DownloadMiniconda ($python_version, $platform_suffix) {
    $webclient = New-Object System.Net.WebClient
    $filename = "Miniconda3-latest-Windows-" + $platform_suffix + ".exe"
    $url = $MINICONDA_URL + $filename

    $basedir = $pwd.Path + "\"
    $filepath = $basedir + $filename
    if (Test-Path $filename) {
        Write-Host "Reusing" $filepath
        return $filepath
    }

    # Download and retry up to 3 times in case of network transient errors.
    Write-Host "Downloading" $filename "from" $url
    $retry_attempts = 2
    for($i=0; $i -lt $retry_attempts; $i++){
        try {
            $webclient.DownloadFile($url, $filepath)
            break
        }
        Catch [Exception]{
            Start-Sleep 1
        }
   }
   if (Test-Path $filepath) {
       Write-Host "File saved at" $filepath
   } else {
       # Retry once to get the error message if any at the last try
       $webclient.DownloadFile($url, $filepath)
   }
   return $filepath
}


function InstallMiniconda ($python_version, $architecture, $python_home) {
    Write-Host "Installing Python" $python_version "for" $architecture "bit architecture to" $python_home
    if (Test-Path $python_home) {
        Write-Host $python_home "already exists, skipping."
        return $false
    }
    if ($architecture -match "32") {
        $platform_suffix = "x86"
    } else {
        $platform_suffix = "x86_64"
    }

    $filepath = DownloadMiniconda $python_version $platform_suffix
    Write-Host "Installing" $filepath "to" $python_home
    $install_log = $python_home + ".log"
    $args = "/S /D=$python_home"
    Write-Host $filepath $args
    Start-Process -FilePath $filepath -ArgumentList $args -Wait -Passthru
    if (Test-Path $python_home) {
        Write-Host "Python $python_version ($architecture) installation complete"
    } else {
        Write-Host "Failed to install Python in $python_home"
        Get-Content -Path $install_log
        Exit 1
    }
}


function InstallCondaPackages ($python_home, $spec) {
    $conda_path = $python_home + "\Scripts\conda.exe"
    $args = "install --yes " + $spec
    Write-Host ("conda " + $args)
    Start-Process -FilePath "$conda_path" -ArgumentList $args -Wait -Passthru
}

function UpdateConda ($python_home) {
    $conda_path = $python_home + "\Scripts\conda.exe"
    Write-Host "Updating conda..."
    $args = "update --yes conda"
    Write-Host $conda_path $args
    Start-Process -FilePath "$conda_path" -ArgumentList $args -Wait -Passthru
}


function main () {
    InstallMiniconda "3.5" $env:PYTHON_ARCH $env:CONDA_ROOT
    UpdateConda $env:CONDA_ROOT
    InstallCondaPackages $env:CONDA_ROOT "conda-build jinja2 anaconda-client"
}

main
