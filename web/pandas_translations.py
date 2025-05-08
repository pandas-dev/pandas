import os
from pathlib import Path
from subprocess import (
    PIPE,
    Popen,
)
import tarfile

import requests
import yaml


def download_translations(url, fname):
    """
    Download the translations from the GitHub repository.
    """
    response = requests.get(url)
    if response.status_code == 200:
        with open(fname, "wb") as f:
            f.write(response.content)
    else:
        raise Exception(f"Failed to download translations: {response.status_code}")


def extract_translations(fpath, dir_name):
    """
    Extract the translations from the tar file.
    """
    with tarfile.open(fpath, "r:gz") as tar:
        tar.extractall(dir_name)
    print(f"Translations extracted to '{dir_name}' directory.")


def load_translations_config(path):
    """
    Load the translations configuration from the YAML file.
    """
    with open(path) as f:
        config = yaml.safe_load(f)
    return config


def load_status_config(path):
    """
    Load the translations configuration from the YAML file.
    """
    with open(path) as f:
        config = yaml.safe_load(f)
    return config


def copy_translations(data, translation_percentage, dir_name, dest_dir):
    """
    Copy the translations to the appropriate directory.
    """
    for lang, value in data.items():
        if value["progress"] >= translation_percentage:
            language_code = lang[:2]
            src_path = (
                f"{dir_name}/pandas-translations-main/web/pandas/{language_code}/"
            )
            dest_path = f"{dest_dir}/{language_code}/"
            cmds = ["rsync", "-av", "--delete", src_path, dest_path]
            p = Popen(cmds, stdout=PIPE, stderr=PIPE)
            stdout, stderr = p.communicate()
            print(stdout.decode())
            print(stderr.decode())


if __name__ == "__main__":
    os.chdir(Path(__file__).parent.parent)
    url = "https://github.com/Scientific-Python-Translations/pandas-translations/archive/refs/heads/main.tar.gz"
    fpath = "web/pandas-translations.tar.gz"
    dir_name = "web/translations"
    dest_dir = "web/pandas"
    config_path = (
        f"{dir_name}/pandas-translations-main/.github/workflows/sync_translations.yml"
    )
    status_path = f"{dir_name}/pandas-translations-main/status.yml"

    download_translations(url, fpath)
    extract_translations(fpath, dir_name)
    config = load_translations_config(config_path)
    variables = config["jobs"]["sync_translations"]["steps"][0]["with"]
    translation_percentage = int(variables["translation-percentage"])
    status = load_status_config(status_path)
    copy_translations(
        status, translation_percentage, dir_name=dir_name, dest_dir=dest_dir
    )
