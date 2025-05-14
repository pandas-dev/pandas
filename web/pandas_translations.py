#!/usr/bin/env python3
"""
Utilities to download and extract translations from the GitHub repository.
"""

import io
import os
import shutil
from subprocess import (
    PIPE,
    Popen,
)
import sys
import tarfile

import requests
import yaml


def get_config(config_fname: str) -> dict:
    """
    Load the config yaml file and return it as a dictionary.
    """
    with open(config_fname, encoding="utf-8") as f:
        context = yaml.safe_load(f)
    return context


def download_and_extract_translations(url: str, dir_name: str) -> None:
    """
    Download the translations from the GitHub repository.
    """
    shutil.rmtree(dir_name, ignore_errors=True)
    response = requests.get(url)
    if response.status_code == 200:
        doc = io.BytesIO(response.content)
        with tarfile.open(None, "r:gz", doc) as tar:
            tar.extractall(dir_name)
    else:
        raise Exception(f"Failed to download translations: {response.status_code}")


def get_languages(source_path: str) -> list[str]:
    """
    Get the list of languages available in the translations directory.
    """
    en_path = f"{source_path}/en/"
    if os.path.exists(en_path):
        shutil.rmtree(en_path)

    paths = os.listdir(source_path)
    return [path for path in paths if os.path.isdir(f"{source_path}/{path}")]


def remove_translations(source_path: str, languages: list[str]) -> None:
    """
    Remove the translations from the source path.
    """
    for language in languages:
        shutil.rmtree(os.path.join(source_path, language), ignore_errors=True)


def copy_translations(source_path: str, target_path: str, languages: list[str]) -> None:
    """
    Copy the translations to the appropriate directory.
    """
    for lang in languages:
        dest = f"{target_path}/{lang}/"
        shutil.rmtree(dest, ignore_errors=True)
        cmds = [
            "rsync",
            "-av",
            "--delete",
            f"{source_path}/{lang}/",
            dest,
        ]
        p = Popen(cmds, stdout=PIPE, stderr=PIPE)
        stdout, stderr = p.communicate()
        sys.stderr.write(f"\nCopying: {lang}...\n\n")
        sys.stderr.write(stdout.decode())
        sys.stderr.write(stderr.decode())


def process_translations(
    config_fname: str, source_path: str
) -> tuple[list[str], list[str]]:
    """
    Process the translations by downloading and extracting them from
    the GitHub repository.
    """
    base_folder = os.path.dirname(__file__)
    config = get_config(os.path.join(source_path, config_fname))
    translations_path = os.path.join(base_folder, f"{config['translations']['folder']}")
    translations_source_path = os.path.join(
        translations_path, config["translations"]["source_path"]
    )
    default_language = config["translations"]["default_language"]

    sys.stderr.write("\nDownloading and extracting translations...\n\n")
    download_and_extract_translations(config["translations"]["url"], translations_path)

    translated_languages = get_languages(translations_source_path)
    remove_translations(source_path, translated_languages)

    languages = [default_language] + translated_languages
    sys.stderr.write("\nCopying translations...\n")
    copy_translations(translations_source_path, source_path, translated_languages)

    return translated_languages, languages
