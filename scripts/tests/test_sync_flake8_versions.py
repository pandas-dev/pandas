import pytest

from ..sync_flake8_versions import get_revisions


def test_wrong_yesqa_flake8(capsys):
    precommit_config = {
        "repos": [
            {
                "repo": "https://gitlab.com/pycqa/flake8",
                "rev": "0.1.1",
                "hooks": [
                    {
                        "id": "flake8",
                    }
                ],
            },
            {
                "repo": "https://github.com/asottile/yesqa",
                "rev": "v1.2.2",
                "hooks": [
                    {
                        "id": "yesqa",
                        "additional_dependencies": [
                            "flake8==0.4.2",
                        ],
                    }
                ],
            },
        ]
    }
    environment = {
        "dependencies": [
            "flake8=0.1.1",
        ]
    }
    with pytest.raises(SystemExit, match=None):
        get_revisions(precommit_config, environment)
    result, _ = capsys.readouterr()
    expected = "flake8 in 'yesqa' does not match in 'flake8' from 'pre-commit'\n"
    assert result == expected


def test_wrong_env_flake8(capsys):
    precommit_config = {
        "repos": [
            {
                "repo": "https://gitlab.com/pycqa/flake8",
                "rev": "0.1.1",
                "hooks": [
                    {
                        "id": "flake8",
                    }
                ],
            },
            {
                "repo": "https://github.com/asottile/yesqa",
                "rev": "v1.2.2",
                "hooks": [
                    {
                        "id": "yesqa",
                        "additional_dependencies": [
                            "flake8==0.4.2",
                        ],
                    }
                ],
            },
        ]
    }
    environment = {
        "dependencies": [
            "flake8=1.5.6",
        ]
    }
    with pytest.raises(SystemExit, match=None):
        get_revisions(precommit_config, environment)
    result, _ = capsys.readouterr()
    expected = (
        "flake8 in 'environment.yml' does not match in 'flake8' from 'pre-commit'\n"
    )
    assert result == expected


def test_wrong_yesqa_add_dep(capsys):
    precommit_config = {
        "repos": [
            {
                "repo": "https://gitlab.com/pycqa/flake8",
                "rev": "0.1.1",
                "hooks": [
                    {
                        "id": "flake8",
                        "additional_dependencies": [
                            "flake8-bugs==1.1.1",
                        ],
                    }
                ],
            },
            {
                "repo": "https://github.com/asottile/yesqa",
                "rev": "v1.2.2",
                "hooks": [
                    {
                        "id": "yesqa",
                        "additional_dependencies": [
                            "flake8==0.4.2",
                            "flake8-bugs>=1.1.1",
                        ],
                    }
                ],
            },
        ]
    }
    environment = {
        "dependencies": [
            "flake8=1.5.6",
            "flake8-bugs=1.1.1",
        ]
    }
    with pytest.raises(SystemExit, match=None):
        get_revisions(precommit_config, environment)
    result, _ = capsys.readouterr()
    expected = (
        "Mismatch of 'flake8-bugs' version between 'flake8' and 'yesqa' in "
        "'.pre-commit-config.yaml'\n"
    )
    assert result == expected


def test_wrong_env_add_dep(capsys):
    precommit_config = {
        "repos": [
            {
                "repo": "https://gitlab.com/pycqa/flake8",
                "rev": "0.1.1",
                "hooks": [
                    {
                        "id": "flake8",
                        "additional_dependencies": [
                            "flake8-bugs==1.1.1",
                        ],
                    }
                ],
            },
            {
                "repo": "https://github.com/asottile/yesqa",
                "rev": "v1.2.2",
                "hooks": [
                    {
                        "id": "yesqa",
                        "additional_dependencies": [
                            "flake8==0.4.2",
                            "flake8-bugs==1.1.1",
                        ],
                    }
                ],
            },
        ]
    }
    environment = {
        "dependencies": [
            "flake8=1.5.6",
            "flake8-bugs=1.1.2",
        ]
    }
    with pytest.raises(SystemExit, match=None):
        get_revisions(precommit_config, environment)
    result, _ = capsys.readouterr()
    expected = (
        "Mismatch of 'flake8-bugs' version between 'enviroment.yml' "
        "and additional dependencies of 'flake8' in '.pre-commit-config.yaml'\n"
    )
    assert result == expected


def test_get_revisions_no_failure(capsys):
    precommit_config = {
        "repos": [
            {
                "repo": "https://gitlab.com/pycqa/flake8",
                "rev": "0.1.1",
                "hooks": [
                    {
                        "id": "flake8",
                        "additional_dependencies": [
                            "pandas-dev-flaker==0.2.0",
                            "flake8-bugs==1.1.1",
                        ],
                    }
                ],
            },
            {
                "repo": "https://github.com/asottile/yesqa",
                "rev": "v1.2.2",
                "hooks": [
                    {
                        "id": "yesqa",
                        "additional_dependencies": [
                            "flake8==0.1.1",
                            "pandas-dev-flaker==0.2.0",
                            "flake8-bugs==1.1.1",
                        ],
                    }
                ],
            },
        ]
    }
    environment = {
        "dependencies": [
            "flake8=0.1.1",
            "flake8-bugs=1.1.1",
            {
                "pip": [
                    "git+https://github.com/pydata/pydata-sphinx-theme.git@master",
                    "pandas-dev-flaker==0.2.0",
                ]
            },
        ]
    }
    # should not raise
    get_revisions(precommit_config, environment)
