import pytest

from ..sync_flake8_versions import get_revisions


def test_wrong_yesqa(capsys):
    precommit_config = {
        "repos": [
            {
                "repo": "https://gitlab.com/pycqa/flake8",
                "rev": "0.1.1",
                "hooks": [
                    {
                        "id": "flake8",
                        "additional_dependencies": ["pandas-dev-flaker==0.2.0"],
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
                            "pandas-dev-flaker==0.2.1",
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


def test_wrong_env(capsys):
    precommit_config = {
        "repos": [
            {
                "repo": "https://gitlab.com/pycqa/flake8",
                "rev": "0.1.1",
                "hooks": [
                    {
                        "id": "flake8",
                        "additional_dependencies": ["pandas-dev-flaker==0.2.0"],
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
                            "pandas-dev-flaker==0.2.1",
                        ],
                    }
                ],
            },
        ]
    }
    environment = {
        "dependencies": [
            "flake8=1.5.6",
            {
                "pip": [
                    "git+https://github.com/pydata/pydata-sphinx-theme.git@master",
                    "pandas-dev-flaker>=0.2.2",
                ]
            },
        ]
    }
    with pytest.raises(SystemExit, match=None):
        get_revisions(precommit_config, environment)
    result, _ = capsys.readouterr()
    expected = (
        "flake8 in 'environment.yml' does not match in 'flake8' from 'pre-commit'\n"
    )
    assert result == expected


def test_get_revisions_yesqa_failure(capsys):
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
                            "flake8==0.4.2",
                            "pandas-dev-flaker==0.2.1",
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
            {
                "pip": [
                    "pandas-dev-flaker>=0.2.2",
                ]
            },
        ]
    }
    with pytest.raises(SystemExit, match=None):
        get_revisions(precommit_config, environment)
    result, _ = capsys.readouterr()
    expected = (
        "Additional depedency of 'flake8' 'flake8-bugs' does not match in 'yesqa'"
    )
    assert result == expected


def test_get_revisions_env_failure(capsys):
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
                            "flake8==0.4.2",
                            "pandas-dev-flaker==0.2.1",
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
            {
                "pip": [
                    "pandas-dev-flaker>=0.2.3",
                ]
            },
        ]
    }
    with pytest.raises(SystemExit, match=None):
        get_revisions(precommit_config, environment)
    result, _ = capsys.readouterr()
    expected = (
        "Additional depedency of 'flake8' 'flake8-bugs' does not "
        "match in 'environment.yml'"
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
                    "pandas-dev-flaker==0.2.0",
                ]
            },
        ]
    }
    # should not raise
    get_revisions(precommit_config, environment)
