from ..sync_flake8_versions import (
    Revisions,
    get_revisions,
)


def test_get_revisions():
    precommit_config = {
        "repos": [
            {
                "repo": "https://gitlab.com/pycqa/flake8",
                "rev": "foo",
                "hooks": [{"id": "flake8"}],
            },
            {
                "repo": "https://github.com/asottile/yesqa",
                "rev": "v1.2.2",
                "hooks": [{"id": "yesqa", "additional_dependencies": ["flake8==bar"]}],
            },
        ]
    }
    environment = {"dependencies": ["flake8=qux"]}
    result = get_revisions(precommit_config, environment)
    expected = Revisions("foo", "bar", "qux")
    assert result == expected
