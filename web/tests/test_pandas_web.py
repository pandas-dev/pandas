import yaml
import os
# Not working in code. the pandas directory is set up for this, but not the web directory
from pandas_web import Preprocessors
import pytest

@pytest.fixture # I think this is fine
def test_home_releases_versions():
    context_path = os.path.join('pandas/config.yml')

    # As it is, the home_add_releases function pull data from the github release page
    # with an http request. After I run the function, I print the data out.
    # In order to run custom data, you would need to modify the function in some way
    with open(context_path, 'r') as context_file:
        context = yaml.safe_load(context_file)

        context["target_path"] = "build"

        Preprocessors.home_add_releases(context)

        for release in context["releases"]:
            print(release)