import pytest


@pytest.fixture(autouse=True)
def clear_font_filehandles(request):
    # https://github.com/matplotlib/matplotlib/issues/22017#issuecomment-998241017
    yield
    import matplotlib

    matplotlib.font_manager._get_font.cache_clear()
