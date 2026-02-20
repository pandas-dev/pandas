"""Finnish search language."""

from __future__ import annotations

import snowballstemmer

from sphinx.search import SearchLanguage
from sphinx.search._stopwords.fi import FINNISH_STOPWORDS


class SearchFinnish(SearchLanguage):
    lang = 'fi'
    language_name = 'Finnish'
    js_stemmer_rawcode = 'finnish-stemmer.js'
    stopwords = FINNISH_STOPWORDS

    def __init__(self, options: dict[str, str]) -> None:
        super().__init__(options)
        self.stemmer = snowballstemmer.stemmer('finnish')

    def stem(self, word: str) -> str:
        return self.stemmer.stemWord(word.lower())
