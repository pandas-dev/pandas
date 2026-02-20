"""Swedish search language."""

from __future__ import annotations

import snowballstemmer

from sphinx.search import SearchLanguage
from sphinx.search._stopwords.sv import SWEDISH_STOPWORDS


class SearchSwedish(SearchLanguage):
    lang = 'sv'
    language_name = 'Swedish'
    js_stemmer_rawcode = 'swedish-stemmer.js'
    stopwords = SWEDISH_STOPWORDS

    def __init__(self, options: dict[str, str]) -> None:
        super().__init__(options)
        self.stemmer = snowballstemmer.stemmer('swedish')

    def stem(self, word: str) -> str:
        return self.stemmer.stemWord(word.lower())
