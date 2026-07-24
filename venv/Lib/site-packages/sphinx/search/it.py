"""Italian search language."""

from __future__ import annotations

import snowballstemmer

from sphinx.search import SearchLanguage
from sphinx.search._stopwords.it import ITALIAN_STOPWORDS


class SearchItalian(SearchLanguage):
    lang = 'it'
    language_name = 'Italian'
    js_stemmer_rawcode = 'italian-stemmer.js'
    stopwords = ITALIAN_STOPWORDS

    def __init__(self, options: dict[str, str]) -> None:
        super().__init__(options)
        self.stemmer = snowballstemmer.stemmer('italian')

    def stem(self, word: str) -> str:
        return self.stemmer.stemWord(word.lower())
