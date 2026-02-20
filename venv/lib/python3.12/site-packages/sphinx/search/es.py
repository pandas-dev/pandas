"""Spanish search language."""

from __future__ import annotations

import snowballstemmer

from sphinx.search import SearchLanguage
from sphinx.search._stopwords.es import SPANISH_STOPWORDS


class SearchSpanish(SearchLanguage):
    lang = 'es'
    language_name = 'Spanish'
    js_stemmer_rawcode = 'spanish-stemmer.js'
    stopwords = SPANISH_STOPWORDS

    def __init__(self, options: dict[str, str]) -> None:
        super().__init__(options)
        self.stemmer = snowballstemmer.stemmer('spanish')

    def stem(self, word: str) -> str:
        return self.stemmer.stemWord(word.lower())
