"""Dutch search language."""

from __future__ import annotations

import snowballstemmer

from sphinx.search import SearchLanguage
from sphinx.search._stopwords.nl import DUTCH_STOPWORDS


class SearchDutch(SearchLanguage):
    lang = 'nl'
    language_name = 'Dutch'
    js_stemmer_rawcode = 'dutch-stemmer.js'
    stopwords = DUTCH_STOPWORDS

    def __init__(self, options: dict[str, str]) -> None:
        super().__init__(options)
        self.stemmer = snowballstemmer.stemmer('dutch')

    def stem(self, word: str) -> str:
        return self.stemmer.stemWord(word.lower())
