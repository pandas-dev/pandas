"""French search language."""

from __future__ import annotations

import snowballstemmer

from sphinx.search import SearchLanguage
from sphinx.search._stopwords.fr import FRENCH_STOPWORDS


class SearchFrench(SearchLanguage):
    lang = 'fr'
    language_name = 'French'
    js_stemmer_rawcode = 'french-stemmer.js'
    stopwords = FRENCH_STOPWORDS

    def __init__(self, options: dict[str, str]) -> None:
        super().__init__(options)
        self.stemmer = snowballstemmer.stemmer('french')

    def stem(self, word: str) -> str:
        return self.stemmer.stemWord(word.lower())
