"""German search language."""

from __future__ import annotations

import snowballstemmer

from sphinx.search import SearchLanguage
from sphinx.search._stopwords.de import GERMAN_STOPWORDS


class SearchGerman(SearchLanguage):
    lang = 'de'
    language_name = 'German'
    js_stemmer_rawcode = 'german-stemmer.js'
    stopwords = GERMAN_STOPWORDS

    def __init__(self, options: dict[str, str]) -> None:
        super().__init__(options)
        self.stemmer = snowballstemmer.stemmer('german')

    def stem(self, word: str) -> str:
        return self.stemmer.stemWord(word.lower())
