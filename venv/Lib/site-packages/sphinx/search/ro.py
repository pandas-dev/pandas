"""Romanian search language."""

from __future__ import annotations

import snowballstemmer

from sphinx.search import SearchLanguage


class SearchRomanian(SearchLanguage):
    lang = 'ro'
    language_name = 'Romanian'
    js_stemmer_rawcode = 'romanian-stemmer.js'
    stopwords = frozenset()

    def __init__(self, options: dict[str, str]) -> None:
        super().__init__(options)
        self.stemmer = snowballstemmer.stemmer('romanian')

    def stem(self, word: str) -> str:
        return self.stemmer.stemWord(word.lower())
