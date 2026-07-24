"""Turkish search language."""

from __future__ import annotations

import snowballstemmer

from sphinx.search import SearchLanguage


class SearchTurkish(SearchLanguage):
    lang = 'tr'
    language_name = 'Turkish'
    js_stemmer_rawcode = 'turkish-stemmer.js'
    stopwords = frozenset()

    def __init__(self, options: dict[str, str]) -> None:
        super().__init__(options)
        self.stemmer = snowballstemmer.stemmer('turkish')

    def stem(self, word: str) -> str:
        return self.stemmer.stemWord(word.lower())
