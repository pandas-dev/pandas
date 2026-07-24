"""Russian search language."""

from __future__ import annotations

import snowballstemmer

from sphinx.search import SearchLanguage
from sphinx.search._stopwords.ru import RUSSIAN_STOPWORDS


class SearchRussian(SearchLanguage):
    lang = 'ru'
    language_name = 'Russian'
    js_stemmer_rawcode = 'russian-stemmer.js'
    stopwords = RUSSIAN_STOPWORDS

    def __init__(self, options: dict[str, str]) -> None:
        super().__init__(options)
        self.stemmer = snowballstemmer.stemmer('russian')

    def stem(self, word: str) -> str:
        return self.stemmer.stemWord(word.lower())
