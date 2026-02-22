"""Portuguese search language."""

from __future__ import annotations

import snowballstemmer

from sphinx.search import SearchLanguage
from sphinx.search._stopwords.pt import PORTUGUESE_STOPWORDS


class SearchPortuguese(SearchLanguage):
    lang = 'pt'
    language_name = 'Portuguese'
    js_stemmer_rawcode = 'portuguese-stemmer.js'
    stopwords = PORTUGUESE_STOPWORDS

    def __init__(self, options: dict[str, str]) -> None:
        super().__init__(options)
        self.stemmer = snowballstemmer.stemmer('portuguese')

    def stem(self, word: str) -> str:
        return self.stemmer.stemWord(word.lower())
