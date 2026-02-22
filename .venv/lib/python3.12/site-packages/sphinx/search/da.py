"""Danish search language."""

from __future__ import annotations

import snowballstemmer

from sphinx.search import SearchLanguage
from sphinx.search._stopwords.da import DANISH_STOPWORDS


class SearchDanish(SearchLanguage):
    lang = 'da'
    language_name = 'Danish'
    js_stemmer_rawcode = 'danish-stemmer.js'
    stopwords = DANISH_STOPWORDS

    def __init__(self, options: dict[str, str]) -> None:
        super().__init__(options)
        self.stemmer = snowballstemmer.stemmer('danish')

    def stem(self, word: str) -> str:
        return self.stemmer.stemWord(word.lower())
