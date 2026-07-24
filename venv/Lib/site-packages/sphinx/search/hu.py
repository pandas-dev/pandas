"""Hungarian search language."""

from __future__ import annotations

import snowballstemmer

from sphinx.search import SearchLanguage
from sphinx.search._stopwords.hu import HUNGARIAN_STOPWORDS


class SearchHungarian(SearchLanguage):
    lang = 'hu'
    language_name = 'Hungarian'
    js_stemmer_rawcode = 'hungarian-stemmer.js'
    stopwords = HUNGARIAN_STOPWORDS

    def __init__(self, options: dict[str, str]) -> None:
        super().__init__(options)
        self.stemmer = snowballstemmer.stemmer('hungarian')

    def stem(self, word: str) -> str:
        return self.stemmer.stemWord(word.lower())
