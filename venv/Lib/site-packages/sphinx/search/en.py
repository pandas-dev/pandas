"""English search language."""

from __future__ import annotations

import snowballstemmer

from sphinx.search import SearchLanguage
from sphinx.search._stopwords.en import ENGLISH_STOPWORDS


class SearchEnglish(SearchLanguage):
    lang = 'en'
    language_name = 'English'
    js_stemmer_rawcode = 'english-stemmer.js'
    stopwords = ENGLISH_STOPWORDS

    def __init__(self, options: dict[str, str]) -> None:
        super().__init__(options)
        self.stemmer = snowballstemmer.stemmer('english')

    def stem(self, word: str) -> str:
        return self.stemmer.stemWord(word.lower())
