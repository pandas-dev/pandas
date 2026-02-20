"""Chinese search language: includes routine to split words."""

from __future__ import annotations

import re
from pathlib import Path
from typing import TYPE_CHECKING

import snowballstemmer

from sphinx.search import SearchLanguage
from sphinx.search._stopwords.en import ENGLISH_STOPWORDS

if TYPE_CHECKING:
    from collections.abc import Iterator

try:
    import jieba  # type: ignore[import-not-found]

    jieba_load_userdict = jieba.load_userdict
    cut_for_search = jieba.cut_for_search
except ImportError:
    JIEBA_DEFAULT_DICT = ''

    def jieba_load_userdict(f: str) -> None:
        pass

    def cut_for_search(sentence: str, HMM: bool = True) -> Iterator[str]:
        yield from ()

else:
    JIEBA_DEFAULT_DICT = (
        Path(jieba.__file__, '..', jieba.DEFAULT_DICT_NAME).resolve().as_posix()
    )
    del jieba


class SearchChinese(SearchLanguage):
    """Chinese search implementation"""

    lang = 'zh'
    language_name = 'Chinese'
    js_stemmer_rawcode = 'english-stemmer.js'
    stopwords = ENGLISH_STOPWORDS
    latin1_letters = re.compile(r'[a-zA-Z0-9_]+')

    def __init__(self, options: dict[str, str]) -> None:
        super().__init__(options)
        self.latin_terms: set[str] = set()
        dict_path = options.get('dict', JIEBA_DEFAULT_DICT)
        if dict_path and Path(dict_path).is_file():
            jieba_load_userdict(str(dict_path))

        self.stemmer = snowballstemmer.stemmer('english')

    def split(self, input: str) -> list[str]:
        chinese: list[str] = list(cut_for_search(input))

        latin1 = [term.strip() for term in self.latin1_letters.findall(input)]
        self.latin_terms.update(latin1)
        return chinese + latin1

    def word_filter(self, stemmed_word: str) -> bool:
        return len(stemmed_word) > 1

    def stem(self, word: str) -> str:
        # Don't stem Latin words that are long enough to be relevant for search
        # if not stemmed, but would be too short after being stemmed
        # avoids some issues with acronyms
        stemmed = self.stemmer.stemWord(word.lower())
        should_not_be_stemmed = (
            len(word) >= 3 > len(stemmed) and word in self.latin_terms
        )
        if should_not_be_stemmed:
            return word.lower()
        return stemmed
