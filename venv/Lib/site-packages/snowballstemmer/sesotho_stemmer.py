# Generated from sesotho.sbl by Snowball 3.1.1 - https://snowballstem.org/

from .basestemmer import BaseStemmer
from .among import Among


class SesothoStemmer(BaseStemmer):
    '''
    This class implements the stemming algorithm defined by a snowball script.
    Generated from sesotho.sbl by Snowball 3.1.1 - https://snowballstem.org/
    '''

    g_v = {"a", "e", "i", "o", "u"}

    I_pV = 0

    def __r_mark_regions(self):
        v_1 = self.cursor
        if not self.go_out_grouping(SesothoStemmer.g_v):
            return False
        self.cursor += 1
        self.I_pV = self.cursor
        self.cursor = v_1
        v_2 = self.cursor
        if self.cursor + 2 > self.limit:
            return False
        self.cursor += 2
        try:
            if self.cursor <= self.I_pV:
                raise lab0()
            self.I_pV = self.cursor
        except lab0: pass
        self.cursor = v_2
        return True

    def __r_remove_noun_prefixes(self):
        self.bra = self.cursor
        if self.find_among(SesothoStemmer.a_0) == 0:
            return False
        self.ket = self.cursor
        v_1 = self.cursor
        if self.cursor >= self.limit:
            return False
        self.cursor += 1
        if self.cursor >= self.limit:
            return False
        self.cursor = v_1
        if not self.go_out_grouping(SesothoStemmer.g_v):
            return False
        self.cursor += 1
        self.slice_del()
        return True

    def __r_remove_verb_suffixes(self):
        if self.cursor < self.I_pV:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_pV
        self.ket = self.cursor
        if self.find_among_b(SesothoStemmer.a_1) == 0:
            self.limit_backward = v_2
            return False
        self.bra = self.cursor
        self.slice_del()
        self.limit_backward = v_2
        return True

    def __r_remove_nominal_suffixes(self):
        if self.cursor < self.I_pV:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_pV
        self.ket = self.cursor
        if self.find_among_b(SesothoStemmer.a_2) == 0:
            self.limit_backward = v_2
            return False
        self.bra = self.cursor
        self.slice_del()
        self.limit_backward = v_2
        return True

    def _stem(self):
        if not self.__r_mark_regions():
            return False
        self.limit_backward = self.cursor
        self.cursor = self.limit
        v_1 = self.limit - self.cursor
        self.__r_remove_nominal_suffixes()
        self.cursor = self.limit - v_1
        v_2 = self.limit - self.cursor
        self.__r_remove_verb_suffixes()
        self.cursor = self.limit - v_2
        self.cursor = self.limit_backward
        v_3 = self.cursor
        self.__r_remove_noun_prefixes()
        self.cursor = v_3
        return True

    a_0 = [
        Among("ba", -1, -1),
        Among("boi", -1, -1),
        Among("le", -1, -1),
        Among("li", -1, -1),
        Among("ma", -1, -1),
        Among("me", -1, -1),
        Among("mo", -1, -1),
        Among("se", -1, -1)
    ]

    a_1 = [
        Among("a", -1, 1),
        Among("ela", 0, 1),
        Among("isa", 0, 1),
        Among("wa", 0, 1),
        Among("ile", -1, 1),
        Among("etse", -1, 1),
        Among("ang", -1, 1),
        Among("eng", -1, 1),
        Among("ong", -1, 1)
    ]

    a_2 = [
        Among("ana", -1, 1),
        Among("nyana", 0, 1),
        Among("oa", -1, 1),
        Among("i", -1, 1),
        Among("ano", -1, 1)
    ]


class lab0(BaseException): pass
