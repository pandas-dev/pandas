# Generated from norwegian.sbl by Snowball 3.1.1 - https://snowballstem.org/

from .basestemmer import BaseStemmer
from .among import Among


class NorwegianStemmer(BaseStemmer):
    '''
    This class implements the stemming algorithm defined by a snowball script.
    Generated from norwegian.sbl by Snowball 3.1.1 - https://snowballstem.org/
    '''

    g_v = {"a", "e", "i", "o", "u", "y", "å", "æ", "ê", "ò", "ó", "ô", "ø"}

    g_s_ending = {"b", "c", "d", "f", "g", "h", "j", "l", "m", "n", "o", "p", "t", "v", "y", "z"}

    I_p1 = 0

    def __r_mark_regions(self):
        self.I_p1 = self.limit
        v_1 = self.cursor
        try:
            while True:
                v_2 = self.cursor
                try:
                    while True:
                        try:
                            if self.cursor == self.limit or self.current[self.cursor] != "'":
                                raise lab2()
                            self.cursor += 1
                            break
                        except lab2: pass
                        if self.cursor >= self.limit:
                            raise lab1()
                        self.cursor += 1
                    break
                except lab1: pass
                self.cursor = v_2
                if not self.go_out_grouping(NorwegianStemmer.g_v):
                    raise lab0()
                self.cursor += 1
                if not self.go_in_grouping(NorwegianStemmer.g_v):
                    raise lab0()
                self.cursor += 1
                break
            self.I_p1 = self.cursor
        except lab0: pass
        self.cursor = v_1
        v_3 = self.cursor
        if self.cursor + 3 > self.limit:
            return False
        self.cursor += 3
        try:
            if self.I_p1 >= self.cursor:
                raise lab0()
            self.I_p1 = self.cursor
        except lab0: pass
        self.cursor = v_3
        return True

    def __r_main_suffix(self):
        if self.cursor < self.I_p1:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_p1
        self.ket = self.cursor
        among_var = self.find_among_b(NorwegianStemmer.a_1)
        if among_var == 0:
            self.limit_backward = v_2
            return False
        self.bra = self.cursor
        self.limit_backward = v_2
        if among_var == 1:
            self.slice_del()
        elif among_var == 2:
            among_var = self.find_among_b(NorwegianStemmer.a_0)
            if among_var == 1:
                self.slice_del()
        elif among_var == 3:
            while True:
                v_3 = self.limit - self.cursor
                try:
                    if not self.in_grouping_b(NorwegianStemmer.g_s_ending):
                        raise lab0()
                    break
                except lab0: pass
                self.cursor = self.limit - v_3
                try:
                    if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "r":
                        raise lab0()
                    self.cursor -= 1
                    try:
                        if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "e":
                            raise lab1()
                        self.cursor -= 1
                        raise lab0()
                    except lab1: pass
                    break
                except lab0: pass
                self.cursor = self.limit - v_3
                if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "k":
                    return False
                self.cursor -= 1
                if not self.out_grouping_b(NorwegianStemmer.g_v):
                    return False
                break
            self.slice_del()
        else:
            self.slice_from("er")
        return True

    def __r_consonant_pair(self):
        v_1 = self.limit - self.cursor
        if self.cursor < self.I_p1:
            return False
        v_3 = self.limit_backward
        self.limit_backward = self.I_p1
        self.ket = self.cursor
        if self.find_among_b(NorwegianStemmer.a_2) == 0:
            self.limit_backward = v_3
            return False
        self.bra = self.cursor
        self.limit_backward = v_3
        self.cursor = self.limit - v_1
        if self.cursor <= self.limit_backward:
            return False
        self.cursor -= 1
        self.bra = self.cursor
        self.slice_del()
        return True

    def __r_other_suffix(self):
        if self.cursor < self.I_p1:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_p1
        self.ket = self.cursor
        if self.find_among_b(NorwegianStemmer.a_3) == 0:
            self.limit_backward = v_2
            return False
        self.bra = self.cursor
        self.limit_backward = v_2
        self.slice_del()
        return True

    def _stem(self):
        if not self.__r_mark_regions():
            return False
        self.limit_backward = self.cursor
        self.cursor = self.limit
        v_1 = self.limit - self.cursor
        self.__r_main_suffix()
        self.cursor = self.limit - v_1
        v_2 = self.limit - self.cursor
        self.__r_consonant_pair()
        self.cursor = self.limit - v_2
        v_3 = self.limit - self.cursor
        self.__r_other_suffix()
        self.cursor = self.limit - v_3
        self.ket = self.cursor
        if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "'":
            return False
        self.cursor -= 1
        self.bra = self.cursor
        self.slice_del()
        self.cursor = self.limit_backward
        return True

    a_0 = [
        Among("", -1, 1),
        Among("ind", 0, -1),
        Among("kk", 0, -1),
        Among("nk", 0, -1),
        Among("amm", 0, -1),
        Among("omm", 0, -1),
        Among("kap", 0, -1),
        Among("skap", 6, 1),
        Among("pp", 0, -1),
        Among("lt", 0, -1),
        Among("ast", 0, -1),
        Among("øst", 0, -1),
        Among("v", 0, -1),
        Among("hav", 12, 1),
        Among("giv", 12, 1)
    ]

    a_1 = [
        Among("a", -1, 1),
        Among("e", -1, 1),
        Among("ede", 1, 1),
        Among("ande", 1, 1),
        Among("ende", 1, 1),
        Among("ane", 1, 1),
        Among("ene", 1, 1),
        Among("hetene", 6, 1),
        Among("erte", 1, 4),
        Among("en", -1, 1),
        Among("heten", 9, 1),
        Among("ar", -1, 1),
        Among("er", -1, 1),
        Among("heter", 12, 1),
        Among("s", -1, 3),
        Among("as", 14, 1),
        Among("es", 14, 1),
        Among("edes", 16, 1),
        Among("endes", 16, 1),
        Among("enes", 16, 1),
        Among("hetenes", 19, 1),
        Among("ens", 14, 1),
        Among("hetens", 21, 1),
        Among("ers", 14, 2),
        Among("ets", 14, 1),
        Among("et", -1, 1),
        Among("het", 25, 1),
        Among("ert", -1, 4),
        Among("ast", -1, 1)
    ]

    a_2 = [
        Among("dt", -1, -1),
        Among("vt", -1, -1)
    ]

    a_3 = [
        Among("leg", -1, 1),
        Among("eleg", 0, 1),
        Among("ig", -1, 1),
        Among("eig", 2, 1),
        Among("lig", 2, 1),
        Among("elig", 4, 1),
        Among("els", -1, 1),
        Among("lov", -1, 1),
        Among("elov", 7, 1),
        Among("slov", 7, 1),
        Among("hetslov", 9, 1)
    ]


class lab0(BaseException): pass


class lab1(BaseException): pass


class lab2(BaseException): pass
