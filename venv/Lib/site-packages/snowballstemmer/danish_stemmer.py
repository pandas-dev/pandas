# Generated from danish.sbl by Snowball 3.1.1 - https://snowballstem.org/

from .basestemmer import BaseStemmer
from .among import Among


class DanishStemmer(BaseStemmer):
    '''
    This class implements the stemming algorithm defined by a snowball script.
    Generated from danish.sbl by Snowball 3.1.1 - https://snowballstem.org/
    '''

    g_undouble_c = {"b", "d", "f", "g", "k", "l", "m", "n", "p", "r", "s", "t"}

    g_v = {"a", "e", "i", "o", "u", "y", "å", "æ", "ø"}

    g_s_ending = {"'", "a", "b", "c", "d", "f", "g", "h", "j", "k", "l", "m", "n", "o", "p", "r", "t", "v", "y", "z", "å"}

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
                if not self.go_out_grouping(DanishStemmer.g_v):
                    raise lab0()
                self.cursor += 1
                if not self.go_in_grouping(DanishStemmer.g_v):
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
        among_var = self.find_among_b(DanishStemmer.a_0)
        if among_var == 0:
            self.limit_backward = v_2
            return False
        self.bra = self.cursor
        self.limit_backward = v_2
        if among_var == 1:
            self.slice_del()
        else:
            if not self.in_grouping_b(DanishStemmer.g_s_ending):
                return False
            self.slice_del()
        return True

    def __r_consonant_pair(self):
        v_1 = self.limit - self.cursor
        if self.cursor < self.I_p1:
            return False
        v_3 = self.limit_backward
        self.limit_backward = self.I_p1
        self.ket = self.cursor
        if self.find_among_b(DanishStemmer.a_1) == 0:
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
        v_1 = self.limit - self.cursor
        try:
            self.ket = self.cursor
            if not self.eq_s_b("st"):
                raise lab0()
            self.bra = self.cursor
            if not self.eq_s_b("ig"):
                raise lab0()
            self.slice_del()
        except lab0: pass
        self.cursor = self.limit - v_1
        if self.cursor < self.I_p1:
            return False
        v_3 = self.limit_backward
        self.limit_backward = self.I_p1
        self.ket = self.cursor
        among_var = self.find_among_b(DanishStemmer.a_2)
        if among_var == 0:
            self.limit_backward = v_3
            return False
        self.bra = self.cursor
        self.limit_backward = v_3
        if among_var == 1:
            self.slice_del()
            v_4 = self.limit - self.cursor
            self.__r_consonant_pair()
            self.cursor = self.limit - v_4
        else:
            self.slice_from("løs")
        return True

    def __r_undouble(self):
        if self.cursor < self.I_p1:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_p1
        self.ket = self.cursor
        if not self.in_grouping_b(DanishStemmer.g_undouble_c):
            self.limit_backward = v_2
            return False
        self.bra = self.cursor
        S_ch = self.slice_to()
        self.limit_backward = v_2
        if not self.eq_s_b(S_ch):
            return False
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
        v_4 = self.limit - self.cursor
        self.__r_undouble()
        self.cursor = self.limit - v_4
        self.ket = self.cursor
        if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "'":
            return False
        self.cursor -= 1
        self.bra = self.cursor
        self.slice_del()
        self.cursor = self.limit_backward
        return True

    a_0 = [
        Among("hed", -1, 1),
        Among("ethed", 0, 1),
        Among("ered", -1, 1),
        Among("e", -1, 1),
        Among("erede", 3, 1),
        Among("ende", 3, 1),
        Among("erende", 5, 1),
        Among("ene", 3, 1),
        Among("erne", 3, 1),
        Among("ere", 3, 1),
        Among("en", -1, 1),
        Among("heden", 10, 1),
        Among("eren", 10, 1),
        Among("er", -1, 1),
        Among("heder", 13, 1),
        Among("erer", 13, 1),
        Among("s", -1, 2),
        Among("heds", 16, 1),
        Among("es", 16, 1),
        Among("endes", 18, 1),
        Among("erendes", 19, 1),
        Among("enes", 18, 1),
        Among("ernes", 18, 1),
        Among("eres", 18, 1),
        Among("ens", 16, 1),
        Among("hedens", 24, 1),
        Among("erens", 24, 1),
        Among("ers", 16, 1),
        Among("ets", 16, 1),
        Among("erets", 28, 1),
        Among("et", -1, 1),
        Among("eret", 30, 1)
    ]

    a_1 = [
        Among("gd", -1, -1),
        Among("dt", -1, -1),
        Among("gt", -1, -1),
        Among("kt", -1, -1)
    ]

    a_2 = [
        Among("ig", -1, 1),
        Among("lig", 0, 1),
        Among("elig", 1, 1),
        Among("els", -1, 1),
        Among("løst", -1, 2)
    ]


class lab0(BaseException): pass


class lab1(BaseException): pass


class lab2(BaseException): pass
