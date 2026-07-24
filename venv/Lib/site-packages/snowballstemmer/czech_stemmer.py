# Generated from czech.sbl by Snowball 3.1.1 - https://snowballstem.org/

from .basestemmer import BaseStemmer
from .among import Among


class CzechStemmer(BaseStemmer):
    '''
    This class implements the stemming algorithm defined by a snowball script.
    Generated from czech.sbl by Snowball 3.1.1 - https://snowballstem.org/
    '''

    g_v = {"a", "e", "i", "o", "u", "y", "á", "é", "í", "ó", "ú", "ý", "ě", "ů"}

    g_v_or_syllabic_c = {"a", "e", "i", "l", "o", "r", "u", "y", "á", "é", "í", "ó", "ú", "ý", "ě", "ů"}

    g_ev_ending = {"h", "k", "n", "r", "t", "z"}

    g_env_ending = {"b", "c", "d", "h", "k", "p", "r", "s", "t", "v", "z", "č", "š", "ž"}

    I_p1 = 0

    def __r_mark_regions(self):
        v_1 = self.cursor
        if self.cursor + 3 > self.limit:
            return False
        self.cursor += 3
        I_x = self.cursor
        self.cursor = v_1
        self.I_p1 = self.limit
        v_2 = self.cursor
        try:
            while True:
                try:
                    if not self.in_grouping(CzechStemmer.g_v):
                        raise lab1()
                    break
                except lab1: pass
                if self.cursor >= self.limit:
                    raise lab0()
                self.cursor += 1
                if not self.go_out_grouping(CzechStemmer.g_v_or_syllabic_c):
                    raise lab0()
                self.cursor += 1
                break
            if not self.go_in_grouping(CzechStemmer.g_v):
                raise lab0()
            self.cursor += 1
            self.I_p1 = self.cursor
            try:
                if self.I_p1 >= I_x:
                    raise lab1()
                self.I_p1 = I_x
            except lab1: pass
        except lab0: pass
        self.cursor = v_2
        return True

    def __r_palatalise_e(self):
        self.ket = self.cursor
        among_var = self.find_among_b(CzechStemmer.a_0)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var > 0:
            self.slice_from(CzechStemmer.as_0[among_var - 1])
        return True

    def __r_palatalise_i(self):
        self.ket = self.cursor
        among_var = self.find_among_b(CzechStemmer.a_1)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var > 0:
            self.slice_from(CzechStemmer.as_1[among_var - 1])
        return True

    def __r_possessive_suffix(self):
        self.ket = self.cursor
        among_var = self.find_among_b(CzechStemmer.a_2)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if self.I_p1 > self.cursor:
            return False
        if among_var == 1:
            self.slice_del()
        else:
            self.slice_del()
            v_1 = self.limit - self.cursor
            try:
                if not self.__r_palatalise_i():
                    self.cursor = self.limit - v_1
                    raise lab0()
            except lab0: pass
        return True

    def __r_case_suffix(self):
        if self.cursor < self.I_p1:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_p1
        self.ket = self.cursor
        among_var = self.find_among_b(CzechStemmer.a_6)
        if among_var == 0:
            self.limit_backward = v_2
            return False
        self.bra = self.cursor
        self.limit_backward = v_2
        if among_var == 1:
            self.slice_del()
        elif among_var == 2:
            self.slice_del()
            v_3 = self.limit - self.cursor
            try:
                if not self.__r_palatalise_e():
                    self.cursor = self.limit - v_3
                    raise lab0()
            except lab0: pass
        elif among_var == 3:
            among_var = self.find_among_b(CzechStemmer.a_3)
            self.slice_from(CzechStemmer.as_3[among_var - 1])
        elif among_var == 4:
            v_4 = self.limit - self.cursor
            if not self.out_grouping_b(CzechStemmer.g_v):
                return False
            self.cursor = self.limit - v_4
            try:
                if not self.eq_s_b("tř"):
                    raise lab0()
                return False
            except lab0: pass
            self.slice_from("b")
        elif among_var == 5:
            v_5 = self.limit - self.cursor
            if not self.out_grouping_b(CzechStemmer.g_v):
                return False
            self.cursor = self.limit - v_5
            self.slice_del()
            self.insert(self.cursor, self.cursor, "c")
            v_6 = self.limit - self.cursor
            try:
                if not self.__r_palatalise_e():
                    self.cursor = self.limit - v_6
                    raise lab0()
            except lab0: pass
        elif among_var == 6:
            v_7 = self.limit - self.cursor
            if not self.out_grouping_b(CzechStemmer.g_v):
                return False
            self.cursor = self.limit - v_7
            v_8 = self.limit - self.cursor
            try:
                if self.find_among_b(CzechStemmer.a_4) == 0:
                    raise lab0()
                return False
            except lab0: pass
            self.cursor = self.limit - v_8
            self.slice_from("k")
        elif among_var == 7:
            if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "n":
                return False
            self.cursor -= 1
            self.bra = self.cursor
            self.slice_from("ňk")
        elif among_var == 8:
            v_9 = self.limit - self.cursor
            if not self.in_grouping_b(CzechStemmer.g_env_ending):
                return False
            self.cursor = self.limit - v_9
            self.slice_from("n")
        elif among_var == 9:
            if self.find_among_b(CzechStemmer.a_5) == 0:
                return False
            self.slice_from("t")
        elif among_var == 10:
            if not self.in_grouping_b(CzechStemmer.g_ev_ending):
                return False
            self.slice_from("v")
        elif among_var == 11:
            self.slice_from("t")
        else:
            self.slice_del()
            v_10 = self.limit - self.cursor
            try:
                if not self.__r_palatalise_i():
                    self.cursor = self.limit - v_10
                    raise lab0()
            except lab0: pass
        return True

    def _stem(self):
        if not self.__r_mark_regions():
            return False
        self.limit_backward = self.cursor
        self.cursor = self.limit
        v_1 = self.limit - self.cursor
        self.__r_case_suffix()
        self.cursor = self.limit - v_1
        v_2 = self.limit - self.cursor
        self.__r_possessive_suffix()
        self.cursor = self.limit - v_2
        self.cursor = self.limit_backward
        return True

    a_0 = [
        Among("c", -1, 1),
        Among("nc", 0, -1),
        Among("ínc", 1, 2),
        Among("avc", 0, -1),
        Among("ovc", 0, -1)
    ]
    as_0 = ("k", "ínk")

    a_1 = [
        Among("c", -1, 1),
        Among("nc", 0, -1),
        Among("ínc", 1, 2),
        Among("avc", 0, -1),
        Among("ovc", 0, -1),
        Among("čt", -1, 3),
        Among("št", -1, 4),
        Among("dešt", 6, -1),
        Among("lešt", 6, -1),
        Among("išt", 6, -1),
        Among("poušt", 6, -1),
        Among("ášt", 6, -1),
        Among("íšt", 6, -1)
    ]
    as_1 = ("k", "ínk", "ck", "sk")

    a_2 = [
        Among("in", -1, 2),
        Among("ov", -1, 1),
        Among("ův", -1, 1)
    ]

    a_3 = [
        Among("", -1, 2),
        Among("l", 0, 1),
        Among("tl", 1, 2),
        Among("s", 0, 1),
        Among("es", 3, 2),
        Among("č", 0, 1),
        Among("eč", 5, 2),
        Among("ř", 0, 1),
        Among("ž", 0, 1)
    ]
    as_3 = ("", "et")

    a_4 = [
        Among("obl", -1, -1),
        Among("sn", -1, -1),
        Among("dot", -1, -1)
    ]

    a_5 = [
        Among("uc", -1, -1),
        Among("h", -1, -1),
        Among("ok", -1, -1),
        Among("kar", -1, -1),
        Among("č", -1, -1)
    ]

    a_6 = [
        Among("a", -1, 1),
        Among("ama", 0, 1),
        Among("ata", 0, 1),
        Among("eb", -1, 4),
        Among("ec", -1, 5),
        Among("e", -1, 2),
        Among("ete", 5, 3),
        Among("ěte", 5, 1),
        Among("ech", -1, 2),
        Among("atech", 8, 1),
        Among("ách", -1, 1),
        Among("ích", -1, 12),
        Among("ých", -1, 1),
        Among("i", -1, 12),
        Among("mi", 13, 1),
        Among("ami", 14, 1),
        Among("emi", 14, 2),
        Among("ími", 14, 12),
        Among("ými", 14, 1),
        Among("ěmi", 14, 1),
        Among("ťmi", 14, 11),
        Among("eti", 13, 3),
        Among("ěti", 13, 1),
        Among("ovi", 13, 1),
        Among("ek", -1, 6),
        Among("ěk", -1, 7),
        Among("em", -1, 2),
        Among("etem", 26, 3),
        Among("ětem", 26, 1),
        Among("ám", -1, 1),
        Among("ém", -1, 1),
        Among("ím", -1, 12),
        Among("ým", -1, 1),
        Among("ěm", -1, 1),
        Among("ům", -1, 1),
        Among("atům", 34, 1),
        Among("o", -1, 1),
        Among("ého", 36, 1),
        Among("ího", 36, 12),
        Among("us", -1, 1),
        Among("at", -1, 1),
        Among("et", -1, 9),
        Among("u", -1, 1),
        Among("ému", 42, 1),
        Among("ímu", 42, 12),
        Among("ou", 42, 1),
        Among("ev", -1, 10),
        Among("y", -1, 1),
        Among("aty", 47, 1),
        Among("á", -1, 1),
        Among("é", -1, 1),
        Among("ové", 50, 1),
        Among("í", -1, 12),
        Among("ý", -1, 1),
        Among("ě", -1, 1),
        Among("eň", -1, 8),
        Among("ť", -1, 11),
        Among("ů", -1, 1)
    ]


class lab0(BaseException): pass


class lab1(BaseException): pass
