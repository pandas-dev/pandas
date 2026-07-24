# Generated from irish.sbl by Snowball 3.1.1 - https://snowballstem.org/

from .basestemmer import BaseStemmer
from .among import Among


class IrishStemmer(BaseStemmer):
    '''
    This class implements the stemming algorithm defined by a snowball script.
    Generated from irish.sbl by Snowball 3.1.1 - https://snowballstem.org/
    '''

    g_v = {"a", "e", "i", "o", "u", "á", "é", "í", "ó", "ú"}

    I_p2 = 0
    I_p1 = 0
    I_pV = 0

    def __r_mark_regions(self):
        self.I_pV = self.limit
        self.I_p1 = self.limit
        self.I_p2 = self.limit
        v_1 = self.cursor
        try:
            if not self.go_out_grouping(IrishStemmer.g_v):
                raise lab0()
            self.cursor += 1
            self.I_pV = self.cursor
            if not self.go_in_grouping(IrishStemmer.g_v):
                raise lab0()
            self.cursor += 1
            self.I_p1 = self.cursor
            if not self.go_out_grouping(IrishStemmer.g_v):
                raise lab0()
            self.cursor += 1
            if not self.go_in_grouping(IrishStemmer.g_v):
                raise lab0()
            self.cursor += 1
            self.I_p2 = self.cursor
        except lab0: pass
        self.cursor = v_1
        return True

    def __r_initial_morph(self):
        self.bra = self.cursor
        among_var = self.find_among(IrishStemmer.a_0)
        if among_var == 0:
            return False
        self.ket = self.cursor
        self.slice_from(IrishStemmer.as_0[among_var - 1])
        return True

    def __r_R1(self):
        return self.I_p1 <= self.cursor

    def __r_R2(self):
        return self.I_p2 <= self.cursor

    def __r_noun_sfx(self):
        self.ket = self.cursor
        among_var = self.find_among_b(IrishStemmer.a_1)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var == 1:
            if not self.__r_R1():
                return False
            self.slice_del()
        else:
            if not self.__r_R2():
                return False
            self.slice_del()
        return True

    def __r_deriv(self):
        self.ket = self.cursor
        among_var = self.find_among_b(IrishStemmer.a_2)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var == 1:
            if not self.__r_R2():
                return False
            self.slice_del()
        elif among_var == 2:
            self.slice_from("arc")
        elif among_var == 3:
            self.slice_from("gin")
        elif among_var == 4:
            self.slice_from("graf")
        elif among_var == 5:
            self.slice_from("paite")
        else:
            self.slice_from("óid")
        return True

    def __r_verb_sfx(self):
        self.ket = self.cursor
        among_var = self.find_among_b(IrishStemmer.a_3)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var == 1:
            if self.I_pV > self.cursor:
                return False
            self.slice_del()
        else:
            if not self.__r_R1():
                return False
            self.slice_del()
        return True

    def _stem(self):
        v_1 = self.cursor
        self.__r_initial_morph()
        self.cursor = v_1
        self.__r_mark_regions()
        self.limit_backward = self.cursor
        self.cursor = self.limit
        v_2 = self.limit - self.cursor
        self.__r_noun_sfx()
        self.cursor = self.limit - v_2
        v_3 = self.limit - self.cursor
        self.__r_deriv()
        self.cursor = self.limit - v_3
        v_4 = self.limit - self.cursor
        self.__r_verb_sfx()
        self.cursor = self.limit - v_4
        self.cursor = self.limit_backward
        return True

    a_0 = [
        Among("b'", -1, 1),
        Among("bh", -1, 4),
        Among("bhf", 1, 2),
        Among("bp", -1, 8),
        Among("ch", -1, 5),
        Among("d'", -1, 1),
        Among("d'fh", 5, 2),
        Among("dh", -1, 6),
        Among("dt", -1, 9),
        Among("fh", -1, 2),
        Among("gc", -1, 5),
        Among("gh", -1, 7),
        Among("h-", -1, 1),
        Among("m'", -1, 1),
        Among("mb", -1, 4),
        Among("mh", -1, 10),
        Among("n-", -1, 1),
        Among("nd", -1, 6),
        Among("ng", -1, 7),
        Among("ph", -1, 8),
        Among("sh", -1, 3),
        Among("t-", -1, 1),
        Among("th", -1, 9),
        Among("ts", -1, 3)
    ]
    as_0 = ("", "f", "s", "b", "c", "d", "g", "p", "t", "m")

    a_1 = [
        Among("íochta", -1, 1),
        Among("aíochta", 0, 1),
        Among("ire", -1, 2),
        Among("aire", 2, 2),
        Among("abh", -1, 1),
        Among("eabh", 4, 1),
        Among("ibh", -1, 1),
        Among("aibh", 6, 1),
        Among("amh", -1, 1),
        Among("eamh", 8, 1),
        Among("imh", -1, 1),
        Among("aimh", 10, 1),
        Among("íocht", -1, 1),
        Among("aíocht", 12, 1),
        Among("irí", -1, 2),
        Among("airí", 14, 2)
    ]

    a_2 = [
        Among("óideacha", -1, 6),
        Among("patacha", -1, 5),
        Among("achta", -1, 1),
        Among("arcachta", 2, 2),
        Among("eachta", 2, 1),
        Among("grafaíochta", -1, 4),
        Among("paite", -1, 5),
        Among("ach", -1, 1),
        Among("each", 7, 1),
        Among("óideach", 8, 6),
        Among("gineach", 8, 3),
        Among("patach", 7, 5),
        Among("grafaíoch", -1, 4),
        Among("pataigh", -1, 5),
        Among("óidigh", -1, 6),
        Among("achtúil", -1, 1),
        Among("eachtúil", 15, 1),
        Among("gineas", -1, 3),
        Among("ginis", -1, 3),
        Among("acht", -1, 1),
        Among("arcacht", 19, 2),
        Among("eacht", 19, 1),
        Among("grafaíocht", -1, 4),
        Among("arcachtaí", -1, 2),
        Among("grafaíochtaí", -1, 4)
    ]

    a_3 = [
        Among("imid", -1, 1),
        Among("aimid", 0, 1),
        Among("ímid", -1, 1),
        Among("aímid", 2, 1),
        Among("adh", -1, 2),
        Among("eadh", 4, 2),
        Among("faidh", -1, 1),
        Among("fidh", -1, 1),
        Among("áil", -1, 2),
        Among("ain", -1, 2),
        Among("tear", -1, 2),
        Among("tar", -1, 2)
    ]


class lab0(BaseException): pass
