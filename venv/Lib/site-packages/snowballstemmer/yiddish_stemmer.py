# Generated from yiddish.sbl by Snowball 3.1.1 - https://snowballstem.org/

from .basestemmer import BaseStemmer
from .among import Among


class YiddishStemmer(BaseStemmer):
    '''
    This class implements the stemming algorithm defined by a snowball script.
    Generated from yiddish.sbl by Snowball 3.1.1 - https://snowballstem.org/
    '''

    g_niked = {"\u05B0", "\u05B1", "\u05B2", "\u05B3", "\u05B4", "\u05B5", "\u05B6", "\u05B7", "\u05B8", "\u05B9", "\u05BB", "\u05BC", "\u05BF", "\u05C1", "\u05C2"}

    g_vowel = {"\u05D0", "\u05D5", "\u05D9", "\u05E2", "\u05F1", "\u05F2"}

    g_consonant = {"\u05D1", "\u05D2", "\u05D3", "\u05D4", "\u05D6", "\u05D7", "\u05D8", "\u05DA", "\u05DB", "\u05DC", "\u05DD", "\u05DE", "\u05DF", "\u05E0", "\u05E1", "\u05E3", "\u05E4", "\u05E5", "\u05E6", "\u05E7", "\u05E8", "\u05E9", "\u05EA", "\u05F0"}

    I_p1 = 0

    def __r_prelude(self):
        v_1 = self.cursor
        try:
            while True:
                v_2 = self.cursor
                try:
                    while True:
                        v_3 = self.cursor
                        try:
                            self.bra = self.cursor
                            among_var = self.find_among(YiddishStemmer.a_0)
                            if among_var == 0:
                                raise lab2()
                            self.ket = self.cursor
                            if among_var == 1:
                                try:
                                    if self.cursor == self.limit or self.current[self.cursor] != "\u05BC":
                                        raise lab3()
                                    self.cursor += 1
                                    raise lab2()
                                except lab3: pass
                                self.slice_from("\u05F0")
                            elif among_var == 2:
                                try:
                                    if self.cursor == self.limit or self.current[self.cursor] != "\u05B4":
                                        raise lab3()
                                    self.cursor += 1
                                    raise lab2()
                                except lab3: pass
                                self.slice_from("\u05F1")
                            elif among_var == 3:
                                try:
                                    if self.cursor == self.limit or self.current[self.cursor] != "\u05B4":
                                        raise lab3()
                                    self.cursor += 1
                                    raise lab2()
                                except lab3: pass
                                self.slice_from("\u05F2")
                            elif among_var == 4:
                                self.slice_from("\u05DB")
                            elif among_var == 5:
                                self.slice_from("\u05DE")
                            elif among_var == 6:
                                self.slice_from("\u05E0")
                            elif among_var == 7:
                                self.slice_from("\u05E4")
                            else:
                                self.slice_from("\u05E6")
                            self.cursor = v_3
                            break
                        except lab2: pass
                        self.cursor = v_3
                        if self.cursor >= self.limit:
                            raise lab1()
                        self.cursor += 1
                    continue
                except lab1: pass
                self.cursor = v_2
                break
        except lab0: pass
        self.cursor = v_1
        v_4 = self.cursor
        try:
            while True:
                v_5 = self.cursor
                try:
                    while True:
                        v_6 = self.cursor
                        try:
                            self.bra = self.cursor
                            if not self.in_grouping(YiddishStemmer.g_niked):
                                raise lab2()
                            self.ket = self.cursor
                            self.slice_del()
                            self.cursor = v_6
                            break
                        except lab2: pass
                        self.cursor = v_6
                        if self.cursor >= self.limit:
                            raise lab1()
                        self.cursor += 1
                    continue
                except lab1: pass
                self.cursor = v_5
                break
        except lab0: pass
        self.cursor = v_4
        return True

    def __r_mark_regions(self):
        self.I_p1 = self.limit
        v_1 = self.cursor
        try:
            self.bra = self.cursor
            if not self.eq_s("\u05D2\u05E2"):
                self.cursor = v_1
                raise lab0()
            self.ket = self.cursor
            v_2 = self.cursor
            try:
                while True:
                    try:
                        if not self.eq_s("\u05DC\u05D8"):
                            raise lab2()
                        break
                    except lab2: pass
                    try:
                        if not self.eq_s("\u05D1\u05E0"):
                            raise lab2()
                        break
                    except lab2: pass
                    if self.cursor < self.limit:
                        raise lab1()
                    break
                self.cursor = v_1
                raise lab0()
            except lab1: pass
            self.cursor = v_2
            self.slice_from("GE")
        except lab0: pass
        v_3 = self.cursor
        try:
            if self.find_among(YiddishStemmer.a_1) == 0:
                self.cursor = v_3
                raise lab0()
            while True:
                v_4 = self.cursor
                try:
                    v_5 = self.cursor
                    while True:
                        try:
                            if not self.eq_s("\u05E6\u05D5\u05D2\u05E0"):
                                raise lab2()
                            break
                        except lab2: pass
                        try:
                            if not self.eq_s("\u05E6\u05D5\u05E7\u05D8"):
                                raise lab2()
                            break
                        except lab2: pass
                        if not self.eq_s("\u05E6\u05D5\u05E7\u05E0"):
                            raise lab1()
                        break
                    if self.cursor < self.limit:
                        raise lab1()
                    self.cursor = v_5
                    break
                except lab1: pass
                self.cursor = v_4
                try:
                    v_6 = self.cursor
                    if not self.eq_s("\u05D2\u05E2\u05D1\u05E0"):
                        raise lab1()
                    self.cursor = v_6
                    break
                except lab1: pass
                self.cursor = v_4
                try:
                    self.bra = self.cursor
                    if not self.eq_s("\u05D2\u05E2"):
                        raise lab1()
                    self.ket = self.cursor
                    self.slice_from("GE")
                    break
                except lab1: pass
                self.cursor = v_4
                self.bra = self.cursor
                if not self.eq_s("\u05E6\u05D5"):
                    self.cursor = v_3
                    raise lab0()
                self.ket = self.cursor
                self.slice_from("TSU")
                break
        except lab0: pass
        v_7 = self.cursor
        if self.cursor + 3 > self.limit:
            return False
        self.cursor += 3
        I_x = self.cursor
        self.cursor = v_7
        v_8 = self.cursor
        try:
            if self.find_among(YiddishStemmer.a_2) == 0:
                self.cursor = v_8
                raise lab0()
        except lab0: pass
        v_9 = self.cursor
        try:
            if not self.in_grouping(YiddishStemmer.g_consonant):
                raise lab0()
            if not self.in_grouping(YiddishStemmer.g_consonant):
                raise lab0()
            if not self.in_grouping(YiddishStemmer.g_consonant):
                raise lab0()
            self.I_p1 = self.cursor
            return False
        except lab0: pass
        self.cursor = v_9
        if not self.go_out_grouping(YiddishStemmer.g_vowel):
            return False
        self.cursor += 1
        if not self.go_in_grouping(YiddishStemmer.g_vowel):
            return False
        self.I_p1 = self.cursor
        try:
            if self.I_p1 >= I_x:
                raise lab0()
            self.I_p1 = I_x
        except lab0: pass
        return True

    def __r_R1(self):
        return self.I_p1 <= self.cursor

    def __r_standard_suffix(self):
        v_1 = self.limit - self.cursor
        try:
            self.ket = self.cursor
            among_var = self.find_among_b(YiddishStemmer.a_4)
            if among_var == 0:
                raise lab0()
            self.bra = self.cursor
            if among_var == 1:
                if not self.__r_R1():
                    raise lab0()
                self.slice_del()
            elif among_var == 2:
                if not self.__r_R1():
                    raise lab0()
                self.slice_from("\u05D9\u05E2")
            elif among_var == 3:
                if not self.__r_R1():
                    raise lab0()
                self.slice_del()
                self.ket = self.cursor
                among_var = self.find_among_b(YiddishStemmer.a_3)
                if among_var == 0:
                    raise lab0()
                self.bra = self.cursor
                self.slice_from(YiddishStemmer.as_3[among_var - 1])
            elif among_var == 4:
                while True:
                    v_2 = self.limit - self.cursor
                    try:
                        if not self.__r_R1():
                            raise lab1()
                        self.slice_del()
                        break
                    except lab1: pass
                    self.cursor = self.limit - v_2
                    self.slice_from("\u05D8")
                    break
                self.ket = self.cursor
                if not self.eq_s_b("\u05D1\u05E8\u05D0\u05DB"):
                    raise lab0()
                v_3 = self.limit - self.cursor
                try:
                    if not self.eq_s_b("\u05D2\u05E2"):
                        self.cursor = self.limit - v_3
                        raise lab1()
                except lab1: pass
                self.bra = self.cursor
                self.slice_from("\u05D1\u05E8\u05E2\u05E0\u05D2")
            elif among_var == 5:
                self.slice_from("\u05D2\u05F2")
            elif among_var == 6:
                self.slice_from("\u05E0\u05E2\u05DE")
            elif among_var == 7:
                self.slice_from("\u05E9\u05E8\u05F2\u05D1")
            elif among_var == 8:
                self.slice_from("\u05DE\u05F2\u05D3")
            elif among_var == 9:
                self.slice_from("\u05D1\u05F2\u05D8")
            elif among_var == 10:
                self.slice_from("\u05D1\u05F2\u05E1")
            elif among_var == 11:
                self.slice_from("\u05F0\u05F2\u05D6")
            elif among_var == 12:
                self.slice_from("\u05D8\u05E8\u05F2\u05D1")
            elif among_var == 13:
                self.slice_from("\u05DC\u05F2\u05D8")
            elif among_var == 14:
                self.slice_from("\u05E7\u05DC\u05F2\u05D1")
            elif among_var == 15:
                self.slice_from("\u05E8\u05F2\u05D1")
            elif among_var == 16:
                self.slice_from("\u05E8\u05F2\u05E1")
            elif among_var == 17:
                self.slice_from("\u05E9\u05F0\u05F2\u05D2")
            elif among_var == 18:
                self.slice_from("\u05E9\u05DE\u05F2\u05E1")
            elif among_var == 19:
                self.slice_from("\u05E9\u05E0\u05F2\u05D3")
            elif among_var == 20:
                self.slice_from("\u05D1\u05D9\u05E0\u05D3")
            elif among_var == 21:
                self.slice_from("\u05F0\u05D9\u05D8\u05E9")
            elif among_var == 22:
                self.slice_from("\u05D6\u05D9\u05E0\u05D2")
            elif among_var == 23:
                self.slice_from("\u05D8\u05E8\u05D9\u05E0\u05E7")
            elif among_var == 24:
                self.slice_from("\u05E6\u05F0\u05D9\u05E0\u05D2")
            elif among_var == 25:
                self.slice_from("\u05E9\u05DC\u05D9\u05E0\u05D2")
            elif among_var == 26:
                self.slice_from("\u05D1\u05F2\u05D2")
            elif among_var == 27:
                self.slice_from("\u05D4\u05F2\u05D1")
            elif among_var == 28:
                self.slice_from("\u05E4\u05D0\u05E8\u05DC\u05D9\u05E8")
            elif among_var == 29:
                self.slice_from("\u05E9\u05D8\u05F2")
            elif among_var == 30:
                self.slice_from("\u05E9\u05F0\u05E2\u05E8")
            elif among_var == 31:
                self.slice_from("\u05D1\u05E8\u05E2\u05E0\u05D2")
            elif among_var == 32:
                if not self.__r_R1():
                    raise lab0()
                self.slice_from("\u05D4")
            elif among_var == 33:
                while True:
                    v_4 = self.limit - self.cursor
                    try:
                        while True:
                            try:
                                if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "\u05D2":
                                    raise lab2()
                                self.cursor -= 1
                                break
                            except lab2: pass
                            if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "\u05E9":
                                raise lab1()
                            self.cursor -= 1
                            break
                        v_5 = self.limit - self.cursor
                        try:
                            if self.I_p1 > (self.cursor + 3):
                                self.cursor = self.limit - v_5
                                raise lab2()
                            self.slice_from("\u05D9\u05E1")
                        except lab2: pass
                        break
                    except lab1: pass
                    self.cursor = self.limit - v_4
                    if not self.__r_R1():
                        raise lab0()
                    self.slice_del()
                    break
        except lab0: pass
        self.cursor = self.limit - v_1
        v_6 = self.limit - self.cursor
        try:
            self.ket = self.cursor
            among_var = self.find_among_b(YiddishStemmer.a_5)
            if among_var == 0:
                raise lab0()
            self.bra = self.cursor
            if among_var == 1:
                if not self.__r_R1():
                    raise lab0()
                self.slice_del()
            else:
                if not self.__r_R1():
                    raise lab0()
                if not self.in_grouping_b(YiddishStemmer.g_consonant):
                    raise lab0()
                self.slice_del()
        except lab0: pass
        self.cursor = self.limit - v_6
        v_7 = self.limit - self.cursor
        try:
            self.ket = self.cursor
            among_var = self.find_among_b(YiddishStemmer.a_6)
            if among_var == 0:
                raise lab0()
            self.bra = self.cursor
            if among_var == 1:
                if not self.__r_R1():
                    raise lab0()
                self.slice_del()
        except lab0: pass
        self.cursor = self.limit - v_7
        v_8 = self.limit - self.cursor
        try:
            while True:
                v_9 = self.limit - self.cursor
                try:
                    while True:
                        v_10 = self.limit - self.cursor
                        try:
                            self.ket = self.cursor
                            while True:
                                try:
                                    if not self.eq_s_b("GE"):
                                        raise lab3()
                                    break
                                except lab3: pass
                                if not self.eq_s_b("TSU"):
                                    raise lab2()
                                break
                            self.bra = self.cursor
                            self.slice_del()
                            self.cursor = self.limit - v_10
                            break
                        except lab2: pass
                        self.cursor = self.limit - v_10
                        if self.cursor <= self.limit_backward:
                            raise lab1()
                        self.cursor -= 1
                    continue
                except lab1: pass
                self.cursor = self.limit - v_9
                break
        except lab0: pass
        self.cursor = self.limit - v_8
        return True

    def _stem(self):
        self.__r_prelude()
        v_1 = self.cursor
        self.__r_mark_regions()
        self.cursor = v_1
        self.limit_backward = self.cursor
        self.cursor = self.limit
        self.__r_standard_suffix()
        self.cursor = self.limit_backward
        return True

    a_0 = [
        Among("\u05D5\u05D5", -1, 1),
        Among("\u05D5\u05D9", -1, 2),
        Among("\u05D9\u05D9", -1, 3),
        Among("\u05DA", -1, 4),
        Among("\u05DD", -1, 5),
        Among("\u05DF", -1, 6),
        Among("\u05E3", -1, 7),
        Among("\u05E5", -1, 8)
    ]

    a_1 = [
        Among("\u05D0\u05D3\u05D5\u05E8\u05DB", -1, 1),
        Among("\u05D0\u05D4\u05D9\u05E0", -1, 1),
        Among("\u05D0\u05D4\u05E2\u05E8", -1, 1),
        Among("\u05D0\u05D4\u05F2\u05DE", -1, 1),
        Among("\u05D0\u05D5\u05DE", -1, 1),
        Among("\u05D0\u05D5\u05E0\u05D8\u05E2\u05E8", -1, 1),
        Among("\u05D0\u05D9\u05D1\u05E2\u05E8", -1, 1),
        Among("\u05D0\u05E0", -1, 1),
        Among("\u05D0\u05E0\u05D8", 7, 1),
        Among("\u05D0\u05E0\u05D8\u05E7\u05E2\u05D2\u05E0", 8, 1),
        Among("\u05D0\u05E0\u05D9\u05D3\u05E2\u05E8", 7, 1),
        Among("\u05D0\u05E4", -1, 1),
        Among("\u05D0\u05E4\u05D9\u05E8", 11, 1),
        Among("\u05D0\u05E7\u05E2\u05D2\u05E0", -1, 1),
        Among("\u05D0\u05E8\u05D0\u05E4", -1, 1),
        Among("\u05D0\u05E8\u05D5\u05DE", -1, 1),
        Among("\u05D0\u05E8\u05D5\u05E0\u05D8\u05E2\u05E8", -1, 1),
        Among("\u05D0\u05E8\u05D9\u05D1\u05E2\u05E8", -1, 1),
        Among("\u05D0\u05E8\u05F1\u05E1", -1, 1),
        Among("\u05D0\u05E8\u05F1\u05E4", -1, 1),
        Among("\u05D0\u05E8\u05F2\u05E0", -1, 1),
        Among("\u05D0\u05F0\u05E2\u05E7", -1, 1),
        Among("\u05D0\u05F1\u05E1", -1, 1),
        Among("\u05D0\u05F1\u05E4", -1, 1),
        Among("\u05D0\u05F2\u05E0", -1, 1),
        Among("\u05D1\u05D0", -1, 1),
        Among("\u05D1\u05F2", -1, 1),
        Among("\u05D3\u05D5\u05E8\u05DB", -1, 1),
        Among("\u05D3\u05E2\u05E8", -1, 1),
        Among("\u05DE\u05D9\u05D8", -1, 1),
        Among("\u05E0\u05D0\u05DB", -1, 1),
        Among("\u05E4\u05D0\u05E8", -1, 1),
        Among("\u05E4\u05D0\u05E8\u05D1\u05F2", 31, 1),
        Among("\u05E4\u05D0\u05E8\u05F1\u05E1", 31, 1),
        Among("\u05E4\u05D5\u05E0\u05D0\u05E0\u05D3\u05E2\u05E8", -1, 1),
        Among("\u05E6\u05D5", -1, 1),
        Among("\u05E6\u05D5\u05D6\u05D0\u05DE\u05E2\u05E0", 35, 1),
        Among("\u05E6\u05D5\u05E0\u05F1\u05E4", 35, 1),
        Among("\u05E6\u05D5\u05E8\u05D9\u05E7", 35, 1),
        Among("\u05E6\u05E2", -1, 1)
    ]

    a_2 = [
        Among("\u05D3\u05D6\u05E9", -1, -1),
        Among("\u05E9\u05D8\u05E8", -1, -1),
        Among("\u05E9\u05D8\u05E9", -1, -1),
        Among("\u05E9\u05E4\u05E8", -1, -1)
    ]

    a_3 = [
        Among("\u05E7\u05DC\u05D9\u05D1", -1, 9),
        Among("\u05E8\u05D9\u05D1", -1, 10),
        Among("\u05D8\u05E8\u05D9\u05D1", 1, 7),
        Among("\u05E9\u05E8\u05D9\u05D1", 1, 15),
        Among("\u05D4\u05F1\u05D1", -1, 23),
        Among("\u05E9\u05F0\u05D9\u05D2", -1, 12),
        Among("\u05D2\u05D0\u05E0\u05D2", -1, 1),
        Among("\u05D6\u05D5\u05E0\u05D2", -1, 18),
        Among("\u05E9\u05DC\u05D5\u05E0\u05D2", -1, 21),
        Among("\u05E6\u05F0\u05D5\u05E0\u05D2", -1, 20),
        Among("\u05D1\u05F1\u05D2", -1, 22),
        Among("\u05D1\u05D5\u05E0\u05D3", -1, 16),
        Among("\u05F0\u05D9\u05D6", -1, 6),
        Among("\u05D1\u05D9\u05D8", -1, 4),
        Among("\u05DC\u05D9\u05D8", -1, 8),
        Among("\u05DE\u05D9\u05D8", -1, 3),
        Among("\u05E9\u05E0\u05D9\u05D8", -1, 14),
        Among("\u05E0\u05D5\u05DE", -1, 2),
        Among("\u05E9\u05D8\u05D0\u05E0", -1, 25),
        Among("\u05D1\u05D9\u05E1", -1, 5),
        Among("\u05E9\u05DE\u05D9\u05E1", -1, 13),
        Among("\u05E8\u05D9\u05E1", -1, 11),
        Among("\u05D8\u05E8\u05D5\u05E0\u05E7", -1, 19),
        Among("\u05E4\u05D0\u05E8\u05DC\u05F1\u05E8", -1, 24),
        Among("\u05E9\u05F0\u05F1\u05E8", -1, 26),
        Among("\u05F0\u05D5\u05D8\u05E9", -1, 17)
    ]
    as_3 = ("\u05D2\u05F2", "\u05E0\u05E2\u05DE", "\u05DE\u05F2\u05D3", "\u05D1\u05F2\u05D8", "\u05D1\u05F2\u05E1", "\u05F0\u05F2\u05D6", "\u05D8\u05E8\u05F2\u05D1", "\u05DC\u05F2\u05D8", "\u05E7\u05DC\u05F2\u05D1", "\u05E8\u05F2\u05D1", "\u05E8\u05F2\u05E1", "\u05E9\u05F0\u05F2\u05D2", "\u05E9\u05DE\u05F2\u05E1", "\u05E9\u05E0\u05F2\u05D3", "\u05E9\u05E8\u05F2\u05D1", "\u05D1\u05D9\u05E0\u05D3", "\u05F0\u05D9\u05D8\u05E9", "\u05D6\u05D9\u05E0\u05D2", "\u05D8\u05E8\u05D9\u05E0\u05E7", "\u05E6\u05F0\u05D9\u05E0\u05D2", "\u05E9\u05DC\u05D9\u05E0\u05D2", "\u05D1\u05F2\u05D2", "\u05D4\u05F2\u05D1", "\u05E4\u05D0\u05E8\u05DC\u05D9\u05E8", "\u05E9\u05D8\u05F2", "\u05E9\u05F0\u05E2\u05E8")

    a_4 = [
        Among("\u05D5\u05E0\u05D2", -1, 1),
        Among("\u05E1\u05D8\u05D5", -1, 1),
        Among("\u05D8", -1, 1),
        Among("\u05D1\u05E8\u05D0\u05DB\u05D8", 2, 31),
        Among("\u05E1\u05D8", 2, 1),
        Among("\u05D9\u05E1\u05D8", 4, 33),
        Among("\u05E2\u05D8", 2, 1),
        Among("\u05E9\u05D0\u05E4\u05D8", 2, 1),
        Among("\u05D4\u05F2\u05D8", 2, 1),
        Among("\u05E7\u05F2\u05D8", 2, 1),
        Among("\u05D9\u05E7\u05F2\u05D8", 9, 1),
        Among("\u05DC\u05E2\u05DB", -1, 1),
        Among("\u05E2\u05DC\u05E2\u05DB", 11, 1),
        Among("\u05D9\u05D6\u05DE", -1, 1),
        Among("\u05D9\u05DE", -1, 1),
        Among("\u05E2\u05DE", -1, 1),
        Among("\u05E2\u05E0\u05E2\u05DE", 15, 3),
        Among("\u05D8\u05E2\u05E0\u05E2\u05DE", 16, 4),
        Among("\u05E0", -1, 1),
        Among("\u05E7\u05DC\u05D9\u05D1\u05E0", 18, 14),
        Among("\u05E8\u05D9\u05D1\u05E0", 18, 15),
        Among("\u05D8\u05E8\u05D9\u05D1\u05E0", 20, 12),
        Among("\u05E9\u05E8\u05D9\u05D1\u05E0", 20, 7),
        Among("\u05D4\u05F1\u05D1\u05E0", 18, 27),
        Among("\u05E9\u05F0\u05D9\u05D2\u05E0", 18, 17),
        Among("\u05D6\u05D5\u05E0\u05D2\u05E0", 18, 22),
        Among("\u05E9\u05DC\u05D5\u05E0\u05D2\u05E0", 18, 25),
        Among("\u05E6\u05F0\u05D5\u05E0\u05D2\u05E0", 18, 24),
        Among("\u05D1\u05F1\u05D2\u05E0", 18, 26),
        Among("\u05D1\u05D5\u05E0\u05D3\u05E0", 18, 20),
        Among("\u05F0\u05D9\u05D6\u05E0", 18, 11),
        Among("\u05D8\u05E0", 18, 4),
        Among("GE\u05D1\u05D9\u05D8\u05E0", 31, 9),
        Among("GE\u05DC\u05D9\u05D8\u05E0", 31, 13),
        Among("GE\u05DE\u05D9\u05D8\u05E0", 31, 8),
        Among("\u05E9\u05E0\u05D9\u05D8\u05E0", 31, 19),
        Among("\u05E1\u05D8\u05E0", 31, 1),
        Among("\u05D9\u05E1\u05D8\u05E0", 36, 1),
        Among("\u05E2\u05D8\u05E0", 31, 1),
        Among("GE\u05D1\u05D9\u05E1\u05E0", 18, 10),
        Among("\u05E9\u05DE\u05D9\u05E1\u05E0", 18, 18),
        Among("GE\u05E8\u05D9\u05E1\u05E0", 18, 16),
        Among("\u05E2\u05E0", 18, 1),
        Among("\u05D2\u05D0\u05E0\u05D2\u05E2\u05E0", 42, 5),
        Among("\u05E2\u05DC\u05E2\u05E0", 42, 1),
        Among("\u05E0\u05D5\u05DE\u05E2\u05E0", 42, 6),
        Among("\u05D9\u05D6\u05DE\u05E2\u05E0", 42, 1),
        Among("\u05E9\u05D8\u05D0\u05E0\u05E2\u05E0", 42, 29),
        Among("\u05D8\u05E8\u05D5\u05E0\u05E7\u05E0", 18, 23),
        Among("\u05E4\u05D0\u05E8\u05DC\u05F1\u05E8\u05E0", 18, 28),
        Among("\u05E9\u05F0\u05F1\u05E8\u05E0", 18, 30),
        Among("\u05F0\u05D5\u05D8\u05E9\u05E0", 18, 21),
        Among("\u05D2\u05F2\u05E0", 18, 5),
        Among("\u05E1", -1, 1),
        Among("\u05D8\u05E1", 53, 4),
        Among("\u05E2\u05D8\u05E1", 54, 1),
        Among("\u05E0\u05E1", 53, 1),
        Among("\u05D8\u05E0\u05E1", 56, 4),
        Among("\u05E2\u05E0\u05E1", 56, 3),
        Among("\u05E2\u05E1", 53, 1),
        Among("\u05D9\u05E2\u05E1", 59, 2),
        Among("\u05E2\u05DC\u05E2\u05E1", 59, 1),
        Among("\u05E2\u05E8\u05E1", 53, 1),
        Among("\u05E2\u05E0\u05E2\u05E8\u05E1", 62, 1),
        Among("\u05E2", -1, 1),
        Among("\u05D8\u05E2", 64, 4),
        Among("\u05E1\u05D8\u05E2", 65, 1),
        Among("\u05E2\u05D8\u05E2", 65, 1),
        Among("\u05D9\u05E2", 64, -1),
        Among("\u05E2\u05DC\u05E2", 64, 1),
        Among("\u05E2\u05E0\u05E2", 64, 3),
        Among("\u05D8\u05E2\u05E0\u05E2", 70, 4),
        Among("\u05E2\u05E8", -1, 1),
        Among("\u05D8\u05E2\u05E8", 72, 4),
        Among("\u05E1\u05D8\u05E2\u05E8", 73, 1),
        Among("\u05E2\u05D8\u05E2\u05E8", 73, 1),
        Among("\u05E2\u05E0\u05E2\u05E8", 72, 3),
        Among("\u05D8\u05E2\u05E0\u05E2\u05E8", 76, 4),
        Among("\u05D5\u05EA", -1, 32)
    ]

    a_5 = [
        Among("\u05D5\u05E0\u05D2", -1, 1),
        Among("\u05E9\u05D0\u05E4\u05D8", -1, 1),
        Among("\u05D4\u05F2\u05D8", -1, 1),
        Among("\u05E7\u05F2\u05D8", -1, 1),
        Among("\u05D9\u05E7\u05F2\u05D8", 3, 1),
        Among("\u05DC", -1, 2)
    ]

    a_6 = [
        Among("\u05D9\u05D2", -1, 1),
        Among("\u05D9\u05E7", -1, 1),
        Among("\u05D3\u05D9\u05E7", 1, 1),
        Among("\u05E0\u05D3\u05D9\u05E7", 2, 1),
        Among("\u05E2\u05E0\u05D3\u05D9\u05E7", 3, 1),
        Among("\u05D1\u05DC\u05D9\u05E7", 1, -1),
        Among("\u05D2\u05DC\u05D9\u05E7", 1, -1),
        Among("\u05E0\u05D9\u05E7", 1, 1),
        Among("\u05D9\u05E9", -1, 1)
    ]


class lab0(BaseException): pass


class lab1(BaseException): pass


class lab2(BaseException): pass


class lab3(BaseException): pass
