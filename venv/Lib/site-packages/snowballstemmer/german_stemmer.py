# Generated from german.sbl by Snowball 3.1.1 - https://snowballstem.org/

from .basestemmer import BaseStemmer
from .among import Among


class GermanStemmer(BaseStemmer):
    '''
    This class implements the stemming algorithm defined by a snowball script.
    Generated from german.sbl by Snowball 3.1.1 - https://snowballstem.org/
    '''

    g_v = {"a", "e", "i", "o", "u", "y", "ä", "ö", "ü"}

    g_et_ending = {"U", "d", "f", "g", "k", "l", "m", "n", "r", "s", "t", "z", "ä"}

    g_s_ending = {"b", "d", "f", "g", "h", "k", "l", "m", "n", "r", "t"}

    g_st_ending = {"b", "d", "f", "g", "h", "k", "l", "m", "n", "t"}

    I_p2 = 0
    I_p1 = 0

    def __r_prelude(self):
        v_1 = self.cursor
        while True:
            v_2 = self.cursor
            try:
                while True:
                    v_3 = self.cursor
                    try:
                        if not self.in_grouping(GermanStemmer.g_v):
                            raise lab1()
                        self.bra = self.cursor
                        while True:
                            v_4 = self.cursor
                            try:
                                if self.cursor == self.limit or self.current[self.cursor] != "u":
                                    raise lab2()
                                self.cursor += 1
                                self.ket = self.cursor
                                if not self.in_grouping(GermanStemmer.g_v):
                                    raise lab2()
                                self.slice_from("U")
                                break
                            except lab2: pass
                            self.cursor = v_4
                            if self.cursor == self.limit or self.current[self.cursor] != "y":
                                raise lab1()
                            self.cursor += 1
                            self.ket = self.cursor
                            if not self.in_grouping(GermanStemmer.g_v):
                                raise lab1()
                            self.slice_from("Y")
                            break
                        self.cursor = v_3
                        break
                    except lab1: pass
                    self.cursor = v_3
                    if self.cursor >= self.limit:
                        raise lab0()
                    self.cursor += 1
                continue
            except lab0: pass
            self.cursor = v_2
            break
        self.cursor = v_1
        while True:
            v_5 = self.cursor
            try:
                self.bra = self.cursor
                among_var = self.find_among(GermanStemmer.a_0)
                self.ket = self.cursor
                if among_var == 1:
                    self.slice_from("ss")
                elif among_var == 2:
                    self.slice_from("ä")
                elif among_var == 3:
                    self.slice_from("ö")
                elif among_var == 4:
                    self.slice_from("ü")
                elif among_var == 5:
                    if self.cursor >= self.limit:
                        raise lab0()
                    self.cursor += 1
                continue
            except lab0: pass
            self.cursor = v_5
            break
        return True

    def __r_mark_regions(self):
        self.I_p1 = self.limit
        self.I_p2 = self.limit
        v_1 = self.cursor
        if self.cursor + 3 > self.limit:
            return False
        self.cursor += 3
        I_x = self.cursor
        self.cursor = v_1
        if not self.go_out_grouping(GermanStemmer.g_v):
            return False
        self.cursor += 1
        if not self.go_in_grouping(GermanStemmer.g_v):
            return False
        self.cursor += 1
        self.I_p1 = self.cursor
        try:
            if self.I_p1 >= I_x:
                raise lab0()
            self.I_p1 = I_x
        except lab0: pass
        if not self.go_out_grouping(GermanStemmer.g_v):
            return False
        self.cursor += 1
        if not self.go_in_grouping(GermanStemmer.g_v):
            return False
        self.cursor += 1
        self.I_p2 = self.cursor
        return True

    def __r_postlude(self):
        while True:
            v_1 = self.cursor
            try:
                self.bra = self.cursor
                among_var = self.find_among(GermanStemmer.a_1)
                self.ket = self.cursor
                if among_var == 1:
                    self.slice_from("y")
                elif among_var == 2:
                    self.slice_from("u")
                elif among_var == 3:
                    self.slice_from("a")
                elif among_var == 4:
                    self.slice_from("o")
                else:
                    if self.cursor >= self.limit:
                        raise lab0()
                    self.cursor += 1
                continue
            except lab0: pass
            self.cursor = v_1
            break
        return True

    def __r_R1(self):
        return self.I_p1 <= self.cursor

    def __r_R2(self):
        return self.I_p2 <= self.cursor

    def __r_standard_suffix(self):
        v_1 = self.limit - self.cursor
        try:
            self.ket = self.cursor
            among_var = self.find_among_b(GermanStemmer.a_2)
            if among_var == 0:
                raise lab0()
            self.bra = self.cursor
            if not self.__r_R1():
                raise lab0()
            if among_var == 1:
                try:
                    if not self.eq_s_b("syst"):
                        raise lab1()
                    raise lab0()
                except lab1: pass
                self.slice_del()
            elif among_var == 2:
                self.slice_del()
            elif among_var == 3:
                self.slice_del()
                v_2 = self.limit - self.cursor
                try:
                    self.ket = self.cursor
                    if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "s":
                        self.cursor = self.limit - v_2
                        raise lab1()
                    self.cursor -= 1
                    self.bra = self.cursor
                    if not self.eq_s_b("nis"):
                        self.cursor = self.limit - v_2
                        raise lab1()
                    self.slice_del()
                except lab1: pass
            elif among_var == 4:
                if not self.in_grouping_b(GermanStemmer.g_s_ending):
                    raise lab0()
                self.slice_del()
            else:
                self.slice_from("l")
        except lab0: pass
        self.cursor = self.limit - v_1
        v_3 = self.limit - self.cursor
        try:
            self.ket = self.cursor
            among_var = self.find_among_b(GermanStemmer.a_4)
            if among_var == 0:
                raise lab0()
            self.bra = self.cursor
            if not self.__r_R1():
                raise lab0()
            if among_var == 1:
                self.slice_del()
            elif among_var == 2:
                if not self.in_grouping_b(GermanStemmer.g_st_ending):
                    raise lab0()
                if self.cursor - 3 < self.limit_backward:
                    raise lab0()
                self.cursor -= 3
                self.slice_del()
            else:
                v_4 = self.limit - self.cursor
                if not self.in_grouping_b(GermanStemmer.g_et_ending):
                    raise lab0()
                self.cursor = self.limit - v_4
                v_5 = self.limit - self.cursor
                try:
                    if self.find_among_b(GermanStemmer.a_3) == 0:
                        raise lab1()
                    raise lab0()
                except lab1: pass
                self.cursor = self.limit - v_5
                self.slice_del()
        except lab0: pass
        self.cursor = self.limit - v_3
        v_6 = self.limit - self.cursor
        try:
            self.ket = self.cursor
            among_var = self.find_among_b(GermanStemmer.a_6)
            if among_var == 0:
                raise lab0()
            self.bra = self.cursor
            if not self.__r_R2():
                raise lab0()
            if among_var == 1:
                self.slice_del()
                v_7 = self.limit - self.cursor
                try:
                    self.ket = self.cursor
                    if not self.eq_s_b("ig"):
                        self.cursor = self.limit - v_7
                        raise lab1()
                    self.bra = self.cursor
                    try:
                        if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "e":
                            raise lab2()
                        self.cursor -= 1
                        self.cursor = self.limit - v_7
                        raise lab1()
                    except lab2: pass
                    if not self.__r_R2():
                        self.cursor = self.limit - v_7
                        raise lab1()
                    self.slice_del()
                except lab1: pass
            elif among_var == 2:
                try:
                    if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "e":
                        raise lab1()
                    self.cursor -= 1
                    raise lab0()
                except lab1: pass
                self.slice_del()
            elif among_var == 3:
                self.slice_del()
                v_8 = self.limit - self.cursor
                try:
                    self.ket = self.cursor
                    while True:
                        try:
                            if not self.eq_s_b("er"):
                                raise lab2()
                            break
                        except lab2: pass
                        if not self.eq_s_b("en"):
                            self.cursor = self.limit - v_8
                            raise lab1()
                        break
                    self.bra = self.cursor
                    if not self.__r_R1():
                        self.cursor = self.limit - v_8
                        raise lab1()
                    self.slice_del()
                except lab1: pass
            else:
                self.slice_del()
                v_9 = self.limit - self.cursor
                try:
                    self.ket = self.cursor
                    if self.find_among_b(GermanStemmer.a_5) == 0:
                        self.cursor = self.limit - v_9
                        raise lab1()
                    self.bra = self.cursor
                    if not self.__r_R2():
                        self.cursor = self.limit - v_9
                        raise lab1()
                    self.slice_del()
                except lab1: pass
        except lab0: pass
        self.cursor = self.limit - v_6
        v_10 = self.limit - self.cursor
        try:
            self.ket = self.cursor
            if self.find_among_b(GermanStemmer.a_7) == 0:
                raise lab0()
            self.bra = self.cursor
            if self.cursor <= self.limit_backward:
                raise lab0()
            self.cursor -= 1
            if self.cursor <= self.limit_backward:
                raise lab0()
            self.slice_del()
        except lab0: pass
        self.cursor = self.limit - v_10
        return True

    def _stem(self):
        v_1 = self.cursor
        self.__r_prelude()
        self.cursor = v_1
        v_2 = self.cursor
        self.__r_mark_regions()
        self.cursor = v_2
        self.limit_backward = self.cursor
        self.cursor = self.limit
        self.__r_standard_suffix()
        self.cursor = self.limit_backward
        v_3 = self.cursor
        self.__r_postlude()
        self.cursor = v_3
        return True

    a_0 = [
        Among("", -1, 5),
        Among("ae", 0, 2),
        Among("oe", 0, 3),
        Among("qu", 0, -1),
        Among("ue", 0, 4),
        Among("ß", 0, 1)
    ]

    a_1 = [
        Among("", -1, 5),
        Among("U", 0, 2),
        Among("Y", 0, 1),
        Among("ä", 0, 3),
        Among("ö", 0, 4),
        Among("ü", 0, 2)
    ]

    a_2 = [
        Among("e", -1, 3),
        Among("em", -1, 1),
        Among("en", -1, 3),
        Among("erinnen", 2, 2),
        Among("erin", -1, 2),
        Among("ln", -1, 5),
        Among("ern", -1, 2),
        Among("er", -1, 2),
        Among("s", -1, 4),
        Among("es", 8, 3),
        Among("lns", 8, 5)
    ]

    a_3 = [
        Among("tick", -1, -1),
        Among("plan", -1, -1),
        Among("geordn", -1, -1),
        Among("intern", -1, -1),
        Among("tr", -1, -1)
    ]

    a_4 = [
        Among("en", -1, 1),
        Among("er", -1, 1),
        Among("et", -1, 3),
        Among("st", -1, 2),
        Among("est", 3, 1)
    ]

    a_5 = [
        Among("ig", -1, 1),
        Among("lich", -1, 1)
    ]

    a_6 = [
        Among("end", -1, 1),
        Among("ig", -1, 2),
        Among("ung", -1, 1),
        Among("lich", -1, 3),
        Among("isch", -1, 2),
        Among("ik", -1, 2),
        Among("heit", -1, 3),
        Among("keit", -1, 4)
    ]

    a_7 = [
        Among("'", -1, 1),
        Among("'sch", -1, 1),
        Among("'s", -1, 1)
    ]


class lab0(BaseException): pass


class lab1(BaseException): pass


class lab2(BaseException): pass
