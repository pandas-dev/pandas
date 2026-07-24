# Generated from dutch_porter.sbl by Snowball 3.1.1 - https://snowballstem.org/

from .basestemmer import BaseStemmer
from .among import Among


class DutchPorterStemmer(BaseStemmer):
    '''
    This class implements the stemming algorithm defined by a snowball script.
    Generated from dutch_porter.sbl by Snowball 3.1.1 - https://snowballstem.org/
    '''

    g_v = {"a", "e", "i", "o", "u", "y", "è"}

    g_v_I = {"I", "a", "e", "i", "o", "u", "y", "è"}

    g_v_j = {"a", "e", "i", "j", "o", "u", "y", "è"}

    I_p2 = 0
    I_p1 = 0
    B_e_found = False

    def __r_prelude(self):
        v_1 = self.cursor
        while True:
            v_2 = self.cursor
            try:
                self.bra = self.cursor
                among_var = self.find_among(DutchPorterStemmer.a_0)
                self.ket = self.cursor
                if among_var == 1:
                    self.slice_from("a")
                elif among_var == 2:
                    self.slice_from("e")
                elif among_var == 3:
                    self.slice_from("i")
                elif among_var == 4:
                    self.slice_from("o")
                elif among_var == 5:
                    self.slice_from("u")
                else:
                    if self.cursor >= self.limit:
                        raise lab0()
                    self.cursor += 1
                continue
            except lab0: pass
            self.cursor = v_2
            break
        self.cursor = v_1
        v_3 = self.cursor
        try:
            self.bra = self.cursor
            if self.cursor == self.limit or self.current[self.cursor] != "y":
                self.cursor = v_3
                raise lab0()
            self.cursor += 1
            self.ket = self.cursor
            self.slice_from("Y")
        except lab0: pass
        while True:
            v_4 = self.cursor
            try:
                if not self.go_out_grouping(DutchPorterStemmer.g_v):
                    raise lab0()
                self.cursor += 1
                v_5 = self.cursor
                try:
                    self.bra = self.cursor
                    while True:
                        v_6 = self.cursor
                        try:
                            if self.cursor == self.limit or self.current[self.cursor] != "i":
                                raise lab2()
                            self.cursor += 1
                            self.ket = self.cursor
                            v_7 = self.cursor
                            try:
                                if not self.in_grouping(DutchPorterStemmer.g_v):
                                    raise lab3()
                                self.slice_from("I")
                            except lab3: pass
                            self.cursor = v_7
                            break
                        except lab2: pass
                        self.cursor = v_6
                        if self.cursor == self.limit or self.current[self.cursor] != "y":
                            self.cursor = v_5
                            raise lab1()
                        self.cursor += 1
                        self.ket = self.cursor
                        self.slice_from("Y")
                        break
                except lab1: pass
                continue
            except lab0: pass
            self.cursor = v_4
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
        if not self.go_out_grouping(DutchPorterStemmer.g_v):
            return False
        self.cursor += 1
        if not self.go_in_grouping(DutchPorterStemmer.g_v):
            return False
        self.cursor += 1
        self.I_p1 = self.cursor
        try:
            if self.I_p1 >= I_x:
                raise lab0()
            self.I_p1 = I_x
        except lab0: pass
        if not self.go_out_grouping(DutchPorterStemmer.g_v):
            return False
        self.cursor += 1
        if not self.go_in_grouping(DutchPorterStemmer.g_v):
            return False
        self.cursor += 1
        self.I_p2 = self.cursor
        return True

    def __r_postlude(self):
        while True:
            v_1 = self.cursor
            try:
                self.bra = self.cursor
                among_var = self.find_among(DutchPorterStemmer.a_1)
                self.ket = self.cursor
                if among_var == 1:
                    self.slice_from("y")
                elif among_var == 2:
                    self.slice_from("i")
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

    def __r_undouble(self):
        v_1 = self.limit - self.cursor
        if self.find_among_b(DutchPorterStemmer.a_2) == 0:
            return False
        self.cursor = self.limit - v_1
        self.ket = self.cursor
        if self.cursor <= self.limit_backward:
            return False
        self.cursor -= 1
        self.bra = self.cursor
        self.slice_del()
        return True

    def __r_e_ending(self):
        self.B_e_found = False
        self.ket = self.cursor
        if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "e":
            return False
        self.cursor -= 1
        self.bra = self.cursor
        if not self.__r_R1():
            return False
        v_1 = self.limit - self.cursor
        if not self.out_grouping_b(DutchPorterStemmer.g_v):
            return False
        self.cursor = self.limit - v_1
        self.slice_del()
        self.B_e_found = True
        return self.__r_undouble()

    def __r_en_ending(self):
        if not self.__r_R1():
            return False
        v_1 = self.limit - self.cursor
        if not self.out_grouping_b(DutchPorterStemmer.g_v):
            return False
        self.cursor = self.limit - v_1
        try:
            if not self.eq_s_b("gem"):
                raise lab0()
            return False
        except lab0: pass
        self.slice_del()
        return self.__r_undouble()

    def __r_standard_suffix(self):
        v_1 = self.limit - self.cursor
        try:
            self.ket = self.cursor
            among_var = self.find_among_b(DutchPorterStemmer.a_3)
            if among_var == 0:
                raise lab0()
            self.bra = self.cursor
            if among_var == 1:
                if not self.__r_R1():
                    raise lab0()
                self.slice_from("heid")
            elif among_var == 2:
                if not self.__r_en_ending():
                    raise lab0()
            else:
                if not self.__r_R1():
                    raise lab0()
                if not self.out_grouping_b(DutchPorterStemmer.g_v_j):
                    raise lab0()
                self.slice_del()
        except lab0: pass
        self.cursor = self.limit - v_1
        v_2 = self.limit - self.cursor
        self.__r_e_ending()
        self.cursor = self.limit - v_2
        v_3 = self.limit - self.cursor
        try:
            self.ket = self.cursor
            if not self.eq_s_b("heid"):
                raise lab0()
            self.bra = self.cursor
            if not self.__r_R2():
                raise lab0()
            try:
                if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "c":
                    raise lab1()
                self.cursor -= 1
                raise lab0()
            except lab1: pass
            self.slice_del()
            self.ket = self.cursor
            if not self.eq_s_b("en"):
                raise lab0()
            self.bra = self.cursor
            if not self.__r_en_ending():
                raise lab0()
        except lab0: pass
        self.cursor = self.limit - v_3
        v_4 = self.limit - self.cursor
        try:
            self.ket = self.cursor
            among_var = self.find_among_b(DutchPorterStemmer.a_4)
            if among_var == 0:
                raise lab0()
            self.bra = self.cursor
            if among_var == 1:
                if not self.__r_R2():
                    raise lab0()
                self.slice_del()
                while True:
                    v_5 = self.limit - self.cursor
                    try:
                        self.ket = self.cursor
                        if not self.eq_s_b("ig"):
                            raise lab1()
                        self.bra = self.cursor
                        if not self.__r_R2():
                            raise lab1()
                        try:
                            if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "e":
                                raise lab2()
                            self.cursor -= 1
                            raise lab1()
                        except lab2: pass
                        self.slice_del()
                        break
                    except lab1: pass
                    self.cursor = self.limit - v_5
                    if not self.__r_undouble():
                        raise lab0()
                    break
            elif among_var == 2:
                if not self.__r_R2():
                    raise lab0()
                try:
                    if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "e":
                        raise lab1()
                    self.cursor -= 1
                    raise lab0()
                except lab1: pass
                self.slice_del()
            elif among_var == 3:
                if not self.__r_R2():
                    raise lab0()
                self.slice_del()
                if not self.__r_e_ending():
                    raise lab0()
            elif among_var == 4:
                if not self.__r_R2():
                    raise lab0()
                self.slice_del()
            else:
                if not self.__r_R2():
                    raise lab0()
                if not self.B_e_found:
                    raise lab0()
                self.slice_del()
        except lab0: pass
        self.cursor = self.limit - v_4
        v_6 = self.limit - self.cursor
        try:
            if not self.out_grouping_b(DutchPorterStemmer.g_v_I):
                raise lab0()
            v_7 = self.limit - self.cursor
            if self.find_among_b(DutchPorterStemmer.a_5) == 0:
                raise lab0()
            if not self.out_grouping_b(DutchPorterStemmer.g_v):
                raise lab0()
            self.cursor = self.limit - v_7
            self.ket = self.cursor
            if self.cursor <= self.limit_backward:
                raise lab0()
            self.cursor -= 1
            self.bra = self.cursor
            self.slice_del()
        except lab0: pass
        self.cursor = self.limit - v_6
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
        Among("", -1, 6),
        Among("á", 0, 1),
        Among("ä", 0, 1),
        Among("é", 0, 2),
        Among("ë", 0, 2),
        Among("í", 0, 3),
        Among("ï", 0, 3),
        Among("ó", 0, 4),
        Among("ö", 0, 4),
        Among("ú", 0, 5),
        Among("ü", 0, 5)
    ]

    a_1 = [
        Among("", -1, 3),
        Among("I", 0, 2),
        Among("Y", 0, 1)
    ]

    a_2 = [
        Among("dd", -1, -1),
        Among("kk", -1, -1),
        Among("tt", -1, -1)
    ]

    a_3 = [
        Among("ene", -1, 2),
        Among("se", -1, 3),
        Among("en", -1, 2),
        Among("heden", 2, 1),
        Among("s", -1, 3)
    ]

    a_4 = [
        Among("end", -1, 1),
        Among("ig", -1, 2),
        Among("ing", -1, 1),
        Among("lijk", -1, 3),
        Among("baar", -1, 4),
        Among("bar", -1, 5)
    ]

    a_5 = [
        Among("aa", -1, -1),
        Among("ee", -1, -1),
        Among("oo", -1, -1),
        Among("uu", -1, -1)
    ]


class lab0(BaseException): pass


class lab1(BaseException): pass


class lab2(BaseException): pass


class lab3(BaseException): pass
