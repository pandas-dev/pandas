# Generated from porter.sbl by Snowball 3.1.1 - https://snowballstem.org/

from .basestemmer import BaseStemmer
from .among import Among


class PorterStemmer(BaseStemmer):
    '''
    This class implements the stemming algorithm defined by a snowball script.
    Generated from porter.sbl by Snowball 3.1.1 - https://snowballstem.org/
    '''

    g_v = {"a", "e", "i", "o", "u", "y"}

    g_v_WXY = {"Y", "a", "e", "i", "o", "u", "w", "x", "y"}

    I_p2 = 0
    I_p1 = 0

    def __r_shortv(self):
        if not self.out_grouping_b(PorterStemmer.g_v_WXY):
            return False
        if not self.in_grouping_b(PorterStemmer.g_v):
            return False
        return self.out_grouping_b(PorterStemmer.g_v)

    def __r_R1(self):
        return self.I_p1 <= self.cursor

    def __r_R2(self):
        return self.I_p2 <= self.cursor

    def __r_Step_1a(self):
        self.ket = self.cursor
        among_var = self.find_among_b(PorterStemmer.a_0)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var > 0:
            self.slice_from(PorterStemmer.as_0[among_var - 1])
        return True

    def __r_Step_1b(self):
        self.ket = self.cursor
        among_var = self.find_among_b(PorterStemmer.a_2)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var == 1:
            if not self.__r_R1():
                return False
            self.slice_from("ee")
        else:
            v_1 = self.limit - self.cursor
            if not self.go_out_grouping_b(PorterStemmer.g_v):
                return False
            self.cursor -= 1
            self.cursor = self.limit - v_1
            self.slice_del()
            v_2 = self.limit - self.cursor
            among_var = self.find_among_b(PorterStemmer.a_1)
            self.cursor = self.limit - v_2
            if among_var == 1:
                c = self.cursor
                self.insert(self.cursor, self.cursor, "e")
                self.cursor = c
            elif among_var == 2:
                self.ket = self.cursor
                if self.cursor <= self.limit_backward:
                    return False
                self.cursor -= 1
                self.bra = self.cursor
                self.slice_del()
            else:
                if self.cursor != self.I_p1:
                    return False
                v_3 = self.limit - self.cursor
                if not self.__r_shortv():
                    return False
                self.cursor = self.limit - v_3
                c = self.cursor
                self.insert(self.cursor, self.cursor, "e")
                self.cursor = c
        return True

    def __r_Step_1c(self):
        self.ket = self.cursor
        while True:
            try:
                if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "y":
                    raise lab0()
                self.cursor -= 1
                break
            except lab0: pass
            if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "Y":
                return False
            self.cursor -= 1
            break
        self.bra = self.cursor
        if not self.go_out_grouping_b(PorterStemmer.g_v):
            return False
        self.cursor -= 1
        self.slice_from("i")
        return True

    def __r_Step_2(self):
        self.ket = self.cursor
        among_var = self.find_among_b(PorterStemmer.a_3)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if not self.__r_R1():
            return False
        self.slice_from(PorterStemmer.as_3[among_var - 1])
        return True

    def __r_Step_3(self):
        self.ket = self.cursor
        among_var = self.find_among_b(PorterStemmer.a_4)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if not self.__r_R1():
            return False
        self.slice_from(PorterStemmer.as_4[among_var - 1])
        return True

    def __r_Step_4(self):
        self.ket = self.cursor
        among_var = self.find_among_b(PorterStemmer.a_5)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if not self.__r_R2():
            return False
        if among_var == 1:
            self.slice_del()
        else:
            while True:
                try:
                    if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "s":
                        raise lab0()
                    self.cursor -= 1
                    break
                except lab0: pass
                if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "t":
                    return False
                self.cursor -= 1
                break
            self.slice_del()
        return True

    def __r_Step_5a(self):
        self.ket = self.cursor
        if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "e":
            return False
        self.cursor -= 1
        self.bra = self.cursor
        while True:
            try:
                if not self.__r_R2():
                    raise lab0()
                break
            except lab0: pass
            if not self.__r_R1():
                return False
            v_1 = self.limit - self.cursor
            try:
                if not self.__r_shortv():
                    raise lab0()
                return False
            except lab0: pass
            self.cursor = self.limit - v_1
            break
        self.slice_del()
        return True

    def __r_Step_5b(self):
        self.ket = self.cursor
        if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "l":
            return False
        self.cursor -= 1
        self.bra = self.cursor
        if not self.__r_R2():
            return False
        if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "l":
            return False
        self.cursor -= 1
        self.slice_del()
        return True

    def _stem(self):
        B_Y_found = False
        v_1 = self.cursor
        try:
            self.bra = self.cursor
            if self.cursor == self.limit or self.current[self.cursor] != "y":
                raise lab0()
            self.cursor += 1
            self.ket = self.cursor
            self.slice_from("Y")
            B_Y_found = True
        except lab0: pass
        self.cursor = v_1
        v_2 = self.cursor
        try:
            while True:
                v_3 = self.cursor
                try:
                    while True:
                        v_4 = self.cursor
                        try:
                            if not self.in_grouping(PorterStemmer.g_v):
                                raise lab2()
                            self.bra = self.cursor
                            if self.cursor == self.limit or self.current[self.cursor] != "y":
                                raise lab2()
                            self.cursor += 1
                            self.ket = self.cursor
                            self.cursor = v_4
                            break
                        except lab2: pass
                        self.cursor = v_4
                        if self.cursor >= self.limit:
                            raise lab1()
                        self.cursor += 1
                    self.slice_from("Y")
                    B_Y_found = True
                    continue
                except lab1: pass
                self.cursor = v_3
                break
        except lab0: pass
        self.cursor = v_2
        self.I_p1 = self.limit
        self.I_p2 = self.limit
        v_5 = self.cursor
        try:
            if not self.go_out_grouping(PorterStemmer.g_v):
                raise lab0()
            self.cursor += 1
            if not self.go_in_grouping(PorterStemmer.g_v):
                raise lab0()
            self.cursor += 1
            self.I_p1 = self.cursor
            if not self.go_out_grouping(PorterStemmer.g_v):
                raise lab0()
            self.cursor += 1
            if not self.go_in_grouping(PorterStemmer.g_v):
                raise lab0()
            self.cursor += 1
            self.I_p2 = self.cursor
        except lab0: pass
        self.cursor = v_5
        self.limit_backward = self.cursor
        self.cursor = self.limit
        v_6 = self.limit - self.cursor
        self.__r_Step_1a()
        self.cursor = self.limit - v_6
        v_7 = self.limit - self.cursor
        self.__r_Step_1b()
        self.cursor = self.limit - v_7
        v_8 = self.limit - self.cursor
        self.__r_Step_1c()
        self.cursor = self.limit - v_8
        v_9 = self.limit - self.cursor
        self.__r_Step_2()
        self.cursor = self.limit - v_9
        v_10 = self.limit - self.cursor
        self.__r_Step_3()
        self.cursor = self.limit - v_10
        v_11 = self.limit - self.cursor
        self.__r_Step_4()
        self.cursor = self.limit - v_11
        v_12 = self.limit - self.cursor
        self.__r_Step_5a()
        self.cursor = self.limit - v_12
        v_13 = self.limit - self.cursor
        self.__r_Step_5b()
        self.cursor = self.limit - v_13
        self.cursor = self.limit_backward
        v_14 = self.cursor
        try:
            if not B_Y_found:
                raise lab0()
            while True:
                v_15 = self.cursor
                try:
                    while True:
                        v_16 = self.cursor
                        try:
                            self.bra = self.cursor
                            if self.cursor == self.limit or self.current[self.cursor] != "Y":
                                raise lab2()
                            self.cursor += 1
                            self.ket = self.cursor
                            self.cursor = v_16
                            break
                        except lab2: pass
                        self.cursor = v_16
                        if self.cursor >= self.limit:
                            raise lab1()
                        self.cursor += 1
                    self.slice_from("y")
                    continue
                except lab1: pass
                self.cursor = v_15
                break
        except lab0: pass
        self.cursor = v_14
        return True

    a_0 = [
        Among("s", -1, 3),
        Among("ies", 0, 2),
        Among("sses", 0, 1),
        Among("ss", 0, -1)
    ]
    as_0 = ("ss", "i", "")

    a_1 = [
        Among("", -1, 3),
        Among("bb", 0, 2),
        Among("dd", 0, 2),
        Among("ff", 0, 2),
        Among("gg", 0, 2),
        Among("bl", 0, 1),
        Among("mm", 0, 2),
        Among("nn", 0, 2),
        Among("pp", 0, 2),
        Among("rr", 0, 2),
        Among("at", 0, 1),
        Among("tt", 0, 2),
        Among("iz", 0, 1)
    ]

    a_2 = [
        Among("ed", -1, 2),
        Among("eed", 0, 1),
        Among("ing", -1, 2)
    ]

    a_3 = [
        Among("anci", -1, 3),
        Among("enci", -1, 2),
        Among("abli", -1, 4),
        Among("eli", -1, 6),
        Among("alli", -1, 9),
        Among("ousli", -1, 11),
        Among("entli", -1, 5),
        Among("aliti", -1, 9),
        Among("biliti", -1, 13),
        Among("iviti", -1, 12),
        Among("tional", -1, 1),
        Among("ational", 10, 8),
        Among("alism", -1, 9),
        Among("ation", -1, 8),
        Among("ization", 13, 7),
        Among("izer", -1, 7),
        Among("ator", -1, 8),
        Among("iveness", -1, 12),
        Among("fulness", -1, 10),
        Among("ousness", -1, 11)
    ]
    as_3 = ("tion", "ence", "ance", "able", "ent", "e", "ize", "ate", "al", "ful", "ous", "ive", "ble")

    a_4 = [
        Among("icate", -1, 2),
        Among("ative", -1, 3),
        Among("alize", -1, 1),
        Among("iciti", -1, 2),
        Among("ical", -1, 2),
        Among("ful", -1, 3),
        Among("ness", -1, 3)
    ]
    as_4 = ("al", "ic", "")

    a_5 = [
        Among("ic", -1, 1),
        Among("ance", -1, 1),
        Among("ence", -1, 1),
        Among("able", -1, 1),
        Among("ible", -1, 1),
        Among("ate", -1, 1),
        Among("ive", -1, 1),
        Among("ize", -1, 1),
        Among("iti", -1, 1),
        Among("al", -1, 1),
        Among("ism", -1, 1),
        Among("ion", -1, 2),
        Among("er", -1, 1),
        Among("ous", -1, 1),
        Among("ant", -1, 1),
        Among("ent", -1, 1),
        Among("ment", 15, 1),
        Among("ement", 16, 1),
        Among("ou", -1, 1)
    ]


class lab0(BaseException): pass


class lab1(BaseException): pass


class lab2(BaseException): pass
