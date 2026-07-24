# Generated from indonesian.sbl by Snowball 3.1.1 - https://snowballstem.org/

from .basestemmer import BaseStemmer
from .among import Among


class IndonesianStemmer(BaseStemmer):
    '''
    This class implements the stemming algorithm defined by a snowball script.
    Generated from indonesian.sbl by Snowball 3.1.1 - https://snowballstem.org/
    '''

    g_vowel = {"a", "e", "i", "o", "u"}

    I_prefix = 0
    I_measure = 0

    def __r_remove_particle(self):
        self.ket = self.cursor
        if self.find_among_b(IndonesianStemmer.a_0) == 0:
            return False
        self.bra = self.cursor
        self.slice_del()
        self.I_measure -= 1
        return True

    def __r_remove_possessive_pronoun(self):
        self.ket = self.cursor
        if self.find_among_b(IndonesianStemmer.a_1) == 0:
            return False
        self.bra = self.cursor
        self.slice_del()
        self.I_measure -= 1
        return True

    def __r_remove_suffix(self):
        self.ket = self.cursor
        among_var = self.find_among_b(IndonesianStemmer.a_2)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var == 1:
            while True:
                v_1 = self.limit - self.cursor
                try:
                    if self.I_prefix == 3:
                        raise lab0()
                    if self.I_prefix == 2:
                        raise lab0()
                    if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "k":
                        raise lab0()
                    self.cursor -= 1
                    self.bra = self.cursor
                    break
                except lab0: pass
                self.cursor = self.limit - v_1
                if self.I_prefix == 1:
                    return False
                break
        else:
            if self.I_prefix > 2:
                return False
            try:
                if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "s":
                    raise lab0()
                self.cursor -= 1
                return False
            except lab0: pass
        self.slice_del()
        self.I_measure -= 1
        return True

    def __r_remove_first_order_prefix(self):
        self.bra = self.cursor
        among_var = self.find_among(IndonesianStemmer.a_3)
        if among_var == 0:
            return False
        self.ket = self.cursor
        if among_var == 1:
            self.slice_del()
            self.I_prefix = 1
            self.I_measure -= 1
        elif among_var == 2:
            while True:
                v_1 = self.cursor
                try:
                    if self.cursor == self.limit or self.current[self.cursor] != "y":
                        raise lab0()
                    self.cursor += 1
                    v_2 = self.cursor
                    if not self.in_grouping(IndonesianStemmer.g_vowel):
                        raise lab0()
                    self.cursor = v_2
                    self.ket = self.cursor
                    self.slice_from("s")
                    self.I_prefix = 1
                    self.I_measure -= 1
                    break
                except lab0: pass
                self.cursor = v_1
                self.slice_del()
                self.I_prefix = 1
                self.I_measure -= 1
                break
        elif among_var == 3:
            self.slice_del()
            self.I_prefix = 3
            self.I_measure -= 1
        elif among_var == 4:
            while True:
                v_3 = self.cursor
                try:
                    if self.cursor == self.limit or self.current[self.cursor] != "y":
                        raise lab0()
                    self.cursor += 1
                    v_4 = self.cursor
                    if not self.in_grouping(IndonesianStemmer.g_vowel):
                        raise lab0()
                    self.cursor = v_4
                    self.ket = self.cursor
                    self.slice_from("s")
                    self.I_prefix = 3
                    self.I_measure -= 1
                    break
                except lab0: pass
                self.cursor = v_3
                self.slice_del()
                self.I_prefix = 3
                self.I_measure -= 1
                break
        elif among_var == 5:
            self.I_prefix = 1
            self.I_measure -= 1
            while True:
                v_5 = self.cursor
                try:
                    v_6 = self.cursor
                    if not self.in_grouping(IndonesianStemmer.g_vowel):
                        raise lab0()
                    self.cursor = v_6
                    self.slice_from("p")
                    break
                except lab0: pass
                self.cursor = v_5
                self.slice_del()
                break
        else:
            self.I_prefix = 3
            self.I_measure -= 1
            while True:
                v_7 = self.cursor
                try:
                    v_8 = self.cursor
                    if not self.in_grouping(IndonesianStemmer.g_vowel):
                        raise lab0()
                    self.cursor = v_8
                    self.slice_from("p")
                    break
                except lab0: pass
                self.cursor = v_7
                self.slice_del()
                break
        return True

    def __r_remove_second_order_prefix(self):
        self.bra = self.cursor
        among_var = self.find_among(IndonesianStemmer.a_4)
        if among_var == 0:
            return False
        if among_var == 1:
            while True:
                v_1 = self.cursor
                try:
                    if self.cursor == self.limit or self.current[self.cursor] != "r":
                        raise lab0()
                    self.cursor += 1
                    self.ket = self.cursor
                    self.I_prefix = 2
                    break
                except lab0: pass
                self.cursor = v_1
                try:
                    if self.cursor == self.limit or self.current[self.cursor] != "l":
                        raise lab0()
                    self.cursor += 1
                    self.ket = self.cursor
                    if not self.eq_s("ajar"):
                        raise lab0()
                    break
                except lab0: pass
                self.cursor = v_1
                self.ket = self.cursor
                self.I_prefix = 2
                break
        else:
            while True:
                v_2 = self.cursor
                try:
                    if self.cursor == self.limit or self.current[self.cursor] != "r":
                        raise lab0()
                    self.cursor += 1
                    self.ket = self.cursor
                    break
                except lab0: pass
                self.cursor = v_2
                try:
                    if self.cursor == self.limit or self.current[self.cursor] != "l":
                        raise lab0()
                    self.cursor += 1
                    self.ket = self.cursor
                    if not self.eq_s("ajar"):
                        raise lab0()
                    break
                except lab0: pass
                self.cursor = v_2
                self.ket = self.cursor
                if not self.out_grouping(IndonesianStemmer.g_vowel):
                    return False
                if not self.eq_s("er"):
                    return False
                break
            self.I_prefix = 4
        self.I_measure -= 1
        self.slice_del()
        return True

    def _stem(self):
        self.I_measure = 0
        v_1 = self.cursor
        try:
            while True:
                v_2 = self.cursor
                try:
                    if not self.go_out_grouping(IndonesianStemmer.g_vowel):
                        raise lab1()
                    self.cursor += 1
                    self.I_measure += 1
                    continue
                except lab1: pass
                self.cursor = v_2
                break
        except lab0: pass
        self.cursor = v_1
        if self.I_measure <= 2:
            return False
        self.I_prefix = 0
        self.limit_backward = self.cursor
        self.cursor = self.limit
        v_3 = self.limit - self.cursor
        self.__r_remove_particle()
        self.cursor = self.limit - v_3
        if self.I_measure <= 2:
            return False
        v_4 = self.limit - self.cursor
        self.__r_remove_possessive_pronoun()
        self.cursor = self.limit - v_4
        self.cursor = self.limit_backward
        if self.I_measure <= 2:
            return False
        while True:
            v_5 = self.cursor
            try:
                v_6 = self.cursor
                if not self.__r_remove_first_order_prefix():
                    raise lab0()
                v_7 = self.cursor
                try:
                    v_8 = self.cursor
                    if self.I_measure <= 2:
                        raise lab1()
                    self.limit_backward = self.cursor
                    self.cursor = self.limit
                    if not self.__r_remove_suffix():
                        raise lab1()
                    self.cursor = self.limit_backward
                    self.cursor = v_8
                    if self.I_measure <= 2:
                        raise lab1()
                    if not self.__r_remove_second_order_prefix():
                        raise lab1()
                except lab1: pass
                self.cursor = v_7
                self.cursor = v_6
                break
            except lab0: pass
            self.cursor = v_5
            v_9 = self.cursor
            self.__r_remove_second_order_prefix()
            self.cursor = v_9
            v_10 = self.cursor
            try:
                if self.I_measure <= 2:
                    raise lab0()
                self.limit_backward = self.cursor
                self.cursor = self.limit
                if not self.__r_remove_suffix():
                    raise lab0()
                self.cursor = self.limit_backward
            except lab0: pass
            self.cursor = v_10
            break
        return True

    a_0 = [
        Among("kah", -1, 1),
        Among("lah", -1, 1),
        Among("pun", -1, 1)
    ]

    a_1 = [
        Among("nya", -1, 1),
        Among("ku", -1, 1),
        Among("mu", -1, 1)
    ]

    a_2 = [
        Among("i", -1, 2),
        Among("an", -1, 1)
    ]

    a_3 = [
        Among("di", -1, 1),
        Among("ke", -1, 3),
        Among("me", -1, 1),
        Among("mem", 2, 5),
        Among("men", 2, 2),
        Among("meng", 4, 1),
        Among("pem", -1, 6),
        Among("pen", -1, 4),
        Among("peng", 7, 3),
        Among("ter", -1, 1)
    ]

    a_4 = [
        Among("be", -1, 2),
        Among("pe", -1, 1)
    ]


class lab0(BaseException): pass


class lab1(BaseException): pass
