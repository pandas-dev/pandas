# Generated from esperanto.sbl by Snowball 3.1.1 - https://snowballstem.org/

from .basestemmer import BaseStemmer
from .among import Among


class EsperantoStemmer(BaseStemmer):
    '''
    This class implements the stemming algorithm defined by a snowball script.
    Generated from esperanto.sbl by Snowball 3.1.1 - https://snowballstem.org/
    '''

    g_vowel = {"a", "e", "i", "o", "u"}

    g_aou = "aou"
    g_digit = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9"}


    def __r_canonical_form(self):
        B_foreign = False
        while True:
            v_1 = self.cursor
            try:
                self.bra = self.cursor
                among_var = self.find_among(EsperantoStemmer.a_0)
                self.ket = self.cursor
                if among_var == 1:
                    self.slice_from("ĉ")
                elif among_var == 2:
                    self.slice_from("ĝ")
                elif among_var == 3:
                    self.slice_from("ĥ")
                elif among_var == 4:
                    self.slice_from("ĵ")
                elif among_var == 5:
                    self.slice_from("ŝ")
                elif among_var == 6:
                    self.slice_from("ŭ")
                elif among_var == 7:
                    self.slice_from("a")
                    B_foreign = True
                elif among_var == 8:
                    self.slice_from("e")
                    B_foreign = True
                elif among_var == 9:
                    self.slice_from("i")
                    B_foreign = True
                elif among_var == 10:
                    self.slice_from("o")
                    B_foreign = True
                elif among_var == 11:
                    self.slice_from("u")
                    B_foreign = True
                elif among_var == 12:
                    B_foreign = True
                elif among_var == 13:
                    B_foreign = False
                else:
                    if self.cursor >= self.limit:
                        raise lab0()
                    self.cursor += 1
                continue
            except lab0: pass
            self.cursor = v_1
            break
        return not B_foreign

    def __r_initial_apostrophe(self):
        self.bra = self.cursor
        if self.cursor == self.limit or self.current[self.cursor] != "'":
            return False
        self.cursor += 1
        self.ket = self.cursor
        if not self.eq_s("st"):
            return False
        if self.find_among(EsperantoStemmer.a_1) == 0:
            return False
        if self.cursor < self.limit:
            return False
        self.slice_from("e")
        return True

    def __r_pronoun(self):
        self.ket = self.cursor
        v_1 = self.limit - self.cursor
        try:
            if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "n":
                self.cursor = self.limit - v_1
                raise lab0()
            self.cursor -= 1
        except lab0: pass
        self.bra = self.cursor
        if self.find_among_b(EsperantoStemmer.a_2) == 0:
            return False
        while True:
            try:
                if self.cursor > self.limit_backward:
                    raise lab0()
                break
            except lab0: pass
            if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "-":
                return False
            self.cursor -= 1
            break
        self.slice_del()
        return True

    def __r_final_apostrophe(self):
        self.ket = self.cursor
        if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "'":
            return False
        self.cursor -= 1
        self.bra = self.cursor
        while True:
            v_1 = self.limit - self.cursor
            try:
                if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "l":
                    raise lab0()
                self.cursor -= 1
                if self.cursor > self.limit_backward:
                    raise lab0()
                self.slice_from("a")
                break
            except lab0: pass
            self.cursor = self.limit - v_1
            try:
                if not self.eq_s_b("un"):
                    raise lab0()
                if self.cursor > self.limit_backward:
                    raise lab0()
                self.slice_from("u")
                break
            except lab0: pass
            self.cursor = self.limit - v_1
            try:
                if self.find_among_b(EsperantoStemmer.a_3) == 0:
                    raise lab0()
                while True:
                    try:
                        if self.cursor > self.limit_backward:
                            raise lab1()
                        break
                    except lab1: pass
                    if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "-":
                        raise lab0()
                    self.cursor -= 1
                    break
                self.slice_from("aŭ")
                break
            except lab0: pass
            self.cursor = self.limit - v_1
            self.slice_from("o")
            break
        return True

    def __r_ujn_suffix(self):
        self.ket = self.cursor
        v_1 = self.limit - self.cursor
        try:
            if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "n":
                self.cursor = self.limit - v_1
                raise lab0()
            self.cursor -= 1
        except lab0: pass
        v_2 = self.limit - self.cursor
        try:
            if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "j":
                self.cursor = self.limit - v_2
                raise lab0()
            self.cursor -= 1
        except lab0: pass
        self.bra = self.cursor
        if self.find_among_b(EsperantoStemmer.a_4) == 0:
            return False
        while True:
            try:
                if self.cursor > self.limit_backward:
                    raise lab0()
                break
            except lab0: pass
            if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "-":
                return False
            self.cursor -= 1
            break
        self.slice_del()
        return True

    def __r_uninflected(self):
        if self.find_among_b(EsperantoStemmer.a_5) == 0:
            return False
        while True:
            try:
                if self.cursor > self.limit_backward:
                    raise lab0()
                break
            except lab0: pass
            if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "-":
                return False
            self.cursor -= 1
            break
        return True

    def __r_merged_numeral(self):
        if self.find_among_b(EsperantoStemmer.a_6) == 0:
            return False
        return self.find_among_b(EsperantoStemmer.a_7) != 0

    def __r_correlative(self):
        self.ket = self.cursor
        self.bra = self.cursor
        v_1 = self.limit - self.cursor
        while True:
            v_2 = self.limit - self.cursor
            try:
                v_3 = self.limit - self.cursor
                try:
                    if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "n":
                        self.cursor = self.limit - v_3
                        raise lab1()
                    self.cursor -= 1
                except lab1: pass
                self.bra = self.cursor
                if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "e":
                    raise lab0()
                self.cursor -= 1
                break
            except lab0: pass
            self.cursor = self.limit - v_2
            v_4 = self.limit - self.cursor
            try:
                if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "n":
                    self.cursor = self.limit - v_4
                    raise lab0()
                self.cursor -= 1
            except lab0: pass
            v_5 = self.limit - self.cursor
            try:
                if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "j":
                    self.cursor = self.limit - v_5
                    raise lab0()
                self.cursor -= 1
            except lab0: pass
            self.bra = self.cursor
            if not self.in_grouping_b(EsperantoStemmer.g_aou):
                return False
            break
        if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "i":
            return False
        self.cursor -= 1
        v_6 = self.limit - self.cursor
        try:
            if self.find_among_b(EsperantoStemmer.a_8) == 0:
                self.cursor = self.limit - v_6
                raise lab0()
        except lab0: pass
        while True:
            try:
                if self.cursor > self.limit_backward:
                    raise lab0()
                break
            except lab0: pass
            if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "-":
                return False
            self.cursor -= 1
            break
        self.cursor = self.limit - v_1
        self.slice_del()
        return True

    def __r_long_word(self):
        while True:
            v_1 = self.limit - self.cursor
            try:
                for _ in 0, 0:
                    if not self.go_out_grouping_b(EsperantoStemmer.g_vowel):
                        raise lab0()
                    self.cursor -= 1
                break
            except lab0: pass
            self.cursor = self.limit - v_1
            try:
                while True:
                    try:
                        if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "-":
                            raise lab1()
                        self.cursor -= 1
                        break
                    except lab1: pass
                    if self.cursor <= self.limit_backward:
                        raise lab0()
                    self.cursor -= 1
                if self.cursor <= self.limit_backward:
                    raise lab0()
                self.cursor -= 1
                break
            except lab0: pass
            self.cursor = self.limit - v_1
            if not self.go_out_grouping_b(EsperantoStemmer.g_digit):
                return False
            self.cursor -= 1
            break
        return True

    def __r_standard_suffix(self):
        self.ket = self.cursor
        among_var = self.find_among_b(EsperantoStemmer.a_9)
        if among_var == 0:
            return False
        if among_var == 1:
            v_1 = self.limit - self.cursor
            while True:
                try:
                    if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "-":
                        raise lab0()
                    self.cursor -= 1
                    break
                except lab0: pass
                if not self.in_grouping_b(EsperantoStemmer.g_digit):
                    return False
                break
            self.cursor = self.limit - v_1
        v_2 = self.limit - self.cursor
        try:
            if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "-":
                self.cursor = self.limit - v_2
                raise lab0()
            self.cursor -= 1
        except lab0: pass
        self.bra = self.cursor
        self.slice_del()
        return True

    def _stem(self):
        v_1 = self.cursor
        if not self.__r_canonical_form():
            return False
        self.cursor = v_1
        v_2 = self.cursor
        self.__r_initial_apostrophe()
        self.cursor = v_2
        self.limit_backward = self.cursor
        self.cursor = self.limit
        v_3 = self.limit - self.cursor
        try:
            if not self.__r_pronoun():
                raise lab0()
            return False
        except lab0: pass
        self.cursor = self.limit - v_3
        v_4 = self.limit - self.cursor
        self.__r_final_apostrophe()
        self.cursor = self.limit - v_4
        v_5 = self.limit - self.cursor
        try:
            if not self.__r_correlative():
                raise lab0()
            return False
        except lab0: pass
        self.cursor = self.limit - v_5
        v_6 = self.limit - self.cursor
        try:
            if not self.__r_uninflected():
                raise lab0()
            return False
        except lab0: pass
        self.cursor = self.limit - v_6
        v_7 = self.limit - self.cursor
        try:
            if not self.__r_merged_numeral():
                raise lab0()
            return False
        except lab0: pass
        self.cursor = self.limit - v_7
        v_8 = self.limit - self.cursor
        try:
            if not self.__r_ujn_suffix():
                raise lab0()
            return False
        except lab0: pass
        self.cursor = self.limit - v_8
        v_9 = self.limit - self.cursor
        if not self.__r_long_word():
            return False
        self.cursor = self.limit - v_9
        if not self.__r_standard_suffix():
            return False
        self.cursor = self.limit_backward
        return True

    a_0 = [
        Among("", -1, 14),
        Among("-", 0, 13),
        Among("cx", 0, 1),
        Among("gx", 0, 2),
        Among("hx", 0, 3),
        Among("jx", 0, 4),
        Among("q", 0, 12),
        Among("sx", 0, 5),
        Among("ux", 0, 6),
        Among("w", 0, 12),
        Among("x", 0, 12),
        Among("y", 0, 12),
        Among("á", 0, 7),
        Among("é", 0, 8),
        Among("í", 0, 9),
        Among("ó", 0, 10),
        Among("ú", 0, 11)
    ]

    a_1 = [
        Among("as", -1, -1),
        Among("i", -1, -1),
        Among("is", 1, -1),
        Among("os", -1, -1),
        Among("u", -1, -1),
        Among("us", 4, -1)
    ]

    a_2 = [
        Among("ci", -1, -1),
        Among("gi", -1, -1),
        Among("hi", -1, -1),
        Among("li", -1, -1),
        Among("ili", 3, -1),
        Among("ŝli", 3, -1),
        Among("mi", -1, -1),
        Among("ni", -1, -1),
        Among("oni", 7, -1),
        Among("ri", -1, -1),
        Among("si", -1, -1),
        Among("vi", -1, -1),
        Among("ivi", 11, -1),
        Among("ĝi", -1, -1),
        Among("ŝi", -1, -1),
        Among("iŝi", 14, -1),
        Among("malŝi", 14, -1)
    ]

    a_3 = [
        Among("amb", -1, -1),
        Among("bald", -1, -1),
        Among("malbald", 1, -1),
        Among("morg", -1, -1),
        Among("postmorg", 3, -1),
        Among("adi", -1, -1),
        Among("hodi", -1, -1),
        Among("ank", -1, -1),
        Among("ĉirk", -1, -1),
        Among("tutĉirk", 8, -1),
        Among("presk", -1, -1),
        Among("almen", -1, -1),
        Among("apen", -1, -1),
        Among("hier", -1, -1),
        Among("antaŭhier", 13, -1),
        Among("malgr", -1, -1),
        Among("ankor", -1, -1),
        Among("kontr", -1, -1),
        Among("anstat", -1, -1),
        Among("kvaz", -1, -1)
    ]

    a_4 = [
        Among("aliu", -1, -1),
        Among("unu", -1, -1)
    ]

    a_5 = [
        Among("aha", -1, -1),
        Among("haha", 0, -1),
        Among("haleluja", -1, -1),
        Among("hola", -1, -1),
        Among("hosana", -1, -1),
        Among("maltra", -1, -1),
        Among("hura", -1, -1),
        Among("ĥaĥa", -1, -1),
        Among("ekde", -1, -1),
        Among("elde", -1, -1),
        Among("disde", -1, -1),
        Among("ehe", -1, -1),
        Among("maltre", -1, -1),
        Among("dirlididi", -1, -1),
        Among("malpli", -1, -1),
        Among("malĉi", -1, -1),
        Among("malkaj", -1, -1),
        Among("amen", -1, -1),
        Among("tamen", 17, -1),
        Among("oho", -1, -1),
        Among("maltro", -1, -1),
        Among("minus", -1, -1),
        Among("uhu", -1, -1),
        Among("muu", -1, -1)
    ]

    a_6 = [
        Among("tri", -1, -1),
        Among("du", -1, -1),
        Among("unu", -1, -1)
    ]

    a_7 = [
        Among("dek", -1, -1),
        Among("cent", -1, -1)
    ]

    a_8 = [
        Among("k", -1, -1),
        Among("kelk", 0, -1),
        Among("nen", -1, -1),
        Among("t", -1, -1),
        Among("mult", 3, -1),
        Among("samt", 3, -1),
        Among("ĉ", -1, -1)
    ]

    a_9 = [
        Among("a", -1, -1),
        Among("e", -1, -1),
        Among("i", -1, -1),
        Among("j", -1, 1),
        Among("aj", 3, -1),
        Among("oj", 3, -1),
        Among("n", -1, 1),
        Among("an", 6, -1),
        Among("en", 6, -1),
        Among("jn", 6, 1),
        Among("ajn", 9, -1),
        Among("ojn", 9, -1),
        Among("on", 6, -1),
        Among("o", -1, -1),
        Among("as", -1, -1),
        Among("is", -1, -1),
        Among("os", -1, -1),
        Among("us", -1, -1),
        Among("u", -1, -1)
    ]


class lab0(BaseException): pass


class lab1(BaseException): pass
