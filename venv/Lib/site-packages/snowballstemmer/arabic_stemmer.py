# Generated from arabic.sbl by Snowball 3.1.1 - https://snowballstem.org/

from .basestemmer import BaseStemmer
from .among import Among


class ArabicStemmer(BaseStemmer):
    '''
    This class implements the stemming algorithm defined by a snowball script.
    Generated from arabic.sbl by Snowball 3.1.1 - https://snowballstem.org/
    '''

    B_is_defined = False
    B_is_verb = False
    B_is_noun = False

    def __r_Normalize_pre(self):
        v_1 = self.cursor
        try:
            while True:
                v_2 = self.cursor
                try:
                    while True:
                        v_3 = self.cursor
                        try:
                            self.bra = self.cursor
                            among_var = self.find_among(ArabicStemmer.a_0)
                            if among_var == 0:
                                raise lab2()
                            self.ket = self.cursor
                            self.slice_from(ArabicStemmer.as_0[among_var - 1])
                            break
                        except lab2: pass
                        self.cursor = v_3
                        if self.cursor >= self.limit:
                            raise lab1()
                        self.cursor += 1
                        break
                    continue
                except lab1: pass
                self.cursor = v_2
                break
        except lab0: pass
        self.cursor = v_1
        return True

    def __r_Normalize_post(self):
        v_1 = self.cursor
        try:
            self.limit_backward = self.cursor
            self.cursor = self.limit
            self.ket = self.cursor
            if self.find_among_b(ArabicStemmer.a_1) == 0:
                raise lab0()
            self.bra = self.cursor
            self.slice_from("\u0621")
            self.cursor = self.limit_backward
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
                            self.bra = self.cursor
                            among_var = self.find_among(ArabicStemmer.a_2)
                            if among_var == 0:
                                raise lab2()
                            self.ket = self.cursor
                            self.slice_from(ArabicStemmer.as_2[among_var - 1])
                            break
                        except lab2: pass
                        self.cursor = v_4
                        if self.cursor >= self.limit:
                            raise lab1()
                        self.cursor += 1
                        break
                    continue
                except lab1: pass
                self.cursor = v_3
                break
        except lab0: pass
        self.cursor = v_2
        return True

    def __r_Checks1(self):
        self.bra = self.cursor
        among_var = self.find_among(ArabicStemmer.a_3)
        if among_var == 0:
            return False
        self.ket = self.cursor
        if among_var == 1:
            if len(self.current) <= 4:
                return False
            self.B_is_noun = True
            self.B_is_verb = False
            self.B_is_defined = True
        else:
            if len(self.current) <= 3:
                return False
            self.B_is_noun = True
            self.B_is_verb = False
            self.B_is_defined = True
        return True

    def __r_Prefix_Step1(self):
        self.bra = self.cursor
        among_var = self.find_among(ArabicStemmer.a_4)
        if among_var == 0:
            return False
        self.ket = self.cursor
        if among_var == 1:
            if len(self.current) <= 3:
                return False
            self.slice_from("\u0623")
        elif among_var == 2:
            if len(self.current) <= 3:
                return False
            self.slice_from("\u0622")
        elif among_var == 3:
            if len(self.current) <= 3:
                return False
            self.slice_from("\u0627")
        else:
            if len(self.current) <= 3:
                return False
            self.slice_from("\u0625")
        return True

    def __r_Prefix_Step2(self):
        self.bra = self.cursor
        if self.find_among(ArabicStemmer.a_5) == 0:
            return False
        self.ket = self.cursor
        if len(self.current) <= 3:
            return False
        try:
            if self.cursor == self.limit or self.current[self.cursor] != "\u0627":
                raise lab0()
            self.cursor += 1
            return False
        except lab0: pass
        self.slice_del()
        return True

    def __r_Prefix_Step3a_Noun(self):
        self.bra = self.cursor
        among_var = self.find_among(ArabicStemmer.a_6)
        if among_var == 0:
            return False
        self.ket = self.cursor
        if among_var == 1:
            if len(self.current) <= 5:
                return False
            self.slice_del()
        else:
            if len(self.current) <= 4:
                return False
            self.slice_del()
        return True

    def __r_Prefix_Step3b_Noun(self):
        self.bra = self.cursor
        among_var = self.find_among(ArabicStemmer.a_7)
        if among_var == 0:
            return False
        self.ket = self.cursor
        if among_var == 1:
            if len(self.current) <= 3:
                return False
            self.slice_del()
        elif among_var == 2:
            if len(self.current) <= 3:
                return False
            self.slice_from("\u0628")
        elif among_var == 3:
            if len(self.current) <= 3:
                return False
            self.slice_from("\u0643")
        return True

    def __r_Prefix_Step3_Verb(self):
        self.bra = self.cursor
        among_var = self.find_among(ArabicStemmer.a_8)
        if among_var == 0:
            return False
        self.ket = self.cursor
        if among_var == 1:
            if len(self.current) <= 4:
                return False
            self.slice_from("\u064A")
        elif among_var == 2:
            if len(self.current) <= 4:
                return False
            self.slice_from("\u062A")
        elif among_var == 3:
            if len(self.current) <= 4:
                return False
            self.slice_from("\u0646")
        else:
            if len(self.current) <= 4:
                return False
            self.slice_from("\u0623")
        return True

    def __r_Prefix_Step4_Verb(self):
        self.bra = self.cursor
        if self.find_among(ArabicStemmer.a_9) == 0:
            return False
        self.ket = self.cursor
        if len(self.current) <= 4:
            return False
        self.B_is_verb = True
        self.B_is_noun = False
        self.slice_from("\u0627\u0633\u062A")
        return True

    def __r_Suffix_Noun_Step1a(self):
        self.ket = self.cursor
        among_var = self.find_among_b(ArabicStemmer.a_10)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var == 1:
            if len(self.current) < 4:
                return False
            self.slice_del()
        elif among_var == 2:
            if len(self.current) < 5:
                return False
            self.slice_del()
        else:
            if len(self.current) < 6:
                return False
            self.slice_del()
        return True

    def __r_Suffix_Noun_Step1b(self):
        self.ket = self.cursor
        if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "\u0646":
            return False
        self.cursor -= 1
        self.bra = self.cursor
        if len(self.current) <= 5:
            return False
        self.slice_del()
        return True

    def __r_Suffix_Noun_Step2a(self):
        self.ket = self.cursor
        if self.find_among_b(ArabicStemmer.a_11) == 0:
            return False
        self.bra = self.cursor
        if len(self.current) <= 4:
            return False
        self.slice_del()
        return True

    def __r_Suffix_Noun_Step2b(self):
        self.ket = self.cursor
        if not self.eq_s_b("\u0627\u062A"):
            return False
        self.bra = self.cursor
        if len(self.current) < 5:
            return False
        self.slice_del()
        return True

    def __r_Suffix_Noun_Step2c1(self):
        self.ket = self.cursor
        if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "\u062A":
            return False
        self.cursor -= 1
        self.bra = self.cursor
        if len(self.current) < 4:
            return False
        self.slice_del()
        return True

    def __r_Suffix_Noun_Step2c2(self):
        self.ket = self.cursor
        if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "\u0629":
            return False
        self.cursor -= 1
        self.bra = self.cursor
        if len(self.current) < 4:
            return False
        self.slice_del()
        return True

    def __r_Suffix_Noun_Step3(self):
        self.ket = self.cursor
        if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "\u064A":
            return False
        self.cursor -= 1
        self.bra = self.cursor
        if len(self.current) < 3:
            return False
        self.slice_del()
        return True

    def __r_Suffix_Verb_Step1(self):
        self.ket = self.cursor
        among_var = self.find_among_b(ArabicStemmer.a_12)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var == 1:
            if len(self.current) < 4:
                return False
            self.slice_del()
        elif among_var == 2:
            if len(self.current) < 5:
                return False
            self.slice_del()
        else:
            if len(self.current) < 6:
                return False
            self.slice_del()
        return True

    def __r_Suffix_Verb_Step2a(self):
        self.ket = self.cursor
        among_var = self.find_among_b(ArabicStemmer.a_13)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var == 1:
            if len(self.current) < 4:
                return False
            self.slice_del()
        elif among_var == 2:
            if len(self.current) < 5:
                return False
            self.slice_del()
        elif among_var == 3:
            if len(self.current) <= 5:
                return False
            self.slice_del()
        else:
            if len(self.current) < 6:
                return False
            self.slice_del()
        return True

    def __r_Suffix_Verb_Step2b(self):
        self.ket = self.cursor
        if self.find_among_b(ArabicStemmer.a_14) == 0:
            return False
        self.bra = self.cursor
        if len(self.current) < 5:
            return False
        self.slice_del()
        return True

    def __r_Suffix_Verb_Step2c(self):
        self.ket = self.cursor
        among_var = self.find_among_b(ArabicStemmer.a_15)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var == 1:
            if len(self.current) < 4:
                return False
            self.slice_del()
        else:
            if len(self.current) < 6:
                return False
            self.slice_del()
        return True

    def __r_Suffix_All_alef_maqsura(self):
        self.ket = self.cursor
        if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "\u0649":
            return False
        self.cursor -= 1
        self.bra = self.cursor
        self.slice_from("\u064A")
        return True

    def _stem(self):
        self.B_is_noun = True
        self.B_is_verb = True
        self.B_is_defined = False
        v_1 = self.cursor
        self.__r_Checks1()
        self.cursor = v_1
        self.__r_Normalize_pre()
        self.limit_backward = self.cursor
        self.cursor = self.limit
        v_2 = self.limit - self.cursor
        try:
            while True:
                v_3 = self.limit - self.cursor
                try:
                    if not self.B_is_verb:
                        raise lab1()
                    while True:
                        v_4 = self.limit - self.cursor
                        try:
                            v_5 = 1
                            while True:
                                v_6 = self.limit - self.cursor
                                try:
                                    if not self.__r_Suffix_Verb_Step1():
                                        raise lab3()
                                    v_5 -= 1
                                    continue
                                except lab3: pass
                                self.cursor = self.limit - v_6
                                break
                            if v_5 > 0:
                                raise lab2()
                            while True:
                                v_7 = self.limit - self.cursor
                                try:
                                    if not self.__r_Suffix_Verb_Step2a():
                                        raise lab3()
                                    break
                                except lab3: pass
                                self.cursor = self.limit - v_7
                                try:
                                    if not self.__r_Suffix_Verb_Step2c():
                                        raise lab3()
                                    break
                                except lab3: pass
                                self.cursor = self.limit - v_7
                                if self.cursor <= self.limit_backward:
                                    raise lab2()
                                self.cursor -= 1
                                break
                            break
                        except lab2: pass
                        self.cursor = self.limit - v_4
                        try:
                            if not self.__r_Suffix_Verb_Step2b():
                                raise lab2()
                            break
                        except lab2: pass
                        self.cursor = self.limit - v_4
                        if not self.__r_Suffix_Verb_Step2a():
                            raise lab1()
                        break
                    break
                except lab1: pass
                self.cursor = self.limit - v_3
                try:
                    if not self.B_is_noun:
                        raise lab1()
                    v_8 = self.limit - self.cursor
                    try:
                        while True:
                            v_9 = self.limit - self.cursor
                            try:
                                if not self.__r_Suffix_Noun_Step2c2():
                                    raise lab3()
                                break
                            except lab3: pass
                            self.cursor = self.limit - v_9
                            try:
                                if self.B_is_defined:
                                    raise lab3()
                                if not self.__r_Suffix_Noun_Step1a():
                                    raise lab3()
                                while True:
                                    v_10 = self.limit - self.cursor
                                    try:
                                        if not self.__r_Suffix_Noun_Step2a():
                                            raise lab4()
                                        break
                                    except lab4: pass
                                    self.cursor = self.limit - v_10
                                    try:
                                        if not self.__r_Suffix_Noun_Step2b():
                                            raise lab4()
                                        break
                                    except lab4: pass
                                    self.cursor = self.limit - v_10
                                    try:
                                        if not self.__r_Suffix_Noun_Step2c1():
                                            raise lab4()
                                        break
                                    except lab4: pass
                                    self.cursor = self.limit - v_10
                                    if self.cursor <= self.limit_backward:
                                        raise lab3()
                                    self.cursor -= 1
                                    break
                                break
                            except lab3: pass
                            self.cursor = self.limit - v_9
                            try:
                                if not self.__r_Suffix_Noun_Step1b():
                                    raise lab3()
                                while True:
                                    v_11 = self.limit - self.cursor
                                    try:
                                        if not self.__r_Suffix_Noun_Step2a():
                                            raise lab4()
                                        break
                                    except lab4: pass
                                    self.cursor = self.limit - v_11
                                    try:
                                        if not self.__r_Suffix_Noun_Step2b():
                                            raise lab4()
                                        break
                                    except lab4: pass
                                    self.cursor = self.limit - v_11
                                    if not self.__r_Suffix_Noun_Step2c1():
                                        raise lab3()
                                    break
                                break
                            except lab3: pass
                            self.cursor = self.limit - v_9
                            try:
                                if self.B_is_defined:
                                    raise lab3()
                                if not self.__r_Suffix_Noun_Step2a():
                                    raise lab3()
                                break
                            except lab3: pass
                            self.cursor = self.limit - v_9
                            if not self.__r_Suffix_Noun_Step2b():
                                self.cursor = self.limit - v_8
                                raise lab2()
                            break
                    except lab2: pass
                    if not self.__r_Suffix_Noun_Step3():
                        raise lab1()
                    break
                except lab1: pass
                self.cursor = self.limit - v_3
                if not self.__r_Suffix_All_alef_maqsura():
                    raise lab0()
                break
        except lab0: pass
        self.cursor = self.limit - v_2
        self.cursor = self.limit_backward
        v_12 = self.cursor
        try:
            v_13 = self.cursor
            try:
                if not self.__r_Prefix_Step1():
                    self.cursor = v_13
                    raise lab1()
            except lab1: pass
            v_14 = self.cursor
            try:
                if not self.__r_Prefix_Step2():
                    self.cursor = v_14
                    raise lab1()
            except lab1: pass
            while True:
                v_15 = self.cursor
                try:
                    if not self.__r_Prefix_Step3a_Noun():
                        raise lab1()
                    break
                except lab1: pass
                self.cursor = v_15
                try:
                    if not self.B_is_noun:
                        raise lab1()
                    if not self.__r_Prefix_Step3b_Noun():
                        raise lab1()
                    break
                except lab1: pass
                self.cursor = v_15
                if not self.B_is_verb:
                    raise lab0()
                v_16 = self.cursor
                try:
                    if not self.__r_Prefix_Step3_Verb():
                        self.cursor = v_16
                        raise lab1()
                except lab1: pass
                if not self.__r_Prefix_Step4_Verb():
                    raise lab0()
                break
        except lab0: pass
        self.cursor = v_12
        self.__r_Normalize_post()
        return True

    a_0 = [
        Among("\u0640", -1, 1),
        Among("\u064B", -1, 1),
        Among("\u064C", -1, 1),
        Among("\u064D", -1, 1),
        Among("\u064E", -1, 1),
        Among("\u064F", -1, 1),
        Among("\u0650", -1, 1),
        Among("\u0651", -1, 1),
        Among("\u0652", -1, 1),
        Among("\u0660", -1, 2),
        Among("\u0661", -1, 3),
        Among("\u0662", -1, 4),
        Among("\u0663", -1, 5),
        Among("\u0664", -1, 6),
        Among("\u0665", -1, 7),
        Among("\u0666", -1, 8),
        Among("\u0667", -1, 9),
        Among("\u0668", -1, 10),
        Among("\u0669", -1, 11),
        Among("\uFE80", -1, 12),
        Among("\uFE81", -1, 16),
        Among("\uFE82", -1, 16),
        Among("\uFE83", -1, 13),
        Among("\uFE84", -1, 13),
        Among("\uFE85", -1, 17),
        Among("\uFE86", -1, 17),
        Among("\uFE87", -1, 14),
        Among("\uFE88", -1, 14),
        Among("\uFE89", -1, 15),
        Among("\uFE8A", -1, 15),
        Among("\uFE8B", -1, 15),
        Among("\uFE8C", -1, 15),
        Among("\uFE8D", -1, 18),
        Among("\uFE8E", -1, 18),
        Among("\uFE8F", -1, 19),
        Among("\uFE90", -1, 19),
        Among("\uFE91", -1, 19),
        Among("\uFE92", -1, 19),
        Among("\uFE93", -1, 20),
        Among("\uFE94", -1, 20),
        Among("\uFE95", -1, 21),
        Among("\uFE96", -1, 21),
        Among("\uFE97", -1, 21),
        Among("\uFE98", -1, 21),
        Among("\uFE99", -1, 22),
        Among("\uFE9A", -1, 22),
        Among("\uFE9B", -1, 22),
        Among("\uFE9C", -1, 22),
        Among("\uFE9D", -1, 23),
        Among("\uFE9E", -1, 23),
        Among("\uFE9F", -1, 23),
        Among("\uFEA0", -1, 23),
        Among("\uFEA1", -1, 24),
        Among("\uFEA2", -1, 24),
        Among("\uFEA3", -1, 24),
        Among("\uFEA4", -1, 24),
        Among("\uFEA5", -1, 25),
        Among("\uFEA6", -1, 25),
        Among("\uFEA7", -1, 25),
        Among("\uFEA8", -1, 25),
        Among("\uFEA9", -1, 26),
        Among("\uFEAA", -1, 26),
        Among("\uFEAB", -1, 27),
        Among("\uFEAC", -1, 27),
        Among("\uFEAD", -1, 28),
        Among("\uFEAE", -1, 28),
        Among("\uFEAF", -1, 29),
        Among("\uFEB0", -1, 29),
        Among("\uFEB1", -1, 30),
        Among("\uFEB2", -1, 30),
        Among("\uFEB3", -1, 30),
        Among("\uFEB4", -1, 30),
        Among("\uFEB5", -1, 31),
        Among("\uFEB6", -1, 31),
        Among("\uFEB7", -1, 31),
        Among("\uFEB8", -1, 31),
        Among("\uFEB9", -1, 32),
        Among("\uFEBA", -1, 32),
        Among("\uFEBB", -1, 32),
        Among("\uFEBC", -1, 32),
        Among("\uFEBD", -1, 33),
        Among("\uFEBE", -1, 33),
        Among("\uFEBF", -1, 33),
        Among("\uFEC0", -1, 33),
        Among("\uFEC1", -1, 34),
        Among("\uFEC2", -1, 34),
        Among("\uFEC3", -1, 34),
        Among("\uFEC4", -1, 34),
        Among("\uFEC5", -1, 35),
        Among("\uFEC6", -1, 35),
        Among("\uFEC7", -1, 35),
        Among("\uFEC8", -1, 35),
        Among("\uFEC9", -1, 36),
        Among("\uFECA", -1, 36),
        Among("\uFECB", -1, 36),
        Among("\uFECC", -1, 36),
        Among("\uFECD", -1, 37),
        Among("\uFECE", -1, 37),
        Among("\uFECF", -1, 37),
        Among("\uFED0", -1, 37),
        Among("\uFED1", -1, 38),
        Among("\uFED2", -1, 38),
        Among("\uFED3", -1, 38),
        Among("\uFED4", -1, 38),
        Among("\uFED5", -1, 39),
        Among("\uFED6", -1, 39),
        Among("\uFED7", -1, 39),
        Among("\uFED8", -1, 39),
        Among("\uFED9", -1, 40),
        Among("\uFEDA", -1, 40),
        Among("\uFEDB", -1, 40),
        Among("\uFEDC", -1, 40),
        Among("\uFEDD", -1, 41),
        Among("\uFEDE", -1, 41),
        Among("\uFEDF", -1, 41),
        Among("\uFEE0", -1, 41),
        Among("\uFEE1", -1, 42),
        Among("\uFEE2", -1, 42),
        Among("\uFEE3", -1, 42),
        Among("\uFEE4", -1, 42),
        Among("\uFEE5", -1, 43),
        Among("\uFEE6", -1, 43),
        Among("\uFEE7", -1, 43),
        Among("\uFEE8", -1, 43),
        Among("\uFEE9", -1, 44),
        Among("\uFEEA", -1, 44),
        Among("\uFEEB", -1, 44),
        Among("\uFEEC", -1, 44),
        Among("\uFEED", -1, 45),
        Among("\uFEEE", -1, 45),
        Among("\uFEEF", -1, 46),
        Among("\uFEF0", -1, 46),
        Among("\uFEF1", -1, 47),
        Among("\uFEF2", -1, 47),
        Among("\uFEF3", -1, 47),
        Among("\uFEF4", -1, 47),
        Among("\uFEF5", -1, 51),
        Among("\uFEF6", -1, 51),
        Among("\uFEF7", -1, 49),
        Among("\uFEF8", -1, 49),
        Among("\uFEF9", -1, 50),
        Among("\uFEFA", -1, 50),
        Among("\uFEFB", -1, 48),
        Among("\uFEFC", -1, 48)
    ]
    as_0 = ("", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "\u0621", "\u0623", "\u0625", "\u0626", "\u0622", "\u0624", "\u0627", "\u0628", "\u0629", "\u062A", "\u062B", "\u062C", "\u062D", "\u062E", "\u062F", "\u0630", "\u0631", "\u0632", "\u0633", "\u0634", "\u0635", "\u0636", "\u0637", "\u0638", "\u0639", "\u063A", "\u0641", "\u0642", "\u0643", "\u0644", "\u0645", "\u0646", "\u0647", "\u0648", "\u0649", "\u064A", "\u0644\u0627", "\u0644\u0623", "\u0644\u0625", "\u0644\u0622")

    a_1 = [
        Among("\u0622", -1, 1),
        Among("\u0623", -1, 1),
        Among("\u0624", -1, 1),
        Among("\u0625", -1, 1),
        Among("\u0626", -1, 1)
    ]

    a_2 = [
        Among("\u0622", -1, 1),
        Among("\u0623", -1, 1),
        Among("\u0624", -1, 2),
        Among("\u0625", -1, 1),
        Among("\u0626", -1, 3)
    ]
    as_2 = ("\u0627", "\u0648", "\u064A")

    a_3 = [
        Among("\u0627\u0644", -1, 2),
        Among("\u0628\u0627\u0644", -1, 1),
        Among("\u0643\u0627\u0644", -1, 1),
        Among("\u0644\u0644", -1, 2)
    ]

    a_4 = [
        Among("\u0623\u0622", -1, 2),
        Among("\u0623\u0623", -1, 1),
        Among("\u0623\u0624", -1, 1),
        Among("\u0623\u0625", -1, 4),
        Among("\u0623\u0627", -1, 3)
    ]

    a_5 = [
        Among("\u0641", -1, 1),
        Among("\u0648", -1, 1)
    ]

    a_6 = [
        Among("\u0627\u0644", -1, 2),
        Among("\u0628\u0627\u0644", -1, 1),
        Among("\u0643\u0627\u0644", -1, 1),
        Among("\u0644\u0644", -1, 2)
    ]

    a_7 = [
        Among("\u0628", -1, 1),
        Among("\u0628\u0627", 0, -1),
        Among("\u0628\u0628", 0, 2),
        Among("\u0643\u0643", -1, 3)
    ]

    a_8 = [
        Among("\u0633\u0623", -1, 4),
        Among("\u0633\u062A", -1, 2),
        Among("\u0633\u0646", -1, 3),
        Among("\u0633\u064A", -1, 1)
    ]

    a_9 = [
        Among("\u062A\u0633\u062A", -1, 1),
        Among("\u0646\u0633\u062A", -1, 1),
        Among("\u064A\u0633\u062A", -1, 1)
    ]

    a_10 = [
        Among("\u0643\u0645\u0627", -1, 3),
        Among("\u0647\u0645\u0627", -1, 3),
        Among("\u0646\u0627", -1, 2),
        Among("\u0647\u0627", -1, 2),
        Among("\u0643", -1, 1),
        Among("\u0643\u0645", -1, 2),
        Among("\u0647\u0645", -1, 2),
        Among("\u0647\u0646", -1, 2),
        Among("\u0647", -1, 1),
        Among("\u064A", -1, 1)
    ]

    a_11 = [
        Among("\u0627", -1, 1),
        Among("\u0648", -1, 1),
        Among("\u064A", -1, 1)
    ]

    a_12 = [
        Among("\u0643\u0645\u0627", -1, 3),
        Among("\u0647\u0645\u0627", -1, 3),
        Among("\u0646\u0627", -1, 2),
        Among("\u0647\u0627", -1, 2),
        Among("\u0643", -1, 1),
        Among("\u0643\u0645", -1, 2),
        Among("\u0647\u0645", -1, 2),
        Among("\u0643\u0646", -1, 2),
        Among("\u0647\u0646", -1, 2),
        Among("\u0647", -1, 1),
        Among("\u0643\u0645\u0648", -1, 3),
        Among("\u0646\u064A", -1, 2)
    ]

    a_13 = [
        Among("\u0627", -1, 1),
        Among("\u062A\u0627", 0, 2),
        Among("\u062A\u0645\u0627", 0, 4),
        Among("\u0646\u0627", 0, 2),
        Among("\u062A", -1, 1),
        Among("\u0646", -1, 1),
        Among("\u0627\u0646", 5, 3),
        Among("\u062A\u0646", 5, 2),
        Among("\u0648\u0646", 5, 3),
        Among("\u064A\u0646", 5, 3),
        Among("\u064A", -1, 1)
    ]

    a_14 = [
        Among("\u0648\u0627", -1, 1),
        Among("\u062A\u0645", -1, 1)
    ]

    a_15 = [
        Among("\u0648", -1, 1),
        Among("\u062A\u0645\u0648", 0, 2)
    ]


class lab0(BaseException): pass


class lab1(BaseException): pass


class lab2(BaseException): pass


class lab3(BaseException): pass


class lab4(BaseException): pass
