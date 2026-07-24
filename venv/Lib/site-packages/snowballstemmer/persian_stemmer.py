# Generated from persian.sbl by Snowball 3.1.1 - https://snowballstem.org/

from .basestemmer import BaseStemmer
from .among import Among


class PersianStemmer(BaseStemmer):
    '''
    This class implements the stemming algorithm defined by a snowball script.
    Generated from persian.sbl by Snowball 3.1.1 - https://snowballstem.org/
    '''

    I_p1 = 0
    B_remove_verb_person_endings = False
    B_saw_present_prefix = False

    def __r_Normalize_Characters(self):
        while True:
            v_1 = self.cursor
            try:
                self.bra = self.cursor
                among_var = self.find_among(PersianStemmer.a_0)
                self.ket = self.cursor
                if among_var == 1:
                    self.slice_from("\u06A9")
                elif among_var == 2:
                    self.slice_from("\u06CC")
                elif among_var == 3:
                    self.slice_from("\u0647")
                elif among_var == 4:
                    self.slice_from("\u0627")
                elif among_var == 5:
                    self.slice_from("\u0648")
                elif among_var == 6:
                    self.slice_del()
                else:
                    if self.cursor >= self.limit:
                        raise lab0()
                    self.cursor += 1
                continue
            except lab0: pass
            self.cursor = v_1
            break
        return True

    def __r_Prefixes(self):
        self.bra = self.cursor
        among_var = self.find_among(PersianStemmer.a_1)
        if among_var == 0:
            return False
        self.ket = self.cursor
        if among_var == 1:
            if self.cursor + 2 > self.limit:
                return False
            self.cursor += 2
            self.B_saw_present_prefix = True
        else:
            if self.cursor + 2 > self.limit:
                return False
            self.cursor += 2
            self.slice_del()
            self.B_saw_present_prefix = True
        return True

    def __r_Delete_ZWNJ(self):
        while True:
            v_1 = self.cursor
            try:
                while True:
                    v_2 = self.cursor
                    try:
                        self.bra = self.cursor
                        if self.cursor == self.limit or self.current[self.cursor] != "\u200C":
                            raise lab1()
                        self.cursor += 1
                        self.ket = self.cursor
                        self.slice_del()
                        self.cursor = v_2
                        break
                    except lab1: pass
                    self.cursor = v_2
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

    def __r_Protect_Lexical_AN(self):
        v_1 = self.limit - self.cursor
        try:
            if not self.__r_AN_Exception():
                raise lab0()
            return False
        except lab0: pass
        self.cursor = self.limit - v_1
        v_2 = self.limit - self.cursor
        try:
            if self.find_among_b(PersianStemmer.a_2) == 0:
                raise lab0()
            return False
        except lab0: pass
        self.cursor = self.limit - v_2
        return True

    def __r_AN_Exception(self):
        if self.find_among_b(PersianStemmer.a_3) == 0:
            return False
        if self.cursor > self.limit_backward:
            return False
        return True

    def __r_Irregular_Noun(self):
        self.ket = self.cursor
        among_var = self.find_among_b(PersianStemmer.a_4)
        if among_var == 0:
            return False
        self.bra = self.cursor
        self.slice_from(PersianStemmer.as_4[among_var - 1])
        return True

    def __r_Stem_Noun_or_Adjective(self):
        while True:
            v_1 = self.limit - self.cursor
            try:
                if not self.__r_Irregular_Noun():
                    raise lab0()
                break
            except lab0: pass
            self.cursor = self.limit - v_1
            if self.cursor < self.I_p1:
                return False
            v_3 = self.limit_backward
            self.limit_backward = self.I_p1
            self.ket = self.cursor
            among_var = self.find_among_b(PersianStemmer.a_5)
            if among_var == 0:
                self.limit_backward = v_3
                return False
            self.bra = self.cursor
            if among_var == 1:
                self.slice_del()
            else:
                if self.cursor <= self.limit_backward:
                    self.limit_backward = v_3
                    return False
                self.slice_del()
            self.limit_backward = v_3
            break
        return True

    def __r_Stem_Verb(self):
        while True:
            v_1 = self.limit - self.cursor
            try:
                self.ket = self.cursor
                if self.find_among_b(PersianStemmer.a_6) == 0:
                    raise lab0()
                self.bra = self.cursor
                if not self.__r_R1():
                    raise lab0()
                self.slice_del()
                break
            except lab0: pass
            self.cursor = self.limit - v_1
            self.ket = self.cursor
            among_var = self.find_among_b(PersianStemmer.a_7)
            if among_var == 0:
                return False
            self.bra = self.cursor
            if among_var == 1:
                if not self.B_remove_verb_person_endings:
                    return False
                if not self.__r_R1():
                    return False
                self.slice_del()
            elif among_var == 2:
                self.slice_from("\u0631\u0641\u062A")
            elif among_var == 3:
                if not self.__r_R1():
                    return False
                self.slice_del()
                self.B_remove_verb_person_endings = True
            elif among_var == 4:
                if self.cursor <= self.limit_backward:
                    return False
                self.slice_from("\u062F")
                self.B_remove_verb_person_endings = True
            else:
                if self.cursor <= self.limit_backward:
                    return False
                self.slice_from("\u062A")
                self.B_remove_verb_person_endings = True
            break
        return True

    def _stem(self):
        self.B_saw_present_prefix = False
        v_1 = self.cursor
        self.__r_Normalize_Characters()
        self.cursor = v_1
        v_2 = self.cursor
        self.__r_Prefixes()
        self.cursor = v_2
        v_3 = self.cursor
        self.__r_Delete_ZWNJ()
        self.cursor = v_3
        self.I_p1 = self.limit
        v_4 = self.cursor
        try:
            if self.cursor + 3 > self.limit:
                raise lab0()
            self.cursor += 3
            self.I_p1 = self.cursor
        except lab0: pass
        self.cursor = v_4
        self.limit_backward = self.cursor
        self.cursor = self.limit
        while True:
            v_5 = self.limit - self.cursor
            try:
                v_6 = self.limit - self.cursor
                self.B_remove_verb_person_endings = False
                try:
                    if not self.B_saw_present_prefix:
                        raise lab1()
                    self.B_remove_verb_person_endings = True
                except lab1: pass
                if not self.__r_Protect_Lexical_AN():
                    raise lab0()
                while True:
                    v_7 = self.limit - self.cursor
                    try:
                        if not self.__r_Stem_Noun_or_Adjective():
                            raise lab1()
                        break
                    except lab1: pass
                    self.cursor = self.limit - v_7
                    if not self.__r_Stem_Verb():
                        raise lab0()
                    break
                self.cursor = self.limit - v_6
                continue
            except lab0: pass
            self.cursor = self.limit - v_5
            break
        self.cursor = self.limit_backward
        return True

    a_0 = [
        Among("", -1, 7),
        Among(" ", 0, 6),
        Among("\u0623", 0, 4),
        Among("\u0624", 0, 5),
        Among("\u0625", 0, 4),
        Among("\u0626", 0, 2),
        Among("\u0629", 0, 3),
        Among("\u0643", 0, 1),
        Among("\u064A", 0, 2),
        Among("\u06C1", 0, 3),
        Among("\u200D", 0, 6)
    ]

    a_1 = [
        Among("\u0645\u06CC\u200C", -1, 2),
        Among("\u0646\u0645\u06CC\u200C", -1, 1)
    ]

    a_2 = [
        Among("\u0633\u062A\u0627\u0646", -1, -1),
        Among("\u0631\u0627\u0646", -1, -1),
        Among("\u0633\u0627\u0646", -1, -1),
        Among("\u0648\u0627\u0646", -1, -1)
    ]

    a_3 = [
        Among("\u0622\u0630\u0631\u0628\u0627\u06CC\u062C\u0627\u0646", -1, 1),
        Among("\u0647\u0645\u062F\u0627\u0646", -1, 1),
        Among("\u062E\u0627\u0646\u062F\u0627\u0646", -1, 1),
        Among("\u0632\u0646\u062F\u0627\u0646", -1, 1),
        Among("\u0645\u06CC\u0632\u0627\u0646", -1, 1),
        Among("\u062F\u0631\u062E\u0634\u0627\u0646", -1, 1),
        Among("\u0622\u062A\u0634\u0641\u0634\u0627\u0646", -1, 1),
        Among("\u0646\u0634\u0627\u0646", -1, 1),
        Among("\u06A9\u0647\u06A9\u0634\u0627\u0646", -1, 1),
        Among("\u0627\u06CC\u0634\u0627\u0646", -1, 1),
        Among("\u067E\u0631\u06CC\u0634\u0627\u0646", -1, 1),
        Among("\u0633\u0644\u0637\u0627\u0646", -1, 1),
        Among("\u06AF\u06CC\u0644\u0627\u0646", -1, 1),
        Among("\u0633\u0627\u062E\u062A\u0645\u0627\u0646", -1, 1),
        Among("\u0631\u0645\u0627\u0646", -1, 1),
        Among("\u062F\u0631\u0645\u0627\u0646", 14, 1),
        Among("\u0642\u0647\u0631\u0645\u0627\u0646", 14, 1),
        Among("\u06A9\u0631\u0645\u0627\u0646", 14, 1),
        Among("\u0633\u0627\u0632\u0645\u0627\u0646", -1, 1),
        Among("\u0647\u0645\u0632\u0645\u0627\u0646", -1, 1),
        Among("\u0622\u0633\u0645\u0627\u0646", -1, 1),
        Among("\u0622\u0644\u0645\u0627\u0646", -1, 1),
        Among("\u0645\u0633\u0644\u0645\u0627\u0646", -1, 1),
        Among("\u0627\u06CC\u0645\u0627\u0646", -1, 1),
        Among("\u0633\u0644\u06CC\u0645\u0627\u0646", -1, 1),
        Among("\u067E\u06CC\u0645\u0627\u0646", -1, 1),
        Among("\u0644\u0628\u0646\u0627\u0646", -1, 1),
        Among("\u06CC\u0648\u0646\u0627\u0646", -1, 1),
        Among("\u0627\u0635\u0641\u0647\u0627\u0646", -1, 1),
        Among("\u0627\u0645\u06A9\u0627\u0646", -1, 1),
        Among("\u067E\u0627\u06CC\u0627\u0646", -1, 1),
        Among("\u0628\u06CC\u0627\u0646", -1, 1),
        Among("\u062C\u0631\u06CC\u0627\u0646", -1, 1)
    ]

    a_4 = [
        Among("\u0627\u0633\u0627\u062A\u06CC\u062F", -1, 2),
        Among("\u0627\u062E\u0628\u0627\u0631", -1, 1)
    ]
    as_4 = ("\u062E\u0628\u0631", "\u0627\u0633\u062A\u0627\u062F")

    a_5 = [
        Among("\u0647\u0627", -1, 1),
        Among("\u0627\u062A", -1, 1),
        Among("\u06CC\u062A", -1, 1),
        Among("\u0645\u0646\u062F", -1, 1),
        Among("\u0648\u0627\u0631", -1, 1),
        Among("\u06AF\u0627\u0631", -1, 1),
        Among("\u062A\u0631", -1, 2),
        Among("\u0627\u0634", -1, 1),
        Among("\u0627\u0645", -1, 1),
        Among("\u0627\u0646", -1, 1),
        Among("\u0628\u0627\u0646", 9, 1),
        Among("\u06AF\u0627\u0646", 9, 1),
        Among("\u06CC\u0627\u0646", 9, 1),
        Among("\u06CC\u0646", -1, 1),
        Among("\u062A\u0631\u06CC\u0646", 13, 1),
        Among("\u06AF\u0627\u0647", -1, 1),
        Among("\u0627\u0646\u0647", -1, 1),
        Among("\u0646\u0627\u06A9", -1, 1),
        Among("\u0647\u0627\u06CC", -1, 1),
        Among("\u0627\u0646\u06CC", -1, 1),
        Among("\u06AF\u06CC", -1, 1),
        Among("\u06CC\u06CC", -1, 1)
    ]

    a_6 = [
        Among("\u0627\u0633\u062A", -1, 1),
        Among("\u0627\u0646\u062F", -1, 1),
        Among("\u06CC\u062F", -1, 1),
        Among("\u0627\u06CC\u062F", 2, 1),
        Among("\u0627\u0633", -1, 1),
        Among("\u06CC\u0645", -1, 1),
        Among("\u0627\u06CC\u0645", 5, 1),
        Among("\u0627\u06CC", -1, 1)
    ]

    a_7 = [
        Among("\u062F", -1, 1),
        Among("\u0627\u0646\u062F", 0, 1),
        Among("\u0631\u0641\u062A\u0627\u0646\u062F", 1, 2),
        Among("\u06CC\u062F", 0, 1),
        Among("\u0631\u0641\u062A\u06CC\u062F", 3, 2),
        Among("\u0645", -1, 1),
        Among("\u0627\u0645", 5, 1),
        Among("\u0631\u0641\u062A\u0645", 5, 2),
        Among("\u06CC\u0645", 5, 1),
        Among("\u0631\u0641\u062A\u06CC\u0645", 8, 2),
        Among("\u0627\u0646", -1, 3),
        Among("\u062A\u0647", -1, 5),
        Among("\u062F\u0647", -1, 4),
        Among("\u0646\u062F\u0647", 12, 3),
        Among("\u0631\u0641\u062A\u06CC", -1, 2)
    ]


class lab0(BaseException): pass


class lab1(BaseException): pass
