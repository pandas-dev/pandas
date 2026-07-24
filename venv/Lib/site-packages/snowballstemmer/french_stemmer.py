# Generated from french.sbl by Snowball 3.1.1 - https://snowballstem.org/

from .basestemmer import BaseStemmer
from .among import Among


class FrenchStemmer(BaseStemmer):
    '''
    This class implements the stemming algorithm defined by a snowball script.
    Generated from french.sbl by Snowball 3.1.1 - https://snowballstem.org/
    '''

    g_v = {"a", "e", "i", "o", "u", "y", "à", "â", "è", "é", "ê", "ë", "î", "ï", "ô", "ù", "û"}

    g_oux_ending = {"b", "h", "j", "l", "n", "p"}

    g_elision_char = {"c", "d", "j", "l", "m", "n", "s", "t"}

    g_keep_with_s = {"a", "i", "o", "s", "u", "è"}

    I_p2 = 0
    I_p1 = 0
    I_pV = 0

    def __r_elisions(self):
        self.bra = self.cursor
        while True:
            try:
                if not self.in_grouping(FrenchStemmer.g_elision_char):
                    raise lab0()
                break
            except lab0: pass
            if not self.eq_s("qu"):
                return False
            break
        if self.cursor == self.limit or self.current[self.cursor] != "'":
            return False
        self.cursor += 1
        self.ket = self.cursor
        if self.cursor >= self.limit:
            return False
        self.slice_del()
        return True

    def __r_prelude(self):
        while True:
            v_1 = self.cursor
            try:
                while True:
                    v_2 = self.cursor
                    try:
                        while True:
                            v_3 = self.cursor
                            try:
                                if not self.in_grouping(FrenchStemmer.g_v):
                                    raise lab2()
                                self.bra = self.cursor
                                while True:
                                    v_4 = self.cursor
                                    try:
                                        if self.cursor == self.limit or self.current[self.cursor] != "u":
                                            raise lab3()
                                        self.cursor += 1
                                        self.ket = self.cursor
                                        if not self.in_grouping(FrenchStemmer.g_v):
                                            raise lab3()
                                        self.slice_from("U")
                                        break
                                    except lab3: pass
                                    self.cursor = v_4
                                    try:
                                        if self.cursor == self.limit or self.current[self.cursor] != "i":
                                            raise lab3()
                                        self.cursor += 1
                                        self.ket = self.cursor
                                        if not self.in_grouping(FrenchStemmer.g_v):
                                            raise lab3()
                                        self.slice_from("I")
                                        break
                                    except lab3: pass
                                    self.cursor = v_4
                                    if self.cursor == self.limit or self.current[self.cursor] != "y":
                                        raise lab2()
                                    self.cursor += 1
                                    self.ket = self.cursor
                                    self.slice_from("Y")
                                    break
                                break
                            except lab2: pass
                            self.cursor = v_3
                            try:
                                self.bra = self.cursor
                                if self.cursor == self.limit or self.current[self.cursor] != "ë":
                                    raise lab2()
                                self.cursor += 1
                                self.ket = self.cursor
                                self.slice_from("He")
                                break
                            except lab2: pass
                            self.cursor = v_3
                            try:
                                self.bra = self.cursor
                                if self.cursor == self.limit or self.current[self.cursor] != "ï":
                                    raise lab2()
                                self.cursor += 1
                                self.ket = self.cursor
                                self.slice_from("Hi")
                                break
                            except lab2: pass
                            self.cursor = v_3
                            try:
                                self.bra = self.cursor
                                if self.cursor == self.limit or self.current[self.cursor] != "y":
                                    raise lab2()
                                self.cursor += 1
                                self.ket = self.cursor
                                if not self.in_grouping(FrenchStemmer.g_v):
                                    raise lab2()
                                self.slice_from("Y")
                                break
                            except lab2: pass
                            self.cursor = v_3
                            if self.cursor == self.limit or self.current[self.cursor] != "q":
                                raise lab1()
                            self.cursor += 1
                            self.bra = self.cursor
                            if self.cursor == self.limit or self.current[self.cursor] != "u":
                                raise lab1()
                            self.cursor += 1
                            self.ket = self.cursor
                            self.slice_from("U")
                            break
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

    def __r_mark_regions(self):
        self.I_pV = self.limit
        self.I_p1 = self.limit
        self.I_p2 = self.limit
        v_1 = self.cursor
        try:
            while True:
                v_2 = self.cursor
                try:
                    if not self.in_grouping(FrenchStemmer.g_v):
                        raise lab1()
                    if not self.in_grouping(FrenchStemmer.g_v):
                        raise lab1()
                    if self.cursor >= self.limit:
                        raise lab1()
                    self.cursor += 1
                    break
                except lab1: pass
                self.cursor = v_2
                try:
                    among_var = self.find_among(FrenchStemmer.a_0)
                    if among_var == 0:
                        raise lab1()
                    if among_var == 1:
                        if not self.in_grouping(FrenchStemmer.g_v):
                            raise lab1()
                    break
                except lab1: pass
                self.cursor = v_2
                if self.cursor >= self.limit:
                    raise lab0()
                self.cursor += 1
                if not self.go_out_grouping(FrenchStemmer.g_v):
                    raise lab0()
                self.cursor += 1
                break
            self.I_pV = self.cursor
        except lab0: pass
        self.cursor = v_1
        v_3 = self.cursor
        try:
            if not self.go_out_grouping(FrenchStemmer.g_v):
                raise lab0()
            self.cursor += 1
            if not self.go_in_grouping(FrenchStemmer.g_v):
                raise lab0()
            self.cursor += 1
            self.I_p1 = self.cursor
            if not self.go_out_grouping(FrenchStemmer.g_v):
                raise lab0()
            self.cursor += 1
            if not self.go_in_grouping(FrenchStemmer.g_v):
                raise lab0()
            self.cursor += 1
            self.I_p2 = self.cursor
        except lab0: pass
        self.cursor = v_3
        return True

    def __r_postlude(self):
        while True:
            v_1 = self.cursor
            try:
                self.bra = self.cursor
                among_var = self.find_among(FrenchStemmer.a_1)
                self.ket = self.cursor
                if among_var == 1:
                    self.slice_from("i")
                elif among_var == 2:
                    self.slice_from("u")
                elif among_var == 3:
                    self.slice_from("y")
                elif among_var == 4:
                    self.slice_from("ë")
                elif among_var == 5:
                    self.slice_from("ï")
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

    def __r_RV(self):
        return self.I_pV <= self.cursor

    def __r_R1(self):
        return self.I_p1 <= self.cursor

    def __r_R2(self):
        return self.I_p2 <= self.cursor

    def __r_standard_suffix(self):
        self.ket = self.cursor
        among_var = self.find_among_b(FrenchStemmer.a_4)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var == 1:
            if not self.__r_R2():
                return False
            self.slice_del()
        elif among_var == 2:
            if not self.__r_R2():
                return False
            self.slice_del()
            v_1 = self.limit - self.cursor
            try:
                self.ket = self.cursor
                if not self.eq_s_b("ic"):
                    self.cursor = self.limit - v_1
                    raise lab0()
                self.bra = self.cursor
                while True:
                    v_2 = self.limit - self.cursor
                    try:
                        if not self.__r_R2():
                            raise lab1()
                        self.slice_del()
                        break
                    except lab1: pass
                    self.cursor = self.limit - v_2
                    self.slice_from("iqU")
                    break
            except lab0: pass
        elif among_var == 3:
            if not self.__r_R2():
                return False
            self.slice_from("log")
        elif among_var == 4:
            if not self.__r_R2():
                return False
            self.slice_from("u")
        elif among_var == 5:
            if not self.__r_R2():
                return False
            self.slice_from("ent")
        elif among_var == 6:
            if not self.__r_RV():
                return False
            self.slice_del()
            v_3 = self.limit - self.cursor
            try:
                self.ket = self.cursor
                among_var = self.find_among_b(FrenchStemmer.a_2)
                if among_var == 0:
                    self.cursor = self.limit - v_3
                    raise lab0()
                self.bra = self.cursor
                if among_var == 1:
                    if not self.__r_R2():
                        self.cursor = self.limit - v_3
                        raise lab0()
                    self.slice_del()
                    self.ket = self.cursor
                    if not self.eq_s_b("at"):
                        self.cursor = self.limit - v_3
                        raise lab0()
                    self.bra = self.cursor
                    if not self.__r_R2():
                        self.cursor = self.limit - v_3
                        raise lab0()
                    self.slice_del()
                elif among_var == 2:
                    while True:
                        v_4 = self.limit - self.cursor
                        try:
                            if not self.__r_R2():
                                raise lab1()
                            self.slice_del()
                            break
                        except lab1: pass
                        self.cursor = self.limit - v_4
                        if not self.__r_R1():
                            self.cursor = self.limit - v_3
                            raise lab0()
                        self.slice_from("eux")
                        break
                elif among_var == 3:
                    if not self.__r_R2():
                        self.cursor = self.limit - v_3
                        raise lab0()
                    self.slice_del()
                else:
                    if not self.__r_RV():
                        self.cursor = self.limit - v_3
                        raise lab0()
                    self.slice_from("i")
            except lab0: pass
        elif among_var == 7:
            if not self.__r_R2():
                return False
            self.slice_del()
            v_5 = self.limit - self.cursor
            try:
                self.ket = self.cursor
                among_var = self.find_among_b(FrenchStemmer.a_3)
                if among_var == 0:
                    self.cursor = self.limit - v_5
                    raise lab0()
                self.bra = self.cursor
                if among_var == 1:
                    while True:
                        v_6 = self.limit - self.cursor
                        try:
                            if not self.__r_R2():
                                raise lab1()
                            self.slice_del()
                            break
                        except lab1: pass
                        self.cursor = self.limit - v_6
                        self.slice_from("abl")
                        break
                elif among_var == 2:
                    while True:
                        v_7 = self.limit - self.cursor
                        try:
                            if not self.__r_R2():
                                raise lab1()
                            self.slice_del()
                            break
                        except lab1: pass
                        self.cursor = self.limit - v_7
                        self.slice_from("iqU")
                        break
                else:
                    if not self.__r_R2():
                        self.cursor = self.limit - v_5
                        raise lab0()
                    self.slice_del()
            except lab0: pass
        elif among_var == 8:
            if not self.__r_R2():
                return False
            self.slice_del()
            v_8 = self.limit - self.cursor
            try:
                self.ket = self.cursor
                if not self.eq_s_b("at"):
                    self.cursor = self.limit - v_8
                    raise lab0()
                self.bra = self.cursor
                if not self.__r_R2():
                    self.cursor = self.limit - v_8
                    raise lab0()
                self.slice_del()
                self.ket = self.cursor
                if not self.eq_s_b("ic"):
                    self.cursor = self.limit - v_8
                    raise lab0()
                self.bra = self.cursor
                while True:
                    v_9 = self.limit - self.cursor
                    try:
                        if not self.__r_R2():
                            raise lab1()
                        self.slice_del()
                        break
                    except lab1: pass
                    self.cursor = self.limit - v_9
                    self.slice_from("iqU")
                    break
            except lab0: pass
        elif among_var == 9:
            self.slice_from("eau")
        elif among_var == 10:
            if not self.__r_R1():
                return False
            self.slice_from("al")
        elif among_var == 11:
            if not self.in_grouping_b(FrenchStemmer.g_oux_ending):
                return False
            self.slice_from("ou")
        elif among_var == 12:
            while True:
                v_10 = self.limit - self.cursor
                try:
                    if not self.__r_R2():
                        raise lab0()
                    self.slice_del()
                    break
                except lab0: pass
                self.cursor = self.limit - v_10
                if not self.__r_R1():
                    return False
                self.slice_from("eux")
                break
        elif among_var == 13:
            if not self.__r_R1():
                return False
            if not self.out_grouping_b(FrenchStemmer.g_v):
                return False
            self.slice_del()
        elif among_var == 14:
            if not self.__r_RV():
                return False
            self.slice_from("ant")
            return False
        elif among_var == 15:
            if not self.__r_RV():
                return False
            self.slice_from("ent")
            return False
        else:
            v_11 = self.limit - self.cursor
            if not self.in_grouping_b(FrenchStemmer.g_v):
                return False
            if not self.__r_RV():
                return False
            self.cursor = self.limit - v_11
            self.slice_del()
            return False
        return True

    def __r_i_verb_suffix(self):
        if self.cursor < self.I_pV:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_pV
        self.ket = self.cursor
        if self.find_among_b(FrenchStemmer.a_5) == 0:
            self.limit_backward = v_2
            return False
        self.bra = self.cursor
        try:
            if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "H":
                raise lab0()
            self.cursor -= 1
            self.limit_backward = v_2
            return False
        except lab0: pass
        if not self.out_grouping_b(FrenchStemmer.g_v):
            self.limit_backward = v_2
            return False
        self.slice_del()
        self.limit_backward = v_2
        return True

    def __r_verb_suffix(self):
        if self.cursor < self.I_pV:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_pV
        self.ket = self.cursor
        among_var = self.find_among_b(FrenchStemmer.a_7)
        if among_var == 0:
            self.limit_backward = v_2
            return False
        self.bra = self.cursor
        self.limit_backward = v_2
        if among_var == 1:
            if not self.__r_R2():
                return False
            self.slice_del()
        elif among_var == 2:
            self.slice_del()
        elif among_var == 3:
            v_3 = self.limit - self.cursor
            try:
                if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "e":
                    self.cursor = self.limit - v_3
                    raise lab0()
                self.cursor -= 1
                if not self.__r_RV():
                    self.cursor = self.limit - v_3
                    raise lab0()
                self.bra = self.cursor
            except lab0: pass
            self.slice_del()
        else:
            v_4 = self.limit - self.cursor
            try:
                among_var = self.find_among_b(FrenchStemmer.a_6)
                if among_var == 0:
                    raise lab0()
                if among_var == 1:
                    if self.cursor <= self.limit_backward:
                        raise lab0()
                    self.cursor -= 1
                    if self.cursor > self.limit_backward:
                        raise lab0()
                return False
            except lab0: pass
            self.cursor = self.limit - v_4
            self.slice_del()
        return True

    def __r_residual_suffix(self):
        v_1 = self.limit - self.cursor
        try:
            self.ket = self.cursor
            if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "s":
                self.cursor = self.limit - v_1
                raise lab0()
            self.cursor -= 1
            self.bra = self.cursor
            v_2 = self.limit - self.cursor
            while True:
                try:
                    if not self.eq_s_b("Hi"):
                        raise lab1()
                    break
                except lab1: pass
                if not self.out_grouping_b(FrenchStemmer.g_keep_with_s):
                    self.cursor = self.limit - v_1
                    raise lab0()
                break
            self.cursor = self.limit - v_2
            self.slice_del()
        except lab0: pass
        if self.cursor < self.I_pV:
            return False
        v_4 = self.limit_backward
        self.limit_backward = self.I_pV
        self.ket = self.cursor
        among_var = self.find_among_b(FrenchStemmer.a_8)
        if among_var == 0:
            self.limit_backward = v_4
            return False
        self.bra = self.cursor
        if among_var == 1:
            if not self.__r_R2():
                self.limit_backward = v_4
                return False
            while True:
                try:
                    if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "s":
                        raise lab0()
                    self.cursor -= 1
                    break
                except lab0: pass
                if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "t":
                    self.limit_backward = v_4
                    return False
                self.cursor -= 1
                break
            self.slice_del()
        elif among_var == 2:
            self.slice_from("i")
        else:
            self.slice_del()
        self.limit_backward = v_4
        return True

    def __r_un_double(self):
        v_1 = self.limit - self.cursor
        if self.find_among_b(FrenchStemmer.a_9) == 0:
            return False
        self.cursor = self.limit - v_1
        self.ket = self.cursor
        if self.cursor <= self.limit_backward:
            return False
        self.cursor -= 1
        self.bra = self.cursor
        self.slice_del()
        return True

    def __r_un_accent(self):
        v_1 = 1
        while True:
            try:
                if not self.out_grouping_b(FrenchStemmer.g_v):
                    raise lab0()
                v_1 -= 1
                continue
            except lab0: pass
            break
        if v_1 > 0:
            return False
        self.ket = self.cursor
        while True:
            try:
                if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "é":
                    raise lab0()
                self.cursor -= 1
                break
            except lab0: pass
            if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "è":
                return False
            self.cursor -= 1
            break
        self.bra = self.cursor
        self.slice_from("e")
        return True

    def _stem(self):
        v_1 = self.cursor
        self.__r_elisions()
        self.cursor = v_1
        v_2 = self.cursor
        self.__r_prelude()
        self.cursor = v_2
        self.__r_mark_regions()
        self.limit_backward = self.cursor
        self.cursor = self.limit
        v_3 = self.limit - self.cursor
        try:
            while True:
                v_4 = self.limit - self.cursor
                try:
                    v_5 = self.limit - self.cursor
                    while True:
                        v_6 = self.limit - self.cursor
                        try:
                            if not self.__r_standard_suffix():
                                raise lab2()
                            break
                        except lab2: pass
                        self.cursor = self.limit - v_6
                        try:
                            if not self.__r_i_verb_suffix():
                                raise lab2()
                            break
                        except lab2: pass
                        self.cursor = self.limit - v_6
                        if not self.__r_verb_suffix():
                            raise lab1()
                        break
                    self.cursor = self.limit - v_5
                    v_7 = self.limit - self.cursor
                    try:
                        self.ket = self.cursor
                        while True:
                            v_8 = self.limit - self.cursor
                            try:
                                if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "Y":
                                    raise lab3()
                                self.cursor -= 1
                                self.bra = self.cursor
                                self.slice_from("i")
                                break
                            except lab3: pass
                            self.cursor = self.limit - v_8
                            if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "ç":
                                self.cursor = self.limit - v_7
                                raise lab2()
                            self.cursor -= 1
                            self.bra = self.cursor
                            self.slice_from("c")
                            break
                    except lab2: pass
                    break
                except lab1: pass
                self.cursor = self.limit - v_4
                if not self.__r_residual_suffix():
                    raise lab0()
                break
        except lab0: pass
        self.cursor = self.limit - v_3
        v_9 = self.limit - self.cursor
        self.__r_un_double()
        self.cursor = self.limit - v_9
        v_10 = self.limit - self.cursor
        self.__r_un_accent()
        self.cursor = self.limit - v_10
        self.cursor = self.limit_backward
        v_11 = self.cursor
        self.__r_postlude()
        self.cursor = v_11
        return True

    a_0 = [
        Among("col", -1, -1),
        Among("ni", -1, 1),
        Among("par", -1, -1),
        Among("tap", -1, -1)
    ]

    a_1 = [
        Among("", -1, 7),
        Among("H", 0, 6),
        Among("He", 1, 4),
        Among("Hi", 1, 5),
        Among("I", 0, 1),
        Among("U", 0, 2),
        Among("Y", 0, 3)
    ]

    a_2 = [
        Among("iqU", -1, 3),
        Among("abl", -1, 3),
        Among("Ièr", -1, 4),
        Among("ièr", -1, 4),
        Among("eus", -1, 2),
        Among("iv", -1, 1)
    ]

    a_3 = [
        Among("ic", -1, 2),
        Among("abil", -1, 1),
        Among("iv", -1, 3)
    ]

    a_4 = [
        Among("iqUe", -1, 1),
        Among("atrice", -1, 2),
        Among("ance", -1, 1),
        Among("ence", -1, 5),
        Among("logie", -1, 3),
        Among("able", -1, 1),
        Among("isme", -1, 1),
        Among("euse", -1, 12),
        Among("iste", -1, 1),
        Among("ive", -1, 8),
        Among("if", -1, 8),
        Among("usion", -1, 4),
        Among("ation", -1, 2),
        Among("ution", -1, 4),
        Among("ateur", -1, 2),
        Among("iqUes", -1, 1),
        Among("atrices", -1, 2),
        Among("ances", -1, 1),
        Among("ences", -1, 5),
        Among("logies", -1, 3),
        Among("ables", -1, 1),
        Among("ismes", -1, 1),
        Among("euses", -1, 12),
        Among("istes", -1, 1),
        Among("ives", -1, 8),
        Among("ifs", -1, 8),
        Among("usions", -1, 4),
        Among("ations", -1, 2),
        Among("utions", -1, 4),
        Among("ateurs", -1, 2),
        Among("ments", -1, 16),
        Among("ements", 30, 6),
        Among("issements", 31, 13),
        Among("ités", -1, 7),
        Among("ment", -1, 16),
        Among("ement", 34, 6),
        Among("issement", 35, 13),
        Among("amment", 34, 14),
        Among("emment", 34, 15),
        Among("aux", -1, 10),
        Among("eaux", 39, 9),
        Among("eux", -1, 1),
        Among("oux", -1, 11),
        Among("ité", -1, 7)
    ]

    a_5 = [
        Among("ira", -1, 1),
        Among("ie", -1, 1),
        Among("isse", -1, 1),
        Among("issante", -1, 1),
        Among("i", -1, 1),
        Among("irai", 4, 1),
        Among("ir", -1, 1),
        Among("iras", -1, 1),
        Among("ies", -1, 1),
        Among("îmes", -1, 1),
        Among("isses", -1, 1),
        Among("issantes", -1, 1),
        Among("îtes", -1, 1),
        Among("is", -1, 1),
        Among("irais", 13, 1),
        Among("issais", 13, 1),
        Among("irions", -1, 1),
        Among("issions", -1, 1),
        Among("irons", -1, 1),
        Among("issons", -1, 1),
        Among("issants", -1, 1),
        Among("it", -1, 1),
        Among("irait", 21, 1),
        Among("issait", 21, 1),
        Among("issant", -1, 1),
        Among("iraIent", -1, 1),
        Among("issaIent", -1, 1),
        Among("irent", -1, 1),
        Among("issent", -1, 1),
        Among("iront", -1, 1),
        Among("ît", -1, 1),
        Among("iriez", -1, 1),
        Among("issiez", -1, 1),
        Among("irez", -1, 1),
        Among("issez", -1, 1)
    ]

    a_6 = [
        Among("al", -1, 1),
        Among("épl", -1, -1),
        Among("auv", -1, -1)
    ]

    a_7 = [
        Among("a", -1, 3),
        Among("era", 0, 2),
        Among("aise", -1, 4),
        Among("asse", -1, 3),
        Among("ante", -1, 3),
        Among("ée", -1, 2),
        Among("ai", -1, 3),
        Among("erai", 6, 2),
        Among("er", -1, 2),
        Among("as", -1, 3),
        Among("eras", 9, 2),
        Among("âmes", -1, 3),
        Among("aises", -1, 4),
        Among("asses", -1, 3),
        Among("antes", -1, 3),
        Among("âtes", -1, 3),
        Among("ées", -1, 2),
        Among("ais", -1, 4),
        Among("eais", 17, 2),
        Among("erais", 17, 2),
        Among("ions", -1, 1),
        Among("erions", 20, 2),
        Among("assions", 20, 3),
        Among("erons", -1, 2),
        Among("ants", -1, 3),
        Among("és", -1, 2),
        Among("ait", -1, 3),
        Among("erait", 26, 2),
        Among("ant", -1, 3),
        Among("aIent", -1, 3),
        Among("eraIent", 29, 2),
        Among("èrent", -1, 2),
        Among("assent", -1, 3),
        Among("eront", -1, 2),
        Among("ât", -1, 3),
        Among("ez", -1, 2),
        Among("iez", 35, 2),
        Among("eriez", 36, 2),
        Among("assiez", 36, 3),
        Among("erez", 35, 2),
        Among("é", -1, 2)
    ]

    a_8 = [
        Among("e", -1, 3),
        Among("Ière", 0, 2),
        Among("ière", 0, 2),
        Among("ion", -1, 1),
        Among("Ier", -1, 2),
        Among("ier", -1, 2)
    ]

    a_9 = [
        Among("ell", -1, -1),
        Among("eill", -1, -1),
        Among("enn", -1, -1),
        Among("onn", -1, -1),
        Among("ett", -1, -1)
    ]


class lab0(BaseException): pass


class lab1(BaseException): pass


class lab2(BaseException): pass


class lab3(BaseException): pass
