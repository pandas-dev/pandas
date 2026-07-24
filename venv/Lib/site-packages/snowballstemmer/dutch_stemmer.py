# Generated from dutch.sbl by Snowball 3.1.1 - https://snowballstem.org/

from .basestemmer import BaseStemmer
from .among import Among


class DutchStemmer(BaseStemmer):
    '''
    This class implements the stemming algorithm defined by a snowball script.
    Generated from dutch.sbl by Snowball 3.1.1 - https://snowballstem.org/
    '''

    g_E = {"e", "è", "é", "ê", "ë"}

    g_AIOU = {"a", "i", "o", "u", "à", "á", "â", "ä", "ì", "í", "î", "ï", "ò", "ó", "ô", "ö", "ù", "ú", "û", "ü"}

    g_AEIOU = {"a", "e", "i", "o", "u", "à", "á", "â", "ä", "è", "é", "ê", "ë", "ì", "í", "î", "ï", "ò", "ó", "ô", "ö", "ù", "ú", "û", "ü"}

    g_v = {"a", "e", "i", "o", "u", "y", "à", "á", "â", "ä", "è", "é", "ê", "ë", "ì", "í", "î", "ï", "ò", "ó", "ô", "ö", "ù", "ú", "û", "ü"}

    g_v_WX = {"a", "e", "i", "o", "u", "w", "x", "y", "à", "á", "â", "ä", "è", "é", "ê", "ë", "ì", "í", "î", "ï", "ò", "ó", "ô", "ö", "ù", "ú", "û", "ü"}

    B_GE_removed = False
    I_p2 = 0
    I_p1 = 0

    def __r_R1(self):
        return self.I_p1 <= self.cursor

    def __r_R2(self):
        return self.I_p2 <= self.cursor

    def __r_V(self):
        v_1 = self.limit - self.cursor
        while True:
            try:
                if not self.in_grouping_b(DutchStemmer.g_v):
                    raise lab0()
                break
            except lab0: pass
            if not self.eq_s_b("ij"):
                return False
            break
        self.cursor = self.limit - v_1
        return True

    def __r_VX(self):
        v_1 = self.limit - self.cursor
        if self.cursor <= self.limit_backward:
            return False
        self.cursor -= 1
        while True:
            try:
                if not self.in_grouping_b(DutchStemmer.g_v):
                    raise lab0()
                break
            except lab0: pass
            if not self.eq_s_b("ij"):
                return False
            break
        self.cursor = self.limit - v_1
        return True

    def __r_C(self):
        v_1 = self.limit - self.cursor
        try:
            if not self.eq_s_b("ij"):
                raise lab0()
            return False
        except lab0: pass
        if not self.out_grouping_b(DutchStemmer.g_v):
            return False
        self.cursor = self.limit - v_1
        return True

    def __r_lengthen_V(self):
        v_1 = self.limit - self.cursor
        try:
            if not self.out_grouping_b(DutchStemmer.g_v_WX):
                raise lab0()
            self.ket = self.cursor
            among_var = self.find_among_b(DutchStemmer.a_0)
            if among_var == 0:
                raise lab0()
            self.bra = self.cursor
            if among_var == 1:
                v_2 = self.limit - self.cursor
                while True:
                    try:
                        if not self.out_grouping_b(DutchStemmer.g_AEIOU):
                            raise lab1()
                        break
                    except lab1: pass
                    if self.cursor > self.limit_backward:
                        raise lab0()
                    break
                self.cursor = self.limit - v_2
                S_ch = self.slice_to()
                c = self.cursor
                self.insert(self.cursor, self.cursor, S_ch)
                self.cursor = c
            elif among_var == 2:
                v_3 = self.limit - self.cursor
                while True:
                    try:
                        if not self.out_grouping_b(DutchStemmer.g_AEIOU):
                            raise lab1()
                        break
                    except lab1: pass
                    if self.cursor > self.limit_backward:
                        raise lab0()
                    break
                v_4 = self.limit - self.cursor
                try:
                    while True:
                        try:
                            if not self.in_grouping_b(DutchStemmer.g_AIOU):
                                raise lab2()
                            break
                        except lab2: pass
                        if not self.in_grouping_b(DutchStemmer.g_E):
                            raise lab1()
                        if self.cursor > self.limit_backward:
                            raise lab1()
                        break
                    raise lab0()
                except lab1: pass
                self.cursor = self.limit - v_4
                v_5 = self.limit - self.cursor
                try:
                    if self.cursor <= self.limit_backward:
                        raise lab1()
                    self.cursor -= 1
                    if not self.in_grouping_b(DutchStemmer.g_AIOU):
                        raise lab1()
                    if not self.out_grouping_b(DutchStemmer.g_AEIOU):
                        raise lab1()
                    raise lab0()
                except lab1: pass
                self.cursor = self.limit - v_5
                self.cursor = self.limit - v_3
                S_ch = self.slice_to()
                c = self.cursor
                self.insert(self.cursor, self.cursor, S_ch)
                self.cursor = c
            elif among_var == 3:
                self.slice_from("eëe")
            else:
                self.slice_from("iee")
        except lab0: pass
        self.cursor = self.limit - v_1
        return True

    def __r_Step_1(self):
        self.ket = self.cursor
        among_var = self.find_among_b(DutchStemmer.a_1)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var == 1:
            self.slice_del()
        elif among_var == 2:
            if not self.__r_R1():
                return False
            v_1 = self.limit - self.cursor
            try:
                if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "t":
                    raise lab0()
                self.cursor -= 1
                if not self.__r_R1():
                    raise lab0()
                return False
            except lab0: pass
            self.cursor = self.limit - v_1
            if not self.__r_C():
                return False
            self.slice_del()
        elif among_var == 3:
            if not self.__r_R1():
                return False
            self.slice_from("ie")
        elif among_var == 4:
            while True:
                v_2 = self.limit - self.cursor
                try:
                    v_3 = self.limit - self.cursor
                    if not self.eq_s_b("ar"):
                        raise lab0()
                    if not self.__r_R1():
                        raise lab0()
                    if not self.__r_C():
                        raise lab0()
                    self.cursor = self.limit - v_3
                    self.slice_del()
                    self.__r_lengthen_V()
                    break
                except lab0: pass
                self.cursor = self.limit - v_2
                try:
                    v_4 = self.limit - self.cursor
                    if not self.eq_s_b("er"):
                        raise lab0()
                    if not self.__r_R1():
                        raise lab0()
                    if not self.__r_C():
                        raise lab0()
                    self.cursor = self.limit - v_4
                    self.slice_del()
                    break
                except lab0: pass
                self.cursor = self.limit - v_2
                if not self.__r_R1():
                    return False
                if not self.__r_C():
                    return False
                self.slice_from("e")
                break
        elif among_var == 5:
            if not self.__r_R1():
                return False
            self.slice_from("é")
        elif among_var == 6:
            if not self.__r_R1():
                return False
            if not self.__r_V():
                return False
            self.slice_from("au")
        elif among_var == 7:
            while True:
                v_5 = self.limit - self.cursor
                try:
                    if not self.eq_s_b("hed"):
                        raise lab0()
                    if not self.__r_R1():
                        raise lab0()
                    self.bra = self.cursor
                    self.slice_from("heid")
                    break
                except lab0: pass
                self.cursor = self.limit - v_5
                try:
                    if not self.eq_s_b("nd"):
                        raise lab0()
                    self.slice_del()
                    break
                except lab0: pass
                self.cursor = self.limit - v_5
                try:
                    if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "d":
                        raise lab0()
                    self.cursor -= 1
                    if not self.__r_R1():
                        raise lab0()
                    if not self.__r_C():
                        raise lab0()
                    self.bra = self.cursor
                    self.slice_del()
                    break
                except lab0: pass
                self.cursor = self.limit - v_5
                try:
                    while True:
                        try:
                            if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "i":
                                raise lab1()
                            self.cursor -= 1
                            break
                        except lab1: pass
                        if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "j":
                            raise lab0()
                        self.cursor -= 1
                        break
                    if not self.__r_V():
                        raise lab0()
                    self.slice_del()
                    break
                except lab0: pass
                self.cursor = self.limit - v_5
                if not self.__r_R1():
                    return False
                if not self.__r_C():
                    return False
                self.slice_del()
                self.__r_lengthen_V()
                break
        else:
            self.slice_from("nd")
        return True

    def __r_Step_2(self):
        self.ket = self.cursor
        among_var = self.find_among_b(DutchStemmer.a_2)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var == 1:
            while True:
                v_1 = self.limit - self.cursor
                try:
                    if not self.eq_s_b("'t"):
                        raise lab0()
                    self.bra = self.cursor
                    self.slice_del()
                    break
                except lab0: pass
                self.cursor = self.limit - v_1
                try:
                    if not self.eq_s_b("et"):
                        raise lab0()
                    self.bra = self.cursor
                    if not self.__r_R1():
                        raise lab0()
                    if not self.__r_C():
                        raise lab0()
                    self.slice_del()
                    break
                except lab0: pass
                self.cursor = self.limit - v_1
                try:
                    if not self.eq_s_b("rnt"):
                        raise lab0()
                    self.bra = self.cursor
                    self.slice_from("rn")
                    break
                except lab0: pass
                self.cursor = self.limit - v_1
                try:
                    if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "t":
                        raise lab0()
                    self.cursor -= 1
                    self.bra = self.cursor
                    if not self.__r_R1():
                        raise lab0()
                    if not self.__r_VX():
                        raise lab0()
                    self.slice_del()
                    break
                except lab0: pass
                self.cursor = self.limit - v_1
                try:
                    if not self.eq_s_b("ink"):
                        raise lab0()
                    self.bra = self.cursor
                    self.slice_from("ing")
                    break
                except lab0: pass
                self.cursor = self.limit - v_1
                try:
                    if not self.eq_s_b("mp"):
                        raise lab0()
                    self.bra = self.cursor
                    self.slice_from("m")
                    break
                except lab0: pass
                self.cursor = self.limit - v_1
                try:
                    if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "'":
                        raise lab0()
                    self.cursor -= 1
                    self.bra = self.cursor
                    if not self.__r_R1():
                        raise lab0()
                    self.slice_del()
                    break
                except lab0: pass
                self.cursor = self.limit - v_1
                self.bra = self.cursor
                if not self.__r_R1():
                    return False
                if not self.__r_C():
                    return False
                self.slice_del()
                break
        elif among_var == 2:
            if not self.__r_R1():
                return False
            self.slice_from("g")
        elif among_var == 3:
            if not self.__r_R1():
                return False
            self.slice_from("lijk")
        elif among_var == 4:
            if not self.__r_R1():
                return False
            self.slice_from("isch")
        elif among_var == 5:
            if not self.__r_R1():
                return False
            if not self.__r_C():
                return False
            self.slice_del()
        elif among_var == 6:
            if not self.__r_R1():
                return False
            self.slice_from("t")
        elif among_var == 7:
            if not self.__r_R1():
                return False
            self.slice_from("s")
        elif among_var == 8:
            if not self.__r_R1():
                return False
            self.slice_from("r")
        elif among_var == 9:
            if not self.__r_R1():
                return False
            self.slice_del()
            self.insert(self.cursor, self.cursor, "l")
            self.__r_lengthen_V()
        elif among_var == 10:
            if not self.__r_R1():
                return False
            if not self.__r_C():
                return False
            self.slice_del()
            self.insert(self.cursor, self.cursor, "en")
            self.__r_lengthen_V()
        else:
            if not self.__r_R1():
                return False
            if not self.__r_C():
                return False
            self.slice_from("ief")
        return True

    def __r_Step_3(self):
        self.ket = self.cursor
        among_var = self.find_among_b(DutchStemmer.a_3)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var == 1:
            if not self.__r_R1():
                return False
            self.slice_from("eer")
        elif among_var == 2:
            if not self.__r_R1():
                return False
            self.slice_del()
            self.__r_lengthen_V()
        elif among_var == 3:
            if not self.__r_R1():
                return False
            self.slice_del()
        elif among_var == 4:
            self.slice_from("r")
        elif among_var == 5:
            while True:
                v_1 = self.limit - self.cursor
                try:
                    if not self.eq_s_b("ild"):
                        raise lab0()
                    self.slice_from("er")
                    break
                except lab0: pass
                self.cursor = self.limit - v_1
                if not self.__r_R1():
                    return False
                self.slice_del()
                self.__r_lengthen_V()
                break
        elif among_var == 6:
            if not self.__r_R1():
                return False
            if not self.__r_C():
                return False
            self.slice_from("aar")
        elif among_var == 7:
            if not self.__r_R2():
                return False
            self.slice_del()
            self.insert(self.cursor, self.cursor, "f")
            self.__r_lengthen_V()
        elif among_var == 8:
            if not self.__r_R2():
                return False
            self.slice_del()
            self.insert(self.cursor, self.cursor, "g")
            self.__r_lengthen_V()
        elif among_var == 9:
            if not self.__r_R1():
                return False
            if not self.__r_C():
                return False
            self.slice_from("t")
        else:
            if not self.__r_R1():
                return False
            if not self.__r_C():
                return False
            self.slice_from("d")
        return True

    def __r_Step_4(self):
        while True:
            v_1 = self.limit - self.cursor
            try:
                self.ket = self.cursor
                among_var = self.find_among_b(DutchStemmer.a_4)
                if among_var == 0:
                    raise lab0()
                self.bra = self.cursor
                if among_var == 1:
                    if not self.__r_R1():
                        raise lab0()
                    self.slice_from("ie")
                elif among_var == 2:
                    if not self.__r_R1():
                        raise lab0()
                    self.slice_from("eer")
                elif among_var == 3:
                    if not self.__r_R1():
                        raise lab0()
                    self.slice_del()
                elif among_var == 4:
                    if not self.__r_R1():
                        raise lab0()
                    if not self.__r_V():
                        raise lab0()
                    self.slice_from("n")
                elif among_var == 5:
                    if not self.__r_R1():
                        raise lab0()
                    if not self.__r_V():
                        raise lab0()
                    self.slice_from("l")
                elif among_var == 6:
                    if not self.__r_R1():
                        raise lab0()
                    if not self.__r_V():
                        raise lab0()
                    self.slice_from("r")
                elif among_var == 7:
                    if not self.__r_R1():
                        raise lab0()
                    self.slice_from("teer")
                elif among_var == 8:
                    if not self.__r_R1():
                        raise lab0()
                    self.slice_from("lijk")
                else:
                    if not self.__r_R1():
                        raise lab0()
                    if not self.__r_C():
                        raise lab0()
                    self.slice_del()
                    self.__r_lengthen_V()
                break
            except lab0: pass
            self.cursor = self.limit - v_1
            self.ket = self.cursor
            if self.find_among_b(DutchStemmer.a_5) == 0:
                return False
            self.bra = self.cursor
            if not self.__r_R1():
                return False
            v_2 = self.limit - self.cursor
            try:
                if not self.eq_s_b("inn"):
                    raise lab0()
                if self.cursor > self.limit_backward:
                    raise lab0()
                return False
            except lab0: pass
            self.cursor = self.limit - v_2
            if not self.__r_C():
                return False
            self.slice_del()
            self.__r_lengthen_V()
            break
        return True

    def __r_Step_7(self):
        self.ket = self.cursor
        among_var = self.find_among_b(DutchStemmer.a_6)
        if among_var == 0:
            return False
        self.bra = self.cursor
        self.slice_from(DutchStemmer.as_6[among_var - 1])
        return True

    def __r_Step_6(self):
        self.ket = self.cursor
        among_var = self.find_among_b(DutchStemmer.a_7)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var == 1:
            self.slice_from("b")
        elif among_var == 2:
            self.slice_from("c")
        elif among_var == 3:
            self.slice_from("d")
        elif among_var == 4:
            self.slice_from("f")
        elif among_var == 5:
            self.slice_from("g")
        elif among_var == 6:
            self.slice_from("h")
        elif among_var == 7:
            self.slice_from("j")
        elif among_var == 8:
            self.slice_from("k")
        elif among_var == 9:
            self.slice_from("l")
        elif among_var == 10:
            self.slice_from("m")
        elif among_var == 11:
            v_1 = self.limit - self.cursor
            try:
                if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "i":
                    raise lab0()
                self.cursor -= 1
                if self.cursor > self.limit_backward:
                    raise lab0()
                return False
            except lab0: pass
            self.cursor = self.limit - v_1
            self.slice_from("n")
        elif among_var == 12:
            self.slice_from("p")
        elif among_var == 13:
            self.slice_from("q")
        elif among_var == 14:
            self.slice_from("r")
        elif among_var == 15:
            self.slice_from("s")
        elif among_var == 16:
            self.slice_from("t")
        elif among_var == 17:
            self.slice_from("v")
        elif among_var == 18:
            self.slice_from("w")
        elif among_var == 19:
            self.slice_from("x")
        else:
            self.slice_from("z")
        return True

    def __r_Step_1c(self):
        self.ket = self.cursor
        among_var = self.find_among_b(DutchStemmer.a_8)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if not self.__r_R1():
            return False
        if not self.__r_C():
            return False
        if among_var == 1:
            v_1 = self.limit - self.cursor
            try:
                if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "n":
                    raise lab0()
                self.cursor -= 1
                if not self.__r_R1():
                    raise lab0()
                return False
            except lab0: pass
            self.cursor = self.limit - v_1
            while True:
                v_2 = self.limit - self.cursor
                try:
                    if not self.eq_s_b("in"):
                        raise lab0()
                    if self.cursor > self.limit_backward:
                        raise lab0()
                    self.slice_from("n")
                    break
                except lab0: pass
                self.cursor = self.limit - v_2
                self.slice_del()
                break
        else:
            v_3 = self.limit - self.cursor
            try:
                if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "h":
                    raise lab0()
                self.cursor -= 1
                if not self.__r_R1():
                    raise lab0()
                return False
            except lab0: pass
            self.cursor = self.limit - v_3
            v_4 = self.limit - self.cursor
            try:
                if not self.eq_s_b("en"):
                    raise lab0()
                if self.cursor > self.limit_backward:
                    raise lab0()
                return False
            except lab0: pass
            self.cursor = self.limit - v_4
            self.slice_del()
        return True

    def __r_Lose_prefix(self):
        self.bra = self.cursor
        if not self.eq_s("ge"):
            return False
        self.ket = self.cursor
        v_1 = self.cursor
        if self.cursor + 3 > self.limit:
            return False
        self.cursor += 3
        self.cursor = v_1
        v_2 = self.cursor
        while True:
            v_3 = self.cursor
            try:
                while True:
                    try:
                        if not self.eq_s("ij"):
                            raise lab1()
                        break
                    except lab1: pass
                    if not self.in_grouping(DutchStemmer.g_v):
                        raise lab0()
                    break
                break
            except lab0: pass
            self.cursor = v_3
            if self.cursor >= self.limit:
                return False
            self.cursor += 1
        while True:
            v_4 = self.cursor
            try:
                while True:
                    try:
                        if not self.eq_s("ij"):
                            raise lab1()
                        break
                    except lab1: pass
                    if not self.in_grouping(DutchStemmer.g_v):
                        raise lab0()
                    break
                continue
            except lab0: pass
            self.cursor = v_4
            break
        if self.cursor >= self.limit:
            return False
        self.cursor = v_2
        among_var = self.find_among(DutchStemmer.a_9)
        if among_var == 1:
            return False
        self.B_GE_removed = True
        self.slice_del()
        v_5 = self.cursor
        try:
            self.bra = self.cursor
            among_var = self.find_among(DutchStemmer.a_10)
            if among_var == 0:
                raise lab0()
            self.ket = self.cursor
            self.slice_from(DutchStemmer.as_10[among_var - 1])
        except lab0: pass
        self.cursor = v_5
        return True

    def __r_Lose_infix(self):
        if self.cursor >= self.limit:
            return False
        self.cursor += 1
        while True:
            try:
                self.bra = self.cursor
                if not self.eq_s("ge"):
                    raise lab0()
                self.ket = self.cursor
                break
            except lab0: pass
            if self.cursor >= self.limit:
                return False
            self.cursor += 1
        v_1 = self.cursor
        if self.cursor + 3 > self.limit:
            return False
        self.cursor += 3
        self.cursor = v_1
        v_2 = self.cursor
        while True:
            v_3 = self.cursor
            try:
                while True:
                    try:
                        if not self.eq_s("ij"):
                            raise lab1()
                        break
                    except lab1: pass
                    if not self.in_grouping(DutchStemmer.g_v):
                        raise lab0()
                    break
                break
            except lab0: pass
            self.cursor = v_3
            if self.cursor >= self.limit:
                return False
            self.cursor += 1
        while True:
            v_4 = self.cursor
            try:
                while True:
                    try:
                        if not self.eq_s("ij"):
                            raise lab1()
                        break
                    except lab1: pass
                    if not self.in_grouping(DutchStemmer.g_v):
                        raise lab0()
                    break
                continue
            except lab0: pass
            self.cursor = v_4
            break
        if self.cursor >= self.limit:
            return False
        self.cursor = v_2
        self.B_GE_removed = True
        self.slice_del()
        v_5 = self.cursor
        try:
            self.bra = self.cursor
            among_var = self.find_among(DutchStemmer.a_11)
            if among_var == 0:
                raise lab0()
            self.ket = self.cursor
            self.slice_from(DutchStemmer.as_11[among_var - 1])
        except lab0: pass
        self.cursor = v_5
        return True

    def __r_measure(self):
        self.I_p1 = self.limit
        self.I_p2 = self.limit
        v_1 = self.cursor
        try:
            while True:
                try:
                    if not self.out_grouping(DutchStemmer.g_v):
                        raise lab1()
                    continue
                except lab1: pass
                break
            v_2 = 1
            while True:
                v_3 = self.cursor
                try:
                    while True:
                        try:
                            if not self.eq_s("ij"):
                                raise lab2()
                            break
                        except lab2: pass
                        if not self.in_grouping(DutchStemmer.g_v):
                            raise lab1()
                        break
                    v_2 -= 1
                    continue
                except lab1: pass
                self.cursor = v_3
                break
            if v_2 > 0:
                raise lab0()
            if not self.out_grouping(DutchStemmer.g_v):
                raise lab0()
            self.I_p1 = self.cursor
            while True:
                try:
                    if not self.out_grouping(DutchStemmer.g_v):
                        raise lab1()
                    continue
                except lab1: pass
                break
            v_4 = 1
            while True:
                v_5 = self.cursor
                try:
                    while True:
                        try:
                            if not self.eq_s("ij"):
                                raise lab2()
                            break
                        except lab2: pass
                        if not self.in_grouping(DutchStemmer.g_v):
                            raise lab1()
                        break
                    v_4 -= 1
                    continue
                except lab1: pass
                self.cursor = v_5
                break
            if v_4 > 0:
                raise lab0()
            if not self.out_grouping(DutchStemmer.g_v):
                raise lab0()
            self.I_p2 = self.cursor
        except lab0: pass
        self.cursor = v_1
        return True

    def _stem(self):
        B_stemmed = False
        self.__r_measure()
        self.limit_backward = self.cursor
        self.cursor = self.limit
        v_1 = self.limit - self.cursor
        try:
            if not self.__r_Step_1():
                raise lab0()
            B_stemmed = True
        except lab0: pass
        self.cursor = self.limit - v_1
        v_2 = self.limit - self.cursor
        try:
            if not self.__r_Step_2():
                raise lab0()
            B_stemmed = True
        except lab0: pass
        self.cursor = self.limit - v_2
        v_3 = self.limit - self.cursor
        try:
            if not self.__r_Step_3():
                raise lab0()
            B_stemmed = True
        except lab0: pass
        self.cursor = self.limit - v_3
        v_4 = self.limit - self.cursor
        try:
            if not self.__r_Step_4():
                raise lab0()
            B_stemmed = True
        except lab0: pass
        self.cursor = self.limit - v_4
        self.cursor = self.limit_backward
        self.B_GE_removed = False
        v_5 = self.cursor
        try:
            v_6 = self.cursor
            if not self.__r_Lose_prefix():
                raise lab0()
            self.cursor = v_6
            self.__r_measure()
        except lab0: pass
        self.cursor = v_5
        self.limit_backward = self.cursor
        self.cursor = self.limit
        v_7 = self.limit - self.cursor
        try:
            if not self.B_GE_removed:
                raise lab0()
            B_stemmed = True
            if not self.__r_Step_1c():
                raise lab0()
        except lab0: pass
        self.cursor = self.limit - v_7
        self.cursor = self.limit_backward
        self.B_GE_removed = False
        v_8 = self.cursor
        try:
            v_9 = self.cursor
            if not self.__r_Lose_infix():
                raise lab0()
            self.cursor = v_9
            self.__r_measure()
        except lab0: pass
        self.cursor = v_8
        self.limit_backward = self.cursor
        self.cursor = self.limit
        v_10 = self.limit - self.cursor
        try:
            if not self.B_GE_removed:
                raise lab0()
            B_stemmed = True
            if not self.__r_Step_1c():
                raise lab0()
        except lab0: pass
        self.cursor = self.limit - v_10
        self.cursor = self.limit_backward
        self.limit_backward = self.cursor
        self.cursor = self.limit
        v_11 = self.limit - self.cursor
        try:
            if not self.__r_Step_7():
                raise lab0()
            B_stemmed = True
        except lab0: pass
        self.cursor = self.limit - v_11
        v_12 = self.limit - self.cursor
        try:
            if not B_stemmed:
                raise lab0()
            if not self.__r_Step_6():
                raise lab0()
        except lab0: pass
        self.cursor = self.limit - v_12
        self.cursor = self.limit_backward
        return True

    a_0 = [
        Among("a", -1, 1),
        Among("e", -1, 2),
        Among("o", -1, 1),
        Among("u", -1, 1),
        Among("à", -1, 1),
        Among("á", -1, 1),
        Among("â", -1, 1),
        Among("ä", -1, 1),
        Among("è", -1, 2),
        Among("é", -1, 2),
        Among("ê", -1, 2),
        Among("eë", -1, 3),
        Among("ië", -1, 4),
        Among("ò", -1, 1),
        Among("ó", -1, 1),
        Among("ô", -1, 1),
        Among("ö", -1, 1),
        Among("ù", -1, 1),
        Among("ú", -1, 1),
        Among("û", -1, 1),
        Among("ü", -1, 1)
    ]

    a_1 = [
        Among("nde", -1, 8),
        Among("en", -1, 7),
        Among("s", -1, 2),
        Among("'s", 2, 1),
        Among("es", 2, 4),
        Among("ies", 4, 3),
        Among("aus", 2, 6),
        Among("és", 2, 5)
    ]

    a_2 = [
        Among("de", -1, 5),
        Among("ge", -1, 2),
        Among("ische", -1, 4),
        Among("je", -1, 1),
        Among("lijke", -1, 3),
        Among("le", -1, 9),
        Among("ene", -1, 10),
        Among("re", -1, 8),
        Among("se", -1, 7),
        Among("te", -1, 6),
        Among("ieve", -1, 11)
    ]

    a_3 = [
        Among("heid", -1, 3),
        Among("fie", -1, 7),
        Among("gie", -1, 8),
        Among("atie", -1, 1),
        Among("isme", -1, 5),
        Among("ing", -1, 5),
        Among("arij", -1, 6),
        Among("erij", -1, 5),
        Among("sel", -1, 3),
        Among("rder", -1, 4),
        Among("ster", -1, 3),
        Among("iteit", -1, 2),
        Among("dst", -1, 10),
        Among("tst", -1, 9)
    ]

    a_4 = [
        Among("end", -1, 9),
        Among("atief", -1, 2),
        Among("erig", -1, 9),
        Among("achtig", -1, 3),
        Among("ioneel", -1, 1),
        Among("baar", -1, 3),
        Among("laar", -1, 5),
        Among("naar", -1, 4),
        Among("raar", -1, 6),
        Among("eriger", -1, 9),
        Among("achtiger", -1, 3),
        Among("lijker", -1, 8),
        Among("tant", -1, 7),
        Among("erigst", -1, 9),
        Among("achtigst", -1, 3),
        Among("lijkst", -1, 8)
    ]

    a_5 = [
        Among("ig", -1, 1),
        Among("iger", -1, 1),
        Among("igst", -1, 1)
    ]

    a_6 = [
        Among("ft", -1, 2),
        Among("kt", -1, 1),
        Among("pt", -1, 3)
    ]
    as_6 = ("k", "f", "p")

    a_7 = [
        Among("bb", -1, 1),
        Among("cc", -1, 2),
        Among("dd", -1, 3),
        Among("ff", -1, 4),
        Among("gg", -1, 5),
        Among("hh", -1, 6),
        Among("jj", -1, 7),
        Among("kk", -1, 8),
        Among("ll", -1, 9),
        Among("mm", -1, 10),
        Among("nn", -1, 11),
        Among("pp", -1, 12),
        Among("qq", -1, 13),
        Among("rr", -1, 14),
        Among("ss", -1, 15),
        Among("tt", -1, 16),
        Among("v", -1, 4),
        Among("vv", 16, 17),
        Among("ww", -1, 18),
        Among("xx", -1, 19),
        Among("z", -1, 15),
        Among("zz", 20, 20)
    ]

    a_8 = [
        Among("d", -1, 1),
        Among("t", -1, 2)
    ]

    a_9 = [
        Among("", -1, -1),
        Among("eft", 0, 1),
        Among("vaa", 0, 1),
        Among("val", 0, 1),
        Among("vali", 3, -1),
        Among("vare", 0, 1)
    ]

    a_10 = [
        Among("ë", -1, 1),
        Among("ï", -1, 2)
    ]
    as_10 = ("e", "i")

    a_11 = [
        Among("ë", -1, 1),
        Among("ï", -1, 2)
    ]
    as_11 = ("e", "i")


class lab0(BaseException): pass


class lab1(BaseException): pass


class lab2(BaseException): pass
