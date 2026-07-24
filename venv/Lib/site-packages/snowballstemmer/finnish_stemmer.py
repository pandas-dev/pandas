# Generated from finnish.sbl by Snowball 3.1.1 - https://snowballstem.org/

from .basestemmer import BaseStemmer
from .among import Among


class FinnishStemmer(BaseStemmer):
    '''
    This class implements the stemming algorithm defined by a snowball script.
    Generated from finnish.sbl by Snowball 3.1.1 - https://snowballstem.org/
    '''

    g_AEI = {"a", "e", "i", "ä"}

    g_C = {"b", "c", "d", "f", "g", "h", "j", "k", "l", "m", "n", "p", "q", "r", "s", "t", "v", "w", "x", "z"}

    g_v = {"a", "e", "i", "o", "u", "y", "ä", "ö"}

    g_particle_end = {"a", "e", "i", "n", "o", "t", "u", "y", "ä", "ö"}

    B_ending_removed = False
    I_p2 = 0
    I_p1 = 0

    def __r_mark_regions(self):
        self.I_p1 = self.limit
        self.I_p2 = self.limit
        if not self.go_out_grouping(FinnishStemmer.g_v):
            return False
        self.cursor += 1
        if not self.go_in_grouping(FinnishStemmer.g_v):
            return False
        self.cursor += 1
        self.I_p1 = self.cursor
        if not self.go_out_grouping(FinnishStemmer.g_v):
            return False
        self.cursor += 1
        if not self.go_in_grouping(FinnishStemmer.g_v):
            return False
        self.cursor += 1
        self.I_p2 = self.cursor
        return True

    def __r_particle_etc(self):
        if self.cursor < self.I_p1:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_p1
        self.ket = self.cursor
        among_var = self.find_among_b(FinnishStemmer.a_0)
        if among_var == 0:
            self.limit_backward = v_2
            return False
        self.bra = self.cursor
        self.limit_backward = v_2
        if among_var == 1:
            if not self.in_grouping_b(FinnishStemmer.g_particle_end):
                return False
        else:
            if self.I_p2 > self.cursor:
                return False
        self.slice_del()
        return True

    def __r_possessive(self):
        if self.cursor < self.I_p1:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_p1
        self.ket = self.cursor
        among_var = self.find_among_b(FinnishStemmer.a_4)
        if among_var == 0:
            self.limit_backward = v_2
            return False
        self.bra = self.cursor
        self.limit_backward = v_2
        if among_var == 1:
            try:
                if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "k":
                    raise lab0()
                self.cursor -= 1
                return False
            except lab0: pass
            self.slice_del()
        elif among_var == 2:
            self.slice_del()
            self.ket = self.cursor
            if not self.eq_s_b("kse"):
                return False
            self.bra = self.cursor
            self.slice_from("ksi")
        elif among_var == 3:
            self.slice_del()
        elif among_var == 4:
            if self.find_among_b(FinnishStemmer.a_1) == 0:
                return False
            self.slice_del()
        elif among_var == 5:
            if self.find_among_b(FinnishStemmer.a_2) == 0:
                return False
            self.slice_del()
        else:
            if self.find_among_b(FinnishStemmer.a_3) == 0:
                return False
            self.slice_del()
        return True

    def __r_LV(self):
        return self.find_among_b(FinnishStemmer.a_5) != 0

    def __r_VI(self):
        return self.find_among_b(FinnishStemmer.a_6) != 0

    def __r_A(self):
        while True:
            try:
                if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "a":
                    raise lab0()
                self.cursor -= 1
                break
            except lab0: pass
            if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "'":
                return False
            self.cursor -= 1
            break
        return True

    def __r_E(self):
        while True:
            try:
                if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "e":
                    raise lab0()
                self.cursor -= 1
                break
            except lab0: pass
            if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "'":
                return False
            self.cursor -= 1
            break
        return True

    def __r_I(self):
        while True:
            try:
                if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "i":
                    raise lab0()
                self.cursor -= 1
                break
            except lab0: pass
            if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "'":
                return False
            self.cursor -= 1
            break
        return True

    def __r_O(self):
        while True:
            try:
                if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "o":
                    raise lab0()
                self.cursor -= 1
                break
            except lab0: pass
            if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "'":
                return False
            self.cursor -= 1
            break
        return True

    def __r_U(self):
        while True:
            try:
                if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "u":
                    raise lab0()
                self.cursor -= 1
                break
            except lab0: pass
            if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "'":
                return False
            self.cursor -= 1
            break
        return True

    def __r_A_(self):
        while True:
            try:
                if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "ä":
                    raise lab0()
                self.cursor -= 1
                break
            except lab0: pass
            if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "'":
                return False
            self.cursor -= 1
            break
        return True

    def __r_O_(self):
        while True:
            try:
                if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "ö":
                    raise lab0()
                self.cursor -= 1
                break
            except lab0: pass
            try:
                if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "ø":
                    raise lab0()
                self.cursor -= 1
                break
            except lab0: pass
            if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "'":
                return False
            self.cursor -= 1
            break
        return True

    def __r_case_ending(self):
        if self.cursor < self.I_p1:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_p1
        self.ket = self.cursor
        among_var = self.find_among_b(FinnishStemmer.a_7)
        if among_var == 0:
            self.limit_backward = v_2
            return False
        self.bra = self.cursor
        self.limit_backward = v_2
        if among_var == 1:
            v_3 = self.limit - self.cursor
            try:
                v_4 = self.limit - self.cursor
                while True:
                    v_5 = self.limit - self.cursor
                    try:
                        if not self.__r_LV():
                            raise lab1()
                        break
                    except lab1: pass
                    self.cursor = self.limit - v_5
                    if not self.eq_s_b("ie"):
                        self.cursor = self.limit - v_3
                        raise lab0()
                    break
                self.cursor = self.limit - v_4
                if self.cursor <= self.limit_backward:
                    self.cursor = self.limit - v_3
                    raise lab0()
                self.cursor -= 1
                self.bra = self.cursor
            except lab0: pass
        elif among_var == 2:
            if not self.in_grouping_b(FinnishStemmer.g_v):
                return False
            if not self.in_grouping_b(FinnishStemmer.g_C):
                return False
        elif among_var == 3:
            if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "e":
                return False
            self.cursor -= 1
        self.slice_del()
        self.B_ending_removed = True
        return True

    def __r_other_endings(self):
        if self.cursor < self.I_p2:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_p2
        self.ket = self.cursor
        among_var = self.find_among_b(FinnishStemmer.a_8)
        if among_var == 0:
            self.limit_backward = v_2
            return False
        self.bra = self.cursor
        self.limit_backward = v_2
        if among_var == 1:
            try:
                if not self.eq_s_b("po"):
                    raise lab0()
                return False
            except lab0: pass
        self.slice_del()
        return True

    def __r_i_plural(self):
        if self.cursor < self.I_p1:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_p1
        self.ket = self.cursor
        if self.find_among_b(FinnishStemmer.a_9) == 0:
            self.limit_backward = v_2
            return False
        self.bra = self.cursor
        self.limit_backward = v_2
        self.slice_del()
        return True

    def __r_t_plural(self):
        if self.cursor < self.I_p1:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_p1
        self.ket = self.cursor
        if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "t":
            self.limit_backward = v_2
            return False
        self.cursor -= 1
        self.bra = self.cursor
        v_3 = self.limit - self.cursor
        if not self.in_grouping_b(FinnishStemmer.g_v):
            self.limit_backward = v_2
            return False
        self.cursor = self.limit - v_3
        self.slice_del()
        self.limit_backward = v_2
        if self.cursor < self.I_p2:
            return False
        v_5 = self.limit_backward
        self.limit_backward = self.I_p2
        self.ket = self.cursor
        among_var = self.find_among_b(FinnishStemmer.a_10)
        if among_var == 0:
            self.limit_backward = v_5
            return False
        self.bra = self.cursor
        self.limit_backward = v_5
        if among_var == 1:
            try:
                if not self.eq_s_b("po"):
                    raise lab0()
                return False
            except lab0: pass
        self.slice_del()
        return True

    def __r_tidy(self):
        if self.cursor < self.I_p1:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_p1
        v_3 = self.limit - self.cursor
        try:
            v_4 = self.limit - self.cursor
            if not self.__r_LV():
                raise lab0()
            self.cursor = self.limit - v_4
            self.ket = self.cursor
            if self.cursor <= self.limit_backward:
                raise lab0()
            self.cursor -= 1
            self.bra = self.cursor
            self.slice_del()
        except lab0: pass
        self.cursor = self.limit - v_3
        v_5 = self.limit - self.cursor
        try:
            self.ket = self.cursor
            if not self.in_grouping_b(FinnishStemmer.g_AEI):
                raise lab0()
            self.bra = self.cursor
            if not self.in_grouping_b(FinnishStemmer.g_C):
                raise lab0()
            self.slice_del()
        except lab0: pass
        self.cursor = self.limit - v_5
        v_6 = self.limit - self.cursor
        try:
            self.ket = self.cursor
            if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "j":
                raise lab0()
            self.cursor -= 1
            self.bra = self.cursor
            while True:
                try:
                    if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "o":
                        raise lab1()
                    self.cursor -= 1
                    break
                except lab1: pass
                if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "u":
                    raise lab0()
                self.cursor -= 1
                break
            self.slice_del()
        except lab0: pass
        self.cursor = self.limit - v_6
        v_7 = self.limit - self.cursor
        try:
            self.ket = self.cursor
            if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "o":
                raise lab0()
            self.cursor -= 1
            self.bra = self.cursor
            if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "j":
                raise lab0()
            self.cursor -= 1
            self.slice_del()
        except lab0: pass
        self.cursor = self.limit - v_7
        self.limit_backward = v_2
        v_8 = self.limit - self.cursor
        try:
            if not self.go_in_grouping_b(FinnishStemmer.g_v):
                raise lab0()
            self.ket = self.cursor
            if not self.in_grouping_b(FinnishStemmer.g_C):
                raise lab0()
            self.bra = self.cursor
            S_x = self.slice_to()
            if not self.eq_s_b(S_x):
                raise lab0()
            self.slice_del()
        except lab0: pass
        self.cursor = self.limit - v_8
        self.ket = self.cursor
        if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "'":
            return False
        self.cursor -= 1
        self.bra = self.cursor
        self.slice_del()
        return True

    def _stem(self):
        v_1 = self.cursor
        self.__r_mark_regions()
        self.cursor = v_1
        self.B_ending_removed = False
        self.limit_backward = self.cursor
        self.cursor = self.limit
        v_2 = self.limit - self.cursor
        self.__r_particle_etc()
        self.cursor = self.limit - v_2
        v_3 = self.limit - self.cursor
        self.__r_possessive()
        self.cursor = self.limit - v_3
        v_4 = self.limit - self.cursor
        self.__r_case_ending()
        self.cursor = self.limit - v_4
        v_5 = self.limit - self.cursor
        self.__r_other_endings()
        self.cursor = self.limit - v_5
        while True:
            try:
                if not self.B_ending_removed:
                    raise lab0()
                v_6 = self.limit - self.cursor
                self.__r_i_plural()
                self.cursor = self.limit - v_6
                break
            except lab0: pass
            v_7 = self.limit - self.cursor
            self.__r_t_plural()
            self.cursor = self.limit - v_7
            break
        v_8 = self.limit - self.cursor
        self.__r_tidy()
        self.cursor = self.limit - v_8
        self.cursor = self.limit_backward
        return True

    a_0 = [
        Among("pa", -1, 1),
        Among("sti", -1, 2),
        Among("kaan", -1, 1),
        Among("han", -1, 1),
        Among("kin", -1, 1),
        Among("hän", -1, 1),
        Among("kään", -1, 1),
        Among("ko", -1, 1),
        Among("pä", -1, 1),
        Among("kö", -1, 1)
    ]

    a_1 = [
        Among("lla", -1, -1),
        Among("na", -1, -1),
        Among("ssa", -1, -1),
        Among("ta", -1, -1),
        Among("lta", 3, -1),
        Among("sta", 3, -1)
    ]

    a_2 = [
        Among("llä", -1, -1),
        Among("nä", -1, -1),
        Among("ssä", -1, -1),
        Among("tä", -1, -1),
        Among("ltä", 3, -1),
        Among("stä", 3, -1)
    ]

    a_3 = [
        Among("lle", -1, -1),
        Among("ine", -1, -1)
    ]

    a_4 = [
        Among("nsa", -1, 3),
        Among("mme", -1, 3),
        Among("nne", -1, 3),
        Among("ni", -1, 2),
        Among("si", -1, 1),
        Among("an", -1, 4),
        Among("en", -1, 6),
        Among("än", -1, 5),
        Among("nsä", -1, 3)
    ]

    a_5 = [
        Among("aa", -1, -1),
        Among("ee", -1, -1),
        Among("ii", -1, -1),
        Among("oo", -1, -1),
        Among("uu", -1, -1),
        Among("ää", -1, -1),
        Among("öö", -1, -1)
    ]

    a_6 = [
        Among("'", -1, -1),
        Among("ai", -1, -1),
        Among("ei", -1, -1),
        Among("ii", -1, -1),
        Among("oi", -1, -1),
        Among("ui", -1, -1),
        Among("äi", -1, -1),
        Among("öi", -1, -1)
    ]

    a_7 = [
        Among("a", -1, 2),
        Among("lla", 0, -1),
        Among("na", 0, -1),
        Among("ssa", 0, -1),
        Among("ta", 0, -1),
        Among("lta", 4, -1),
        Among("sta", 4, -1),
        Among("tta", 4, 3),
        Among("lle", -1, -1),
        Among("ine", -1, -1),
        Among("ksi", -1, -1),
        Among("n", -1, 1),
        Among("han", 11, -1, __r_A),
        Among("den", 11, -1, __r_VI),
        Among("seen", 11, -1, __r_LV),
        Among("hen", 11, -1, __r_E),
        Among("tten", 11, -1, __r_VI),
        Among("hin", 11, -1, __r_I),
        Among("siin", 11, -1, __r_VI),
        Among("hon", 11, -1, __r_O),
        Among("hun", 11, -1, __r_U),
        Among("hän", 11, -1, __r_A_),
        Among("hön", 11, -1, __r_O_),
        Among("ä", -1, 2),
        Among("llä", 23, -1),
        Among("nä", 23, -1),
        Among("ssä", 23, -1),
        Among("tä", 23, -1),
        Among("ltä", 27, -1),
        Among("stä", 27, -1),
        Among("ttä", 27, 3)
    ]

    a_8 = [
        Among("eja", -1, -1),
        Among("mma", -1, 1),
        Among("imma", 1, -1),
        Among("mpa", -1, 1),
        Among("impa", 3, -1),
        Among("mmi", -1, 1),
        Among("immi", 5, -1),
        Among("mpi", -1, 1),
        Among("impi", 7, -1),
        Among("ejä", -1, -1),
        Among("mmä", -1, 1),
        Among("immä", 10, -1),
        Among("mpä", -1, 1),
        Among("impä", 12, -1)
    ]

    a_9 = [
        Among("i", -1, -1),
        Among("j", -1, -1)
    ]

    a_10 = [
        Among("mma", -1, 1),
        Among("imma", 0, -1)
    ]


class lab0(BaseException): pass


class lab1(BaseException): pass
