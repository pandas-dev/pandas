# Generated from polish.sbl by Snowball 3.1.1 - https://snowballstem.org/

from .basestemmer import BaseStemmer
from .among import Among


class PolishStemmer(BaseStemmer):
    '''
    This class implements the stemming algorithm defined by a snowball script.
    Generated from polish.sbl by Snowball 3.1.1 - https://snowballstem.org/
    '''

    g_v = {"a", "e", "i", "o", "u", "y", "ó", "ą", "ę"}

    I_p1 = 0

    def __r_mark_regions(self):
        self.I_p1 = self.limit
        if not self.go_out_grouping(PolishStemmer.g_v):
            return False
        self.cursor += 1
        if not self.go_in_grouping(PolishStemmer.g_v):
            return False
        self.cursor += 1
        self.I_p1 = self.cursor
        return True

    def __r_R1(self):
        return self.I_p1 <= self.cursor

    def __r_remove_endings(self):
        v_1 = self.limit - self.cursor
        try:
            if self.cursor < self.I_p1:
                raise lab0()
            v_3 = self.limit_backward
            self.limit_backward = self.I_p1
            self.ket = self.cursor
            if self.find_among_b(PolishStemmer.a_0) == 0:
                self.limit_backward = v_3
                raise lab0()
            self.bra = self.cursor
            self.limit_backward = v_3
            self.slice_del()
        except lab0: pass
        self.cursor = self.limit - v_1
        self.ket = self.cursor
        among_var = self.find_among_b(PolishStemmer.a_2)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var == 1:
            self.slice_del()
        elif among_var == 2:
            self.slice_from("s")
        elif among_var == 3:
            while True:
                v_4 = self.limit - self.cursor
                try:
                    if not self.__r_R1():
                        raise lab0()
                    self.slice_del()
                    break
                except lab0: pass
                self.cursor = self.limit - v_4
                self.slice_from("s")
                break
        elif among_var == 4:
            self.slice_from("ł")
        else:
            self.slice_del()
            v_5 = self.limit - self.cursor
            try:
                self.ket = self.cursor
                among_var = self.find_among_b(PolishStemmer.a_1)
                if among_var == 0:
                    self.cursor = self.limit - v_5
                    raise lab0()
                self.bra = self.cursor
                self.slice_from(PolishStemmer.as_1[among_var - 1])
            except lab0: pass
        v_6 = self.limit - self.cursor
        try:
            self.ket = self.cursor
            if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "'":
                self.cursor = self.limit - v_6
                raise lab0()
            self.cursor -= 1
            self.bra = self.cursor
            self.slice_del()
        except lab0: pass
        return True

    def __r_normalize_consonant(self):
        self.ket = self.cursor
        among_var = self.find_among_b(PolishStemmer.a_3)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if self.cursor <= self.limit_backward:
            return False
        self.slice_from(PolishStemmer.as_3[among_var - 1])
        return True

    def _stem(self):
        v_1 = self.cursor
        self.__r_mark_regions()
        self.cursor = v_1
        while True:
            v_2 = self.cursor
            try:
                if self.cursor + 2 > self.limit:
                    raise lab0()
                self.cursor += 2
                self.limit_backward = self.cursor
                self.cursor = self.limit
                if not self.__r_remove_endings():
                    raise lab0()
                self.cursor = self.limit_backward
                break
            except lab0: pass
            self.cursor = v_2
            self.limit_backward = self.cursor
            self.cursor = self.limit
            if not self.__r_normalize_consonant():
                return False
            self.cursor = self.limit_backward
            break
        return True

    a_0 = [
        Among("byście", -1, 1),
        Among("bym", -1, 1),
        Among("by", -1, 1),
        Among("byśmy", -1, 1),
        Among("byś", -1, 1)
    ]

    a_1 = [
        Among("ąc", -1, 1),
        Among("ając", 0, 1),
        Among("sząc", 0, 2),
        Among("sz", -1, 1),
        Among("iejsz", 3, 1)
    ]
    as_1 = ("", "s")

    a_2 = [
        Among("a", -1, 1, __r_R1),
        Among("ąca", 0, 1),
        Among("ająca", 1, 1),
        Among("sząca", 1, 2),
        Among("ia", 0, 1, __r_R1),
        Among("sza", 0, 1),
        Among("iejsza", 5, 1),
        Among("ała", 0, 1),
        Among("iała", 7, 1),
        Among("iła", 0, 1),
        Among("ąc", -1, 1),
        Among("ając", 10, 1),
        Among("e", -1, 1, __r_R1),
        Among("ące", 12, 1),
        Among("ające", 13, 1),
        Among("szące", 13, 2),
        Among("ie", 12, 1, __r_R1),
        Among("cie", 16, 1),
        Among("acie", 17, 1),
        Among("ecie", 17, 1),
        Among("icie", 17, 1),
        Among("ajcie", 17, 1),
        Among("liście", 17, 4),
        Among("aliście", 22, 1),
        Among("ieliście", 22, 1),
        Among("iliście", 22, 1),
        Among("łyście", 17, 4),
        Among("ałyście", 26, 1),
        Among("iałyście", 27, 1),
        Among("iłyście", 26, 1),
        Among("sze", 12, 1),
        Among("iejsze", 30, 1),
        Among("ach", -1, 1, __r_R1),
        Among("iach", 32, 1, __r_R1),
        Among("ich", -1, 5),
        Among("ych", -1, 5),
        Among("i", -1, 1, __r_R1),
        Among("ali", 36, 1),
        Among("ieli", 36, 1),
        Among("ili", 36, 1),
        Among("ami", 36, 1, __r_R1),
        Among("iami", 40, 1, __r_R1),
        Among("imi", 36, 5),
        Among("ymi", 36, 5),
        Among("owi", 36, 1, __r_R1),
        Among("iowi", 44, 1, __r_R1),
        Among("aj", -1, 1),
        Among("ej", -1, 5),
        Among("iej", 47, 5),
        Among("am", -1, 1),
        Among("ałam", 49, 1),
        Among("iałam", 50, 1),
        Among("iłam", 49, 1),
        Among("em", -1, 1, __r_R1),
        Among("iem", 53, 1, __r_R1),
        Among("ałem", 53, 1),
        Among("iałem", 55, 1),
        Among("iłem", 53, 1),
        Among("im", -1, 5),
        Among("om", -1, 1, __r_R1),
        Among("iom", 59, 1, __r_R1),
        Among("ym", -1, 5),
        Among("o", -1, 1, __r_R1),
        Among("ego", 62, 5),
        Among("iego", 63, 5),
        Among("ało", 62, 1),
        Among("iało", 65, 1),
        Among("iło", 62, 1),
        Among("u", -1, 1, __r_R1),
        Among("iu", 68, 1, __r_R1),
        Among("emu", 68, 5),
        Among("iemu", 70, 5),
        Among("ów", -1, 1, __r_R1),
        Among("y", -1, 5),
        Among("amy", 73, 1),
        Among("emy", 73, 1),
        Among("imy", 73, 1),
        Among("liśmy", 73, 4),
        Among("aliśmy", 77, 1),
        Among("ieliśmy", 77, 1),
        Among("iliśmy", 77, 1),
        Among("łyśmy", 73, 4),
        Among("ałyśmy", 81, 1),
        Among("iałyśmy", 82, 1),
        Among("iłyśmy", 81, 1),
        Among("ały", 73, 1),
        Among("iały", 85, 1),
        Among("iły", 73, 1),
        Among("asz", -1, 1),
        Among("esz", -1, 1),
        Among("isz", -1, 1),
        Among("ą", -1, 1, __r_R1),
        Among("ącą", 91, 1),
        Among("ającą", 92, 1),
        Among("szącą", 92, 2),
        Among("ią", 91, 1, __r_R1),
        Among("ają", 91, 1),
        Among("szą", 91, 3),
        Among("iejszą", 97, 1),
        Among("ać", -1, 1),
        Among("ieć", -1, 1),
        Among("ić", -1, 1),
        Among("ąć", -1, 1),
        Among("aść", -1, 1),
        Among("eść", -1, 1),
        Among("ę", -1, 1),
        Among("szę", 105, 2),
        Among("ał", -1, 1),
        Among("iał", 107, 1),
        Among("ił", -1, 1),
        Among("łaś", -1, 4),
        Among("ałaś", 110, 1),
        Among("iałaś", 111, 1),
        Among("iłaś", 110, 1),
        Among("łeś", -1, 4),
        Among("ałeś", 114, 1),
        Among("iałeś", 115, 1),
        Among("iłeś", 114, 1)
    ]

    a_3 = [
        Among("ć", -1, 1),
        Among("ń", -1, 2),
        Among("ś", -1, 3),
        Among("ź", -1, 4)
    ]
    as_3 = ("c", "n", "s", "z")


class lab0(BaseException): pass
