# Generated from russian.sbl by Snowball 3.1.1 - https://snowballstem.org/

from .basestemmer import BaseStemmer
from .among import Among


class RussianStemmer(BaseStemmer):
    '''
    This class implements the stemming algorithm defined by a snowball script.
    Generated from russian.sbl by Snowball 3.1.1 - https://snowballstem.org/
    '''

    g_v = {"а", "е", "и", "о", "у", "ы", "э", "ю", "я"}

    I_p2 = 0
    I_pV = 0

    def __r_mark_regions(self):
        self.I_pV = self.limit
        self.I_p2 = self.limit
        v_1 = self.cursor
        try:
            if not self.go_out_grouping(RussianStemmer.g_v):
                raise lab0()
            self.cursor += 1
            self.I_pV = self.cursor
            if not self.go_in_grouping(RussianStemmer.g_v):
                raise lab0()
            self.cursor += 1
            if not self.go_out_grouping(RussianStemmer.g_v):
                raise lab0()
            self.cursor += 1
            if not self.go_in_grouping(RussianStemmer.g_v):
                raise lab0()
            self.cursor += 1
            self.I_p2 = self.cursor
        except lab0: pass
        self.cursor = v_1
        return True

    def __r_perfective_gerund(self):
        self.ket = self.cursor
        among_var = self.find_among_b(RussianStemmer.a_0)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var == 1:
            while True:
                try:
                    if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "а":
                        raise lab0()
                    self.cursor -= 1
                    break
                except lab0: pass
                if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "я":
                    return False
                self.cursor -= 1
                break
            self.slice_del()
        else:
            self.slice_del()
        return True

    def __r_adjective(self):
        self.ket = self.cursor
        if self.find_among_b(RussianStemmer.a_1) == 0:
            return False
        self.bra = self.cursor
        self.slice_del()
        return True

    def __r_adjectival(self):
        if not self.__r_adjective():
            return False
        v_1 = self.limit - self.cursor
        try:
            self.ket = self.cursor
            among_var = self.find_among_b(RussianStemmer.a_2)
            if among_var == 0:
                self.cursor = self.limit - v_1
                raise lab0()
            self.bra = self.cursor
            if among_var == 1:
                while True:
                    try:
                        if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "а":
                            raise lab1()
                        self.cursor -= 1
                        break
                    except lab1: pass
                    if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "я":
                        self.cursor = self.limit - v_1
                        raise lab0()
                    self.cursor -= 1
                    break
                self.slice_del()
            else:
                self.slice_del()
        except lab0: pass
        return True

    def __r_reflexive(self):
        self.ket = self.cursor
        if self.find_among_b(RussianStemmer.a_3) == 0:
            return False
        self.bra = self.cursor
        self.slice_del()
        return True

    def __r_verb(self):
        self.ket = self.cursor
        among_var = self.find_among_b(RussianStemmer.a_4)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var == 1:
            while True:
                try:
                    if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "а":
                        raise lab0()
                    self.cursor -= 1
                    break
                except lab0: pass
                if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "я":
                    return False
                self.cursor -= 1
                break
            self.slice_del()
        else:
            self.slice_del()
        return True

    def __r_noun(self):
        self.ket = self.cursor
        if self.find_among_b(RussianStemmer.a_5) == 0:
            return False
        self.bra = self.cursor
        self.slice_del()
        return True

    def __r_derivational(self):
        self.ket = self.cursor
        if self.find_among_b(RussianStemmer.a_6) == 0:
            return False
        self.bra = self.cursor
        if self.I_p2 > self.cursor:
            return False
        self.slice_del()
        return True

    def __r_tidy_up(self):
        self.ket = self.cursor
        among_var = self.find_among_b(RussianStemmer.a_7)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var == 1:
            self.slice_del()
            self.ket = self.cursor
            if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "н":
                return False
            self.cursor -= 1
            self.bra = self.cursor
            if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "н":
                return False
            self.cursor -= 1
            self.slice_del()
        elif among_var == 2:
            if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "н":
                return False
            self.cursor -= 1
            self.slice_del()
        else:
            self.slice_del()
        return True

    def _stem(self):
        v_1 = self.cursor
        try:
            while True:
                v_2 = self.cursor
                try:
                    while True:
                        v_3 = self.cursor
                        try:
                            self.bra = self.cursor
                            if self.cursor == self.limit or self.current[self.cursor] != "ё":
                                raise lab2()
                            self.cursor += 1
                            self.ket = self.cursor
                            self.cursor = v_3
                            break
                        except lab2: pass
                        self.cursor = v_3
                        if self.cursor >= self.limit:
                            raise lab1()
                        self.cursor += 1
                    self.slice_from("е")
                    continue
                except lab1: pass
                self.cursor = v_2
                break
        except lab0: pass
        self.cursor = v_1
        self.__r_mark_regions()
        self.limit_backward = self.cursor
        self.cursor = self.limit
        if self.cursor < self.I_pV:
            return False
        v_5 = self.limit_backward
        self.limit_backward = self.I_pV
        v_6 = self.limit - self.cursor
        try:
            while True:
                v_7 = self.limit - self.cursor
                try:
                    if not self.__r_perfective_gerund():
                        raise lab1()
                    break
                except lab1: pass
                self.cursor = self.limit - v_7
                v_8 = self.limit - self.cursor
                try:
                    if not self.__r_reflexive():
                        self.cursor = self.limit - v_8
                        raise lab1()
                except lab1: pass
                while True:
                    v_9 = self.limit - self.cursor
                    try:
                        if not self.__r_adjectival():
                            raise lab1()
                        break
                    except lab1: pass
                    self.cursor = self.limit - v_9
                    try:
                        if not self.__r_verb():
                            raise lab1()
                        break
                    except lab1: pass
                    self.cursor = self.limit - v_9
                    if not self.__r_noun():
                        raise lab0()
                    break
                break
        except lab0: pass
        self.cursor = self.limit - v_6
        v_10 = self.limit - self.cursor
        try:
            self.ket = self.cursor
            if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "и":
                self.cursor = self.limit - v_10
                raise lab0()
            self.cursor -= 1
            self.bra = self.cursor
            self.slice_del()
        except lab0: pass
        v_11 = self.limit - self.cursor
        self.__r_derivational()
        self.cursor = self.limit - v_11
        v_12 = self.limit - self.cursor
        self.__r_tidy_up()
        self.cursor = self.limit - v_12
        self.limit_backward = v_5
        self.cursor = self.limit_backward
        return True

    a_0 = [
        Among("в", -1, 1),
        Among("ив", 0, 2),
        Among("ыв", 0, 2),
        Among("вши", -1, 1),
        Among("ивши", 3, 2),
        Among("ывши", 3, 2),
        Among("вшись", -1, 1),
        Among("ившись", 6, 2),
        Among("ывшись", 6, 2)
    ]

    a_1 = [
        Among("ее", -1, 1),
        Among("ие", -1, 1),
        Among("ое", -1, 1),
        Among("ые", -1, 1),
        Among("ими", -1, 1),
        Among("ыми", -1, 1),
        Among("ей", -1, 1),
        Among("ий", -1, 1),
        Among("ой", -1, 1),
        Among("ый", -1, 1),
        Among("ем", -1, 1),
        Among("им", -1, 1),
        Among("ом", -1, 1),
        Among("ым", -1, 1),
        Among("его", -1, 1),
        Among("ого", -1, 1),
        Among("ему", -1, 1),
        Among("ому", -1, 1),
        Among("их", -1, 1),
        Among("ых", -1, 1),
        Among("ею", -1, 1),
        Among("ою", -1, 1),
        Among("ую", -1, 1),
        Among("юю", -1, 1),
        Among("ая", -1, 1),
        Among("яя", -1, 1)
    ]

    a_2 = [
        Among("ем", -1, 1),
        Among("нн", -1, 1),
        Among("вш", -1, 1),
        Among("ивш", 2, 2),
        Among("ывш", 2, 2),
        Among("щ", -1, 1),
        Among("ющ", 5, 1),
        Among("ующ", 6, 2)
    ]

    a_3 = [
        Among("сь", -1, 1),
        Among("ся", -1, 1)
    ]

    a_4 = [
        Among("ла", -1, 1),
        Among("ила", 0, 2),
        Among("ыла", 0, 2),
        Among("на", -1, 1),
        Among("ена", 3, 2),
        Among("ете", -1, 1),
        Among("ите", -1, 2),
        Among("йте", -1, 1),
        Among("ейте", 7, 2),
        Among("уйте", 7, 2),
        Among("ли", -1, 1),
        Among("или", 10, 2),
        Among("ыли", 10, 2),
        Among("й", -1, 1),
        Among("ей", 13, 2),
        Among("уй", 13, 2),
        Among("л", -1, 1),
        Among("ил", 16, 2),
        Among("ыл", 16, 2),
        Among("ем", -1, 1),
        Among("им", -1, 2),
        Among("ым", -1, 2),
        Among("н", -1, 1),
        Among("ен", 22, 2),
        Among("ло", -1, 1),
        Among("ило", 24, 2),
        Among("ыло", 24, 2),
        Among("но", -1, 1),
        Among("ено", 27, 2),
        Among("нно", 27, 1),
        Among("ет", -1, 1),
        Among("ует", 30, 2),
        Among("ит", -1, 2),
        Among("ыт", -1, 2),
        Among("ют", -1, 1),
        Among("уют", 34, 2),
        Among("ят", -1, 2),
        Among("ны", -1, 1),
        Among("ены", 37, 2),
        Among("ть", -1, 1),
        Among("ить", 39, 2),
        Among("ыть", 39, 2),
        Among("ешь", -1, 1),
        Among("ишь", -1, 2),
        Among("ю", -1, 2),
        Among("ую", 44, 2)
    ]

    a_5 = [
        Among("а", -1, 1),
        Among("ев", -1, 1),
        Among("ов", -1, 1),
        Among("е", -1, 1),
        Among("ие", 3, 1),
        Among("ье", 3, 1),
        Among("и", -1, 1),
        Among("еи", 6, 1),
        Among("ии", 6, 1),
        Among("ами", 6, 1),
        Among("ями", 6, 1),
        Among("иями", 10, 1),
        Among("й", -1, 1),
        Among("ей", 12, 1),
        Among("ией", 13, 1),
        Among("ий", 12, 1),
        Among("ой", 12, 1),
        Among("ам", -1, 1),
        Among("ем", -1, 1),
        Among("ием", 18, 1),
        Among("ом", -1, 1),
        Among("ям", -1, 1),
        Among("иям", 21, 1),
        Among("о", -1, 1),
        Among("у", -1, 1),
        Among("ах", -1, 1),
        Among("ях", -1, 1),
        Among("иях", 26, 1),
        Among("ы", -1, 1),
        Among("ь", -1, 1),
        Among("ю", -1, 1),
        Among("ию", 30, 1),
        Among("ью", 30, 1),
        Among("я", -1, 1),
        Among("ия", 33, 1),
        Among("ья", 33, 1)
    ]

    a_6 = [
        Among("ост", -1, 1),
        Among("ость", -1, 1)
    ]

    a_7 = [
        Among("ейше", -1, 1),
        Among("н", -1, 2),
        Among("ейш", -1, 1),
        Among("ь", -1, 3)
    ]


class lab0(BaseException): pass


class lab1(BaseException): pass


class lab2(BaseException): pass
