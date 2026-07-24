# Generated from romanian.sbl by Snowball 3.1.1 - https://snowballstem.org/

from .basestemmer import BaseStemmer
from .among import Among


class RomanianStemmer(BaseStemmer):
    '''
    This class implements the stemming algorithm defined by a snowball script.
    Generated from romanian.sbl by Snowball 3.1.1 - https://snowballstem.org/
    '''

    g_v = {"a", "e", "i", "o", "u", "â", "î", "ă"}

    B_standard_suffix_removed = False
    I_p2 = 0
    I_p1 = 0
    I_pV = 0

    def __r_norm(self):
        v_1 = self.cursor
        try:
            while True:
                v_2 = self.cursor
                try:
                    while True:
                        v_3 = self.cursor
                        try:
                            self.bra = self.cursor
                            among_var = self.find_among(RomanianStemmer.a_0)
                            if among_var == 0:
                                raise lab2()
                            self.ket = self.cursor
                            self.slice_from(RomanianStemmer.as_0[among_var - 1])
                            self.cursor = v_3
                            break
                        except lab2: pass
                        self.cursor = v_3
                        if self.cursor >= self.limit:
                            raise lab1()
                        self.cursor += 1
                    continue
                except lab1: pass
                self.cursor = v_2
                break
        except lab0: pass
        self.cursor = v_1
        return True

    def __r_prelude(self):
        while True:
            v_1 = self.cursor
            try:
                while True:
                    v_2 = self.cursor
                    try:
                        if not self.in_grouping(RomanianStemmer.g_v):
                            raise lab1()
                        self.bra = self.cursor
                        while True:
                            v_3 = self.cursor
                            try:
                                if self.cursor == self.limit or self.current[self.cursor] != "u":
                                    raise lab2()
                                self.cursor += 1
                                self.ket = self.cursor
                                if not self.in_grouping(RomanianStemmer.g_v):
                                    raise lab2()
                                self.slice_from("U")
                                break
                            except lab2: pass
                            self.cursor = v_3
                            if self.cursor == self.limit or self.current[self.cursor] != "i":
                                raise lab1()
                            self.cursor += 1
                            self.ket = self.cursor
                            if not self.in_grouping(RomanianStemmer.g_v):
                                raise lab1()
                            self.slice_from("I")
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
                    if not self.in_grouping(RomanianStemmer.g_v):
                        raise lab1()
                    while True:
                        v_3 = self.cursor
                        try:
                            if not self.out_grouping(RomanianStemmer.g_v):
                                raise lab2()
                            if not self.go_out_grouping(RomanianStemmer.g_v):
                                raise lab2()
                            self.cursor += 1
                            break
                        except lab2: pass
                        self.cursor = v_3
                        if not self.in_grouping(RomanianStemmer.g_v):
                            raise lab1()
                        if not self.go_in_grouping(RomanianStemmer.g_v):
                            raise lab1()
                        self.cursor += 1
                        break
                    break
                except lab1: pass
                self.cursor = v_2
                if not self.out_grouping(RomanianStemmer.g_v):
                    raise lab0()
                while True:
                    v_4 = self.cursor
                    try:
                        if not self.out_grouping(RomanianStemmer.g_v):
                            raise lab1()
                        if not self.go_out_grouping(RomanianStemmer.g_v):
                            raise lab1()
                        self.cursor += 1
                        break
                    except lab1: pass
                    self.cursor = v_4
                    if not self.in_grouping(RomanianStemmer.g_v):
                        raise lab0()
                    if self.cursor >= self.limit:
                        raise lab0()
                    self.cursor += 1
                    break
                break
            self.I_pV = self.cursor
        except lab0: pass
        self.cursor = v_1
        v_5 = self.cursor
        try:
            if not self.go_out_grouping(RomanianStemmer.g_v):
                raise lab0()
            self.cursor += 1
            if not self.go_in_grouping(RomanianStemmer.g_v):
                raise lab0()
            self.cursor += 1
            self.I_p1 = self.cursor
            if not self.go_out_grouping(RomanianStemmer.g_v):
                raise lab0()
            self.cursor += 1
            if not self.go_in_grouping(RomanianStemmer.g_v):
                raise lab0()
            self.cursor += 1
            self.I_p2 = self.cursor
        except lab0: pass
        self.cursor = v_5
        return True

    def __r_postlude(self):
        while True:
            v_1 = self.cursor
            try:
                self.bra = self.cursor
                among_var = self.find_among(RomanianStemmer.a_1)
                self.ket = self.cursor
                if among_var == 1:
                    self.slice_from("i")
                elif among_var == 2:
                    self.slice_from("u")
                else:
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

    def __r_step_0(self):
        self.ket = self.cursor
        among_var = self.find_among_b(RomanianStemmer.a_2)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if not self.__r_R1():
            return False
        if among_var == 1:
            self.slice_del()
        elif among_var == 2:
            self.slice_from("a")
        elif among_var == 3:
            self.slice_from("e")
        elif among_var == 4:
            self.slice_from("i")
        elif among_var == 5:
            try:
                if not self.eq_s_b("ab"):
                    raise lab0()
                return False
            except lab0: pass
            self.slice_from("i")
        elif among_var == 6:
            self.slice_from("at")
        else:
            self.slice_from("ați")
        return True

    def __r_combo_suffix(self):
        v_1 = self.limit - self.cursor
        self.ket = self.cursor
        among_var = self.find_among_b(RomanianStemmer.a_3)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if not self.__r_R1():
            return False
        self.slice_from(RomanianStemmer.as_3[among_var - 1])
        self.B_standard_suffix_removed = True
        self.cursor = self.limit - v_1
        return True

    def __r_standard_suffix(self):
        self.B_standard_suffix_removed = False
        while True:
            v_1 = self.limit - self.cursor
            try:
                if not self.__r_combo_suffix():
                    raise lab0()
                continue
            except lab0: pass
            self.cursor = self.limit - v_1
            break
        self.ket = self.cursor
        among_var = self.find_among_b(RomanianStemmer.a_4)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if self.I_p2 > self.cursor:
            return False
        if among_var == 1:
            self.slice_del()
        elif among_var == 2:
            if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "ț":
                return False
            self.cursor -= 1
            self.bra = self.cursor
            self.slice_from("t")
        else:
            self.slice_from("ist")
        self.B_standard_suffix_removed = True
        return True

    def __r_verb_suffix(self):
        if self.cursor < self.I_pV:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_pV
        self.ket = self.cursor
        among_var = self.find_among_b(RomanianStemmer.a_5)
        if among_var == 0:
            self.limit_backward = v_2
            return False
        self.bra = self.cursor
        if among_var == 1:
            while True:
                try:
                    if not self.out_grouping_b(RomanianStemmer.g_v):
                        raise lab0()
                    break
                except lab0: pass
                if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "u":
                    self.limit_backward = v_2
                    return False
                self.cursor -= 1
                break
            self.slice_del()
        else:
            self.slice_del()
        self.limit_backward = v_2
        return True

    def __r_vowel_suffix(self):
        self.ket = self.cursor
        if self.find_among_b(RomanianStemmer.a_6) == 0:
            return False
        self.bra = self.cursor
        if self.I_pV > self.cursor:
            return False
        self.slice_del()
        return True

    def _stem(self):
        self.__r_norm()
        v_1 = self.cursor
        self.__r_prelude()
        self.cursor = v_1
        self.__r_mark_regions()
        self.limit_backward = self.cursor
        self.cursor = self.limit
        v_2 = self.limit - self.cursor
        self.__r_step_0()
        self.cursor = self.limit - v_2
        v_3 = self.limit - self.cursor
        self.__r_standard_suffix()
        self.cursor = self.limit - v_3
        v_4 = self.limit - self.cursor
        try:
            while True:
                try:
                    if not self.B_standard_suffix_removed:
                        raise lab1()
                    break
                except lab1: pass
                if not self.__r_verb_suffix():
                    raise lab0()
                break
        except lab0: pass
        self.cursor = self.limit - v_4
        v_5 = self.limit - self.cursor
        self.__r_vowel_suffix()
        self.cursor = self.limit - v_5
        self.cursor = self.limit_backward
        v_6 = self.cursor
        self.__r_postlude()
        self.cursor = v_6
        return True

    a_0 = [
        Among("ş", -1, 1),
        Among("ţ", -1, 2)
    ]
    as_0 = ("ș", "ț")

    a_1 = [
        Among("", -1, 3),
        Among("I", 0, 1),
        Among("U", 0, 2)
    ]

    a_2 = [
        Among("ea", -1, 3),
        Among("ația", -1, 7),
        Among("aua", -1, 2),
        Among("iua", -1, 4),
        Among("ație", -1, 7),
        Among("ele", -1, 3),
        Among("ile", -1, 5),
        Among("iile", 6, 4),
        Among("iei", -1, 4),
        Among("atei", -1, 6),
        Among("ii", -1, 4),
        Among("ului", -1, 1),
        Among("ul", -1, 1),
        Among("elor", -1, 3),
        Among("ilor", -1, 4),
        Among("iilor", 14, 4)
    ]

    a_3 = [
        Among("icala", -1, 4),
        Among("iciva", -1, 4),
        Among("ativa", -1, 5),
        Among("itiva", -1, 6),
        Among("icale", -1, 4),
        Among("ațiune", -1, 5),
        Among("ițiune", -1, 6),
        Among("atoare", -1, 5),
        Among("itoare", -1, 6),
        Among("ătoare", -1, 5),
        Among("icitate", -1, 4),
        Among("abilitate", -1, 1),
        Among("ibilitate", -1, 2),
        Among("ivitate", -1, 3),
        Among("icive", -1, 4),
        Among("ative", -1, 5),
        Among("itive", -1, 6),
        Among("icali", -1, 4),
        Among("atori", -1, 5),
        Among("icatori", 18, 4),
        Among("itori", -1, 6),
        Among("ători", -1, 5),
        Among("icitati", -1, 4),
        Among("abilitati", -1, 1),
        Among("ivitati", -1, 3),
        Among("icivi", -1, 4),
        Among("ativi", -1, 5),
        Among("itivi", -1, 6),
        Among("icităi", -1, 4),
        Among("abilităi", -1, 1),
        Among("ivităi", -1, 3),
        Among("icități", -1, 4),
        Among("abilități", -1, 1),
        Among("ivități", -1, 3),
        Among("ical", -1, 4),
        Among("ator", -1, 5),
        Among("icator", 35, 4),
        Among("itor", -1, 6),
        Among("ător", -1, 5),
        Among("iciv", -1, 4),
        Among("ativ", -1, 5),
        Among("itiv", -1, 6),
        Among("icală", -1, 4),
        Among("icivă", -1, 4),
        Among("ativă", -1, 5),
        Among("itivă", -1, 6)
    ]
    as_3 = ("abil", "ibil", "iv", "ic", "at", "it")

    a_4 = [
        Among("ica", -1, 1),
        Among("abila", -1, 1),
        Among("ibila", -1, 1),
        Among("oasa", -1, 1),
        Among("ata", -1, 1),
        Among("ita", -1, 1),
        Among("anta", -1, 1),
        Among("ista", -1, 3),
        Among("uta", -1, 1),
        Among("iva", -1, 1),
        Among("ic", -1, 1),
        Among("ice", -1, 1),
        Among("abile", -1, 1),
        Among("ibile", -1, 1),
        Among("isme", -1, 3),
        Among("iune", -1, 2),
        Among("oase", -1, 1),
        Among("ate", -1, 1),
        Among("itate", 17, 1),
        Among("ite", -1, 1),
        Among("ante", -1, 1),
        Among("iste", -1, 3),
        Among("ute", -1, 1),
        Among("ive", -1, 1),
        Among("ici", -1, 1),
        Among("abili", -1, 1),
        Among("ibili", -1, 1),
        Among("iuni", -1, 2),
        Among("atori", -1, 1),
        Among("osi", -1, 1),
        Among("ati", -1, 1),
        Among("itati", 30, 1),
        Among("iti", -1, 1),
        Among("anti", -1, 1),
        Among("isti", -1, 3),
        Among("uti", -1, 1),
        Among("iști", -1, 3),
        Among("ivi", -1, 1),
        Among("ităi", -1, 1),
        Among("oși", -1, 1),
        Among("ități", -1, 1),
        Among("abil", -1, 1),
        Among("ibil", -1, 1),
        Among("ism", -1, 3),
        Among("ator", -1, 1),
        Among("os", -1, 1),
        Among("at", -1, 1),
        Among("it", -1, 1),
        Among("ant", -1, 1),
        Among("ist", -1, 3),
        Among("ut", -1, 1),
        Among("iv", -1, 1),
        Among("ică", -1, 1),
        Among("abilă", -1, 1),
        Among("ibilă", -1, 1),
        Among("oasă", -1, 1),
        Among("ată", -1, 1),
        Among("ită", -1, 1),
        Among("antă", -1, 1),
        Among("istă", -1, 3),
        Among("ută", -1, 1),
        Among("ivă", -1, 1)
    ]

    a_5 = [
        Among("ea", -1, 1),
        Among("ia", -1, 1),
        Among("esc", -1, 1),
        Among("ăsc", -1, 1),
        Among("ind", -1, 1),
        Among("ând", -1, 1),
        Among("are", -1, 1),
        Among("ere", -1, 1),
        Among("ire", -1, 1),
        Among("âre", -1, 1),
        Among("se", -1, 2),
        Among("ase", 10, 1),
        Among("sese", 10, 2),
        Among("ise", 10, 1),
        Among("use", 10, 1),
        Among("âse", 10, 1),
        Among("ește", -1, 1),
        Among("ăște", -1, 1),
        Among("eze", -1, 1),
        Among("ai", -1, 1),
        Among("eai", 19, 1),
        Among("iai", 19, 1),
        Among("sei", -1, 2),
        Among("ești", -1, 1),
        Among("ăști", -1, 1),
        Among("ui", -1, 1),
        Among("ezi", -1, 1),
        Among("âi", -1, 1),
        Among("ași", -1, 1),
        Among("seși", -1, 2),
        Among("aseși", 29, 1),
        Among("seseși", 29, 2),
        Among("iseși", 29, 1),
        Among("useși", 29, 1),
        Among("âseși", 29, 1),
        Among("iși", -1, 1),
        Among("uși", -1, 1),
        Among("âși", -1, 1),
        Among("ați", -1, 2),
        Among("eați", 38, 1),
        Among("iați", 38, 1),
        Among("eți", -1, 2),
        Among("iți", -1, 2),
        Among("âți", -1, 2),
        Among("arăți", -1, 1),
        Among("serăți", -1, 2),
        Among("aserăți", 45, 1),
        Among("seserăți", 45, 2),
        Among("iserăți", 45, 1),
        Among("userăți", 45, 1),
        Among("âserăți", 45, 1),
        Among("irăți", -1, 1),
        Among("urăți", -1, 1),
        Among("ârăți", -1, 1),
        Among("am", -1, 1),
        Among("eam", 54, 1),
        Among("iam", 54, 1),
        Among("em", -1, 2),
        Among("asem", 57, 1),
        Among("sesem", 57, 2),
        Among("isem", 57, 1),
        Among("usem", 57, 1),
        Among("âsem", 57, 1),
        Among("im", -1, 2),
        Among("âm", -1, 2),
        Among("ăm", -1, 2),
        Among("arăm", 65, 1),
        Among("serăm", 65, 2),
        Among("aserăm", 67, 1),
        Among("seserăm", 67, 2),
        Among("iserăm", 67, 1),
        Among("userăm", 67, 1),
        Among("âserăm", 67, 1),
        Among("irăm", 65, 1),
        Among("urăm", 65, 1),
        Among("ârăm", 65, 1),
        Among("au", -1, 1),
        Among("eau", 76, 1),
        Among("iau", 76, 1),
        Among("indu", -1, 1),
        Among("ându", -1, 1),
        Among("ez", -1, 1),
        Among("ească", -1, 1),
        Among("ară", -1, 1),
        Among("seră", -1, 2),
        Among("aseră", 84, 1),
        Among("seseră", 84, 2),
        Among("iseră", 84, 1),
        Among("useră", 84, 1),
        Among("âseră", 84, 1),
        Among("iră", -1, 1),
        Among("ură", -1, 1),
        Among("âră", -1, 1),
        Among("ează", -1, 1)
    ]

    a_6 = [
        Among("a", -1, 1),
        Among("e", -1, 1),
        Among("ie", 1, 1),
        Among("i", -1, 1),
        Among("ă", -1, 1)
    ]


class lab0(BaseException): pass


class lab1(BaseException): pass


class lab2(BaseException): pass
