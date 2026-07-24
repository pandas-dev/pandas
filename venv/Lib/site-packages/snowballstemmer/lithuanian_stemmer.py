# Generated from lithuanian.sbl by Snowball 3.1.1 - https://snowballstem.org/

from .basestemmer import BaseStemmer
from .among import Among


class LithuanianStemmer(BaseStemmer):
    '''
    This class implements the stemming algorithm defined by a snowball script.
    Generated from lithuanian.sbl by Snowball 3.1.1 - https://snowballstem.org/
    '''

    g_v = {"a", "e", "i", "o", "u", "y", "ą", "ė", "ę", "į", "ū", "ų"}

    I_p1 = 0

    def __r_step1(self):
        if self.cursor < self.I_p1:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_p1
        self.ket = self.cursor
        if self.find_among_b(LithuanianStemmer.a_0) == 0:
            self.limit_backward = v_2
            return False
        self.bra = self.cursor
        self.limit_backward = v_2
        self.slice_del()
        return True

    def __r_step2(self):
        while True:
            v_1 = self.limit - self.cursor
            try:
                if self.cursor < self.I_p1:
                    raise lab0()
                v_3 = self.limit_backward
                self.limit_backward = self.I_p1
                self.ket = self.cursor
                if self.find_among_b(LithuanianStemmer.a_1) == 0:
                    self.limit_backward = v_3
                    raise lab0()
                self.bra = self.cursor
                self.limit_backward = v_3
                self.slice_del()
                continue
            except lab0: pass
            self.cursor = self.limit - v_1
            break
        return True

    def __r_fix_conflicts(self):
        self.ket = self.cursor
        among_var = self.find_among_b(LithuanianStemmer.a_2)
        if among_var == 0:
            return False
        self.bra = self.cursor
        self.slice_from(LithuanianStemmer.as_2[among_var - 1])
        return True

    def __r_fix_chdz(self):
        self.ket = self.cursor
        among_var = self.find_among_b(LithuanianStemmer.a_3)
        if among_var == 0:
            return False
        self.bra = self.cursor
        self.slice_from(LithuanianStemmer.as_3[among_var - 1])
        return True

    def __r_fix_gd(self):
        self.ket = self.cursor
        if not self.eq_s_b("gd"):
            return False
        self.bra = self.cursor
        self.slice_from("g")
        return True

    def _stem(self):
        self.I_p1 = self.limit
        v_1 = self.cursor
        try:
            v_2 = self.cursor
            try:
                if self.cursor == self.limit or self.current[self.cursor] != "a":
                    self.cursor = v_2
                    raise lab1()
                self.cursor += 1
                if len(self.current) <= 6:
                    self.cursor = v_2
                    raise lab1()
            except lab1: pass
            if not self.go_out_grouping(LithuanianStemmer.g_v):
                raise lab0()
            self.cursor += 1
            if not self.go_in_grouping(LithuanianStemmer.g_v):
                raise lab0()
            self.cursor += 1
            self.I_p1 = self.cursor
        except lab0: pass
        self.cursor = v_1
        self.limit_backward = self.cursor
        self.cursor = self.limit
        v_3 = self.limit - self.cursor
        self.__r_fix_conflicts()
        self.cursor = self.limit - v_3
        v_4 = self.limit - self.cursor
        self.__r_step1()
        self.cursor = self.limit - v_4
        v_5 = self.limit - self.cursor
        self.__r_fix_chdz()
        self.cursor = self.limit - v_5
        v_6 = self.limit - self.cursor
        self.__r_step2()
        self.cursor = self.limit - v_6
        v_7 = self.limit - self.cursor
        self.__r_fix_chdz()
        self.cursor = self.limit - v_7
        v_8 = self.limit - self.cursor
        self.__r_fix_gd()
        self.cursor = self.limit - v_8
        self.ket = self.cursor
        if self.cursor <= self.limit_backward or self.current[self.cursor - 1] != "'":
            return False
        self.cursor -= 1
        self.bra = self.cursor
        self.slice_del()
        self.cursor = self.limit_backward
        return True

    a_0 = [
        Among("a", -1, -1),
        Among("ia", 0, -1),
        Among("osna", 0, -1),
        Among("iosna", 2, -1),
        Among("uosna", 2, -1),
        Among("iuosna", 4, -1),
        Among("ysna", 0, -1),
        Among("ėsna", 0, -1),
        Among("e", -1, -1),
        Among("ie", 8, -1),
        Among("enie", 9, -1),
        Among("oje", 8, -1),
        Among("ioje", 11, -1),
        Among("uje", 8, -1),
        Among("iuje", 13, -1),
        Among("yje", 8, -1),
        Among("enyje", 15, -1),
        Among("ėje", 8, -1),
        Among("ame", 8, -1),
        Among("iame", 18, -1),
        Among("sime", 8, -1),
        Among("ome", 8, -1),
        Among("ėme", 8, -1),
        Among("tumėme", 22, -1),
        Among("ose", 8, -1),
        Among("iose", 24, -1),
        Among("uose", 24, -1),
        Among("iuose", 26, -1),
        Among("yse", 8, -1),
        Among("enyse", 28, -1),
        Among("ėse", 8, -1),
        Among("ate", 8, -1),
        Among("iate", 31, -1),
        Among("ite", 8, -1),
        Among("kite", 33, -1),
        Among("site", 33, -1),
        Among("ote", 8, -1),
        Among("tute", 8, -1),
        Among("ėte", 8, -1),
        Among("tumėte", 38, -1),
        Among("i", -1, -1),
        Among("ai", 40, -1),
        Among("iai", 41, -1),
        Among("ei", 40, -1),
        Among("tumei", 43, -1),
        Among("ki", 40, -1),
        Among("imi", 40, -1),
        Among("umi", 40, -1),
        Among("iumi", 47, -1),
        Among("si", 40, -1),
        Among("asi", 49, -1),
        Among("iasi", 50, -1),
        Among("esi", 49, -1),
        Among("iesi", 52, -1),
        Among("siesi", 53, -1),
        Among("isi", 49, -1),
        Among("aisi", 55, -1),
        Among("eisi", 55, -1),
        Among("tumeisi", 57, -1),
        Among("uisi", 55, -1),
        Among("osi", 49, -1),
        Among("ėjosi", 60, -1),
        Among("uosi", 60, -1),
        Among("iuosi", 62, -1),
        Among("siuosi", 63, -1),
        Among("usi", 49, -1),
        Among("ausi", 65, -1),
        Among("čiausi", 66, -1),
        Among("ąsi", 49, -1),
        Among("ėsi", 49, -1),
        Among("ųsi", 49, -1),
        Among("tųsi", 70, -1),
        Among("ti", 40, -1),
        Among("enti", 72, -1),
        Among("inti", 72, -1),
        Among("oti", 72, -1),
        Among("ioti", 75, -1),
        Among("uoti", 75, -1),
        Among("iuoti", 77, -1),
        Among("auti", 72, -1),
        Among("iauti", 79, -1),
        Among("yti", 72, -1),
        Among("ėti", 72, -1),
        Among("telėti", 82, -1),
        Among("inėti", 82, -1),
        Among("terėti", 82, -1),
        Among("ui", 40, -1),
        Among("iui", 86, -1),
        Among("eniui", 87, -1),
        Among("oj", -1, -1),
        Among("ėj", -1, -1),
        Among("k", -1, -1),
        Among("am", -1, -1),
        Among("iam", 92, -1),
        Among("iem", -1, -1),
        Among("im", -1, -1),
        Among("sim", 95, -1),
        Among("om", -1, -1),
        Among("tum", -1, -1),
        Among("ėm", -1, -1),
        Among("tumėm", 99, -1),
        Among("an", -1, -1),
        Among("on", -1, -1),
        Among("ion", 102, -1),
        Among("un", -1, -1),
        Among("iun", 104, -1),
        Among("ėn", -1, -1),
        Among("o", -1, -1),
        Among("io", 107, -1),
        Among("enio", 108, -1),
        Among("ėjo", 107, -1),
        Among("uo", 107, -1),
        Among("s", -1, -1),
        Among("as", 112, -1),
        Among("ias", 113, -1),
        Among("es", 112, -1),
        Among("ies", 115, -1),
        Among("is", 112, -1),
        Among("ais", 117, -1),
        Among("iais", 118, -1),
        Among("tumeis", 117, -1),
        Among("imis", 117, -1),
        Among("enimis", 121, -1),
        Among("omis", 117, -1),
        Among("iomis", 123, -1),
        Among("umis", 117, -1),
        Among("ėmis", 117, -1),
        Among("enis", 117, -1),
        Among("asis", 117, -1),
        Among("ysis", 117, -1),
        Among("ams", 112, -1),
        Among("iams", 130, -1),
        Among("iems", 112, -1),
        Among("ims", 112, -1),
        Among("enims", 133, -1),
        Among("oms", 112, -1),
        Among("ioms", 135, -1),
        Among("ums", 112, -1),
        Among("ėms", 112, -1),
        Among("ens", 112, -1),
        Among("os", 112, -1),
        Among("ios", 140, -1),
        Among("uos", 140, -1),
        Among("iuos", 142, -1),
        Among("us", 112, -1),
        Among("aus", 144, -1),
        Among("iaus", 145, -1),
        Among("ius", 144, -1),
        Among("ys", 112, -1),
        Among("enys", 148, -1),
        Among("ąs", 112, -1),
        Among("iąs", 150, -1),
        Among("ės", 112, -1),
        Among("amės", 152, -1),
        Among("iamės", 153, -1),
        Among("imės", 152, -1),
        Among("kimės", 155, -1),
        Among("simės", 155, -1),
        Among("omės", 152, -1),
        Among("ėmės", 152, -1),
        Among("tumėmės", 159, -1),
        Among("atės", 152, -1),
        Among("iatės", 161, -1),
        Among("sitės", 152, -1),
        Among("otės", 152, -1),
        Among("ėtės", 152, -1),
        Among("tumėtės", 165, -1),
        Among("įs", 112, -1),
        Among("ūs", 112, -1),
        Among("tųs", 112, -1),
        Among("at", -1, -1),
        Among("iat", 170, -1),
        Among("it", -1, -1),
        Among("sit", 172, -1),
        Among("ot", -1, -1),
        Among("ėt", -1, -1),
        Among("tumėt", 175, -1),
        Among("u", -1, -1),
        Among("au", 177, -1),
        Among("iau", 178, -1),
        Among("čiau", 179, -1),
        Among("iu", 177, -1),
        Among("eniu", 181, -1),
        Among("siu", 181, -1),
        Among("y", -1, -1),
        Among("ą", -1, -1),
        Among("ią", 185, -1),
        Among("ė", -1, -1),
        Among("ę", -1, -1),
        Among("į", -1, -1),
        Among("enį", 189, -1),
        Among("ų", -1, -1),
        Among("ių", 191, -1)
    ]

    a_1 = [
        Among("ing", -1, -1),
        Among("aj", -1, -1),
        Among("iaj", 1, -1),
        Among("iej", -1, -1),
        Among("oj", -1, -1),
        Among("ioj", 4, -1),
        Among("uoj", 4, -1),
        Among("iuoj", 6, -1),
        Among("auj", -1, -1),
        Among("ąj", -1, -1),
        Among("iąj", 9, -1),
        Among("ėj", -1, -1),
        Among("ųj", -1, -1),
        Among("iųj", 12, -1),
        Among("ok", -1, -1),
        Among("iok", 14, -1),
        Among("iuk", -1, -1),
        Among("uliuk", 16, -1),
        Among("učiuk", 16, -1),
        Among("išk", -1, -1),
        Among("iul", -1, -1),
        Among("yl", -1, -1),
        Among("ėl", -1, -1),
        Among("am", -1, -1),
        Among("dam", 23, -1),
        Among("jam", 23, -1),
        Among("zgan", -1, -1),
        Among("ain", -1, -1),
        Among("esn", -1, -1),
        Among("op", -1, -1),
        Among("iop", 29, -1),
        Among("ias", -1, -1),
        Among("ies", -1, -1),
        Among("ais", -1, -1),
        Among("iais", 33, -1),
        Among("os", -1, -1),
        Among("ios", 35, -1),
        Among("uos", 35, -1),
        Among("iuos", 37, -1),
        Among("aus", -1, -1),
        Among("iaus", 39, -1),
        Among("ąs", -1, -1),
        Among("iąs", 41, -1),
        Among("ęs", -1, -1),
        Among("utėait", -1, -1),
        Among("ant", -1, -1),
        Among("iant", 45, -1),
        Among("siant", 46, -1),
        Among("int", -1, -1),
        Among("ot", -1, -1),
        Among("uot", 49, -1),
        Among("iuot", 50, -1),
        Among("yt", -1, -1),
        Among("ėt", -1, -1),
        Among("ykšt", -1, -1),
        Among("iau", -1, -1),
        Among("dav", -1, -1),
        Among("sv", -1, -1),
        Among("šv", -1, -1),
        Among("ykšč", -1, -1),
        Among("ę", -1, -1),
        Among("ėję", 60, -1)
    ]

    a_2 = [
        Among("ojime", -1, 7),
        Among("ėjime", -1, 3),
        Among("avime", -1, 6),
        Among("okate", -1, 8),
        Among("aite", -1, 1),
        Among("uote", -1, 2),
        Among("asius", -1, 5),
        Among("okatės", -1, 8),
        Among("aitės", -1, 1),
        Among("uotės", -1, 2),
        Among("esiu", -1, 4)
    ]
    as_2 = ("aitė", "uotė", "ėjimas", "esys", "asys", "avimas", "ojimas", "okatė")

    a_3 = [
        Among("č", -1, 1),
        Among("dž", -1, 2)
    ]
    as_3 = ("t", "d")


class lab0(BaseException): pass


class lab1(BaseException): pass
