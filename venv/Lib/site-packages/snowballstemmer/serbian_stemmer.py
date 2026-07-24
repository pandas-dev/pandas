# Generated from serbian.sbl by Snowball 3.1.1 - https://snowballstem.org/

from .basestemmer import BaseStemmer
from .among import Among


class SerbianStemmer(BaseStemmer):
    '''
    This class implements the stemming algorithm defined by a snowball script.
    Generated from serbian.sbl by Snowball 3.1.1 - https://snowballstem.org/
    '''

    g_v = {"a", "e", "i", "o", "u"}

    g_sa = {"ć", "č", "đ", "š", "ž"}

    g_ca = {"b", "c", "d", "f", "g", "h", "j", "k", "l", "m", "n", "p", "r", "s", "t", "v", "z", "ć", "č", "đ", "š", "ž"}

    g_rg = "r"
    I_p1 = 0
    B_no_diacritics = False

    def __r_cyr_to_lat(self):
        v_1 = self.cursor
        try:
            while True:
                v_2 = self.cursor
                try:
                    while True:
                        v_3 = self.cursor
                        try:
                            self.bra = self.cursor
                            among_var = self.find_among(SerbianStemmer.a_0)
                            if among_var == 0:
                                raise lab2()
                            self.ket = self.cursor
                            self.slice_from(SerbianStemmer.as_0[among_var - 1])
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
        v_1 = self.cursor
        try:
            while True:
                v_2 = self.cursor
                try:
                    while True:
                        v_3 = self.cursor
                        try:
                            if not self.in_grouping(SerbianStemmer.g_ca):
                                raise lab2()
                            self.bra = self.cursor
                            if not self.eq_s("ije"):
                                raise lab2()
                            self.ket = self.cursor
                            if not self.in_grouping(SerbianStemmer.g_ca):
                                raise lab2()
                            self.slice_from("e")
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
        v_4 = self.cursor
        try:
            while True:
                v_5 = self.cursor
                try:
                    while True:
                        v_6 = self.cursor
                        try:
                            if not self.in_grouping(SerbianStemmer.g_ca):
                                raise lab2()
                            self.bra = self.cursor
                            if not self.eq_s("je"):
                                raise lab2()
                            self.ket = self.cursor
                            if not self.in_grouping(SerbianStemmer.g_ca):
                                raise lab2()
                            self.slice_from("e")
                            self.cursor = v_6
                            break
                        except lab2: pass
                        self.cursor = v_6
                        if self.cursor >= self.limit:
                            raise lab1()
                        self.cursor += 1
                    continue
                except lab1: pass
                self.cursor = v_5
                break
        except lab0: pass
        self.cursor = v_4
        v_7 = self.cursor
        try:
            while True:
                v_8 = self.cursor
                try:
                    while True:
                        v_9 = self.cursor
                        try:
                            self.bra = self.cursor
                            if not self.eq_s("dj"):
                                raise lab2()
                            self.ket = self.cursor
                            self.slice_from("đ")
                            self.cursor = v_9
                            break
                        except lab2: pass
                        self.cursor = v_9
                        if self.cursor >= self.limit:
                            raise lab1()
                        self.cursor += 1
                    continue
                except lab1: pass
                self.cursor = v_8
                break
        except lab0: pass
        self.cursor = v_7
        return True

    def __r_mark_regions(self):
        self.B_no_diacritics = True
        v_1 = self.cursor
        try:
            if not self.go_out_grouping(SerbianStemmer.g_sa):
                raise lab0()
            self.cursor += 1
            self.B_no_diacritics = False
        except lab0: pass
        self.cursor = v_1
        self.I_p1 = self.limit
        v_2 = self.cursor
        try:
            if not self.go_out_grouping(SerbianStemmer.g_v):
                raise lab0()
            self.cursor += 1
            self.I_p1 = self.cursor
            if self.I_p1 >= 2:
                raise lab0()
            if not self.go_in_grouping(SerbianStemmer.g_v):
                raise lab0()
            self.cursor += 1
            self.I_p1 = self.cursor
        except lab0: pass
        self.cursor = v_2
        v_3 = self.cursor
        try:
            while True:
                try:
                    if self.cursor == self.limit or self.current[self.cursor] != "r":
                        raise lab1()
                    self.cursor += 1
                    break
                except lab1: pass
                if self.cursor >= self.limit:
                    raise lab0()
                self.cursor += 1
            while True:
                try:
                    if self.cursor < 2:
                        raise lab1()
                    break
                except lab1: pass
                if not self.go_in_grouping(SerbianStemmer.g_rg):
                    raise lab0()
                self.cursor += 1
                break
            if (self.I_p1 - self.cursor) <= 1:
                raise lab0()
            self.I_p1 = self.cursor
        except lab0: pass
        self.cursor = v_3
        return True

    def __r_R1(self):
        return self.I_p1 <= self.cursor

    def __r_Step_1(self):
        self.ket = self.cursor
        among_var = self.find_among_b(SerbianStemmer.a_1)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var == 1:
            self.slice_from("loga")
        elif among_var == 2:
            self.slice_from("peh")
        elif among_var == 3:
            self.slice_from("vojka")
        elif among_var == 4:
            self.slice_from("bojka")
        elif among_var == 5:
            self.slice_from("jak")
        elif among_var == 6:
            self.slice_from("čajni")
        elif among_var == 7:
            if not self.B_no_diacritics:
                return False
            self.slice_from("cajni")
        elif among_var == 8:
            self.slice_from("erni")
        elif among_var == 9:
            self.slice_from("larni")
        elif among_var == 10:
            self.slice_from("esni")
        elif among_var == 11:
            self.slice_from("anjca")
        elif among_var == 12:
            self.slice_from("ajca")
        elif among_var == 13:
            self.slice_from("ljca")
        elif among_var == 14:
            self.slice_from("ejca")
        elif among_var == 15:
            self.slice_from("ojca")
        elif among_var == 16:
            self.slice_from("ajka")
        elif among_var == 17:
            self.slice_from("ojka")
        elif among_var == 18:
            self.slice_from("šca")
        elif among_var == 19:
            self.slice_from("ing")
        elif among_var == 20:
            self.slice_from("tvenik")
        elif among_var == 21:
            self.slice_from("tetika")
        elif among_var == 22:
            self.slice_from("nstva")
        elif among_var == 23:
            self.slice_from("nik")
        elif among_var == 24:
            self.slice_from("tik")
        elif among_var == 25:
            self.slice_from("zik")
        elif among_var == 26:
            self.slice_from("snik")
        elif among_var == 27:
            self.slice_from("kusi")
        elif among_var == 28:
            self.slice_from("kusni")
        elif among_var == 29:
            self.slice_from("kustva")
        elif among_var == 30:
            self.slice_from("dušni")
        elif among_var == 31:
            if not self.B_no_diacritics:
                return False
            self.slice_from("dusni")
        elif among_var == 32:
            self.slice_from("antni")
        elif among_var == 33:
            self.slice_from("bilni")
        elif among_var == 34:
            self.slice_from("tilni")
        elif among_var == 35:
            self.slice_from("avilni")
        elif among_var == 36:
            self.slice_from("silni")
        elif among_var == 37:
            self.slice_from("gilni")
        elif among_var == 38:
            self.slice_from("rilni")
        elif among_var == 39:
            self.slice_from("nilni")
        elif among_var == 40:
            self.slice_from("alni")
        elif among_var == 41:
            self.slice_from("ozni")
        elif among_var == 42:
            self.slice_from("ravi")
        elif among_var == 43:
            self.slice_from("stavni")
        elif among_var == 44:
            self.slice_from("pravni")
        elif among_var == 45:
            self.slice_from("tivni")
        elif among_var == 46:
            self.slice_from("sivni")
        elif among_var == 47:
            self.slice_from("atni")
        elif among_var == 48:
            self.slice_from("enta")
        elif among_var == 49:
            self.slice_from("tetni")
        elif among_var == 50:
            self.slice_from("pletni")
        elif among_var == 51:
            self.slice_from("šavi")
        elif among_var == 52:
            if not self.B_no_diacritics:
                return False
            self.slice_from("savi")
        elif among_var == 53:
            self.slice_from("anta")
        elif among_var == 54:
            self.slice_from("ačka")
        elif among_var == 55:
            if not self.B_no_diacritics:
                return False
            self.slice_from("acka")
        elif among_var == 56:
            self.slice_from("uška")
        elif among_var == 57:
            if not self.B_no_diacritics:
                return False
            self.slice_from("uska")
        elif among_var == 58:
            self.slice_from("atka")
        elif among_var == 59:
            self.slice_from("etka")
        elif among_var == 60:
            self.slice_from("itka")
        elif among_var == 61:
            self.slice_from("otka")
        elif among_var == 62:
            self.slice_from("utka")
        elif among_var == 63:
            self.slice_from("eskna")
        elif among_var == 64:
            self.slice_from("tični")
        elif among_var == 65:
            if not self.B_no_diacritics:
                return False
            self.slice_from("ticni")
        elif among_var == 66:
            self.slice_from("ojska")
        elif among_var == 67:
            self.slice_from("esma")
        elif among_var == 68:
            self.slice_from("metra")
        elif among_var == 69:
            self.slice_from("centra")
        elif among_var == 70:
            self.slice_from("istra")
        elif among_var == 71:
            self.slice_from("osti")
        elif among_var == 72:
            if not self.B_no_diacritics:
                return False
            self.slice_from("osti")
        elif among_var == 73:
            self.slice_from("dba")
        elif among_var == 74:
            self.slice_from("čka")
        elif among_var == 75:
            self.slice_from("mca")
        elif among_var == 76:
            self.slice_from("nca")
        elif among_var == 77:
            self.slice_from("voljni")
        elif among_var == 78:
            self.slice_from("anki")
        elif among_var == 79:
            self.slice_from("vca")
        elif among_var == 80:
            self.slice_from("sca")
        elif among_var == 81:
            self.slice_from("rca")
        elif among_var == 82:
            self.slice_from("alca")
        elif among_var == 83:
            self.slice_from("elca")
        elif among_var == 84:
            self.slice_from("olca")
        elif among_var == 85:
            self.slice_from("njca")
        elif among_var == 86:
            self.slice_from("ekta")
        elif among_var == 87:
            self.slice_from("izma")
        elif among_var == 88:
            self.slice_from("jebi")
        elif among_var == 89:
            self.slice_from("baci")
        elif among_var == 90:
            self.slice_from("ašni")
        else:
            if not self.B_no_diacritics:
                return False
            self.slice_from("asni")
        return True

    def __r_Step_2(self):
        self.ket = self.cursor
        among_var = self.find_among_b(SerbianStemmer.a_2)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if not self.__r_R1():
            return False
        if among_var == 1:
            self.slice_from("sk")
        elif among_var == 2:
            self.slice_from("šk")
        elif among_var == 3:
            self.slice_from("stv")
        elif among_var == 4:
            self.slice_from("štv")
        elif among_var == 5:
            self.slice_from("tanij")
        elif among_var == 6:
            self.slice_from("manij")
        elif among_var == 7:
            self.slice_from("panij")
        elif among_var == 8:
            self.slice_from("ranij")
        elif among_var == 9:
            self.slice_from("ganij")
        elif among_var == 10:
            self.slice_from("an")
        elif among_var == 11:
            self.slice_from("in")
        elif among_var == 12:
            self.slice_from("on")
        elif among_var == 13:
            self.slice_from("n")
        elif among_var == 14:
            self.slice_from("ać")
        elif among_var == 15:
            self.slice_from("eć")
        elif among_var == 16:
            self.slice_from("uć")
        elif among_var == 17:
            self.slice_from("ugov")
        elif among_var == 18:
            self.slice_from("ug")
        elif among_var == 19:
            self.slice_from("log")
        elif among_var == 20:
            self.slice_from("g")
        elif among_var == 21:
            self.slice_from("rari")
        elif among_var == 22:
            self.slice_from("oti")
        elif among_var == 23:
            self.slice_from("si")
        elif among_var == 24:
            self.slice_from("li")
        elif among_var == 25:
            self.slice_from("uj")
        elif among_var == 26:
            self.slice_from("caj")
        elif among_var == 27:
            self.slice_from("čaj")
        elif among_var == 28:
            self.slice_from("ćaj")
        elif among_var == 29:
            self.slice_from("đaj")
        elif among_var == 30:
            self.slice_from("laj")
        elif among_var == 31:
            self.slice_from("raj")
        elif among_var == 32:
            self.slice_from("bij")
        elif among_var == 33:
            self.slice_from("cij")
        elif among_var == 34:
            self.slice_from("dij")
        elif among_var == 35:
            self.slice_from("lij")
        elif among_var == 36:
            self.slice_from("nij")
        elif among_var == 37:
            self.slice_from("mij")
        elif among_var == 38:
            self.slice_from("žij")
        elif among_var == 39:
            self.slice_from("gij")
        elif among_var == 40:
            self.slice_from("fij")
        elif among_var == 41:
            self.slice_from("pij")
        elif among_var == 42:
            self.slice_from("rij")
        elif among_var == 43:
            self.slice_from("sij")
        elif among_var == 44:
            self.slice_from("tij")
        elif among_var == 45:
            self.slice_from("zij")
        elif among_var == 46:
            self.slice_from("nal")
        elif among_var == 47:
            self.slice_from("ijal")
        elif among_var == 48:
            self.slice_from("ozil")
        elif among_var == 49:
            self.slice_from("olov")
        elif among_var == 50:
            self.slice_from("ol")
        elif among_var == 51:
            self.slice_from("lem")
        elif among_var == 52:
            self.slice_from("ram")
        elif among_var == 53:
            self.slice_from("ar")
        elif among_var == 54:
            self.slice_from("dr")
        elif among_var == 55:
            self.slice_from("er")
        elif among_var == 56:
            self.slice_from("or")
        elif among_var == 57:
            self.slice_from("es")
        elif among_var == 58:
            self.slice_from("is")
        elif among_var == 59:
            self.slice_from("taš")
        elif among_var == 60:
            self.slice_from("naš")
        elif among_var == 61:
            self.slice_from("jaš")
        elif among_var == 62:
            self.slice_from("kaš")
        elif among_var == 63:
            self.slice_from("baš")
        elif among_var == 64:
            self.slice_from("gaš")
        elif among_var == 65:
            self.slice_from("vaš")
        elif among_var == 66:
            self.slice_from("eš")
        elif among_var == 67:
            self.slice_from("iš")
        elif among_var == 68:
            self.slice_from("ikat")
        elif among_var == 69:
            self.slice_from("lat")
        elif among_var == 70:
            self.slice_from("et")
        elif among_var == 71:
            self.slice_from("est")
        elif among_var == 72:
            self.slice_from("ist")
        elif among_var == 73:
            self.slice_from("kst")
        elif among_var == 74:
            self.slice_from("ost")
        elif among_var == 75:
            self.slice_from("išt")
        elif among_var == 76:
            self.slice_from("ova")
        elif among_var == 77:
            self.slice_from("av")
        elif among_var == 78:
            self.slice_from("ev")
        elif among_var == 79:
            self.slice_from("iv")
        elif among_var == 80:
            self.slice_from("ov")
        elif among_var == 81:
            self.slice_from("mov")
        elif among_var == 82:
            self.slice_from("lov")
        elif among_var == 83:
            self.slice_from("el")
        elif among_var == 84:
            self.slice_from("anj")
        elif among_var == 85:
            self.slice_from("enj")
        elif among_var == 86:
            self.slice_from("šnj")
        elif among_var == 87:
            self.slice_from("en")
        elif among_var == 88:
            self.slice_from("šn")
        elif among_var == 89:
            self.slice_from("čin")
        elif among_var == 90:
            self.slice_from("roši")
        elif among_var == 91:
            self.slice_from("oš")
        elif among_var == 92:
            self.slice_from("evit")
        elif among_var == 93:
            self.slice_from("ovit")
        elif among_var == 94:
            self.slice_from("ast")
        elif among_var == 95:
            self.slice_from("k")
        elif among_var == 96:
            self.slice_from("eva")
        elif among_var == 97:
            self.slice_from("ava")
        elif among_var == 98:
            self.slice_from("iva")
        elif among_var == 99:
            self.slice_from("uva")
        elif among_var == 100:
            self.slice_from("ir")
        elif among_var == 101:
            self.slice_from("ač")
        elif among_var == 102:
            self.slice_from("ača")
        elif among_var == 103:
            self.slice_from("ni")
        elif among_var == 104:
            self.slice_from("a")
        elif among_var == 105:
            self.slice_from("ur")
        elif among_var == 106:
            self.slice_from("astaj")
        elif among_var == 107:
            self.slice_from("istaj")
        elif among_var == 108:
            self.slice_from("ostaj")
        elif among_var == 109:
            self.slice_from("aj")
        elif among_var == 110:
            self.slice_from("asta")
        elif among_var == 111:
            self.slice_from("ista")
        elif among_var == 112:
            self.slice_from("osta")
        elif among_var == 113:
            self.slice_from("ta")
        elif among_var == 114:
            self.slice_from("inj")
        elif among_var == 115:
            self.slice_from("as")
        elif among_var == 116:
            self.slice_from("i")
        elif among_var == 117:
            self.slice_from("luč")
        elif among_var == 118:
            self.slice_from("jeti")
        elif among_var == 119:
            self.slice_from("e")
        elif among_var == 120:
            self.slice_from("at")
        elif among_var == 121:
            if not self.B_no_diacritics:
                return False
            self.slice_from("luc")
        elif among_var == 122:
            if not self.B_no_diacritics:
                return False
            self.slice_from("snj")
        elif among_var == 123:
            if not self.B_no_diacritics:
                return False
            self.slice_from("os")
        elif among_var == 124:
            if not self.B_no_diacritics:
                return False
            self.slice_from("ac")
        elif among_var == 125:
            if not self.B_no_diacritics:
                return False
            self.slice_from("ec")
        elif among_var == 126:
            if not self.B_no_diacritics:
                return False
            self.slice_from("uc")
        elif among_var == 127:
            if not self.B_no_diacritics:
                return False
            self.slice_from("rosi")
        elif among_var == 128:
            if not self.B_no_diacritics:
                return False
            self.slice_from("aca")
        elif among_var == 129:
            if not self.B_no_diacritics:
                return False
            self.slice_from("jas")
        elif among_var == 130:
            if not self.B_no_diacritics:
                return False
            self.slice_from("tas")
        elif among_var == 131:
            if not self.B_no_diacritics:
                return False
            self.slice_from("gas")
        elif among_var == 132:
            if not self.B_no_diacritics:
                return False
            self.slice_from("nas")
        elif among_var == 133:
            if not self.B_no_diacritics:
                return False
            self.slice_from("kas")
        elif among_var == 134:
            if not self.B_no_diacritics:
                return False
            self.slice_from("vas")
        elif among_var == 135:
            if not self.B_no_diacritics:
                return False
            self.slice_from("bas")
        elif among_var == 136:
            if not self.B_no_diacritics:
                return False
            self.slice_from("as")
        elif among_var == 137:
            if not self.B_no_diacritics:
                return False
            self.slice_from("cin")
        elif among_var == 138:
            if not self.B_no_diacritics:
                return False
            self.slice_from("astaj")
        elif among_var == 139:
            if not self.B_no_diacritics:
                return False
            self.slice_from("istaj")
        elif among_var == 140:
            if not self.B_no_diacritics:
                return False
            self.slice_from("ostaj")
        elif among_var == 141:
            if not self.B_no_diacritics:
                return False
            self.slice_from("asta")
        elif among_var == 142:
            if not self.B_no_diacritics:
                return False
            self.slice_from("ista")
        elif among_var == 143:
            if not self.B_no_diacritics:
                return False
            self.slice_from("osta")
        elif among_var == 144:
            if not self.B_no_diacritics:
                return False
            self.slice_from("ava")
        elif among_var == 145:
            if not self.B_no_diacritics:
                return False
            self.slice_from("eva")
        elif among_var == 146:
            if not self.B_no_diacritics:
                return False
            self.slice_from("iva")
        elif among_var == 147:
            if not self.B_no_diacritics:
                return False
            self.slice_from("uva")
        elif among_var == 148:
            if not self.B_no_diacritics:
                return False
            self.slice_from("ova")
        elif among_var == 149:
            if not self.B_no_diacritics:
                return False
            self.slice_from("jeti")
        elif among_var == 150:
            if not self.B_no_diacritics:
                return False
            self.slice_from("inj")
        elif among_var == 151:
            if not self.B_no_diacritics:
                return False
            self.slice_from("ist")
        elif among_var == 152:
            if not self.B_no_diacritics:
                return False
            self.slice_from("es")
        elif among_var == 153:
            if not self.B_no_diacritics:
                return False
            self.slice_from("et")
        elif among_var == 154:
            if not self.B_no_diacritics:
                return False
            self.slice_from("is")
        elif among_var == 155:
            if not self.B_no_diacritics:
                return False
            self.slice_from("ir")
        elif among_var == 156:
            if not self.B_no_diacritics:
                return False
            self.slice_from("ur")
        elif among_var == 157:
            if not self.B_no_diacritics:
                return False
            self.slice_from("uj")
        elif among_var == 158:
            if not self.B_no_diacritics:
                return False
            self.slice_from("ni")
        elif among_var == 159:
            if not self.B_no_diacritics:
                return False
            self.slice_from("sn")
        elif among_var == 160:
            if not self.B_no_diacritics:
                return False
            self.slice_from("ta")
        elif among_var == 161:
            if not self.B_no_diacritics:
                return False
            self.slice_from("a")
        elif among_var == 162:
            if not self.B_no_diacritics:
                return False
            self.slice_from("i")
        elif among_var == 163:
            if not self.B_no_diacritics:
                return False
            self.slice_from("e")
        else:
            if not self.B_no_diacritics:
                return False
            self.slice_from("n")
        return True

    def __r_Step_3(self):
        self.ket = self.cursor
        if self.find_among_b(SerbianStemmer.a_3) == 0:
            return False
        self.bra = self.cursor
        if not self.__r_R1():
            return False
        self.slice_del()
        return True

    def _stem(self):
        self.__r_cyr_to_lat()
        self.__r_prelude()
        self.__r_mark_regions()
        self.limit_backward = self.cursor
        self.cursor = self.limit
        v_1 = self.limit - self.cursor
        self.__r_Step_1()
        self.cursor = self.limit - v_1
        v_2 = self.limit - self.cursor
        try:
            while True:
                v_3 = self.limit - self.cursor
                try:
                    if not self.__r_Step_2():
                        raise lab1()
                    break
                except lab1: pass
                self.cursor = self.limit - v_3
                if not self.__r_Step_3():
                    raise lab0()
                break
        except lab0: pass
        self.cursor = self.limit - v_2
        self.cursor = self.limit_backward
        return True

    a_0 = [
        Among("а", -1, 1),
        Among("б", -1, 2),
        Among("в", -1, 3),
        Among("г", -1, 4),
        Among("д", -1, 5),
        Among("е", -1, 7),
        Among("ж", -1, 8),
        Among("з", -1, 9),
        Among("и", -1, 10),
        Among("к", -1, 12),
        Among("л", -1, 13),
        Among("м", -1, 15),
        Among("н", -1, 16),
        Among("о", -1, 18),
        Among("п", -1, 19),
        Among("р", -1, 20),
        Among("с", -1, 21),
        Among("т", -1, 22),
        Among("у", -1, 24),
        Among("ф", -1, 25),
        Among("х", -1, 26),
        Among("ц", -1, 27),
        Among("ч", -1, 28),
        Among("ш", -1, 30),
        Among("ђ", -1, 6),
        Among("ј", -1, 11),
        Among("љ", -1, 14),
        Among("њ", -1, 17),
        Among("ћ", -1, 23),
        Among("џ", -1, 29)
    ]
    as_0 = ("a", "b", "v", "g", "d", "đ", "e", "ž", "z", "i", "j", "k", "l", "lj", "m", "n", "nj", "o", "p", "r", "s", "t", "ć", "u", "f", "h", "c", "č", "dž", "š")

    a_1 = [
        Among("daba", -1, 73),
        Among("ajaca", -1, 12),
        Among("ejaca", -1, 14),
        Among("ljaca", -1, 13),
        Among("njaca", -1, 85),
        Among("ojaca", -1, 15),
        Among("alaca", -1, 82),
        Among("elaca", -1, 83),
        Among("olaca", -1, 84),
        Among("maca", -1, 75),
        Among("naca", -1, 76),
        Among("raca", -1, 81),
        Among("saca", -1, 80),
        Among("vaca", -1, 79),
        Among("šaca", -1, 18),
        Among("aoca", -1, 82),
        Among("acaka", -1, 55),
        Among("ajaka", -1, 16),
        Among("ojaka", -1, 17),
        Among("anaka", -1, 78),
        Among("ataka", -1, 58),
        Among("etaka", -1, 59),
        Among("itaka", -1, 60),
        Among("otaka", -1, 61),
        Among("utaka", -1, 62),
        Among("ačaka", -1, 54),
        Among("esama", -1, 67),
        Among("izama", -1, 87),
        Among("jacima", -1, 5),
        Among("nicima", -1, 23),
        Among("ticima", -1, 24),
        Among("teticima", 30, 21),
        Among("zicima", -1, 25),
        Among("atcima", -1, 58),
        Among("utcima", -1, 62),
        Among("čcima", -1, 74),
        Among("pesima", -1, 2),
        Among("inzima", -1, 19),
        Among("lozima", -1, 1),
        Among("metara", -1, 68),
        Among("centara", -1, 69),
        Among("istara", -1, 70),
        Among("ekata", -1, 86),
        Among("anata", -1, 53),
        Among("nstava", -1, 22),
        Among("kustava", -1, 29),
        Among("ajac", -1, 12),
        Among("ejac", -1, 14),
        Among("ljac", -1, 13),
        Among("njac", -1, 85),
        Among("anjac", 49, 11),
        Among("ojac", -1, 15),
        Among("alac", -1, 82),
        Among("elac", -1, 83),
        Among("olac", -1, 84),
        Among("mac", -1, 75),
        Among("nac", -1, 76),
        Among("rac", -1, 81),
        Among("sac", -1, 80),
        Among("vac", -1, 79),
        Among("šac", -1, 18),
        Among("jebe", -1, 88),
        Among("olce", -1, 84),
        Among("kuse", -1, 27),
        Among("rave", -1, 42),
        Among("save", -1, 52),
        Among("šave", -1, 51),
        Among("baci", -1, 89),
        Among("jaci", -1, 5),
        Among("tvenici", -1, 20),
        Among("snici", -1, 26),
        Among("tetici", -1, 21),
        Among("bojci", -1, 4),
        Among("vojci", -1, 3),
        Among("ojsci", -1, 66),
        Among("atci", -1, 58),
        Among("itci", -1, 60),
        Among("utci", -1, 62),
        Among("čci", -1, 74),
        Among("pesi", -1, 2),
        Among("inzi", -1, 19),
        Among("lozi", -1, 1),
        Among("acak", -1, 55),
        Among("usak", -1, 57),
        Among("atak", -1, 58),
        Among("etak", -1, 59),
        Among("itak", -1, 60),
        Among("otak", -1, 61),
        Among("utak", -1, 62),
        Among("ačak", -1, 54),
        Among("ušak", -1, 56),
        Among("izam", -1, 87),
        Among("tican", -1, 65),
        Among("cajan", -1, 7),
        Among("čajan", -1, 6),
        Among("voljan", -1, 77),
        Among("eskan", -1, 63),
        Among("alan", -1, 40),
        Among("bilan", -1, 33),
        Among("gilan", -1, 37),
        Among("nilan", -1, 39),
        Among("rilan", -1, 38),
        Among("silan", -1, 36),
        Among("tilan", -1, 34),
        Among("avilan", -1, 35),
        Among("laran", -1, 9),
        Among("eran", -1, 8),
        Among("asan", -1, 91),
        Among("esan", -1, 10),
        Among("dusan", -1, 31),
        Among("kusan", -1, 28),
        Among("atan", -1, 47),
        Among("pletan", -1, 50),
        Among("tetan", -1, 49),
        Among("antan", -1, 32),
        Among("pravan", -1, 44),
        Among("stavan", -1, 43),
        Among("sivan", -1, 46),
        Among("tivan", -1, 45),
        Among("ozan", -1, 41),
        Among("tičan", -1, 64),
        Among("ašan", -1, 90),
        Among("dušan", -1, 30),
        Among("metar", -1, 68),
        Among("centar", -1, 69),
        Among("istar", -1, 70),
        Among("ekat", -1, 86),
        Among("enat", -1, 48),
        Among("oscu", -1, 72),
        Among("ošću", -1, 71)
    ]

    a_2 = [
        Among("aca", -1, 124),
        Among("eca", -1, 125),
        Among("uca", -1, 126),
        Among("ga", -1, 20),
        Among("acega", 3, 124),
        Among("ecega", 3, 125),
        Among("ucega", 3, 126),
        Among("anjijega", 3, 84),
        Among("enjijega", 3, 85),
        Among("snjijega", 3, 122),
        Among("šnjijega", 3, 86),
        Among("kijega", 3, 95),
        Among("skijega", 11, 1),
        Among("škijega", 11, 2),
        Among("elijega", 3, 83),
        Among("nijega", 3, 13),
        Among("osijega", 3, 123),
        Among("atijega", 3, 120),
        Among("evitijega", 3, 92),
        Among("ovitijega", 3, 93),
        Among("astijega", 3, 94),
        Among("avijega", 3, 77),
        Among("evijega", 3, 78),
        Among("ivijega", 3, 79),
        Among("ovijega", 3, 80),
        Among("ošijega", 3, 91),
        Among("anjega", 3, 84),
        Among("enjega", 3, 85),
        Among("snjega", 3, 122),
        Among("šnjega", 3, 86),
        Among("kega", 3, 95),
        Among("skega", 30, 1),
        Among("škega", 30, 2),
        Among("elega", 3, 83),
        Among("nega", 3, 13),
        Among("anega", 34, 10),
        Among("enega", 34, 87),
        Among("snega", 34, 159),
        Among("šnega", 34, 88),
        Among("osega", 3, 123),
        Among("atega", 3, 120),
        Among("evitega", 3, 92),
        Among("ovitega", 3, 93),
        Among("astega", 3, 94),
        Among("avega", 3, 77),
        Among("evega", 3, 78),
        Among("ivega", 3, 79),
        Among("ovega", 3, 80),
        Among("aćega", 3, 14),
        Among("ećega", 3, 15),
        Among("ućega", 3, 16),
        Among("ošega", 3, 91),
        Among("acoga", 3, 124),
        Among("ecoga", 3, 125),
        Among("ucoga", 3, 126),
        Among("anjoga", 3, 84),
        Among("enjoga", 3, 85),
        Among("snjoga", 3, 122),
        Among("šnjoga", 3, 86),
        Among("koga", 3, 95),
        Among("skoga", 59, 1),
        Among("škoga", 59, 2),
        Among("loga", 3, 19),
        Among("eloga", 62, 83),
        Among("noga", 3, 13),
        Among("cinoga", 64, 137),
        Among("činoga", 64, 89),
        Among("osoga", 3, 123),
        Among("atoga", 3, 120),
        Among("evitoga", 3, 92),
        Among("ovitoga", 3, 93),
        Among("astoga", 3, 94),
        Among("avoga", 3, 77),
        Among("evoga", 3, 78),
        Among("ivoga", 3, 79),
        Among("ovoga", 3, 80),
        Among("aćoga", 3, 14),
        Among("ećoga", 3, 15),
        Among("ućoga", 3, 16),
        Among("ošoga", 3, 91),
        Among("uga", 3, 18),
        Among("aja", -1, 109),
        Among("caja", 81, 26),
        Among("laja", 81, 30),
        Among("raja", 81, 31),
        Among("ćaja", 81, 28),
        Among("čaja", 81, 27),
        Among("đaja", 81, 29),
        Among("bija", -1, 32),
        Among("cija", -1, 33),
        Among("dija", -1, 34),
        Among("fija", -1, 40),
        Among("gija", -1, 39),
        Among("anjija", -1, 84),
        Among("enjija", -1, 85),
        Among("snjija", -1, 122),
        Among("šnjija", -1, 86),
        Among("kija", -1, 95),
        Among("skija", 97, 1),
        Among("škija", 97, 2),
        Among("lija", -1, 24),
        Among("elija", 100, 83),
        Among("mija", -1, 37),
        Among("nija", -1, 13),
        Among("ganija", 103, 9),
        Among("manija", 103, 6),
        Among("panija", 103, 7),
        Among("ranija", 103, 8),
        Among("tanija", 103, 5),
        Among("pija", -1, 41),
        Among("rija", -1, 42),
        Among("rarija", 110, 21),
        Among("sija", -1, 23),
        Among("osija", 112, 123),
        Among("tija", -1, 44),
        Among("atija", 114, 120),
        Among("evitija", 114, 92),
        Among("ovitija", 114, 93),
        Among("otija", 114, 22),
        Among("astija", 114, 94),
        Among("avija", -1, 77),
        Among("evija", -1, 78),
        Among("ivija", -1, 79),
        Among("ovija", -1, 80),
        Among("zija", -1, 45),
        Among("ošija", -1, 91),
        Among("žija", -1, 38),
        Among("anja", -1, 84),
        Among("enja", -1, 85),
        Among("snja", -1, 122),
        Among("šnja", -1, 86),
        Among("ka", -1, 95),
        Among("ska", 131, 1),
        Among("ška", 131, 2),
        Among("ala", -1, 104),
        Among("acala", 134, 128),
        Among("astajala", 134, 106),
        Among("istajala", 134, 107),
        Among("ostajala", 134, 108),
        Among("ijala", 134, 47),
        Among("injala", 134, 114),
        Among("nala", 134, 46),
        Among("irala", 134, 100),
        Among("urala", 134, 105),
        Among("tala", 134, 113),
        Among("astala", 144, 110),
        Among("istala", 144, 111),
        Among("ostala", 144, 112),
        Among("avala", 134, 97),
        Among("evala", 134, 96),
        Among("ivala", 134, 98),
        Among("ovala", 134, 76),
        Among("uvala", 134, 99),
        Among("ačala", 134, 102),
        Among("ela", -1, 83),
        Among("ila", -1, 116),
        Among("acila", 155, 124),
        Among("lucila", 155, 121),
        Among("nila", 155, 103),
        Among("astanila", 158, 110),
        Among("istanila", 158, 111),
        Among("ostanila", 158, 112),
        Among("rosila", 155, 127),
        Among("jetila", 155, 118),
        Among("ozila", 155, 48),
        Among("ačila", 155, 101),
        Among("lučila", 155, 117),
        Among("rošila", 155, 90),
        Among("ola", -1, 50),
        Among("asla", -1, 115),
        Among("nula", -1, 13),
        Among("gama", -1, 20),
        Among("logama", 171, 19),
        Among("ugama", 171, 18),
        Among("ajama", -1, 109),
        Among("cajama", 174, 26),
        Among("lajama", 174, 30),
        Among("rajama", 174, 31),
        Among("ćajama", 174, 28),
        Among("čajama", 174, 27),
        Among("đajama", 174, 29),
        Among("bijama", -1, 32),
        Among("cijama", -1, 33),
        Among("dijama", -1, 34),
        Among("fijama", -1, 40),
        Among("gijama", -1, 39),
        Among("lijama", -1, 35),
        Among("mijama", -1, 37),
        Among("nijama", -1, 36),
        Among("ganijama", 188, 9),
        Among("manijama", 188, 6),
        Among("panijama", 188, 7),
        Among("ranijama", 188, 8),
        Among("tanijama", 188, 5),
        Among("pijama", -1, 41),
        Among("rijama", -1, 42),
        Among("sijama", -1, 43),
        Among("tijama", -1, 44),
        Among("zijama", -1, 45),
        Among("žijama", -1, 38),
        Among("alama", -1, 104),
        Among("ijalama", 200, 47),
        Among("nalama", 200, 46),
        Among("elama", -1, 119),
        Among("ilama", -1, 116),
        Among("ramama", -1, 52),
        Among("lemama", -1, 51),
        Among("inama", -1, 11),
        Among("cinama", 207, 137),
        Among("činama", 207, 89),
        Among("rama", -1, 52),
        Among("arama", 210, 53),
        Among("drama", 210, 54),
        Among("erama", 210, 55),
        Among("orama", 210, 56),
        Among("basama", -1, 135),
        Among("gasama", -1, 131),
        Among("jasama", -1, 129),
        Among("kasama", -1, 133),
        Among("nasama", -1, 132),
        Among("tasama", -1, 130),
        Among("vasama", -1, 134),
        Among("esama", -1, 152),
        Among("isama", -1, 154),
        Among("etama", -1, 70),
        Among("estama", -1, 71),
        Among("istama", -1, 72),
        Among("kstama", -1, 73),
        Among("ostama", -1, 74),
        Among("avama", -1, 77),
        Among("evama", -1, 78),
        Among("ivama", -1, 79),
        Among("bašama", -1, 63),
        Among("gašama", -1, 64),
        Among("jašama", -1, 61),
        Among("kašama", -1, 62),
        Among("našama", -1, 60),
        Among("tašama", -1, 59),
        Among("vašama", -1, 65),
        Among("ešama", -1, 66),
        Among("išama", -1, 67),
        Among("lema", -1, 51),
        Among("acima", -1, 124),
        Among("ecima", -1, 125),
        Among("ucima", -1, 126),
        Among("ajima", -1, 109),
        Among("cajima", 245, 26),
        Among("lajima", 245, 30),
        Among("rajima", 245, 31),
        Among("ćajima", 245, 28),
        Among("čajima", 245, 27),
        Among("đajima", 245, 29),
        Among("bijima", -1, 32),
        Among("cijima", -1, 33),
        Among("dijima", -1, 34),
        Among("fijima", -1, 40),
        Among("gijima", -1, 39),
        Among("anjijima", -1, 84),
        Among("enjijima", -1, 85),
        Among("snjijima", -1, 122),
        Among("šnjijima", -1, 86),
        Among("kijima", -1, 95),
        Among("skijima", 261, 1),
        Among("škijima", 261, 2),
        Among("lijima", -1, 35),
        Among("elijima", 264, 83),
        Among("mijima", -1, 37),
        Among("nijima", -1, 13),
        Among("ganijima", 267, 9),
        Among("manijima", 267, 6),
        Among("panijima", 267, 7),
        Among("ranijima", 267, 8),
        Among("tanijima", 267, 5),
        Among("pijima", -1, 41),
        Among("rijima", -1, 42),
        Among("sijima", -1, 43),
        Among("osijima", 275, 123),
        Among("tijima", -1, 44),
        Among("atijima", 277, 120),
        Among("evitijima", 277, 92),
        Among("ovitijima", 277, 93),
        Among("astijima", 277, 94),
        Among("avijima", -1, 77),
        Among("evijima", -1, 78),
        Among("ivijima", -1, 79),
        Among("ovijima", -1, 80),
        Among("zijima", -1, 45),
        Among("ošijima", -1, 91),
        Among("žijima", -1, 38),
        Among("anjima", -1, 84),
        Among("enjima", -1, 85),
        Among("snjima", -1, 122),
        Among("šnjima", -1, 86),
        Among("kima", -1, 95),
        Among("skima", 293, 1),
        Among("škima", 293, 2),
        Among("alima", -1, 104),
        Among("ijalima", 296, 47),
        Among("nalima", 296, 46),
        Among("elima", -1, 83),
        Among("ilima", -1, 116),
        Among("ozilima", 300, 48),
        Among("olima", -1, 50),
        Among("lemima", -1, 51),
        Among("nima", -1, 13),
        Among("anima", 304, 10),
        Among("inima", 304, 11),
        Among("cinima", 306, 137),
        Among("činima", 306, 89),
        Among("onima", 304, 12),
        Among("arima", -1, 53),
        Among("drima", -1, 54),
        Among("erima", -1, 55),
        Among("orima", -1, 56),
        Among("basima", -1, 135),
        Among("gasima", -1, 131),
        Among("jasima", -1, 129),
        Among("kasima", -1, 133),
        Among("nasima", -1, 132),
        Among("tasima", -1, 130),
        Among("vasima", -1, 134),
        Among("esima", -1, 57),
        Among("isima", -1, 58),
        Among("osima", -1, 123),
        Among("atima", -1, 120),
        Among("ikatima", 324, 68),
        Among("latima", 324, 69),
        Among("etima", -1, 70),
        Among("evitima", -1, 92),
        Among("ovitima", -1, 93),
        Among("astima", -1, 94),
        Among("estima", -1, 71),
        Among("istima", -1, 72),
        Among("kstima", -1, 73),
        Among("ostima", -1, 74),
        Among("ištima", -1, 75),
        Among("avima", -1, 77),
        Among("evima", -1, 78),
        Among("ajevima", 337, 109),
        Among("cajevima", 338, 26),
        Among("lajevima", 338, 30),
        Among("rajevima", 338, 31),
        Among("ćajevima", 338, 28),
        Among("čajevima", 338, 27),
        Among("đajevima", 338, 29),
        Among("ivima", -1, 79),
        Among("ovima", -1, 80),
        Among("govima", 346, 20),
        Among("ugovima", 347, 17),
        Among("lovima", 346, 82),
        Among("olovima", 349, 49),
        Among("movima", 346, 81),
        Among("onovima", 346, 12),
        Among("stvima", -1, 3),
        Among("štvima", -1, 4),
        Among("aćima", -1, 14),
        Among("ećima", -1, 15),
        Among("ućima", -1, 16),
        Among("bašima", -1, 63),
        Among("gašima", -1, 64),
        Among("jašima", -1, 61),
        Among("kašima", -1, 62),
        Among("našima", -1, 60),
        Among("tašima", -1, 59),
        Among("vašima", -1, 65),
        Among("ešima", -1, 66),
        Among("išima", -1, 67),
        Among("ošima", -1, 91),
        Among("na", -1, 13),
        Among("ana", 368, 10),
        Among("acana", 369, 128),
        Among("urana", 369, 105),
        Among("tana", 369, 113),
        Among("avana", 369, 97),
        Among("evana", 369, 96),
        Among("ivana", 369, 98),
        Among("uvana", 369, 99),
        Among("ačana", 369, 102),
        Among("acena", 368, 124),
        Among("lucena", 368, 121),
        Among("ačena", 368, 101),
        Among("lučena", 368, 117),
        Among("ina", 368, 11),
        Among("cina", 382, 137),
        Among("anina", 382, 10),
        Among("čina", 382, 89),
        Among("ona", 368, 12),
        Among("ara", -1, 53),
        Among("dra", -1, 54),
        Among("era", -1, 55),
        Among("ora", -1, 56),
        Among("basa", -1, 135),
        Among("gasa", -1, 131),
        Among("jasa", -1, 129),
        Among("kasa", -1, 133),
        Among("nasa", -1, 132),
        Among("tasa", -1, 130),
        Among("vasa", -1, 134),
        Among("esa", -1, 57),
        Among("isa", -1, 58),
        Among("osa", -1, 123),
        Among("ata", -1, 120),
        Among("ikata", 401, 68),
        Among("lata", 401, 69),
        Among("eta", -1, 70),
        Among("evita", -1, 92),
        Among("ovita", -1, 93),
        Among("asta", -1, 94),
        Among("esta", -1, 71),
        Among("ista", -1, 72),
        Among("ksta", -1, 73),
        Among("osta", -1, 74),
        Among("nuta", -1, 13),
        Among("išta", -1, 75),
        Among("ava", -1, 77),
        Among("eva", -1, 78),
        Among("ajeva", 415, 109),
        Among("cajeva", 416, 26),
        Among("lajeva", 416, 30),
        Among("rajeva", 416, 31),
        Among("ćajeva", 416, 28),
        Among("čajeva", 416, 27),
        Among("đajeva", 416, 29),
        Among("iva", -1, 79),
        Among("ova", -1, 80),
        Among("gova", 424, 20),
        Among("ugova", 425, 17),
        Among("lova", 424, 82),
        Among("olova", 427, 49),
        Among("mova", 424, 81),
        Among("onova", 424, 12),
        Among("stva", -1, 3),
        Among("štva", -1, 4),
        Among("aća", -1, 14),
        Among("eća", -1, 15),
        Among("uća", -1, 16),
        Among("baša", -1, 63),
        Among("gaša", -1, 64),
        Among("jaša", -1, 61),
        Among("kaša", -1, 62),
        Among("naša", -1, 60),
        Among("taša", -1, 59),
        Among("vaša", -1, 65),
        Among("eša", -1, 66),
        Among("iša", -1, 67),
        Among("oša", -1, 91),
        Among("ace", -1, 124),
        Among("ece", -1, 125),
        Among("uce", -1, 126),
        Among("luce", 448, 121),
        Among("astade", -1, 110),
        Among("istade", -1, 111),
        Among("ostade", -1, 112),
        Among("ge", -1, 20),
        Among("loge", 453, 19),
        Among("uge", 453, 18),
        Among("aje", -1, 104),
        Among("caje", 456, 26),
        Among("laje", 456, 30),
        Among("raje", 456, 31),
        Among("astaje", 456, 106),
        Among("istaje", 456, 107),
        Among("ostaje", 456, 108),
        Among("ćaje", 456, 28),
        Among("čaje", 456, 27),
        Among("đaje", 456, 29),
        Among("ije", -1, 116),
        Among("bije", 466, 32),
        Among("cije", 466, 33),
        Among("dije", 466, 34),
        Among("fije", 466, 40),
        Among("gije", 466, 39),
        Among("anjije", 466, 84),
        Among("enjije", 466, 85),
        Among("snjije", 466, 122),
        Among("šnjije", 466, 86),
        Among("kije", 466, 95),
        Among("skije", 476, 1),
        Among("škije", 476, 2),
        Among("lije", 466, 35),
        Among("elije", 479, 83),
        Among("mije", 466, 37),
        Among("nije", 466, 13),
        Among("ganije", 482, 9),
        Among("manije", 482, 6),
        Among("panije", 482, 7),
        Among("ranije", 482, 8),
        Among("tanije", 482, 5),
        Among("pije", 466, 41),
        Among("rije", 466, 42),
        Among("sije", 466, 43),
        Among("osije", 490, 123),
        Among("tije", 466, 44),
        Among("atije", 492, 120),
        Among("evitije", 492, 92),
        Among("ovitije", 492, 93),
        Among("astije", 492, 94),
        Among("avije", 466, 77),
        Among("evije", 466, 78),
        Among("ivije", 466, 79),
        Among("ovije", 466, 80),
        Among("zije", 466, 45),
        Among("ošije", 466, 91),
        Among("žije", 466, 38),
        Among("anje", -1, 84),
        Among("enje", -1, 85),
        Among("snje", -1, 122),
        Among("šnje", -1, 86),
        Among("uje", -1, 25),
        Among("lucuje", 508, 121),
        Among("iruje", 508, 100),
        Among("lučuje", 508, 117),
        Among("ke", -1, 95),
        Among("ske", 512, 1),
        Among("ške", 512, 2),
        Among("ale", -1, 104),
        Among("acale", 515, 128),
        Among("astajale", 515, 106),
        Among("istajale", 515, 107),
        Among("ostajale", 515, 108),
        Among("ijale", 515, 47),
        Among("injale", 515, 114),
        Among("nale", 515, 46),
        Among("irale", 515, 100),
        Among("urale", 515, 105),
        Among("tale", 515, 113),
        Among("astale", 525, 110),
        Among("istale", 525, 111),
        Among("ostale", 525, 112),
        Among("avale", 515, 97),
        Among("evale", 515, 96),
        Among("ivale", 515, 98),
        Among("ovale", 515, 76),
        Among("uvale", 515, 99),
        Among("ačale", 515, 102),
        Among("ele", -1, 83),
        Among("ile", -1, 116),
        Among("acile", 536, 124),
        Among("lucile", 536, 121),
        Among("nile", 536, 103),
        Among("rosile", 536, 127),
        Among("jetile", 536, 118),
        Among("ozile", 536, 48),
        Among("ačile", 536, 101),
        Among("lučile", 536, 117),
        Among("rošile", 536, 90),
        Among("ole", -1, 50),
        Among("asle", -1, 115),
        Among("nule", -1, 13),
        Among("rame", -1, 52),
        Among("leme", -1, 51),
        Among("acome", -1, 124),
        Among("ecome", -1, 125),
        Among("ucome", -1, 126),
        Among("anjome", -1, 84),
        Among("enjome", -1, 85),
        Among("snjome", -1, 122),
        Among("šnjome", -1, 86),
        Among("kome", -1, 95),
        Among("skome", 558, 1),
        Among("škome", 558, 2),
        Among("elome", -1, 83),
        Among("nome", -1, 13),
        Among("cinome", 562, 137),
        Among("činome", 562, 89),
        Among("osome", -1, 123),
        Among("atome", -1, 120),
        Among("evitome", -1, 92),
        Among("ovitome", -1, 93),
        Among("astome", -1, 94),
        Among("avome", -1, 77),
        Among("evome", -1, 78),
        Among("ivome", -1, 79),
        Among("ovome", -1, 80),
        Among("aćome", -1, 14),
        Among("ećome", -1, 15),
        Among("ućome", -1, 16),
        Among("ošome", -1, 91),
        Among("ne", -1, 13),
        Among("ane", 578, 10),
        Among("acane", 579, 128),
        Among("urane", 579, 105),
        Among("tane", 579, 113),
        Among("astane", 582, 110),
        Among("istane", 582, 111),
        Among("ostane", 582, 112),
        Among("avane", 579, 97),
        Among("evane", 579, 96),
        Among("ivane", 579, 98),
        Among("uvane", 579, 99),
        Among("ačane", 579, 102),
        Among("acene", 578, 124),
        Among("lucene", 578, 121),
        Among("ačene", 578, 101),
        Among("lučene", 578, 117),
        Among("ine", 578, 11),
        Among("cine", 595, 137),
        Among("anine", 595, 10),
        Among("čine", 595, 89),
        Among("one", 578, 12),
        Among("are", -1, 53),
        Among("dre", -1, 54),
        Among("ere", -1, 55),
        Among("ore", -1, 56),
        Among("ase", -1, 161),
        Among("base", 604, 135),
        Among("acase", 604, 128),
        Among("gase", 604, 131),
        Among("jase", 604, 129),
        Among("astajase", 608, 138),
        Among("istajase", 608, 139),
        Among("ostajase", 608, 140),
        Among("injase", 608, 150),
        Among("kase", 604, 133),
        Among("nase", 604, 132),
        Among("irase", 604, 155),
        Among("urase", 604, 156),
        Among("tase", 604, 130),
        Among("vase", 604, 134),
        Among("avase", 618, 144),
        Among("evase", 618, 145),
        Among("ivase", 618, 146),
        Among("ovase", 618, 148),
        Among("uvase", 618, 147),
        Among("ese", -1, 57),
        Among("ise", -1, 58),
        Among("acise", 625, 124),
        Among("lucise", 625, 121),
        Among("rosise", 625, 127),
        Among("jetise", 625, 149),
        Among("ose", -1, 123),
        Among("astadose", 630, 141),
        Among("istadose", 630, 142),
        Among("ostadose", 630, 143),
        Among("ate", -1, 104),
        Among("acate", 634, 128),
        Among("ikate", 634, 68),
        Among("late", 634, 69),
        Among("irate", 634, 100),
        Among("urate", 634, 105),
        Among("tate", 634, 113),
        Among("avate", 634, 97),
        Among("evate", 634, 96),
        Among("ivate", 634, 98),
        Among("uvate", 634, 99),
        Among("ačate", 634, 102),
        Among("ete", -1, 70),
        Among("astadete", 646, 110),
        Among("istadete", 646, 111),
        Among("ostadete", 646, 112),
        Among("astajete", 646, 106),
        Among("istajete", 646, 107),
        Among("ostajete", 646, 108),
        Among("ijete", 646, 116),
        Among("injete", 646, 114),
        Among("ujete", 646, 25),
        Among("lucujete", 655, 121),
        Among("irujete", 655, 100),
        Among("lučujete", 655, 117),
        Among("nete", 646, 13),
        Among("astanete", 659, 110),
        Among("istanete", 659, 111),
        Among("ostanete", 659, 112),
        Among("astete", 646, 115),
        Among("ite", -1, 116),
        Among("acite", 664, 124),
        Among("lucite", 664, 121),
        Among("nite", 664, 13),
        Among("astanite", 667, 110),
        Among("istanite", 667, 111),
        Among("ostanite", 667, 112),
        Among("rosite", 664, 127),
        Among("jetite", 664, 118),
        Among("astite", 664, 115),
        Among("evite", 664, 92),
        Among("ovite", 664, 93),
        Among("ačite", 664, 101),
        Among("lučite", 664, 117),
        Among("rošite", 664, 90),
        Among("ajte", -1, 104),
        Among("urajte", 679, 105),
        Among("tajte", 679, 113),
        Among("astajte", 681, 106),
        Among("istajte", 681, 107),
        Among("ostajte", 681, 108),
        Among("avajte", 679, 97),
        Among("evajte", 679, 96),
        Among("ivajte", 679, 98),
        Among("uvajte", 679, 99),
        Among("ijte", -1, 116),
        Among("lucujte", -1, 121),
        Among("irujte", -1, 100),
        Among("lučujte", -1, 117),
        Among("aste", -1, 94),
        Among("acaste", 693, 128),
        Among("astajaste", 693, 106),
        Among("istajaste", 693, 107),
        Among("ostajaste", 693, 108),
        Among("injaste", 693, 114),
        Among("iraste", 693, 100),
        Among("uraste", 693, 105),
        Among("taste", 693, 113),
        Among("avaste", 693, 97),
        Among("evaste", 693, 96),
        Among("ivaste", 693, 98),
        Among("ovaste", 693, 76),
        Among("uvaste", 693, 99),
        Among("ačaste", 693, 102),
        Among("este", -1, 71),
        Among("iste", -1, 72),
        Among("aciste", 709, 124),
        Among("luciste", 709, 121),
        Among("niste", 709, 103),
        Among("rosiste", 709, 127),
        Among("jetiste", 709, 118),
        Among("ačiste", 709, 101),
        Among("lučiste", 709, 117),
        Among("rošiste", 709, 90),
        Among("kste", -1, 73),
        Among("oste", -1, 74),
        Among("astadoste", 719, 110),
        Among("istadoste", 719, 111),
        Among("ostadoste", 719, 112),
        Among("nuste", -1, 13),
        Among("ište", -1, 75),
        Among("ave", -1, 77),
        Among("eve", -1, 78),
        Among("ajeve", 726, 109),
        Among("cajeve", 727, 26),
        Among("lajeve", 727, 30),
        Among("rajeve", 727, 31),
        Among("ćajeve", 727, 28),
        Among("čajeve", 727, 27),
        Among("đajeve", 727, 29),
        Among("ive", -1, 79),
        Among("ove", -1, 80),
        Among("gove", 735, 20),
        Among("ugove", 736, 17),
        Among("love", 735, 82),
        Among("olove", 738, 49),
        Among("move", 735, 81),
        Among("onove", 735, 12),
        Among("aće", -1, 14),
        Among("eće", -1, 15),
        Among("uće", -1, 16),
        Among("ače", -1, 101),
        Among("luče", -1, 117),
        Among("aše", -1, 104),
        Among("baše", 747, 63),
        Among("gaše", 747, 64),
        Among("jaše", 747, 61),
        Among("astajaše", 750, 106),
        Among("istajaše", 750, 107),
        Among("ostajaše", 750, 108),
        Among("injaše", 750, 114),
        Among("kaše", 747, 62),
        Among("naše", 747, 60),
        Among("iraše", 747, 100),
        Among("uraše", 747, 105),
        Among("taše", 747, 59),
        Among("vaše", 747, 65),
        Among("avaše", 760, 97),
        Among("evaše", 760, 96),
        Among("ivaše", 760, 98),
        Among("ovaše", 760, 76),
        Among("uvaše", 760, 99),
        Among("ačaše", 747, 102),
        Among("eše", -1, 66),
        Among("iše", -1, 67),
        Among("jetiše", 768, 118),
        Among("ačiše", 768, 101),
        Among("lučiše", 768, 117),
        Among("rošiše", 768, 90),
        Among("oše", -1, 91),
        Among("astadoše", 773, 110),
        Among("istadoše", 773, 111),
        Among("ostadoše", 773, 112),
        Among("aceg", -1, 124),
        Among("eceg", -1, 125),
        Among("uceg", -1, 126),
        Among("anjijeg", -1, 84),
        Among("enjijeg", -1, 85),
        Among("snjijeg", -1, 122),
        Among("šnjijeg", -1, 86),
        Among("kijeg", -1, 95),
        Among("skijeg", 784, 1),
        Among("škijeg", 784, 2),
        Among("elijeg", -1, 83),
        Among("nijeg", -1, 13),
        Among("osijeg", -1, 123),
        Among("atijeg", -1, 120),
        Among("evitijeg", -1, 92),
        Among("ovitijeg", -1, 93),
        Among("astijeg", -1, 94),
        Among("avijeg", -1, 77),
        Among("evijeg", -1, 78),
        Among("ivijeg", -1, 79),
        Among("ovijeg", -1, 80),
        Among("ošijeg", -1, 91),
        Among("anjeg", -1, 84),
        Among("enjeg", -1, 85),
        Among("snjeg", -1, 122),
        Among("šnjeg", -1, 86),
        Among("keg", -1, 95),
        Among("eleg", -1, 83),
        Among("neg", -1, 13),
        Among("aneg", 805, 10),
        Among("eneg", 805, 87),
        Among("sneg", 805, 159),
        Among("šneg", 805, 88),
        Among("oseg", -1, 123),
        Among("ateg", -1, 120),
        Among("aveg", -1, 77),
        Among("eveg", -1, 78),
        Among("iveg", -1, 79),
        Among("oveg", -1, 80),
        Among("aćeg", -1, 14),
        Among("ećeg", -1, 15),
        Among("ućeg", -1, 16),
        Among("ošeg", -1, 91),
        Among("acog", -1, 124),
        Among("ecog", -1, 125),
        Among("ucog", -1, 126),
        Among("anjog", -1, 84),
        Among("enjog", -1, 85),
        Among("snjog", -1, 122),
        Among("šnjog", -1, 86),
        Among("kog", -1, 95),
        Among("skog", 827, 1),
        Among("škog", 827, 2),
        Among("elog", -1, 83),
        Among("nog", -1, 13),
        Among("cinog", 831, 137),
        Among("činog", 831, 89),
        Among("osog", -1, 123),
        Among("atog", -1, 120),
        Among("evitog", -1, 92),
        Among("ovitog", -1, 93),
        Among("astog", -1, 94),
        Among("avog", -1, 77),
        Among("evog", -1, 78),
        Among("ivog", -1, 79),
        Among("ovog", -1, 80),
        Among("aćog", -1, 14),
        Among("ećog", -1, 15),
        Among("ućog", -1, 16),
        Among("ošog", -1, 91),
        Among("ah", -1, 104),
        Among("acah", 847, 128),
        Among("astajah", 847, 106),
        Among("istajah", 847, 107),
        Among("ostajah", 847, 108),
        Among("injah", 847, 114),
        Among("irah", 847, 100),
        Among("urah", 847, 105),
        Among("tah", 847, 113),
        Among("avah", 847, 97),
        Among("evah", 847, 96),
        Among("ivah", 847, 98),
        Among("ovah", 847, 76),
        Among("uvah", 847, 99),
        Among("ačah", 847, 102),
        Among("ih", -1, 116),
        Among("acih", 862, 124),
        Among("ecih", 862, 125),
        Among("ucih", 862, 126),
        Among("lucih", 865, 121),
        Among("anjijih", 862, 84),
        Among("enjijih", 862, 85),
        Among("snjijih", 862, 122),
        Among("šnjijih", 862, 86),
        Among("kijih", 862, 95),
        Among("skijih", 871, 1),
        Among("škijih", 871, 2),
        Among("elijih", 862, 83),
        Among("nijih", 862, 13),
        Among("osijih", 862, 123),
        Among("atijih", 862, 120),
        Among("evitijih", 862, 92),
        Among("ovitijih", 862, 93),
        Among("astijih", 862, 94),
        Among("avijih", 862, 77),
        Among("evijih", 862, 78),
        Among("ivijih", 862, 79),
        Among("ovijih", 862, 80),
        Among("ošijih", 862, 91),
        Among("anjih", 862, 84),
        Among("enjih", 862, 85),
        Among("snjih", 862, 122),
        Among("šnjih", 862, 86),
        Among("kih", 862, 95),
        Among("skih", 890, 1),
        Among("ških", 890, 2),
        Among("elih", 862, 83),
        Among("nih", 862, 13),
        Among("cinih", 894, 137),
        Among("činih", 894, 89),
        Among("osih", 862, 123),
        Among("rosih", 897, 127),
        Among("atih", 862, 120),
        Among("jetih", 862, 118),
        Among("evitih", 862, 92),
        Among("ovitih", 862, 93),
        Among("astih", 862, 94),
        Among("avih", 862, 77),
        Among("evih", 862, 78),
        Among("ivih", 862, 79),
        Among("ovih", 862, 80),
        Among("aćih", 862, 14),
        Among("ećih", 862, 15),
        Among("ućih", 862, 16),
        Among("ačih", 862, 101),
        Among("lučih", 862, 117),
        Among("oših", 862, 91),
        Among("roših", 913, 90),
        Among("astadoh", -1, 110),
        Among("istadoh", -1, 111),
        Among("ostadoh", -1, 112),
        Among("acuh", -1, 124),
        Among("ecuh", -1, 125),
        Among("ucuh", -1, 126),
        Among("aćuh", -1, 14),
        Among("ećuh", -1, 15),
        Among("ućuh", -1, 16),
        Among("aci", -1, 124),
        Among("aceci", -1, 124),
        Among("ieci", -1, 162),
        Among("ajuci", -1, 161),
        Among("irajuci", 927, 155),
        Among("urajuci", 927, 156),
        Among("astajuci", 927, 138),
        Among("istajuci", 927, 139),
        Among("ostajuci", 927, 140),
        Among("avajuci", 927, 144),
        Among("evajuci", 927, 145),
        Among("ivajuci", 927, 146),
        Among("uvajuci", 927, 147),
        Among("ujuci", -1, 157),
        Among("lucujuci", 937, 121),
        Among("irujuci", 937, 155),
        Among("luci", -1, 121),
        Among("nuci", -1, 164),
        Among("etuci", -1, 153),
        Among("astuci", -1, 136),
        Among("gi", -1, 20),
        Among("ugi", 944, 18),
        Among("aji", -1, 109),
        Among("caji", 946, 26),
        Among("laji", 946, 30),
        Among("raji", 946, 31),
        Among("ćaji", 946, 28),
        Among("čaji", 946, 27),
        Among("đaji", 946, 29),
        Among("biji", -1, 32),
        Among("ciji", -1, 33),
        Among("diji", -1, 34),
        Among("fiji", -1, 40),
        Among("giji", -1, 39),
        Among("anjiji", -1, 84),
        Among("enjiji", -1, 85),
        Among("snjiji", -1, 122),
        Among("šnjiji", -1, 86),
        Among("kiji", -1, 95),
        Among("skiji", 962, 1),
        Among("škiji", 962, 2),
        Among("liji", -1, 35),
        Among("eliji", 965, 83),
        Among("miji", -1, 37),
        Among("niji", -1, 13),
        Among("ganiji", 968, 9),
        Among("maniji", 968, 6),
        Among("paniji", 968, 7),
        Among("raniji", 968, 8),
        Among("taniji", 968, 5),
        Among("piji", -1, 41),
        Among("riji", -1, 42),
        Among("siji", -1, 43),
        Among("osiji", 976, 123),
        Among("tiji", -1, 44),
        Among("atiji", 978, 120),
        Among("evitiji", 978, 92),
        Among("ovitiji", 978, 93),
        Among("astiji", 978, 94),
        Among("aviji", -1, 77),
        Among("eviji", -1, 78),
        Among("iviji", -1, 79),
        Among("oviji", -1, 80),
        Among("ziji", -1, 45),
        Among("ošiji", -1, 91),
        Among("žiji", -1, 38),
        Among("anji", -1, 84),
        Among("enji", -1, 85),
        Among("snji", -1, 122),
        Among("šnji", -1, 86),
        Among("ki", -1, 95),
        Among("ski", 994, 1),
        Among("ški", 994, 2),
        Among("ali", -1, 104),
        Among("acali", 997, 128),
        Among("astajali", 997, 106),
        Among("istajali", 997, 107),
        Among("ostajali", 997, 108),
        Among("ijali", 997, 47),
        Among("injali", 997, 114),
        Among("nali", 997, 46),
        Among("irali", 997, 100),
        Among("urali", 997, 105),
        Among("tali", 997, 113),
        Among("astali", 1007, 110),
        Among("istali", 1007, 111),
        Among("ostali", 1007, 112),
        Among("avali", 997, 97),
        Among("evali", 997, 96),
        Among("ivali", 997, 98),
        Among("ovali", 997, 76),
        Among("uvali", 997, 99),
        Among("ačali", 997, 102),
        Among("eli", -1, 83),
        Among("ili", -1, 116),
        Among("acili", 1018, 124),
        Among("lucili", 1018, 121),
        Among("nili", 1018, 103),
        Among("rosili", 1018, 127),
        Among("jetili", 1018, 118),
        Among("ozili", 1018, 48),
        Among("ačili", 1018, 101),
        Among("lučili", 1018, 117),
        Among("rošili", 1018, 90),
        Among("oli", -1, 50),
        Among("asli", -1, 115),
        Among("nuli", -1, 13),
        Among("rami", -1, 52),
        Among("lemi", -1, 51),
        Among("ni", -1, 13),
        Among("ani", 1033, 10),
        Among("acani", 1034, 128),
        Among("urani", 1034, 105),
        Among("tani", 1034, 113),
        Among("avani", 1034, 97),
        Among("evani", 1034, 96),
        Among("ivani", 1034, 98),
        Among("uvani", 1034, 99),
        Among("ačani", 1034, 102),
        Among("aceni", 1033, 124),
        Among("luceni", 1033, 121),
        Among("ačeni", 1033, 101),
        Among("lučeni", 1033, 117),
        Among("ini", 1033, 11),
        Among("cini", 1047, 137),
        Among("čini", 1047, 89),
        Among("oni", 1033, 12),
        Among("ari", -1, 53),
        Among("dri", -1, 54),
        Among("eri", -1, 55),
        Among("ori", -1, 56),
        Among("basi", -1, 135),
        Among("gasi", -1, 131),
        Among("jasi", -1, 129),
        Among("kasi", -1, 133),
        Among("nasi", -1, 132),
        Among("tasi", -1, 130),
        Among("vasi", -1, 134),
        Among("esi", -1, 152),
        Among("isi", -1, 154),
        Among("osi", -1, 123),
        Among("avsi", -1, 161),
        Among("acavsi", 1065, 128),
        Among("iravsi", 1065, 155),
        Among("tavsi", 1065, 160),
        Among("etavsi", 1068, 153),
        Among("astavsi", 1068, 141),
        Among("istavsi", 1068, 142),
        Among("ostavsi", 1068, 143),
        Among("ivsi", -1, 162),
        Among("nivsi", 1073, 158),
        Among("rosivsi", 1073, 127),
        Among("nuvsi", -1, 164),
        Among("ati", -1, 104),
        Among("acati", 1077, 128),
        Among("astajati", 1077, 106),
        Among("istajati", 1077, 107),
        Among("ostajati", 1077, 108),
        Among("injati", 1077, 114),
        Among("ikati", 1077, 68),
        Among("lati", 1077, 69),
        Among("irati", 1077, 100),
        Among("urati", 1077, 105),
        Among("tati", 1077, 113),
        Among("astati", 1087, 110),
        Among("istati", 1087, 111),
        Among("ostati", 1087, 112),
        Among("avati", 1077, 97),
        Among("evati", 1077, 96),
        Among("ivati", 1077, 98),
        Among("ovati", 1077, 76),
        Among("uvati", 1077, 99),
        Among("ačati", 1077, 102),
        Among("eti", -1, 70),
        Among("iti", -1, 116),
        Among("aciti", 1098, 124),
        Among("luciti", 1098, 121),
        Among("niti", 1098, 103),
        Among("rositi", 1098, 127),
        Among("jetiti", 1098, 118),
        Among("eviti", 1098, 92),
        Among("oviti", 1098, 93),
        Among("ačiti", 1098, 101),
        Among("lučiti", 1098, 117),
        Among("rošiti", 1098, 90),
        Among("asti", -1, 94),
        Among("esti", -1, 71),
        Among("isti", -1, 72),
        Among("ksti", -1, 73),
        Among("osti", -1, 74),
        Among("nuti", -1, 13),
        Among("avi", -1, 77),
        Among("evi", -1, 78),
        Among("ajevi", 1116, 109),
        Among("cajevi", 1117, 26),
        Among("lajevi", 1117, 30),
        Among("rajevi", 1117, 31),
        Among("ćajevi", 1117, 28),
        Among("čajevi", 1117, 27),
        Among("đajevi", 1117, 29),
        Among("ivi", -1, 79),
        Among("ovi", -1, 80),
        Among("govi", 1125, 20),
        Among("ugovi", 1126, 17),
        Among("lovi", 1125, 82),
        Among("olovi", 1128, 49),
        Among("movi", 1125, 81),
        Among("onovi", 1125, 12),
        Among("ieći", -1, 116),
        Among("ačeći", -1, 101),
        Among("ajući", -1, 104),
        Among("irajući", 1134, 100),
        Among("urajući", 1134, 105),
        Among("astajući", 1134, 106),
        Among("istajući", 1134, 107),
        Among("ostajući", 1134, 108),
        Among("avajući", 1134, 97),
        Among("evajući", 1134, 96),
        Among("ivajući", 1134, 98),
        Among("uvajući", 1134, 99),
        Among("ujući", -1, 25),
        Among("irujući", 1144, 100),
        Among("lučujući", 1144, 117),
        Among("nući", -1, 13),
        Among("etući", -1, 70),
        Among("astući", -1, 115),
        Among("ači", -1, 101),
        Among("luči", -1, 117),
        Among("baši", -1, 63),
        Among("gaši", -1, 64),
        Among("jaši", -1, 61),
        Among("kaši", -1, 62),
        Among("naši", -1, 60),
        Among("taši", -1, 59),
        Among("vaši", -1, 65),
        Among("eši", -1, 66),
        Among("iši", -1, 67),
        Among("oši", -1, 91),
        Among("avši", -1, 104),
        Among("iravši", 1162, 100),
        Among("tavši", 1162, 113),
        Among("etavši", 1164, 70),
        Among("astavši", 1164, 110),
        Among("istavši", 1164, 111),
        Among("ostavši", 1164, 112),
        Among("ačavši", 1162, 102),
        Among("ivši", -1, 116),
        Among("nivši", 1170, 103),
        Among("rošivši", 1170, 90),
        Among("nuvši", -1, 13),
        Among("aj", -1, 104),
        Among("uraj", 1174, 105),
        Among("taj", 1174, 113),
        Among("avaj", 1174, 97),
        Among("evaj", 1174, 96),
        Among("ivaj", 1174, 98),
        Among("uvaj", 1174, 99),
        Among("ij", -1, 116),
        Among("acoj", -1, 124),
        Among("ecoj", -1, 125),
        Among("ucoj", -1, 126),
        Among("anjijoj", -1, 84),
        Among("enjijoj", -1, 85),
        Among("snjijoj", -1, 122),
        Among("šnjijoj", -1, 86),
        Among("kijoj", -1, 95),
        Among("skijoj", 1189, 1),
        Among("škijoj", 1189, 2),
        Among("elijoj", -1, 83),
        Among("nijoj", -1, 13),
        Among("osijoj", -1, 123),
        Among("evitijoj", -1, 92),
        Among("ovitijoj", -1, 93),
        Among("astijoj", -1, 94),
        Among("avijoj", -1, 77),
        Among("evijoj", -1, 78),
        Among("ivijoj", -1, 79),
        Among("ovijoj", -1, 80),
        Among("ošijoj", -1, 91),
        Among("anjoj", -1, 84),
        Among("enjoj", -1, 85),
        Among("snjoj", -1, 122),
        Among("šnjoj", -1, 86),
        Among("koj", -1, 95),
        Among("skoj", 1207, 1),
        Among("škoj", 1207, 2),
        Among("aloj", -1, 104),
        Among("eloj", -1, 83),
        Among("noj", -1, 13),
        Among("cinoj", 1212, 137),
        Among("činoj", 1212, 89),
        Among("osoj", -1, 123),
        Among("atoj", -1, 120),
        Among("evitoj", -1, 92),
        Among("ovitoj", -1, 93),
        Among("astoj", -1, 94),
        Among("avoj", -1, 77),
        Among("evoj", -1, 78),
        Among("ivoj", -1, 79),
        Among("ovoj", -1, 80),
        Among("aćoj", -1, 14),
        Among("ećoj", -1, 15),
        Among("ućoj", -1, 16),
        Among("ošoj", -1, 91),
        Among("lucuj", -1, 121),
        Among("iruj", -1, 100),
        Among("lučuj", -1, 117),
        Among("al", -1, 104),
        Among("iral", 1231, 100),
        Among("ural", 1231, 105),
        Among("el", -1, 119),
        Among("il", -1, 116),
        Among("am", -1, 104),
        Among("acam", 1236, 128),
        Among("iram", 1236, 100),
        Among("uram", 1236, 105),
        Among("tam", 1236, 113),
        Among("avam", 1236, 97),
        Among("evam", 1236, 96),
        Among("ivam", 1236, 98),
        Among("uvam", 1236, 99),
        Among("ačam", 1236, 102),
        Among("em", -1, 119),
        Among("acem", 1246, 124),
        Among("ecem", 1246, 125),
        Among("ucem", 1246, 126),
        Among("astadem", 1246, 110),
        Among("istadem", 1246, 111),
        Among("ostadem", 1246, 112),
        Among("ajem", 1246, 104),
        Among("cajem", 1253, 26),
        Among("lajem", 1253, 30),
        Among("rajem", 1253, 31),
        Among("astajem", 1253, 106),
        Among("istajem", 1253, 107),
        Among("ostajem", 1253, 108),
        Among("ćajem", 1253, 28),
        Among("čajem", 1253, 27),
        Among("đajem", 1253, 29),
        Among("ijem", 1246, 116),
        Among("anjijem", 1263, 84),
        Among("enjijem", 1263, 85),
        Among("snjijem", 1263, 123),
        Among("šnjijem", 1263, 86),
        Among("kijem", 1263, 95),
        Among("skijem", 1268, 1),
        Among("škijem", 1268, 2),
        Among("lijem", 1263, 24),
        Among("elijem", 1271, 83),
        Among("nijem", 1263, 13),
        Among("rarijem", 1263, 21),
        Among("sijem", 1263, 23),
        Among("osijem", 1275, 123),
        Among("atijem", 1263, 120),
        Among("evitijem", 1263, 92),
        Among("ovitijem", 1263, 93),
        Among("otijem", 1263, 22),
        Among("astijem", 1263, 94),
        Among("avijem", 1263, 77),
        Among("evijem", 1263, 78),
        Among("ivijem", 1263, 79),
        Among("ovijem", 1263, 80),
        Among("ošijem", 1263, 91),
        Among("anjem", 1246, 84),
        Among("enjem", 1246, 85),
        Among("injem", 1246, 114),
        Among("snjem", 1246, 122),
        Among("šnjem", 1246, 86),
        Among("ujem", 1246, 25),
        Among("lucujem", 1292, 121),
        Among("irujem", 1292, 100),
        Among("lučujem", 1292, 117),
        Among("kem", 1246, 95),
        Among("skem", 1296, 1),
        Among("škem", 1296, 2),
        Among("elem", 1246, 83),
        Among("nem", 1246, 13),
        Among("anem", 1300, 10),
        Among("astanem", 1301, 110),
        Among("istanem", 1301, 111),
        Among("ostanem", 1301, 112),
        Among("enem", 1300, 87),
        Among("snem", 1300, 159),
        Among("šnem", 1300, 88),
        Among("basem", 1246, 135),
        Among("gasem", 1246, 131),
        Among("jasem", 1246, 129),
        Among("kasem", 1246, 133),
        Among("nasem", 1246, 132),
        Among("tasem", 1246, 130),
        Among("vasem", 1246, 134),
        Among("esem", 1246, 152),
        Among("isem", 1246, 154),
        Among("osem", 1246, 123),
        Among("atem", 1246, 120),
        Among("etem", 1246, 70),
        Among("evitem", 1246, 92),
        Among("ovitem", 1246, 93),
        Among("astem", 1246, 94),
        Among("istem", 1246, 151),
        Among("ištem", 1246, 75),
        Among("avem", 1246, 77),
        Among("evem", 1246, 78),
        Among("ivem", 1246, 79),
        Among("aćem", 1246, 14),
        Among("ećem", 1246, 15),
        Among("ućem", 1246, 16),
        Among("bašem", 1246, 63),
        Among("gašem", 1246, 64),
        Among("jašem", 1246, 61),
        Among("kašem", 1246, 62),
        Among("našem", 1246, 60),
        Among("tašem", 1246, 59),
        Among("vašem", 1246, 65),
        Among("ešem", 1246, 66),
        Among("išem", 1246, 67),
        Among("ošem", 1246, 91),
        Among("im", -1, 116),
        Among("acim", 1341, 124),
        Among("ecim", 1341, 125),
        Among("ucim", 1341, 126),
        Among("lucim", 1344, 121),
        Among("anjijim", 1341, 84),
        Among("enjijim", 1341, 85),
        Among("snjijim", 1341, 122),
        Among("šnjijim", 1341, 86),
        Among("kijim", 1341, 95),
        Among("skijim", 1350, 1),
        Among("škijim", 1350, 2),
        Among("elijim", 1341, 83),
        Among("nijim", 1341, 13),
        Among("osijim", 1341, 123),
        Among("atijim", 1341, 120),
        Among("evitijim", 1341, 92),
        Among("ovitijim", 1341, 93),
        Among("astijim", 1341, 94),
        Among("avijim", 1341, 77),
        Among("evijim", 1341, 78),
        Among("ivijim", 1341, 79),
        Among("ovijim", 1341, 80),
        Among("ošijim", 1341, 91),
        Among("anjim", 1341, 84),
        Among("enjim", 1341, 85),
        Among("snjim", 1341, 122),
        Among("šnjim", 1341, 86),
        Among("kim", 1341, 95),
        Among("skim", 1369, 1),
        Among("škim", 1369, 2),
        Among("elim", 1341, 83),
        Among("nim", 1341, 13),
        Among("cinim", 1373, 137),
        Among("činim", 1373, 89),
        Among("osim", 1341, 123),
        Among("rosim", 1376, 127),
        Among("atim", 1341, 120),
        Among("jetim", 1341, 118),
        Among("evitim", 1341, 92),
        Among("ovitim", 1341, 93),
        Among("astim", 1341, 94),
        Among("avim", 1341, 77),
        Among("evim", 1341, 78),
        Among("ivim", 1341, 79),
        Among("ovim", 1341, 80),
        Among("aćim", 1341, 14),
        Among("ećim", 1341, 15),
        Among("ućim", 1341, 16),
        Among("ačim", 1341, 101),
        Among("lučim", 1341, 117),
        Among("ošim", 1341, 91),
        Among("rošim", 1392, 90),
        Among("acom", -1, 124),
        Among("ecom", -1, 125),
        Among("ucom", -1, 126),
        Among("gom", -1, 20),
        Among("logom", 1397, 19),
        Among("ugom", 1397, 18),
        Among("bijom", -1, 32),
        Among("cijom", -1, 33),
        Among("dijom", -1, 34),
        Among("fijom", -1, 40),
        Among("gijom", -1, 39),
        Among("lijom", -1, 35),
        Among("mijom", -1, 37),
        Among("nijom", -1, 36),
        Among("ganijom", 1407, 9),
        Among("manijom", 1407, 6),
        Among("panijom", 1407, 7),
        Among("ranijom", 1407, 8),
        Among("tanijom", 1407, 5),
        Among("pijom", -1, 41),
        Among("rijom", -1, 42),
        Among("sijom", -1, 43),
        Among("tijom", -1, 44),
        Among("zijom", -1, 45),
        Among("žijom", -1, 38),
        Among("anjom", -1, 84),
        Among("enjom", -1, 85),
        Among("snjom", -1, 122),
        Among("šnjom", -1, 86),
        Among("kom", -1, 95),
        Among("skom", 1423, 1),
        Among("škom", 1423, 2),
        Among("alom", -1, 104),
        Among("ijalom", 1426, 47),
        Among("nalom", 1426, 46),
        Among("elom", -1, 83),
        Among("ilom", -1, 116),
        Among("ozilom", 1430, 48),
        Among("olom", -1, 50),
        Among("ramom", -1, 52),
        Among("lemom", -1, 51),
        Among("nom", -1, 13),
        Among("anom", 1435, 10),
        Among("inom", 1435, 11),
        Among("cinom", 1437, 137),
        Among("aninom", 1437, 10),
        Among("činom", 1437, 89),
        Among("onom", 1435, 12),
        Among("arom", -1, 53),
        Among("drom", -1, 54),
        Among("erom", -1, 55),
        Among("orom", -1, 56),
        Among("basom", -1, 135),
        Among("gasom", -1, 131),
        Among("jasom", -1, 129),
        Among("kasom", -1, 133),
        Among("nasom", -1, 132),
        Among("tasom", -1, 130),
        Among("vasom", -1, 134),
        Among("esom", -1, 57),
        Among("isom", -1, 58),
        Among("osom", -1, 123),
        Among("atom", -1, 120),
        Among("ikatom", 1456, 68),
        Among("latom", 1456, 69),
        Among("etom", -1, 70),
        Among("evitom", -1, 92),
        Among("ovitom", -1, 93),
        Among("astom", -1, 94),
        Among("estom", -1, 71),
        Among("istom", -1, 72),
        Among("kstom", -1, 73),
        Among("ostom", -1, 74),
        Among("avom", -1, 77),
        Among("evom", -1, 78),
        Among("ivom", -1, 79),
        Among("ovom", -1, 80),
        Among("lovom", 1470, 82),
        Among("movom", 1470, 81),
        Among("stvom", -1, 3),
        Among("štvom", -1, 4),
        Among("aćom", -1, 14),
        Among("ećom", -1, 15),
        Among("ućom", -1, 16),
        Among("bašom", -1, 63),
        Among("gašom", -1, 64),
        Among("jašom", -1, 61),
        Among("kašom", -1, 62),
        Among("našom", -1, 60),
        Among("tašom", -1, 59),
        Among("vašom", -1, 65),
        Among("ešom", -1, 66),
        Among("išom", -1, 67),
        Among("ošom", -1, 91),
        Among("an", -1, 104),
        Among("acan", 1488, 128),
        Among("iran", 1488, 100),
        Among("uran", 1488, 105),
        Among("tan", 1488, 113),
        Among("avan", 1488, 97),
        Among("evan", 1488, 96),
        Among("ivan", 1488, 98),
        Among("uvan", 1488, 99),
        Among("ačan", 1488, 102),
        Among("acen", -1, 124),
        Among("lucen", -1, 121),
        Among("ačen", -1, 101),
        Among("lučen", -1, 117),
        Among("anin", -1, 10),
        Among("ao", -1, 104),
        Among("acao", 1503, 128),
        Among("astajao", 1503, 106),
        Among("istajao", 1503, 107),
        Among("ostajao", 1503, 108),
        Among("injao", 1503, 114),
        Among("irao", 1503, 100),
        Among("urao", 1503, 105),
        Among("tao", 1503, 113),
        Among("astao", 1511, 110),
        Among("istao", 1511, 111),
        Among("ostao", 1511, 112),
        Among("avao", 1503, 97),
        Among("evao", 1503, 96),
        Among("ivao", 1503, 98),
        Among("ovao", 1503, 76),
        Among("uvao", 1503, 99),
        Among("ačao", 1503, 102),
        Among("go", -1, 20),
        Among("ugo", 1521, 18),
        Among("io", -1, 116),
        Among("acio", 1523, 124),
        Among("lucio", 1523, 121),
        Among("lio", 1523, 24),
        Among("nio", 1523, 103),
        Among("rario", 1523, 21),
        Among("sio", 1523, 23),
        Among("rosio", 1529, 127),
        Among("jetio", 1523, 118),
        Among("otio", 1523, 22),
        Among("ačio", 1523, 101),
        Among("lučio", 1523, 117),
        Among("rošio", 1523, 90),
        Among("bijo", -1, 32),
        Among("cijo", -1, 33),
        Among("dijo", -1, 34),
        Among("fijo", -1, 40),
        Among("gijo", -1, 39),
        Among("lijo", -1, 35),
        Among("mijo", -1, 37),
        Among("nijo", -1, 36),
        Among("pijo", -1, 41),
        Among("rijo", -1, 42),
        Among("sijo", -1, 43),
        Among("tijo", -1, 44),
        Among("zijo", -1, 45),
        Among("žijo", -1, 38),
        Among("anjo", -1, 84),
        Among("enjo", -1, 85),
        Among("snjo", -1, 122),
        Among("šnjo", -1, 86),
        Among("ko", -1, 95),
        Among("sko", 1554, 1),
        Among("ško", 1554, 2),
        Among("alo", -1, 104),
        Among("acalo", 1557, 128),
        Among("astajalo", 1557, 106),
        Among("istajalo", 1557, 107),
        Among("ostajalo", 1557, 108),
        Among("ijalo", 1557, 47),
        Among("injalo", 1557, 114),
        Among("nalo", 1557, 46),
        Among("iralo", 1557, 100),
        Among("uralo", 1557, 105),
        Among("talo", 1557, 113),
        Among("astalo", 1567, 110),
        Among("istalo", 1567, 111),
        Among("ostalo", 1567, 112),
        Among("avalo", 1557, 97),
        Among("evalo", 1557, 96),
        Among("ivalo", 1557, 98),
        Among("ovalo", 1557, 76),
        Among("uvalo", 1557, 99),
        Among("ačalo", 1557, 102),
        Among("elo", -1, 83),
        Among("ilo", -1, 116),
        Among("acilo", 1578, 124),
        Among("lucilo", 1578, 121),
        Among("nilo", 1578, 103),
        Among("rosilo", 1578, 127),
        Among("jetilo", 1578, 118),
        Among("ačilo", 1578, 101),
        Among("lučilo", 1578, 117),
        Among("rošilo", 1578, 90),
        Among("aslo", -1, 115),
        Among("nulo", -1, 13),
        Among("amo", -1, 104),
        Among("acamo", 1589, 128),
        Among("ramo", 1589, 52),
        Among("iramo", 1591, 100),
        Among("uramo", 1591, 105),
        Among("tamo", 1589, 113),
        Among("avamo", 1589, 97),
        Among("evamo", 1589, 96),
        Among("ivamo", 1589, 98),
        Among("uvamo", 1589, 99),
        Among("ačamo", 1589, 102),
        Among("emo", -1, 119),
        Among("astademo", 1600, 110),
        Among("istademo", 1600, 111),
        Among("ostademo", 1600, 112),
        Among("astajemo", 1600, 106),
        Among("istajemo", 1600, 107),
        Among("ostajemo", 1600, 108),
        Among("ijemo", 1600, 116),
        Among("injemo", 1600, 114),
        Among("ujemo", 1600, 25),
        Among("lucujemo", 1609, 121),
        Among("irujemo", 1609, 100),
        Among("lučujemo", 1609, 117),
        Among("lemo", 1600, 51),
        Among("nemo", 1600, 13),
        Among("astanemo", 1614, 110),
        Among("istanemo", 1614, 111),
        Among("ostanemo", 1614, 112),
        Among("etemo", 1600, 70),
        Among("astemo", 1600, 115),
        Among("imo", -1, 116),
        Among("acimo", 1620, 124),
        Among("lucimo", 1620, 121),
        Among("nimo", 1620, 13),
        Among("astanimo", 1623, 110),
        Among("istanimo", 1623, 111),
        Among("ostanimo", 1623, 112),
        Among("rosimo", 1620, 127),
        Among("etimo", 1620, 70),
        Among("jetimo", 1628, 118),
        Among("astimo", 1620, 115),
        Among("ačimo", 1620, 101),
        Among("lučimo", 1620, 117),
        Among("rošimo", 1620, 90),
        Among("ajmo", -1, 104),
        Among("urajmo", 1634, 105),
        Among("tajmo", 1634, 113),
        Among("astajmo", 1636, 106),
        Among("istajmo", 1636, 107),
        Among("ostajmo", 1636, 108),
        Among("avajmo", 1634, 97),
        Among("evajmo", 1634, 96),
        Among("ivajmo", 1634, 98),
        Among("uvajmo", 1634, 99),
        Among("ijmo", -1, 116),
        Among("ujmo", -1, 25),
        Among("lucujmo", 1645, 121),
        Among("irujmo", 1645, 100),
        Among("lučujmo", 1645, 117),
        Among("asmo", -1, 104),
        Among("acasmo", 1649, 128),
        Among("astajasmo", 1649, 106),
        Among("istajasmo", 1649, 107),
        Among("ostajasmo", 1649, 108),
        Among("injasmo", 1649, 114),
        Among("irasmo", 1649, 100),
        Among("urasmo", 1649, 105),
        Among("tasmo", 1649, 113),
        Among("avasmo", 1649, 97),
        Among("evasmo", 1649, 96),
        Among("ivasmo", 1649, 98),
        Among("ovasmo", 1649, 76),
        Among("uvasmo", 1649, 99),
        Among("ačasmo", 1649, 102),
        Among("ismo", -1, 116),
        Among("acismo", 1664, 124),
        Among("lucismo", 1664, 121),
        Among("nismo", 1664, 103),
        Among("rosismo", 1664, 127),
        Among("jetismo", 1664, 118),
        Among("ačismo", 1664, 101),
        Among("lučismo", 1664, 117),
        Among("rošismo", 1664, 90),
        Among("astadosmo", -1, 110),
        Among("istadosmo", -1, 111),
        Among("ostadosmo", -1, 112),
        Among("nusmo", -1, 13),
        Among("no", -1, 13),
        Among("ano", 1677, 104),
        Among("acano", 1678, 128),
        Among("urano", 1678, 105),
        Among("tano", 1678, 113),
        Among("avano", 1678, 97),
        Among("evano", 1678, 96),
        Among("ivano", 1678, 98),
        Among("uvano", 1678, 99),
        Among("ačano", 1678, 102),
        Among("aceno", 1677, 124),
        Among("luceno", 1677, 121),
        Among("ačeno", 1677, 101),
        Among("lučeno", 1677, 117),
        Among("ino", 1677, 11),
        Among("cino", 1691, 137),
        Among("čino", 1691, 89),
        Among("ato", -1, 120),
        Among("ikato", 1694, 68),
        Among("lato", 1694, 69),
        Among("eto", -1, 70),
        Among("evito", -1, 92),
        Among("ovito", -1, 93),
        Among("asto", -1, 94),
        Among("esto", -1, 71),
        Among("isto", -1, 72),
        Among("ksto", -1, 73),
        Among("osto", -1, 74),
        Among("nuto", -1, 13),
        Among("nuo", -1, 13),
        Among("avo", -1, 77),
        Among("evo", -1, 78),
        Among("ivo", -1, 79),
        Among("ovo", -1, 80),
        Among("stvo", -1, 3),
        Among("štvo", -1, 4),
        Among("as", -1, 161),
        Among("acas", 1713, 128),
        Among("iras", 1713, 155),
        Among("uras", 1713, 156),
        Among("tas", 1713, 160),
        Among("avas", 1713, 144),
        Among("evas", 1713, 145),
        Among("ivas", 1713, 146),
        Among("uvas", 1713, 147),
        Among("es", -1, 163),
        Among("astades", 1722, 141),
        Among("istades", 1722, 142),
        Among("ostades", 1722, 143),
        Among("astajes", 1722, 138),
        Among("istajes", 1722, 139),
        Among("ostajes", 1722, 140),
        Among("ijes", 1722, 162),
        Among("injes", 1722, 150),
        Among("ujes", 1722, 157),
        Among("lucujes", 1731, 121),
        Among("irujes", 1731, 155),
        Among("nes", 1722, 164),
        Among("astanes", 1734, 141),
        Among("istanes", 1734, 142),
        Among("ostanes", 1734, 143),
        Among("etes", 1722, 153),
        Among("astes", 1722, 136),
        Among("is", -1, 162),
        Among("acis", 1740, 124),
        Among("lucis", 1740, 121),
        Among("nis", 1740, 158),
        Among("rosis", 1740, 127),
        Among("jetis", 1740, 149),
        Among("at", -1, 104),
        Among("acat", 1746, 128),
        Among("astajat", 1746, 106),
        Among("istajat", 1746, 107),
        Among("ostajat", 1746, 108),
        Among("injat", 1746, 114),
        Among("irat", 1746, 100),
        Among("urat", 1746, 105),
        Among("tat", 1746, 113),
        Among("astat", 1754, 110),
        Among("istat", 1754, 111),
        Among("ostat", 1754, 112),
        Among("avat", 1746, 97),
        Among("evat", 1746, 96),
        Among("ivat", 1746, 98),
        Among("irivat", 1760, 100),
        Among("ovat", 1746, 76),
        Among("uvat", 1746, 99),
        Among("ačat", 1746, 102),
        Among("it", -1, 116),
        Among("acit", 1765, 124),
        Among("lucit", 1765, 121),
        Among("rosit", 1765, 127),
        Among("jetit", 1765, 118),
        Among("ačit", 1765, 101),
        Among("lučit", 1765, 117),
        Among("rošit", 1765, 90),
        Among("nut", -1, 13),
        Among("astadu", -1, 110),
        Among("istadu", -1, 111),
        Among("ostadu", -1, 112),
        Among("gu", -1, 20),
        Among("logu", 1777, 19),
        Among("ugu", 1777, 18),
        Among("ahu", -1, 104),
        Among("acahu", 1780, 128),
        Among("astajahu", 1780, 106),
        Among("istajahu", 1780, 107),
        Among("ostajahu", 1780, 108),
        Among("injahu", 1780, 114),
        Among("irahu", 1780, 100),
        Among("urahu", 1780, 105),
        Among("avahu", 1780, 97),
        Among("evahu", 1780, 96),
        Among("ivahu", 1780, 98),
        Among("ovahu", 1780, 76),
        Among("uvahu", 1780, 99),
        Among("ačahu", 1780, 102),
        Among("aju", -1, 104),
        Among("caju", 1794, 26),
        Among("acaju", 1795, 128),
        Among("laju", 1794, 30),
        Among("raju", 1794, 31),
        Among("iraju", 1798, 100),
        Among("uraju", 1798, 105),
        Among("taju", 1794, 113),
        Among("astaju", 1801, 106),
        Among("istaju", 1801, 107),
        Among("ostaju", 1801, 108),
        Among("avaju", 1794, 97),
        Among("evaju", 1794, 96),
        Among("ivaju", 1794, 98),
        Among("uvaju", 1794, 99),
        Among("ćaju", 1794, 28),
        Among("čaju", 1794, 27),
        Among("ačaju", 1810, 102),
        Among("đaju", 1794, 29),
        Among("iju", -1, 116),
        Among("biju", 1813, 32),
        Among("ciju", 1813, 33),
        Among("diju", 1813, 34),
        Among("fiju", 1813, 40),
        Among("giju", 1813, 39),
        Among("anjiju", 1813, 84),
        Among("enjiju", 1813, 85),
        Among("snjiju", 1813, 122),
        Among("šnjiju", 1813, 86),
        Among("kiju", 1813, 95),
        Among("liju", 1813, 24),
        Among("eliju", 1824, 83),
        Among("miju", 1813, 37),
        Among("niju", 1813, 13),
        Among("ganiju", 1827, 9),
        Among("maniju", 1827, 6),
        Among("paniju", 1827, 7),
        Among("raniju", 1827, 8),
        Among("taniju", 1827, 5),
        Among("piju", 1813, 41),
        Among("riju", 1813, 42),
        Among("rariju", 1834, 21),
        Among("siju", 1813, 23),
        Among("osiju", 1836, 123),
        Among("tiju", 1813, 44),
        Among("atiju", 1838, 120),
        Among("otiju", 1838, 22),
        Among("aviju", 1813, 77),
        Among("eviju", 1813, 78),
        Among("iviju", 1813, 79),
        Among("oviju", 1813, 80),
        Among("ziju", 1813, 45),
        Among("ošiju", 1813, 91),
        Among("žiju", 1813, 38),
        Among("anju", -1, 84),
        Among("enju", -1, 85),
        Among("snju", -1, 122),
        Among("šnju", -1, 86),
        Among("uju", -1, 25),
        Among("lucuju", 1852, 121),
        Among("iruju", 1852, 100),
        Among("lučuju", 1852, 117),
        Among("ku", -1, 95),
        Among("sku", 1856, 1),
        Among("šku", 1856, 2),
        Among("alu", -1, 104),
        Among("ijalu", 1859, 47),
        Among("nalu", 1859, 46),
        Among("elu", -1, 83),
        Among("ilu", -1, 116),
        Among("ozilu", 1863, 48),
        Among("olu", -1, 50),
        Among("ramu", -1, 52),
        Among("acemu", -1, 124),
        Among("ecemu", -1, 125),
        Among("ucemu", -1, 126),
        Among("anjijemu", -1, 84),
        Among("enjijemu", -1, 85),
        Among("snjijemu", -1, 122),
        Among("šnjijemu", -1, 86),
        Among("kijemu", -1, 95),
        Among("skijemu", 1874, 1),
        Among("škijemu", 1874, 2),
        Among("elijemu", -1, 83),
        Among("nijemu", -1, 13),
        Among("osijemu", -1, 123),
        Among("atijemu", -1, 120),
        Among("evitijemu", -1, 92),
        Among("ovitijemu", -1, 93),
        Among("astijemu", -1, 94),
        Among("avijemu", -1, 77),
        Among("evijemu", -1, 78),
        Among("ivijemu", -1, 79),
        Among("ovijemu", -1, 80),
        Among("ošijemu", -1, 91),
        Among("anjemu", -1, 84),
        Among("enjemu", -1, 85),
        Among("snjemu", -1, 122),
        Among("šnjemu", -1, 86),
        Among("kemu", -1, 95),
        Among("skemu", 1893, 1),
        Among("škemu", 1893, 2),
        Among("lemu", -1, 51),
        Among("elemu", 1896, 83),
        Among("nemu", -1, 13),
        Among("anemu", 1898, 10),
        Among("enemu", 1898, 87),
        Among("snemu", 1898, 159),
        Among("šnemu", 1898, 88),
        Among("osemu", -1, 123),
        Among("atemu", -1, 120),
        Among("evitemu", -1, 92),
        Among("ovitemu", -1, 93),
        Among("astemu", -1, 94),
        Among("avemu", -1, 77),
        Among("evemu", -1, 78),
        Among("ivemu", -1, 79),
        Among("ovemu", -1, 80),
        Among("aćemu", -1, 14),
        Among("ećemu", -1, 15),
        Among("ućemu", -1, 16),
        Among("ošemu", -1, 91),
        Among("acomu", -1, 124),
        Among("ecomu", -1, 125),
        Among("ucomu", -1, 126),
        Among("anjomu", -1, 84),
        Among("enjomu", -1, 85),
        Among("snjomu", -1, 122),
        Among("šnjomu", -1, 86),
        Among("komu", -1, 95),
        Among("skomu", 1923, 1),
        Among("škomu", 1923, 2),
        Among("elomu", -1, 83),
        Among("nomu", -1, 13),
        Among("cinomu", 1927, 137),
        Among("činomu", 1927, 89),
        Among("osomu", -1, 123),
        Among("atomu", -1, 120),
        Among("evitomu", -1, 92),
        Among("ovitomu", -1, 93),
        Among("astomu", -1, 94),
        Among("avomu", -1, 77),
        Among("evomu", -1, 78),
        Among("ivomu", -1, 79),
        Among("ovomu", -1, 80),
        Among("aćomu", -1, 14),
        Among("ećomu", -1, 15),
        Among("ućomu", -1, 16),
        Among("ošomu", -1, 91),
        Among("nu", -1, 13),
        Among("anu", 1943, 10),
        Among("astanu", 1944, 110),
        Among("istanu", 1944, 111),
        Among("ostanu", 1944, 112),
        Among("inu", 1943, 11),
        Among("cinu", 1948, 137),
        Among("aninu", 1948, 10),
        Among("činu", 1948, 89),
        Among("onu", 1943, 12),
        Among("aru", -1, 53),
        Among("dru", -1, 54),
        Among("eru", -1, 55),
        Among("oru", -1, 56),
        Among("basu", -1, 135),
        Among("gasu", -1, 131),
        Among("jasu", -1, 129),
        Among("kasu", -1, 133),
        Among("nasu", -1, 132),
        Among("tasu", -1, 130),
        Among("vasu", -1, 134),
        Among("esu", -1, 57),
        Among("isu", -1, 58),
        Among("osu", -1, 123),
        Among("atu", -1, 120),
        Among("ikatu", 1967, 68),
        Among("latu", 1967, 69),
        Among("etu", -1, 70),
        Among("evitu", -1, 92),
        Among("ovitu", -1, 93),
        Among("astu", -1, 94),
        Among("estu", -1, 71),
        Among("istu", -1, 72),
        Among("kstu", -1, 73),
        Among("ostu", -1, 74),
        Among("ištu", -1, 75),
        Among("avu", -1, 77),
        Among("evu", -1, 78),
        Among("ivu", -1, 79),
        Among("ovu", -1, 80),
        Among("lovu", 1982, 82),
        Among("movu", 1982, 81),
        Among("stvu", -1, 3),
        Among("štvu", -1, 4),
        Among("bašu", -1, 63),
        Among("gašu", -1, 64),
        Among("jašu", -1, 61),
        Among("kašu", -1, 62),
        Among("našu", -1, 60),
        Among("tašu", -1, 59),
        Among("vašu", -1, 65),
        Among("ešu", -1, 66),
        Among("išu", -1, 67),
        Among("ošu", -1, 91),
        Among("avav", -1, 97),
        Among("evav", -1, 96),
        Among("ivav", -1, 98),
        Among("uvav", -1, 99),
        Among("kov", -1, 95),
        Among("aš", -1, 104),
        Among("iraš", 2002, 100),
        Among("uraš", 2002, 105),
        Among("taš", 2002, 113),
        Among("avaš", 2002, 97),
        Among("evaš", 2002, 96),
        Among("ivaš", 2002, 98),
        Among("uvaš", 2002, 99),
        Among("ačaš", 2002, 102),
        Among("eš", -1, 119),
        Among("astadeš", 2011, 110),
        Among("istadeš", 2011, 111),
        Among("ostadeš", 2011, 112),
        Among("astaješ", 2011, 106),
        Among("istaješ", 2011, 107),
        Among("ostaješ", 2011, 108),
        Among("iješ", 2011, 116),
        Among("inješ", 2011, 114),
        Among("uješ", 2011, 25),
        Among("iruješ", 2020, 100),
        Among("lučuješ", 2020, 117),
        Among("neš", 2011, 13),
        Among("astaneš", 2023, 110),
        Among("istaneš", 2023, 111),
        Among("ostaneš", 2023, 112),
        Among("eteš", 2011, 70),
        Among("asteš", 2011, 115),
        Among("iš", -1, 116),
        Among("niš", 2029, 103),
        Among("jetiš", 2029, 118),
        Among("ačiš", 2029, 101),
        Among("lučiš", 2029, 117),
        Among("rošiš", 2029, 90)
    ]

    a_3 = [
        Among("a", -1, 1),
        Among("oga", 0, 1),
        Among("ama", 0, 1),
        Among("ima", 0, 1),
        Among("ena", 0, 1),
        Among("e", -1, 1),
        Among("og", -1, 1),
        Among("anog", 6, 1),
        Among("enog", 6, 1),
        Among("anih", -1, 1),
        Among("enih", -1, 1),
        Among("i", -1, 1),
        Among("ani", 11, 1),
        Among("eni", 11, 1),
        Among("anoj", -1, 1),
        Among("enoj", -1, 1),
        Among("anim", -1, 1),
        Among("enim", -1, 1),
        Among("om", -1, 1),
        Among("enom", 18, 1),
        Among("o", -1, 1),
        Among("ano", 20, 1),
        Among("eno", 20, 1),
        Among("ost", -1, 1),
        Among("u", -1, 1),
        Among("enu", 24, 1)
    ]


class lab0(BaseException): pass


class lab1(BaseException): pass


class lab2(BaseException): pass
