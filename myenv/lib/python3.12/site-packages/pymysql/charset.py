# Internal use only. Do not use directly.

MBLENGTH = {8: 1, 33: 3, 88: 2, 91: 2}


class Charset:
    def __init__(self, id, name, collation, is_default=False):
        self.id, self.name, self.collation = id, name, collation
        self.is_default = is_default

    def __repr__(self):
        return (
            f"Charset(id={self.id}, name={self.name!r}, collation={self.collation!r})"
        )

    @property
    def encoding(self):
        name = self.name
        if name in ("utf8mb4", "utf8mb3"):
            return "utf8"
        if name == "latin1":
            return "cp1252"
        if name == "koi8r":
            return "koi8_r"
        if name == "koi8u":
            return "koi8_u"
        return name

    @property
    def is_binary(self):
        return self.id == 63


class Charsets:
    def __init__(self):
        self._by_id = {}
        self._by_name = {}

    def add(self, c):
        self._by_id[c.id] = c
        if c.is_default:
            self._by_name[c.name] = c

    def by_id(self, id):
        return self._by_id[id]

    def by_name(self, name):
        if name == "utf8":
            name = "utf8mb4"
        return self._by_name.get(name.lower())


_charsets = Charsets()
charset_by_name = _charsets.by_name
charset_by_id = _charsets.by_id

"""
TODO: update this script.

Generated with:

mysql -N -s -e "select id, character_set_name, collation_name, is_default
from information_schema.collations order by id;" | python -c "import sys
for l in sys.stdin.readlines():
    id, name, collation, is_default  = l.split(chr(9))
    if is_default.strip() == "Yes":
        print('_charsets.add(Charset(%s, \'%s\', \'%s\', True))' \
              % (id, name, collation))
    else:
        print('_charsets.add(Charset(%s, \'%s\', \'%s\'))' \
              % (id, name, collation, bool(is_default.strip()))
"""

_charsets.add(Charset(1, "big5", "big5_chinese_ci", True))
_charsets.add(Charset(2, "latin2", "latin2_czech_cs"))
_charsets.add(Charset(3, "dec8", "dec8_swedish_ci", True))
_charsets.add(Charset(4, "cp850", "cp850_general_ci", True))
_charsets.add(Charset(5, "latin1", "latin1_german1_ci"))
_charsets.add(Charset(6, "hp8", "hp8_english_ci", True))
_charsets.add(Charset(7, "koi8r", "koi8r_general_ci", True))
_charsets.add(Charset(8, "latin1", "latin1_swedish_ci", True))
_charsets.add(Charset(9, "latin2", "latin2_general_ci", True))
_charsets.add(Charset(10, "swe7", "swe7_swedish_ci", True))
_charsets.add(Charset(11, "ascii", "ascii_general_ci", True))
_charsets.add(Charset(12, "ujis", "ujis_japanese_ci", True))
_charsets.add(Charset(13, "sjis", "sjis_japanese_ci", True))
_charsets.add(Charset(14, "cp1251", "cp1251_bulgarian_ci"))
_charsets.add(Charset(15, "latin1", "latin1_danish_ci"))
_charsets.add(Charset(16, "hebrew", "hebrew_general_ci", True))
_charsets.add(Charset(18, "tis620", "tis620_thai_ci", True))
_charsets.add(Charset(19, "euckr", "euckr_korean_ci", True))
_charsets.add(Charset(20, "latin7", "latin7_estonian_cs"))
_charsets.add(Charset(21, "latin2", "latin2_hungarian_ci"))
_charsets.add(Charset(22, "koi8u", "koi8u_general_ci", True))
_charsets.add(Charset(23, "cp1251", "cp1251_ukrainian_ci"))
_charsets.add(Charset(24, "gb2312", "gb2312_chinese_ci", True))
_charsets.add(Charset(25, "greek", "greek_general_ci", True))
_charsets.add(Charset(26, "cp1250", "cp1250_general_ci", True))
_charsets.add(Charset(27, "latin2", "latin2_croatian_ci"))
_charsets.add(Charset(28, "gbk", "gbk_chinese_ci", True))
_charsets.add(Charset(29, "cp1257", "cp1257_lithuanian_ci"))
_charsets.add(Charset(30, "latin5", "latin5_turkish_ci", True))
_charsets.add(Charset(31, "latin1", "latin1_german2_ci"))
_charsets.add(Charset(32, "armscii8", "armscii8_general_ci", True))
_charsets.add(Charset(33, "utf8mb3", "utf8mb3_general_ci", True))
_charsets.add(Charset(34, "cp1250", "cp1250_czech_cs"))
_charsets.add(Charset(36, "cp866", "cp866_general_ci", True))
_charsets.add(Charset(37, "keybcs2", "keybcs2_general_ci", True))
_charsets.add(Charset(38, "macce", "macce_general_ci", True))
_charsets.add(Charset(39, "macroman", "macroman_general_ci", True))
_charsets.add(Charset(40, "cp852", "cp852_general_ci", True))
_charsets.add(Charset(41, "latin7", "latin7_general_ci", True))
_charsets.add(Charset(42, "latin7", "latin7_general_cs"))
_charsets.add(Charset(43, "macce", "macce_bin"))
_charsets.add(Charset(44, "cp1250", "cp1250_croatian_ci"))
_charsets.add(Charset(45, "utf8mb4", "utf8mb4_general_ci", True))
_charsets.add(Charset(46, "utf8mb4", "utf8mb4_bin"))
_charsets.add(Charset(47, "latin1", "latin1_bin"))
_charsets.add(Charset(48, "latin1", "latin1_general_ci"))
_charsets.add(Charset(49, "latin1", "latin1_general_cs"))
_charsets.add(Charset(50, "cp1251", "cp1251_bin"))
_charsets.add(Charset(51, "cp1251", "cp1251_general_ci", True))
_charsets.add(Charset(52, "cp1251", "cp1251_general_cs"))
_charsets.add(Charset(53, "macroman", "macroman_bin"))
_charsets.add(Charset(57, "cp1256", "cp1256_general_ci", True))
_charsets.add(Charset(58, "cp1257", "cp1257_bin"))
_charsets.add(Charset(59, "cp1257", "cp1257_general_ci", True))
_charsets.add(Charset(63, "binary", "binary", True))
_charsets.add(Charset(64, "armscii8", "armscii8_bin"))
_charsets.add(Charset(65, "ascii", "ascii_bin"))
_charsets.add(Charset(66, "cp1250", "cp1250_bin"))
_charsets.add(Charset(67, "cp1256", "cp1256_bin"))
_charsets.add(Charset(68, "cp866", "cp866_bin"))
_charsets.add(Charset(69, "dec8", "dec8_bin"))
_charsets.add(Charset(70, "greek", "greek_bin"))
_charsets.add(Charset(71, "hebrew", "hebrew_bin"))
_charsets.add(Charset(72, "hp8", "hp8_bin"))
_charsets.add(Charset(73, "keybcs2", "keybcs2_bin"))
_charsets.add(Charset(74, "koi8r", "koi8r_bin"))
_charsets.add(Charset(75, "koi8u", "koi8u_bin"))
_charsets.add(Charset(76, "utf8mb3", "utf8mb3_tolower_ci"))
_charsets.add(Charset(77, "latin2", "latin2_bin"))
_charsets.add(Charset(78, "latin5", "latin5_bin"))
_charsets.add(Charset(79, "latin7", "latin7_bin"))
_charsets.add(Charset(80, "cp850", "cp850_bin"))
_charsets.add(Charset(81, "cp852", "cp852_bin"))
_charsets.add(Charset(82, "swe7", "swe7_bin"))
_charsets.add(Charset(83, "utf8mb3", "utf8mb3_bin"))
_charsets.add(Charset(84, "big5", "big5_bin"))
_charsets.add(Charset(85, "euckr", "euckr_bin"))
_charsets.add(Charset(86, "gb2312", "gb2312_bin"))
_charsets.add(Charset(87, "gbk", "gbk_bin"))
_charsets.add(Charset(88, "sjis", "sjis_bin"))
_charsets.add(Charset(89, "tis620", "tis620_bin"))
_charsets.add(Charset(91, "ujis", "ujis_bin"))
_charsets.add(Charset(92, "geostd8", "geostd8_general_ci", True))
_charsets.add(Charset(93, "geostd8", "geostd8_bin"))
_charsets.add(Charset(94, "latin1", "latin1_spanish_ci"))
_charsets.add(Charset(95, "cp932", "cp932_japanese_ci", True))
_charsets.add(Charset(96, "cp932", "cp932_bin"))
_charsets.add(Charset(97, "eucjpms", "eucjpms_japanese_ci", True))
_charsets.add(Charset(98, "eucjpms", "eucjpms_bin"))
_charsets.add(Charset(99, "cp1250", "cp1250_polish_ci"))
_charsets.add(Charset(192, "utf8mb3", "utf8mb3_unicode_ci"))
_charsets.add(Charset(193, "utf8mb3", "utf8mb3_icelandic_ci"))
_charsets.add(Charset(194, "utf8mb3", "utf8mb3_latvian_ci"))
_charsets.add(Charset(195, "utf8mb3", "utf8mb3_romanian_ci"))
_charsets.add(Charset(196, "utf8mb3", "utf8mb3_slovenian_ci"))
_charsets.add(Charset(197, "utf8mb3", "utf8mb3_polish_ci"))
_charsets.add(Charset(198, "utf8mb3", "utf8mb3_estonian_ci"))
_charsets.add(Charset(199, "utf8mb3", "utf8mb3_spanish_ci"))
_charsets.add(Charset(200, "utf8mb3", "utf8mb3_swedish_ci"))
_charsets.add(Charset(201, "utf8mb3", "utf8mb3_turkish_ci"))
_charsets.add(Charset(202, "utf8mb3", "utf8mb3_czech_ci"))
_charsets.add(Charset(203, "utf8mb3", "utf8mb3_danish_ci"))
_charsets.add(Charset(204, "utf8mb3", "utf8mb3_lithuanian_ci"))
_charsets.add(Charset(205, "utf8mb3", "utf8mb3_slovak_ci"))
_charsets.add(Charset(206, "utf8mb3", "utf8mb3_spanish2_ci"))
_charsets.add(Charset(207, "utf8mb3", "utf8mb3_roman_ci"))
_charsets.add(Charset(208, "utf8mb3", "utf8mb3_persian_ci"))
_charsets.add(Charset(209, "utf8mb3", "utf8mb3_esperanto_ci"))
_charsets.add(Charset(210, "utf8mb3", "utf8mb3_hungarian_ci"))
_charsets.add(Charset(211, "utf8mb3", "utf8mb3_sinhala_ci"))
_charsets.add(Charset(212, "utf8mb3", "utf8mb3_german2_ci"))
_charsets.add(Charset(213, "utf8mb3", "utf8mb3_croatian_ci"))
_charsets.add(Charset(214, "utf8mb3", "utf8mb3_unicode_520_ci"))
_charsets.add(Charset(215, "utf8mb3", "utf8mb3_vietnamese_ci"))
_charsets.add(Charset(223, "utf8mb3", "utf8mb3_general_mysql500_ci"))
_charsets.add(Charset(224, "utf8mb4", "utf8mb4_unicode_ci"))
_charsets.add(Charset(225, "utf8mb4", "utf8mb4_icelandic_ci"))
_charsets.add(Charset(226, "utf8mb4", "utf8mb4_latvian_ci"))
_charsets.add(Charset(227, "utf8mb4", "utf8mb4_romanian_ci"))
_charsets.add(Charset(228, "utf8mb4", "utf8mb4_slovenian_ci"))
_charsets.add(Charset(229, "utf8mb4", "utf8mb4_polish_ci"))
_charsets.add(Charset(230, "utf8mb4", "utf8mb4_estonian_ci"))
_charsets.add(Charset(231, "utf8mb4", "utf8mb4_spanish_ci"))
_charsets.add(Charset(232, "utf8mb4", "utf8mb4_swedish_ci"))
_charsets.add(Charset(233, "utf8mb4", "utf8mb4_turkish_ci"))
_charsets.add(Charset(234, "utf8mb4", "utf8mb4_czech_ci"))
_charsets.add(Charset(235, "utf8mb4", "utf8mb4_danish_ci"))
_charsets.add(Charset(236, "utf8mb4", "utf8mb4_lithuanian_ci"))
_charsets.add(Charset(237, "utf8mb4", "utf8mb4_slovak_ci"))
_charsets.add(Charset(238, "utf8mb4", "utf8mb4_spanish2_ci"))
_charsets.add(Charset(239, "utf8mb4", "utf8mb4_roman_ci"))
_charsets.add(Charset(240, "utf8mb4", "utf8mb4_persian_ci"))
_charsets.add(Charset(241, "utf8mb4", "utf8mb4_esperanto_ci"))
_charsets.add(Charset(242, "utf8mb4", "utf8mb4_hungarian_ci"))
_charsets.add(Charset(243, "utf8mb4", "utf8mb4_sinhala_ci"))
_charsets.add(Charset(244, "utf8mb4", "utf8mb4_german2_ci"))
_charsets.add(Charset(245, "utf8mb4", "utf8mb4_croatian_ci"))
_charsets.add(Charset(246, "utf8mb4", "utf8mb4_unicode_520_ci"))
_charsets.add(Charset(247, "utf8mb4", "utf8mb4_vietnamese_ci"))
_charsets.add(Charset(248, "gb18030", "gb18030_chinese_ci", True))
_charsets.add(Charset(249, "gb18030", "gb18030_bin"))
_charsets.add(Charset(250, "gb18030", "gb18030_unicode_520_ci"))
_charsets.add(Charset(255, "utf8mb4", "utf8mb4_0900_ai_ci"))
