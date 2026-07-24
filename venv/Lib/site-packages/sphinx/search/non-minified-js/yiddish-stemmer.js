// Generated from yiddish.sbl by Snowball 3.0.1 - https://snowballstem.org/

/**@constructor*/
var YiddishStemmer = function() {
    var base = new BaseStemmer();

    /** @const */ var a_0 = [
        ["\u05D5\u05D5", -1, 1],
        ["\u05D5\u05D9", -1, 2],
        ["\u05D9\u05D9", -1, 3],
        ["\u05DA", -1, 4],
        ["\u05DD", -1, 5],
        ["\u05DF", -1, 6],
        ["\u05E3", -1, 7],
        ["\u05E5", -1, 8]
    ];

    /** @const */ var a_1 = [
        ["\u05D0\u05D3\u05D5\u05E8\u05DB", -1, 1],
        ["\u05D0\u05D4\u05D9\u05E0", -1, 1],
        ["\u05D0\u05D4\u05E2\u05E8", -1, 1],
        ["\u05D0\u05D4\u05F2\u05DE", -1, 1],
        ["\u05D0\u05D5\u05DE", -1, 1],
        ["\u05D0\u05D5\u05E0\u05D8\u05E2\u05E8", -1, 1],
        ["\u05D0\u05D9\u05D1\u05E2\u05E8", -1, 1],
        ["\u05D0\u05E0", -1, 1],
        ["\u05D0\u05E0\u05D8", 7, 1],
        ["\u05D0\u05E0\u05D8\u05E7\u05E2\u05D2\u05E0", 8, 1],
        ["\u05D0\u05E0\u05D9\u05D3\u05E2\u05E8", 7, 1],
        ["\u05D0\u05E4", -1, 1],
        ["\u05D0\u05E4\u05D9\u05E8", 11, 1],
        ["\u05D0\u05E7\u05E2\u05D2\u05E0", -1, 1],
        ["\u05D0\u05E8\u05D0\u05E4", -1, 1],
        ["\u05D0\u05E8\u05D5\u05DE", -1, 1],
        ["\u05D0\u05E8\u05D5\u05E0\u05D8\u05E2\u05E8", -1, 1],
        ["\u05D0\u05E8\u05D9\u05D1\u05E2\u05E8", -1, 1],
        ["\u05D0\u05E8\u05F1\u05E1", -1, 1],
        ["\u05D0\u05E8\u05F1\u05E4", -1, 1],
        ["\u05D0\u05E8\u05F2\u05E0", -1, 1],
        ["\u05D0\u05F0\u05E2\u05E7", -1, 1],
        ["\u05D0\u05F1\u05E1", -1, 1],
        ["\u05D0\u05F1\u05E4", -1, 1],
        ["\u05D0\u05F2\u05E0", -1, 1],
        ["\u05D1\u05D0", -1, 1],
        ["\u05D1\u05F2", -1, 1],
        ["\u05D3\u05D5\u05E8\u05DB", -1, 1],
        ["\u05D3\u05E2\u05E8", -1, 1],
        ["\u05DE\u05D9\u05D8", -1, 1],
        ["\u05E0\u05D0\u05DB", -1, 1],
        ["\u05E4\u05D0\u05E8", -1, 1],
        ["\u05E4\u05D0\u05E8\u05D1\u05F2", 31, 1],
        ["\u05E4\u05D0\u05E8\u05F1\u05E1", 31, 1],
        ["\u05E4\u05D5\u05E0\u05D0\u05E0\u05D3\u05E2\u05E8", -1, 1],
        ["\u05E6\u05D5", -1, 1],
        ["\u05E6\u05D5\u05D6\u05D0\u05DE\u05E2\u05E0", 35, 1],
        ["\u05E6\u05D5\u05E0\u05F1\u05E4", 35, 1],
        ["\u05E6\u05D5\u05E8\u05D9\u05E7", 35, 1],
        ["\u05E6\u05E2", -1, 1]
    ];

    /** @const */ var a_2 = [
        ["\u05D3\u05D6\u05E9", -1, -1],
        ["\u05E9\u05D8\u05E8", -1, -1],
        ["\u05E9\u05D8\u05E9", -1, -1],
        ["\u05E9\u05E4\u05E8", -1, -1]
    ];

    /** @const */ var a_3 = [
        ["\u05E7\u05DC\u05D9\u05D1", -1, 9],
        ["\u05E8\u05D9\u05D1", -1, 10],
        ["\u05D8\u05E8\u05D9\u05D1", 1, 7],
        ["\u05E9\u05E8\u05D9\u05D1", 1, 15],
        ["\u05D4\u05F1\u05D1", -1, 23],
        ["\u05E9\u05F0\u05D9\u05D2", -1, 12],
        ["\u05D2\u05D0\u05E0\u05D2", -1, 1],
        ["\u05D6\u05D5\u05E0\u05D2", -1, 18],
        ["\u05E9\u05DC\u05D5\u05E0\u05D2", -1, 21],
        ["\u05E6\u05F0\u05D5\u05E0\u05D2", -1, 20],
        ["\u05D1\u05F1\u05D2", -1, 22],
        ["\u05D1\u05D5\u05E0\u05D3", -1, 16],
        ["\u05F0\u05D9\u05D6", -1, 6],
        ["\u05D1\u05D9\u05D8", -1, 4],
        ["\u05DC\u05D9\u05D8", -1, 8],
        ["\u05DE\u05D9\u05D8", -1, 3],
        ["\u05E9\u05E0\u05D9\u05D8", -1, 14],
        ["\u05E0\u05D5\u05DE", -1, 2],
        ["\u05E9\u05D8\u05D0\u05E0", -1, 25],
        ["\u05D1\u05D9\u05E1", -1, 5],
        ["\u05E9\u05DE\u05D9\u05E1", -1, 13],
        ["\u05E8\u05D9\u05E1", -1, 11],
        ["\u05D8\u05E8\u05D5\u05E0\u05E7", -1, 19],
        ["\u05E4\u05D0\u05E8\u05DC\u05F1\u05E8", -1, 24],
        ["\u05E9\u05F0\u05F1\u05E8", -1, 26],
        ["\u05F0\u05D5\u05D8\u05E9", -1, 17]
    ];

    /** @const */ var a_4 = [
        ["\u05D5\u05E0\u05D2", -1, 1],
        ["\u05E1\u05D8\u05D5", -1, 1],
        ["\u05D8", -1, 1],
        ["\u05D1\u05E8\u05D0\u05DB\u05D8", 2, 31],
        ["\u05E1\u05D8", 2, 1],
        ["\u05D9\u05E1\u05D8", 4, 33],
        ["\u05E2\u05D8", 2, 1],
        ["\u05E9\u05D0\u05E4\u05D8", 2, 1],
        ["\u05D4\u05F2\u05D8", 2, 1],
        ["\u05E7\u05F2\u05D8", 2, 1],
        ["\u05D9\u05E7\u05F2\u05D8", 9, 1],
        ["\u05DC\u05E2\u05DB", -1, 1],
        ["\u05E2\u05DC\u05E2\u05DB", 11, 1],
        ["\u05D9\u05D6\u05DE", -1, 1],
        ["\u05D9\u05DE", -1, 1],
        ["\u05E2\u05DE", -1, 1],
        ["\u05E2\u05E0\u05E2\u05DE", 15, 3],
        ["\u05D8\u05E2\u05E0\u05E2\u05DE", 16, 4],
        ["\u05E0", -1, 1],
        ["\u05E7\u05DC\u05D9\u05D1\u05E0", 18, 14],
        ["\u05E8\u05D9\u05D1\u05E0", 18, 15],
        ["\u05D8\u05E8\u05D9\u05D1\u05E0", 20, 12],
        ["\u05E9\u05E8\u05D9\u05D1\u05E0", 20, 7],
        ["\u05D4\u05F1\u05D1\u05E0", 18, 27],
        ["\u05E9\u05F0\u05D9\u05D2\u05E0", 18, 17],
        ["\u05D6\u05D5\u05E0\u05D2\u05E0", 18, 22],
        ["\u05E9\u05DC\u05D5\u05E0\u05D2\u05E0", 18, 25],
        ["\u05E6\u05F0\u05D5\u05E0\u05D2\u05E0", 18, 24],
        ["\u05D1\u05F1\u05D2\u05E0", 18, 26],
        ["\u05D1\u05D5\u05E0\u05D3\u05E0", 18, 20],
        ["\u05F0\u05D9\u05D6\u05E0", 18, 11],
        ["\u05D8\u05E0", 18, 4],
        ["GE\u05D1\u05D9\u05D8\u05E0", 31, 9],
        ["GE\u05DC\u05D9\u05D8\u05E0", 31, 13],
        ["GE\u05DE\u05D9\u05D8\u05E0", 31, 8],
        ["\u05E9\u05E0\u05D9\u05D8\u05E0", 31, 19],
        ["\u05E1\u05D8\u05E0", 31, 1],
        ["\u05D9\u05E1\u05D8\u05E0", 36, 1],
        ["\u05E2\u05D8\u05E0", 31, 1],
        ["GE\u05D1\u05D9\u05E1\u05E0", 18, 10],
        ["\u05E9\u05DE\u05D9\u05E1\u05E0", 18, 18],
        ["GE\u05E8\u05D9\u05E1\u05E0", 18, 16],
        ["\u05E2\u05E0", 18, 1],
        ["\u05D2\u05D0\u05E0\u05D2\u05E2\u05E0", 42, 5],
        ["\u05E2\u05DC\u05E2\u05E0", 42, 1],
        ["\u05E0\u05D5\u05DE\u05E2\u05E0", 42, 6],
        ["\u05D9\u05D6\u05DE\u05E2\u05E0", 42, 1],
        ["\u05E9\u05D8\u05D0\u05E0\u05E2\u05E0", 42, 29],
        ["\u05D8\u05E8\u05D5\u05E0\u05E7\u05E0", 18, 23],
        ["\u05E4\u05D0\u05E8\u05DC\u05F1\u05E8\u05E0", 18, 28],
        ["\u05E9\u05F0\u05F1\u05E8\u05E0", 18, 30],
        ["\u05F0\u05D5\u05D8\u05E9\u05E0", 18, 21],
        ["\u05D2\u05F2\u05E0", 18, 5],
        ["\u05E1", -1, 1],
        ["\u05D8\u05E1", 53, 4],
        ["\u05E2\u05D8\u05E1", 54, 1],
        ["\u05E0\u05E1", 53, 1],
        ["\u05D8\u05E0\u05E1", 56, 4],
        ["\u05E2\u05E0\u05E1", 56, 3],
        ["\u05E2\u05E1", 53, 1],
        ["\u05D9\u05E2\u05E1", 59, 2],
        ["\u05E2\u05DC\u05E2\u05E1", 59, 1],
        ["\u05E2\u05E8\u05E1", 53, 1],
        ["\u05E2\u05E0\u05E2\u05E8\u05E1", 62, 1],
        ["\u05E2", -1, 1],
        ["\u05D8\u05E2", 64, 4],
        ["\u05E1\u05D8\u05E2", 65, 1],
        ["\u05E2\u05D8\u05E2", 65, 1],
        ["\u05D9\u05E2", 64, -1],
        ["\u05E2\u05DC\u05E2", 64, 1],
        ["\u05E2\u05E0\u05E2", 64, 3],
        ["\u05D8\u05E2\u05E0\u05E2", 70, 4],
        ["\u05E2\u05E8", -1, 1],
        ["\u05D8\u05E2\u05E8", 72, 4],
        ["\u05E1\u05D8\u05E2\u05E8", 73, 1],
        ["\u05E2\u05D8\u05E2\u05E8", 73, 1],
        ["\u05E2\u05E0\u05E2\u05E8", 72, 3],
        ["\u05D8\u05E2\u05E0\u05E2\u05E8", 76, 4],
        ["\u05D5\u05EA", -1, 32]
    ];

    /** @const */ var a_5 = [
        ["\u05D5\u05E0\u05D2", -1, 1],
        ["\u05E9\u05D0\u05E4\u05D8", -1, 1],
        ["\u05D4\u05F2\u05D8", -1, 1],
        ["\u05E7\u05F2\u05D8", -1, 1],
        ["\u05D9\u05E7\u05F2\u05D8", 3, 1],
        ["\u05DC", -1, 2]
    ];

    /** @const */ var a_6 = [
        ["\u05D9\u05D2", -1, 1],
        ["\u05D9\u05E7", -1, 1],
        ["\u05D3\u05D9\u05E7", 1, 1],
        ["\u05E0\u05D3\u05D9\u05E7", 2, 1],
        ["\u05E2\u05E0\u05D3\u05D9\u05E7", 3, 1],
        ["\u05D1\u05DC\u05D9\u05E7", 1, -1],
        ["\u05D2\u05DC\u05D9\u05E7", 1, -1],
        ["\u05E0\u05D9\u05E7", 1, 1],
        ["\u05D9\u05E9", -1, 1]
    ];

    /** @const */ var /** Array<int> */ g_niked = [255, 155, 6];

    /** @const */ var /** Array<int> */ g_vowel = [33, 2, 4, 0, 6];

    /** @const */ var /** Array<int> */ g_consonant = [239, 254, 253, 131];

    var /** number */ I_x = 0;
    var /** number */ I_p1 = 0;


    /** @return {boolean} */
    function r_prelude() {
        var /** number */ among_var;
        /** @const */ var /** number */ v_1 = base.cursor;
        lab0: {
            while(true)
            {
                /** @const */ var /** number */ v_2 = base.cursor;
                lab1: {
                    golab2: while(true)
                    {
                        /** @const */ var /** number */ v_3 = base.cursor;
                        lab3: {
                            base.bra = base.cursor;
                            among_var = base.find_among(a_0);
                            if (among_var == 0)
                            {
                                break lab3;
                            }
                            base.ket = base.cursor;
                            switch (among_var) {
                                case 1:
                                    {
                                        /** @const */ var /** number */ v_4 = base.cursor;
                                        lab4: {
                                            if (!(base.eq_s("\u05BC")))
                                            {
                                                break lab4;
                                            }
                                            break lab3;
                                        }
                                        base.cursor = v_4;
                                    }
                                    if (!base.slice_from("\u05F0"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 2:
                                    {
                                        /** @const */ var /** number */ v_5 = base.cursor;
                                        lab5: {
                                            if (!(base.eq_s("\u05B4")))
                                            {
                                                break lab5;
                                            }
                                            break lab3;
                                        }
                                        base.cursor = v_5;
                                    }
                                    if (!base.slice_from("\u05F1"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 3:
                                    {
                                        /** @const */ var /** number */ v_6 = base.cursor;
                                        lab6: {
                                            if (!(base.eq_s("\u05B4")))
                                            {
                                                break lab6;
                                            }
                                            break lab3;
                                        }
                                        base.cursor = v_6;
                                    }
                                    if (!base.slice_from("\u05F2"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 4:
                                    if (!base.slice_from("\u05DB"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 5:
                                    if (!base.slice_from("\u05DE"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 6:
                                    if (!base.slice_from("\u05E0"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 7:
                                    if (!base.slice_from("\u05E4"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 8:
                                    if (!base.slice_from("\u05E6"))
                                    {
                                        return false;
                                    }
                                    break;
                            }
                            base.cursor = v_3;
                            break golab2;
                        }
                        base.cursor = v_3;
                        if (base.cursor >= base.limit)
                        {
                            break lab1;
                        }
                        base.cursor++;
                    }
                    continue;
                }
                base.cursor = v_2;
                break;
            }
        }
        base.cursor = v_1;
        /** @const */ var /** number */ v_7 = base.cursor;
        lab7: {
            while(true)
            {
                /** @const */ var /** number */ v_8 = base.cursor;
                lab8: {
                    golab9: while(true)
                    {
                        /** @const */ var /** number */ v_9 = base.cursor;
                        lab10: {
                            base.bra = base.cursor;
                            if (!(base.in_grouping(g_niked, 1456, 1474)))
                            {
                                break lab10;
                            }
                            base.ket = base.cursor;
                            if (!base.slice_del())
                            {
                                return false;
                            }
                            base.cursor = v_9;
                            break golab9;
                        }
                        base.cursor = v_9;
                        if (base.cursor >= base.limit)
                        {
                            break lab8;
                        }
                        base.cursor++;
                    }
                    continue;
                }
                base.cursor = v_8;
                break;
            }
        }
        base.cursor = v_7;
        return true;
    };

    /** @return {boolean} */
    function r_mark_regions() {
        I_p1 = base.limit;
        /** @const */ var /** number */ v_1 = base.cursor;
        lab0: {
            base.bra = base.cursor;
            if (!(base.eq_s("\u05D2\u05E2")))
            {
                base.cursor = v_1;
                break lab0;
            }
            base.ket = base.cursor;
            {
                /** @const */ var /** number */ v_2 = base.cursor;
                lab1: {
                    lab2: {
                        /** @const */ var /** number */ v_3 = base.cursor;
                        lab3: {
                            if (!(base.eq_s("\u05DC\u05D8")))
                            {
                                break lab3;
                            }
                            break lab2;
                        }
                        base.cursor = v_3;
                        lab4: {
                            if (!(base.eq_s("\u05D1\u05E0")))
                            {
                                break lab4;
                            }
                            break lab2;
                        }
                        base.cursor = v_3;
                        if (base.cursor < base.limit)
                        {
                            break lab1;
                        }
                    }
                    base.cursor = v_1;
                    break lab0;
                }
                base.cursor = v_2;
            }
            if (!base.slice_from("GE"))
            {
                return false;
            }
        }
        /** @const */ var /** number */ v_4 = base.cursor;
        lab5: {
            if (base.find_among(a_1) == 0)
            {
                base.cursor = v_4;
                break lab5;
            }
            lab6: {
                /** @const */ var /** number */ v_5 = base.cursor;
                lab7: {
                    /** @const */ var /** number */ v_6 = base.cursor;
                    lab8: {
                        /** @const */ var /** number */ v_7 = base.cursor;
                        lab9: {
                            if (!(base.eq_s("\u05E6\u05D5\u05D2\u05E0")))
                            {
                                break lab9;
                            }
                            break lab8;
                        }
                        base.cursor = v_7;
                        lab10: {
                            if (!(base.eq_s("\u05E6\u05D5\u05E7\u05D8")))
                            {
                                break lab10;
                            }
                            break lab8;
                        }
                        base.cursor = v_7;
                        if (!(base.eq_s("\u05E6\u05D5\u05E7\u05E0")))
                        {
                            break lab7;
                        }
                    }
                    if (base.cursor < base.limit)
                    {
                        break lab7;
                    }
                    base.cursor = v_6;
                    break lab6;
                }
                base.cursor = v_5;
                lab11: {
                    /** @const */ var /** number */ v_8 = base.cursor;
                    if (!(base.eq_s("\u05D2\u05E2\u05D1\u05E0")))
                    {
                        break lab11;
                    }
                    base.cursor = v_8;
                    break lab6;
                }
                base.cursor = v_5;
                lab12: {
                    base.bra = base.cursor;
                    if (!(base.eq_s("\u05D2\u05E2")))
                    {
                        break lab12;
                    }
                    base.ket = base.cursor;
                    if (!base.slice_from("GE"))
                    {
                        return false;
                    }
                    break lab6;
                }
                base.cursor = v_5;
                base.bra = base.cursor;
                if (!(base.eq_s("\u05E6\u05D5")))
                {
                    base.cursor = v_4;
                    break lab5;
                }
                base.ket = base.cursor;
                if (!base.slice_from("TSU"))
                {
                    return false;
                }
            }
        }
        /** @const */ var /** number */ v_9 = base.cursor;
        {
            /** @const */ var /** number */ c1 = base.cursor + 3;
            if (c1 > base.limit)
            {
                return false;
            }
            base.cursor = c1;
        }
        I_x = base.cursor;
        base.cursor = v_9;
        /** @const */ var /** number */ v_10 = base.cursor;
        lab13: {
            if (base.find_among(a_2) == 0)
            {
                base.cursor = v_10;
                break lab13;
            }
        }
        {
            /** @const */ var /** number */ v_11 = base.cursor;
            lab14: {
                if (!(base.in_grouping(g_consonant, 1489, 1520)))
                {
                    break lab14;
                }
                if (!(base.in_grouping(g_consonant, 1489, 1520)))
                {
                    break lab14;
                }
                if (!(base.in_grouping(g_consonant, 1489, 1520)))
                {
                    break lab14;
                }
                I_p1 = base.cursor;
                return false;
            }
            base.cursor = v_11;
        }
        if (!base.go_out_grouping(g_vowel, 1488, 1522))
        {
            return false;
        }
        base.cursor++;
        if (!base.go_in_grouping(g_vowel, 1488, 1522))
        {
            return false;
        }
        I_p1 = base.cursor;
        lab15: {
            if (I_p1 >= I_x)
            {
                break lab15;
            }
            I_p1 = I_x;
        }
        return true;
    };

    /** @return {boolean} */
    function r_R1() {
        return I_p1 <= base.cursor;
    };

    /** @return {boolean} */
    function r_R1plus3() {
        return I_p1 <= (base.cursor + 3);
    };

    /** @return {boolean} */
    function r_standard_suffix() {
        var /** number */ among_var;
        /** @const */ var /** number */ v_1 = base.limit - base.cursor;
        lab0: {
            base.ket = base.cursor;
            among_var = base.find_among_b(a_4);
            if (among_var == 0)
            {
                break lab0;
            }
            base.bra = base.cursor;
            switch (among_var) {
                case 1:
                    if (!r_R1())
                    {
                        break lab0;
                    }
                    if (!base.slice_del())
                    {
                        return false;
                    }
                    break;
                case 2:
                    if (!r_R1())
                    {
                        break lab0;
                    }
                    if (!base.slice_from("\u05D9\u05E2"))
                    {
                        return false;
                    }
                    break;
                case 3:
                    if (!r_R1())
                    {
                        break lab0;
                    }
                    if (!base.slice_del())
                    {
                        return false;
                    }
                    base.ket = base.cursor;
                    among_var = base.find_among_b(a_3);
                    if (among_var == 0)
                    {
                        break lab0;
                    }
                    base.bra = base.cursor;
                    switch (among_var) {
                        case 1:
                            if (!base.slice_from("\u05D2\u05F2"))
                            {
                                return false;
                            }
                            break;
                        case 2:
                            if (!base.slice_from("\u05E0\u05E2\u05DE"))
                            {
                                return false;
                            }
                            break;
                        case 3:
                            if (!base.slice_from("\u05DE\u05F2\u05D3"))
                            {
                                return false;
                            }
                            break;
                        case 4:
                            if (!base.slice_from("\u05D1\u05F2\u05D8"))
                            {
                                return false;
                            }
                            break;
                        case 5:
                            if (!base.slice_from("\u05D1\u05F2\u05E1"))
                            {
                                return false;
                            }
                            break;
                        case 6:
                            if (!base.slice_from("\u05F0\u05F2\u05D6"))
                            {
                                return false;
                            }
                            break;
                        case 7:
                            if (!base.slice_from("\u05D8\u05E8\u05F2\u05D1"))
                            {
                                return false;
                            }
                            break;
                        case 8:
                            if (!base.slice_from("\u05DC\u05F2\u05D8"))
                            {
                                return false;
                            }
                            break;
                        case 9:
                            if (!base.slice_from("\u05E7\u05DC\u05F2\u05D1"))
                            {
                                return false;
                            }
                            break;
                        case 10:
                            if (!base.slice_from("\u05E8\u05F2\u05D1"))
                            {
                                return false;
                            }
                            break;
                        case 11:
                            if (!base.slice_from("\u05E8\u05F2\u05E1"))
                            {
                                return false;
                            }
                            break;
                        case 12:
                            if (!base.slice_from("\u05E9\u05F0\u05F2\u05D2"))
                            {
                                return false;
                            }
                            break;
                        case 13:
                            if (!base.slice_from("\u05E9\u05DE\u05F2\u05E1"))
                            {
                                return false;
                            }
                            break;
                        case 14:
                            if (!base.slice_from("\u05E9\u05E0\u05F2\u05D3"))
                            {
                                return false;
                            }
                            break;
                        case 15:
                            if (!base.slice_from("\u05E9\u05E8\u05F2\u05D1"))
                            {
                                return false;
                            }
                            break;
                        case 16:
                            if (!base.slice_from("\u05D1\u05D9\u05E0\u05D3"))
                            {
                                return false;
                            }
                            break;
                        case 17:
                            if (!base.slice_from("\u05F0\u05D9\u05D8\u05E9"))
                            {
                                return false;
                            }
                            break;
                        case 18:
                            if (!base.slice_from("\u05D6\u05D9\u05E0\u05D2"))
                            {
                                return false;
                            }
                            break;
                        case 19:
                            if (!base.slice_from("\u05D8\u05E8\u05D9\u05E0\u05E7"))
                            {
                                return false;
                            }
                            break;
                        case 20:
                            if (!base.slice_from("\u05E6\u05F0\u05D9\u05E0\u05D2"))
                            {
                                return false;
                            }
                            break;
                        case 21:
                            if (!base.slice_from("\u05E9\u05DC\u05D9\u05E0\u05D2"))
                            {
                                return false;
                            }
                            break;
                        case 22:
                            if (!base.slice_from("\u05D1\u05F2\u05D2"))
                            {
                                return false;
                            }
                            break;
                        case 23:
                            if (!base.slice_from("\u05D4\u05F2\u05D1"))
                            {
                                return false;
                            }
                            break;
                        case 24:
                            if (!base.slice_from("\u05E4\u05D0\u05E8\u05DC\u05D9\u05E8"))
                            {
                                return false;
                            }
                            break;
                        case 25:
                            if (!base.slice_from("\u05E9\u05D8\u05F2"))
                            {
                                return false;
                            }
                            break;
                        case 26:
                            if (!base.slice_from("\u05E9\u05F0\u05E2\u05E8"))
                            {
                                return false;
                            }
                            break;
                    }
                    break;
                case 4:
                    lab1: {
                        /** @const */ var /** number */ v_2 = base.limit - base.cursor;
                        lab2: {
                            if (!r_R1())
                            {
                                break lab2;
                            }
                            if (!base.slice_del())
                            {
                                return false;
                            }
                            break lab1;
                        }
                        base.cursor = base.limit - v_2;
                        if (!base.slice_from("\u05D8"))
                        {
                            return false;
                        }
                    }
                    base.ket = base.cursor;
                    if (!(base.eq_s_b("\u05D1\u05E8\u05D0\u05DB")))
                    {
                        break lab0;
                    }
                    /** @const */ var /** number */ v_3 = base.limit - base.cursor;
                    lab3: {
                        if (!(base.eq_s_b("\u05D2\u05E2")))
                        {
                            base.cursor = base.limit - v_3;
                            break lab3;
                        }
                    }
                    base.bra = base.cursor;
                    if (!base.slice_from("\u05D1\u05E8\u05E2\u05E0\u05D2"))
                    {
                        return false;
                    }
                    break;
                case 5:
                    if (!base.slice_from("\u05D2\u05F2"))
                    {
                        return false;
                    }
                    break;
                case 6:
                    if (!base.slice_from("\u05E0\u05E2\u05DE"))
                    {
                        return false;
                    }
                    break;
                case 7:
                    if (!base.slice_from("\u05E9\u05E8\u05F2\u05D1"))
                    {
                        return false;
                    }
                    break;
                case 8:
                    if (!base.slice_from("\u05DE\u05F2\u05D3"))
                    {
                        return false;
                    }
                    break;
                case 9:
                    if (!base.slice_from("\u05D1\u05F2\u05D8"))
                    {
                        return false;
                    }
                    break;
                case 10:
                    if (!base.slice_from("\u05D1\u05F2\u05E1"))
                    {
                        return false;
                    }
                    break;
                case 11:
                    if (!base.slice_from("\u05F0\u05F2\u05D6"))
                    {
                        return false;
                    }
                    break;
                case 12:
                    if (!base.slice_from("\u05D8\u05E8\u05F2\u05D1"))
                    {
                        return false;
                    }
                    break;
                case 13:
                    if (!base.slice_from("\u05DC\u05F2\u05D8"))
                    {
                        return false;
                    }
                    break;
                case 14:
                    if (!base.slice_from("\u05E7\u05DC\u05F2\u05D1"))
                    {
                        return false;
                    }
                    break;
                case 15:
                    if (!base.slice_from("\u05E8\u05F2\u05D1"))
                    {
                        return false;
                    }
                    break;
                case 16:
                    if (!base.slice_from("\u05E8\u05F2\u05E1"))
                    {
                        return false;
                    }
                    break;
                case 17:
                    if (!base.slice_from("\u05E9\u05F0\u05F2\u05D2"))
                    {
                        return false;
                    }
                    break;
                case 18:
                    if (!base.slice_from("\u05E9\u05DE\u05F2\u05E1"))
                    {
                        return false;
                    }
                    break;
                case 19:
                    if (!base.slice_from("\u05E9\u05E0\u05F2\u05D3"))
                    {
                        return false;
                    }
                    break;
                case 20:
                    if (!base.slice_from("\u05D1\u05D9\u05E0\u05D3"))
                    {
                        return false;
                    }
                    break;
                case 21:
                    if (!base.slice_from("\u05F0\u05D9\u05D8\u05E9"))
                    {
                        return false;
                    }
                    break;
                case 22:
                    if (!base.slice_from("\u05D6\u05D9\u05E0\u05D2"))
                    {
                        return false;
                    }
                    break;
                case 23:
                    if (!base.slice_from("\u05D8\u05E8\u05D9\u05E0\u05E7"))
                    {
                        return false;
                    }
                    break;
                case 24:
                    if (!base.slice_from("\u05E6\u05F0\u05D9\u05E0\u05D2"))
                    {
                        return false;
                    }
                    break;
                case 25:
                    if (!base.slice_from("\u05E9\u05DC\u05D9\u05E0\u05D2"))
                    {
                        return false;
                    }
                    break;
                case 26:
                    if (!base.slice_from("\u05D1\u05F2\u05D2"))
                    {
                        return false;
                    }
                    break;
                case 27:
                    if (!base.slice_from("\u05D4\u05F2\u05D1"))
                    {
                        return false;
                    }
                    break;
                case 28:
                    if (!base.slice_from("\u05E4\u05D0\u05E8\u05DC\u05D9\u05E8"))
                    {
                        return false;
                    }
                    break;
                case 29:
                    if (!base.slice_from("\u05E9\u05D8\u05F2"))
                    {
                        return false;
                    }
                    break;
                case 30:
                    if (!base.slice_from("\u05E9\u05F0\u05E2\u05E8"))
                    {
                        return false;
                    }
                    break;
                case 31:
                    if (!base.slice_from("\u05D1\u05E8\u05E2\u05E0\u05D2"))
                    {
                        return false;
                    }
                    break;
                case 32:
                    if (!r_R1())
                    {
                        break lab0;
                    }
                    if (!base.slice_from("\u05D4"))
                    {
                        return false;
                    }
                    break;
                case 33:
                    lab4: {
                        /** @const */ var /** number */ v_4 = base.limit - base.cursor;
                        lab5: {
                            lab6: {
                                /** @const */ var /** number */ v_5 = base.limit - base.cursor;
                                lab7: {
                                    if (!(base.eq_s_b("\u05D2")))
                                    {
                                        break lab7;
                                    }
                                    break lab6;
                                }
                                base.cursor = base.limit - v_5;
                                if (!(base.eq_s_b("\u05E9")))
                                {
                                    break lab5;
                                }
                            }
                            /** @const */ var /** number */ v_6 = base.limit - base.cursor;
                            lab8: {
                                if (!r_R1plus3())
                                {
                                    base.cursor = base.limit - v_6;
                                    break lab8;
                                }
                                if (!base.slice_from("\u05D9\u05E1"))
                                {
                                    return false;
                                }
                            }
                            break lab4;
                        }
                        base.cursor = base.limit - v_4;
                        if (!r_R1())
                        {
                            break lab0;
                        }
                        if (!base.slice_del())
                        {
                            return false;
                        }
                    }
                    break;
            }
        }
        base.cursor = base.limit - v_1;
        /** @const */ var /** number */ v_7 = base.limit - base.cursor;
        lab9: {
            base.ket = base.cursor;
            among_var = base.find_among_b(a_5);
            if (among_var == 0)
            {
                break lab9;
            }
            base.bra = base.cursor;
            switch (among_var) {
                case 1:
                    if (!r_R1())
                    {
                        break lab9;
                    }
                    if (!base.slice_del())
                    {
                        return false;
                    }
                    break;
                case 2:
                    if (!r_R1())
                    {
                        break lab9;
                    }
                    if (!(base.in_grouping_b(g_consonant, 1489, 1520)))
                    {
                        break lab9;
                    }
                    if (!base.slice_del())
                    {
                        return false;
                    }
                    break;
            }
        }
        base.cursor = base.limit - v_7;
        /** @const */ var /** number */ v_8 = base.limit - base.cursor;
        lab10: {
            base.ket = base.cursor;
            among_var = base.find_among_b(a_6);
            if (among_var == 0)
            {
                break lab10;
            }
            base.bra = base.cursor;
            switch (among_var) {
                case 1:
                    if (!r_R1())
                    {
                        break lab10;
                    }
                    if (!base.slice_del())
                    {
                        return false;
                    }
                    break;
            }
        }
        base.cursor = base.limit - v_8;
        /** @const */ var /** number */ v_9 = base.limit - base.cursor;
        lab11: {
            while(true)
            {
                /** @const */ var /** number */ v_10 = base.limit - base.cursor;
                lab12: {
                    golab13: while(true)
                    {
                        /** @const */ var /** number */ v_11 = base.limit - base.cursor;
                        lab14: {
                            base.ket = base.cursor;
                            lab15: {
                                /** @const */ var /** number */ v_12 = base.limit - base.cursor;
                                lab16: {
                                    if (!(base.eq_s_b("GE")))
                                    {
                                        break lab16;
                                    }
                                    break lab15;
                                }
                                base.cursor = base.limit - v_12;
                                if (!(base.eq_s_b("TSU")))
                                {
                                    break lab14;
                                }
                            }
                            base.bra = base.cursor;
                            if (!base.slice_del())
                            {
                                return false;
                            }
                            base.cursor = base.limit - v_11;
                            break golab13;
                        }
                        base.cursor = base.limit - v_11;
                        if (base.cursor <= base.limit_backward)
                        {
                            break lab12;
                        }
                        base.cursor--;
                    }
                    continue;
                }
                base.cursor = base.limit - v_10;
                break;
            }
        }
        base.cursor = base.limit - v_9;
        return true;
    };

    this.stem = /** @return {boolean} */ function() {
        r_prelude();
        /** @const */ var /** number */ v_1 = base.cursor;
        r_mark_regions();
        base.cursor = v_1;
        base.limit_backward = base.cursor; base.cursor = base.limit;
        r_standard_suffix();
        base.cursor = base.limit_backward;
        return true;
    };

    /**@return{string}*/
    this['stemWord'] = function(/**string*/word) {
        base.setCurrent(word);
        this.stem();
        return base.getCurrent();
    };
};
