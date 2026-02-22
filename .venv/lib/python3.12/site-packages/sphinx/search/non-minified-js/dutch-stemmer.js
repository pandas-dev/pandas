// Generated from dutch.sbl by Snowball 3.0.1 - https://snowballstem.org/

/**@constructor*/
var DutchStemmer = function() {
    var base = new BaseStemmer();

    /** @const */ var a_0 = [
        ["a", -1, 1],
        ["e", -1, 2],
        ["o", -1, 1],
        ["u", -1, 1],
        ["\u00E0", -1, 1],
        ["\u00E1", -1, 1],
        ["\u00E2", -1, 1],
        ["\u00E4", -1, 1],
        ["\u00E8", -1, 2],
        ["\u00E9", -1, 2],
        ["\u00EA", -1, 2],
        ["e\u00EB", -1, 3],
        ["i\u00EB", -1, 4],
        ["\u00F2", -1, 1],
        ["\u00F3", -1, 1],
        ["\u00F4", -1, 1],
        ["\u00F6", -1, 1],
        ["\u00F9", -1, 1],
        ["\u00FA", -1, 1],
        ["\u00FB", -1, 1],
        ["\u00FC", -1, 1]
    ];

    /** @const */ var a_1 = [
        ["nde", -1, 8],
        ["en", -1, 7],
        ["s", -1, 2],
        ["'s", 2, 1],
        ["es", 2, 4],
        ["ies", 4, 3],
        ["aus", 2, 6],
        ["\u00E9s", 2, 5]
    ];

    /** @const */ var a_2 = [
        ["de", -1, 5],
        ["ge", -1, 2],
        ["ische", -1, 4],
        ["je", -1, 1],
        ["lijke", -1, 3],
        ["le", -1, 9],
        ["ene", -1, 10],
        ["re", -1, 8],
        ["se", -1, 7],
        ["te", -1, 6],
        ["ieve", -1, 11]
    ];

    /** @const */ var a_3 = [
        ["heid", -1, 3],
        ["fie", -1, 7],
        ["gie", -1, 8],
        ["atie", -1, 1],
        ["isme", -1, 5],
        ["ing", -1, 5],
        ["arij", -1, 6],
        ["erij", -1, 5],
        ["sel", -1, 3],
        ["rder", -1, 4],
        ["ster", -1, 3],
        ["iteit", -1, 2],
        ["dst", -1, 10],
        ["tst", -1, 9]
    ];

    /** @const */ var a_4 = [
        ["end", -1, 9],
        ["atief", -1, 2],
        ["erig", -1, 9],
        ["achtig", -1, 3],
        ["ioneel", -1, 1],
        ["baar", -1, 3],
        ["laar", -1, 5],
        ["naar", -1, 4],
        ["raar", -1, 6],
        ["eriger", -1, 9],
        ["achtiger", -1, 3],
        ["lijker", -1, 8],
        ["tant", -1, 7],
        ["erigst", -1, 9],
        ["achtigst", -1, 3],
        ["lijkst", -1, 8]
    ];

    /** @const */ var a_5 = [
        ["ig", -1, 1],
        ["iger", -1, 1],
        ["igst", -1, 1]
    ];

    /** @const */ var a_6 = [
        ["ft", -1, 2],
        ["kt", -1, 1],
        ["pt", -1, 3]
    ];

    /** @const */ var a_7 = [
        ["bb", -1, 1],
        ["cc", -1, 2],
        ["dd", -1, 3],
        ["ff", -1, 4],
        ["gg", -1, 5],
        ["hh", -1, 6],
        ["jj", -1, 7],
        ["kk", -1, 8],
        ["ll", -1, 9],
        ["mm", -1, 10],
        ["nn", -1, 11],
        ["pp", -1, 12],
        ["qq", -1, 13],
        ["rr", -1, 14],
        ["ss", -1, 15],
        ["tt", -1, 16],
        ["v", -1, 4],
        ["vv", 16, 17],
        ["ww", -1, 18],
        ["xx", -1, 19],
        ["z", -1, 15],
        ["zz", 20, 20]
    ];

    /** @const */ var a_8 = [
        ["d", -1, 1],
        ["t", -1, 2]
    ];

    /** @const */ var a_9 = [
        ["", -1, -1],
        ["eft", 0, 1],
        ["vaa", 0, 1],
        ["val", 0, 1],
        ["vali", 3, -1],
        ["vare", 0, 1]
    ];

    /** @const */ var a_10 = [
        ["\u00EB", -1, 1],
        ["\u00EF", -1, 2]
    ];

    /** @const */ var a_11 = [
        ["\u00EB", -1, 1],
        ["\u00EF", -1, 2]
    ];

    /** @const */ var /** Array<int> */ g_E = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 120];

    /** @const */ var /** Array<int> */ g_AIOU = [1, 65, 16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 128, 11, 120, 46, 15];

    /** @const */ var /** Array<int> */ g_AEIOU = [17, 65, 16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 128, 139, 127, 46, 15];

    /** @const */ var /** Array<int> */ g_v = [17, 65, 16, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 128, 139, 127, 46, 15];

    /** @const */ var /** Array<int> */ g_v_WX = [17, 65, 208, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 128, 139, 127, 46, 15];

    var /** boolean */ B_GE_removed = false;
    var /** boolean */ B_stemmed = false;
    var /** number */ I_p2 = 0;
    var /** number */ I_p1 = 0;
    var /** string */ S_ch = '';


    /** @return {boolean} */
    function r_R1() {
        return I_p1 <= base.cursor;
    };

    /** @return {boolean} */
    function r_R2() {
        return I_p2 <= base.cursor;
    };

    /** @return {boolean} */
    function r_V() {
        /** @const */ var /** number */ v_1 = base.limit - base.cursor;
        lab0: {
            /** @const */ var /** number */ v_2 = base.limit - base.cursor;
            lab1: {
                if (!(base.in_grouping_b(g_v, 97, 252)))
                {
                    break lab1;
                }
                break lab0;
            }
            base.cursor = base.limit - v_2;
            if (!(base.eq_s_b("ij")))
            {
                return false;
            }
        }
        base.cursor = base.limit - v_1;
        return true;
    };

    /** @return {boolean} */
    function r_VX() {
        /** @const */ var /** number */ v_1 = base.limit - base.cursor;
        if (base.cursor <= base.limit_backward)
        {
            return false;
        }
        base.cursor--;
        lab0: {
            /** @const */ var /** number */ v_2 = base.limit - base.cursor;
            lab1: {
                if (!(base.in_grouping_b(g_v, 97, 252)))
                {
                    break lab1;
                }
                break lab0;
            }
            base.cursor = base.limit - v_2;
            if (!(base.eq_s_b("ij")))
            {
                return false;
            }
        }
        base.cursor = base.limit - v_1;
        return true;
    };

    /** @return {boolean} */
    function r_C() {
        /** @const */ var /** number */ v_1 = base.limit - base.cursor;
        {
            /** @const */ var /** number */ v_2 = base.limit - base.cursor;
            lab0: {
                if (!(base.eq_s_b("ij")))
                {
                    break lab0;
                }
                return false;
            }
            base.cursor = base.limit - v_2;
        }
        if (!(base.out_grouping_b(g_v, 97, 252)))
        {
            return false;
        }
        base.cursor = base.limit - v_1;
        return true;
    };

    /** @return {boolean} */
    function r_lengthen_V() {
        var /** number */ among_var;
        /** @const */ var /** number */ v_1 = base.limit - base.cursor;
        lab0: {
            if (!(base.out_grouping_b(g_v_WX, 97, 252)))
            {
                break lab0;
            }
            base.ket = base.cursor;
            among_var = base.find_among_b(a_0);
            if (among_var == 0)
            {
                break lab0;
            }
            base.bra = base.cursor;
            switch (among_var) {
                case 1:
                    /** @const */ var /** number */ v_2 = base.limit - base.cursor;
                    lab1: {
                        /** @const */ var /** number */ v_3 = base.limit - base.cursor;
                        lab2: {
                            if (!(base.out_grouping_b(g_AEIOU, 97, 252)))
                            {
                                break lab2;
                            }
                            break lab1;
                        }
                        base.cursor = base.limit - v_3;
                        if (base.cursor > base.limit_backward)
                        {
                            break lab0;
                        }
                    }
                    base.cursor = base.limit - v_2;
                    S_ch = base.slice_to();
                    if (S_ch == '')
                    {
                        return false;
                    }
                    {
                        /** @const */ var /** number */ c1 = base.cursor;
                        base.insert(base.cursor, base.cursor, S_ch);
                        base.cursor = c1;
                    }
                    break;
                case 2:
                    /** @const */ var /** number */ v_4 = base.limit - base.cursor;
                    lab3: {
                        /** @const */ var /** number */ v_5 = base.limit - base.cursor;
                        lab4: {
                            if (!(base.out_grouping_b(g_AEIOU, 97, 252)))
                            {
                                break lab4;
                            }
                            break lab3;
                        }
                        base.cursor = base.limit - v_5;
                        if (base.cursor > base.limit_backward)
                        {
                            break lab0;
                        }
                    }
                    {
                        /** @const */ var /** number */ v_6 = base.limit - base.cursor;
                        lab5: {
                            lab6: {
                                /** @const */ var /** number */ v_7 = base.limit - base.cursor;
                                lab7: {
                                    if (!(base.in_grouping_b(g_AIOU, 97, 252)))
                                    {
                                        break lab7;
                                    }
                                    break lab6;
                                }
                                base.cursor = base.limit - v_7;
                                if (!(base.in_grouping_b(g_E, 101, 235)))
                                {
                                    break lab5;
                                }
                                if (base.cursor > base.limit_backward)
                                {
                                    break lab5;
                                }
                            }
                            break lab0;
                        }
                        base.cursor = base.limit - v_6;
                    }
                    {
                        /** @const */ var /** number */ v_8 = base.limit - base.cursor;
                        lab8: {
                            if (base.cursor <= base.limit_backward)
                            {
                                break lab8;
                            }
                            base.cursor--;
                            if (!(base.in_grouping_b(g_AIOU, 97, 252)))
                            {
                                break lab8;
                            }
                            if (!(base.out_grouping_b(g_AEIOU, 97, 252)))
                            {
                                break lab8;
                            }
                            break lab0;
                        }
                        base.cursor = base.limit - v_8;
                    }
                    base.cursor = base.limit - v_4;
                    S_ch = base.slice_to();
                    if (S_ch == '')
                    {
                        return false;
                    }
                    {
                        /** @const */ var /** number */ c2 = base.cursor;
                        base.insert(base.cursor, base.cursor, S_ch);
                        base.cursor = c2;
                    }
                    break;
                case 3:
                    if (!base.slice_from("e\u00EBe"))
                    {
                        return false;
                    }
                    break;
                case 4:
                    if (!base.slice_from("iee"))
                    {
                        return false;
                    }
                    break;
            }
        }
        base.cursor = base.limit - v_1;
        return true;
    };

    /** @return {boolean} */
    function r_Step_1() {
        var /** number */ among_var;
        base.ket = base.cursor;
        among_var = base.find_among_b(a_1);
        if (among_var == 0)
        {
            return false;
        }
        base.bra = base.cursor;
        switch (among_var) {
            case 1:
                if (!base.slice_del())
                {
                    return false;
                }
                break;
            case 2:
                if (!r_R1())
                {
                    return false;
                }
                {
                    /** @const */ var /** number */ v_1 = base.limit - base.cursor;
                    lab0: {
                        if (!(base.eq_s_b("t")))
                        {
                            break lab0;
                        }
                        if (!r_R1())
                        {
                            break lab0;
                        }
                        return false;
                    }
                    base.cursor = base.limit - v_1;
                }
                if (!r_C())
                {
                    return false;
                }
                if (!base.slice_del())
                {
                    return false;
                }
                break;
            case 3:
                if (!r_R1())
                {
                    return false;
                }
                if (!base.slice_from("ie"))
                {
                    return false;
                }
                break;
            case 4:
                lab1: {
                    /** @const */ var /** number */ v_2 = base.limit - base.cursor;
                    lab2: {
                        /** @const */ var /** number */ v_3 = base.limit - base.cursor;
                        if (!(base.eq_s_b("ar")))
                        {
                            break lab2;
                        }
                        if (!r_R1())
                        {
                            break lab2;
                        }
                        if (!r_C())
                        {
                            break lab2;
                        }
                        base.cursor = base.limit - v_3;
                        if (!base.slice_del())
                        {
                            return false;
                        }
                        r_lengthen_V();
                        break lab1;
                    }
                    base.cursor = base.limit - v_2;
                    lab3: {
                        /** @const */ var /** number */ v_4 = base.limit - base.cursor;
                        if (!(base.eq_s_b("er")))
                        {
                            break lab3;
                        }
                        if (!r_R1())
                        {
                            break lab3;
                        }
                        if (!r_C())
                        {
                            break lab3;
                        }
                        base.cursor = base.limit - v_4;
                        if (!base.slice_del())
                        {
                            return false;
                        }
                        break lab1;
                    }
                    base.cursor = base.limit - v_2;
                    if (!r_R1())
                    {
                        return false;
                    }
                    if (!r_C())
                    {
                        return false;
                    }
                    if (!base.slice_from("e"))
                    {
                        return false;
                    }
                }
                break;
            case 5:
                if (!r_R1())
                {
                    return false;
                }
                if (!base.slice_from("\u00E9"))
                {
                    return false;
                }
                break;
            case 6:
                if (!r_R1())
                {
                    return false;
                }
                if (!r_V())
                {
                    return false;
                }
                if (!base.slice_from("au"))
                {
                    return false;
                }
                break;
            case 7:
                lab4: {
                    /** @const */ var /** number */ v_5 = base.limit - base.cursor;
                    lab5: {
                        if (!(base.eq_s_b("hed")))
                        {
                            break lab5;
                        }
                        if (!r_R1())
                        {
                            break lab5;
                        }
                        base.bra = base.cursor;
                        if (!base.slice_from("heid"))
                        {
                            return false;
                        }
                        break lab4;
                    }
                    base.cursor = base.limit - v_5;
                    lab6: {
                        if (!(base.eq_s_b("nd")))
                        {
                            break lab6;
                        }
                        if (!base.slice_del())
                        {
                            return false;
                        }
                        break lab4;
                    }
                    base.cursor = base.limit - v_5;
                    lab7: {
                        if (!(base.eq_s_b("d")))
                        {
                            break lab7;
                        }
                        if (!r_R1())
                        {
                            break lab7;
                        }
                        if (!r_C())
                        {
                            break lab7;
                        }
                        base.bra = base.cursor;
                        if (!base.slice_del())
                        {
                            return false;
                        }
                        break lab4;
                    }
                    base.cursor = base.limit - v_5;
                    lab8: {
                        lab9: {
                            /** @const */ var /** number */ v_6 = base.limit - base.cursor;
                            lab10: {
                                if (!(base.eq_s_b("i")))
                                {
                                    break lab10;
                                }
                                break lab9;
                            }
                            base.cursor = base.limit - v_6;
                            if (!(base.eq_s_b("j")))
                            {
                                break lab8;
                            }
                        }
                        if (!r_V())
                        {
                            break lab8;
                        }
                        if (!base.slice_del())
                        {
                            return false;
                        }
                        break lab4;
                    }
                    base.cursor = base.limit - v_5;
                    if (!r_R1())
                    {
                        return false;
                    }
                    if (!r_C())
                    {
                        return false;
                    }
                    if (!base.slice_del())
                    {
                        return false;
                    }
                    r_lengthen_V();
                }
                break;
            case 8:
                if (!base.slice_from("nd"))
                {
                    return false;
                }
                break;
        }
        return true;
    };

    /** @return {boolean} */
    function r_Step_2() {
        var /** number */ among_var;
        base.ket = base.cursor;
        among_var = base.find_among_b(a_2);
        if (among_var == 0)
        {
            return false;
        }
        base.bra = base.cursor;
        switch (among_var) {
            case 1:
                lab0: {
                    /** @const */ var /** number */ v_1 = base.limit - base.cursor;
                    lab1: {
                        if (!(base.eq_s_b("'t")))
                        {
                            break lab1;
                        }
                        base.bra = base.cursor;
                        if (!base.slice_del())
                        {
                            return false;
                        }
                        break lab0;
                    }
                    base.cursor = base.limit - v_1;
                    lab2: {
                        if (!(base.eq_s_b("et")))
                        {
                            break lab2;
                        }
                        base.bra = base.cursor;
                        if (!r_R1())
                        {
                            break lab2;
                        }
                        if (!r_C())
                        {
                            break lab2;
                        }
                        if (!base.slice_del())
                        {
                            return false;
                        }
                        break lab0;
                    }
                    base.cursor = base.limit - v_1;
                    lab3: {
                        if (!(base.eq_s_b("rnt")))
                        {
                            break lab3;
                        }
                        base.bra = base.cursor;
                        if (!base.slice_from("rn"))
                        {
                            return false;
                        }
                        break lab0;
                    }
                    base.cursor = base.limit - v_1;
                    lab4: {
                        if (!(base.eq_s_b("t")))
                        {
                            break lab4;
                        }
                        base.bra = base.cursor;
                        if (!r_R1())
                        {
                            break lab4;
                        }
                        if (!r_VX())
                        {
                            break lab4;
                        }
                        if (!base.slice_del())
                        {
                            return false;
                        }
                        break lab0;
                    }
                    base.cursor = base.limit - v_1;
                    lab5: {
                        if (!(base.eq_s_b("ink")))
                        {
                            break lab5;
                        }
                        base.bra = base.cursor;
                        if (!base.slice_from("ing"))
                        {
                            return false;
                        }
                        break lab0;
                    }
                    base.cursor = base.limit - v_1;
                    lab6: {
                        if (!(base.eq_s_b("mp")))
                        {
                            break lab6;
                        }
                        base.bra = base.cursor;
                        if (!base.slice_from("m"))
                        {
                            return false;
                        }
                        break lab0;
                    }
                    base.cursor = base.limit - v_1;
                    lab7: {
                        if (!(base.eq_s_b("'")))
                        {
                            break lab7;
                        }
                        base.bra = base.cursor;
                        if (!r_R1())
                        {
                            break lab7;
                        }
                        if (!base.slice_del())
                        {
                            return false;
                        }
                        break lab0;
                    }
                    base.cursor = base.limit - v_1;
                    base.bra = base.cursor;
                    if (!r_R1())
                    {
                        return false;
                    }
                    if (!r_C())
                    {
                        return false;
                    }
                    if (!base.slice_del())
                    {
                        return false;
                    }
                }
                break;
            case 2:
                if (!r_R1())
                {
                    return false;
                }
                if (!base.slice_from("g"))
                {
                    return false;
                }
                break;
            case 3:
                if (!r_R1())
                {
                    return false;
                }
                if (!base.slice_from("lijk"))
                {
                    return false;
                }
                break;
            case 4:
                if (!r_R1())
                {
                    return false;
                }
                if (!base.slice_from("isch"))
                {
                    return false;
                }
                break;
            case 5:
                if (!r_R1())
                {
                    return false;
                }
                if (!r_C())
                {
                    return false;
                }
                if (!base.slice_del())
                {
                    return false;
                }
                break;
            case 6:
                if (!r_R1())
                {
                    return false;
                }
                if (!base.slice_from("t"))
                {
                    return false;
                }
                break;
            case 7:
                if (!r_R1())
                {
                    return false;
                }
                if (!base.slice_from("s"))
                {
                    return false;
                }
                break;
            case 8:
                if (!r_R1())
                {
                    return false;
                }
                if (!base.slice_from("r"))
                {
                    return false;
                }
                break;
            case 9:
                if (!r_R1())
                {
                    return false;
                }
                if (!base.slice_del())
                {
                    return false;
                }
                base.insert(base.cursor, base.cursor, "l");
                r_lengthen_V();
                break;
            case 10:
                if (!r_R1())
                {
                    return false;
                }
                if (!r_C())
                {
                    return false;
                }
                if (!base.slice_del())
                {
                    return false;
                }
                base.insert(base.cursor, base.cursor, "en");
                r_lengthen_V();
                break;
            case 11:
                if (!r_R1())
                {
                    return false;
                }
                if (!r_C())
                {
                    return false;
                }
                if (!base.slice_from("ief"))
                {
                    return false;
                }
                break;
        }
        return true;
    };

    /** @return {boolean} */
    function r_Step_3() {
        var /** number */ among_var;
        base.ket = base.cursor;
        among_var = base.find_among_b(a_3);
        if (among_var == 0)
        {
            return false;
        }
        base.bra = base.cursor;
        switch (among_var) {
            case 1:
                if (!r_R1())
                {
                    return false;
                }
                if (!base.slice_from("eer"))
                {
                    return false;
                }
                break;
            case 2:
                if (!r_R1())
                {
                    return false;
                }
                if (!base.slice_del())
                {
                    return false;
                }
                r_lengthen_V();
                break;
            case 3:
                if (!r_R1())
                {
                    return false;
                }
                if (!base.slice_del())
                {
                    return false;
                }
                break;
            case 4:
                if (!base.slice_from("r"))
                {
                    return false;
                }
                break;
            case 5:
                lab0: {
                    /** @const */ var /** number */ v_1 = base.limit - base.cursor;
                    lab1: {
                        if (!(base.eq_s_b("ild")))
                        {
                            break lab1;
                        }
                        if (!base.slice_from("er"))
                        {
                            return false;
                        }
                        break lab0;
                    }
                    base.cursor = base.limit - v_1;
                    if (!r_R1())
                    {
                        return false;
                    }
                    if (!base.slice_del())
                    {
                        return false;
                    }
                    r_lengthen_V();
                }
                break;
            case 6:
                if (!r_R1())
                {
                    return false;
                }
                if (!r_C())
                {
                    return false;
                }
                if (!base.slice_from("aar"))
                {
                    return false;
                }
                break;
            case 7:
                if (!r_R2())
                {
                    return false;
                }
                if (!base.slice_del())
                {
                    return false;
                }
                base.insert(base.cursor, base.cursor, "f");
                r_lengthen_V();
                break;
            case 8:
                if (!r_R2())
                {
                    return false;
                }
                if (!base.slice_del())
                {
                    return false;
                }
                base.insert(base.cursor, base.cursor, "g");
                r_lengthen_V();
                break;
            case 9:
                if (!r_R1())
                {
                    return false;
                }
                if (!r_C())
                {
                    return false;
                }
                if (!base.slice_from("t"))
                {
                    return false;
                }
                break;
            case 10:
                if (!r_R1())
                {
                    return false;
                }
                if (!r_C())
                {
                    return false;
                }
                if (!base.slice_from("d"))
                {
                    return false;
                }
                break;
        }
        return true;
    };

    /** @return {boolean} */
    function r_Step_4() {
        var /** number */ among_var;
        lab0: {
            /** @const */ var /** number */ v_1 = base.limit - base.cursor;
            lab1: {
                base.ket = base.cursor;
                among_var = base.find_among_b(a_4);
                if (among_var == 0)
                {
                    break lab1;
                }
                base.bra = base.cursor;
                switch (among_var) {
                    case 1:
                        if (!r_R1())
                        {
                            break lab1;
                        }
                        if (!base.slice_from("ie"))
                        {
                            return false;
                        }
                        break;
                    case 2:
                        if (!r_R1())
                        {
                            break lab1;
                        }
                        if (!base.slice_from("eer"))
                        {
                            return false;
                        }
                        break;
                    case 3:
                        if (!r_R1())
                        {
                            break lab1;
                        }
                        if (!base.slice_del())
                        {
                            return false;
                        }
                        break;
                    case 4:
                        if (!r_R1())
                        {
                            break lab1;
                        }
                        if (!r_V())
                        {
                            break lab1;
                        }
                        if (!base.slice_from("n"))
                        {
                            return false;
                        }
                        break;
                    case 5:
                        if (!r_R1())
                        {
                            break lab1;
                        }
                        if (!r_V())
                        {
                            break lab1;
                        }
                        if (!base.slice_from("l"))
                        {
                            return false;
                        }
                        break;
                    case 6:
                        if (!r_R1())
                        {
                            break lab1;
                        }
                        if (!r_V())
                        {
                            break lab1;
                        }
                        if (!base.slice_from("r"))
                        {
                            return false;
                        }
                        break;
                    case 7:
                        if (!r_R1())
                        {
                            break lab1;
                        }
                        if (!base.slice_from("teer"))
                        {
                            return false;
                        }
                        break;
                    case 8:
                        if (!r_R1())
                        {
                            break lab1;
                        }
                        if (!base.slice_from("lijk"))
                        {
                            return false;
                        }
                        break;
                    case 9:
                        if (!r_R1())
                        {
                            break lab1;
                        }
                        if (!r_C())
                        {
                            break lab1;
                        }
                        if (!base.slice_del())
                        {
                            return false;
                        }
                        r_lengthen_V();
                        break;
                }
                break lab0;
            }
            base.cursor = base.limit - v_1;
            base.ket = base.cursor;
            if (base.find_among_b(a_5) == 0)
            {
                return false;
            }
            base.bra = base.cursor;
            if (!r_R1())
            {
                return false;
            }
            {
                /** @const */ var /** number */ v_2 = base.limit - base.cursor;
                lab2: {
                    if (!(base.eq_s_b("inn")))
                    {
                        break lab2;
                    }
                    if (base.cursor > base.limit_backward)
                    {
                        break lab2;
                    }
                    return false;
                }
                base.cursor = base.limit - v_2;
            }
            if (!r_C())
            {
                return false;
            }
            if (!base.slice_del())
            {
                return false;
            }
            r_lengthen_V();
        }
        return true;
    };

    /** @return {boolean} */
    function r_Step_7() {
        var /** number */ among_var;
        base.ket = base.cursor;
        among_var = base.find_among_b(a_6);
        if (among_var == 0)
        {
            return false;
        }
        base.bra = base.cursor;
        switch (among_var) {
            case 1:
                if (!base.slice_from("k"))
                {
                    return false;
                }
                break;
            case 2:
                if (!base.slice_from("f"))
                {
                    return false;
                }
                break;
            case 3:
                if (!base.slice_from("p"))
                {
                    return false;
                }
                break;
        }
        return true;
    };

    /** @return {boolean} */
    function r_Step_6() {
        var /** number */ among_var;
        base.ket = base.cursor;
        among_var = base.find_among_b(a_7);
        if (among_var == 0)
        {
            return false;
        }
        base.bra = base.cursor;
        switch (among_var) {
            case 1:
                if (!base.slice_from("b"))
                {
                    return false;
                }
                break;
            case 2:
                if (!base.slice_from("c"))
                {
                    return false;
                }
                break;
            case 3:
                if (!base.slice_from("d"))
                {
                    return false;
                }
                break;
            case 4:
                if (!base.slice_from("f"))
                {
                    return false;
                }
                break;
            case 5:
                if (!base.slice_from("g"))
                {
                    return false;
                }
                break;
            case 6:
                if (!base.slice_from("h"))
                {
                    return false;
                }
                break;
            case 7:
                if (!base.slice_from("j"))
                {
                    return false;
                }
                break;
            case 8:
                if (!base.slice_from("k"))
                {
                    return false;
                }
                break;
            case 9:
                if (!base.slice_from("l"))
                {
                    return false;
                }
                break;
            case 10:
                if (!base.slice_from("m"))
                {
                    return false;
                }
                break;
            case 11:
                {
                    /** @const */ var /** number */ v_1 = base.limit - base.cursor;
                    lab0: {
                        if (!(base.eq_s_b("i")))
                        {
                            break lab0;
                        }
                        if (base.cursor > base.limit_backward)
                        {
                            break lab0;
                        }
                        return false;
                    }
                    base.cursor = base.limit - v_1;
                }
                if (!base.slice_from("n"))
                {
                    return false;
                }
                break;
            case 12:
                if (!base.slice_from("p"))
                {
                    return false;
                }
                break;
            case 13:
                if (!base.slice_from("q"))
                {
                    return false;
                }
                break;
            case 14:
                if (!base.slice_from("r"))
                {
                    return false;
                }
                break;
            case 15:
                if (!base.slice_from("s"))
                {
                    return false;
                }
                break;
            case 16:
                if (!base.slice_from("t"))
                {
                    return false;
                }
                break;
            case 17:
                if (!base.slice_from("v"))
                {
                    return false;
                }
                break;
            case 18:
                if (!base.slice_from("w"))
                {
                    return false;
                }
                break;
            case 19:
                if (!base.slice_from("x"))
                {
                    return false;
                }
                break;
            case 20:
                if (!base.slice_from("z"))
                {
                    return false;
                }
                break;
        }
        return true;
    };

    /** @return {boolean} */
    function r_Step_1c() {
        var /** number */ among_var;
        base.ket = base.cursor;
        among_var = base.find_among_b(a_8);
        if (among_var == 0)
        {
            return false;
        }
        base.bra = base.cursor;
        if (!r_R1())
        {
            return false;
        }
        if (!r_C())
        {
            return false;
        }
        switch (among_var) {
            case 1:
                {
                    /** @const */ var /** number */ v_1 = base.limit - base.cursor;
                    lab0: {
                        if (!(base.eq_s_b("n")))
                        {
                            break lab0;
                        }
                        if (!r_R1())
                        {
                            break lab0;
                        }
                        return false;
                    }
                    base.cursor = base.limit - v_1;
                }
                lab1: {
                    /** @const */ var /** number */ v_2 = base.limit - base.cursor;
                    lab2: {
                        if (!(base.eq_s_b("in")))
                        {
                            break lab2;
                        }
                        if (base.cursor > base.limit_backward)
                        {
                            break lab2;
                        }
                        if (!base.slice_from("n"))
                        {
                            return false;
                        }
                        break lab1;
                    }
                    base.cursor = base.limit - v_2;
                    if (!base.slice_del())
                    {
                        return false;
                    }
                }
                break;
            case 2:
                {
                    /** @const */ var /** number */ v_3 = base.limit - base.cursor;
                    lab3: {
                        if (!(base.eq_s_b("h")))
                        {
                            break lab3;
                        }
                        if (!r_R1())
                        {
                            break lab3;
                        }
                        return false;
                    }
                    base.cursor = base.limit - v_3;
                }
                {
                    /** @const */ var /** number */ v_4 = base.limit - base.cursor;
                    lab4: {
                        if (!(base.eq_s_b("en")))
                        {
                            break lab4;
                        }
                        if (base.cursor > base.limit_backward)
                        {
                            break lab4;
                        }
                        return false;
                    }
                    base.cursor = base.limit - v_4;
                }
                if (!base.slice_del())
                {
                    return false;
                }
                break;
        }
        return true;
    };

    /** @return {boolean} */
    function r_Lose_prefix() {
        var /** number */ among_var;
        base.bra = base.cursor;
        if (!(base.eq_s("ge")))
        {
            return false;
        }
        base.ket = base.cursor;
        /** @const */ var /** number */ v_1 = base.cursor;
        {
            /** @const */ var /** number */ c1 = base.cursor + 3;
            if (c1 > base.limit)
            {
                return false;
            }
            base.cursor = c1;
        }
        base.cursor = v_1;
        /** @const */ var /** number */ v_2 = base.cursor;
        golab0: while(true)
        {
            /** @const */ var /** number */ v_3 = base.cursor;
            lab1: {
                lab2: {
                    /** @const */ var /** number */ v_4 = base.cursor;
                    lab3: {
                        if (!(base.eq_s("ij")))
                        {
                            break lab3;
                        }
                        break lab2;
                    }
                    base.cursor = v_4;
                    if (!(base.in_grouping(g_v, 97, 252)))
                    {
                        break lab1;
                    }
                }
                break golab0;
            }
            base.cursor = v_3;
            if (base.cursor >= base.limit)
            {
                return false;
            }
            base.cursor++;
        }
        while(true)
        {
            /** @const */ var /** number */ v_5 = base.cursor;
            lab4: {
                lab5: {
                    /** @const */ var /** number */ v_6 = base.cursor;
                    lab6: {
                        if (!(base.eq_s("ij")))
                        {
                            break lab6;
                        }
                        break lab5;
                    }
                    base.cursor = v_6;
                    if (!(base.in_grouping(g_v, 97, 252)))
                    {
                        break lab4;
                    }
                }
                continue;
            }
            base.cursor = v_5;
            break;
        }
        lab7: {
            if (base.cursor < base.limit)
            {
                break lab7;
            }
            return false;
        }
        base.cursor = v_2;
        among_var = base.find_among(a_9);
        switch (among_var) {
            case 1:
                return false;
        }
        B_GE_removed = true;
        if (!base.slice_del())
        {
            return false;
        }
        /** @const */ var /** number */ v_7 = base.cursor;
        lab8: {
            base.bra = base.cursor;
            among_var = base.find_among(a_10);
            if (among_var == 0)
            {
                break lab8;
            }
            base.ket = base.cursor;
            switch (among_var) {
                case 1:
                    if (!base.slice_from("e"))
                    {
                        return false;
                    }
                    break;
                case 2:
                    if (!base.slice_from("i"))
                    {
                        return false;
                    }
                    break;
            }
        }
        base.cursor = v_7;
        return true;
    };

    /** @return {boolean} */
    function r_Lose_infix() {
        var /** number */ among_var;
        if (base.cursor >= base.limit)
        {
            return false;
        }
        base.cursor++;
        golab0: while(true)
        {
            lab1: {
                base.bra = base.cursor;
                if (!(base.eq_s("ge")))
                {
                    break lab1;
                }
                base.ket = base.cursor;
                break golab0;
            }
            if (base.cursor >= base.limit)
            {
                return false;
            }
            base.cursor++;
        }
        /** @const */ var /** number */ v_1 = base.cursor;
        {
            /** @const */ var /** number */ c1 = base.cursor + 3;
            if (c1 > base.limit)
            {
                return false;
            }
            base.cursor = c1;
        }
        base.cursor = v_1;
        /** @const */ var /** number */ v_2 = base.cursor;
        golab2: while(true)
        {
            /** @const */ var /** number */ v_3 = base.cursor;
            lab3: {
                lab4: {
                    /** @const */ var /** number */ v_4 = base.cursor;
                    lab5: {
                        if (!(base.eq_s("ij")))
                        {
                            break lab5;
                        }
                        break lab4;
                    }
                    base.cursor = v_4;
                    if (!(base.in_grouping(g_v, 97, 252)))
                    {
                        break lab3;
                    }
                }
                break golab2;
            }
            base.cursor = v_3;
            if (base.cursor >= base.limit)
            {
                return false;
            }
            base.cursor++;
        }
        while(true)
        {
            /** @const */ var /** number */ v_5 = base.cursor;
            lab6: {
                lab7: {
                    /** @const */ var /** number */ v_6 = base.cursor;
                    lab8: {
                        if (!(base.eq_s("ij")))
                        {
                            break lab8;
                        }
                        break lab7;
                    }
                    base.cursor = v_6;
                    if (!(base.in_grouping(g_v, 97, 252)))
                    {
                        break lab6;
                    }
                }
                continue;
            }
            base.cursor = v_5;
            break;
        }
        lab9: {
            if (base.cursor < base.limit)
            {
                break lab9;
            }
            return false;
        }
        base.cursor = v_2;
        B_GE_removed = true;
        if (!base.slice_del())
        {
            return false;
        }
        /** @const */ var /** number */ v_7 = base.cursor;
        lab10: {
            base.bra = base.cursor;
            among_var = base.find_among(a_11);
            if (among_var == 0)
            {
                break lab10;
            }
            base.ket = base.cursor;
            switch (among_var) {
                case 1:
                    if (!base.slice_from("e"))
                    {
                        return false;
                    }
                    break;
                case 2:
                    if (!base.slice_from("i"))
                    {
                        return false;
                    }
                    break;
            }
        }
        base.cursor = v_7;
        return true;
    };

    /** @return {boolean} */
    function r_measure() {
        I_p1 = base.limit;
        I_p2 = base.limit;
        /** @const */ var /** number */ v_1 = base.cursor;
        lab0: {
            while(true)
            {
                lab1: {
                    if (!(base.out_grouping(g_v, 97, 252)))
                    {
                        break lab1;
                    }
                    continue;
                }
                break;
            }
            {
                var v_2 = 1;
                while(true)
                {
                    /** @const */ var /** number */ v_3 = base.cursor;
                    lab2: {
                        lab3: {
                            /** @const */ var /** number */ v_4 = base.cursor;
                            lab4: {
                                if (!(base.eq_s("ij")))
                                {
                                    break lab4;
                                }
                                break lab3;
                            }
                            base.cursor = v_4;
                            if (!(base.in_grouping(g_v, 97, 252)))
                            {
                                break lab2;
                            }
                        }
                        v_2--;
                        continue;
                    }
                    base.cursor = v_3;
                    break;
                }
                if (v_2 > 0)
                {
                    break lab0;
                }
            }
            if (!(base.out_grouping(g_v, 97, 252)))
            {
                break lab0;
            }
            I_p1 = base.cursor;
            while(true)
            {
                lab5: {
                    if (!(base.out_grouping(g_v, 97, 252)))
                    {
                        break lab5;
                    }
                    continue;
                }
                break;
            }
            {
                var v_5 = 1;
                while(true)
                {
                    /** @const */ var /** number */ v_6 = base.cursor;
                    lab6: {
                        lab7: {
                            /** @const */ var /** number */ v_7 = base.cursor;
                            lab8: {
                                if (!(base.eq_s("ij")))
                                {
                                    break lab8;
                                }
                                break lab7;
                            }
                            base.cursor = v_7;
                            if (!(base.in_grouping(g_v, 97, 252)))
                            {
                                break lab6;
                            }
                        }
                        v_5--;
                        continue;
                    }
                    base.cursor = v_6;
                    break;
                }
                if (v_5 > 0)
                {
                    break lab0;
                }
            }
            if (!(base.out_grouping(g_v, 97, 252)))
            {
                break lab0;
            }
            I_p2 = base.cursor;
        }
        base.cursor = v_1;
        return true;
    };

    this.stem = /** @return {boolean} */ function() {
        B_stemmed = false;
        r_measure();
        base.limit_backward = base.cursor; base.cursor = base.limit;
        /** @const */ var /** number */ v_1 = base.limit - base.cursor;
        lab0: {
            if (!r_Step_1())
            {
                break lab0;
            }
            B_stemmed = true;
        }
        base.cursor = base.limit - v_1;
        /** @const */ var /** number */ v_2 = base.limit - base.cursor;
        lab1: {
            if (!r_Step_2())
            {
                break lab1;
            }
            B_stemmed = true;
        }
        base.cursor = base.limit - v_2;
        /** @const */ var /** number */ v_3 = base.limit - base.cursor;
        lab2: {
            if (!r_Step_3())
            {
                break lab2;
            }
            B_stemmed = true;
        }
        base.cursor = base.limit - v_3;
        /** @const */ var /** number */ v_4 = base.limit - base.cursor;
        lab3: {
            if (!r_Step_4())
            {
                break lab3;
            }
            B_stemmed = true;
        }
        base.cursor = base.limit - v_4;
        base.cursor = base.limit_backward;
        B_GE_removed = false;
        /** @const */ var /** number */ v_5 = base.cursor;
        lab4: {
            /** @const */ var /** number */ v_6 = base.cursor;
            if (!r_Lose_prefix())
            {
                break lab4;
            }
            base.cursor = v_6;
            r_measure();
        }
        base.cursor = v_5;
        base.limit_backward = base.cursor; base.cursor = base.limit;
        /** @const */ var /** number */ v_7 = base.limit - base.cursor;
        lab5: {
            if (!B_GE_removed)
            {
                break lab5;
            }
            B_stemmed = true;
            if (!r_Step_1c())
            {
                break lab5;
            }
        }
        base.cursor = base.limit - v_7;
        base.cursor = base.limit_backward;
        B_GE_removed = false;
        /** @const */ var /** number */ v_8 = base.cursor;
        lab6: {
            /** @const */ var /** number */ v_9 = base.cursor;
            if (!r_Lose_infix())
            {
                break lab6;
            }
            base.cursor = v_9;
            r_measure();
        }
        base.cursor = v_8;
        base.limit_backward = base.cursor; base.cursor = base.limit;
        /** @const */ var /** number */ v_10 = base.limit - base.cursor;
        lab7: {
            if (!B_GE_removed)
            {
                break lab7;
            }
            B_stemmed = true;
            if (!r_Step_1c())
            {
                break lab7;
            }
        }
        base.cursor = base.limit - v_10;
        base.cursor = base.limit_backward;
        base.limit_backward = base.cursor; base.cursor = base.limit;
        /** @const */ var /** number */ v_11 = base.limit - base.cursor;
        lab8: {
            if (!r_Step_7())
            {
                break lab8;
            }
            B_stemmed = true;
        }
        base.cursor = base.limit - v_11;
        /** @const */ var /** number */ v_12 = base.limit - base.cursor;
        lab9: {
            if (!B_stemmed)
            {
                break lab9;
            }
            if (!r_Step_6())
            {
                break lab9;
            }
        }
        base.cursor = base.limit - v_12;
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
