// Generated from german.sbl by Snowball 3.0.1 - https://snowballstem.org/

/**@constructor*/
var GermanStemmer = function() {
    var base = new BaseStemmer();

    /** @const */ var a_0 = [
        ["", -1, 5],
        ["ae", 0, 2],
        ["oe", 0, 3],
        ["qu", 0, -1],
        ["ue", 0, 4],
        ["\u00DF", 0, 1]
    ];

    /** @const */ var a_1 = [
        ["", -1, 5],
        ["U", 0, 2],
        ["Y", 0, 1],
        ["\u00E4", 0, 3],
        ["\u00F6", 0, 4],
        ["\u00FC", 0, 2]
    ];

    /** @const */ var a_2 = [
        ["e", -1, 3],
        ["em", -1, 1],
        ["en", -1, 3],
        ["erinnen", 2, 2],
        ["erin", -1, 2],
        ["ln", -1, 5],
        ["ern", -1, 2],
        ["er", -1, 2],
        ["s", -1, 4],
        ["es", 8, 3],
        ["lns", 8, 5]
    ];

    /** @const */ var a_3 = [
        ["tick", -1, -1],
        ["plan", -1, -1],
        ["geordn", -1, -1],
        ["intern", -1, -1],
        ["tr", -1, -1]
    ];

    /** @const */ var a_4 = [
        ["en", -1, 1],
        ["er", -1, 1],
        ["et", -1, 3],
        ["st", -1, 2],
        ["est", 3, 1]
    ];

    /** @const */ var a_5 = [
        ["ig", -1, 1],
        ["lich", -1, 1]
    ];

    /** @const */ var a_6 = [
        ["end", -1, 1],
        ["ig", -1, 2],
        ["ung", -1, 1],
        ["lich", -1, 3],
        ["isch", -1, 2],
        ["ik", -1, 2],
        ["heit", -1, 3],
        ["keit", -1, 4]
    ];

    /** @const */ var /** Array<int> */ g_v = [17, 65, 16, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 0, 32, 8];

    /** @const */ var /** Array<int> */ g_et_ending = [1, 128, 198, 227, 32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 128];

    /** @const */ var /** Array<int> */ g_s_ending = [117, 30, 5];

    /** @const */ var /** Array<int> */ g_st_ending = [117, 30, 4];

    var /** number */ I_x = 0;
    var /** number */ I_p2 = 0;
    var /** number */ I_p1 = 0;


    /** @return {boolean} */
    function r_prelude() {
        var /** number */ among_var;
        /** @const */ var /** number */ v_1 = base.cursor;
        while(true)
        {
            /** @const */ var /** number */ v_2 = base.cursor;
            lab0: {
                golab1: while(true)
                {
                    /** @const */ var /** number */ v_3 = base.cursor;
                    lab2: {
                        if (!(base.in_grouping(g_v, 97, 252)))
                        {
                            break lab2;
                        }
                        base.bra = base.cursor;
                        lab3: {
                            /** @const */ var /** number */ v_4 = base.cursor;
                            lab4: {
                                if (!(base.eq_s("u")))
                                {
                                    break lab4;
                                }
                                base.ket = base.cursor;
                                if (!(base.in_grouping(g_v, 97, 252)))
                                {
                                    break lab4;
                                }
                                if (!base.slice_from("U"))
                                {
                                    return false;
                                }
                                break lab3;
                            }
                            base.cursor = v_4;
                            if (!(base.eq_s("y")))
                            {
                                break lab2;
                            }
                            base.ket = base.cursor;
                            if (!(base.in_grouping(g_v, 97, 252)))
                            {
                                break lab2;
                            }
                            if (!base.slice_from("Y"))
                            {
                                return false;
                            }
                        }
                        base.cursor = v_3;
                        break golab1;
                    }
                    base.cursor = v_3;
                    if (base.cursor >= base.limit)
                    {
                        break lab0;
                    }
                    base.cursor++;
                }
                continue;
            }
            base.cursor = v_2;
            break;
        }
        base.cursor = v_1;
        while(true)
        {
            /** @const */ var /** number */ v_5 = base.cursor;
            lab5: {
                base.bra = base.cursor;
                among_var = base.find_among(a_0);
                base.ket = base.cursor;
                switch (among_var) {
                    case 1:
                        if (!base.slice_from("ss"))
                        {
                            return false;
                        }
                        break;
                    case 2:
                        if (!base.slice_from("\u00E4"))
                        {
                            return false;
                        }
                        break;
                    case 3:
                        if (!base.slice_from("\u00F6"))
                        {
                            return false;
                        }
                        break;
                    case 4:
                        if (!base.slice_from("\u00FC"))
                        {
                            return false;
                        }
                        break;
                    case 5:
                        if (base.cursor >= base.limit)
                        {
                            break lab5;
                        }
                        base.cursor++;
                        break;
                }
                continue;
            }
            base.cursor = v_5;
            break;
        }
        return true;
    };

    /** @return {boolean} */
    function r_mark_regions() {
        I_p1 = base.limit;
        I_p2 = base.limit;
        /** @const */ var /** number */ v_1 = base.cursor;
        {
            /** @const */ var /** number */ c1 = base.cursor + 3;
            if (c1 > base.limit)
            {
                return false;
            }
            base.cursor = c1;
        }
        I_x = base.cursor;
        base.cursor = v_1;
        if (!base.go_out_grouping(g_v, 97, 252))
        {
            return false;
        }
        base.cursor++;
        if (!base.go_in_grouping(g_v, 97, 252))
        {
            return false;
        }
        base.cursor++;
        I_p1 = base.cursor;
        lab0: {
            if (I_p1 >= I_x)
            {
                break lab0;
            }
            I_p1 = I_x;
        }
        if (!base.go_out_grouping(g_v, 97, 252))
        {
            return false;
        }
        base.cursor++;
        if (!base.go_in_grouping(g_v, 97, 252))
        {
            return false;
        }
        base.cursor++;
        I_p2 = base.cursor;
        return true;
    };

    /** @return {boolean} */
    function r_postlude() {
        var /** number */ among_var;
        while(true)
        {
            /** @const */ var /** number */ v_1 = base.cursor;
            lab0: {
                base.bra = base.cursor;
                among_var = base.find_among(a_1);
                base.ket = base.cursor;
                switch (among_var) {
                    case 1:
                        if (!base.slice_from("y"))
                        {
                            return false;
                        }
                        break;
                    case 2:
                        if (!base.slice_from("u"))
                        {
                            return false;
                        }
                        break;
                    case 3:
                        if (!base.slice_from("a"))
                        {
                            return false;
                        }
                        break;
                    case 4:
                        if (!base.slice_from("o"))
                        {
                            return false;
                        }
                        break;
                    case 5:
                        if (base.cursor >= base.limit)
                        {
                            break lab0;
                        }
                        base.cursor++;
                        break;
                }
                continue;
            }
            base.cursor = v_1;
            break;
        }
        return true;
    };

    /** @return {boolean} */
    function r_R1() {
        return I_p1 <= base.cursor;
    };

    /** @return {boolean} */
    function r_R2() {
        return I_p2 <= base.cursor;
    };

    /** @return {boolean} */
    function r_standard_suffix() {
        var /** number */ among_var;
        /** @const */ var /** number */ v_1 = base.limit - base.cursor;
        lab0: {
            base.ket = base.cursor;
            among_var = base.find_among_b(a_2);
            if (among_var == 0)
            {
                break lab0;
            }
            base.bra = base.cursor;
            if (!r_R1())
            {
                break lab0;
            }
            switch (among_var) {
                case 1:
                    {
                        /** @const */ var /** number */ v_2 = base.limit - base.cursor;
                        lab1: {
                            if (!(base.eq_s_b("syst")))
                            {
                                break lab1;
                            }
                            break lab0;
                        }
                        base.cursor = base.limit - v_2;
                    }
                    if (!base.slice_del())
                    {
                        return false;
                    }
                    break;
                case 2:
                    if (!base.slice_del())
                    {
                        return false;
                    }
                    break;
                case 3:
                    if (!base.slice_del())
                    {
                        return false;
                    }
                    /** @const */ var /** number */ v_3 = base.limit - base.cursor;
                    lab2: {
                        base.ket = base.cursor;
                        if (!(base.eq_s_b("s")))
                        {
                            base.cursor = base.limit - v_3;
                            break lab2;
                        }
                        base.bra = base.cursor;
                        if (!(base.eq_s_b("nis")))
                        {
                            base.cursor = base.limit - v_3;
                            break lab2;
                        }
                        if (!base.slice_del())
                        {
                            return false;
                        }
                    }
                    break;
                case 4:
                    if (!(base.in_grouping_b(g_s_ending, 98, 116)))
                    {
                        break lab0;
                    }
                    if (!base.slice_del())
                    {
                        return false;
                    }
                    break;
                case 5:
                    if (!base.slice_from("l"))
                    {
                        return false;
                    }
                    break;
            }
        }
        base.cursor = base.limit - v_1;
        /** @const */ var /** number */ v_4 = base.limit - base.cursor;
        lab3: {
            base.ket = base.cursor;
            among_var = base.find_among_b(a_4);
            if (among_var == 0)
            {
                break lab3;
            }
            base.bra = base.cursor;
            if (!r_R1())
            {
                break lab3;
            }
            switch (among_var) {
                case 1:
                    if (!base.slice_del())
                    {
                        return false;
                    }
                    break;
                case 2:
                    if (!(base.in_grouping_b(g_st_ending, 98, 116)))
                    {
                        break lab3;
                    }
                    {
                        /** @const */ var /** number */ c1 = base.cursor - 3;
                        if (c1 < base.limit_backward)
                        {
                            break lab3;
                        }
                        base.cursor = c1;
                    }
                    if (!base.slice_del())
                    {
                        return false;
                    }
                    break;
                case 3:
                    /** @const */ var /** number */ v_5 = base.limit - base.cursor;
                    if (!(base.in_grouping_b(g_et_ending, 85, 228)))
                    {
                        break lab3;
                    }
                    base.cursor = base.limit - v_5;
                    {
                        /** @const */ var /** number */ v_6 = base.limit - base.cursor;
                        lab4: {
                            if (base.find_among_b(a_3) == 0)
                            {
                                break lab4;
                            }
                            break lab3;
                        }
                        base.cursor = base.limit - v_6;
                    }
                    if (!base.slice_del())
                    {
                        return false;
                    }
                    break;
            }
        }
        base.cursor = base.limit - v_4;
        /** @const */ var /** number */ v_7 = base.limit - base.cursor;
        lab5: {
            base.ket = base.cursor;
            among_var = base.find_among_b(a_6);
            if (among_var == 0)
            {
                break lab5;
            }
            base.bra = base.cursor;
            if (!r_R2())
            {
                break lab5;
            }
            switch (among_var) {
                case 1:
                    if (!base.slice_del())
                    {
                        return false;
                    }
                    /** @const */ var /** number */ v_8 = base.limit - base.cursor;
                    lab6: {
                        base.ket = base.cursor;
                        if (!(base.eq_s_b("ig")))
                        {
                            base.cursor = base.limit - v_8;
                            break lab6;
                        }
                        base.bra = base.cursor;
                        {
                            /** @const */ var /** number */ v_9 = base.limit - base.cursor;
                            lab7: {
                                if (!(base.eq_s_b("e")))
                                {
                                    break lab7;
                                }
                                base.cursor = base.limit - v_8;
                                break lab6;
                            }
                            base.cursor = base.limit - v_9;
                        }
                        if (!r_R2())
                        {
                            base.cursor = base.limit - v_8;
                            break lab6;
                        }
                        if (!base.slice_del())
                        {
                            return false;
                        }
                    }
                    break;
                case 2:
                    {
                        /** @const */ var /** number */ v_10 = base.limit - base.cursor;
                        lab8: {
                            if (!(base.eq_s_b("e")))
                            {
                                break lab8;
                            }
                            break lab5;
                        }
                        base.cursor = base.limit - v_10;
                    }
                    if (!base.slice_del())
                    {
                        return false;
                    }
                    break;
                case 3:
                    if (!base.slice_del())
                    {
                        return false;
                    }
                    /** @const */ var /** number */ v_11 = base.limit - base.cursor;
                    lab9: {
                        base.ket = base.cursor;
                        lab10: {
                            /** @const */ var /** number */ v_12 = base.limit - base.cursor;
                            lab11: {
                                if (!(base.eq_s_b("er")))
                                {
                                    break lab11;
                                }
                                break lab10;
                            }
                            base.cursor = base.limit - v_12;
                            if (!(base.eq_s_b("en")))
                            {
                                base.cursor = base.limit - v_11;
                                break lab9;
                            }
                        }
                        base.bra = base.cursor;
                        if (!r_R1())
                        {
                            base.cursor = base.limit - v_11;
                            break lab9;
                        }
                        if (!base.slice_del())
                        {
                            return false;
                        }
                    }
                    break;
                case 4:
                    if (!base.slice_del())
                    {
                        return false;
                    }
                    /** @const */ var /** number */ v_13 = base.limit - base.cursor;
                    lab12: {
                        base.ket = base.cursor;
                        if (base.find_among_b(a_5) == 0)
                        {
                            base.cursor = base.limit - v_13;
                            break lab12;
                        }
                        base.bra = base.cursor;
                        if (!r_R2())
                        {
                            base.cursor = base.limit - v_13;
                            break lab12;
                        }
                        if (!base.slice_del())
                        {
                            return false;
                        }
                    }
                    break;
            }
        }
        base.cursor = base.limit - v_7;
        return true;
    };

    this.stem = /** @return {boolean} */ function() {
        /** @const */ var /** number */ v_1 = base.cursor;
        r_prelude();
        base.cursor = v_1;
        /** @const */ var /** number */ v_2 = base.cursor;
        r_mark_regions();
        base.cursor = v_2;
        base.limit_backward = base.cursor; base.cursor = base.limit;
        r_standard_suffix();
        base.cursor = base.limit_backward;
        /** @const */ var /** number */ v_3 = base.cursor;
        r_postlude();
        base.cursor = v_3;
        return true;
    };

    /**@return{string}*/
    this['stemWord'] = function(/**string*/word) {
        base.setCurrent(word);
        this.stem();
        return base.getCurrent();
    };
};
