// Generated from arabic.sbl by Snowball 3.0.1 - https://snowballstem.org/

/**@constructor*/
var ArabicStemmer = function() {
    var base = new BaseStemmer();

    /** @const */ var a_0 = [
        ["\u0640", -1, 1],
        ["\u064B", -1, 1],
        ["\u064C", -1, 1],
        ["\u064D", -1, 1],
        ["\u064E", -1, 1],
        ["\u064F", -1, 1],
        ["\u0650", -1, 1],
        ["\u0651", -1, 1],
        ["\u0652", -1, 1],
        ["\u0660", -1, 2],
        ["\u0661", -1, 3],
        ["\u0662", -1, 4],
        ["\u0663", -1, 5],
        ["\u0664", -1, 6],
        ["\u0665", -1, 7],
        ["\u0666", -1, 8],
        ["\u0667", -1, 9],
        ["\u0668", -1, 10],
        ["\u0669", -1, 11],
        ["\uFE80", -1, 12],
        ["\uFE81", -1, 16],
        ["\uFE82", -1, 16],
        ["\uFE83", -1, 13],
        ["\uFE84", -1, 13],
        ["\uFE85", -1, 17],
        ["\uFE86", -1, 17],
        ["\uFE87", -1, 14],
        ["\uFE88", -1, 14],
        ["\uFE89", -1, 15],
        ["\uFE8A", -1, 15],
        ["\uFE8B", -1, 15],
        ["\uFE8C", -1, 15],
        ["\uFE8D", -1, 18],
        ["\uFE8E", -1, 18],
        ["\uFE8F", -1, 19],
        ["\uFE90", -1, 19],
        ["\uFE91", -1, 19],
        ["\uFE92", -1, 19],
        ["\uFE93", -1, 20],
        ["\uFE94", -1, 20],
        ["\uFE95", -1, 21],
        ["\uFE96", -1, 21],
        ["\uFE97", -1, 21],
        ["\uFE98", -1, 21],
        ["\uFE99", -1, 22],
        ["\uFE9A", -1, 22],
        ["\uFE9B", -1, 22],
        ["\uFE9C", -1, 22],
        ["\uFE9D", -1, 23],
        ["\uFE9E", -1, 23],
        ["\uFE9F", -1, 23],
        ["\uFEA0", -1, 23],
        ["\uFEA1", -1, 24],
        ["\uFEA2", -1, 24],
        ["\uFEA3", -1, 24],
        ["\uFEA4", -1, 24],
        ["\uFEA5", -1, 25],
        ["\uFEA6", -1, 25],
        ["\uFEA7", -1, 25],
        ["\uFEA8", -1, 25],
        ["\uFEA9", -1, 26],
        ["\uFEAA", -1, 26],
        ["\uFEAB", -1, 27],
        ["\uFEAC", -1, 27],
        ["\uFEAD", -1, 28],
        ["\uFEAE", -1, 28],
        ["\uFEAF", -1, 29],
        ["\uFEB0", -1, 29],
        ["\uFEB1", -1, 30],
        ["\uFEB2", -1, 30],
        ["\uFEB3", -1, 30],
        ["\uFEB4", -1, 30],
        ["\uFEB5", -1, 31],
        ["\uFEB6", -1, 31],
        ["\uFEB7", -1, 31],
        ["\uFEB8", -1, 31],
        ["\uFEB9", -1, 32],
        ["\uFEBA", -1, 32],
        ["\uFEBB", -1, 32],
        ["\uFEBC", -1, 32],
        ["\uFEBD", -1, 33],
        ["\uFEBE", -1, 33],
        ["\uFEBF", -1, 33],
        ["\uFEC0", -1, 33],
        ["\uFEC1", -1, 34],
        ["\uFEC2", -1, 34],
        ["\uFEC3", -1, 34],
        ["\uFEC4", -1, 34],
        ["\uFEC5", -1, 35],
        ["\uFEC6", -1, 35],
        ["\uFEC7", -1, 35],
        ["\uFEC8", -1, 35],
        ["\uFEC9", -1, 36],
        ["\uFECA", -1, 36],
        ["\uFECB", -1, 36],
        ["\uFECC", -1, 36],
        ["\uFECD", -1, 37],
        ["\uFECE", -1, 37],
        ["\uFECF", -1, 37],
        ["\uFED0", -1, 37],
        ["\uFED1", -1, 38],
        ["\uFED2", -1, 38],
        ["\uFED3", -1, 38],
        ["\uFED4", -1, 38],
        ["\uFED5", -1, 39],
        ["\uFED6", -1, 39],
        ["\uFED7", -1, 39],
        ["\uFED8", -1, 39],
        ["\uFED9", -1, 40],
        ["\uFEDA", -1, 40],
        ["\uFEDB", -1, 40],
        ["\uFEDC", -1, 40],
        ["\uFEDD", -1, 41],
        ["\uFEDE", -1, 41],
        ["\uFEDF", -1, 41],
        ["\uFEE0", -1, 41],
        ["\uFEE1", -1, 42],
        ["\uFEE2", -1, 42],
        ["\uFEE3", -1, 42],
        ["\uFEE4", -1, 42],
        ["\uFEE5", -1, 43],
        ["\uFEE6", -1, 43],
        ["\uFEE7", -1, 43],
        ["\uFEE8", -1, 43],
        ["\uFEE9", -1, 44],
        ["\uFEEA", -1, 44],
        ["\uFEEB", -1, 44],
        ["\uFEEC", -1, 44],
        ["\uFEED", -1, 45],
        ["\uFEEE", -1, 45],
        ["\uFEEF", -1, 46],
        ["\uFEF0", -1, 46],
        ["\uFEF1", -1, 47],
        ["\uFEF2", -1, 47],
        ["\uFEF3", -1, 47],
        ["\uFEF4", -1, 47],
        ["\uFEF5", -1, 51],
        ["\uFEF6", -1, 51],
        ["\uFEF7", -1, 49],
        ["\uFEF8", -1, 49],
        ["\uFEF9", -1, 50],
        ["\uFEFA", -1, 50],
        ["\uFEFB", -1, 48],
        ["\uFEFC", -1, 48]
    ];

    /** @const */ var a_1 = [
        ["\u0622", -1, 1],
        ["\u0623", -1, 1],
        ["\u0624", -1, 1],
        ["\u0625", -1, 1],
        ["\u0626", -1, 1]
    ];

    /** @const */ var a_2 = [
        ["\u0622", -1, 1],
        ["\u0623", -1, 1],
        ["\u0624", -1, 2],
        ["\u0625", -1, 1],
        ["\u0626", -1, 3]
    ];

    /** @const */ var a_3 = [
        ["\u0627\u0644", -1, 2],
        ["\u0628\u0627\u0644", -1, 1],
        ["\u0643\u0627\u0644", -1, 1],
        ["\u0644\u0644", -1, 2]
    ];

    /** @const */ var a_4 = [
        ["\u0623\u0622", -1, 2],
        ["\u0623\u0623", -1, 1],
        ["\u0623\u0624", -1, 1],
        ["\u0623\u0625", -1, 4],
        ["\u0623\u0627", -1, 3]
    ];

    /** @const */ var a_5 = [
        ["\u0641", -1, 1],
        ["\u0648", -1, 1]
    ];

    /** @const */ var a_6 = [
        ["\u0627\u0644", -1, 2],
        ["\u0628\u0627\u0644", -1, 1],
        ["\u0643\u0627\u0644", -1, 1],
        ["\u0644\u0644", -1, 2]
    ];

    /** @const */ var a_7 = [
        ["\u0628", -1, 1],
        ["\u0628\u0627", 0, -1],
        ["\u0628\u0628", 0, 2],
        ["\u0643\u0643", -1, 3]
    ];

    /** @const */ var a_8 = [
        ["\u0633\u0623", -1, 4],
        ["\u0633\u062A", -1, 2],
        ["\u0633\u0646", -1, 3],
        ["\u0633\u064A", -1, 1]
    ];

    /** @const */ var a_9 = [
        ["\u062A\u0633\u062A", -1, 1],
        ["\u0646\u0633\u062A", -1, 1],
        ["\u064A\u0633\u062A", -1, 1]
    ];

    /** @const */ var a_10 = [
        ["\u0643\u0645\u0627", -1, 3],
        ["\u0647\u0645\u0627", -1, 3],
        ["\u0646\u0627", -1, 2],
        ["\u0647\u0627", -1, 2],
        ["\u0643", -1, 1],
        ["\u0643\u0645", -1, 2],
        ["\u0647\u0645", -1, 2],
        ["\u0647\u0646", -1, 2],
        ["\u0647", -1, 1],
        ["\u064A", -1, 1]
    ];

    /** @const */ var a_11 = [
        ["\u0646", -1, 1]
    ];

    /** @const */ var a_12 = [
        ["\u0627", -1, 1],
        ["\u0648", -1, 1],
        ["\u064A", -1, 1]
    ];

    /** @const */ var a_13 = [
        ["\u0627\u062A", -1, 1]
    ];

    /** @const */ var a_14 = [
        ["\u062A", -1, 1]
    ];

    /** @const */ var a_15 = [
        ["\u0629", -1, 1]
    ];

    /** @const */ var a_16 = [
        ["\u064A", -1, 1]
    ];

    /** @const */ var a_17 = [
        ["\u0643\u0645\u0627", -1, 3],
        ["\u0647\u0645\u0627", -1, 3],
        ["\u0646\u0627", -1, 2],
        ["\u0647\u0627", -1, 2],
        ["\u0643", -1, 1],
        ["\u0643\u0645", -1, 2],
        ["\u0647\u0645", -1, 2],
        ["\u0643\u0646", -1, 2],
        ["\u0647\u0646", -1, 2],
        ["\u0647", -1, 1],
        ["\u0643\u0645\u0648", -1, 3],
        ["\u0646\u064A", -1, 2]
    ];

    /** @const */ var a_18 = [
        ["\u0627", -1, 1],
        ["\u062A\u0627", 0, 2],
        ["\u062A\u0645\u0627", 0, 4],
        ["\u0646\u0627", 0, 2],
        ["\u062A", -1, 1],
        ["\u0646", -1, 1],
        ["\u0627\u0646", 5, 3],
        ["\u062A\u0646", 5, 2],
        ["\u0648\u0646", 5, 3],
        ["\u064A\u0646", 5, 3],
        ["\u064A", -1, 1]
    ];

    /** @const */ var a_19 = [
        ["\u0648\u0627", -1, 1],
        ["\u062A\u0645", -1, 1]
    ];

    /** @const */ var a_20 = [
        ["\u0648", -1, 1],
        ["\u062A\u0645\u0648", 0, 2]
    ];

    /** @const */ var a_21 = [
        ["\u0649", -1, 1]
    ];

    var /** boolean */ B_is_defined = false;
    var /** boolean */ B_is_verb = false;
    var /** boolean */ B_is_noun = false;


    /** @return {boolean} */
    function r_Normalize_pre() {
        var /** number */ among_var;
        /** @const */ var /** number */ v_1 = base.cursor;
        lab0: {
            while(true)
            {
                /** @const */ var /** number */ v_2 = base.cursor;
                lab1: {
                    lab2: {
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
                                    if (!base.slice_del())
                                    {
                                        return false;
                                    }
                                    break;
                                case 2:
                                    if (!base.slice_from("0"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 3:
                                    if (!base.slice_from("1"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 4:
                                    if (!base.slice_from("2"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 5:
                                    if (!base.slice_from("3"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 6:
                                    if (!base.slice_from("4"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 7:
                                    if (!base.slice_from("5"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 8:
                                    if (!base.slice_from("6"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 9:
                                    if (!base.slice_from("7"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 10:
                                    if (!base.slice_from("8"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 11:
                                    if (!base.slice_from("9"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 12:
                                    if (!base.slice_from("\u0621"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 13:
                                    if (!base.slice_from("\u0623"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 14:
                                    if (!base.slice_from("\u0625"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 15:
                                    if (!base.slice_from("\u0626"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 16:
                                    if (!base.slice_from("\u0622"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 17:
                                    if (!base.slice_from("\u0624"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 18:
                                    if (!base.slice_from("\u0627"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 19:
                                    if (!base.slice_from("\u0628"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 20:
                                    if (!base.slice_from("\u0629"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 21:
                                    if (!base.slice_from("\u062A"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 22:
                                    if (!base.slice_from("\u062B"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 23:
                                    if (!base.slice_from("\u062C"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 24:
                                    if (!base.slice_from("\u062D"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 25:
                                    if (!base.slice_from("\u062E"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 26:
                                    if (!base.slice_from("\u062F"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 27:
                                    if (!base.slice_from("\u0630"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 28:
                                    if (!base.slice_from("\u0631"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 29:
                                    if (!base.slice_from("\u0632"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 30:
                                    if (!base.slice_from("\u0633"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 31:
                                    if (!base.slice_from("\u0634"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 32:
                                    if (!base.slice_from("\u0635"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 33:
                                    if (!base.slice_from("\u0636"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 34:
                                    if (!base.slice_from("\u0637"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 35:
                                    if (!base.slice_from("\u0638"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 36:
                                    if (!base.slice_from("\u0639"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 37:
                                    if (!base.slice_from("\u063A"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 38:
                                    if (!base.slice_from("\u0641"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 39:
                                    if (!base.slice_from("\u0642"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 40:
                                    if (!base.slice_from("\u0643"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 41:
                                    if (!base.slice_from("\u0644"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 42:
                                    if (!base.slice_from("\u0645"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 43:
                                    if (!base.slice_from("\u0646"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 44:
                                    if (!base.slice_from("\u0647"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 45:
                                    if (!base.slice_from("\u0648"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 46:
                                    if (!base.slice_from("\u0649"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 47:
                                    if (!base.slice_from("\u064A"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 48:
                                    if (!base.slice_from("\u0644\u0627"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 49:
                                    if (!base.slice_from("\u0644\u0623"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 50:
                                    if (!base.slice_from("\u0644\u0625"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 51:
                                    if (!base.slice_from("\u0644\u0622"))
                                    {
                                        return false;
                                    }
                                    break;
                            }
                            break lab2;
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
        return true;
    };

    /** @return {boolean} */
    function r_Normalize_post() {
        var /** number */ among_var;
        /** @const */ var /** number */ v_1 = base.cursor;
        lab0: {
            base.limit_backward = base.cursor; base.cursor = base.limit;
            base.ket = base.cursor;
            if (base.find_among_b(a_1) == 0)
            {
                break lab0;
            }
            base.bra = base.cursor;
            if (!base.slice_from("\u0621"))
            {
                return false;
            }
            base.cursor = base.limit_backward;
        }
        base.cursor = v_1;
        /** @const */ var /** number */ v_2 = base.cursor;
        lab1: {
            while(true)
            {
                /** @const */ var /** number */ v_3 = base.cursor;
                lab2: {
                    lab3: {
                        /** @const */ var /** number */ v_4 = base.cursor;
                        lab4: {
                            base.bra = base.cursor;
                            among_var = base.find_among(a_2);
                            if (among_var == 0)
                            {
                                break lab4;
                            }
                            base.ket = base.cursor;
                            switch (among_var) {
                                case 1:
                                    if (!base.slice_from("\u0627"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 2:
                                    if (!base.slice_from("\u0648"))
                                    {
                                        return false;
                                    }
                                    break;
                                case 3:
                                    if (!base.slice_from("\u064A"))
                                    {
                                        return false;
                                    }
                                    break;
                            }
                            break lab3;
                        }
                        base.cursor = v_4;
                        if (base.cursor >= base.limit)
                        {
                            break lab2;
                        }
                        base.cursor++;
                    }
                    continue;
                }
                base.cursor = v_3;
                break;
            }
        }
        base.cursor = v_2;
        return true;
    };

    /** @return {boolean} */
    function r_Checks1() {
        var /** number */ among_var;
        base.bra = base.cursor;
        among_var = base.find_among(a_3);
        if (among_var == 0)
        {
            return false;
        }
        base.ket = base.cursor;
        switch (among_var) {
            case 1:
                if (base.current.length <= 4)
                {
                    return false;
                }
                B_is_noun = true;
                B_is_verb = false;
                B_is_defined = true;
                break;
            case 2:
                if (base.current.length <= 3)
                {
                    return false;
                }
                B_is_noun = true;
                B_is_verb = false;
                B_is_defined = true;
                break;
        }
        return true;
    };

    /** @return {boolean} */
    function r_Prefix_Step1() {
        var /** number */ among_var;
        base.bra = base.cursor;
        among_var = base.find_among(a_4);
        if (among_var == 0)
        {
            return false;
        }
        base.ket = base.cursor;
        switch (among_var) {
            case 1:
                if (base.current.length <= 3)
                {
                    return false;
                }
                if (!base.slice_from("\u0623"))
                {
                    return false;
                }
                break;
            case 2:
                if (base.current.length <= 3)
                {
                    return false;
                }
                if (!base.slice_from("\u0622"))
                {
                    return false;
                }
                break;
            case 3:
                if (base.current.length <= 3)
                {
                    return false;
                }
                if (!base.slice_from("\u0627"))
                {
                    return false;
                }
                break;
            case 4:
                if (base.current.length <= 3)
                {
                    return false;
                }
                if (!base.slice_from("\u0625"))
                {
                    return false;
                }
                break;
        }
        return true;
    };

    /** @return {boolean} */
    function r_Prefix_Step2() {
        base.bra = base.cursor;
        if (base.find_among(a_5) == 0)
        {
            return false;
        }
        base.ket = base.cursor;
        if (base.current.length <= 3)
        {
            return false;
        }
        {
            /** @const */ var /** number */ v_1 = base.cursor;
            lab0: {
                if (!(base.eq_s("\u0627")))
                {
                    break lab0;
                }
                return false;
            }
            base.cursor = v_1;
        }
        if (!base.slice_del())
        {
            return false;
        }
        return true;
    };

    /** @return {boolean} */
    function r_Prefix_Step3a_Noun() {
        var /** number */ among_var;
        base.bra = base.cursor;
        among_var = base.find_among(a_6);
        if (among_var == 0)
        {
            return false;
        }
        base.ket = base.cursor;
        switch (among_var) {
            case 1:
                if (base.current.length <= 5)
                {
                    return false;
                }
                if (!base.slice_del())
                {
                    return false;
                }
                break;
            case 2:
                if (base.current.length <= 4)
                {
                    return false;
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
    function r_Prefix_Step3b_Noun() {
        var /** number */ among_var;
        base.bra = base.cursor;
        among_var = base.find_among(a_7);
        if (among_var == 0)
        {
            return false;
        }
        base.ket = base.cursor;
        switch (among_var) {
            case 1:
                if (base.current.length <= 3)
                {
                    return false;
                }
                if (!base.slice_del())
                {
                    return false;
                }
                break;
            case 2:
                if (base.current.length <= 3)
                {
                    return false;
                }
                if (!base.slice_from("\u0628"))
                {
                    return false;
                }
                break;
            case 3:
                if (base.current.length <= 3)
                {
                    return false;
                }
                if (!base.slice_from("\u0643"))
                {
                    return false;
                }
                break;
        }
        return true;
    };

    /** @return {boolean} */
    function r_Prefix_Step3_Verb() {
        var /** number */ among_var;
        base.bra = base.cursor;
        among_var = base.find_among(a_8);
        if (among_var == 0)
        {
            return false;
        }
        base.ket = base.cursor;
        switch (among_var) {
            case 1:
                if (base.current.length <= 4)
                {
                    return false;
                }
                if (!base.slice_from("\u064A"))
                {
                    return false;
                }
                break;
            case 2:
                if (base.current.length <= 4)
                {
                    return false;
                }
                if (!base.slice_from("\u062A"))
                {
                    return false;
                }
                break;
            case 3:
                if (base.current.length <= 4)
                {
                    return false;
                }
                if (!base.slice_from("\u0646"))
                {
                    return false;
                }
                break;
            case 4:
                if (base.current.length <= 4)
                {
                    return false;
                }
                if (!base.slice_from("\u0623"))
                {
                    return false;
                }
                break;
        }
        return true;
    };

    /** @return {boolean} */
    function r_Prefix_Step4_Verb() {
        base.bra = base.cursor;
        if (base.find_among(a_9) == 0)
        {
            return false;
        }
        base.ket = base.cursor;
        if (base.current.length <= 4)
        {
            return false;
        }
        B_is_verb = true;
        B_is_noun = false;
        if (!base.slice_from("\u0627\u0633\u062A"))
        {
            return false;
        }
        return true;
    };

    /** @return {boolean} */
    function r_Suffix_Noun_Step1a() {
        var /** number */ among_var;
        base.ket = base.cursor;
        among_var = base.find_among_b(a_10);
        if (among_var == 0)
        {
            return false;
        }
        base.bra = base.cursor;
        switch (among_var) {
            case 1:
                if (base.current.length < 4)
                {
                    return false;
                }
                if (!base.slice_del())
                {
                    return false;
                }
                break;
            case 2:
                if (base.current.length < 5)
                {
                    return false;
                }
                if (!base.slice_del())
                {
                    return false;
                }
                break;
            case 3:
                if (base.current.length < 6)
                {
                    return false;
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
    function r_Suffix_Noun_Step1b() {
        base.ket = base.cursor;
        if (base.find_among_b(a_11) == 0)
        {
            return false;
        }
        base.bra = base.cursor;
        if (base.current.length <= 5)
        {
            return false;
        }
        if (!base.slice_del())
        {
            return false;
        }
        return true;
    };

    /** @return {boolean} */
    function r_Suffix_Noun_Step2a() {
        base.ket = base.cursor;
        if (base.find_among_b(a_12) == 0)
        {
            return false;
        }
        base.bra = base.cursor;
        if (base.current.length <= 4)
        {
            return false;
        }
        if (!base.slice_del())
        {
            return false;
        }
        return true;
    };

    /** @return {boolean} */
    function r_Suffix_Noun_Step2b() {
        base.ket = base.cursor;
        if (base.find_among_b(a_13) == 0)
        {
            return false;
        }
        base.bra = base.cursor;
        if (base.current.length < 5)
        {
            return false;
        }
        if (!base.slice_del())
        {
            return false;
        }
        return true;
    };

    /** @return {boolean} */
    function r_Suffix_Noun_Step2c1() {
        base.ket = base.cursor;
        if (base.find_among_b(a_14) == 0)
        {
            return false;
        }
        base.bra = base.cursor;
        if (base.current.length < 4)
        {
            return false;
        }
        if (!base.slice_del())
        {
            return false;
        }
        return true;
    };

    /** @return {boolean} */
    function r_Suffix_Noun_Step2c2() {
        base.ket = base.cursor;
        if (base.find_among_b(a_15) == 0)
        {
            return false;
        }
        base.bra = base.cursor;
        if (base.current.length < 4)
        {
            return false;
        }
        if (!base.slice_del())
        {
            return false;
        }
        return true;
    };

    /** @return {boolean} */
    function r_Suffix_Noun_Step3() {
        base.ket = base.cursor;
        if (base.find_among_b(a_16) == 0)
        {
            return false;
        }
        base.bra = base.cursor;
        if (base.current.length < 3)
        {
            return false;
        }
        if (!base.slice_del())
        {
            return false;
        }
        return true;
    };

    /** @return {boolean} */
    function r_Suffix_Verb_Step1() {
        var /** number */ among_var;
        base.ket = base.cursor;
        among_var = base.find_among_b(a_17);
        if (among_var == 0)
        {
            return false;
        }
        base.bra = base.cursor;
        switch (among_var) {
            case 1:
                if (base.current.length < 4)
                {
                    return false;
                }
                if (!base.slice_del())
                {
                    return false;
                }
                break;
            case 2:
                if (base.current.length < 5)
                {
                    return false;
                }
                if (!base.slice_del())
                {
                    return false;
                }
                break;
            case 3:
                if (base.current.length < 6)
                {
                    return false;
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
    function r_Suffix_Verb_Step2a() {
        var /** number */ among_var;
        base.ket = base.cursor;
        among_var = base.find_among_b(a_18);
        if (among_var == 0)
        {
            return false;
        }
        base.bra = base.cursor;
        switch (among_var) {
            case 1:
                if (base.current.length < 4)
                {
                    return false;
                }
                if (!base.slice_del())
                {
                    return false;
                }
                break;
            case 2:
                if (base.current.length < 5)
                {
                    return false;
                }
                if (!base.slice_del())
                {
                    return false;
                }
                break;
            case 3:
                if (base.current.length <= 5)
                {
                    return false;
                }
                if (!base.slice_del())
                {
                    return false;
                }
                break;
            case 4:
                if (base.current.length < 6)
                {
                    return false;
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
    function r_Suffix_Verb_Step2b() {
        base.ket = base.cursor;
        if (base.find_among_b(a_19) == 0)
        {
            return false;
        }
        base.bra = base.cursor;
        if (base.current.length < 5)
        {
            return false;
        }
        if (!base.slice_del())
        {
            return false;
        }
        return true;
    };

    /** @return {boolean} */
    function r_Suffix_Verb_Step2c() {
        var /** number */ among_var;
        base.ket = base.cursor;
        among_var = base.find_among_b(a_20);
        if (among_var == 0)
        {
            return false;
        }
        base.bra = base.cursor;
        switch (among_var) {
            case 1:
                if (base.current.length < 4)
                {
                    return false;
                }
                if (!base.slice_del())
                {
                    return false;
                }
                break;
            case 2:
                if (base.current.length < 6)
                {
                    return false;
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
    function r_Suffix_All_alef_maqsura() {
        base.ket = base.cursor;
        if (base.find_among_b(a_21) == 0)
        {
            return false;
        }
        base.bra = base.cursor;
        if (!base.slice_from("\u064A"))
        {
            return false;
        }
        return true;
    };

    this.stem = /** @return {boolean} */ function() {
        B_is_noun = true;
        B_is_verb = true;
        B_is_defined = false;
        /** @const */ var /** number */ v_1 = base.cursor;
        r_Checks1();
        base.cursor = v_1;
        r_Normalize_pre();
        base.limit_backward = base.cursor; base.cursor = base.limit;
        /** @const */ var /** number */ v_2 = base.limit - base.cursor;
        lab0: {
            lab1: {
                /** @const */ var /** number */ v_3 = base.limit - base.cursor;
                lab2: {
                    if (!B_is_verb)
                    {
                        break lab2;
                    }
                    lab3: {
                        /** @const */ var /** number */ v_4 = base.limit - base.cursor;
                        lab4: {
                            {
                                var v_5 = 1;
                                while(true)
                                {
                                    /** @const */ var /** number */ v_6 = base.limit - base.cursor;
                                    lab5: {
                                        if (!r_Suffix_Verb_Step1())
                                        {
                                            break lab5;
                                        }
                                        v_5--;
                                        continue;
                                    }
                                    base.cursor = base.limit - v_6;
                                    break;
                                }
                                if (v_5 > 0)
                                {
                                    break lab4;
                                }
                            }
                            lab6: {
                                /** @const */ var /** number */ v_7 = base.limit - base.cursor;
                                lab7: {
                                    if (!r_Suffix_Verb_Step2a())
                                    {
                                        break lab7;
                                    }
                                    break lab6;
                                }
                                base.cursor = base.limit - v_7;
                                lab8: {
                                    if (!r_Suffix_Verb_Step2c())
                                    {
                                        break lab8;
                                    }
                                    break lab6;
                                }
                                base.cursor = base.limit - v_7;
                                if (base.cursor <= base.limit_backward)
                                {
                                    break lab4;
                                }
                                base.cursor--;
                            }
                            break lab3;
                        }
                        base.cursor = base.limit - v_4;
                        lab9: {
                            if (!r_Suffix_Verb_Step2b())
                            {
                                break lab9;
                            }
                            break lab3;
                        }
                        base.cursor = base.limit - v_4;
                        if (!r_Suffix_Verb_Step2a())
                        {
                            break lab2;
                        }
                    }
                    break lab1;
                }
                base.cursor = base.limit - v_3;
                lab10: {
                    if (!B_is_noun)
                    {
                        break lab10;
                    }
                    /** @const */ var /** number */ v_8 = base.limit - base.cursor;
                    lab11: {
                        lab12: {
                            /** @const */ var /** number */ v_9 = base.limit - base.cursor;
                            lab13: {
                                if (!r_Suffix_Noun_Step2c2())
                                {
                                    break lab13;
                                }
                                break lab12;
                            }
                            base.cursor = base.limit - v_9;
                            lab14: {
                                lab15: {
                                    if (!B_is_defined)
                                    {
                                        break lab15;
                                    }
                                    break lab14;
                                }
                                if (!r_Suffix_Noun_Step1a())
                                {
                                    break lab14;
                                }
                                lab16: {
                                    /** @const */ var /** number */ v_10 = base.limit - base.cursor;
                                    lab17: {
                                        if (!r_Suffix_Noun_Step2a())
                                        {
                                            break lab17;
                                        }
                                        break lab16;
                                    }
                                    base.cursor = base.limit - v_10;
                                    lab18: {
                                        if (!r_Suffix_Noun_Step2b())
                                        {
                                            break lab18;
                                        }
                                        break lab16;
                                    }
                                    base.cursor = base.limit - v_10;
                                    lab19: {
                                        if (!r_Suffix_Noun_Step2c1())
                                        {
                                            break lab19;
                                        }
                                        break lab16;
                                    }
                                    base.cursor = base.limit - v_10;
                                    if (base.cursor <= base.limit_backward)
                                    {
                                        break lab14;
                                    }
                                    base.cursor--;
                                }
                                break lab12;
                            }
                            base.cursor = base.limit - v_9;
                            lab20: {
                                if (!r_Suffix_Noun_Step1b())
                                {
                                    break lab20;
                                }
                                lab21: {
                                    /** @const */ var /** number */ v_11 = base.limit - base.cursor;
                                    lab22: {
                                        if (!r_Suffix_Noun_Step2a())
                                        {
                                            break lab22;
                                        }
                                        break lab21;
                                    }
                                    base.cursor = base.limit - v_11;
                                    lab23: {
                                        if (!r_Suffix_Noun_Step2b())
                                        {
                                            break lab23;
                                        }
                                        break lab21;
                                    }
                                    base.cursor = base.limit - v_11;
                                    if (!r_Suffix_Noun_Step2c1())
                                    {
                                        break lab20;
                                    }
                                }
                                break lab12;
                            }
                            base.cursor = base.limit - v_9;
                            lab24: {
                                lab25: {
                                    if (!B_is_defined)
                                    {
                                        break lab25;
                                    }
                                    break lab24;
                                }
                                if (!r_Suffix_Noun_Step2a())
                                {
                                    break lab24;
                                }
                                break lab12;
                            }
                            base.cursor = base.limit - v_9;
                            if (!r_Suffix_Noun_Step2b())
                            {
                                base.cursor = base.limit - v_8;
                                break lab11;
                            }
                        }
                    }
                    if (!r_Suffix_Noun_Step3())
                    {
                        break lab10;
                    }
                    break lab1;
                }
                base.cursor = base.limit - v_3;
                if (!r_Suffix_All_alef_maqsura())
                {
                    break lab0;
                }
            }
        }
        base.cursor = base.limit - v_2;
        base.cursor = base.limit_backward;
        /** @const */ var /** number */ v_12 = base.cursor;
        lab26: {
            /** @const */ var /** number */ v_13 = base.cursor;
            lab27: {
                if (!r_Prefix_Step1())
                {
                    base.cursor = v_13;
                    break lab27;
                }
            }
            /** @const */ var /** number */ v_14 = base.cursor;
            lab28: {
                if (!r_Prefix_Step2())
                {
                    base.cursor = v_14;
                    break lab28;
                }
            }
            lab29: {
                /** @const */ var /** number */ v_15 = base.cursor;
                lab30: {
                    if (!r_Prefix_Step3a_Noun())
                    {
                        break lab30;
                    }
                    break lab29;
                }
                base.cursor = v_15;
                lab31: {
                    if (!B_is_noun)
                    {
                        break lab31;
                    }
                    if (!r_Prefix_Step3b_Noun())
                    {
                        break lab31;
                    }
                    break lab29;
                }
                base.cursor = v_15;
                if (!B_is_verb)
                {
                    break lab26;
                }
                /** @const */ var /** number */ v_16 = base.cursor;
                lab32: {
                    if (!r_Prefix_Step3_Verb())
                    {
                        base.cursor = v_16;
                        break lab32;
                    }
                }
                if (!r_Prefix_Step4_Verb())
                {
                    break lab26;
                }
            }
        }
        base.cursor = v_12;
        r_Normalize_post();
        return true;
    };

    /**@return{string}*/
    this['stemWord'] = function(/**string*/word) {
        base.setCurrent(word);
        this.stem();
        return base.getCurrent();
    };
};
