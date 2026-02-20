// Generated from french.sbl by Snowball 3.0.1 - https://snowballstem.org/

/**@constructor*/
var FrenchStemmer = function() {
    var base = new BaseStemmer();

    /** @const */ var a_0 = [
        ["col", -1, -1],
        ["ni", -1, 1],
        ["par", -1, -1],
        ["tap", -1, -1]
    ];

    /** @const */ var a_1 = [
        ["", -1, 7],
        ["H", 0, 6],
        ["He", 1, 4],
        ["Hi", 1, 5],
        ["I", 0, 1],
        ["U", 0, 2],
        ["Y", 0, 3]
    ];

    /** @const */ var a_2 = [
        ["iqU", -1, 3],
        ["abl", -1, 3],
        ["I\u00E8r", -1, 4],
        ["i\u00E8r", -1, 4],
        ["eus", -1, 2],
        ["iv", -1, 1]
    ];

    /** @const */ var a_3 = [
        ["ic", -1, 2],
        ["abil", -1, 1],
        ["iv", -1, 3]
    ];

    /** @const */ var a_4 = [
        ["iqUe", -1, 1],
        ["atrice", -1, 2],
        ["ance", -1, 1],
        ["ence", -1, 5],
        ["logie", -1, 3],
        ["able", -1, 1],
        ["isme", -1, 1],
        ["euse", -1, 12],
        ["iste", -1, 1],
        ["ive", -1, 8],
        ["if", -1, 8],
        ["usion", -1, 4],
        ["ation", -1, 2],
        ["ution", -1, 4],
        ["ateur", -1, 2],
        ["iqUes", -1, 1],
        ["atrices", -1, 2],
        ["ances", -1, 1],
        ["ences", -1, 5],
        ["logies", -1, 3],
        ["ables", -1, 1],
        ["ismes", -1, 1],
        ["euses", -1, 12],
        ["istes", -1, 1],
        ["ives", -1, 8],
        ["ifs", -1, 8],
        ["usions", -1, 4],
        ["ations", -1, 2],
        ["utions", -1, 4],
        ["ateurs", -1, 2],
        ["ments", -1, 16],
        ["ements", 30, 6],
        ["issements", 31, 13],
        ["it\u00E9s", -1, 7],
        ["ment", -1, 16],
        ["ement", 34, 6],
        ["issement", 35, 13],
        ["amment", 34, 14],
        ["emment", 34, 15],
        ["aux", -1, 10],
        ["eaux", 39, 9],
        ["eux", -1, 1],
        ["oux", -1, 11],
        ["it\u00E9", -1, 7]
    ];

    /** @const */ var a_5 = [
        ["ira", -1, 1],
        ["ie", -1, 1],
        ["isse", -1, 1],
        ["issante", -1, 1],
        ["i", -1, 1],
        ["irai", 4, 1],
        ["ir", -1, 1],
        ["iras", -1, 1],
        ["ies", -1, 1],
        ["\u00EEmes", -1, 1],
        ["isses", -1, 1],
        ["issantes", -1, 1],
        ["\u00EEtes", -1, 1],
        ["is", -1, 1],
        ["irais", 13, 1],
        ["issais", 13, 1],
        ["irions", -1, 1],
        ["issions", -1, 1],
        ["irons", -1, 1],
        ["issons", -1, 1],
        ["issants", -1, 1],
        ["it", -1, 1],
        ["irait", 21, 1],
        ["issait", 21, 1],
        ["issant", -1, 1],
        ["iraIent", -1, 1],
        ["issaIent", -1, 1],
        ["irent", -1, 1],
        ["issent", -1, 1],
        ["iront", -1, 1],
        ["\u00EEt", -1, 1],
        ["iriez", -1, 1],
        ["issiez", -1, 1],
        ["irez", -1, 1],
        ["issez", -1, 1]
    ];

    /** @const */ var a_6 = [
        ["al", -1, 1],
        ["\u00E9pl", -1, -1],
        ["auv", -1, -1]
    ];

    /** @const */ var a_7 = [
        ["a", -1, 3],
        ["era", 0, 2],
        ["aise", -1, 4],
        ["asse", -1, 3],
        ["ante", -1, 3],
        ["\u00E9e", -1, 2],
        ["ai", -1, 3],
        ["erai", 6, 2],
        ["er", -1, 2],
        ["as", -1, 3],
        ["eras", 9, 2],
        ["\u00E2mes", -1, 3],
        ["aises", -1, 4],
        ["asses", -1, 3],
        ["antes", -1, 3],
        ["\u00E2tes", -1, 3],
        ["\u00E9es", -1, 2],
        ["ais", -1, 4],
        ["eais", 17, 2],
        ["erais", 17, 2],
        ["ions", -1, 1],
        ["erions", 20, 2],
        ["assions", 20, 3],
        ["erons", -1, 2],
        ["ants", -1, 3],
        ["\u00E9s", -1, 2],
        ["ait", -1, 3],
        ["erait", 26, 2],
        ["ant", -1, 3],
        ["aIent", -1, 3],
        ["eraIent", 29, 2],
        ["\u00E8rent", -1, 2],
        ["assent", -1, 3],
        ["eront", -1, 2],
        ["\u00E2t", -1, 3],
        ["ez", -1, 2],
        ["iez", 35, 2],
        ["eriez", 36, 2],
        ["assiez", 36, 3],
        ["erez", 35, 2],
        ["\u00E9", -1, 2]
    ];

    /** @const */ var a_8 = [
        ["e", -1, 3],
        ["I\u00E8re", 0, 2],
        ["i\u00E8re", 0, 2],
        ["ion", -1, 1],
        ["Ier", -1, 2],
        ["ier", -1, 2]
    ];

    /** @const */ var a_9 = [
        ["ell", -1, -1],
        ["eill", -1, -1],
        ["enn", -1, -1],
        ["onn", -1, -1],
        ["ett", -1, -1]
    ];

    /** @const */ var /** Array<int> */ g_v = [17, 65, 16, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 128, 130, 103, 8, 5];

    /** @const */ var /** Array<int> */ g_oux_ending = [65, 85];

    /** @const */ var /** Array<int> */ g_elision_char = [131, 14, 3];

    /** @const */ var /** Array<int> */ g_keep_with_s = [1, 65, 20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 128];

    var /** number */ I_p2 = 0;
    var /** number */ I_p1 = 0;
    var /** number */ I_pV = 0;


    /** @return {boolean} */
    function r_elisions() {
        base.bra = base.cursor;
        lab0: {
            /** @const */ var /** number */ v_1 = base.cursor;
            lab1: {
                if (!(base.in_grouping(g_elision_char, 99, 116)))
                {
                    break lab1;
                }
                break lab0;
            }
            base.cursor = v_1;
            if (!(base.eq_s("qu")))
            {
                return false;
            }
        }
        if (!(base.eq_s("'")))
        {
            return false;
        }
        base.ket = base.cursor;
        lab2: {
            if (base.cursor < base.limit)
            {
                break lab2;
            }
            return false;
        }
        if (!base.slice_del())
        {
            return false;
        }
        return true;
    };

    /** @return {boolean} */
    function r_prelude() {
        while(true)
        {
            /** @const */ var /** number */ v_1 = base.cursor;
            lab0: {
                golab1: while(true)
                {
                    /** @const */ var /** number */ v_2 = base.cursor;
                    lab2: {
                        lab3: {
                            /** @const */ var /** number */ v_3 = base.cursor;
                            lab4: {
                                if (!(base.in_grouping(g_v, 97, 251)))
                                {
                                    break lab4;
                                }
                                base.bra = base.cursor;
                                lab5: {
                                    /** @const */ var /** number */ v_4 = base.cursor;
                                    lab6: {
                                        if (!(base.eq_s("u")))
                                        {
                                            break lab6;
                                        }
                                        base.ket = base.cursor;
                                        if (!(base.in_grouping(g_v, 97, 251)))
                                        {
                                            break lab6;
                                        }
                                        if (!base.slice_from("U"))
                                        {
                                            return false;
                                        }
                                        break lab5;
                                    }
                                    base.cursor = v_4;
                                    lab7: {
                                        if (!(base.eq_s("i")))
                                        {
                                            break lab7;
                                        }
                                        base.ket = base.cursor;
                                        if (!(base.in_grouping(g_v, 97, 251)))
                                        {
                                            break lab7;
                                        }
                                        if (!base.slice_from("I"))
                                        {
                                            return false;
                                        }
                                        break lab5;
                                    }
                                    base.cursor = v_4;
                                    if (!(base.eq_s("y")))
                                    {
                                        break lab4;
                                    }
                                    base.ket = base.cursor;
                                    if (!base.slice_from("Y"))
                                    {
                                        return false;
                                    }
                                }
                                break lab3;
                            }
                            base.cursor = v_3;
                            lab8: {
                                base.bra = base.cursor;
                                if (!(base.eq_s("\u00EB")))
                                {
                                    break lab8;
                                }
                                base.ket = base.cursor;
                                if (!base.slice_from("He"))
                                {
                                    return false;
                                }
                                break lab3;
                            }
                            base.cursor = v_3;
                            lab9: {
                                base.bra = base.cursor;
                                if (!(base.eq_s("\u00EF")))
                                {
                                    break lab9;
                                }
                                base.ket = base.cursor;
                                if (!base.slice_from("Hi"))
                                {
                                    return false;
                                }
                                break lab3;
                            }
                            base.cursor = v_3;
                            lab10: {
                                base.bra = base.cursor;
                                if (!(base.eq_s("y")))
                                {
                                    break lab10;
                                }
                                base.ket = base.cursor;
                                if (!(base.in_grouping(g_v, 97, 251)))
                                {
                                    break lab10;
                                }
                                if (!base.slice_from("Y"))
                                {
                                    return false;
                                }
                                break lab3;
                            }
                            base.cursor = v_3;
                            if (!(base.eq_s("q")))
                            {
                                break lab2;
                            }
                            base.bra = base.cursor;
                            if (!(base.eq_s("u")))
                            {
                                break lab2;
                            }
                            base.ket = base.cursor;
                            if (!base.slice_from("U"))
                            {
                                return false;
                            }
                        }
                        base.cursor = v_2;
                        break golab1;
                    }
                    base.cursor = v_2;
                    if (base.cursor >= base.limit)
                    {
                        break lab0;
                    }
                    base.cursor++;
                }
                continue;
            }
            base.cursor = v_1;
            break;
        }
        return true;
    };

    /** @return {boolean} */
    function r_mark_regions() {
        var /** number */ among_var;
        I_pV = base.limit;
        I_p1 = base.limit;
        I_p2 = base.limit;
        /** @const */ var /** number */ v_1 = base.cursor;
        lab0: {
            lab1: {
                /** @const */ var /** number */ v_2 = base.cursor;
                lab2: {
                    if (!(base.in_grouping(g_v, 97, 251)))
                    {
                        break lab2;
                    }
                    if (!(base.in_grouping(g_v, 97, 251)))
                    {
                        break lab2;
                    }
                    if (base.cursor >= base.limit)
                    {
                        break lab2;
                    }
                    base.cursor++;
                    break lab1;
                }
                base.cursor = v_2;
                lab3: {
                    among_var = base.find_among(a_0);
                    if (among_var == 0)
                    {
                        break lab3;
                    }
                    switch (among_var) {
                        case 1:
                            if (!(base.in_grouping(g_v, 97, 251)))
                            {
                                break lab3;
                            }
                            break;
                    }
                    break lab1;
                }
                base.cursor = v_2;
                if (base.cursor >= base.limit)
                {
                    break lab0;
                }
                base.cursor++;
                if (!base.go_out_grouping(g_v, 97, 251))
                {
                    break lab0;
                }
                base.cursor++;
            }
            I_pV = base.cursor;
        }
        base.cursor = v_1;
        /** @const */ var /** number */ v_3 = base.cursor;
        lab4: {
            if (!base.go_out_grouping(g_v, 97, 251))
            {
                break lab4;
            }
            base.cursor++;
            if (!base.go_in_grouping(g_v, 97, 251))
            {
                break lab4;
            }
            base.cursor++;
            I_p1 = base.cursor;
            if (!base.go_out_grouping(g_v, 97, 251))
            {
                break lab4;
            }
            base.cursor++;
            if (!base.go_in_grouping(g_v, 97, 251))
            {
                break lab4;
            }
            base.cursor++;
            I_p2 = base.cursor;
        }
        base.cursor = v_3;
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
                        if (!base.slice_from("i"))
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
                        if (!base.slice_from("y"))
                        {
                            return false;
                        }
                        break;
                    case 4:
                        if (!base.slice_from("\u00EB"))
                        {
                            return false;
                        }
                        break;
                    case 5:
                        if (!base.slice_from("\u00EF"))
                        {
                            return false;
                        }
                        break;
                    case 6:
                        if (!base.slice_del())
                        {
                            return false;
                        }
                        break;
                    case 7:
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
    function r_RV() {
        return I_pV <= base.cursor;
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
        base.ket = base.cursor;
        among_var = base.find_among_b(a_4);
        if (among_var == 0)
        {
            return false;
        }
        base.bra = base.cursor;
        switch (among_var) {
            case 1:
                if (!r_R2())
                {
                    return false;
                }
                if (!base.slice_del())
                {
                    return false;
                }
                break;
            case 2:
                if (!r_R2())
                {
                    return false;
                }
                if (!base.slice_del())
                {
                    return false;
                }
                /** @const */ var /** number */ v_1 = base.limit - base.cursor;
                lab0: {
                    base.ket = base.cursor;
                    if (!(base.eq_s_b("ic")))
                    {
                        base.cursor = base.limit - v_1;
                        break lab0;
                    }
                    base.bra = base.cursor;
                    lab1: {
                        /** @const */ var /** number */ v_2 = base.limit - base.cursor;
                        lab2: {
                            if (!r_R2())
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
                        if (!base.slice_from("iqU"))
                        {
                            return false;
                        }
                    }
                }
                break;
            case 3:
                if (!r_R2())
                {
                    return false;
                }
                if (!base.slice_from("log"))
                {
                    return false;
                }
                break;
            case 4:
                if (!r_R2())
                {
                    return false;
                }
                if (!base.slice_from("u"))
                {
                    return false;
                }
                break;
            case 5:
                if (!r_R2())
                {
                    return false;
                }
                if (!base.slice_from("ent"))
                {
                    return false;
                }
                break;
            case 6:
                if (!r_RV())
                {
                    return false;
                }
                if (!base.slice_del())
                {
                    return false;
                }
                /** @const */ var /** number */ v_3 = base.limit - base.cursor;
                lab3: {
                    base.ket = base.cursor;
                    among_var = base.find_among_b(a_2);
                    if (among_var == 0)
                    {
                        base.cursor = base.limit - v_3;
                        break lab3;
                    }
                    base.bra = base.cursor;
                    switch (among_var) {
                        case 1:
                            if (!r_R2())
                            {
                                base.cursor = base.limit - v_3;
                                break lab3;
                            }
                            if (!base.slice_del())
                            {
                                return false;
                            }
                            base.ket = base.cursor;
                            if (!(base.eq_s_b("at")))
                            {
                                base.cursor = base.limit - v_3;
                                break lab3;
                            }
                            base.bra = base.cursor;
                            if (!r_R2())
                            {
                                base.cursor = base.limit - v_3;
                                break lab3;
                            }
                            if (!base.slice_del())
                            {
                                return false;
                            }
                            break;
                        case 2:
                            lab4: {
                                /** @const */ var /** number */ v_4 = base.limit - base.cursor;
                                lab5: {
                                    if (!r_R2())
                                    {
                                        break lab5;
                                    }
                                    if (!base.slice_del())
                                    {
                                        return false;
                                    }
                                    break lab4;
                                }
                                base.cursor = base.limit - v_4;
                                if (!r_R1())
                                {
                                    base.cursor = base.limit - v_3;
                                    break lab3;
                                }
                                if (!base.slice_from("eux"))
                                {
                                    return false;
                                }
                            }
                            break;
                        case 3:
                            if (!r_R2())
                            {
                                base.cursor = base.limit - v_3;
                                break lab3;
                            }
                            if (!base.slice_del())
                            {
                                return false;
                            }
                            break;
                        case 4:
                            if (!r_RV())
                            {
                                base.cursor = base.limit - v_3;
                                break lab3;
                            }
                            if (!base.slice_from("i"))
                            {
                                return false;
                            }
                            break;
                    }
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
                /** @const */ var /** number */ v_5 = base.limit - base.cursor;
                lab6: {
                    base.ket = base.cursor;
                    among_var = base.find_among_b(a_3);
                    if (among_var == 0)
                    {
                        base.cursor = base.limit - v_5;
                        break lab6;
                    }
                    base.bra = base.cursor;
                    switch (among_var) {
                        case 1:
                            lab7: {
                                /** @const */ var /** number */ v_6 = base.limit - base.cursor;
                                lab8: {
                                    if (!r_R2())
                                    {
                                        break lab8;
                                    }
                                    if (!base.slice_del())
                                    {
                                        return false;
                                    }
                                    break lab7;
                                }
                                base.cursor = base.limit - v_6;
                                if (!base.slice_from("abl"))
                                {
                                    return false;
                                }
                            }
                            break;
                        case 2:
                            lab9: {
                                /** @const */ var /** number */ v_7 = base.limit - base.cursor;
                                lab10: {
                                    if (!r_R2())
                                    {
                                        break lab10;
                                    }
                                    if (!base.slice_del())
                                    {
                                        return false;
                                    }
                                    break lab9;
                                }
                                base.cursor = base.limit - v_7;
                                if (!base.slice_from("iqU"))
                                {
                                    return false;
                                }
                            }
                            break;
                        case 3:
                            if (!r_R2())
                            {
                                base.cursor = base.limit - v_5;
                                break lab6;
                            }
                            if (!base.slice_del())
                            {
                                return false;
                            }
                            break;
                    }
                }
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
                /** @const */ var /** number */ v_8 = base.limit - base.cursor;
                lab11: {
                    base.ket = base.cursor;
                    if (!(base.eq_s_b("at")))
                    {
                        base.cursor = base.limit - v_8;
                        break lab11;
                    }
                    base.bra = base.cursor;
                    if (!r_R2())
                    {
                        base.cursor = base.limit - v_8;
                        break lab11;
                    }
                    if (!base.slice_del())
                    {
                        return false;
                    }
                    base.ket = base.cursor;
                    if (!(base.eq_s_b("ic")))
                    {
                        base.cursor = base.limit - v_8;
                        break lab11;
                    }
                    base.bra = base.cursor;
                    lab12: {
                        /** @const */ var /** number */ v_9 = base.limit - base.cursor;
                        lab13: {
                            if (!r_R2())
                            {
                                break lab13;
                            }
                            if (!base.slice_del())
                            {
                                return false;
                            }
                            break lab12;
                        }
                        base.cursor = base.limit - v_9;
                        if (!base.slice_from("iqU"))
                        {
                            return false;
                        }
                    }
                }
                break;
            case 9:
                if (!base.slice_from("eau"))
                {
                    return false;
                }
                break;
            case 10:
                if (!r_R1())
                {
                    return false;
                }
                if (!base.slice_from("al"))
                {
                    return false;
                }
                break;
            case 11:
                if (!(base.in_grouping_b(g_oux_ending, 98, 112)))
                {
                    return false;
                }
                if (!base.slice_from("ou"))
                {
                    return false;
                }
                break;
            case 12:
                lab14: {
                    /** @const */ var /** number */ v_10 = base.limit - base.cursor;
                    lab15: {
                        if (!r_R2())
                        {
                            break lab15;
                        }
                        if (!base.slice_del())
                        {
                            return false;
                        }
                        break lab14;
                    }
                    base.cursor = base.limit - v_10;
                    if (!r_R1())
                    {
                        return false;
                    }
                    if (!base.slice_from("eux"))
                    {
                        return false;
                    }
                }
                break;
            case 13:
                if (!r_R1())
                {
                    return false;
                }
                if (!(base.out_grouping_b(g_v, 97, 251)))
                {
                    return false;
                }
                if (!base.slice_del())
                {
                    return false;
                }
                break;
            case 14:
                if (!r_RV())
                {
                    return false;
                }
                if (!base.slice_from("ant"))
                {
                    return false;
                }
                return false;
            case 15:
                if (!r_RV())
                {
                    return false;
                }
                if (!base.slice_from("ent"))
                {
                    return false;
                }
                return false;
            case 16:
                /** @const */ var /** number */ v_11 = base.limit - base.cursor;
                if (!(base.in_grouping_b(g_v, 97, 251)))
                {
                    return false;
                }
                if (!r_RV())
                {
                    return false;
                }
                base.cursor = base.limit - v_11;
                if (!base.slice_del())
                {
                    return false;
                }
                return false;
        }
        return true;
    };

    /** @return {boolean} */
    function r_i_verb_suffix() {
        if (base.cursor < I_pV)
        {
            return false;
        }
        /** @const */ var /** number */ v_1 = base.limit_backward;
        base.limit_backward = I_pV;
        base.ket = base.cursor;
        if (base.find_among_b(a_5) == 0)
        {
            base.limit_backward = v_1;
            return false;
        }
        base.bra = base.cursor;
        {
            /** @const */ var /** number */ v_2 = base.limit - base.cursor;
            lab0: {
                if (!(base.eq_s_b("H")))
                {
                    break lab0;
                }
                base.limit_backward = v_1;
                return false;
            }
            base.cursor = base.limit - v_2;
        }
        if (!(base.out_grouping_b(g_v, 97, 251)))
        {
            base.limit_backward = v_1;
            return false;
        }
        if (!base.slice_del())
        {
            return false;
        }
        base.limit_backward = v_1;
        return true;
    };

    /** @return {boolean} */
    function r_verb_suffix() {
        var /** number */ among_var;
        if (base.cursor < I_pV)
        {
            return false;
        }
        /** @const */ var /** number */ v_1 = base.limit_backward;
        base.limit_backward = I_pV;
        base.ket = base.cursor;
        among_var = base.find_among_b(a_7);
        if (among_var == 0)
        {
            base.limit_backward = v_1;
            return false;
        }
        base.bra = base.cursor;
        base.limit_backward = v_1;
        switch (among_var) {
            case 1:
                if (!r_R2())
                {
                    return false;
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
                /** @const */ var /** number */ v_2 = base.limit - base.cursor;
                lab0: {
                    if (!(base.eq_s_b("e")))
                    {
                        base.cursor = base.limit - v_2;
                        break lab0;
                    }
                    if (!r_RV())
                    {
                        base.cursor = base.limit - v_2;
                        break lab0;
                    }
                    base.bra = base.cursor;
                }
                if (!base.slice_del())
                {
                    return false;
                }
                break;
            case 4:
                {
                    /** @const */ var /** number */ v_3 = base.limit - base.cursor;
                    lab1: {
                        among_var = base.find_among_b(a_6);
                        if (among_var == 0)
                        {
                            break lab1;
                        }
                        switch (among_var) {
                            case 1:
                                if (base.cursor <= base.limit_backward)
                                {
                                    break lab1;
                                }
                                base.cursor--;
                                if (base.cursor > base.limit_backward)
                                {
                                    break lab1;
                                }
                                break;
                        }
                        return false;
                    }
                    base.cursor = base.limit - v_3;
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
    function r_residual_suffix() {
        var /** number */ among_var;
        /** @const */ var /** number */ v_1 = base.limit - base.cursor;
        lab0: {
            base.ket = base.cursor;
            if (!(base.eq_s_b("s")))
            {
                base.cursor = base.limit - v_1;
                break lab0;
            }
            base.bra = base.cursor;
            /** @const */ var /** number */ v_2 = base.limit - base.cursor;
            lab1: {
                /** @const */ var /** number */ v_3 = base.limit - base.cursor;
                lab2: {
                    if (!(base.eq_s_b("Hi")))
                    {
                        break lab2;
                    }
                    break lab1;
                }
                base.cursor = base.limit - v_3;
                if (!(base.out_grouping_b(g_keep_with_s, 97, 232)))
                {
                    base.cursor = base.limit - v_1;
                    break lab0;
                }
            }
            base.cursor = base.limit - v_2;
            if (!base.slice_del())
            {
                return false;
            }
        }
        if (base.cursor < I_pV)
        {
            return false;
        }
        /** @const */ var /** number */ v_4 = base.limit_backward;
        base.limit_backward = I_pV;
        base.ket = base.cursor;
        among_var = base.find_among_b(a_8);
        if (among_var == 0)
        {
            base.limit_backward = v_4;
            return false;
        }
        base.bra = base.cursor;
        switch (among_var) {
            case 1:
                if (!r_R2())
                {
                    base.limit_backward = v_4;
                    return false;
                }
                lab3: {
                    /** @const */ var /** number */ v_5 = base.limit - base.cursor;
                    lab4: {
                        if (!(base.eq_s_b("s")))
                        {
                            break lab4;
                        }
                        break lab3;
                    }
                    base.cursor = base.limit - v_5;
                    if (!(base.eq_s_b("t")))
                    {
                        base.limit_backward = v_4;
                        return false;
                    }
                }
                if (!base.slice_del())
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
            case 3:
                if (!base.slice_del())
                {
                    return false;
                }
                break;
        }
        base.limit_backward = v_4;
        return true;
    };

    /** @return {boolean} */
    function r_un_double() {
        /** @const */ var /** number */ v_1 = base.limit - base.cursor;
        if (base.find_among_b(a_9) == 0)
        {
            return false;
        }
        base.cursor = base.limit - v_1;
        base.ket = base.cursor;
        if (base.cursor <= base.limit_backward)
        {
            return false;
        }
        base.cursor--;
        base.bra = base.cursor;
        if (!base.slice_del())
        {
            return false;
        }
        return true;
    };

    /** @return {boolean} */
    function r_un_accent() {
        {
            var v_1 = 1;
            while(true)
            {
                lab0: {
                    if (!(base.out_grouping_b(g_v, 97, 251)))
                    {
                        break lab0;
                    }
                    v_1--;
                    continue;
                }
                break;
            }
            if (v_1 > 0)
            {
                return false;
            }
        }
        base.ket = base.cursor;
        lab1: {
            /** @const */ var /** number */ v_2 = base.limit - base.cursor;
            lab2: {
                if (!(base.eq_s_b("\u00E9")))
                {
                    break lab2;
                }
                break lab1;
            }
            base.cursor = base.limit - v_2;
            if (!(base.eq_s_b("\u00E8")))
            {
                return false;
            }
        }
        base.bra = base.cursor;
        if (!base.slice_from("e"))
        {
            return false;
        }
        return true;
    };

    this.stem = /** @return {boolean} */ function() {
        /** @const */ var /** number */ v_1 = base.cursor;
        r_elisions();
        base.cursor = v_1;
        /** @const */ var /** number */ v_2 = base.cursor;
        r_prelude();
        base.cursor = v_2;
        r_mark_regions();
        base.limit_backward = base.cursor; base.cursor = base.limit;
        /** @const */ var /** number */ v_3 = base.limit - base.cursor;
        lab0: {
            lab1: {
                /** @const */ var /** number */ v_4 = base.limit - base.cursor;
                lab2: {
                    /** @const */ var /** number */ v_5 = base.limit - base.cursor;
                    lab3: {
                        /** @const */ var /** number */ v_6 = base.limit - base.cursor;
                        lab4: {
                            if (!r_standard_suffix())
                            {
                                break lab4;
                            }
                            break lab3;
                        }
                        base.cursor = base.limit - v_6;
                        lab5: {
                            if (!r_i_verb_suffix())
                            {
                                break lab5;
                            }
                            break lab3;
                        }
                        base.cursor = base.limit - v_6;
                        if (!r_verb_suffix())
                        {
                            break lab2;
                        }
                    }
                    base.cursor = base.limit - v_5;
                    /** @const */ var /** number */ v_7 = base.limit - base.cursor;
                    lab6: {
                        base.ket = base.cursor;
                        lab7: {
                            /** @const */ var /** number */ v_8 = base.limit - base.cursor;
                            lab8: {
                                if (!(base.eq_s_b("Y")))
                                {
                                    break lab8;
                                }
                                base.bra = base.cursor;
                                if (!base.slice_from("i"))
                                {
                                    return false;
                                }
                                break lab7;
                            }
                            base.cursor = base.limit - v_8;
                            if (!(base.eq_s_b("\u00E7")))
                            {
                                base.cursor = base.limit - v_7;
                                break lab6;
                            }
                            base.bra = base.cursor;
                            if (!base.slice_from("c"))
                            {
                                return false;
                            }
                        }
                    }
                    break lab1;
                }
                base.cursor = base.limit - v_4;
                if (!r_residual_suffix())
                {
                    break lab0;
                }
            }
        }
        base.cursor = base.limit - v_3;
        /** @const */ var /** number */ v_9 = base.limit - base.cursor;
        r_un_double();
        base.cursor = base.limit - v_9;
        /** @const */ var /** number */ v_10 = base.limit - base.cursor;
        r_un_accent();
        base.cursor = base.limit - v_10;
        base.cursor = base.limit_backward;
        /** @const */ var /** number */ v_11 = base.cursor;
        r_postlude();
        base.cursor = v_11;
        return true;
    };

    /**@return{string}*/
    this['stemWord'] = function(/**string*/word) {
        base.setCurrent(word);
        this.stem();
        return base.getCurrent();
    };
};
