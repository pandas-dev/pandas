// Generated from estonian.sbl by Snowball 3.0.1 - https://snowballstem.org/

/**@constructor*/
var EstonianStemmer = function() {
    var base = new BaseStemmer();

    /** @const */ var a_0 = [
        ["gi", -1, 1],
        ["ki", -1, 2]
    ];

    /** @const */ var a_1 = [
        ["da", -1, 3],
        ["mata", -1, 1],
        ["b", -1, 3],
        ["ksid", -1, 1],
        ["nuksid", 3, 1],
        ["me", -1, 3],
        ["sime", 5, 1],
        ["ksime", 6, 1],
        ["nuksime", 7, 1],
        ["akse", -1, 2],
        ["dakse", 9, 1],
        ["takse", 9, 1],
        ["site", -1, 1],
        ["ksite", 12, 1],
        ["nuksite", 13, 1],
        ["n", -1, 3],
        ["sin", 15, 1],
        ["ksin", 16, 1],
        ["nuksin", 17, 1],
        ["daks", -1, 1],
        ["taks", -1, 1]
    ];

    /** @const */ var a_2 = [
        ["aa", -1, -1],
        ["ee", -1, -1],
        ["ii", -1, -1],
        ["oo", -1, -1],
        ["uu", -1, -1],
        ["\u00E4\u00E4", -1, -1],
        ["\u00F5\u00F5", -1, -1],
        ["\u00F6\u00F6", -1, -1],
        ["\u00FC\u00FC", -1, -1]
    ];

    /** @const */ var a_3 = [
        ["i", -1, 1]
    ];

    /** @const */ var a_4 = [
        ["lane", -1, 1],
        ["line", -1, 3],
        ["mine", -1, 2],
        ["lasse", -1, 1],
        ["lisse", -1, 3],
        ["misse", -1, 2],
        ["lasi", -1, 1],
        ["lisi", -1, 3],
        ["misi", -1, 2],
        ["last", -1, 1],
        ["list", -1, 3],
        ["mist", -1, 2]
    ];

    /** @const */ var a_5 = [
        ["ga", -1, 1],
        ["ta", -1, 1],
        ["le", -1, 1],
        ["sse", -1, 1],
        ["l", -1, 1],
        ["s", -1, 1],
        ["ks", 5, 1],
        ["t", -1, 2],
        ["lt", 7, 1],
        ["st", 7, 1]
    ];

    /** @const */ var a_6 = [
        ["", -1, 2],
        ["las", 0, 1],
        ["lis", 0, 1],
        ["mis", 0, 1],
        ["t", 0, -1]
    ];

    /** @const */ var a_7 = [
        ["d", -1, 4],
        ["sid", 0, 2],
        ["de", -1, 4],
        ["ikkude", 2, 1],
        ["ike", -1, 1],
        ["ikke", -1, 1],
        ["te", -1, 3]
    ];

    /** @const */ var a_8 = [
        ["va", -1, -1],
        ["du", -1, -1],
        ["nu", -1, -1],
        ["tu", -1, -1]
    ];

    /** @const */ var a_9 = [
        ["kk", -1, 1],
        ["pp", -1, 2],
        ["tt", -1, 3]
    ];

    /** @const */ var a_10 = [
        ["ma", -1, 2],
        ["mai", -1, 1],
        ["m", -1, 1]
    ];

    /** @const */ var a_11 = [
        ["joob", -1, 1],
        ["jood", -1, 1],
        ["joodakse", 1, 1],
        ["jooma", -1, 1],
        ["joomata", 3, 1],
        ["joome", -1, 1],
        ["joon", -1, 1],
        ["joote", -1, 1],
        ["joovad", -1, 1],
        ["juua", -1, 1],
        ["juuakse", 9, 1],
        ["j\u00E4i", -1, 12],
        ["j\u00E4id", 11, 12],
        ["j\u00E4ime", 11, 12],
        ["j\u00E4in", 11, 12],
        ["j\u00E4ite", 11, 12],
        ["j\u00E4\u00E4b", -1, 12],
        ["j\u00E4\u00E4d", -1, 12],
        ["j\u00E4\u00E4da", 17, 12],
        ["j\u00E4\u00E4dakse", 18, 12],
        ["j\u00E4\u00E4di", 17, 12],
        ["j\u00E4\u00E4ks", -1, 12],
        ["j\u00E4\u00E4ksid", 21, 12],
        ["j\u00E4\u00E4ksime", 21, 12],
        ["j\u00E4\u00E4ksin", 21, 12],
        ["j\u00E4\u00E4ksite", 21, 12],
        ["j\u00E4\u00E4ma", -1, 12],
        ["j\u00E4\u00E4mata", 26, 12],
        ["j\u00E4\u00E4me", -1, 12],
        ["j\u00E4\u00E4n", -1, 12],
        ["j\u00E4\u00E4te", -1, 12],
        ["j\u00E4\u00E4vad", -1, 12],
        ["j\u00F5i", -1, 1],
        ["j\u00F5id", 32, 1],
        ["j\u00F5ime", 32, 1],
        ["j\u00F5in", 32, 1],
        ["j\u00F5ite", 32, 1],
        ["keeb", -1, 4],
        ["keed", -1, 4],
        ["keedakse", 38, 4],
        ["keeks", -1, 4],
        ["keeksid", 40, 4],
        ["keeksime", 40, 4],
        ["keeksin", 40, 4],
        ["keeksite", 40, 4],
        ["keema", -1, 4],
        ["keemata", 45, 4],
        ["keeme", -1, 4],
        ["keen", -1, 4],
        ["kees", -1, 4],
        ["keeta", -1, 4],
        ["keete", -1, 4],
        ["keevad", -1, 4],
        ["k\u00E4ia", -1, 8],
        ["k\u00E4iakse", 53, 8],
        ["k\u00E4ib", -1, 8],
        ["k\u00E4id", -1, 8],
        ["k\u00E4idi", 56, 8],
        ["k\u00E4iks", -1, 8],
        ["k\u00E4iksid", 58, 8],
        ["k\u00E4iksime", 58, 8],
        ["k\u00E4iksin", 58, 8],
        ["k\u00E4iksite", 58, 8],
        ["k\u00E4ima", -1, 8],
        ["k\u00E4imata", 63, 8],
        ["k\u00E4ime", -1, 8],
        ["k\u00E4in", -1, 8],
        ["k\u00E4is", -1, 8],
        ["k\u00E4ite", -1, 8],
        ["k\u00E4ivad", -1, 8],
        ["laob", -1, 16],
        ["laod", -1, 16],
        ["laoks", -1, 16],
        ["laoksid", 72, 16],
        ["laoksime", 72, 16],
        ["laoksin", 72, 16],
        ["laoksite", 72, 16],
        ["laome", -1, 16],
        ["laon", -1, 16],
        ["laote", -1, 16],
        ["laovad", -1, 16],
        ["loeb", -1, 14],
        ["loed", -1, 14],
        ["loeks", -1, 14],
        ["loeksid", 83, 14],
        ["loeksime", 83, 14],
        ["loeksin", 83, 14],
        ["loeksite", 83, 14],
        ["loeme", -1, 14],
        ["loen", -1, 14],
        ["loete", -1, 14],
        ["loevad", -1, 14],
        ["loob", -1, 7],
        ["lood", -1, 7],
        ["loodi", 93, 7],
        ["looks", -1, 7],
        ["looksid", 95, 7],
        ["looksime", 95, 7],
        ["looksin", 95, 7],
        ["looksite", 95, 7],
        ["looma", -1, 7],
        ["loomata", 100, 7],
        ["loome", -1, 7],
        ["loon", -1, 7],
        ["loote", -1, 7],
        ["loovad", -1, 7],
        ["luua", -1, 7],
        ["luuakse", 106, 7],
        ["l\u00F5i", -1, 6],
        ["l\u00F5id", 108, 6],
        ["l\u00F5ime", 108, 6],
        ["l\u00F5in", 108, 6],
        ["l\u00F5ite", 108, 6],
        ["l\u00F6\u00F6b", -1, 5],
        ["l\u00F6\u00F6d", -1, 5],
        ["l\u00F6\u00F6dakse", 114, 5],
        ["l\u00F6\u00F6di", 114, 5],
        ["l\u00F6\u00F6ks", -1, 5],
        ["l\u00F6\u00F6ksid", 117, 5],
        ["l\u00F6\u00F6ksime", 117, 5],
        ["l\u00F6\u00F6ksin", 117, 5],
        ["l\u00F6\u00F6ksite", 117, 5],
        ["l\u00F6\u00F6ma", -1, 5],
        ["l\u00F6\u00F6mata", 122, 5],
        ["l\u00F6\u00F6me", -1, 5],
        ["l\u00F6\u00F6n", -1, 5],
        ["l\u00F6\u00F6te", -1, 5],
        ["l\u00F6\u00F6vad", -1, 5],
        ["l\u00FC\u00FCa", -1, 5],
        ["l\u00FC\u00FCakse", 128, 5],
        ["m\u00FC\u00FCa", -1, 13],
        ["m\u00FC\u00FCakse", 130, 13],
        ["m\u00FC\u00FCb", -1, 13],
        ["m\u00FC\u00FCd", -1, 13],
        ["m\u00FC\u00FCdi", 133, 13],
        ["m\u00FC\u00FCks", -1, 13],
        ["m\u00FC\u00FCksid", 135, 13],
        ["m\u00FC\u00FCksime", 135, 13],
        ["m\u00FC\u00FCksin", 135, 13],
        ["m\u00FC\u00FCksite", 135, 13],
        ["m\u00FC\u00FCma", -1, 13],
        ["m\u00FC\u00FCmata", 140, 13],
        ["m\u00FC\u00FCme", -1, 13],
        ["m\u00FC\u00FCn", -1, 13],
        ["m\u00FC\u00FCs", -1, 13],
        ["m\u00FC\u00FCte", -1, 13],
        ["m\u00FC\u00FCvad", -1, 13],
        ["n\u00E4eb", -1, 18],
        ["n\u00E4ed", -1, 18],
        ["n\u00E4eks", -1, 18],
        ["n\u00E4eksid", 149, 18],
        ["n\u00E4eksime", 149, 18],
        ["n\u00E4eksin", 149, 18],
        ["n\u00E4eksite", 149, 18],
        ["n\u00E4eme", -1, 18],
        ["n\u00E4en", -1, 18],
        ["n\u00E4ete", -1, 18],
        ["n\u00E4evad", -1, 18],
        ["n\u00E4gema", -1, 18],
        ["n\u00E4gemata", 158, 18],
        ["n\u00E4ha", -1, 18],
        ["n\u00E4hakse", 160, 18],
        ["n\u00E4hti", -1, 18],
        ["p\u00F5eb", -1, 15],
        ["p\u00F5ed", -1, 15],
        ["p\u00F5eks", -1, 15],
        ["p\u00F5eksid", 165, 15],
        ["p\u00F5eksime", 165, 15],
        ["p\u00F5eksin", 165, 15],
        ["p\u00F5eksite", 165, 15],
        ["p\u00F5eme", -1, 15],
        ["p\u00F5en", -1, 15],
        ["p\u00F5ete", -1, 15],
        ["p\u00F5evad", -1, 15],
        ["saab", -1, 2],
        ["saad", -1, 2],
        ["saada", 175, 2],
        ["saadakse", 176, 2],
        ["saadi", 175, 2],
        ["saaks", -1, 2],
        ["saaksid", 179, 2],
        ["saaksime", 179, 2],
        ["saaksin", 179, 2],
        ["saaksite", 179, 2],
        ["saama", -1, 2],
        ["saamata", 184, 2],
        ["saame", -1, 2],
        ["saan", -1, 2],
        ["saate", -1, 2],
        ["saavad", -1, 2],
        ["sai", -1, 2],
        ["said", 190, 2],
        ["saime", 190, 2],
        ["sain", 190, 2],
        ["saite", 190, 2],
        ["s\u00F5i", -1, 9],
        ["s\u00F5id", 195, 9],
        ["s\u00F5ime", 195, 9],
        ["s\u00F5in", 195, 9],
        ["s\u00F5ite", 195, 9],
        ["s\u00F6\u00F6b", -1, 9],
        ["s\u00F6\u00F6d", -1, 9],
        ["s\u00F6\u00F6dakse", 201, 9],
        ["s\u00F6\u00F6di", 201, 9],
        ["s\u00F6\u00F6ks", -1, 9],
        ["s\u00F6\u00F6ksid", 204, 9],
        ["s\u00F6\u00F6ksime", 204, 9],
        ["s\u00F6\u00F6ksin", 204, 9],
        ["s\u00F6\u00F6ksite", 204, 9],
        ["s\u00F6\u00F6ma", -1, 9],
        ["s\u00F6\u00F6mata", 209, 9],
        ["s\u00F6\u00F6me", -1, 9],
        ["s\u00F6\u00F6n", -1, 9],
        ["s\u00F6\u00F6te", -1, 9],
        ["s\u00F6\u00F6vad", -1, 9],
        ["s\u00FC\u00FCa", -1, 9],
        ["s\u00FC\u00FCakse", 215, 9],
        ["teeb", -1, 17],
        ["teed", -1, 17],
        ["teeks", -1, 17],
        ["teeksid", 219, 17],
        ["teeksime", 219, 17],
        ["teeksin", 219, 17],
        ["teeksite", 219, 17],
        ["teeme", -1, 17],
        ["teen", -1, 17],
        ["teete", -1, 17],
        ["teevad", -1, 17],
        ["tegema", -1, 17],
        ["tegemata", 228, 17],
        ["teha", -1, 17],
        ["tehakse", 230, 17],
        ["tehti", -1, 17],
        ["toob", -1, 10],
        ["tood", -1, 10],
        ["toodi", 234, 10],
        ["tooks", -1, 10],
        ["tooksid", 236, 10],
        ["tooksime", 236, 10],
        ["tooksin", 236, 10],
        ["tooksite", 236, 10],
        ["tooma", -1, 10],
        ["toomata", 241, 10],
        ["toome", -1, 10],
        ["toon", -1, 10],
        ["toote", -1, 10],
        ["toovad", -1, 10],
        ["tuua", -1, 10],
        ["tuuakse", 247, 10],
        ["t\u00F5i", -1, 10],
        ["t\u00F5id", 249, 10],
        ["t\u00F5ime", 249, 10],
        ["t\u00F5in", 249, 10],
        ["t\u00F5ite", 249, 10],
        ["viia", -1, 3],
        ["viiakse", 254, 3],
        ["viib", -1, 3],
        ["viid", -1, 3],
        ["viidi", 257, 3],
        ["viiks", -1, 3],
        ["viiksid", 259, 3],
        ["viiksime", 259, 3],
        ["viiksin", 259, 3],
        ["viiksite", 259, 3],
        ["viima", -1, 3],
        ["viimata", 264, 3],
        ["viime", -1, 3],
        ["viin", -1, 3],
        ["viisime", -1, 3],
        ["viisin", -1, 3],
        ["viisite", -1, 3],
        ["viite", -1, 3],
        ["viivad", -1, 3],
        ["v\u00F5ib", -1, 11],
        ["v\u00F5id", -1, 11],
        ["v\u00F5ida", 274, 11],
        ["v\u00F5idakse", 275, 11],
        ["v\u00F5idi", 274, 11],
        ["v\u00F5iks", -1, 11],
        ["v\u00F5iksid", 278, 11],
        ["v\u00F5iksime", 278, 11],
        ["v\u00F5iksin", 278, 11],
        ["v\u00F5iksite", 278, 11],
        ["v\u00F5ima", -1, 11],
        ["v\u00F5imata", 283, 11],
        ["v\u00F5ime", -1, 11],
        ["v\u00F5in", -1, 11],
        ["v\u00F5is", -1, 11],
        ["v\u00F5ite", -1, 11],
        ["v\u00F5ivad", -1, 11]
    ];

    /** @const */ var /** Array<int> */ g_V1 = [17, 65, 16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 0, 48, 8];

    /** @const */ var /** Array<int> */ g_RV = [17, 65, 16];

    /** @const */ var /** Array<int> */ g_KI = [117, 66, 6, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 128, 0, 0, 0, 16];

    /** @const */ var /** Array<int> */ g_GI = [21, 123, 243, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 0, 48, 8];

    var /** number */ I_p1 = 0;


    /** @return {boolean} */
    function r_mark_regions() {
        I_p1 = base.limit;
        if (!base.go_out_grouping(g_V1, 97, 252))
        {
            return false;
        }
        base.cursor++;
        if (!base.go_in_grouping(g_V1, 97, 252))
        {
            return false;
        }
        base.cursor++;
        I_p1 = base.cursor;
        return true;
    };

    /** @return {boolean} */
    function r_emphasis() {
        var /** number */ among_var;
        if (base.cursor < I_p1)
        {
            return false;
        }
        /** @const */ var /** number */ v_1 = base.limit_backward;
        base.limit_backward = I_p1;
        base.ket = base.cursor;
        among_var = base.find_among_b(a_0);
        if (among_var == 0)
        {
            base.limit_backward = v_1;
            return false;
        }
        base.bra = base.cursor;
        base.limit_backward = v_1;
        /** @const */ var /** number */ v_2 = base.limit - base.cursor;
        {
            /** @const */ var /** number */ c1 = base.cursor - 4;
            if (c1 < base.limit_backward)
            {
                return false;
            }
            base.cursor = c1;
        }
        base.cursor = base.limit - v_2;
        switch (among_var) {
            case 1:
                /** @const */ var /** number */ v_3 = base.limit - base.cursor;
                if (!(base.in_grouping_b(g_GI, 97, 252)))
                {
                    return false;
                }
                base.cursor = base.limit - v_3;
                {
                    /** @const */ var /** number */ v_4 = base.limit - base.cursor;
                    lab0: {
                        if (!r_LONGV())
                        {
                            break lab0;
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
            case 2:
                if (!(base.in_grouping_b(g_KI, 98, 382)))
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
    function r_verb() {
        var /** number */ among_var;
        if (base.cursor < I_p1)
        {
            return false;
        }
        /** @const */ var /** number */ v_1 = base.limit_backward;
        base.limit_backward = I_p1;
        base.ket = base.cursor;
        among_var = base.find_among_b(a_1);
        if (among_var == 0)
        {
            base.limit_backward = v_1;
            return false;
        }
        base.bra = base.cursor;
        base.limit_backward = v_1;
        switch (among_var) {
            case 1:
                if (!base.slice_del())
                {
                    return false;
                }
                break;
            case 2:
                if (!base.slice_from("a"))
                {
                    return false;
                }
                break;
            case 3:
                if (!(base.in_grouping_b(g_V1, 97, 252)))
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
    function r_LONGV() {
        if (base.find_among_b(a_2) == 0)
        {
            return false;
        }
        return true;
    };

    /** @return {boolean} */
    function r_i_plural() {
        if (base.cursor < I_p1)
        {
            return false;
        }
        /** @const */ var /** number */ v_1 = base.limit_backward;
        base.limit_backward = I_p1;
        base.ket = base.cursor;
        if (base.find_among_b(a_3) == 0)
        {
            base.limit_backward = v_1;
            return false;
        }
        base.bra = base.cursor;
        base.limit_backward = v_1;
        if (!(base.in_grouping_b(g_RV, 97, 117)))
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
    function r_special_noun_endings() {
        var /** number */ among_var;
        if (base.cursor < I_p1)
        {
            return false;
        }
        /** @const */ var /** number */ v_1 = base.limit_backward;
        base.limit_backward = I_p1;
        base.ket = base.cursor;
        among_var = base.find_among_b(a_4);
        if (among_var == 0)
        {
            base.limit_backward = v_1;
            return false;
        }
        base.bra = base.cursor;
        base.limit_backward = v_1;
        switch (among_var) {
            case 1:
                if (!base.slice_from("lase"))
                {
                    return false;
                }
                break;
            case 2:
                if (!base.slice_from("mise"))
                {
                    return false;
                }
                break;
            case 3:
                if (!base.slice_from("lise"))
                {
                    return false;
                }
                break;
        }
        return true;
    };

    /** @return {boolean} */
    function r_case_ending() {
        var /** number */ among_var;
        if (base.cursor < I_p1)
        {
            return false;
        }
        /** @const */ var /** number */ v_1 = base.limit_backward;
        base.limit_backward = I_p1;
        base.ket = base.cursor;
        among_var = base.find_among_b(a_5);
        if (among_var == 0)
        {
            base.limit_backward = v_1;
            return false;
        }
        base.bra = base.cursor;
        base.limit_backward = v_1;
        switch (among_var) {
            case 1:
                lab0: {
                    /** @const */ var /** number */ v_2 = base.limit - base.cursor;
                    lab1: {
                        if (!(base.in_grouping_b(g_RV, 97, 117)))
                        {
                            break lab1;
                        }
                        break lab0;
                    }
                    base.cursor = base.limit - v_2;
                    if (!r_LONGV())
                    {
                        return false;
                    }
                }
                break;
            case 2:
                /** @const */ var /** number */ v_3 = base.limit - base.cursor;
                {
                    /** @const */ var /** number */ c1 = base.cursor - 4;
                    if (c1 < base.limit_backward)
                    {
                        return false;
                    }
                    base.cursor = c1;
                }
                base.cursor = base.limit - v_3;
                break;
        }
        if (!base.slice_del())
        {
            return false;
        }
        return true;
    };

    /** @return {boolean} */
    function r_plural_three_first_cases() {
        var /** number */ among_var;
        if (base.cursor < I_p1)
        {
            return false;
        }
        /** @const */ var /** number */ v_1 = base.limit_backward;
        base.limit_backward = I_p1;
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
                if (!base.slice_from("iku"))
                {
                    return false;
                }
                break;
            case 2:
                {
                    /** @const */ var /** number */ v_2 = base.limit - base.cursor;
                    lab0: {
                        if (!r_LONGV())
                        {
                            break lab0;
                        }
                        return false;
                    }
                    base.cursor = base.limit - v_2;
                }
                if (!base.slice_del())
                {
                    return false;
                }
                break;
            case 3:
                lab1: {
                    /** @const */ var /** number */ v_3 = base.limit - base.cursor;
                    lab2: {
                        /** @const */ var /** number */ v_4 = base.limit - base.cursor;
                        {
                            /** @const */ var /** number */ c1 = base.cursor - 4;
                            if (c1 < base.limit_backward)
                            {
                                break lab2;
                            }
                            base.cursor = c1;
                        }
                        base.cursor = base.limit - v_4;
                        among_var = base.find_among_b(a_6);
                        switch (among_var) {
                            case 1:
                                if (!base.slice_from("e"))
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
                        }
                        break lab1;
                    }
                    base.cursor = base.limit - v_3;
                    if (!base.slice_from("t"))
                    {
                        return false;
                    }
                }
                break;
            case 4:
                lab3: {
                    /** @const */ var /** number */ v_5 = base.limit - base.cursor;
                    lab4: {
                        if (!(base.in_grouping_b(g_RV, 97, 117)))
                        {
                            break lab4;
                        }
                        break lab3;
                    }
                    base.cursor = base.limit - v_5;
                    if (!r_LONGV())
                    {
                        return false;
                    }
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
    function r_nu() {
        if (base.cursor < I_p1)
        {
            return false;
        }
        /** @const */ var /** number */ v_1 = base.limit_backward;
        base.limit_backward = I_p1;
        base.ket = base.cursor;
        if (base.find_among_b(a_8) == 0)
        {
            base.limit_backward = v_1;
            return false;
        }
        base.bra = base.cursor;
        base.limit_backward = v_1;
        if (!base.slice_del())
        {
            return false;
        }
        return true;
    };

    /** @return {boolean} */
    function r_undouble_kpt() {
        var /** number */ among_var;
        if (!(base.in_grouping_b(g_V1, 97, 252)))
        {
            return false;
        }
        if (I_p1 > base.cursor)
        {
            return false;
        }
        base.ket = base.cursor;
        among_var = base.find_among_b(a_9);
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
                if (!base.slice_from("p"))
                {
                    return false;
                }
                break;
            case 3:
                if (!base.slice_from("t"))
                {
                    return false;
                }
                break;
        }
        return true;
    };

    /** @return {boolean} */
    function r_degrees() {
        var /** number */ among_var;
        if (base.cursor < I_p1)
        {
            return false;
        }
        /** @const */ var /** number */ v_1 = base.limit_backward;
        base.limit_backward = I_p1;
        base.ket = base.cursor;
        among_var = base.find_among_b(a_10);
        if (among_var == 0)
        {
            base.limit_backward = v_1;
            return false;
        }
        base.bra = base.cursor;
        base.limit_backward = v_1;
        switch (among_var) {
            case 1:
                if (!(base.in_grouping_b(g_RV, 97, 117)))
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
        }
        return true;
    };

    /** @return {boolean} */
    function r_substantive() {
        /** @const */ var /** number */ v_1 = base.limit - base.cursor;
        r_special_noun_endings();
        base.cursor = base.limit - v_1;
        /** @const */ var /** number */ v_2 = base.limit - base.cursor;
        r_case_ending();
        base.cursor = base.limit - v_2;
        /** @const */ var /** number */ v_3 = base.limit - base.cursor;
        r_plural_three_first_cases();
        base.cursor = base.limit - v_3;
        /** @const */ var /** number */ v_4 = base.limit - base.cursor;
        r_degrees();
        base.cursor = base.limit - v_4;
        /** @const */ var /** number */ v_5 = base.limit - base.cursor;
        r_i_plural();
        base.cursor = base.limit - v_5;
        /** @const */ var /** number */ v_6 = base.limit - base.cursor;
        r_nu();
        base.cursor = base.limit - v_6;
        return true;
    };

    /** @return {boolean} */
    function r_verb_exceptions() {
        var /** number */ among_var;
        base.bra = base.cursor;
        among_var = base.find_among(a_11);
        if (among_var == 0)
        {
            return false;
        }
        base.ket = base.cursor;
        if (base.cursor < base.limit)
        {
            return false;
        }
        switch (among_var) {
            case 1:
                if (!base.slice_from("joo"))
                {
                    return false;
                }
                break;
            case 2:
                if (!base.slice_from("saa"))
                {
                    return false;
                }
                break;
            case 3:
                if (!base.slice_from("viima"))
                {
                    return false;
                }
                break;
            case 4:
                if (!base.slice_from("keesi"))
                {
                    return false;
                }
                break;
            case 5:
                if (!base.slice_from("l\u00F6\u00F6"))
                {
                    return false;
                }
                break;
            case 6:
                if (!base.slice_from("l\u00F5i"))
                {
                    return false;
                }
                break;
            case 7:
                if (!base.slice_from("loo"))
                {
                    return false;
                }
                break;
            case 8:
                if (!base.slice_from("k\u00E4isi"))
                {
                    return false;
                }
                break;
            case 9:
                if (!base.slice_from("s\u00F6\u00F6"))
                {
                    return false;
                }
                break;
            case 10:
                if (!base.slice_from("too"))
                {
                    return false;
                }
                break;
            case 11:
                if (!base.slice_from("v\u00F5isi"))
                {
                    return false;
                }
                break;
            case 12:
                if (!base.slice_from("j\u00E4\u00E4ma"))
                {
                    return false;
                }
                break;
            case 13:
                if (!base.slice_from("m\u00FC\u00FCsi"))
                {
                    return false;
                }
                break;
            case 14:
                if (!base.slice_from("luge"))
                {
                    return false;
                }
                break;
            case 15:
                if (!base.slice_from("p\u00F5de"))
                {
                    return false;
                }
                break;
            case 16:
                if (!base.slice_from("ladu"))
                {
                    return false;
                }
                break;
            case 17:
                if (!base.slice_from("tegi"))
                {
                    return false;
                }
                break;
            case 18:
                if (!base.slice_from("n\u00E4gi"))
                {
                    return false;
                }
                break;
        }
        return true;
    };

    this.stem = /** @return {boolean} */ function() {
        {
            /** @const */ var /** number */ v_1 = base.cursor;
            lab0: {
                if (!r_verb_exceptions())
                {
                    break lab0;
                }
                return false;
            }
            base.cursor = v_1;
        }
        /** @const */ var /** number */ v_2 = base.cursor;
        r_mark_regions();
        base.cursor = v_2;
        base.limit_backward = base.cursor; base.cursor = base.limit;
        /** @const */ var /** number */ v_3 = base.limit - base.cursor;
        r_emphasis();
        base.cursor = base.limit - v_3;
        /** @const */ var /** number */ v_4 = base.limit - base.cursor;
        lab1: {
            lab2: {
                /** @const */ var /** number */ v_5 = base.limit - base.cursor;
                lab3: {
                    if (!r_verb())
                    {
                        break lab3;
                    }
                    break lab2;
                }
                base.cursor = base.limit - v_5;
                r_substantive();
            }
        }
        base.cursor = base.limit - v_4;
        /** @const */ var /** number */ v_6 = base.limit - base.cursor;
        r_undouble_kpt();
        base.cursor = base.limit - v_6;
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
