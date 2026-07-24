// Generated from irish.sbl by Snowball 3.0.1 - https://snowballstem.org/

/**@constructor*/
var IrishStemmer = function() {
    var base = new BaseStemmer();

    /** @const */ var a_0 = [
        ["b'", -1, 1],
        ["bh", -1, 4],
        ["bhf", 1, 2],
        ["bp", -1, 8],
        ["ch", -1, 5],
        ["d'", -1, 1],
        ["d'fh", 5, 2],
        ["dh", -1, 6],
        ["dt", -1, 9],
        ["fh", -1, 2],
        ["gc", -1, 5],
        ["gh", -1, 7],
        ["h-", -1, 1],
        ["m'", -1, 1],
        ["mb", -1, 4],
        ["mh", -1, 10],
        ["n-", -1, 1],
        ["nd", -1, 6],
        ["ng", -1, 7],
        ["ph", -1, 8],
        ["sh", -1, 3],
        ["t-", -1, 1],
        ["th", -1, 9],
        ["ts", -1, 3]
    ];

    /** @const */ var a_1 = [
        ["\u00EDochta", -1, 1],
        ["a\u00EDochta", 0, 1],
        ["ire", -1, 2],
        ["aire", 2, 2],
        ["abh", -1, 1],
        ["eabh", 4, 1],
        ["ibh", -1, 1],
        ["aibh", 6, 1],
        ["amh", -1, 1],
        ["eamh", 8, 1],
        ["imh", -1, 1],
        ["aimh", 10, 1],
        ["\u00EDocht", -1, 1],
        ["a\u00EDocht", 12, 1],
        ["ir\u00ED", -1, 2],
        ["air\u00ED", 14, 2]
    ];

    /** @const */ var a_2 = [
        ["\u00F3ideacha", -1, 6],
        ["patacha", -1, 5],
        ["achta", -1, 1],
        ["arcachta", 2, 2],
        ["eachta", 2, 1],
        ["grafa\u00EDochta", -1, 4],
        ["paite", -1, 5],
        ["ach", -1, 1],
        ["each", 7, 1],
        ["\u00F3ideach", 8, 6],
        ["gineach", 8, 3],
        ["patach", 7, 5],
        ["grafa\u00EDoch", -1, 4],
        ["pataigh", -1, 5],
        ["\u00F3idigh", -1, 6],
        ["acht\u00FAil", -1, 1],
        ["eacht\u00FAil", 15, 1],
        ["gineas", -1, 3],
        ["ginis", -1, 3],
        ["acht", -1, 1],
        ["arcacht", 19, 2],
        ["eacht", 19, 1],
        ["grafa\u00EDocht", -1, 4],
        ["arcachta\u00ED", -1, 2],
        ["grafa\u00EDochta\u00ED", -1, 4]
    ];

    /** @const */ var a_3 = [
        ["imid", -1, 1],
        ["aimid", 0, 1],
        ["\u00EDmid", -1, 1],
        ["a\u00EDmid", 2, 1],
        ["adh", -1, 2],
        ["eadh", 4, 2],
        ["faidh", -1, 1],
        ["fidh", -1, 1],
        ["\u00E1il", -1, 2],
        ["ain", -1, 2],
        ["tear", -1, 2],
        ["tar", -1, 2]
    ];

    /** @const */ var /** Array<int> */ g_v = [17, 65, 16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 17, 4, 2];

    var /** number */ I_p2 = 0;
    var /** number */ I_p1 = 0;
    var /** number */ I_pV = 0;


    /** @return {boolean} */
    function r_mark_regions() {
        I_pV = base.limit;
        I_p1 = base.limit;
        I_p2 = base.limit;
        /** @const */ var /** number */ v_1 = base.cursor;
        lab0: {
            if (!base.go_out_grouping(g_v, 97, 250))
            {
                break lab0;
            }
            base.cursor++;
            I_pV = base.cursor;
            if (!base.go_in_grouping(g_v, 97, 250))
            {
                break lab0;
            }
            base.cursor++;
            I_p1 = base.cursor;
            if (!base.go_out_grouping(g_v, 97, 250))
            {
                break lab0;
            }
            base.cursor++;
            if (!base.go_in_grouping(g_v, 97, 250))
            {
                break lab0;
            }
            base.cursor++;
            I_p2 = base.cursor;
        }
        base.cursor = v_1;
        return true;
    };

    /** @return {boolean} */
    function r_initial_morph() {
        var /** number */ among_var;
        base.bra = base.cursor;
        among_var = base.find_among(a_0);
        if (among_var == 0)
        {
            return false;
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
                if (!base.slice_from("f"))
                {
                    return false;
                }
                break;
            case 3:
                if (!base.slice_from("s"))
                {
                    return false;
                }
                break;
            case 4:
                if (!base.slice_from("b"))
                {
                    return false;
                }
                break;
            case 5:
                if (!base.slice_from("c"))
                {
                    return false;
                }
                break;
            case 6:
                if (!base.slice_from("d"))
                {
                    return false;
                }
                break;
            case 7:
                if (!base.slice_from("g"))
                {
                    return false;
                }
                break;
            case 8:
                if (!base.slice_from("p"))
                {
                    return false;
                }
                break;
            case 9:
                if (!base.slice_from("t"))
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
    function r_noun_sfx() {
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
                if (!r_R1())
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
                break;
        }
        return true;
    };

    /** @return {boolean} */
    function r_deriv() {
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
                if (!base.slice_from("arc"))
                {
                    return false;
                }
                break;
            case 3:
                if (!base.slice_from("gin"))
                {
                    return false;
                }
                break;
            case 4:
                if (!base.slice_from("graf"))
                {
                    return false;
                }
                break;
            case 5:
                if (!base.slice_from("paite"))
                {
                    return false;
                }
                break;
            case 6:
                if (!base.slice_from("\u00F3id"))
                {
                    return false;
                }
                break;
        }
        return true;
    };

    /** @return {boolean} */
    function r_verb_sfx() {
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
                if (!r_RV())
                {
                    return false;
                }
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
                if (!base.slice_del())
                {
                    return false;
                }
                break;
        }
        return true;
    };

    this.stem = /** @return {boolean} */ function() {
        /** @const */ var /** number */ v_1 = base.cursor;
        r_initial_morph();
        base.cursor = v_1;
        r_mark_regions();
        base.limit_backward = base.cursor; base.cursor = base.limit;
        /** @const */ var /** number */ v_2 = base.limit - base.cursor;
        r_noun_sfx();
        base.cursor = base.limit - v_2;
        /** @const */ var /** number */ v_3 = base.limit - base.cursor;
        r_deriv();
        base.cursor = base.limit - v_3;
        /** @const */ var /** number */ v_4 = base.limit - base.cursor;
        r_verb_sfx();
        base.cursor = base.limit - v_4;
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
