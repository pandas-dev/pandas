// Generated from swedish.sbl by Snowball 3.0.1 - https://snowballstem.org/

/**@constructor*/
var SwedishStemmer = function() {
    var base = new BaseStemmer();

    /** @const */ var a_0 = [
        ["fab", -1, -1],
        ["h", -1, -1],
        ["pak", -1, -1],
        ["rak", -1, -1],
        ["stak", -1, -1],
        ["kom", -1, -1],
        ["iet", -1, -1],
        ["cit", -1, -1],
        ["dit", -1, -1],
        ["alit", -1, -1],
        ["ilit", -1, -1],
        ["mit", -1, -1],
        ["nit", -1, -1],
        ["pit", -1, -1],
        ["rit", -1, -1],
        ["sit", -1, -1],
        ["tit", -1, -1],
        ["uit", -1, -1],
        ["ivit", -1, -1],
        ["kvit", -1, -1],
        ["xit", -1, -1]
    ];

    /** @const */ var a_1 = [
        ["a", -1, 1],
        ["arna", 0, 1],
        ["erna", 0, 1],
        ["heterna", 2, 1],
        ["orna", 0, 1],
        ["ad", -1, 1],
        ["e", -1, 1],
        ["ade", 6, 1],
        ["ande", 6, 1],
        ["arne", 6, 1],
        ["are", 6, 1],
        ["aste", 6, 1],
        ["en", -1, 1],
        ["anden", 12, 1],
        ["aren", 12, 1],
        ["heten", 12, 1],
        ["ern", -1, 1],
        ["ar", -1, 1],
        ["er", -1, 1],
        ["heter", 18, 1],
        ["or", -1, 1],
        ["s", -1, 2],
        ["as", 21, 1],
        ["arnas", 22, 1],
        ["ernas", 22, 1],
        ["ornas", 22, 1],
        ["es", 21, 1],
        ["ades", 26, 1],
        ["andes", 26, 1],
        ["ens", 21, 1],
        ["arens", 29, 1],
        ["hetens", 29, 1],
        ["erns", 21, 1],
        ["at", -1, 1],
        ["et", -1, 3],
        ["andet", 34, 1],
        ["het", 34, 1],
        ["ast", -1, 1]
    ];

    /** @const */ var a_2 = [
        ["dd", -1, -1],
        ["gd", -1, -1],
        ["nn", -1, -1],
        ["dt", -1, -1],
        ["gt", -1, -1],
        ["kt", -1, -1],
        ["tt", -1, -1]
    ];

    /** @const */ var a_3 = [
        ["ig", -1, 1],
        ["lig", 0, 1],
        ["els", -1, 1],
        ["fullt", -1, 3],
        ["\u00F6st", -1, 2]
    ];

    /** @const */ var /** Array<int> */ g_v = [17, 65, 16, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 24, 0, 32];

    /** @const */ var /** Array<int> */ g_s_ending = [119, 127, 149];

    /** @const */ var /** Array<int> */ g_ost_ending = [173, 58];

    var /** number */ I_x = 0;
    var /** number */ I_p1 = 0;


    /** @return {boolean} */
    function r_mark_regions() {
        I_p1 = base.limit;
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
        if (!base.go_out_grouping(g_v, 97, 246))
        {
            return false;
        }
        base.cursor++;
        if (!base.go_in_grouping(g_v, 97, 246))
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
        return true;
    };

    /** @return {boolean} */
    function r_et_condition() {
        /** @const */ var /** number */ v_1 = base.limit - base.cursor;
        if (!(base.out_grouping_b(g_v, 97, 246)))
        {
            return false;
        }
        if (!(base.in_grouping_b(g_v, 97, 246)))
        {
            return false;
        }
        lab0: {
            if (base.cursor > base.limit_backward)
            {
                break lab0;
            }
            return false;
        }
        base.cursor = base.limit - v_1;
        {
            /** @const */ var /** number */ v_2 = base.limit - base.cursor;
            lab1: {
                if (base.find_among_b(a_0) == 0)
                {
                    break lab1;
                }
                return false;
            }
            base.cursor = base.limit - v_2;
        }
        return true;
    };

    /** @return {boolean} */
    function r_main_suffix() {
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
                lab0: {
                    /** @const */ var /** number */ v_2 = base.limit - base.cursor;
                    lab1: {
                        if (!(base.eq_s_b("et")))
                        {
                            break lab1;
                        }
                        if (!r_et_condition())
                        {
                            break lab1;
                        }
                        base.bra = base.cursor;
                        break lab0;
                    }
                    base.cursor = base.limit - v_2;
                    if (!(base.in_grouping_b(g_s_ending, 98, 121)))
                    {
                        return false;
                    }
                }
                if (!base.slice_del())
                {
                    return false;
                }
                break;
            case 3:
                if (!r_et_condition())
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
    function r_consonant_pair() {
        if (base.cursor < I_p1)
        {
            return false;
        }
        /** @const */ var /** number */ v_1 = base.limit_backward;
        base.limit_backward = I_p1;
        /** @const */ var /** number */ v_2 = base.limit - base.cursor;
        if (base.find_among_b(a_2) == 0)
        {
            base.limit_backward = v_1;
            return false;
        }
        base.cursor = base.limit - v_2;
        base.ket = base.cursor;
        if (base.cursor <= base.limit_backward)
        {
            base.limit_backward = v_1;
            return false;
        }
        base.cursor--;
        base.bra = base.cursor;
        if (!base.slice_del())
        {
            return false;
        }
        base.limit_backward = v_1;
        return true;
    };

    /** @return {boolean} */
    function r_other_suffix() {
        var /** number */ among_var;
        if (base.cursor < I_p1)
        {
            return false;
        }
        /** @const */ var /** number */ v_1 = base.limit_backward;
        base.limit_backward = I_p1;
        base.ket = base.cursor;
        among_var = base.find_among_b(a_3);
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
                if (!(base.in_grouping_b(g_ost_ending, 105, 118)))
                {
                    return false;
                }
                if (!base.slice_from("\u00F6s"))
                {
                    return false;
                }
                break;
            case 3:
                if (!base.slice_from("full"))
                {
                    return false;
                }
                break;
        }
        return true;
    };

    this.stem = /** @return {boolean} */ function() {
        /** @const */ var /** number */ v_1 = base.cursor;
        r_mark_regions();
        base.cursor = v_1;
        base.limit_backward = base.cursor; base.cursor = base.limit;
        /** @const */ var /** number */ v_2 = base.limit - base.cursor;
        r_main_suffix();
        base.cursor = base.limit - v_2;
        /** @const */ var /** number */ v_3 = base.limit - base.cursor;
        r_consonant_pair();
        base.cursor = base.limit - v_3;
        /** @const */ var /** number */ v_4 = base.limit - base.cursor;
        r_other_suffix();
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
