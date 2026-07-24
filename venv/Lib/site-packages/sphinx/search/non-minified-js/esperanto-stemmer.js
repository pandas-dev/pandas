// Generated from esperanto.sbl by Snowball 3.0.1 - https://snowballstem.org/

/**@constructor*/
var EsperantoStemmer = function() {
    var base = new BaseStemmer();

    /** @const */ var a_0 = [
        ["", -1, 14],
        ["-", 0, 13],
        ["cx", 0, 1],
        ["gx", 0, 2],
        ["hx", 0, 3],
        ["jx", 0, 4],
        ["q", 0, 12],
        ["sx", 0, 5],
        ["ux", 0, 6],
        ["w", 0, 12],
        ["x", 0, 12],
        ["y", 0, 12],
        ["\u00E1", 0, 7],
        ["\u00E9", 0, 8],
        ["\u00ED", 0, 9],
        ["\u00F3", 0, 10],
        ["\u00FA", 0, 11]
    ];

    /** @const */ var a_1 = [
        ["as", -1, -1],
        ["i", -1, -1],
        ["is", 1, -1],
        ["os", -1, -1],
        ["u", -1, -1],
        ["us", 4, -1]
    ];

    /** @const */ var a_2 = [
        ["ci", -1, -1],
        ["gi", -1, -1],
        ["hi", -1, -1],
        ["li", -1, -1],
        ["ili", 3, -1],
        ["\u015Dli", 3, -1],
        ["mi", -1, -1],
        ["ni", -1, -1],
        ["oni", 7, -1],
        ["ri", -1, -1],
        ["si", -1, -1],
        ["vi", -1, -1],
        ["ivi", 11, -1],
        ["\u011Di", -1, -1],
        ["\u015Di", -1, -1],
        ["i\u015Di", 14, -1],
        ["mal\u015Di", 14, -1]
    ];

    /** @const */ var a_3 = [
        ["amb", -1, -1],
        ["bald", -1, -1],
        ["malbald", 1, -1],
        ["morg", -1, -1],
        ["postmorg", 3, -1],
        ["adi", -1, -1],
        ["hodi", -1, -1],
        ["ank", -1, -1],
        ["\u0109irk", -1, -1],
        ["tut\u0109irk", 8, -1],
        ["presk", -1, -1],
        ["almen", -1, -1],
        ["apen", -1, -1],
        ["hier", -1, -1],
        ["anta\u016Dhier", 13, -1],
        ["malgr", -1, -1],
        ["ankor", -1, -1],
        ["kontr", -1, -1],
        ["anstat", -1, -1],
        ["kvaz", -1, -1]
    ];

    /** @const */ var a_4 = [
        ["aliu", -1, -1],
        ["unu", -1, -1]
    ];

    /** @const */ var a_5 = [
        ["aha", -1, -1],
        ["haha", 0, -1],
        ["haleluja", -1, -1],
        ["hola", -1, -1],
        ["hosana", -1, -1],
        ["maltra", -1, -1],
        ["hura", -1, -1],
        ["\u0125a\u0125a", -1, -1],
        ["ekde", -1, -1],
        ["elde", -1, -1],
        ["disde", -1, -1],
        ["ehe", -1, -1],
        ["maltre", -1, -1],
        ["dirlididi", -1, -1],
        ["malpli", -1, -1],
        ["mal\u0109i", -1, -1],
        ["malkaj", -1, -1],
        ["amen", -1, -1],
        ["tamen", 17, -1],
        ["oho", -1, -1],
        ["maltro", -1, -1],
        ["minus", -1, -1],
        ["uhu", -1, -1],
        ["muu", -1, -1]
    ];

    /** @const */ var a_6 = [
        ["tri", -1, -1],
        ["du", -1, -1],
        ["unu", -1, -1]
    ];

    /** @const */ var a_7 = [
        ["dek", -1, -1],
        ["cent", -1, -1]
    ];

    /** @const */ var a_8 = [
        ["k", -1, -1],
        ["kelk", 0, -1],
        ["nen", -1, -1],
        ["t", -1, -1],
        ["mult", 3, -1],
        ["samt", 3, -1],
        ["\u0109", -1, -1]
    ];

    /** @const */ var a_9 = [
        ["a", -1, -1],
        ["e", -1, -1],
        ["i", -1, -1],
        ["j", -1, -1, r_not_after_letter],
        ["aj", 3, -1],
        ["oj", 3, -1],
        ["n", -1, -1, r_not_after_letter],
        ["an", 6, -1],
        ["en", 6, -1],
        ["jn", 6, -1, r_not_after_letter],
        ["ajn", 9, -1],
        ["ojn", 9, -1],
        ["on", 6, -1],
        ["o", -1, -1],
        ["as", -1, -1],
        ["is", -1, -1],
        ["os", -1, -1],
        ["us", -1, -1],
        ["u", -1, -1]
    ];

    /** @const */ var /** Array<int> */ g_vowel = [17, 65, 16];

    /** @const */ var /** Array<int> */ g_aou = [1, 64, 16];

    /** @const */ var /** Array<int> */ g_digit = [255, 3];

    var /** boolean */ B_foreign = false;


    /** @return {boolean} */
    function r_canonical_form() {
        var /** number */ among_var;
        B_foreign = false;
        while(true)
        {
            /** @const */ var /** number */ v_1 = base.cursor;
            lab0: {
                base.bra = base.cursor;
                among_var = base.find_among(a_0);
                base.ket = base.cursor;
                switch (among_var) {
                    case 1:
                        if (!base.slice_from("\u0109"))
                        {
                            return false;
                        }
                        break;
                    case 2:
                        if (!base.slice_from("\u011D"))
                        {
                            return false;
                        }
                        break;
                    case 3:
                        if (!base.slice_from("\u0125"))
                        {
                            return false;
                        }
                        break;
                    case 4:
                        if (!base.slice_from("\u0135"))
                        {
                            return false;
                        }
                        break;
                    case 5:
                        if (!base.slice_from("\u015D"))
                        {
                            return false;
                        }
                        break;
                    case 6:
                        if (!base.slice_from("\u016D"))
                        {
                            return false;
                        }
                        break;
                    case 7:
                        if (!base.slice_from("a"))
                        {
                            return false;
                        }
                        B_foreign = true;
                        break;
                    case 8:
                        if (!base.slice_from("e"))
                        {
                            return false;
                        }
                        B_foreign = true;
                        break;
                    case 9:
                        if (!base.slice_from("i"))
                        {
                            return false;
                        }
                        B_foreign = true;
                        break;
                    case 10:
                        if (!base.slice_from("o"))
                        {
                            return false;
                        }
                        B_foreign = true;
                        break;
                    case 11:
                        if (!base.slice_from("u"))
                        {
                            return false;
                        }
                        B_foreign = true;
                        break;
                    case 12:
                        B_foreign = true;
                        break;
                    case 13:
                        B_foreign = false;
                        break;
                    case 14:
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
        lab1: {
            if (!B_foreign)
            {
                break lab1;
            }
            return false;
        }
        return true;
    };

    /** @return {boolean} */
    function r_initial_apostrophe() {
        base.bra = base.cursor;
        if (!(base.eq_s("'")))
        {
            return false;
        }
        base.ket = base.cursor;
        if (!(base.eq_s("st")))
        {
            return false;
        }
        if (base.find_among(a_1) == 0)
        {
            return false;
        }
        if (base.cursor < base.limit)
        {
            return false;
        }
        if (!base.slice_from("e"))
        {
            return false;
        }
        return true;
    };

    /** @return {boolean} */
    function r_pronoun() {
        base.ket = base.cursor;
        /** @const */ var /** number */ v_1 = base.limit - base.cursor;
        lab0: {
            if (!(base.eq_s_b("n")))
            {
                base.cursor = base.limit - v_1;
                break lab0;
            }
        }
        base.bra = base.cursor;
        if (base.find_among_b(a_2) == 0)
        {
            return false;
        }
        lab1: {
            /** @const */ var /** number */ v_2 = base.limit - base.cursor;
            lab2: {
                if (base.cursor > base.limit_backward)
                {
                    break lab2;
                }
                break lab1;
            }
            base.cursor = base.limit - v_2;
            if (!(base.eq_s_b("-")))
            {
                return false;
            }
        }
        if (!base.slice_del())
        {
            return false;
        }
        return true;
    };

    /** @return {boolean} */
    function r_final_apostrophe() {
        base.ket = base.cursor;
        if (!(base.eq_s_b("'")))
        {
            return false;
        }
        base.bra = base.cursor;
        lab0: {
            /** @const */ var /** number */ v_1 = base.limit - base.cursor;
            lab1: {
                if (!(base.eq_s_b("l")))
                {
                    break lab1;
                }
                if (base.cursor > base.limit_backward)
                {
                    break lab1;
                }
                if (!base.slice_from("a"))
                {
                    return false;
                }
                break lab0;
            }
            base.cursor = base.limit - v_1;
            lab2: {
                if (!(base.eq_s_b("un")))
                {
                    break lab2;
                }
                if (base.cursor > base.limit_backward)
                {
                    break lab2;
                }
                if (!base.slice_from("u"))
                {
                    return false;
                }
                break lab0;
            }
            base.cursor = base.limit - v_1;
            lab3: {
                if (base.find_among_b(a_3) == 0)
                {
                    break lab3;
                }
                lab4: {
                    /** @const */ var /** number */ v_2 = base.limit - base.cursor;
                    lab5: {
                        if (base.cursor > base.limit_backward)
                        {
                            break lab5;
                        }
                        break lab4;
                    }
                    base.cursor = base.limit - v_2;
                    if (!(base.eq_s_b("-")))
                    {
                        break lab3;
                    }
                }
                if (!base.slice_from("a\u016D"))
                {
                    return false;
                }
                break lab0;
            }
            base.cursor = base.limit - v_1;
            if (!base.slice_from("o"))
            {
                return false;
            }
        }
        return true;
    };

    /** @return {boolean} */
    function r_ujn_suffix() {
        base.ket = base.cursor;
        /** @const */ var /** number */ v_1 = base.limit - base.cursor;
        lab0: {
            if (!(base.eq_s_b("n")))
            {
                base.cursor = base.limit - v_1;
                break lab0;
            }
        }
        /** @const */ var /** number */ v_2 = base.limit - base.cursor;
        lab1: {
            if (!(base.eq_s_b("j")))
            {
                base.cursor = base.limit - v_2;
                break lab1;
            }
        }
        base.bra = base.cursor;
        if (base.find_among_b(a_4) == 0)
        {
            return false;
        }
        lab2: {
            /** @const */ var /** number */ v_3 = base.limit - base.cursor;
            lab3: {
                if (base.cursor > base.limit_backward)
                {
                    break lab3;
                }
                break lab2;
            }
            base.cursor = base.limit - v_3;
            if (!(base.eq_s_b("-")))
            {
                return false;
            }
        }
        if (!base.slice_del())
        {
            return false;
        }
        return true;
    };

    /** @return {boolean} */
    function r_uninflected() {
        if (base.find_among_b(a_5) == 0)
        {
            return false;
        }
        lab0: {
            /** @const */ var /** number */ v_1 = base.limit - base.cursor;
            lab1: {
                if (base.cursor > base.limit_backward)
                {
                    break lab1;
                }
                break lab0;
            }
            base.cursor = base.limit - v_1;
            if (!(base.eq_s_b("-")))
            {
                return false;
            }
        }
        return true;
    };

    /** @return {boolean} */
    function r_merged_numeral() {
        if (base.find_among_b(a_6) == 0)
        {
            return false;
        }
        if (base.find_among_b(a_7) == 0)
        {
            return false;
        }
        return true;
    };

    /** @return {boolean} */
    function r_correlative() {
        base.ket = base.cursor;
        base.bra = base.cursor;
        /** @const */ var /** number */ v_1 = base.limit - base.cursor;
        lab0: {
            /** @const */ var /** number */ v_2 = base.limit - base.cursor;
            lab1: {
                /** @const */ var /** number */ v_3 = base.limit - base.cursor;
                lab2: {
                    if (!(base.eq_s_b("n")))
                    {
                        base.cursor = base.limit - v_3;
                        break lab2;
                    }
                }
                base.bra = base.cursor;
                if (!(base.eq_s_b("e")))
                {
                    break lab1;
                }
                break lab0;
            }
            base.cursor = base.limit - v_2;
            /** @const */ var /** number */ v_4 = base.limit - base.cursor;
            lab3: {
                if (!(base.eq_s_b("n")))
                {
                    base.cursor = base.limit - v_4;
                    break lab3;
                }
            }
            /** @const */ var /** number */ v_5 = base.limit - base.cursor;
            lab4: {
                if (!(base.eq_s_b("j")))
                {
                    base.cursor = base.limit - v_5;
                    break lab4;
                }
            }
            base.bra = base.cursor;
            if (!(base.in_grouping_b(g_aou, 97, 117)))
            {
                return false;
            }
        }
        if (!(base.eq_s_b("i")))
        {
            return false;
        }
        /** @const */ var /** number */ v_6 = base.limit - base.cursor;
        lab5: {
            if (base.find_among_b(a_8) == 0)
            {
                base.cursor = base.limit - v_6;
                break lab5;
            }
        }
        lab6: {
            /** @const */ var /** number */ v_7 = base.limit - base.cursor;
            lab7: {
                if (base.cursor > base.limit_backward)
                {
                    break lab7;
                }
                break lab6;
            }
            base.cursor = base.limit - v_7;
            if (!(base.eq_s_b("-")))
            {
                return false;
            }
        }
        base.cursor = base.limit - v_1;
        if (!base.slice_del())
        {
            return false;
        }
        return true;
    };

    /** @return {boolean} */
    function r_long_word() {
        lab0: {
            /** @const */ var /** number */ v_1 = base.limit - base.cursor;
            lab1: {
                for (var /** number */ v_2 = 2; v_2 > 0; v_2--)
                {
                    if (!base.go_out_grouping_b(g_vowel, 97, 117))
                    {
                        break lab1;
                    }
                    base.cursor--;
                }
                break lab0;
            }
            base.cursor = base.limit - v_1;
            lab2: {
                golab3: while(true)
                {
                    lab4: {
                        if (!(base.eq_s_b("-")))
                        {
                            break lab4;
                        }
                        break golab3;
                    }
                    if (base.cursor <= base.limit_backward)
                    {
                        break lab2;
                    }
                    base.cursor--;
                }
                if (base.cursor <= base.limit_backward)
                {
                    break lab2;
                }
                base.cursor--;
                break lab0;
            }
            base.cursor = base.limit - v_1;
            if (!base.go_out_grouping_b(g_digit, 48, 57))
            {
                return false;
            }
            base.cursor--;
        }
        return true;
    };

    /** @return {boolean} */
    function r_not_after_letter() {
        lab0: {
            /** @const */ var /** number */ v_1 = base.limit - base.cursor;
            lab1: {
                if (!(base.eq_s_b("-")))
                {
                    break lab1;
                }
                break lab0;
            }
            base.cursor = base.limit - v_1;
            if (!(base.in_grouping_b(g_digit, 48, 57)))
            {
                return false;
            }
        }
        return true;
    };

    /** @return {boolean} */
    function r_standard_suffix() {
        base.ket = base.cursor;
        if (base.find_among_b(a_9) == 0)
        {
            return false;
        }
        /** @const */ var /** number */ v_1 = base.limit - base.cursor;
        lab0: {
            if (!(base.eq_s_b("-")))
            {
                base.cursor = base.limit - v_1;
                break lab0;
            }
        }
        base.bra = base.cursor;
        if (!base.slice_del())
        {
            return false;
        }
        return true;
    };

    this.stem = /** @return {boolean} */ function() {
        /** @const */ var /** number */ v_1 = base.cursor;
        if (!r_canonical_form())
        {
            return false;
        }
        base.cursor = v_1;
        /** @const */ var /** number */ v_2 = base.cursor;
        r_initial_apostrophe();
        base.cursor = v_2;
        base.limit_backward = base.cursor; base.cursor = base.limit;
        {
            /** @const */ var /** number */ v_3 = base.limit - base.cursor;
            lab0: {
                if (!r_pronoun())
                {
                    break lab0;
                }
                return false;
            }
            base.cursor = base.limit - v_3;
        }
        /** @const */ var /** number */ v_4 = base.limit - base.cursor;
        r_final_apostrophe();
        base.cursor = base.limit - v_4;
        {
            /** @const */ var /** number */ v_5 = base.limit - base.cursor;
            lab1: {
                if (!r_correlative())
                {
                    break lab1;
                }
                return false;
            }
            base.cursor = base.limit - v_5;
        }
        {
            /** @const */ var /** number */ v_6 = base.limit - base.cursor;
            lab2: {
                if (!r_uninflected())
                {
                    break lab2;
                }
                return false;
            }
            base.cursor = base.limit - v_6;
        }
        {
            /** @const */ var /** number */ v_7 = base.limit - base.cursor;
            lab3: {
                if (!r_merged_numeral())
                {
                    break lab3;
                }
                return false;
            }
            base.cursor = base.limit - v_7;
        }
        {
            /** @const */ var /** number */ v_8 = base.limit - base.cursor;
            lab4: {
                if (!r_ujn_suffix())
                {
                    break lab4;
                }
                return false;
            }
            base.cursor = base.limit - v_8;
        }
        /** @const */ var /** number */ v_9 = base.limit - base.cursor;
        if (!r_long_word())
        {
            return false;
        }
        base.cursor = base.limit - v_9;
        if (!r_standard_suffix())
        {
            return false;
        }
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
