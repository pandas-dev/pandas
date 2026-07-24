// Generated from indonesian.sbl by Snowball 3.0.1 - https://snowballstem.org/

/**@constructor*/
var IndonesianStemmer = function() {
    var base = new BaseStemmer();

    /** @const */ var a_0 = [
        ["kah", -1, 1],
        ["lah", -1, 1],
        ["pun", -1, 1]
    ];

    /** @const */ var a_1 = [
        ["nya", -1, 1],
        ["ku", -1, 1],
        ["mu", -1, 1]
    ];

    /** @const */ var a_2 = [
        ["i", -1, 1, r_SUFFIX_I_OK],
        ["an", -1, 1, r_SUFFIX_AN_OK],
        ["kan", 1, 1, r_SUFFIX_KAN_OK]
    ];

    /** @const */ var a_3 = [
        ["di", -1, 1],
        ["ke", -1, 2],
        ["me", -1, 1],
        ["mem", 2, 5],
        ["men", 2, 1],
        ["meng", 4, 1],
        ["meny", 4, 3, r_VOWEL],
        ["pem", -1, 6],
        ["pen", -1, 2],
        ["peng", 8, 2],
        ["peny", 8, 4, r_VOWEL],
        ["ter", -1, 1]
    ];

    /** @const */ var a_4 = [
        ["be", -1, 3, r_KER],
        ["belajar", 0, 4],
        ["ber", 0, 3],
        ["pe", -1, 1],
        ["pelajar", 3, 2],
        ["per", 3, 1]
    ];

    /** @const */ var /** Array<int> */ g_vowel = [17, 65, 16];

    var /** number */ I_prefix = 0;
    var /** number */ I_measure = 0;


    /** @return {boolean} */
    function r_remove_particle() {
        base.ket = base.cursor;
        if (base.find_among_b(a_0) == 0)
        {
            return false;
        }
        base.bra = base.cursor;
        if (!base.slice_del())
        {
            return false;
        }
        I_measure -= 1;
        return true;
    };

    /** @return {boolean} */
    function r_remove_possessive_pronoun() {
        base.ket = base.cursor;
        if (base.find_among_b(a_1) == 0)
        {
            return false;
        }
        base.bra = base.cursor;
        if (!base.slice_del())
        {
            return false;
        }
        I_measure -= 1;
        return true;
    };

    /** @return {boolean} */
    function r_SUFFIX_KAN_OK() {
        if (I_prefix == 3)
        {
            return false;
        }
        if (I_prefix == 2)
        {
            return false;
        }
        return true;
    };

    /** @return {boolean} */
    function r_SUFFIX_AN_OK() {
        return I_prefix != 1;
    };

    /** @return {boolean} */
    function r_SUFFIX_I_OK() {
        if (I_prefix > 2)
        {
            return false;
        }
        {
            /** @const */ var /** number */ v_1 = base.limit - base.cursor;
            lab0: {
                if (!(base.eq_s_b("s")))
                {
                    break lab0;
                }
                return false;
            }
            base.cursor = base.limit - v_1;
        }
        return true;
    };

    /** @return {boolean} */
    function r_remove_suffix() {
        base.ket = base.cursor;
        if (base.find_among_b(a_2) == 0)
        {
            return false;
        }
        base.bra = base.cursor;
        if (!base.slice_del())
        {
            return false;
        }
        I_measure -= 1;
        return true;
    };

    /** @return {boolean} */
    function r_VOWEL() {
        if (!(base.in_grouping(g_vowel, 97, 117)))
        {
            return false;
        }
        return true;
    };

    /** @return {boolean} */
    function r_KER() {
        if (!(base.out_grouping(g_vowel, 97, 117)))
        {
            return false;
        }
        if (!(base.eq_s("er")))
        {
            return false;
        }
        return true;
    };

    /** @return {boolean} */
    function r_remove_first_order_prefix() {
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
                if (!base.slice_del())
                {
                    return false;
                }
                I_prefix = 1;
                I_measure -= 1;
                break;
            case 2:
                if (!base.slice_del())
                {
                    return false;
                }
                I_prefix = 3;
                I_measure -= 1;
                break;
            case 3:
                I_prefix = 1;
                if (!base.slice_from("s"))
                {
                    return false;
                }
                I_measure -= 1;
                break;
            case 4:
                I_prefix = 3;
                if (!base.slice_from("s"))
                {
                    return false;
                }
                I_measure -= 1;
                break;
            case 5:
                I_prefix = 1;
                I_measure -= 1;
                lab0: {
                    /** @const */ var /** number */ v_1 = base.cursor;
                    lab1: {
                        /** @const */ var /** number */ v_2 = base.cursor;
                        if (!(base.in_grouping(g_vowel, 97, 117)))
                        {
                            break lab1;
                        }
                        base.cursor = v_2;
                        if (!base.slice_from("p"))
                        {
                            return false;
                        }
                        break lab0;
                    }
                    base.cursor = v_1;
                    if (!base.slice_del())
                    {
                        return false;
                    }
                }
                break;
            case 6:
                I_prefix = 3;
                I_measure -= 1;
                lab2: {
                    /** @const */ var /** number */ v_3 = base.cursor;
                    lab3: {
                        /** @const */ var /** number */ v_4 = base.cursor;
                        if (!(base.in_grouping(g_vowel, 97, 117)))
                        {
                            break lab3;
                        }
                        base.cursor = v_4;
                        if (!base.slice_from("p"))
                        {
                            return false;
                        }
                        break lab2;
                    }
                    base.cursor = v_3;
                    if (!base.slice_del())
                    {
                        return false;
                    }
                }
                break;
        }
        return true;
    };

    /** @return {boolean} */
    function r_remove_second_order_prefix() {
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
                if (!base.slice_del())
                {
                    return false;
                }
                I_prefix = 2;
                I_measure -= 1;
                break;
            case 2:
                if (!base.slice_from("ajar"))
                {
                    return false;
                }
                I_measure -= 1;
                break;
            case 3:
                if (!base.slice_del())
                {
                    return false;
                }
                I_prefix = 4;
                I_measure -= 1;
                break;
            case 4:
                if (!base.slice_from("ajar"))
                {
                    return false;
                }
                I_prefix = 4;
                I_measure -= 1;
                break;
        }
        return true;
    };

    this.stem = /** @return {boolean} */ function() {
        I_measure = 0;
        /** @const */ var /** number */ v_1 = base.cursor;
        lab0: {
            while(true)
            {
                lab1: {
                    if (!base.go_out_grouping(g_vowel, 97, 117))
                    {
                        break lab1;
                    }
                    base.cursor++;
                    I_measure += 1;
                    continue;
                }
                break;
            }
        }
        base.cursor = v_1;
        if (I_measure <= 2)
        {
            return false;
        }
        I_prefix = 0;
        base.limit_backward = base.cursor; base.cursor = base.limit;
        /** @const */ var /** number */ v_2 = base.limit - base.cursor;
        r_remove_particle();
        base.cursor = base.limit - v_2;
        if (I_measure <= 2)
        {
            return false;
        }
        /** @const */ var /** number */ v_3 = base.limit - base.cursor;
        r_remove_possessive_pronoun();
        base.cursor = base.limit - v_3;
        base.cursor = base.limit_backward;
        if (I_measure <= 2)
        {
            return false;
        }
        lab2: {
            /** @const */ var /** number */ v_4 = base.cursor;
            lab3: {
                /** @const */ var /** number */ v_5 = base.cursor;
                if (!r_remove_first_order_prefix())
                {
                    break lab3;
                }
                /** @const */ var /** number */ v_6 = base.cursor;
                lab4: {
                    /** @const */ var /** number */ v_7 = base.cursor;
                    if (I_measure <= 2)
                    {
                        break lab4;
                    }
                    base.limit_backward = base.cursor; base.cursor = base.limit;
                    if (!r_remove_suffix())
                    {
                        break lab4;
                    }
                    base.cursor = base.limit_backward;
                    base.cursor = v_7;
                    if (I_measure <= 2)
                    {
                        break lab4;
                    }
                    if (!r_remove_second_order_prefix())
                    {
                        break lab4;
                    }
                }
                base.cursor = v_6;
                base.cursor = v_5;
                break lab2;
            }
            base.cursor = v_4;
            /** @const */ var /** number */ v_8 = base.cursor;
            r_remove_second_order_prefix();
            base.cursor = v_8;
            /** @const */ var /** number */ v_9 = base.cursor;
            lab5: {
                if (I_measure <= 2)
                {
                    break lab5;
                }
                base.limit_backward = base.cursor; base.cursor = base.limit;
                if (!r_remove_suffix())
                {
                    break lab5;
                }
                base.cursor = base.limit_backward;
            }
            base.cursor = v_9;
        }
        return true;
    };

    /**@return{string}*/
    this['stemWord'] = function(/**string*/word) {
        base.setCurrent(word);
        this.stem();
        return base.getCurrent();
    };
};
