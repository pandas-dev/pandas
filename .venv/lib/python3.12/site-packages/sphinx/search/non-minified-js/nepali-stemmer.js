// Generated from nepali.sbl by Snowball 3.0.1 - https://snowballstem.org/

/**@constructor*/
var NepaliStemmer = function() {
    var base = new BaseStemmer();

    /** @const */ var a_0 = [
        ["\u0932\u093E\u0907", -1, 1],
        ["\u0932\u093E\u0908", -1, 1],
        ["\u0938\u0901\u0917", -1, 1],
        ["\u0938\u0902\u0917", -1, 1],
        ["\u092E\u093E\u0930\u094D\u092B\u0924", -1, 1],
        ["\u0930\u0924", -1, 1],
        ["\u0915\u093E", -1, 2],
        ["\u092E\u093E", -1, 1],
        ["\u0926\u094D\u0935\u093E\u0930\u093E", -1, 1],
        ["\u0915\u093F", -1, 2],
        ["\u092A\u091B\u093F", -1, 1],
        ["\u0915\u0940", -1, 2],
        ["\u0932\u0947", -1, 1],
        ["\u0915\u0948", -1, 2],
        ["\u0938\u0901\u0917\u0948", -1, 1],
        ["\u092E\u0948", -1, 1],
        ["\u0915\u094B", -1, 2]
    ];

    /** @const */ var a_1 = [
        ["\u0901", -1, 1],
        ["\u0902", -1, 1],
        ["\u0948", -1, 2]
    ];

    /** @const */ var a_2 = [
        ["\u0925\u093F\u090F", -1, 1],
        ["\u091B", -1, 1],
        ["\u0907\u091B", 1, 1],
        ["\u090F\u091B", 1, 1],
        ["\u093F\u091B", 1, 1],
        ["\u0947\u091B", 1, 1],
        ["\u0928\u0947\u091B", 5, 1],
        ["\u0939\u0941\u0928\u0947\u091B", 6, 1],
        ["\u0907\u0928\u094D\u091B", 1, 1],
        ["\u093F\u0928\u094D\u091B", 1, 1],
        ["\u0939\u0941\u0928\u094D\u091B", 1, 1],
        ["\u090F\u0915\u093E", -1, 1],
        ["\u0907\u090F\u0915\u093E", 11, 1],
        ["\u093F\u090F\u0915\u093E", 11, 1],
        ["\u0947\u0915\u093E", -1, 1],
        ["\u0928\u0947\u0915\u093E", 14, 1],
        ["\u0926\u093E", -1, 1],
        ["\u0907\u0926\u093E", 16, 1],
        ["\u093F\u0926\u093E", 16, 1],
        ["\u0926\u0947\u0916\u093F", -1, 1],
        ["\u092E\u093E\u0925\u093F", -1, 1],
        ["\u090F\u0915\u0940", -1, 1],
        ["\u0907\u090F\u0915\u0940", 21, 1],
        ["\u093F\u090F\u0915\u0940", 21, 1],
        ["\u0947\u0915\u0940", -1, 1],
        ["\u0926\u0947\u0916\u0940", -1, 1],
        ["\u0925\u0940", -1, 1],
        ["\u0926\u0940", -1, 1],
        ["\u091B\u0941", -1, 1],
        ["\u090F\u091B\u0941", 28, 1],
        ["\u0947\u091B\u0941", 28, 1],
        ["\u0928\u0947\u091B\u0941", 30, 1],
        ["\u0928\u0941", -1, 1],
        ["\u0939\u0930\u0941", -1, 1],
        ["\u0939\u0930\u0942", -1, 1],
        ["\u091B\u0947", -1, 1],
        ["\u0925\u0947", -1, 1],
        ["\u0928\u0947", -1, 1],
        ["\u090F\u0915\u0948", -1, 1],
        ["\u0947\u0915\u0948", -1, 1],
        ["\u0928\u0947\u0915\u0948", 39, 1],
        ["\u0926\u0948", -1, 1],
        ["\u0907\u0926\u0948", 41, 1],
        ["\u093F\u0926\u0948", 41, 1],
        ["\u090F\u0915\u094B", -1, 1],
        ["\u0907\u090F\u0915\u094B", 44, 1],
        ["\u093F\u090F\u0915\u094B", 44, 1],
        ["\u0947\u0915\u094B", -1, 1],
        ["\u0928\u0947\u0915\u094B", 47, 1],
        ["\u0926\u094B", -1, 1],
        ["\u0907\u0926\u094B", 49, 1],
        ["\u093F\u0926\u094B", 49, 1],
        ["\u092F\u094B", -1, 1],
        ["\u0907\u092F\u094B", 52, 1],
        ["\u092D\u092F\u094B", 52, 1],
        ["\u093F\u092F\u094B", 52, 1],
        ["\u0925\u093F\u092F\u094B", 55, 1],
        ["\u0926\u093F\u092F\u094B", 55, 1],
        ["\u0925\u094D\u092F\u094B", 52, 1],
        ["\u091B\u094C", -1, 1],
        ["\u0907\u091B\u094C", 59, 1],
        ["\u090F\u091B\u094C", 59, 1],
        ["\u093F\u091B\u094C", 59, 1],
        ["\u0947\u091B\u094C", 59, 1],
        ["\u0928\u0947\u091B\u094C", 63, 1],
        ["\u092F\u094C", -1, 1],
        ["\u0925\u093F\u092F\u094C", 65, 1],
        ["\u091B\u094D\u092F\u094C", 65, 1],
        ["\u0925\u094D\u092F\u094C", 65, 1],
        ["\u091B\u0928\u094D", -1, 1],
        ["\u0907\u091B\u0928\u094D", 69, 1],
        ["\u090F\u091B\u0928\u094D", 69, 1],
        ["\u093F\u091B\u0928\u094D", 69, 1],
        ["\u0947\u091B\u0928\u094D", 69, 1],
        ["\u0928\u0947\u091B\u0928\u094D", 73, 1],
        ["\u0932\u093E\u0928\u094D", -1, 1],
        ["\u091B\u093F\u0928\u094D", -1, 1],
        ["\u0925\u093F\u0928\u094D", -1, 1],
        ["\u092A\u0930\u094D", -1, 1],
        ["\u0907\u0938\u094D", -1, 1],
        ["\u0925\u093F\u0907\u0938\u094D", 79, 1],
        ["\u091B\u0938\u094D", -1, 1],
        ["\u0907\u091B\u0938\u094D", 81, 1],
        ["\u090F\u091B\u0938\u094D", 81, 1],
        ["\u093F\u091B\u0938\u094D", 81, 1],
        ["\u0947\u091B\u0938\u094D", 81, 1],
        ["\u0928\u0947\u091B\u0938\u094D", 85, 1],
        ["\u093F\u0938\u094D", -1, 1],
        ["\u0925\u093F\u0938\u094D", 87, 1],
        ["\u091B\u0947\u0938\u094D", -1, 1],
        ["\u0939\u094B\u0938\u094D", -1, 1]
    ];


    /** @return {boolean} */
    function r_remove_category_1() {
        var /** number */ among_var;
        base.ket = base.cursor;
        among_var = base.find_among_b(a_0);
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
                lab0: {
                    /** @const */ var /** number */ v_1 = base.limit - base.cursor;
                    lab1: {
                        if (!(base.eq_s_b("\u090F")))
                        {
                            break lab1;
                        }
                        break lab0;
                    }
                    base.cursor = base.limit - v_1;
                    lab2: {
                        if (!(base.eq_s_b("\u0947")))
                        {
                            break lab2;
                        }
                        break lab0;
                    }
                    base.cursor = base.limit - v_1;
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
    function r_remove_category_2() {
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
                lab0: {
                    /** @const */ var /** number */ v_1 = base.limit - base.cursor;
                    lab1: {
                        if (!(base.eq_s_b("\u092F\u094C")))
                        {
                            break lab1;
                        }
                        break lab0;
                    }
                    base.cursor = base.limit - v_1;
                    lab2: {
                        if (!(base.eq_s_b("\u091B\u094C")))
                        {
                            break lab2;
                        }
                        break lab0;
                    }
                    base.cursor = base.limit - v_1;
                    lab3: {
                        if (!(base.eq_s_b("\u0928\u094C")))
                        {
                            break lab3;
                        }
                        break lab0;
                    }
                    base.cursor = base.limit - v_1;
                    if (!(base.eq_s_b("\u0925\u0947")))
                    {
                        return false;
                    }
                }
                if (!base.slice_del())
                {
                    return false;
                }
                break;
            case 2:
                if (!(base.eq_s_b("\u0924\u094D\u0930")))
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
    function r_remove_category_3() {
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
        return true;
    };

    this.stem = /** @return {boolean} */ function() {
        base.limit_backward = base.cursor; base.cursor = base.limit;
        /** @const */ var /** number */ v_1 = base.limit - base.cursor;
        r_remove_category_1();
        base.cursor = base.limit - v_1;
        while(true)
        {
            /** @const */ var /** number */ v_2 = base.limit - base.cursor;
            lab0: {
                /** @const */ var /** number */ v_3 = base.limit - base.cursor;
                r_remove_category_2();
                base.cursor = base.limit - v_3;
                if (!r_remove_category_3())
                {
                    break lab0;
                }
                continue;
            }
            base.cursor = base.limit - v_2;
            break;
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
