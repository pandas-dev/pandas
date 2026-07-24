// Generated from tamil.sbl by Snowball 3.0.1 - https://snowballstem.org/

/**@constructor*/
var TamilStemmer = function() {
    var base = new BaseStemmer();

    /** @const */ var a_0 = [
        ["\u0BB5\u0BC1", -1, 3],
        ["\u0BB5\u0BC2", -1, 4],
        ["\u0BB5\u0BCA", -1, 2],
        ["\u0BB5\u0BCB", -1, 1]
    ];

    /** @const */ var a_1 = [
        ["\u0B95", -1, -1],
        ["\u0B99", -1, -1],
        ["\u0B9A", -1, -1],
        ["\u0B9E", -1, -1],
        ["\u0BA4", -1, -1],
        ["\u0BA8", -1, -1],
        ["\u0BAA", -1, -1],
        ["\u0BAE", -1, -1],
        ["\u0BAF", -1, -1],
        ["\u0BB5", -1, -1]
    ];

    /** @const */ var a_2 = [
        ["\u0BBF", -1, -1],
        ["\u0BC0", -1, -1],
        ["\u0BC8", -1, -1]
    ];

    /** @const */ var a_3 = [
        ["\u0BBE", -1, -1],
        ["\u0BBF", -1, -1],
        ["\u0BC0", -1, -1],
        ["\u0BC1", -1, -1],
        ["\u0BC2", -1, -1],
        ["\u0BC6", -1, -1],
        ["\u0BC7", -1, -1],
        ["\u0BC8", -1, -1]
    ];

    /** @const */ var a_4 = [
        ["", -1, 2],
        ["\u0BC8", 0, 1],
        ["\u0BCD", 0, 1]
    ];

    /** @const */ var a_5 = [
        ["\u0BA8\u0BCD\u0BA4", -1, 1],
        ["\u0BAF", -1, 1],
        ["\u0BB5", -1, 1],
        ["\u0BA9\u0BC1", -1, 8],
        ["\u0BC1\u0B95\u0BCD", -1, 7],
        ["\u0BC1\u0B95\u0BCD\u0B95\u0BCD", -1, 7],
        ["\u0B9F\u0BCD\u0B95\u0BCD", -1, 3],
        ["\u0BB1\u0BCD\u0B95\u0BCD", -1, 4],
        ["\u0B99\u0BCD", -1, 9],
        ["\u0B9F\u0BCD\u0B9F\u0BCD", -1, 5],
        ["\u0BA4\u0BCD\u0BA4\u0BCD", -1, 6],
        ["\u0BA8\u0BCD\u0BA4\u0BCD", -1, 1],
        ["\u0BA8\u0BCD", -1, 1],
        ["\u0B9F\u0BCD\u0BAA\u0BCD", -1, 3],
        ["\u0BAF\u0BCD", -1, 2],
        ["\u0BA9\u0BCD\u0BB1\u0BCD", -1, 4],
        ["\u0BB5\u0BCD", -1, 1]
    ];

    /** @const */ var a_6 = [
        ["\u0B95", -1, -1],
        ["\u0B9A", -1, -1],
        ["\u0B9F", -1, -1],
        ["\u0BA4", -1, -1],
        ["\u0BAA", -1, -1],
        ["\u0BB1", -1, -1]
    ];

    /** @const */ var a_7 = [
        ["\u0B95", -1, -1],
        ["\u0B9A", -1, -1],
        ["\u0B9F", -1, -1],
        ["\u0BA4", -1, -1],
        ["\u0BAA", -1, -1],
        ["\u0BB1", -1, -1]
    ];

    /** @const */ var a_8 = [
        ["\u0B9E", -1, -1],
        ["\u0BA3", -1, -1],
        ["\u0BA8", -1, -1],
        ["\u0BA9", -1, -1],
        ["\u0BAE", -1, -1],
        ["\u0BAF", -1, -1],
        ["\u0BB0", -1, -1],
        ["\u0BB2", -1, -1],
        ["\u0BB3", -1, -1],
        ["\u0BB4", -1, -1],
        ["\u0BB5", -1, -1]
    ];

    /** @const */ var a_9 = [
        ["\u0BBE", -1, -1],
        ["\u0BBF", -1, -1],
        ["\u0BC0", -1, -1],
        ["\u0BC1", -1, -1],
        ["\u0BC2", -1, -1],
        ["\u0BC6", -1, -1],
        ["\u0BC7", -1, -1],
        ["\u0BC8", -1, -1],
        ["\u0BCD", -1, -1]
    ];

    /** @const */ var a_10 = [
        ["\u0B85", -1, -1],
        ["\u0B87", -1, -1],
        ["\u0B89", -1, -1]
    ];

    /** @const */ var a_11 = [
        ["\u0B95", -1, -1],
        ["\u0B99", -1, -1],
        ["\u0B9A", -1, -1],
        ["\u0B9E", -1, -1],
        ["\u0BA4", -1, -1],
        ["\u0BA8", -1, -1],
        ["\u0BAA", -1, -1],
        ["\u0BAE", -1, -1],
        ["\u0BAF", -1, -1],
        ["\u0BB5", -1, -1]
    ];

    /** @const */ var a_12 = [
        ["\u0B95", -1, -1],
        ["\u0B9A", -1, -1],
        ["\u0B9F", -1, -1],
        ["\u0BA4", -1, -1],
        ["\u0BAA", -1, -1],
        ["\u0BB1", -1, -1]
    ];

    /** @const */ var a_13 = [
        ["\u0B95\u0BB3\u0BCD", -1, 4],
        ["\u0BC1\u0B99\u0BCD\u0B95\u0BB3\u0BCD", 0, 1],
        ["\u0B9F\u0BCD\u0B95\u0BB3\u0BCD", 0, 3],
        ["\u0BB1\u0BCD\u0B95\u0BB3\u0BCD", 0, 2]
    ];

    /** @const */ var a_14 = [
        ["\u0BBE", -1, -1],
        ["\u0BC7", -1, -1],
        ["\u0BCB", -1, -1]
    ];

    /** @const */ var a_15 = [
        ["\u0BAA\u0BBF", -1, -1],
        ["\u0BB5\u0BBF", -1, -1]
    ];

    /** @const */ var a_16 = [
        ["\u0BBE", -1, -1],
        ["\u0BBF", -1, -1],
        ["\u0BC0", -1, -1],
        ["\u0BC1", -1, -1],
        ["\u0BC2", -1, -1],
        ["\u0BC6", -1, -1],
        ["\u0BC7", -1, -1],
        ["\u0BC8", -1, -1]
    ];

    /** @const */ var a_17 = [
        ["\u0BAA\u0B9F\u0BCD\u0B9F", -1, 3],
        ["\u0BAA\u0B9F\u0BCD\u0B9F\u0BA3", -1, 3],
        ["\u0BA4\u0BBE\u0BA9", -1, 3],
        ["\u0BAA\u0B9F\u0BBF\u0BA4\u0BBE\u0BA9", 2, 3],
        ["\u0BC6\u0BA9", -1, 1],
        ["\u0BBE\u0B95\u0BBF\u0BAF", -1, 1],
        ["\u0B95\u0BC1\u0BB0\u0BBF\u0BAF", -1, 3],
        ["\u0BC1\u0B9F\u0BC8\u0BAF", -1, 1],
        ["\u0BB2\u0BCD\u0BB2", -1, 2],
        ["\u0BC1\u0BB3\u0BCD\u0BB3", -1, 1],
        ["\u0BBE\u0B95\u0BBF", -1, 1],
        ["\u0BAA\u0B9F\u0BBF", -1, 3],
        ["\u0BBF\u0BA9\u0BCD\u0BB1\u0BBF", -1, 1],
        ["\u0BAA\u0BB1\u0BCD\u0BB1\u0BBF", -1, 3],
        ["\u0BAA\u0B9F\u0BC1", -1, 3],
        ["\u0BB5\u0BBF\u0B9F\u0BC1", -1, 3],
        ["\u0BAA\u0B9F\u0BCD\u0B9F\u0BC1", -1, 3],
        ["\u0BB5\u0BBF\u0B9F\u0BCD\u0B9F\u0BC1", -1, 3],
        ["\u0BAA\u0B9F\u0BCD\u0B9F\u0BA4\u0BC1", -1, 3],
        ["\u0BC6\u0BA9\u0BCD\u0BB1\u0BC1", -1, 1],
        ["\u0BC1\u0B9F\u0BC8", -1, 1],
        ["\u0BBF\u0BB2\u0BCD\u0BB2\u0BC8", -1, 1],
        ["\u0BC1\u0B9F\u0BA9\u0BCD", -1, 1],
        ["\u0BBF\u0B9F\u0BAE\u0BCD", -1, 1],
        ["\u0BC6\u0BB2\u0BCD\u0BB2\u0BBE\u0BAE\u0BCD", -1, 3],
        ["\u0BC6\u0BA9\u0BC1\u0BAE\u0BCD", -1, 1]
    ];

    /** @const */ var a_18 = [
        ["\u0BBE", -1, -1],
        ["\u0BBF", -1, -1],
        ["\u0BC0", -1, -1],
        ["\u0BC1", -1, -1],
        ["\u0BC2", -1, -1],
        ["\u0BC6", -1, -1],
        ["\u0BC7", -1, -1],
        ["\u0BC8", -1, -1]
    ];

    /** @const */ var a_19 = [
        ["\u0BBE", -1, -1],
        ["\u0BBF", -1, -1],
        ["\u0BC0", -1, -1],
        ["\u0BC1", -1, -1],
        ["\u0BC2", -1, -1],
        ["\u0BC6", -1, -1],
        ["\u0BC7", -1, -1],
        ["\u0BC8", -1, -1]
    ];

    /** @const */ var a_20 = [
        ["\u0BB5\u0BBF\u0B9F", -1, 2],
        ["\u0BC0", -1, 7],
        ["\u0BCA\u0B9F\u0BC1", -1, 2],
        ["\u0BCB\u0B9F\u0BC1", -1, 2],
        ["\u0BA4\u0BC1", -1, 6],
        ["\u0BBF\u0BB0\u0BC1\u0BA8\u0BCD\u0BA4\u0BC1", 4, 2],
        ["\u0BBF\u0BA9\u0BCD\u0BB1\u0BC1", -1, 2],
        ["\u0BC1\u0B9F\u0BC8", -1, 2],
        ["\u0BA9\u0BC8", -1, 1],
        ["\u0B95\u0BA3\u0BCD", -1, 1],
        ["\u0BBF\u0BA9\u0BCD", -1, 3],
        ["\u0BAE\u0BC1\u0BA9\u0BCD", -1, 1],
        ["\u0BBF\u0B9F\u0BAE\u0BCD", -1, 4],
        ["\u0BBF\u0BB1\u0BCD", -1, 2],
        ["\u0BAE\u0BC7\u0BB1\u0BCD", -1, 1],
        ["\u0BB2\u0BCD", -1, 5],
        ["\u0BBE\u0BAE\u0BB2\u0BCD", 15, 2],
        ["\u0BBE\u0BB2\u0BCD", 15, 2],
        ["\u0BBF\u0BB2\u0BCD", 15, 2],
        ["\u0BAE\u0BC7\u0BB2\u0BCD", 15, 1],
        ["\u0BC1\u0BB3\u0BCD", -1, 2],
        ["\u0B95\u0BC0\u0BB4\u0BCD", -1, 1]
    ];

    /** @const */ var a_21 = [
        ["\u0B95", -1, -1],
        ["\u0B9A", -1, -1],
        ["\u0B9F", -1, -1],
        ["\u0BA4", -1, -1],
        ["\u0BAA", -1, -1],
        ["\u0BB1", -1, -1]
    ];

    /** @const */ var a_22 = [
        ["\u0B95", -1, -1],
        ["\u0B9A", -1, -1],
        ["\u0B9F", -1, -1],
        ["\u0BA4", -1, -1],
        ["\u0BAA", -1, -1],
        ["\u0BB1", -1, -1]
    ];

    /** @const */ var a_23 = [
        ["\u0B85", -1, -1],
        ["\u0B86", -1, -1],
        ["\u0B87", -1, -1],
        ["\u0B88", -1, -1],
        ["\u0B89", -1, -1],
        ["\u0B8A", -1, -1],
        ["\u0B8E", -1, -1],
        ["\u0B8F", -1, -1],
        ["\u0B90", -1, -1],
        ["\u0B92", -1, -1],
        ["\u0B93", -1, -1],
        ["\u0B94", -1, -1]
    ];

    /** @const */ var a_24 = [
        ["\u0BBE", -1, -1],
        ["\u0BBF", -1, -1],
        ["\u0BC0", -1, -1],
        ["\u0BC1", -1, -1],
        ["\u0BC2", -1, -1],
        ["\u0BC6", -1, -1],
        ["\u0BC7", -1, -1],
        ["\u0BC8", -1, -1]
    ];

    /** @const */ var a_25 = [
        ["\u0B95", -1, 1],
        ["\u0BA4", -1, 1],
        ["\u0BA9", -1, 1],
        ["\u0BAA", -1, 1],
        ["\u0BAF", -1, 1],
        ["\u0BBE", -1, 5],
        ["\u0B95\u0BC1", -1, 6],
        ["\u0BAA\u0B9F\u0BC1", -1, 1],
        ["\u0BA4\u0BC1", -1, 3],
        ["\u0BBF\u0BB1\u0BCD\u0BB1\u0BC1", -1, 1],
        ["\u0BA9\u0BC8", -1, 1],
        ["\u0BB5\u0BC8", -1, 1],
        ["\u0BA9\u0BA9\u0BCD", -1, 1],
        ["\u0BAA\u0BA9\u0BCD", -1, 1],
        ["\u0BB5\u0BA9\u0BCD", -1, 2],
        ["\u0BBE\u0BA9\u0BCD", -1, 4],
        ["\u0BA9\u0BBE\u0BA9\u0BCD", 15, 1],
        ["\u0BAE\u0BBF\u0BA9\u0BCD", -1, 1],
        ["\u0BA9\u0BC6\u0BA9\u0BCD", -1, 1],
        ["\u0BC7\u0BA9\u0BCD", -1, 5],
        ["\u0BA9\u0BAE\u0BCD", -1, 1],
        ["\u0BAA\u0BAE\u0BCD", -1, 1],
        ["\u0BBE\u0BAE\u0BCD", -1, 5],
        ["\u0B95\u0BC1\u0BAE\u0BCD", -1, 1],
        ["\u0B9F\u0BC1\u0BAE\u0BCD", -1, 5],
        ["\u0BA4\u0BC1\u0BAE\u0BCD", -1, 1],
        ["\u0BB1\u0BC1\u0BAE\u0BCD", -1, 1],
        ["\u0BC6\u0BAE\u0BCD", -1, 5],
        ["\u0BC7\u0BAE\u0BCD", -1, 5],
        ["\u0BCB\u0BAE\u0BCD", -1, 5],
        ["\u0BBE\u0BAF\u0BCD", -1, 5],
        ["\u0BA9\u0BB0\u0BCD", -1, 1],
        ["\u0BAA\u0BB0\u0BCD", -1, 1],
        ["\u0BC0\u0BAF\u0BB0\u0BCD", -1, 5],
        ["\u0BB5\u0BB0\u0BCD", -1, 1],
        ["\u0BBE\u0BB0\u0BCD", -1, 5],
        ["\u0BA9\u0BBE\u0BB0\u0BCD", 35, 1],
        ["\u0BAE\u0BBE\u0BB0\u0BCD", 35, 1],
        ["\u0B95\u0BCA\u0BA3\u0BCD\u0B9F\u0BBF\u0BB0\u0BCD", -1, 1],
        ["\u0BA9\u0BBF\u0BB0\u0BCD", -1, 5],
        ["\u0BC0\u0BB0\u0BCD", -1, 5],
        ["\u0BA9\u0BB3\u0BCD", -1, 1],
        ["\u0BAA\u0BB3\u0BCD", -1, 1],
        ["\u0BB5\u0BB3\u0BCD", -1, 1],
        ["\u0BBE\u0BB3\u0BCD", -1, 5],
        ["\u0BA9\u0BBE\u0BB3\u0BCD", 44, 1]
    ];

    /** @const */ var a_26 = [
        ["\u0B95\u0BBF\u0BB1", -1, -1],
        ["\u0B95\u0BBF\u0BA9\u0BCD\u0BB1", -1, -1],
        ["\u0BBE\u0BA8\u0BBF\u0BA9\u0BCD\u0BB1", -1, -1],
        ["\u0B95\u0BBF\u0BB1\u0BCD", -1, -1],
        ["\u0B95\u0BBF\u0BA9\u0BCD\u0BB1\u0BCD", -1, -1],
        ["\u0BBE\u0BA8\u0BBF\u0BA9\u0BCD\u0BB1\u0BCD", -1, -1]
    ];

    var /** boolean */ B_found_vetrumai_urupu = false;
    var /** boolean */ B_found_a_match = false;


    /** @return {boolean} */
    function r_has_min_length() {
        return base.current.length > 4;
    };

    /** @return {boolean} */
    function r_fix_va_start() {
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
                if (!base.slice_from("\u0B93"))
                {
                    return false;
                }
                break;
            case 2:
                if (!base.slice_from("\u0B92"))
                {
                    return false;
                }
                break;
            case 3:
                if (!base.slice_from("\u0B89"))
                {
                    return false;
                }
                break;
            case 4:
                if (!base.slice_from("\u0B8A"))
                {
                    return false;
                }
                break;
        }
        return true;
    };

    /** @return {boolean} */
    function r_fix_endings() {
        /** @const */ var /** number */ v_1 = base.cursor;
        lab0: {
            while(true)
            {
                /** @const */ var /** number */ v_2 = base.cursor;
                lab1: {
                    if (!r_fix_ending())
                    {
                        break lab1;
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
    function r_remove_question_prefixes() {
        base.bra = base.cursor;
        if (!(base.eq_s("\u0B8E")))
        {
            return false;
        }
        if (base.find_among(a_1) == 0)
        {
            return false;
        }
        if (!(base.eq_s("\u0BCD")))
        {
            return false;
        }
        base.ket = base.cursor;
        if (!base.slice_del())
        {
            return false;
        }
        /** @const */ var /** number */ v_1 = base.cursor;
        r_fix_va_start();
        base.cursor = v_1;
        return true;
    };

    /** @return {boolean} */
    function r_fix_ending() {
        var /** number */ among_var;
        if (base.current.length <= 3)
        {
            return false;
        }
        base.limit_backward = base.cursor; base.cursor = base.limit;
        lab0: {
            /** @const */ var /** number */ v_1 = base.limit - base.cursor;
            lab1: {
                base.ket = base.cursor;
                among_var = base.find_among_b(a_5);
                if (among_var == 0)
                {
                    break lab1;
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
                        /** @const */ var /** number */ v_2 = base.limit - base.cursor;
                        if (base.find_among_b(a_2) == 0)
                        {
                            break lab1;
                        }
                        base.cursor = base.limit - v_2;
                        if (!base.slice_del())
                        {
                            return false;
                        }
                        break;
                    case 3:
                        if (!base.slice_from("\u0BB3\u0BCD"))
                        {
                            return false;
                        }
                        break;
                    case 4:
                        if (!base.slice_from("\u0BB2\u0BCD"))
                        {
                            return false;
                        }
                        break;
                    case 5:
                        if (!base.slice_from("\u0B9F\u0BC1"))
                        {
                            return false;
                        }
                        break;
                    case 6:
                        if (!B_found_vetrumai_urupu)
                        {
                            break lab1;
                        }
                        {
                            /** @const */ var /** number */ v_3 = base.limit - base.cursor;
                            lab2: {
                                if (!(base.eq_s_b("\u0BC8")))
                                {
                                    break lab2;
                                }
                                break lab1;
                            }
                            base.cursor = base.limit - v_3;
                        }
                        if (!base.slice_from("\u0BAE\u0BCD"))
                        {
                            return false;
                        }
                        break;
                    case 7:
                        if (!base.slice_from("\u0BCD"))
                        {
                            return false;
                        }
                        break;
                    case 8:
                        {
                            /** @const */ var /** number */ v_4 = base.limit - base.cursor;
                            lab3: {
                                if (base.find_among_b(a_3) == 0)
                                {
                                    break lab3;
                                }
                                break lab1;
                            }
                            base.cursor = base.limit - v_4;
                        }
                        if (!base.slice_del())
                        {
                            return false;
                        }
                        break;
                    case 9:
                        among_var = base.find_among_b(a_4);
                        switch (among_var) {
                            case 1:
                                if (!base.slice_del())
                                {
                                    return false;
                                }
                                break;
                            case 2:
                                if (!base.slice_from("\u0BAE\u0BCD"))
                                {
                                    return false;
                                }
                                break;
                        }
                        break;
                }
                break lab0;
            }
            base.cursor = base.limit - v_1;
            base.ket = base.cursor;
            if (!(base.eq_s_b("\u0BCD")))
            {
                return false;
            }
            lab4: {
                /** @const */ var /** number */ v_5 = base.limit - base.cursor;
                lab5: {
                    if (base.find_among_b(a_6) == 0)
                    {
                        break lab5;
                    }
                    /** @const */ var /** number */ v_6 = base.limit - base.cursor;
                    lab6: {
                        if (!(base.eq_s_b("\u0BCD")))
                        {
                            base.cursor = base.limit - v_6;
                            break lab6;
                        }
                        if (base.find_among_b(a_7) == 0)
                        {
                            base.cursor = base.limit - v_6;
                            break lab6;
                        }
                    }
                    base.bra = base.cursor;
                    if (!base.slice_del())
                    {
                        return false;
                    }
                    break lab4;
                }
                base.cursor = base.limit - v_5;
                lab7: {
                    if (base.find_among_b(a_8) == 0)
                    {
                        break lab7;
                    }
                    base.bra = base.cursor;
                    if (!(base.eq_s_b("\u0BCD")))
                    {
                        break lab7;
                    }
                    if (!base.slice_del())
                    {
                        return false;
                    }
                    break lab4;
                }
                base.cursor = base.limit - v_5;
                /** @const */ var /** number */ v_7 = base.limit - base.cursor;
                if (base.find_among_b(a_9) == 0)
                {
                    return false;
                }
                base.cursor = base.limit - v_7;
                base.bra = base.cursor;
                if (!base.slice_del())
                {
                    return false;
                }
            }
        }
        base.cursor = base.limit_backward;
        return true;
    };

    /** @return {boolean} */
    function r_remove_pronoun_prefixes() {
        base.bra = base.cursor;
        if (base.find_among(a_10) == 0)
        {
            return false;
        }
        if (base.find_among(a_11) == 0)
        {
            return false;
        }
        if (!(base.eq_s("\u0BCD")))
        {
            return false;
        }
        base.ket = base.cursor;
        if (!base.slice_del())
        {
            return false;
        }
        /** @const */ var /** number */ v_1 = base.cursor;
        r_fix_va_start();
        base.cursor = v_1;
        return true;
    };

    /** @return {boolean} */
    function r_remove_plural_suffix() {
        var /** number */ among_var;
        base.limit_backward = base.cursor; base.cursor = base.limit;
        base.ket = base.cursor;
        among_var = base.find_among_b(a_13);
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
                        if (base.find_among_b(a_12) == 0)
                        {
                            break lab1;
                        }
                        if (!base.slice_from("\u0BC1\u0B99\u0BCD"))
                        {
                            return false;
                        }
                        break lab0;
                    }
                    base.cursor = base.limit - v_1;
                    if (!base.slice_from("\u0BCD"))
                    {
                        return false;
                    }
                }
                break;
            case 2:
                if (!base.slice_from("\u0BB2\u0BCD"))
                {
                    return false;
                }
                break;
            case 3:
                if (!base.slice_from("\u0BB3\u0BCD"))
                {
                    return false;
                }
                break;
            case 4:
                if (!base.slice_del())
                {
                    return false;
                }
                break;
        }
        base.cursor = base.limit_backward;
        return true;
    };

    /** @return {boolean} */
    function r_remove_question_suffixes() {
        if (!r_has_min_length())
        {
            return false;
        }
        base.limit_backward = base.cursor; base.cursor = base.limit;
        /** @const */ var /** number */ v_1 = base.limit - base.cursor;
        lab0: {
            base.ket = base.cursor;
            if (base.find_among_b(a_14) == 0)
            {
                break lab0;
            }
            base.bra = base.cursor;
            if (!base.slice_from("\u0BCD"))
            {
                return false;
            }
        }
        base.cursor = base.limit - v_1;
        base.cursor = base.limit_backward;
        r_fix_endings();
        return true;
    };

    /** @return {boolean} */
    function r_remove_command_suffixes() {
        if (!r_has_min_length())
        {
            return false;
        }
        base.limit_backward = base.cursor; base.cursor = base.limit;
        base.ket = base.cursor;
        if (base.find_among_b(a_15) == 0)
        {
            return false;
        }
        base.bra = base.cursor;
        if (!base.slice_del())
        {
            return false;
        }
        base.cursor = base.limit_backward;
        return true;
    };

    /** @return {boolean} */
    function r_remove_um() {
        if (!r_has_min_length())
        {
            return false;
        }
        base.limit_backward = base.cursor; base.cursor = base.limit;
        base.ket = base.cursor;
        if (!(base.eq_s_b("\u0BC1\u0BAE\u0BCD")))
        {
            return false;
        }
        base.bra = base.cursor;
        if (!base.slice_from("\u0BCD"))
        {
            return false;
        }
        base.cursor = base.limit_backward;
        /** @const */ var /** number */ v_1 = base.cursor;
        r_fix_ending();
        base.cursor = v_1;
        return true;
    };

    /** @return {boolean} */
    function r_remove_common_word_endings() {
        var /** number */ among_var;
        if (!r_has_min_length())
        {
            return false;
        }
        base.limit_backward = base.cursor; base.cursor = base.limit;
        base.ket = base.cursor;
        among_var = base.find_among_b(a_17);
        if (among_var == 0)
        {
            return false;
        }
        base.bra = base.cursor;
        switch (among_var) {
            case 1:
                if (!base.slice_from("\u0BCD"))
                {
                    return false;
                }
                break;
            case 2:
                {
                    /** @const */ var /** number */ v_1 = base.limit - base.cursor;
                    lab0: {
                        if (base.find_among_b(a_16) == 0)
                        {
                            break lab0;
                        }
                        return false;
                    }
                    base.cursor = base.limit - v_1;
                }
                if (!base.slice_from("\u0BCD"))
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
        base.cursor = base.limit_backward;
        r_fix_endings();
        return true;
    };

    /** @return {boolean} */
    function r_remove_vetrumai_urupukal() {
        var /** number */ among_var;
        B_found_vetrumai_urupu = false;
        if (!r_has_min_length())
        {
            return false;
        }
        base.limit_backward = base.cursor; base.cursor = base.limit;
        lab0: {
            /** @const */ var /** number */ v_1 = base.limit - base.cursor;
            lab1: {
                /** @const */ var /** number */ v_2 = base.limit - base.cursor;
                base.ket = base.cursor;
                among_var = base.find_among_b(a_20);
                if (among_var == 0)
                {
                    break lab1;
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
                        if (!base.slice_from("\u0BCD"))
                        {
                            return false;
                        }
                        break;
                    case 3:
                        {
                            /** @const */ var /** number */ v_3 = base.limit - base.cursor;
                            lab2: {
                                if (!(base.eq_s_b("\u0BAE")))
                                {
                                    break lab2;
                                }
                                break lab1;
                            }
                            base.cursor = base.limit - v_3;
                        }
                        if (!base.slice_from("\u0BCD"))
                        {
                            return false;
                        }
                        break;
                    case 4:
                        if (base.current.length < 7)
                        {
                            break lab1;
                        }
                        if (!base.slice_from("\u0BCD"))
                        {
                            return false;
                        }
                        break;
                    case 5:
                        {
                            /** @const */ var /** number */ v_4 = base.limit - base.cursor;
                            lab3: {
                                if (base.find_among_b(a_18) == 0)
                                {
                                    break lab3;
                                }
                                break lab1;
                            }
                            base.cursor = base.limit - v_4;
                        }
                        if (!base.slice_from("\u0BCD"))
                        {
                            return false;
                        }
                        break;
                    case 6:
                        {
                            /** @const */ var /** number */ v_5 = base.limit - base.cursor;
                            lab4: {
                                if (base.find_among_b(a_19) == 0)
                                {
                                    break lab4;
                                }
                                break lab1;
                            }
                            base.cursor = base.limit - v_5;
                        }
                        if (!base.slice_del())
                        {
                            return false;
                        }
                        break;
                    case 7:
                        if (!base.slice_from("\u0BBF"))
                        {
                            return false;
                        }
                        break;
                }
                base.cursor = base.limit - v_2;
                break lab0;
            }
            base.cursor = base.limit - v_1;
            /** @const */ var /** number */ v_6 = base.limit - base.cursor;
            base.ket = base.cursor;
            if (!(base.eq_s_b("\u0BC8")))
            {
                return false;
            }
            lab5: {
                /** @const */ var /** number */ v_7 = base.limit - base.cursor;
                lab6: {
                    {
                        /** @const */ var /** number */ v_8 = base.limit - base.cursor;
                        lab7: {
                            if (base.find_among_b(a_21) == 0)
                            {
                                break lab7;
                            }
                            break lab6;
                        }
                        base.cursor = base.limit - v_8;
                    }
                    break lab5;
                }
                base.cursor = base.limit - v_7;
                /** @const */ var /** number */ v_9 = base.limit - base.cursor;
                if (base.find_among_b(a_22) == 0)
                {
                    return false;
                }
                if (!(base.eq_s_b("\u0BCD")))
                {
                    return false;
                }
                base.cursor = base.limit - v_9;
            }
            base.bra = base.cursor;
            if (!base.slice_from("\u0BCD"))
            {
                return false;
            }
            base.cursor = base.limit - v_6;
        }
        B_found_vetrumai_urupu = true;
        /** @const */ var /** number */ v_10 = base.limit - base.cursor;
        lab8: {
            base.ket = base.cursor;
            if (!(base.eq_s_b("\u0BBF\u0BA9\u0BCD")))
            {
                break lab8;
            }
            base.bra = base.cursor;
            if (!base.slice_from("\u0BCD"))
            {
                return false;
            }
        }
        base.cursor = base.limit - v_10;
        base.cursor = base.limit_backward;
        r_fix_endings();
        return true;
    };

    /** @return {boolean} */
    function r_remove_tense_suffixes() {
        while(true)
        {
            /** @const */ var /** number */ v_1 = base.cursor;
            lab0: {
                if (!r_remove_tense_suffix())
                {
                    break lab0;
                }
                continue;
            }
            base.cursor = v_1;
            break;
        }
        return true;
    };

    /** @return {boolean} */
    function r_remove_tense_suffix() {
        var /** number */ among_var;
        B_found_a_match = false;
        if (!r_has_min_length())
        {
            return false;
        }
        base.limit_backward = base.cursor; base.cursor = base.limit;
        /** @const */ var /** number */ v_1 = base.limit - base.cursor;
        lab0: {
            /** @const */ var /** number */ v_2 = base.limit - base.cursor;
            base.ket = base.cursor;
            among_var = base.find_among_b(a_25);
            if (among_var == 0)
            {
                break lab0;
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
                    {
                        /** @const */ var /** number */ v_3 = base.limit - base.cursor;
                        lab1: {
                            if (base.find_among_b(a_23) == 0)
                            {
                                break lab1;
                            }
                            break lab0;
                        }
                        base.cursor = base.limit - v_3;
                    }
                    if (!base.slice_del())
                    {
                        return false;
                    }
                    break;
                case 3:
                    {
                        /** @const */ var /** number */ v_4 = base.limit - base.cursor;
                        lab2: {
                            if (base.find_among_b(a_24) == 0)
                            {
                                break lab2;
                            }
                            break lab0;
                        }
                        base.cursor = base.limit - v_4;
                    }
                    if (!base.slice_del())
                    {
                        return false;
                    }
                    break;
                case 4:
                    {
                        /** @const */ var /** number */ v_5 = base.limit - base.cursor;
                        lab3: {
                            if (!(base.eq_s_b("\u0B9A")))
                            {
                                break lab3;
                            }
                            break lab0;
                        }
                        base.cursor = base.limit - v_5;
                    }
                    if (!base.slice_from("\u0BCD"))
                    {
                        return false;
                    }
                    break;
                case 5:
                    if (!base.slice_from("\u0BCD"))
                    {
                        return false;
                    }
                    break;
                case 6:
                    /** @const */ var /** number */ v_6 = base.limit - base.cursor;
                    if (!(base.eq_s_b("\u0BCD")))
                    {
                        break lab0;
                    }
                    base.cursor = base.limit - v_6;
                    if (!base.slice_del())
                    {
                        return false;
                    }
                    break;
            }
            B_found_a_match = true;
            base.cursor = base.limit - v_2;
        }
        base.cursor = base.limit - v_1;
        /** @const */ var /** number */ v_7 = base.limit - base.cursor;
        lab4: {
            base.ket = base.cursor;
            if (base.find_among_b(a_26) == 0)
            {
                break lab4;
            }
            base.bra = base.cursor;
            if (!base.slice_del())
            {
                return false;
            }
            B_found_a_match = true;
        }
        base.cursor = base.limit - v_7;
        base.cursor = base.limit_backward;
        r_fix_endings();
        if (!B_found_a_match)
        {
            return false;
        }
        return true;
    };

    this.stem = /** @return {boolean} */ function() {
        B_found_vetrumai_urupu = false;
        /** @const */ var /** number */ v_1 = base.cursor;
        r_fix_ending();
        base.cursor = v_1;
        if (!r_has_min_length())
        {
            return false;
        }
        /** @const */ var /** number */ v_2 = base.cursor;
        r_remove_question_prefixes();
        base.cursor = v_2;
        /** @const */ var /** number */ v_3 = base.cursor;
        r_remove_pronoun_prefixes();
        base.cursor = v_3;
        /** @const */ var /** number */ v_4 = base.cursor;
        r_remove_question_suffixes();
        base.cursor = v_4;
        /** @const */ var /** number */ v_5 = base.cursor;
        r_remove_um();
        base.cursor = v_5;
        /** @const */ var /** number */ v_6 = base.cursor;
        r_remove_common_word_endings();
        base.cursor = v_6;
        /** @const */ var /** number */ v_7 = base.cursor;
        r_remove_vetrumai_urupukal();
        base.cursor = v_7;
        /** @const */ var /** number */ v_8 = base.cursor;
        r_remove_plural_suffix();
        base.cursor = v_8;
        /** @const */ var /** number */ v_9 = base.cursor;
        r_remove_command_suffixes();
        base.cursor = v_9;
        /** @const */ var /** number */ v_10 = base.cursor;
        r_remove_tense_suffixes();
        base.cursor = v_10;
        return true;
    };

    /**@return{string}*/
    this['stemWord'] = function(/**string*/word) {
        base.setCurrent(word);
        this.stem();
        return base.getCurrent();
    };
};
