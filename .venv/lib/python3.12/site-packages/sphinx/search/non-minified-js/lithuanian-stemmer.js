// Generated from lithuanian.sbl by Snowball 3.0.1 - https://snowballstem.org/

/**@constructor*/
var LithuanianStemmer = function() {
    var base = new BaseStemmer();

    /** @const */ var a_0 = [
        ["a", -1, -1],
        ["ia", 0, -1],
        ["eria", 1, -1],
        ["osna", 0, -1],
        ["iosna", 3, -1],
        ["uosna", 3, -1],
        ["iuosna", 5, -1],
        ["ysna", 0, -1],
        ["\u0117sna", 0, -1],
        ["e", -1, -1],
        ["ie", 9, -1],
        ["enie", 10, -1],
        ["erie", 10, -1],
        ["oje", 9, -1],
        ["ioje", 13, -1],
        ["uje", 9, -1],
        ["iuje", 15, -1],
        ["yje", 9, -1],
        ["enyje", 17, -1],
        ["eryje", 17, -1],
        ["\u0117je", 9, -1],
        ["ame", 9, -1],
        ["iame", 21, -1],
        ["sime", 9, -1],
        ["ome", 9, -1],
        ["\u0117me", 9, -1],
        ["tum\u0117me", 25, -1],
        ["ose", 9, -1],
        ["iose", 27, -1],
        ["uose", 27, -1],
        ["iuose", 29, -1],
        ["yse", 9, -1],
        ["enyse", 31, -1],
        ["eryse", 31, -1],
        ["\u0117se", 9, -1],
        ["ate", 9, -1],
        ["iate", 35, -1],
        ["ite", 9, -1],
        ["kite", 37, -1],
        ["site", 37, -1],
        ["ote", 9, -1],
        ["tute", 9, -1],
        ["\u0117te", 9, -1],
        ["tum\u0117te", 42, -1],
        ["i", -1, -1],
        ["ai", 44, -1],
        ["iai", 45, -1],
        ["eriai", 46, -1],
        ["ei", 44, -1],
        ["tumei", 48, -1],
        ["ki", 44, -1],
        ["imi", 44, -1],
        ["erimi", 51, -1],
        ["umi", 44, -1],
        ["iumi", 53, -1],
        ["si", 44, -1],
        ["asi", 55, -1],
        ["iasi", 56, -1],
        ["esi", 55, -1],
        ["iesi", 58, -1],
        ["siesi", 59, -1],
        ["isi", 55, -1],
        ["aisi", 61, -1],
        ["eisi", 61, -1],
        ["tumeisi", 63, -1],
        ["uisi", 61, -1],
        ["osi", 55, -1],
        ["\u0117josi", 66, -1],
        ["uosi", 66, -1],
        ["iuosi", 68, -1],
        ["siuosi", 69, -1],
        ["usi", 55, -1],
        ["ausi", 71, -1],
        ["\u010Diausi", 72, -1],
        ["\u0105si", 55, -1],
        ["\u0117si", 55, -1],
        ["\u0173si", 55, -1],
        ["t\u0173si", 76, -1],
        ["ti", 44, -1],
        ["enti", 78, -1],
        ["inti", 78, -1],
        ["oti", 78, -1],
        ["ioti", 81, -1],
        ["uoti", 81, -1],
        ["iuoti", 83, -1],
        ["auti", 78, -1],
        ["iauti", 85, -1],
        ["yti", 78, -1],
        ["\u0117ti", 78, -1],
        ["tel\u0117ti", 88, -1],
        ["in\u0117ti", 88, -1],
        ["ter\u0117ti", 88, -1],
        ["ui", 44, -1],
        ["iui", 92, -1],
        ["eniui", 93, -1],
        ["oj", -1, -1],
        ["\u0117j", -1, -1],
        ["k", -1, -1],
        ["am", -1, -1],
        ["iam", 98, -1],
        ["iem", -1, -1],
        ["im", -1, -1],
        ["sim", 101, -1],
        ["om", -1, -1],
        ["tum", -1, -1],
        ["\u0117m", -1, -1],
        ["tum\u0117m", 105, -1],
        ["an", -1, -1],
        ["on", -1, -1],
        ["ion", 108, -1],
        ["un", -1, -1],
        ["iun", 110, -1],
        ["\u0117n", -1, -1],
        ["o", -1, -1],
        ["io", 113, -1],
        ["enio", 114, -1],
        ["\u0117jo", 113, -1],
        ["uo", 113, -1],
        ["s", -1, -1],
        ["as", 118, -1],
        ["ias", 119, -1],
        ["es", 118, -1],
        ["ies", 121, -1],
        ["is", 118, -1],
        ["ais", 123, -1],
        ["iais", 124, -1],
        ["tumeis", 123, -1],
        ["imis", 123, -1],
        ["enimis", 127, -1],
        ["omis", 123, -1],
        ["iomis", 129, -1],
        ["umis", 123, -1],
        ["\u0117mis", 123, -1],
        ["enis", 123, -1],
        ["asis", 123, -1],
        ["ysis", 123, -1],
        ["ams", 118, -1],
        ["iams", 136, -1],
        ["iems", 118, -1],
        ["ims", 118, -1],
        ["enims", 139, -1],
        ["erims", 139, -1],
        ["oms", 118, -1],
        ["ioms", 142, -1],
        ["ums", 118, -1],
        ["\u0117ms", 118, -1],
        ["ens", 118, -1],
        ["os", 118, -1],
        ["ios", 147, -1],
        ["uos", 147, -1],
        ["iuos", 149, -1],
        ["ers", 118, -1],
        ["us", 118, -1],
        ["aus", 152, -1],
        ["iaus", 153, -1],
        ["ius", 152, -1],
        ["ys", 118, -1],
        ["enys", 156, -1],
        ["erys", 156, -1],
        ["\u0105s", 118, -1],
        ["i\u0105s", 159, -1],
        ["\u0117s", 118, -1],
        ["am\u0117s", 161, -1],
        ["iam\u0117s", 162, -1],
        ["im\u0117s", 161, -1],
        ["kim\u0117s", 164, -1],
        ["sim\u0117s", 164, -1],
        ["om\u0117s", 161, -1],
        ["\u0117m\u0117s", 161, -1],
        ["tum\u0117m\u0117s", 168, -1],
        ["at\u0117s", 161, -1],
        ["iat\u0117s", 170, -1],
        ["sit\u0117s", 161, -1],
        ["ot\u0117s", 161, -1],
        ["\u0117t\u0117s", 161, -1],
        ["tum\u0117t\u0117s", 174, -1],
        ["\u012Fs", 118, -1],
        ["\u016Bs", 118, -1],
        ["t\u0173s", 118, -1],
        ["at", -1, -1],
        ["iat", 179, -1],
        ["it", -1, -1],
        ["sit", 181, -1],
        ["ot", -1, -1],
        ["\u0117t", -1, -1],
        ["tum\u0117t", 184, -1],
        ["u", -1, -1],
        ["au", 186, -1],
        ["iau", 187, -1],
        ["\u010Diau", 188, -1],
        ["iu", 186, -1],
        ["eniu", 190, -1],
        ["siu", 190, -1],
        ["y", -1, -1],
        ["\u0105", -1, -1],
        ["i\u0105", 194, -1],
        ["\u0117", -1, -1],
        ["\u0119", -1, -1],
        ["\u012F", -1, -1],
        ["en\u012F", 198, -1],
        ["er\u012F", 198, -1],
        ["\u0173", -1, -1],
        ["i\u0173", 201, -1],
        ["er\u0173", 201, -1]
    ];

    /** @const */ var a_1 = [
        ["ing", -1, -1],
        ["aj", -1, -1],
        ["iaj", 1, -1],
        ["iej", -1, -1],
        ["oj", -1, -1],
        ["ioj", 4, -1],
        ["uoj", 4, -1],
        ["iuoj", 6, -1],
        ["auj", -1, -1],
        ["\u0105j", -1, -1],
        ["i\u0105j", 9, -1],
        ["\u0117j", -1, -1],
        ["\u0173j", -1, -1],
        ["i\u0173j", 12, -1],
        ["ok", -1, -1],
        ["iok", 14, -1],
        ["iuk", -1, -1],
        ["uliuk", 16, -1],
        ["u\u010Diuk", 16, -1],
        ["i\u0161k", -1, -1],
        ["iul", -1, -1],
        ["yl", -1, -1],
        ["\u0117l", -1, -1],
        ["am", -1, -1],
        ["dam", 23, -1],
        ["jam", 23, -1],
        ["zgan", -1, -1],
        ["ain", -1, -1],
        ["esn", -1, -1],
        ["op", -1, -1],
        ["iop", 29, -1],
        ["ias", -1, -1],
        ["ies", -1, -1],
        ["ais", -1, -1],
        ["iais", 33, -1],
        ["os", -1, -1],
        ["ios", 35, -1],
        ["uos", 35, -1],
        ["iuos", 37, -1],
        ["aus", -1, -1],
        ["iaus", 39, -1],
        ["\u0105s", -1, -1],
        ["i\u0105s", 41, -1],
        ["\u0119s", -1, -1],
        ["ut\u0117ait", -1, -1],
        ["ant", -1, -1],
        ["iant", 45, -1],
        ["siant", 46, -1],
        ["int", -1, -1],
        ["ot", -1, -1],
        ["uot", 49, -1],
        ["iuot", 50, -1],
        ["yt", -1, -1],
        ["\u0117t", -1, -1],
        ["yk\u0161t", -1, -1],
        ["iau", -1, -1],
        ["dav", -1, -1],
        ["sv", -1, -1],
        ["\u0161v", -1, -1],
        ["yk\u0161\u010D", -1, -1],
        ["\u0119", -1, -1],
        ["\u0117j\u0119", 60, -1]
    ];

    /** @const */ var a_2 = [
        ["ojime", -1, 7],
        ["\u0117jime", -1, 3],
        ["avime", -1, 6],
        ["okate", -1, 8],
        ["aite", -1, 1],
        ["uote", -1, 2],
        ["asius", -1, 5],
        ["okat\u0117s", -1, 8],
        ["ait\u0117s", -1, 1],
        ["uot\u0117s", -1, 2],
        ["esiu", -1, 4]
    ];

    /** @const */ var a_3 = [
        ["\u010D", -1, 1],
        ["d\u017E", -1, 2]
    ];

    /** @const */ var a_4 = [
        ["gd", -1, 1]
    ];

    /** @const */ var /** Array<int> */ g_v = [17, 65, 16, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16, 0, 64, 1, 0, 64, 0, 0, 0, 0, 0, 0, 0, 4, 4];

    var /** number */ I_p1 = 0;


    /** @return {boolean} */
    function r_step1() {
        if (base.cursor < I_p1)
        {
            return false;
        }
        /** @const */ var /** number */ v_1 = base.limit_backward;
        base.limit_backward = I_p1;
        base.ket = base.cursor;
        if (base.find_among_b(a_0) == 0)
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
    function r_step2() {
        while(true)
        {
            /** @const */ var /** number */ v_1 = base.limit - base.cursor;
            lab0: {
                if (base.cursor < I_p1)
                {
                    break lab0;
                }
                /** @const */ var /** number */ v_2 = base.limit_backward;
                base.limit_backward = I_p1;
                base.ket = base.cursor;
                if (base.find_among_b(a_1) == 0)
                {
                    base.limit_backward = v_2;
                    break lab0;
                }
                base.bra = base.cursor;
                base.limit_backward = v_2;
                if (!base.slice_del())
                {
                    return false;
                }
                continue;
            }
            base.cursor = base.limit - v_1;
            break;
        }
        return true;
    };

    /** @return {boolean} */
    function r_fix_conflicts() {
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
                if (!base.slice_from("ait\u0117"))
                {
                    return false;
                }
                break;
            case 2:
                if (!base.slice_from("uot\u0117"))
                {
                    return false;
                }
                break;
            case 3:
                if (!base.slice_from("\u0117jimas"))
                {
                    return false;
                }
                break;
            case 4:
                if (!base.slice_from("esys"))
                {
                    return false;
                }
                break;
            case 5:
                if (!base.slice_from("asys"))
                {
                    return false;
                }
                break;
            case 6:
                if (!base.slice_from("avimas"))
                {
                    return false;
                }
                break;
            case 7:
                if (!base.slice_from("ojimas"))
                {
                    return false;
                }
                break;
            case 8:
                if (!base.slice_from("okat\u0117"))
                {
                    return false;
                }
                break;
        }
        return true;
    };

    /** @return {boolean} */
    function r_fix_chdz() {
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
                if (!base.slice_from("t"))
                {
                    return false;
                }
                break;
            case 2:
                if (!base.slice_from("d"))
                {
                    return false;
                }
                break;
        }
        return true;
    };

    /** @return {boolean} */
    function r_fix_gd() {
        base.ket = base.cursor;
        if (base.find_among_b(a_4) == 0)
        {
            return false;
        }
        base.bra = base.cursor;
        if (!base.slice_from("g"))
        {
            return false;
        }
        return true;
    };

    this.stem = /** @return {boolean} */ function() {
        I_p1 = base.limit;
        /** @const */ var /** number */ v_1 = base.cursor;
        lab0: {
            /** @const */ var /** number */ v_2 = base.cursor;
            lab1: {
                /** @const */ var /** number */ v_3 = base.cursor;
                if (!(base.eq_s("a")))
                {
                    base.cursor = v_2;
                    break lab1;
                }
                base.cursor = v_3;
                if (base.current.length <= 6)
                {
                    base.cursor = v_2;
                    break lab1;
                }
                if (base.cursor >= base.limit)
                {
                    base.cursor = v_2;
                    break lab1;
                }
                base.cursor++;
            }
            if (!base.go_out_grouping(g_v, 97, 371))
            {
                break lab0;
            }
            base.cursor++;
            if (!base.go_in_grouping(g_v, 97, 371))
            {
                break lab0;
            }
            base.cursor++;
            I_p1 = base.cursor;
        }
        base.cursor = v_1;
        base.limit_backward = base.cursor; base.cursor = base.limit;
        /** @const */ var /** number */ v_4 = base.limit - base.cursor;
        r_fix_conflicts();
        base.cursor = base.limit - v_4;
        /** @const */ var /** number */ v_5 = base.limit - base.cursor;
        r_step1();
        base.cursor = base.limit - v_5;
        /** @const */ var /** number */ v_6 = base.limit - base.cursor;
        r_fix_chdz();
        base.cursor = base.limit - v_6;
        /** @const */ var /** number */ v_7 = base.limit - base.cursor;
        r_step2();
        base.cursor = base.limit - v_7;
        /** @const */ var /** number */ v_8 = base.limit - base.cursor;
        r_fix_chdz();
        base.cursor = base.limit - v_8;
        /** @const */ var /** number */ v_9 = base.limit - base.cursor;
        r_fix_gd();
        base.cursor = base.limit - v_9;
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
