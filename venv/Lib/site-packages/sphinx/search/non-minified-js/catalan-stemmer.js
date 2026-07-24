// Generated from catalan.sbl by Snowball 3.0.1 - https://snowballstem.org/

/**@constructor*/
var CatalanStemmer = function() {
    var base = new BaseStemmer();

    /** @const */ var a_0 = [
        ["", -1, 7],
        ["\u00B7", 0, 6],
        ["\u00E0", 0, 1],
        ["\u00E1", 0, 1],
        ["\u00E8", 0, 2],
        ["\u00E9", 0, 2],
        ["\u00EC", 0, 3],
        ["\u00ED", 0, 3],
        ["\u00EF", 0, 3],
        ["\u00F2", 0, 4],
        ["\u00F3", 0, 4],
        ["\u00FA", 0, 5],
        ["\u00FC", 0, 5]
    ];

    /** @const */ var a_1 = [
        ["la", -1, 1],
        ["-la", 0, 1],
        ["sela", 0, 1],
        ["le", -1, 1],
        ["me", -1, 1],
        ["-me", 4, 1],
        ["se", -1, 1],
        ["-te", -1, 1],
        ["hi", -1, 1],
        ["'hi", 8, 1],
        ["li", -1, 1],
        ["-li", 10, 1],
        ["'l", -1, 1],
        ["'m", -1, 1],
        ["-m", -1, 1],
        ["'n", -1, 1],
        ["-n", -1, 1],
        ["ho", -1, 1],
        ["'ho", 17, 1],
        ["lo", -1, 1],
        ["selo", 19, 1],
        ["'s", -1, 1],
        ["las", -1, 1],
        ["selas", 22, 1],
        ["les", -1, 1],
        ["-les", 24, 1],
        ["'ls", -1, 1],
        ["-ls", -1, 1],
        ["'ns", -1, 1],
        ["-ns", -1, 1],
        ["ens", -1, 1],
        ["los", -1, 1],
        ["selos", 31, 1],
        ["nos", -1, 1],
        ["-nos", 33, 1],
        ["vos", -1, 1],
        ["us", -1, 1],
        ["-us", 36, 1],
        ["'t", -1, 1]
    ];

    /** @const */ var a_2 = [
        ["ica", -1, 4],
        ["l\u00F3gica", 0, 3],
        ["enca", -1, 1],
        ["ada", -1, 2],
        ["ancia", -1, 1],
        ["encia", -1, 1],
        ["\u00E8ncia", -1, 1],
        ["\u00EDcia", -1, 1],
        ["logia", -1, 3],
        ["inia", -1, 1],
        ["\u00EDinia", 9, 1],
        ["eria", -1, 1],
        ["\u00E0ria", -1, 1],
        ["at\u00F2ria", -1, 1],
        ["alla", -1, 1],
        ["ella", -1, 1],
        ["\u00EDvola", -1, 1],
        ["ima", -1, 1],
        ["\u00EDssima", 17, 1],
        ["qu\u00EDssima", 18, 5],
        ["ana", -1, 1],
        ["ina", -1, 1],
        ["era", -1, 1],
        ["sfera", 22, 1],
        ["ora", -1, 1],
        ["dora", 24, 1],
        ["adora", 25, 1],
        ["adura", -1, 1],
        ["esa", -1, 1],
        ["osa", -1, 1],
        ["assa", -1, 1],
        ["essa", -1, 1],
        ["issa", -1, 1],
        ["eta", -1, 1],
        ["ita", -1, 1],
        ["ota", -1, 1],
        ["ista", -1, 1],
        ["ialista", 36, 1],
        ["ionista", 36, 1],
        ["iva", -1, 1],
        ["ativa", 39, 1],
        ["n\u00E7a", -1, 1],
        ["log\u00EDa", -1, 3],
        ["ic", -1, 4],
        ["\u00EDstic", 43, 1],
        ["enc", -1, 1],
        ["esc", -1, 1],
        ["ud", -1, 1],
        ["atge", -1, 1],
        ["ble", -1, 1],
        ["able", 49, 1],
        ["ible", 49, 1],
        ["isme", -1, 1],
        ["ialisme", 52, 1],
        ["ionisme", 52, 1],
        ["ivisme", 52, 1],
        ["aire", -1, 1],
        ["icte", -1, 1],
        ["iste", -1, 1],
        ["ici", -1, 1],
        ["\u00EDci", -1, 1],
        ["logi", -1, 3],
        ["ari", -1, 1],
        ["tori", -1, 1],
        ["al", -1, 1],
        ["il", -1, 1],
        ["all", -1, 1],
        ["ell", -1, 1],
        ["\u00EDvol", -1, 1],
        ["isam", -1, 1],
        ["issem", -1, 1],
        ["\u00ECssem", -1, 1],
        ["\u00EDssem", -1, 1],
        ["\u00EDssim", -1, 1],
        ["qu\u00EDssim", 73, 5],
        ["amen", -1, 1],
        ["\u00ECssin", -1, 1],
        ["ar", -1, 1],
        ["ificar", 77, 1],
        ["egar", 77, 1],
        ["ejar", 77, 1],
        ["itar", 77, 1],
        ["itzar", 77, 1],
        ["fer", -1, 1],
        ["or", -1, 1],
        ["dor", 84, 1],
        ["dur", -1, 1],
        ["doras", -1, 1],
        ["ics", -1, 4],
        ["l\u00F3gics", 88, 3],
        ["uds", -1, 1],
        ["nces", -1, 1],
        ["ades", -1, 2],
        ["ancies", -1, 1],
        ["encies", -1, 1],
        ["\u00E8ncies", -1, 1],
        ["\u00EDcies", -1, 1],
        ["logies", -1, 3],
        ["inies", -1, 1],
        ["\u00EDnies", -1, 1],
        ["eries", -1, 1],
        ["\u00E0ries", -1, 1],
        ["at\u00F2ries", -1, 1],
        ["bles", -1, 1],
        ["ables", 103, 1],
        ["ibles", 103, 1],
        ["imes", -1, 1],
        ["\u00EDssimes", 106, 1],
        ["qu\u00EDssimes", 107, 5],
        ["formes", -1, 1],
        ["ismes", -1, 1],
        ["ialismes", 110, 1],
        ["ines", -1, 1],
        ["eres", -1, 1],
        ["ores", -1, 1],
        ["dores", 114, 1],
        ["idores", 115, 1],
        ["dures", -1, 1],
        ["eses", -1, 1],
        ["oses", -1, 1],
        ["asses", -1, 1],
        ["ictes", -1, 1],
        ["ites", -1, 1],
        ["otes", -1, 1],
        ["istes", -1, 1],
        ["ialistes", 124, 1],
        ["ionistes", 124, 1],
        ["iques", -1, 4],
        ["l\u00F3giques", 127, 3],
        ["ives", -1, 1],
        ["atives", 129, 1],
        ["log\u00EDes", -1, 3],
        ["alleng\u00FCes", -1, 1],
        ["icis", -1, 1],
        ["\u00EDcis", -1, 1],
        ["logis", -1, 3],
        ["aris", -1, 1],
        ["toris", -1, 1],
        ["ls", -1, 1],
        ["als", 138, 1],
        ["ells", 138, 1],
        ["ims", -1, 1],
        ["\u00EDssims", 141, 1],
        ["qu\u00EDssims", 142, 5],
        ["ions", -1, 1],
        ["cions", 144, 1],
        ["acions", 145, 2],
        ["esos", -1, 1],
        ["osos", -1, 1],
        ["assos", -1, 1],
        ["issos", -1, 1],
        ["ers", -1, 1],
        ["ors", -1, 1],
        ["dors", 152, 1],
        ["adors", 153, 1],
        ["idors", 153, 1],
        ["ats", -1, 1],
        ["itats", 156, 1],
        ["bilitats", 157, 1],
        ["ivitats", 157, 1],
        ["ativitats", 159, 1],
        ["\u00EFtats", 156, 1],
        ["ets", -1, 1],
        ["ants", -1, 1],
        ["ents", -1, 1],
        ["ments", 164, 1],
        ["aments", 165, 1],
        ["ots", -1, 1],
        ["uts", -1, 1],
        ["ius", -1, 1],
        ["trius", 169, 1],
        ["atius", 169, 1],
        ["\u00E8s", -1, 1],
        ["\u00E9s", -1, 1],
        ["\u00EDs", -1, 1],
        ["d\u00EDs", 174, 1],
        ["\u00F3s", -1, 1],
        ["itat", -1, 1],
        ["bilitat", 177, 1],
        ["ivitat", 177, 1],
        ["ativitat", 179, 1],
        ["\u00EFtat", -1, 1],
        ["et", -1, 1],
        ["ant", -1, 1],
        ["ent", -1, 1],
        ["ient", 184, 1],
        ["ment", 184, 1],
        ["ament", 186, 1],
        ["isament", 187, 1],
        ["ot", -1, 1],
        ["isseu", -1, 1],
        ["\u00ECsseu", -1, 1],
        ["\u00EDsseu", -1, 1],
        ["triu", -1, 1],
        ["\u00EDssiu", -1, 1],
        ["atiu", -1, 1],
        ["\u00F3", -1, 1],
        ["i\u00F3", 196, 1],
        ["ci\u00F3", 197, 1],
        ["aci\u00F3", 198, 1]
    ];

    /** @const */ var a_3 = [
        ["aba", -1, 1],
        ["esca", -1, 1],
        ["isca", -1, 1],
        ["\u00EFsca", -1, 1],
        ["ada", -1, 1],
        ["ida", -1, 1],
        ["uda", -1, 1],
        ["\u00EFda", -1, 1],
        ["ia", -1, 1],
        ["aria", 8, 1],
        ["iria", 8, 1],
        ["ara", -1, 1],
        ["iera", -1, 1],
        ["ira", -1, 1],
        ["adora", -1, 1],
        ["\u00EFra", -1, 1],
        ["ava", -1, 1],
        ["ixa", -1, 1],
        ["itza", -1, 1],
        ["\u00EDa", -1, 1],
        ["ar\u00EDa", 19, 1],
        ["er\u00EDa", 19, 1],
        ["ir\u00EDa", 19, 1],
        ["\u00EFa", -1, 1],
        ["isc", -1, 1],
        ["\u00EFsc", -1, 1],
        ["ad", -1, 1],
        ["ed", -1, 1],
        ["id", -1, 1],
        ["ie", -1, 1],
        ["re", -1, 1],
        ["dre", 30, 1],
        ["ase", -1, 1],
        ["iese", -1, 1],
        ["aste", -1, 1],
        ["iste", -1, 1],
        ["ii", -1, 1],
        ["ini", -1, 1],
        ["esqui", -1, 1],
        ["eixi", -1, 1],
        ["itzi", -1, 1],
        ["am", -1, 1],
        ["em", -1, 1],
        ["arem", 42, 1],
        ["irem", 42, 1],
        ["\u00E0rem", 42, 1],
        ["\u00EDrem", 42, 1],
        ["\u00E0ssem", 42, 1],
        ["\u00E9ssem", 42, 1],
        ["iguem", 42, 1],
        ["\u00EFguem", 42, 1],
        ["avem", 42, 1],
        ["\u00E0vem", 42, 1],
        ["\u00E1vem", 42, 1],
        ["ir\u00ECem", 42, 1],
        ["\u00EDem", 42, 1],
        ["ar\u00EDem", 55, 1],
        ["ir\u00EDem", 55, 1],
        ["assim", -1, 1],
        ["essim", -1, 1],
        ["issim", -1, 1],
        ["\u00E0ssim", -1, 1],
        ["\u00E8ssim", -1, 1],
        ["\u00E9ssim", -1, 1],
        ["\u00EDssim", -1, 1],
        ["\u00EFm", -1, 1],
        ["an", -1, 1],
        ["aban", 66, 1],
        ["arian", 66, 1],
        ["aran", 66, 1],
        ["ieran", 66, 1],
        ["iran", 66, 1],
        ["\u00EDan", 66, 1],
        ["ar\u00EDan", 72, 1],
        ["er\u00EDan", 72, 1],
        ["ir\u00EDan", 72, 1],
        ["en", -1, 1],
        ["ien", 76, 1],
        ["arien", 77, 1],
        ["irien", 77, 1],
        ["aren", 76, 1],
        ["eren", 76, 1],
        ["iren", 76, 1],
        ["\u00E0ren", 76, 1],
        ["\u00EFren", 76, 1],
        ["asen", 76, 1],
        ["iesen", 76, 1],
        ["assen", 76, 1],
        ["essen", 76, 1],
        ["issen", 76, 1],
        ["\u00E9ssen", 76, 1],
        ["\u00EFssen", 76, 1],
        ["esquen", 76, 1],
        ["isquen", 76, 1],
        ["\u00EFsquen", 76, 1],
        ["aven", 76, 1],
        ["ixen", 76, 1],
        ["eixen", 96, 1],
        ["\u00EFxen", 76, 1],
        ["\u00EFen", 76, 1],
        ["in", -1, 1],
        ["inin", 100, 1],
        ["sin", 100, 1],
        ["isin", 102, 1],
        ["assin", 102, 1],
        ["essin", 102, 1],
        ["issin", 102, 1],
        ["\u00EFssin", 102, 1],
        ["esquin", 100, 1],
        ["eixin", 100, 1],
        ["aron", -1, 1],
        ["ieron", -1, 1],
        ["ar\u00E1n", -1, 1],
        ["er\u00E1n", -1, 1],
        ["ir\u00E1n", -1, 1],
        ["i\u00EFn", -1, 1],
        ["ado", -1, 1],
        ["ido", -1, 1],
        ["ando", -1, 2],
        ["iendo", -1, 1],
        ["io", -1, 1],
        ["ixo", -1, 1],
        ["eixo", 121, 1],
        ["\u00EFxo", -1, 1],
        ["itzo", -1, 1],
        ["ar", -1, 1],
        ["tzar", 125, 1],
        ["er", -1, 1],
        ["eixer", 127, 1],
        ["ir", -1, 1],
        ["ador", -1, 1],
        ["as", -1, 1],
        ["abas", 131, 1],
        ["adas", 131, 1],
        ["idas", 131, 1],
        ["aras", 131, 1],
        ["ieras", 131, 1],
        ["\u00EDas", 131, 1],
        ["ar\u00EDas", 137, 1],
        ["er\u00EDas", 137, 1],
        ["ir\u00EDas", 137, 1],
        ["ids", -1, 1],
        ["es", -1, 1],
        ["ades", 142, 1],
        ["ides", 142, 1],
        ["udes", 142, 1],
        ["\u00EFdes", 142, 1],
        ["atges", 142, 1],
        ["ies", 142, 1],
        ["aries", 148, 1],
        ["iries", 148, 1],
        ["ares", 142, 1],
        ["ires", 142, 1],
        ["adores", 142, 1],
        ["\u00EFres", 142, 1],
        ["ases", 142, 1],
        ["ieses", 142, 1],
        ["asses", 142, 1],
        ["esses", 142, 1],
        ["isses", 142, 1],
        ["\u00EFsses", 142, 1],
        ["ques", 142, 1],
        ["esques", 161, 1],
        ["\u00EFsques", 161, 1],
        ["aves", 142, 1],
        ["ixes", 142, 1],
        ["eixes", 165, 1],
        ["\u00EFxes", 142, 1],
        ["\u00EFes", 142, 1],
        ["abais", -1, 1],
        ["arais", -1, 1],
        ["ierais", -1, 1],
        ["\u00EDais", -1, 1],
        ["ar\u00EDais", 172, 1],
        ["er\u00EDais", 172, 1],
        ["ir\u00EDais", 172, 1],
        ["aseis", -1, 1],
        ["ieseis", -1, 1],
        ["asteis", -1, 1],
        ["isteis", -1, 1],
        ["inis", -1, 1],
        ["sis", -1, 1],
        ["isis", 181, 1],
        ["assis", 181, 1],
        ["essis", 181, 1],
        ["issis", 181, 1],
        ["\u00EFssis", 181, 1],
        ["esquis", -1, 1],
        ["eixis", -1, 1],
        ["itzis", -1, 1],
        ["\u00E1is", -1, 1],
        ["ar\u00E9is", -1, 1],
        ["er\u00E9is", -1, 1],
        ["ir\u00E9is", -1, 1],
        ["ams", -1, 1],
        ["ados", -1, 1],
        ["idos", -1, 1],
        ["amos", -1, 1],
        ["\u00E1bamos", 197, 1],
        ["\u00E1ramos", 197, 1],
        ["i\u00E9ramos", 197, 1],
        ["\u00EDamos", 197, 1],
        ["ar\u00EDamos", 201, 1],
        ["er\u00EDamos", 201, 1],
        ["ir\u00EDamos", 201, 1],
        ["aremos", -1, 1],
        ["eremos", -1, 1],
        ["iremos", -1, 1],
        ["\u00E1semos", -1, 1],
        ["i\u00E9semos", -1, 1],
        ["imos", -1, 1],
        ["adors", -1, 1],
        ["ass", -1, 1],
        ["erass", 212, 1],
        ["ess", -1, 1],
        ["ats", -1, 1],
        ["its", -1, 1],
        ["ents", -1, 1],
        ["\u00E0s", -1, 1],
        ["ar\u00E0s", 218, 1],
        ["ir\u00E0s", 218, 1],
        ["ar\u00E1s", -1, 1],
        ["er\u00E1s", -1, 1],
        ["ir\u00E1s", -1, 1],
        ["\u00E9s", -1, 1],
        ["ar\u00E9s", 224, 1],
        ["\u00EDs", -1, 1],
        ["i\u00EFs", -1, 1],
        ["at", -1, 1],
        ["it", -1, 1],
        ["ant", -1, 1],
        ["ent", -1, 1],
        ["int", -1, 1],
        ["ut", -1, 1],
        ["\u00EFt", -1, 1],
        ["au", -1, 1],
        ["erau", 235, 1],
        ["ieu", -1, 1],
        ["ineu", -1, 1],
        ["areu", -1, 1],
        ["ireu", -1, 1],
        ["\u00E0reu", -1, 1],
        ["\u00EDreu", -1, 1],
        ["asseu", -1, 1],
        ["esseu", -1, 1],
        ["eresseu", 244, 1],
        ["\u00E0sseu", -1, 1],
        ["\u00E9sseu", -1, 1],
        ["igueu", -1, 1],
        ["\u00EFgueu", -1, 1],
        ["\u00E0veu", -1, 1],
        ["\u00E1veu", -1, 1],
        ["itzeu", -1, 1],
        ["\u00ECeu", -1, 1],
        ["ir\u00ECeu", 253, 1],
        ["\u00EDeu", -1, 1],
        ["ar\u00EDeu", 255, 1],
        ["ir\u00EDeu", 255, 1],
        ["assiu", -1, 1],
        ["issiu", -1, 1],
        ["\u00E0ssiu", -1, 1],
        ["\u00E8ssiu", -1, 1],
        ["\u00E9ssiu", -1, 1],
        ["\u00EDssiu", -1, 1],
        ["\u00EFu", -1, 1],
        ["ix", -1, 1],
        ["eix", 265, 1],
        ["\u00EFx", -1, 1],
        ["itz", -1, 1],
        ["i\u00E0", -1, 1],
        ["ar\u00E0", -1, 1],
        ["ir\u00E0", -1, 1],
        ["itz\u00E0", -1, 1],
        ["ar\u00E1", -1, 1],
        ["er\u00E1", -1, 1],
        ["ir\u00E1", -1, 1],
        ["ir\u00E8", -1, 1],
        ["ar\u00E9", -1, 1],
        ["er\u00E9", -1, 1],
        ["ir\u00E9", -1, 1],
        ["\u00ED", -1, 1],
        ["i\u00EF", -1, 1],
        ["i\u00F3", -1, 1]
    ];

    /** @const */ var a_4 = [
        ["a", -1, 1],
        ["e", -1, 1],
        ["i", -1, 1],
        ["\u00EFn", -1, 1],
        ["o", -1, 1],
        ["ir", -1, 1],
        ["s", -1, 1],
        ["is", 6, 1],
        ["os", 6, 1],
        ["\u00EFs", 6, 1],
        ["it", -1, 1],
        ["eu", -1, 1],
        ["iu", -1, 1],
        ["iqu", -1, 2],
        ["itz", -1, 1],
        ["\u00E0", -1, 1],
        ["\u00E1", -1, 1],
        ["\u00E9", -1, 1],
        ["\u00EC", -1, 1],
        ["\u00ED", -1, 1],
        ["\u00EF", -1, 1],
        ["\u00F3", -1, 1]
    ];

    /** @const */ var /** Array<int> */ g_v = [17, 65, 16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 128, 129, 81, 6, 10];

    var /** number */ I_p2 = 0;
    var /** number */ I_p1 = 0;


    /** @return {boolean} */
    function r_mark_regions() {
        I_p1 = base.limit;
        I_p2 = base.limit;
        /** @const */ var /** number */ v_1 = base.cursor;
        lab0: {
            if (!base.go_out_grouping(g_v, 97, 252))
            {
                break lab0;
            }
            base.cursor++;
            if (!base.go_in_grouping(g_v, 97, 252))
            {
                break lab0;
            }
            base.cursor++;
            I_p1 = base.cursor;
            if (!base.go_out_grouping(g_v, 97, 252))
            {
                break lab0;
            }
            base.cursor++;
            if (!base.go_in_grouping(g_v, 97, 252))
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
    function r_cleaning() {
        var /** number */ among_var;
        while(true)
        {
            /** @const */ var /** number */ v_1 = base.cursor;
            lab0: {
                base.bra = base.cursor;
                among_var = base.find_among(a_0);
                base.ket = base.cursor;
                switch (among_var) {
                    case 1:
                        if (!base.slice_from("a"))
                        {
                            return false;
                        }
                        break;
                    case 2:
                        if (!base.slice_from("e"))
                        {
                            return false;
                        }
                        break;
                    case 3:
                        if (!base.slice_from("i"))
                        {
                            return false;
                        }
                        break;
                    case 4:
                        if (!base.slice_from("o"))
                        {
                            return false;
                        }
                        break;
                    case 5:
                        if (!base.slice_from("u"))
                        {
                            return false;
                        }
                        break;
                    case 6:
                        if (!base.slice_from("."))
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
    function r_R1() {
        return I_p1 <= base.cursor;
    };

    /** @return {boolean} */
    function r_R2() {
        return I_p2 <= base.cursor;
    };

    /** @return {boolean} */
    function r_attached_pronoun() {
        base.ket = base.cursor;
        if (base.find_among_b(a_1) == 0)
        {
            return false;
        }
        base.bra = base.cursor;
        if (!r_R1())
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
    function r_standard_suffix() {
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
                if (!base.slice_from("ic"))
                {
                    return false;
                }
                break;
            case 5:
                if (!r_R1())
                {
                    return false;
                }
                if (!base.slice_from("c"))
                {
                    return false;
                }
                break;
        }
        return true;
    };

    /** @return {boolean} */
    function r_verb_suffix() {
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
    function r_residual_suffix() {
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
                if (!r_R1())
                {
                    return false;
                }
                if (!base.slice_from("ic"))
                {
                    return false;
                }
                break;
        }
        return true;
    };

    this.stem = /** @return {boolean} */ function() {
        r_mark_regions();
        base.limit_backward = base.cursor; base.cursor = base.limit;
        /** @const */ var /** number */ v_1 = base.limit - base.cursor;
        r_attached_pronoun();
        base.cursor = base.limit - v_1;
        /** @const */ var /** number */ v_2 = base.limit - base.cursor;
        lab0: {
            lab1: {
                /** @const */ var /** number */ v_3 = base.limit - base.cursor;
                lab2: {
                    if (!r_standard_suffix())
                    {
                        break lab2;
                    }
                    break lab1;
                }
                base.cursor = base.limit - v_3;
                if (!r_verb_suffix())
                {
                    break lab0;
                }
            }
        }
        base.cursor = base.limit - v_2;
        /** @const */ var /** number */ v_4 = base.limit - base.cursor;
        r_residual_suffix();
        base.cursor = base.limit - v_4;
        base.cursor = base.limit_backward;
        /** @const */ var /** number */ v_5 = base.cursor;
        r_cleaning();
        base.cursor = v_5;
        return true;
    };

    /**@return{string}*/
    this['stemWord'] = function(/**string*/word) {
        base.setCurrent(word);
        this.stem();
        return base.getCurrent();
    };
};
