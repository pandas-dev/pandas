"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[7817],{

/***/ 77817:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   forth: () => (/* binding */ forth)
/* harmony export */ });
function toWordList(words) {
  var ret = [];
  words.split(' ').forEach(function(e){
    ret.push({name: e});
  });
  return ret;
}

var coreWordList = toWordList(
  'INVERT AND OR XOR\
 2* 2/ LSHIFT RSHIFT\
 0= = 0< < > U< MIN MAX\
 2DROP 2DUP 2OVER 2SWAP ?DUP DEPTH DROP DUP OVER ROT SWAP\
 >R R> R@\
 + - 1+ 1- ABS NEGATE\
 S>D * M* UM*\
 FM/MOD SM/REM UM/MOD */ */MOD / /MOD MOD\
 HERE , @ ! CELL+ CELLS C, C@ C! CHARS 2@ 2!\
 ALIGN ALIGNED +! ALLOT\
 CHAR [CHAR] [ ] BL\
 FIND EXECUTE IMMEDIATE COUNT LITERAL STATE\
 ; DOES> >BODY\
 EVALUATE\
 SOURCE >IN\
 <# # #S #> HOLD SIGN BASE >NUMBER HEX DECIMAL\
 FILL MOVE\
 . CR EMIT SPACE SPACES TYPE U. .R U.R\
 ACCEPT\
 TRUE FALSE\
 <> U> 0<> 0>\
 NIP TUCK ROLL PICK\
 2>R 2R@ 2R>\
 WITHIN UNUSED MARKER\
 I J\
 TO\
 COMPILE, [COMPILE]\
 SAVE-INPUT RESTORE-INPUT\
 PAD ERASE\
 2LITERAL DNEGATE\
 D- D+ D0< D0= D2* D2/ D< D= DMAX DMIN D>S DABS\
 M+ M*/ D. D.R 2ROT DU<\
 CATCH THROW\
 FREE RESIZE ALLOCATE\
 CS-PICK CS-ROLL\
 GET-CURRENT SET-CURRENT FORTH-WORDLIST GET-ORDER SET-ORDER\
 PREVIOUS SEARCH-WORDLIST WORDLIST FIND ALSO ONLY FORTH DEFINITIONS ORDER\
 -TRAILING /STRING SEARCH COMPARE CMOVE CMOVE> BLANK SLITERAL');

var immediateWordList = toWordList('IF ELSE THEN BEGIN WHILE REPEAT UNTIL RECURSE [IF] [ELSE] [THEN] ?DO DO LOOP +LOOP UNLOOP LEAVE EXIT AGAIN CASE OF ENDOF ENDCASE');

function searchWordList (wordList, word) {
  var i;
  for (i = wordList.length - 1; i >= 0; i--) {
    if (wordList[i].name === word.toUpperCase()) {
      return wordList[i];
    }
  }
  return undefined;
}
const forth = {
  name: "forth",
  startState: function() {
    return {
      state: '',
      base: 10,
      coreWordList: coreWordList,
      immediateWordList: immediateWordList,
      wordList: []
    };
  },
  token: function (stream, stt) {
    var mat;
    if (stream.eatSpace()) {
      return null;
    }
    if (stt.state === '') { // interpretation
      if (stream.match(/^(\]|:NONAME)(\s|$)/i)) {
        stt.state = ' compilation';
        return 'builtin';
      }
      mat = stream.match(/^(\:)\s+(\S+)(\s|$)+/);
      if (mat) {
        stt.wordList.push({name: mat[2].toUpperCase()});
        stt.state = ' compilation';
        return 'def';
      }
      mat = stream.match(/^(VARIABLE|2VARIABLE|CONSTANT|2CONSTANT|CREATE|POSTPONE|VALUE|WORD)\s+(\S+)(\s|$)+/i);
      if (mat) {
        stt.wordList.push({name: mat[2].toUpperCase()});
        return 'def';
      }
      mat = stream.match(/^(\'|\[\'\])\s+(\S+)(\s|$)+/);
      if (mat) {
        return 'builtin'
      }
    } else { // compilation
      // ; [
      if (stream.match(/^(\;|\[)(\s)/)) {
        stt.state = '';
        stream.backUp(1);
        return 'builtin';
      }
      if (stream.match(/^(\;|\[)($)/)) {
        stt.state = '';
        return 'builtin';
      }
      if (stream.match(/^(POSTPONE)\s+\S+(\s|$)+/)) {
        return 'builtin';
      }
    }

    // dynamic wordlist
    mat = stream.match(/^(\S+)(\s+|$)/);
    if (mat) {
      if (searchWordList(stt.wordList, mat[1]) !== undefined) {
        return 'variable';
      }

      // comments
      if (mat[1] === '\\') {
        stream.skipToEnd();
        return 'comment';
      }

      // core words
      if (searchWordList(stt.coreWordList, mat[1]) !== undefined) {
        return 'builtin';
      }
      if (searchWordList(stt.immediateWordList, mat[1]) !== undefined) {
        return 'keyword';
      }

      if (mat[1] === '(') {
        stream.eatWhile(function (s) { return s !== ')'; });
        stream.eat(')');
        return 'comment';
      }

      // // strings
      if (mat[1] === '.(') {
        stream.eatWhile(function (s) { return s !== ')'; });
        stream.eat(')');
        return 'string';
      }
      if (mat[1] === 'S"' || mat[1] === '."' || mat[1] === 'C"') {
        stream.eatWhile(function (s) { return s !== '"'; });
        stream.eat('"');
        return 'string';
      }

      // numbers
      if (mat[1] - 0xfffffffff) {
        return 'number';
      }
      // if (mat[1].match(/^[-+]?[0-9]+\.[0-9]*/)) {
      //     return 'number';
      // }

      return 'atom';
    }
  }
};


/***/ })

}]);
//# sourceMappingURL=7817.74b742c39300a07a9efa.js.map?v=74b742c39300a07a9efa