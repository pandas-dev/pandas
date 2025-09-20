"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[3211],{

/***/ 3211:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   liveScript: () => (/* binding */ liveScript)
/* harmony export */ });
var tokenBase = function(stream, state) {
  var next_rule = state.next || "start";
  if (next_rule) {
    state.next = state.next;
    var nr = Rules[next_rule];
    if (nr.splice) {
      for (var i$ = 0; i$ < nr.length; ++i$) {
        var r = nr[i$];
        if (r.regex && stream.match(r.regex)) {
          state.next = r.next || state.next;
          return r.token;
        }
      }
      stream.next();
      return 'error';
    }
    if (stream.match(r = Rules[next_rule])) {
      if (r.regex && stream.match(r.regex)) {
        state.next = r.next;
        return r.token;
      } else {
        stream.next();
        return 'error';
      }
    }
  }
  stream.next();
  return 'error';
};

var identifier = '(?![\\d\\s])[$\\w\\xAA-\\uFFDC](?:(?!\\s)[$\\w\\xAA-\\uFFDC]|-[A-Za-z])*';
var indenter = RegExp('(?:[({[=:]|[-~]>|\\b(?:e(?:lse|xport)|d(?:o|efault)|t(?:ry|hen)|finally|import(?:\\s*all)?|const|var|let|new|catch(?:\\s*' + identifier + ')?))\\s*$');
var keywordend = '(?![$\\w]|-[A-Za-z]|\\s*:(?![:=]))';
var stringfill = {
  token: 'string',
  regex: '.+'
};
var Rules = {
  start: [
    {
      token: 'docComment',
      regex: '/\\*',
      next: 'comment'
    }, {
      token: 'comment',
      regex: '#.*'
    }, {
      token: 'keyword',
      regex: '(?:t(?:h(?:is|row|en)|ry|ypeof!?)|c(?:on(?:tinue|st)|a(?:se|tch)|lass)|i(?:n(?:stanceof)?|mp(?:ort(?:\\s+all)?|lements)|[fs])|d(?:e(?:fault|lete|bugger)|o)|f(?:or(?:\\s+own)?|inally|unction)|s(?:uper|witch)|e(?:lse|x(?:tends|port)|val)|a(?:nd|rguments)|n(?:ew|ot)|un(?:less|til)|w(?:hile|ith)|o[fr]|return|break|let|var|loop)' + keywordend
    }, {
      token: 'atom',
      regex: '(?:true|false|yes|no|on|off|null|void|undefined)' + keywordend
    }, {
      token: 'invalid',
      regex: '(?:p(?:ackage|r(?:ivate|otected)|ublic)|i(?:mplements|nterface)|enum|static|yield)' + keywordend
    }, {
      token: 'className.standard',
      regex: '(?:R(?:e(?:gExp|ferenceError)|angeError)|S(?:tring|yntaxError)|E(?:rror|valError)|Array|Boolean|Date|Function|Number|Object|TypeError|URIError)' + keywordend
    }, {
      token: 'variableName.function.standard',
      regex: '(?:is(?:NaN|Finite)|parse(?:Int|Float)|Math|JSON|(?:en|de)codeURI(?:Component)?)' + keywordend
    }, {
      token: 'variableName.standard',
      regex: '(?:t(?:hat|il|o)|f(?:rom|allthrough)|it|by|e)' + keywordend
    }, {
      token: 'variableName',
      regex: identifier + '\\s*:(?![:=])'
    }, {
      token: 'variableName',
      regex: identifier
    }, {
      token: 'operatorKeyword',
      regex: '(?:\\.{3}|\\s+\\?)'
    }, {
      token: 'keyword',
      regex: '(?:@+|::|\\.\\.)',
      next: 'key'
    }, {
      token: 'operatorKeyword',
      regex: '\\.\\s*',
      next: 'key'
    }, {
      token: 'string',
      regex: '\\\\\\S[^\\s,;)}\\]]*'
    }, {
      token: 'docString',
      regex: '\'\'\'',
      next: 'qdoc'
    }, {
      token: 'docString',
      regex: '"""',
      next: 'qqdoc'
    }, {
      token: 'string',
      regex: '\'',
      next: 'qstring'
    }, {
      token: 'string',
      regex: '"',
      next: 'qqstring'
    }, {
      token: 'string',
      regex: '`',
      next: 'js'
    }, {
      token: 'string',
      regex: '<\\[',
      next: 'words'
    }, {
      token: 'regexp',
      regex: '//',
      next: 'heregex'
    }, {
      token: 'regexp',
      regex: '\\/(?:[^[\\/\\n\\\\]*(?:(?:\\\\.|\\[[^\\]\\n\\\\]*(?:\\\\.[^\\]\\n\\\\]*)*\\])[^[\\/\\n\\\\]*)*)\\/[gimy$]{0,4}',
      next: 'key'
    }, {
      token: 'number',
      regex: '(?:0x[\\da-fA-F][\\da-fA-F_]*|(?:[2-9]|[12]\\d|3[0-6])r[\\da-zA-Z][\\da-zA-Z_]*|(?:\\d[\\d_]*(?:\\.\\d[\\d_]*)?|\\.\\d[\\d_]*)(?:e[+-]?\\d[\\d_]*)?[\\w$]*)'
    }, {
      token: 'paren',
      regex: '[({[]'
    }, {
      token: 'paren',
      regex: '[)}\\]]',
      next: 'key'
    }, {
      token: 'operatorKeyword',
      regex: '\\S+'
    }, {
      token: 'content',
      regex: '\\s+'
    }
  ],
  heregex: [
    {
      token: 'regexp',
      regex: '.*?//[gimy$?]{0,4}',
      next: 'start'
    }, {
      token: 'regexp',
      regex: '\\s*#{'
    }, {
      token: 'comment',
      regex: '\\s+(?:#.*)?'
    }, {
      token: 'regexp',
      regex: '\\S+'
    }
  ],
  key: [
    {
      token: 'operatorKeyword',
      regex: '[.?@!]+'
    }, {
      token: 'variableName',
      regex: identifier,
      next: 'start'
    }, {
      token: 'content',
      regex: '',
      next: 'start'
    }
  ],
  comment: [
    {
      token: 'docComment',
      regex: '.*?\\*/',
      next: 'start'
    }, {
      token: 'docComment',
      regex: '.+'
    }
  ],
  qdoc: [
    {
      token: 'string',
      regex: ".*?'''",
      next: 'key'
    }, stringfill
  ],
  qqdoc: [
    {
      token: 'string',
      regex: '.*?"""',
      next: 'key'
    }, stringfill
  ],
  qstring: [
    {
      token: 'string',
      regex: '[^\\\\\']*(?:\\\\.[^\\\\\']*)*\'',
      next: 'key'
    }, stringfill
  ],
  qqstring: [
    {
      token: 'string',
      regex: '[^\\\\"]*(?:\\\\.[^\\\\"]*)*"',
      next: 'key'
    }, stringfill
  ],
  js: [
    {
      token: 'string',
      regex: '[^\\\\`]*(?:\\\\.[^\\\\`]*)*`',
      next: 'key'
    }, stringfill
  ],
  words: [
    {
      token: 'string',
      regex: '.*?\\]>',
      next: 'key'
    }, stringfill
  ]
};
for (var idx in Rules) {
  var r = Rules[idx];
  if (r.splice) {
    for (var i = 0, len = r.length; i < len; ++i) {
      var rr = r[i];
      if (typeof rr.regex === 'string') {
        Rules[idx][i].regex = new RegExp('^' + rr.regex);
      }
    }
  } else if (typeof rr.regex === 'string') {
    Rules[idx].regex = new RegExp('^' + r.regex);
  }
}

const liveScript = {
  name: "livescript",
  startState: function(){
    return {
      next: 'start',
      lastToken: {style: null, indent: 0, content: ""}
    };
  },
  token: function(stream, state){
    while (stream.pos == stream.start)
      var style = tokenBase(stream, state);
    state.lastToken = {
      style: style,
      indent: stream.indentation(),
      content: stream.current()
    };
    return style.replace(/\./g, ' ');
  },
  indent: function(state){
    var indentation = state.lastToken.indent;
    if (state.lastToken.content.match(indenter)) {
      indentation += 2;
    }
    return indentation;
  }
};


/***/ })

}]);
//# sourceMappingURL=3211.2e93fd406e5c4e53774f.js.map?v=2e93fd406e5c4e53774f