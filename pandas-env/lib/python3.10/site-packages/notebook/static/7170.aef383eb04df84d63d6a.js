"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[7170],{

/***/ 27170:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   shell: () => (/* binding */ shell)
/* harmony export */ });
var words = {};
function define(style, dict) {
  for(var i = 0; i < dict.length; i++) {
    words[dict[i]] = style;
  }
};

var commonAtoms = ["true", "false"];
var commonKeywords = ["if", "then", "do", "else", "elif", "while", "until", "for", "in", "esac", "fi",
                      "fin", "fil", "done", "exit", "set", "unset", "export", "function"];
var commonCommands = ["ab", "awk", "bash", "beep", "cat", "cc", "cd", "chown", "chmod", "chroot", "clear",
                      "cp", "curl", "cut", "diff", "echo", "find", "gawk", "gcc", "get", "git", "grep", "hg", "kill", "killall",
                      "ln", "ls", "make", "mkdir", "openssl", "mv", "nc", "nl", "node", "npm", "ping", "ps", "restart", "rm",
                      "rmdir", "sed", "service", "sh", "shopt", "shred", "source", "sort", "sleep", "ssh", "start", "stop",
                      "su", "sudo", "svn", "tee", "telnet", "top", "touch", "vi", "vim", "wall", "wc", "wget", "who", "write",
                      "yes", "zsh"];

define('atom', commonAtoms);
define('keyword', commonKeywords);
define('builtin', commonCommands);

function tokenBase(stream, state) {
  if (stream.eatSpace()) return null;

  var sol = stream.sol();
  var ch = stream.next();

  if (ch === '\\') {
    stream.next();
    return null;
  }
  if (ch === '\'' || ch === '"' || ch === '`') {
    state.tokens.unshift(tokenString(ch, ch === "`" ? "quote" : "string"));
    return tokenize(stream, state);
  }
  if (ch === '#') {
    if (sol && stream.eat('!')) {
      stream.skipToEnd();
      return 'meta'; // 'comment'?
    }
    stream.skipToEnd();
    return 'comment';
  }
  if (ch === '$') {
    state.tokens.unshift(tokenDollar);
    return tokenize(stream, state);
  }
  if (ch === '+' || ch === '=') {
    return 'operator';
  }
  if (ch === '-') {
    stream.eat('-');
    stream.eatWhile(/\w/);
    return 'attribute';
  }
  if (ch == "<") {
    if (stream.match("<<")) return "operator"
    var heredoc = stream.match(/^<-?\s*(?:['"]([^'"]*)['"]|([^'"\s]*))/)
    if (heredoc) {
      state.tokens.unshift(tokenHeredoc(heredoc[1] || heredoc[2]))
      return 'string.special'
    }
  }
  if (/\d/.test(ch)) {
    stream.eatWhile(/\d/);
    if(stream.eol() || !/\w/.test(stream.peek())) {
      return 'number';
    }
  }
  stream.eatWhile(/[\w-]/);
  var cur = stream.current();
  if (stream.peek() === '=' && /\w+/.test(cur)) return 'def';
  return words.hasOwnProperty(cur) ? words[cur] : null;
}

function tokenString(quote, style) {
  var close = quote == "(" ? ")" : quote == "{" ? "}" : quote
  return function(stream, state) {
    var next, escaped = false;
    while ((next = stream.next()) != null) {
      if (next === close && !escaped) {
        state.tokens.shift();
        break;
      } else if (next === '$' && !escaped && quote !== "'" && stream.peek() != close) {
        escaped = true;
        stream.backUp(1);
        state.tokens.unshift(tokenDollar);
        break;
      } else if (!escaped && quote !== close && next === quote) {
        state.tokens.unshift(tokenString(quote, style))
        return tokenize(stream, state)
      } else if (!escaped && /['"]/.test(next) && !/['"]/.test(quote)) {
        state.tokens.unshift(tokenStringStart(next, "string"));
        stream.backUp(1);
        break;
      }
      escaped = !escaped && next === '\\';
    }
    return style;
  };
};

function tokenStringStart(quote, style) {
  return function(stream, state) {
    state.tokens[0] = tokenString(quote, style)
    stream.next()
    return tokenize(stream, state)
  }
}

var tokenDollar = function(stream, state) {
  if (state.tokens.length > 1) stream.eat('$');
  var ch = stream.next()
  if (/['"({]/.test(ch)) {
    state.tokens[0] = tokenString(ch, ch == "(" ? "quote" : ch == "{" ? "def" : "string");
    return tokenize(stream, state);
  }
  if (!/\d/.test(ch)) stream.eatWhile(/\w/);
  state.tokens.shift();
  return 'def';
};

function tokenHeredoc(delim) {
  return function(stream, state) {
    if (stream.sol() && stream.string == delim) state.tokens.shift()
    stream.skipToEnd()
    return "string.special"
  }
}

function tokenize(stream, state) {
  return (state.tokens[0] || tokenBase) (stream, state);
};

const shell = {
  name: "shell",
  startState: function() {return {tokens:[]};},
  token: function(stream, state) {
    return tokenize(stream, state);
  },
  languageData: {
    autocomplete: commonAtoms.concat(commonKeywords, commonCommands),
    closeBrackets: {brackets: ["(", "[", "{", "'", '"', "`"]},
    commentTokens: {line: "#"}
  }
};


/***/ })

}]);
//# sourceMappingURL=7170.aef383eb04df84d63d6a.js.map?v=aef383eb04df84d63d6a