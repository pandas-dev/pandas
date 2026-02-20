"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[1088],{

/***/ 31088:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   q: () => (/* binding */ q)
/* harmony export */ });
var curPunc,
    keywords=buildRE(["abs","acos","aj","aj0","all","and","any","asc","asin","asof","atan","attr","avg","avgs","bin","by","ceiling","cols","cor","cos","count","cov","cross","csv","cut","delete","deltas","desc","dev","differ","distinct","div","do","each","ej","enlist","eval","except","exec","exit","exp","fby","fills","first","fkeys","flip","floor","from","get","getenv","group","gtime","hclose","hcount","hdel","hopen","hsym","iasc","idesc","if","ij","in","insert","inter","inv","key","keys","last","like","list","lj","load","log","lower","lsq","ltime","ltrim","mavg","max","maxs","mcount","md5","mdev","med","meta","min","mins","mmax","mmin","mmu","mod","msum","neg","next","not","null","or","over","parse","peach","pj","plist","prd","prds","prev","prior","rand","rank","ratios","raze","read0","read1","reciprocal","reverse","rload","rotate","rsave","rtrim","save","scan","select","set","setenv","show","signum","sin","sqrt","ss","ssr","string","sublist","sum","sums","sv","system","tables","tan","til","trim","txf","type","uj","ungroup","union","update","upper","upsert","value","var","view","views","vs","wavg","where","where","while","within","wj","wj1","wsum","xasc","xbar","xcol","xcols","xdesc","xexp","xgroup","xkey","xlog","xprev","xrank"]),
    E=/[|/&^!+:\\\-*%$=~#;@><,?_\'\"\[\(\]\)\s{}]/;
function buildRE(w){return new RegExp("^("+w.join("|")+")$");}
function tokenBase(stream,state){
  var sol=stream.sol(),c=stream.next();
  curPunc=null;
  if(sol)
    if(c=="/")
      return(state.tokenize=tokenLineComment)(stream,state);
  else if(c=="\\"){
    if(stream.eol()||/\s/.test(stream.peek()))
      return stream.skipToEnd(),/^\\\s*$/.test(stream.current())?(state.tokenize=tokenCommentToEOF)(stream):state.tokenize=tokenBase,"comment";
    else
      return state.tokenize=tokenBase,"builtin";
  }
  if(/\s/.test(c))
    return stream.peek()=="/"?(stream.skipToEnd(),"comment"):"null";
  if(c=='"')
    return(state.tokenize=tokenString)(stream,state);
  if(c=='`')
    return stream.eatWhile(/[A-Za-z\d_:\/.]/),"macroName";
  if(("."==c&&/\d/.test(stream.peek()))||/\d/.test(c)){
    var t=null;
    stream.backUp(1);
    if(stream.match(/^\d{4}\.\d{2}(m|\.\d{2}([DT](\d{2}(:\d{2}(:\d{2}(\.\d{1,9})?)?)?)?)?)/)
       || stream.match(/^\d+D(\d{2}(:\d{2}(:\d{2}(\.\d{1,9})?)?)?)/)
       || stream.match(/^\d{2}:\d{2}(:\d{2}(\.\d{1,9})?)?/)
       || stream.match(/^\d+[ptuv]{1}/))
      t="temporal";
    else if(stream.match(/^0[NwW]{1}/)
            || stream.match(/^0x[\da-fA-F]*/)
            || stream.match(/^[01]+[b]{1}/)
            || stream.match(/^\d+[chijn]{1}/)
            || stream.match(/-?\d*(\.\d*)?(e[+\-]?\d+)?(e|f)?/))
      t="number";
    return(t&&(!(c=stream.peek())||E.test(c)))?t:(stream.next(),"error");
  }
  if(/[A-Za-z]|\./.test(c))
    return stream.eatWhile(/[A-Za-z._\d]/),keywords.test(stream.current())?"keyword":"variable";
  if(/[|/&^!+:\\\-*%$=~#;@><\.,?_\']/.test(c))
    return null;
  if(/[{}\(\[\]\)]/.test(c))
    return null;
  return"error";
}
function tokenLineComment(stream,state){
  return stream.skipToEnd(),/^\/\s*$/.test(stream.current())?(state.tokenize=tokenBlockComment)(stream,state):(state.tokenize=tokenBase),"comment";
}
function tokenBlockComment(stream,state){
  var f=stream.sol()&&stream.peek()=="\\";
  stream.skipToEnd();
  if(f&&/^\\\s*$/.test(stream.current()))
    state.tokenize=tokenBase;
  return"comment";
}
function tokenCommentToEOF(stream){return stream.skipToEnd(),"comment";}
function tokenString(stream,state){
  var escaped=false,next,end=false;
  while((next=stream.next())){
    if(next=="\""&&!escaped){end=true;break;}
    escaped=!escaped&&next=="\\";
  }
  if(end)state.tokenize=tokenBase;
  return"string";
}
function pushContext(state,type,col){state.context={prev:state.context,indent:state.indent,col:col,type:type};}
function popContext(state){state.indent=state.context.indent;state.context=state.context.prev;}
const q = {
  name: "q",
  startState:function(){
    return{tokenize:tokenBase,
           context:null,
           indent:0,
           col:0};
  },
  token:function(stream,state){
    if(stream.sol()){
      if(state.context&&state.context.align==null)
        state.context.align=false;
      state.indent=stream.indentation();
    }
    //if (stream.eatSpace()) return null;
    var style=state.tokenize(stream,state);
    if(style!="comment"&&state.context&&state.context.align==null&&state.context.type!="pattern"){
      state.context.align=true;
    }
    if(curPunc=="(")pushContext(state,")",stream.column());
    else if(curPunc=="[")pushContext(state,"]",stream.column());
    else if(curPunc=="{")pushContext(state,"}",stream.column());
    else if(/[\]\}\)]/.test(curPunc)){
      while(state.context&&state.context.type=="pattern")popContext(state);
      if(state.context&&curPunc==state.context.type)popContext(state);
    }
    else if(curPunc=="."&&state.context&&state.context.type=="pattern")popContext(state);
    else if(/atom|string|variable/.test(style)&&state.context){
      if(/[\}\]]/.test(state.context.type))
        pushContext(state,"pattern",stream.column());
      else if(state.context.type=="pattern"&&!state.context.align){
        state.context.align=true;
        state.context.col=stream.column();
      }
    }
    return style;
  },
  indent:function(state,textAfter,cx){
    var firstChar=textAfter&&textAfter.charAt(0);
    var context=state.context;
    if(/[\]\}]/.test(firstChar))
      while (context&&context.type=="pattern")context=context.prev;
    var closing=context&&firstChar==context.type;
    if(!context)
      return 0;
    else if(context.type=="pattern")
      return context.col;
    else if(context.align)
      return context.col+(closing?0:1);
    else
      return context.indent+(closing?0:cx.unit);
  },
  languageData: {
    commentTokens: { line: "/" },
  },
};


/***/ })

}]);
//# sourceMappingURL=1088.f26c568e858d1f160276.js.map?v=f26c568e858d1f160276