"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[403],{

/***/ 70403:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   asterisk: () => (/* binding */ asterisk)
/* harmony export */ });
var atoms    = ["exten", "same", "include","ignorepat","switch"],
    dpcmd    = ["#include","#exec"],
    apps     = [
      "addqueuemember","adsiprog","aelsub","agentlogin","agentmonitoroutgoing","agi",
      "alarmreceiver","amd","answer","authenticate","background","backgrounddetect",
      "bridge","busy","callcompletioncancel","callcompletionrequest","celgenuserevent",
      "changemonitor","chanisavail","channelredirect","chanspy","clearhash","confbridge",
      "congestion","continuewhile","controlplayback","dahdiacceptr2call","dahdibarge",
      "dahdiras","dahdiscan","dahdisendcallreroutingfacility","dahdisendkeypadfacility",
      "datetime","dbdel","dbdeltree","deadagi","dial","dictate","directory","disa",
      "dumpchan","eagi","echo","endwhile","exec","execif","execiftime","exitwhile","extenspy",
      "externalivr","festival","flash","followme","forkcdr","getcpeid","gosub","gosubif",
      "goto","gotoif","gotoiftime","hangup","iax2provision","ices","importvar","incomplete",
      "ivrdemo","jabberjoin","jabberleave","jabbersend","jabbersendgroup","jabberstatus",
      "jack","log","macro","macroexclusive","macroexit","macroif","mailboxexists","meetme",
      "meetmeadmin","meetmechanneladmin","meetmecount","milliwatt","minivmaccmess","minivmdelete",
      "minivmgreet","minivmmwi","minivmnotify","minivmrecord","mixmonitor","monitor","morsecode",
      "mp3player","mset","musiconhold","nbscat","nocdr","noop","odbc","odbc","odbcfinish",
      "originate","ospauth","ospfinish","osplookup","ospnext","page","park","parkandannounce",
      "parkedcall","pausemonitor","pausequeuemember","pickup","pickupchan","playback","playtones",
      "privacymanager","proceeding","progress","queue","queuelog","raiseexception","read","readexten",
      "readfile","receivefax","receivefax","receivefax","record","removequeuemember",
      "resetcdr","retrydial","return","ringing","sayalpha","saycountedadj","saycountednoun",
      "saycountpl","saydigits","saynumber","sayphonetic","sayunixtime","senddtmf","sendfax",
      "sendfax","sendfax","sendimage","sendtext","sendurl","set","setamaflags",
      "setcallerpres","setmusiconhold","sipaddheader","sipdtmfmode","sipremoveheader","skel",
      "slastation","slatrunk","sms","softhangup","speechactivategrammar","speechbackground",
      "speechcreate","speechdeactivategrammar","speechdestroy","speechloadgrammar","speechprocessingsound",
      "speechstart","speechunloadgrammar","stackpop","startmusiconhold","stopmixmonitor","stopmonitor",
      "stopmusiconhold","stopplaytones","system","testclient","testserver","transfer","tryexec",
      "trysystem","unpausemonitor","unpausequeuemember","userevent","verbose","vmauthenticate",
      "vmsayname","voicemail","voicemailmain","wait","waitexten","waitfornoise","waitforring",
      "waitforsilence","waitmusiconhold","waituntil","while","zapateller"
    ];

function basicToken(stream,state){
  var cur = '';
  var ch = stream.next();
  // comment
  if (state.blockComment) {
    if (ch == "-" && stream.match("-;", true)) {
      state.blockComment = false;
    } else if (stream.skipTo("--;")) {
      stream.next();
      stream.next();
      stream.next();
      state.blockComment = false;
    } else {
      stream.skipToEnd();
    }
    return "comment";
  }
  if(ch == ";") {
    if (stream.match("--", true)) {
      if (!stream.match("-", false)) {  // Except ;--- is not a block comment
        state.blockComment = true;
        return "comment";
      }
    }
    stream.skipToEnd();
    return "comment";
  }
  // context
  if(ch == '[') {
    stream.skipTo(']');
    stream.eat(']');
    return "header";
  }
  // string
  if(ch == '"') {
    stream.skipTo('"');
    return "string";
  }
  if(ch == "'") {
    stream.skipTo("'");
    return "string.special";
  }
  // dialplan commands
  if(ch == '#') {
    stream.eatWhile(/\w/);
    cur = stream.current();
    if(dpcmd.indexOf(cur) !== -1) {
      stream.skipToEnd();
      return "strong";
    }
  }
  // application args
  if(ch == '$'){
    var ch1 = stream.peek();
    if(ch1 == '{'){
      stream.skipTo('}');
      stream.eat('}');
      return "variableName.special";
    }
  }
  // extension
  stream.eatWhile(/\w/);
  cur = stream.current();
  if(atoms.indexOf(cur) !== -1) {
    state.extenStart = true;
    switch(cur) {
    case 'same': state.extenSame = true; break;
    case 'include':
    case 'switch':
    case 'ignorepat':
      state.extenInclude = true;break;
    default:break;
    }
    return "atom";
  }
}

const asterisk = {
  name: "asterisk",
  startState: function() {
    return {
      blockComment: false,
      extenStart: false,
      extenSame:  false,
      extenInclude: false,
      extenExten: false,
      extenPriority: false,
      extenApplication: false
    };
  },
  token: function(stream, state) {

    var cur = '';
    if(stream.eatSpace()) return null;
    // extension started
    if(state.extenStart){
      stream.eatWhile(/[^\s]/);
      cur = stream.current();
      if(/^=>?$/.test(cur)){
        state.extenExten = true;
        state.extenStart = false;
        return "strong";
      } else {
        state.extenStart = false;
        stream.skipToEnd();
        return "error";
      }
    } else if(state.extenExten) {
      // set exten and priority
      state.extenExten = false;
      state.extenPriority = true;
      stream.eatWhile(/[^,]/);
      if(state.extenInclude) {
        stream.skipToEnd();
        state.extenPriority = false;
        state.extenInclude = false;
      }
      if(state.extenSame) {
        state.extenPriority = false;
        state.extenSame = false;
        state.extenApplication = true;
      }
      return "tag";
    } else if(state.extenPriority) {
      state.extenPriority = false;
      state.extenApplication = true;
      stream.next(); // get comma
      if(state.extenSame) return null;
      stream.eatWhile(/[^,]/);
      return "number";
    } else if(state.extenApplication) {
      stream.eatWhile(/,/);
      cur = stream.current();
      if(cur === ',') return null;
      stream.eatWhile(/\w/);
      cur = stream.current().toLowerCase();
      state.extenApplication = false;
      if(apps.indexOf(cur) !== -1){
        return "def";
      }
    } else{
      return basicToken(stream,state);
    }

    return null;
  },

  languageData: {
    commentTokens: {line: ";", block: {open: ";--", close: "--;"}}
  }
};


/***/ })

}]);
//# sourceMappingURL=403.270ca5cf44874182bd4d.js.map?v=270ca5cf44874182bd4d