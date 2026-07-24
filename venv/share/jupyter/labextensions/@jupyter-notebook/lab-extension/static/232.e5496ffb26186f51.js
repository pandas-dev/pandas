"use strict";(self.rspackChunk_jupyter_notebook_lab_extension=self.rspackChunk_jupyter_notebook_lab_extension||[]).push([[232],{4296(e,t,n){n.d(t,{A:()=>c});var r=n(6758),o=n.n(r),a=n(935),i=n.n(a)()(o());i.push([e.id,`/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

.jp-InterfaceSwitcher {
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
}

.jp-InterfaceSwitcher .lm-MenuBar-itemIcon svg {
  vertical-align: sub;
}

.jp-nb-interface-switcher-button > .jp-ToolbarButtonComponent::part(content) {
  flex-direction: row-reverse;
}

.jp-nb-interface-switcher-button > .jp-ToolbarButtonComponent > svg {
  padding-left: 3px;
}
`,""]);let c=i},935(e){e.exports=function(e){var t=[];return t.toString=function(){return this.map(function(t){var n="",r=void 0!==t[5];return t[4]&&(n+="@supports (".concat(t[4],") {")),t[2]&&(n+="@media ".concat(t[2]," {")),r&&(n+="@layer".concat(t[5].length>0?" ".concat(t[5]):""," {")),n+=e(t),r&&(n+="}"),t[2]&&(n+="}"),t[4]&&(n+="}"),n}).join("")},t.i=function(e,n,r,o,a){"string"==typeof e&&(e=[[null,e,void 0]]);var i={};if(r)for(var c=0;c<this.length;c++){var s=this[c][0];null!=s&&(i[s]=!0)}for(var u=0;u<e.length;u++){var p=[].concat(e[u]);r&&i[p[0]]||(void 0!==a&&(void 0===p[5]||(p[1]="@layer".concat(p[5].length>0?" ".concat(p[5]):""," {").concat(p[1],"}")),p[5]=a),n&&(p[2]&&(p[1]="@media ".concat(p[2]," {").concat(p[1],"}")),p[2]=n),o&&(p[4]?(p[1]="@supports (".concat(p[4],") {").concat(p[1],"}"),p[4]=o):p[4]="".concat(o)),t.push(p))}},t}},6758(e){e.exports=function(e){return e[1]}},2591(e){var t=[];function n(e){for(var n=-1,r=0;r<t.length;r++)if(t[r].identifier===e){n=r;break}return n}function r(e,r){for(var o={},a=[],i=0;i<e.length;i++){var c=e[i],s=r.base?c[0]+r.base:c[0],u=o[s]||0,p="".concat(s," ").concat(u);o[s]=u+1;var l=n(p),f={css:c[1],media:c[2],sourceMap:c[3],supports:c[4],layer:c[5]};if(-1!==l)t[l].references++,t[l].updater(f);else{var d=function(e,t){var n=t.domAPI(t);return n.update(e),function(t){t?(t.css!==e.css||t.media!==e.media||t.sourceMap!==e.sourceMap||t.supports!==e.supports||t.layer!==e.layer)&&n.update(e=t):n.remove()}}(f,r);r.byIndex=i,t.splice(i,0,{identifier:p,updater:d,references:1})}a.push(p)}return a}e.exports=function(e,o){var a=r(e=e||[],o=o||{});return function(e){e=e||[];for(var i=0;i<a.length;i++){var c=n(a[i]);t[c].references--}for(var s=r(e,o),u=0;u<a.length;u++){var p=n(a[u]);0===t[p].references&&(t[p].updater(),t.splice(p,1))}a=s}}},8128(e){var t={};e.exports=function(e,n){var r=function(e){if(void 0===t[e]){var n=document.querySelector(e);if(window.HTMLIFrameElement&&n instanceof window.HTMLIFrameElement)try{n=n.contentDocument.head}catch(e){n=null}t[e]=n}return t[e]}(e);if(!r)throw Error("Couldn't find a style target. This probably means that the value for the 'insert' parameter is invalid.");r.appendChild(n)}},3051(e){e.exports=function(e){var t=document.createElement("style");return e.setAttributes(t,e.attributes),e.insert(t,e.options),t}},855(e,t,n){e.exports=function(e){var t=n.nc;t&&e.setAttribute("nonce",t)}},1740(e){e.exports=function(e){if("u"<typeof document)return{update:function(){},remove:function(){}};var t=e.insertStyleElement(e);return{update:function(n){var r,o,a;r="",n.supports&&(r+="@supports (".concat(n.supports,") {")),n.media&&(r+="@media ".concat(n.media," {")),(o=void 0!==n.layer)&&(r+="@layer".concat(n.layer.length>0?" ".concat(n.layer):""," {")),r+=n.css,o&&(r+="}"),n.media&&(r+="}"),n.supports&&(r+="}"),(a=n.sourceMap)&&"u">typeof btoa&&(r+="\n/*# sourceMappingURL=data:application/json;base64,".concat(btoa(unescape(encodeURIComponent(JSON.stringify(a))))," */")),e.styleTagTransform(r,t,e.options)},remove:function(){var e;null===(e=t).parentNode||e.parentNode.removeChild(e)}}}},3656(e){e.exports=function(e,t){if(t.styleSheet)t.styleSheet.cssText=e;else{for(;t.firstChild;)t.removeChild(t.firstChild);t.appendChild(document.createTextNode(e))}}},8579(e,t,n){var r=n(2591),o=n.n(r),a=n(1740),i=n.n(a),c=n(8128),s=n.n(c),u=n(855),p=n.n(u),l=n(3051),f=n.n(l),d=n(3656),v=n.n(d),m=n(4296),h={};h.styleTagTransform=v(),h.setAttributes=p(),h.insert=s().bind(null,"head"),h.domAPI=i(),h.insertStyleElement=f(),o()(m.A,h),m.A&&m.A.locals&&m.A.locals}}]);