(self.rspackChunk_JUPYTERLAB_CORE_OUTPUT=self.rspackChunk_JUPYTERLAB_CORE_OUTPUT||[]).push([[2343],{71439(){},80458(e,t,o){"use strict";o.r(t);var n=o(93673),a=o(83044);async function r(e,t){await new Promise((t,o)=>{let n=document.createElement("script");n.onerror=o,n.onload=t,n.async=!0,document.head.appendChild(n),n.src=e}),await o.I("default");let n=window._JUPYTERLAB[t];await n.init(o.S.default)}async function i(e,t){try{let o=(await window._JUPYTERLAB[e].get(t))();return o.__scope__=e,o}catch(o){throw console.warn(`Failed to create module: package: ${e}; module: ${t}`),o}}async function l(){let e=[o(38537),o(38789),o(52331),o(18655)],t=await Promise.all(e),l=[o(70241),o(363),o(3481),o(83989),o(73941),o(81121),o(13093),o(41387).default.filter(({id:e})=>["@jupyterlab/application-extension:context-menu","@jupyterlab/application-extension:faviconbusy","@jupyterlab/application-extension:router","@jupyterlab/application-extension:top-bar","@jupyterlab/application-extension:top-spacer"].includes(e)),o(19689).default.filter(({id:e})=>["@jupyterlab/apputils-extension:kernels-settings","@jupyterlab/apputils-extension:palette","@jupyterlab/apputils-extension:notification","@jupyterlab/apputils-extension:sanitizer","@jupyterlab/apputils-extension:sessionDialogs","@jupyterlab/apputils-extension:settings","@jupyterlab/apputils-extension:state","@jupyterlab/apputils-extension:themes","@jupyterlab/apputils-extension:themes-palette-menu","@jupyterlab/apputils-extension:toolbar-registry","@jupyterlab/apputils-extension:utilityCommands"].includes(e)),o(53739),o(11057),o(40805).default.filter(({id:e})=>["@jupyterlab/completer-extension:base-service","@jupyterlab/completer-extension:inline-completer","@jupyterlab/completer-extension:inline-completer-factory","@jupyterlab/completer-extension:inline-history","@jupyterlab/completer-extension:manager"].includes(e)),o(30081).default.filter(({id:e})=>["@jupyterlab/console-extension:cell-executor","@jupyterlab/console-extension:completer","@jupyterlab/console-extension:factory","@jupyterlab/console-extension:foreign","@jupyterlab/console-extension:tracker"].includes(e)),o(27787),o(9881).default.filter(({id:e})=>["@jupyterlab/docmanager-extension:plugin","@jupyterlab/docmanager-extension:download","@jupyterlab/docmanager-extension:contexts","@jupyterlab/docmanager-extension:manager"].includes(e)),o(96581).default.filter(({id:e})=>["@jupyterlab/documentsearch-extension:plugin"].includes(e)),o(34803).default.filter(({id:e})=>["@jupyterlab/filebrowser-extension:factory","@jupyterlab/filebrowser-extension:default-file-browser"].includes(e)),o(18433).default.filter(({id:e})=>["@jupyterlab/fileeditor-extension:plugin","@jupyterlab/fileeditor-extension:widget-factory"].includes(e)),o(79713).default.filter(({id:e})=>["@jupyterlab/help-extension:resources"].includes(e)),o(22705),o(6925),o(79297),o(13357),o(17457).default.filter(({id:e})=>["@jupyterlab/mainmenu-extension:plugin"].includes(e)),o(12993),o(1297),o(35437),o(94885).default.filter(({id:e})=>["@jupyterlab/notebook-extension:cell-executor","@jupyterlab/notebook-extension:code-console","@jupyterlab/notebook-extension:export","@jupyterlab/notebook-extension:factory","@jupyterlab/notebook-extension:tracker","@jupyterlab/notebook-extension:widget-factory"].includes(e)),o(16747),o(78201),o(9762),o(25441),o(29611),o(50769),o(18121),o(45413),o(78273),o(80401)];switch(`/${n.PageConfig.getOption("notebookPage")}`){case"/tree":l=l.concat([o(41387).default.filter(({id:e})=>["@jupyterlab/application-extension:commands"].includes(e)),o(37565),o(61289),o(34803).default.filter(({id:e})=>["@jupyterlab/filebrowser-extension:browser","@jupyterlab/filebrowser-extension:create-new-language-file","@jupyterlab/filebrowser-extension:download","@jupyterlab/filebrowser-extension:file-upload-status","@jupyterlab/filebrowser-extension:open-with","@jupyterlab/filebrowser-extension:search","@jupyterlab/filebrowser-extension:share-file"].includes(e)),o(20021),o(8869).default.filter(({id:e})=>["@jupyterlab/running-extension:plugin"].includes(e)),o(83481)]);break;case"/notebooks":l=l.concat([o(61729),o(37565),o(42025).default.filter(({id:e})=>["@jupyterlab/debugger-extension:completions","@jupyterlab/debugger-extension:config","@jupyterlab/debugger-extension:debug-console","@jupyterlab/debugger-extension:main","@jupyterlab/debugger-extension:notebooks","@jupyterlab/debugger-extension:service","@jupyterlab/debugger-extension:sidebar","@jupyterlab/debugger-extension:sources","@jupyterlab/debugger-extension:display-registry"].includes(e)),o(24449),o(92673),o(94885).default.filter(({id:e})=>["@jupyterlab/notebook-extension:active-cell-tool","@jupyterlab/notebook-extension:completer","@jupyterlab/notebook-extension:copy-output","@jupyterlab/notebook-extension:metadata-editor","@jupyterlab/notebook-extension:search","@jupyterlab/notebook-extension:toc","@jupyterlab/notebook-extension:tools","@jupyterlab/notebook-extension:update-raw-mimetype"].includes(e)),o(37225).default.filter(({id:e})=>["@jupyterlab/toc-extension:registry","@jupyterlab/toc-extension:tracker"].includes(e)),o(93617).default.filter(({id:e})=>["@jupyterlab/tooltip-extension:manager","@jupyterlab/tooltip-extension:notebooks"].includes(e))]);break;case"/consoles":l=l.concat([o(93617).default.filter(({id:e})=>["@jupyterlab/tooltip-extension:manager","@jupyterlab/tooltip-extension:consoles"].includes(e))]);break;case"/edit":l=l.concat([o(18433).default.filter(({id:e})=>["@jupyterlab/fileeditor-extension:completer","@jupyterlab/fileeditor-extension:search"].includes(e)),o(72097)])}let p=[],s=[];function*u(e){let t;for(let o of Array.isArray(t=Object.prototype.hasOwnProperty.call(e,"__esModule")?e.default:e)?t:[t]){let t=n.PageConfig.Extension.isDisabled(o.id);if(s.push({id:o.id,description:o.description,requires:o.requires??[],optional:o.optional??[],provides:o.provides??null,autoStart:o.autoStart,enabled:!t,extension:e.__scope__}),t){p.push(o.id);continue}yield o}}let b=JSON.parse(n.PageConfig.getOption("federated_extensions")),d=[],c=[],j=[],y=[];(await Promise.allSettled(b.map(async e=>(await r(`${n.URLExt.join(n.PageConfig.getOption("fullLabextensionsUrl"),e.name,e.load)}`,e.name),e)))).forEach(e=>{if("rejected"===e.status)return void console.error(e.reason);let t=e.value;t.extension&&c.push(i(t.name,t.extension)),t.mimeExtension&&j.push(i(t.name,t.mimeExtension)),t.style&&!n.PageConfig.Extension.isDisabled(t.name)&&y.push(i(t.name,t.style))}),(await Promise.all(l)).forEach(e=>{for(let t of u(e))d.push(t)}),(await Promise.allSettled(c)).forEach(e=>{if("fulfilled"===e.status)for(let t of u(e.value))d.push(t);else console.error(e.reason)}),(await Promise.allSettled(j)).forEach(e=>{if("fulfilled"===e.status)for(let o of u(e.value))t.push(o);else console.error(e.reason)}),(await Promise.allSettled(y)).filter(({status:e})=>"rejected"===e).forEach(({reason:e})=>{console.error(e)}),n.PageConfig.setOption("allPlugins",'{"/":{"@jupyter-notebook/application-extension":true,"@jupyter-notebook/console-extension":true,"@jupyter-notebook/docmanager-extension":true,"@jupyter-notebook/documentsearch-extension":true,"@jupyter-notebook/help-extension":true,"@jupyter-notebook/notebook-extension":true,"@jupyter-notebook/terminal-extension":true,"@jupyterlab/application-extension":["@jupyterlab/application-extension:context-menu","@jupyterlab/application-extension:faviconbusy","@jupyterlab/application-extension:router","@jupyterlab/application-extension:top-bar","@jupyterlab/application-extension:top-spacer"],"@jupyterlab/apputils-extension":["@jupyterlab/apputils-extension:kernels-settings","@jupyterlab/apputils-extension:palette","@jupyterlab/apputils-extension:notification","@jupyterlab/apputils-extension:sanitizer","@jupyterlab/apputils-extension:sessionDialogs","@jupyterlab/apputils-extension:settings","@jupyterlab/apputils-extension:state","@jupyterlab/apputils-extension:themes","@jupyterlab/apputils-extension:themes-palette-menu","@jupyterlab/apputils-extension:toolbar-registry","@jupyterlab/apputils-extension:utilityCommands"],"@jupyterlab/audio-extension":true,"@jupyterlab/codemirror-extension":true,"@jupyterlab/completer-extension":["@jupyterlab/completer-extension:base-service","@jupyterlab/completer-extension:inline-completer","@jupyterlab/completer-extension:inline-completer-factory","@jupyterlab/completer-extension:inline-history","@jupyterlab/completer-extension:manager"],"@jupyterlab/console-extension":["@jupyterlab/console-extension:cell-executor","@jupyterlab/console-extension:completer","@jupyterlab/console-extension:factory","@jupyterlab/console-extension:foreign","@jupyterlab/console-extension:tracker"],"@jupyterlab/csvviewer-extension":true,"@jupyterlab/docmanager-extension":["@jupyterlab/docmanager-extension:plugin","@jupyterlab/docmanager-extension:download","@jupyterlab/docmanager-extension:contexts","@jupyterlab/docmanager-extension:manager"],"@jupyterlab/documentsearch-extension":["@jupyterlab/documentsearch-extension:plugin"],"@jupyterlab/filebrowser-extension":["@jupyterlab/filebrowser-extension:factory","@jupyterlab/filebrowser-extension:default-file-browser"],"@jupyterlab/fileeditor-extension":["@jupyterlab/fileeditor-extension:plugin","@jupyterlab/fileeditor-extension:widget-factory"],"@jupyterlab/help-extension":["@jupyterlab/help-extension:resources"],"@jupyterlab/htmlviewer-extension":true,"@jupyterlab/hub-extension":true,"@jupyterlab/imageviewer-extension":true,"@jupyterlab/lsp-extension":true,"@jupyterlab/mainmenu-extension":["@jupyterlab/mainmenu-extension:plugin"],"@jupyterlab/markedparser-extension":true,"@jupyterlab/mathjax-extension":true,"@jupyterlab/mermaid-extension":true,"@jupyterlab/notebook-extension":["@jupyterlab/notebook-extension:cell-executor","@jupyterlab/notebook-extension:code-console","@jupyterlab/notebook-extension:export","@jupyterlab/notebook-extension:factory","@jupyterlab/notebook-extension:tracker","@jupyterlab/notebook-extension:widget-factory"],"@jupyterlab/pluginmanager-extension":true,"@jupyterlab/services-extension":true,"@jupyterlab/shortcuts-extension":true,"@jupyterlab/terminal-extension":true,"@jupyterlab/theme-light-extension":true,"@jupyterlab/theme-dark-extension":true,"@jupyterlab/theme-dark-high-contrast-extension":true,"@jupyterlab/translation-extension":true,"@jupyterlab/ui-components-extension":true,"@jupyterlab/video-extension":true},"/tree":{"@jupyterlab/application-extension":["@jupyterlab/application-extension:commands"],"@jupyterlab/cell-toolbar-extension":true,"@jupyterlab/extensionmanager-extension":true,"@jupyterlab/filebrowser-extension":["@jupyterlab/filebrowser-extension:browser","@jupyterlab/filebrowser-extension:create-new-language-file","@jupyterlab/filebrowser-extension:download","@jupyterlab/filebrowser-extension:file-upload-status","@jupyterlab/filebrowser-extension:open-with","@jupyterlab/filebrowser-extension:search","@jupyterlab/filebrowser-extension:share-file"],"@jupyter-notebook/tree-extension":true,"@jupyterlab/running-extension":["@jupyterlab/running-extension:plugin"],"@jupyterlab/settingeditor-extension":true},"/notebooks":{"@jupyterlab/celltags-extension":true,"@jupyterlab/cell-toolbar-extension":true,"@jupyterlab/debugger-extension":["@jupyterlab/debugger-extension:completions","@jupyterlab/debugger-extension:config","@jupyterlab/debugger-extension:debug-console","@jupyterlab/debugger-extension:main","@jupyterlab/debugger-extension:notebooks","@jupyterlab/debugger-extension:service","@jupyterlab/debugger-extension:sidebar","@jupyterlab/debugger-extension:sources","@jupyterlab/debugger-extension:display-registry"],"@jupyterlab/logconsole-extension":true,"@jupyterlab/metadataform-extension":true,"@jupyterlab/notebook-extension":["@jupyterlab/notebook-extension:active-cell-tool","@jupyterlab/notebook-extension:completer","@jupyterlab/notebook-extension:copy-output","@jupyterlab/notebook-extension:metadata-editor","@jupyterlab/notebook-extension:search","@jupyterlab/notebook-extension:toc","@jupyterlab/notebook-extension:tools","@jupyterlab/notebook-extension:update-raw-mimetype"],"@jupyterlab/toc-extension":["@jupyterlab/toc-extension:registry","@jupyterlab/toc-extension:tracker"],"@jupyterlab/tooltip-extension":["@jupyterlab/tooltip-extension:manager","@jupyterlab/tooltip-extension:notebooks"]},"/consoles":{"@jupyterlab/tooltip-extension":["@jupyterlab/tooltip-extension:manager","@jupyterlab/tooltip-extension:consoles"]},"/edit":{"@jupyterlab/fileeditor-extension":["@jupyterlab/fileeditor-extension:completer","@jupyterlab/fileeditor-extension:search"],"@jupyterlab/markdownviewer-extension":true}}');let m=new a.PluginRegistry,x=o(33177).NotebookApp;m.registerPlugins(d);let g=o(19420).IServiceManager,h=await m.resolveRequiredService(g),f=new x({pluginRegistry:m,serviceManager:h,mimeExtensions:t,availablePlugins:s});"true"===(n.PageConfig.getOption("exposeAppInBrowser")||"").toLowerCase()&&(window.jupyterapp=f),await f.start()}o(77570),o(71439),window.addEventListener("load",l)},77570(e,t,o){"use strict";o(50595),o(75200),o(41680);var n=o(95292),a=o.n(n),r=o(49893),i=o.n(r),l=o(9383),p=o.n(l),s=o(56884),u=o.n(s),b=o(99088),d=o.n(b),c=o(27997),j=o.n(c),y=o(5296),m={};m.styleTagTransform=j(),m.setAttributes=u(),m.insert=p().bind(null,"head"),m.domAPI=i(),m.insertStyleElement=d(),a()(y.A,m),y.A&&y.A.locals&&y.A.locals;var x=o(83146),g={};g.styleTagTransform=j(),g.setAttributes=u(),g.insert=p().bind(null,"head"),g.domAPI=i(),g.insertStyleElement=d(),a()(x.A,g),x.A&&x.A.locals&&x.A.locals,o(48095);var h=o(98384),f={};f.styleTagTransform=j(),f.setAttributes=u(),f.insert=p().bind(null,"head"),f.domAPI=i(),f.insertStyleElement=d(),a()(h.A,f),h.A&&h.A.locals&&h.A.locals;var k=o(9755),w={};w.styleTagTransform=j(),w.setAttributes=u(),w.insert=p().bind(null,"head"),w.domAPI=i(),w.insertStyleElement=d(),a()(k.A,w),k.A&&k.A.locals&&k.A.locals;var v=o(99433),A={};A.styleTagTransform=j(),A.setAttributes=u(),A.insert=p().bind(null,"head"),A.domAPI=i(),A.insertStyleElement=d(),a()(v.A,A),v.A&&v.A.locals&&v.A.locals;var T=o(11781),P={};P.styleTagTransform=j(),P.setAttributes=u(),P.insert=p().bind(null,"head"),P.domAPI=i(),P.insertStyleElement=d(),a()(T.A,P),T.A&&T.A.locals&&T.A.locals;var B=o(75131),S={};S.styleTagTransform=j(),S.setAttributes=u(),S.insert=p().bind(null,"head"),S.domAPI=i(),S.insertStyleElement=d(),a()(B.A,S),B.A&&B.A.locals&&B.A.locals;var C=o(7965),N={};N.styleTagTransform=j(),N.setAttributes=u(),N.insert=p().bind(null,"head"),N.domAPI=i(),N.insertStyleElement=d(),a()(C.A,N),C.A&&C.A.locals&&C.A.locals;var E=o(29404),D={};D.styleTagTransform=j(),D.setAttributes=u(),D.insert=p().bind(null,"head"),D.domAPI=i(),D.insertStyleElement=d(),a()(E.A,D),E.A&&E.A.locals&&E.A.locals,o(68591);var L=o(95296),z={};z.styleTagTransform=j(),z.setAttributes=u(),z.insert=p().bind(null,"head"),z.domAPI=i(),z.insertStyleElement=d(),a()(L.A,z),L.A&&L.A.locals&&L.A.locals;var F=o(58688),M={};M.styleTagTransform=j(),M.setAttributes=u(),M.insert=p().bind(null,"head"),M.domAPI=i(),M.insertStyleElement=d(),a()(F.A,M),F.A&&F.A.locals&&F.A.locals;var _=o(25455),I={};I.styleTagTransform=j(),I.setAttributes=u(),I.insert=p().bind(null,"head"),I.domAPI=i(),I.insertStyleElement=d(),a()(_.A,I),_.A&&_.A.locals&&_.A.locals,o(19894),o(89277),o(59344),o(47568),o(20550),o(14640),o(2645),o(66156),o(71462),o(44349),o(48175),o(5207),o(42442),o(50301),o(92045),o(93712),o(88551),o(51357),o(58119),o(71670),o(66029),o(85317),o(69471),o(17704),o(98419),o(39601),o(63562),o(38257),o(26021),o(53555),o(91749),o(69632),o(89821),o(83199),o(93992),o(59158),o(21462),o(43684),o(313),o(41315),o(28114),o(52400)},98384(e,t,o){"use strict";o.d(t,{A:()=>l});var n=o(8645),a=o.n(n),r=o(60278),i=o.n(r)()(a());i.push([e.id,`/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
|
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

.jp-NotebookSpacer {
  flex-grow: 1;
  flex-shrink: 1;
}

.jp-MainAreaWidget {
  height: 100%;
}

.jp-Toolbar > .jp-Toolbar-item {
  height: unset;
}

#jp-UserMenu {
  flex: 0 0 auto;
  display: flex;
  text-align: center;
  margin-top: 8px;
}

.jp-MimeDocument .jp-RenderedJSON {
  background: var(--jp-layout-color0);
}

/* Hide the stub toolbar that appears above terminals and documents */

.jp-MainAreaWidget > .jp-Toolbar-micro {
  display: none;
}

#jp-NotebookLogo {
  /* bring logo to the front so it is selectable by tab*/
  z-index: 10;
}

/* Hide the notification status item */
.jp-Notification-Status {
  display: none;
}
`,""]);let l=i},5296(e,t,o){"use strict";o.d(t,{A:()=>l});var n=o(8645),a=o.n(n),r=o(60278),i=o.n(r)()(a());i.push([e.id,`/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

:root {
  --jp-private-topbar-height: 28px;
  /* Override the layout-2 color for the dark theme */
  --md-grey-800: #323232;
  --jp-notebook-max-width: 1200px;
}

/*
  Override the default background
  See https://github.com/jupyterlab/jupyterlab/pull/16519 for more information
*/
body.jp-ThemedContainer {
  margin: 0;
  padding: 0;
  background: var(--jp-layout-color2);
}

#main.jp-ThemedContainer {
  position: absolute;
  top: 0;
  left: 0;
  right: 0;
  bottom: 0;
  background: var(--jp-layout-color2);
}

#top-panel-wrapper {
  min-height: calc(1.5 * var(--jp-private-topbar-height));
  border-bottom: var(--jp-border-width) solid var(--jp-border-color0);
  background: var(--jp-layout-color1);
}

#top-panel {
  display: flex;
  min-height: calc(1.5 * var(--jp-private-topbar-height));
  padding-left: 5px;
  padding-right: 5px;
  margin-left: auto;
  margin-right: auto;
  max-width: 1200px;
}

#menu-panel-wrapper {
  min-height: var(--jp-private-topbar-height);
  background: var(--jp-layout-color1);
  border-bottom: var(--jp-border-width) solid var(--jp-border-color0);
  box-shadow: var(--jp-elevation-z1);
}

#menu-panel {
  display: flex;
  min-height: var(--jp-private-topbar-height);
  background: var(--jp-layout-color1);
  padding-left: 5px;
  padding-right: 5px;
  margin-left: auto;
  margin-right: auto;
  max-width: var(--jp-notebook-max-width);
}

#main-panel {
  margin-left: auto;
  margin-right: auto;
  max-width: var(--jp-notebook-max-width);
}

#spacer-widget-top {
  min-height: 16px;
}

/* Only edit pages should have a bottom space */

body[data-notebook='edit'] #spacer-widget-bottom {
  min-height: 16px;
}

/* Special case notebooks as document oriented pages */

[data-notebook]:not(body[data-notebook='notebooks']) #main-panel {
  box-shadow: var(--jp-elevation-z4);
}

.jp-TreePanel > .lm-TabPanel-stackedPanel {
  box-shadow: var(--jp-elevation-z4);
}

body[data-notebook='notebooks'] #main-panel {
  margin-left: unset;
  margin-right: unset;
  max-width: unset;
}

body[data-notebook='notebooks'] #spacer-widget-top {
  min-height: unset;
}

#main-panel > .jp-TreePanel {
  padding: 0px 5px;
}

@media only screen and (max-width: 760px) {
  #main-panel > .jp-TreePanel {
    margin: 0px -5px;
  }
}
`,""]);let l=i},83146(e,t,o){"use strict";o.d(t,{A:()=>l});var n=o(8645),a=o.n(n),r=o(60278),i=o.n(r)()(a());i.push([e.id,`/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|
| Adapted from JupyterLab's packages/application/style/sidepanel.css.
|----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
| Variables
|----------------------------------------------------------------------------*/

:root {
  --jp-private-sidebar-tab-width: 32px;
}

/*-----------------------------------------------------------------------------
| SideBar
|----------------------------------------------------------------------------*/

/* Stack panels */

#jp-right-stack,
#jp-left-stack {
  display: flex;
  flex-direction: column;
  min-width: var(--jp-sidebar-min-width);
}

#jp-left-stack .jp-SidePanel-collapse,
#jp-right-stack .jp-SidePanel-collapse {
  display: flex;
  flex: 0 0 auto;
  min-height: 0;
  padding: 0;
}

#jp-left-stack .jp-SidePanel-collapse {
  justify-content: right;
}

#jp-right-stack .jp-SidePanel-collapse {
  justify-content: left;
}

#jp-left-stack .lm-StackedPanel,
#jp-right-stack .lm-StackedPanel {
  flex: 1 1 auto;
}
`,""]);let l=i},9755(e,t,o){"use strict";o.d(t,{A:()=>l});var n=o(8645),a=o.n(n),r=o(60278),i=o.n(r)()(a());i.push([e.id,"",""]);let l=i},99433(e,t,o){"use strict";o.d(t,{A:()=>l});var n=o(8645),a=o.n(n),r=o(60278),i=o.n(r)()(a());i.push([e.id,`.jp-Document {
  height: 100%;
}
`,""]);let l=i},11781(e,t,o){"use strict";o.d(t,{A:()=>l});var n=o(8645),a=o.n(n),r=o(60278),i=o.n(r)()(a());i.push([e.id,"",""]);let l=i},75131(e,t,o){"use strict";o.d(t,{A:()=>l});var n=o(8645),a=o.n(n),r=o(60278),i=o.n(r)()(a());i.push([e.id,`.jp-AboutNotebook .jp-Dialog-header {
  justify-content: center;
  padding: 0;
}

.jp-AboutNotebook-header {
  display: flex;
  flex-direction: row;
  align-items: center;
  padding: var(--jp-flat-button-padding);
}

.jp-AboutNotebook-header-text {
  margin-left: 16px;
}

.jp-AboutNotebook-version {
  color: var(--jp-ui-font-color1);
  font-size: var(--jp-ui-font-size1);
  padding-bottom: 30px;
  font-weight: 400;
  letter-spacing: 0.4px;
  line-height: 1.12;
  min-width: 360px;
  text-align: center;
}

.jp-AboutNotebook-body {
  display: flex;
  font-size: var(--jp-ui-font-size2);
  padding: var(--jp-flat-button-padding);
  color: var(--jp-ui-font-color1);
  text-align: center;
  flex-direction: column;
  min-width: 360px;
  overflow: hidden;
}

.jp-AboutNotebook-about-body pre {
  white-space: pre-wrap;
}

.jp-AboutNotebook-about-externalLinks {
  display: flex;
  flex-direction: column;
  justify-content: flex-start;
  align-items: flex-start;
  color: var(--jp-warn-color0);
}

.jp-AboutNotebook-about-copyright {
  padding-top: 25px;
}
`,""]);let l=i},7965(e,t,o){"use strict";o.d(t,{A:()=>s});var n=o(8645),a=o.n(n),r=o(60278),i=o.n(r),l=o(16623),p=i()(a());p.i(l.A),p.push([e.id,`/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
|
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

/**
  Document oriented look for the notebook.
  This includes changes to the look and feel of the JupyterLab Notebook
  component like:
  - scrollbar to the right of the page
  - drop shadow on the notebook
  - smaller empty space at the bottom of the notebook
  - compact view on mobile
*/

/* Make the notebook take up the full width of the page when jp-mod-fullwidth is set */

body[data-notebook='notebooks']
  .jp-NotebookPanel.jp-mod-fullwidth
  .jp-WindowedPanel-outer {
  padding-left: unset;
  padding-right: unset !important;
  width: unset;
}

/* Keep the notebook centered on the page */

body[data-notebook='notebooks'] .jp-NotebookPanel-toolbar {
  padding-left: calc(calc(100% - var(--jp-notebook-max-width)) * 0.5);
  padding-right: calc(calc(100% - var(--jp-notebook-max-width)) * 0.5);
}

body[data-notebook='notebooks'] .jp-WindowedPanel-outer {
  width: unset !important;
  padding-top: unset;
  padding-left: calc(calc(100% - var(--jp-notebook-max-width)) * 0.5);
  padding-right: calc(
    calc(
        100% - var(--jp-notebook-max-width) - var(--jp-notebook-padding-offset)
      ) * 0.5
  ) !important;
  background: var(--jp-layout-color2);
}

body[data-notebook='notebooks'] .jp-WindowedPanel-inner {
  margin-top: var(--jp-notebook-toolbar-margin-bottom);
  /* Adjustments for the extra top and bottom notebook padding */
  margin-bottom: calc(4 * var(--jp-notebook-padding));
}

body[data-notebook='notebooks'] .jp-Notebook-cell {
  background: var(--jp-layout-color0);
}

/* Empty space at the bottom of the notebook (similar to classic) */
body[data-notebook='notebooks']
  .jp-Notebook.jp-mod-scrollPastEnd
  .jp-WindowedPanel-outer::after {
  min-height: 100px;
}

/* Fix background colors */

body[data-notebook='notebooks'] .jp-WindowedPanel-outer > * {
  background: var(--jp-layout-color0);
}

body[data-notebook='notebooks']
  .jp-Notebook.jp-mod-commandMode
  .jp-Cell.jp-mod-active.jp-mod-selected:not(.jp-mod-multiSelected) {
  background: var(--jp-layout-color0) !important;
}

body[data-notebook='notebooks']
  .jp-Notebook
  .jp-Notebook-cell:not(:first-child)::before {
  content: ' ';
  height: 100%;
  position: absolute;
  top: 0;
  width: 11px;
}

/* Cell toolbar adjustments */

body[data-notebook='notebooks'] .jp-cell-toolbar {
  background: unset;
  box-shadow: unset;
}

/** first code cell on mobile
    (keep the selector above the media query)
*/
body[data-notebook='notebooks']
  .jp-CodeCell[data-windowed-list-index='0']
  .jp-cell-toolbar {
  top: unset;
}

@media only screen and (max-width: 760px) {
  /* first code cell on mobile */
  body[data-notebook='notebooks']
    .jp-CodeCell[data-windowed-list-index='0']
    .jp-cell-toolbar {
    top: var(--jp-notebook-padding);
  }

  body[data-notebook='notebooks'] .jp-MarkdownCell .jp-cell-toolbar,
  body[data-notebook='notebooks'] .jp-RawCell .jp-cell-toolbar {
    top: calc(0.5 * var(--jp-notebook-padding));
  }
}

/* Tweak the notebook footer (to add a new cell) */
body[data-notebook='notebooks'] .jp-Notebook-footer {
  background: unset;
  width: 100%;
  margin-left: unset;
}

/* Mobile View */

body[data-format='mobile'] .jp-NotebookCheckpoint {
  display: none;
}

body[data-format='mobile'] .jp-WindowedPanel-outer > *:first-child {
  margin-top: 0;
}

body[data-format='mobile'] .jp-ToolbarButton .jp-DebuggerBugButton {
  display: none;
}

body[data-notebook='notebooks'] .jp-WindowedPanel-viewport {
  background: var(--jp-layout-color0);
  box-shadow: var(--jp-elevation-z4);

  /* Extra padding at the top and bottom so the notebook looks nicer */
  padding-top: calc(2 * var(--jp-notebook-padding));
  padding-bottom: calc(2 * var(--jp-notebook-padding));
}

/* Notebook box shadow */

body[data-notebook='notebooks']
  .jp-Notebook
  > *:first-child:last-child::before {
  content: '';
  position: absolute;
  top: 0;
  bottom: 0;
  left: 0;
  right: 0;
  box-shadow: 0px 0px 12px 1px var(--jp-shadow-umbra-color);
}

/* Additional customizations of the components on the notebook page */

.jp-NotebookKernelLogo {
  flex: 0 0 auto;
  display: flex;
  align-items: center;
  text-align: center;
  margin-right: 8px;
}

.jp-NotebookKernelLogo img {
  max-width: 28px;
  max-height: 28px;
  display: flex;
}

.jp-NotebookKernelStatus {
  margin: 0;
  font-weight: normal;
  font-size: var(--jp-ui-font-size1);
  color: var(--jp-ui-font-color0);
  font-family: var(--jp-ui-font-family);
  line-height: var(--jp-private-title-panel-height);
  padding-left: var(--jp-kernel-status-padding);
  padding-right: var(--jp-kernel-status-padding);
}

.jp-NotebookKernelStatus-error {
  background-color: var(--jp-error-color0);
}

.jp-NotebookKernelStatus-warn {
  background-color: var(--jp-warn-color0);
}

.jp-NotebookKernelStatus-info {
  background-color: var(--jp-info-color0);
}

.jp-NotebookKernelStatus-fade {
  animation: 0.5s fade-out forwards;
}

.jp-NotebookTrustedStatus {
  background: var(--jp-layout-color1);
  color: var(--jp-ui-font-color1);
  margin-top: 4px;
  margin-bottom: 4px;
  border: solid 1px var(--jp-border-color2);
  cursor: help;
}

.jp-NotebookTrustedStatus-not-trusted {
  cursor: pointer;
}

@keyframes fade-out {
  0% {
    opacity: 1;
  }
  100% {
    opacity: 0;
  }
}

#jp-title h1 {
  cursor: pointer;
  font-size: 18px;
  margin: 0;
  font-weight: normal;
  color: var(--jp-ui-font-color0);
  font-family: var(--jp-ui-font-family);
  line-height: calc(1.5 * var(--jp-private-title-panel-height));
  text-overflow: ellipsis;
  overflow: hidden;
  white-space: nowrap;
}

#jp-title h1:hover {
  background: var(--jp-layout-color2);
}

.jp-NotebookCheckpoint {
  font-size: 14px;
  margin-left: 5px;
  margin-right: 5px;
  font-weight: normal;
  color: var(--jp-ui-font-color0);
  font-family: var(--jp-ui-font-family);
  line-height: calc(1.5 * var(--jp-private-title-panel-height));
  text-overflow: ellipsis;
  overflow: hidden;
  white-space: nowrap;
}

.jp-skiplink {
  position: absolute;
  top: -100em;
}

.jp-skiplink:focus-within {
  position: absolute;
  z-index: 10000;
  top: 0;
  left: 46%;
  margin: 0 auto;
  padding: 1em;
  width: 15%;
  box-shadow: var(--jp-elevation-z4);
  border-radius: 4px;
  background: var(--jp-layout-color0);
  text-align: center;
}

.jp-skiplink:focus-within a {
  text-decoration: underline;
  color: var(--jp-content-link-color);
}
`,""]);let s=p},16623(e,t,o){"use strict";o.d(t,{A:()=>l});var n=o(8645),a=o.n(n),r=o(60278),i=o.n(r)()(a());i.push([e.id,`:root {
  --jp-notebook-toolbar-margin-bottom: 20px;
  --jp-notebook-padding-offset: 20px;

  --jp-kernel-status-padding: 5px;
}
`,""]);let l=i},29404(e,t,o){"use strict";o.d(t,{A:()=>l});var n=o(8645),a=o.n(n),r=o(60278),i=o.n(r)()(a());i.push([e.id,"",""]);let l=i},58688(e,t,o){"use strict";o.d(t,{A:()=>l});var n=o(8645),a=o.n(n),r=o(60278),i=o.n(r)()(a());i.push([e.id,`/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
|
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

.jp-FileBrowser-toolbar .jp-Toolbar-item.jp-DropdownMenu,
.jp-FileBrowser-toolbar .jp-Toolbar-item.jp-ToolbarButton,
.jp-FileBrowser-toolbar .jp-Toolbar-item.jp-CommandToolbarButton {
  border: solid 1px var(--jp-border-color2);
  margin: 1px;
  padding: 0px;
}

.jp-FileBrowser-toolbar > .jp-Toolbar-item.jp-ToolbarButton:hover,
.jp-FileBrowser-toolbar > .jp-Toolbar-item.jp-CommandToolbarButton:hover,
.jp-FileBrowser-toolbar > .jp-Toolbar-item.jp-DropdownMenu:hover {
  background: var(--neutral-fill-stealth-hover);
}

.jp-FileBrowser-toolbar .lm-MenuBar-item {
  height: var(--jp-private-toolbar-height);
  display: inline-flex;
  align-items: center;
}

.jp-FileBrowser-toolbar .jp-ToolbarButtonComponent {
  height: var(--jp-flat-button-height);
}

.jp-FileBrowser-toolbar jp-button.jp-ToolbarButtonComponent:hover {
  background: inherit;
}

.jp-DirListing-content .jp-DirListing-checkboxWrapper {
  visibility: visible;
}

/* Action buttons */

.jp-FileBrowser-toolbar > .jp-FileAction > .jp-ToolbarButtonComponent > svg {
  display: none;
}

.jp-FileBrowser-toolbar > #fileAction-delete {
  background-color: var(--jp-error-color1);
}

.jp-FileBrowser-toolbar
  .jp-ToolbarButtonComponent[data-command='filebrowser:delete']
  .jp-ToolbarButtonComponent-label {
  color: var(--jp-ui-inverse-font-color1);
}

.jp-FileBrowser-toolbar .jp-FileAction {
  border: solid 1px var(--jp-border-color2);
  margin: 1px;
  min-height: var(--jp-private-toolbar-height);
}

body[data-format='mobile'] #fileAction-placeholder {
  display: none;
}
`,""]);let l=i},95296(e,t,o){"use strict";o.d(t,{A:()=>l});var n=o(8645),a=o.n(n),r=o(60278),i=o.n(r)()(a());i.push([e.id,`.jp-FileBrowser {
  height: 100%;
}

.lm-TabPanel {
  height: 100%;
}

.jp-TreePanel .lm-TabPanel-tabBar {
  overflow: visible;
  min-height: 32px;
  border-bottom: unset;
  height: var(--jp-private-toolbar-height);
}

.jp-TreePanel .lm-TabBar-content {
  height: 100%;
}

.jp-TreePanel .lm-TabBar-tab {
  flex: 0 1 auto;
  color: var(--jp-ui-font-color0);
  font-size: var(--jp-ui-font-size1);
  height: 100%;
}

.jp-TreePanel .lm-TabBar-tabLabel {
  padding-left: 5px;
  padding-right: 5px;
}

.jp-FileBrowser-toolbar.jp-Toolbar .jp-ToolbarButtonComponent {
  width: unset;
}

.jp-FileBrowser-toolbar > .jp-Toolbar-item {
  flex-direction: column;
  justify-content: center;
}

.jp-DropdownMenu .lm-MenuBar-itemIcon svg {
  vertical-align: sub;
}

jp-button[data-command='filebrowser:refresh'] .jp-ToolbarButtonComponent-label {
  display: none;
}

.jp-TreePanel .lm-TabBar-tabIcon svg {
  vertical-align: sub;
}
`,""]);let l=i},25455(e,t,o){"use strict";o.d(t,{A:()=>l});var n=o(8645),a=o.n(n),r=o(60278),i=o.n(r)()(a());i.push([e.id,"",""]);let l=i}}]);