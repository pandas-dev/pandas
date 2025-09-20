/*! For license information please see 232.5419cbec68e3fd0cf431.js.LICENSE.txt */
"use strict";(self["webpackChunk_jupyterlab_application_top"]=self["webpackChunk_jupyterlab_application_top"]||[]).push([[232],{50232:(e,t,i)=>{i.r(t);i.d(t,{ARIAGlobalStatesAndProperties:()=>je,Accordion:()=>qe,AccordionExpandMode:()=>Be,AccordionItem:()=>He,Anchor:()=>_e,AnchoredRegion:()=>Os,Avatar:()=>Ms,Badge:()=>Hs,BaseProgress:()=>zr,Breadcrumb:()=>Bs,BreadcrumbItem:()=>zs,Button:()=>Qs,Calendar:()=>eo,CalendarTitleTemplate:()=>ho,Card:()=>go,CheckableFormAssociated:()=>Gs,Checkbox:()=>wo,Combobox:()=>Zo,ComboboxAutocomplete:()=>Qo,ComponentPresentation:()=>Se,Container:()=>U,ContainerConfiguration:()=>V,ContainerImpl:()=>ge,DI:()=>q,DataGrid:()=>lo,DataGridCell:()=>ro,DataGridCellTypes:()=>io,DataGridRow:()=>ao,DataGridRowTypes:()=>so,DateFormatter:()=>Js,DefaultComponentPresentation:()=>Ae,DefaultResolver:()=>H,DelegatesARIAButton:()=>Zs,DelegatesARIACombobox:()=>Jo,DelegatesARIALink:()=>Ke,DelegatesARIAListbox:()=>Wo,DelegatesARIAListboxOption:()=>_o,DelegatesARIASearch:()=>ea,DelegatesARIASelect:()=>oa,DelegatesARIATextbox:()=>Fr,DelegatesARIAToolbar:()=>Ha,DesignSystem:()=>Pn,DesignToken:()=>Dn,Dialog:()=>qn,Disclosure:()=>Gn,Divider:()=>Zn,DividerRole:()=>Qn,ElementDisambiguation:()=>Sn,FactoryImpl:()=>de,Flipper:()=>tr,FlipperDirection:()=>Jn,FlyoutPosBottom:()=>Rs,FlyoutPosBottomFill:()=>As,FlyoutPosTallest:()=>Ds,FlyoutPosTallestFill:()=>Fs,FlyoutPosTop:()=>Es,FlyoutPosTopFill:()=>Ss,FormAssociated:()=>Ws,FoundationElement:()=>Fe,FoundationElementRegistry:()=>Me,GenerateHeaderOptions:()=>to,HorizontalScroll:()=>Wr,Listbox:()=>Ko,ListboxElement:()=>sr,ListboxOption:()=>jo,MatchMediaBehavior:()=>Ka,MatchMediaStyleSheetBehavior:()=>Wa,Menu:()=>Tr,MenuItem:()=>kr,MenuItemRole:()=>wr,NumberField:()=>Pr,Picker:()=>br,PickerList:()=>lr,PickerListItem:()=>hr,PickerMenu:()=>nr,PickerMenuOption:()=>ar,PropertyStyleSheetBehavior:()=>Qa,Radio:()=>Kr,RadioGroup:()=>qr,Registration:()=>xe,ResolverBuilder:()=>M,ResolverImpl:()=>re,Search:()=>Jr,Select:()=>sa,SelectPosition:()=>Go,ServiceLocator:()=>j,Skeleton:()=>aa,Slider:()=>va,SliderLabel:()=>ca,SliderMode:()=>ma,StartEnd:()=>o,Switch:()=>Ca,Tab:()=>Ia,TabPanel:()=>wa,Tabs:()=>Ta,TabsOrientation:()=>Oa,TextArea:()=>Aa,TextAreaResize:()=>Ea,TextField:()=>Ar,TextFieldType:()=>Sr,Toolbar:()=>Pa,Tooltip:()=>Na,TooltipPosition:()=>za,TreeItem:()=>Ua,TreeView:()=>_a,accordionItemTemplate:()=>d,accordionTemplate:()=>Ve,all:()=>J,anchorTemplate:()=>Ue,anchoredRegionTemplate:()=>We,applyMixins:()=>Pe,avatarTemplate:()=>Ls,badgeTemplate:()=>Ps,breadcrumbItemTemplate:()=>Vs,breadcrumbTemplate:()=>Ns,buttonTemplate:()=>qs,calendarCellTemplate:()=>uo,calendarRowTemplate:()=>po,calendarTemplate:()=>vo,calendarWeekdayTemplate:()=>co,cardTemplate:()=>bo,checkboxTemplate:()=>yo,comboboxTemplate:()=>en,composedContains:()=>dn,composedParent:()=>ln,darkModeStylesheetBehavior:()=>Xa,dataGridCellTemplate:()=>an,dataGridRowTemplate:()=>rn,dataGridTemplate:()=>sn,dialogTemplate:()=>Nn,disabledCursor:()=>Za,disclosureTemplate:()=>Wn,display:()=>el,dividerTemplate:()=>Xn,endSlotTemplate:()=>n,endTemplate:()=>a,flipperTemplate:()=>er,focusVisible:()=>tl,forcedColorsStylesheetBehavior:()=>Ga,getDirection:()=>Is,hidden:()=>Ja,horizontalScrollTemplate:()=>Gr,ignore:()=>ie,inject:()=>K,interactiveCalendarGridTemplate:()=>fo,isListboxOption:()=>Uo,isTreeItemElement:()=>qa,lazy:()=>ee,lightModeStylesheetBehavior:()=>Ya,listboxOptionTemplate:()=>ir,listboxTemplate:()=>or,menuItemTemplate:()=>Ir,menuTemplate:()=>Or,newInstanceForScope:()=>se,newInstanceOf:()=>oe,noninteractiveCalendarTemplate:()=>mo,numberFieldTemplate:()=>Er,optional:()=>te,pickerListItemTemplate:()=>xr,pickerListTemplate:()=>Cr,pickerMenuOptionTemplate:()=>yr,pickerMenuTemplate:()=>gr,pickerTemplate:()=>pr,progressRingTemplate:()=>Vr,progressTemplate:()=>Nr,radioGroupTemplate:()=>Br,radioTemplate:()=>Ur,reflectAttributes:()=>Kn,roleForMenuItem:()=>$r,searchTemplate:()=>Yr,selectTemplate:()=>na,singleton:()=>Q,skeletonTemplate:()=>ra,sliderLabelTemplate:()=>la,sliderTemplate:()=>ua,startSlotTemplate:()=>r,startTemplate:()=>l,supportsElementInternals:()=>_s,switchTemplate:()=>ba,tabPanelTemplate:()=>xa,tabTemplate:()=>$a,tabsTemplate:()=>ka,textAreaTemplate:()=>Ra,textFieldTemplate:()=>Fa,toolbarTemplate:()=>La,tooltipTemplate:()=>Va,transient:()=>G,treeItemTemplate:()=>Ba,treeViewTemplate:()=>ja,validateKey:()=>we,whitespaceFilter:()=>Xr});var s=i(29690);class o{handleStartContentChange(){this.startContainer.classList.toggle("start",this.start.assignedNodes().length>0)}handleEndContentChange(){this.endContainer.classList.toggle("end",this.end.assignedNodes().length>0)}}const n=(e,t)=>(0,s.html)`
    <span
        part="end"
        ${(0,s.ref)("endContainer")}
        class=${e=>t.end?"end":void 0}
    >
        <slot name="end" ${(0,s.ref)("end")} @slotchange="${e=>e.handleEndContentChange()}">
            ${t.end||""}
        </slot>
    </span>
`;const r=(e,t)=>(0,s.html)`
    <span
        part="start"
        ${(0,s.ref)("startContainer")}
        class="${e=>t.start?"start":void 0}"
    >
        <slot
            name="start"
            ${(0,s.ref)("start")}
            @slotchange="${e=>e.handleStartContentChange()}"
        >
            ${t.start||""}
        </slot>
    </span>
`;const a=(0,s.html)`
    <span part="end" ${(0,s.ref)("endContainer")}>
        <slot
            name="end"
            ${(0,s.ref)("end")}
            @slotchange="${e=>e.handleEndContentChange()}"
        ></slot>
    </span>
`;const l=(0,s.html)`
    <span part="start" ${(0,s.ref)("startContainer")}>
        <slot
            name="start"
            ${(0,s.ref)("start")}
            @slotchange="${e=>e.handleStartContentChange()}"
        ></slot>
    </span>
`;const d=(e,t)=>(0,s.html)`
    <template class="${e=>e.expanded?"expanded":""}">
        <div
            class="heading"
            part="heading"
            role="heading"
            aria-level="${e=>e.headinglevel}"
        >
            <button
                class="button"
                part="button"
                ${(0,s.ref)("expandbutton")}
                aria-expanded="${e=>e.expanded}"
                aria-controls="${e=>e.id}-panel"
                id="${e=>e.id}"
                @click="${(e,t)=>e.clickHandler(t.event)}"
            >
                <span class="heading-content" part="heading-content">
                    <slot name="heading"></slot>
                </span>
            </button>
            ${r(e,t)}
            ${n(e,t)}
            <span class="icon" part="icon" aria-hidden="true">
                <slot name="expanded-icon" part="expanded-icon">
                    ${t.expandedIcon||""}
                </slot>
                <slot name="collapsed-icon" part="collapsed-icon">
                    ${t.collapsedIcon||""}
                </slot>
            <span>
        </div>
        <div
            class="region"
            part="region"
            id="${e=>e.id}-panel"
            role="region"
            aria-labelledby="${e=>e.id}"
        >
            <slot></slot>
        </div>
    </template>
`;var h=function(e,t){h=Object.setPrototypeOf||{__proto__:[]}instanceof Array&&function(e,t){e.__proto__=t}||function(e,t){for(var i in t)if(t.hasOwnProperty(i))e[i]=t[i]};return h(e,t)};function c(e,t){h(e,t);function i(){this.constructor=e}e.prototype=t===null?Object.create(t):(i.prototype=t.prototype,new i)}var u=function(){u=Object.assign||function e(t){for(var i,s=1,o=arguments.length;s<o;s++){i=arguments[s];for(var n in i)if(Object.prototype.hasOwnProperty.call(i,n))t[n]=i[n]}return t};return u.apply(this,arguments)};function p(e,t){var i={};for(var s in e)if(Object.prototype.hasOwnProperty.call(e,s)&&t.indexOf(s)<0)i[s]=e[s];if(e!=null&&typeof Object.getOwnPropertySymbols==="function")for(var o=0,s=Object.getOwnPropertySymbols(e);o<s.length;o++){if(t.indexOf(s[o])<0&&Object.prototype.propertyIsEnumerable.call(e,s[o]))i[s[o]]=e[s[o]]}return i}function f(e,t,i,s){var o=arguments.length,n=o<3?t:s===null?s=Object.getOwnPropertyDescriptor(t,i):s,r;if(typeof Reflect==="object"&&typeof Reflect.decorate==="function")n=Reflect.decorate(e,t,i,s);else for(var a=e.length-1;a>=0;a--)if(r=e[a])n=(o<3?r(n):o>3?r(t,i,n):r(t,i))||n;return o>3&&n&&Object.defineProperty(t,i,n),n}function m(e,t){return function(i,s){t(i,s,e)}}function v(e,t){if(typeof Reflect==="object"&&typeof Reflect.metadata==="function")return Reflect.metadata(e,t)}function b(e,t,i,s){function o(e){return e instanceof i?e:new i((function(t){t(e)}))}return new(i||(i=Promise))((function(i,n){function r(e){try{l(s.next(e))}catch(t){n(t)}}function a(e){try{l(s["throw"](e))}catch(t){n(t)}}function l(e){e.done?i(e.value):o(e.value).then(r,a)}l((s=s.apply(e,t||[])).next())}))}function g(e,t){var i={label:0,sent:function(){if(n[0]&1)throw n[1];return n[1]},trys:[],ops:[]},s,o,n,r;return r={next:a(0),throw:a(1),return:a(2)},typeof Symbol==="function"&&(r[Symbol.iterator]=function(){return this}),r;function a(e){return function(t){return l([e,t])}}function l(r){if(s)throw new TypeError("Generator is already executing.");while(i)try{if(s=1,o&&(n=r[0]&2?o["return"]:r[0]?o["throw"]||((n=o["return"])&&n.call(o),0):o.next)&&!(n=n.call(o,r[1])).done)return n;if(o=0,n)r=[r[0]&2,n.value];switch(r[0]){case 0:case 1:n=r;break;case 4:i.label++;return{value:r[1],done:false};case 5:i.label++;o=r[1];r=[0];continue;case 7:r=i.ops.pop();i.trys.pop();continue;default:if(!(n=i.trys,n=n.length>0&&n[n.length-1])&&(r[0]===6||r[0]===2)){i=0;continue}if(r[0]===3&&(!n||r[1]>n[0]&&r[1]<n[3])){i.label=r[1];break}if(r[0]===6&&i.label<n[1]){i.label=n[1];n=r;break}if(n&&i.label<n[2]){i.label=n[2];i.ops.push(r);break}if(n[2])i.ops.pop();i.trys.pop();continue}r=t.call(e,i)}catch(a){r=[6,a];o=0}finally{s=n=0}if(r[0]&5)throw r[1];return{value:r[0]?r[1]:void 0,done:true}}}function y(e,t,i,s){if(s===undefined)s=i;e[s]=t[i]}function C(e,t){for(var i in e)if(i!=="default"&&!t.hasOwnProperty(i))t[i]=e[i]}function x(e){var t=typeof Symbol==="function"&&Symbol.iterator,i=t&&e[t],s=0;if(i)return i.call(e);if(e&&typeof e.length==="number")return{next:function(){if(e&&s>=e.length)e=void 0;return{value:e&&e[s++],done:!e}}};throw new TypeError(t?"Object is not iterable.":"Symbol.iterator is not defined.")}function w(e,t){var i=typeof Symbol==="function"&&e[Symbol.iterator];if(!i)return e;var s=i.call(e),o,n=[],r;try{while((t===void 0||t-- >0)&&!(o=s.next()).done)n.push(o.value)}catch(a){r={error:a}}finally{try{if(o&&!o.done&&(i=s["return"]))i.call(s)}finally{if(r)throw r.error}}return n}function $(){for(var e=[],t=0;t<arguments.length;t++)e=e.concat(w(arguments[t]));return e}function I(){for(var e=0,t=0,i=arguments.length;t<i;t++)e+=arguments[t].length;for(var s=Array(e),o=0,t=0;t<i;t++)for(var n=arguments[t],r=0,a=n.length;r<a;r++,o++)s[o]=n[r];return s}function k(e){return this instanceof k?(this.v=e,this):new k(e)}function O(e,t,i){if(!Symbol.asyncIterator)throw new TypeError("Symbol.asyncIterator is not defined.");var s=i.apply(e,t||[]),o,n=[];return o={},r("next"),r("throw"),r("return"),o[Symbol.asyncIterator]=function(){return this},o;function r(e){if(s[e])o[e]=function(t){return new Promise((function(i,s){n.push([e,t,i,s])>1||a(e,t)}))}}function a(e,t){try{l(s[e](t))}catch(i){c(n[0][3],i)}}function l(e){e.value instanceof k?Promise.resolve(e.value.v).then(d,h):c(n[0][2],e)}function d(e){a("next",e)}function h(e){a("throw",e)}function c(e,t){if(e(t),n.shift(),n.length)a(n[0][0],n[0][1])}}function T(e){var t,i;return t={},s("next"),s("throw",(function(e){throw e})),s("return"),t[Symbol.iterator]=function(){return this},t;function s(s,o){t[s]=e[s]?function(t){return(i=!i)?{value:k(e[s](t)),done:s==="return"}:o?o(t):t}:o}}function E(e){if(!Symbol.asyncIterator)throw new TypeError("Symbol.asyncIterator is not defined.");var t=e[Symbol.asyncIterator],i;return t?t.call(e):(e=typeof x==="function"?x(e):e[Symbol.iterator](),i={},s("next"),s("throw"),s("return"),i[Symbol.asyncIterator]=function(){return this},i);function s(t){i[t]=e[t]&&function(i){return new Promise((function(s,n){i=e[t](i),o(s,n,i.done,i.value)}))}}function o(e,t,i,s){Promise.resolve(s).then((function(t){e({value:t,done:i})}),t)}}function R(e,t){if(Object.defineProperty){Object.defineProperty(e,"raw",{value:t})}else{e.raw=t}return e}function D(e){if(e&&e.__esModule)return e;var t={};if(e!=null)for(var i in e)if(Object.hasOwnProperty.call(e,i))t[i]=e[i];t.default=e;return t}function S(e){return e&&e.__esModule?e:{default:e}}function A(e,t){if(!t.has(e)){throw new TypeError("attempted to get private field on non-instance")}return t.get(e)}function F(e,t,i){if(!t.has(e)){throw new TypeError("attempted to set private field on non-instance")}t.set(e,i);return i}const L=new Map;if(!("metadata"in Reflect)){Reflect.metadata=function(e,t){return function(i){Reflect.defineMetadata(e,t,i)}};Reflect.defineMetadata=function(e,t,i){let s=L.get(i);if(s===void 0){L.set(i,s=new Map)}s.set(e,t)};Reflect.getOwnMetadata=function(e,t){const i=L.get(t);if(i!==void 0){return i.get(e)}return void 0}}class M{constructor(e,t){this.container=e;this.key=t}instance(e){return this.registerResolver(0,e)}singleton(e){return this.registerResolver(1,e)}transient(e){return this.registerResolver(2,e)}callback(e){return this.registerResolver(3,e)}cachedCallback(e){return this.registerResolver(3,Ce(e))}aliasTo(e){return this.registerResolver(5,e)}registerResolver(e,t){const{container:i,key:s}=this;this.container=this.key=void 0;return i.registerResolver(s,new re(s,e,t))}}function P(e){const t=e.slice();const i=Object.keys(e);const s=i.length;let o;for(let n=0;n<s;++n){o=i[n];if(!Ee(o)){t[o]=e[o]}}return t}const H=Object.freeze({none(e){throw Error(`${e.toString()} not registered, did you forget to add @singleton()?`)},singleton(e){return new re(e,1,e)},transient(e){return new re(e,2,e)}});const V=Object.freeze({default:Object.freeze({parentLocator:()=>null,responsibleForOwnerRequests:false,defaultResolver:H.singleton})});const z=new Map;function N(e){return t=>Reflect.getOwnMetadata(e,t)}let B=null;const q=Object.freeze({createContainer(e){return new ge(null,Object.assign({},V.default,e))},findResponsibleContainer(e){const t=e.$$container$$;if(t&&t.responsibleForOwnerRequests){return t}return q.findParentContainer(e)},findParentContainer(e){const t=new CustomEvent(ve,{bubbles:true,composed:true,cancelable:true,detail:{container:void 0}});e.dispatchEvent(t);return t.detail.container||q.getOrCreateDOMContainer()},getOrCreateDOMContainer(e,t){if(!e){return B||(B=new ge(null,Object.assign({},V.default,t,{parentLocator:()=>null})))}return e.$$container$$||new ge(e,Object.assign({},V.default,t,{parentLocator:q.findParentContainer}))},getDesignParamtypes:N("design:paramtypes"),getAnnotationParamtypes:N("di:paramtypes"),getOrCreateAnnotationParamTypes(e){let t=this.getAnnotationParamtypes(e);if(t===void 0){Reflect.defineMetadata("di:paramtypes",t=[],e)}return t},getDependencies(e){let t=z.get(e);if(t===void 0){const i=e.inject;if(i===void 0){const i=q.getDesignParamtypes(e);const s=q.getAnnotationParamtypes(e);if(i===void 0){if(s===void 0){const i=Object.getPrototypeOf(e);if(typeof i==="function"&&i!==Function.prototype){t=P(q.getDependencies(i))}else{t=[]}}else{t=P(s)}}else if(s===void 0){t=P(i)}else{t=P(i);let e=s.length;let o;for(let i=0;i<e;++i){o=s[i];if(o!==void 0){t[i]=o}}const n=Object.keys(s);e=n.length;let r;for(let i=0;i<e;++i){r=n[i];if(!Ee(r)){t[r]=s[r]}}}}else{t=P(i)}z.set(e,t)}return t},defineProperty(e,t,i,o=false){const n=`$di_${t}`;Reflect.defineProperty(e,t,{get:function(){let e=this[n];if(e===void 0){const r=this instanceof HTMLElement?q.findResponsibleContainer(this):q.getOrCreateDOMContainer();e=r.get(i);this[n]=e;if(o&&this instanceof s.FASTElement){const s=this.$fastController;const o=()=>{const o=q.findResponsibleContainer(this);const r=o.get(i);const a=this[n];if(r!==a){this[n]=e;s.notify(t)}};s.subscribe({handleChange:o},"isConnected")}}return e}})},createInterface(e,t){const i=typeof e==="function"?e:t;const s=typeof e==="string"?e:e&&"friendlyName"in e?e.friendlyName||Ie:Ie;const o=typeof e==="string"?false:e&&"respectConnection"in e?e.respectConnection||false:false;const n=function(e,t,i){if(e==null||new.target!==undefined){throw new Error(`No registration for interface: '${n.friendlyName}'`)}if(t){q.defineProperty(e,t,n,o)}else{const t=q.getOrCreateAnnotationParamTypes(e);t[i]=n}};n.$isInterface=true;n.friendlyName=s==null?"(anonymous)":s;if(i!=null){n.register=function(e,t){return i(new M(e,t!==null&&t!==void 0?t:n))}}n.toString=function e(){return`InterfaceSymbol<${n.friendlyName}>`};return n},inject(...e){return function(t,i,s){if(typeof s==="number"){const i=q.getOrCreateAnnotationParamTypes(t);const o=e[0];if(o!==void 0){i[s]=o}}else if(i){q.defineProperty(t,i,e[0])}else{const i=s?q.getOrCreateAnnotationParamTypes(s.value):q.getOrCreateAnnotationParamTypes(t);let o;for(let t=0;t<e.length;++t){o=e[t];if(o!==void 0){i[t]=o}}}}},transient(e){e.register=function t(i){const s=xe.transient(e,e);return s.register(i)};e.registerInRequestor=false;return e},singleton(e,t=X){e.register=function t(i){const s=xe.singleton(e,e);return s.register(i)};e.registerInRequestor=t.scoped;return e}});const U=q.createInterface("Container");const j=U;function _(e){return function(t){const i=function(e,t,s){q.inject(i)(e,t,s)};i.$isResolver=true;i.resolve=function(i,s){return e(t,i,s)};return i}}const K=q.inject;function W(e){return q.transient(e)}function G(e){return e==null?W:W(e)}const X={scoped:false};function Y(e){return q.singleton(e)}function Q(e){if(typeof e==="function"){return q.singleton(e)}return function(t){return q.singleton(t,e)}}function Z(e){return function(t,i){i=!!i;const s=function(e,t,i){q.inject(s)(e,t,i)};s.$isResolver=true;s.resolve=function(s,o){return e(t,s,o,i)};return s}}const J=Z(((e,t,i,s)=>i.getAll(e,s)));const ee=_(((e,t,i)=>()=>i.get(e)));const te=_(((e,t,i)=>{if(i.has(e,true)){return i.get(e)}else{return undefined}}));function ie(e,t,i){q.inject(ie)(e,t,i)}ie.$isResolver=true;ie.resolve=()=>undefined;const se=_(((e,t,i)=>{const s=ne(e,t);const o=new re(e,0,s);i.registerResolver(e,o);return s}));const oe=_(((e,t,i)=>ne(e,t)));function ne(e,t){return t.getFactory(e).construct(t)}class re{constructor(e,t,i){this.key=e;this.strategy=t;this.state=i;this.resolving=false}get $isResolver(){return true}register(e){return e.registerResolver(this.key,this)}resolve(e,t){switch(this.strategy){case 0:return this.state;case 1:{if(this.resolving){throw new Error(`Cyclic dependency found: ${this.state.name}`)}this.resolving=true;this.state=e.getFactory(this.state).construct(t);this.strategy=0;this.resolving=false;return this.state}case 2:{const i=e.getFactory(this.state);if(i===null){throw new Error(`Resolver for ${String(this.key)} returned a null factory`)}return i.construct(t)}case 3:return this.state(e,t,this);case 4:return this.state[0].resolve(e,t);case 5:return t.get(this.state);default:throw new Error(`Invalid resolver strategy specified: ${this.strategy}.`)}}getFactory(e){var t,i,s;switch(this.strategy){case 1:case 2:return e.getFactory(this.state);case 5:return(s=(i=(t=e.getResolver(this.state))===null||t===void 0?void 0:t.getFactory)===null||i===void 0?void 0:i.call(t,e))!==null&&s!==void 0?s:null;default:return null}}}function ae(e){return this.get(e)}function le(e,t){return t(e)}class de{constructor(e,t){this.Type=e;this.dependencies=t;this.transformers=null}construct(e,t){let i;if(t===void 0){i=new this.Type(...this.dependencies.map(ae,e))}else{i=new this.Type(...this.dependencies.map(ae,e),...t)}if(this.transformers==null){return i}return this.transformers.reduce(le,i)}registerTransformer(e){(this.transformers||(this.transformers=[])).push(e)}}const he={$isResolver:true,resolve(e,t){return t}};function ce(e){return typeof e.register==="function"}function ue(e){return ce(e)&&typeof e.registerInRequestor==="boolean"}function pe(e){return ue(e)&&e.registerInRequestor}function fe(e){return e.prototype!==void 0}const me=new Set(["Array","ArrayBuffer","Boolean","DataView","Date","Error","EvalError","Float32Array","Float64Array","Function","Int8Array","Int16Array","Int32Array","Map","Number","Object","Promise","RangeError","ReferenceError","RegExp","Set","SharedArrayBuffer","String","SyntaxError","TypeError","Uint8Array","Uint8ClampedArray","Uint16Array","Uint32Array","URIError","WeakMap","WeakSet"]);const ve="__DI_LOCATE_PARENT__";const be=new Map;class ge{constructor(e,t){this.owner=e;this.config=t;this._parent=void 0;this.registerDepth=0;this.context=null;if(e!==null){e.$$container$$=this}this.resolvers=new Map;this.resolvers.set(U,he);if(e instanceof Node){e.addEventListener(ve,(e=>{if(e.composedPath()[0]!==this.owner){e.detail.container=this;e.stopImmediatePropagation()}}))}}get parent(){if(this._parent===void 0){this._parent=this.config.parentLocator(this.owner)}return this._parent}get depth(){return this.parent===null?0:this.parent.depth+1}get responsibleForOwnerRequests(){return this.config.responsibleForOwnerRequests}registerWithContext(e,...t){this.context=e;this.register(...t);this.context=null;return this}register(...e){if(++this.registerDepth===100){throw new Error("Unable to autoregister dependency")}let t;let i;let s;let o;let n;const r=this.context;for(let a=0,l=e.length;a<l;++a){t=e[a];if(!ke(t)){continue}if(ce(t)){t.register(this,r)}else if(fe(t)){xe.singleton(t,t).register(this)}else{i=Object.keys(t);o=0;n=i.length;for(;o<n;++o){s=t[i[o]];if(!ke(s)){continue}if(ce(s)){s.register(this,r)}else{this.register(s)}}}}--this.registerDepth;return this}registerResolver(e,t){we(e);const i=this.resolvers;const s=i.get(e);if(s==null){i.set(e,t)}else if(s instanceof re&&s.strategy===4){s.state.push(t)}else{i.set(e,new re(e,4,[s,t]))}return t}registerTransformer(e,t){const i=this.getResolver(e);if(i==null){return false}if(i.getFactory){const e=i.getFactory(this);if(e==null){return false}e.registerTransformer(t);return true}return false}getResolver(e,t=true){we(e);if(e.resolve!==void 0){return e}let i=this;let s;while(i!=null){s=i.resolvers.get(e);if(s==null){if(i.parent==null){const s=pe(e)?this:i;return t?this.jitRegister(e,s):null}i=i.parent}else{return s}}return null}has(e,t=false){return this.resolvers.has(e)?true:t&&this.parent!=null?this.parent.has(e,true):false}get(e){we(e);if(e.$isResolver){return e.resolve(this,this)}let t=this;let i;while(t!=null){i=t.resolvers.get(e);if(i==null){if(t.parent==null){const s=pe(e)?this:t;i=this.jitRegister(e,s);return i.resolve(t,this)}t=t.parent}else{return i.resolve(t,this)}}throw new Error(`Unable to resolve key: ${String(e)}`)}getAll(e,t=false){we(e);const i=this;let o=i;let n;if(t){let t=s.emptyArray;while(o!=null){n=o.resolvers.get(e);if(n!=null){t=t.concat($e(n,o,i))}o=o.parent}return t}else{while(o!=null){n=o.resolvers.get(e);if(n==null){o=o.parent;if(o==null){return s.emptyArray}}else{return $e(n,o,i)}}}return s.emptyArray}getFactory(e){let t=be.get(e);if(t===void 0){if(Oe(e)){throw new Error(`${e.name} is a native function and therefore cannot be safely constructed by DI. If this is intentional, please use a callback or cachedCallback resolver.`)}be.set(e,t=new de(e,q.getDependencies(e)))}return t}registerFactory(e,t){be.set(e,t)}createChild(e){return new ge(null,Object.assign({},this.config,e,{parentLocator:()=>this}))}jitRegister(e,t){if(typeof e!=="function"){throw new Error(`Attempted to jitRegister something that is not a constructor: '${e}'. Did you forget to register this dependency?`)}if(me.has(e.name)){throw new Error(`Attempted to jitRegister an intrinsic type: ${e.name}. Did you forget to add @inject(Key)`)}if(ce(e)){const i=e.register(t);if(!(i instanceof Object)||i.resolve==null){const i=t.resolvers.get(e);if(i!=void 0){return i}throw new Error("A valid resolver was not returned from the static register method")}return i}else if(e.$isInterface){throw new Error(`Attempted to jitRegister an interface: ${e.friendlyName}`)}else{const i=this.config.defaultResolver(e,t);t.resolvers.set(e,i);return i}}}const ye=new WeakMap;function Ce(e){return function(t,i,s){if(ye.has(s)){return ye.get(s)}const o=e(t,i,s);ye.set(s,o);return o}}const xe=Object.freeze({instance(e,t){return new re(e,0,t)},singleton(e,t){return new re(e,1,t)},transient(e,t){return new re(e,2,t)},callback(e,t){return new re(e,3,t)},cachedCallback(e,t){return new re(e,3,Ce(t))},aliasTo(e,t){return new re(t,5,e)}});function we(e){if(e===null||e===void 0){throw new Error("key/value cannot be null or undefined. Are you trying to inject/register something that doesn't exist with DI?")}}function $e(e,t,i){if(e instanceof re&&e.strategy===4){const s=e.state;let o=s.length;const n=new Array(o);while(o--){n[o]=s[o].resolve(t,i)}return n}return[e.resolve(t,i)]}const Ie="(anonymous)";function ke(e){return typeof e==="object"&&e!==null||typeof e==="function"}const Oe=function(){const e=new WeakMap;let t=false;let i="";let s=0;return function(o){t=e.get(o);if(t===void 0){i=o.toString();s=i.length;t=s>=29&&s<=100&&i.charCodeAt(s-1)===125&&i.charCodeAt(s-2)<=32&&i.charCodeAt(s-3)===93&&i.charCodeAt(s-4)===101&&i.charCodeAt(s-5)===100&&i.charCodeAt(s-6)===111&&i.charCodeAt(s-7)===99&&i.charCodeAt(s-8)===32&&i.charCodeAt(s-9)===101&&i.charCodeAt(s-10)===118&&i.charCodeAt(s-11)===105&&i.charCodeAt(s-12)===116&&i.charCodeAt(s-13)===97&&i.charCodeAt(s-14)===110&&i.charCodeAt(s-15)===88;e.set(o,t)}return t}}();const Te={};function Ee(e){switch(typeof e){case"number":return e>=0&&(e|0)===e;case"string":{const t=Te[e];if(t!==void 0){return t}const i=e.length;if(i===0){return Te[e]=false}let s=0;for(let o=0;o<i;++o){s=e.charCodeAt(o);if(o===0&&s===48&&i>1||s<48||s>57){return Te[e]=false}}return Te[e]=true}default:return false}}function Re(e){return`${e.toLowerCase()}:presentation`}const De=new Map;const Se=Object.freeze({define(e,t,i){const s=Re(e);const o=De.get(s);if(o===void 0){De.set(s,t)}else{De.set(s,false)}i.register(xe.instance(s,t))},forTag(e,t){const i=Re(e);const s=De.get(i);if(s===false){const e=q.findResponsibleContainer(t);return e.get(i)}return s||null}});class Ae{constructor(e,t){this.template=e||null;this.styles=t===void 0?null:Array.isArray(t)?s.ElementStyles.create(t):t instanceof s.ElementStyles?t:s.ElementStyles.create([t])}applyTo(e){const t=e.$fastController;if(t.template===null){t.template=this.template}if(t.styles===null){t.styles=this.styles}}}class Fe extends s.FASTElement{constructor(){super(...arguments);this._presentation=void 0}get $presentation(){if(this._presentation===void 0){this._presentation=Se.forTag(this.tagName,this)}return this._presentation}templateChanged(){if(this.template!==undefined){this.$fastController.template=this.template}}stylesChanged(){if(this.styles!==undefined){this.$fastController.styles=this.styles}}connectedCallback(){if(this.$presentation!==null){this.$presentation.applyTo(this)}super.connectedCallback()}static compose(e){return(t={})=>new Me(this===Fe?class extends Fe{}:this,e,t)}}f([s.observable],Fe.prototype,"template",void 0);f([s.observable],Fe.prototype,"styles",void 0);function Le(e,t,i){if(typeof e==="function"){return e(t,i)}return e}class Me{constructor(e,t,i){this.type=e;this.elementDefinition=t;this.overrideDefinition=i;this.definition=Object.assign(Object.assign({},this.elementDefinition),this.overrideDefinition)}register(e,t){const i=this.definition;const s=this.overrideDefinition;const o=i.prefix||t.elementPrefix;const n=`${o}-${i.baseName}`;t.tryDefineElement({name:n,type:this.type,baseClass:this.elementDefinition.baseClass,callback:e=>{const t=new Ae(Le(i.template,e,i),Le(i.styles,e,i));e.definePresentation(t);let o=Le(i.shadowOptions,e,i);if(e.shadowRootMode){if(o){if(!s.shadowOptions){o.mode=e.shadowRootMode}}else if(o!==null){o={mode:e.shadowRootMode}}}e.defineElement({elementOptions:Le(i.elementOptions,e,i),shadowOptions:o,attributes:Le(i.attributes,e,i)})}})}}function Pe(e,...t){const i=s.AttributeConfiguration.locate(e);t.forEach((t=>{Object.getOwnPropertyNames(t.prototype).forEach((i=>{if(i!=="constructor"){Object.defineProperty(e.prototype,i,Object.getOwnPropertyDescriptor(t.prototype,i))}}));const o=s.AttributeConfiguration.locate(t);o.forEach((e=>i.push(e)))}))}class He extends Fe{constructor(){super(...arguments);this.headinglevel=2;this.expanded=false;this.clickHandler=e=>{this.expanded=!this.expanded;this.change()};this.change=()=>{this.$emit("change")}}}f([(0,s.attr)({attribute:"heading-level",mode:"fromView",converter:s.nullableNumberConverter})],He.prototype,"headinglevel",void 0);f([(0,s.attr)({mode:"boolean"})],He.prototype,"expanded",void 0);f([s.attr],He.prototype,"id",void 0);Pe(He,o);const Ve=(e,t)=>(0,s.html)`
    <template>
        <slot ${(0,s.slotted)({property:"accordionItems",filter:(0,s.elements)()})}></slot>
        <slot name="item" part="item" ${(0,s.slotted)("accordionItems")}></slot>
    </template>
`;var ze=i(74291);var Ne=i(83021);const Be={single:"single",multi:"multi"};class qe extends Fe{constructor(){super(...arguments);this.expandmode=Be.multi;this.activeItemIndex=0;this.change=()=>{this.$emit("change",this.activeid)};this.setItems=()=>{var e;if(this.accordionItems.length===0){return}this.accordionIds=this.getItemIds();this.accordionItems.forEach(((e,t)=>{if(e instanceof He){e.addEventListener("change",this.activeItemChange);if(this.isSingleExpandMode()){this.activeItemIndex!==t?e.expanded=false:e.expanded=true}}const i=this.accordionIds[t];e.setAttribute("id",typeof i!=="string"?`accordion-${t+1}`:i);this.activeid=this.accordionIds[this.activeItemIndex];e.addEventListener("keydown",this.handleItemKeyDown);e.addEventListener("focus",this.handleItemFocus)}));if(this.isSingleExpandMode()){const t=(e=this.findExpandedItem())!==null&&e!==void 0?e:this.accordionItems[0];t.setAttribute("aria-disabled","true")}};this.removeItemListeners=e=>{e.forEach(((e,t)=>{e.removeEventListener("change",this.activeItemChange);e.removeEventListener("keydown",this.handleItemKeyDown);e.removeEventListener("focus",this.handleItemFocus)}))};this.activeItemChange=e=>{if(e.defaultPrevented||e.target!==e.currentTarget){return}e.preventDefault();const t=e.target;this.activeid=t.getAttribute("id");if(this.isSingleExpandMode()){this.resetItems();t.expanded=true;t.setAttribute("aria-disabled","true");this.accordionItems.forEach((e=>{if(!e.hasAttribute("disabled")&&e.id!==this.activeid){e.removeAttribute("aria-disabled")}}))}this.activeItemIndex=Array.from(this.accordionItems).indexOf(t);this.change()};this.handleItemKeyDown=e=>{if(e.target!==e.currentTarget){return}this.accordionIds=this.getItemIds();switch(e.key){case ze.I5:e.preventDefault();this.adjust(-1);break;case ze.HX:e.preventDefault();this.adjust(1);break;case ze.Tg:this.activeItemIndex=0;this.focusItem();break;case ze.FM:this.activeItemIndex=this.accordionItems.length-1;this.focusItem();break}};this.handleItemFocus=e=>{if(e.target===e.currentTarget){const t=e.target;const i=this.activeItemIndex=Array.from(this.accordionItems).indexOf(t);if(this.activeItemIndex!==i&&i!==-1){this.activeItemIndex=i;this.activeid=this.accordionIds[this.activeItemIndex]}}}}accordionItemsChanged(e,t){if(this.$fastController.isConnected){this.removeItemListeners(e);this.setItems()}}findExpandedItem(){for(let e=0;e<this.accordionItems.length;e++){if(this.accordionItems[e].getAttribute("expanded")==="true"){return this.accordionItems[e]}}return null}resetItems(){this.accordionItems.forEach(((e,t)=>{e.expanded=false}))}getItemIds(){return this.accordionItems.map((e=>e.getAttribute("id")))}isSingleExpandMode(){return this.expandmode===Be.single}adjust(e){this.activeItemIndex=(0,Ne.Vf)(0,this.accordionItems.length-1,this.activeItemIndex+e);this.focusItem()}focusItem(){const e=this.accordionItems[this.activeItemIndex];if(e instanceof He){e.expandbutton.focus()}}}f([(0,s.attr)({attribute:"expand-mode"})],qe.prototype,"expandmode",void 0);f([s.observable],qe.prototype,"accordionItems",void 0);const Ue=(e,t)=>(0,s.html)`
    <a
        class="control"
        part="control"
        download="${e=>e.download}"
        href="${e=>e.href}"
        hreflang="${e=>e.hreflang}"
        ping="${e=>e.ping}"
        referrerpolicy="${e=>e.referrerpolicy}"
        rel="${e=>e.rel}"
        target="${e=>e.target}"
        type="${e=>e.type}"
        aria-atomic="${e=>e.ariaAtomic}"
        aria-busy="${e=>e.ariaBusy}"
        aria-controls="${e=>e.ariaControls}"
        aria-current="${e=>e.ariaCurrent}"
        aria-describedby="${e=>e.ariaDescribedby}"
        aria-details="${e=>e.ariaDetails}"
        aria-disabled="${e=>e.ariaDisabled}"
        aria-errormessage="${e=>e.ariaErrormessage}"
        aria-expanded="${e=>e.ariaExpanded}"
        aria-flowto="${e=>e.ariaFlowto}"
        aria-haspopup="${e=>e.ariaHaspopup}"
        aria-hidden="${e=>e.ariaHidden}"
        aria-invalid="${e=>e.ariaInvalid}"
        aria-keyshortcuts="${e=>e.ariaKeyshortcuts}"
        aria-label="${e=>e.ariaLabel}"
        aria-labelledby="${e=>e.ariaLabelledby}"
        aria-live="${e=>e.ariaLive}"
        aria-owns="${e=>e.ariaOwns}"
        aria-relevant="${e=>e.ariaRelevant}"
        aria-roledescription="${e=>e.ariaRoledescription}"
        ${(0,s.ref)("control")}
    >
        ${r(e,t)}
        <span class="content" part="content">
            <slot ${(0,s.slotted)("defaultSlottedContent")}></slot>
        </span>
        ${n(e,t)}
    </a>
`;class je{}f([(0,s.attr)({attribute:"aria-atomic"})],je.prototype,"ariaAtomic",void 0);f([(0,s.attr)({attribute:"aria-busy"})],je.prototype,"ariaBusy",void 0);f([(0,s.attr)({attribute:"aria-controls"})],je.prototype,"ariaControls",void 0);f([(0,s.attr)({attribute:"aria-current"})],je.prototype,"ariaCurrent",void 0);f([(0,s.attr)({attribute:"aria-describedby"})],je.prototype,"ariaDescribedby",void 0);f([(0,s.attr)({attribute:"aria-details"})],je.prototype,"ariaDetails",void 0);f([(0,s.attr)({attribute:"aria-disabled"})],je.prototype,"ariaDisabled",void 0);f([(0,s.attr)({attribute:"aria-errormessage"})],je.prototype,"ariaErrormessage",void 0);f([(0,s.attr)({attribute:"aria-flowto"})],je.prototype,"ariaFlowto",void 0);f([(0,s.attr)({attribute:"aria-haspopup"})],je.prototype,"ariaHaspopup",void 0);f([(0,s.attr)({attribute:"aria-hidden"})],je.prototype,"ariaHidden",void 0);f([(0,s.attr)({attribute:"aria-invalid"})],je.prototype,"ariaInvalid",void 0);f([(0,s.attr)({attribute:"aria-keyshortcuts"})],je.prototype,"ariaKeyshortcuts",void 0);f([(0,s.attr)({attribute:"aria-label"})],je.prototype,"ariaLabel",void 0);f([(0,s.attr)({attribute:"aria-labelledby"})],je.prototype,"ariaLabelledby",void 0);f([(0,s.attr)({attribute:"aria-live"})],je.prototype,"ariaLive",void 0);f([(0,s.attr)({attribute:"aria-owns"})],je.prototype,"ariaOwns",void 0);f([(0,s.attr)({attribute:"aria-relevant"})],je.prototype,"ariaRelevant",void 0);f([(0,s.attr)({attribute:"aria-roledescription"})],je.prototype,"ariaRoledescription",void 0);class _e extends Fe{constructor(){super(...arguments);this.handleUnsupportedDelegatesFocus=()=>{var e;if(window.ShadowRoot&&!window.ShadowRoot.prototype.hasOwnProperty("delegatesFocus")&&((e=this.$fastController.definition.shadowOptions)===null||e===void 0?void 0:e.delegatesFocus)){this.focus=()=>{var e;(e=this.control)===null||e===void 0?void 0:e.focus()}}}}connectedCallback(){super.connectedCallback();this.handleUnsupportedDelegatesFocus()}}f([s.attr],_e.prototype,"download",void 0);f([s.attr],_e.prototype,"href",void 0);f([s.attr],_e.prototype,"hreflang",void 0);f([s.attr],_e.prototype,"ping",void 0);f([s.attr],_e.prototype,"referrerpolicy",void 0);f([s.attr],_e.prototype,"rel",void 0);f([s.attr],_e.prototype,"target",void 0);f([s.attr],_e.prototype,"type",void 0);f([s.observable],_e.prototype,"defaultSlottedContent",void 0);class Ke{}f([(0,s.attr)({attribute:"aria-expanded"})],Ke.prototype,"ariaExpanded",void 0);Pe(Ke,je);Pe(_e,o,Ke);const We=(e,t)=>(0,s.html)`
    <template class="${e=>e.initialLayoutComplete?"loaded":""}">
        ${(0,s.when)((e=>e.initialLayoutComplete),(0,s.html)`
                <slot></slot>
            `)}
    </template>
`;var Ge=i(30086);const Xe="abort";const Ye="afterprint";const Qe="animationcancel";const Ze="animationend";const Je="animationiteration";const et="animationstart";const tt="appinstalled";const it="beforeprint";const st="beforeunload";const ot="beginEvent";const nt="blocked";const rt="blur";const at="canplay";const lt="canplaythrough";const dt="change";const ht="chargingchange";const ct="chargingtimechange";const ut="click";const pt="close";const ft="complete";const mt="compositionend";const vt="compositionstart";const bt="compositionupdate";const gt="contextmenu";const yt="copy";const Ct="cut";const xt="dblclick";const wt="devicechange";const $t="devicemotion";const It="deviceorientation";const kt="dischargingtimechange";const Ot="drag";const Tt="dragend";const Et="dragenter";const Rt="dragleave";const Dt="dragover";const St="dragstart";const At="drop";const Ft="durationchange";const Lt="emptied";const Mt="ended";const Pt="endevent";const Ht="error";const Vt="focus";const zt="focusin";const Nt="focusout";const Bt="fullscreenchange";const qt="fullscreenerror";const Ut="gamepadconnected";const jt="gamepaddisconnected";const _t="gotpointercapture";const Kt="hashchange";const Wt="lostpointercapture";const Gt="input";const Xt="invalid";const Yt="keydown";const Qt="keyup";const Zt="levelchange";const Jt="load";const ei="loadeddata";const ti="loadedmetadata";const ii="loadend";const si="loadstart";const oi="message";const ni="messageerror";const ri="mousedown";const ai="mouseenter";const li="mouseleave";const di="mousemove";const hi="mouseout";const ci="mouseover";const ui="mouseup";const pi="notificationclick";const fi="offline";const mi="online";const vi="open";const bi="orientationchange";const gi="pagehide";const yi="pageshow";const Ci="paste";const xi="pause";const wi="pointercancel";const $i="pointerdown";const Ii="pointerenter";const ki="pointerleave";const Oi="pointerlockchange";const Ti="pointerlockerror";const Ei="pointermove";const Ri="pointerout";const Di="pointerover";const Si="pointerup";const Ai="play";const Fi="playing";const Li="popstate";const Mi="progress";const Pi="push";const Hi="pushsubscriptionchange";const Vi="ratechange";const zi="readystatechange";const Ni="repeatevent";const Bi="reset";const qi="resize";const Ui="resourcetimingbufferfull";const ji="scroll";const _i="seeked";const Ki="seeking";const Wi="select";const Gi="show";const Xi="slotchange";const Yi="stalled";const Qi="start";const Zi="storage";const Ji="submit";const es="success";const ts="suspend";const is="SVGAbort";const ss="SVGError";const os="SVGLoad";const ns="SVGResize";const rs="SVGScroll";const as="SVGUnload";const ls="SVGZoom";const ds="timeout";const hs="timeupdate";const cs="touchcancel";const us="touchend";const ps="touchmove";const fs="touchstart";const ms="transitionend";const vs="unload";const bs="upgradeneeded";const gs="userproximity";const ys="versionchange";const Cs="visibilitychange";const xs="volumechange";const ws="waiting";const $s="wheel";const Is=e=>{const t=e.closest("[dir]");return t!==null&&t.dir==="rtl"?Ge.O.rtl:Ge.O.ltr};class ks{constructor(){this.intersectionDetector=null;this.observedElements=new Map;this.requestPosition=(e,t)=>{var i;if(this.intersectionDetector===null){return}if(this.observedElements.has(e)){(i=this.observedElements.get(e))===null||i===void 0?void 0:i.push(t);return}this.observedElements.set(e,[t]);this.intersectionDetector.observe(e)};this.cancelRequestPosition=(e,t)=>{const i=this.observedElements.get(e);if(i!==undefined){const e=i.indexOf(t);if(e!==-1){i.splice(e,1)}}};this.initializeIntersectionDetector=()=>{if(!s.$global.IntersectionObserver){return}this.intersectionDetector=new IntersectionObserver(this.handleIntersection,{root:null,rootMargin:"0px",threshold:[0,1]})};this.handleIntersection=e=>{if(this.intersectionDetector===null){return}const t=[];const i=[];e.forEach((e=>{var s;(s=this.intersectionDetector)===null||s===void 0?void 0:s.unobserve(e.target);const o=this.observedElements.get(e.target);if(o!==undefined){o.forEach((s=>{let o=t.indexOf(s);if(o===-1){o=t.length;t.push(s);i.push([])}i[o].push(e)}));this.observedElements.delete(e.target)}}));t.forEach(((e,t)=>{e(i[t])}))};this.initializeIntersectionDetector()}}class Os extends Fe{constructor(){super(...arguments);this.anchor="";this.viewport="";this.horizontalPositioningMode="uncontrolled";this.horizontalDefaultPosition="unset";this.horizontalViewportLock=false;this.horizontalInset=false;this.horizontalScaling="content";this.verticalPositioningMode="uncontrolled";this.verticalDefaultPosition="unset";this.verticalViewportLock=false;this.verticalInset=false;this.verticalScaling="content";this.fixedPlacement=false;this.autoUpdateMode="anchor";this.anchorElement=null;this.viewportElement=null;this.initialLayoutComplete=false;this.resizeDetector=null;this.baseHorizontalOffset=0;this.baseVerticalOffset=0;this.pendingPositioningUpdate=false;this.pendingReset=false;this.currentDirection=Ge.O.ltr;this.regionVisible=false;this.forceUpdate=false;this.updateThreshold=.5;this.update=()=>{if(!this.pendingPositioningUpdate){this.requestPositionUpdates()}};this.startObservers=()=>{this.stopObservers();if(this.anchorElement===null){return}this.requestPositionUpdates();if(this.resizeDetector!==null){this.resizeDetector.observe(this.anchorElement);this.resizeDetector.observe(this)}};this.requestPositionUpdates=()=>{if(this.anchorElement===null||this.pendingPositioningUpdate){return}Os.intersectionService.requestPosition(this,this.handleIntersection);Os.intersectionService.requestPosition(this.anchorElement,this.handleIntersection);if(this.viewportElement!==null){Os.intersectionService.requestPosition(this.viewportElement,this.handleIntersection)}this.pendingPositioningUpdate=true};this.stopObservers=()=>{if(this.pendingPositioningUpdate){this.pendingPositioningUpdate=false;Os.intersectionService.cancelRequestPosition(this,this.handleIntersection);if(this.anchorElement!==null){Os.intersectionService.cancelRequestPosition(this.anchorElement,this.handleIntersection)}if(this.viewportElement!==null){Os.intersectionService.cancelRequestPosition(this.viewportElement,this.handleIntersection)}}if(this.resizeDetector!==null){this.resizeDetector.disconnect()}};this.getViewport=()=>{if(typeof this.viewport!=="string"||this.viewport===""){return document.documentElement}return document.getElementById(this.viewport)};this.getAnchor=()=>document.getElementById(this.anchor);this.handleIntersection=e=>{if(!this.pendingPositioningUpdate){return}this.pendingPositioningUpdate=false;if(!this.applyIntersectionEntries(e)){return}this.updateLayout()};this.applyIntersectionEntries=e=>{const t=e.find((e=>e.target===this));const i=e.find((e=>e.target===this.anchorElement));const s=e.find((e=>e.target===this.viewportElement));if(t===undefined||s===undefined||i===undefined){return false}if(!this.regionVisible||this.forceUpdate||this.regionRect===undefined||this.anchorRect===undefined||this.viewportRect===undefined||this.isRectDifferent(this.anchorRect,i.boundingClientRect)||this.isRectDifferent(this.viewportRect,s.boundingClientRect)||this.isRectDifferent(this.regionRect,t.boundingClientRect)){this.regionRect=t.boundingClientRect;this.anchorRect=i.boundingClientRect;if(this.viewportElement===document.documentElement){this.viewportRect=new DOMRectReadOnly(s.boundingClientRect.x+document.documentElement.scrollLeft,s.boundingClientRect.y+document.documentElement.scrollTop,s.boundingClientRect.width,s.boundingClientRect.height)}else{this.viewportRect=s.boundingClientRect}this.updateRegionOffset();this.forceUpdate=false;return true}return false};this.updateRegionOffset=()=>{if(this.anchorRect&&this.regionRect){this.baseHorizontalOffset=this.baseHorizontalOffset+(this.anchorRect.left-this.regionRect.left)+(this.translateX-this.baseHorizontalOffset);this.baseVerticalOffset=this.baseVerticalOffset+(this.anchorRect.top-this.regionRect.top)+(this.translateY-this.baseVerticalOffset)}};this.isRectDifferent=(e,t)=>{if(Math.abs(e.top-t.top)>this.updateThreshold||Math.abs(e.right-t.right)>this.updateThreshold||Math.abs(e.bottom-t.bottom)>this.updateThreshold||Math.abs(e.left-t.left)>this.updateThreshold){return true}return false};this.handleResize=e=>{this.update()};this.reset=()=>{if(!this.pendingReset){return}this.pendingReset=false;if(this.anchorElement===null){this.anchorElement=this.getAnchor()}if(this.viewportElement===null){this.viewportElement=this.getViewport()}this.currentDirection=Is(this);this.startObservers()};this.updateLayout=()=>{let e=undefined;let t=undefined;if(this.horizontalPositioningMode!=="uncontrolled"){const e=this.getPositioningOptions(this.horizontalInset);if(this.horizontalDefaultPosition==="center"){t="center"}else if(this.horizontalDefaultPosition!=="unset"){let e=this.horizontalDefaultPosition;if(e==="start"||e==="end"){const t=Is(this);if(t!==this.currentDirection){this.currentDirection=t;this.initialize();return}if(this.currentDirection===Ge.O.ltr){e=e==="start"?"left":"right"}else{e=e==="start"?"right":"left"}}switch(e){case"left":t=this.horizontalInset?"insetStart":"start";break;case"right":t=this.horizontalInset?"insetEnd":"end";break}}const i=this.horizontalThreshold!==undefined?this.horizontalThreshold:this.regionRect!==undefined?this.regionRect.width:0;const s=this.anchorRect!==undefined?this.anchorRect.left:0;const o=this.anchorRect!==undefined?this.anchorRect.right:0;const n=this.anchorRect!==undefined?this.anchorRect.width:0;const r=this.viewportRect!==undefined?this.viewportRect.left:0;const a=this.viewportRect!==undefined?this.viewportRect.right:0;if(t===undefined||!(this.horizontalPositioningMode==="locktodefault")&&this.getAvailableSpace(t,s,o,n,r,a)<i){t=this.getAvailableSpace(e[0],s,o,n,r,a)>this.getAvailableSpace(e[1],s,o,n,r,a)?e[0]:e[1]}}if(this.verticalPositioningMode!=="uncontrolled"){const t=this.getPositioningOptions(this.verticalInset);if(this.verticalDefaultPosition==="center"){e="center"}else if(this.verticalDefaultPosition!=="unset"){switch(this.verticalDefaultPosition){case"top":e=this.verticalInset?"insetStart":"start";break;case"bottom":e=this.verticalInset?"insetEnd":"end";break}}const i=this.verticalThreshold!==undefined?this.verticalThreshold:this.regionRect!==undefined?this.regionRect.height:0;const s=this.anchorRect!==undefined?this.anchorRect.top:0;const o=this.anchorRect!==undefined?this.anchorRect.bottom:0;const n=this.anchorRect!==undefined?this.anchorRect.height:0;const r=this.viewportRect!==undefined?this.viewportRect.top:0;const a=this.viewportRect!==undefined?this.viewportRect.bottom:0;if(e===undefined||!(this.verticalPositioningMode==="locktodefault")&&this.getAvailableSpace(e,s,o,n,r,a)<i){e=this.getAvailableSpace(t[0],s,o,n,r,a)>this.getAvailableSpace(t[1],s,o,n,r,a)?t[0]:t[1]}}const i=this.getNextRegionDimension(t,e);const s=this.horizontalPosition!==t||this.verticalPosition!==e;this.setHorizontalPosition(t,i);this.setVerticalPosition(e,i);this.updateRegionStyle();if(!this.initialLayoutComplete){this.initialLayoutComplete=true;this.requestPositionUpdates();return}if(!this.regionVisible){this.regionVisible=true;this.style.removeProperty("pointer-events");this.style.removeProperty("opacity");this.classList.toggle("loaded",true);this.$emit("loaded",this,{bubbles:false})}this.updatePositionClasses();if(s){this.$emit("positionchange",this,{bubbles:false})}};this.updateRegionStyle=()=>{this.style.width=this.regionWidth;this.style.height=this.regionHeight;this.style.transform=`translate(${this.translateX}px, ${this.translateY}px)`};this.updatePositionClasses=()=>{this.classList.toggle("top",this.verticalPosition==="start");this.classList.toggle("bottom",this.verticalPosition==="end");this.classList.toggle("inset-top",this.verticalPosition==="insetStart");this.classList.toggle("inset-bottom",this.verticalPosition==="insetEnd");this.classList.toggle("vertical-center",this.verticalPosition==="center");this.classList.toggle("left",this.horizontalPosition==="start");this.classList.toggle("right",this.horizontalPosition==="end");this.classList.toggle("inset-left",this.horizontalPosition==="insetStart");this.classList.toggle("inset-right",this.horizontalPosition==="insetEnd");this.classList.toggle("horizontal-center",this.horizontalPosition==="center")};this.setHorizontalPosition=(e,t)=>{if(e===undefined||this.regionRect===undefined||this.anchorRect===undefined||this.viewportRect===undefined){return}let i=0;switch(this.horizontalScaling){case"anchor":case"fill":i=this.horizontalViewportLock?this.viewportRect.width:t.width;this.regionWidth=`${i}px`;break;case"content":i=this.regionRect.width;this.regionWidth="unset";break}let s=0;switch(e){case"start":this.translateX=this.baseHorizontalOffset-i;if(this.horizontalViewportLock&&this.anchorRect.left>this.viewportRect.right){this.translateX=this.translateX-(this.anchorRect.left-this.viewportRect.right)}break;case"insetStart":this.translateX=this.baseHorizontalOffset-i+this.anchorRect.width;if(this.horizontalViewportLock&&this.anchorRect.right>this.viewportRect.right){this.translateX=this.translateX-(this.anchorRect.right-this.viewportRect.right)}break;case"insetEnd":this.translateX=this.baseHorizontalOffset;if(this.horizontalViewportLock&&this.anchorRect.left<this.viewportRect.left){this.translateX=this.translateX-(this.anchorRect.left-this.viewportRect.left)}break;case"end":this.translateX=this.baseHorizontalOffset+this.anchorRect.width;if(this.horizontalViewportLock&&this.anchorRect.right<this.viewportRect.left){this.translateX=this.translateX-(this.anchorRect.right-this.viewportRect.left)}break;case"center":s=(this.anchorRect.width-i)/2;this.translateX=this.baseHorizontalOffset+s;if(this.horizontalViewportLock){const e=this.anchorRect.left+s;const t=this.anchorRect.right-s;if(e<this.viewportRect.left&&!(t>this.viewportRect.right)){this.translateX=this.translateX-(e-this.viewportRect.left)}else if(t>this.viewportRect.right&&!(e<this.viewportRect.left)){this.translateX=this.translateX-(t-this.viewportRect.right)}}break}this.horizontalPosition=e};this.setVerticalPosition=(e,t)=>{if(e===undefined||this.regionRect===undefined||this.anchorRect===undefined||this.viewportRect===undefined){return}let i=0;switch(this.verticalScaling){case"anchor":case"fill":i=this.verticalViewportLock?this.viewportRect.height:t.height;this.regionHeight=`${i}px`;break;case"content":i=this.regionRect.height;this.regionHeight="unset";break}let s=0;switch(e){case"start":this.translateY=this.baseVerticalOffset-i;if(this.verticalViewportLock&&this.anchorRect.top>this.viewportRect.bottom){this.translateY=this.translateY-(this.anchorRect.top-this.viewportRect.bottom)}break;case"insetStart":this.translateY=this.baseVerticalOffset-i+this.anchorRect.height;if(this.verticalViewportLock&&this.anchorRect.bottom>this.viewportRect.bottom){this.translateY=this.translateY-(this.anchorRect.bottom-this.viewportRect.bottom)}break;case"insetEnd":this.translateY=this.baseVerticalOffset;if(this.verticalViewportLock&&this.anchorRect.top<this.viewportRect.top){this.translateY=this.translateY-(this.anchorRect.top-this.viewportRect.top)}break;case"end":this.translateY=this.baseVerticalOffset+this.anchorRect.height;if(this.verticalViewportLock&&this.anchorRect.bottom<this.viewportRect.top){this.translateY=this.translateY-(this.anchorRect.bottom-this.viewportRect.top)}break;case"center":s=(this.anchorRect.height-i)/2;this.translateY=this.baseVerticalOffset+s;if(this.verticalViewportLock){const e=this.anchorRect.top+s;const t=this.anchorRect.bottom-s;if(e<this.viewportRect.top&&!(t>this.viewportRect.bottom)){this.translateY=this.translateY-(e-this.viewportRect.top)}else if(t>this.viewportRect.bottom&&!(e<this.viewportRect.top)){this.translateY=this.translateY-(t-this.viewportRect.bottom)}}}this.verticalPosition=e};this.getPositioningOptions=e=>{if(e){return["insetStart","insetEnd"]}return["start","end"]};this.getAvailableSpace=(e,t,i,s,o,n)=>{const r=t-o;const a=n-(t+s);switch(e){case"start":return r;case"insetStart":return r+s;case"insetEnd":return a+s;case"end":return a;case"center":return Math.min(r,a)*2+s}};this.getNextRegionDimension=(e,t)=>{const i={height:this.regionRect!==undefined?this.regionRect.height:0,width:this.regionRect!==undefined?this.regionRect.width:0};if(e!==undefined&&this.horizontalScaling==="fill"){i.width=this.getAvailableSpace(e,this.anchorRect!==undefined?this.anchorRect.left:0,this.anchorRect!==undefined?this.anchorRect.right:0,this.anchorRect!==undefined?this.anchorRect.width:0,this.viewportRect!==undefined?this.viewportRect.left:0,this.viewportRect!==undefined?this.viewportRect.right:0)}else if(this.horizontalScaling==="anchor"){i.width=this.anchorRect!==undefined?this.anchorRect.width:0}if(t!==undefined&&this.verticalScaling==="fill"){i.height=this.getAvailableSpace(t,this.anchorRect!==undefined?this.anchorRect.top:0,this.anchorRect!==undefined?this.anchorRect.bottom:0,this.anchorRect!==undefined?this.anchorRect.height:0,this.viewportRect!==undefined?this.viewportRect.top:0,this.viewportRect!==undefined?this.viewportRect.bottom:0)}else if(this.verticalScaling==="anchor"){i.height=this.anchorRect!==undefined?this.anchorRect.height:0}return i};this.startAutoUpdateEventListeners=()=>{window.addEventListener(qi,this.update,{passive:true});window.addEventListener(ji,this.update,{passive:true,capture:true});if(this.resizeDetector!==null&&this.viewportElement!==null){this.resizeDetector.observe(this.viewportElement)}};this.stopAutoUpdateEventListeners=()=>{window.removeEventListener(qi,this.update);window.removeEventListener(ji,this.update);if(this.resizeDetector!==null&&this.viewportElement!==null){this.resizeDetector.unobserve(this.viewportElement)}}}anchorChanged(){if(this.initialLayoutComplete){this.anchorElement=this.getAnchor()}}viewportChanged(){if(this.initialLayoutComplete){this.viewportElement=this.getViewport()}}horizontalPositioningModeChanged(){this.requestReset()}horizontalDefaultPositionChanged(){this.updateForAttributeChange()}horizontalViewportLockChanged(){this.updateForAttributeChange()}horizontalInsetChanged(){this.updateForAttributeChange()}horizontalThresholdChanged(){this.updateForAttributeChange()}horizontalScalingChanged(){this.updateForAttributeChange()}verticalPositioningModeChanged(){this.requestReset()}verticalDefaultPositionChanged(){this.updateForAttributeChange()}verticalViewportLockChanged(){this.updateForAttributeChange()}verticalInsetChanged(){this.updateForAttributeChange()}verticalThresholdChanged(){this.updateForAttributeChange()}verticalScalingChanged(){this.updateForAttributeChange()}fixedPlacementChanged(){if(this.$fastController.isConnected&&this.initialLayoutComplete){this.initialize()}}autoUpdateModeChanged(e,t){if(this.$fastController.isConnected&&this.initialLayoutComplete){if(e==="auto"){this.stopAutoUpdateEventListeners()}if(t==="auto"){this.startAutoUpdateEventListeners()}}}anchorElementChanged(){this.requestReset()}viewportElementChanged(){if(this.$fastController.isConnected&&this.initialLayoutComplete){this.initialize()}}connectedCallback(){super.connectedCallback();if(this.autoUpdateMode==="auto"){this.startAutoUpdateEventListeners()}this.initialize()}disconnectedCallback(){super.disconnectedCallback();if(this.autoUpdateMode==="auto"){this.stopAutoUpdateEventListeners()}this.stopObservers();this.disconnectResizeDetector()}adoptedCallback(){this.initialize()}disconnectResizeDetector(){if(this.resizeDetector!==null){this.resizeDetector.disconnect();this.resizeDetector=null}}initializeResizeDetector(){this.disconnectResizeDetector();this.resizeDetector=new window.ResizeObserver(this.handleResize)}updateForAttributeChange(){if(this.$fastController.isConnected&&this.initialLayoutComplete){this.forceUpdate=true;this.update()}}initialize(){this.initializeResizeDetector();if(this.anchorElement===null){this.anchorElement=this.getAnchor()}this.requestReset()}requestReset(){if(this.$fastController.isConnected&&this.pendingReset===false){this.setInitialState();s.DOM.queueUpdate((()=>this.reset()));this.pendingReset=true}}setInitialState(){this.initialLayoutComplete=false;this.regionVisible=false;this.translateX=0;this.translateY=0;this.baseHorizontalOffset=0;this.baseVerticalOffset=0;this.viewportRect=undefined;this.regionRect=undefined;this.anchorRect=undefined;this.verticalPosition=undefined;this.horizontalPosition=undefined;this.style.opacity="0";this.style.pointerEvents="none";this.forceUpdate=false;this.style.position=this.fixedPlacement?"fixed":"absolute";this.updatePositionClasses();this.updateRegionStyle()}}Os.intersectionService=new ks;f([s.attr],Os.prototype,"anchor",void 0);f([s.attr],Os.prototype,"viewport",void 0);f([(0,s.attr)({attribute:"horizontal-positioning-mode"})],Os.prototype,"horizontalPositioningMode",void 0);f([(0,s.attr)({attribute:"horizontal-default-position"})],Os.prototype,"horizontalDefaultPosition",void 0);f([(0,s.attr)({attribute:"horizontal-viewport-lock",mode:"boolean"})],Os.prototype,"horizontalViewportLock",void 0);f([(0,s.attr)({attribute:"horizontal-inset",mode:"boolean"})],Os.prototype,"horizontalInset",void 0);f([(0,s.attr)({attribute:"horizontal-threshold"})],Os.prototype,"horizontalThreshold",void 0);f([(0,s.attr)({attribute:"horizontal-scaling"})],Os.prototype,"horizontalScaling",void 0);f([(0,s.attr)({attribute:"vertical-positioning-mode"})],Os.prototype,"verticalPositioningMode",void 0);f([(0,s.attr)({attribute:"vertical-default-position"})],Os.prototype,"verticalDefaultPosition",void 0);f([(0,s.attr)({attribute:"vertical-viewport-lock",mode:"boolean"})],Os.prototype,"verticalViewportLock",void 0);f([(0,s.attr)({attribute:"vertical-inset",mode:"boolean"})],Os.prototype,"verticalInset",void 0);f([(0,s.attr)({attribute:"vertical-threshold"})],Os.prototype,"verticalThreshold",void 0);f([(0,s.attr)({attribute:"vertical-scaling"})],Os.prototype,"verticalScaling",void 0);f([(0,s.attr)({attribute:"fixed-placement",mode:"boolean"})],Os.prototype,"fixedPlacement",void 0);f([(0,s.attr)({attribute:"auto-update-mode"})],Os.prototype,"autoUpdateMode",void 0);f([s.observable],Os.prototype,"anchorElement",void 0);f([s.observable],Os.prototype,"viewportElement",void 0);f([s.observable],Os.prototype,"initialLayoutComplete",void 0);const Ts={horizontalDefaultPosition:"center",horizontalPositioningMode:"locktodefault",horizontalInset:false,horizontalScaling:"anchor"};const Es=Object.assign(Object.assign({},Ts),{verticalDefaultPosition:"top",verticalPositioningMode:"locktodefault",verticalInset:false,verticalScaling:"content"});const Rs=Object.assign(Object.assign({},Ts),{verticalDefaultPosition:"bottom",verticalPositioningMode:"locktodefault",verticalInset:false,verticalScaling:"content"});const Ds=Object.assign(Object.assign({},Ts),{verticalPositioningMode:"dynamic",verticalInset:false,verticalScaling:"content"});const Ss=Object.assign(Object.assign({},Es),{verticalScaling:"fill"});const As=Object.assign(Object.assign({},Rs),{verticalScaling:"fill"});const Fs=Object.assign(Object.assign({},Ds),{verticalScaling:"fill"});const Ls=(e,t)=>(0,s.html)`
    <div
        class="backplate ${e=>e.shape}"
        part="backplate"
        style="${e=>e.fill?`background-color: var(--avatar-fill-${e.fill});`:void 0}"
    >
        <a
            class="link"
            part="link"
            href="${e=>e.link?e.link:void 0}"
            style="${e=>e.color?`color: var(--avatar-color-${e.color});`:void 0}"
        >
            <slot name="media" part="media">${t.media||""}</slot>
            <slot class="content" part="content"><slot>
        </a>
    </div>
    <slot name="badge" part="badge"></slot>
`;class Ms extends Fe{connectedCallback(){super.connectedCallback();if(!this.shape){this.shape="circle"}}}f([s.attr],Ms.prototype,"fill",void 0);f([s.attr],Ms.prototype,"color",void 0);f([s.attr],Ms.prototype,"link",void 0);f([s.attr],Ms.prototype,"shape",void 0);const Ps=(e,t)=>(0,s.html)`
    <template class="${e=>e.circular?"circular":""}">
        <div class="control" part="control" style="${e=>e.generateBadgeStyle()}">
            <slot></slot>
        </div>
    </template>
`;class Hs extends Fe{constructor(){super(...arguments);this.generateBadgeStyle=()=>{if(!this.fill&&!this.color){return}const e=`background-color: var(--badge-fill-${this.fill});`;const t=`color: var(--badge-color-${this.color});`;if(this.fill&&!this.color){return e}else if(this.color&&!this.fill){return t}else{return`${t} ${e}`}}}}f([(0,s.attr)({attribute:"fill"})],Hs.prototype,"fill",void 0);f([(0,s.attr)({attribute:"color"})],Hs.prototype,"color",void 0);f([(0,s.attr)({mode:"boolean"})],Hs.prototype,"circular",void 0);const Vs=(e,t)=>(0,s.html)`
    <div role="listitem" class="listitem" part="listitem">
        ${(0,s.when)((e=>e.href&&e.href.length>0),(0,s.html)`
                ${Ue(e,t)}
            `)}
        ${(0,s.when)((e=>!e.href),(0,s.html)`
                ${r(e,t)}
                <slot></slot>
                ${n(e,t)}
            `)}
        ${(0,s.when)((e=>e.separator),(0,s.html)`
                <span class="separator" part="separator" aria-hidden="true">
                    <slot name="separator">${t.separator||""}</slot>
                </span>
            `)}
    </div>
`;class zs extends _e{constructor(){super(...arguments);this.separator=true}}f([s.observable],zs.prototype,"separator",void 0);Pe(zs,o,Ke);const Ns=(e,t)=>(0,s.html)`
    <template role="navigation">
        <div role="list" class="list" part="list">
            <slot
                ${(0,s.slotted)({property:"slottedBreadcrumbItems",filter:(0,s.elements)()})}
            ></slot>
        </div>
    </template>
`;class Bs extends Fe{slottedBreadcrumbItemsChanged(){if(this.$fastController.isConnected){if(this.slottedBreadcrumbItems===undefined||this.slottedBreadcrumbItems.length===0){return}const e=this.slottedBreadcrumbItems[this.slottedBreadcrumbItems.length-1];this.slottedBreadcrumbItems.forEach((t=>{const i=t===e;this.setItemSeparator(t,i);this.setAriaCurrent(t,i)}))}}setItemSeparator(e,t){if(e instanceof zs){e.separator=!t}}findChildWithHref(e){var t,i;if(e.childElementCount>0){return e.querySelector("a[href]")}else if((t=e.shadowRoot)===null||t===void 0?void 0:t.childElementCount){return(i=e.shadowRoot)===null||i===void 0?void 0:i.querySelector("a[href]")}else return null}setAriaCurrent(e,t){const i=this.findChildWithHref(e);if(i===null&&e.hasAttribute("href")&&e instanceof zs){t?e.setAttribute("aria-current","page"):e.removeAttribute("aria-current")}else if(i!==null){t?i.setAttribute("aria-current","page"):i.removeAttribute("aria-current")}}}f([s.observable],Bs.prototype,"slottedBreadcrumbItems",void 0);const qs=(e,t)=>(0,s.html)`
    <button
        class="control"
        part="control"
        ?autofocus="${e=>e.autofocus}"
        ?disabled="${e=>e.disabled}"
        form="${e=>e.formId}"
        formaction="${e=>e.formaction}"
        formenctype="${e=>e.formenctype}"
        formmethod="${e=>e.formmethod}"
        formnovalidate="${e=>e.formnovalidate}"
        formtarget="${e=>e.formtarget}"
        name="${e=>e.name}"
        type="${e=>e.type}"
        value="${e=>e.value}"
        aria-atomic="${e=>e.ariaAtomic}"
        aria-busy="${e=>e.ariaBusy}"
        aria-controls="${e=>e.ariaControls}"
        aria-current="${e=>e.ariaCurrent}"
        aria-describedby="${e=>e.ariaDescribedby}"
        aria-details="${e=>e.ariaDetails}"
        aria-disabled="${e=>e.ariaDisabled}"
        aria-errormessage="${e=>e.ariaErrormessage}"
        aria-expanded="${e=>e.ariaExpanded}"
        aria-flowto="${e=>e.ariaFlowto}"
        aria-haspopup="${e=>e.ariaHaspopup}"
        aria-hidden="${e=>e.ariaHidden}"
        aria-invalid="${e=>e.ariaInvalid}"
        aria-keyshortcuts="${e=>e.ariaKeyshortcuts}"
        aria-label="${e=>e.ariaLabel}"
        aria-labelledby="${e=>e.ariaLabelledby}"
        aria-live="${e=>e.ariaLive}"
        aria-owns="${e=>e.ariaOwns}"
        aria-pressed="${e=>e.ariaPressed}"
        aria-relevant="${e=>e.ariaRelevant}"
        aria-roledescription="${e=>e.ariaRoledescription}"
        ${(0,s.ref)("control")}
    >
        ${r(e,t)}
        <span class="content" part="content">
            <slot ${(0,s.slotted)("defaultSlottedContent")}></slot>
        </span>
        ${n(e,t)}
    </button>
`;const Us="form-associated-proxy";const js="ElementInternals";const _s=js in window&&"setFormValue"in window[js].prototype;const Ks=new WeakMap;function Ws(e){const t=class extends e{constructor(...e){super(...e);this.dirtyValue=false;this.disabled=false;this.proxyEventsToBlock=["change","click"];this.proxyInitialized=false;this.required=false;this.initialValue=this.initialValue||"";if(!this.elementInternals){this.formResetCallback=this.formResetCallback.bind(this)}}static get formAssociated(){return _s}get validity(){return this.elementInternals?this.elementInternals.validity:this.proxy.validity}get form(){return this.elementInternals?this.elementInternals.form:this.proxy.form}get validationMessage(){return this.elementInternals?this.elementInternals.validationMessage:this.proxy.validationMessage}get willValidate(){return this.elementInternals?this.elementInternals.willValidate:this.proxy.willValidate}get labels(){if(this.elementInternals){return Object.freeze(Array.from(this.elementInternals.labels))}else if(this.proxy instanceof HTMLElement&&this.proxy.ownerDocument&&this.id){const e=this.proxy.labels;const t=Array.from(this.proxy.getRootNode().querySelectorAll(`[for='${this.id}']`));const i=e?t.concat(Array.from(e)):t;return Object.freeze(i)}else{return s.emptyArray}}valueChanged(e,t){this.dirtyValue=true;if(this.proxy instanceof HTMLElement){this.proxy.value=this.value}this.currentValue=this.value;this.setFormValue(this.value);this.validate()}currentValueChanged(){this.value=this.currentValue}initialValueChanged(e,t){if(!this.dirtyValue){this.value=this.initialValue;this.dirtyValue=false}}disabledChanged(e,t){if(this.proxy instanceof HTMLElement){this.proxy.disabled=this.disabled}s.DOM.queueUpdate((()=>this.classList.toggle("disabled",this.disabled)))}nameChanged(e,t){if(this.proxy instanceof HTMLElement){this.proxy.name=this.name}}requiredChanged(e,t){if(this.proxy instanceof HTMLElement){this.proxy.required=this.required}s.DOM.queueUpdate((()=>this.classList.toggle("required",this.required)));this.validate()}get elementInternals(){if(!_s){return null}let e=Ks.get(this);if(!e){e=this.attachInternals();Ks.set(this,e)}return e}connectedCallback(){super.connectedCallback();this.addEventListener("keypress",this._keypressHandler);if(!this.value){this.value=this.initialValue;this.dirtyValue=false}if(!this.elementInternals){this.attachProxy();if(this.form){this.form.addEventListener("reset",this.formResetCallback)}}}disconnectedCallback(){super.disconnectedCallback();this.proxyEventsToBlock.forEach((e=>this.proxy.removeEventListener(e,this.stopPropagation)));if(!this.elementInternals&&this.form){this.form.removeEventListener("reset",this.formResetCallback)}}checkValidity(){return this.elementInternals?this.elementInternals.checkValidity():this.proxy.checkValidity()}reportValidity(){return this.elementInternals?this.elementInternals.reportValidity():this.proxy.reportValidity()}setValidity(e,t,i){if(this.elementInternals){this.elementInternals.setValidity(e,t,i)}else if(typeof t==="string"){this.proxy.setCustomValidity(t)}}formDisabledCallback(e){this.disabled=e}formResetCallback(){this.value=this.initialValue;this.dirtyValue=false}attachProxy(){var e;if(!this.proxyInitialized){this.proxyInitialized=true;this.proxy.style.display="none";this.proxyEventsToBlock.forEach((e=>this.proxy.addEventListener(e,this.stopPropagation)));this.proxy.disabled=this.disabled;this.proxy.required=this.required;if(typeof this.name==="string"){this.proxy.name=this.name}if(typeof this.value==="string"){this.proxy.value=this.value}this.proxy.setAttribute("slot",Us);this.proxySlot=document.createElement("slot");this.proxySlot.setAttribute("name",Us)}(e=this.shadowRoot)===null||e===void 0?void 0:e.appendChild(this.proxySlot);this.appendChild(this.proxy)}detachProxy(){var e;this.removeChild(this.proxy);(e=this.shadowRoot)===null||e===void 0?void 0:e.removeChild(this.proxySlot)}validate(e){if(this.proxy instanceof HTMLElement){this.setValidity(this.proxy.validity,this.proxy.validationMessage,e)}}setFormValue(e,t){if(this.elementInternals){this.elementInternals.setFormValue(e,t||e)}}_keypressHandler(e){switch(e.key){case ze.Mm:if(this.form instanceof HTMLFormElement){const e=this.form.querySelector("[type=submit]");e===null||e===void 0?void 0:e.click()}break}}stopPropagation(e){e.stopPropagation()}};(0,s.attr)({mode:"boolean"})(t.prototype,"disabled");(0,s.attr)({mode:"fromView",attribute:"value"})(t.prototype,"initialValue");(0,s.attr)({attribute:"current-value"})(t.prototype,"currentValue");(0,s.attr)(t.prototype,"name");(0,s.attr)({mode:"boolean"})(t.prototype,"required");(0,s.observable)(t.prototype,"value");return t}function Gs(e){class t extends(Ws(e)){}class i extends t{constructor(...e){super(e);this.dirtyChecked=false;this.checkedAttribute=false;this.checked=false;this.dirtyChecked=false}checkedAttributeChanged(){this.defaultChecked=this.checkedAttribute}defaultCheckedChanged(){if(!this.dirtyChecked){this.checked=this.defaultChecked;this.dirtyChecked=false}}checkedChanged(e,t){if(!this.dirtyChecked){this.dirtyChecked=true}this.currentChecked=this.checked;this.updateForm();if(this.proxy instanceof HTMLInputElement){this.proxy.checked=this.checked}if(e!==undefined){this.$emit("change")}this.validate()}currentCheckedChanged(e,t){this.checked=this.currentChecked}updateForm(){const e=this.checked?this.value:null;this.setFormValue(e,e)}connectedCallback(){super.connectedCallback();this.updateForm()}formResetCallback(){super.formResetCallback();this.checked=!!this.checkedAttribute;this.dirtyChecked=false}}(0,s.attr)({attribute:"checked",mode:"boolean"})(i.prototype,"checkedAttribute");(0,s.attr)({attribute:"current-checked",converter:s.booleanConverter})(i.prototype,"currentChecked");(0,s.observable)(i.prototype,"defaultChecked");(0,s.observable)(i.prototype,"checked");return i}class Xs extends Fe{}class Ys extends(Ws(Xs)){constructor(){super(...arguments);this.proxy=document.createElement("input")}}class Qs extends Ys{constructor(){super(...arguments);this.handleClick=e=>{var t;if(this.disabled&&((t=this.defaultSlottedContent)===null||t===void 0?void 0:t.length)<=1){e.stopPropagation()}};this.handleSubmission=()=>{if(!this.form){return}const e=this.proxy.isConnected;if(!e){this.attachProxy()}typeof this.form.requestSubmit==="function"?this.form.requestSubmit(this.proxy):this.proxy.click();if(!e){this.detachProxy()}};this.handleFormReset=()=>{var e;(e=this.form)===null||e===void 0?void 0:e.reset()};this.handleUnsupportedDelegatesFocus=()=>{var e;if(window.ShadowRoot&&!window.ShadowRoot.prototype.hasOwnProperty("delegatesFocus")&&((e=this.$fastController.definition.shadowOptions)===null||e===void 0?void 0:e.delegatesFocus)){this.focus=()=>{this.control.focus()}}}}formactionChanged(){if(this.proxy instanceof HTMLInputElement){this.proxy.formAction=this.formaction}}formenctypeChanged(){if(this.proxy instanceof HTMLInputElement){this.proxy.formEnctype=this.formenctype}}formmethodChanged(){if(this.proxy instanceof HTMLInputElement){this.proxy.formMethod=this.formmethod}}formnovalidateChanged(){if(this.proxy instanceof HTMLInputElement){this.proxy.formNoValidate=this.formnovalidate}}formtargetChanged(){if(this.proxy instanceof HTMLInputElement){this.proxy.formTarget=this.formtarget}}typeChanged(e,t){if(this.proxy instanceof HTMLInputElement){this.proxy.type=this.type}t==="submit"&&this.addEventListener("click",this.handleSubmission);e==="submit"&&this.removeEventListener("click",this.handleSubmission);t==="reset"&&this.addEventListener("click",this.handleFormReset);e==="reset"&&this.removeEventListener("click",this.handleFormReset)}validate(){super.validate(this.control)}connectedCallback(){var e;super.connectedCallback();this.proxy.setAttribute("type",this.type);this.handleUnsupportedDelegatesFocus();const t=Array.from((e=this.control)===null||e===void 0?void 0:e.children);if(t){t.forEach((e=>{e.addEventListener("click",this.handleClick)}))}}disconnectedCallback(){var e;super.disconnectedCallback();const t=Array.from((e=this.control)===null||e===void 0?void 0:e.children);if(t){t.forEach((e=>{e.removeEventListener("click",this.handleClick)}))}}}f([(0,s.attr)({mode:"boolean"})],Qs.prototype,"autofocus",void 0);f([(0,s.attr)({attribute:"form"})],Qs.prototype,"formId",void 0);f([s.attr],Qs.prototype,"formaction",void 0);f([s.attr],Qs.prototype,"formenctype",void 0);f([s.attr],Qs.prototype,"formmethod",void 0);f([(0,s.attr)({mode:"boolean"})],Qs.prototype,"formnovalidate",void 0);f([s.attr],Qs.prototype,"formtarget",void 0);f([s.attr],Qs.prototype,"type",void 0);f([s.observable],Qs.prototype,"defaultSlottedContent",void 0);class Zs{}f([(0,s.attr)({attribute:"aria-expanded"})],Zs.prototype,"ariaExpanded",void 0);f([(0,s.attr)({attribute:"aria-pressed"})],Zs.prototype,"ariaPressed",void 0);Pe(Zs,je);Pe(Qs,o,Zs);class Js{constructor(e){this.dayFormat="numeric";this.weekdayFormat="long";this.monthFormat="long";this.yearFormat="numeric";this.date=new Date;if(e){for(const t in e){const i=e[t];if(t==="date"){this.date=this.getDateObject(i)}else{this[t]=i}}}}getDateObject(e){if(typeof e==="string"){const t=e.split(/[/-]/);if(t.length<3){return new Date}return new Date(parseInt(t[2],10),parseInt(t[0],10)-1,parseInt(t[1],10))}else if("day"in e&&"month"in e&&"year"in e){const{day:t,month:i,year:s}=e;return new Date(s,i-1,t)}return e}getDate(e=this.date,t={weekday:this.weekdayFormat,month:this.monthFormat,day:this.dayFormat,year:this.yearFormat},i=this.locale){const s=this.getDateObject(e);if(!s.getTime()){return""}const o=Object.assign({timeZone:Intl.DateTimeFormat().resolvedOptions().timeZone},t);return new Intl.DateTimeFormat(i,o).format(s)}getDay(e=this.date.getDate(),t=this.dayFormat,i=this.locale){return this.getDate({month:1,day:e,year:2020},{day:t},i)}getMonth(e=this.date.getMonth()+1,t=this.monthFormat,i=this.locale){return this.getDate({month:e,day:2,year:2020},{month:t},i)}getYear(e=this.date.getFullYear(),t=this.yearFormat,i=this.locale){return this.getDate({month:2,day:2,year:e},{year:t},i)}getWeekday(e=0,t=this.weekdayFormat,i=this.locale){const s=`1-${e+1}-2017`;return this.getDate(s,{weekday:t},i)}getWeekdays(e=this.weekdayFormat,t=this.locale){return Array(7).fill(null).map(((i,s)=>this.getWeekday(s,e,t)))}}class eo extends Fe{constructor(){super(...arguments);this.dateFormatter=new Js;this.readonly=false;this.locale="en-US";this.month=(new Date).getMonth()+1;this.year=(new Date).getFullYear();this.dayFormat="numeric";this.weekdayFormat="short";this.monthFormat="long";this.yearFormat="numeric";this.minWeeks=0;this.disabledDates="";this.selectedDates="";this.oneDayInMs=864e5}localeChanged(){this.dateFormatter.locale=this.locale}dayFormatChanged(){this.dateFormatter.dayFormat=this.dayFormat}weekdayFormatChanged(){this.dateFormatter.weekdayFormat=this.weekdayFormat}monthFormatChanged(){this.dateFormatter.monthFormat=this.monthFormat}yearFormatChanged(){this.dateFormatter.yearFormat=this.yearFormat}getMonthInfo(e=this.month,t=this.year){const i=e=>new Date(e.getFullYear(),e.getMonth(),1).getDay();const s=e=>{const t=new Date(e.getFullYear(),e.getMonth()+1,1);return new Date(t.getTime()-this.oneDayInMs).getDate()};const o=new Date(t,e-1);const n=new Date(t,e);const r=new Date(t,e-2);return{length:s(o),month:e,start:i(o),year:t,previous:{length:s(r),month:r.getMonth()+1,start:i(r),year:r.getFullYear()},next:{length:s(n),month:n.getMonth()+1,start:i(n),year:n.getFullYear()}}}getDays(e=this.getMonthInfo(),t=this.minWeeks){t=t>10?10:t;const{start:i,length:s,previous:o,next:n}=e;const r=[];let a=1-i;while(a<s+1||r.length<t||r[r.length-1].length%7!==0){const{month:t,year:i}=a<1?o:a>s?n:e;const l=a<1?o.length+a:a>s?a-s:a;const d=`${t}-${l}-${i}`;const h=this.dateInString(d,this.disabledDates);const c=this.dateInString(d,this.selectedDates);const u={day:l,month:t,year:i,disabled:h,selected:c};const p=r[r.length-1];if(r.length===0||p.length%7===0){r.push([u])}else{p.push(u)}a++}return r}dateInString(e,t){const i=t.split(",").map((e=>e.trim()));e=typeof e==="string"?e:`${e.getMonth()+1}-${e.getDate()}-${e.getFullYear()}`;return i.some((t=>t===e))}getDayClassNames(e,t){const{day:i,month:s,year:o,disabled:n,selected:r}=e;const a=t===`${s}-${i}-${o}`;const l=this.month!==s;return["day",a&&"today",l&&"inactive",n&&"disabled",r&&"selected"].filter(Boolean).join(" ")}getWeekdayText(){const e=this.dateFormatter.getWeekdays().map((e=>({text:e})));if(this.weekdayFormat!=="long"){const t=this.dateFormatter.getWeekdays("long");e.forEach(((e,i)=>{e.abbr=t[i]}))}return e}handleDateSelect(e,t){e.preventDefault;this.$emit("dateselected",t)}handleKeydown(e,t){if(e.key===ze.Mm){this.handleDateSelect(e,t)}return true}}f([(0,s.attr)({mode:"boolean"})],eo.prototype,"readonly",void 0);f([s.attr],eo.prototype,"locale",void 0);f([(0,s.attr)({converter:s.nullableNumberConverter})],eo.prototype,"month",void 0);f([(0,s.attr)({converter:s.nullableNumberConverter})],eo.prototype,"year",void 0);f([(0,s.attr)({attribute:"day-format",mode:"fromView"})],eo.prototype,"dayFormat",void 0);f([(0,s.attr)({attribute:"weekday-format",mode:"fromView"})],eo.prototype,"weekdayFormat",void 0);f([(0,s.attr)({attribute:"month-format",mode:"fromView"})],eo.prototype,"monthFormat",void 0);f([(0,s.attr)({attribute:"year-format",mode:"fromView"})],eo.prototype,"yearFormat",void 0);f([(0,s.attr)({attribute:"min-weeks",converter:s.nullableNumberConverter})],eo.prototype,"minWeeks",void 0);f([(0,s.attr)({attribute:"disabled-dates"})],eo.prototype,"disabledDates",void 0);f([(0,s.attr)({attribute:"selected-dates"})],eo.prototype,"selectedDates",void 0);const to={none:"none",default:"default",sticky:"sticky"};const io={default:"default",columnHeader:"columnheader",rowHeader:"rowheader"};const so={default:"default",header:"header",stickyHeader:"sticky-header"};const oo=(0,s.html)`
    <template>
        ${e=>e.rowData===null||e.columnDefinition===null||e.columnDefinition.columnDataKey===null?null:e.rowData[e.columnDefinition.columnDataKey]}
    </template>
`;const no=(0,s.html)`
    <template>
        ${e=>e.columnDefinition===null?null:e.columnDefinition.title===undefined?e.columnDefinition.columnDataKey:e.columnDefinition.title}
    </template>
`;class ro extends Fe{constructor(){super(...arguments);this.cellType=io.default;this.rowData=null;this.columnDefinition=null;this.isActiveCell=false;this.customCellView=null;this.updateCellStyle=()=>{this.style.gridColumn=this.gridColumn}}cellTypeChanged(){if(this.$fastController.isConnected){this.updateCellView()}}gridColumnChanged(){if(this.$fastController.isConnected){this.updateCellStyle()}}columnDefinitionChanged(e,t){if(this.$fastController.isConnected){this.updateCellView()}}connectedCallback(){var e;super.connectedCallback();this.addEventListener(zt,this.handleFocusin);this.addEventListener(Nt,this.handleFocusout);this.addEventListener(Yt,this.handleKeydown);this.style.gridColumn=`${((e=this.columnDefinition)===null||e===void 0?void 0:e.gridColumn)===undefined?0:this.columnDefinition.gridColumn}`;this.updateCellView();this.updateCellStyle()}disconnectedCallback(){super.disconnectedCallback();this.removeEventListener(zt,this.handleFocusin);this.removeEventListener(Nt,this.handleFocusout);this.removeEventListener(Yt,this.handleKeydown);this.disconnectCellView()}handleFocusin(e){if(this.isActiveCell){return}this.isActiveCell=true;switch(this.cellType){case io.columnHeader:if(this.columnDefinition!==null&&this.columnDefinition.headerCellInternalFocusQueue!==true&&typeof this.columnDefinition.headerCellFocusTargetCallback==="function"){const e=this.columnDefinition.headerCellFocusTargetCallback(this);if(e!==null){e.focus()}}break;default:if(this.columnDefinition!==null&&this.columnDefinition.cellInternalFocusQueue!==true&&typeof this.columnDefinition.cellFocusTargetCallback==="function"){const e=this.columnDefinition.cellFocusTargetCallback(this);if(e!==null){e.focus()}}break}this.$emit("cell-focused",this)}handleFocusout(e){if(this!==document.activeElement&&!this.contains(document.activeElement)){this.isActiveCell=false}}handleKeydown(e){if(e.defaultPrevented||this.columnDefinition===null||this.cellType===io.default&&this.columnDefinition.cellInternalFocusQueue!==true||this.cellType===io.columnHeader&&this.columnDefinition.headerCellInternalFocusQueue!==true){return}switch(e.key){case ze.Mm:case ze.Ac:if(this.contains(document.activeElement)&&document.activeElement!==this){return}switch(this.cellType){case io.columnHeader:if(this.columnDefinition.headerCellFocusTargetCallback!==undefined){const t=this.columnDefinition.headerCellFocusTargetCallback(this);if(t!==null){t.focus()}e.preventDefault()}break;default:if(this.columnDefinition.cellFocusTargetCallback!==undefined){const t=this.columnDefinition.cellFocusTargetCallback(this);if(t!==null){t.focus()}e.preventDefault()}break}break;case ze.F9:if(this.contains(document.activeElement)&&document.activeElement!==this){this.focus();e.preventDefault()}break}}updateCellView(){this.disconnectCellView();if(this.columnDefinition===null){return}switch(this.cellType){case io.columnHeader:if(this.columnDefinition.headerCellTemplate!==undefined){this.customCellView=this.columnDefinition.headerCellTemplate.render(this,this)}else{this.customCellView=no.render(this,this)}break;case undefined:case io.rowHeader:case io.default:if(this.columnDefinition.cellTemplate!==undefined){this.customCellView=this.columnDefinition.cellTemplate.render(this,this)}else{this.customCellView=oo.render(this,this)}break}}disconnectCellView(){if(this.customCellView!==null){this.customCellView.dispose();this.customCellView=null}}}f([(0,s.attr)({attribute:"cell-type"})],ro.prototype,"cellType",void 0);f([(0,s.attr)({attribute:"grid-column"})],ro.prototype,"gridColumn",void 0);f([s.observable],ro.prototype,"rowData",void 0);f([s.observable],ro.prototype,"columnDefinition",void 0);class ao extends Fe{constructor(){super(...arguments);this.rowType=so.default;this.rowData=null;this.columnDefinitions=null;this.isActiveRow=false;this.cellsRepeatBehavior=null;this.cellsPlaceholder=null;this.focusColumnIndex=0;this.refocusOnLoad=false;this.updateRowStyle=()=>{this.style.gridTemplateColumns=this.gridTemplateColumns}}gridTemplateColumnsChanged(){if(this.$fastController.isConnected){this.updateRowStyle()}}rowTypeChanged(){if(this.$fastController.isConnected){this.updateItemTemplate()}}rowDataChanged(){if(this.rowData!==null&&this.isActiveRow){this.refocusOnLoad=true;return}}cellItemTemplateChanged(){this.updateItemTemplate()}headerCellItemTemplateChanged(){this.updateItemTemplate()}connectedCallback(){super.connectedCallback();if(this.cellsRepeatBehavior===null){this.cellsPlaceholder=document.createComment("");this.appendChild(this.cellsPlaceholder);this.updateItemTemplate();this.cellsRepeatBehavior=new s.RepeatDirective((e=>e.columnDefinitions),(e=>e.activeCellItemTemplate),{positioning:true}).createBehavior(this.cellsPlaceholder);this.$fastController.addBehaviors([this.cellsRepeatBehavior])}this.addEventListener("cell-focused",this.handleCellFocus);this.addEventListener(Nt,this.handleFocusout);this.addEventListener(Yt,this.handleKeydown);this.updateRowStyle();if(this.refocusOnLoad){this.refocusOnLoad=false;if(this.cellElements.length>this.focusColumnIndex){this.cellElements[this.focusColumnIndex].focus()}}}disconnectedCallback(){super.disconnectedCallback();this.removeEventListener("cell-focused",this.handleCellFocus);this.removeEventListener(Nt,this.handleFocusout);this.removeEventListener(Yt,this.handleKeydown)}handleFocusout(e){if(!this.contains(e.target)){this.isActiveRow=false;this.focusColumnIndex=0}}handleCellFocus(e){this.isActiveRow=true;this.focusColumnIndex=this.cellElements.indexOf(e.target);this.$emit("row-focused",this)}handleKeydown(e){if(e.defaultPrevented){return}let t=0;switch(e.key){case ze.kT:t=Math.max(0,this.focusColumnIndex-1);this.cellElements[t].focus();e.preventDefault();break;case ze.bb:t=Math.min(this.cellElements.length-1,this.focusColumnIndex+1);this.cellElements[t].focus();e.preventDefault();break;case ze.Tg:if(!e.ctrlKey){this.cellElements[0].focus();e.preventDefault()}break;case ze.FM:if(!e.ctrlKey){this.cellElements[this.cellElements.length-1].focus();e.preventDefault()}break}}updateItemTemplate(){this.activeCellItemTemplate=this.rowType===so.default&&this.cellItemTemplate!==undefined?this.cellItemTemplate:this.rowType===so.default&&this.cellItemTemplate===undefined?this.defaultCellItemTemplate:this.headerCellItemTemplate!==undefined?this.headerCellItemTemplate:this.defaultHeaderCellItemTemplate}}f([(0,s.attr)({attribute:"grid-template-columns"})],ao.prototype,"gridTemplateColumns",void 0);f([(0,s.attr)({attribute:"row-type"})],ao.prototype,"rowType",void 0);f([s.observable],ao.prototype,"rowData",void 0);f([s.observable],ao.prototype,"columnDefinitions",void 0);f([s.observable],ao.prototype,"cellItemTemplate",void 0);f([s.observable],ao.prototype,"headerCellItemTemplate",void 0);f([s.observable],ao.prototype,"rowIndex",void 0);f([s.observable],ao.prototype,"isActiveRow",void 0);f([s.observable],ao.prototype,"activeCellItemTemplate",void 0);f([s.observable],ao.prototype,"defaultCellItemTemplate",void 0);f([s.observable],ao.prototype,"defaultHeaderCellItemTemplate",void 0);f([s.observable],ao.prototype,"cellElements",void 0);class lo extends Fe{constructor(){super();this.noTabbing=false;this.generateHeader=to.default;this.rowsData=[];this.columnDefinitions=null;this.focusRowIndex=0;this.focusColumnIndex=0;this.rowsPlaceholder=null;this.generatedHeader=null;this.isUpdatingFocus=false;this.pendingFocusUpdate=false;this.rowindexUpdateQueued=false;this.columnDefinitionsStale=true;this.generatedGridTemplateColumns="";this.focusOnCell=(e,t,i)=>{if(this.rowElements.length===0){this.focusRowIndex=0;this.focusColumnIndex=0;return}const s=Math.max(0,Math.min(this.rowElements.length-1,e));const o=this.rowElements[s];const n=o.querySelectorAll('[role="cell"], [role="gridcell"], [role="columnheader"], [role="rowheader"]');const r=Math.max(0,Math.min(n.length-1,t));const a=n[r];if(i&&this.scrollHeight!==this.clientHeight&&(s<this.focusRowIndex&&this.scrollTop>0||s>this.focusRowIndex&&this.scrollTop<this.scrollHeight-this.clientHeight)){a.scrollIntoView({block:"center",inline:"center"})}a.focus()};this.onChildListChange=(e,t)=>{if(e&&e.length){e.forEach((e=>{e.addedNodes.forEach((e=>{if(e.nodeType===1&&e.getAttribute("role")==="row"){e.columnDefinitions=this.columnDefinitions}}))}));this.queueRowIndexUpdate()}};this.queueRowIndexUpdate=()=>{if(!this.rowindexUpdateQueued){this.rowindexUpdateQueued=true;s.DOM.queueUpdate(this.updateRowIndexes)}};this.updateRowIndexes=()=>{let e=this.gridTemplateColumns;if(e===undefined){if(this.generatedGridTemplateColumns===""&&this.rowElements.length>0){const e=this.rowElements[0];this.generatedGridTemplateColumns=new Array(e.cellElements.length).fill("1fr").join(" ")}e=this.generatedGridTemplateColumns}this.rowElements.forEach(((t,i)=>{const s=t;s.rowIndex=i;s.gridTemplateColumns=e;if(this.columnDefinitionsStale){s.columnDefinitions=this.columnDefinitions}}));this.rowindexUpdateQueued=false;this.columnDefinitionsStale=false}}static generateTemplateColumns(e){let t="";e.forEach((e=>{t=`${t}${t===""?"":" "}${"1fr"}`}));return t}noTabbingChanged(){if(this.$fastController.isConnected){if(this.noTabbing){this.setAttribute("tabIndex","-1")}else{this.setAttribute("tabIndex",this.contains(document.activeElement)||this===document.activeElement?"-1":"0")}}}generateHeaderChanged(){if(this.$fastController.isConnected){this.toggleGeneratedHeader()}}gridTemplateColumnsChanged(){if(this.$fastController.isConnected){this.updateRowIndexes()}}rowsDataChanged(){if(this.columnDefinitions===null&&this.rowsData.length>0){this.columnDefinitions=lo.generateColumns(this.rowsData[0])}if(this.$fastController.isConnected){this.toggleGeneratedHeader()}}columnDefinitionsChanged(){if(this.columnDefinitions===null){this.generatedGridTemplateColumns="";return}this.generatedGridTemplateColumns=lo.generateTemplateColumns(this.columnDefinitions);if(this.$fastController.isConnected){this.columnDefinitionsStale=true;this.queueRowIndexUpdate()}}headerCellItemTemplateChanged(){if(this.$fastController.isConnected){if(this.generatedHeader!==null){this.generatedHeader.headerCellItemTemplate=this.headerCellItemTemplate}}}focusRowIndexChanged(){if(this.$fastController.isConnected){this.queueFocusUpdate()}}focusColumnIndexChanged(){if(this.$fastController.isConnected){this.queueFocusUpdate()}}connectedCallback(){super.connectedCallback();if(this.rowItemTemplate===undefined){this.rowItemTemplate=this.defaultRowItemTemplate}this.rowsPlaceholder=document.createComment("");this.appendChild(this.rowsPlaceholder);this.toggleGeneratedHeader();this.rowsRepeatBehavior=new s.RepeatDirective((e=>e.rowsData),(e=>e.rowItemTemplate),{positioning:true}).createBehavior(this.rowsPlaceholder);this.$fastController.addBehaviors([this.rowsRepeatBehavior]);this.addEventListener("row-focused",this.handleRowFocus);this.addEventListener(Vt,this.handleFocus);this.addEventListener(Yt,this.handleKeydown);this.addEventListener(Nt,this.handleFocusOut);this.observer=new MutationObserver(this.onChildListChange);this.observer.observe(this,{childList:true});if(this.noTabbing){this.setAttribute("tabindex","-1")}s.DOM.queueUpdate(this.queueRowIndexUpdate)}disconnectedCallback(){super.disconnectedCallback();this.removeEventListener("row-focused",this.handleRowFocus);this.removeEventListener(Vt,this.handleFocus);this.removeEventListener(Yt,this.handleKeydown);this.removeEventListener(Nt,this.handleFocusOut);this.observer.disconnect();this.rowsPlaceholder=null;this.generatedHeader=null}handleRowFocus(e){this.isUpdatingFocus=true;const t=e.target;this.focusRowIndex=this.rowElements.indexOf(t);this.focusColumnIndex=t.focusColumnIndex;this.setAttribute("tabIndex","-1");this.isUpdatingFocus=false}handleFocus(e){this.focusOnCell(this.focusRowIndex,this.focusColumnIndex,true)}handleFocusOut(e){if(e.relatedTarget===null||!this.contains(e.relatedTarget)){this.setAttribute("tabIndex",this.noTabbing?"-1":"0")}}handleKeydown(e){if(e.defaultPrevented){return}let t;const i=this.rowElements.length-1;const s=this.offsetHeight+this.scrollTop;const o=this.rowElements[i];switch(e.key){case ze.I5:e.preventDefault();this.focusOnCell(this.focusRowIndex-1,this.focusColumnIndex,true);break;case ze.HX:e.preventDefault();this.focusOnCell(this.focusRowIndex+1,this.focusColumnIndex,true);break;case ze.oK:e.preventDefault();if(this.rowElements.length===0){this.focusOnCell(0,0,false);break}if(this.focusRowIndex===0){this.focusOnCell(0,this.focusColumnIndex,false);return}t=this.focusRowIndex-1;for(t;t>=0;t--){const e=this.rowElements[t];if(e.offsetTop<this.scrollTop){this.scrollTop=e.offsetTop+e.clientHeight-this.clientHeight;break}}this.focusOnCell(t,this.focusColumnIndex,false);break;case ze.f_:e.preventDefault();if(this.rowElements.length===0){this.focusOnCell(0,0,false);break}if(this.focusRowIndex>=i||o.offsetTop+o.offsetHeight<=s){this.focusOnCell(i,this.focusColumnIndex,false);return}t=this.focusRowIndex+1;for(t;t<=i;t++){const e=this.rowElements[t];if(e.offsetTop+e.offsetHeight>s){let t=0;if(this.generateHeader===to.sticky&&this.generatedHeader!==null){t=this.generatedHeader.clientHeight}this.scrollTop=e.offsetTop-t;break}}this.focusOnCell(t,this.focusColumnIndex,false);break;case ze.Tg:if(e.ctrlKey){e.preventDefault();this.focusOnCell(0,0,true)}break;case ze.FM:if(e.ctrlKey&&this.columnDefinitions!==null){e.preventDefault();this.focusOnCell(this.rowElements.length-1,this.columnDefinitions.length-1,true)}break}}queueFocusUpdate(){if(this.isUpdatingFocus&&(this.contains(document.activeElement)||this===document.activeElement)){return}if(this.pendingFocusUpdate===false){this.pendingFocusUpdate=true;s.DOM.queueUpdate((()=>this.updateFocus()))}}updateFocus(){this.pendingFocusUpdate=false;this.focusOnCell(this.focusRowIndex,this.focusColumnIndex,true)}toggleGeneratedHeader(){if(this.generatedHeader!==null){this.removeChild(this.generatedHeader);this.generatedHeader=null}if(this.generateHeader!==to.none&&this.rowsData.length>0){const e=document.createElement(this.rowElementTag);this.generatedHeader=e;this.generatedHeader.columnDefinitions=this.columnDefinitions;this.generatedHeader.gridTemplateColumns=this.gridTemplateColumns;this.generatedHeader.rowType=this.generateHeader===to.sticky?so.stickyHeader:so.header;if(this.firstChild!==null||this.rowsPlaceholder!==null){this.insertBefore(e,this.firstChild!==null?this.firstChild:this.rowsPlaceholder)}return}}}lo.generateColumns=e=>Object.getOwnPropertyNames(e).map(((e,t)=>({columnDataKey:e,gridColumn:`${t}`})));f([(0,s.attr)({attribute:"no-tabbing",mode:"boolean"})],lo.prototype,"noTabbing",void 0);f([(0,s.attr)({attribute:"generate-header"})],lo.prototype,"generateHeader",void 0);f([(0,s.attr)({attribute:"grid-template-columns"})],lo.prototype,"gridTemplateColumns",void 0);f([s.observable],lo.prototype,"rowsData",void 0);f([s.observable],lo.prototype,"columnDefinitions",void 0);f([s.observable],lo.prototype,"rowItemTemplate",void 0);f([s.observable],lo.prototype,"cellItemTemplate",void 0);f([s.observable],lo.prototype,"headerCellItemTemplate",void 0);f([s.observable],lo.prototype,"focusRowIndex",void 0);f([s.observable],lo.prototype,"focusColumnIndex",void 0);f([s.observable],lo.prototype,"defaultRowItemTemplate",void 0);f([s.observable],lo.prototype,"rowElementTag",void 0);f([s.observable],lo.prototype,"rowElements",void 0);const ho=(0,s.html)`
    <div
        class="title"
        part="title"
        aria-label="${e=>e.dateFormatter.getDate(`${e.month}-2-${e.year}`,{month:"long",year:"numeric"})}"
    >
        <span part="month">
            ${e=>e.dateFormatter.getMonth(e.month)}
        </span>
        <span part="year">${e=>e.dateFormatter.getYear(e.year)}</span>
    </div>
`;const co=e=>{const t=e.tagFor(ro);return(0,s.html)`
        <${t}
            class="week-day"
            part="week-day"
            tabindex="-1"
            grid-column="${(e,t)=>t.index+1}"
            abbr="${e=>e.abbr}"
        >
            ${e=>e.text}
        </${t}>
    `};const uo=(e,t)=>{const i=e.tagFor(ro);return(0,s.html)`
        <${i}
            class="${(e,i)=>i.parentContext.parent.getDayClassNames(e,t)}"
            part="day"
            tabindex="-1"
            role="gridcell"
            grid-column="${(e,t)=>t.index+1}"
            @click="${(e,t)=>t.parentContext.parent.handleDateSelect(t.event,e)}"
            @keydown="${(e,t)=>t.parentContext.parent.handleKeydown(t.event,e)}"
            aria-label="${(e,t)=>t.parentContext.parent.dateFormatter.getDate(`${e.month}-${e.day}-${e.year}`,{month:"long",day:"numeric"})}"
        >
            <div
                class="date"
                part="${e=>t===`${e.month}-${e.day}-${e.year}`?"today":"date"}"
            >
                ${(e,t)=>t.parentContext.parent.dateFormatter.getDay(e.day)}
            </div>
            <slot name="${e=>e.month}-${e=>e.day}-${e=>e.year}"></slot>
        </${i}>
    `};const po=(e,t)=>{const i=e.tagFor(ao);return(0,s.html)`
        <${i}
            class="week"
            part="week"
            role="row"
            role-type="default"
            grid-template-columns="1fr 1fr 1fr 1fr 1fr 1fr 1fr"
        >
        ${(0,s.repeat)((e=>e),uo(e,t),{positioning:true})}
        </${i}>
    `};const fo=(e,t)=>{const i=e.tagFor(lo);const o=e.tagFor(ao);return(0,s.html)`
    <${i} class="days interact" part="days" generate-header="none">
        <${o}
            class="week-days"
            part="week-days"
            role="row"
            row-type="header"
            grid-template-columns="1fr 1fr 1fr 1fr 1fr 1fr 1fr"
        >
            ${(0,s.repeat)((e=>e.getWeekdayText()),co(e),{positioning:true})}
        </${o}>
        ${(0,s.repeat)((e=>e.getDays()),po(e,t))}
    </${i}>
`};const mo=e=>(0,s.html)`
        <div class="days" part="days">
            <div class="week-days" part="week-days">
                ${(0,s.repeat)((e=>e.getWeekdayText()),(0,s.html)`
                        <div class="week-day" part="week-day" abbr="${e=>e.abbr}">
                            ${e=>e.text}
                        </div>
                    `)}
            </div>
            ${(0,s.repeat)((e=>e.getDays()),(0,s.html)`
                    <div class="week">
                        ${(0,s.repeat)((e=>e),(0,s.html)`
                                <div
                                    class="${(t,i)=>i.parentContext.parent.getDayClassNames(t,e)}"
                                    part="day"
                                    aria-label="${(e,t)=>t.parentContext.parent.dateFormatter.getDate(`${e.month}-${e.day}-${e.year}`,{month:"long",day:"numeric"})}"
                                >
                                    <div
                                        class="date"
                                        part="${t=>e===`${t.month}-${t.day}-${t.year}`?"today":"date"}"
                                    >
                                        ${(e,t)=>t.parentContext.parent.dateFormatter.getDay(e.day)}
                                    </div>
                                    <slot
                                        name="${e=>e.month}-${e=>e.day}-${e=>e.year}"
                                    ></slot>
                                </div>
                            `)}
                    </div>
                `)}
        </div>
    `;const vo=(e,t)=>{var i;const o=new Date;const n=`${o.getMonth()+1}-${o.getDate()}-${o.getFullYear()}`;return(0,s.html)`
        <template>
            ${l}
            ${t.title instanceof Function?t.title(e,t):(i=t.title)!==null&&i!==void 0?i:""}
            <slot></slot>
            ${(0,s.when)((e=>e.readonly),mo(n),fo(e,n))}
            ${a}
        </template>
    `};const bo=(e,t)=>(0,s.html)`
    <slot></slot>
`;class go extends Fe{}const yo=(e,t)=>(0,s.html)`
    <template
        role="checkbox"
        aria-checked="${e=>e.checked}"
        aria-required="${e=>e.required}"
        aria-disabled="${e=>e.disabled}"
        aria-readonly="${e=>e.readOnly}"
        tabindex="${e=>e.disabled?null:0}"
        @keypress="${(e,t)=>e.keypressHandler(t.event)}"
        @click="${(e,t)=>e.clickHandler(t.event)}"
        class="${e=>e.readOnly?"readonly":""} ${e=>e.checked?"checked":""} ${e=>e.indeterminate?"indeterminate":""}"
    >
        <div part="control" class="control">
            <slot name="checked-indicator">
                ${t.checkedIndicator||""}
            </slot>
            <slot name="indeterminate-indicator">
                ${t.indeterminateIndicator||""}
            </slot>
        </div>
        <label
            part="label"
            class="${e=>e.defaultSlottedNodes&&e.defaultSlottedNodes.length?"label":"label label__hidden"}"
        >
            <slot ${(0,s.slotted)("defaultSlottedNodes")}></slot>
        </label>
    </template>
`;class Co extends Fe{}class xo extends(Gs(Co)){constructor(){super(...arguments);this.proxy=document.createElement("input")}}class wo extends xo{constructor(){super();this.initialValue="on";this.indeterminate=false;this.keypressHandler=e=>{if(this.readOnly){return}switch(e.key){case ze.gG:if(this.indeterminate){this.indeterminate=false}this.checked=!this.checked;break}};this.clickHandler=e=>{if(!this.disabled&&!this.readOnly){if(this.indeterminate){this.indeterminate=false}this.checked=!this.checked}};this.proxy.setAttribute("type","checkbox")}readOnlyChanged(){if(this.proxy instanceof HTMLInputElement){this.proxy.readOnly=this.readOnly}}}f([(0,s.attr)({attribute:"readonly",mode:"boolean"})],wo.prototype,"readOnly",void 0);f([s.observable],wo.prototype,"defaultSlottedNodes",void 0);f([s.observable],wo.prototype,"indeterminate",void 0);let $o=0;function Io(e=""){return`${e}${$o++}`}function ko(e,...t){return e.replace(/{(\d+)}/g,(function(e,i){if(i>=t.length){return e}const s=t[i];if(typeof s!=="number"&&!s){return""}return s}))}function Oo(e,t,i=0){if(!e||!t){return false}return e.substr(i,t.length)===t}function To(e){return!e||!e.trim()}function Eo(e){let t=`${e}`.replace(new RegExp(/[-_]+/,"g")," ").replace(new RegExp(/[^\w\s]/,"g"),"").replace(/^\s+|\s+$|\s+(?=\s)/g,"").replace(new RegExp(/\s+(.)(\w*)/,"g"),((e,t,i)=>`${t.toUpperCase()+i.toLowerCase()}`)).replace(new RegExp(/\w/),(e=>e.toUpperCase()));let i=0;for(let s=0;s<t.length;s++){const e=t.charAt(s);if(e==e.toLowerCase()){i=s;break}}if(i>1){t=`${t.charAt(0).toUpperCase()}${t.slice(1,i-1).toLowerCase()}`+t.slice(i-1)}return t}function Ro(e){const t=`${e.charAt(0).toLowerCase()}${e.slice(1)}`;return t.replace(/([A-Z]|[0-9])/g,(function(e,t){return`-${t.toLowerCase()}`}))}function Do(e,t){let i=e.length;while(i--){if(t(e[i],i,e)){return i}}return-1}function So(){return!!(typeof window!=="undefined"&&window.document&&window.document.createElement)}function Ao(...e){return e.every((e=>e instanceof HTMLElement))}function Fo(e,t){if(!e||!t||!Ao(e)){return}const i=Array.from(e.querySelectorAll(t));return i.filter((e=>e.offsetParent!==null))}function Lo(e){return e===null?null:e.which||e.keyCode||e.charCode}function Mo(){const e=document.querySelector('meta[property="csp-nonce"]');if(e){return e.getAttribute("content")}else{return null}}let Po;function Ho(){if(typeof Po==="boolean"){return Po}if(!So()){Po=false;return Po}const e=document.createElement("style");const t=Mo();if(t!==null){e.setAttribute("nonce",t)}document.head.appendChild(e);try{e.sheet.insertRule("foo:focus-visible {color:inherit}",0);Po=true}catch(i){Po=false}finally{document.head.removeChild(e)}return Po}let Vo;function zo(){if(typeof Vo==="boolean"){return Vo}try{Vo=CSS.supports("display","grid")}catch(e){Vo=false}return Vo}function No(){return canUseDOM()&&(window.matchMedia("(forced-colors: none)").matches||window.matchMedia("(forced-colors: active)").matches)}function Bo(){Vo=undefined;Po=undefined}const qo=null&&No;function Uo(e){return Ao(e)&&(e.getAttribute("role")==="option"||e instanceof HTMLOptionElement)}class jo extends Fe{constructor(e,t,i,s){super();this.defaultSelected=false;this.dirtySelected=false;this.selected=this.defaultSelected;this.dirtyValue=false;if(e){this.textContent=e}if(t){this.initialValue=t}if(i){this.defaultSelected=i}if(s){this.selected=s}this.proxy=new Option(`${this.textContent}`,this.initialValue,this.defaultSelected,this.selected);this.proxy.disabled=this.disabled}checkedChanged(e,t){if(typeof t==="boolean"){this.ariaChecked=t?"true":"false";return}this.ariaChecked=null}contentChanged(e,t){if(this.proxy instanceof HTMLOptionElement){this.proxy.textContent=this.textContent}this.$emit("contentchange",null,{bubbles:true})}defaultSelectedChanged(){if(!this.dirtySelected){this.selected=this.defaultSelected;if(this.proxy instanceof HTMLOptionElement){this.proxy.selected=this.defaultSelected}}}disabledChanged(e,t){this.ariaDisabled=this.disabled?"true":"false";if(this.proxy instanceof HTMLOptionElement){this.proxy.disabled=this.disabled}}selectedAttributeChanged(){this.defaultSelected=this.selectedAttribute;if(this.proxy instanceof HTMLOptionElement){this.proxy.defaultSelected=this.defaultSelected}}selectedChanged(){this.ariaSelected=this.selected?"true":"false";if(!this.dirtySelected){this.dirtySelected=true}if(this.proxy instanceof HTMLOptionElement){this.proxy.selected=this.selected}}initialValueChanged(e,t){if(!this.dirtyValue){this.value=this.initialValue;this.dirtyValue=false}}get label(){var e;return(e=this.value)!==null&&e!==void 0?e:this.text}get text(){var e,t;return(t=(e=this.textContent)===null||e===void 0?void 0:e.replace(/\s+/g," ").trim())!==null&&t!==void 0?t:""}set value(e){const t=`${e!==null&&e!==void 0?e:""}`;this._value=t;this.dirtyValue=true;if(this.proxy instanceof HTMLOptionElement){this.proxy.value=t}s.Observable.notify(this,"value")}get value(){var e;s.Observable.track(this,"value");return(e=this._value)!==null&&e!==void 0?e:this.text}get form(){return this.proxy?this.proxy.form:null}}f([s.observable],jo.prototype,"checked",void 0);f([s.observable],jo.prototype,"content",void 0);f([s.observable],jo.prototype,"defaultSelected",void 0);f([(0,s.attr)({mode:"boolean"})],jo.prototype,"disabled",void 0);f([(0,s.attr)({attribute:"selected",mode:"boolean"})],jo.prototype,"selectedAttribute",void 0);f([s.observable],jo.prototype,"selected",void 0);f([(0,s.attr)({attribute:"value",mode:"fromView"})],jo.prototype,"initialValue",void 0);class _o{}f([s.observable],_o.prototype,"ariaChecked",void 0);f([s.observable],_o.prototype,"ariaPosInSet",void 0);f([s.observable],_o.prototype,"ariaSelected",void 0);f([s.observable],_o.prototype,"ariaSetSize",void 0);Pe(_o,je);Pe(jo,o,_o);class Ko extends Fe{constructor(){super(...arguments);this._options=[];this.selectedIndex=-1;this.selectedOptions=[];this.shouldSkipFocus=false;this.typeaheadBuffer="";this.typeaheadExpired=true;this.typeaheadTimeout=-1}get firstSelectedOption(){var e;return(e=this.selectedOptions[0])!==null&&e!==void 0?e:null}get hasSelectableOptions(){return this.options.length>0&&!this.options.every((e=>e.disabled))}get length(){var e,t;return(t=(e=this.options)===null||e===void 0?void 0:e.length)!==null&&t!==void 0?t:0}get options(){s.Observable.track(this,"options");return this._options}set options(e){this._options=e;s.Observable.notify(this,"options")}get typeAheadExpired(){return this.typeaheadExpired}set typeAheadExpired(e){this.typeaheadExpired=e}clickHandler(e){const t=e.target.closest(`option,[role=option]`);if(t&&!t.disabled){this.selectedIndex=this.options.indexOf(t);return true}}focusAndScrollOptionIntoView(e=this.firstSelectedOption){if(this.contains(document.activeElement)&&e!==null){e.focus();requestAnimationFrame((()=>{e.scrollIntoView({block:"nearest"})}))}}focusinHandler(e){if(!this.shouldSkipFocus&&e.target===e.currentTarget){this.setSelectedOptions();this.focusAndScrollOptionIntoView()}this.shouldSkipFocus=false}getTypeaheadMatches(){const e=this.typeaheadBuffer.replace(/[.*+\-?^${}()|[\]\\]/g,"\\$&");const t=new RegExp(`^${e}`,"gi");return this.options.filter((e=>e.text.trim().match(t)))}getSelectableIndex(e=this.selectedIndex,t){const i=e>t?-1:e<t?1:0;const s=e+i;let o=null;switch(i){case-1:{o=this.options.reduceRight(((e,t,i)=>!e&&!t.disabled&&i<s?t:e),o);break}case 1:{o=this.options.reduce(((e,t,i)=>!e&&!t.disabled&&i>s?t:e),o);break}}return this.options.indexOf(o)}handleChange(e,t){switch(t){case"selected":{if(Ko.slottedOptionFilter(e)){this.selectedIndex=this.options.indexOf(e)}this.setSelectedOptions();break}}}handleTypeAhead(e){if(this.typeaheadTimeout){window.clearTimeout(this.typeaheadTimeout)}this.typeaheadTimeout=window.setTimeout((()=>this.typeaheadExpired=true),Ko.TYPE_AHEAD_TIMEOUT_MS);if(e.length>1){return}this.typeaheadBuffer=`${this.typeaheadExpired?"":this.typeaheadBuffer}${e}`}keydownHandler(e){if(this.disabled){return true}this.shouldSkipFocus=false;const t=e.key;switch(t){case ze.Tg:{if(!e.shiftKey){e.preventDefault();this.selectFirstOption()}break}case ze.HX:{if(!e.shiftKey){e.preventDefault();this.selectNextOption()}break}case ze.I5:{if(!e.shiftKey){e.preventDefault();this.selectPreviousOption()}break}case ze.FM:{e.preventDefault();this.selectLastOption();break}case ze.J9:{this.focusAndScrollOptionIntoView();return true}case ze.Mm:case ze.F9:{return true}case ze.gG:{if(this.typeaheadExpired){return true}}default:{if(t.length===1){this.handleTypeAhead(`${t}`)}return true}}}mousedownHandler(e){this.shouldSkipFocus=!this.contains(document.activeElement);return true}multipleChanged(e,t){this.ariaMultiSelectable=t?"true":null}selectedIndexChanged(e,t){var i;if(!this.hasSelectableOptions){this.selectedIndex=-1;return}if(((i=this.options[this.selectedIndex])===null||i===void 0?void 0:i.disabled)&&typeof e==="number"){const i=this.getSelectableIndex(e,t);const s=i>-1?i:e;this.selectedIndex=s;if(t===s){this.selectedIndexChanged(t,s)}return}this.setSelectedOptions()}selectedOptionsChanged(e,t){var i;const o=t.filter(Ko.slottedOptionFilter);(i=this.options)===null||i===void 0?void 0:i.forEach((e=>{const t=s.Observable.getNotifier(e);t.unsubscribe(this,"selected");e.selected=o.includes(e);t.subscribe(this,"selected")}))}selectFirstOption(){var e,t;if(!this.disabled){this.selectedIndex=(t=(e=this.options)===null||e===void 0?void 0:e.findIndex((e=>!e.disabled)))!==null&&t!==void 0?t:-1}}selectLastOption(){if(!this.disabled){this.selectedIndex=Do(this.options,(e=>!e.disabled))}}selectNextOption(){if(!this.disabled&&this.selectedIndex<this.options.length-1){this.selectedIndex+=1}}selectPreviousOption(){if(!this.disabled&&this.selectedIndex>0){this.selectedIndex=this.selectedIndex-1}}setDefaultSelectedOption(){var e,t;this.selectedIndex=(t=(e=this.options)===null||e===void 0?void 0:e.findIndex((e=>e.defaultSelected)))!==null&&t!==void 0?t:-1}setSelectedOptions(){var e,t,i;if((e=this.options)===null||e===void 0?void 0:e.length){this.selectedOptions=[this.options[this.selectedIndex]];this.ariaActiveDescendant=(i=(t=this.firstSelectedOption)===null||t===void 0?void 0:t.id)!==null&&i!==void 0?i:"";this.focusAndScrollOptionIntoView()}}slottedOptionsChanged(e,t){this.options=t.reduce(((e,t)=>{if(Uo(t)){e.push(t)}return e}),[]);const i=`${this.options.length}`;this.options.forEach(((e,t)=>{if(!e.id){e.id=Io("option-")}e.ariaPosInSet=`${t+1}`;e.ariaSetSize=i}));if(this.$fastController.isConnected){this.setSelectedOptions();this.setDefaultSelectedOption()}}typeaheadBufferChanged(e,t){if(this.$fastController.isConnected){const e=this.getTypeaheadMatches();if(e.length){const t=this.options.indexOf(e[0]);if(t>-1){this.selectedIndex=t}}this.typeaheadExpired=false}}}Ko.slottedOptionFilter=e=>Uo(e)&&!e.hidden;Ko.TYPE_AHEAD_TIMEOUT_MS=1e3;f([(0,s.attr)({mode:"boolean"})],Ko.prototype,"disabled",void 0);f([s.observable],Ko.prototype,"selectedIndex",void 0);f([s.observable],Ko.prototype,"selectedOptions",void 0);f([s.observable],Ko.prototype,"slottedOptions",void 0);f([s.observable],Ko.prototype,"typeaheadBuffer",void 0);class Wo{}f([s.observable],Wo.prototype,"ariaActiveDescendant",void 0);f([s.observable],Wo.prototype,"ariaDisabled",void 0);f([s.observable],Wo.prototype,"ariaExpanded",void 0);f([s.observable],Wo.prototype,"ariaMultiSelectable",void 0);Pe(Wo,je);Pe(Ko,Wo);const Go={above:"above",below:"below"};class Xo extends Ko{}class Yo extends(Ws(Xo)){constructor(){super(...arguments);this.proxy=document.createElement("input")}}const Qo={inline:"inline",list:"list",both:"both",none:"none"};class Zo extends Yo{constructor(){super(...arguments);this._value="";this.filteredOptions=[];this.filter="";this.forcedPosition=false;this.listboxId=Io("listbox-");this.maxHeight=0;this.open=false}formResetCallback(){super.formResetCallback();this.setDefaultSelectedOption();this.updateValue()}validate(){super.validate(this.control)}get isAutocompleteInline(){return this.autocomplete===Qo.inline||this.isAutocompleteBoth}get isAutocompleteList(){return this.autocomplete===Qo.list||this.isAutocompleteBoth}get isAutocompleteBoth(){return this.autocomplete===Qo.both}openChanged(){if(this.open){this.ariaControls=this.listboxId;this.ariaExpanded="true";this.setPositioning();this.focusAndScrollOptionIntoView();s.DOM.queueUpdate((()=>this.focus()));return}this.ariaControls="";this.ariaExpanded="false"}get options(){s.Observable.track(this,"options");return this.filteredOptions.length?this.filteredOptions:this._options}set options(e){this._options=e;s.Observable.notify(this,"options")}placeholderChanged(){if(this.proxy instanceof HTMLInputElement){this.proxy.placeholder=this.placeholder}}positionChanged(e,t){this.positionAttribute=t;this.setPositioning()}get value(){s.Observable.track(this,"value");return this._value}set value(e){var t,i,o;const n=`${this._value}`;if(this.$fastController.isConnected&&this.options){const s=this.options.findIndex((t=>t.text.toLowerCase()===e.toLowerCase()));const n=(t=this.options[this.selectedIndex])===null||t===void 0?void 0:t.text;const r=(i=this.options[s])===null||i===void 0?void 0:i.text;this.selectedIndex=n!==r?s:this.selectedIndex;e=((o=this.firstSelectedOption)===null||o===void 0?void 0:o.text)||e}if(n!==e){this._value=e;super.valueChanged(n,e);s.Observable.notify(this,"value")}}clickHandler(e){if(this.disabled){return}if(this.open){const t=e.target.closest(`option,[role=option]`);if(!t||t.disabled){return}this.selectedOptions=[t];this.control.value=t.text;this.clearSelectionRange();this.updateValue(true)}this.open=!this.open;if(this.open){this.control.focus()}return true}connectedCallback(){super.connectedCallback();this.forcedPosition=!!this.positionAttribute;if(this.value){this.initialValue=this.value}}disabledChanged(e,t){if(super.disabledChanged){super.disabledChanged(e,t)}this.ariaDisabled=this.disabled?"true":"false"}filterOptions(){if(!this.autocomplete||this.autocomplete===Qo.none){this.filter=""}const e=this.filter.toLowerCase();this.filteredOptions=this._options.filter((e=>e.text.toLowerCase().startsWith(this.filter.toLowerCase())));if(this.isAutocompleteList){if(!this.filteredOptions.length&&!e){this.filteredOptions=this._options}this._options.forEach((e=>{e.hidden=!this.filteredOptions.includes(e)}))}}focusAndScrollOptionIntoView(){if(this.contains(document.activeElement)){this.control.focus();if(this.firstSelectedOption){requestAnimationFrame((()=>{var e;(e=this.firstSelectedOption)===null||e===void 0?void 0:e.scrollIntoView({block:"nearest"})}))}}}focusoutHandler(e){this.syncValue();if(!this.open){return true}const t=e.relatedTarget;if(this.isSameNode(t)){this.focus();return}if(!this.options||!this.options.includes(t)){this.open=false}}inputHandler(e){this.filter=this.control.value;this.filterOptions();if(!this.isAutocompleteInline){this.selectedIndex=this.options.map((e=>e.text)).indexOf(this.control.value)}if(e.inputType.includes("deleteContent")||!this.filter.length){return true}if(this.isAutocompleteList&&!this.open){this.open=true}if(this.isAutocompleteInline){if(this.filteredOptions.length){this.selectedOptions=[this.filteredOptions[0]];this.selectedIndex=this.options.indexOf(this.firstSelectedOption);this.setInlineSelection()}else{this.selectedIndex=-1}}return}keydownHandler(e){const t=e.key;if(e.ctrlKey||e.shiftKey){return true}switch(t){case"Enter":{this.syncValue();if(this.isAutocompleteInline){this.filter=this.value}this.open=false;this.clearSelectionRange();break}case"Escape":{if(!this.isAutocompleteInline){this.selectedIndex=-1}if(this.open){this.open=false;break}this.value="";this.control.value="";this.filter="";this.filterOptions();break}case"Tab":{this.setInputToSelection();if(!this.open){return true}e.preventDefault();this.open=false;break}case"ArrowUp":case"ArrowDown":{this.filterOptions();if(!this.open){this.open=true;break}if(this.filteredOptions.length>0){super.keydownHandler(e)}if(this.isAutocompleteInline){this.setInlineSelection()}break}default:{return true}}}keyupHandler(e){const t=e.key;switch(t){case"ArrowLeft":case"ArrowRight":case"Backspace":case"Delete":case"Home":case"End":{this.filter=this.control.value;this.selectedIndex=-1;this.filterOptions();break}}}selectedIndexChanged(e,t){if(this.$fastController.isConnected){t=(0,Ne.AB)(-1,this.options.length-1,t);if(t!==this.selectedIndex){this.selectedIndex=t;return}super.selectedIndexChanged(e,t)}}selectPreviousOption(){if(!this.disabled&&this.selectedIndex>=0){this.selectedIndex=this.selectedIndex-1}}setDefaultSelectedOption(){if(this.$fastController.isConnected&&this.options){const e=this.options.findIndex((e=>e.getAttribute("selected")!==null||e.selected));this.selectedIndex=e;if(!this.dirtyValue&&this.firstSelectedOption){this.value=this.firstSelectedOption.text}this.setSelectedOptions()}}setInputToSelection(){if(this.firstSelectedOption){this.control.value=this.firstSelectedOption.text;this.control.focus()}}setInlineSelection(){if(this.firstSelectedOption){this.setInputToSelection();this.control.setSelectionRange(this.filter.length,this.control.value.length,"backward")}}syncValue(){var e;const t=this.selectedIndex>-1?(e=this.firstSelectedOption)===null||e===void 0?void 0:e.text:this.control.value;this.updateValue(this.value!==t)}setPositioning(){const e=this.getBoundingClientRect();const t=window.innerHeight;const i=t-e.bottom;this.position=this.forcedPosition?this.positionAttribute:e.top>i?Go.above:Go.below;this.positionAttribute=this.forcedPosition?this.positionAttribute:this.position;this.maxHeight=this.position===Go.above?~~e.top:~~i}selectedOptionsChanged(e,t){if(this.$fastController.isConnected){this._options.forEach((e=>{e.selected=t.includes(e)}))}}slottedOptionsChanged(e,t){super.slottedOptionsChanged(e,t);this.updateValue()}updateValue(e){var t;if(this.$fastController.isConnected){this.value=((t=this.firstSelectedOption)===null||t===void 0?void 0:t.text)||this.control.value;this.control.value=this.value}if(e){this.$emit("change")}}clearSelectionRange(){const e=this.control.value.length;this.control.setSelectionRange(e,e)}}f([(0,s.attr)({attribute:"autocomplete",mode:"fromView"})],Zo.prototype,"autocomplete",void 0);f([s.observable],Zo.prototype,"maxHeight",void 0);f([(0,s.attr)({attribute:"open",mode:"boolean"})],Zo.prototype,"open",void 0);f([s.attr],Zo.prototype,"placeholder",void 0);f([(0,s.attr)({attribute:"position"})],Zo.prototype,"positionAttribute",void 0);f([s.observable],Zo.prototype,"position",void 0);class Jo{}f([s.observable],Jo.prototype,"ariaAutoComplete",void 0);f([s.observable],Jo.prototype,"ariaControls",void 0);Pe(Jo,Wo);Pe(Zo,o,Jo);const en=(e,t)=>(0,s.html)`
    <template
        aria-disabled="${e=>e.ariaDisabled}"
        autocomplete="${e=>e.autocomplete}"
        class="${e=>e.open?"open":""} ${e=>e.disabled?"disabled":""} ${e=>e.position}"
        ?open="${e=>e.open}"
        tabindex="${e=>!e.disabled?"0":null}"
        @click="${(e,t)=>e.clickHandler(t.event)}"
        @focusout="${(e,t)=>e.focusoutHandler(t.event)}"
        @keydown="${(e,t)=>e.keydownHandler(t.event)}"
    >
        <div class="control" part="control">
            ${r(e,t)}
            <slot name="control">
                <input
                    aria-activedescendant="${e=>e.open?e.ariaActiveDescendant:null}"
                    aria-autocomplete="${e=>e.ariaAutoComplete}"
                    aria-controls="${e=>e.ariaControls}"
                    aria-disabled="${e=>e.ariaDisabled}"
                    aria-expanded="${e=>e.ariaExpanded}"
                    aria-haspopup="listbox"
                    class="selected-value"
                    part="selected-value"
                    placeholder="${e=>e.placeholder}"
                    role="combobox"
                    type="text"
                    ?disabled="${e=>e.disabled}"
                    :value="${e=>e.value}"
                    @input="${(e,t)=>e.inputHandler(t.event)}"
                    @keyup="${(e,t)=>e.keyupHandler(t.event)}"
                    ${(0,s.ref)("control")}
                />
                <div class="indicator" part="indicator" aria-hidden="true">
                    <slot name="indicator">
                        ${t.indicator||""}
                    </slot>
                </div>
            </slot>
            ${n(e,t)}
        </div>
        <div
            class="listbox"
            id="${e=>e.listboxId}"
            part="listbox"
            role="listbox"
            ?disabled="${e=>e.disabled}"
            ?hidden="${e=>!e.open}"
            ${(0,s.ref)("listbox")}
        >
            <slot
                ${(0,s.slotted)({filter:Ko.slottedOptionFilter,flatten:true,property:"slottedOptions"})}
            ></slot>
        </div>
    </template>
`;function tn(e){const t=e.tagFor(ao);return(0,s.html)`
    <${t}
        :rowData="${e=>e}"
        :cellItemTemplate="${(e,t)=>t.parent.cellItemTemplate}"
        :headerCellItemTemplate="${(e,t)=>t.parent.headerCellItemTemplate}"
    ></${t}>
`}const sn=(e,t)=>{const i=tn(e);const o=e.tagFor(ao);return(0,s.html)`
        <template
            role="grid"
            tabindex="0"
            :rowElementTag="${()=>o}"
            :defaultRowItemTemplate="${i}"
            ${(0,s.children)({property:"rowElements",filter:(0,s.elements)("[role=row]")})}
        >
            <slot></slot>
        </template>
    `};function on(e){const t=e.tagFor(ro);return(0,s.html)`
    <${t}
        cell-type="${e=>e.isRowHeader?"rowheader":undefined}"
        grid-column="${(e,t)=>t.index+1}"
        :rowData="${(e,t)=>t.parent.rowData}"
        :columnDefinition="${e=>e}"
    ></${t}>
`}function nn(e){const t=e.tagFor(ro);return(0,s.html)`
    <${t}
        cell-type="columnheader"
        grid-column="${(e,t)=>t.index+1}"
        :columnDefinition="${e=>e}"
    ></${t}>
`}const rn=(e,t)=>{const i=on(e);const o=nn(e);return(0,s.html)`
        <template
            role="row"
            class="${e=>e.rowType!=="default"?e.rowType:""}"
            :defaultCellItemTemplate="${i}"
            :defaultHeaderCellItemTemplate="${o}"
            ${(0,s.children)({property:"cellElements",filter:(0,s.elements)('[role="cell"],[role="gridcell"],[role="columnheader"],[role="rowheader"]')})}
        >
            <slot ${(0,s.slotted)("slottedCellElements")}></slot>
        </template>
    `};const an=(e,t)=>(0,s.html)`
        <template
            tabindex="-1"
            role="${e=>!e.cellType||e.cellType==="default"?"gridcell":e.cellType}"
            class="
            ${e=>e.cellType==="columnheader"?"column-header":e.cellType==="rowheader"?"row-header":""}
            "
        >
            <slot></slot>
        </template>
    `;function ln(e){const t=e.parentElement;if(t){return t}else{const t=e.getRootNode();if(t.host instanceof HTMLElement){return t.host}}return null}function dn(e,t){let i=t;while(i!==null){if(i===e){return true}i=ln(i)}return false}const hn=document.createElement("div");function cn(e){return e instanceof s.FASTElement}class un{setProperty(e,t){s.DOM.queueUpdate((()=>this.target.setProperty(e,t)))}removeProperty(e){s.DOM.queueUpdate((()=>this.target.removeProperty(e)))}}class pn extends un{constructor(e){super();const t=new CSSStyleSheet;this.target=t.cssRules[t.insertRule(":host{}")].style;e.$fastController.addStyles(s.ElementStyles.create([t]))}}class fn extends un{constructor(){super();const e=new CSSStyleSheet;this.target=e.cssRules[e.insertRule(":root{}")].style;document.adoptedStyleSheets=[...document.adoptedStyleSheets,e]}}class mn extends un{constructor(){super();this.style=document.createElement("style");document.head.appendChild(this.style);const{sheet:e}=this.style;if(e){const t=e.insertRule(":root{}",e.cssRules.length);this.target=e.cssRules[t].style}}}class vn{constructor(e){this.store=new Map;this.target=null;const t=e.$fastController;this.style=document.createElement("style");t.addStyles(this.style);s.Observable.getNotifier(t).subscribe(this,"isConnected");this.handleChange(t,"isConnected")}targetChanged(){if(this.target!==null){for(const[e,t]of this.store.entries()){this.target.setProperty(e,t)}}}setProperty(e,t){this.store.set(e,t);s.DOM.queueUpdate((()=>{if(this.target!==null){this.target.setProperty(e,t)}}))}removeProperty(e){this.store.delete(e);s.DOM.queueUpdate((()=>{if(this.target!==null){this.target.removeProperty(e)}}))}handleChange(e,t){const{sheet:i}=this.style;if(i){const e=i.insertRule(":host{}",i.cssRules.length);this.target=i.cssRules[e].style}else{this.target=null}}}f([s.observable],vn.prototype,"target",void 0);class bn{constructor(e){this.target=e.style}setProperty(e,t){s.DOM.queueUpdate((()=>this.target.setProperty(e,t)))}removeProperty(e){s.DOM.queueUpdate((()=>this.target.removeProperty(e)))}}class gn{setProperty(e,t){gn.properties[e]=t;for(const i of gn.roots.values()){xn.getOrCreate(gn.normalizeRoot(i)).setProperty(e,t)}}removeProperty(e){delete gn.properties[e];for(const t of gn.roots.values()){xn.getOrCreate(gn.normalizeRoot(t)).removeProperty(e)}}static registerRoot(e){const{roots:t}=gn;if(!t.has(e)){t.add(e);const i=xn.getOrCreate(this.normalizeRoot(e));for(const e in gn.properties){i.setProperty(e,gn.properties[e])}}}static unregisterRoot(e){const{roots:t}=gn;if(t.has(e)){t.delete(e);const i=xn.getOrCreate(gn.normalizeRoot(e));for(const e in gn.properties){i.removeProperty(e)}}}static normalizeRoot(e){return e===hn?document:e}}gn.roots=new Set;gn.properties={};const yn=new WeakMap;const Cn=s.DOM.supportsAdoptedStyleSheets?pn:vn;const xn=Object.freeze({getOrCreate(e){if(yn.has(e)){return yn.get(e)}let t;if(e===hn){t=new gn}else if(e instanceof Document){t=s.DOM.supportsAdoptedStyleSheets?new fn:new mn}else if(cn(e)){t=new Cn(e)}else{t=new bn(e)}yn.set(e,t);return t}});class wn extends s.CSSDirective{constructor(e){super();this.subscribers=new WeakMap;this._appliedTo=new Set;this.name=e.name;if(e.cssCustomPropertyName!==null){this.cssCustomProperty=`--${e.cssCustomPropertyName}`;this.cssVar=`var(${this.cssCustomProperty})`}this.id=wn.uniqueId();wn.tokensById.set(this.id,this)}get appliedTo(){return[...this._appliedTo]}static from(e){return new wn({name:typeof e==="string"?e:e.name,cssCustomPropertyName:typeof e==="string"?e:e.cssCustomPropertyName===void 0?e.name:e.cssCustomPropertyName})}static isCSSDesignToken(e){return typeof e.cssCustomProperty==="string"}static isDerivedDesignTokenValue(e){return typeof e==="function"}static getTokenById(e){return wn.tokensById.get(e)}getOrCreateSubscriberSet(e=this){return this.subscribers.get(e)||this.subscribers.set(e,new Set)&&this.subscribers.get(e)}createCSS(){return this.cssVar||""}getValueFor(e){const t=En.getOrCreate(e).get(this);if(t!==undefined){return t}throw new Error(`Value could not be retrieved for token named "${this.name}". Ensure the value is set for ${e} or an ancestor of ${e}.`)}setValueFor(e,t){this._appliedTo.add(e);if(t instanceof wn){t=this.alias(t)}En.getOrCreate(e).set(this,t);return this}deleteValueFor(e){this._appliedTo.delete(e);if(En.existsFor(e)){En.getOrCreate(e).delete(this)}return this}withDefault(e){this.setValueFor(hn,e);return this}subscribe(e,t){const i=this.getOrCreateSubscriberSet(t);if(t&&!En.existsFor(t)){En.getOrCreate(t)}if(!i.has(e)){i.add(e)}}unsubscribe(e,t){const i=this.subscribers.get(t||this);if(i&&i.has(e)){i.delete(e)}}notify(e){const t=Object.freeze({token:this,target:e});if(this.subscribers.has(this)){this.subscribers.get(this).forEach((e=>e.handleChange(t)))}if(this.subscribers.has(e)){this.subscribers.get(e).forEach((e=>e.handleChange(t)))}}alias(e){return t=>e.getValueFor(t)}}wn.uniqueId=(()=>{let e=0;return()=>{e++;return e.toString(16)}})();wn.tokensById=new Map;class $n{startReflection(e,t){e.subscribe(this,t);this.handleChange({token:e,target:t})}stopReflection(e,t){e.unsubscribe(this,t);this.remove(e,t)}handleChange(e){const{token:t,target:i}=e;this.add(t,i)}add(e,t){xn.getOrCreate(t).setProperty(e.cssCustomProperty,this.resolveCSSValue(En.getOrCreate(t).get(e)))}remove(e,t){xn.getOrCreate(t).removeProperty(e.cssCustomProperty)}resolveCSSValue(e){return e&&typeof e.createCSS==="function"?e.createCSS():e}}class In{constructor(e,t,i){this.source=e;this.token=t;this.node=i;this.dependencies=new Set;this.observer=s.Observable.binding(e,this,false);this.observer.handleChange=this.observer.call;this.handleChange()}disconnect(){this.observer.disconnect()}handleChange(){this.node.store.set(this.token,this.observer.observe(this.node.target,s.defaultExecutionContext))}}class kn{constructor(){this.values=new Map}set(e,t){if(this.values.get(e)!==t){this.values.set(e,t);s.Observable.getNotifier(this).notify(e.id)}}get(e){s.Observable.track(this,e.id);return this.values.get(e)}delete(e){this.values.delete(e)}all(){return this.values.entries()}}const On=new WeakMap;const Tn=new WeakMap;class En{constructor(e){this.target=e;this.store=new kn;this.children=[];this.assignedValues=new Map;this.reflecting=new Set;this.bindingObservers=new Map;this.tokenValueChangeHandler={handleChange:(e,t)=>{const i=wn.getTokenById(t);if(i){i.notify(this.target);if(wn.isCSSDesignToken(i)){const t=this.parent;const s=this.isReflecting(i);if(t){const o=t.get(i);const n=e.get(i);if(o!==n&&!s){this.reflectToCSS(i)}else if(o===n&&s){this.stopReflectToCSS(i)}}else if(!s){this.reflectToCSS(i)}}}}};On.set(e,this);s.Observable.getNotifier(this.store).subscribe(this.tokenValueChangeHandler);if(e instanceof s.FASTElement){e.$fastController.addBehaviors([this])}else if(e.isConnected){this.bind()}}static getOrCreate(e){return On.get(e)||new En(e)}static existsFor(e){return On.has(e)}static findParent(e){if(!(hn===e.target)){let t=ln(e.target);while(t!==null){if(On.has(t)){return On.get(t)}t=ln(t)}return En.getOrCreate(hn)}return null}static findClosestAssignedNode(e,t){let i=t;do{if(i.has(e)){return i}i=i.parent?i.parent:i.target!==hn?En.getOrCreate(hn):null}while(i!==null);return null}get parent(){return Tn.get(this)||null}has(e){return this.assignedValues.has(e)}get(e){const t=this.store.get(e);if(t!==undefined){return t}const i=this.getRaw(e);if(i!==undefined){this.hydrate(e,i);return this.get(e)}}getRaw(e){var t;if(this.assignedValues.has(e)){return this.assignedValues.get(e)}return(t=En.findClosestAssignedNode(e,this))===null||t===void 0?void 0:t.getRaw(e)}set(e,t){if(wn.isDerivedDesignTokenValue(this.assignedValues.get(e))){this.tearDownBindingObserver(e)}this.assignedValues.set(e,t);if(wn.isDerivedDesignTokenValue(t)){this.setupBindingObserver(e,t)}else{this.store.set(e,t)}}delete(e){this.assignedValues.delete(e);this.tearDownBindingObserver(e);const t=this.getRaw(e);if(t){this.hydrate(e,t)}else{this.store.delete(e)}}bind(){const e=En.findParent(this);if(e){e.appendChild(this)}for(const t of this.assignedValues.keys()){t.notify(this.target)}}unbind(){if(this.parent){const e=Tn.get(this);e.removeChild(this)}}appendChild(e){if(e.parent){Tn.get(e).removeChild(e)}const t=this.children.filter((t=>e.contains(t)));Tn.set(e,this);this.children.push(e);t.forEach((t=>e.appendChild(t)));s.Observable.getNotifier(this.store).subscribe(e);for(const[i,s]of this.store.all()){e.hydrate(i,this.bindingObservers.has(i)?this.getRaw(i):s)}}removeChild(e){const t=this.children.indexOf(e);if(t!==-1){this.children.splice(t,1)}s.Observable.getNotifier(this.store).unsubscribe(e);return e.parent===this?Tn.delete(e):false}contains(e){return dn(this.target,e.target)}reflectToCSS(e){if(!this.isReflecting(e)){this.reflecting.add(e);En.cssCustomPropertyReflector.startReflection(e,this.target)}}stopReflectToCSS(e){if(this.isReflecting(e)){this.reflecting.delete(e);En.cssCustomPropertyReflector.stopReflection(e,this.target)}}isReflecting(e){return this.reflecting.has(e)}handleChange(e,t){const i=wn.getTokenById(t);if(!i){return}this.hydrate(i,this.getRaw(i))}hydrate(e,t){if(!this.has(e)){const i=this.bindingObservers.get(e);if(wn.isDerivedDesignTokenValue(t)){if(i){if(i.source!==t){this.tearDownBindingObserver(e);this.setupBindingObserver(e,t)}}else{this.setupBindingObserver(e,t)}}else{if(i){this.tearDownBindingObserver(e)}this.store.set(e,t)}}}setupBindingObserver(e,t){const i=new In(t,e,this);this.bindingObservers.set(e,i);return i}tearDownBindingObserver(e){if(this.bindingObservers.has(e)){this.bindingObservers.get(e).disconnect();this.bindingObservers.delete(e);return true}return false}}En.cssCustomPropertyReflector=new $n;f([s.observable],En.prototype,"children",void 0);function Rn(e){return wn.from(e)}const Dn=Object.freeze({create:Rn,notifyConnection(e){if(!e.isConnected||!En.existsFor(e)){return false}En.getOrCreate(e).bind();return true},notifyDisconnection(e){if(e.isConnected||!En.existsFor(e)){return false}En.getOrCreate(e).unbind();return true},registerRoot(e=hn){gn.registerRoot(e)},unregisterRoot(e=hn){gn.unregisterRoot(e)}});const Sn=Object.freeze({definitionCallbackOnly:null,ignoreDuplicate:Symbol()});const An=new Map;const Fn=new Map;let Ln=null;const Mn=q.createInterface((e=>e.cachedCallback((e=>{if(Ln===null){Ln=new Vn(null,e)}return Ln}))));const Pn=Object.freeze({tagFor(e){return Fn.get(e)},responsibleFor(e){const t=e.$$designSystem$$;if(t){return t}const i=q.findResponsibleContainer(e);return i.get(Mn)},getOrCreate(e){if(!e){if(Ln===null){Ln=q.getOrCreateDOMContainer().get(Mn)}return Ln}const t=e.$$designSystem$$;if(t){return t}const i=q.getOrCreateDOMContainer(e);if(i.has(Mn,false)){return i.get(Mn)}else{const t=new Vn(e,i);i.register(xe.instance(Mn,t));return t}}});function Hn(e,t,i){if(typeof e==="string"){return{name:e,type:t,callback:i}}else{return e}}class Vn{constructor(e,t){this.owner=e;this.container=t;this.designTokensInitialized=false;this.prefix="fast";this.shadowRootMode=undefined;this.disambiguate=()=>Sn.definitionCallbackOnly;if(e!==null){e.$$designSystem$$=this}}withPrefix(e){this.prefix=e;return this}withShadowRootMode(e){this.shadowRootMode=e;return this}withElementDisambiguation(e){this.disambiguate=e;return this}withDesignTokenRoot(e){this.designTokenRoot=e;return this}register(...e){const t=this.container;const i=[];const s=this.disambiguate;const o=this.shadowRootMode;const n={elementPrefix:this.prefix,tryDefineElement(e,n,r){const a=Hn(e,n,r);const{name:l,callback:d,baseClass:h}=a;let{type:c}=a;let u=l;let p=An.get(u);let f=true;while(p){const e=s(u,c,p);switch(e){case Sn.ignoreDuplicate:return;case Sn.definitionCallbackOnly:f=false;p=void 0;break;default:u=e;p=An.get(u);break}}if(f){if(Fn.has(c)||c===Fe){c=class extends c{}}An.set(u,c);Fn.set(c,u);if(h){Fn.set(h,u)}}i.push(new zn(t,u,c,o,d,f))}};if(!this.designTokensInitialized){this.designTokensInitialized=true;if(this.designTokenRoot!==null){Dn.registerRoot(this.designTokenRoot)}}t.registerWithContext(n,...e);for(const r of i){r.callback(r);if(r.willDefine&&r.definition!==null){r.definition.define()}}return this}}class zn{constructor(e,t,i,s,o,n){this.container=e;this.name=t;this.type=i;this.shadowRootMode=s;this.callback=o;this.willDefine=n;this.definition=null}definePresentation(e){Se.define(this.name,e,this.container)}defineElement(e){this.definition=new s.FASTElementDefinition(this.type,Object.assign(Object.assign({},e),{name:this.name}))}tagFor(e){return Pn.tagFor(e)}}const Nn=(e,t)=>(0,s.html)`
    <div class="positioning-region" part="positioning-region">
        ${(0,s.when)((e=>e.modal),(0,s.html)`
                <div
                    class="overlay"
                    part="overlay"
                    role="presentation"
                    @click="${e=>e.dismiss()}"
                ></div>
            `)}
        <div
            role="dialog"
            tabindex="-1"
            class="control"
            part="control"
            aria-modal="${e=>e.modal}"
            aria-describedby="${e=>e.ariaDescribedby}"
            aria-labelledby="${e=>e.ariaLabelledby}"
            aria-label="${e=>e.ariaLabel}"
            ${(0,s.ref)("dialog")}
        >
            <slot></slot>
        </div>
    </div>
`;var Bn=i(49054);class qn extends Fe{constructor(){super(...arguments);this.modal=true;this.hidden=false;this.trapFocus=true;this.trapFocusChanged=()=>{if(this.$fastController.isConnected){this.updateTrapFocus()}};this.isTrappingFocus=false;this.handleDocumentKeydown=e=>{if(!e.defaultPrevented&&!this.hidden){switch(e.key){case ze.F9:this.dismiss();e.preventDefault();break;case ze.J9:this.handleTabKeyDown(e);break}}};this.handleDocumentFocus=e=>{if(!e.defaultPrevented&&this.shouldForceFocus(e.target)){this.focusFirstElement();e.preventDefault()}};this.handleTabKeyDown=e=>{if(!this.trapFocus||this.hidden){return}const t=this.getTabQueueBounds();if(t.length===0){return}if(t.length===1){t[0].focus();e.preventDefault();return}if(e.shiftKey&&e.target===t[0]){t[t.length-1].focus();e.preventDefault()}else if(!e.shiftKey&&e.target===t[t.length-1]){t[0].focus();e.preventDefault()}return};this.getTabQueueBounds=()=>{const e=[];return qn.reduceTabbableItems(e,this)};this.focusFirstElement=()=>{const e=this.getTabQueueBounds();if(e.length>0){e[0].focus()}else{if(this.dialog instanceof HTMLElement){this.dialog.focus()}}};this.shouldForceFocus=e=>this.isTrappingFocus&&!this.contains(e);this.shouldTrapFocus=()=>this.trapFocus&&!this.hidden;this.updateTrapFocus=e=>{const t=e===undefined?this.shouldTrapFocus():e;if(t&&!this.isTrappingFocus){this.isTrappingFocus=true;document.addEventListener("focusin",this.handleDocumentFocus);s.DOM.queueUpdate((()=>{if(this.shouldForceFocus(document.activeElement)){this.focusFirstElement()}}))}else if(!t&&this.isTrappingFocus){this.isTrappingFocus=false;document.removeEventListener("focusin",this.handleDocumentFocus)}}}dismiss(){this.$emit("dismiss");this.$emit("cancel")}show(){this.hidden=false}hide(){this.hidden=true;this.$emit("close")}connectedCallback(){super.connectedCallback();document.addEventListener("keydown",this.handleDocumentKeydown);this.notifier=s.Observable.getNotifier(this);this.notifier.subscribe(this,"hidden");this.updateTrapFocus()}disconnectedCallback(){super.disconnectedCallback();document.removeEventListener("keydown",this.handleDocumentKeydown);this.updateTrapFocus(false);this.notifier.unsubscribe(this,"hidden")}handleChange(e,t){switch(t){case"hidden":this.updateTrapFocus();break;default:break}}static reduceTabbableItems(e,t){if(t.getAttribute("tabindex")==="-1"){return e}if((0,Bn.AO)(t)||qn.isFocusableFastElement(t)&&qn.hasTabbableShadow(t)){e.push(t);return e}if(t.childElementCount){return e.concat(Array.from(t.children).reduce(qn.reduceTabbableItems,[]))}return e}static isFocusableFastElement(e){var t,i;return!!((i=(t=e.$fastController)===null||t===void 0?void 0:t.definition.shadowOptions)===null||i===void 0?void 0:i.delegatesFocus)}static hasTabbableShadow(e){var t,i;return Array.from((i=(t=e.shadowRoot)===null||t===void 0?void 0:t.querySelectorAll("*"))!==null&&i!==void 0?i:[]).some((e=>(0,Bn.AO)(e)))}}f([(0,s.attr)({mode:"boolean"})],qn.prototype,"modal",void 0);f([(0,s.attr)({mode:"boolean"})],qn.prototype,"hidden",void 0);f([(0,s.attr)({attribute:"trap-focus",mode:"boolean"})],qn.prototype,"trapFocus",void 0);f([(0,s.attr)({attribute:"aria-describedby"})],qn.prototype,"ariaDescribedby",void 0);f([(0,s.attr)({attribute:"aria-labelledby"})],qn.prototype,"ariaLabelledby",void 0);f([(0,s.attr)({attribute:"aria-label"})],qn.prototype,"ariaLabel",void 0);const Un=new MutationObserver((e=>{for(const t of e){jn.getOrCreateFor(t.target).notify(t.attributeName)}}));class jn extends s.SubscriberSet{constructor(e){super(e);this.watchedAttributes=new Set;jn.subscriberCache.set(e,this)}subscribe(e){super.subscribe(e);if(!this.watchedAttributes.has(e.attributes)){this.watchedAttributes.add(e.attributes);this.observe()}}unsubscribe(e){super.unsubscribe(e);if(this.watchedAttributes.has(e.attributes)){this.watchedAttributes.delete(e.attributes);this.observe()}}static getOrCreateFor(e){return this.subscriberCache.get(e)||new jn(e)}observe(){const e=[];for(const t of this.watchedAttributes.values()){for(let i=0;i<t.length;i++){e.push(t[i])}}Un.observe(this.source,{attributeFilter:e})}}jn.subscriberCache=new WeakMap;class _n{constructor(e,t){this.target=e;this.attributes=Object.freeze(t)}bind(e){jn.getOrCreateFor(e).subscribe(this);if(e.hasAttributes()){for(let t=0;t<e.attributes.length;t++){this.handleChange(e,e.attributes[t].name)}}}unbind(e){jn.getOrCreateFor(e).unsubscribe(this)}handleChange(e,t){if(this.attributes.includes(t)){s.DOM.setAttribute(this.target,t,e.getAttribute(t))}}}function Kn(...e){return new s.AttachedBehaviorHTMLDirective("fast-reflect-attr",_n,e)}const Wn=(e,t)=>(0,s.html)`
    <details class="disclosure" ${(0,s.ref)("details")}>
        <summary
            class="invoker"
            role="button"
            aria-controls="disclosure-content"
            aria-expanded="${e=>e.expanded}"
        >
            <slot name="start"></slot>
            <slot name="title">${e=>e.title}</slot>
            <slot name="end"></slot>
        </summary>
        <div id="disclosure-content"><slot></slot></div>
    </details>
`;class Gn extends Fe{connectedCallback(){super.connectedCallback();this.setup()}disconnectedCallback(){super.disconnectedCallback();this.details.removeEventListener("toggle",this.onToggle)}show(){this.details.open=true}hide(){this.details.open=false}toggle(){this.details.open=!this.details.open}setup(){this.onToggle=this.onToggle.bind(this);this.details.addEventListener("toggle",this.onToggle);if(this.expanded){this.show()}}onToggle(){this.expanded=this.details.open;this.$emit("toggle")}}f([(0,s.attr)({mode:"boolean"})],Gn.prototype,"expanded",void 0);f([s.attr],Gn.prototype,"title",void 0);const Xn=(e,t)=>(0,s.html)`
    <template role="${e=>e.role}" aria-orientation="${e=>e.orientation}"></template>
`;var Yn=i(67002);const Qn={separator:"separator",presentation:"presentation"};class Zn extends Fe{constructor(){super(...arguments);this.role=Qn.separator;this.orientation=Yn.t.horizontal}}f([s.attr],Zn.prototype,"role",void 0);f([s.attr],Zn.prototype,"orientation",void 0);const Jn={next:"next",previous:"previous"};const er=(e,t)=>(0,s.html)`
    <template
        role="button"
        aria-disabled="${e=>e.disabled?true:void 0}"
        tabindex="${e=>e.hiddenFromAT?-1:0}"
        class="${e=>e.direction} ${e=>e.disabled?"disabled":""}"
        @keyup="${(e,t)=>e.keyupHandler(t.event)}"
    >
        ${(0,s.when)((e=>e.direction===Jn.next),(0,s.html)`
                <span part="next" class="next">
                    <slot name="next">
                        ${t.next||""}
                    </slot>
                </span>
            `)}
        ${(0,s.when)((e=>e.direction===Jn.previous),(0,s.html)`
                <span part="previous" class="previous">
                    <slot name="previous">
                        ${t.previous||""}
                    </slot>
                </span>
            `)}
    </template>
`;class tr extends Fe{constructor(){super(...arguments);this.hiddenFromAT=true;this.direction=Jn.next}keyupHandler(e){if(!this.hiddenFromAT){const t=e.key;if(t==="Enter"||t==="Space"){this.$emit("click",e)}if(t==="Escape"){this.blur()}}}}f([(0,s.attr)({mode:"boolean"})],tr.prototype,"disabled",void 0);f([(0,s.attr)({attribute:"aria-hidden",converter:s.booleanConverter})],tr.prototype,"hiddenFromAT",void 0);f([s.attr],tr.prototype,"direction",void 0);const ir=(e,t)=>(0,s.html)`
    <template
        aria-checked="${e=>e.ariaChecked}"
        aria-disabled="${e=>e.ariaDisabled}"
        aria-posinset="${e=>e.ariaPosInSet}"
        aria-selected="${e=>e.ariaSelected}"
        aria-setsize="${e=>e.ariaSetSize}"
        class="${e=>[e.checked&&"checked",e.selected&&"selected",e.disabled&&"disabled"].filter(Boolean).join(" ")}"
        role="option"
    >
        ${r(e,t)}
        <span class="content" part="content">
            <slot ${(0,s.slotted)("content")}></slot>
        </span>
        ${n(e,t)}
    </template>
`;class sr extends Ko{constructor(){super(...arguments);this.activeIndex=-1;this.rangeStartIndex=-1}get activeOption(){return this.options[this.activeIndex]}get checkedOptions(){var e;return(e=this.options)===null||e===void 0?void 0:e.filter((e=>e.checked))}get firstSelectedOptionIndex(){return this.options.indexOf(this.firstSelectedOption)}activeIndexChanged(e,t){var i,s;this.ariaActiveDescendant=(s=(i=this.options[t])===null||i===void 0?void 0:i.id)!==null&&s!==void 0?s:"";this.focusAndScrollOptionIntoView()}checkActiveIndex(){if(!this.multiple){return}const e=this.activeOption;if(e){e.checked=true}}checkFirstOption(e=false){if(e){if(this.rangeStartIndex===-1){this.rangeStartIndex=this.activeIndex+1}this.options.forEach(((e,t)=>{e.checked=(0,Ne.r4)(t,this.rangeStartIndex)}))}else{this.uncheckAllOptions()}this.activeIndex=0;this.checkActiveIndex()}checkLastOption(e=false){if(e){if(this.rangeStartIndex===-1){this.rangeStartIndex=this.activeIndex}this.options.forEach(((e,t)=>{e.checked=(0,Ne.r4)(t,this.rangeStartIndex,this.options.length)}))}else{this.uncheckAllOptions()}this.activeIndex=this.options.length-1;this.checkActiveIndex()}connectedCallback(){super.connectedCallback();this.addEventListener("focusout",this.focusoutHandler)}disconnectedCallback(){this.removeEventListener("focusout",this.focusoutHandler);super.disconnectedCallback()}checkNextOption(e=false){if(e){if(this.rangeStartIndex===-1){this.rangeStartIndex=this.activeIndex}this.options.forEach(((e,t)=>{e.checked=(0,Ne.r4)(t,this.rangeStartIndex,this.activeIndex+1)}))}else{this.uncheckAllOptions()}this.activeIndex+=this.activeIndex<this.options.length-1?1:0;this.checkActiveIndex()}checkPreviousOption(e=false){if(e){if(this.rangeStartIndex===-1){this.rangeStartIndex=this.activeIndex}if(this.checkedOptions.length===1){this.rangeStartIndex+=1}this.options.forEach(((e,t)=>{e.checked=(0,Ne.r4)(t,this.activeIndex,this.rangeStartIndex)}))}else{this.uncheckAllOptions()}this.activeIndex-=this.activeIndex>0?1:0;this.checkActiveIndex()}clickHandler(e){var t;if(!this.multiple){return super.clickHandler(e)}const i=(t=e.target)===null||t===void 0?void 0:t.closest(`[role=option]`);if(!i||i.disabled){return}this.uncheckAllOptions();this.activeIndex=this.options.indexOf(i);this.checkActiveIndex();this.toggleSelectedForAllCheckedOptions();return true}focusAndScrollOptionIntoView(){super.focusAndScrollOptionIntoView(this.activeOption)}focusinHandler(e){if(!this.multiple){return super.focusinHandler(e)}if(!this.shouldSkipFocus&&e.target===e.currentTarget){this.uncheckAllOptions();if(this.activeIndex===-1){this.activeIndex=this.firstSelectedOptionIndex!==-1?this.firstSelectedOptionIndex:0}this.checkActiveIndex();this.setSelectedOptions();this.focusAndScrollOptionIntoView()}this.shouldSkipFocus=false}focusoutHandler(e){if(this.multiple){this.uncheckAllOptions()}}keydownHandler(e){if(!this.multiple){return super.keydownHandler(e)}if(this.disabled){return true}const{key:t,shiftKey:i}=e;this.shouldSkipFocus=false;switch(t){case ze.Tg:{this.checkFirstOption(i);return}case ze.HX:{this.checkNextOption(i);return}case ze.I5:{this.checkPreviousOption(i);return}case ze.FM:{this.checkLastOption(i);return}case ze.J9:{this.focusAndScrollOptionIntoView();return true}case ze.F9:{this.uncheckAllOptions();this.checkActiveIndex();return true}case ze.gG:{e.preventDefault();if(this.typeAheadExpired){this.toggleSelectedForAllCheckedOptions();return}}default:{if(t.length===1){this.handleTypeAhead(`${t}`)}return true}}}mousedownHandler(e){if(e.offsetX>=0&&e.offsetX<=this.scrollWidth){return super.mousedownHandler(e)}}multipleChanged(e,t){var i;this.ariaMultiSelectable=t?"true":null;(i=this.options)===null||i===void 0?void 0:i.forEach((e=>{e.checked=t?false:undefined}));this.setSelectedOptions()}setSelectedOptions(){if(!this.multiple){super.setSelectedOptions();return}if(this.$fastController.isConnected&&this.options){this.selectedOptions=this.options.filter((e=>e.selected));this.focusAndScrollOptionIntoView()}}sizeChanged(e,t){var i;const o=Math.max(0,parseInt((i=t===null||t===void 0?void 0:t.toFixed())!==null&&i!==void 0?i:"",10));if(o!==t){s.DOM.queueUpdate((()=>{this.size=o}))}}toggleSelectedForAllCheckedOptions(){const e=this.checkedOptions.filter((e=>!e.disabled));const t=!e.every((e=>e.selected));e.forEach((e=>e.selected=t));this.selectedIndex=this.options.indexOf(e[e.length-1]);this.setSelectedOptions()}typeaheadBufferChanged(e,t){if(!this.multiple){super.typeaheadBufferChanged(e,t);return}if(this.$fastController.isConnected){const e=this.getTypeaheadMatches();const t=this.options.indexOf(e[0]);if(t>-1){this.activeIndex=t;this.uncheckAllOptions();this.checkActiveIndex()}this.typeAheadExpired=false}}uncheckAllOptions(e=false){this.options.forEach((e=>e.checked=this.multiple?false:undefined));if(!e){this.rangeStartIndex=-1}}}f([s.observable],sr.prototype,"activeIndex",void 0);f([(0,s.attr)({mode:"boolean"})],sr.prototype,"multiple",void 0);f([(0,s.attr)({converter:s.nullableNumberConverter})],sr.prototype,"size",void 0);const or=(e,t)=>(0,s.html)`
    <template
        aria-activedescendant="${e=>e.ariaActiveDescendant}"
        aria-multiselectable="${e=>e.ariaMultiSelectable}"
        class="listbox"
        role="listbox"
        tabindex="${e=>!e.disabled?"0":null}"
        @click="${(e,t)=>e.clickHandler(t.event)}"
        @focusin="${(e,t)=>e.focusinHandler(t.event)}"
        @keydown="${(e,t)=>e.keydownHandler(t.event)}"
        @mousedown="${(e,t)=>e.mousedownHandler(t.event)}"
    >
        <slot
            ${(0,s.slotted)({filter:sr.slottedOptionFilter,flatten:true,property:"slottedOptions"})}
        ></slot>
    </template>
`;class nr extends Fe{constructor(){super(...arguments);this.optionElements=[]}menuElementsChanged(){this.updateOptions()}headerElementsChanged(){this.updateOptions()}footerElementsChanged(){this.updateOptions()}updateOptions(){this.optionElements.splice(0,this.optionElements.length);this.addSlottedListItems(this.headerElements);this.addSlottedListItems(this.menuElements);this.addSlottedListItems(this.footerElements);this.$emit("optionsupdated",{bubbles:false})}addSlottedListItems(e){if(e===undefined){return}e.forEach((e=>{if(e.nodeType===1&&e.getAttribute("role")==="listitem"){e.id=e.id||Io("option-");this.optionElements.push(e)}}))}}f([s.observable],nr.prototype,"menuElements",void 0);f([s.observable],nr.prototype,"headerElements",void 0);f([s.observable],nr.prototype,"footerElements",void 0);f([s.observable],nr.prototype,"suggestionsAvailableText",void 0);const rr=(0,s.html)`
    <template>
        ${e=>e.value}
    </template>
`;class ar extends Fe{contentsTemplateChanged(){if(this.$fastController.isConnected){this.updateView()}}connectedCallback(){super.connectedCallback();this.updateView()}disconnectedCallback(){super.disconnectedCallback();this.disconnectView()}handleClick(e){if(e.defaultPrevented){return false}this.handleInvoked();return false}handleInvoked(){this.$emit("pickeroptioninvoked")}updateView(){var e,t;this.disconnectView();this.customView=(t=(e=this.contentsTemplate)===null||e===void 0?void 0:e.render(this,this))!==null&&t!==void 0?t:rr.render(this,this)}disconnectView(){var e;(e=this.customView)===null||e===void 0?void 0:e.dispose();this.customView=undefined}}f([(0,s.attr)({attribute:"value"})],ar.prototype,"value",void 0);f([s.observable],ar.prototype,"contentsTemplate",void 0);class lr extends Fe{}const dr=(0,s.html)`
    <template>
        ${e=>e.value}
    </template>
`;class hr extends Fe{contentsTemplateChanged(){if(this.$fastController.isConnected){this.updateView()}}connectedCallback(){super.connectedCallback();this.updateView()}disconnectedCallback(){this.disconnectView();super.disconnectedCallback()}handleKeyDown(e){if(e.defaultPrevented){return false}if(e.key===ze.Mm){this.handleInvoke();return false}return true}handleClick(e){if(!e.defaultPrevented){this.handleInvoke()}return false}handleInvoke(){this.$emit("pickeriteminvoked")}updateView(){var e,t;this.disconnectView();this.customView=(t=(e=this.contentsTemplate)===null||e===void 0?void 0:e.render(this,this))!==null&&t!==void 0?t:dr.render(this,this)}disconnectView(){var e;(e=this.customView)===null||e===void 0?void 0:e.dispose();this.customView=undefined}}f([(0,s.attr)({attribute:"value"})],hr.prototype,"value",void 0);f([s.observable],hr.prototype,"contentsTemplate",void 0);function cr(e){const t=e.tagFor(hr);return(0,s.html)`
    <${t}
        value="${e=>e}"
        :contentsTemplate="${(e,t)=>t.parent.listItemContentsTemplate}"
    >
    </${t}>
    `}function ur(e){const t=e.tagFor(ar);return(0,s.html)`
    <${t}
        value="${e=>e}"
        :contentsTemplate="${(e,t)=>t.parent.menuOptionContentsTemplate}"
    >
    </${t}>
    `}const pr=(e,t)=>{const i=e.tagFor(Os);const o=e.tagFor(nr);const n=e.tagFor(lr);const r=e.tagFor(lr);const a=cr(e);const l=ur(e);return(0,s.html)`
        <template
            :selectedListTag="${()=>n}"
            :menuTag="${()=>o}"
            :defaultListItemTemplate="${a}"
            :defaultMenuOptionTemplate="${l}"
            @focusin="${(e,t)=>e.handleFocusIn(t.event)}"
            @focusout="${(e,t)=>e.handleFocusOut(t.event)}"
            @keydown="${(e,t)=>e.handleKeyDown(t.event)}"
            @pickeriteminvoked="${(e,t)=>e.handleItemInvoke(t.event)}"
            @pickeroptioninvoked="${(e,t)=>e.handleOptionInvoke(t.event)}"
        >
            <slot name="list-region"></slot>

            ${(0,s.when)((e=>e.flyoutOpen),(0,s.html)`
                <${i}
                    class="region"
                    part="region"
                    auto-update-mode="${e=>e.menuConfig.autoUpdateMode}"
                    fixed-placement="${e=>e.menuConfig.fixedPlacement}"
                    vertical-positioning-mode="${e=>e.menuConfig.verticalPositioningMode}"
                    vertical-default-position="${e=>e.menuConfig.verticalDefaultPosition}"
                    vertical-scaling="${e=>e.menuConfig.verticalScaling}"
                    vertical-inset="${e=>e.menuConfig.verticalInset}"
                    vertical-viewport-lock="${e=>e.menuConfig.verticalViewportLock}"
                    horizontal-positioning-mode="${e=>e.menuConfig.horizontalPositioningMode}"
                    horizontal-default-position="${e=>e.menuConfig.horizontalDefaultPosition}"
                    horizontal-scaling="${e=>e.menuConfig.horizontalScaling}"
                    horizontal-inset="${e=>e.menuConfig.horizontalInset}"
                    horizontal-viewport-lock="${e=>e.menuConfig.horizontalViewportLock}"
                    @loaded="${(e,t)=>e.handleRegionLoaded(t.event)}"
                    ${(0,s.ref)("region")}
                >
                    ${(0,s.when)((e=>!e.showNoOptions&&!e.showLoading),(0,s.html)`
                            <slot name="menu-region"></slot>
                        `)}
                    ${(0,s.when)((e=>e.showNoOptions&&!e.showLoading),(0,s.html)`
                            <div class="no-options-display" part="no-options-display">
                                <slot name="no-options-region">
                                    ${e=>e.noSuggestionsText}
                                </slot>
                            </div>
                        `)}
                    ${(0,s.when)((e=>e.showLoading),(0,s.html)`
                            <div class="loading-display" part="loading-display">
                                <slot name="loading-region">
                                    <${r}
                                        part="loading-progress"
                                        class="loading-progress
                                        slot="loading-region"
                                    ></${r}>
                                        ${e=>e.loadingText}
                                </slot>
                            </div>
                        `)}
                </${i}>
            `)}
        </template>
    `};class fr extends Fe{}class mr extends(Ws(fr)){constructor(){super(...arguments);this.proxy=document.createElement("input")}}const vr=(0,s.html)`
    <input
        slot="input-region"
        role="combobox"
        type="text"
        autocapitalize="off"
        autocomplete="off"
        haspopup="list"
        aria-label="${e=>e.label}"
        aria-labelledby="${e=>e.labelledBy}"
        placeholder="${e=>e.placeholder}"
        ${(0,s.ref)("inputElement")}
    ></input>
`;class br extends mr{constructor(){super(...arguments);this.selection="";this.filterSelected=true;this.filterQuery=true;this.noSuggestionsText="No suggestions available";this.suggestionsAvailableText="Suggestions available";this.loadingText="Loading suggestions";this.menuPlacement="bottom-fill";this.showLoading=false;this.optionsList=[];this.filteredOptionsList=[];this.flyoutOpen=false;this.menuFocusIndex=-1;this.showNoOptions=false;this.selectedItems=[];this.inputElementView=null;this.handleTextInput=e=>{this.query=this.inputElement.value};this.handleInputClick=e=>{e.preventDefault();this.toggleFlyout(true)};this.setRegionProps=()=>{if(!this.flyoutOpen){return}if(this.region===null||this.region===undefined){s.DOM.queueUpdate(this.setRegionProps);return}this.region.anchorElement=this.inputElement};this.configLookup={top:Es,bottom:Rs,tallest:Ds,"top-fill":Ss,"bottom-fill":As,"tallest-fill":Fs}}selectionChanged(){if(this.$fastController.isConnected){this.handleSelectionChange();if(this.proxy instanceof HTMLInputElement){this.proxy.value=this.selection;this.validate()}}}optionsChanged(){this.optionsList=this.options.split(",").map((e=>e.trim())).filter((e=>e!==""))}menuPlacementChanged(){if(this.$fastController.isConnected){this.updateMenuConfig()}}showLoadingChanged(){if(this.$fastController.isConnected){s.DOM.queueUpdate((()=>{this.setFocusedOption(0)}))}}listItemTemplateChanged(){this.updateListItemTemplate()}defaultListItemTemplateChanged(){this.updateListItemTemplate()}menuOptionTemplateChanged(){this.updateOptionTemplate()}defaultMenuOptionTemplateChanged(){this.updateOptionTemplate()}optionsListChanged(){this.updateFilteredOptions()}queryChanged(){if(this.$fastController.isConnected){if(this.inputElement.value!==this.query){this.inputElement.value=this.query}this.updateFilteredOptions();this.$emit("querychange",{bubbles:false})}}filteredOptionsListChanged(){if(this.$fastController.isConnected){this.showNoOptions=this.filteredOptionsList.length===0&&this.menuElement.querySelectorAll('[role="listitem"]').length===0;this.setFocusedOption(this.showNoOptions?-1:0)}}flyoutOpenChanged(){if(this.flyoutOpen){s.DOM.queueUpdate(this.setRegionProps);this.$emit("menuopening",{bubbles:false})}else{this.$emit("menuclosing",{bubbles:false})}}showNoOptionsChanged(){if(this.$fastController.isConnected){s.DOM.queueUpdate((()=>{this.setFocusedOption(0)}))}}connectedCallback(){super.connectedCallback();this.listElement=document.createElement(this.selectedListTag);this.appendChild(this.listElement);this.itemsPlaceholderElement=document.createComment("");this.listElement.append(this.itemsPlaceholderElement);this.inputElementView=vr.render(this,this.listElement);const e=this.menuTag.toUpperCase();this.menuElement=Array.from(this.children).find((t=>t.tagName===e));if(this.menuElement===undefined){this.menuElement=document.createElement(this.menuTag);this.appendChild(this.menuElement)}if(this.menuElement.id===""){this.menuElement.id=Io("listbox-")}this.menuId=this.menuElement.id;this.optionsPlaceholder=document.createComment("");this.menuElement.append(this.optionsPlaceholder);this.updateMenuConfig();s.DOM.queueUpdate((()=>this.initialize()))}disconnectedCallback(){super.disconnectedCallback();this.toggleFlyout(false);this.inputElement.removeEventListener("input",this.handleTextInput);this.inputElement.removeEventListener("click",this.handleInputClick);if(this.inputElementView!==null){this.inputElementView.dispose();this.inputElementView=null}}focus(){this.inputElement.focus()}initialize(){this.updateListItemTemplate();this.updateOptionTemplate();this.itemsRepeatBehavior=new s.RepeatDirective((e=>e.selectedItems),(e=>e.activeListItemTemplate),{positioning:true}).createBehavior(this.itemsPlaceholderElement);this.inputElement.addEventListener("input",this.handleTextInput);this.inputElement.addEventListener("click",this.handleInputClick);this.$fastController.addBehaviors([this.itemsRepeatBehavior]);this.menuElement.suggestionsAvailableText=this.suggestionsAvailableText;this.menuElement.addEventListener("optionsupdated",this.handleMenuOptionsUpdated);this.optionsRepeatBehavior=new s.RepeatDirective((e=>e.filteredOptionsList),(e=>e.activeMenuOptionTemplate),{positioning:true}).createBehavior(this.optionsPlaceholder);this.$fastController.addBehaviors([this.optionsRepeatBehavior]);this.handleSelectionChange()}toggleFlyout(e){if(this.flyoutOpen===e){return}if(e&&document.activeElement===this.inputElement){this.flyoutOpen=e;s.DOM.queueUpdate((()=>{if(this.menuElement!==undefined){this.setFocusedOption(0)}else{this.disableMenu()}}));return}this.flyoutOpen=false;this.disableMenu();return}handleMenuOptionsUpdated(e){e.preventDefault();if(this.flyoutOpen){this.setFocusedOption(0)}}handleKeyDown(e){if(e.defaultPrevented){return false}switch(e.key){case ze.HX:{if(!this.flyoutOpen){this.toggleFlyout(true)}else{const e=this.flyoutOpen?Math.min(this.menuFocusIndex+1,this.menuElement.optionElements.length-1):0;this.setFocusedOption(e)}return false}case ze.I5:{if(!this.flyoutOpen){this.toggleFlyout(true)}else{const e=this.flyoutOpen?Math.max(this.menuFocusIndex-1,0):0;this.setFocusedOption(e)}return false}case ze.F9:{this.toggleFlyout(false);return false}case ze.Mm:{if(this.menuFocusIndex!==-1&&this.menuElement.optionElements.length>this.menuFocusIndex){this.menuElement.optionElements[this.menuFocusIndex].click()}return false}case ze.bb:{if(document.activeElement!==this.inputElement){this.incrementFocusedItem(1);return false}return true}case ze.kT:{if(this.inputElement.selectionStart===0){this.incrementFocusedItem(-1);return false}return true}case ze.De:case ze.R9:{if(document.activeElement===null){return true}if(document.activeElement===this.inputElement){if(this.inputElement.selectionStart===0){this.selection=this.selectedItems.slice(0,this.selectedItems.length-1).toString();this.toggleFlyout(false);return false}return true}const e=Array.from(this.listElement.children);const t=e.indexOf(document.activeElement);if(t>-1){this.selection=this.selectedItems.splice(t,1).toString();s.DOM.queueUpdate((()=>{e[Math.min(e.length,t)].focus()}));return false}return true}}this.toggleFlyout(true);return true}handleFocusIn(e){return false}handleFocusOut(e){if(this.menuElement===undefined||!this.menuElement.contains(e.relatedTarget)){this.toggleFlyout(false)}return false}handleSelectionChange(){if(this.selectedItems.toString()===this.selection){return}this.selectedItems=this.selection===""?[]:this.selection.split(",");this.updateFilteredOptions();s.DOM.queueUpdate((()=>{this.checkMaxItems()}));this.$emit("selectionchange",{bubbles:false})}handleRegionLoaded(e){s.DOM.queueUpdate((()=>{this.setFocusedOption(0);this.$emit("menuloaded",{bubbles:false})}))}checkMaxItems(){if(this.inputElement===undefined){return}if(this.maxSelected!==undefined&&this.selectedItems.length>=this.maxSelected){if(document.activeElement===this.inputElement){const e=Array.from(this.listElement.querySelectorAll("[role='listitem']"));e[e.length-1].focus()}this.inputElement.hidden=true}else{this.inputElement.hidden=false}}handleItemInvoke(e){if(e.defaultPrevented){return false}if(e.target instanceof hr){const t=Array.from(this.listElement.querySelectorAll("[role='listitem']"));const i=t.indexOf(e.target);if(i!==-1){const e=this.selectedItems.slice();e.splice(i,1);this.selection=e.toString();s.DOM.queueUpdate((()=>this.incrementFocusedItem(0)))}return false}return true}handleOptionInvoke(e){if(e.defaultPrevented){return false}if(e.target instanceof ar){if(e.target.value!==undefined){this.selection=`${this.selection}${this.selection===""?"":","}${e.target.value}`}this.inputElement.value="";this.query="";this.inputElement.focus();this.toggleFlyout(false);return false}return true}incrementFocusedItem(e){if(this.selectedItems.length===0){this.inputElement.focus();return}const t=Array.from(this.listElement.querySelectorAll("[role='listitem']"));if(document.activeElement!==null){let i=t.indexOf(document.activeElement);if(i===-1){i=t.length}const s=Math.min(t.length,Math.max(0,i+e));if(s===t.length){if(this.maxSelected!==undefined&&this.selectedItems.length>=this.maxSelected){t[s-1].focus()}else{this.inputElement.focus()}}else{t[s].focus()}}}disableMenu(){var e,t,i;this.menuFocusIndex=-1;this.menuFocusOptionId=undefined;(e=this.inputElement)===null||e===void 0?void 0:e.removeAttribute("aria-activedescendant");(t=this.inputElement)===null||t===void 0?void 0:t.removeAttribute("aria-owns");(i=this.inputElement)===null||i===void 0?void 0:i.removeAttribute("aria-expanded")}setFocusedOption(e){if(!this.flyoutOpen||e===-1||this.showNoOptions||this.showLoading){this.disableMenu();return}if(this.menuElement.optionElements.length===0){return}this.menuElement.optionElements.forEach((e=>{e.setAttribute("aria-selected","false")}));this.menuFocusIndex=e;if(this.menuFocusIndex>this.menuElement.optionElements.length-1){this.menuFocusIndex=this.menuElement.optionElements.length-1}this.menuFocusOptionId=this.menuElement.optionElements[this.menuFocusIndex].id;this.inputElement.setAttribute("aria-owns",this.menuId);this.inputElement.setAttribute("aria-expanded","true");this.inputElement.setAttribute("aria-activedescendant",this.menuFocusOptionId);const t=this.menuElement.optionElements[this.menuFocusIndex];t.setAttribute("aria-selected","true");this.menuElement.scrollTo(0,t.offsetTop)}updateListItemTemplate(){var e;this.activeListItemTemplate=(e=this.listItemTemplate)!==null&&e!==void 0?e:this.defaultListItemTemplate}updateOptionTemplate(){var e;this.activeMenuOptionTemplate=(e=this.menuOptionTemplate)!==null&&e!==void 0?e:this.defaultMenuOptionTemplate}updateFilteredOptions(){this.filteredOptionsList=this.optionsList.slice(0);if(this.filterSelected){this.filteredOptionsList=this.filteredOptionsList.filter((e=>this.selectedItems.indexOf(e)===-1))}if(this.filterQuery&&this.query!==""&&this.query!==undefined){this.filteredOptionsList=this.filteredOptionsList.filter((e=>e.indexOf(this.query)!==-1))}}updateMenuConfig(){let e=this.configLookup[this.menuPlacement];if(e===null){e=As}this.menuConfig=Object.assign(Object.assign({},e),{autoUpdateMode:"auto",fixedPlacement:true,horizontalViewportLock:false,verticalViewportLock:false})}}f([(0,s.attr)({attribute:"selection"})],br.prototype,"selection",void 0);f([(0,s.attr)({attribute:"options"})],br.prototype,"options",void 0);f([(0,s.attr)({attribute:"filter-selected",mode:"boolean"})],br.prototype,"filterSelected",void 0);f([(0,s.attr)({attribute:"filter-query",mode:"boolean"})],br.prototype,"filterQuery",void 0);f([(0,s.attr)({attribute:"max-selected"})],br.prototype,"maxSelected",void 0);f([(0,s.attr)({attribute:"no-suggestions-text"})],br.prototype,"noSuggestionsText",void 0);f([(0,s.attr)({attribute:"suggestions-available-text"})],br.prototype,"suggestionsAvailableText",void 0);f([(0,s.attr)({attribute:"loading-text"})],br.prototype,"loadingText",void 0);f([(0,s.attr)({attribute:"label"})],br.prototype,"label",void 0);f([(0,s.attr)({attribute:"labelledby"})],br.prototype,"labelledBy",void 0);f([(0,s.attr)({attribute:"placeholder"})],br.prototype,"placeholder",void 0);f([(0,s.attr)({attribute:"menu-placement"})],br.prototype,"menuPlacement",void 0);f([s.observable],br.prototype,"showLoading",void 0);f([s.observable],br.prototype,"listItemTemplate",void 0);f([s.observable],br.prototype,"defaultListItemTemplate",void 0);f([s.observable],br.prototype,"activeListItemTemplate",void 0);f([s.observable],br.prototype,"menuOptionTemplate",void 0);f([s.observable],br.prototype,"defaultMenuOptionTemplate",void 0);f([s.observable],br.prototype,"activeMenuOptionTemplate",void 0);f([s.observable],br.prototype,"listItemContentsTemplate",void 0);f([s.observable],br.prototype,"menuOptionContentsTemplate",void 0);f([s.observable],br.prototype,"optionsList",void 0);f([s.observable],br.prototype,"query",void 0);f([s.observable],br.prototype,"filteredOptionsList",void 0);f([s.observable],br.prototype,"flyoutOpen",void 0);f([s.observable],br.prototype,"menuId",void 0);f([s.observable],br.prototype,"selectedListTag",void 0);f([s.observable],br.prototype,"menuTag",void 0);f([s.observable],br.prototype,"menuFocusIndex",void 0);f([s.observable],br.prototype,"menuFocusOptionId",void 0);f([s.observable],br.prototype,"showNoOptions",void 0);f([s.observable],br.prototype,"menuConfig",void 0);f([s.observable],br.prototype,"selectedItems",void 0);const gr=(e,t)=>(0,s.html)`
        <template role="list" slot="menu-region">
            <div class="options-display" part="options-display">
                <div class="header-region" part="header-region">
                    <slot name="header-region" ${(0,s.slotted)("headerElements")}></slot>
                </div>

                <slot ${(0,s.slotted)("menuElements")}></slot>
                <div class="footer-region" part="footer-region">
                    <slot name="footer-region" ${(0,s.slotted)("footerElements")}></slot>
                </div>
                <div
                    role="alert"
                    aria-live="polite"
                    part="suggestions-available-alert"
                    class="suggestions-available-alert"
                >
                    ${e=>e.suggestionsAvailableText}
                </div>
            </div>
        </template>
    `;const yr=(e,t)=>(0,s.html)`
        <template
            role="listitem"
            tabindex="-1"
            @click="${(e,t)=>e.handleClick(t.event)}"
        >
            <slot></slot>
        </template>
    `;const Cr=(e,t)=>(0,s.html)`
        <template slot="list-region" role="list" class="picker-list">
            <slot></slot>
            <slot name="input-region"></slot>
        </template>
    `;const xr=(e,t)=>(0,s.html)`
        <template
            role="listitem"
            tabindex="0"
            @click="${(e,t)=>e.handleClick(t.event)}"
            @keydown="${(e,t)=>e.handleKeyDown(t.event)}"
        >
            <slot></slot>
        </template>
    `;const wr={menuitem:"menuitem",menuitemcheckbox:"menuitemcheckbox",menuitemradio:"menuitemradio"};const $r={[wr.menuitem]:"menuitem",[wr.menuitemcheckbox]:"menuitemcheckbox",[wr.menuitemradio]:"menuitemradio"};const Ir=(e,t)=>(0,s.html)`
    <template
        role="${e=>e.role}"
        aria-haspopup="${e=>e.hasSubmenu?"menu":void 0}"
        aria-checked="${e=>e.role!==wr.menuitem?e.checked:void 0}"
        aria-disabled="${e=>e.disabled}"
        aria-expanded="${e=>e.expanded}"
        @keydown="${(e,t)=>e.handleMenuItemKeyDown(t.event)}"
        @click="${(e,t)=>e.handleMenuItemClick(t.event)}"
        @mouseover="${(e,t)=>e.handleMouseOver(t.event)}"
        @mouseout="${(e,t)=>e.handleMouseOut(t.event)}"
        class="${e=>e.disabled?"disabled":""} ${e=>e.expanded?"expanded":""} ${e=>`indent-${e.startColumnCount}`}"
    >
            ${(0,s.when)((e=>e.role===wr.menuitemcheckbox),(0,s.html)`
                    <div part="input-container" class="input-container">
                        <span part="checkbox" class="checkbox">
                            <slot name="checkbox-indicator">
                                ${t.checkboxIndicator||""}
                            </slot>
                        </span>
                    </div>
                `)}
            ${(0,s.when)((e=>e.role===wr.menuitemradio),(0,s.html)`
                    <div part="input-container" class="input-container">
                        <span part="radio" class="radio">
                            <slot name="radio-indicator">
                                ${t.radioIndicator||""}
                            </slot>
                        </span>
                    </div>
                `)}
        </div>
        ${r(e,t)}
        <span class="content" part="content">
            <slot></slot>
        </span>
        ${n(e,t)}
        ${(0,s.when)((e=>e.hasSubmenu),(0,s.html)`
                <div
                    part="expand-collapse-glyph-container"
                    class="expand-collapse-glyph-container"
                >
                    <span part="expand-collapse" class="expand-collapse">
                        <slot name="expand-collapse-indicator">
                            ${t.expandCollapseGlyph||""}
                        </slot>
                    </span>
                </div>
            `)}
        ${(0,s.when)((e=>e.expanded),(0,s.html)`
                <${e.tagFor(Os)}
                    :anchorElement="${e=>e}"
                    vertical-positioning-mode="dynamic"
                    vertical-default-position="bottom"
                    vertical-inset="true"
                    horizontal-positioning-mode="dynamic"
                    horizontal-default-position="end"
                    class="submenu-region"
                    dir="${e=>e.currentDirection}"
                    @loaded="${e=>e.submenuLoaded()}"
                    ${(0,s.ref)("submenuRegion")}
                    part="submenu-region"
                >
                    <slot name="submenu"></slot>
                </${e.tagFor(Os)}>
            `)}
    </template>
`;class kr extends Fe{constructor(){super(...arguments);this.role=wr.menuitem;this.hasSubmenu=false;this.currentDirection=Ge.O.ltr;this.focusSubmenuOnLoad=false;this.handleMenuItemKeyDown=e=>{if(e.defaultPrevented){return false}switch(e.key){case ze.Mm:case ze.gG:this.invoke();return false;case ze.bb:this.expandAndFocus();return false;case ze.kT:if(this.expanded){this.expanded=false;this.focus();return false}}return true};this.handleMenuItemClick=e=>{if(e.defaultPrevented||this.disabled){return false}this.invoke();return false};this.submenuLoaded=()=>{if(!this.focusSubmenuOnLoad){return}this.focusSubmenuOnLoad=false;if(this.hasSubmenu){this.submenu.focus();this.setAttribute("tabindex","-1")}};this.handleMouseOver=e=>{if(this.disabled||!this.hasSubmenu||this.expanded){return false}this.expanded=true;return false};this.handleMouseOut=e=>{if(!this.expanded||this.contains(document.activeElement)){return false}this.expanded=false;return false};this.expandAndFocus=()=>{if(!this.hasSubmenu){return}this.focusSubmenuOnLoad=true;this.expanded=true};this.invoke=()=>{if(this.disabled){return}switch(this.role){case wr.menuitemcheckbox:this.checked=!this.checked;break;case wr.menuitem:this.updateSubmenu();if(this.hasSubmenu){this.expandAndFocus()}else{this.$emit("change")}break;case wr.menuitemradio:if(!this.checked){this.checked=true}break}};this.updateSubmenu=()=>{this.submenu=this.domChildren().find((e=>e.getAttribute("role")==="menu"));this.hasSubmenu=this.submenu===undefined?false:true}}expandedChanged(e){if(this.$fastController.isConnected){if(this.submenu===undefined){return}if(this.expanded===false){this.submenu.collapseExpandedItem()}else{this.currentDirection=Is(this)}this.$emit("expanded-change",this,{bubbles:false})}}checkedChanged(e,t){if(this.$fastController.isConnected){this.$emit("change")}}connectedCallback(){super.connectedCallback();s.DOM.queueUpdate((()=>{this.updateSubmenu()}));if(!this.startColumnCount){this.startColumnCount=1}this.observer=new MutationObserver(this.updateSubmenu)}disconnectedCallback(){super.disconnectedCallback();this.submenu=undefined;if(this.observer!==undefined){this.observer.disconnect();this.observer=undefined}}domChildren(){return Array.from(this.children).filter((e=>!e.hasAttribute("hidden")))}}f([(0,s.attr)({mode:"boolean"})],kr.prototype,"disabled",void 0);f([(0,s.attr)({mode:"boolean"})],kr.prototype,"expanded",void 0);f([s.observable],kr.prototype,"startColumnCount",void 0);f([s.attr],kr.prototype,"role",void 0);f([(0,s.attr)({mode:"boolean"})],kr.prototype,"checked",void 0);f([s.observable],kr.prototype,"submenuRegion",void 0);f([s.observable],kr.prototype,"hasSubmenu",void 0);f([s.observable],kr.prototype,"currentDirection",void 0);f([s.observable],kr.prototype,"submenu",void 0);Pe(kr,o);const Or=(e,t)=>(0,s.html)`
    <template
        slot="${e=>e.slot?e.slot:e.isNestedMenu()?"submenu":void 0}"
        role="menu"
        @keydown="${(e,t)=>e.handleMenuKeyDown(t.event)}"
        @focusout="${(e,t)=>e.handleFocusOut(t.event)}"
    >
        <slot ${(0,s.slotted)("items")}></slot>
    </template>
`;class Tr extends Fe{constructor(){super(...arguments);this.expandedItem=null;this.focusIndex=-1;this.isNestedMenu=()=>this.parentElement!==null&&Ao(this.parentElement)&&this.parentElement.getAttribute("role")==="menuitem";this.handleFocusOut=e=>{if(!this.contains(e.relatedTarget)&&this.menuItems!==undefined){this.collapseExpandedItem();const e=this.menuItems.findIndex(this.isFocusableElement);this.menuItems[this.focusIndex].setAttribute("tabindex","-1");this.menuItems[e].setAttribute("tabindex","0");this.focusIndex=e}};this.handleItemFocus=e=>{const t=e.target;if(this.menuItems!==undefined&&t!==this.menuItems[this.focusIndex]){this.menuItems[this.focusIndex].setAttribute("tabindex","-1");this.focusIndex=this.menuItems.indexOf(t);t.setAttribute("tabindex","0")}};this.handleExpandedChanged=e=>{if(e.defaultPrevented||e.target===null||this.menuItems===undefined||this.menuItems.indexOf(e.target)<0){return}e.preventDefault();const t=e.target;if(this.expandedItem!==null&&t===this.expandedItem&&t.expanded===false){this.expandedItem=null;return}if(t.expanded){if(this.expandedItem!==null&&this.expandedItem!==t){this.expandedItem.expanded=false}this.menuItems[this.focusIndex].setAttribute("tabindex","-1");this.expandedItem=t;this.focusIndex=this.menuItems.indexOf(t);t.setAttribute("tabindex","0")}};this.removeItemListeners=()=>{if(this.menuItems!==undefined){this.menuItems.forEach((e=>{e.removeEventListener("expanded-change",this.handleExpandedChanged);e.removeEventListener("focus",this.handleItemFocus)}))}};this.setItems=()=>{const e=this.domChildren();this.removeItemListeners();this.menuItems=e;const t=this.menuItems.filter(this.isMenuItemElement);if(t.length){this.focusIndex=0}function i(e){const t=e.getAttribute("role");const i=e.querySelector("[slot=start]");if(t!==wr.menuitem&&i===null){return 1}else if(t===wr.menuitem&&i!==null){return 1}else if(t!==wr.menuitem&&i!==null){return 2}else{return 0}}const s=t.reduce(((e,t)=>{const s=i(t);return e>s?e:s}),0);t.forEach(((e,t)=>{e.setAttribute("tabindex",t===0?"0":"-1");e.addEventListener("expanded-change",this.handleExpandedChanged);e.addEventListener("focus",this.handleItemFocus);if(e instanceof kr){e.startColumnCount=s}}))};this.changeHandler=e=>{if(this.menuItems===undefined){return}const t=e.target;const i=this.menuItems.indexOf(t);if(i===-1){return}if(t.role==="menuitemradio"&&t.checked===true){for(let t=i-1;t>=0;--t){const e=this.menuItems[t];const i=e.getAttribute("role");if(i===wr.menuitemradio){e.checked=false}if(i==="separator"){break}}const e=this.menuItems.length-1;for(let t=i+1;t<=e;++t){const e=this.menuItems[t];const i=e.getAttribute("role");if(i===wr.menuitemradio){e.checked=false}if(i==="separator"){break}}}};this.isMenuItemElement=e=>Ao(e)&&Tr.focusableElementRoles.hasOwnProperty(e.getAttribute("role"));this.isFocusableElement=e=>this.isMenuItemElement(e)}itemsChanged(e,t){if(this.$fastController.isConnected&&this.menuItems!==undefined){this.setItems()}}connectedCallback(){super.connectedCallback();s.DOM.queueUpdate((()=>{this.setItems()}));this.addEventListener("change",this.changeHandler)}disconnectedCallback(){super.disconnectedCallback();this.removeItemListeners();this.menuItems=undefined;this.removeEventListener("change",this.changeHandler)}focus(){this.setFocus(0,1)}collapseExpandedItem(){if(this.expandedItem!==null){this.expandedItem.expanded=false;this.expandedItem=null}}handleMenuKeyDown(e){if(e.defaultPrevented||this.menuItems===undefined){return}switch(e.key){case ze.HX:this.setFocus(this.focusIndex+1,1);return;case ze.I5:this.setFocus(this.focusIndex-1,-1);return;case ze.FM:this.setFocus(this.menuItems.length-1,-1);return;case ze.Tg:this.setFocus(0,1);return;default:return true}}domChildren(){return Array.from(this.children).filter((e=>!e.hasAttribute("hidden")))}setFocus(e,t){if(this.menuItems===undefined){return}while(e>=0&&e<this.menuItems.length){const i=this.menuItems[e];if(this.isFocusableElement(i)){if(this.focusIndex>-1&&this.menuItems.length>=this.focusIndex-1){this.menuItems[this.focusIndex].setAttribute("tabindex","-1")}this.focusIndex=e;i.setAttribute("tabindex","0");i.focus();break}e+=t}}}Tr.focusableElementRoles=$r;f([s.observable],Tr.prototype,"items",void 0);const Er=(e,t)=>(0,s.html)`
    <template class="${e=>e.readOnly?"readonly":""}">
        <label
            part="label"
            for="control"
            class="${e=>e.defaultSlottedNodes&&e.defaultSlottedNodes.length?"label":"label label__hidden"}"
        >
            <slot ${(0,s.slotted)("defaultSlottedNodes")}></slot>
        </label>
        <div class="root" part="root">
            ${r(e,t)}
            <input
                class="control"
                part="control"
                id="control"
                @input="${e=>e.handleTextInput()}"
                @change="${e=>e.handleChange()}"
                @keydown="${(e,t)=>e.handleKeyDown(t.event)}"
                @blur="${(e,t)=>e.handleBlur()}"
                ?autofocus="${e=>e.autofocus}"
                ?disabled="${e=>e.disabled}"
                list="${e=>e.list}"
                maxlength="${e=>e.maxlength}"
                minlength="${e=>e.minlength}"
                placeholder="${e=>e.placeholder}"
                ?readonly="${e=>e.readOnly}"
                ?required="${e=>e.required}"
                size="${e=>e.size}"
                type="text"
                inputmode="numeric"
                min="${e=>e.min}"
                max="${e=>e.max}"
                step="${e=>e.step}"
                aria-atomic="${e=>e.ariaAtomic}"
                aria-busy="${e=>e.ariaBusy}"
                aria-controls="${e=>e.ariaControls}"
                aria-current="${e=>e.ariaCurrent}"
                aria-describedby="${e=>e.ariaDescribedby}"
                aria-details="${e=>e.ariaDetails}"
                aria-disabled="${e=>e.ariaDisabled}"
                aria-errormessage="${e=>e.ariaErrormessage}"
                aria-flowto="${e=>e.ariaFlowto}"
                aria-haspopup="${e=>e.ariaHaspopup}"
                aria-hidden="${e=>e.ariaHidden}"
                aria-invalid="${e=>e.ariaInvalid}"
                aria-keyshortcuts="${e=>e.ariaKeyshortcuts}"
                aria-label="${e=>e.ariaLabel}"
                aria-labelledby="${e=>e.ariaLabelledby}"
                aria-live="${e=>e.ariaLive}"
                aria-owns="${e=>e.ariaOwns}"
                aria-relevant="${e=>e.ariaRelevant}"
                aria-roledescription="${e=>e.ariaRoledescription}"
                ${(0,s.ref)("control")}
            />
            ${(0,s.when)((e=>!e.hideStep&&!e.readOnly&&!e.disabled),(0,s.html)`
                    <div class="controls" part="controls">
                        <div class="step-up" part="step-up" @click="${e=>e.stepUp()}">
                            <slot name="step-up-glyph">
                                ${t.stepUpGlyph||""}
                            </slot>
                        </div>
                        <div
                            class="step-down"
                            part="step-down"
                            @click="${e=>e.stepDown()}"
                        >
                            <slot name="step-down-glyph">
                                ${t.stepDownGlyph||""}
                            </slot>
                        </div>
                    </div>
                `)}
            ${n(e,t)}
        </div>
    </template>
`;class Rr extends Fe{}class Dr extends(Ws(Rr)){constructor(){super(...arguments);this.proxy=document.createElement("input")}}const Sr={email:"email",password:"password",tel:"tel",text:"text",url:"url"};class Ar extends Dr{constructor(){super(...arguments);this.type=Sr.text}readOnlyChanged(){if(this.proxy instanceof HTMLInputElement){this.proxy.readOnly=this.readOnly;this.validate()}}autofocusChanged(){if(this.proxy instanceof HTMLInputElement){this.proxy.autofocus=this.autofocus;this.validate()}}placeholderChanged(){if(this.proxy instanceof HTMLInputElement){this.proxy.placeholder=this.placeholder}}typeChanged(){if(this.proxy instanceof HTMLInputElement){this.proxy.type=this.type;this.validate()}}listChanged(){if(this.proxy instanceof HTMLInputElement){this.proxy.setAttribute("list",this.list);this.validate()}}maxlengthChanged(){if(this.proxy instanceof HTMLInputElement){this.proxy.maxLength=this.maxlength;this.validate()}}minlengthChanged(){if(this.proxy instanceof HTMLInputElement){this.proxy.minLength=this.minlength;this.validate()}}patternChanged(){if(this.proxy instanceof HTMLInputElement){this.proxy.pattern=this.pattern;this.validate()}}sizeChanged(){if(this.proxy instanceof HTMLInputElement){this.proxy.size=this.size}}spellcheckChanged(){if(this.proxy instanceof HTMLInputElement){this.proxy.spellcheck=this.spellcheck}}connectedCallback(){super.connectedCallback();this.proxy.setAttribute("type",this.type);this.validate();if(this.autofocus){s.DOM.queueUpdate((()=>{this.focus()}))}}select(){this.control.select();this.$emit("select")}handleTextInput(){this.value=this.control.value}handleChange(){this.$emit("change")}validate(){super.validate(this.control)}}f([(0,s.attr)({attribute:"readonly",mode:"boolean"})],Ar.prototype,"readOnly",void 0);f([(0,s.attr)({mode:"boolean"})],Ar.prototype,"autofocus",void 0);f([s.attr],Ar.prototype,"placeholder",void 0);f([s.attr],Ar.prototype,"type",void 0);f([s.attr],Ar.prototype,"list",void 0);f([(0,s.attr)({converter:s.nullableNumberConverter})],Ar.prototype,"maxlength",void 0);f([(0,s.attr)({converter:s.nullableNumberConverter})],Ar.prototype,"minlength",void 0);f([s.attr],Ar.prototype,"pattern",void 0);f([(0,s.attr)({converter:s.nullableNumberConverter})],Ar.prototype,"size",void 0);f([(0,s.attr)({mode:"boolean"})],Ar.prototype,"spellcheck",void 0);f([s.observable],Ar.prototype,"defaultSlottedNodes",void 0);class Fr{}Pe(Fr,je);Pe(Ar,o,Fr);class Lr extends Fe{}class Mr extends(Ws(Lr)){constructor(){super(...arguments);this.proxy=document.createElement("input")}}class Pr extends Mr{constructor(){super(...arguments);this.hideStep=false;this.step=1;this.isUserInput=false}maxChanged(e,t){var i;this.max=Math.max(t,(i=this.min)!==null&&i!==void 0?i:t);const s=Math.min(this.min,this.max);if(this.min!==undefined&&this.min!==s){this.min=s}this.value=this.getValidValue(this.value)}minChanged(e,t){var i;this.min=Math.min(t,(i=this.max)!==null&&i!==void 0?i:t);const s=Math.max(this.min,this.max);if(this.max!==undefined&&this.max!==s){this.max=s}this.value=this.getValidValue(this.value)}get valueAsNumber(){return parseFloat(super.value)}set valueAsNumber(e){this.value=e.toString()}valueChanged(e,t){this.value=this.getValidValue(t);if(t!==this.value){return}if(this.control&&!this.isUserInput){this.control.value=this.value}super.valueChanged(e,this.value);if(e!==undefined&&!this.isUserInput){this.$emit("input");this.$emit("change")}this.isUserInput=false}validate(){super.validate(this.control)}getValidValue(e){var t,i;let s=parseFloat(parseFloat(e).toPrecision(12));if(isNaN(s)){s=""}else{s=Math.min(s,(t=this.max)!==null&&t!==void 0?t:s);s=Math.max(s,(i=this.min)!==null&&i!==void 0?i:s).toString()}return s}stepUp(){const e=parseFloat(this.value);const t=!isNaN(e)?e+this.step:this.min>0?this.min:this.max<0?this.max:!this.min?this.step:0;this.value=t.toString()}stepDown(){const e=parseFloat(this.value);const t=!isNaN(e)?e-this.step:this.min>0?this.min:this.max<0?this.max:!this.min?0-this.step:0;this.value=t.toString()}connectedCallback(){super.connectedCallback();this.proxy.setAttribute("type","number");this.validate();this.control.value=this.value;if(this.autofocus){s.DOM.queueUpdate((()=>{this.focus()}))}}select(){this.control.select();this.$emit("select")}handleTextInput(){this.control.value=this.control.value.replace(/[^0-9\-+e.]/g,"");this.isUserInput=true;this.value=this.control.value}handleChange(){this.$emit("change")}handleKeyDown(e){const t=e.key;switch(t){case ze.I5:this.stepUp();return false;case ze.HX:this.stepDown();return false}return true}handleBlur(){this.control.value=this.value}}f([(0,s.attr)({attribute:"readonly",mode:"boolean"})],Pr.prototype,"readOnly",void 0);f([(0,s.attr)({mode:"boolean"})],Pr.prototype,"autofocus",void 0);f([(0,s.attr)({attribute:"hide-step",mode:"boolean"})],Pr.prototype,"hideStep",void 0);f([s.attr],Pr.prototype,"placeholder",void 0);f([s.attr],Pr.prototype,"list",void 0);f([(0,s.attr)({converter:s.nullableNumberConverter})],Pr.prototype,"maxlength",void 0);f([(0,s.attr)({converter:s.nullableNumberConverter})],Pr.prototype,"minlength",void 0);f([(0,s.attr)({converter:s.nullableNumberConverter})],Pr.prototype,"size",void 0);f([(0,s.attr)({converter:s.nullableNumberConverter})],Pr.prototype,"step",void 0);f([(0,s.attr)({converter:s.nullableNumberConverter})],Pr.prototype,"max",void 0);f([(0,s.attr)({converter:s.nullableNumberConverter})],Pr.prototype,"min",void 0);f([s.observable],Pr.prototype,"defaultSlottedNodes",void 0);Pe(Pr,o,Fr);const Hr=44;const Vr=(e,t)=>(0,s.html)`
    <template
        role="progressbar"
        aria-valuenow="${e=>e.value}"
        aria-valuemin="${e=>e.min}"
        aria-valuemax="${e=>e.max}"
        class="${e=>e.paused?"paused":""}"
    >
        ${(0,s.when)((e=>typeof e.value==="number"),(0,s.html)`
                <svg
                    class="progress"
                    part="progress"
                    viewBox="0 0 16 16"
                    slot="determinate"
                >
                    <circle
                        class="background"
                        part="background"
                        cx="8px"
                        cy="8px"
                        r="7px"
                    ></circle>
                    <circle
                        class="determinate"
                        part="determinate"
                        style="stroke-dasharray: ${e=>Hr*e.percentComplete/100}px ${Hr}px"
                        cx="8px"
                        cy="8px"
                        r="7px"
                    ></circle>
                </svg>
            `,(0,s.html)`
                <slot name="indeterminate" slot="indeterminate">
                    ${t.indeterminateIndicator||""}
                </slot>
            `)}
    </template>
`;class zr extends Fe{constructor(){super(...arguments);this.percentComplete=0}valueChanged(){if(this.$fastController.isConnected){this.updatePercentComplete()}}minChanged(){if(this.$fastController.isConnected){this.updatePercentComplete()}}maxChanged(){if(this.$fastController.isConnected){this.updatePercentComplete()}}connectedCallback(){super.connectedCallback();this.updatePercentComplete()}updatePercentComplete(){const e=typeof this.min==="number"?this.min:0;const t=typeof this.max==="number"?this.max:100;const i=typeof this.value==="number"?this.value:0;const s=t-e;this.percentComplete=s===0?0:Math.fround((i-e)/s*100)}}f([(0,s.attr)({converter:s.nullableNumberConverter})],zr.prototype,"value",void 0);f([(0,s.attr)({converter:s.nullableNumberConverter})],zr.prototype,"min",void 0);f([(0,s.attr)({converter:s.nullableNumberConverter})],zr.prototype,"max",void 0);f([(0,s.attr)({mode:"boolean"})],zr.prototype,"paused",void 0);f([s.observable],zr.prototype,"percentComplete",void 0);const Nr=(e,t)=>(0,s.html)`
    <template
        role="progressbar"
        aria-valuenow="${e=>e.value}"
        aria-valuemin="${e=>e.min}"
        aria-valuemax="${e=>e.max}"
        class="${e=>e.paused?"paused":""}"
    >
        ${(0,s.when)((e=>typeof e.value==="number"),(0,s.html)`
                <div class="progress" part="progress" slot="determinate">
                    <div
                        class="determinate"
                        part="determinate"
                        style="width: ${e=>e.percentComplete}%"
                    ></div>
                </div>
            `,(0,s.html)`
                <div class="progress" part="progress" slot="indeterminate">
                    <slot class="indeterminate" name="indeterminate">
                        ${t.indeterminateIndicator1||""}
                        ${t.indeterminateIndicator2||""}
                    </slot>
                </div>
            `)}
    </template>
`;const Br=(e,t)=>(0,s.html)`
    <template
        role="radiogroup"
        aria-disabled="${e=>e.disabled}"
        aria-readonly="${e=>e.readOnly}"
        @click="${(e,t)=>e.clickHandler(t.event)}"
        @keydown="${(e,t)=>e.keydownHandler(t.event)}"
        @focusout="${(e,t)=>e.focusOutHandler(t.event)}"
    >
        <slot name="label"></slot>
        <div
            class="positioning-region ${e=>e.orientation===Yn.t.horizontal?"horizontal":"vertical"}"
            part="positioning-region"
        >
            <slot
                ${(0,s.slotted)({property:"slottedRadioButtons",filter:(0,s.elements)("[role=radio]")})}
            ></slot>
        </div>
    </template>
`;class qr extends Fe{constructor(){super(...arguments);this.orientation=Yn.t.horizontal;this.radioChangeHandler=e=>{const t=e.target;if(t.checked){this.slottedRadioButtons.forEach((e=>{if(e!==t){e.checked=false;if(!this.isInsideFoundationToolbar){e.setAttribute("tabindex","-1")}}}));this.selectedRadio=t;this.value=t.value;t.setAttribute("tabindex","0");this.focusedRadio=t}e.stopPropagation()};this.moveToRadioByIndex=(e,t)=>{const i=e[t];if(!this.isInsideToolbar){i.setAttribute("tabindex","0");if(i.readOnly){this.slottedRadioButtons.forEach((e=>{if(e!==i){e.setAttribute("tabindex","-1")}}))}else{i.checked=true;this.selectedRadio=i}}this.focusedRadio=i;i.focus()};this.moveRightOffGroup=()=>{var e;(e=this.nextElementSibling)===null||e===void 0?void 0:e.focus()};this.moveLeftOffGroup=()=>{var e;(e=this.previousElementSibling)===null||e===void 0?void 0:e.focus()};this.focusOutHandler=e=>{const t=this.slottedRadioButtons;const i=e.target;const s=i!==null?t.indexOf(i):0;const o=this.focusedRadio?t.indexOf(this.focusedRadio):-1;if(o===0&&s===o||o===t.length-1&&o===s){if(!this.selectedRadio){this.focusedRadio=t[0];this.focusedRadio.setAttribute("tabindex","0");t.forEach((e=>{if(e!==this.focusedRadio){e.setAttribute("tabindex","-1")}}))}else{this.focusedRadio=this.selectedRadio;if(!this.isInsideFoundationToolbar){this.selectedRadio.setAttribute("tabindex","0");t.forEach((e=>{if(e!==this.selectedRadio){e.setAttribute("tabindex","-1")}}))}}}return true};this.clickHandler=e=>{const t=e.target;if(t){const e=this.slottedRadioButtons;if(t.checked||e.indexOf(t)===0){t.setAttribute("tabindex","0");this.selectedRadio=t}else{t.setAttribute("tabindex","-1");this.selectedRadio=null}this.focusedRadio=t}e.preventDefault()};this.shouldMoveOffGroupToTheRight=(e,t,i)=>e===t.length&&this.isInsideToolbar&&i===ze.bb;this.shouldMoveOffGroupToTheLeft=(e,t)=>{const i=this.focusedRadio?e.indexOf(this.focusedRadio)-1:0;return i<0&&this.isInsideToolbar&&t===ze.kT};this.checkFocusedRadio=()=>{if(this.focusedRadio!==null&&!this.focusedRadio.readOnly&&!this.focusedRadio.checked){this.focusedRadio.checked=true;this.focusedRadio.setAttribute("tabindex","0");this.focusedRadio.focus();this.selectedRadio=this.focusedRadio}};this.moveRight=e=>{const t=this.slottedRadioButtons;let i=0;i=this.focusedRadio?t.indexOf(this.focusedRadio)+1:1;if(this.shouldMoveOffGroupToTheRight(i,t,e.key)){this.moveRightOffGroup();return}else if(i===t.length){i=0}while(i<t.length&&t.length>1){if(!t[i].disabled){this.moveToRadioByIndex(t,i);break}else if(this.focusedRadio&&i===t.indexOf(this.focusedRadio)){break}else if(i+1>=t.length){if(this.isInsideToolbar){break}else{i=0}}else{i+=1}}};this.moveLeft=e=>{const t=this.slottedRadioButtons;let i=0;i=this.focusedRadio?t.indexOf(this.focusedRadio)-1:0;i=i<0?t.length-1:i;if(this.shouldMoveOffGroupToTheLeft(t,e.key)){this.moveLeftOffGroup();return}while(i>=0&&t.length>1){if(!t[i].disabled){this.moveToRadioByIndex(t,i);break}else if(this.focusedRadio&&i===t.indexOf(this.focusedRadio)){break}else if(i-1<0){i=t.length-1}else{i-=1}}};this.keydownHandler=e=>{const t=e.key;if(t in ze.Is&&this.isInsideFoundationToolbar){return true}switch(t){case ze.Mm:{this.checkFocusedRadio();break}case ze.bb:case ze.HX:{if(this.direction===Ge.O.ltr){this.moveRight(e)}else{this.moveLeft(e)}break}case ze.kT:case ze.I5:{if(this.direction===Ge.O.ltr){this.moveLeft(e)}else{this.moveRight(e)}break}default:{return true}}}}readOnlyChanged(){if(this.slottedRadioButtons!==undefined){this.slottedRadioButtons.forEach((e=>{if(this.readOnly){e.readOnly=true}else{e.readOnly=false}}))}}disabledChanged(){if(this.slottedRadioButtons!==undefined){this.slottedRadioButtons.forEach((e=>{if(this.disabled){e.disabled=true}else{e.disabled=false}}))}}nameChanged(){if(this.slottedRadioButtons){this.slottedRadioButtons.forEach((e=>{e.setAttribute("name",this.name)}))}}valueChanged(){if(this.slottedRadioButtons){this.slottedRadioButtons.forEach((e=>{if(e.value===this.value){e.checked=true;this.selectedRadio=e}}))}this.$emit("change")}slottedRadioButtonsChanged(e,t){if(this.slottedRadioButtons&&this.slottedRadioButtons.length>0){this.setupRadioButtons()}}get parentToolbar(){return this.closest('[role="toolbar"]')}get isInsideToolbar(){var e;return(e=this.parentToolbar)!==null&&e!==void 0?e:false}get isInsideFoundationToolbar(){var e;return!!((e=this.parentToolbar)===null||e===void 0?void 0:e["$fastController"])}connectedCallback(){super.connectedCallback();this.direction=Is(this);this.setupRadioButtons()}disconnectedCallback(){this.slottedRadioButtons.forEach((e=>{e.removeEventListener("change",this.radioChangeHandler)}))}setupRadioButtons(){const e=this.slottedRadioButtons.filter((e=>e.hasAttribute("checked")));const t=e?e.length:0;if(t>1){const i=e[t-1];i.checked=true}let i=false;this.slottedRadioButtons.forEach((e=>{if(this.name!==undefined){e.setAttribute("name",this.name)}if(this.disabled){e.disabled=true}if(this.readOnly){e.readOnly=true}if(this.value&&this.value===e.value){this.selectedRadio=e;this.focusedRadio=e;e.checked=true;e.setAttribute("tabindex","0");i=true}else{if(!this.isInsideFoundationToolbar){e.setAttribute("tabindex","-1")}e.checked=false}e.addEventListener("change",this.radioChangeHandler)}));if(this.value===undefined&&this.slottedRadioButtons.length>0){const e=this.slottedRadioButtons.filter((e=>e.hasAttribute("checked")));const t=e!==null?e.length:0;if(t>0&&!i){const i=e[t-1];i.checked=true;this.focusedRadio=i;i.setAttribute("tabindex","0")}else{this.slottedRadioButtons[0].setAttribute("tabindex","0");this.focusedRadio=this.slottedRadioButtons[0]}}}}f([(0,s.attr)({attribute:"readonly",mode:"boolean"})],qr.prototype,"readOnly",void 0);f([(0,s.attr)({attribute:"disabled",mode:"boolean"})],qr.prototype,"disabled",void 0);f([s.attr],qr.prototype,"name",void 0);f([s.attr],qr.prototype,"value",void 0);f([s.attr],qr.prototype,"orientation",void 0);f([s.observable],qr.prototype,"childItems",void 0);f([s.observable],qr.prototype,"slottedRadioButtons",void 0);const Ur=(e,t)=>(0,s.html)`
    <template
        role="radio"
        class="${e=>e.checked?"checked":""} ${e=>e.readOnly?"readonly":""}"
        aria-checked="${e=>e.checked}"
        aria-required="${e=>e.required}"
        aria-disabled="${e=>e.disabled}"
        aria-readonly="${e=>e.readOnly}"
        @keypress="${(e,t)=>e.keypressHandler(t.event)}"
        @click="${(e,t)=>e.clickHandler(t.event)}"
    >
        <div part="control" class="control">
            <slot name="checked-indicator">
                ${t.checkedIndicator||""}
            </slot>
        </div>
        <label
            part="label"
            class="${e=>e.defaultSlottedNodes&&e.defaultSlottedNodes.length?"label":"label label__hidden"}"
        >
            <slot ${(0,s.slotted)("defaultSlottedNodes")}></slot>
        </label>
    </template>
`;class jr extends Fe{}class _r extends(Gs(jr)){constructor(){super(...arguments);this.proxy=document.createElement("input")}}class Kr extends _r{constructor(){super();this.initialValue="on";this.keypressHandler=e=>{switch(e.key){case ze.gG:if(!this.checked&&!this.readOnly){this.checked=true}return}return true};this.proxy.setAttribute("type","radio")}readOnlyChanged(){if(this.proxy instanceof HTMLInputElement){this.proxy.readOnly=this.readOnly}}defaultCheckedChanged(){var e;if(this.$fastController.isConnected&&!this.dirtyChecked){if(!this.isInsideRadioGroup()){this.checked=(e=this.defaultChecked)!==null&&e!==void 0?e:false;this.dirtyChecked=false}}}connectedCallback(){var e,t;super.connectedCallback();this.validate();if(((e=this.parentElement)===null||e===void 0?void 0:e.getAttribute("role"))!=="radiogroup"&&this.getAttribute("tabindex")===null){if(!this.disabled){this.setAttribute("tabindex","0")}}if(this.checkedAttribute){if(!this.dirtyChecked){if(!this.isInsideRadioGroup()){this.checked=(t=this.defaultChecked)!==null&&t!==void 0?t:false;this.dirtyChecked=false}}}}isInsideRadioGroup(){const e=this.closest("[role=radiogroup]");return e!==null}clickHandler(e){if(!this.disabled&&!this.readOnly&&!this.checked){this.checked=true}}}f([(0,s.attr)({attribute:"readonly",mode:"boolean"})],Kr.prototype,"readOnly",void 0);f([s.observable],Kr.prototype,"name",void 0);f([s.observable],Kr.prototype,"defaultSlottedNodes",void 0);class Wr extends Fe{constructor(){super(...arguments);this.framesPerSecond=60;this.updatingItems=false;this.speed=600;this.easing="ease-in-out";this.flippersHiddenFromAT=false;this.scrolling=false;this.resizeDetector=null}get frameTime(){return 1e3/this.framesPerSecond}scrollingChanged(e,t){if(this.scrollContainer){const e=this.scrolling==true?"scrollstart":"scrollend";this.$emit(e,this.scrollContainer.scrollLeft)}}get isRtl(){return this.scrollItems.length>1&&this.scrollItems[0].offsetLeft>this.scrollItems[1].offsetLeft}connectedCallback(){super.connectedCallback();this.initializeResizeDetector()}disconnectedCallback(){this.disconnectResizeDetector();super.disconnectedCallback()}scrollItemsChanged(e,t){if(t&&!this.updatingItems){s.DOM.queueUpdate((()=>this.setStops()))}}disconnectResizeDetector(){if(this.resizeDetector){this.resizeDetector.disconnect();this.resizeDetector=null}}initializeResizeDetector(){this.disconnectResizeDetector();this.resizeDetector=new window.ResizeObserver(this.resized.bind(this));this.resizeDetector.observe(this)}updateScrollStops(){this.updatingItems=true;const e=this.scrollItems.reduce(((e,t)=>{if(t instanceof HTMLSlotElement){return e.concat(t.assignedElements())}e.push(t);return e}),[]);this.scrollItems=e;this.updatingItems=false}setStops(){this.updateScrollStops();const{scrollContainer:e}=this;const{scrollLeft:t}=e;const{width:i,left:s}=e.getBoundingClientRect();this.width=i;let o=0;let n=this.scrollItems.map(((e,i)=>{const{left:n,width:r}=e.getBoundingClientRect();const a=Math.round(n+t-s);const l=Math.round(a+r);if(this.isRtl){return-l}o=l;return i===0?0:a})).concat(o);n=this.fixScrollMisalign(n);n.sort(((e,t)=>Math.abs(e)-Math.abs(t)));this.scrollStops=n;this.setFlippers()}validateStops(e=true){const t=()=>!!this.scrollStops.find((e=>e>0));if(!t()&&e){this.setStops()}return t()}fixScrollMisalign(e){if(this.isRtl&&e.some((e=>e>0))){e.sort(((e,t)=>t-e));const t=e[0];e=e.map((e=>e-t))}return e}setFlippers(){var e,t;const i=this.scrollContainer.scrollLeft;(e=this.previousFlipperContainer)===null||e===void 0?void 0:e.classList.toggle("disabled",i===0);if(this.scrollStops){const e=Math.abs(this.scrollStops[this.scrollStops.length-1]);(t=this.nextFlipperContainer)===null||t===void 0?void 0:t.classList.toggle("disabled",this.validateStops(false)&&Math.abs(i)+this.width>=e)}}scrollInView(e,t=0,i){var s;if(typeof e!=="number"&&e){e=this.scrollItems.findIndex((t=>t===e||t.contains(e)))}if(e!==undefined){i=i!==null&&i!==void 0?i:t;const{scrollContainer:o,scrollStops:n,scrollItems:r}=this;const{scrollLeft:a}=this.scrollContainer;const{width:l}=o.getBoundingClientRect();const d=n[e];const{width:h}=r[e].getBoundingClientRect();const c=d+h;const u=a+t>d;if(u||a+l-i<c){const e=[...n].sort(((e,t)=>u?t-e:e-t));const o=(s=e.find((e=>u?e+t<d:e+l-(i!==null&&i!==void 0?i:0)>c)))!==null&&s!==void 0?s:0;this.scrollToPosition(o)}}}keyupHandler(e){const t=e.key;switch(t){case"ArrowLeft":this.scrollToPrevious();break;case"ArrowRight":this.scrollToNext();break}}scrollToPrevious(){this.validateStops();const e=this.scrollContainer.scrollLeft;const t=this.scrollStops.findIndex(((t,i)=>t>=e&&(this.isRtl||i===this.scrollStops.length-1||this.scrollStops[i+1]>e)));const i=Math.abs(this.scrollStops[t+1]);let s=this.scrollStops.findIndex((e=>Math.abs(e)+this.width>i));if(s>=t||s===-1){s=t>0?t-1:0}this.scrollToPosition(this.scrollStops[s],e)}scrollToNext(){this.validateStops();const e=this.scrollContainer.scrollLeft;const t=this.scrollStops.findIndex((t=>Math.abs(t)>=Math.abs(e)));const i=this.scrollStops.findIndex((t=>Math.abs(e)+this.width<=Math.abs(t)));let s=t;if(i>t+2){s=i-2}else if(t<this.scrollStops.length-2){s=t+1}this.scrollToPosition(this.scrollStops[s],e)}scrollToPosition(e,t=this.scrollContainer.scrollLeft){var i;if(this.scrolling){return}this.scrolling=true;const s=(i=this.duration)!==null&&i!==void 0?i:`${Math.abs(e-t)/this.speed}s`;this.content.style.setProperty("transition-duration",s);const o=parseFloat(getComputedStyle(this.content).getPropertyValue("transition-duration"));const n=t=>{if(t&&t.target!==t.currentTarget){return}this.content.style.setProperty("transition-duration","0s");this.content.style.removeProperty("transform");this.scrollContainer.style.setProperty("scroll-behavior","auto");this.scrollContainer.scrollLeft=e;this.setFlippers();this.content.removeEventListener("transitionend",n);this.scrolling=false};if(o===0){n();return}this.content.addEventListener("transitionend",n);const r=this.scrollContainer.scrollWidth-this.scrollContainer.clientWidth;let a=this.scrollContainer.scrollLeft-Math.min(e,r);if(this.isRtl){a=this.scrollContainer.scrollLeft+Math.min(Math.abs(e),r)}this.content.style.setProperty("transition-property","transform");this.content.style.setProperty("transition-timing-function",this.easing);this.content.style.setProperty("transform",`translateX(${a}px)`)}resized(){if(this.resizeTimeout){this.resizeTimeout=clearTimeout(this.resizeTimeout)}this.resizeTimeout=setTimeout((()=>{this.width=this.scrollContainer.offsetWidth;this.setFlippers()}),this.frameTime)}scrolled(){if(this.scrollTimeout){this.scrollTimeout=clearTimeout(this.scrollTimeout)}this.scrollTimeout=setTimeout((()=>{this.setFlippers()}),this.frameTime)}}f([(0,s.attr)({converter:s.nullableNumberConverter})],Wr.prototype,"speed",void 0);f([s.attr],Wr.prototype,"duration",void 0);f([s.attr],Wr.prototype,"easing",void 0);f([(0,s.attr)({attribute:"flippers-hidden-from-at",converter:s.booleanConverter})],Wr.prototype,"flippersHiddenFromAT",void 0);f([s.observable],Wr.prototype,"scrolling",void 0);f([s.observable],Wr.prototype,"scrollItems",void 0);f([(0,s.attr)({attribute:"view"})],Wr.prototype,"view",void 0);const Gr=(e,t)=>{var i,o;return(0,s.html)`
    <template
        class="horizontal-scroll"
        @keyup="${(e,t)=>e.keyupHandler(t.event)}"
    >
        ${r(e,t)}
        <div class="scroll-area" part="scroll-area">
            <div
                class="scroll-view"
                part="scroll-view"
                @scroll="${e=>e.scrolled()}"
                ${(0,s.ref)("scrollContainer")}
            >
                <div class="content-container" part="content-container" ${(0,s.ref)("content")}>
                    <slot
                        ${(0,s.slotted)({property:"scrollItems",filter:(0,s.elements)()})}
                    ></slot>
                </div>
            </div>
            ${(0,s.when)((e=>e.view!=="mobile"),(0,s.html)`
                    <div
                        class="scroll scroll-prev"
                        part="scroll-prev"
                        ${(0,s.ref)("previousFlipperContainer")}
                    >
                        <div class="scroll-action" part="scroll-action-previous">
                            <slot name="previous-flipper">
                                ${t.previousFlipper instanceof Function?t.previousFlipper(e,t):(i=t.previousFlipper)!==null&&i!==void 0?i:""}
                            </slot>
                        </div>
                    </div>
                    <div
                        class="scroll scroll-next"
                        part="scroll-next"
                        ${(0,s.ref)("nextFlipperContainer")}
                    >
                        <div class="scroll-action" part="scroll-action-next">
                            <slot name="next-flipper">
                                ${t.nextFlipper instanceof Function?t.nextFlipper(e,t):(o=t.nextFlipper)!==null&&o!==void 0?o:""}
                            </slot>
                        </div>
                    </div>
                `)}
        </div>
        ${n(e,t)}
    </template>
`};function Xr(e,t,i){return e.nodeType!==Node.TEXT_NODE?true:typeof e.nodeValue==="string"&&!!e.nodeValue.trim().length}const Yr=(e,t)=>(0,s.html)`
    <template
        class="
            ${e=>e.readOnly?"readonly":""}
        "
    >
        <label
            part="label"
            for="control"
            class="${e=>e.defaultSlottedNodes&&e.defaultSlottedNodes.length?"label":"label label__hidden"}"
        >
            <slot
                ${(0,s.slotted)({property:"defaultSlottedNodes",filter:Xr})}
            ></slot>
        </label>
        <div class="root" part="root" ${(0,s.ref)("root")}>
            ${r(e,t)}
            <div class="input-wrapper" part="input-wrapper">
                <input
                    class="control"
                    part="control"
                    id="control"
                    @input="${e=>e.handleTextInput()}"
                    @change="${e=>e.handleChange()}"
                    ?autofocus="${e=>e.autofocus}"
                    ?disabled="${e=>e.disabled}"
                    list="${e=>e.list}"
                    maxlength="${e=>e.maxlength}"
                    minlength="${e=>e.minlength}"
                    pattern="${e=>e.pattern}"
                    placeholder="${e=>e.placeholder}"
                    ?readonly="${e=>e.readOnly}"
                    ?required="${e=>e.required}"
                    size="${e=>e.size}"
                    ?spellcheck="${e=>e.spellcheck}"
                    :value="${e=>e.value}"
                    type="search"
                    aria-atomic="${e=>e.ariaAtomic}"
                    aria-busy="${e=>e.ariaBusy}"
                    aria-controls="${e=>e.ariaControls}"
                    aria-current="${e=>e.ariaCurrent}"
                    aria-describedby="${e=>e.ariaDescribedby}"
                    aria-details="${e=>e.ariaDetails}"
                    aria-disabled="${e=>e.ariaDisabled}"
                    aria-errormessage="${e=>e.ariaErrormessage}"
                    aria-flowto="${e=>e.ariaFlowto}"
                    aria-haspopup="${e=>e.ariaHaspopup}"
                    aria-hidden="${e=>e.ariaHidden}"
                    aria-invalid="${e=>e.ariaInvalid}"
                    aria-keyshortcuts="${e=>e.ariaKeyshortcuts}"
                    aria-label="${e=>e.ariaLabel}"
                    aria-labelledby="${e=>e.ariaLabelledby}"
                    aria-live="${e=>e.ariaLive}"
                    aria-owns="${e=>e.ariaOwns}"
                    aria-relevant="${e=>e.ariaRelevant}"
                    aria-roledescription="${e=>e.ariaRoledescription}"
                    ${(0,s.ref)("control")}
                />
                <slot name="close-button">
                    <button
                        class="clear-button ${e=>e.value?"":"clear-button__hidden"}"
                        part="clear-button"
                        tabindex="-1"
                        @click=${e=>e.handleClearInput()}
                    >
                        <slot name="close-glyph">
                            <svg
                                width="9"
                                height="9"
                                viewBox="0 0 9 9"
                                xmlns="http://www.w3.org/2000/svg"
                            >
                                <path
                                    d="M0.146447 0.146447C0.338683 -0.0478972 0.645911 -0.0270359 0.853553 0.146447L4.5 3.793L8.14645 0.146447C8.34171 -0.0488155 8.65829 -0.0488155 8.85355 0.146447C9.04882 0.341709 9.04882 0.658291 8.85355 0.853553L5.207 4.5L8.85355 8.14645C9.05934 8.35223 9.03129 8.67582 8.85355 8.85355C8.67582 9.03129 8.35409 9.02703 8.14645 8.85355L4.5 5.207L0.853553 8.85355C0.658291 9.04882 0.341709 9.04882 0.146447 8.85355C-0.0488155 8.65829 -0.0488155 8.34171 0.146447 8.14645L3.793 4.5L0.146447 0.853553C-0.0268697 0.680237 -0.0457894 0.34079 0.146447 0.146447Z"
                                />
                            </svg>
                        </slot>
                    </button>
                </slot>
            </div>
            ${n(e,t)}
        </div>
    </template>
`;class Qr extends Fe{}class Zr extends(Ws(Qr)){constructor(){super(...arguments);this.proxy=document.createElement("input")}}class Jr extends Zr{readOnlyChanged(){if(this.proxy instanceof HTMLInputElement){this.proxy.readOnly=this.readOnly;this.validate()}}autofocusChanged(){if(this.proxy instanceof HTMLInputElement){this.proxy.autofocus=this.autofocus;this.validate()}}placeholderChanged(){if(this.proxy instanceof HTMLInputElement){this.proxy.placeholder=this.placeholder}}listChanged(){if(this.proxy instanceof HTMLInputElement){this.proxy.setAttribute("list",this.list);this.validate()}}maxlengthChanged(){if(this.proxy instanceof HTMLInputElement){this.proxy.maxLength=this.maxlength;this.validate()}}minlengthChanged(){if(this.proxy instanceof HTMLInputElement){this.proxy.minLength=this.minlength;this.validate()}}patternChanged(){if(this.proxy instanceof HTMLInputElement){this.proxy.pattern=this.pattern;this.validate()}}sizeChanged(){if(this.proxy instanceof HTMLInputElement){this.proxy.size=this.size}}spellcheckChanged(){if(this.proxy instanceof HTMLInputElement){this.proxy.spellcheck=this.spellcheck}}connectedCallback(){super.connectedCallback();this.validate();if(this.autofocus){s.DOM.queueUpdate((()=>{this.focus()}))}}validate(){super.validate(this.control)}handleTextInput(){this.value=this.control.value}handleClearInput(){this.value="";this.control.focus();this.handleChange()}handleChange(){this.$emit("change")}}f([(0,s.attr)({attribute:"readonly",mode:"boolean"})],Jr.prototype,"readOnly",void 0);f([(0,s.attr)({mode:"boolean"})],Jr.prototype,"autofocus",void 0);f([s.attr],Jr.prototype,"placeholder",void 0);f([s.attr],Jr.prototype,"list",void 0);f([(0,s.attr)({converter:s.nullableNumberConverter})],Jr.prototype,"maxlength",void 0);f([(0,s.attr)({converter:s.nullableNumberConverter})],Jr.prototype,"minlength",void 0);f([s.attr],Jr.prototype,"pattern",void 0);f([(0,s.attr)({converter:s.nullableNumberConverter})],Jr.prototype,"size",void 0);f([(0,s.attr)({mode:"boolean"})],Jr.prototype,"spellcheck",void 0);f([s.observable],Jr.prototype,"defaultSlottedNodes",void 0);class ea{}Pe(ea,je);Pe(Jr,o,ea);class ta extends sr{}class ia extends(Ws(ta)){constructor(){super(...arguments);this.proxy=document.createElement("select")}}class sa extends ia{constructor(){super(...arguments);this.open=false;this.forcedPosition=false;this.listboxId=Io("listbox-");this.maxHeight=0}openChanged(e,t){if(!this.collapsible){return}if(this.open){this.ariaControls=this.listboxId;this.ariaExpanded="true";this.setPositioning();this.focusAndScrollOptionIntoView();this.indexWhenOpened=this.selectedIndex;s.DOM.queueUpdate((()=>this.focus()));return}this.ariaControls="";this.ariaExpanded="false"}get collapsible(){return!(this.multiple||typeof this.size==="number")}get value(){s.Observable.track(this,"value");return this._value}set value(e){var t,i,o,n,r,a,l;const d=`${this._value}`;if((t=this._options)===null||t===void 0?void 0:t.length){const t=this._options.findIndex((t=>t.value===e));const s=(o=(i=this._options[this.selectedIndex])===null||i===void 0?void 0:i.value)!==null&&o!==void 0?o:null;const d=(r=(n=this._options[t])===null||n===void 0?void 0:n.value)!==null&&r!==void 0?r:null;if(t===-1||s!==d){e="";this.selectedIndex=t}e=(l=(a=this.firstSelectedOption)===null||a===void 0?void 0:a.value)!==null&&l!==void 0?l:e}if(d!==e){this._value=e;super.valueChanged(d,e);s.Observable.notify(this,"value");this.updateDisplayValue()}}updateValue(e){var t,i;if(this.$fastController.isConnected){this.value=(i=(t=this.firstSelectedOption)===null||t===void 0?void 0:t.value)!==null&&i!==void 0?i:""}if(e){this.$emit("input");this.$emit("change",this,{bubbles:true,composed:undefined})}}selectedIndexChanged(e,t){super.selectedIndexChanged(e,t);this.updateValue()}positionChanged(e,t){this.positionAttribute=t;this.setPositioning()}setPositioning(){const e=this.getBoundingClientRect();const t=window.innerHeight;const i=t-e.bottom;this.position=this.forcedPosition?this.positionAttribute:e.top>i?Go.above:Go.below;this.positionAttribute=this.forcedPosition?this.positionAttribute:this.position;this.maxHeight=this.position===Go.above?~~e.top:~~i}get displayValue(){var e,t;s.Observable.track(this,"displayValue");return(t=(e=this.firstSelectedOption)===null||e===void 0?void 0:e.text)!==null&&t!==void 0?t:""}disabledChanged(e,t){if(super.disabledChanged){super.disabledChanged(e,t)}this.ariaDisabled=this.disabled?"true":"false"}formResetCallback(){this.setProxyOptions();super.setDefaultSelectedOption();if(this.selectedIndex===-1){this.selectedIndex=0}}clickHandler(e){if(this.disabled){return}if(this.open){const t=e.target.closest(`option,[role=option]`);if(t&&t.disabled){return}}super.clickHandler(e);this.open=this.collapsible&&!this.open;if(!this.open&&this.indexWhenOpened!==this.selectedIndex){this.updateValue(true)}return true}focusoutHandler(e){var t;super.focusoutHandler(e);if(!this.open){return true}const i=e.relatedTarget;if(this.isSameNode(i)){this.focus();return}if(!((t=this.options)===null||t===void 0?void 0:t.includes(i))){this.open=false;if(this.indexWhenOpened!==this.selectedIndex){this.updateValue(true)}}}handleChange(e,t){super.handleChange(e,t);if(t==="value"){this.updateValue()}}slottedOptionsChanged(e,t){this.options.forEach((e=>{const t=s.Observable.getNotifier(e);t.unsubscribe(this,"value")}));super.slottedOptionsChanged(e,t);this.options.forEach((e=>{const t=s.Observable.getNotifier(e);t.subscribe(this,"value")}));this.setProxyOptions();this.updateValue()}mousedownHandler(e){var t;if(e.offsetX>=0&&e.offsetX<=((t=this.listbox)===null||t===void 0?void 0:t.scrollWidth)){return super.mousedownHandler(e)}return this.collapsible}multipleChanged(e,t){super.multipleChanged(e,t);if(this.proxy){this.proxy.multiple=t}}selectedOptionsChanged(e,t){var i;super.selectedOptionsChanged(e,t);(i=this.options)===null||i===void 0?void 0:i.forEach(((e,t)=>{var i;const s=(i=this.proxy)===null||i===void 0?void 0:i.options.item(t);if(s){s.selected=e.selected}}))}setDefaultSelectedOption(){var e;const t=(e=this.options)!==null&&e!==void 0?e:Array.from(this.children).filter(Ko.slottedOptionFilter);const i=t===null||t===void 0?void 0:t.findIndex((e=>e.hasAttribute("selected")||e.selected||e.value===this.value));if(i!==-1){this.selectedIndex=i;return}this.selectedIndex=0}setProxyOptions(){if(this.proxy instanceof HTMLSelectElement&&this.options){this.proxy.options.length=0;this.options.forEach((e=>{const t=e.proxy||(e instanceof HTMLOptionElement?e.cloneNode():null);if(t){this.proxy.options.add(t)}}))}}keydownHandler(e){super.keydownHandler(e);const t=e.key||e.key.charCodeAt(0);switch(t){case ze.gG:{e.preventDefault();if(this.collapsible&&this.typeAheadExpired){this.open=!this.open}break}case ze.Tg:case ze.FM:{e.preventDefault();break}case ze.Mm:{e.preventDefault();this.open=!this.open;break}case ze.F9:{if(this.collapsible&&this.open){e.preventDefault();this.open=false}break}case ze.J9:{if(this.collapsible&&this.open){e.preventDefault();this.open=false}return true}}if(!this.open&&this.indexWhenOpened!==this.selectedIndex){this.updateValue(true);this.indexWhenOpened=this.selectedIndex}return!(t===ze.HX||t===ze.I5)}connectedCallback(){super.connectedCallback();this.forcedPosition=!!this.positionAttribute;this.addEventListener("contentchange",this.updateDisplayValue)}disconnectedCallback(){this.removeEventListener("contentchange",this.updateDisplayValue);super.disconnectedCallback()}sizeChanged(e,t){super.sizeChanged(e,t);if(this.proxy){this.proxy.size=t}}updateDisplayValue(){if(this.collapsible){s.Observable.notify(this,"displayValue")}}}f([(0,s.attr)({attribute:"open",mode:"boolean"})],sa.prototype,"open",void 0);f([s.volatile],sa.prototype,"collapsible",null);f([s.observable],sa.prototype,"control",void 0);f([(0,s.attr)({attribute:"position"})],sa.prototype,"positionAttribute",void 0);f([s.observable],sa.prototype,"position",void 0);f([s.observable],sa.prototype,"maxHeight",void 0);class oa{}f([s.observable],oa.prototype,"ariaControls",void 0);Pe(oa,Wo);Pe(sa,o,oa);const na=(e,t)=>(0,s.html)`
    <template
        class="${e=>[e.collapsible&&"collapsible",e.collapsible&&e.open&&"open",e.disabled&&"disabled",e.collapsible&&e.position].filter(Boolean).join(" ")}"
        aria-activedescendant="${e=>e.ariaActiveDescendant}"
        aria-controls="${e=>e.ariaControls}"
        aria-disabled="${e=>e.ariaDisabled}"
        aria-expanded="${e=>e.ariaExpanded}"
        aria-haspopup="${e=>e.collapsible?"listbox":null}"
        aria-multiselectable="${e=>e.ariaMultiSelectable}"
        ?open="${e=>e.open}"
        role="combobox"
        tabindex="${e=>!e.disabled?"0":null}"
        @click="${(e,t)=>e.clickHandler(t.event)}"
        @focusin="${(e,t)=>e.focusinHandler(t.event)}"
        @focusout="${(e,t)=>e.focusoutHandler(t.event)}"
        @keydown="${(e,t)=>e.keydownHandler(t.event)}"
        @mousedown="${(e,t)=>e.mousedownHandler(t.event)}"
    >
        ${(0,s.when)((e=>e.collapsible),(0,s.html)`
                <div
                    class="control"
                    part="control"
                    ?disabled="${e=>e.disabled}"
                    ${(0,s.ref)("control")}
                >
                    ${r(e,t)}
                    <slot name="button-container">
                        <div class="selected-value" part="selected-value">
                            <slot name="selected-value">${e=>e.displayValue}</slot>
                        </div>
                        <div aria-hidden="true" class="indicator" part="indicator">
                            <slot name="indicator">
                                ${t.indicator||""}
                            </slot>
                        </div>
                    </slot>
                    ${n(e,t)}
                </div>
            `)}
        <div
            class="listbox"
            id="${e=>e.listboxId}"
            part="listbox"
            role="listbox"
            ?disabled="${e=>e.disabled}"
            ?hidden="${e=>e.collapsible?!e.open:false}"
            ${(0,s.ref)("listbox")}
        >
            <slot
                ${(0,s.slotted)({filter:Ko.slottedOptionFilter,flatten:true,property:"slottedOptions"})}
            ></slot>
        </div>
    </template>
`;const ra=(e,t)=>(0,s.html)`
    <template
        class="${e=>e.shape==="circle"?"circle":"rect"}"
        pattern="${e=>e.pattern}"
        ?shimmer="${e=>e.shimmer}"
    >
        ${(0,s.when)((e=>e.shimmer===true),(0,s.html)`
                <span class="shimmer"></span>
            `)}
        <object type="image/svg+xml" data="${e=>e.pattern}" role="presentation">
            <img class="pattern" src="${e=>e.pattern}" />
        </object>
        <slot></slot>
    </template>
`;class aa extends Fe{constructor(){super(...arguments);this.shape="rect"}}f([s.attr],aa.prototype,"fill",void 0);f([s.attr],aa.prototype,"shape",void 0);f([s.attr],aa.prototype,"pattern",void 0);f([(0,s.attr)({mode:"boolean"})],aa.prototype,"shimmer",void 0);const la=(e,t)=>(0,s.html)`
    <template
        aria-disabled="${e=>e.disabled}"
        class="${e=>e.sliderOrientation||Yn.t.horizontal}
            ${e=>e.disabled?"disabled":""}"
    >
        <div ${(0,s.ref)("root")} part="root" class="root" style="${e=>e.positionStyle}">
            <div class="container">
                ${(0,s.when)((e=>!e.hideMark),(0,s.html)`
                        <div class="mark"></div>
                    `)}
                <div class="label">
                    <slot></slot>
                </div>
            </div>
        </div>
    </template>
`;function da(e,t,i,s){let o=(0,Ne.AB)(0,1,(e-t)/(i-t));if(s===Ge.O.rtl){o=1-o}return o}const ha={min:0,max:0,direction:Ge.O.ltr,orientation:Yn.t.horizontal,disabled:false};class ca extends Fe{constructor(){super(...arguments);this.hideMark=false;this.sliderDirection=Ge.O.ltr;this.getSliderConfiguration=()=>{if(!this.isSliderConfig(this.parentNode)){this.sliderDirection=ha.direction||Ge.O.ltr;this.sliderOrientation=ha.orientation||Yn.t.horizontal;this.sliderMaxPosition=ha.max;this.sliderMinPosition=ha.min}else{const e=this.parentNode;const{min:t,max:i,direction:s,orientation:o,disabled:n}=e;if(n!==undefined){this.disabled=n}this.sliderDirection=s||Ge.O.ltr;this.sliderOrientation=o||Yn.t.horizontal;this.sliderMaxPosition=i;this.sliderMinPosition=t}};this.positionAsStyle=()=>{const e=this.sliderDirection?this.sliderDirection:Ge.O.ltr;const t=da(Number(this.position),Number(this.sliderMinPosition),Number(this.sliderMaxPosition));let i=Math.round((1-t)*100);let s=Math.round(t*100);if(Number.isNaN(s)&&Number.isNaN(i)){i=50;s=50}if(this.sliderOrientation===Yn.t.horizontal){return e===Ge.O.rtl?`right: ${s}%; left: ${i}%;`:`left: ${s}%; right: ${i}%;`}else{return`top: ${s}%; bottom: ${i}%;`}}}positionChanged(){this.positionStyle=this.positionAsStyle()}sliderOrientationChanged(){void 0}connectedCallback(){super.connectedCallback();this.getSliderConfiguration();this.positionStyle=this.positionAsStyle();this.notifier=s.Observable.getNotifier(this.parentNode);this.notifier.subscribe(this,"orientation");this.notifier.subscribe(this,"direction");this.notifier.subscribe(this,"max");this.notifier.subscribe(this,"min")}disconnectedCallback(){super.disconnectedCallback();this.notifier.unsubscribe(this,"orientation");this.notifier.unsubscribe(this,"direction");this.notifier.unsubscribe(this,"max");this.notifier.unsubscribe(this,"min")}handleChange(e,t){switch(t){case"direction":this.sliderDirection=e.direction;break;case"orientation":this.sliderOrientation=e.orientation;break;case"max":this.sliderMaxPosition=e.max;break;case"min":this.sliderMinPosition=e.min;break;default:break}this.positionStyle=this.positionAsStyle()}isSliderConfig(e){return e.max!==undefined&&e.min!==undefined}}f([s.observable],ca.prototype,"positionStyle",void 0);f([s.attr],ca.prototype,"position",void 0);f([(0,s.attr)({attribute:"hide-mark",mode:"boolean"})],ca.prototype,"hideMark",void 0);f([(0,s.attr)({attribute:"disabled",mode:"boolean"})],ca.prototype,"disabled",void 0);f([s.observable],ca.prototype,"sliderOrientation",void 0);f([s.observable],ca.prototype,"sliderMinPosition",void 0);f([s.observable],ca.prototype,"sliderMaxPosition",void 0);f([s.observable],ca.prototype,"sliderDirection",void 0);const ua=(e,t)=>(0,s.html)`
    <template
        role="slider"
        class="${e=>e.readOnly?"readonly":""}
        ${e=>e.orientation||Yn.t.horizontal}"
        tabindex="${e=>e.disabled?null:0}"
        aria-valuetext="${e=>e.valueTextFormatter(e.value)}"
        aria-valuenow="${e=>e.value}"
        aria-valuemin="${e=>e.min}"
        aria-valuemax="${e=>e.max}"
        aria-disabled="${e=>e.disabled?true:void 0}"
        aria-readonly="${e=>e.readOnly?true:void 0}"
        aria-orientation="${e=>e.orientation}"
        class="${e=>e.orientation}"
    >
        <div part="positioning-region" class="positioning-region">
            <div ${(0,s.ref)("track")} part="track-container" class="track">
                <slot name="track"></slot>
                <div part="track-start" class="track-start" style="${e=>e.position}">
                    <slot name="track-start"></slot>
                </div>
            </div>
            <slot></slot>
            <div
                ${(0,s.ref)("thumb")}
                part="thumb-container"
                class="thumb-container"
                style="${e=>e.position}"
            >
                <slot name="thumb">${t.thumb||""}</slot>
            </div>
        </div>
    </template>
`;class pa extends Fe{}class fa extends(Ws(pa)){constructor(){super(...arguments);this.proxy=document.createElement("input")}}const ma={singleValue:"single-value"};class va extends fa{constructor(){super(...arguments);this.direction=Ge.O.ltr;this.isDragging=false;this.trackWidth=0;this.trackMinWidth=0;this.trackHeight=0;this.trackLeft=0;this.trackMinHeight=0;this.valueTextFormatter=()=>null;this.min=0;this.max=10;this.step=1;this.orientation=Yn.t.horizontal;this.mode=ma.singleValue;this.keypressHandler=e=>{if(this.readOnly){return}if(e.key===ze.Tg){e.preventDefault();this.value=`${this.min}`}else if(e.key===ze.FM){e.preventDefault();this.value=`${this.max}`}else if(!e.shiftKey){switch(e.key){case ze.bb:case ze.I5:e.preventDefault();this.increment();break;case ze.kT:case ze.HX:e.preventDefault();this.decrement();break}}};this.setupTrackConstraints=()=>{const e=this.track.getBoundingClientRect();this.trackWidth=this.track.clientWidth;this.trackMinWidth=this.track.clientLeft;this.trackHeight=e.bottom;this.trackMinHeight=e.top;this.trackLeft=this.getBoundingClientRect().left;if(this.trackWidth===0){this.trackWidth=1}};this.setupListeners=(e=false)=>{const t=`${e?"remove":"add"}EventListener`;this[t]("keydown",this.keypressHandler);this[t]("mousedown",this.handleMouseDown);this.thumb[t]("mousedown",this.handleThumbMouseDown,{passive:true});this.thumb[t]("touchstart",this.handleThumbMouseDown,{passive:true});if(e){this.handleMouseDown(null);this.handleThumbMouseDown(null)}};this.initialValue="";this.handleThumbMouseDown=e=>{if(e){if(this.readOnly||this.disabled||e.defaultPrevented){return}e.target.focus()}const t=`${e!==null?"add":"remove"}EventListener`;window[t]("mouseup",this.handleWindowMouseUp);window[t]("mousemove",this.handleMouseMove,{passive:true});window[t]("touchmove",this.handleMouseMove,{passive:true});window[t]("touchend",this.handleWindowMouseUp);this.isDragging=e!==null};this.handleMouseMove=e=>{if(this.readOnly||this.disabled||e.defaultPrevented){return}const t=window.TouchEvent&&e instanceof TouchEvent?e.touches[0]:e;const i=this.orientation===Yn.t.horizontal?t.pageX-document.documentElement.scrollLeft-this.trackLeft:t.pageY-document.documentElement.scrollTop;this.value=`${this.calculateNewValue(i)}`};this.calculateNewValue=e=>{const t=da(e,this.orientation===Yn.t.horizontal?this.trackMinWidth:this.trackMinHeight,this.orientation===Yn.t.horizontal?this.trackWidth:this.trackHeight,this.direction);const i=(this.max-this.min)*t+this.min;return this.convertToConstrainedValue(i)};this.handleWindowMouseUp=e=>{this.stopDragging()};this.stopDragging=()=>{this.isDragging=false;this.handleMouseDown(null);this.handleThumbMouseDown(null)};this.handleMouseDown=e=>{const t=`${e!==null?"add":"remove"}EventListener`;if(e===null||!this.disabled&&!this.readOnly){window[t]("mouseup",this.handleWindowMouseUp);window.document[t]("mouseleave",this.handleWindowMouseUp);window[t]("mousemove",this.handleMouseMove);if(e){e.preventDefault();this.setupTrackConstraints();e.target.focus();const t=this.orientation===Yn.t.horizontal?e.pageX-document.documentElement.scrollLeft-this.trackLeft:e.pageY-document.documentElement.scrollTop;this.value=`${this.calculateNewValue(t)}`}}};this.convertToConstrainedValue=e=>{if(isNaN(e)){e=this.min}let t=e-this.min;const i=Math.round(t/this.step);const s=t-i*(this.stepMultiplier*this.step)/this.stepMultiplier;t=s>=Number(this.step)/2?t-s+Number(this.step):t-s;return t+this.min}}readOnlyChanged(){if(this.proxy instanceof HTMLInputElement){this.proxy.readOnly=this.readOnly}}get valueAsNumber(){return parseFloat(super.value)}set valueAsNumber(e){this.value=e.toString()}valueChanged(e,t){super.valueChanged(e,t);if(this.$fastController.isConnected){this.setThumbPositionForOrientation(this.direction)}this.$emit("change")}minChanged(){if(this.proxy instanceof HTMLInputElement){this.proxy.min=`${this.min}`}this.validate()}maxChanged(){if(this.proxy instanceof HTMLInputElement){this.proxy.max=`${this.max}`}this.validate()}stepChanged(){if(this.proxy instanceof HTMLInputElement){this.proxy.step=`${this.step}`}this.updateStepMultiplier();this.validate()}orientationChanged(){if(this.$fastController.isConnected){this.setThumbPositionForOrientation(this.direction)}}connectedCallback(){super.connectedCallback();this.proxy.setAttribute("type","range");this.direction=Is(this);this.updateStepMultiplier();this.setupTrackConstraints();this.setupListeners();this.setupDefaultValue();this.setThumbPositionForOrientation(this.direction)}disconnectedCallback(){this.setupListeners(true)}increment(){const e=this.direction!==Ge.O.rtl&&this.orientation!==Yn.t.vertical?Number(this.value)+Number(this.step):Number(this.value)-Number(this.step);const t=this.convertToConstrainedValue(e);const i=t<Number(this.max)?`${t}`:`${this.max}`;this.value=i}decrement(){const e=this.direction!==Ge.O.rtl&&this.orientation!==Yn.t.vertical?Number(this.value)-Number(this.step):Number(this.value)+Number(this.step);const t=this.convertToConstrainedValue(e);const i=t>Number(this.min)?`${t}`:`${this.min}`;this.value=i}setThumbPositionForOrientation(e){const t=da(Number(this.value),Number(this.min),Number(this.max),e);const i=(1-t)*100;if(this.orientation===Yn.t.horizontal){this.position=this.isDragging?`right: ${i}%; transition: none;`:`right: ${i}%; transition: all 0.2s ease;`}else{this.position=this.isDragging?`bottom: ${i}%; transition: none;`:`bottom: ${i}%; transition: all 0.2s ease;`}}updateStepMultiplier(){const e=this.step+"";const t=!!(this.step%1)?e.length-e.indexOf(".")-1:0;this.stepMultiplier=Math.pow(10,t)}get midpoint(){return`${this.convertToConstrainedValue((this.max+this.min)/2)}`}setupDefaultValue(){if(typeof this.value==="string"){if(this.value.length===0){this.initialValue=this.midpoint}else{const e=parseFloat(this.value);if(!Number.isNaN(e)&&(e<this.min||e>this.max)){this.value=this.midpoint}}}}}f([(0,s.attr)({attribute:"readonly",mode:"boolean"})],va.prototype,"readOnly",void 0);f([s.observable],va.prototype,"direction",void 0);f([s.observable],va.prototype,"isDragging",void 0);f([s.observable],va.prototype,"position",void 0);f([s.observable],va.prototype,"trackWidth",void 0);f([s.observable],va.prototype,"trackMinWidth",void 0);f([s.observable],va.prototype,"trackHeight",void 0);f([s.observable],va.prototype,"trackLeft",void 0);f([s.observable],va.prototype,"trackMinHeight",void 0);f([s.observable],va.prototype,"valueTextFormatter",void 0);f([(0,s.attr)({converter:s.nullableNumberConverter})],va.prototype,"min",void 0);f([(0,s.attr)({converter:s.nullableNumberConverter})],va.prototype,"max",void 0);f([(0,s.attr)({converter:s.nullableNumberConverter})],va.prototype,"step",void 0);f([s.attr],va.prototype,"orientation",void 0);f([s.attr],va.prototype,"mode",void 0);const ba=(e,t)=>(0,s.html)`
    <template
        role="switch"
        aria-checked="${e=>e.checked}"
        aria-disabled="${e=>e.disabled}"
        aria-readonly="${e=>e.readOnly}"
        tabindex="${e=>e.disabled?null:0}"
        @keypress="${(e,t)=>e.keypressHandler(t.event)}"
        @click="${(e,t)=>e.clickHandler(t.event)}"
        class="${e=>e.checked?"checked":""}"
    >
        <label
            part="label"
            class="${e=>e.defaultSlottedNodes&&e.defaultSlottedNodes.length?"label":"label label__hidden"}"
        >
            <slot ${(0,s.slotted)("defaultSlottedNodes")}></slot>
        </label>
        <div part="switch" class="switch">
            <slot name="switch">${t.switch||""}</slot>
        </div>
        <span class="status-message" part="status-message">
            <span class="checked-message" part="checked-message">
                <slot name="checked-message"></slot>
            </span>
            <span class="unchecked-message" part="unchecked-message">
                <slot name="unchecked-message"></slot>
            </span>
        </span>
    </template>
`;class ga extends Fe{}class ya extends(Gs(ga)){constructor(){super(...arguments);this.proxy=document.createElement("input")}}class Ca extends ya{constructor(){super();this.initialValue="on";this.keypressHandler=e=>{if(this.readOnly){return}switch(e.key){case ze.Mm:case ze.gG:this.checked=!this.checked;break}};this.clickHandler=e=>{if(!this.disabled&&!this.readOnly){this.checked=!this.checked}};this.proxy.setAttribute("type","checkbox")}readOnlyChanged(){if(this.proxy instanceof HTMLInputElement){this.proxy.readOnly=this.readOnly}this.readOnly?this.classList.add("readonly"):this.classList.remove("readonly")}checkedChanged(e,t){super.checkedChanged(e,t);this.checked?this.classList.add("checked"):this.classList.remove("checked")}}f([(0,s.attr)({attribute:"readonly",mode:"boolean"})],Ca.prototype,"readOnly",void 0);f([s.observable],Ca.prototype,"defaultSlottedNodes",void 0);const xa=(e,t)=>(0,s.html)`
    <template slot="tabpanel" role="tabpanel">
        <slot></slot>
    </template>
`;class wa extends Fe{}const $a=(e,t)=>(0,s.html)`
    <template slot="tab" role="tab" aria-disabled="${e=>e.disabled}">
        <slot></slot>
    </template>
`;class Ia extends Fe{}f([(0,s.attr)({mode:"boolean"})],Ia.prototype,"disabled",void 0);const ka=(e,t)=>(0,s.html)`
    <template class="${e=>e.orientation}">
        ${r(e,t)}
        <div class="tablist" part="tablist" role="tablist">
            <slot class="tab" name="tab" part="tab" ${(0,s.slotted)("tabs")}></slot>

            ${(0,s.when)((e=>e.showActiveIndicator),(0,s.html)`
                    <div
                        ${(0,s.ref)("activeIndicatorRef")}
                        class="activeIndicator"
                        part="activeIndicator"
                    ></div>
                `)}
        </div>
        ${n(e,t)}
        <div class="tabpanel" part="tabpanel">
            <slot name="tabpanel" ${(0,s.slotted)("tabpanels")}></slot>
        </div>
    </template>
`;const Oa={vertical:"vertical",horizontal:"horizontal"};class Ta extends Fe{constructor(){super(...arguments);this.orientation=Oa.horizontal;this.activeindicator=true;this.showActiveIndicator=true;this.prevActiveTabIndex=0;this.activeTabIndex=0;this.ticking=false;this.change=()=>{this.$emit("change",this.activetab)};this.isDisabledElement=e=>e.getAttribute("aria-disabled")==="true";this.isHiddenElement=e=>e.hasAttribute("hidden");this.isFocusableElement=e=>!this.isDisabledElement(e)&&!this.isHiddenElement(e);this.setTabs=()=>{const e="gridColumn";const t="gridRow";const i=this.isHorizontal()?e:t;this.activeTabIndex=this.getActiveIndex();this.showActiveIndicator=false;this.tabs.forEach(((s,o)=>{if(s.slot==="tab"){const e=this.activeTabIndex===o&&this.isFocusableElement(s);if(this.activeindicator&&this.isFocusableElement(s)){this.showActiveIndicator=true}const t=this.tabIds[o];const i=this.tabpanelIds[o];s.setAttribute("id",t);s.setAttribute("aria-selected",e?"true":"false");s.setAttribute("aria-controls",i);s.addEventListener("click",this.handleTabClick);s.addEventListener("keydown",this.handleTabKeyDown);s.setAttribute("tabindex",e?"0":"-1");if(e){this.activetab=s;this.activeid=t}}s.style[e]="";s.style[t]="";s.style[i]=`${o+1}`;!this.isHorizontal()?s.classList.add("vertical"):s.classList.remove("vertical")}))};this.setTabPanels=()=>{this.tabpanels.forEach(((e,t)=>{const i=this.tabIds[t];const s=this.tabpanelIds[t];e.setAttribute("id",s);e.setAttribute("aria-labelledby",i);this.activeTabIndex!==t?e.setAttribute("hidden",""):e.removeAttribute("hidden")}))};this.handleTabClick=e=>{const t=e.currentTarget;if(t.nodeType===1&&this.isFocusableElement(t)){this.prevActiveTabIndex=this.activeTabIndex;this.activeTabIndex=this.tabs.indexOf(t);this.setComponent()}};this.handleTabKeyDown=e=>{if(this.isHorizontal()){switch(e.key){case ze.kT:e.preventDefault();this.adjustBackward(e);break;case ze.bb:e.preventDefault();this.adjustForward(e);break}}else{switch(e.key){case ze.I5:e.preventDefault();this.adjustBackward(e);break;case ze.HX:e.preventDefault();this.adjustForward(e);break}}switch(e.key){case ze.Tg:e.preventDefault();this.adjust(-this.activeTabIndex);break;case ze.FM:e.preventDefault();this.adjust(this.tabs.length-this.activeTabIndex-1);break}};this.adjustForward=e=>{const t=this.tabs;let i=0;i=this.activetab?t.indexOf(this.activetab)+1:1;if(i===t.length){i=0}while(i<t.length&&t.length>1){if(this.isFocusableElement(t[i])){this.moveToTabByIndex(t,i);break}else if(this.activetab&&i===t.indexOf(this.activetab)){break}else if(i+1>=t.length){i=0}else{i+=1}}};this.adjustBackward=e=>{const t=this.tabs;let i=0;i=this.activetab?t.indexOf(this.activetab)-1:0;i=i<0?t.length-1:i;while(i>=0&&t.length>1){if(this.isFocusableElement(t[i])){this.moveToTabByIndex(t,i);break}else if(i-1<0){i=t.length-1}else{i-=1}}};this.moveToTabByIndex=(e,t)=>{const i=e[t];this.activetab=i;this.prevActiveTabIndex=this.activeTabIndex;this.activeTabIndex=t;i.focus();this.setComponent()}}orientationChanged(){if(this.$fastController.isConnected){this.setTabs();this.setTabPanels();this.handleActiveIndicatorPosition()}}activeidChanged(e,t){if(this.$fastController.isConnected&&this.tabs.length<=this.tabpanels.length){this.prevActiveTabIndex=this.tabs.findIndex((t=>t.id===e));this.setTabs();this.setTabPanels();this.handleActiveIndicatorPosition()}}tabsChanged(){if(this.$fastController.isConnected&&this.tabs.length<=this.tabpanels.length){this.tabIds=this.getTabIds();this.tabpanelIds=this.getTabPanelIds();this.setTabs();this.setTabPanels();this.handleActiveIndicatorPosition()}}tabpanelsChanged(){if(this.$fastController.isConnected&&this.tabpanels.length<=this.tabs.length){this.tabIds=this.getTabIds();this.tabpanelIds=this.getTabPanelIds();this.setTabs();this.setTabPanels();this.handleActiveIndicatorPosition()}}getActiveIndex(){const e=this.activeid;if(e!==undefined){return this.tabIds.indexOf(this.activeid)===-1?0:this.tabIds.indexOf(this.activeid)}else{return 0}}getTabIds(){return this.tabs.map((e=>{var t;return(t=e.getAttribute("id"))!==null&&t!==void 0?t:`tab-${Io()}`}))}getTabPanelIds(){return this.tabpanels.map((e=>{var t;return(t=e.getAttribute("id"))!==null&&t!==void 0?t:`panel-${Io()}`}))}setComponent(){if(this.activeTabIndex!==this.prevActiveTabIndex){this.activeid=this.tabIds[this.activeTabIndex];this.focusTab();this.change()}}isHorizontal(){return this.orientation===Oa.horizontal}handleActiveIndicatorPosition(){if(this.showActiveIndicator&&this.activeindicator&&this.activeTabIndex!==this.prevActiveTabIndex){if(this.ticking){this.ticking=false}else{this.ticking=true;this.animateActiveIndicator()}}}animateActiveIndicator(){this.ticking=true;const e=this.isHorizontal()?"gridColumn":"gridRow";const t=this.isHorizontal()?"translateX":"translateY";const i=this.isHorizontal()?"offsetLeft":"offsetTop";const s=this.activeIndicatorRef[i];this.activeIndicatorRef.style[e]=`${this.activeTabIndex+1}`;const o=this.activeIndicatorRef[i];this.activeIndicatorRef.style[e]=`${this.prevActiveTabIndex+1}`;const n=o-s;this.activeIndicatorRef.style.transform=`${t}(${n}px)`;this.activeIndicatorRef.classList.add("activeIndicatorTransition");this.activeIndicatorRef.addEventListener("transitionend",(()=>{this.ticking=false;this.activeIndicatorRef.style[e]=`${this.activeTabIndex+1}`;this.activeIndicatorRef.style.transform=`${t}(0px)`;this.activeIndicatorRef.classList.remove("activeIndicatorTransition")}))}adjust(e){const t=this.tabs.filter((e=>this.isFocusableElement(e)));const i=t.indexOf(this.activetab);const s=(0,Ne.AB)(0,t.length-1,i+e);const o=this.tabs.indexOf(t[s]);if(o>-1){this.moveToTabByIndex(this.tabs,o)}}focusTab(){this.tabs[this.activeTabIndex].focus()}connectedCallback(){super.connectedCallback();this.tabIds=this.getTabIds();this.tabpanelIds=this.getTabPanelIds();this.activeTabIndex=this.getActiveIndex()}}f([s.attr],Ta.prototype,"orientation",void 0);f([s.attr],Ta.prototype,"activeid",void 0);f([s.observable],Ta.prototype,"tabs",void 0);f([s.observable],Ta.prototype,"tabpanels",void 0);f([(0,s.attr)({mode:"boolean"})],Ta.prototype,"activeindicator",void 0);f([s.observable],Ta.prototype,"activeIndicatorRef",void 0);f([s.observable],Ta.prototype,"showActiveIndicator",void 0);Pe(Ta,o);const Ea={none:"none",both:"both",horizontal:"horizontal",vertical:"vertical"};const Ra=(e,t)=>(0,s.html)`
    <template
        class="
            ${e=>e.readOnly?"readonly":""}
            ${e=>e.resize!==Ea.none?`resize-${e.resize}`:""}"
    >
        <label
            part="label"
            for="control"
            class="${e=>e.defaultSlottedNodes&&e.defaultSlottedNodes.length?"label":"label label__hidden"}"
        >
            <slot ${(0,s.slotted)("defaultSlottedNodes")}></slot>
        </label>
        <textarea
            part="control"
            class="control"
            id="control"
            ?autofocus="${e=>e.autofocus}"
            cols="${e=>e.cols}"
            ?disabled="${e=>e.disabled}"
            form="${e=>e.form}"
            list="${e=>e.list}"
            maxlength="${e=>e.maxlength}"
            minlength="${e=>e.minlength}"
            name="${e=>e.name}"
            placeholder="${e=>e.placeholder}"
            ?readonly="${e=>e.readOnly}"
            ?required="${e=>e.required}"
            rows="${e=>e.rows}"
            ?spellcheck="${e=>e.spellcheck}"
            :value="${e=>e.value}"
            aria-atomic="${e=>e.ariaAtomic}"
            aria-busy="${e=>e.ariaBusy}"
            aria-controls="${e=>e.ariaControls}"
            aria-current="${e=>e.ariaCurrent}"
            aria-describedby="${e=>e.ariaDescribedby}"
            aria-details="${e=>e.ariaDetails}"
            aria-disabled="${e=>e.ariaDisabled}"
            aria-errormessage="${e=>e.ariaErrormessage}"
            aria-flowto="${e=>e.ariaFlowto}"
            aria-haspopup="${e=>e.ariaHaspopup}"
            aria-hidden="${e=>e.ariaHidden}"
            aria-invalid="${e=>e.ariaInvalid}"
            aria-keyshortcuts="${e=>e.ariaKeyshortcuts}"
            aria-label="${e=>e.ariaLabel}"
            aria-labelledby="${e=>e.ariaLabelledby}"
            aria-live="${e=>e.ariaLive}"
            aria-owns="${e=>e.ariaOwns}"
            aria-relevant="${e=>e.ariaRelevant}"
            aria-roledescription="${e=>e.ariaRoledescription}"
            @input="${(e,t)=>e.handleTextInput()}"
            @change="${e=>e.handleChange()}"
            ${(0,s.ref)("control")}
        ></textarea>
    </template>
`;class Da extends Fe{}class Sa extends(Ws(Da)){constructor(){super(...arguments);this.proxy=document.createElement("textarea")}}class Aa extends Sa{constructor(){super(...arguments);this.resize=Ea.none;this.cols=20;this.handleTextInput=()=>{this.value=this.control.value}}readOnlyChanged(){if(this.proxy instanceof HTMLTextAreaElement){this.proxy.readOnly=this.readOnly}}autofocusChanged(){if(this.proxy instanceof HTMLTextAreaElement){this.proxy.autofocus=this.autofocus}}listChanged(){if(this.proxy instanceof HTMLTextAreaElement){this.proxy.setAttribute("list",this.list)}}maxlengthChanged(){if(this.proxy instanceof HTMLTextAreaElement){this.proxy.maxLength=this.maxlength}}minlengthChanged(){if(this.proxy instanceof HTMLTextAreaElement){this.proxy.minLength=this.minlength}}spellcheckChanged(){if(this.proxy instanceof HTMLTextAreaElement){this.proxy.spellcheck=this.spellcheck}}select(){this.control.select();this.$emit("select")}handleChange(){this.$emit("change")}validate(){super.validate(this.control)}}f([(0,s.attr)({mode:"boolean"})],Aa.prototype,"readOnly",void 0);f([s.attr],Aa.prototype,"resize",void 0);f([(0,s.attr)({mode:"boolean"})],Aa.prototype,"autofocus",void 0);f([(0,s.attr)({attribute:"form"})],Aa.prototype,"formId",void 0);f([s.attr],Aa.prototype,"list",void 0);f([(0,s.attr)({converter:s.nullableNumberConverter})],Aa.prototype,"maxlength",void 0);f([(0,s.attr)({converter:s.nullableNumberConverter})],Aa.prototype,"minlength",void 0);f([s.attr],Aa.prototype,"name",void 0);f([s.attr],Aa.prototype,"placeholder",void 0);f([(0,s.attr)({converter:s.nullableNumberConverter,mode:"fromView"})],Aa.prototype,"cols",void 0);f([(0,s.attr)({converter:s.nullableNumberConverter,mode:"fromView"})],Aa.prototype,"rows",void 0);f([(0,s.attr)({mode:"boolean"})],Aa.prototype,"spellcheck",void 0);f([s.observable],Aa.prototype,"defaultSlottedNodes",void 0);Pe(Aa,Fr);const Fa=(e,t)=>(0,s.html)`
    <template
        class="
            ${e=>e.readOnly?"readonly":""}
        "
    >
        <label
            part="label"
            for="control"
            class="${e=>e.defaultSlottedNodes&&e.defaultSlottedNodes.length?"label":"label label__hidden"}"
        >
            <slot
                ${(0,s.slotted)({property:"defaultSlottedNodes",filter:Xr})}
            ></slot>
        </label>
        <div class="root" part="root">
            ${r(e,t)}
            <input
                class="control"
                part="control"
                id="control"
                @input="${e=>e.handleTextInput()}"
                @change="${e=>e.handleChange()}"
                ?autofocus="${e=>e.autofocus}"
                ?disabled="${e=>e.disabled}"
                list="${e=>e.list}"
                maxlength="${e=>e.maxlength}"
                minlength="${e=>e.minlength}"
                pattern="${e=>e.pattern}"
                placeholder="${e=>e.placeholder}"
                ?readonly="${e=>e.readOnly}"
                ?required="${e=>e.required}"
                size="${e=>e.size}"
                ?spellcheck="${e=>e.spellcheck}"
                :value="${e=>e.value}"
                type="${e=>e.type}"
                aria-atomic="${e=>e.ariaAtomic}"
                aria-busy="${e=>e.ariaBusy}"
                aria-controls="${e=>e.ariaControls}"
                aria-current="${e=>e.ariaCurrent}"
                aria-describedby="${e=>e.ariaDescribedby}"
                aria-details="${e=>e.ariaDetails}"
                aria-disabled="${e=>e.ariaDisabled}"
                aria-errormessage="${e=>e.ariaErrormessage}"
                aria-flowto="${e=>e.ariaFlowto}"
                aria-haspopup="${e=>e.ariaHaspopup}"
                aria-hidden="${e=>e.ariaHidden}"
                aria-invalid="${e=>e.ariaInvalid}"
                aria-keyshortcuts="${e=>e.ariaKeyshortcuts}"
                aria-label="${e=>e.ariaLabel}"
                aria-labelledby="${e=>e.ariaLabelledby}"
                aria-live="${e=>e.ariaLive}"
                aria-owns="${e=>e.ariaOwns}"
                aria-relevant="${e=>e.ariaRelevant}"
                aria-roledescription="${e=>e.ariaRoledescription}"
                ${(0,s.ref)("control")}
            />
            ${n(e,t)}
        </div>
    </template>
`;const La=(e,t)=>(0,s.html)`
    <template
        aria-label="${e=>e.ariaLabel}"
        aria-labelledby="${e=>e.ariaLabelledby}"
        aria-orientation="${e=>e.orientation}"
        orientation="${e=>e.orientation}"
        role="toolbar"
        @mousedown="${(e,t)=>e.mouseDownHandler(t.event)}"
        @focusin="${(e,t)=>e.focusinHandler(t.event)}"
        @keydown="${(e,t)=>e.keydownHandler(t.event)}"
        ${(0,s.children)({property:"childItems",attributeFilter:["disabled","hidden"],filter:(0,s.elements)(),subtree:true})}
    >
        <slot name="label"></slot>
        <div class="positioning-region" part="positioning-region">
            ${r(e,t)}
            <slot
                ${(0,s.slotted)({filter:(0,s.elements)(),property:"slottedItems"})}
            ></slot>
            ${n(e,t)}
        </div>
    </template>
`;const Ma=Object.freeze({[ze.Is.ArrowUp]:{[Yn.t.vertical]:-1},[ze.Is.ArrowDown]:{[Yn.t.vertical]:1},[ze.Is.ArrowLeft]:{[Yn.t.horizontal]:{[Ge.O.ltr]:-1,[Ge.O.rtl]:1}},[ze.Is.ArrowRight]:{[Yn.t.horizontal]:{[Ge.O.ltr]:1,[Ge.O.rtl]:-1}}});class Pa extends Fe{constructor(){super(...arguments);this._activeIndex=0;this.direction=Ge.O.ltr;this.orientation=Yn.t.horizontal}get activeIndex(){s.Observable.track(this,"activeIndex");return this._activeIndex}set activeIndex(e){if(this.$fastController.isConnected){this._activeIndex=(0,Ne.AB)(0,this.focusableElements.length-1,e);s.Observable.notify(this,"activeIndex")}}slottedItemsChanged(){if(this.$fastController.isConnected){this.reduceFocusableElements()}}mouseDownHandler(e){var t;const i=(t=this.focusableElements)===null||t===void 0?void 0:t.findIndex((t=>t.contains(e.target)));if(i>-1&&this.activeIndex!==i){this.setFocusedElement(i)}return true}childItemsChanged(e,t){if(this.$fastController.isConnected){this.reduceFocusableElements()}}connectedCallback(){super.connectedCallback();this.direction=Is(this)}focusinHandler(e){const t=e.relatedTarget;if(!t||this.contains(t)){return}this.setFocusedElement()}getDirectionalIncrementer(e){var t,i,s,o,n;return(n=(s=(i=(t=Ma[e])===null||t===void 0?void 0:t[this.orientation])===null||i===void 0?void 0:i[this.direction])!==null&&s!==void 0?s:(o=Ma[e])===null||o===void 0?void 0:o[this.orientation])!==null&&n!==void 0?n:0}keydownHandler(e){const t=e.key;if(!(t in ze.Is)||e.defaultPrevented||e.shiftKey){return true}const i=this.getDirectionalIncrementer(t);if(!i){return!e.target.closest("[role=radiogroup]")}const s=this.activeIndex+i;if(this.focusableElements[s]){e.preventDefault()}this.setFocusedElement(s);return true}get allSlottedItems(){return[...this.start.assignedElements(),...this.slottedItems,...this.end.assignedElements()]}reduceFocusableElements(){var e;const t=(e=this.focusableElements)===null||e===void 0?void 0:e[this.activeIndex];this.focusableElements=this.allSlottedItems.reduce(Pa.reduceFocusableItems,[]);const i=this.focusableElements.indexOf(t);this.activeIndex=Math.max(0,i);this.setFocusableElements()}setFocusedElement(e=this.activeIndex){var t;this.activeIndex=e;this.setFocusableElements();(t=this.focusableElements[this.activeIndex])===null||t===void 0?void 0:t.focus()}static reduceFocusableItems(e,t){var i,s,o,n;const r=t.getAttribute("role")==="radio";const a=(s=(i=t.$fastController)===null||i===void 0?void 0:i.definition.shadowOptions)===null||s===void 0?void 0:s.delegatesFocus;const l=Array.from((n=(o=t.shadowRoot)===null||o===void 0?void 0:o.querySelectorAll("*"))!==null&&n!==void 0?n:[]).some((e=>(0,Bn.tp)(e)));if(!t.hasAttribute("disabled")&&!t.hasAttribute("hidden")&&((0,Bn.tp)(t)||r||a||l)){e.push(t);return e}if(t.childElementCount){return e.concat(Array.from(t.children).reduce(Pa.reduceFocusableItems,[]))}return e}setFocusableElements(){if(this.$fastController.isConnected&&this.focusableElements.length>0){this.focusableElements.forEach(((e,t)=>{e.tabIndex=this.activeIndex===t?0:-1}))}}}f([s.observable],Pa.prototype,"direction",void 0);f([s.attr],Pa.prototype,"orientation",void 0);f([s.observable],Pa.prototype,"slottedItems",void 0);f([s.observable],Pa.prototype,"slottedLabel",void 0);f([s.observable],Pa.prototype,"childItems",void 0);class Ha{}f([(0,s.attr)({attribute:"aria-labelledby"})],Ha.prototype,"ariaLabelledby",void 0);f([(0,s.attr)({attribute:"aria-label"})],Ha.prototype,"ariaLabel",void 0);Pe(Ha,je);Pe(Pa,o,Ha);const Va=(e,t)=>(0,s.html)`
        ${(0,s.when)((e=>e.tooltipVisible),(0,s.html)`
            <${e.tagFor(Os)}
                fixed-placement="true"
                auto-update-mode="${e=>e.autoUpdateMode}"
                vertical-positioning-mode="${e=>e.verticalPositioningMode}"
                vertical-default-position="${e=>e.verticalDefaultPosition}"
                vertical-inset="${e=>e.verticalInset}"
                vertical-scaling="${e=>e.verticalScaling}"
                horizontal-positioning-mode="${e=>e.horizontalPositioningMode}"
                horizontal-default-position="${e=>e.horizontalDefaultPosition}"
                horizontal-scaling="${e=>e.horizontalScaling}"
                horizontal-inset="${e=>e.horizontalInset}"
                vertical-viewport-lock="${e=>e.horizontalViewportLock}"
                horizontal-viewport-lock="${e=>e.verticalViewportLock}"
                dir="${e=>e.currentDirection}"
                ${(0,s.ref)("region")}
            >
                <div class="tooltip" part="tooltip" role="tooltip">
                    <slot></slot>
                </div>
            </${e.tagFor(Os)}>
        `)}
    `;const za={top:"top",right:"right",bottom:"bottom",left:"left",start:"start",end:"end",topLeft:"top-left",topRight:"top-right",bottomLeft:"bottom-left",bottomRight:"bottom-right",topStart:"top-start",topEnd:"top-end",bottomStart:"bottom-start",bottomEnd:"bottom-end"};class Na extends Fe{constructor(){super(...arguments);this.anchor="";this.delay=300;this.autoUpdateMode="anchor";this.anchorElement=null;this.viewportElement=null;this.verticalPositioningMode="dynamic";this.horizontalPositioningMode="dynamic";this.horizontalInset="false";this.verticalInset="false";this.horizontalScaling="content";this.verticalScaling="content";this.verticalDefaultPosition=undefined;this.horizontalDefaultPosition=undefined;this.tooltipVisible=false;this.currentDirection=Ge.O.ltr;this.showDelayTimer=null;this.hideDelayTimer=null;this.isAnchorHoveredFocused=false;this.isRegionHovered=false;this.handlePositionChange=e=>{this.classList.toggle("top",this.region.verticalPosition==="start");this.classList.toggle("bottom",this.region.verticalPosition==="end");this.classList.toggle("inset-top",this.region.verticalPosition==="insetStart");this.classList.toggle("inset-bottom",this.region.verticalPosition==="insetEnd");this.classList.toggle("center-vertical",this.region.verticalPosition==="center");this.classList.toggle("left",this.region.horizontalPosition==="start");this.classList.toggle("right",this.region.horizontalPosition==="end");this.classList.toggle("inset-left",this.region.horizontalPosition==="insetStart");this.classList.toggle("inset-right",this.region.horizontalPosition==="insetEnd");this.classList.toggle("center-horizontal",this.region.horizontalPosition==="center")};this.handleRegionMouseOver=e=>{this.isRegionHovered=true};this.handleRegionMouseOut=e=>{this.isRegionHovered=false;this.startHideDelayTimer()};this.handleAnchorMouseOver=e=>{if(this.tooltipVisible){this.isAnchorHoveredFocused=true;return}this.startShowDelayTimer()};this.handleAnchorMouseOut=e=>{this.isAnchorHoveredFocused=false;this.clearShowDelayTimer();this.startHideDelayTimer()};this.handleAnchorFocusIn=e=>{this.startShowDelayTimer()};this.handleAnchorFocusOut=e=>{this.isAnchorHoveredFocused=false;this.clearShowDelayTimer();this.startHideDelayTimer()};this.startHideDelayTimer=()=>{this.clearHideDelayTimer();if(!this.tooltipVisible){return}this.hideDelayTimer=window.setTimeout((()=>{this.updateTooltipVisibility()}),60)};this.clearHideDelayTimer=()=>{if(this.hideDelayTimer!==null){clearTimeout(this.hideDelayTimer);this.hideDelayTimer=null}};this.startShowDelayTimer=()=>{if(this.isAnchorHoveredFocused){return}if(this.delay>1){if(this.showDelayTimer===null)this.showDelayTimer=window.setTimeout((()=>{this.startHover()}),this.delay);return}this.startHover()};this.startHover=()=>{this.isAnchorHoveredFocused=true;this.updateTooltipVisibility()};this.clearShowDelayTimer=()=>{if(this.showDelayTimer!==null){clearTimeout(this.showDelayTimer);this.showDelayTimer=null}};this.getAnchor=()=>{const e=this.getRootNode();if(e instanceof ShadowRoot){return e.getElementById(this.anchor)}return document.getElementById(this.anchor)};this.handleDocumentKeydown=e=>{if(!e.defaultPrevented&&this.tooltipVisible){switch(e.key){case ze.F9:this.isAnchorHoveredFocused=false;this.updateTooltipVisibility();this.$emit("dismiss");break}}};this.updateTooltipVisibility=()=>{if(this.visible===false){this.hideTooltip()}else if(this.visible===true){this.showTooltip();return}else{if(this.isAnchorHoveredFocused||this.isRegionHovered){this.showTooltip();return}this.hideTooltip()}};this.showTooltip=()=>{if(this.tooltipVisible){return}this.currentDirection=Is(this);this.tooltipVisible=true;document.addEventListener("keydown",this.handleDocumentKeydown);s.DOM.queueUpdate(this.setRegionProps)};this.hideTooltip=()=>{if(!this.tooltipVisible){return}this.clearHideDelayTimer();if(this.region!==null&&this.region!==undefined){this.region.removeEventListener("positionchange",this.handlePositionChange);this.region.viewportElement=null;this.region.anchorElement=null;this.region.removeEventListener("mouseover",this.handleRegionMouseOver);this.region.removeEventListener("mouseout",this.handleRegionMouseOut)}document.removeEventListener("keydown",this.handleDocumentKeydown);this.tooltipVisible=false};this.setRegionProps=()=>{if(!this.tooltipVisible){return}this.region.viewportElement=this.viewportElement;this.region.anchorElement=this.anchorElement;this.region.addEventListener("positionchange",this.handlePositionChange);this.region.addEventListener("mouseover",this.handleRegionMouseOver,{passive:true});this.region.addEventListener("mouseout",this.handleRegionMouseOut,{passive:true})}}visibleChanged(){if(this.$fastController.isConnected){this.updateTooltipVisibility();this.updateLayout()}}anchorChanged(){if(this.$fastController.isConnected){this.anchorElement=this.getAnchor()}}positionChanged(){if(this.$fastController.isConnected){this.updateLayout()}}anchorElementChanged(e){if(this.$fastController.isConnected){if(e!==null&&e!==undefined){e.removeEventListener("mouseover",this.handleAnchorMouseOver);e.removeEventListener("mouseout",this.handleAnchorMouseOut);e.removeEventListener("focusin",this.handleAnchorFocusIn);e.removeEventListener("focusout",this.handleAnchorFocusOut)}if(this.anchorElement!==null&&this.anchorElement!==undefined){this.anchorElement.addEventListener("mouseover",this.handleAnchorMouseOver,{passive:true});this.anchorElement.addEventListener("mouseout",this.handleAnchorMouseOut,{passive:true});this.anchorElement.addEventListener("focusin",this.handleAnchorFocusIn,{passive:true});this.anchorElement.addEventListener("focusout",this.handleAnchorFocusOut,{passive:true});const e=this.anchorElement.id;if(this.anchorElement.parentElement!==null){this.anchorElement.parentElement.querySelectorAll(":hover").forEach((t=>{if(t.id===e){this.startShowDelayTimer()}}))}}if(this.region!==null&&this.region!==undefined&&this.tooltipVisible){this.region.anchorElement=this.anchorElement}this.updateLayout()}}viewportElementChanged(){if(this.region!==null&&this.region!==undefined){this.region.viewportElement=this.viewportElement}this.updateLayout()}connectedCallback(){super.connectedCallback();this.anchorElement=this.getAnchor();this.updateTooltipVisibility()}disconnectedCallback(){this.hideTooltip();this.clearShowDelayTimer();this.clearHideDelayTimer();super.disconnectedCallback()}updateLayout(){this.verticalPositioningMode="locktodefault";this.horizontalPositioningMode="locktodefault";switch(this.position){case za.top:case za.bottom:this.verticalDefaultPosition=this.position;this.horizontalDefaultPosition="center";break;case za.right:case za.left:case za.start:case za.end:this.verticalDefaultPosition="center";this.horizontalDefaultPosition=this.position;break;case za.topLeft:this.verticalDefaultPosition="top";this.horizontalDefaultPosition="left";break;case za.topRight:this.verticalDefaultPosition="top";this.horizontalDefaultPosition="right";break;case za.bottomLeft:this.verticalDefaultPosition="bottom";this.horizontalDefaultPosition="left";break;case za.bottomRight:this.verticalDefaultPosition="bottom";this.horizontalDefaultPosition="right";break;case za.topStart:this.verticalDefaultPosition="top";this.horizontalDefaultPosition="start";break;case za.topEnd:this.verticalDefaultPosition="top";this.horizontalDefaultPosition="end";break;case za.bottomStart:this.verticalDefaultPosition="bottom";this.horizontalDefaultPosition="start";break;case za.bottomEnd:this.verticalDefaultPosition="bottom";this.horizontalDefaultPosition="end";break;default:this.verticalPositioningMode="dynamic";this.horizontalPositioningMode="dynamic";this.verticalDefaultPosition=void 0;this.horizontalDefaultPosition="center";break}}}f([(0,s.attr)({mode:"boolean"})],Na.prototype,"visible",void 0);f([s.attr],Na.prototype,"anchor",void 0);f([s.attr],Na.prototype,"delay",void 0);f([s.attr],Na.prototype,"position",void 0);f([(0,s.attr)({attribute:"auto-update-mode"})],Na.prototype,"autoUpdateMode",void 0);f([(0,s.attr)({attribute:"horizontal-viewport-lock"})],Na.prototype,"horizontalViewportLock",void 0);f([(0,s.attr)({attribute:"vertical-viewport-lock"})],Na.prototype,"verticalViewportLock",void 0);f([s.observable],Na.prototype,"anchorElement",void 0);f([s.observable],Na.prototype,"viewportElement",void 0);f([s.observable],Na.prototype,"verticalPositioningMode",void 0);f([s.observable],Na.prototype,"horizontalPositioningMode",void 0);f([s.observable],Na.prototype,"horizontalInset",void 0);f([s.observable],Na.prototype,"verticalInset",void 0);f([s.observable],Na.prototype,"horizontalScaling",void 0);f([s.observable],Na.prototype,"verticalScaling",void 0);f([s.observable],Na.prototype,"verticalDefaultPosition",void 0);f([s.observable],Na.prototype,"horizontalDefaultPosition",void 0);f([s.observable],Na.prototype,"tooltipVisible",void 0);f([s.observable],Na.prototype,"currentDirection",void 0);const Ba=(e,t)=>(0,s.html)`
    <template
        role="treeitem"
        slot="${e=>e.isNestedItem()?"item":void 0}"
        tabindex="-1"
        class="${e=>e.expanded?"expanded":""} ${e=>e.selected?"selected":""} ${e=>e.nested?"nested":""}
            ${e=>e.disabled?"disabled":""}"
        aria-expanded="${e=>e.childItems&&e.childItemLength()>0?e.expanded:void 0}"
        aria-selected="${e=>e.selected}"
        aria-disabled="${e=>e.disabled}"
        @focusin="${(e,t)=>e.handleFocus(t.event)}"
        @focusout="${(e,t)=>e.handleBlur(t.event)}"
        ${(0,s.children)({property:"childItems",filter:(0,s.elements)()})}
    >
        <div class="positioning-region" part="positioning-region">
            <div class="content-region" part="content-region">
                ${(0,s.when)((e=>e.childItems&&e.childItemLength()>0),(0,s.html)`
                        <div
                            aria-hidden="true"
                            class="expand-collapse-button"
                            part="expand-collapse-button"
                            @click="${(e,t)=>e.handleExpandCollapseButtonClick(t.event)}"
                            ${(0,s.ref)("expandCollapseButton")}
                        >
                            <slot name="expand-collapse-glyph">
                                ${t.expandCollapseGlyph||""}
                            </slot>
                        </div>
                    `)}
                ${r(e,t)}
                <slot></slot>
                ${n(e,t)}
            </div>
        </div>
        ${(0,s.when)((e=>e.childItems&&e.childItemLength()>0&&(e.expanded||e.renderCollapsedChildren)),(0,s.html)`
                <div role="group" class="items" part="items">
                    <slot name="item" ${(0,s.slotted)("items")}></slot>
                </div>
            `)}
    </template>
`;function qa(e){return Ao(e)&&e.getAttribute("role")==="treeitem"}class Ua extends Fe{constructor(){super(...arguments);this.expanded=false;this.focusable=false;this.isNestedItem=()=>qa(this.parentElement);this.handleExpandCollapseButtonClick=e=>{if(!this.disabled&&!e.defaultPrevented){this.expanded=!this.expanded}};this.handleFocus=e=>{this.setAttribute("tabindex","0")};this.handleBlur=e=>{this.setAttribute("tabindex","-1")}}expandedChanged(){if(this.$fastController.isConnected){this.$emit("expanded-change",this)}}selectedChanged(){if(this.$fastController.isConnected){this.$emit("selected-change",this)}}itemsChanged(e,t){if(this.$fastController.isConnected){this.items.forEach((e=>{if(qa(e)){e.nested=true}}))}}static focusItem(e){e.focusable=true;e.focus()}childItemLength(){const e=this.childItems.filter((e=>qa(e)));return e?e.length:0}}f([(0,s.attr)({mode:"boolean"})],Ua.prototype,"expanded",void 0);f([(0,s.attr)({mode:"boolean"})],Ua.prototype,"selected",void 0);f([(0,s.attr)({mode:"boolean"})],Ua.prototype,"disabled",void 0);f([s.observable],Ua.prototype,"focusable",void 0);f([s.observable],Ua.prototype,"childItems",void 0);f([s.observable],Ua.prototype,"items",void 0);f([s.observable],Ua.prototype,"nested",void 0);f([s.observable],Ua.prototype,"renderCollapsedChildren",void 0);Pe(Ua,o);const ja=(e,t)=>(0,s.html)`
    <template
        role="tree"
        ${(0,s.ref)("treeView")}
        @keydown="${(e,t)=>e.handleKeyDown(t.event)}"
        @focusin="${(e,t)=>e.handleFocus(t.event)}"
        @focusout="${(e,t)=>e.handleBlur(t.event)}"
        @click="${(e,t)=>e.handleClick(t.event)}"
        @selected-change="${(e,t)=>e.handleSelectedChange(t.event)}"
    >
        <slot ${(0,s.slotted)("slottedTreeItems")}></slot>
    </template>
`;class _a extends Fe{constructor(){super(...arguments);this.currentFocused=null;this.handleFocus=e=>{if(this.slottedTreeItems.length<1){return}if(e.target===this){if(this.currentFocused===null){this.currentFocused=this.getValidFocusableItem()}if(this.currentFocused!==null){Ua.focusItem(this.currentFocused)}return}if(this.contains(e.target)){this.setAttribute("tabindex","-1");this.currentFocused=e.target}};this.handleBlur=e=>{if(e.target instanceof HTMLElement&&(e.relatedTarget===null||!this.contains(e.relatedTarget))){this.setAttribute("tabindex","0")}};this.handleKeyDown=e=>{if(e.defaultPrevented){return}if(this.slottedTreeItems.length<1){return true}const t=this.getVisibleNodes();switch(e.key){case ze.Tg:if(t.length){Ua.focusItem(t[0])}return;case ze.FM:if(t.length){Ua.focusItem(t[t.length-1])}return;case ze.kT:if(e.target&&this.isFocusableElement(e.target)){const t=e.target;if(t instanceof Ua&&t.childItemLength()>0&&t.expanded){t.expanded=false}else if(t instanceof Ua&&t.parentElement instanceof Ua){Ua.focusItem(t.parentElement)}}return false;case ze.bb:if(e.target&&this.isFocusableElement(e.target)){const t=e.target;if(t instanceof Ua&&t.childItemLength()>0&&!t.expanded){t.expanded=true}else if(t instanceof Ua&&t.childItemLength()>0){this.focusNextNode(1,e.target)}}return;case ze.HX:if(e.target&&this.isFocusableElement(e.target)){this.focusNextNode(1,e.target)}return;case ze.I5:if(e.target&&this.isFocusableElement(e.target)){this.focusNextNode(-1,e.target)}return;case ze.Mm:this.handleClick(e);return}return true};this.handleSelectedChange=e=>{if(e.defaultPrevented){return}if(!(e.target instanceof Element)||!qa(e.target)){return true}const t=e.target;if(t.selected){if(this.currentSelected&&this.currentSelected!==t){this.currentSelected.selected=false}this.currentSelected=t}else if(!t.selected&&this.currentSelected===t){this.currentSelected=null}return};this.setItems=()=>{const e=this.treeView.querySelector("[aria-selected='true']");this.currentSelected=e;if(this.currentFocused===null||!this.contains(this.currentFocused)){this.currentFocused=this.getValidFocusableItem()}this.nested=this.checkForNestedItems();const t=this.getVisibleNodes();t.forEach((e=>{if(qa(e)){e.nested=this.nested}}))};this.isFocusableElement=e=>qa(e);this.isSelectedElement=e=>e.selected}slottedTreeItemsChanged(){if(this.$fastController.isConnected){this.setItems()}}connectedCallback(){super.connectedCallback();this.setAttribute("tabindex","0");s.DOM.queueUpdate((()=>{this.setItems()}))}handleClick(e){if(e.defaultPrevented){return}if(!(e.target instanceof Element)||!qa(e.target)){return true}const t=e.target;if(!t.disabled){t.selected=!t.selected}return}focusNextNode(e,t){const i=this.getVisibleNodes();if(!i){return}const s=i[i.indexOf(t)+e];if(Ao(s)){Ua.focusItem(s)}}getValidFocusableItem(){const e=this.getVisibleNodes();let t=e.findIndex(this.isSelectedElement);if(t===-1){t=e.findIndex(this.isFocusableElement)}if(t!==-1){return e[t]}return null}checkForNestedItems(){return this.slottedTreeItems.some((e=>qa(e)&&e.querySelector("[role='treeitem']")))}getVisibleNodes(){return Fo(this,"[role='treeitem']")||[]}}f([(0,s.attr)({attribute:"render-collapsed-nodes"})],_a.prototype,"renderCollapsedNodes",void 0);f([s.observable],_a.prototype,"currentSelected",void 0);f([s.observable],_a.prototype,"slottedTreeItems",void 0);class Ka{constructor(e){this.listenerCache=new WeakMap;this.query=e}bind(e){const{query:t}=this;const i=this.constructListener(e);i.bind(t)();t.addListener(i);this.listenerCache.set(e,i)}unbind(e){const t=this.listenerCache.get(e);if(t){this.query.removeListener(t);this.listenerCache.delete(e)}}}class Wa extends Ka{constructor(e,t){super(e);this.styles=t}static with(e){return t=>new Wa(e,t)}constructListener(e){let t=false;const i=this.styles;return function s(){const{matches:o}=this;if(o&&!t){e.$fastController.addStyles(i);t=o}else if(!o&&t){e.$fastController.removeStyles(i);t=o}}}unbind(e){super.unbind(e);e.$fastController.removeStyles(this.styles)}}const Ga=Wa.with(window.matchMedia("(forced-colors)"));const Xa=Wa.with(window.matchMedia("(prefers-color-scheme: dark)"));const Ya=Wa.with(window.matchMedia("(prefers-color-scheme: light)"));class Qa{constructor(e,t,i){this.propertyName=e;this.value=t;this.styles=i}bind(e){s.Observable.getNotifier(e).subscribe(this,this.propertyName);this.handleChange(e,this.propertyName)}unbind(e){s.Observable.getNotifier(e).unsubscribe(this,this.propertyName);e.$fastController.removeStyles(this.styles)}handleChange(e,t){if(e[t]===this.value){e.$fastController.addStyles(this.styles)}else{e.$fastController.removeStyles(this.styles)}}}const Za="not-allowed";const Ja=`:host([hidden]){display:none}`;function el(e){return`${Ja}:host{display:${e}}`}const tl=Ho()?"focus-visible":"focus"}}]);