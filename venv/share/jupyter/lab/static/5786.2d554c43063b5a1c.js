"use strict";(self.rspackChunk_jupyterlab_application_top=self.rspackChunk_jupyterlab_application_top||[]).push([[5786],{92225(e,t,i){let s,o,n,a,r,l;i.r(t),i.d(t,{ARIAGlobalStatesAndProperties:()=>eI,Accordion:()=>e$,AccordionExpandMode:()=>ex,AccordionItem:()=>ef,Anchor:()=>ek,AnchoredRegion:()=>eD,Avatar:()=>eV,Badge:()=>eB,BaseProgress:()=>iZ,Breadcrumb:()=>e_,BreadcrumbItem:()=>eU,Button:()=>e1,Calendar:()=>e8,CalendarTitleTemplate:()=>ta,Card:()=>tm,CheckableFormAssociated:()=>eZ,Checkbox:()=>tg,Combobox:()=>tS,ComboboxAutocomplete:()=>tD,ComponentPresentation:()=>ec,Container:()=>O,ContainerConfiguration:()=>$,ContainerImpl:()=>Z,DI:()=>E,DataGrid:()=>tn,DataGridCell:()=>ts,DataGridCellTypes:()=>e7,DataGridRow:()=>to,DataGridRowTypes:()=>te,DateFormatter:()=>e4,DefaultComponentPresentation:()=>eu,DefaultResolver:()=>x,DelegatesARIAButton:()=>e5,DelegatesARIACombobox:()=>tA,DelegatesARIALink:()=>eE,DelegatesARIAListbox:()=>tE,DelegatesARIAListboxOption:()=>tI,DelegatesARIASearch:()=>ss,DelegatesARIASelect:()=>sr,DelegatesARIATextbox:()=>iW,DelegatesARIAToolbar:()=>sN,DesignSystem:()=>t9,DesignToken:()=>t1,Dialog:()=>ii,Disclosure:()=>ih,Divider:()=>ip,DividerRole:()=>iu,ElementDisambiguation:()=>t5,FactoryImpl:()=>_,Flipper:()=>ib,FlipperDirection:()=>im,FlyoutPosBottom:()=>eF,FlyoutPosBottomFill:()=>eP,FlyoutPosTallest:()=>eL,FlyoutPosTallestFill:()=>eH,FlyoutPosTop:()=>eA,FlyoutPosTopFill:()=>eM,FormAssociated:()=>eQ,FoundationElement:()=>ep,FoundationElementRegistry:()=>ev,GenerateHeaderOptions:()=>e6,HorizontalScroll:()=>i3,Listbox:()=>tk,ListboxElement:()=>iy,ListboxOption:()=>tw,MatchMediaBehavior:()=>sG,MatchMediaStyleSheetBehavior:()=>sY,Menu:()=>iB,MenuItem:()=>iV,MenuItemRole:()=>iP,NumberField:()=>iY,Picker:()=>iS,PickerList:()=>iI,PickerListItem:()=>iE,PickerMenu:()=>ix,PickerMenuOption:()=>iw,PropertyStyleSheetBehavior:()=>s0,Radio:()=>i2,RadioGroup:()=>i1,Registration:()=>et,ResolverBuilder:()=>y,ResolverImpl:()=>q,Search:()=>si,Select:()=>sa,SelectPosition:()=>tO,ServiceLocator:()=>T,Skeleton:()=>sd,Slider:()=>sy,SliderLabel:()=>sm,SliderMode:()=>sg,StartEnd:()=>c,Switch:()=>sw,Tab:()=>sO,TabPanel:()=>sk,Tabs:()=>sD,TabsOrientation:()=>sR,TextArea:()=>sM,TextAreaResize:()=>sS,TextField:()=>iK,TextFieldType:()=>i_,Toolbar:()=>sV,Tooltip:()=>sU,TooltipPosition:()=>sq,TreeItem:()=>sK,TreeView:()=>sX,accordionItemTemplate:()=>b,accordionTemplate:()=>eg,all:()=>M,anchorTemplate:()=>ew,anchoredRegionTemplate:()=>eO,applyMixins:()=>eb,avatarTemplate:()=>ez,badgeTemplate:()=>eN,breadcrumbItemTemplate:()=>eq,breadcrumbTemplate:()=>ej,buttonTemplate:()=>eK,calendarCellTemplate:()=>tl,calendarRowTemplate:()=>th,calendarTemplate:()=>tu,calendarWeekdayTemplate:()=>tr,cardTemplate:()=>tp,checkboxTemplate:()=>tv,comboboxTemplate:()=>tF,composedContains:()=>tz,composedParent:()=>tH,darkModeStylesheetBehavior:()=>sZ,dataGridCellTemplate:()=>tP,dataGridRowTemplate:()=>tM,dataGridTemplate:()=>tL,dialogTemplate:()=>ie,disabledCursor:()=>s1,disclosureTemplate:()=>il,display:()=>s4,dividerTemplate:()=>id,endSlotTemplate:()=>u,endTemplate:()=>m,flipperTemplate:()=>iv,focusVisible:()=>s8,forcedColorsStylesheetBehavior:()=>sQ,getDirection:()=>eR,hidden:()=>s5,horizontalScrollTemplate:()=>i9,ignore:()=>z,inject:()=>D,interactiveCalendarGridTemplate:()=>td,isListboxOption:()=>t$,isTreeItemElement:()=>s_,lazy:()=>P,lightModeStylesheetBehavior:()=>sJ,listboxOptionTemplate:()=>ig,listboxTemplate:()=>iC,menuItemTemplate:()=>iz,menuTemplate:()=>iN,newInstanceForScope:()=>V,newInstanceOf:()=>N,noninteractiveCalendarTemplate:()=>tc,numberFieldTemplate:()=>iq,optional:()=>H,pickerListItemTemplate:()=>iM,pickerListTemplate:()=>iL,pickerMenuOptionTemplate:()=>iF,pickerMenuTemplate:()=>iA,pickerTemplate:()=>iO,progressRingTemplate:()=>iQ,progressTemplate:()=>iJ,radioGroupTemplate:()=>i0,radioTemplate:()=>i5,reflectAttributes:()=>ir,roleForMenuItem:()=>iH,searchTemplate:()=>i7,selectTemplate:()=>sl,singleton:()=>L,skeletonTemplate:()=>sh,sliderLabelTemplate:()=>sc,sliderTemplate:()=>sv,startSlotTemplate:()=>p,startTemplate:()=>v,supportsElementInternals:()=>eG,switchTemplate:()=>sC,tabPanelTemplate:()=>sI,tabTemplate:()=>sE,tabsTemplate:()=>sT,textAreaTemplate:()=>sA,textFieldTemplate:()=>sP,toolbarTemplate:()=>sH,tooltipTemplate:()=>sB,transient:()=>A,treeItemTemplate:()=>sj,treeViewTemplate:()=>sW,validateKey:()=>ei,whitespaceFilter:()=>i6});var h,d=i(23453);class c{handleStartContentChange(){this.startContainer.classList.toggle("start",this.start.assignedNodes().length>0)}handleEndContentChange(){this.endContainer.classList.toggle("end",this.end.assignedNodes().length>0)}}let u=(e,t)=>(0,d.html)`
    <span
        part="end"
        ${(0,d.ref)("endContainer")}
        class=${e=>t.end?"end":void 0}
    >
        <slot name="end" ${(0,d.ref)("end")} @slotchange="${e=>e.handleEndContentChange()}">
            ${t.end||""}
        </slot>
    </span>
`,p=(e,t)=>(0,d.html)`
    <span
        part="start"
        ${(0,d.ref)("startContainer")}
        class="${e=>t.start?"start":void 0}"
    >
        <slot
            name="start"
            ${(0,d.ref)("start")}
            @slotchange="${e=>e.handleStartContentChange()}"
        >
            ${t.start||""}
        </slot>
    </span>
`,m=(0,d.html)`
    <span part="end" ${(0,d.ref)("endContainer")}>
        <slot
            name="end"
            ${(0,d.ref)("end")}
            @slotchange="${e=>e.handleEndContentChange()}"
        ></slot>
    </span>
`,v=(0,d.html)`
    <span part="start" ${(0,d.ref)("startContainer")}>
        <slot
            name="start"
            ${(0,d.ref)("start")}
            @slotchange="${e=>e.handleStartContentChange()}"
        ></slot>
    </span>
`,b=(e,t)=>(0,d.html)`
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
                ${(0,d.ref)("expandbutton")}
                aria-expanded="${e=>e.expanded}"
                aria-controls="${e=>e.id}-panel"
                id="${e=>e.id}"
                @click="${(e,t)=>e.clickHandler(t.event)}"
            >
                <span class="heading-content" part="heading-content">
                    <slot name="heading"></slot>
                </span>
            </button>
            ${p(e,t)}
            ${u(e,t)}
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
`;function f(e,t,i,s){var o,n=arguments.length,a=n<3?t:null===s?s=Object.getOwnPropertyDescriptor(t,i):s;if("object"==typeof Reflect&&"function"==typeof Reflect.decorate)a=Reflect.decorate(e,t,i,s);else for(var r=e.length-1;r>=0;r--)(o=e[r])&&(a=(n<3?o(a):n>3?o(t,i,a):o(t,i))||a);return n>3&&a&&Object.defineProperty(t,i,a),a}let g=new Map;"metadata"in Reflect||(Reflect.metadata=function(e,t){return function(i){Reflect.defineMetadata(e,t,i)}},Reflect.defineMetadata=function(e,t,i){let s=g.get(i);void 0===s&&g.set(i,s=new Map),s.set(e,t)},Reflect.getOwnMetadata=function(e,t){let i=g.get(t);if(void 0!==i)return i.get(e)});class y{constructor(e,t){this.container=e,this.key=t}instance(e){return this.registerResolver(0,e)}singleton(e){return this.registerResolver(1,e)}transient(e){return this.registerResolver(2,e)}callback(e){return this.registerResolver(3,e)}cachedCallback(e){return this.registerResolver(3,ee(e))}aliasTo(e){return this.registerResolver(5,e)}registerResolver(e,t){let{container:i,key:s}=this;return this.container=this.key=void 0,i.registerResolver(s,new q(s,e,t))}}function C(e){let t,i=e.slice(),s=Object.keys(e),o=s.length;for(let n=0;n<o;++n)el(t=s[n])||(i[t]=e[t]);return i}let x=Object.freeze({none(e){throw Error(`${e.toString()} not registered, did you forget to add @singleton()?`)},singleton:e=>new q(e,1,e),transient:e=>new q(e,2,e)}),$=Object.freeze({default:Object.freeze({parentLocator:()=>null,responsibleForOwnerRequests:!1,defaultResolver:x.singleton})}),w=new Map;function I(e){return t=>Reflect.getOwnMetadata(e,t)}let k=null,E=Object.freeze({createContainer:e=>new Z(null,Object.assign({},$.default,e)),findResponsibleContainer(e){let t=e.$$container$$;return t&&t.responsibleForOwnerRequests?t:E.findParentContainer(e)},findParentContainer(e){let t=new CustomEvent(Y,{bubbles:!0,composed:!0,cancelable:!0,detail:{container:void 0}});return e.dispatchEvent(t),t.detail.container||E.getOrCreateDOMContainer()},getOrCreateDOMContainer:(e,t)=>e?e.$$container$$||new Z(e,Object.assign({},$.default,t,{parentLocator:E.findParentContainer})):k||(k=new Z(null,Object.assign({},$.default,t,{parentLocator:()=>null}))),getDesignParamtypes:I("design:paramtypes"),getAnnotationParamtypes:I("di:paramtypes"),getOrCreateAnnotationParamTypes(e){let t=this.getAnnotationParamtypes(e);return void 0===t&&Reflect.defineMetadata("di:paramtypes",t=[],e),t},getDependencies(e){let t=w.get(e);if(void 0===t){let i=e.inject;if(void 0===i){let i=E.getDesignParamtypes(e),s=E.getAnnotationParamtypes(e);if(void 0===i)if(void 0===s){let i=Object.getPrototypeOf(e);t="function"==typeof i&&i!==Function.prototype?C(E.getDependencies(i)):[]}else t=C(s);else if(void 0===s)t=C(i);else{let e,o;t=C(i);let n=s.length;for(let i=0;i<n;++i)void 0!==(e=s[i])&&(t[i]=e);let a=Object.keys(s);n=a.length;for(let e=0;e<n;++e)el(o=a[e])||(t[o]=s[o])}}else t=C(i);w.set(e,t)}return t},defineProperty(e,t,i,s=!1){let o=`$di_${t}`;Reflect.defineProperty(e,t,{get:function(){let e=this[o];if(void 0===e&&(e=(this instanceof HTMLElement?E.findResponsibleContainer(this):E.getOrCreateDOMContainer()).get(i),this[o]=e,s&&this instanceof d.FASTElement)){let s=this.$fastController,n=()=>{E.findResponsibleContainer(this).get(i)!==this[o]&&(this[o]=e,s.notify(t))};s.subscribe({handleChange:n},"isConnected")}return e}})},createInterface(e,t){let i="function"==typeof e?e:t,s="string"==typeof e?e:e&&"friendlyName"in e&&e.friendlyName||eo,o="string"!=typeof e&&!!e&&"respectConnection"in e&&(e.respectConnection||!1),n=function(e,t,i){if(null==e||new.target!==void 0)throw Error(`No registration for interface: '${n.friendlyName}'`);t?E.defineProperty(e,t,n,o):E.getOrCreateAnnotationParamTypes(e)[i]=n};return n.$isInterface=!0,n.friendlyName=null==s?"(anonymous)":s,null!=i&&(n.register=function(e,t){return i(new y(e,null!=t?t:n))}),n.toString=function(){return`InterfaceSymbol<${n.friendlyName}>`},n},inject:(...e)=>function(t,i,s){if("number"==typeof s){let i=E.getOrCreateAnnotationParamTypes(t),o=e[0];void 0!==o&&(i[s]=o)}else if(i)E.defineProperty(t,i,e[0]);else{let i,o=s?E.getOrCreateAnnotationParamTypes(s.value):E.getOrCreateAnnotationParamTypes(t);for(let t=0;t<e.length;++t)void 0!==(i=e[t])&&(o[t]=i)}},transient:e=>(e.register=function(t){return et.transient(e,e).register(t)},e.registerInRequestor=!1,e),singleton:(e,t=F)=>(e.register=function(t){return et.singleton(e,e).register(t)},e.registerInRequestor=t.scoped,e)}),O=E.createInterface("Container"),T=O;function R(e){return function(t){let i=function(e,t,s){E.inject(i)(e,t,s)};return i.$isResolver=!0,i.resolve=function(i,s){return e(t,i,s)},i}}let D=E.inject;function S(e){return E.transient(e)}function A(e){return null==e?S:S(e)}let F={scoped:!1};function L(e){return"function"==typeof e?E.singleton(e):function(t){return E.singleton(t,e)}}let M=(h=(e,t,i,s)=>i.getAll(e,s),function(e,t){t=!!t;let i=function(e,t,s){E.inject(i)(e,t,s)};return i.$isResolver=!0,i.resolve=function(i,s){return h(e,i,s,t)},i}),P=R((e,t,i)=>()=>i.get(e)),H=R((e,t,i)=>i.has(e,!0)?i.get(e):void 0);function z(e,t,i){E.inject(z)(e,t,i)}z.$isResolver=!0,z.resolve=()=>void 0;let V=R((e,t,i)=>{let s=B(e,t),o=new q(e,0,s);return i.registerResolver(e,o),s}),N=R((e,t,i)=>B(e,t));function B(e,t){return t.getFactory(e).construct(t)}class q{constructor(e,t,i){this.key=e,this.strategy=t,this.state=i,this.resolving=!1}get $isResolver(){return!0}register(e){return e.registerResolver(this.key,this)}resolve(e,t){switch(this.strategy){case 0:return this.state;case 1:if(this.resolving)throw Error(`Cyclic dependency found: ${this.state.name}`);return this.resolving=!0,this.state=e.getFactory(this.state).construct(t),this.strategy=0,this.resolving=!1,this.state;case 2:{let i=e.getFactory(this.state);if(null===i)throw Error(`Resolver for ${String(this.key)} returned a null factory`);return i.construct(t)}case 3:return this.state(e,t,this);case 4:return this.state[0].resolve(e,t);case 5:return t.get(this.state);default:throw Error(`Invalid resolver strategy specified: ${this.strategy}.`)}}getFactory(e){var t,i,s;switch(this.strategy){case 1:case 2:return e.getFactory(this.state);case 5:return null!=(s=null==(i=null==(t=e.getResolver(this.state))?void 0:t.getFactory)?void 0:i.call(t,e))?s:null;default:return null}}}function U(e){return this.get(e)}function j(e,t){return t(e)}class _{constructor(e,t){this.Type=e,this.dependencies=t,this.transformers=null}construct(e,t){let i;return(i=void 0===t?new this.Type(...this.dependencies.map(U,e)):new this.Type(...this.dependencies.map(U,e),...t),null==this.transformers)?i:this.transformers.reduce(j,i)}registerTransformer(e){(this.transformers||(this.transformers=[])).push(e)}}let K={$isResolver:!0,resolve:(e,t)=>t};function W(e){return"function"==typeof e.register}function X(e){return W(e)&&"boolean"==typeof e.registerInRequestor&&e.registerInRequestor}let G=new Set(["Array","ArrayBuffer","Boolean","DataView","Date","Error","EvalError","Float32Array","Float64Array","Function","Int8Array","Int16Array","Int32Array","Map","Number","Object","Promise","RangeError","ReferenceError","RegExp","Set","SharedArrayBuffer","String","SyntaxError","TypeError","Uint8Array","Uint8ClampedArray","Uint16Array","Uint32Array","URIError","WeakMap","WeakSet"]),Y="__DI_LOCATE_PARENT__",Q=new Map;class Z{constructor(e,t){this.owner=e,this.config=t,this._parent=void 0,this.registerDepth=0,this.context=null,null!==e&&(e.$$container$$=this),this.resolvers=new Map,this.resolvers.set(O,K),e instanceof Node&&e.addEventListener(Y,e=>{e.composedPath()[0]!==this.owner&&(e.detail.container=this,e.stopImmediatePropagation())})}get parent(){return void 0===this._parent&&(this._parent=this.config.parentLocator(this.owner)),this._parent}get depth(){return null===this.parent?0:this.parent.depth+1}get responsibleForOwnerRequests(){return this.config.responsibleForOwnerRequests}registerWithContext(e,...t){return this.context=e,this.register(...t),this.context=null,this}register(...e){let t,i,s,o,n;if(100==++this.registerDepth)throw Error("Unable to autoregister dependency");let a=this.context;for(let r=0,l=e.length;r<l;++r)if(en(t=e[r]))if(W(t))t.register(this,a);else if(void 0!==t.prototype)et.singleton(t,t).register(this);else for(i=Object.keys(t),o=0,n=i.length;o<n;++o)en(s=t[i[o]])&&(W(s)?s.register(this,a):this.register(s));return--this.registerDepth,this}registerResolver(e,t){ei(e);let i=this.resolvers,s=i.get(e);return null==s?i.set(e,t):s instanceof q&&4===s.strategy?s.state.push(t):i.set(e,new q(e,4,[s,t])),t}registerTransformer(e,t){let i=this.getResolver(e);if(null==i)return!1;if(i.getFactory){let e=i.getFactory(this);return null!=e&&(e.registerTransformer(t),!0)}return!1}getResolver(e,t=!0){let i;if(ei(e),void 0!==e.resolve)return e;let s=this;for(;null!=s;){if(null!=(i=s.resolvers.get(e)))return i;if(null==s.parent){let i=X(e)?this:s;return t?this.jitRegister(e,i):null}s=s.parent}return null}has(e,t=!1){return!!this.resolvers.has(e)||!!t&&null!=this.parent&&this.parent.has(e,!0)}get(e){let t;if(ei(e),e.$isResolver)return e.resolve(this,this);let i=this;for(;null!=i;){if(null!=(t=i.resolvers.get(e)))return t.resolve(i,this);if(null==i.parent){let s=X(e)?this:i;return(t=this.jitRegister(e,s)).resolve(i,this)}i=i.parent}throw Error(`Unable to resolve key: ${String(e)}`)}getAll(e,t=!1){let i;ei(e);let s=this;if(t){let t=d.emptyArray;for(;null!=s;)null!=(i=s.resolvers.get(e))&&(t=t.concat(es(i,s,this))),s=s.parent;return t}for(;null!=s;){if(null!=(i=s.resolvers.get(e)))return es(i,s,this);if(null==(s=s.parent))break}return d.emptyArray}getFactory(e){let t=Q.get(e);if(void 0===t){if(ea(e))throw Error(`${e.name} is a native function and therefore cannot be safely constructed by DI. If this is intentional, please use a callback or cachedCallback resolver.`);Q.set(e,t=new _(e,E.getDependencies(e)))}return t}registerFactory(e,t){Q.set(e,t)}createChild(e){return new Z(null,Object.assign({},this.config,e,{parentLocator:()=>this}))}jitRegister(e,t){if("function"!=typeof e)throw Error(`Attempted to jitRegister something that is not a constructor: '${e}'. Did you forget to register this dependency?`);if(G.has(e.name))throw Error(`Attempted to jitRegister an intrinsic type: ${e.name}. Did you forget to add @inject(Key)`);if(W(e)){let i=e.register(t);if(!(i instanceof Object)||null==i.resolve){let i=t.resolvers.get(e);if(void 0!=i)return i;throw Error("A valid resolver was not returned from the static register method")}return i}if(e.$isInterface)throw Error(`Attempted to jitRegister an interface: ${e.friendlyName}`);{let i=this.config.defaultResolver(e,t);return t.resolvers.set(e,i),i}}}let J=new WeakMap;function ee(e){return function(t,i,s){if(J.has(s))return J.get(s);let o=e(t,i,s);return J.set(s,o),o}}let et=Object.freeze({instance:(e,t)=>new q(e,0,t),singleton:(e,t)=>new q(e,1,t),transient:(e,t)=>new q(e,2,t),callback:(e,t)=>new q(e,3,t),cachedCallback:(e,t)=>new q(e,3,ee(t)),aliasTo:(e,t)=>new q(t,5,e)});function ei(e){if(null==e)throw Error("key/value cannot be null or undefined. Are you trying to inject/register something that doesn't exist with DI?")}function es(e,t,i){if(e instanceof q&&4===e.strategy){let s=e.state,o=s.length,n=Array(o);for(;o--;)n[o]=s[o].resolve(t,i);return n}return[e.resolve(t,i)]}let eo="(anonymous)";function en(e){return"object"==typeof e&&null!==e||"function"==typeof e}let ea=(o=new WeakMap,n=!1,a="",r=0,function(e){return void 0===(n=o.get(e))&&(n=(r=(a=e.toString()).length)>=29&&r<=100&&125===a.charCodeAt(r-1)&&32>=a.charCodeAt(r-2)&&93===a.charCodeAt(r-3)&&101===a.charCodeAt(r-4)&&100===a.charCodeAt(r-5)&&111===a.charCodeAt(r-6)&&99===a.charCodeAt(r-7)&&32===a.charCodeAt(r-8)&&101===a.charCodeAt(r-9)&&118===a.charCodeAt(r-10)&&105===a.charCodeAt(r-11)&&116===a.charCodeAt(r-12)&&97===a.charCodeAt(r-13)&&110===a.charCodeAt(r-14)&&88===a.charCodeAt(r-15),o.set(e,n)),n}),er={};function el(e){switch(typeof e){case"number":return e>=0&&(0|e)===e;case"string":{let t=er[e];if(void 0!==t)return t;let i=e.length;if(0===i)return er[e]=!1;let s=0;for(let t=0;t<i;++t)if(s=e.charCodeAt(t),0===t&&48===s&&i>1||s<48||s>57)return er[e]=!1;return er[e]=!0}default:return!1}}function eh(e){return`${e.toLowerCase()}:presentation`}let ed=new Map,ec=Object.freeze({define(e,t,i){let s=eh(e);void 0===ed.get(s)?ed.set(s,t):ed.set(s,!1),i.register(et.instance(s,t))},forTag(e,t){let i=eh(e),s=ed.get(i);return!1===s?E.findResponsibleContainer(t).get(i):s||null}});class eu{constructor(e,t){this.template=e||null,this.styles=void 0===t?null:Array.isArray(t)?d.ElementStyles.create(t):t instanceof d.ElementStyles?t:d.ElementStyles.create([t])}applyTo(e){let t=e.$fastController;null===t.template&&(t.template=this.template),null===t.styles&&(t.styles=this.styles)}}class ep extends d.FASTElement{constructor(){super(...arguments),this._presentation=void 0}get $presentation(){return void 0===this._presentation&&(this._presentation=ec.forTag(this.tagName,this)),this._presentation}templateChanged(){void 0!==this.template&&(this.$fastController.template=this.template)}stylesChanged(){void 0!==this.styles&&(this.$fastController.styles=this.styles)}connectedCallback(){null!==this.$presentation&&this.$presentation.applyTo(this),super.connectedCallback()}static compose(e){return (t={})=>new ev(this===ep?class extends ep{}:this,e,t)}}function em(e,t,i){return"function"==typeof e?e(t,i):e}f([d.observable],ep.prototype,"template",void 0),f([d.observable],ep.prototype,"styles",void 0);class ev{constructor(e,t,i){this.type=e,this.elementDefinition=t,this.overrideDefinition=i,this.definition=Object.assign(Object.assign({},this.elementDefinition),this.overrideDefinition)}register(e,t){let i=this.definition,s=this.overrideDefinition,o=i.prefix||t.elementPrefix,n=`${o}-${i.baseName}`;t.tryDefineElement({name:n,type:this.type,baseClass:this.elementDefinition.baseClass,callback:e=>{let t=new eu(em(i.template,e,i),em(i.styles,e,i));e.definePresentation(t);let o=em(i.shadowOptions,e,i);e.shadowRootMode&&(o?s.shadowOptions||(o.mode=e.shadowRootMode):null!==o&&(o={mode:e.shadowRootMode})),e.defineElement({elementOptions:em(i.elementOptions,e,i),shadowOptions:o,attributes:em(i.attributes,e,i)})}})}}function eb(e,...t){let i=d.AttributeConfiguration.locate(e);t.forEach(t=>{Object.getOwnPropertyNames(t.prototype).forEach(i=>{"constructor"!==i&&Object.defineProperty(e.prototype,i,Object.getOwnPropertyDescriptor(t.prototype,i))}),d.AttributeConfiguration.locate(t).forEach(e=>i.push(e))})}class ef extends ep{constructor(){super(...arguments),this.headinglevel=2,this.expanded=!1,this.clickHandler=e=>{this.expanded=!this.expanded,this.change()},this.change=()=>{this.$emit("change")}}}f([(0,d.attr)({attribute:"heading-level",mode:"fromView",converter:d.nullableNumberConverter})],ef.prototype,"headinglevel",void 0),f([(0,d.attr)({mode:"boolean"})],ef.prototype,"expanded",void 0),f([d.attr],ef.prototype,"id",void 0),eb(ef,c);let eg=(e,t)=>(0,d.html)`
    <template>
        <slot ${(0,d.slotted)({property:"accordionItems",filter:(0,d.elements)()})}></slot>
        <slot name="item" part="item" ${(0,d.slotted)("accordionItems")}></slot>
    </template>
`;var ey=i(74291),eC=i(83021);let ex={single:"single",multi:"multi"};class e$ extends ep{constructor(){super(...arguments),this.expandmode=ex.multi,this.activeItemIndex=0,this.change=()=>{this.$emit("change",this.activeid)},this.setItems=()=>{var e;0!==this.accordionItems.length&&(this.accordionIds=this.getItemIds(),this.accordionItems.forEach((e,t)=>{e instanceof ef&&(e.addEventListener("change",this.activeItemChange),this.isSingleExpandMode()&&(this.activeItemIndex!==t?e.expanded=!1:e.expanded=!0));let i=this.accordionIds[t];e.setAttribute("id","string"!=typeof i?`accordion-${t+1}`:i),this.activeid=this.accordionIds[this.activeItemIndex],e.addEventListener("keydown",this.handleItemKeyDown),e.addEventListener("focus",this.handleItemFocus)}),this.isSingleExpandMode()&&(null!=(e=this.findExpandedItem())?e:this.accordionItems[0]).setAttribute("aria-disabled","true"))},this.removeItemListeners=e=>{e.forEach((e,t)=>{e.removeEventListener("change",this.activeItemChange),e.removeEventListener("keydown",this.handleItemKeyDown),e.removeEventListener("focus",this.handleItemFocus)})},this.activeItemChange=e=>{if(e.defaultPrevented||e.target!==e.currentTarget)return;e.preventDefault();let t=e.target;this.activeid=t.getAttribute("id"),this.isSingleExpandMode()&&(this.resetItems(),t.expanded=!0,t.setAttribute("aria-disabled","true"),this.accordionItems.forEach(e=>{e.hasAttribute("disabled")||e.id===this.activeid||e.removeAttribute("aria-disabled")})),this.activeItemIndex=Array.from(this.accordionItems).indexOf(t),this.change()},this.handleItemKeyDown=e=>{if(e.target===e.currentTarget)switch(this.accordionIds=this.getItemIds(),e.key){case ey.I5:e.preventDefault(),this.adjust(-1);break;case ey.HX:e.preventDefault(),this.adjust(1);break;case"Home":this.activeItemIndex=0,this.focusItem();break;case"End":this.activeItemIndex=this.accordionItems.length-1,this.focusItem()}},this.handleItemFocus=e=>{if(e.target===e.currentTarget){let t=e.target,i=this.activeItemIndex=Array.from(this.accordionItems).indexOf(t);this.activeItemIndex!==i&&-1!==i&&(this.activeItemIndex=i,this.activeid=this.accordionIds[this.activeItemIndex])}}}accordionItemsChanged(e,t){this.$fastController.isConnected&&(this.removeItemListeners(e),this.setItems())}findExpandedItem(){for(let e=0;e<this.accordionItems.length;e++)if("true"===this.accordionItems[e].getAttribute("expanded"))return this.accordionItems[e];return null}resetItems(){this.accordionItems.forEach((e,t)=>{e.expanded=!1})}getItemIds(){return this.accordionItems.map(e=>e.getAttribute("id"))}isSingleExpandMode(){return this.expandmode===ex.single}adjust(e){this.activeItemIndex=(0,eC.Vf)(0,this.accordionItems.length-1,this.activeItemIndex+e),this.focusItem()}focusItem(){let e=this.accordionItems[this.activeItemIndex];e instanceof ef&&e.expandbutton.focus()}}f([(0,d.attr)({attribute:"expand-mode"})],e$.prototype,"expandmode",void 0),f([d.observable],e$.prototype,"accordionItems",void 0);let ew=(e,t)=>(0,d.html)`
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
        ${(0,d.ref)("control")}
    >
        ${p(e,t)}
        <span class="content" part="content">
            <slot ${(0,d.slotted)("defaultSlottedContent")}></slot>
        </span>
        ${u(e,t)}
    </a>
`;class eI{}f([(0,d.attr)({attribute:"aria-atomic"})],eI.prototype,"ariaAtomic",void 0),f([(0,d.attr)({attribute:"aria-busy"})],eI.prototype,"ariaBusy",void 0),f([(0,d.attr)({attribute:"aria-controls"})],eI.prototype,"ariaControls",void 0),f([(0,d.attr)({attribute:"aria-current"})],eI.prototype,"ariaCurrent",void 0),f([(0,d.attr)({attribute:"aria-describedby"})],eI.prototype,"ariaDescribedby",void 0),f([(0,d.attr)({attribute:"aria-details"})],eI.prototype,"ariaDetails",void 0),f([(0,d.attr)({attribute:"aria-disabled"})],eI.prototype,"ariaDisabled",void 0),f([(0,d.attr)({attribute:"aria-errormessage"})],eI.prototype,"ariaErrormessage",void 0),f([(0,d.attr)({attribute:"aria-flowto"})],eI.prototype,"ariaFlowto",void 0),f([(0,d.attr)({attribute:"aria-haspopup"})],eI.prototype,"ariaHaspopup",void 0),f([(0,d.attr)({attribute:"aria-hidden"})],eI.prototype,"ariaHidden",void 0),f([(0,d.attr)({attribute:"aria-invalid"})],eI.prototype,"ariaInvalid",void 0),f([(0,d.attr)({attribute:"aria-keyshortcuts"})],eI.prototype,"ariaKeyshortcuts",void 0),f([(0,d.attr)({attribute:"aria-label"})],eI.prototype,"ariaLabel",void 0),f([(0,d.attr)({attribute:"aria-labelledby"})],eI.prototype,"ariaLabelledby",void 0),f([(0,d.attr)({attribute:"aria-live"})],eI.prototype,"ariaLive",void 0),f([(0,d.attr)({attribute:"aria-owns"})],eI.prototype,"ariaOwns",void 0),f([(0,d.attr)({attribute:"aria-relevant"})],eI.prototype,"ariaRelevant",void 0),f([(0,d.attr)({attribute:"aria-roledescription"})],eI.prototype,"ariaRoledescription",void 0);class ek extends ep{constructor(){super(...arguments),this.handleUnsupportedDelegatesFocus=()=>{var e;window.ShadowRoot&&!window.ShadowRoot.prototype.hasOwnProperty("delegatesFocus")&&(null==(e=this.$fastController.definition.shadowOptions)?void 0:e.delegatesFocus)&&(this.focus=()=>{var e;null==(e=this.control)||e.focus()})}}connectedCallback(){super.connectedCallback(),this.handleUnsupportedDelegatesFocus()}}f([d.attr],ek.prototype,"download",void 0),f([d.attr],ek.prototype,"href",void 0),f([d.attr],ek.prototype,"hreflang",void 0),f([d.attr],ek.prototype,"ping",void 0),f([d.attr],ek.prototype,"referrerpolicy",void 0),f([d.attr],ek.prototype,"rel",void 0),f([d.attr],ek.prototype,"target",void 0),f([d.attr],ek.prototype,"type",void 0),f([d.observable],ek.prototype,"defaultSlottedContent",void 0);class eE{}f([(0,d.attr)({attribute:"aria-expanded"})],eE.prototype,"ariaExpanded",void 0),eb(eE,eI),eb(ek,c,eE);let eO=(e,t)=>(0,d.html)`
    <template class="${e=>e.initialLayoutComplete?"loaded":""}">
        ${(0,d.when)(e=>e.initialLayoutComplete,(0,d.html)`
                <slot></slot>
            `)}
    </template>
`;var eT=i(30086);let eR=e=>{let t=e.closest("[dir]");return null!==t&&"rtl"===t.dir?eT.O.rtl:eT.O.ltr};class eD extends ep{constructor(){super(...arguments),this.anchor="",this.viewport="",this.horizontalPositioningMode="uncontrolled",this.horizontalDefaultPosition="unset",this.horizontalViewportLock=!1,this.horizontalInset=!1,this.horizontalScaling="content",this.verticalPositioningMode="uncontrolled",this.verticalDefaultPosition="unset",this.verticalViewportLock=!1,this.verticalInset=!1,this.verticalScaling="content",this.fixedPlacement=!1,this.autoUpdateMode="anchor",this.anchorElement=null,this.viewportElement=null,this.initialLayoutComplete=!1,this.resizeDetector=null,this.baseHorizontalOffset=0,this.baseVerticalOffset=0,this.pendingPositioningUpdate=!1,this.pendingReset=!1,this.currentDirection=eT.O.ltr,this.regionVisible=!1,this.forceUpdate=!1,this.updateThreshold=.5,this.update=()=>{this.pendingPositioningUpdate||this.requestPositionUpdates()},this.startObservers=()=>{this.stopObservers(),null!==this.anchorElement&&(this.requestPositionUpdates(),null!==this.resizeDetector&&(this.resizeDetector.observe(this.anchorElement),this.resizeDetector.observe(this)))},this.requestPositionUpdates=()=>{null===this.anchorElement||this.pendingPositioningUpdate||(eD.intersectionService.requestPosition(this,this.handleIntersection),eD.intersectionService.requestPosition(this.anchorElement,this.handleIntersection),null!==this.viewportElement&&eD.intersectionService.requestPosition(this.viewportElement,this.handleIntersection),this.pendingPositioningUpdate=!0)},this.stopObservers=()=>{this.pendingPositioningUpdate&&(this.pendingPositioningUpdate=!1,eD.intersectionService.cancelRequestPosition(this,this.handleIntersection),null!==this.anchorElement&&eD.intersectionService.cancelRequestPosition(this.anchorElement,this.handleIntersection),null!==this.viewportElement&&eD.intersectionService.cancelRequestPosition(this.viewportElement,this.handleIntersection)),null!==this.resizeDetector&&this.resizeDetector.disconnect()},this.getViewport=()=>"string"!=typeof this.viewport||""===this.viewport?document.documentElement:document.getElementById(this.viewport),this.getAnchor=()=>document.getElementById(this.anchor),this.handleIntersection=e=>{!this.pendingPositioningUpdate||(this.pendingPositioningUpdate=!1,this.applyIntersectionEntries(e)&&this.updateLayout())},this.applyIntersectionEntries=e=>{let t=e.find(e=>e.target===this),i=e.find(e=>e.target===this.anchorElement),s=e.find(e=>e.target===this.viewportElement);return void 0!==t&&void 0!==s&&void 0!==i&&!!(!this.regionVisible||this.forceUpdate||void 0===this.regionRect||void 0===this.anchorRect||void 0===this.viewportRect||this.isRectDifferent(this.anchorRect,i.boundingClientRect)||this.isRectDifferent(this.viewportRect,s.boundingClientRect)||this.isRectDifferent(this.regionRect,t.boundingClientRect))&&(this.regionRect=t.boundingClientRect,this.anchorRect=i.boundingClientRect,this.viewportElement===document.documentElement?this.viewportRect=new DOMRectReadOnly(s.boundingClientRect.x+document.documentElement.scrollLeft,s.boundingClientRect.y+document.documentElement.scrollTop,s.boundingClientRect.width,s.boundingClientRect.height):this.viewportRect=s.boundingClientRect,this.updateRegionOffset(),this.forceUpdate=!1,!0)},this.updateRegionOffset=()=>{this.anchorRect&&this.regionRect&&(this.baseHorizontalOffset=this.baseHorizontalOffset+(this.anchorRect.left-this.regionRect.left)+(this.translateX-this.baseHorizontalOffset),this.baseVerticalOffset=this.baseVerticalOffset+(this.anchorRect.top-this.regionRect.top)+(this.translateY-this.baseVerticalOffset))},this.isRectDifferent=(e,t)=>!!(Math.abs(e.top-t.top)>this.updateThreshold||Math.abs(e.right-t.right)>this.updateThreshold||Math.abs(e.bottom-t.bottom)>this.updateThreshold||Math.abs(e.left-t.left)>this.updateThreshold),this.handleResize=e=>{this.update()},this.reset=()=>{this.pendingReset&&(this.pendingReset=!1,null===this.anchorElement&&(this.anchorElement=this.getAnchor()),null===this.viewportElement&&(this.viewportElement=this.getViewport()),this.currentDirection=eR(this),this.startObservers())},this.updateLayout=()=>{let e,t;if("uncontrolled"!==this.horizontalPositioningMode){let e=this.getPositioningOptions(this.horizontalInset);if("center"===this.horizontalDefaultPosition)t="center";else if("unset"!==this.horizontalDefaultPosition){let e=this.horizontalDefaultPosition;if("start"===e||"end"===e){let t=eR(this);if(t!==this.currentDirection){this.currentDirection=t,this.initialize();return}e=this.currentDirection===eT.O.ltr?"start"===e?"left":"right":"start"===e?"right":"left"}switch(e){case"left":t=this.horizontalInset?"insetStart":"start";break;case"right":t=this.horizontalInset?"insetEnd":"end"}}let i=void 0!==this.horizontalThreshold?this.horizontalThreshold:void 0!==this.regionRect?this.regionRect.width:0,s=void 0!==this.anchorRect?this.anchorRect.left:0,o=void 0!==this.anchorRect?this.anchorRect.right:0,n=void 0!==this.anchorRect?this.anchorRect.width:0,a=void 0!==this.viewportRect?this.viewportRect.left:0,r=void 0!==this.viewportRect?this.viewportRect.right:0;(void 0===t||"locktodefault"!==this.horizontalPositioningMode&&this.getAvailableSpace(t,s,o,n,a,r)<i)&&(t=this.getAvailableSpace(e[0],s,o,n,a,r)>this.getAvailableSpace(e[1],s,o,n,a,r)?e[0]:e[1])}if("uncontrolled"!==this.verticalPositioningMode){let t=this.getPositioningOptions(this.verticalInset);if("center"===this.verticalDefaultPosition)e="center";else if("unset"!==this.verticalDefaultPosition)switch(this.verticalDefaultPosition){case"top":e=this.verticalInset?"insetStart":"start";break;case"bottom":e=this.verticalInset?"insetEnd":"end"}let i=void 0!==this.verticalThreshold?this.verticalThreshold:void 0!==this.regionRect?this.regionRect.height:0,s=void 0!==this.anchorRect?this.anchorRect.top:0,o=void 0!==this.anchorRect?this.anchorRect.bottom:0,n=void 0!==this.anchorRect?this.anchorRect.height:0,a=void 0!==this.viewportRect?this.viewportRect.top:0,r=void 0!==this.viewportRect?this.viewportRect.bottom:0;(void 0===e||"locktodefault"!==this.verticalPositioningMode&&this.getAvailableSpace(e,s,o,n,a,r)<i)&&(e=this.getAvailableSpace(t[0],s,o,n,a,r)>this.getAvailableSpace(t[1],s,o,n,a,r)?t[0]:t[1])}let i=this.getNextRegionDimension(t,e),s=this.horizontalPosition!==t||this.verticalPosition!==e;if(this.setHorizontalPosition(t,i),this.setVerticalPosition(e,i),this.updateRegionStyle(),!this.initialLayoutComplete){this.initialLayoutComplete=!0,this.requestPositionUpdates();return}this.regionVisible||(this.regionVisible=!0,this.style.removeProperty("pointer-events"),this.style.removeProperty("opacity"),this.classList.toggle("loaded",!0),this.$emit("loaded",this,{bubbles:!1})),this.updatePositionClasses(),s&&this.$emit("positionchange",this,{bubbles:!1})},this.updateRegionStyle=()=>{this.style.width=this.regionWidth,this.style.height=this.regionHeight,this.style.transform=`translate(${this.translateX}px, ${this.translateY}px)`},this.updatePositionClasses=()=>{this.classList.toggle("top","start"===this.verticalPosition),this.classList.toggle("bottom","end"===this.verticalPosition),this.classList.toggle("inset-top","insetStart"===this.verticalPosition),this.classList.toggle("inset-bottom","insetEnd"===this.verticalPosition),this.classList.toggle("vertical-center","center"===this.verticalPosition),this.classList.toggle("left","start"===this.horizontalPosition),this.classList.toggle("right","end"===this.horizontalPosition),this.classList.toggle("inset-left","insetStart"===this.horizontalPosition),this.classList.toggle("inset-right","insetEnd"===this.horizontalPosition),this.classList.toggle("horizontal-center","center"===this.horizontalPosition)},this.setHorizontalPosition=(e,t)=>{if(void 0===e||void 0===this.regionRect||void 0===this.anchorRect||void 0===this.viewportRect)return;let i=0;switch(this.horizontalScaling){case"anchor":case"fill":i=this.horizontalViewportLock?this.viewportRect.width:t.width,this.regionWidth=`${i}px`;break;case"content":i=this.regionRect.width,this.regionWidth="unset"}let s=0;switch(e){case"start":this.translateX=this.baseHorizontalOffset-i,this.horizontalViewportLock&&this.anchorRect.left>this.viewportRect.right&&(this.translateX=this.translateX-(this.anchorRect.left-this.viewportRect.right));break;case"insetStart":this.translateX=this.baseHorizontalOffset-i+this.anchorRect.width,this.horizontalViewportLock&&this.anchorRect.right>this.viewportRect.right&&(this.translateX=this.translateX-(this.anchorRect.right-this.viewportRect.right));break;case"insetEnd":this.translateX=this.baseHorizontalOffset,this.horizontalViewportLock&&this.anchorRect.left<this.viewportRect.left&&(this.translateX=this.translateX-(this.anchorRect.left-this.viewportRect.left));break;case"end":this.translateX=this.baseHorizontalOffset+this.anchorRect.width,this.horizontalViewportLock&&this.anchorRect.right<this.viewportRect.left&&(this.translateX=this.translateX-(this.anchorRect.right-this.viewportRect.left));break;case"center":if(s=(this.anchorRect.width-i)/2,this.translateX=this.baseHorizontalOffset+s,this.horizontalViewportLock){let e=this.anchorRect.left+s,t=this.anchorRect.right-s;e<this.viewportRect.left&&!(t>this.viewportRect.right)?this.translateX=this.translateX-(e-this.viewportRect.left):t>this.viewportRect.right&&!(e<this.viewportRect.left)&&(this.translateX=this.translateX-(t-this.viewportRect.right))}}this.horizontalPosition=e},this.setVerticalPosition=(e,t)=>{if(void 0===e||void 0===this.regionRect||void 0===this.anchorRect||void 0===this.viewportRect)return;let i=0;switch(this.verticalScaling){case"anchor":case"fill":i=this.verticalViewportLock?this.viewportRect.height:t.height,this.regionHeight=`${i}px`;break;case"content":i=this.regionRect.height,this.regionHeight="unset"}let s=0;switch(e){case"start":this.translateY=this.baseVerticalOffset-i,this.verticalViewportLock&&this.anchorRect.top>this.viewportRect.bottom&&(this.translateY=this.translateY-(this.anchorRect.top-this.viewportRect.bottom));break;case"insetStart":this.translateY=this.baseVerticalOffset-i+this.anchorRect.height,this.verticalViewportLock&&this.anchorRect.bottom>this.viewportRect.bottom&&(this.translateY=this.translateY-(this.anchorRect.bottom-this.viewportRect.bottom));break;case"insetEnd":this.translateY=this.baseVerticalOffset,this.verticalViewportLock&&this.anchorRect.top<this.viewportRect.top&&(this.translateY=this.translateY-(this.anchorRect.top-this.viewportRect.top));break;case"end":this.translateY=this.baseVerticalOffset+this.anchorRect.height,this.verticalViewportLock&&this.anchorRect.bottom<this.viewportRect.top&&(this.translateY=this.translateY-(this.anchorRect.bottom-this.viewportRect.top));break;case"center":if(s=(this.anchorRect.height-i)/2,this.translateY=this.baseVerticalOffset+s,this.verticalViewportLock){let e=this.anchorRect.top+s,t=this.anchorRect.bottom-s;e<this.viewportRect.top&&!(t>this.viewportRect.bottom)?this.translateY=this.translateY-(e-this.viewportRect.top):t>this.viewportRect.bottom&&!(e<this.viewportRect.top)&&(this.translateY=this.translateY-(t-this.viewportRect.bottom))}}this.verticalPosition=e},this.getPositioningOptions=e=>e?["insetStart","insetEnd"]:["start","end"],this.getAvailableSpace=(e,t,i,s,o,n)=>{let a=t-o,r=n-(t+s);switch(e){case"start":return a;case"insetStart":return a+s;case"insetEnd":return r+s;case"end":return r;case"center":return 2*Math.min(a,r)+s}},this.getNextRegionDimension=(e,t)=>{let i={height:void 0!==this.regionRect?this.regionRect.height:0,width:void 0!==this.regionRect?this.regionRect.width:0};return void 0!==e&&"fill"===this.horizontalScaling?i.width=this.getAvailableSpace(e,void 0!==this.anchorRect?this.anchorRect.left:0,void 0!==this.anchorRect?this.anchorRect.right:0,void 0!==this.anchorRect?this.anchorRect.width:0,void 0!==this.viewportRect?this.viewportRect.left:0,void 0!==this.viewportRect?this.viewportRect.right:0):"anchor"===this.horizontalScaling&&(i.width=void 0!==this.anchorRect?this.anchorRect.width:0),void 0!==t&&"fill"===this.verticalScaling?i.height=this.getAvailableSpace(t,void 0!==this.anchorRect?this.anchorRect.top:0,void 0!==this.anchorRect?this.anchorRect.bottom:0,void 0!==this.anchorRect?this.anchorRect.height:0,void 0!==this.viewportRect?this.viewportRect.top:0,void 0!==this.viewportRect?this.viewportRect.bottom:0):"anchor"===this.verticalScaling&&(i.height=void 0!==this.anchorRect?this.anchorRect.height:0),i},this.startAutoUpdateEventListeners=()=>{window.addEventListener("resize",this.update,{passive:!0}),window.addEventListener("scroll",this.update,{passive:!0,capture:!0}),null!==this.resizeDetector&&null!==this.viewportElement&&this.resizeDetector.observe(this.viewportElement)},this.stopAutoUpdateEventListeners=()=>{window.removeEventListener("resize",this.update),window.removeEventListener("scroll",this.update),null!==this.resizeDetector&&null!==this.viewportElement&&this.resizeDetector.unobserve(this.viewportElement)}}anchorChanged(){this.initialLayoutComplete&&(this.anchorElement=this.getAnchor())}viewportChanged(){this.initialLayoutComplete&&(this.viewportElement=this.getViewport())}horizontalPositioningModeChanged(){this.requestReset()}horizontalDefaultPositionChanged(){this.updateForAttributeChange()}horizontalViewportLockChanged(){this.updateForAttributeChange()}horizontalInsetChanged(){this.updateForAttributeChange()}horizontalThresholdChanged(){this.updateForAttributeChange()}horizontalScalingChanged(){this.updateForAttributeChange()}verticalPositioningModeChanged(){this.requestReset()}verticalDefaultPositionChanged(){this.updateForAttributeChange()}verticalViewportLockChanged(){this.updateForAttributeChange()}verticalInsetChanged(){this.updateForAttributeChange()}verticalThresholdChanged(){this.updateForAttributeChange()}verticalScalingChanged(){this.updateForAttributeChange()}fixedPlacementChanged(){this.$fastController.isConnected&&this.initialLayoutComplete&&this.initialize()}autoUpdateModeChanged(e,t){this.$fastController.isConnected&&this.initialLayoutComplete&&("auto"===e&&this.stopAutoUpdateEventListeners(),"auto"===t&&this.startAutoUpdateEventListeners())}anchorElementChanged(){this.requestReset()}viewportElementChanged(){this.$fastController.isConnected&&this.initialLayoutComplete&&this.initialize()}connectedCallback(){super.connectedCallback(),"auto"===this.autoUpdateMode&&this.startAutoUpdateEventListeners(),this.initialize()}disconnectedCallback(){super.disconnectedCallback(),"auto"===this.autoUpdateMode&&this.stopAutoUpdateEventListeners(),this.stopObservers(),this.disconnectResizeDetector()}adoptedCallback(){this.initialize()}disconnectResizeDetector(){null!==this.resizeDetector&&(this.resizeDetector.disconnect(),this.resizeDetector=null)}initializeResizeDetector(){this.disconnectResizeDetector(),this.resizeDetector=new window.ResizeObserver(this.handleResize)}updateForAttributeChange(){this.$fastController.isConnected&&this.initialLayoutComplete&&(this.forceUpdate=!0,this.update())}initialize(){this.initializeResizeDetector(),null===this.anchorElement&&(this.anchorElement=this.getAnchor()),this.requestReset()}requestReset(){this.$fastController.isConnected&&!1===this.pendingReset&&(this.setInitialState(),d.DOM.queueUpdate(()=>this.reset()),this.pendingReset=!0)}setInitialState(){this.initialLayoutComplete=!1,this.regionVisible=!1,this.translateX=0,this.translateY=0,this.baseHorizontalOffset=0,this.baseVerticalOffset=0,this.viewportRect=void 0,this.regionRect=void 0,this.anchorRect=void 0,this.verticalPosition=void 0,this.horizontalPosition=void 0,this.style.opacity="0",this.style.pointerEvents="none",this.forceUpdate=!1,this.style.position=this.fixedPlacement?"fixed":"absolute",this.updatePositionClasses(),this.updateRegionStyle()}}eD.intersectionService=new class{constructor(){this.intersectionDetector=null,this.observedElements=new Map,this.requestPosition=(e,t)=>{var i;if(null!==this.intersectionDetector){if(this.observedElements.has(e)){null==(i=this.observedElements.get(e))||i.push(t);return}this.observedElements.set(e,[t]),this.intersectionDetector.observe(e)}},this.cancelRequestPosition=(e,t)=>{let i=this.observedElements.get(e);if(void 0!==i){let e=i.indexOf(t);-1!==e&&i.splice(e,1)}},this.initializeIntersectionDetector=()=>{d.$global.IntersectionObserver&&(this.intersectionDetector=new IntersectionObserver(this.handleIntersection,{root:null,rootMargin:"0px",threshold:[0,1]}))},this.handleIntersection=e=>{if(null===this.intersectionDetector)return;let t=[],i=[];e.forEach(e=>{var s;null==(s=this.intersectionDetector)||s.unobserve(e.target);let o=this.observedElements.get(e.target);void 0!==o&&(o.forEach(s=>{let o=t.indexOf(s);-1===o&&(o=t.length,t.push(s),i.push([])),i[o].push(e)}),this.observedElements.delete(e.target))}),t.forEach((e,t)=>{e(i[t])})},this.initializeIntersectionDetector()}},f([d.attr],eD.prototype,"anchor",void 0),f([d.attr],eD.prototype,"viewport",void 0),f([(0,d.attr)({attribute:"horizontal-positioning-mode"})],eD.prototype,"horizontalPositioningMode",void 0),f([(0,d.attr)({attribute:"horizontal-default-position"})],eD.prototype,"horizontalDefaultPosition",void 0),f([(0,d.attr)({attribute:"horizontal-viewport-lock",mode:"boolean"})],eD.prototype,"horizontalViewportLock",void 0),f([(0,d.attr)({attribute:"horizontal-inset",mode:"boolean"})],eD.prototype,"horizontalInset",void 0),f([(0,d.attr)({attribute:"horizontal-threshold"})],eD.prototype,"horizontalThreshold",void 0),f([(0,d.attr)({attribute:"horizontal-scaling"})],eD.prototype,"horizontalScaling",void 0),f([(0,d.attr)({attribute:"vertical-positioning-mode"})],eD.prototype,"verticalPositioningMode",void 0),f([(0,d.attr)({attribute:"vertical-default-position"})],eD.prototype,"verticalDefaultPosition",void 0),f([(0,d.attr)({attribute:"vertical-viewport-lock",mode:"boolean"})],eD.prototype,"verticalViewportLock",void 0),f([(0,d.attr)({attribute:"vertical-inset",mode:"boolean"})],eD.prototype,"verticalInset",void 0),f([(0,d.attr)({attribute:"vertical-threshold"})],eD.prototype,"verticalThreshold",void 0),f([(0,d.attr)({attribute:"vertical-scaling"})],eD.prototype,"verticalScaling",void 0),f([(0,d.attr)({attribute:"fixed-placement",mode:"boolean"})],eD.prototype,"fixedPlacement",void 0),f([(0,d.attr)({attribute:"auto-update-mode"})],eD.prototype,"autoUpdateMode",void 0),f([d.observable],eD.prototype,"anchorElement",void 0),f([d.observable],eD.prototype,"viewportElement",void 0),f([d.observable],eD.prototype,"initialLayoutComplete",void 0);let eS={horizontalDefaultPosition:"center",horizontalPositioningMode:"locktodefault",horizontalInset:!1,horizontalScaling:"anchor"},eA=Object.assign(Object.assign({},eS),{verticalDefaultPosition:"top",verticalPositioningMode:"locktodefault",verticalInset:!1,verticalScaling:"content"}),eF=Object.assign(Object.assign({},eS),{verticalDefaultPosition:"bottom",verticalPositioningMode:"locktodefault",verticalInset:!1,verticalScaling:"content"}),eL=Object.assign(Object.assign({},eS),{verticalPositioningMode:"dynamic",verticalInset:!1,verticalScaling:"content"}),eM=Object.assign(Object.assign({},eA),{verticalScaling:"fill"}),eP=Object.assign(Object.assign({},eF),{verticalScaling:"fill"}),eH=Object.assign(Object.assign({},eL),{verticalScaling:"fill"}),ez=(e,t)=>(0,d.html)`
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
`;class eV extends ep{connectedCallback(){super.connectedCallback(),this.shape||(this.shape="circle")}}f([d.attr],eV.prototype,"fill",void 0),f([d.attr],eV.prototype,"color",void 0),f([d.attr],eV.prototype,"link",void 0),f([d.attr],eV.prototype,"shape",void 0);let eN=(e,t)=>(0,d.html)`
    <template class="${e=>e.circular?"circular":""}">
        <div class="control" part="control" style="${e=>e.generateBadgeStyle()}">
            <slot></slot>
        </div>
    </template>
`;class eB extends ep{constructor(){super(...arguments),this.generateBadgeStyle=()=>{if(!this.fill&&!this.color)return;let e=`background-color: var(--badge-fill-${this.fill});`,t=`color: var(--badge-color-${this.color});`;return this.fill&&!this.color?e:this.color&&!this.fill?t:`${t} ${e}`}}}f([(0,d.attr)({attribute:"fill"})],eB.prototype,"fill",void 0),f([(0,d.attr)({attribute:"color"})],eB.prototype,"color",void 0),f([(0,d.attr)({mode:"boolean"})],eB.prototype,"circular",void 0);let eq=(e,t)=>(0,d.html)`
    <div role="listitem" class="listitem" part="listitem">
        ${(0,d.when)(e=>e.href&&e.href.length>0,(0,d.html)`
                ${ew(e,t)}
            `)}
        ${(0,d.when)(e=>!e.href,(0,d.html)`
                ${p(e,t)}
                <slot></slot>
                ${u(e,t)}
            `)}
        ${(0,d.when)(e=>e.separator,(0,d.html)`
                <span class="separator" part="separator" aria-hidden="true">
                    <slot name="separator">${t.separator||""}</slot>
                </span>
            `)}
    </div>
`;class eU extends ek{constructor(){super(...arguments),this.separator=!0}}f([d.observable],eU.prototype,"separator",void 0),eb(eU,c,eE);let ej=(e,t)=>(0,d.html)`
    <template role="navigation">
        <div role="list" class="list" part="list">
            <slot
                ${(0,d.slotted)({property:"slottedBreadcrumbItems",filter:(0,d.elements)()})}
            ></slot>
        </div>
    </template>
`;class e_ extends ep{slottedBreadcrumbItemsChanged(){if(this.$fastController.isConnected){if(void 0===this.slottedBreadcrumbItems||0===this.slottedBreadcrumbItems.length)return;let e=this.slottedBreadcrumbItems[this.slottedBreadcrumbItems.length-1];this.slottedBreadcrumbItems.forEach(t=>{let i=t===e;this.setItemSeparator(t,i),this.setAriaCurrent(t,i)})}}setItemSeparator(e,t){e instanceof eU&&(e.separator=!t)}findChildWithHref(e){var t,i;return e.childElementCount>0?e.querySelector("a[href]"):(null==(t=e.shadowRoot)?void 0:t.childElementCount)?null==(i=e.shadowRoot)?void 0:i.querySelector("a[href]"):null}setAriaCurrent(e,t){let i=this.findChildWithHref(e);null===i&&e.hasAttribute("href")&&e instanceof eU?t?e.setAttribute("aria-current","page"):e.removeAttribute("aria-current"):null!==i&&(t?i.setAttribute("aria-current","page"):i.removeAttribute("aria-current"))}}f([d.observable],e_.prototype,"slottedBreadcrumbItems",void 0);let eK=(e,t)=>(0,d.html)`
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
        ${(0,d.ref)("control")}
    >
        ${p(e,t)}
        <span class="content" part="content">
            <slot ${(0,d.slotted)("defaultSlottedContent")}></slot>
        </span>
        ${u(e,t)}
    </button>
`,eW="form-associated-proxy",eX="ElementInternals",eG=eX in window&&"setFormValue"in window[eX].prototype,eY=new WeakMap;function eQ(e){let t=class extends e{constructor(...e){super(...e),this.dirtyValue=!1,this.disabled=!1,this.proxyEventsToBlock=["change","click"],this.proxyInitialized=!1,this.required=!1,this.initialValue=this.initialValue||"",this.elementInternals||(this.formResetCallback=this.formResetCallback.bind(this))}static get formAssociated(){return eG}get validity(){return this.elementInternals?this.elementInternals.validity:this.proxy.validity}get form(){return this.elementInternals?this.elementInternals.form:this.proxy.form}get validationMessage(){return this.elementInternals?this.elementInternals.validationMessage:this.proxy.validationMessage}get willValidate(){return this.elementInternals?this.elementInternals.willValidate:this.proxy.willValidate}get labels(){if(this.elementInternals)return Object.freeze(Array.from(this.elementInternals.labels));if(!(this.proxy instanceof HTMLElement)||!this.proxy.ownerDocument||!this.id)return d.emptyArray;{let e=this.proxy.labels,t=Array.from(this.proxy.getRootNode().querySelectorAll(`[for='${this.id}']`));return Object.freeze(e?t.concat(Array.from(e)):t)}}valueChanged(e,t){this.dirtyValue=!0,this.proxy instanceof HTMLElement&&(this.proxy.value=this.value),this.currentValue=this.value,this.setFormValue(this.value),this.validate()}currentValueChanged(){this.value=this.currentValue}initialValueChanged(e,t){this.dirtyValue||(this.value=this.initialValue,this.dirtyValue=!1)}disabledChanged(e,t){this.proxy instanceof HTMLElement&&(this.proxy.disabled=this.disabled),d.DOM.queueUpdate(()=>this.classList.toggle("disabled",this.disabled))}nameChanged(e,t){this.proxy instanceof HTMLElement&&(this.proxy.name=this.name)}requiredChanged(e,t){this.proxy instanceof HTMLElement&&(this.proxy.required=this.required),d.DOM.queueUpdate(()=>this.classList.toggle("required",this.required)),this.validate()}get elementInternals(){if(!eG)return null;let e=eY.get(this);return e||(e=this.attachInternals(),eY.set(this,e)),e}connectedCallback(){super.connectedCallback(),this.addEventListener("keypress",this._keypressHandler),this.value||(this.value=this.initialValue,this.dirtyValue=!1),!this.elementInternals&&(this.attachProxy(),this.form&&this.form.addEventListener("reset",this.formResetCallback))}disconnectedCallback(){super.disconnectedCallback(),this.proxyEventsToBlock.forEach(e=>this.proxy.removeEventListener(e,this.stopPropagation)),!this.elementInternals&&this.form&&this.form.removeEventListener("reset",this.formResetCallback)}checkValidity(){return this.elementInternals?this.elementInternals.checkValidity():this.proxy.checkValidity()}reportValidity(){return this.elementInternals?this.elementInternals.reportValidity():this.proxy.reportValidity()}setValidity(e,t,i){this.elementInternals?this.elementInternals.setValidity(e,t,i):"string"==typeof t&&this.proxy.setCustomValidity(t)}formDisabledCallback(e){this.disabled=e}formResetCallback(){this.value=this.initialValue,this.dirtyValue=!1}attachProxy(){var e;this.proxyInitialized||(this.proxyInitialized=!0,this.proxy.style.display="none",this.proxyEventsToBlock.forEach(e=>this.proxy.addEventListener(e,this.stopPropagation)),this.proxy.disabled=this.disabled,this.proxy.required=this.required,"string"==typeof this.name&&(this.proxy.name=this.name),"string"==typeof this.value&&(this.proxy.value=this.value),this.proxy.setAttribute("slot",eW),this.proxySlot=document.createElement("slot"),this.proxySlot.setAttribute("name",eW)),null==(e=this.shadowRoot)||e.appendChild(this.proxySlot),this.appendChild(this.proxy)}detachProxy(){var e;this.removeChild(this.proxy),null==(e=this.shadowRoot)||e.removeChild(this.proxySlot)}validate(e){this.proxy instanceof HTMLElement&&this.setValidity(this.proxy.validity,this.proxy.validationMessage,e)}setFormValue(e,t){this.elementInternals&&this.elementInternals.setFormValue(e,t||e)}_keypressHandler(e){if("Enter"===e.key&&this.form instanceof HTMLFormElement){let e=this.form.querySelector("[type=submit]");null==e||e.click()}}stopPropagation(e){e.stopPropagation()}};return(0,d.attr)({mode:"boolean"})(t.prototype,"disabled"),(0,d.attr)({mode:"fromView",attribute:"value"})(t.prototype,"initialValue"),(0,d.attr)({attribute:"current-value"})(t.prototype,"currentValue"),(0,d.attr)(t.prototype,"name"),(0,d.attr)({mode:"boolean"})(t.prototype,"required"),(0,d.observable)(t.prototype,"value"),t}function eZ(e){class t extends eQ(e){}class i extends t{constructor(...e){super(e),this.dirtyChecked=!1,this.checkedAttribute=!1,this.checked=!1,this.dirtyChecked=!1}checkedAttributeChanged(){this.defaultChecked=this.checkedAttribute}defaultCheckedChanged(){this.dirtyChecked||(this.checked=this.defaultChecked,this.dirtyChecked=!1)}checkedChanged(e,t){this.dirtyChecked||(this.dirtyChecked=!0),this.currentChecked=this.checked,this.updateForm(),this.proxy instanceof HTMLInputElement&&(this.proxy.checked=this.checked),void 0!==e&&this.$emit("change"),this.validate()}currentCheckedChanged(e,t){this.checked=this.currentChecked}updateForm(){let e=this.checked?this.value:null;this.setFormValue(e,e)}connectedCallback(){super.connectedCallback(),this.updateForm()}formResetCallback(){super.formResetCallback(),this.checked=!!this.checkedAttribute,this.dirtyChecked=!1}}return(0,d.attr)({attribute:"checked",mode:"boolean"})(i.prototype,"checkedAttribute"),(0,d.attr)({attribute:"current-checked",converter:d.booleanConverter})(i.prototype,"currentChecked"),(0,d.observable)(i.prototype,"defaultChecked"),(0,d.observable)(i.prototype,"checked"),i}class eJ extends ep{}class e0 extends eQ(eJ){constructor(){super(...arguments),this.proxy=document.createElement("input")}}class e1 extends e0{constructor(){super(...arguments),this.handleClick=e=>{var t;this.disabled&&(null==(t=this.defaultSlottedContent)?void 0:t.length)<=1&&e.stopPropagation()},this.handleSubmission=()=>{if(!this.form)return;let e=this.proxy.isConnected;e||this.attachProxy(),"function"==typeof this.form.requestSubmit?this.form.requestSubmit(this.proxy):this.proxy.click(),e||this.detachProxy()},this.handleFormReset=()=>{var e;null==(e=this.form)||e.reset()},this.handleUnsupportedDelegatesFocus=()=>{var e;window.ShadowRoot&&!window.ShadowRoot.prototype.hasOwnProperty("delegatesFocus")&&(null==(e=this.$fastController.definition.shadowOptions)?void 0:e.delegatesFocus)&&(this.focus=()=>{this.control.focus()})}}formactionChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.formAction=this.formaction)}formenctypeChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.formEnctype=this.formenctype)}formmethodChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.formMethod=this.formmethod)}formnovalidateChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.formNoValidate=this.formnovalidate)}formtargetChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.formTarget=this.formtarget)}typeChanged(e,t){this.proxy instanceof HTMLInputElement&&(this.proxy.type=this.type),"submit"===t&&this.addEventListener("click",this.handleSubmission),"submit"===e&&this.removeEventListener("click",this.handleSubmission),"reset"===t&&this.addEventListener("click",this.handleFormReset),"reset"===e&&this.removeEventListener("click",this.handleFormReset)}validate(){super.validate(this.control)}connectedCallback(){var e;super.connectedCallback(),this.proxy.setAttribute("type",this.type),this.handleUnsupportedDelegatesFocus();let t=Array.from(null==(e=this.control)?void 0:e.children);t&&t.forEach(e=>{e.addEventListener("click",this.handleClick)})}disconnectedCallback(){var e;super.disconnectedCallback();let t=Array.from(null==(e=this.control)?void 0:e.children);t&&t.forEach(e=>{e.removeEventListener("click",this.handleClick)})}}f([(0,d.attr)({mode:"boolean"})],e1.prototype,"autofocus",void 0),f([(0,d.attr)({attribute:"form"})],e1.prototype,"formId",void 0),f([d.attr],e1.prototype,"formaction",void 0),f([d.attr],e1.prototype,"formenctype",void 0),f([d.attr],e1.prototype,"formmethod",void 0),f([(0,d.attr)({mode:"boolean"})],e1.prototype,"formnovalidate",void 0),f([d.attr],e1.prototype,"formtarget",void 0),f([d.attr],e1.prototype,"type",void 0),f([d.observable],e1.prototype,"defaultSlottedContent",void 0);class e5{}f([(0,d.attr)({attribute:"aria-expanded"})],e5.prototype,"ariaExpanded",void 0),f([(0,d.attr)({attribute:"aria-pressed"})],e5.prototype,"ariaPressed",void 0),eb(e5,eI),eb(e1,c,e5);class e4{constructor(e){if(this.dayFormat="numeric",this.weekdayFormat="long",this.monthFormat="long",this.yearFormat="numeric",this.date=new Date,e)for(const t in e){const i=e[t];"date"===t?this.date=this.getDateObject(i):this[t]=i}}getDateObject(e){if("string"==typeof e){let t=e.split(/[/-]/);return t.length<3?new Date:new Date(parseInt(t[2],10),parseInt(t[0],10)-1,parseInt(t[1],10))}if("day"in e&&"month"in e&&"year"in e){let{day:t,month:i,year:s}=e;return new Date(s,i-1,t)}return e}getDate(e=this.date,t={weekday:this.weekdayFormat,month:this.monthFormat,day:this.dayFormat,year:this.yearFormat},i=this.locale){let s=this.getDateObject(e);if(!s.getTime())return"";let o=Object.assign({timeZone:Intl.DateTimeFormat().resolvedOptions().timeZone},t);return new Intl.DateTimeFormat(i,o).format(s)}getDay(e=this.date.getDate(),t=this.dayFormat,i=this.locale){return this.getDate({month:1,day:e,year:2020},{day:t},i)}getMonth(e=this.date.getMonth()+1,t=this.monthFormat,i=this.locale){return this.getDate({month:e,day:2,year:2020},{month:t},i)}getYear(e=this.date.getFullYear(),t=this.yearFormat,i=this.locale){return this.getDate({month:2,day:2,year:e},{year:t},i)}getWeekday(e=0,t=this.weekdayFormat,i=this.locale){let s=`1-${e+1}-2017`;return this.getDate(s,{weekday:t},i)}getWeekdays(e=this.weekdayFormat,t=this.locale){return Array(7).fill(null).map((i,s)=>this.getWeekday(s,e,t))}}class e8 extends ep{constructor(){super(...arguments),this.dateFormatter=new e4,this.readonly=!1,this.locale="en-US",this.month=new Date().getMonth()+1,this.year=new Date().getFullYear(),this.dayFormat="numeric",this.weekdayFormat="short",this.monthFormat="long",this.yearFormat="numeric",this.minWeeks=0,this.disabledDates="",this.selectedDates="",this.oneDayInMs=864e5}localeChanged(){this.dateFormatter.locale=this.locale}dayFormatChanged(){this.dateFormatter.dayFormat=this.dayFormat}weekdayFormatChanged(){this.dateFormatter.weekdayFormat=this.weekdayFormat}monthFormatChanged(){this.dateFormatter.monthFormat=this.monthFormat}yearFormatChanged(){this.dateFormatter.yearFormat=this.yearFormat}getMonthInfo(e=this.month,t=this.year){let i=e=>new Date(e.getFullYear(),e.getMonth(),1).getDay(),s=e=>new Date(new Date(e.getFullYear(),e.getMonth()+1,1).getTime()-this.oneDayInMs).getDate(),o=new Date(t,e-1),n=new Date(t,e),a=new Date(t,e-2);return{length:s(o),month:e,start:i(o),year:t,previous:{length:s(a),month:a.getMonth()+1,start:i(a),year:a.getFullYear()},next:{length:s(n),month:n.getMonth()+1,start:i(n),year:n.getFullYear()}}}getDays(e=this.getMonthInfo(),t=this.minWeeks){t=t>10?10:t;let{start:i,length:s,previous:o,next:n}=e,a=[],r=1-i;for(;r<s+1||a.length<t||a[a.length-1].length%7!=0;){let{month:t,year:i}=r<1?o:r>s?n:e,l=r<1?o.length+r:r>s?r-s:r,h=`${t}-${l}-${i}`,d={day:l,month:t,year:i,disabled:this.dateInString(h,this.disabledDates),selected:this.dateInString(h,this.selectedDates)},c=a[a.length-1];0===a.length||c.length%7==0?a.push([d]):c.push(d),r++}return a}dateInString(e,t){let i=t.split(",").map(e=>e.trim());return e="string"==typeof e?e:`${e.getMonth()+1}-${e.getDate()}-${e.getFullYear()}`,i.some(t=>t===e)}getDayClassNames(e,t){let{day:i,month:s,year:o,disabled:n,selected:a}=e;return["day",t===`${s}-${i}-${o}`&&"today",this.month!==s&&"inactive",n&&"disabled",a&&"selected"].filter(Boolean).join(" ")}getWeekdayText(){let e=this.dateFormatter.getWeekdays().map(e=>({text:e}));if("long"!==this.weekdayFormat){let t=this.dateFormatter.getWeekdays("long");e.forEach((e,i)=>{e.abbr=t[i]})}return e}handleDateSelect(e,t){e.preventDefault,this.$emit("dateselected",t)}handleKeydown(e,t){return"Enter"===e.key&&this.handleDateSelect(e,t),!0}}f([(0,d.attr)({mode:"boolean"})],e8.prototype,"readonly",void 0),f([d.attr],e8.prototype,"locale",void 0),f([(0,d.attr)({converter:d.nullableNumberConverter})],e8.prototype,"month",void 0),f([(0,d.attr)({converter:d.nullableNumberConverter})],e8.prototype,"year",void 0),f([(0,d.attr)({attribute:"day-format",mode:"fromView"})],e8.prototype,"dayFormat",void 0),f([(0,d.attr)({attribute:"weekday-format",mode:"fromView"})],e8.prototype,"weekdayFormat",void 0),f([(0,d.attr)({attribute:"month-format",mode:"fromView"})],e8.prototype,"monthFormat",void 0),f([(0,d.attr)({attribute:"year-format",mode:"fromView"})],e8.prototype,"yearFormat",void 0),f([(0,d.attr)({attribute:"min-weeks",converter:d.nullableNumberConverter})],e8.prototype,"minWeeks",void 0),f([(0,d.attr)({attribute:"disabled-dates"})],e8.prototype,"disabledDates",void 0),f([(0,d.attr)({attribute:"selected-dates"})],e8.prototype,"selectedDates",void 0);let e2="focusin",e3="focusout",e9="keydown",e6={none:"none",default:"default",sticky:"sticky"},e7={default:"default",columnHeader:"columnheader",rowHeader:"rowheader"},te={default:"default",header:"header",stickyHeader:"sticky-header"},tt=(0,d.html)`
    <template>
        ${e=>null===e.rowData||null===e.columnDefinition||null===e.columnDefinition.columnDataKey?null:e.rowData[e.columnDefinition.columnDataKey]}
    </template>
`,ti=(0,d.html)`
    <template>
        ${e=>null===e.columnDefinition?null:void 0===e.columnDefinition.title?e.columnDefinition.columnDataKey:e.columnDefinition.title}
    </template>
`;class ts extends ep{constructor(){super(...arguments),this.cellType=e7.default,this.rowData=null,this.columnDefinition=null,this.isActiveCell=!1,this.customCellView=null,this.updateCellStyle=()=>{this.style.gridColumn=this.gridColumn}}cellTypeChanged(){this.$fastController.isConnected&&this.updateCellView()}gridColumnChanged(){this.$fastController.isConnected&&this.updateCellStyle()}columnDefinitionChanged(e,t){this.$fastController.isConnected&&this.updateCellView()}connectedCallback(){var e;super.connectedCallback(),this.addEventListener(e2,this.handleFocusin),this.addEventListener(e3,this.handleFocusout),this.addEventListener(e9,this.handleKeydown),this.style.gridColumn=`${(null==(e=this.columnDefinition)?void 0:e.gridColumn)===void 0?0:this.columnDefinition.gridColumn}`,this.updateCellView(),this.updateCellStyle()}disconnectedCallback(){super.disconnectedCallback(),this.removeEventListener(e2,this.handleFocusin),this.removeEventListener(e3,this.handleFocusout),this.removeEventListener(e9,this.handleKeydown),this.disconnectCellView()}handleFocusin(e){if(!this.isActiveCell){if(this.isActiveCell=!0,this.cellType===e7.columnHeader){if(null!==this.columnDefinition&&!0!==this.columnDefinition.headerCellInternalFocusQueue&&"function"==typeof this.columnDefinition.headerCellFocusTargetCallback){let e=this.columnDefinition.headerCellFocusTargetCallback(this);null!==e&&e.focus()}}else if(null!==this.columnDefinition&&!0!==this.columnDefinition.cellInternalFocusQueue&&"function"==typeof this.columnDefinition.cellFocusTargetCallback){let e=this.columnDefinition.cellFocusTargetCallback(this);null!==e&&e.focus()}this.$emit("cell-focused",this)}}handleFocusout(e){this===document.activeElement||this.contains(document.activeElement)||(this.isActiveCell=!1)}handleKeydown(e){if(!e.defaultPrevented&&null!==this.columnDefinition&&(this.cellType!==e7.default||!0===this.columnDefinition.cellInternalFocusQueue)&&(this.cellType!==e7.columnHeader||!0===this.columnDefinition.headerCellInternalFocusQueue))switch(e.key){case"Enter":case"F2":if(this.contains(document.activeElement)&&document.activeElement!==this)return;if(this.cellType===e7.columnHeader){if(void 0!==this.columnDefinition.headerCellFocusTargetCallback){let t=this.columnDefinition.headerCellFocusTargetCallback(this);null!==t&&t.focus(),e.preventDefault()}}else if(void 0!==this.columnDefinition.cellFocusTargetCallback){let t=this.columnDefinition.cellFocusTargetCallback(this);null!==t&&t.focus(),e.preventDefault()}break;case"Escape":this.contains(document.activeElement)&&document.activeElement!==this&&(this.focus(),e.preventDefault())}}updateCellView(){if(this.disconnectCellView(),null!==this.columnDefinition)switch(this.cellType){case e7.columnHeader:void 0!==this.columnDefinition.headerCellTemplate?this.customCellView=this.columnDefinition.headerCellTemplate.render(this,this):this.customCellView=ti.render(this,this);break;case void 0:case e7.rowHeader:case e7.default:void 0!==this.columnDefinition.cellTemplate?this.customCellView=this.columnDefinition.cellTemplate.render(this,this):this.customCellView=tt.render(this,this)}}disconnectCellView(){null!==this.customCellView&&(this.customCellView.dispose(),this.customCellView=null)}}f([(0,d.attr)({attribute:"cell-type"})],ts.prototype,"cellType",void 0),f([(0,d.attr)({attribute:"grid-column"})],ts.prototype,"gridColumn",void 0),f([d.observable],ts.prototype,"rowData",void 0),f([d.observable],ts.prototype,"columnDefinition",void 0);class to extends ep{constructor(){super(...arguments),this.rowType=te.default,this.rowData=null,this.columnDefinitions=null,this.isActiveRow=!1,this.cellsRepeatBehavior=null,this.cellsPlaceholder=null,this.focusColumnIndex=0,this.refocusOnLoad=!1,this.updateRowStyle=()=>{this.style.gridTemplateColumns=this.gridTemplateColumns}}gridTemplateColumnsChanged(){this.$fastController.isConnected&&this.updateRowStyle()}rowTypeChanged(){this.$fastController.isConnected&&this.updateItemTemplate()}rowDataChanged(){if(null!==this.rowData&&this.isActiveRow){this.refocusOnLoad=!0;return}}cellItemTemplateChanged(){this.updateItemTemplate()}headerCellItemTemplateChanged(){this.updateItemTemplate()}connectedCallback(){super.connectedCallback(),null===this.cellsRepeatBehavior&&(this.cellsPlaceholder=document.createComment(""),this.appendChild(this.cellsPlaceholder),this.updateItemTemplate(),this.cellsRepeatBehavior=new d.RepeatDirective(e=>e.columnDefinitions,e=>e.activeCellItemTemplate,{positioning:!0}).createBehavior(this.cellsPlaceholder),this.$fastController.addBehaviors([this.cellsRepeatBehavior])),this.addEventListener("cell-focused",this.handleCellFocus),this.addEventListener(e3,this.handleFocusout),this.addEventListener(e9,this.handleKeydown),this.updateRowStyle(),this.refocusOnLoad&&(this.refocusOnLoad=!1,this.cellElements.length>this.focusColumnIndex&&this.cellElements[this.focusColumnIndex].focus())}disconnectedCallback(){super.disconnectedCallback(),this.removeEventListener("cell-focused",this.handleCellFocus),this.removeEventListener(e3,this.handleFocusout),this.removeEventListener(e9,this.handleKeydown)}handleFocusout(e){this.contains(e.target)||(this.isActiveRow=!1,this.focusColumnIndex=0)}handleCellFocus(e){this.isActiveRow=!0,this.focusColumnIndex=this.cellElements.indexOf(e.target),this.$emit("row-focused",this)}handleKeydown(e){if(e.defaultPrevented)return;let t=0;switch(e.key){case ey.kT:t=Math.max(0,this.focusColumnIndex-1),this.cellElements[t].focus(),e.preventDefault();break;case ey.bb:t=Math.min(this.cellElements.length-1,this.focusColumnIndex+1),this.cellElements[t].focus(),e.preventDefault();break;case"Home":e.ctrlKey||(this.cellElements[0].focus(),e.preventDefault());break;case"End":e.ctrlKey||(this.cellElements[this.cellElements.length-1].focus(),e.preventDefault())}}updateItemTemplate(){this.activeCellItemTemplate=this.rowType===te.default&&void 0!==this.cellItemTemplate?this.cellItemTemplate:this.rowType===te.default&&void 0===this.cellItemTemplate?this.defaultCellItemTemplate:void 0!==this.headerCellItemTemplate?this.headerCellItemTemplate:this.defaultHeaderCellItemTemplate}}f([(0,d.attr)({attribute:"grid-template-columns"})],to.prototype,"gridTemplateColumns",void 0),f([(0,d.attr)({attribute:"row-type"})],to.prototype,"rowType",void 0),f([d.observable],to.prototype,"rowData",void 0),f([d.observable],to.prototype,"columnDefinitions",void 0),f([d.observable],to.prototype,"cellItemTemplate",void 0),f([d.observable],to.prototype,"headerCellItemTemplate",void 0),f([d.observable],to.prototype,"rowIndex",void 0),f([d.observable],to.prototype,"isActiveRow",void 0),f([d.observable],to.prototype,"activeCellItemTemplate",void 0),f([d.observable],to.prototype,"defaultCellItemTemplate",void 0),f([d.observable],to.prototype,"defaultHeaderCellItemTemplate",void 0),f([d.observable],to.prototype,"cellElements",void 0);class tn extends ep{constructor(){super(),this.noTabbing=!1,this.generateHeader=e6.default,this.rowsData=[],this.columnDefinitions=null,this.focusRowIndex=0,this.focusColumnIndex=0,this.rowsPlaceholder=null,this.generatedHeader=null,this.isUpdatingFocus=!1,this.pendingFocusUpdate=!1,this.rowindexUpdateQueued=!1,this.columnDefinitionsStale=!0,this.generatedGridTemplateColumns="",this.focusOnCell=(e,t,i)=>{if(0===this.rowElements.length){this.focusRowIndex=0,this.focusColumnIndex=0;return}let s=Math.max(0,Math.min(this.rowElements.length-1,e)),o=this.rowElements[s].querySelectorAll('[role="cell"], [role="gridcell"], [role="columnheader"], [role="rowheader"]'),n=Math.max(0,Math.min(o.length-1,t)),a=o[n];i&&this.scrollHeight!==this.clientHeight&&(s<this.focusRowIndex&&this.scrollTop>0||s>this.focusRowIndex&&this.scrollTop<this.scrollHeight-this.clientHeight)&&a.scrollIntoView({block:"center",inline:"center"}),a.focus()},this.onChildListChange=(e,t)=>{e&&e.length&&(e.forEach(e=>{e.addedNodes.forEach(e=>{1===e.nodeType&&"row"===e.getAttribute("role")&&(e.columnDefinitions=this.columnDefinitions)})}),this.queueRowIndexUpdate())},this.queueRowIndexUpdate=()=>{this.rowindexUpdateQueued||(this.rowindexUpdateQueued=!0,d.DOM.queueUpdate(this.updateRowIndexes))},this.updateRowIndexes=()=>{let e=this.gridTemplateColumns;if(void 0===e){if(""===this.generatedGridTemplateColumns&&this.rowElements.length>0){let e=this.rowElements[0];this.generatedGridTemplateColumns=Array(e.cellElements.length).fill("1fr").join(" ")}e=this.generatedGridTemplateColumns}this.rowElements.forEach((t,i)=>{t.rowIndex=i,t.gridTemplateColumns=e,this.columnDefinitionsStale&&(t.columnDefinitions=this.columnDefinitions)}),this.rowindexUpdateQueued=!1,this.columnDefinitionsStale=!1}}static generateTemplateColumns(e){let t="";return e.forEach(e=>{t=`${t}${""===t?"":" "}1fr`}),t}noTabbingChanged(){this.$fastController.isConnected&&(this.noTabbing?this.setAttribute("tabIndex","-1"):this.setAttribute("tabIndex",this.contains(document.activeElement)||this===document.activeElement?"-1":"0"))}generateHeaderChanged(){this.$fastController.isConnected&&this.toggleGeneratedHeader()}gridTemplateColumnsChanged(){this.$fastController.isConnected&&this.updateRowIndexes()}rowsDataChanged(){null===this.columnDefinitions&&this.rowsData.length>0&&(this.columnDefinitions=tn.generateColumns(this.rowsData[0])),this.$fastController.isConnected&&this.toggleGeneratedHeader()}columnDefinitionsChanged(){if(null===this.columnDefinitions){this.generatedGridTemplateColumns="";return}this.generatedGridTemplateColumns=tn.generateTemplateColumns(this.columnDefinitions),this.$fastController.isConnected&&(this.columnDefinitionsStale=!0,this.queueRowIndexUpdate())}headerCellItemTemplateChanged(){this.$fastController.isConnected&&null!==this.generatedHeader&&(this.generatedHeader.headerCellItemTemplate=this.headerCellItemTemplate)}focusRowIndexChanged(){this.$fastController.isConnected&&this.queueFocusUpdate()}focusColumnIndexChanged(){this.$fastController.isConnected&&this.queueFocusUpdate()}connectedCallback(){super.connectedCallback(),void 0===this.rowItemTemplate&&(this.rowItemTemplate=this.defaultRowItemTemplate),this.rowsPlaceholder=document.createComment(""),this.appendChild(this.rowsPlaceholder),this.toggleGeneratedHeader(),this.rowsRepeatBehavior=new d.RepeatDirective(e=>e.rowsData,e=>e.rowItemTemplate,{positioning:!0}).createBehavior(this.rowsPlaceholder),this.$fastController.addBehaviors([this.rowsRepeatBehavior]),this.addEventListener("row-focused",this.handleRowFocus),this.addEventListener("focus",this.handleFocus),this.addEventListener(e9,this.handleKeydown),this.addEventListener(e3,this.handleFocusOut),this.observer=new MutationObserver(this.onChildListChange),this.observer.observe(this,{childList:!0}),this.noTabbing&&this.setAttribute("tabindex","-1"),d.DOM.queueUpdate(this.queueRowIndexUpdate)}disconnectedCallback(){super.disconnectedCallback(),this.removeEventListener("row-focused",this.handleRowFocus),this.removeEventListener("focus",this.handleFocus),this.removeEventListener(e9,this.handleKeydown),this.removeEventListener(e3,this.handleFocusOut),this.observer.disconnect(),this.rowsPlaceholder=null,this.generatedHeader=null}handleRowFocus(e){this.isUpdatingFocus=!0;let t=e.target;this.focusRowIndex=this.rowElements.indexOf(t),this.focusColumnIndex=t.focusColumnIndex,this.setAttribute("tabIndex","-1"),this.isUpdatingFocus=!1}handleFocus(e){this.focusOnCell(this.focusRowIndex,this.focusColumnIndex,!0)}handleFocusOut(e){null!==e.relatedTarget&&this.contains(e.relatedTarget)||this.setAttribute("tabIndex",this.noTabbing?"-1":"0")}handleKeydown(e){let t;if(e.defaultPrevented)return;let i=this.rowElements.length-1,s=this.offsetHeight+this.scrollTop,o=this.rowElements[i];switch(e.key){case ey.I5:e.preventDefault(),this.focusOnCell(this.focusRowIndex-1,this.focusColumnIndex,!0);break;case ey.HX:e.preventDefault(),this.focusOnCell(this.focusRowIndex+1,this.focusColumnIndex,!0);break;case"PageUp":if(e.preventDefault(),0===this.rowElements.length){this.focusOnCell(0,0,!1);break}if(0===this.focusRowIndex)return void this.focusOnCell(0,this.focusColumnIndex,!1);for(t=this.focusRowIndex-1;t>=0;t--){let e=this.rowElements[t];if(e.offsetTop<this.scrollTop){this.scrollTop=e.offsetTop+e.clientHeight-this.clientHeight;break}}this.focusOnCell(t,this.focusColumnIndex,!1);break;case ey.f_:if(e.preventDefault(),0===this.rowElements.length){this.focusOnCell(0,0,!1);break}if(this.focusRowIndex>=i||o.offsetTop+o.offsetHeight<=s)return void this.focusOnCell(i,this.focusColumnIndex,!1);for(t=this.focusRowIndex+1;t<=i;t++){let e=this.rowElements[t];if(e.offsetTop+e.offsetHeight>s){let t=0;this.generateHeader===e6.sticky&&null!==this.generatedHeader&&(t=this.generatedHeader.clientHeight),this.scrollTop=e.offsetTop-t;break}}this.focusOnCell(t,this.focusColumnIndex,!1);break;case"Home":e.ctrlKey&&(e.preventDefault(),this.focusOnCell(0,0,!0));break;case"End":e.ctrlKey&&null!==this.columnDefinitions&&(e.preventDefault(),this.focusOnCell(this.rowElements.length-1,this.columnDefinitions.length-1,!0))}}queueFocusUpdate(){this.isUpdatingFocus&&(this.contains(document.activeElement)||this===document.activeElement)||!1===this.pendingFocusUpdate&&(this.pendingFocusUpdate=!0,d.DOM.queueUpdate(()=>this.updateFocus()))}updateFocus(){this.pendingFocusUpdate=!1,this.focusOnCell(this.focusRowIndex,this.focusColumnIndex,!0)}toggleGeneratedHeader(){if(null!==this.generatedHeader&&(this.removeChild(this.generatedHeader),this.generatedHeader=null),this.generateHeader!==e6.none&&this.rowsData.length>0){let e=document.createElement(this.rowElementTag);this.generatedHeader=e,this.generatedHeader.columnDefinitions=this.columnDefinitions,this.generatedHeader.gridTemplateColumns=this.gridTemplateColumns,this.generatedHeader.rowType=this.generateHeader===e6.sticky?te.stickyHeader:te.header,(null!==this.firstChild||null!==this.rowsPlaceholder)&&this.insertBefore(e,null!==this.firstChild?this.firstChild:this.rowsPlaceholder);return}}}tn.generateColumns=e=>Object.getOwnPropertyNames(e).map((e,t)=>({columnDataKey:e,gridColumn:`${t}`})),f([(0,d.attr)({attribute:"no-tabbing",mode:"boolean"})],tn.prototype,"noTabbing",void 0),f([(0,d.attr)({attribute:"generate-header"})],tn.prototype,"generateHeader",void 0),f([(0,d.attr)({attribute:"grid-template-columns"})],tn.prototype,"gridTemplateColumns",void 0),f([d.observable],tn.prototype,"rowsData",void 0),f([d.observable],tn.prototype,"columnDefinitions",void 0),f([d.observable],tn.prototype,"rowItemTemplate",void 0),f([d.observable],tn.prototype,"cellItemTemplate",void 0),f([d.observable],tn.prototype,"headerCellItemTemplate",void 0),f([d.observable],tn.prototype,"focusRowIndex",void 0),f([d.observable],tn.prototype,"focusColumnIndex",void 0),f([d.observable],tn.prototype,"defaultRowItemTemplate",void 0),f([d.observable],tn.prototype,"rowElementTag",void 0),f([d.observable],tn.prototype,"rowElements",void 0);let ta=(0,d.html)`
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
`,tr=e=>{let t=e.tagFor(ts);return(0,d.html)`
        <${t}
            class="week-day"
            part="week-day"
            tabindex="-1"
            grid-column="${(e,t)=>t.index+1}"
            abbr="${e=>e.abbr}"
        >
            ${e=>e.text}
        </${t}>
    `},tl=(e,t)=>{let i=e.tagFor(ts);return(0,d.html)`
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
    `},th=(e,t)=>{let i=e.tagFor(to);return(0,d.html)`
        <${i}
            class="week"
            part="week"
            role="row"
            role-type="default"
            grid-template-columns="1fr 1fr 1fr 1fr 1fr 1fr 1fr"
        >
        ${(0,d.repeat)(e=>e,tl(e,t),{positioning:!0})}
        </${i}>
    `},td=(e,t)=>{let i=e.tagFor(tn),s=e.tagFor(to);return(0,d.html)`
    <${i} class="days interact" part="days" generate-header="none">
        <${s}
            class="week-days"
            part="week-days"
            role="row"
            row-type="header"
            grid-template-columns="1fr 1fr 1fr 1fr 1fr 1fr 1fr"
        >
            ${(0,d.repeat)(e=>e.getWeekdayText(),tr(e),{positioning:!0})}
        </${s}>
        ${(0,d.repeat)(e=>e.getDays(),th(e,t))}
    </${i}>
`},tc=e=>(0,d.html)`
        <div class="days" part="days">
            <div class="week-days" part="week-days">
                ${(0,d.repeat)(e=>e.getWeekdayText(),(0,d.html)`
                        <div class="week-day" part="week-day" abbr="${e=>e.abbr}">
                            ${e=>e.text}
                        </div>
                    `)}
            </div>
            ${(0,d.repeat)(e=>e.getDays(),(0,d.html)`
                    <div class="week">
                        ${(0,d.repeat)(e=>e,(0,d.html)`
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
    `,tu=(e,t)=>{var i;let s=new Date,o=`${s.getMonth()+1}-${s.getDate()}-${s.getFullYear()}`;return(0,d.html)`
        <template>
            ${v}
            ${t.title instanceof Function?t.title(e,t):null!=(i=t.title)?i:""}
            <slot></slot>
            ${(0,d.when)(e=>e.readonly,tc(o),td(e,o))}
            ${m}
        </template>
    `},tp=(e,t)=>(0,d.html)`
    <slot></slot>
`;class tm extends ep{}let tv=(e,t)=>(0,d.html)`
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
            <slot ${(0,d.slotted)("defaultSlottedNodes")}></slot>
        </label>
    </template>
`;class tb extends ep{}class tf extends eZ(tb){constructor(){super(...arguments),this.proxy=document.createElement("input")}}class tg extends tf{constructor(){super(),this.initialValue="on",this.indeterminate=!1,this.keypressHandler=e=>{this.readOnly||" "===e.key&&(this.indeterminate&&(this.indeterminate=!1),this.checked=!this.checked)},this.clickHandler=e=>{this.disabled||this.readOnly||(this.indeterminate&&(this.indeterminate=!1),this.checked=!this.checked)},this.proxy.setAttribute("type","checkbox")}readOnlyChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.readOnly=this.readOnly)}}f([(0,d.attr)({attribute:"readonly",mode:"boolean"})],tg.prototype,"readOnly",void 0),f([d.observable],tg.prototype,"defaultSlottedNodes",void 0),f([d.observable],tg.prototype,"indeterminate",void 0);let ty=0;function tC(e=""){return`${e}${ty++}`}function tx(...e){return e.every(e=>e instanceof HTMLElement)}function t$(e){return tx(e)&&("option"===e.getAttribute("role")||e instanceof HTMLOptionElement)}class tw extends ep{constructor(e,t,i,s){super(),this.defaultSelected=!1,this.dirtySelected=!1,this.selected=this.defaultSelected,this.dirtyValue=!1,e&&(this.textContent=e),t&&(this.initialValue=t),i&&(this.defaultSelected=i),s&&(this.selected=s),this.proxy=new Option(`${this.textContent}`,this.initialValue,this.defaultSelected,this.selected),this.proxy.disabled=this.disabled}checkedChanged(e,t){if("boolean"==typeof t){this.ariaChecked=t?"true":"false";return}this.ariaChecked=null}contentChanged(e,t){this.proxy instanceof HTMLOptionElement&&(this.proxy.textContent=this.textContent),this.$emit("contentchange",null,{bubbles:!0})}defaultSelectedChanged(){!this.dirtySelected&&(this.selected=this.defaultSelected,this.proxy instanceof HTMLOptionElement&&(this.proxy.selected=this.defaultSelected))}disabledChanged(e,t){this.ariaDisabled=this.disabled?"true":"false",this.proxy instanceof HTMLOptionElement&&(this.proxy.disabled=this.disabled)}selectedAttributeChanged(){this.defaultSelected=this.selectedAttribute,this.proxy instanceof HTMLOptionElement&&(this.proxy.defaultSelected=this.defaultSelected)}selectedChanged(){this.ariaSelected=this.selected?"true":"false",this.dirtySelected||(this.dirtySelected=!0),this.proxy instanceof HTMLOptionElement&&(this.proxy.selected=this.selected)}initialValueChanged(e,t){this.dirtyValue||(this.value=this.initialValue,this.dirtyValue=!1)}get label(){var e;return null!=(e=this.value)?e:this.text}get text(){var e,t;return null!=(t=null==(e=this.textContent)?void 0:e.replace(/\s+/g," ").trim())?t:""}set value(e){let t=`${null!=e?e:""}`;this._value=t,this.dirtyValue=!0,this.proxy instanceof HTMLOptionElement&&(this.proxy.value=t),d.Observable.notify(this,"value")}get value(){var e;return d.Observable.track(this,"value"),null!=(e=this._value)?e:this.text}get form(){return this.proxy?this.proxy.form:null}}f([d.observable],tw.prototype,"checked",void 0),f([d.observable],tw.prototype,"content",void 0),f([d.observable],tw.prototype,"defaultSelected",void 0),f([(0,d.attr)({mode:"boolean"})],tw.prototype,"disabled",void 0),f([(0,d.attr)({attribute:"selected",mode:"boolean"})],tw.prototype,"selectedAttribute",void 0),f([d.observable],tw.prototype,"selected",void 0),f([(0,d.attr)({attribute:"value",mode:"fromView"})],tw.prototype,"initialValue",void 0);class tI{}f([d.observable],tI.prototype,"ariaChecked",void 0),f([d.observable],tI.prototype,"ariaPosInSet",void 0),f([d.observable],tI.prototype,"ariaSelected",void 0),f([d.observable],tI.prototype,"ariaSetSize",void 0),eb(tI,eI),eb(tw,c,tI);class tk extends ep{constructor(){super(...arguments),this._options=[],this.selectedIndex=-1,this.selectedOptions=[],this.shouldSkipFocus=!1,this.typeaheadBuffer="",this.typeaheadExpired=!0,this.typeaheadTimeout=-1}get firstSelectedOption(){var e;return null!=(e=this.selectedOptions[0])?e:null}get hasSelectableOptions(){return this.options.length>0&&!this.options.every(e=>e.disabled)}get length(){var e,t;return null!=(t=null==(e=this.options)?void 0:e.length)?t:0}get options(){return d.Observable.track(this,"options"),this._options}set options(e){this._options=e,d.Observable.notify(this,"options")}get typeAheadExpired(){return this.typeaheadExpired}set typeAheadExpired(e){this.typeaheadExpired=e}clickHandler(e){let t=e.target.closest("option,[role=option]");if(t&&!t.disabled)return this.selectedIndex=this.options.indexOf(t),!0}focusAndScrollOptionIntoView(e=this.firstSelectedOption){this.contains(document.activeElement)&&null!==e&&(e.focus(),requestAnimationFrame(()=>{e.scrollIntoView({block:"nearest"})}))}focusinHandler(e){this.shouldSkipFocus||e.target!==e.currentTarget||(this.setSelectedOptions(),this.focusAndScrollOptionIntoView()),this.shouldSkipFocus=!1}getTypeaheadMatches(){let e=this.typeaheadBuffer.replace(/[.*+\-?^${}()|[\]\\]/g,"\\$&"),t=RegExp(`^${e}`,"gi");return this.options.filter(e=>e.text.trim().match(t))}getSelectableIndex(e=this.selectedIndex,t){let i=e>t?-1:+(e<t),s=e+i,o=null;switch(i){case -1:o=this.options.reduceRight((e,t,i)=>e||t.disabled||!(i<s)?e:t,o);break;case 1:o=this.options.reduce((e,t,i)=>e||t.disabled||!(i>s)?e:t,o)}return this.options.indexOf(o)}handleChange(e,t){"selected"===t&&(tk.slottedOptionFilter(e)&&(this.selectedIndex=this.options.indexOf(e)),this.setSelectedOptions())}handleTypeAhead(e){this.typeaheadTimeout&&window.clearTimeout(this.typeaheadTimeout),this.typeaheadTimeout=window.setTimeout(()=>this.typeaheadExpired=!0,tk.TYPE_AHEAD_TIMEOUT_MS),e.length>1||(this.typeaheadBuffer=`${this.typeaheadExpired?"":this.typeaheadBuffer}${e}`)}keydownHandler(e){if(this.disabled)return!0;this.shouldSkipFocus=!1;let t=e.key;switch(t){case"Home":e.shiftKey||(e.preventDefault(),this.selectFirstOption());break;case ey.HX:e.shiftKey||(e.preventDefault(),this.selectNextOption());break;case ey.I5:e.shiftKey||(e.preventDefault(),this.selectPreviousOption());break;case"End":e.preventDefault(),this.selectLastOption();break;case"Tab":return this.focusAndScrollOptionIntoView(),!0;case"Enter":case"Escape":return!0;case" ":if(this.typeaheadExpired)return!0;default:return 1===t.length&&this.handleTypeAhead(`${t}`),!0}}mousedownHandler(e){return this.shouldSkipFocus=!this.contains(document.activeElement),!0}multipleChanged(e,t){this.ariaMultiSelectable=t?"true":null}selectedIndexChanged(e,t){var i;if(!this.hasSelectableOptions){this.selectedIndex=-1;return}if((null==(i=this.options[this.selectedIndex])?void 0:i.disabled)&&"number"==typeof e){let i=this.getSelectableIndex(e,t),s=i>-1?i:e;this.selectedIndex=s,t===s&&this.selectedIndexChanged(t,s);return}this.setSelectedOptions()}selectedOptionsChanged(e,t){var i;let s=t.filter(tk.slottedOptionFilter);null==(i=this.options)||i.forEach(e=>{let t=d.Observable.getNotifier(e);t.unsubscribe(this,"selected"),e.selected=s.includes(e),t.subscribe(this,"selected")})}selectFirstOption(){var e,t;this.disabled||(this.selectedIndex=null!=(t=null==(e=this.options)?void 0:e.findIndex(e=>!e.disabled))?t:-1)}selectLastOption(){this.disabled||(this.selectedIndex=function(e,t){let i=e.length;for(;i--;)if(t(e[i],i,e))return i;return -1}(this.options,e=>!e.disabled))}selectNextOption(){!this.disabled&&this.selectedIndex<this.options.length-1&&(this.selectedIndex+=1)}selectPreviousOption(){!this.disabled&&this.selectedIndex>0&&(this.selectedIndex=this.selectedIndex-1)}setDefaultSelectedOption(){var e,t;this.selectedIndex=null!=(t=null==(e=this.options)?void 0:e.findIndex(e=>e.defaultSelected))?t:-1}setSelectedOptions(){var e,t,i;(null==(e=this.options)?void 0:e.length)&&(this.selectedOptions=[this.options[this.selectedIndex]],this.ariaActiveDescendant=null!=(i=null==(t=this.firstSelectedOption)?void 0:t.id)?i:"",this.focusAndScrollOptionIntoView())}slottedOptionsChanged(e,t){this.options=t.reduce((e,t)=>(t$(t)&&e.push(t),e),[]);let i=`${this.options.length}`;this.options.forEach((e,t)=>{e.id||(e.id=tC("option-")),e.ariaPosInSet=`${t+1}`,e.ariaSetSize=i}),this.$fastController.isConnected&&(this.setSelectedOptions(),this.setDefaultSelectedOption())}typeaheadBufferChanged(e,t){if(this.$fastController.isConnected){let e=this.getTypeaheadMatches();if(e.length){let t=this.options.indexOf(e[0]);t>-1&&(this.selectedIndex=t)}this.typeaheadExpired=!1}}}tk.slottedOptionFilter=e=>t$(e)&&!e.hidden,tk.TYPE_AHEAD_TIMEOUT_MS=1e3,f([(0,d.attr)({mode:"boolean"})],tk.prototype,"disabled",void 0),f([d.observable],tk.prototype,"selectedIndex",void 0),f([d.observable],tk.prototype,"selectedOptions",void 0),f([d.observable],tk.prototype,"slottedOptions",void 0),f([d.observable],tk.prototype,"typeaheadBuffer",void 0);class tE{}f([d.observable],tE.prototype,"ariaActiveDescendant",void 0),f([d.observable],tE.prototype,"ariaDisabled",void 0),f([d.observable],tE.prototype,"ariaExpanded",void 0),f([d.observable],tE.prototype,"ariaMultiSelectable",void 0),eb(tE,eI),eb(tk,tE);let tO={above:"above",below:"below"};class tT extends tk{}class tR extends eQ(tT){constructor(){super(...arguments),this.proxy=document.createElement("input")}}let tD={inline:"inline",list:"list",both:"both",none:"none"};class tS extends tR{constructor(){super(...arguments),this._value="",this.filteredOptions=[],this.filter="",this.forcedPosition=!1,this.listboxId=tC("listbox-"),this.maxHeight=0,this.open=!1}formResetCallback(){super.formResetCallback(),this.setDefaultSelectedOption(),this.updateValue()}validate(){super.validate(this.control)}get isAutocompleteInline(){return this.autocomplete===tD.inline||this.isAutocompleteBoth}get isAutocompleteList(){return this.autocomplete===tD.list||this.isAutocompleteBoth}get isAutocompleteBoth(){return this.autocomplete===tD.both}openChanged(){if(this.open){this.ariaControls=this.listboxId,this.ariaExpanded="true",this.setPositioning(),this.focusAndScrollOptionIntoView(),d.DOM.queueUpdate(()=>this.focus());return}this.ariaControls="",this.ariaExpanded="false"}get options(){return d.Observable.track(this,"options"),this.filteredOptions.length?this.filteredOptions:this._options}set options(e){this._options=e,d.Observable.notify(this,"options")}placeholderChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.placeholder=this.placeholder)}positionChanged(e,t){this.positionAttribute=t,this.setPositioning()}get value(){return d.Observable.track(this,"value"),this._value}set value(e){var t,i,s;let o=`${this._value}`;if(this.$fastController.isConnected&&this.options){let o=this.options.findIndex(t=>t.text.toLowerCase()===e.toLowerCase()),n=null==(t=this.options[this.selectedIndex])?void 0:t.text,a=null==(i=this.options[o])?void 0:i.text;this.selectedIndex=n!==a?o:this.selectedIndex,e=(null==(s=this.firstSelectedOption)?void 0:s.text)||e}o!==e&&(this._value=e,super.valueChanged(o,e),d.Observable.notify(this,"value"))}clickHandler(e){if(!this.disabled){if(this.open){let t=e.target.closest("option,[role=option]");if(!t||t.disabled)return;this.selectedOptions=[t],this.control.value=t.text,this.clearSelectionRange(),this.updateValue(!0)}return this.open=!this.open,this.open&&this.control.focus(),!0}}connectedCallback(){super.connectedCallback(),this.forcedPosition=!!this.positionAttribute,this.value&&(this.initialValue=this.value)}disabledChanged(e,t){super.disabledChanged&&super.disabledChanged(e,t),this.ariaDisabled=this.disabled?"true":"false"}filterOptions(){this.autocomplete&&this.autocomplete!==tD.none||(this.filter="");let e=this.filter.toLowerCase();this.filteredOptions=this._options.filter(e=>e.text.toLowerCase().startsWith(this.filter.toLowerCase())),this.isAutocompleteList&&(this.filteredOptions.length||e||(this.filteredOptions=this._options),this._options.forEach(e=>{e.hidden=!this.filteredOptions.includes(e)}))}focusAndScrollOptionIntoView(){this.contains(document.activeElement)&&(this.control.focus(),this.firstSelectedOption&&requestAnimationFrame(()=>{var e;null==(e=this.firstSelectedOption)||e.scrollIntoView({block:"nearest"})}))}focusoutHandler(e){if(this.syncValue(),!this.open)return!0;let t=e.relatedTarget;this.isSameNode(t)?this.focus():this.options&&this.options.includes(t)||(this.open=!1)}inputHandler(e){if(this.filter=this.control.value,this.filterOptions(),this.isAutocompleteInline||(this.selectedIndex=this.options.map(e=>e.text).indexOf(this.control.value)),e.inputType.includes("deleteContent")||!this.filter.length)return!0;this.isAutocompleteList&&!this.open&&(this.open=!0),this.isAutocompleteInline&&(this.filteredOptions.length?(this.selectedOptions=[this.filteredOptions[0]],this.selectedIndex=this.options.indexOf(this.firstSelectedOption),this.setInlineSelection()):this.selectedIndex=-1)}keydownHandler(e){let t=e.key;if(e.ctrlKey||e.shiftKey)return!0;switch(t){case"Enter":this.syncValue(),this.isAutocompleteInline&&(this.filter=this.value),this.open=!1,this.clearSelectionRange();break;case"Escape":if(this.isAutocompleteInline||(this.selectedIndex=-1),this.open){this.open=!1;break}this.value="",this.control.value="",this.filter="",this.filterOptions();break;case"Tab":if(this.setInputToSelection(),!this.open)return!0;e.preventDefault(),this.open=!1;break;case"ArrowUp":case"ArrowDown":if(this.filterOptions(),!this.open){this.open=!0;break}this.filteredOptions.length>0&&super.keydownHandler(e),this.isAutocompleteInline&&this.setInlineSelection();break;default:return!0}}keyupHandler(e){switch(e.key){case"ArrowLeft":case"ArrowRight":case"Backspace":case"Delete":case"Home":case"End":this.filter=this.control.value,this.selectedIndex=-1,this.filterOptions()}}selectedIndexChanged(e,t){if(this.$fastController.isConnected){if((t=(0,eC.AB)(-1,this.options.length-1,t))!==this.selectedIndex){this.selectedIndex=t;return}super.selectedIndexChanged(e,t)}}selectPreviousOption(){!this.disabled&&this.selectedIndex>=0&&(this.selectedIndex=this.selectedIndex-1)}setDefaultSelectedOption(){if(this.$fastController.isConnected&&this.options){let e=this.options.findIndex(e=>null!==e.getAttribute("selected")||e.selected);this.selectedIndex=e,!this.dirtyValue&&this.firstSelectedOption&&(this.value=this.firstSelectedOption.text),this.setSelectedOptions()}}setInputToSelection(){this.firstSelectedOption&&(this.control.value=this.firstSelectedOption.text,this.control.focus())}setInlineSelection(){this.firstSelectedOption&&(this.setInputToSelection(),this.control.setSelectionRange(this.filter.length,this.control.value.length,"backward"))}syncValue(){var e;let t=this.selectedIndex>-1?null==(e=this.firstSelectedOption)?void 0:e.text:this.control.value;this.updateValue(this.value!==t)}setPositioning(){let e=this.getBoundingClientRect(),t=window.innerHeight-e.bottom;this.position=this.forcedPosition?this.positionAttribute:e.top>t?tO.above:tO.below,this.positionAttribute=this.forcedPosition?this.positionAttribute:this.position,this.maxHeight=this.position===tO.above?~~e.top:~~t}selectedOptionsChanged(e,t){this.$fastController.isConnected&&this._options.forEach(e=>{e.selected=t.includes(e)})}slottedOptionsChanged(e,t){super.slottedOptionsChanged(e,t),this.updateValue()}updateValue(e){var t;this.$fastController.isConnected&&(this.value=(null==(t=this.firstSelectedOption)?void 0:t.text)||this.control.value,this.control.value=this.value),e&&this.$emit("change")}clearSelectionRange(){let e=this.control.value.length;this.control.setSelectionRange(e,e)}}f([(0,d.attr)({attribute:"autocomplete",mode:"fromView"})],tS.prototype,"autocomplete",void 0),f([d.observable],tS.prototype,"maxHeight",void 0),f([(0,d.attr)({attribute:"open",mode:"boolean"})],tS.prototype,"open",void 0),f([d.attr],tS.prototype,"placeholder",void 0),f([(0,d.attr)({attribute:"position"})],tS.prototype,"positionAttribute",void 0),f([d.observable],tS.prototype,"position",void 0);class tA{}f([d.observable],tA.prototype,"ariaAutoComplete",void 0),f([d.observable],tA.prototype,"ariaControls",void 0),eb(tA,tE),eb(tS,c,tA);let tF=(e,t)=>(0,d.html)`
    <template
        aria-disabled="${e=>e.ariaDisabled}"
        autocomplete="${e=>e.autocomplete}"
        class="${e=>e.open?"open":""} ${e=>e.disabled?"disabled":""} ${e=>e.position}"
        ?open="${e=>e.open}"
        tabindex="${e=>e.disabled?null:"0"}"
        @click="${(e,t)=>e.clickHandler(t.event)}"
        @focusout="${(e,t)=>e.focusoutHandler(t.event)}"
        @keydown="${(e,t)=>e.keydownHandler(t.event)}"
    >
        <div class="control" part="control">
            ${p(e,t)}
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
                    ${(0,d.ref)("control")}
                />
                <div class="indicator" part="indicator" aria-hidden="true">
                    <slot name="indicator">
                        ${t.indicator||""}
                    </slot>
                </div>
            </slot>
            ${u(e,t)}
        </div>
        <div
            class="listbox"
            id="${e=>e.listboxId}"
            part="listbox"
            role="listbox"
            ?disabled="${e=>e.disabled}"
            ?hidden="${e=>!e.open}"
            ${(0,d.ref)("listbox")}
        >
            <slot
                ${(0,d.slotted)({filter:tk.slottedOptionFilter,flatten:!0,property:"slottedOptions"})}
            ></slot>
        </div>
    </template>
`,tL=(e,t)=>{let i,s=(i=e.tagFor(to),(0,d.html)`
    <${i}
        :rowData="${e=>e}"
        :cellItemTemplate="${(e,t)=>t.parent.cellItemTemplate}"
        :headerCellItemTemplate="${(e,t)=>t.parent.headerCellItemTemplate}"
    ></${i}>
`),o=e.tagFor(to);return(0,d.html)`
        <template
            role="grid"
            tabindex="0"
            :rowElementTag="${()=>o}"
            :defaultRowItemTemplate="${s}"
            ${(0,d.children)({property:"rowElements",filter:(0,d.elements)("[role=row]")})}
        >
            <slot></slot>
        </template>
    `},tM=(e,t)=>{let i,s,o=(i=e.tagFor(ts),(0,d.html)`
    <${i}
        cell-type="${e=>e.isRowHeader?"rowheader":void 0}"
        grid-column="${(e,t)=>t.index+1}"
        :rowData="${(e,t)=>t.parent.rowData}"
        :columnDefinition="${e=>e}"
    ></${i}>
`),n=(s=e.tagFor(ts),(0,d.html)`
    <${s}
        cell-type="columnheader"
        grid-column="${(e,t)=>t.index+1}"
        :columnDefinition="${e=>e}"
    ></${s}>
`);return(0,d.html)`
        <template
            role="row"
            class="${e=>"default"!==e.rowType?e.rowType:""}"
            :defaultCellItemTemplate="${o}"
            :defaultHeaderCellItemTemplate="${n}"
            ${(0,d.children)({property:"cellElements",filter:(0,d.elements)('[role="cell"],[role="gridcell"],[role="columnheader"],[role="rowheader"]')})}
        >
            <slot ${(0,d.slotted)("slottedCellElements")}></slot>
        </template>
    `},tP=(e,t)=>(0,d.html)`
        <template
            tabindex="-1"
            role="${e=>e.cellType&&"default"!==e.cellType?e.cellType:"gridcell"}"
            class="
            ${e=>"columnheader"===e.cellType?"column-header":"rowheader"===e.cellType?"row-header":""}
            "
        >
            <slot></slot>
        </template>
    `;function tH(e){let t=e.parentElement;if(t)return t;{let t=e.getRootNode();if(t.host instanceof HTMLElement)return t.host}return null}function tz(e,t){let i=t;for(;null!==i;){if(i===e)return!0;i=tH(i)}return!1}let tV=document.createElement("div");class tN{setProperty(e,t){d.DOM.queueUpdate(()=>this.target.setProperty(e,t))}removeProperty(e){d.DOM.queueUpdate(()=>this.target.removeProperty(e))}}class tB extends tN{constructor(){super();const e=new CSSStyleSheet;this.target=e.cssRules[e.insertRule(":root{}")].style,document.adoptedStyleSheets=[...document.adoptedStyleSheets,e]}}class tq extends tN{constructor(){super(),this.style=document.createElement("style"),document.head.appendChild(this.style);const{sheet:e}=this.style;if(e){const t=e.insertRule(":root{}",e.cssRules.length);this.target=e.cssRules[t].style}}}class tU{constructor(e){this.store=new Map,this.target=null;const t=e.$fastController;this.style=document.createElement("style"),t.addStyles(this.style),d.Observable.getNotifier(t).subscribe(this,"isConnected"),this.handleChange(t,"isConnected")}targetChanged(){if(null!==this.target)for(let[e,t]of this.store.entries())this.target.setProperty(e,t)}setProperty(e,t){this.store.set(e,t),d.DOM.queueUpdate(()=>{null!==this.target&&this.target.setProperty(e,t)})}removeProperty(e){this.store.delete(e),d.DOM.queueUpdate(()=>{null!==this.target&&this.target.removeProperty(e)})}handleChange(e,t){let{sheet:i}=this.style;if(i){let e=i.insertRule(":host{}",i.cssRules.length);this.target=i.cssRules[e].style}else this.target=null}}f([d.observable],tU.prototype,"target",void 0);class tj{constructor(e){this.target=e.style}setProperty(e,t){d.DOM.queueUpdate(()=>this.target.setProperty(e,t))}removeProperty(e){d.DOM.queueUpdate(()=>this.target.removeProperty(e))}}class t_{setProperty(e,t){for(let i of(t_.properties[e]=t,t_.roots.values()))tX.getOrCreate(t_.normalizeRoot(i)).setProperty(e,t)}removeProperty(e){for(let t of(delete t_.properties[e],t_.roots.values()))tX.getOrCreate(t_.normalizeRoot(t)).removeProperty(e)}static registerRoot(e){let{roots:t}=t_;if(!t.has(e)){t.add(e);let i=tX.getOrCreate(this.normalizeRoot(e));for(let e in t_.properties)i.setProperty(e,t_.properties[e])}}static unregisterRoot(e){let{roots:t}=t_;if(t.has(e)){t.delete(e);let i=tX.getOrCreate(t_.normalizeRoot(e));for(let e in t_.properties)i.removeProperty(e)}}static normalizeRoot(e){return e===tV?document:e}}t_.roots=new Set,t_.properties={};let tK=new WeakMap,tW=d.DOM.supportsAdoptedStyleSheets?class extends tN{constructor(e){super();const t=new CSSStyleSheet;this.target=t.cssRules[t.insertRule(":host{}")].style,e.$fastController.addStyles(d.ElementStyles.create([t]))}}:tU,tX=Object.freeze({getOrCreate(e){let t;return tK.has(e)?tK.get(e):(t=e===tV?new t_:e instanceof Document?d.DOM.supportsAdoptedStyleSheets?new tB:new tq:e instanceof d.FASTElement?new tW(e):new tj(e),tK.set(e,t),t)}});class tG extends d.CSSDirective{constructor(e){super(),this.subscribers=new WeakMap,this._appliedTo=new Set,this.name=e.name,null!==e.cssCustomPropertyName&&(this.cssCustomProperty=`--${e.cssCustomPropertyName}`,this.cssVar=`var(${this.cssCustomProperty})`),this.id=tG.uniqueId(),tG.tokensById.set(this.id,this)}get appliedTo(){return[...this._appliedTo]}static from(e){return new tG({name:"string"==typeof e?e:e.name,cssCustomPropertyName:"string"==typeof e?e:void 0===e.cssCustomPropertyName?e.name:e.cssCustomPropertyName})}static isCSSDesignToken(e){return"string"==typeof e.cssCustomProperty}static isDerivedDesignTokenValue(e){return"function"==typeof e}static getTokenById(e){return tG.tokensById.get(e)}getOrCreateSubscriberSet(e=this){return this.subscribers.get(e)||this.subscribers.set(e,new Set)&&this.subscribers.get(e)}createCSS(){return this.cssVar||""}getValueFor(e){let t=t0.getOrCreate(e).get(this);if(void 0!==t)return t;throw Error(`Value could not be retrieved for token named "${this.name}". Ensure the value is set for ${e} or an ancestor of ${e}.`)}setValueFor(e,t){return this._appliedTo.add(e),t instanceof tG&&(t=this.alias(t)),t0.getOrCreate(e).set(this,t),this}deleteValueFor(e){return this._appliedTo.delete(e),t0.existsFor(e)&&t0.getOrCreate(e).delete(this),this}withDefault(e){return this.setValueFor(tV,e),this}subscribe(e,t){let i=this.getOrCreateSubscriberSet(t);t&&!t0.existsFor(t)&&t0.getOrCreate(t),i.has(e)||i.add(e)}unsubscribe(e,t){let i=this.subscribers.get(t||this);i&&i.has(e)&&i.delete(e)}notify(e){let t=Object.freeze({token:this,target:e});this.subscribers.has(this)&&this.subscribers.get(this).forEach(e=>e.handleChange(t)),this.subscribers.has(e)&&this.subscribers.get(e).forEach(e=>e.handleChange(t))}alias(e){return t=>e.getValueFor(t)}}l=0,tG.uniqueId=()=>(l++,l.toString(16)),tG.tokensById=new Map;class tY{constructor(e,t,i){this.source=e,this.token=t,this.node=i,this.dependencies=new Set,this.observer=d.Observable.binding(e,this,!1),this.observer.handleChange=this.observer.call,this.handleChange()}disconnect(){this.observer.disconnect()}handleChange(){this.node.store.set(this.token,this.observer.observe(this.node.target,d.defaultExecutionContext))}}class tQ{constructor(){this.values=new Map}set(e,t){this.values.get(e)!==t&&(this.values.set(e,t),d.Observable.getNotifier(this).notify(e.id))}get(e){return d.Observable.track(this,e.id),this.values.get(e)}delete(e){this.values.delete(e)}all(){return this.values.entries()}}let tZ=new WeakMap,tJ=new WeakMap;class t0{constructor(e){this.target=e,this.store=new tQ,this.children=[],this.assignedValues=new Map,this.reflecting=new Set,this.bindingObservers=new Map,this.tokenValueChangeHandler={handleChange:(e,t)=>{let i=tG.getTokenById(t);if(i&&(i.notify(this.target),tG.isCSSDesignToken(i))){let t=this.parent,s=this.isReflecting(i);if(t){let o=t.get(i),n=e.get(i);o===n||s?o===n&&s&&this.stopReflectToCSS(i):this.reflectToCSS(i)}else s||this.reflectToCSS(i)}}},tZ.set(e,this),d.Observable.getNotifier(this.store).subscribe(this.tokenValueChangeHandler),e instanceof d.FASTElement?e.$fastController.addBehaviors([this]):e.isConnected&&this.bind()}static getOrCreate(e){return tZ.get(e)||new t0(e)}static existsFor(e){return tZ.has(e)}static findParent(e){if(tV!==e.target){let t=tH(e.target);for(;null!==t;){if(tZ.has(t))return tZ.get(t);t=tH(t)}return t0.getOrCreate(tV)}return null}static findClosestAssignedNode(e,t){let i=t;do{if(i.has(e))return i;i=i.parent?i.parent:i.target!==tV?t0.getOrCreate(tV):null}while(null!==i);return null}get parent(){return tJ.get(this)||null}has(e){return this.assignedValues.has(e)}get(e){let t=this.store.get(e);if(void 0!==t)return t;let i=this.getRaw(e);if(void 0!==i)return this.hydrate(e,i),this.get(e)}getRaw(e){var t;return this.assignedValues.has(e)?this.assignedValues.get(e):null==(t=t0.findClosestAssignedNode(e,this))?void 0:t.getRaw(e)}set(e,t){tG.isDerivedDesignTokenValue(this.assignedValues.get(e))&&this.tearDownBindingObserver(e),this.assignedValues.set(e,t),tG.isDerivedDesignTokenValue(t)?this.setupBindingObserver(e,t):this.store.set(e,t)}delete(e){this.assignedValues.delete(e),this.tearDownBindingObserver(e);let t=this.getRaw(e);t?this.hydrate(e,t):this.store.delete(e)}bind(){let e=t0.findParent(this);for(let t of(e&&e.appendChild(this),this.assignedValues.keys()))t.notify(this.target)}unbind(){this.parent&&tJ.get(this).removeChild(this)}appendChild(e){e.parent&&tJ.get(e).removeChild(e);let t=this.children.filter(t=>e.contains(t));for(let[i,s]of(tJ.set(e,this),this.children.push(e),t.forEach(t=>e.appendChild(t)),d.Observable.getNotifier(this.store).subscribe(e),this.store.all()))e.hydrate(i,this.bindingObservers.has(i)?this.getRaw(i):s)}removeChild(e){let t=this.children.indexOf(e);return -1!==t&&this.children.splice(t,1),d.Observable.getNotifier(this.store).unsubscribe(e),e.parent===this&&tJ.delete(e)}contains(e){return tz(this.target,e.target)}reflectToCSS(e){this.isReflecting(e)||(this.reflecting.add(e),t0.cssCustomPropertyReflector.startReflection(e,this.target))}stopReflectToCSS(e){this.isReflecting(e)&&(this.reflecting.delete(e),t0.cssCustomPropertyReflector.stopReflection(e,this.target))}isReflecting(e){return this.reflecting.has(e)}handleChange(e,t){let i=tG.getTokenById(t);i&&this.hydrate(i,this.getRaw(i))}hydrate(e,t){if(!this.has(e)){let i=this.bindingObservers.get(e);tG.isDerivedDesignTokenValue(t)?i?i.source!==t&&(this.tearDownBindingObserver(e),this.setupBindingObserver(e,t)):this.setupBindingObserver(e,t):(i&&this.tearDownBindingObserver(e),this.store.set(e,t))}}setupBindingObserver(e,t){let i=new tY(t,e,this);return this.bindingObservers.set(e,i),i}tearDownBindingObserver(e){return!!this.bindingObservers.has(e)&&(this.bindingObservers.get(e).disconnect(),this.bindingObservers.delete(e),!0)}}t0.cssCustomPropertyReflector=new class{startReflection(e,t){e.subscribe(this,t),this.handleChange({token:e,target:t})}stopReflection(e,t){e.unsubscribe(this,t),this.remove(e,t)}handleChange(e){let{token:t,target:i}=e;this.add(t,i)}add(e,t){tX.getOrCreate(t).setProperty(e.cssCustomProperty,this.resolveCSSValue(t0.getOrCreate(t).get(e)))}remove(e,t){tX.getOrCreate(t).removeProperty(e.cssCustomProperty)}resolveCSSValue(e){return e&&"function"==typeof e.createCSS?e.createCSS():e}},f([d.observable],t0.prototype,"children",void 0);let t1=Object.freeze({create:function(e){return tG.from(e)},notifyConnection:e=>!!e.isConnected&&!!t0.existsFor(e)&&(t0.getOrCreate(e).bind(),!0),notifyDisconnection:e=>!e.isConnected&&!!t0.existsFor(e)&&(t0.getOrCreate(e).unbind(),!0),registerRoot(e=tV){t_.registerRoot(e)},unregisterRoot(e=tV){t_.unregisterRoot(e)}}),t5=Object.freeze({definitionCallbackOnly:null,ignoreDuplicate:Symbol()}),t4=new Map,t8=new Map,t2=null,t3=E.createInterface(e=>e.cachedCallback(e=>(null===t2&&(t2=new t6(null,e)),t2))),t9=Object.freeze({tagFor:e=>t8.get(e),responsibleFor(e){let t=e.$$designSystem$$;return t||E.findResponsibleContainer(e).get(t3)},getOrCreate(e){if(!e)return null===t2&&(t2=E.getOrCreateDOMContainer().get(t3)),t2;let t=e.$$designSystem$$;if(t)return t;let i=E.getOrCreateDOMContainer(e);if(i.has(t3,!1))return i.get(t3);{let t=new t6(e,i);return i.register(et.instance(t3,t)),t}}});class t6{constructor(e,t){this.owner=e,this.container=t,this.designTokensInitialized=!1,this.prefix="fast",this.shadowRootMode=void 0,this.disambiguate=()=>t5.definitionCallbackOnly,null!==e&&(e.$$designSystem$$=this)}withPrefix(e){return this.prefix=e,this}withShadowRootMode(e){return this.shadowRootMode=e,this}withElementDisambiguation(e){return this.disambiguate=e,this}withDesignTokenRoot(e){return this.designTokenRoot=e,this}register(...e){let t=this.container,i=[],s=this.disambiguate,o=this.shadowRootMode,n={elementPrefix:this.prefix,tryDefineElement(e,n,a){let r="string"==typeof e?{name:e,type:n,callback:a}:e,{name:l,callback:h,baseClass:d}=r,{type:c}=r,u=l,p=t4.get(u),m=!0;for(;p;){let e=s(u,c,p);switch(e){case t5.ignoreDuplicate:return;case t5.definitionCallbackOnly:m=!1,p=void 0;break;default:u=e,p=t4.get(u)}}m&&((t8.has(c)||c===ep)&&(c=class extends c{}),t4.set(u,c),t8.set(c,u),d&&t8.set(d,u)),i.push(new t7(t,u,c,o,h,m))}};for(let s of(this.designTokensInitialized||(this.designTokensInitialized=!0,null!==this.designTokenRoot&&t1.registerRoot(this.designTokenRoot)),t.registerWithContext(n,...e),i))s.callback(s),s.willDefine&&null!==s.definition&&s.definition.define();return this}}class t7{constructor(e,t,i,s,o,n){this.container=e,this.name=t,this.type=i,this.shadowRootMode=s,this.callback=o,this.willDefine=n,this.definition=null}definePresentation(e){ec.define(this.name,e,this.container)}defineElement(e){this.definition=new d.FASTElementDefinition(this.type,Object.assign(Object.assign({},e),{name:this.name}))}tagFor(e){return t9.tagFor(e)}}let ie=(e,t)=>(0,d.html)`
    <div class="positioning-region" part="positioning-region">
        ${(0,d.when)(e=>e.modal,(0,d.html)`
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
            ${(0,d.ref)("dialog")}
        >
            <slot></slot>
        </div>
    </div>
`;var it=i(49054);class ii extends ep{constructor(){super(...arguments),this.modal=!0,this.hidden=!1,this.trapFocus=!0,this.trapFocusChanged=()=>{this.$fastController.isConnected&&this.updateTrapFocus()},this.isTrappingFocus=!1,this.handleDocumentKeydown=e=>{if(!e.defaultPrevented&&!this.hidden)switch(e.key){case"Escape":this.dismiss(),e.preventDefault();break;case"Tab":this.handleTabKeyDown(e)}},this.handleDocumentFocus=e=>{!e.defaultPrevented&&this.shouldForceFocus(e.target)&&(this.focusFirstElement(),e.preventDefault())},this.handleTabKeyDown=e=>{if(!this.trapFocus||this.hidden)return;let t=this.getTabQueueBounds();if(0!==t.length){if(1===t.length){t[0].focus(),e.preventDefault();return}e.shiftKey&&e.target===t[0]?(t[t.length-1].focus(),e.preventDefault()):e.shiftKey||e.target!==t[t.length-1]||(t[0].focus(),e.preventDefault())}},this.getTabQueueBounds=()=>ii.reduceTabbableItems([],this),this.focusFirstElement=()=>{let e=this.getTabQueueBounds();e.length>0?e[0].focus():this.dialog instanceof HTMLElement&&this.dialog.focus()},this.shouldForceFocus=e=>this.isTrappingFocus&&!this.contains(e),this.shouldTrapFocus=()=>this.trapFocus&&!this.hidden,this.updateTrapFocus=e=>{let t=void 0===e?this.shouldTrapFocus():e;t&&!this.isTrappingFocus?(this.isTrappingFocus=!0,document.addEventListener("focusin",this.handleDocumentFocus),d.DOM.queueUpdate(()=>{this.shouldForceFocus(document.activeElement)&&this.focusFirstElement()})):!t&&this.isTrappingFocus&&(this.isTrappingFocus=!1,document.removeEventListener("focusin",this.handleDocumentFocus))}}dismiss(){this.$emit("dismiss"),this.$emit("cancel")}show(){this.hidden=!1}hide(){this.hidden=!0,this.$emit("close")}connectedCallback(){super.connectedCallback(),document.addEventListener("keydown",this.handleDocumentKeydown),this.notifier=d.Observable.getNotifier(this),this.notifier.subscribe(this,"hidden"),this.updateTrapFocus()}disconnectedCallback(){super.disconnectedCallback(),document.removeEventListener("keydown",this.handleDocumentKeydown),this.updateTrapFocus(!1),this.notifier.unsubscribe(this,"hidden")}handleChange(e,t){"hidden"===t&&this.updateTrapFocus()}static reduceTabbableItems(e,t){return"-1"===t.getAttribute("tabindex")?e:(0,it.AO)(t)||ii.isFocusableFastElement(t)&&ii.hasTabbableShadow(t)?(e.push(t),e):t.childElementCount?e.concat(Array.from(t.children).reduce(ii.reduceTabbableItems,[])):e}static isFocusableFastElement(e){var t,i;return!!(null==(i=null==(t=e.$fastController)?void 0:t.definition.shadowOptions)?void 0:i.delegatesFocus)}static hasTabbableShadow(e){var t,i;return Array.from(null!=(i=null==(t=e.shadowRoot)?void 0:t.querySelectorAll("*"))?i:[]).some(e=>(0,it.AO)(e))}}f([(0,d.attr)({mode:"boolean"})],ii.prototype,"modal",void 0),f([(0,d.attr)({mode:"boolean"})],ii.prototype,"hidden",void 0),f([(0,d.attr)({attribute:"trap-focus",mode:"boolean"})],ii.prototype,"trapFocus",void 0),f([(0,d.attr)({attribute:"aria-describedby"})],ii.prototype,"ariaDescribedby",void 0),f([(0,d.attr)({attribute:"aria-labelledby"})],ii.prototype,"ariaLabelledby",void 0),f([(0,d.attr)({attribute:"aria-label"})],ii.prototype,"ariaLabel",void 0);let is=new MutationObserver(e=>{for(let t of e)io.getOrCreateFor(t.target).notify(t.attributeName)});class io extends d.SubscriberSet{constructor(e){super(e),this.watchedAttributes=new Set,io.subscriberCache.set(e,this)}subscribe(e){super.subscribe(e),this.watchedAttributes.has(e.attributes)||(this.watchedAttributes.add(e.attributes),this.observe())}unsubscribe(e){super.unsubscribe(e),this.watchedAttributes.has(e.attributes)&&(this.watchedAttributes.delete(e.attributes),this.observe())}static getOrCreateFor(e){return this.subscriberCache.get(e)||new io(e)}observe(){let e=[];for(let t of this.watchedAttributes.values())for(let i=0;i<t.length;i++)e.push(t[i]);is.observe(this.source,{attributeFilter:e})}}io.subscriberCache=new WeakMap;class ia{constructor(e,t){this.target=e,this.attributes=Object.freeze(t)}bind(e){if(io.getOrCreateFor(e).subscribe(this),e.hasAttributes())for(let t=0;t<e.attributes.length;t++)this.handleChange(e,e.attributes[t].name)}unbind(e){io.getOrCreateFor(e).unsubscribe(this)}handleChange(e,t){this.attributes.includes(t)&&d.DOM.setAttribute(this.target,t,e.getAttribute(t))}}function ir(...e){return new d.AttachedBehaviorHTMLDirective("fast-reflect-attr",ia,e)}let il=(e,t)=>(0,d.html)`
    <details class="disclosure" ${(0,d.ref)("details")}>
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
`;class ih extends ep{connectedCallback(){super.connectedCallback(),this.setup()}disconnectedCallback(){super.disconnectedCallback(),this.details.removeEventListener("toggle",this.onToggle)}show(){this.details.open=!0}hide(){this.details.open=!1}toggle(){this.details.open=!this.details.open}setup(){this.onToggle=this.onToggle.bind(this),this.details.addEventListener("toggle",this.onToggle),this.expanded&&this.show()}onToggle(){this.expanded=this.details.open,this.$emit("toggle")}}f([(0,d.attr)({mode:"boolean"})],ih.prototype,"expanded",void 0),f([d.attr],ih.prototype,"title",void 0);let id=(e,t)=>(0,d.html)`
    <template role="${e=>e.role}" aria-orientation="${e=>e.orientation}"></template>
`;var ic=i(67002);let iu={separator:"separator",presentation:"presentation"};class ip extends ep{constructor(){super(...arguments),this.role=iu.separator,this.orientation=ic.t.horizontal}}f([d.attr],ip.prototype,"role",void 0),f([d.attr],ip.prototype,"orientation",void 0);let im={next:"next",previous:"previous"},iv=(e,t)=>(0,d.html)`
    <template
        role="button"
        aria-disabled="${e=>!!e.disabled||void 0}"
        tabindex="${e=>e.hiddenFromAT?-1:0}"
        class="${e=>e.direction} ${e=>e.disabled?"disabled":""}"
        @keyup="${(e,t)=>e.keyupHandler(t.event)}"
    >
        ${(0,d.when)(e=>e.direction===im.next,(0,d.html)`
                <span part="next" class="next">
                    <slot name="next">
                        ${t.next||""}
                    </slot>
                </span>
            `)}
        ${(0,d.when)(e=>e.direction===im.previous,(0,d.html)`
                <span part="previous" class="previous">
                    <slot name="previous">
                        ${t.previous||""}
                    </slot>
                </span>
            `)}
    </template>
`;class ib extends ep{constructor(){super(...arguments),this.hiddenFromAT=!0,this.direction=im.next}keyupHandler(e){if(!this.hiddenFromAT){let t=e.key;("Enter"===t||"Space"===t)&&this.$emit("click",e),"Escape"===t&&this.blur()}}}f([(0,d.attr)({mode:"boolean"})],ib.prototype,"disabled",void 0),f([(0,d.attr)({attribute:"aria-hidden",converter:d.booleanConverter})],ib.prototype,"hiddenFromAT",void 0),f([d.attr],ib.prototype,"direction",void 0);let ig=(e,t)=>(0,d.html)`
    <template
        aria-checked="${e=>e.ariaChecked}"
        aria-disabled="${e=>e.ariaDisabled}"
        aria-posinset="${e=>e.ariaPosInSet}"
        aria-selected="${e=>e.ariaSelected}"
        aria-setsize="${e=>e.ariaSetSize}"
        class="${e=>[e.checked&&"checked",e.selected&&"selected",e.disabled&&"disabled"].filter(Boolean).join(" ")}"
        role="option"
    >
        ${p(e,t)}
        <span class="content" part="content">
            <slot ${(0,d.slotted)("content")}></slot>
        </span>
        ${u(e,t)}
    </template>
`;class iy extends tk{constructor(){super(...arguments),this.activeIndex=-1,this.rangeStartIndex=-1}get activeOption(){return this.options[this.activeIndex]}get checkedOptions(){var e;return null==(e=this.options)?void 0:e.filter(e=>e.checked)}get firstSelectedOptionIndex(){return this.options.indexOf(this.firstSelectedOption)}activeIndexChanged(e,t){var i,s;this.ariaActiveDescendant=null!=(s=null==(i=this.options[t])?void 0:i.id)?s:"",this.focusAndScrollOptionIntoView()}checkActiveIndex(){if(!this.multiple)return;let e=this.activeOption;e&&(e.checked=!0)}checkFirstOption(e=!1){e?(-1===this.rangeStartIndex&&(this.rangeStartIndex=this.activeIndex+1),this.options.forEach((e,t)=>{e.checked=(0,eC.r4)(t,this.rangeStartIndex)})):this.uncheckAllOptions(),this.activeIndex=0,this.checkActiveIndex()}checkLastOption(e=!1){e?(-1===this.rangeStartIndex&&(this.rangeStartIndex=this.activeIndex),this.options.forEach((e,t)=>{e.checked=(0,eC.r4)(t,this.rangeStartIndex,this.options.length)})):this.uncheckAllOptions(),this.activeIndex=this.options.length-1,this.checkActiveIndex()}connectedCallback(){super.connectedCallback(),this.addEventListener("focusout",this.focusoutHandler)}disconnectedCallback(){this.removeEventListener("focusout",this.focusoutHandler),super.disconnectedCallback()}checkNextOption(e=!1){e?(-1===this.rangeStartIndex&&(this.rangeStartIndex=this.activeIndex),this.options.forEach((e,t)=>{e.checked=(0,eC.r4)(t,this.rangeStartIndex,this.activeIndex+1)})):this.uncheckAllOptions(),this.activeIndex+=+(this.activeIndex<this.options.length-1),this.checkActiveIndex()}checkPreviousOption(e=!1){e?(-1===this.rangeStartIndex&&(this.rangeStartIndex=this.activeIndex),1===this.checkedOptions.length&&(this.rangeStartIndex+=1),this.options.forEach((e,t)=>{e.checked=(0,eC.r4)(t,this.activeIndex,this.rangeStartIndex)})):this.uncheckAllOptions(),this.activeIndex-=this.activeIndex>0,this.checkActiveIndex()}clickHandler(e){var t;if(!this.multiple)return super.clickHandler(e);let i=null==(t=e.target)?void 0:t.closest("[role=option]");if(i&&!i.disabled)return this.uncheckAllOptions(),this.activeIndex=this.options.indexOf(i),this.checkActiveIndex(),this.toggleSelectedForAllCheckedOptions(),!0}focusAndScrollOptionIntoView(){super.focusAndScrollOptionIntoView(this.activeOption)}focusinHandler(e){if(!this.multiple)return super.focusinHandler(e);this.shouldSkipFocus||e.target!==e.currentTarget||(this.uncheckAllOptions(),-1===this.activeIndex&&(this.activeIndex=-1!==this.firstSelectedOptionIndex?this.firstSelectedOptionIndex:0),this.checkActiveIndex(),this.setSelectedOptions(),this.focusAndScrollOptionIntoView()),this.shouldSkipFocus=!1}focusoutHandler(e){this.multiple&&this.uncheckAllOptions()}keydownHandler(e){if(!this.multiple)return super.keydownHandler(e);if(this.disabled)return!0;let{key:t,shiftKey:i}=e;switch(this.shouldSkipFocus=!1,t){case"Home":return void this.checkFirstOption(i);case ey.HX:return void this.checkNextOption(i);case ey.I5:return void this.checkPreviousOption(i);case"End":return void this.checkLastOption(i);case"Tab":return this.focusAndScrollOptionIntoView(),!0;case"Escape":return this.uncheckAllOptions(),this.checkActiveIndex(),!0;case" ":if(e.preventDefault(),this.typeAheadExpired)return void this.toggleSelectedForAllCheckedOptions();default:return 1===t.length&&this.handleTypeAhead(`${t}`),!0}}mousedownHandler(e){if(e.offsetX>=0&&e.offsetX<=this.scrollWidth)return super.mousedownHandler(e)}multipleChanged(e,t){var i;this.ariaMultiSelectable=t?"true":null,null==(i=this.options)||i.forEach(e=>{e.checked=!t&&void 0}),this.setSelectedOptions()}setSelectedOptions(){this.multiple?this.$fastController.isConnected&&this.options&&(this.selectedOptions=this.options.filter(e=>e.selected),this.focusAndScrollOptionIntoView()):super.setSelectedOptions()}sizeChanged(e,t){var i;let s=Math.max(0,parseInt(null!=(i=null==t?void 0:t.toFixed())?i:"",10));s!==t&&d.DOM.queueUpdate(()=>{this.size=s})}toggleSelectedForAllCheckedOptions(){let e=this.checkedOptions.filter(e=>!e.disabled),t=!e.every(e=>e.selected);e.forEach(e=>e.selected=t),this.selectedIndex=this.options.indexOf(e[e.length-1]),this.setSelectedOptions()}typeaheadBufferChanged(e,t){if(!this.multiple)return void super.typeaheadBufferChanged(e,t);if(this.$fastController.isConnected){let e=this.getTypeaheadMatches(),t=this.options.indexOf(e[0]);t>-1&&(this.activeIndex=t,this.uncheckAllOptions(),this.checkActiveIndex()),this.typeAheadExpired=!1}}uncheckAllOptions(e=!1){this.options.forEach(e=>e.checked=!this.multiple&&void 0),e||(this.rangeStartIndex=-1)}}f([d.observable],iy.prototype,"activeIndex",void 0),f([(0,d.attr)({mode:"boolean"})],iy.prototype,"multiple",void 0),f([(0,d.attr)({converter:d.nullableNumberConverter})],iy.prototype,"size",void 0);let iC=(e,t)=>(0,d.html)`
    <template
        aria-activedescendant="${e=>e.ariaActiveDescendant}"
        aria-multiselectable="${e=>e.ariaMultiSelectable}"
        class="listbox"
        role="listbox"
        tabindex="${e=>e.disabled?null:"0"}"
        @click="${(e,t)=>e.clickHandler(t.event)}"
        @focusin="${(e,t)=>e.focusinHandler(t.event)}"
        @keydown="${(e,t)=>e.keydownHandler(t.event)}"
        @mousedown="${(e,t)=>e.mousedownHandler(t.event)}"
    >
        <slot
            ${(0,d.slotted)({filter:iy.slottedOptionFilter,flatten:!0,property:"slottedOptions"})}
        ></slot>
    </template>
`;class ix extends ep{constructor(){super(...arguments),this.optionElements=[]}menuElementsChanged(){this.updateOptions()}headerElementsChanged(){this.updateOptions()}footerElementsChanged(){this.updateOptions()}updateOptions(){this.optionElements.splice(0,this.optionElements.length),this.addSlottedListItems(this.headerElements),this.addSlottedListItems(this.menuElements),this.addSlottedListItems(this.footerElements),this.$emit("optionsupdated",{bubbles:!1})}addSlottedListItems(e){void 0!==e&&e.forEach(e=>{1===e.nodeType&&"listitem"===e.getAttribute("role")&&(e.id=e.id||tC("option-"),this.optionElements.push(e))})}}f([d.observable],ix.prototype,"menuElements",void 0),f([d.observable],ix.prototype,"headerElements",void 0),f([d.observable],ix.prototype,"footerElements",void 0),f([d.observable],ix.prototype,"suggestionsAvailableText",void 0);let i$=(0,d.html)`
    <template>
        ${e=>e.value}
    </template>
`;class iw extends ep{contentsTemplateChanged(){this.$fastController.isConnected&&this.updateView()}connectedCallback(){super.connectedCallback(),this.updateView()}disconnectedCallback(){super.disconnectedCallback(),this.disconnectView()}handleClick(e){return!e.defaultPrevented&&(this.handleInvoked(),!1)}handleInvoked(){this.$emit("pickeroptioninvoked")}updateView(){var e,t;this.disconnectView(),this.customView=null!=(t=null==(e=this.contentsTemplate)?void 0:e.render(this,this))?t:i$.render(this,this)}disconnectView(){var e;null==(e=this.customView)||e.dispose(),this.customView=void 0}}f([(0,d.attr)({attribute:"value"})],iw.prototype,"value",void 0),f([d.observable],iw.prototype,"contentsTemplate",void 0);class iI extends ep{}let ik=(0,d.html)`
    <template>
        ${e=>e.value}
    </template>
`;class iE extends ep{contentsTemplateChanged(){this.$fastController.isConnected&&this.updateView()}connectedCallback(){super.connectedCallback(),this.updateView()}disconnectedCallback(){this.disconnectView(),super.disconnectedCallback()}handleKeyDown(e){return!e.defaultPrevented&&("Enter"!==e.key||(this.handleInvoke(),!1))}handleClick(e){return e.defaultPrevented||this.handleInvoke(),!1}handleInvoke(){this.$emit("pickeriteminvoked")}updateView(){var e,t;this.disconnectView(),this.customView=null!=(t=null==(e=this.contentsTemplate)?void 0:e.render(this,this))?t:ik.render(this,this)}disconnectView(){var e;null==(e=this.customView)||e.dispose(),this.customView=void 0}}f([(0,d.attr)({attribute:"value"})],iE.prototype,"value",void 0),f([d.observable],iE.prototype,"contentsTemplate",void 0);let iO=(e,t)=>{let i,s,o=e.tagFor(eD),n=e.tagFor(ix),a=e.tagFor(iI),r=e.tagFor(iI),l=(i=e.tagFor(iE),(0,d.html)`
    <${i}
        value="${e=>e}"
        :contentsTemplate="${(e,t)=>t.parent.listItemContentsTemplate}"
    >
    </${i}>
    `),h=(s=e.tagFor(iw),(0,d.html)`
    <${s}
        value="${e=>e}"
        :contentsTemplate="${(e,t)=>t.parent.menuOptionContentsTemplate}"
    >
    </${s}>
    `);return(0,d.html)`
        <template
            :selectedListTag="${()=>a}"
            :menuTag="${()=>n}"
            :defaultListItemTemplate="${l}"
            :defaultMenuOptionTemplate="${h}"
            @focusin="${(e,t)=>e.handleFocusIn(t.event)}"
            @focusout="${(e,t)=>e.handleFocusOut(t.event)}"
            @keydown="${(e,t)=>e.handleKeyDown(t.event)}"
            @pickeriteminvoked="${(e,t)=>e.handleItemInvoke(t.event)}"
            @pickeroptioninvoked="${(e,t)=>e.handleOptionInvoke(t.event)}"
        >
            <slot name="list-region"></slot>

            ${(0,d.when)(e=>e.flyoutOpen,(0,d.html)`
                <${o}
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
                    ${(0,d.ref)("region")}
                >
                    ${(0,d.when)(e=>!e.showNoOptions&&!e.showLoading,(0,d.html)`
                            <slot name="menu-region"></slot>
                        `)}
                    ${(0,d.when)(e=>e.showNoOptions&&!e.showLoading,(0,d.html)`
                            <div class="no-options-display" part="no-options-display">
                                <slot name="no-options-region">
                                    ${e=>e.noSuggestionsText}
                                </slot>
                            </div>
                        `)}
                    ${(0,d.when)(e=>e.showLoading,(0,d.html)`
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
                </${o}>
            `)}
        </template>
    `};class iT extends ep{}class iR extends eQ(iT){constructor(){super(...arguments),this.proxy=document.createElement("input")}}let iD=(0,d.html)`
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
        ${(0,d.ref)("inputElement")}
    ></input>
`;class iS extends iR{constructor(){super(...arguments),this.selection="",this.filterSelected=!0,this.filterQuery=!0,this.noSuggestionsText="No suggestions available",this.suggestionsAvailableText="Suggestions available",this.loadingText="Loading suggestions",this.menuPlacement="bottom-fill",this.showLoading=!1,this.optionsList=[],this.filteredOptionsList=[],this.flyoutOpen=!1,this.menuFocusIndex=-1,this.showNoOptions=!1,this.selectedItems=[],this.inputElementView=null,this.handleTextInput=e=>{this.query=this.inputElement.value},this.handleInputClick=e=>{e.preventDefault(),this.toggleFlyout(!0)},this.setRegionProps=()=>{if(this.flyoutOpen){if(null===this.region||void 0===this.region)return void d.DOM.queueUpdate(this.setRegionProps);this.region.anchorElement=this.inputElement}},this.configLookup={top:eA,bottom:eF,tallest:eL,"top-fill":eM,"bottom-fill":eP,"tallest-fill":eH}}selectionChanged(){this.$fastController.isConnected&&(this.handleSelectionChange(),this.proxy instanceof HTMLInputElement&&(this.proxy.value=this.selection,this.validate()))}optionsChanged(){this.optionsList=this.options.split(",").map(e=>e.trim()).filter(e=>""!==e)}menuPlacementChanged(){this.$fastController.isConnected&&this.updateMenuConfig()}showLoadingChanged(){this.$fastController.isConnected&&d.DOM.queueUpdate(()=>{this.setFocusedOption(0)})}listItemTemplateChanged(){this.updateListItemTemplate()}defaultListItemTemplateChanged(){this.updateListItemTemplate()}menuOptionTemplateChanged(){this.updateOptionTemplate()}defaultMenuOptionTemplateChanged(){this.updateOptionTemplate()}optionsListChanged(){this.updateFilteredOptions()}queryChanged(){this.$fastController.isConnected&&(this.inputElement.value!==this.query&&(this.inputElement.value=this.query),this.updateFilteredOptions(),this.$emit("querychange",{bubbles:!1}))}filteredOptionsListChanged(){this.$fastController.isConnected&&(this.showNoOptions=0===this.filteredOptionsList.length&&0===this.menuElement.querySelectorAll('[role="listitem"]').length,this.setFocusedOption(this.showNoOptions?-1:0))}flyoutOpenChanged(){this.flyoutOpen?(d.DOM.queueUpdate(this.setRegionProps),this.$emit("menuopening",{bubbles:!1})):this.$emit("menuclosing",{bubbles:!1})}showNoOptionsChanged(){this.$fastController.isConnected&&d.DOM.queueUpdate(()=>{this.setFocusedOption(0)})}connectedCallback(){super.connectedCallback(),this.listElement=document.createElement(this.selectedListTag),this.appendChild(this.listElement),this.itemsPlaceholderElement=document.createComment(""),this.listElement.append(this.itemsPlaceholderElement),this.inputElementView=iD.render(this,this.listElement);let e=this.menuTag.toUpperCase();this.menuElement=Array.from(this.children).find(t=>t.tagName===e),void 0===this.menuElement&&(this.menuElement=document.createElement(this.menuTag),this.appendChild(this.menuElement)),""===this.menuElement.id&&(this.menuElement.id=tC("listbox-")),this.menuId=this.menuElement.id,this.optionsPlaceholder=document.createComment(""),this.menuElement.append(this.optionsPlaceholder),this.updateMenuConfig(),d.DOM.queueUpdate(()=>this.initialize())}disconnectedCallback(){super.disconnectedCallback(),this.toggleFlyout(!1),this.inputElement.removeEventListener("input",this.handleTextInput),this.inputElement.removeEventListener("click",this.handleInputClick),null!==this.inputElementView&&(this.inputElementView.dispose(),this.inputElementView=null)}focus(){this.inputElement.focus()}initialize(){this.updateListItemTemplate(),this.updateOptionTemplate(),this.itemsRepeatBehavior=new d.RepeatDirective(e=>e.selectedItems,e=>e.activeListItemTemplate,{positioning:!0}).createBehavior(this.itemsPlaceholderElement),this.inputElement.addEventListener("input",this.handleTextInput),this.inputElement.addEventListener("click",this.handleInputClick),this.$fastController.addBehaviors([this.itemsRepeatBehavior]),this.menuElement.suggestionsAvailableText=this.suggestionsAvailableText,this.menuElement.addEventListener("optionsupdated",this.handleMenuOptionsUpdated),this.optionsRepeatBehavior=new d.RepeatDirective(e=>e.filteredOptionsList,e=>e.activeMenuOptionTemplate,{positioning:!0}).createBehavior(this.optionsPlaceholder),this.$fastController.addBehaviors([this.optionsRepeatBehavior]),this.handleSelectionChange()}toggleFlyout(e){if(this.flyoutOpen!==e){if(e&&document.activeElement===this.inputElement){this.flyoutOpen=e,d.DOM.queueUpdate(()=>{void 0!==this.menuElement?this.setFocusedOption(0):this.disableMenu()});return}this.flyoutOpen=!1,this.disableMenu()}}handleMenuOptionsUpdated(e){e.preventDefault(),this.flyoutOpen&&this.setFocusedOption(0)}handleKeyDown(e){if(e.defaultPrevented)return!1;switch(e.key){case ey.HX:if(this.flyoutOpen){let e=this.flyoutOpen?Math.min(this.menuFocusIndex+1,this.menuElement.optionElements.length-1):0;this.setFocusedOption(e)}else this.toggleFlyout(!0);return!1;case ey.I5:if(this.flyoutOpen){let e=this.flyoutOpen?Math.max(this.menuFocusIndex-1,0):0;this.setFocusedOption(e)}else this.toggleFlyout(!0);return!1;case"Escape":return this.toggleFlyout(!1),!1;case"Enter":return -1!==this.menuFocusIndex&&this.menuElement.optionElements.length>this.menuFocusIndex&&this.menuElement.optionElements[this.menuFocusIndex].click(),!1;case ey.bb:if(document.activeElement!==this.inputElement)return this.incrementFocusedItem(1),!1;return!0;case ey.kT:if(0===this.inputElement.selectionStart)return this.incrementFocusedItem(-1),!1;return!0;case"Delete":case ey.R9:{if(null===document.activeElement)return!0;if(document.activeElement===this.inputElement){if(0===this.inputElement.selectionStart)return this.selection=this.selectedItems.slice(0,this.selectedItems.length-1).toString(),this.toggleFlyout(!1),!1;return!0}let e=Array.from(this.listElement.children),t=e.indexOf(document.activeElement);if(t>-1)return this.selection=this.selectedItems.splice(t,1).toString(),d.DOM.queueUpdate(()=>{e[Math.min(e.length,t)].focus()}),!1;return!0}}return this.toggleFlyout(!0),!0}handleFocusIn(e){return!1}handleFocusOut(e){return void 0!==this.menuElement&&this.menuElement.contains(e.relatedTarget)||this.toggleFlyout(!1),!1}handleSelectionChange(){this.selectedItems.toString()!==this.selection&&(this.selectedItems=""===this.selection?[]:this.selection.split(","),this.updateFilteredOptions(),d.DOM.queueUpdate(()=>{this.checkMaxItems()}),this.$emit("selectionchange",{bubbles:!1}))}handleRegionLoaded(e){d.DOM.queueUpdate(()=>{this.setFocusedOption(0),this.$emit("menuloaded",{bubbles:!1})})}checkMaxItems(){if(void 0!==this.inputElement)if(void 0!==this.maxSelected&&this.selectedItems.length>=this.maxSelected){if(document.activeElement===this.inputElement){let e=Array.from(this.listElement.querySelectorAll("[role='listitem']"));e[e.length-1].focus()}this.inputElement.hidden=!0}else this.inputElement.hidden=!1}handleItemInvoke(e){if(e.defaultPrevented)return!1;if(e.target instanceof iE){let t=Array.from(this.listElement.querySelectorAll("[role='listitem']")).indexOf(e.target);if(-1!==t){let e=this.selectedItems.slice();e.splice(t,1),this.selection=e.toString(),d.DOM.queueUpdate(()=>this.incrementFocusedItem(0))}return!1}return!0}handleOptionInvoke(e){return!e.defaultPrevented&&(!(e.target instanceof iw)||(void 0!==e.target.value&&(this.selection=`${this.selection}${""===this.selection?"":","}${e.target.value}`),this.inputElement.value="",this.query="",this.inputElement.focus(),this.toggleFlyout(!1),!1))}incrementFocusedItem(e){if(0===this.selectedItems.length)return void this.inputElement.focus();let t=Array.from(this.listElement.querySelectorAll("[role='listitem']"));if(null!==document.activeElement){let i=t.indexOf(document.activeElement);-1===i&&(i=t.length);let s=Math.min(t.length,Math.max(0,i+e));s===t.length?void 0!==this.maxSelected&&this.selectedItems.length>=this.maxSelected?t[s-1].focus():this.inputElement.focus():t[s].focus()}}disableMenu(){var e,t,i;this.menuFocusIndex=-1,this.menuFocusOptionId=void 0,null==(e=this.inputElement)||e.removeAttribute("aria-activedescendant"),null==(t=this.inputElement)||t.removeAttribute("aria-owns"),null==(i=this.inputElement)||i.removeAttribute("aria-expanded")}setFocusedOption(e){if(!this.flyoutOpen||-1===e||this.showNoOptions||this.showLoading)return void this.disableMenu();if(0===this.menuElement.optionElements.length)return;this.menuElement.optionElements.forEach(e=>{e.setAttribute("aria-selected","false")}),this.menuFocusIndex=e,this.menuFocusIndex>this.menuElement.optionElements.length-1&&(this.menuFocusIndex=this.menuElement.optionElements.length-1),this.menuFocusOptionId=this.menuElement.optionElements[this.menuFocusIndex].id,this.inputElement.setAttribute("aria-owns",this.menuId),this.inputElement.setAttribute("aria-expanded","true"),this.inputElement.setAttribute("aria-activedescendant",this.menuFocusOptionId);let t=this.menuElement.optionElements[this.menuFocusIndex];t.setAttribute("aria-selected","true"),this.menuElement.scrollTo(0,t.offsetTop)}updateListItemTemplate(){var e;this.activeListItemTemplate=null!=(e=this.listItemTemplate)?e:this.defaultListItemTemplate}updateOptionTemplate(){var e;this.activeMenuOptionTemplate=null!=(e=this.menuOptionTemplate)?e:this.defaultMenuOptionTemplate}updateFilteredOptions(){this.filteredOptionsList=this.optionsList.slice(0),this.filterSelected&&(this.filteredOptionsList=this.filteredOptionsList.filter(e=>-1===this.selectedItems.indexOf(e))),this.filterQuery&&""!==this.query&&void 0!==this.query&&(this.filteredOptionsList=this.filteredOptionsList.filter(e=>-1!==e.indexOf(this.query)))}updateMenuConfig(){let e=this.configLookup[this.menuPlacement];null===e&&(e=eP),this.menuConfig=Object.assign(Object.assign({},e),{autoUpdateMode:"auto",fixedPlacement:!0,horizontalViewportLock:!1,verticalViewportLock:!1})}}f([(0,d.attr)({attribute:"selection"})],iS.prototype,"selection",void 0),f([(0,d.attr)({attribute:"options"})],iS.prototype,"options",void 0),f([(0,d.attr)({attribute:"filter-selected",mode:"boolean"})],iS.prototype,"filterSelected",void 0),f([(0,d.attr)({attribute:"filter-query",mode:"boolean"})],iS.prototype,"filterQuery",void 0),f([(0,d.attr)({attribute:"max-selected"})],iS.prototype,"maxSelected",void 0),f([(0,d.attr)({attribute:"no-suggestions-text"})],iS.prototype,"noSuggestionsText",void 0),f([(0,d.attr)({attribute:"suggestions-available-text"})],iS.prototype,"suggestionsAvailableText",void 0),f([(0,d.attr)({attribute:"loading-text"})],iS.prototype,"loadingText",void 0),f([(0,d.attr)({attribute:"label"})],iS.prototype,"label",void 0),f([(0,d.attr)({attribute:"labelledby"})],iS.prototype,"labelledBy",void 0),f([(0,d.attr)({attribute:"placeholder"})],iS.prototype,"placeholder",void 0),f([(0,d.attr)({attribute:"menu-placement"})],iS.prototype,"menuPlacement",void 0),f([d.observable],iS.prototype,"showLoading",void 0),f([d.observable],iS.prototype,"listItemTemplate",void 0),f([d.observable],iS.prototype,"defaultListItemTemplate",void 0),f([d.observable],iS.prototype,"activeListItemTemplate",void 0),f([d.observable],iS.prototype,"menuOptionTemplate",void 0),f([d.observable],iS.prototype,"defaultMenuOptionTemplate",void 0),f([d.observable],iS.prototype,"activeMenuOptionTemplate",void 0),f([d.observable],iS.prototype,"listItemContentsTemplate",void 0),f([d.observable],iS.prototype,"menuOptionContentsTemplate",void 0),f([d.observable],iS.prototype,"optionsList",void 0),f([d.observable],iS.prototype,"query",void 0),f([d.observable],iS.prototype,"filteredOptionsList",void 0),f([d.observable],iS.prototype,"flyoutOpen",void 0),f([d.observable],iS.prototype,"menuId",void 0),f([d.observable],iS.prototype,"selectedListTag",void 0),f([d.observable],iS.prototype,"menuTag",void 0),f([d.observable],iS.prototype,"menuFocusIndex",void 0),f([d.observable],iS.prototype,"menuFocusOptionId",void 0),f([d.observable],iS.prototype,"showNoOptions",void 0),f([d.observable],iS.prototype,"menuConfig",void 0),f([d.observable],iS.prototype,"selectedItems",void 0);let iA=(e,t)=>(0,d.html)`
        <template role="list" slot="menu-region">
            <div class="options-display" part="options-display">
                <div class="header-region" part="header-region">
                    <slot name="header-region" ${(0,d.slotted)("headerElements")}></slot>
                </div>

                <slot ${(0,d.slotted)("menuElements")}></slot>
                <div class="footer-region" part="footer-region">
                    <slot name="footer-region" ${(0,d.slotted)("footerElements")}></slot>
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
    `,iF=(e,t)=>(0,d.html)`
        <template
            role="listitem"
            tabindex="-1"
            @click="${(e,t)=>e.handleClick(t.event)}"
        >
            <slot></slot>
        </template>
    `,iL=(e,t)=>(0,d.html)`
        <template slot="list-region" role="list" class="picker-list">
            <slot></slot>
            <slot name="input-region"></slot>
        </template>
    `,iM=(e,t)=>(0,d.html)`
        <template
            role="listitem"
            tabindex="0"
            @click="${(e,t)=>e.handleClick(t.event)}"
            @keydown="${(e,t)=>e.handleKeyDown(t.event)}"
        >
            <slot></slot>
        </template>
    `,iP={menuitem:"menuitem",menuitemcheckbox:"menuitemcheckbox",menuitemradio:"menuitemradio"},iH={[iP.menuitem]:"menuitem",[iP.menuitemcheckbox]:"menuitemcheckbox",[iP.menuitemradio]:"menuitemradio"},iz=(e,t)=>(0,d.html)`
    <template
        role="${e=>e.role}"
        aria-haspopup="${e=>e.hasSubmenu?"menu":void 0}"
        aria-checked="${e=>e.role!==iP.menuitem?e.checked:void 0}"
        aria-disabled="${e=>e.disabled}"
        aria-expanded="${e=>e.expanded}"
        @keydown="${(e,t)=>e.handleMenuItemKeyDown(t.event)}"
        @click="${(e,t)=>e.handleMenuItemClick(t.event)}"
        @mouseover="${(e,t)=>e.handleMouseOver(t.event)}"
        @mouseout="${(e,t)=>e.handleMouseOut(t.event)}"
        class="${e=>e.disabled?"disabled":""} ${e=>e.expanded?"expanded":""} ${e=>`indent-${e.startColumnCount}`}"
    >
            ${(0,d.when)(e=>e.role===iP.menuitemcheckbox,(0,d.html)`
                    <div part="input-container" class="input-container">
                        <span part="checkbox" class="checkbox">
                            <slot name="checkbox-indicator">
                                ${t.checkboxIndicator||""}
                            </slot>
                        </span>
                    </div>
                `)}
            ${(0,d.when)(e=>e.role===iP.menuitemradio,(0,d.html)`
                    <div part="input-container" class="input-container">
                        <span part="radio" class="radio">
                            <slot name="radio-indicator">
                                ${t.radioIndicator||""}
                            </slot>
                        </span>
                    </div>
                `)}
        </div>
        ${p(e,t)}
        <span class="content" part="content">
            <slot></slot>
        </span>
        ${u(e,t)}
        ${(0,d.when)(e=>e.hasSubmenu,(0,d.html)`
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
        ${(0,d.when)(e=>e.expanded,(0,d.html)`
                <${e.tagFor(eD)}
                    :anchorElement="${e=>e}"
                    vertical-positioning-mode="dynamic"
                    vertical-default-position="bottom"
                    vertical-inset="true"
                    horizontal-positioning-mode="dynamic"
                    horizontal-default-position="end"
                    class="submenu-region"
                    dir="${e=>e.currentDirection}"
                    @loaded="${e=>e.submenuLoaded()}"
                    ${(0,d.ref)("submenuRegion")}
                    part="submenu-region"
                >
                    <slot name="submenu"></slot>
                </${e.tagFor(eD)}>
            `)}
    </template>
`;class iV extends ep{constructor(){super(...arguments),this.role=iP.menuitem,this.hasSubmenu=!1,this.currentDirection=eT.O.ltr,this.focusSubmenuOnLoad=!1,this.handleMenuItemKeyDown=e=>{if(e.defaultPrevented)return!1;switch(e.key){case"Enter":case" ":return this.invoke(),!1;case ey.bb:return this.expandAndFocus(),!1;case ey.kT:if(this.expanded)return this.expanded=!1,this.focus(),!1}return!0},this.handleMenuItemClick=e=>!e.defaultPrevented&&!this.disabled&&(this.invoke(),!1),this.submenuLoaded=()=>{this.focusSubmenuOnLoad&&(this.focusSubmenuOnLoad=!1,this.hasSubmenu&&(this.submenu.focus(),this.setAttribute("tabindex","-1")))},this.handleMouseOver=e=>!this.disabled&&!!this.hasSubmenu&&!this.expanded&&(this.expanded=!0,!1),this.handleMouseOut=e=>!(!this.expanded||this.contains(document.activeElement))&&(this.expanded=!1,!1),this.expandAndFocus=()=>{this.hasSubmenu&&(this.focusSubmenuOnLoad=!0,this.expanded=!0)},this.invoke=()=>{if(!this.disabled)switch(this.role){case iP.menuitemcheckbox:this.checked=!this.checked;break;case iP.menuitem:this.updateSubmenu(),this.hasSubmenu?this.expandAndFocus():this.$emit("change");break;case iP.menuitemradio:this.checked||(this.checked=!0)}},this.updateSubmenu=()=>{this.submenu=this.domChildren().find(e=>"menu"===e.getAttribute("role")),this.hasSubmenu=void 0!==this.submenu}}expandedChanged(e){this.$fastController.isConnected&&void 0!==this.submenu&&(!1===this.expanded?this.submenu.collapseExpandedItem():this.currentDirection=eR(this),this.$emit("expanded-change",this,{bubbles:!1}))}checkedChanged(e,t){this.$fastController.isConnected&&this.$emit("change")}connectedCallback(){super.connectedCallback(),d.DOM.queueUpdate(()=>{this.updateSubmenu()}),this.startColumnCount||(this.startColumnCount=1),this.observer=new MutationObserver(this.updateSubmenu)}disconnectedCallback(){super.disconnectedCallback(),this.submenu=void 0,void 0!==this.observer&&(this.observer.disconnect(),this.observer=void 0)}domChildren(){return Array.from(this.children).filter(e=>!e.hasAttribute("hidden"))}}f([(0,d.attr)({mode:"boolean"})],iV.prototype,"disabled",void 0),f([(0,d.attr)({mode:"boolean"})],iV.prototype,"expanded",void 0),f([d.observable],iV.prototype,"startColumnCount",void 0),f([d.attr],iV.prototype,"role",void 0),f([(0,d.attr)({mode:"boolean"})],iV.prototype,"checked",void 0),f([d.observable],iV.prototype,"submenuRegion",void 0),f([d.observable],iV.prototype,"hasSubmenu",void 0),f([d.observable],iV.prototype,"currentDirection",void 0),f([d.observable],iV.prototype,"submenu",void 0),eb(iV,c);let iN=(e,t)=>(0,d.html)`
    <template
        slot="${e=>e.slot?e.slot:e.isNestedMenu()?"submenu":void 0}"
        role="menu"
        @keydown="${(e,t)=>e.handleMenuKeyDown(t.event)}"
        @focusout="${(e,t)=>e.handleFocusOut(t.event)}"
    >
        <slot ${(0,d.slotted)("items")}></slot>
    </template>
`;class iB extends ep{constructor(){super(...arguments),this.expandedItem=null,this.focusIndex=-1,this.isNestedMenu=()=>null!==this.parentElement&&tx(this.parentElement)&&"menuitem"===this.parentElement.getAttribute("role"),this.handleFocusOut=e=>{if(!this.contains(e.relatedTarget)&&void 0!==this.menuItems){this.collapseExpandedItem();let e=this.menuItems.findIndex(this.isFocusableElement);this.menuItems[this.focusIndex].setAttribute("tabindex","-1"),this.menuItems[e].setAttribute("tabindex","0"),this.focusIndex=e}},this.handleItemFocus=e=>{let t=e.target;void 0!==this.menuItems&&t!==this.menuItems[this.focusIndex]&&(this.menuItems[this.focusIndex].setAttribute("tabindex","-1"),this.focusIndex=this.menuItems.indexOf(t),t.setAttribute("tabindex","0"))},this.handleExpandedChanged=e=>{if(e.defaultPrevented||null===e.target||void 0===this.menuItems||0>this.menuItems.indexOf(e.target))return;e.preventDefault();let t=e.target;if(null!==this.expandedItem&&t===this.expandedItem&&!1===t.expanded){this.expandedItem=null;return}t.expanded&&(null!==this.expandedItem&&this.expandedItem!==t&&(this.expandedItem.expanded=!1),this.menuItems[this.focusIndex].setAttribute("tabindex","-1"),this.expandedItem=t,this.focusIndex=this.menuItems.indexOf(t),t.setAttribute("tabindex","0"))},this.removeItemListeners=()=>{void 0!==this.menuItems&&this.menuItems.forEach(e=>{e.removeEventListener("expanded-change",this.handleExpandedChanged),e.removeEventListener("focus",this.handleItemFocus)})},this.setItems=()=>{let e=this.domChildren();this.removeItemListeners(),this.menuItems=e;let t=this.menuItems.filter(this.isMenuItemElement);t.length&&(this.focusIndex=0);let i=t.reduce((e,t)=>{var i;let s,o,n=(s=(i=t).getAttribute("role"),o=i.querySelector("[slot=start]"),s!==iP.menuitem&&null===o||s===iP.menuitem&&null!==o?1:2*(s!==iP.menuitem&&null!==o));return e>n?e:n},0);t.forEach((e,t)=>{e.setAttribute("tabindex",0===t?"0":"-1"),e.addEventListener("expanded-change",this.handleExpandedChanged),e.addEventListener("focus",this.handleItemFocus),e instanceof iV&&(e.startColumnCount=i)})},this.changeHandler=e=>{if(void 0===this.menuItems)return;let t=e.target,i=this.menuItems.indexOf(t);if(-1!==i&&"menuitemradio"===t.role&&!0===t.checked){for(let e=i-1;e>=0;--e){let t=this.menuItems[e],i=t.getAttribute("role");if(i===iP.menuitemradio&&(t.checked=!1),"separator"===i)break}let e=this.menuItems.length-1;for(let t=i+1;t<=e;++t){let e=this.menuItems[t],i=e.getAttribute("role");if(i===iP.menuitemradio&&(e.checked=!1),"separator"===i)break}}},this.isMenuItemElement=e=>tx(e)&&iB.focusableElementRoles.hasOwnProperty(e.getAttribute("role")),this.isFocusableElement=e=>this.isMenuItemElement(e)}itemsChanged(e,t){this.$fastController.isConnected&&void 0!==this.menuItems&&this.setItems()}connectedCallback(){super.connectedCallback(),d.DOM.queueUpdate(()=>{this.setItems()}),this.addEventListener("change",this.changeHandler)}disconnectedCallback(){super.disconnectedCallback(),this.removeItemListeners(),this.menuItems=void 0,this.removeEventListener("change",this.changeHandler)}focus(){this.setFocus(0,1)}collapseExpandedItem(){null!==this.expandedItem&&(this.expandedItem.expanded=!1,this.expandedItem=null)}handleMenuKeyDown(e){if(!e.defaultPrevented&&void 0!==this.menuItems)switch(e.key){case ey.HX:this.setFocus(this.focusIndex+1,1);return;case ey.I5:this.setFocus(this.focusIndex-1,-1);return;case"End":this.setFocus(this.menuItems.length-1,-1);return;case"Home":this.setFocus(0,1);return;default:return!0}}domChildren(){return Array.from(this.children).filter(e=>!e.hasAttribute("hidden"))}setFocus(e,t){if(void 0!==this.menuItems)for(;e>=0&&e<this.menuItems.length;){let i=this.menuItems[e];if(this.isFocusableElement(i)){this.focusIndex>-1&&this.menuItems.length>=this.focusIndex-1&&this.menuItems[this.focusIndex].setAttribute("tabindex","-1"),this.focusIndex=e,i.setAttribute("tabindex","0"),i.focus();break}e+=t}}}iB.focusableElementRoles=iH,f([d.observable],iB.prototype,"items",void 0);let iq=(e,t)=>(0,d.html)`
    <template class="${e=>e.readOnly?"readonly":""}">
        <label
            part="label"
            for="control"
            class="${e=>e.defaultSlottedNodes&&e.defaultSlottedNodes.length?"label":"label label__hidden"}"
        >
            <slot ${(0,d.slotted)("defaultSlottedNodes")}></slot>
        </label>
        <div class="root" part="root">
            ${p(e,t)}
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
                ${(0,d.ref)("control")}
            />
            ${(0,d.when)(e=>!e.hideStep&&!e.readOnly&&!e.disabled,(0,d.html)`
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
            ${u(e,t)}
        </div>
    </template>
`;class iU extends ep{}class ij extends eQ(iU){constructor(){super(...arguments),this.proxy=document.createElement("input")}}let i_={email:"email",password:"password",tel:"tel",text:"text",url:"url"};class iK extends ij{constructor(){super(...arguments),this.type=i_.text}readOnlyChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.readOnly=this.readOnly,this.validate())}autofocusChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.autofocus=this.autofocus,this.validate())}placeholderChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.placeholder=this.placeholder)}typeChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.type=this.type,this.validate())}listChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.setAttribute("list",this.list),this.validate())}maxlengthChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.maxLength=this.maxlength,this.validate())}minlengthChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.minLength=this.minlength,this.validate())}patternChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.pattern=this.pattern,this.validate())}sizeChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.size=this.size)}spellcheckChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.spellcheck=this.spellcheck)}connectedCallback(){super.connectedCallback(),this.proxy.setAttribute("type",this.type),this.validate(),this.autofocus&&d.DOM.queueUpdate(()=>{this.focus()})}select(){this.control.select(),this.$emit("select")}handleTextInput(){this.value=this.control.value}handleChange(){this.$emit("change")}validate(){super.validate(this.control)}}f([(0,d.attr)({attribute:"readonly",mode:"boolean"})],iK.prototype,"readOnly",void 0),f([(0,d.attr)({mode:"boolean"})],iK.prototype,"autofocus",void 0),f([d.attr],iK.prototype,"placeholder",void 0),f([d.attr],iK.prototype,"type",void 0),f([d.attr],iK.prototype,"list",void 0),f([(0,d.attr)({converter:d.nullableNumberConverter})],iK.prototype,"maxlength",void 0),f([(0,d.attr)({converter:d.nullableNumberConverter})],iK.prototype,"minlength",void 0),f([d.attr],iK.prototype,"pattern",void 0),f([(0,d.attr)({converter:d.nullableNumberConverter})],iK.prototype,"size",void 0),f([(0,d.attr)({mode:"boolean"})],iK.prototype,"spellcheck",void 0),f([d.observable],iK.prototype,"defaultSlottedNodes",void 0);class iW{}eb(iW,eI),eb(iK,c,iW);class iX extends ep{}class iG extends eQ(iX){constructor(){super(...arguments),this.proxy=document.createElement("input")}}class iY extends iG{constructor(){super(...arguments),this.hideStep=!1,this.step=1,this.isUserInput=!1}maxChanged(e,t){var i;this.max=Math.max(t,null!=(i=this.min)?i:t);let s=Math.min(this.min,this.max);void 0!==this.min&&this.min!==s&&(this.min=s),this.value=this.getValidValue(this.value)}minChanged(e,t){var i;this.min=Math.min(t,null!=(i=this.max)?i:t);let s=Math.max(this.min,this.max);void 0!==this.max&&this.max!==s&&(this.max=s),this.value=this.getValidValue(this.value)}get valueAsNumber(){return parseFloat(super.value)}set valueAsNumber(e){this.value=e.toString()}valueChanged(e,t){this.value=this.getValidValue(t),t===this.value&&(this.control&&!this.isUserInput&&(this.control.value=this.value),super.valueChanged(e,this.value),void 0===e||this.isUserInput||(this.$emit("input"),this.$emit("change")),this.isUserInput=!1)}validate(){super.validate(this.control)}getValidValue(e){var t,i;let s=parseFloat(parseFloat(e).toPrecision(12));return isNaN(s)?"":Math.max(s=Math.min(s,null!=(t=this.max)?t:s),null!=(i=this.min)?i:s).toString()}stepUp(){let e=parseFloat(this.value),t=isNaN(e)?this.min>0?this.min:this.max<0?this.max:this.min?0:this.step:e+this.step;this.value=t.toString()}stepDown(){let e=parseFloat(this.value),t=isNaN(e)?this.min>0?this.min:this.max<0?this.max:this.min?0:0-this.step:e-this.step;this.value=t.toString()}connectedCallback(){super.connectedCallback(),this.proxy.setAttribute("type","number"),this.validate(),this.control.value=this.value,this.autofocus&&d.DOM.queueUpdate(()=>{this.focus()})}select(){this.control.select(),this.$emit("select")}handleTextInput(){this.control.value=this.control.value.replace(/[^0-9\-+e.]/g,""),this.isUserInput=!0,this.value=this.control.value}handleChange(){this.$emit("change")}handleKeyDown(e){switch(e.key){case ey.I5:return this.stepUp(),!1;case ey.HX:return this.stepDown(),!1}return!0}handleBlur(){this.control.value=this.value}}f([(0,d.attr)({attribute:"readonly",mode:"boolean"})],iY.prototype,"readOnly",void 0),f([(0,d.attr)({mode:"boolean"})],iY.prototype,"autofocus",void 0),f([(0,d.attr)({attribute:"hide-step",mode:"boolean"})],iY.prototype,"hideStep",void 0),f([d.attr],iY.prototype,"placeholder",void 0),f([d.attr],iY.prototype,"list",void 0),f([(0,d.attr)({converter:d.nullableNumberConverter})],iY.prototype,"maxlength",void 0),f([(0,d.attr)({converter:d.nullableNumberConverter})],iY.prototype,"minlength",void 0),f([(0,d.attr)({converter:d.nullableNumberConverter})],iY.prototype,"size",void 0),f([(0,d.attr)({converter:d.nullableNumberConverter})],iY.prototype,"step",void 0),f([(0,d.attr)({converter:d.nullableNumberConverter})],iY.prototype,"max",void 0),f([(0,d.attr)({converter:d.nullableNumberConverter})],iY.prototype,"min",void 0),f([d.observable],iY.prototype,"defaultSlottedNodes",void 0),eb(iY,c,iW);let iQ=(e,t)=>(0,d.html)`
    <template
        role="progressbar"
        aria-valuenow="${e=>e.value}"
        aria-valuemin="${e=>e.min}"
        aria-valuemax="${e=>e.max}"
        class="${e=>e.paused?"paused":""}"
    >
        ${(0,d.when)(e=>"number"==typeof e.value,(0,d.html)`
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
                        style="stroke-dasharray: ${e=>44*e.percentComplete/100}px ${44}px"
                        cx="8px"
                        cy="8px"
                        r="7px"
                    ></circle>
                </svg>
            `,(0,d.html)`
                <slot name="indeterminate" slot="indeterminate">
                    ${t.indeterminateIndicator||""}
                </slot>
            `)}
    </template>
`;class iZ extends ep{constructor(){super(...arguments),this.percentComplete=0}valueChanged(){this.$fastController.isConnected&&this.updatePercentComplete()}minChanged(){this.$fastController.isConnected&&this.updatePercentComplete()}maxChanged(){this.$fastController.isConnected&&this.updatePercentComplete()}connectedCallback(){super.connectedCallback(),this.updatePercentComplete()}updatePercentComplete(){let e="number"==typeof this.min?this.min:0,t="number"==typeof this.max?this.max:100,i="number"==typeof this.value?this.value:0,s=t-e;this.percentComplete=0===s?0:Math.fround((i-e)/s*100)}}f([(0,d.attr)({converter:d.nullableNumberConverter})],iZ.prototype,"value",void 0),f([(0,d.attr)({converter:d.nullableNumberConverter})],iZ.prototype,"min",void 0),f([(0,d.attr)({converter:d.nullableNumberConverter})],iZ.prototype,"max",void 0),f([(0,d.attr)({mode:"boolean"})],iZ.prototype,"paused",void 0),f([d.observable],iZ.prototype,"percentComplete",void 0);let iJ=(e,t)=>(0,d.html)`
    <template
        role="progressbar"
        aria-valuenow="${e=>e.value}"
        aria-valuemin="${e=>e.min}"
        aria-valuemax="${e=>e.max}"
        class="${e=>e.paused?"paused":""}"
    >
        ${(0,d.when)(e=>"number"==typeof e.value,(0,d.html)`
                <div class="progress" part="progress" slot="determinate">
                    <div
                        class="determinate"
                        part="determinate"
                        style="width: ${e=>e.percentComplete}%"
                    ></div>
                </div>
            `,(0,d.html)`
                <div class="progress" part="progress" slot="indeterminate">
                    <slot class="indeterminate" name="indeterminate">
                        ${t.indeterminateIndicator1||""}
                        ${t.indeterminateIndicator2||""}
                    </slot>
                </div>
            `)}
    </template>
`,i0=(e,t)=>(0,d.html)`
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
            class="positioning-region ${e=>e.orientation===ic.t.horizontal?"horizontal":"vertical"}"
            part="positioning-region"
        >
            <slot
                ${(0,d.slotted)({property:"slottedRadioButtons",filter:(0,d.elements)("[role=radio]")})}
            ></slot>
        </div>
    </template>
`;class i1 extends ep{constructor(){super(...arguments),this.orientation=ic.t.horizontal,this.radioChangeHandler=e=>{let t=e.target;t.checked&&(this.slottedRadioButtons.forEach(e=>{e!==t&&(e.checked=!1,this.isInsideFoundationToolbar||e.setAttribute("tabindex","-1"))}),this.selectedRadio=t,this.value=t.value,t.setAttribute("tabindex","0"),this.focusedRadio=t),e.stopPropagation()},this.moveToRadioByIndex=(e,t)=>{let i=e[t];this.isInsideToolbar||(i.setAttribute("tabindex","0"),i.readOnly?this.slottedRadioButtons.forEach(e=>{e!==i&&e.setAttribute("tabindex","-1")}):(i.checked=!0,this.selectedRadio=i)),this.focusedRadio=i,i.focus()},this.moveRightOffGroup=()=>{var e;null==(e=this.nextElementSibling)||e.focus()},this.moveLeftOffGroup=()=>{var e;null==(e=this.previousElementSibling)||e.focus()},this.focusOutHandler=e=>{let t=this.slottedRadioButtons,i=e.target,s=null!==i?t.indexOf(i):0,o=this.focusedRadio?t.indexOf(this.focusedRadio):-1;return(0===o&&s===o||o===t.length-1&&o===s)&&(this.selectedRadio?(this.focusedRadio=this.selectedRadio,this.isInsideFoundationToolbar||(this.selectedRadio.setAttribute("tabindex","0"),t.forEach(e=>{e!==this.selectedRadio&&e.setAttribute("tabindex","-1")}))):(this.focusedRadio=t[0],this.focusedRadio.setAttribute("tabindex","0"),t.forEach(e=>{e!==this.focusedRadio&&e.setAttribute("tabindex","-1")}))),!0},this.clickHandler=e=>{let t=e.target;if(t){let e=this.slottedRadioButtons;t.checked||0===e.indexOf(t)?(t.setAttribute("tabindex","0"),this.selectedRadio=t):(t.setAttribute("tabindex","-1"),this.selectedRadio=null),this.focusedRadio=t}e.preventDefault()},this.shouldMoveOffGroupToTheRight=(e,t,i)=>e===t.length&&this.isInsideToolbar&&i===ey.bb,this.shouldMoveOffGroupToTheLeft=(e,t)=>(this.focusedRadio?e.indexOf(this.focusedRadio)-1:0)<0&&this.isInsideToolbar&&t===ey.kT,this.checkFocusedRadio=()=>{null===this.focusedRadio||this.focusedRadio.readOnly||this.focusedRadio.checked||(this.focusedRadio.checked=!0,this.focusedRadio.setAttribute("tabindex","0"),this.focusedRadio.focus(),this.selectedRadio=this.focusedRadio)},this.moveRight=e=>{let t=this.slottedRadioButtons,i=0;if(i=this.focusedRadio?t.indexOf(this.focusedRadio)+1:1,this.shouldMoveOffGroupToTheRight(i,t,e.key))return void this.moveRightOffGroup();for(i===t.length&&(i=0);i<t.length&&t.length>1;)if(t[i].disabled)if(this.focusedRadio&&i===t.indexOf(this.focusedRadio))break;else if(i+1>=t.length)if(this.isInsideToolbar)break;else i=0;else i+=1;else{this.moveToRadioByIndex(t,i);break}},this.moveLeft=e=>{let t=this.slottedRadioButtons,i=0;if(i=(i=this.focusedRadio?t.indexOf(this.focusedRadio)-1:0)<0?t.length-1:i,this.shouldMoveOffGroupToTheLeft(t,e.key))return void this.moveLeftOffGroup();for(;i>=0&&t.length>1;)if(t[i].disabled)if(this.focusedRadio&&i===t.indexOf(this.focusedRadio))break;else i-1<0?i=t.length-1:i-=1;else{this.moveToRadioByIndex(t,i);break}},this.keydownHandler=e=>{let t=e.key;if(t in ey.Is&&this.isInsideFoundationToolbar)return!0;switch(t){case"Enter":this.checkFocusedRadio();break;case ey.bb:case ey.HX:this.direction===eT.O.ltr?this.moveRight(e):this.moveLeft(e);break;case ey.kT:case ey.I5:this.direction===eT.O.ltr?this.moveLeft(e):this.moveRight(e);break;default:return!0}}}readOnlyChanged(){void 0!==this.slottedRadioButtons&&this.slottedRadioButtons.forEach(e=>{this.readOnly?e.readOnly=!0:e.readOnly=!1})}disabledChanged(){void 0!==this.slottedRadioButtons&&this.slottedRadioButtons.forEach(e=>{this.disabled?e.disabled=!0:e.disabled=!1})}nameChanged(){this.slottedRadioButtons&&this.slottedRadioButtons.forEach(e=>{e.setAttribute("name",this.name)})}valueChanged(){this.slottedRadioButtons&&this.slottedRadioButtons.forEach(e=>{e.value===this.value&&(e.checked=!0,this.selectedRadio=e)}),this.$emit("change")}slottedRadioButtonsChanged(e,t){this.slottedRadioButtons&&this.slottedRadioButtons.length>0&&this.setupRadioButtons()}get parentToolbar(){return this.closest('[role="toolbar"]')}get isInsideToolbar(){var e;return null!=(e=this.parentToolbar)&&e}get isInsideFoundationToolbar(){var e;return!!(null==(e=this.parentToolbar)?void 0:e.$fastController)}connectedCallback(){super.connectedCallback(),this.direction=eR(this),this.setupRadioButtons()}disconnectedCallback(){this.slottedRadioButtons.forEach(e=>{e.removeEventListener("change",this.radioChangeHandler)})}setupRadioButtons(){let e=this.slottedRadioButtons.filter(e=>e.hasAttribute("checked")),t=e?e.length:0;t>1&&(e[t-1].checked=!0);let i=!1;if(this.slottedRadioButtons.forEach(e=>{void 0!==this.name&&e.setAttribute("name",this.name),this.disabled&&(e.disabled=!0),this.readOnly&&(e.readOnly=!0),this.value&&this.value===e.value?(this.selectedRadio=e,this.focusedRadio=e,e.checked=!0,e.setAttribute("tabindex","0"),i=!0):(this.isInsideFoundationToolbar||e.setAttribute("tabindex","-1"),e.checked=!1),e.addEventListener("change",this.radioChangeHandler)}),void 0===this.value&&this.slottedRadioButtons.length>0){let e=this.slottedRadioButtons.filter(e=>e.hasAttribute("checked")),t=null!==e?e.length:0;if(t>0&&!i){let i=e[t-1];i.checked=!0,this.focusedRadio=i,i.setAttribute("tabindex","0")}else this.slottedRadioButtons[0].setAttribute("tabindex","0"),this.focusedRadio=this.slottedRadioButtons[0]}}}f([(0,d.attr)({attribute:"readonly",mode:"boolean"})],i1.prototype,"readOnly",void 0),f([(0,d.attr)({attribute:"disabled",mode:"boolean"})],i1.prototype,"disabled",void 0),f([d.attr],i1.prototype,"name",void 0),f([d.attr],i1.prototype,"value",void 0),f([d.attr],i1.prototype,"orientation",void 0),f([d.observable],i1.prototype,"childItems",void 0),f([d.observable],i1.prototype,"slottedRadioButtons",void 0);let i5=(e,t)=>(0,d.html)`
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
            <slot ${(0,d.slotted)("defaultSlottedNodes")}></slot>
        </label>
    </template>
`;class i4 extends ep{}class i8 extends eZ(i4){constructor(){super(...arguments),this.proxy=document.createElement("input")}}class i2 extends i8{constructor(){super(),this.initialValue="on",this.keypressHandler=e=>{if(" "===e.key){this.checked||this.readOnly||(this.checked=!0);return}return!0},this.proxy.setAttribute("type","radio")}readOnlyChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.readOnly=this.readOnly)}defaultCheckedChanged(){var e;!this.$fastController.isConnected||this.dirtyChecked||this.isInsideRadioGroup()||(this.checked=null!=(e=this.defaultChecked)&&e,this.dirtyChecked=!1)}connectedCallback(){var e,t;super.connectedCallback(),this.validate(),(null==(e=this.parentElement)?void 0:e.getAttribute("role"))==="radiogroup"||null!==this.getAttribute("tabindex")||this.disabled||this.setAttribute("tabindex","0"),!this.checkedAttribute||this.dirtyChecked||this.isInsideRadioGroup()||(this.checked=null!=(t=this.defaultChecked)&&t,this.dirtyChecked=!1)}isInsideRadioGroup(){return null!==this.closest("[role=radiogroup]")}clickHandler(e){this.disabled||this.readOnly||this.checked||(this.checked=!0)}}f([(0,d.attr)({attribute:"readonly",mode:"boolean"})],i2.prototype,"readOnly",void 0),f([d.observable],i2.prototype,"name",void 0),f([d.observable],i2.prototype,"defaultSlottedNodes",void 0);class i3 extends ep{constructor(){super(...arguments),this.framesPerSecond=60,this.updatingItems=!1,this.speed=600,this.easing="ease-in-out",this.flippersHiddenFromAT=!1,this.scrolling=!1,this.resizeDetector=null}get frameTime(){return 1e3/this.framesPerSecond}scrollingChanged(e,t){if(this.scrollContainer){let e=!0==this.scrolling?"scrollstart":"scrollend";this.$emit(e,this.scrollContainer.scrollLeft)}}get isRtl(){return this.scrollItems.length>1&&this.scrollItems[0].offsetLeft>this.scrollItems[1].offsetLeft}connectedCallback(){super.connectedCallback(),this.initializeResizeDetector()}disconnectedCallback(){this.disconnectResizeDetector(),super.disconnectedCallback()}scrollItemsChanged(e,t){t&&!this.updatingItems&&d.DOM.queueUpdate(()=>this.setStops())}disconnectResizeDetector(){this.resizeDetector&&(this.resizeDetector.disconnect(),this.resizeDetector=null)}initializeResizeDetector(){this.disconnectResizeDetector(),this.resizeDetector=new window.ResizeObserver(this.resized.bind(this)),this.resizeDetector.observe(this)}updateScrollStops(){this.updatingItems=!0;let e=this.scrollItems.reduce((e,t)=>t instanceof HTMLSlotElement?e.concat(t.assignedElements()):(e.push(t),e),[]);this.scrollItems=e,this.updatingItems=!1}setStops(){this.updateScrollStops();let{scrollContainer:e}=this,{scrollLeft:t}=e,{width:i,left:s}=e.getBoundingClientRect();this.width=i;let o=0,n=this.scrollItems.map((e,i)=>{let{left:n,width:a}=e.getBoundingClientRect(),r=Math.round(n+t-s),l=Math.round(r+a);return this.isRtl?-l:(o=l,0===i?0:r)}).concat(o);(n=this.fixScrollMisalign(n)).sort((e,t)=>Math.abs(e)-Math.abs(t)),this.scrollStops=n,this.setFlippers()}validateStops(e=!0){let t=()=>!!this.scrollStops.find(e=>e>0);return!t()&&e&&this.setStops(),t()}fixScrollMisalign(e){if(this.isRtl&&e.some(e=>e>0)){e.sort((e,t)=>t-e);let t=e[0];e=e.map(e=>e-t)}return e}setFlippers(){var e,t;let i=this.scrollContainer.scrollLeft;if(null==(e=this.previousFlipperContainer)||e.classList.toggle("disabled",0===i),this.scrollStops){let e=Math.abs(this.scrollStops[this.scrollStops.length-1]);null==(t=this.nextFlipperContainer)||t.classList.toggle("disabled",this.validateStops(!1)&&Math.abs(i)+this.width>=e)}}scrollInView(e,t=0,i){var s;if("number"!=typeof e&&e&&(e=this.scrollItems.findIndex(t=>t===e||t.contains(e))),void 0!==e){i=null!=i?i:t;let{scrollContainer:o,scrollStops:n,scrollItems:a}=this,{scrollLeft:r}=this.scrollContainer,{width:l}=o.getBoundingClientRect(),h=n[e],{width:d}=a[e].getBoundingClientRect(),c=h+d,u=r+t>h;if(u||r+l-i<c){let e=null!=(s=[...n].sort((e,t)=>u?t-e:e-t).find(e=>u?e+t<h:e+l-(null!=i?i:0)>c))?s:0;this.scrollToPosition(e)}}}keyupHandler(e){switch(e.key){case"ArrowLeft":this.scrollToPrevious();break;case"ArrowRight":this.scrollToNext()}}scrollToPrevious(){this.validateStops();let e=this.scrollContainer.scrollLeft,t=this.scrollStops.findIndex((t,i)=>t>=e&&(this.isRtl||i===this.scrollStops.length-1||this.scrollStops[i+1]>e)),i=Math.abs(this.scrollStops[t+1]),s=this.scrollStops.findIndex(e=>Math.abs(e)+this.width>i);(s>=t||-1===s)&&(s=t>0?t-1:0),this.scrollToPosition(this.scrollStops[s],e)}scrollToNext(){this.validateStops();let e=this.scrollContainer.scrollLeft,t=this.scrollStops.findIndex(t=>Math.abs(t)>=Math.abs(e)),i=this.scrollStops.findIndex(t=>Math.abs(e)+this.width<=Math.abs(t)),s=t;i>t+2?s=i-2:t<this.scrollStops.length-2&&(s=t+1),this.scrollToPosition(this.scrollStops[s],e)}scrollToPosition(e,t=this.scrollContainer.scrollLeft){var i;if(this.scrolling)return;this.scrolling=!0;let s=null!=(i=this.duration)?i:`${Math.abs(e-t)/this.speed}s`;this.content.style.setProperty("transition-duration",s);let o=parseFloat(getComputedStyle(this.content).getPropertyValue("transition-duration")),n=t=>{t&&t.target!==t.currentTarget||(this.content.style.setProperty("transition-duration","0s"),this.content.style.removeProperty("transform"),this.scrollContainer.style.setProperty("scroll-behavior","auto"),this.scrollContainer.scrollLeft=e,this.setFlippers(),this.content.removeEventListener("transitionend",n),this.scrolling=!1)};if(0===o)return void n();this.content.addEventListener("transitionend",n);let a=this.scrollContainer.scrollWidth-this.scrollContainer.clientWidth,r=this.scrollContainer.scrollLeft-Math.min(e,a);this.isRtl&&(r=this.scrollContainer.scrollLeft+Math.min(Math.abs(e),a)),this.content.style.setProperty("transition-property","transform"),this.content.style.setProperty("transition-timing-function",this.easing),this.content.style.setProperty("transform",`translateX(${r}px)`)}resized(){this.resizeTimeout&&(this.resizeTimeout=clearTimeout(this.resizeTimeout)),this.resizeTimeout=setTimeout(()=>{this.width=this.scrollContainer.offsetWidth,this.setFlippers()},this.frameTime)}scrolled(){this.scrollTimeout&&(this.scrollTimeout=clearTimeout(this.scrollTimeout)),this.scrollTimeout=setTimeout(()=>{this.setFlippers()},this.frameTime)}}f([(0,d.attr)({converter:d.nullableNumberConverter})],i3.prototype,"speed",void 0),f([d.attr],i3.prototype,"duration",void 0),f([d.attr],i3.prototype,"easing",void 0),f([(0,d.attr)({attribute:"flippers-hidden-from-at",converter:d.booleanConverter})],i3.prototype,"flippersHiddenFromAT",void 0),f([d.observable],i3.prototype,"scrolling",void 0),f([d.observable],i3.prototype,"scrollItems",void 0),f([(0,d.attr)({attribute:"view"})],i3.prototype,"view",void 0);let i9=(e,t)=>{var i,s;return(0,d.html)`
    <template
        class="horizontal-scroll"
        @keyup="${(e,t)=>e.keyupHandler(t.event)}"
    >
        ${p(e,t)}
        <div class="scroll-area" part="scroll-area">
            <div
                class="scroll-view"
                part="scroll-view"
                @scroll="${e=>e.scrolled()}"
                ${(0,d.ref)("scrollContainer")}
            >
                <div class="content-container" part="content-container" ${(0,d.ref)("content")}>
                    <slot
                        ${(0,d.slotted)({property:"scrollItems",filter:(0,d.elements)()})}
                    ></slot>
                </div>
            </div>
            ${(0,d.when)(e=>"mobile"!==e.view,(0,d.html)`
                    <div
                        class="scroll scroll-prev"
                        part="scroll-prev"
                        ${(0,d.ref)("previousFlipperContainer")}
                    >
                        <div class="scroll-action" part="scroll-action-previous">
                            <slot name="previous-flipper">
                                ${t.previousFlipper instanceof Function?t.previousFlipper(e,t):null!=(i=t.previousFlipper)?i:""}
                            </slot>
                        </div>
                    </div>
                    <div
                        class="scroll scroll-next"
                        part="scroll-next"
                        ${(0,d.ref)("nextFlipperContainer")}
                    >
                        <div class="scroll-action" part="scroll-action-next">
                            <slot name="next-flipper">
                                ${t.nextFlipper instanceof Function?t.nextFlipper(e,t):null!=(s=t.nextFlipper)?s:""}
                            </slot>
                        </div>
                    </div>
                `)}
        </div>
        ${u(e,t)}
    </template>
`};function i6(e,t,i){return e.nodeType!==Node.TEXT_NODE||"string"==typeof e.nodeValue&&!!e.nodeValue.trim().length}let i7=(e,t)=>(0,d.html)`
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
                ${(0,d.slotted)({property:"defaultSlottedNodes",filter:i6})}
            ></slot>
        </label>
        <div class="root" part="root" ${(0,d.ref)("root")}>
            ${p(e,t)}
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
                    ${(0,d.ref)("control")}
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
            ${u(e,t)}
        </div>
    </template>
`;class se extends ep{}class st extends eQ(se){constructor(){super(...arguments),this.proxy=document.createElement("input")}}class si extends st{readOnlyChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.readOnly=this.readOnly,this.validate())}autofocusChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.autofocus=this.autofocus,this.validate())}placeholderChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.placeholder=this.placeholder)}listChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.setAttribute("list",this.list),this.validate())}maxlengthChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.maxLength=this.maxlength,this.validate())}minlengthChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.minLength=this.minlength,this.validate())}patternChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.pattern=this.pattern,this.validate())}sizeChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.size=this.size)}spellcheckChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.spellcheck=this.spellcheck)}connectedCallback(){super.connectedCallback(),this.validate(),this.autofocus&&d.DOM.queueUpdate(()=>{this.focus()})}validate(){super.validate(this.control)}handleTextInput(){this.value=this.control.value}handleClearInput(){this.value="",this.control.focus(),this.handleChange()}handleChange(){this.$emit("change")}}f([(0,d.attr)({attribute:"readonly",mode:"boolean"})],si.prototype,"readOnly",void 0),f([(0,d.attr)({mode:"boolean"})],si.prototype,"autofocus",void 0),f([d.attr],si.prototype,"placeholder",void 0),f([d.attr],si.prototype,"list",void 0),f([(0,d.attr)({converter:d.nullableNumberConverter})],si.prototype,"maxlength",void 0),f([(0,d.attr)({converter:d.nullableNumberConverter})],si.prototype,"minlength",void 0),f([d.attr],si.prototype,"pattern",void 0),f([(0,d.attr)({converter:d.nullableNumberConverter})],si.prototype,"size",void 0),f([(0,d.attr)({mode:"boolean"})],si.prototype,"spellcheck",void 0),f([d.observable],si.prototype,"defaultSlottedNodes",void 0);class ss{}eb(ss,eI),eb(si,c,ss);class so extends iy{}class sn extends eQ(so){constructor(){super(...arguments),this.proxy=document.createElement("select")}}class sa extends sn{constructor(){super(...arguments),this.open=!1,this.forcedPosition=!1,this.listboxId=tC("listbox-"),this.maxHeight=0}openChanged(e,t){if(this.collapsible){if(this.open){this.ariaControls=this.listboxId,this.ariaExpanded="true",this.setPositioning(),this.focusAndScrollOptionIntoView(),this.indexWhenOpened=this.selectedIndex,d.DOM.queueUpdate(()=>this.focus());return}this.ariaControls="",this.ariaExpanded="false"}}get collapsible(){return!(this.multiple||"number"==typeof this.size)}get value(){return d.Observable.track(this,"value"),this._value}set value(e){var t,i,s,o,n,a,r;let l=`${this._value}`;if(null==(t=this._options)?void 0:t.length){let t=this._options.findIndex(t=>t.value===e),l=null!=(s=null==(i=this._options[this.selectedIndex])?void 0:i.value)?s:null,h=null!=(n=null==(o=this._options[t])?void 0:o.value)?n:null;(-1===t||l!==h)&&(e="",this.selectedIndex=t),e=null!=(r=null==(a=this.firstSelectedOption)?void 0:a.value)?r:e}l!==e&&(this._value=e,super.valueChanged(l,e),d.Observable.notify(this,"value"),this.updateDisplayValue())}updateValue(e){var t,i;this.$fastController.isConnected&&(this.value=null!=(i=null==(t=this.firstSelectedOption)?void 0:t.value)?i:""),e&&(this.$emit("input"),this.$emit("change",this,{bubbles:!0,composed:void 0}))}selectedIndexChanged(e,t){super.selectedIndexChanged(e,t),this.updateValue()}positionChanged(e,t){this.positionAttribute=t,this.setPositioning()}setPositioning(){let e=this.getBoundingClientRect(),t=window.innerHeight-e.bottom;this.position=this.forcedPosition?this.positionAttribute:e.top>t?tO.above:tO.below,this.positionAttribute=this.forcedPosition?this.positionAttribute:this.position,this.maxHeight=this.position===tO.above?~~e.top:~~t}get displayValue(){var e,t;return d.Observable.track(this,"displayValue"),null!=(t=null==(e=this.firstSelectedOption)?void 0:e.text)?t:""}disabledChanged(e,t){super.disabledChanged&&super.disabledChanged(e,t),this.ariaDisabled=this.disabled?"true":"false"}formResetCallback(){this.setProxyOptions(),super.setDefaultSelectedOption(),-1===this.selectedIndex&&(this.selectedIndex=0)}clickHandler(e){if(!this.disabled){if(this.open){let t=e.target.closest("option,[role=option]");if(t&&t.disabled)return}return super.clickHandler(e),this.open=this.collapsible&&!this.open,this.open||this.indexWhenOpened===this.selectedIndex||this.updateValue(!0),!0}}focusoutHandler(e){var t;if(super.focusoutHandler(e),!this.open)return!0;let i=e.relatedTarget;this.isSameNode(i)?this.focus():(null==(t=this.options)?void 0:t.includes(i))||(this.open=!1,this.indexWhenOpened!==this.selectedIndex&&this.updateValue(!0))}handleChange(e,t){super.handleChange(e,t),"value"===t&&this.updateValue()}slottedOptionsChanged(e,t){this.options.forEach(e=>{d.Observable.getNotifier(e).unsubscribe(this,"value")}),super.slottedOptionsChanged(e,t),this.options.forEach(e=>{d.Observable.getNotifier(e).subscribe(this,"value")}),this.setProxyOptions(),this.updateValue()}mousedownHandler(e){var t;return e.offsetX>=0&&e.offsetX<=(null==(t=this.listbox)?void 0:t.scrollWidth)?super.mousedownHandler(e):this.collapsible}multipleChanged(e,t){super.multipleChanged(e,t),this.proxy&&(this.proxy.multiple=t)}selectedOptionsChanged(e,t){var i;super.selectedOptionsChanged(e,t),null==(i=this.options)||i.forEach((e,t)=>{var i;let s=null==(i=this.proxy)?void 0:i.options.item(t);s&&(s.selected=e.selected)})}setDefaultSelectedOption(){var e;let t=null!=(e=this.options)?e:Array.from(this.children).filter(tk.slottedOptionFilter),i=null==t?void 0:t.findIndex(e=>e.hasAttribute("selected")||e.selected||e.value===this.value);if(-1!==i){this.selectedIndex=i;return}this.selectedIndex=0}setProxyOptions(){this.proxy instanceof HTMLSelectElement&&this.options&&(this.proxy.options.length=0,this.options.forEach(e=>{let t=e.proxy||(e instanceof HTMLOptionElement?e.cloneNode():null);t&&this.proxy.options.add(t)}))}keydownHandler(e){super.keydownHandler(e);let t=e.key||e.key.charCodeAt(0);switch(t){case" ":e.preventDefault(),this.collapsible&&this.typeAheadExpired&&(this.open=!this.open);break;case"Home":case"End":e.preventDefault();break;case"Enter":e.preventDefault(),this.open=!this.open;break;case"Escape":this.collapsible&&this.open&&(e.preventDefault(),this.open=!1);break;case"Tab":return this.collapsible&&this.open&&(e.preventDefault(),this.open=!1),!0}return this.open||this.indexWhenOpened===this.selectedIndex||(this.updateValue(!0),this.indexWhenOpened=this.selectedIndex),t!==ey.HX&&t!==ey.I5}connectedCallback(){super.connectedCallback(),this.forcedPosition=!!this.positionAttribute,this.addEventListener("contentchange",this.updateDisplayValue)}disconnectedCallback(){this.removeEventListener("contentchange",this.updateDisplayValue),super.disconnectedCallback()}sizeChanged(e,t){super.sizeChanged(e,t),this.proxy&&(this.proxy.size=t)}updateDisplayValue(){this.collapsible&&d.Observable.notify(this,"displayValue")}}f([(0,d.attr)({attribute:"open",mode:"boolean"})],sa.prototype,"open",void 0),f([d.volatile],sa.prototype,"collapsible",null),f([d.observable],sa.prototype,"control",void 0),f([(0,d.attr)({attribute:"position"})],sa.prototype,"positionAttribute",void 0),f([d.observable],sa.prototype,"position",void 0),f([d.observable],sa.prototype,"maxHeight",void 0);class sr{}f([d.observable],sr.prototype,"ariaControls",void 0),eb(sr,tE),eb(sa,c,sr);let sl=(e,t)=>(0,d.html)`
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
        tabindex="${e=>e.disabled?null:"0"}"
        @click="${(e,t)=>e.clickHandler(t.event)}"
        @focusin="${(e,t)=>e.focusinHandler(t.event)}"
        @focusout="${(e,t)=>e.focusoutHandler(t.event)}"
        @keydown="${(e,t)=>e.keydownHandler(t.event)}"
        @mousedown="${(e,t)=>e.mousedownHandler(t.event)}"
    >
        ${(0,d.when)(e=>e.collapsible,(0,d.html)`
                <div
                    class="control"
                    part="control"
                    ?disabled="${e=>e.disabled}"
                    ${(0,d.ref)("control")}
                >
                    ${p(e,t)}
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
                    ${u(e,t)}
                </div>
            `)}
        <div
            class="listbox"
            id="${e=>e.listboxId}"
            part="listbox"
            role="listbox"
            ?disabled="${e=>e.disabled}"
            ?hidden="${e=>!!e.collapsible&&!e.open}"
            ${(0,d.ref)("listbox")}
        >
            <slot
                ${(0,d.slotted)({filter:tk.slottedOptionFilter,flatten:!0,property:"slottedOptions"})}
            ></slot>
        </div>
    </template>
`,sh=(e,t)=>(0,d.html)`
    <template
        class="${e=>"circle"===e.shape?"circle":"rect"}"
        pattern="${e=>e.pattern}"
        ?shimmer="${e=>e.shimmer}"
    >
        ${(0,d.when)(e=>!0===e.shimmer,(0,d.html)`
                <span class="shimmer"></span>
            `)}
        <object type="image/svg+xml" data="${e=>e.pattern}" role="presentation">
            <img class="pattern" src="${e=>e.pattern}" />
        </object>
        <slot></slot>
    </template>
`;class sd extends ep{constructor(){super(...arguments),this.shape="rect"}}f([d.attr],sd.prototype,"fill",void 0),f([d.attr],sd.prototype,"shape",void 0),f([d.attr],sd.prototype,"pattern",void 0),f([(0,d.attr)({mode:"boolean"})],sd.prototype,"shimmer",void 0);let sc=(e,t)=>(0,d.html)`
    <template
        aria-disabled="${e=>e.disabled}"
        class="${e=>e.sliderOrientation||ic.t.horizontal}
            ${e=>e.disabled?"disabled":""}"
    >
        <div ${(0,d.ref)("root")} part="root" class="root" style="${e=>e.positionStyle}">
            <div class="container">
                ${(0,d.when)(e=>!e.hideMark,(0,d.html)`
                        <div class="mark"></div>
                    `)}
                <div class="label">
                    <slot></slot>
                </div>
            </div>
        </div>
    </template>
`;function su(e,t,i,s){let o=(0,eC.AB)(0,1,(e-t)/(i-t));return s===eT.O.rtl&&(o=1-o),o}let sp={min:0,max:0,direction:eT.O.ltr,orientation:ic.t.horizontal,disabled:!1};class sm extends ep{constructor(){super(...arguments),this.hideMark=!1,this.sliderDirection=eT.O.ltr,this.getSliderConfiguration=()=>{if(this.isSliderConfig(this.parentNode)){let{min:e,max:t,direction:i,orientation:s,disabled:o}=this.parentNode;void 0!==o&&(this.disabled=o),this.sliderDirection=i||eT.O.ltr,this.sliderOrientation=s||ic.t.horizontal,this.sliderMaxPosition=t,this.sliderMinPosition=e}else this.sliderDirection=sp.direction||eT.O.ltr,this.sliderOrientation=sp.orientation||ic.t.horizontal,this.sliderMaxPosition=sp.max,this.sliderMinPosition=sp.min},this.positionAsStyle=()=>{let e=this.sliderDirection?this.sliderDirection:eT.O.ltr,t=su(Number(this.position),Number(this.sliderMinPosition),Number(this.sliderMaxPosition)),i=Math.round((1-t)*100),s=Math.round(100*t);return(Number.isNaN(s)&&Number.isNaN(i)&&(i=50,s=50),this.sliderOrientation===ic.t.horizontal)?e===eT.O.rtl?`right: ${s}%; left: ${i}%;`:`left: ${s}%; right: ${i}%;`:`top: ${s}%; bottom: ${i}%;`}}positionChanged(){this.positionStyle=this.positionAsStyle()}sliderOrientationChanged(){}connectedCallback(){super.connectedCallback(),this.getSliderConfiguration(),this.positionStyle=this.positionAsStyle(),this.notifier=d.Observable.getNotifier(this.parentNode),this.notifier.subscribe(this,"orientation"),this.notifier.subscribe(this,"direction"),this.notifier.subscribe(this,"max"),this.notifier.subscribe(this,"min")}disconnectedCallback(){super.disconnectedCallback(),this.notifier.unsubscribe(this,"orientation"),this.notifier.unsubscribe(this,"direction"),this.notifier.unsubscribe(this,"max"),this.notifier.unsubscribe(this,"min")}handleChange(e,t){switch(t){case"direction":this.sliderDirection=e.direction;break;case"orientation":this.sliderOrientation=e.orientation;break;case"max":this.sliderMaxPosition=e.max;break;case"min":this.sliderMinPosition=e.min}this.positionStyle=this.positionAsStyle()}isSliderConfig(e){return void 0!==e.max&&void 0!==e.min}}f([d.observable],sm.prototype,"positionStyle",void 0),f([d.attr],sm.prototype,"position",void 0),f([(0,d.attr)({attribute:"hide-mark",mode:"boolean"})],sm.prototype,"hideMark",void 0),f([(0,d.attr)({attribute:"disabled",mode:"boolean"})],sm.prototype,"disabled",void 0),f([d.observable],sm.prototype,"sliderOrientation",void 0),f([d.observable],sm.prototype,"sliderMinPosition",void 0),f([d.observable],sm.prototype,"sliderMaxPosition",void 0),f([d.observable],sm.prototype,"sliderDirection",void 0);let sv=(e,t)=>(0,d.html)`
    <template
        role="slider"
        class="${e=>e.readOnly?"readonly":""}
        ${e=>e.orientation||ic.t.horizontal}"
        tabindex="${e=>e.disabled?null:0}"
        aria-valuetext="${e=>e.valueTextFormatter(e.value)}"
        aria-valuenow="${e=>e.value}"
        aria-valuemin="${e=>e.min}"
        aria-valuemax="${e=>e.max}"
        aria-disabled="${e=>!!e.disabled||void 0}"
        aria-readonly="${e=>!!e.readOnly||void 0}"
        aria-orientation="${e=>e.orientation}"
        class="${e=>e.orientation}"
    >
        <div part="positioning-region" class="positioning-region">
            <div ${(0,d.ref)("track")} part="track-container" class="track">
                <slot name="track"></slot>
                <div part="track-start" class="track-start" style="${e=>e.position}">
                    <slot name="track-start"></slot>
                </div>
            </div>
            <slot></slot>
            <div
                ${(0,d.ref)("thumb")}
                part="thumb-container"
                class="thumb-container"
                style="${e=>e.position}"
            >
                <slot name="thumb">${t.thumb||""}</slot>
            </div>
        </div>
    </template>
`;class sb extends ep{}class sf extends eQ(sb){constructor(){super(...arguments),this.proxy=document.createElement("input")}}let sg={singleValue:"single-value"};class sy extends sf{constructor(){super(...arguments),this.direction=eT.O.ltr,this.isDragging=!1,this.trackWidth=0,this.trackMinWidth=0,this.trackHeight=0,this.trackLeft=0,this.trackMinHeight=0,this.valueTextFormatter=()=>null,this.min=0,this.max=10,this.step=1,this.orientation=ic.t.horizontal,this.mode=sg.singleValue,this.keypressHandler=e=>{if(!this.readOnly){if("Home"===e.key)e.preventDefault(),this.value=`${this.min}`;else if("End"===e.key)e.preventDefault(),this.value=`${this.max}`;else if(!e.shiftKey)switch(e.key){case ey.bb:case ey.I5:e.preventDefault(),this.increment();break;case ey.kT:case ey.HX:e.preventDefault(),this.decrement()}}},this.setupTrackConstraints=()=>{let e=this.track.getBoundingClientRect();this.trackWidth=this.track.clientWidth,this.trackMinWidth=this.track.clientLeft,this.trackHeight=e.bottom,this.trackMinHeight=e.top,this.trackLeft=this.getBoundingClientRect().left,0===this.trackWidth&&(this.trackWidth=1)},this.setupListeners=(e=!1)=>{let t=`${e?"remove":"add"}EventListener`;this[t]("keydown",this.keypressHandler),this[t]("mousedown",this.handleMouseDown),this.thumb[t]("mousedown",this.handleThumbMouseDown,{passive:!0}),this.thumb[t]("touchstart",this.handleThumbMouseDown,{passive:!0}),e&&(this.handleMouseDown(null),this.handleThumbMouseDown(null))},this.initialValue="",this.handleThumbMouseDown=e=>{if(e){if(this.readOnly||this.disabled||e.defaultPrevented)return;e.target.focus()}let t=`${null!==e?"add":"remove"}EventListener`;window[t]("mouseup",this.handleWindowMouseUp),window[t]("mousemove",this.handleMouseMove,{passive:!0}),window[t]("touchmove",this.handleMouseMove,{passive:!0}),window[t]("touchend",this.handleWindowMouseUp),this.isDragging=null!==e},this.handleMouseMove=e=>{if(this.readOnly||this.disabled||e.defaultPrevented)return;let t=window.TouchEvent&&e instanceof TouchEvent?e.touches[0]:e,i=this.orientation===ic.t.horizontal?t.pageX-document.documentElement.scrollLeft-this.trackLeft:t.pageY-document.documentElement.scrollTop;this.value=`${this.calculateNewValue(i)}`},this.calculateNewValue=e=>{let t=su(e,this.orientation===ic.t.horizontal?this.trackMinWidth:this.trackMinHeight,this.orientation===ic.t.horizontal?this.trackWidth:this.trackHeight,this.direction),i=(this.max-this.min)*t+this.min;return this.convertToConstrainedValue(i)},this.handleWindowMouseUp=e=>{this.stopDragging()},this.stopDragging=()=>{this.isDragging=!1,this.handleMouseDown(null),this.handleThumbMouseDown(null)},this.handleMouseDown=e=>{let t=`${null!==e?"add":"remove"}EventListener`;if((null===e||!this.disabled&&!this.readOnly)&&(window[t]("mouseup",this.handleWindowMouseUp),window.document[t]("mouseleave",this.handleWindowMouseUp),window[t]("mousemove",this.handleMouseMove),e)){e.preventDefault(),this.setupTrackConstraints(),e.target.focus();let t=this.orientation===ic.t.horizontal?e.pageX-document.documentElement.scrollLeft-this.trackLeft:e.pageY-document.documentElement.scrollTop;this.value=`${this.calculateNewValue(t)}`}},this.convertToConstrainedValue=e=>{isNaN(e)&&(e=this.min);let t=e-this.min,i=Math.round(t/this.step),s=t-i*(this.stepMultiplier*this.step)/this.stepMultiplier;return(t=s>=Number(this.step)/2?t-s+Number(this.step):t-s)+this.min}}readOnlyChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.readOnly=this.readOnly)}get valueAsNumber(){return parseFloat(super.value)}set valueAsNumber(e){this.value=e.toString()}valueChanged(e,t){super.valueChanged(e,t),this.$fastController.isConnected&&this.setThumbPositionForOrientation(this.direction),this.$emit("change")}minChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.min=`${this.min}`),this.validate()}maxChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.max=`${this.max}`),this.validate()}stepChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.step=`${this.step}`),this.updateStepMultiplier(),this.validate()}orientationChanged(){this.$fastController.isConnected&&this.setThumbPositionForOrientation(this.direction)}connectedCallback(){super.connectedCallback(),this.proxy.setAttribute("type","range"),this.direction=eR(this),this.updateStepMultiplier(),this.setupTrackConstraints(),this.setupListeners(),this.setupDefaultValue(),this.setThumbPositionForOrientation(this.direction)}disconnectedCallback(){this.setupListeners(!0)}increment(){let e=this.direction!==eT.O.rtl&&this.orientation!==ic.t.vertical?Number(this.value)+Number(this.step):Number(this.value)-Number(this.step),t=this.convertToConstrainedValue(e),i=t<Number(this.max)?`${t}`:`${this.max}`;this.value=i}decrement(){let e=this.direction!==eT.O.rtl&&this.orientation!==ic.t.vertical?Number(this.value)-Number(this.step):Number(this.value)+Number(this.step),t=this.convertToConstrainedValue(e),i=t>Number(this.min)?`${t}`:`${this.min}`;this.value=i}setThumbPositionForOrientation(e){let t=(1-su(Number(this.value),Number(this.min),Number(this.max),e))*100;this.orientation===ic.t.horizontal?this.position=this.isDragging?`right: ${t}%; transition: none;`:`right: ${t}%; transition: all 0.2s ease;`:this.position=this.isDragging?`bottom: ${t}%; transition: none;`:`bottom: ${t}%; transition: all 0.2s ease;`}updateStepMultiplier(){let e=this.step+"",t=this.step%1?e.length-e.indexOf(".")-1:0;this.stepMultiplier=Math.pow(10,t)}get midpoint(){return`${this.convertToConstrainedValue((this.max+this.min)/2)}`}setupDefaultValue(){if("string"==typeof this.value)if(0===this.value.length)this.initialValue=this.midpoint;else{let e=parseFloat(this.value);!Number.isNaN(e)&&(e<this.min||e>this.max)&&(this.value=this.midpoint)}}}f([(0,d.attr)({attribute:"readonly",mode:"boolean"})],sy.prototype,"readOnly",void 0),f([d.observable],sy.prototype,"direction",void 0),f([d.observable],sy.prototype,"isDragging",void 0),f([d.observable],sy.prototype,"position",void 0),f([d.observable],sy.prototype,"trackWidth",void 0),f([d.observable],sy.prototype,"trackMinWidth",void 0),f([d.observable],sy.prototype,"trackHeight",void 0),f([d.observable],sy.prototype,"trackLeft",void 0),f([d.observable],sy.prototype,"trackMinHeight",void 0),f([d.observable],sy.prototype,"valueTextFormatter",void 0),f([(0,d.attr)({converter:d.nullableNumberConverter})],sy.prototype,"min",void 0),f([(0,d.attr)({converter:d.nullableNumberConverter})],sy.prototype,"max",void 0),f([(0,d.attr)({converter:d.nullableNumberConverter})],sy.prototype,"step",void 0),f([d.attr],sy.prototype,"orientation",void 0),f([d.attr],sy.prototype,"mode",void 0);let sC=(e,t)=>(0,d.html)`
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
            <slot ${(0,d.slotted)("defaultSlottedNodes")}></slot>
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
`;class sx extends ep{}class s$ extends eZ(sx){constructor(){super(...arguments),this.proxy=document.createElement("input")}}class sw extends s${constructor(){super(),this.initialValue="on",this.keypressHandler=e=>{if(!this.readOnly)switch(e.key){case"Enter":case" ":this.checked=!this.checked}},this.clickHandler=e=>{this.disabled||this.readOnly||(this.checked=!this.checked)},this.proxy.setAttribute("type","checkbox")}readOnlyChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.readOnly=this.readOnly),this.readOnly?this.classList.add("readonly"):this.classList.remove("readonly")}checkedChanged(e,t){super.checkedChanged(e,t),this.checked?this.classList.add("checked"):this.classList.remove("checked")}}f([(0,d.attr)({attribute:"readonly",mode:"boolean"})],sw.prototype,"readOnly",void 0),f([d.observable],sw.prototype,"defaultSlottedNodes",void 0);let sI=(e,t)=>(0,d.html)`
    <template slot="tabpanel" role="tabpanel">
        <slot></slot>
    </template>
`;class sk extends ep{}let sE=(e,t)=>(0,d.html)`
    <template slot="tab" role="tab" aria-disabled="${e=>e.disabled}">
        <slot></slot>
    </template>
`;class sO extends ep{}f([(0,d.attr)({mode:"boolean"})],sO.prototype,"disabled",void 0);let sT=(e,t)=>(0,d.html)`
    <template class="${e=>e.orientation}">
        ${p(e,t)}
        <div class="tablist" part="tablist" role="tablist">
            <slot class="tab" name="tab" part="tab" ${(0,d.slotted)("tabs")}></slot>

            ${(0,d.when)(e=>e.showActiveIndicator,(0,d.html)`
                    <div
                        ${(0,d.ref)("activeIndicatorRef")}
                        class="activeIndicator"
                        part="activeIndicator"
                    ></div>
                `)}
        </div>
        ${u(e,t)}
        <div class="tabpanel" part="tabpanel">
            <slot name="tabpanel" ${(0,d.slotted)("tabpanels")}></slot>
        </div>
    </template>
`,sR={vertical:"vertical",horizontal:"horizontal"};class sD extends ep{constructor(){super(...arguments),this.orientation=sR.horizontal,this.activeindicator=!0,this.showActiveIndicator=!0,this.prevActiveTabIndex=0,this.activeTabIndex=0,this.ticking=!1,this.change=()=>{this.$emit("change",this.activetab)},this.isDisabledElement=e=>"true"===e.getAttribute("aria-disabled"),this.isHiddenElement=e=>e.hasAttribute("hidden"),this.isFocusableElement=e=>!this.isDisabledElement(e)&&!this.isHiddenElement(e),this.setTabs=()=>{let e="gridColumn",t="gridRow",i=this.isHorizontal()?e:t;this.activeTabIndex=this.getActiveIndex(),this.showActiveIndicator=!1,this.tabs.forEach((s,o)=>{if("tab"===s.slot){let e=this.activeTabIndex===o&&this.isFocusableElement(s);this.activeindicator&&this.isFocusableElement(s)&&(this.showActiveIndicator=!0);let t=this.tabIds[o],i=this.tabpanelIds[o];s.setAttribute("id",t),s.setAttribute("aria-selected",e?"true":"false"),s.setAttribute("aria-controls",i),s.addEventListener("click",this.handleTabClick),s.addEventListener("keydown",this.handleTabKeyDown),s.setAttribute("tabindex",e?"0":"-1"),e&&(this.activetab=s,this.activeid=t)}s.style[e]="",s.style[t]="",s.style[i]=`${o+1}`,this.isHorizontal()?s.classList.remove("vertical"):s.classList.add("vertical")})},this.setTabPanels=()=>{this.tabpanels.forEach((e,t)=>{let i=this.tabIds[t],s=this.tabpanelIds[t];e.setAttribute("id",s),e.setAttribute("aria-labelledby",i),this.activeTabIndex!==t?e.setAttribute("hidden",""):e.removeAttribute("hidden")})},this.handleTabClick=e=>{let t=e.currentTarget;1===t.nodeType&&this.isFocusableElement(t)&&(this.prevActiveTabIndex=this.activeTabIndex,this.activeTabIndex=this.tabs.indexOf(t),this.setComponent())},this.handleTabKeyDown=e=>{if(this.isHorizontal())switch(e.key){case ey.kT:e.preventDefault(),this.adjustBackward(e);break;case ey.bb:e.preventDefault(),this.adjustForward(e)}else switch(e.key){case ey.I5:e.preventDefault(),this.adjustBackward(e);break;case ey.HX:e.preventDefault(),this.adjustForward(e)}switch(e.key){case"Home":e.preventDefault(),this.adjust(-this.activeTabIndex);break;case"End":e.preventDefault(),this.adjust(this.tabs.length-this.activeTabIndex-1)}},this.adjustForward=e=>{let t=this.tabs,i=0;for((i=this.activetab?t.indexOf(this.activetab)+1:1)===t.length&&(i=0);i<t.length&&t.length>1;)if(this.isFocusableElement(t[i])){this.moveToTabByIndex(t,i);break}else if(this.activetab&&i===t.indexOf(this.activetab))break;else i+1>=t.length?i=0:i+=1},this.adjustBackward=e=>{let t=this.tabs,i=0;for(i=(i=this.activetab?t.indexOf(this.activetab)-1:0)<0?t.length-1:i;i>=0&&t.length>1;)if(this.isFocusableElement(t[i])){this.moveToTabByIndex(t,i);break}else i-1<0?i=t.length-1:i-=1},this.moveToTabByIndex=(e,t)=>{let i=e[t];this.activetab=i,this.prevActiveTabIndex=this.activeTabIndex,this.activeTabIndex=t,i.focus(),this.setComponent()}}orientationChanged(){this.$fastController.isConnected&&(this.setTabs(),this.setTabPanels(),this.handleActiveIndicatorPosition())}activeidChanged(e,t){this.$fastController.isConnected&&this.tabs.length<=this.tabpanels.length&&(this.prevActiveTabIndex=this.tabs.findIndex(t=>t.id===e),this.setTabs(),this.setTabPanels(),this.handleActiveIndicatorPosition())}tabsChanged(){this.$fastController.isConnected&&this.tabs.length<=this.tabpanels.length&&(this.tabIds=this.getTabIds(),this.tabpanelIds=this.getTabPanelIds(),this.setTabs(),this.setTabPanels(),this.handleActiveIndicatorPosition())}tabpanelsChanged(){this.$fastController.isConnected&&this.tabpanels.length<=this.tabs.length&&(this.tabIds=this.getTabIds(),this.tabpanelIds=this.getTabPanelIds(),this.setTabs(),this.setTabPanels(),this.handleActiveIndicatorPosition())}getActiveIndex(){return void 0!==this.activeid?-1===this.tabIds.indexOf(this.activeid)?0:this.tabIds.indexOf(this.activeid):0}getTabIds(){return this.tabs.map(e=>{var t;return null!=(t=e.getAttribute("id"))?t:`tab-${tC()}`})}getTabPanelIds(){return this.tabpanels.map(e=>{var t;return null!=(t=e.getAttribute("id"))?t:`panel-${tC()}`})}setComponent(){this.activeTabIndex!==this.prevActiveTabIndex&&(this.activeid=this.tabIds[this.activeTabIndex],this.focusTab(),this.change())}isHorizontal(){return this.orientation===sR.horizontal}handleActiveIndicatorPosition(){this.showActiveIndicator&&this.activeindicator&&this.activeTabIndex!==this.prevActiveTabIndex&&(this.ticking?this.ticking=!1:(this.ticking=!0,this.animateActiveIndicator()))}animateActiveIndicator(){this.ticking=!0;let e=this.isHorizontal()?"gridColumn":"gridRow",t=this.isHorizontal()?"translateX":"translateY",i=this.isHorizontal()?"offsetLeft":"offsetTop",s=this.activeIndicatorRef[i];this.activeIndicatorRef.style[e]=`${this.activeTabIndex+1}`;let o=this.activeIndicatorRef[i];this.activeIndicatorRef.style[e]=`${this.prevActiveTabIndex+1}`,this.activeIndicatorRef.style.transform=`${t}(${o-s}px)`,this.activeIndicatorRef.classList.add("activeIndicatorTransition"),this.activeIndicatorRef.addEventListener("transitionend",()=>{this.ticking=!1,this.activeIndicatorRef.style[e]=`${this.activeTabIndex+1}`,this.activeIndicatorRef.style.transform=`${t}(0px)`,this.activeIndicatorRef.classList.remove("activeIndicatorTransition")})}adjust(e){let t=this.tabs.filter(e=>this.isFocusableElement(e)),i=t.indexOf(this.activetab),s=(0,eC.AB)(0,t.length-1,i+e),o=this.tabs.indexOf(t[s]);o>-1&&this.moveToTabByIndex(this.tabs,o)}focusTab(){this.tabs[this.activeTabIndex].focus()}connectedCallback(){super.connectedCallback(),this.tabIds=this.getTabIds(),this.tabpanelIds=this.getTabPanelIds(),this.activeTabIndex=this.getActiveIndex()}}f([d.attr],sD.prototype,"orientation",void 0),f([d.attr],sD.prototype,"activeid",void 0),f([d.observable],sD.prototype,"tabs",void 0),f([d.observable],sD.prototype,"tabpanels",void 0),f([(0,d.attr)({mode:"boolean"})],sD.prototype,"activeindicator",void 0),f([d.observable],sD.prototype,"activeIndicatorRef",void 0),f([d.observable],sD.prototype,"showActiveIndicator",void 0),eb(sD,c);let sS={none:"none",both:"both",horizontal:"horizontal",vertical:"vertical"},sA=(e,t)=>(0,d.html)`
    <template
        class="
            ${e=>e.readOnly?"readonly":""}
            ${e=>e.resize!==sS.none?`resize-${e.resize}`:""}"
    >
        <label
            part="label"
            for="control"
            class="${e=>e.defaultSlottedNodes&&e.defaultSlottedNodes.length?"label":"label label__hidden"}"
        >
            <slot ${(0,d.slotted)("defaultSlottedNodes")}></slot>
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
            ${(0,d.ref)("control")}
        ></textarea>
    </template>
`;class sF extends ep{}class sL extends eQ(sF){constructor(){super(...arguments),this.proxy=document.createElement("textarea")}}class sM extends sL{constructor(){super(...arguments),this.resize=sS.none,this.cols=20,this.handleTextInput=()=>{this.value=this.control.value}}readOnlyChanged(){this.proxy instanceof HTMLTextAreaElement&&(this.proxy.readOnly=this.readOnly)}autofocusChanged(){this.proxy instanceof HTMLTextAreaElement&&(this.proxy.autofocus=this.autofocus)}listChanged(){this.proxy instanceof HTMLTextAreaElement&&this.proxy.setAttribute("list",this.list)}maxlengthChanged(){this.proxy instanceof HTMLTextAreaElement&&(this.proxy.maxLength=this.maxlength)}minlengthChanged(){this.proxy instanceof HTMLTextAreaElement&&(this.proxy.minLength=this.minlength)}spellcheckChanged(){this.proxy instanceof HTMLTextAreaElement&&(this.proxy.spellcheck=this.spellcheck)}select(){this.control.select(),this.$emit("select")}handleChange(){this.$emit("change")}validate(){super.validate(this.control)}}f([(0,d.attr)({mode:"boolean"})],sM.prototype,"readOnly",void 0),f([d.attr],sM.prototype,"resize",void 0),f([(0,d.attr)({mode:"boolean"})],sM.prototype,"autofocus",void 0),f([(0,d.attr)({attribute:"form"})],sM.prototype,"formId",void 0),f([d.attr],sM.prototype,"list",void 0),f([(0,d.attr)({converter:d.nullableNumberConverter})],sM.prototype,"maxlength",void 0),f([(0,d.attr)({converter:d.nullableNumberConverter})],sM.prototype,"minlength",void 0),f([d.attr],sM.prototype,"name",void 0),f([d.attr],sM.prototype,"placeholder",void 0),f([(0,d.attr)({converter:d.nullableNumberConverter,mode:"fromView"})],sM.prototype,"cols",void 0),f([(0,d.attr)({converter:d.nullableNumberConverter,mode:"fromView"})],sM.prototype,"rows",void 0),f([(0,d.attr)({mode:"boolean"})],sM.prototype,"spellcheck",void 0),f([d.observable],sM.prototype,"defaultSlottedNodes",void 0),eb(sM,iW);let sP=(e,t)=>(0,d.html)`
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
                ${(0,d.slotted)({property:"defaultSlottedNodes",filter:i6})}
            ></slot>
        </label>
        <div class="root" part="root">
            ${p(e,t)}
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
                ${(0,d.ref)("control")}
            />
            ${u(e,t)}
        </div>
    </template>
`,sH=(e,t)=>(0,d.html)`
    <template
        aria-label="${e=>e.ariaLabel}"
        aria-labelledby="${e=>e.ariaLabelledby}"
        aria-orientation="${e=>e.orientation}"
        orientation="${e=>e.orientation}"
        role="toolbar"
        @mousedown="${(e,t)=>e.mouseDownHandler(t.event)}"
        @focusin="${(e,t)=>e.focusinHandler(t.event)}"
        @keydown="${(e,t)=>e.keydownHandler(t.event)}"
        ${(0,d.children)({property:"childItems",attributeFilter:["disabled","hidden"],filter:(0,d.elements)(),subtree:!0})}
    >
        <slot name="label"></slot>
        <div class="positioning-region" part="positioning-region">
            ${p(e,t)}
            <slot
                ${(0,d.slotted)({filter:(0,d.elements)(),property:"slottedItems"})}
            ></slot>
            ${u(e,t)}
        </div>
    </template>
`,sz=Object.freeze({[ey.Is.ArrowUp]:{[ic.t.vertical]:-1},[ey.Is.ArrowDown]:{[ic.t.vertical]:1},[ey.Is.ArrowLeft]:{[ic.t.horizontal]:{[eT.O.ltr]:-1,[eT.O.rtl]:1}},[ey.Is.ArrowRight]:{[ic.t.horizontal]:{[eT.O.ltr]:1,[eT.O.rtl]:-1}}});class sV extends ep{constructor(){super(...arguments),this._activeIndex=0,this.direction=eT.O.ltr,this.orientation=ic.t.horizontal}get activeIndex(){return d.Observable.track(this,"activeIndex"),this._activeIndex}set activeIndex(e){this.$fastController.isConnected&&(this._activeIndex=(0,eC.AB)(0,this.focusableElements.length-1,e),d.Observable.notify(this,"activeIndex"))}slottedItemsChanged(){this.$fastController.isConnected&&this.reduceFocusableElements()}mouseDownHandler(e){var t;let i=null==(t=this.focusableElements)?void 0:t.findIndex(t=>t.contains(e.target));return i>-1&&this.activeIndex!==i&&this.setFocusedElement(i),!0}childItemsChanged(e,t){this.$fastController.isConnected&&this.reduceFocusableElements()}connectedCallback(){super.connectedCallback(),this.direction=eR(this)}focusinHandler(e){let t=e.relatedTarget;!t||this.contains(t)||this.setFocusedElement()}getDirectionalIncrementer(e){var t,i,s,o,n;return null!=(n=null!=(s=null==(i=null==(t=sz[e])?void 0:t[this.orientation])?void 0:i[this.direction])?s:null==(o=sz[e])?void 0:o[this.orientation])?n:0}keydownHandler(e){let t=e.key;if(!(t in ey.Is)||e.defaultPrevented||e.shiftKey)return!0;let i=this.getDirectionalIncrementer(t);if(!i)return!e.target.closest("[role=radiogroup]");let s=this.activeIndex+i;return this.focusableElements[s]&&e.preventDefault(),this.setFocusedElement(s),!0}get allSlottedItems(){return[...this.start.assignedElements(),...this.slottedItems,...this.end.assignedElements()]}reduceFocusableElements(){var e;let t=null==(e=this.focusableElements)?void 0:e[this.activeIndex];this.focusableElements=this.allSlottedItems.reduce(sV.reduceFocusableItems,[]);let i=this.focusableElements.indexOf(t);this.activeIndex=Math.max(0,i),this.setFocusableElements()}setFocusedElement(e=this.activeIndex){var t;this.activeIndex=e,this.setFocusableElements(),null==(t=this.focusableElements[this.activeIndex])||t.focus()}static reduceFocusableItems(e,t){var i,s,o,n;let a="radio"===t.getAttribute("role"),r=null==(s=null==(i=t.$fastController)?void 0:i.definition.shadowOptions)?void 0:s.delegatesFocus,l=Array.from(null!=(n=null==(o=t.shadowRoot)?void 0:o.querySelectorAll("*"))?n:[]).some(e=>(0,it.tp)(e));return!t.hasAttribute("disabled")&&!t.hasAttribute("hidden")&&((0,it.tp)(t)||a||r||l)?(e.push(t),e):t.childElementCount?e.concat(Array.from(t.children).reduce(sV.reduceFocusableItems,[])):e}setFocusableElements(){this.$fastController.isConnected&&this.focusableElements.length>0&&this.focusableElements.forEach((e,t)=>{e.tabIndex=this.activeIndex===t?0:-1})}}f([d.observable],sV.prototype,"direction",void 0),f([d.attr],sV.prototype,"orientation",void 0),f([d.observable],sV.prototype,"slottedItems",void 0),f([d.observable],sV.prototype,"slottedLabel",void 0),f([d.observable],sV.prototype,"childItems",void 0);class sN{}f([(0,d.attr)({attribute:"aria-labelledby"})],sN.prototype,"ariaLabelledby",void 0),f([(0,d.attr)({attribute:"aria-label"})],sN.prototype,"ariaLabel",void 0),eb(sN,eI),eb(sV,c,sN);let sB=(e,t)=>(0,d.html)`
        ${(0,d.when)(e=>e.tooltipVisible,(0,d.html)`
            <${e.tagFor(eD)}
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
                ${(0,d.ref)("region")}
            >
                <div class="tooltip" part="tooltip" role="tooltip">
                    <slot></slot>
                </div>
            </${e.tagFor(eD)}>
        `)}
    `,sq={top:"top",right:"right",bottom:"bottom",left:"left",start:"start",end:"end",topLeft:"top-left",topRight:"top-right",bottomLeft:"bottom-left",bottomRight:"bottom-right",topStart:"top-start",topEnd:"top-end",bottomStart:"bottom-start",bottomEnd:"bottom-end"};class sU extends ep{constructor(){super(...arguments),this.anchor="",this.delay=300,this.autoUpdateMode="anchor",this.anchorElement=null,this.viewportElement=null,this.verticalPositioningMode="dynamic",this.horizontalPositioningMode="dynamic",this.horizontalInset="false",this.verticalInset="false",this.horizontalScaling="content",this.verticalScaling="content",this.verticalDefaultPosition=void 0,this.horizontalDefaultPosition=void 0,this.tooltipVisible=!1,this.currentDirection=eT.O.ltr,this.showDelayTimer=null,this.hideDelayTimer=null,this.isAnchorHoveredFocused=!1,this.isRegionHovered=!1,this.handlePositionChange=e=>{this.classList.toggle("top","start"===this.region.verticalPosition),this.classList.toggle("bottom","end"===this.region.verticalPosition),this.classList.toggle("inset-top","insetStart"===this.region.verticalPosition),this.classList.toggle("inset-bottom","insetEnd"===this.region.verticalPosition),this.classList.toggle("center-vertical","center"===this.region.verticalPosition),this.classList.toggle("left","start"===this.region.horizontalPosition),this.classList.toggle("right","end"===this.region.horizontalPosition),this.classList.toggle("inset-left","insetStart"===this.region.horizontalPosition),this.classList.toggle("inset-right","insetEnd"===this.region.horizontalPosition),this.classList.toggle("center-horizontal","center"===this.region.horizontalPosition)},this.handleRegionMouseOver=e=>{this.isRegionHovered=!0},this.handleRegionMouseOut=e=>{this.isRegionHovered=!1,this.startHideDelayTimer()},this.handleAnchorMouseOver=e=>{if(this.tooltipVisible){this.isAnchorHoveredFocused=!0;return}this.startShowDelayTimer()},this.handleAnchorMouseOut=e=>{this.isAnchorHoveredFocused=!1,this.clearShowDelayTimer(),this.startHideDelayTimer()},this.handleAnchorFocusIn=e=>{this.startShowDelayTimer()},this.handleAnchorFocusOut=e=>{this.isAnchorHoveredFocused=!1,this.clearShowDelayTimer(),this.startHideDelayTimer()},this.startHideDelayTimer=()=>{this.clearHideDelayTimer(),this.tooltipVisible&&(this.hideDelayTimer=window.setTimeout(()=>{this.updateTooltipVisibility()},60))},this.clearHideDelayTimer=()=>{null!==this.hideDelayTimer&&(clearTimeout(this.hideDelayTimer),this.hideDelayTimer=null)},this.startShowDelayTimer=()=>{if(!this.isAnchorHoveredFocused){if(this.delay>1){null===this.showDelayTimer&&(this.showDelayTimer=window.setTimeout(()=>{this.startHover()},this.delay));return}this.startHover()}},this.startHover=()=>{this.isAnchorHoveredFocused=!0,this.updateTooltipVisibility()},this.clearShowDelayTimer=()=>{null!==this.showDelayTimer&&(clearTimeout(this.showDelayTimer),this.showDelayTimer=null)},this.getAnchor=()=>{let e=this.getRootNode();return e instanceof ShadowRoot?e.getElementById(this.anchor):document.getElementById(this.anchor)},this.handleDocumentKeydown=e=>{!e.defaultPrevented&&this.tooltipVisible&&"Escape"===e.key&&(this.isAnchorHoveredFocused=!1,this.updateTooltipVisibility(),this.$emit("dismiss"))},this.updateTooltipVisibility=()=>{if(!1===this.visible)this.hideTooltip();else{if(!0===this.visible||this.isAnchorHoveredFocused||this.isRegionHovered)return void this.showTooltip();this.hideTooltip()}},this.showTooltip=()=>{this.tooltipVisible||(this.currentDirection=eR(this),this.tooltipVisible=!0,document.addEventListener("keydown",this.handleDocumentKeydown),d.DOM.queueUpdate(this.setRegionProps))},this.hideTooltip=()=>{this.tooltipVisible&&(this.clearHideDelayTimer(),null!==this.region&&void 0!==this.region&&(this.region.removeEventListener("positionchange",this.handlePositionChange),this.region.viewportElement=null,this.region.anchorElement=null,this.region.removeEventListener("mouseover",this.handleRegionMouseOver),this.region.removeEventListener("mouseout",this.handleRegionMouseOut)),document.removeEventListener("keydown",this.handleDocumentKeydown),this.tooltipVisible=!1)},this.setRegionProps=()=>{this.tooltipVisible&&(this.region.viewportElement=this.viewportElement,this.region.anchorElement=this.anchorElement,this.region.addEventListener("positionchange",this.handlePositionChange),this.region.addEventListener("mouseover",this.handleRegionMouseOver,{passive:!0}),this.region.addEventListener("mouseout",this.handleRegionMouseOut,{passive:!0}))}}visibleChanged(){this.$fastController.isConnected&&(this.updateTooltipVisibility(),this.updateLayout())}anchorChanged(){this.$fastController.isConnected&&(this.anchorElement=this.getAnchor())}positionChanged(){this.$fastController.isConnected&&this.updateLayout()}anchorElementChanged(e){if(this.$fastController.isConnected){if(null!=e&&(e.removeEventListener("mouseover",this.handleAnchorMouseOver),e.removeEventListener("mouseout",this.handleAnchorMouseOut),e.removeEventListener("focusin",this.handleAnchorFocusIn),e.removeEventListener("focusout",this.handleAnchorFocusOut)),null!==this.anchorElement&&void 0!==this.anchorElement){this.anchorElement.addEventListener("mouseover",this.handleAnchorMouseOver,{passive:!0}),this.anchorElement.addEventListener("mouseout",this.handleAnchorMouseOut,{passive:!0}),this.anchorElement.addEventListener("focusin",this.handleAnchorFocusIn,{passive:!0}),this.anchorElement.addEventListener("focusout",this.handleAnchorFocusOut,{passive:!0});let e=this.anchorElement.id;null!==this.anchorElement.parentElement&&this.anchorElement.parentElement.querySelectorAll(":hover").forEach(t=>{t.id===e&&this.startShowDelayTimer()})}null!==this.region&&void 0!==this.region&&this.tooltipVisible&&(this.region.anchorElement=this.anchorElement),this.updateLayout()}}viewportElementChanged(){null!==this.region&&void 0!==this.region&&(this.region.viewportElement=this.viewportElement),this.updateLayout()}connectedCallback(){super.connectedCallback(),this.anchorElement=this.getAnchor(),this.updateTooltipVisibility()}disconnectedCallback(){this.hideTooltip(),this.clearShowDelayTimer(),this.clearHideDelayTimer(),super.disconnectedCallback()}updateLayout(){switch(this.verticalPositioningMode="locktodefault",this.horizontalPositioningMode="locktodefault",this.position){case sq.top:case sq.bottom:this.verticalDefaultPosition=this.position,this.horizontalDefaultPosition="center";break;case sq.right:case sq.left:case sq.start:case sq.end:this.verticalDefaultPosition="center",this.horizontalDefaultPosition=this.position;break;case sq.topLeft:this.verticalDefaultPosition="top",this.horizontalDefaultPosition="left";break;case sq.topRight:this.verticalDefaultPosition="top",this.horizontalDefaultPosition="right";break;case sq.bottomLeft:this.verticalDefaultPosition="bottom",this.horizontalDefaultPosition="left";break;case sq.bottomRight:this.verticalDefaultPosition="bottom",this.horizontalDefaultPosition="right";break;case sq.topStart:this.verticalDefaultPosition="top",this.horizontalDefaultPosition="start";break;case sq.topEnd:this.verticalDefaultPosition="top",this.horizontalDefaultPosition="end";break;case sq.bottomStart:this.verticalDefaultPosition="bottom",this.horizontalDefaultPosition="start";break;case sq.bottomEnd:this.verticalDefaultPosition="bottom",this.horizontalDefaultPosition="end";break;default:this.verticalPositioningMode="dynamic",this.horizontalPositioningMode="dynamic",this.verticalDefaultPosition=void 0,this.horizontalDefaultPosition="center"}}}f([(0,d.attr)({mode:"boolean"})],sU.prototype,"visible",void 0),f([d.attr],sU.prototype,"anchor",void 0),f([d.attr],sU.prototype,"delay",void 0),f([d.attr],sU.prototype,"position",void 0),f([(0,d.attr)({attribute:"auto-update-mode"})],sU.prototype,"autoUpdateMode",void 0),f([(0,d.attr)({attribute:"horizontal-viewport-lock"})],sU.prototype,"horizontalViewportLock",void 0),f([(0,d.attr)({attribute:"vertical-viewport-lock"})],sU.prototype,"verticalViewportLock",void 0),f([d.observable],sU.prototype,"anchorElement",void 0),f([d.observable],sU.prototype,"viewportElement",void 0),f([d.observable],sU.prototype,"verticalPositioningMode",void 0),f([d.observable],sU.prototype,"horizontalPositioningMode",void 0),f([d.observable],sU.prototype,"horizontalInset",void 0),f([d.observable],sU.prototype,"verticalInset",void 0),f([d.observable],sU.prototype,"horizontalScaling",void 0),f([d.observable],sU.prototype,"verticalScaling",void 0),f([d.observable],sU.prototype,"verticalDefaultPosition",void 0),f([d.observable],sU.prototype,"horizontalDefaultPosition",void 0),f([d.observable],sU.prototype,"tooltipVisible",void 0),f([d.observable],sU.prototype,"currentDirection",void 0);let sj=(e,t)=>(0,d.html)`
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
        ${(0,d.children)({property:"childItems",filter:(0,d.elements)()})}
    >
        <div class="positioning-region" part="positioning-region">
            <div class="content-region" part="content-region">
                ${(0,d.when)(e=>e.childItems&&e.childItemLength()>0,(0,d.html)`
                        <div
                            aria-hidden="true"
                            class="expand-collapse-button"
                            part="expand-collapse-button"
                            @click="${(e,t)=>e.handleExpandCollapseButtonClick(t.event)}"
                            ${(0,d.ref)("expandCollapseButton")}
                        >
                            <slot name="expand-collapse-glyph">
                                ${t.expandCollapseGlyph||""}
                            </slot>
                        </div>
                    `)}
                ${p(e,t)}
                <slot></slot>
                ${u(e,t)}
            </div>
        </div>
        ${(0,d.when)(e=>e.childItems&&e.childItemLength()>0&&(e.expanded||e.renderCollapsedChildren),(0,d.html)`
                <div role="group" class="items" part="items">
                    <slot name="item" ${(0,d.slotted)("items")}></slot>
                </div>
            `)}
    </template>
`;function s_(e){return tx(e)&&"treeitem"===e.getAttribute("role")}class sK extends ep{constructor(){super(...arguments),this.expanded=!1,this.focusable=!1,this.isNestedItem=()=>s_(this.parentElement),this.handleExpandCollapseButtonClick=e=>{this.disabled||e.defaultPrevented||(this.expanded=!this.expanded)},this.handleFocus=e=>{this.setAttribute("tabindex","0")},this.handleBlur=e=>{this.setAttribute("tabindex","-1")}}expandedChanged(){this.$fastController.isConnected&&this.$emit("expanded-change",this)}selectedChanged(){this.$fastController.isConnected&&this.$emit("selected-change",this)}itemsChanged(e,t){this.$fastController.isConnected&&this.items.forEach(e=>{s_(e)&&(e.nested=!0)})}static focusItem(e){e.focusable=!0,e.focus()}childItemLength(){let e=this.childItems.filter(e=>s_(e));return e?e.length:0}}f([(0,d.attr)({mode:"boolean"})],sK.prototype,"expanded",void 0),f([(0,d.attr)({mode:"boolean"})],sK.prototype,"selected",void 0),f([(0,d.attr)({mode:"boolean"})],sK.prototype,"disabled",void 0),f([d.observable],sK.prototype,"focusable",void 0),f([d.observable],sK.prototype,"childItems",void 0),f([d.observable],sK.prototype,"items",void 0),f([d.observable],sK.prototype,"nested",void 0),f([d.observable],sK.prototype,"renderCollapsedChildren",void 0),eb(sK,c);let sW=(e,t)=>(0,d.html)`
    <template
        role="tree"
        ${(0,d.ref)("treeView")}
        @keydown="${(e,t)=>e.handleKeyDown(t.event)}"
        @focusin="${(e,t)=>e.handleFocus(t.event)}"
        @focusout="${(e,t)=>e.handleBlur(t.event)}"
        @click="${(e,t)=>e.handleClick(t.event)}"
        @selected-change="${(e,t)=>e.handleSelectedChange(t.event)}"
    >
        <slot ${(0,d.slotted)("slottedTreeItems")}></slot>
    </template>
`;class sX extends ep{constructor(){super(...arguments),this.currentFocused=null,this.handleFocus=e=>{if(!(this.slottedTreeItems.length<1)){if(e.target===this){null===this.currentFocused&&(this.currentFocused=this.getValidFocusableItem()),null!==this.currentFocused&&sK.focusItem(this.currentFocused);return}this.contains(e.target)&&(this.setAttribute("tabindex","-1"),this.currentFocused=e.target)}},this.handleBlur=e=>{e.target instanceof HTMLElement&&(null===e.relatedTarget||!this.contains(e.relatedTarget))&&this.setAttribute("tabindex","0")},this.handleKeyDown=e=>{if(e.defaultPrevented)return;if(this.slottedTreeItems.length<1)return!0;let t=this.getVisibleNodes();switch(e.key){case"Home":t.length&&sK.focusItem(t[0]);return;case"End":t.length&&sK.focusItem(t[t.length-1]);return;case ey.kT:if(e.target&&this.isFocusableElement(e.target)){let t=e.target;t instanceof sK&&t.childItemLength()>0&&t.expanded?t.expanded=!1:t instanceof sK&&t.parentElement instanceof sK&&sK.focusItem(t.parentElement)}return!1;case ey.bb:if(e.target&&this.isFocusableElement(e.target)){let t=e.target;t instanceof sK&&t.childItemLength()>0&&!t.expanded?t.expanded=!0:t instanceof sK&&t.childItemLength()>0&&this.focusNextNode(1,e.target)}return;case ey.HX:e.target&&this.isFocusableElement(e.target)&&this.focusNextNode(1,e.target);return;case ey.I5:e.target&&this.isFocusableElement(e.target)&&this.focusNextNode(-1,e.target);return;case"Enter":this.handleClick(e);return}return!0},this.handleSelectedChange=e=>{if(e.defaultPrevented)return;if(!(e.target instanceof Element)||!s_(e.target))return!0;let t=e.target;t.selected?(this.currentSelected&&this.currentSelected!==t&&(this.currentSelected.selected=!1),this.currentSelected=t):t.selected||this.currentSelected!==t||(this.currentSelected=null)},this.setItems=()=>{let e=this.treeView.querySelector("[aria-selected='true']");this.currentSelected=e,null!==this.currentFocused&&this.contains(this.currentFocused)||(this.currentFocused=this.getValidFocusableItem()),this.nested=this.checkForNestedItems(),this.getVisibleNodes().forEach(e=>{s_(e)&&(e.nested=this.nested)})},this.isFocusableElement=e=>s_(e),this.isSelectedElement=e=>e.selected}slottedTreeItemsChanged(){this.$fastController.isConnected&&this.setItems()}connectedCallback(){super.connectedCallback(),this.setAttribute("tabindex","0"),d.DOM.queueUpdate(()=>{this.setItems()})}handleClick(e){if(e.defaultPrevented)return;if(!(e.target instanceof Element)||!s_(e.target))return!0;let t=e.target;t.disabled||(t.selected=!t.selected)}focusNextNode(e,t){let i=this.getVisibleNodes();if(!i)return;let s=i[i.indexOf(t)+e];tx(s)&&sK.focusItem(s)}getValidFocusableItem(){let e=this.getVisibleNodes(),t=e.findIndex(this.isSelectedElement);return(-1===t&&(t=e.findIndex(this.isFocusableElement)),-1!==t)?e[t]:null}checkForNestedItems(){return this.slottedTreeItems.some(e=>s_(e)&&e.querySelector("[role='treeitem']"))}getVisibleNodes(){return function(e,t){if(e&&t&&tx(e))return Array.from(e.querySelectorAll(t)).filter(e=>null!==e.offsetParent)}(this,"[role='treeitem']")||[]}}f([(0,d.attr)({attribute:"render-collapsed-nodes"})],sX.prototype,"renderCollapsedNodes",void 0),f([d.observable],sX.prototype,"currentSelected",void 0),f([d.observable],sX.prototype,"slottedTreeItems",void 0);class sG{constructor(e){this.listenerCache=new WeakMap,this.query=e}bind(e){let{query:t}=this,i=this.constructListener(e);i.bind(t)(),t.addListener(i),this.listenerCache.set(e,i)}unbind(e){let t=this.listenerCache.get(e);t&&(this.query.removeListener(t),this.listenerCache.delete(e))}}class sY extends sG{constructor(e,t){super(e),this.styles=t}static with(e){return t=>new sY(e,t)}constructListener(e){let t=!1,i=this.styles;return function(){let{matches:s}=this;s&&!t?(e.$fastController.addStyles(i),t=s):!s&&t&&(e.$fastController.removeStyles(i),t=s)}}unbind(e){super.unbind(e),e.$fastController.removeStyles(this.styles)}}let sQ=sY.with(window.matchMedia("(forced-colors)")),sZ=sY.with(window.matchMedia("(prefers-color-scheme: dark)")),sJ=sY.with(window.matchMedia("(prefers-color-scheme: light)"));class s0{constructor(e,t,i){this.propertyName=e,this.value=t,this.styles=i}bind(e){d.Observable.getNotifier(e).subscribe(this,this.propertyName),this.handleChange(e,this.propertyName)}unbind(e){d.Observable.getNotifier(e).unsubscribe(this,this.propertyName),e.$fastController.removeStyles(this.styles)}handleChange(e,t){e[t]===this.value?e.$fastController.addStyles(this.styles):e.$fastController.removeStyles(this.styles)}}let s1="not-allowed",s5=":host([hidden]){display:none}";function s4(e){return`${s5}:host{display:${e}}`}let s8=!function(){let e;if("boolean"==typeof s)return s;if(!("u">typeof window&&window.document&&window.document.createElement))return s=!1;let t=document.createElement("style"),i=(e=document.querySelector('meta[property="csp-nonce"]'))?e.getAttribute("content"):null;null!==i&&t.setAttribute("nonce",i),document.head.appendChild(t);try{t.sheet.insertRule("foo:focus-visible {color:inherit}",0),s=!0}catch(e){s=!1}finally{document.head.removeChild(t)}return s}()?"focus":"focus-visible"}}]);