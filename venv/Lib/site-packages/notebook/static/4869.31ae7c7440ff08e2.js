"use strict";(self.rspackChunk_JUPYTERLAB_CORE_OUTPUT=self.rspackChunk_JUPYTERLAB_CORE_OUTPUT||[]).push([[4869],{7716(e,t,i){let o,s,r,a,n,l;function h(e,t,i){return isNaN(e)||e<=t?t:e>=i?i:e}function d(e,t,i){return isNaN(e)||e<=t?0:e>=i?1:e/(i-t)}function c(e,t,i){return isNaN(e)?t:t+e*(i-t)}function u(e){return Math.PI/180*e}function p(e,t,i){return isNaN(e)||e<=0?t:e>=1?i:t+e*(i-t)}function g(e,t,i){if(e<=0)return t%360;if(e>=1)return i%360;let o=(t-i+360)%360;return o<=(i-t+360)%360?(t-o*e+360)%360:(t+o*e+360)%360}function m(e,t){let i=Math.pow(10,t);return Math.round(e*i)/i}i.r(t),i.d(t,{Accordion:()=>rk,AccordionItem:()=>r$,Anchor:()=>rN,AnchoredRegion:()=>rK,Avatar:()=>r1,Badge:()=>r4,Breadcrumb:()=>ae,BreadcrumbItem:()=>ao,Button:()=>af,Card:()=>aw,Checkbox:()=>aO,Combobox:()=>aZ,DataGrid:()=>n$,DataGridCell:()=>nf,DataGridRow:()=>ny,DateField:()=>nE,DelegatesARIAToolbar:()=>h8,DesignSystemProvider:()=>nN,Dialog:()=>lc,DirectionalStyleSheetBehavior:()=>rZ,Disclosure:()=>lm,Divider:()=>l$,FoundationToolbar:()=>h9,Listbox:()=>lk,Menu:()=>lE,MenuItem:()=>lA,NumberField:()=>lN,Option:()=>lq,PaletteRGB:()=>ti,Picker:()=>dB,PickerList:()=>dW,PickerListItem:()=>dX,PickerMenu:()=>dq,PickerMenuOption:()=>d_,Progress:()=>lW,ProgressRing:()=>lY,Radio:()=>l5,RadioGroup:()=>l9,Search:()=>ha,Select:()=>hl,Skeleton:()=>hu,Slider:()=>h$,SliderLabel:()=>hS,StandardLuminance:()=>G,SwatchRGB:()=>q,Switch:()=>hA,Tab:()=>hj,TabPanel:()=>hz,Tabs:()=>hW,TextArea:()=>h0,TextField:()=>h5,Toolbar:()=>h7,Tooltip:()=>ds,TreeItem:()=>dg,TreeView:()=>dv,accentColor:()=>ic,accentFillActive:()=>iD,accentFillActiveDelta:()=>tU,accentFillFocus:()=>iE,accentFillFocusDelta:()=>t_,accentFillHover:()=>iO,accentFillHoverDelta:()=>tq,accentFillRecipe:()=>iF,accentFillRest:()=>iS,accentFillRestDelta:()=>tj,accentForegroundActive:()=>iG,accentForegroundActiveDelta:()=>tK,accentForegroundFocus:()=>iW,accentForegroundFocusDelta:()=>tX,accentForegroundHover:()=>i_,accentForegroundHoverDelta:()=>tW,accentForegroundRecipe:()=>iq,accentForegroundRest:()=>iU,accentForegroundRestDelta:()=>tG,accentPalette:()=>iu,accordionItemStyles:()=>rx,accordionStyles:()=>rp,addJupyterLabThemeChangeListener:()=>oW,allComponents:()=>dQ,anchorStyles:()=>rM,anchoredRegionStyles:()=>rW,applyJupyterTheme:()=>oY,avatarStyles:()=>r0,badgeStyles:()=>r3,baseHeightMultiplier:()=>tp,baseHorizontalSpacingMultiplier:()=>tg,baseLayerLuminance:()=>tm,bodyFont:()=>tu,breadcrumbItemStyles:()=>ai,breadcrumbStyles:()=>r7,buttonStyles:()=>ab,cardStyles:()=>a$,checkboxStyles:()=>aF,checkboxTemplate:()=>aS,comboboxStyles:()=>aQ,controlCornerRadius:()=>tb,dataGridCellStyles:()=>nb,dataGridRowStyles:()=>nm,dataGridStyles:()=>ng,dateFieldStyles:()=>nL,dateFieldTemplate:()=>nV,density:()=>tf,designSystemProviderStyles:()=>nj,designSystemProviderTemplate:()=>nB,designUnit:()=>tv,dialogStyles:()=>ld,direction:()=>tx,disabledOpacity:()=>t$,disclosureStyles:()=>lg,dividerStyles:()=>lx,elementScale:()=>ty,errorColor:()=>ow,errorFillActive:()=>oF,errorFillFocus:()=>oS,errorFillHover:()=>oT,errorFillRecipe:()=>oC,errorFillRest:()=>oI,errorForegroundActive:()=>oq,errorForegroundFocus:()=>oU,errorForegroundHover:()=>oj,errorForegroundRecipe:()=>oN,errorForegroundRest:()=>oB,errorPalette:()=>ok,fillColor:()=>iT,focusStrokeInner:()=>oh,focusStrokeInnerRecipe:()=>ol,focusStrokeOuter:()=>on,focusStrokeOuterRecipe:()=>oa,focusStrokeWidth:()=>tk,foregroundOnAccentActive:()=>iP,foregroundOnAccentActiveLarge:()=>iB,foregroundOnAccentFocus:()=>iH,foregroundOnAccentFocusLarge:()=>ij,foregroundOnAccentHover:()=>iV,foregroundOnAccentHoverLarge:()=>iN,foregroundOnAccentLargeRecipe:()=>iz,foregroundOnAccentRecipe:()=>iL,foregroundOnAccentRest:()=>iA,foregroundOnAccentRestLarge:()=>iM,foregroundOnErrorActive:()=>oL,foregroundOnErrorActiveLarge:()=>oz,foregroundOnErrorFocus:()=>oA,foregroundOnErrorFocusLarge:()=>oM,foregroundOnErrorHover:()=>oR,foregroundOnErrorHoverLarge:()=>oH,foregroundOnErrorLargeRecipe:()=>oV,foregroundOnErrorRecipe:()=>oD,foregroundOnErrorRest:()=>oE,foregroundOnErrorRestLarge:()=>oP,heightNumberAsToken:()=>o$,horizontalSliderLabelStyles:()=>hI,imgTemplate:()=>r2,isDark:()=>K,jpAccordion:()=>rC,jpAccordionItem:()=>rw,jpAnchor:()=>rB,jpAnchoredRegion:()=>rX,jpAvatar:()=>r5,jpBadge:()=>r6,jpBreadcrumb:()=>at,jpBreadcrumbItem:()=>as,jpButton:()=>av,jpCard:()=>ak,jpCheckbox:()=>aD,jpCombobox:()=>aJ,jpDataGrid:()=>nw,jpDataGridCell:()=>nv,jpDataGridRow:()=>nx,jpDateField:()=>nP,jpDesignSystemProvider:()=>nq,jpDialog:()=>lu,jpDisclosure:()=>lb,jpDivider:()=>lw,jpListbox:()=>lC,jpMenu:()=>lR,jpMenuItem:()=>lV,jpNumberField:()=>lB,jpOption:()=>lU,jpPicker:()=>dj,jpPickerList:()=>dK,jpPickerListItem:()=>dY,jpPickerMenu:()=>dU,jpPickerMenuOption:()=>dG,jpProgress:()=>lK,jpProgressRing:()=>lQ,jpRadio:()=>l3,jpRadioGroup:()=>l8,jpSearch:()=>hn,jpSelect:()=>hh,jpSkeleton:()=>hp,jpSlider:()=>hw,jpSliderLabel:()=>hO,jpSwitch:()=>hV,jpTab:()=>hq,jpTabPanel:()=>hM,jpTabs:()=>hK,jpTextArea:()=>h1,jpTextField:()=>h3,jpToolbar:()=>de,jpTooltip:()=>dr,jpTreeItem:()=>dm,jpTreeView:()=>dy,listboxStyles:()=>aX,menuItemStyles:()=>lL,menuStyles:()=>lD,neutralColor:()=>ih,neutralFillActive:()=>iQ,neutralFillActiveDelta:()=>tZ,neutralFillFocus:()=>iZ,neutralFillFocusDelta:()=>tJ,neutralFillHover:()=>iY,neutralFillHoverDelta:()=>tQ,neutralFillInputActive:()=>i2,neutralFillInputActiveDelta:()=>t2,neutralFillInputFocus:()=>i5,neutralFillInputFocusDelta:()=>t5,neutralFillInputHover:()=>i1,neutralFillInputHoverDelta:()=>t1,neutralFillInputRecipe:()=>iJ,neutralFillInputRest:()=>i0,neutralFillInputRestDelta:()=>t0,neutralFillLayerRecipe:()=>os,neutralFillLayerRest:()=>or,neutralFillLayerRestDelta:()=>ii,neutralFillRecipe:()=>iK,neutralFillRest:()=>iX,neutralFillRestDelta:()=>tY,neutralFillStealthActive:()=>i9,neutralFillStealthActiveDelta:()=>t6,neutralFillStealthFocus:()=>i8,neutralFillStealthFocusDelta:()=>t9,neutralFillStealthHover:()=>i6,neutralFillStealthHoverDelta:()=>t4,neutralFillStealthRecipe:()=>i3,neutralFillStealthRest:()=>i4,neutralFillStealthRestDelta:()=>t3,neutralFillStrongActive:()=>oi,neutralFillStrongActiveDelta:()=>ie,neutralFillStrongFocus:()=>oo,neutralFillStrongFocusDelta:()=>it,neutralFillStrongHover:()=>ot,neutralFillStrongHoverDelta:()=>t7,neutralFillStrongRecipe:()=>i7,neutralFillStrongRest:()=>oe,neutralFillStrongRestDelta:()=>t8,neutralForegroundHint:()=>oc,neutralForegroundHintRecipe:()=>od,neutralForegroundRecipe:()=>ou,neutralForegroundRest:()=>op,neutralLayer1:()=>iy,neutralLayer1Recipe:()=>iv,neutralLayer2:()=>i$,neutralLayer2Recipe:()=>ix,neutralLayer3:()=>ik,neutralLayer3Recipe:()=>iw,neutralLayer4:()=>iI,neutralLayer4Recipe:()=>iC,neutralLayerCardContainer:()=>ig,neutralLayerCardContainerRecipe:()=>ip,neutralLayerFloating:()=>ib,neutralLayerFloatingRecipe:()=>im,neutralPalette:()=>id,neutralStrokeActive:()=>of,neutralStrokeActiveDelta:()=>ir,neutralStrokeDividerRecipe:()=>oy,neutralStrokeDividerRest:()=>ox,neutralStrokeDividerRestDelta:()=>il,neutralStrokeFocus:()=>ov,neutralStrokeFocusDelta:()=>ia,neutralStrokeHover:()=>ob,neutralStrokeHoverDelta:()=>is,neutralStrokeRecipe:()=>og,neutralStrokeRest:()=>om,neutralStrokeRestDelta:()=>io,numberFieldStyles:()=>lM,optionStyles:()=>lj,pickerListItemStyles:()=>dN,pickerMenuOptionStyles:()=>dM,pickerMenuStyles:()=>dz,pickerStyles:()=>dH,progressRingStyles:()=>lX,progressStyles:()=>lG,provideJupyterDesignSystem:()=>d6,radioGroupStyles:()=>l6,radioStyles:()=>l1,radioTemplate:()=>l2,searchStyles:()=>hr,selectStyles:()=>aY,skeletonStyles:()=>hc,sliderLabelStyles:()=>hF,sliderStyles:()=>hx,strokeWidth:()=>tw,switchStyles:()=>hL,tabPanelStyles:()=>hH,tabStyles:()=>hB,tabsStyles:()=>hG,textAreaStyles:()=>hJ,textFieldStyles:()=>h2,toolbarStyles:()=>h4,tooltipStyles:()=>di,treeItemStyles:()=>dp,treeViewStyles:()=>df,typeRampBaseFontSize:()=>tC,typeRampBaseLineHeight:()=>tI,typeRampMinus1FontSize:()=>tT,typeRampMinus1LineHeight:()=>tF,typeRampMinus2FontSize:()=>tS,typeRampMinus2LineHeight:()=>tO,typeRampPlus1FontSize:()=>tD,typeRampPlus1LineHeight:()=>tE,typeRampPlus2FontSize:()=>tR,typeRampPlus2LineHeight:()=>tL,typeRampPlus3FontSize:()=>tA,typeRampPlus3LineHeight:()=>tV,typeRampPlus4FontSize:()=>tP,typeRampPlus4LineHeight:()=>tH,typeRampPlus5FontSize:()=>tz,typeRampPlus5LineHeight:()=>tM,typeRampPlus6FontSize:()=>tN,typeRampPlus6LineHeight:()=>tB,verticalSliderLabelStyles:()=>hT});class b{constructor(e,t,i,o){this.r=e,this.g=t,this.b=i,this.a="number"!=typeof o||isNaN(o)?1:o}static fromObject(e){return!e||isNaN(e.r)||isNaN(e.g)||isNaN(e.b)?null:new b(e.r,e.g,e.b,e.a)}equalValue(e){return this.r===e.r&&this.g===e.g&&this.b===e.b&&this.a===e.a}toStringHexRGB(){return"#"+[this.r,this.g,this.b].map(this.formatHexValue).join("")}toStringHexRGBA(){return this.toStringHexRGB()+this.formatHexValue(this.a)}toStringHexARGB(){return"#"+[this.a,this.r,this.g,this.b].map(this.formatHexValue).join("")}toStringWebRGB(){return`rgb(${Math.round(c(this.r,0,255))},${Math.round(c(this.g,0,255))},${Math.round(c(this.b,0,255))})`}toStringWebRGBA(){return`rgba(${Math.round(c(this.r,0,255))},${Math.round(c(this.g,0,255))},${Math.round(c(this.b,0,255))},${h(this.a,0,1)})`}roundToPrecision(e){return new b(m(this.r,e),m(this.g,e),m(this.b,e),m(this.a,e))}clamp(){return new b(h(this.r,0,1),h(this.g,0,1),h(this.b,0,1),h(this.a,0,1))}toObject(){return{r:this.r,g:this.g,b:this.b,a:this.a}}formatHexValue(e){let t;return 1===(t=Math.round(h(c(e,0,255),0,255)).toString(16)).length?"0"+t:t}}let f={aliceblue:{r:.941176,g:.972549,b:1},antiquewhite:{r:.980392,g:.921569,b:.843137},aqua:{r:0,g:1,b:1},aquamarine:{r:.498039,g:1,b:.831373},azure:{r:.941176,g:1,b:1},beige:{r:.960784,g:.960784,b:.862745},bisque:{r:1,g:.894118,b:.768627},black:{r:0,g:0,b:0},blanchedalmond:{r:1,g:.921569,b:.803922},blue:{r:0,g:0,b:1},blueviolet:{r:.541176,g:.168627,b:.886275},brown:{r:.647059,g:.164706,b:.164706},burlywood:{r:.870588,g:.721569,b:.529412},cadetblue:{r:.372549,g:.619608,b:.627451},chartreuse:{r:.498039,g:1,b:0},chocolate:{r:.823529,g:.411765,b:.117647},coral:{r:1,g:.498039,b:.313725},cornflowerblue:{r:.392157,g:.584314,b:.929412},cornsilk:{r:1,g:.972549,b:.862745},crimson:{r:.862745,g:.078431,b:.235294},cyan:{r:0,g:1,b:1},darkblue:{r:0,g:0,b:.545098},darkcyan:{r:0,g:.545098,b:.545098},darkgoldenrod:{r:.721569,g:.52549,b:.043137},darkgray:{r:.662745,g:.662745,b:.662745},darkgreen:{r:0,g:.392157,b:0},darkgrey:{r:.662745,g:.662745,b:.662745},darkkhaki:{r:.741176,g:.717647,b:.419608},darkmagenta:{r:.545098,g:0,b:.545098},darkolivegreen:{r:.333333,g:.419608,b:.184314},darkorange:{r:1,g:.54902,b:0},darkorchid:{r:.6,g:.196078,b:.8},darkred:{r:.545098,g:0,b:0},darksalmon:{r:.913725,g:.588235,b:.478431},darkseagreen:{r:.560784,g:.737255,b:.560784},darkslateblue:{r:.282353,g:.239216,b:.545098},darkslategray:{r:.184314,g:.309804,b:.309804},darkslategrey:{r:.184314,g:.309804,b:.309804},darkturquoise:{r:0,g:.807843,b:.819608},darkviolet:{r:.580392,g:0,b:.827451},deeppink:{r:1,g:.078431,b:.576471},deepskyblue:{r:0,g:.74902,b:1},dimgray:{r:.411765,g:.411765,b:.411765},dimgrey:{r:.411765,g:.411765,b:.411765},dodgerblue:{r:.117647,g:.564706,b:1},firebrick:{r:.698039,g:.133333,b:.133333},floralwhite:{r:1,g:.980392,b:.941176},forestgreen:{r:.133333,g:.545098,b:.133333},fuchsia:{r:1,g:0,b:1},gainsboro:{r:.862745,g:.862745,b:.862745},ghostwhite:{r:.972549,g:.972549,b:1},gold:{r:1,g:.843137,b:0},goldenrod:{r:.854902,g:.647059,b:.12549},gray:{r:.501961,g:.501961,b:.501961},green:{r:0,g:.501961,b:0},greenyellow:{r:.678431,g:1,b:.184314},grey:{r:.501961,g:.501961,b:.501961},honeydew:{r:.941176,g:1,b:.941176},hotpink:{r:1,g:.411765,b:.705882},indianred:{r:.803922,g:.360784,b:.360784},indigo:{r:.294118,g:0,b:.509804},ivory:{r:1,g:1,b:.941176},khaki:{r:.941176,g:.901961,b:.54902},lavender:{r:.901961,g:.901961,b:.980392},lavenderblush:{r:1,g:.941176,b:.960784},lawngreen:{r:.486275,g:.988235,b:0},lemonchiffon:{r:1,g:.980392,b:.803922},lightblue:{r:.678431,g:.847059,b:.901961},lightcoral:{r:.941176,g:.501961,b:.501961},lightcyan:{r:.878431,g:1,b:1},lightgoldenrodyellow:{r:.980392,g:.980392,b:.823529},lightgray:{r:.827451,g:.827451,b:.827451},lightgreen:{r:.564706,g:.933333,b:.564706},lightgrey:{r:.827451,g:.827451,b:.827451},lightpink:{r:1,g:.713725,b:.756863},lightsalmon:{r:1,g:.627451,b:.478431},lightseagreen:{r:.12549,g:.698039,b:.666667},lightskyblue:{r:.529412,g:.807843,b:.980392},lightslategray:{r:.466667,g:.533333,b:.6},lightslategrey:{r:.466667,g:.533333,b:.6},lightsteelblue:{r:.690196,g:.768627,b:.870588},lightyellow:{r:1,g:1,b:.878431},lime:{r:0,g:1,b:0},limegreen:{r:.196078,g:.803922,b:.196078},linen:{r:.980392,g:.941176,b:.901961},magenta:{r:1,g:0,b:1},maroon:{r:.501961,g:0,b:0},mediumaquamarine:{r:.4,g:.803922,b:.666667},mediumblue:{r:0,g:0,b:.803922},mediumorchid:{r:.729412,g:.333333,b:.827451},mediumpurple:{r:.576471,g:.439216,b:.858824},mediumseagreen:{r:.235294,g:.701961,b:.443137},mediumslateblue:{r:.482353,g:.407843,b:.933333},mediumspringgreen:{r:0,g:.980392,b:.603922},mediumturquoise:{r:.282353,g:.819608,b:.8},mediumvioletred:{r:.780392,g:.082353,b:.521569},midnightblue:{r:.098039,g:.098039,b:.439216},mintcream:{r:.960784,g:1,b:.980392},mistyrose:{r:1,g:.894118,b:.882353},moccasin:{r:1,g:.894118,b:.709804},navajowhite:{r:1,g:.870588,b:.678431},navy:{r:0,g:0,b:.501961},oldlace:{r:.992157,g:.960784,b:.901961},olive:{r:.501961,g:.501961,b:0},olivedrab:{r:.419608,g:.556863,b:.137255},orange:{r:1,g:.647059,b:0},orangered:{r:1,g:.270588,b:0},orchid:{r:.854902,g:.439216,b:.839216},palegoldenrod:{r:.933333,g:.909804,b:.666667},palegreen:{r:.596078,g:.984314,b:.596078},paleturquoise:{r:.686275,g:.933333,b:.933333},palevioletred:{r:.858824,g:.439216,b:.576471},papayawhip:{r:1,g:.937255,b:.835294},peachpuff:{r:1,g:.854902,b:.72549},peru:{r:.803922,g:.521569,b:.247059},pink:{r:1,g:.752941,b:.796078},plum:{r:.866667,g:.627451,b:.866667},powderblue:{r:.690196,g:.878431,b:.901961},purple:{r:.501961,g:0,b:.501961},red:{r:1,g:0,b:0},rosybrown:{r:.737255,g:.560784,b:.560784},royalblue:{r:.254902,g:.411765,b:.882353},saddlebrown:{r:.545098,g:.270588,b:.07451},salmon:{r:.980392,g:.501961,b:.447059},sandybrown:{r:.956863,g:.643137,b:.376471},seagreen:{r:.180392,g:.545098,b:.341176},seashell:{r:1,g:.960784,b:.933333},sienna:{r:.627451,g:.321569,b:.176471},silver:{r:.752941,g:.752941,b:.752941},skyblue:{r:.529412,g:.807843,b:.921569},slateblue:{r:.415686,g:.352941,b:.803922},slategray:{r:.439216,g:.501961,b:.564706},slategrey:{r:.439216,g:.501961,b:.564706},snow:{r:1,g:.980392,b:.980392},springgreen:{r:0,g:1,b:.498039},steelblue:{r:.27451,g:.509804,b:.705882},tan:{r:.823529,g:.705882,b:.54902},teal:{r:0,g:.501961,b:.501961},thistle:{r:.847059,g:.74902,b:.847059},tomato:{r:1,g:.388235,b:.278431},transparent:{r:0,g:0,b:0,a:0},turquoise:{r:.25098,g:.878431,b:.815686},violet:{r:.933333,g:.509804,b:.933333},wheat:{r:.960784,g:.870588,b:.701961},white:{r:1,g:1,b:1},whitesmoke:{r:.960784,g:.960784,b:.960784},yellow:{r:1,g:1,b:0},yellowgreen:{r:.603922,g:.803922,b:.196078}},v=/^rgb\(\s*((?:(?:25[0-5]|2[0-4]\d|1\d\d|\d{1,2})\s*,\s*){2}(?:25[0-5]|2[0-4]\d|1\d\d|\d{1,2})\s*)\)$/i,y=/^rgba\(\s*((?:(?:25[0-5]|2[0-4]\d|1\d\d|\d{1,2})\s*,\s*){3}(?:0|1|0?\.\d*)\s*)\)$/i,x=/^#((?:[0-9a-f]{6}|[0-9a-f]{3}))$/i,$=/^#((?:[0-9a-f]{8}|[0-9a-f]{4}))$/i;function w(e){let t=x.exec(e);if(null===t)return null;let i=t[1];if(3===i.length){let e=i.charAt(0),t=i.charAt(1),o=i.charAt(2);i=e.concat(e,t,t,o,o)}let o=parseInt(i,16);return isNaN(o)?null:new b(d((0xff0000&o)>>>16,0,255),d((65280&o)>>>8,0,255),d(255&o,0,255),1)}function k(e){let t,i=e.toLowerCase();return x.test(i)?w(i):$.test(i)?function(e){let t=$.exec(e);if(null===t)return null;let i=t[1];if(4===i.length){let e=i.charAt(0),t=i.charAt(1),o=i.charAt(2),s=i.charAt(3);i=e.concat(e,t,t,o,o,s,s)}let o=parseInt(i,16);return isNaN(o)?null:new b(d((0xff0000&o)>>>16,0,255),d((65280&o)>>>8,0,255),d(255&o,0,255),d((0xff000000&o)>>>24,0,255))}(i):v.test(i)?function(e){let t=v.exec(e);if(null===t)return null;let i=t[1].split(",");return new b(d(Number(i[0]),0,255),d(Number(i[1]),0,255),d(Number(i[2]),0,255),1)}(i):y.test(i)?function(e){let t=y.exec(e);if(null===t)return null;let i=t[1].split(",");return 4===i.length?new b(d(Number(i[0]),0,255),d(Number(i[1]),0,255),d(Number(i[2]),0,255),Number(i[3])):null}(i):f.hasOwnProperty(i)&&(t=f[i.toLowerCase()])?new b(t.r,t.g,t.b,t.hasOwnProperty("a")?t.a:void 0):null}class C{constructor(e,t,i){this.h=e,this.s=t,this.l=i}static fromObject(e){return!e||isNaN(e.h)||isNaN(e.s)||isNaN(e.l)?null:new C(e.h,e.s,e.l)}equalValue(e){return this.h===e.h&&this.s===e.s&&this.l===e.l}roundToPrecision(e){return new C(m(this.h,e),m(this.s,e),m(this.l,e))}toObject(){return{h:this.h,s:this.s,l:this.l}}}class I{constructor(e,t,i){this.h=e,this.s=t,this.v=i}static fromObject(e){return!e||isNaN(e.h)||isNaN(e.s)||isNaN(e.v)?null:new I(e.h,e.s,e.v)}equalValue(e){return this.h===e.h&&this.s===e.s&&this.v===e.v}roundToPrecision(e){return new I(m(this.h,e),m(this.s,e),m(this.v,e))}toObject(){return{h:this.h,s:this.s,v:this.v}}}class T{constructor(e,t,i){this.l=e,this.a=t,this.b=i}static fromObject(e){return!e||isNaN(e.l)||isNaN(e.a)||isNaN(e.b)?null:new T(e.l,e.a,e.b)}equalValue(e){return this.l===e.l&&this.a===e.a&&this.b===e.b}roundToPrecision(e){return new T(m(this.l,e),m(this.a,e),m(this.b,e))}toObject(){return{l:this.l,a:this.a,b:this.b}}}T.epsilon=216/24389,T.kappa=24389/27;class F{constructor(e,t,i){this.l=e,this.c=t,this.h=i}static fromObject(e){return!e||isNaN(e.l)||isNaN(e.c)||isNaN(e.h)?null:new F(e.l,e.c,e.h)}equalValue(e){return this.l===e.l&&this.c===e.c&&this.h===e.h}roundToPrecision(e){return new F(m(this.l,e),m(this.c,e),m(this.h,e))}toObject(){return{l:this.l,c:this.c,h:this.h}}}class S{constructor(e,t,i){this.x=e,this.y=t,this.z=i}static fromObject(e){return!e||isNaN(e.x)||isNaN(e.y)||isNaN(e.z)?null:new S(e.x,e.y,e.z)}equalValue(e){return this.x===e.x&&this.y===e.y&&this.z===e.z}roundToPrecision(e){return new S(m(this.x,e),m(this.y,e),m(this.z,e))}toObject(){return{x:this.x,y:this.y,z:this.z}}}function O(e){return .2126*e.r+.7152*e.g+.0722*e.b}function D(e){function t(e){return e<=.03928?e/12.92:Math.pow((e+.055)/1.055,2.4)}return O(new b(t(e.r),t(e.g),t(e.b),1))}S.whitePoint=new S(.95047,1,1.08883);let E=(e,t)=>(e+.05)/(t+.05);function R(e,t){let i=D(e),o=D(t);return i>o?E(i,o):E(o,i)}function L(e){let t=Math.max(e.r,e.g,e.b),i=Math.min(e.r,e.g,e.b),o=t-i,s=0;0!==o&&(s=t===e.r?60*((e.g-e.b)/o%6):t===e.g?60*((e.b-e.r)/o+2):60*((e.r-e.g)/o+4)),s<0&&(s+=360);let r=(t+i)/2,a=0;return 0!==o&&(a=o/(1-Math.abs(2*r-1))),new C(s,a,r)}function A(e,t=1){let i=(1-Math.abs(2*e.l-1))*e.s,o=i*(1-Math.abs(e.h/60%2-1)),s=e.l-i/2,r=0,a=0,n=0;return e.h<60?(r=i,a=o,n=0):e.h<120?(r=o,a=i,n=0):e.h<180?(r=0,a=i,n=o):e.h<240?(r=0,a=o,n=i):e.h<300?(r=o,a=0,n=i):e.h<360&&(r=i,a=0,n=o),new b(r+s,a+s,n+s,t)}function V(e){let t=Math.max(e.r,e.g,e.b),i=t-Math.min(e.r,e.g,e.b),o=0;0!==i&&(o=t===e.r?60*((e.g-e.b)/i%6):t===e.g?60*((e.b-e.r)/i+2):60*((e.r-e.g)/i+4)),o<0&&(o+=360);let s=0;return 0!==t&&(s=i/t),new I(o,s,t)}function P(e){function t(e){return e<=.04045?e/12.92:Math.pow((e+.055)/1.055,2.4)}let i=t(e.r),o=t(e.g),s=t(e.b);return new S(.4124564*i+.3575761*o+.1804375*s,.2126729*i+.7151522*o+.072175*s,.0193339*i+.119192*o+.9503041*s)}function H(e,t=1){function i(e){return e<=.0031308?12.92*e:1.055*Math.pow(e,1/2.4)-.055}return new b(i(3.2404542*e.x-1.5371385*e.y-.4985314*e.z),i(-.969266*e.x+1.8760108*e.y+.041556*e.z),i(.0556434*e.x-.2040259*e.y+1.0572252*e.z),t)}function z(e){return function(e){function t(e){return e>T.epsilon?Math.pow(e,1/3):(T.kappa*e+16)/116}let i=t(e.x/S.whitePoint.x),o=t(e.y/S.whitePoint.y);return new T(116*o-16,500*(i-o),200*(o-t(e.z/S.whitePoint.z)))}(P(e))}function M(e,t=1){let i,o,s,r,a,n,l,h,d;return H((o=(i=(e.l+16)/116)+e.a/500,s=i-e.b/200,r=Math.pow(o,3),a=Math.pow(i,3),n=Math.pow(s,3),l=0,l=r>T.epsilon?r:(116*o-16)/T.kappa,h=0,h=e.l>T.epsilon*T.kappa?a:e.l/T.kappa,d=0,d=n>T.epsilon?n:(116*s-16)/T.kappa,l=S.whitePoint.x*l,h=S.whitePoint.y*h,d=S.whitePoint.z*d,new S(l,h,d)),t)}function N(e){var t=z(e);let i=0;(Math.abs(t.b)>.001||Math.abs(t.a)>.001)&&(i=180/Math.PI*Math.atan2(t.b,t.a)),i<0&&(i+=360);let o=Math.sqrt(t.a*t.a+t.b*t.b);return new F(t.l,o,i)}function B(e,t=1){let i,o;return M((i=0,o=0,0!==e.h&&(i=Math.cos(u(e.h))*e.c,o=Math.sin(u(e.h))*e.c),new T(e.l,i,o)),t)}function j(e,t){let i=e.relativeLuminance>t.relativeLuminance?e:t,o=e.relativeLuminance>t.relativeLuminance?t:e;return(i.relativeLuminance+.05)/(o.relativeLuminance+.05)}let q=Object.freeze({create:(e,t,i)=>new U(e,t,i),from:e=>new U(e.r,e.g,e.b)});class U extends b{constructor(e,t,i){super(e,t,i,1),this.toColorString=this.toStringHexRGB,this.contrast=j.bind(null,this),this.createCSS=this.toColorString,this.relativeLuminance=D(this)}static fromObject(e){return new U(e.r,e.g,e.b)}}function _(e){return q.create(e,e,e)}let G={LightMode:1,DarkMode:.23},W=(-.1+Math.sqrt(.21))/2;function K(e){return e.relativeLuminance<=W}function X(e,t,i,o){var s,r=arguments.length,a=r<3?t:null===o?o=Object.getOwnPropertyDescriptor(t,i):o;if("object"==typeof Reflect&&"function"==typeof Reflect.decorate)a=Reflect.decorate(e,t,i,o);else for(var n=e.length-1;n>=0;n--)(s=e[n])&&(a=(r<3?s(a):r>3?s(t,i,a):s(t,i))||a);return r>3&&a&&Object.defineProperty(t,i,a),a}class Y{createCSS(){return""}createBehavior(){}}let Q=function(){if("u">typeof globalThis)return globalThis;if(void 0!==i.g)return i.g;if("u">typeof self)return self;if("u">typeof window)return window;try{return Function("return this")()}catch(e){return{}}}();void 0===Q.trustedTypes&&(Q.trustedTypes={createPolicy:(e,t)=>t});let Z={configurable:!1,enumerable:!1,writable:!1};void 0===Q.FAST&&Reflect.defineProperty(Q,"FAST",Object.assign({value:Object.create(null)},Z));let J=Q.FAST;if(void 0===J.getById){let e=Object.create(null);Reflect.defineProperty(J,"getById",Object.assign({value(t,i){let o=e[t];return void 0===o&&(o=i?e[t]=i():null),o}},Z))}let ee=Object.freeze([]);function et(){let e=new WeakMap;return function(t){let i=e.get(t);if(void 0===i){let o=Reflect.getPrototypeOf(t);for(;void 0===i&&null!==o;)i=e.get(o),o=Reflect.getPrototypeOf(o);i=void 0===i?[]:i.slice(0),e.set(t,i)}return i}}let ei=Q.FAST.getById(1,()=>{let e=[],t=[];function i(){if(t.length)throw t.shift()}function o(){let o=0;for(;o<e.length;){var s=e[o];try{s.call()}catch(e){t.push(e),setTimeout(i,0)}if(++o>1024){for(let t=0,i=e.length-o;t<i;t++)e[t]=e[t+o];e.length-=o,o=0}}e.length=0}return Object.freeze({enqueue:function(t){e.length<1&&Q.requestAnimationFrame(o),e.push(t)},process:o})}),eo=Q.trustedTypes.createPolicy("fast-html",{createHTML:e=>e}),es=eo,er=`fast-${Math.random().toString(36).substring(2,8)}`,ea=`${er}{`,en=`}${er}`,el=Object.freeze({supportsAdoptedStyleSheets:Array.isArray(document.adoptedStyleSheets)&&"replace"in CSSStyleSheet.prototype,setHTMLPolicy(e){if(es!==eo)throw Error("The HTML policy can only be set once.");es=e},createHTML:e=>es.createHTML(e),isMarker:e=>e&&8===e.nodeType&&e.data.startsWith(er),extractDirectiveIndexFromMarker:e=>parseInt(e.data.replace(`${er}:`,"")),createInterpolationPlaceholder:e=>`${ea}${e}${en}`,createCustomAttributePlaceholder(e,t){return`${e}="${this.createInterpolationPlaceholder(t)}"`},createBlockPlaceholder:e=>`<!--${er}:${e}-->`,queueUpdate:ei.enqueue,processUpdates:ei.process,nextUpdate:()=>new Promise(ei.enqueue),setAttribute(e,t,i){null==i?e.removeAttribute(t):e.setAttribute(t,i)},setBooleanAttribute(e,t,i){i?e.setAttribute(t,""):e.removeAttribute(t)},removeChildNodes(e){for(let t=e.firstChild;null!==t;t=e.firstChild)e.removeChild(t)},createTemplateWalker:e=>document.createTreeWalker(e,133,null,!1)});class eh{constructor(e,t){this.sub1=void 0,this.sub2=void 0,this.spillover=void 0,this.source=e,this.sub1=t}has(e){return void 0===this.spillover?this.sub1===e||this.sub2===e:-1!==this.spillover.indexOf(e)}subscribe(e){let t=this.spillover;if(void 0===t){if(this.has(e))return;if(void 0===this.sub1){this.sub1=e;return}if(void 0===this.sub2){this.sub2=e;return}this.spillover=[this.sub1,this.sub2,e],this.sub1=void 0,this.sub2=void 0}else -1===t.indexOf(e)&&t.push(e)}unsubscribe(e){let t=this.spillover;if(void 0===t)this.sub1===e?this.sub1=void 0:this.sub2===e&&(this.sub2=void 0);else{let i=t.indexOf(e);-1!==i&&t.splice(i,1)}}notify(e){let t=this.spillover,i=this.source;if(void 0===t){let t=this.sub1,o=this.sub2;void 0!==t&&t.handleChange(i,e),void 0!==o&&o.handleChange(i,e)}else for(let o=0,s=t.length;o<s;++o)t[o].handleChange(i,e)}}class ed{constructor(e){this.subscribers={},this.sourceSubscribers=null,this.source=e}notify(e){var t;let i=this.subscribers[e];void 0!==i&&i.notify(e),null==(t=this.sourceSubscribers)||t.notify(e)}subscribe(e,t){var i;if(t){let i=this.subscribers[t];void 0===i&&(this.subscribers[t]=i=new eh(this.source)),i.subscribe(e)}else this.sourceSubscribers=null!=(i=this.sourceSubscribers)?i:new eh(this.source),this.sourceSubscribers.subscribe(e)}unsubscribe(e,t){var i;if(t){let i=this.subscribers[t];void 0!==i&&i.unsubscribe(e)}else null==(i=this.sourceSubscribers)||i.unsubscribe(e)}}let ec=J.getById(2,()=>{let e,t=/(:|&&|\|\||if)/,i=new WeakMap,o=el.queueUpdate,s=e=>{throw Error("Must call enableArrayObservation before observing arrays.")};function r(e){let t=e.$fastController||i.get(e);return void 0===t&&(Array.isArray(e)?t=s(e):i.set(e,t=new ed(e))),t}let a=et();class n{constructor(e){this.name=e,this.field=`_${e}`,this.callback=`${e}Changed`}getValue(t){return void 0!==e&&e.watch(t,this.name),t[this.field]}setValue(e,t){let i=this.field,o=e[i];if(o!==t){e[i]=t;let s=e[this.callback];"function"==typeof s&&s.call(e,o,t),r(e).notify(this.name)}}}class l extends eh{constructor(e,t,i=!1){super(e,t),this.binding=e,this.isVolatileBinding=i,this.needsRefresh=!0,this.needsQueue=!0,this.first=this,this.last=null,this.propertySource=void 0,this.propertyName=void 0,this.notifier=void 0,this.next=void 0}observe(t,i){this.needsRefresh&&null!==this.last&&this.disconnect();let o=e;e=this.needsRefresh?this:void 0,this.needsRefresh=this.isVolatileBinding;let s=this.binding(t,i);return e=o,s}disconnect(){if(null!==this.last){let e=this.first;for(;void 0!==e;)e.notifier.unsubscribe(this,e.propertyName),e=e.next;this.last=null,this.needsRefresh=this.needsQueue=!0}}watch(t,i){let o=this.last,s=r(t),a=null===o?this.first:{};if(a.propertySource=t,a.propertyName=i,a.notifier=s,s.subscribe(this,i),null!==o){if(!this.needsRefresh){let i;e=void 0,i=o.propertySource[o.propertyName],e=this,t===i&&(this.needsRefresh=!0)}o.next=a}this.last=a}handleChange(){this.needsQueue&&(this.needsQueue=!1,o(this))}call(){null!==this.last&&(this.needsQueue=!0,this.notify(this))}records(){let e=this.first;return{next:()=>{let t=e;return void 0===t?{value:void 0,done:!0}:(e=e.next,{value:t,done:!1})},[Symbol.iterator]:function(){return this}}}}return Object.freeze({setArrayObserverFactory(e){s=e},getNotifier:r,track(t,i){void 0!==e&&e.watch(t,i)},trackVolatile(){void 0!==e&&(e.needsRefresh=!0)},notify(e,t){r(e).notify(t)},defineProperty(e,t){"string"==typeof t&&(t=new n(t)),a(e).push(t),Reflect.defineProperty(e,t.name,{enumerable:!0,get:function(){return t.getValue(this)},set:function(e){t.setValue(this,e)}})},getAccessors:a,binding(e,t,i=this.isVolatileBinding(e)){return new l(e,t,i)},isVolatileBinding:e=>t.test(e.toString())})});function eu(e,t){ec.defineProperty(e,t)}let ep=J.getById(3,()=>{let e=null;return{get:()=>e,set(t){e=t}}});class eg{constructor(){this.index=0,this.length=0,this.parent=null,this.parentContext=null}get event(){return ep.get()}get isEven(){return this.index%2==0}get isOdd(){return this.index%2!=0}get isFirst(){return 0===this.index}get isInMiddle(){return!this.isFirst&&!this.isLast}get isLast(){return this.index===this.length-1}static setEvent(e){ep.set(e)}}ec.defineProperty(eg.prototype,"index"),ec.defineProperty(eg.prototype,"length");let em=Object.seal(new eg);class eb{constructor(){this.targets=new WeakSet}addStylesTo(e){this.targets.add(e)}removeStylesFrom(e){this.targets.delete(e)}isAttachedTo(e){return this.targets.has(e)}withBehaviors(...e){return this.behaviors=null===this.behaviors?e:this.behaviors.concat(e),this}}function ef(e){return e.map(e=>e instanceof eb?ef(e.styles):[e]).reduce((e,t)=>e.concat(t),[])}function ev(e){return e.map(e=>e instanceof eb?e.behaviors:null).reduce((e,t)=>null===t?e:(null===e&&(e=[]),e.concat(t)),null)}eb.create=(()=>{if(el.supportsAdoptedStyleSheets){let e=new Map;return t=>new e$(t,e)}return e=>new ek(e)})();let ey=(e,t)=>{e.adoptedStyleSheets=[...e.adoptedStyleSheets,...t]},ex=(e,t)=>{e.adoptedStyleSheets=e.adoptedStyleSheets.filter(e=>-1===t.indexOf(e))};if(el.supportsAdoptedStyleSheets)try{document.adoptedStyleSheets.push(),document.adoptedStyleSheets.splice(),ey=(e,t)=>{e.adoptedStyleSheets.push(...t)},ex=(e,t)=>{for(let i of t){let t=e.adoptedStyleSheets.indexOf(i);-1!==t&&e.adoptedStyleSheets.splice(t,1)}}}catch(e){}class e$ extends eb{constructor(e,t){super(),this.styles=e,this.styleSheetCache=t,this._styleSheets=void 0,this.behaviors=ev(e)}get styleSheets(){if(void 0===this._styleSheets){let e=this.styles,t=this.styleSheetCache;this._styleSheets=ef(e).map(e=>{if(e instanceof CSSStyleSheet)return e;let i=t.get(e);return void 0===i&&((i=new CSSStyleSheet).replaceSync(e),t.set(e,i)),i})}return this._styleSheets}addStylesTo(e){ey(e,this.styleSheets),super.addStylesTo(e)}removeStylesFrom(e){ex(e,this.styleSheets),super.removeStylesFrom(e)}}let ew=0;class ek extends eb{constructor(e){super(),this.styles=e,this.behaviors=null,this.behaviors=ev(e),this.styleSheets=ef(e),this.styleClass=`fast-style-class-${++ew}`}addStylesTo(e){let t=this.styleSheets,i=this.styleClass;e=this.normalizeTarget(e);for(let o=0;o<t.length;o++){let s=document.createElement("style");s.innerHTML=t[o],s.className=i,e.append(s)}super.addStylesTo(e)}removeStylesFrom(e){let t=(e=this.normalizeTarget(e)).querySelectorAll(`.${this.styleClass}`);for(let i=0,o=t.length;i<o;++i)e.removeChild(t[i]);super.removeStylesFrom(e)}isAttachedTo(e){return super.isAttachedTo(this.normalizeTarget(e))}normalizeTarget(e){return e===document?document.body:e}}let eC=Object.freeze({locate:et()}),eI={toView:e=>e?"true":"false",fromView:e=>null!=e&&"false"!==e&&!1!==e&&0!==e},eT={toView(e){if(null==e)return null;let t=+e;return isNaN(t)?null:t.toString()},fromView(e){if(null==e)return null;let t=+e;return isNaN(t)?null:t}};class eF{constructor(e,t,i=t.toLowerCase(),o="reflect",s){this.guards=new Set,this.Owner=e,this.name=t,this.attribute=i,this.mode=o,this.converter=s,this.fieldName=`_${t}`,this.callbackName=`${t}Changed`,this.hasCallback=this.callbackName in e.prototype,"boolean"===o&&void 0===s&&(this.converter=eI)}setValue(e,t){let i=e[this.fieldName],o=this.converter;void 0!==o&&(t=o.fromView(t)),i!==t&&(e[this.fieldName]=t,this.tryReflectToAttribute(e),this.hasCallback&&e[this.callbackName](i,t),e.$fastController.notify(this.name))}getValue(e){return ec.track(e,this.name),e[this.fieldName]}onAttributeChangedCallback(e,t){this.guards.has(e)||(this.guards.add(e),this.setValue(e,t),this.guards.delete(e))}tryReflectToAttribute(e){let t=this.mode,i=this.guards;i.has(e)||"fromView"===t||el.queueUpdate(()=>{i.add(e);let o=e[this.fieldName];switch(t){case"reflect":let s=this.converter;el.setAttribute(e,this.attribute,void 0!==s?s.toView(o):o);break;case"boolean":el.setBooleanAttribute(e,this.attribute,o)}i.delete(e)})}static collect(e,...t){let i=[];t.push(eC.locate(e));for(let o=0,s=t.length;o<s;++o){let s=t[o];if(void 0!==s)for(let t=0,o=s.length;t<o;++t){let o=s[t];"string"==typeof o?i.push(new eF(e,o)):i.push(new eF(e,o.property,o.attribute,o.mode,o.converter))}}return i}}function eS(e,t){let i;function o(e,t){arguments.length>1&&(i.property=t),eC.locate(e.constructor).push(i)}if(arguments.length>1){i={},o(e,t);return}return i=void 0===e?{}:e,o}let eO={mode:"open"},eD={},eE=J.getById(4,()=>{let e=new Map;return Object.freeze({register:t=>!e.has(t.type)&&(e.set(t.type,t),!0),getByType:t=>e.get(t)})});class eR{constructor(e,t=e.definition){"string"==typeof t&&(t={name:t}),this.type=e,this.name=t.name,this.template=t.template;const i=eF.collect(e,t.attributes),o=Array(i.length),s={},r={};for(let e=0,t=i.length;e<t;++e){const t=i[e];o[e]=t.attribute,s[t.name]=t,r[t.attribute]=t}this.attributes=i,this.observedAttributes=o,this.propertyLookup=s,this.attributeLookup=r,this.shadowOptions=void 0===t.shadowOptions?eO:null===t.shadowOptions?void 0:Object.assign(Object.assign({},eO),t.shadowOptions),this.elementOptions=void 0===t.elementOptions?eD:Object.assign(Object.assign({},eD),t.elementOptions),this.styles=void 0===t.styles?void 0:Array.isArray(t.styles)?eb.create(t.styles):t.styles instanceof eb?t.styles:eb.create([t.styles])}get isDefined(){return!!eE.getByType(this.type)}define(e=customElements){let t=this.type;if(eE.register(this)){let e=this.attributes,i=t.prototype;for(let t=0,o=e.length;t<o;++t)ec.defineProperty(i,e[t]);Reflect.defineProperty(t,"observedAttributes",{value:this.observedAttributes,enumerable:!0})}return e.get(this.name)||e.define(this.name,t,this.elementOptions),this}}eR.forType=eE.getByType;let eL=new WeakMap,eA={bubbles:!0,composed:!0,cancelable:!0};function eV(e){return e.shadowRoot||eL.get(e)||null}class eP extends ed{constructor(e,t){super(e),this.boundObservables=null,this.behaviors=null,this.needsInitialization=!0,this._template=null,this._styles=null,this._isConnected=!1,this.$fastController=this,this.view=null,this.element=e,this.definition=t;const i=t.shadowOptions;if(void 0!==i){const t=e.attachShadow(i);"closed"===i.mode&&eL.set(e,t)}const o=ec.getAccessors(e);if(o.length>0){const t=this.boundObservables=Object.create(null);for(let i=0,s=o.length;i<s;++i){const s=o[i].name,r=e[s];void 0!==r&&(delete e[s],t[s]=r)}}}get isConnected(){return ec.track(this,"isConnected"),this._isConnected}setIsConnected(e){this._isConnected=e,ec.notify(this,"isConnected")}get template(){return this._template}set template(e){this._template!==e&&(this._template=e,this.needsInitialization||this.renderTemplate(e))}get styles(){return this._styles}set styles(e){this._styles!==e&&(null!==this._styles&&this.removeStyles(this._styles),this._styles=e,this.needsInitialization||null===e||this.addStyles(e))}addStyles(e){let t=eV(this.element)||this.element.getRootNode();if(e instanceof HTMLStyleElement)t.append(e);else if(!e.isAttachedTo(t)){let i=e.behaviors;e.addStylesTo(t),null!==i&&this.addBehaviors(i)}}removeStyles(e){let t=eV(this.element)||this.element.getRootNode();if(e instanceof HTMLStyleElement)t.removeChild(e);else if(e.isAttachedTo(t)){let i=e.behaviors;e.removeStylesFrom(t),null!==i&&this.removeBehaviors(i)}}addBehaviors(e){let t=this.behaviors||(this.behaviors=new Map),i=e.length,o=[];for(let s=0;s<i;++s){let i=e[s];t.has(i)?t.set(i,t.get(i)+1):(t.set(i,1),o.push(i))}if(this._isConnected){let e=this.element;for(let t=0;t<o.length;++t)o[t].bind(e,em)}}removeBehaviors(e,t=!1){let i=this.behaviors;if(null===i)return;let o=e.length,s=[];for(let r=0;r<o;++r){let o=e[r];if(i.has(o)){let e=i.get(o)-1;0===e||t?i.delete(o)&&s.push(o):i.set(o,e)}}if(this._isConnected){let e=this.element;for(let t=0;t<s.length;++t)s[t].unbind(e)}}onConnectedCallback(){if(this._isConnected)return;let e=this.element;this.needsInitialization?this.finishInitialization():null!==this.view&&this.view.bind(e,em);let t=this.behaviors;if(null!==t)for(let[i]of t)i.bind(e,em);this.setIsConnected(!0)}onDisconnectedCallback(){if(!this._isConnected)return;this.setIsConnected(!1);let e=this.view;null!==e&&e.unbind();let t=this.behaviors;if(null!==t){let e=this.element;for(let[i]of t)i.unbind(e)}}onAttributeChangedCallback(e,t,i){let o=this.definition.attributeLookup[e];void 0!==o&&o.onAttributeChangedCallback(this.element,i)}emit(e,t,i){return!!this._isConnected&&this.element.dispatchEvent(new CustomEvent(e,Object.assign(Object.assign({detail:t},eA),i)))}finishInitialization(){let e=this.element,t=this.boundObservables;if(null!==t){let i=Object.keys(t);for(let o=0,s=i.length;o<s;++o){let s=i[o];e[s]=t[s]}this.boundObservables=null}let i=this.definition;null===this._template&&(this.element.resolveTemplate?this._template=this.element.resolveTemplate():i.template&&(this._template=i.template||null)),null!==this._template&&this.renderTemplate(this._template),null===this._styles&&(this.element.resolveStyles?this._styles=this.element.resolveStyles():i.styles&&(this._styles=i.styles||null)),null!==this._styles&&this.addStyles(this._styles),this.needsInitialization=!1}renderTemplate(e){let t=this.element,i=eV(t)||t;null!==this.view?(this.view.dispose(),this.view=null):this.needsInitialization||el.removeChildNodes(i),e&&(this.view=e.render(t,i,t))}static forCustomElement(e){let t=e.$fastController;if(void 0!==t)return t;let i=eR.forType(e.constructor);if(void 0===i)throw Error("Missing FASTElement definition.");return e.$fastController=new eP(e,i)}}function eH(e){return class extends e{constructor(){super(),eP.forCustomElement(this)}$emit(e,t,i){return this.$fastController.emit(e,t,i)}connectedCallback(){this.$fastController.onConnectedCallback()}disconnectedCallback(){this.$fastController.onDisconnectedCallback()}attributeChangedCallback(e,t,i){this.$fastController.onAttributeChangedCallback(e,t,i)}}}let ez=Object.assign(eH(HTMLElement),{from:e=>eH(e),define:(e,t)=>new eR(e,t).define().type});function eM(e){let t=e.parentElement;if(t)return t;{let t=e.getRootNode();if(t.host instanceof HTMLElement)return t.host}return null}let eN=document.createElement("div");class eB{setProperty(e,t){el.queueUpdate(()=>this.target.setProperty(e,t))}removeProperty(e){el.queueUpdate(()=>this.target.removeProperty(e))}}class ej extends eB{constructor(){super();const e=new CSSStyleSheet;this.target=e.cssRules[e.insertRule(":root{}")].style,document.adoptedStyleSheets=[...document.adoptedStyleSheets,e]}}class eq extends eB{constructor(){super(),this.style=document.createElement("style"),document.head.appendChild(this.style);const{sheet:e}=this.style;if(e){const t=e.insertRule(":root{}",e.cssRules.length);this.target=e.cssRules[t].style}}}class eU{constructor(e){this.store=new Map,this.target=null;const t=e.$fastController;this.style=document.createElement("style"),t.addStyles(this.style),ec.getNotifier(t).subscribe(this,"isConnected"),this.handleChange(t,"isConnected")}targetChanged(){if(null!==this.target)for(let[e,t]of this.store.entries())this.target.setProperty(e,t)}setProperty(e,t){this.store.set(e,t),el.queueUpdate(()=>{null!==this.target&&this.target.setProperty(e,t)})}removeProperty(e){this.store.delete(e),el.queueUpdate(()=>{null!==this.target&&this.target.removeProperty(e)})}handleChange(e,t){let{sheet:i}=this.style;if(i){let e=i.insertRule(":host{}",i.cssRules.length);this.target=i.cssRules[e].style}else this.target=null}}X([eu],eU.prototype,"target",void 0);class e_{constructor(e){this.target=e.style}setProperty(e,t){el.queueUpdate(()=>this.target.setProperty(e,t))}removeProperty(e){el.queueUpdate(()=>this.target.removeProperty(e))}}class eG{setProperty(e,t){for(let i of(eG.properties[e]=t,eG.roots.values()))eX.getOrCreate(eG.normalizeRoot(i)).setProperty(e,t)}removeProperty(e){for(let t of(delete eG.properties[e],eG.roots.values()))eX.getOrCreate(eG.normalizeRoot(t)).removeProperty(e)}static registerRoot(e){let{roots:t}=eG;if(!t.has(e)){t.add(e);let i=eX.getOrCreate(this.normalizeRoot(e));for(let e in eG.properties)i.setProperty(e,eG.properties[e])}}static unregisterRoot(e){let{roots:t}=eG;if(t.has(e)){t.delete(e);let i=eX.getOrCreate(eG.normalizeRoot(e));for(let e in eG.properties)i.removeProperty(e)}}static normalizeRoot(e){return e===eN?document:e}}eG.roots=new Set,eG.properties={};let eW=new WeakMap,eK=el.supportsAdoptedStyleSheets?class extends eB{constructor(e){super();const t=new CSSStyleSheet;this.target=t.cssRules[t.insertRule(":host{}")].style,e.$fastController.addStyles(eb.create([t]))}}:eU,eX=Object.freeze({getOrCreate(e){let t;return eW.has(e)?eW.get(e):(t=e===eN?new eG:e instanceof Document?el.supportsAdoptedStyleSheets?new ej:new eq:e instanceof ez?new eK(e):new e_(e),eW.set(e,t),t)}});class eY extends Y{constructor(e){super(),this.subscribers=new WeakMap,this._appliedTo=new Set,this.name=e.name,null!==e.cssCustomPropertyName&&(this.cssCustomProperty=`--${e.cssCustomPropertyName}`,this.cssVar=`var(${this.cssCustomProperty})`),this.id=eY.uniqueId(),eY.tokensById.set(this.id,this)}get appliedTo(){return[...this._appliedTo]}static from(e){return new eY({name:"string"==typeof e?e:e.name,cssCustomPropertyName:"string"==typeof e?e:void 0===e.cssCustomPropertyName?e.name:e.cssCustomPropertyName})}static isCSSDesignToken(e){return"string"==typeof e.cssCustomProperty}static isDerivedDesignTokenValue(e){return"function"==typeof e}static getTokenById(e){return eY.tokensById.get(e)}getOrCreateSubscriberSet(e=this){return this.subscribers.get(e)||this.subscribers.set(e,new Set)&&this.subscribers.get(e)}createCSS(){return this.cssVar||""}getValueFor(e){let t=e1.getOrCreate(e).get(this);if(void 0!==t)return t;throw Error(`Value could not be retrieved for token named "${this.name}". Ensure the value is set for ${e} or an ancestor of ${e}.`)}setValueFor(e,t){return this._appliedTo.add(e),t instanceof eY&&(t=this.alias(t)),e1.getOrCreate(e).set(this,t),this}deleteValueFor(e){return this._appliedTo.delete(e),e1.existsFor(e)&&e1.getOrCreate(e).delete(this),this}withDefault(e){return this.setValueFor(eN,e),this}subscribe(e,t){let i=this.getOrCreateSubscriberSet(t);t&&!e1.existsFor(t)&&e1.getOrCreate(t),i.has(e)||i.add(e)}unsubscribe(e,t){let i=this.subscribers.get(t||this);i&&i.has(e)&&i.delete(e)}notify(e){let t=Object.freeze({token:this,target:e});this.subscribers.has(this)&&this.subscribers.get(this).forEach(e=>e.handleChange(t)),this.subscribers.has(e)&&this.subscribers.get(e).forEach(e=>e.handleChange(t))}alias(e){return t=>e.getValueFor(t)}}s=0,eY.uniqueId=()=>(s++,s.toString(16)),eY.tokensById=new Map;class eQ{constructor(e,t,i){this.source=e,this.token=t,this.node=i,this.dependencies=new Set,this.observer=ec.binding(e,this,!1),this.observer.handleChange=this.observer.call,this.handleChange()}disconnect(){this.observer.disconnect()}handleChange(){this.node.store.set(this.token,this.observer.observe(this.node.target,em))}}class eZ{constructor(){this.values=new Map}set(e,t){this.values.get(e)!==t&&(this.values.set(e,t),ec.getNotifier(this).notify(e.id))}get(e){return ec.track(this,e.id),this.values.get(e)}delete(e){this.values.delete(e)}all(){return this.values.entries()}}let eJ=new WeakMap,e0=new WeakMap;class e1{constructor(e){this.target=e,this.store=new eZ,this.children=[],this.assignedValues=new Map,this.reflecting=new Set,this.bindingObservers=new Map,this.tokenValueChangeHandler={handleChange:(e,t)=>{let i=eY.getTokenById(t);i&&(i.notify(this.target),this.updateCSSTokenReflection(e,i))}},eJ.set(e,this),ec.getNotifier(this.store).subscribe(this.tokenValueChangeHandler),e instanceof ez?e.$fastController.addBehaviors([this]):e.isConnected&&this.bind()}static getOrCreate(e){return eJ.get(e)||new e1(e)}static existsFor(e){return eJ.has(e)}static findParent(e){if(eN!==e.target){let t=eM(e.target);for(;null!==t;){if(eJ.has(t))return eJ.get(t);t=eM(t)}return e1.getOrCreate(eN)}return null}static findClosestAssignedNode(e,t){let i=t;do{if(i.has(e))return i;i=i.parent?i.parent:i.target!==eN?e1.getOrCreate(eN):null}while(null!==i);return null}get parent(){return e0.get(this)||null}updateCSSTokenReflection(e,t){if(eY.isCSSDesignToken(t)){let i=this.parent,o=this.isReflecting(t);if(i){let s=i.get(t),r=e.get(t);s===r||o?s===r&&o&&this.stopReflectToCSS(t):this.reflectToCSS(t)}else o||this.reflectToCSS(t)}}has(e){return this.assignedValues.has(e)}get(e){let t=this.store.get(e);if(void 0!==t)return t;let i=this.getRaw(e);if(void 0!==i)return this.hydrate(e,i),this.get(e)}getRaw(e){var t;return this.assignedValues.has(e)?this.assignedValues.get(e):null==(t=e1.findClosestAssignedNode(e,this))?void 0:t.getRaw(e)}set(e,t){eY.isDerivedDesignTokenValue(this.assignedValues.get(e))&&this.tearDownBindingObserver(e),this.assignedValues.set(e,t),eY.isDerivedDesignTokenValue(t)?this.setupBindingObserver(e,t):this.store.set(e,t)}delete(e){this.assignedValues.delete(e),this.tearDownBindingObserver(e);let t=this.getRaw(e);t?this.hydrate(e,t):this.store.delete(e)}bind(){let e=e1.findParent(this);for(let t of(e&&e.appendChild(this),this.assignedValues.keys()))t.notify(this.target)}unbind(){this.parent&&e0.get(this).removeChild(this)}appendChild(e){e.parent&&e0.get(e).removeChild(e);let t=this.children.filter(t=>e.contains(t));for(let[i,o]of(e0.set(e,this),this.children.push(e),t.forEach(t=>e.appendChild(t)),ec.getNotifier(this.store).subscribe(e),this.store.all()))e.hydrate(i,this.bindingObservers.has(i)?this.getRaw(i):o)}removeChild(e){let t=this.children.indexOf(e);return -1!==t&&this.children.splice(t,1),ec.getNotifier(this.store).unsubscribe(e),e.parent===this&&e0.delete(e)}contains(e){return function(e,t){let i=t;for(;null!==i;){if(i===e)return!0;i=eM(i)}return!1}(this.target,e.target)}reflectToCSS(e){this.isReflecting(e)||(this.reflecting.add(e),e1.cssCustomPropertyReflector.startReflection(e,this.target))}stopReflectToCSS(e){this.isReflecting(e)&&(this.reflecting.delete(e),e1.cssCustomPropertyReflector.stopReflection(e,this.target))}isReflecting(e){return this.reflecting.has(e)}handleChange(e,t){let i=eY.getTokenById(t);i&&(this.hydrate(i,this.getRaw(i)),this.updateCSSTokenReflection(this.store,i))}hydrate(e,t){if(!this.has(e)){let i=this.bindingObservers.get(e);eY.isDerivedDesignTokenValue(t)?i?i.source!==t&&(this.tearDownBindingObserver(e),this.setupBindingObserver(e,t)):this.setupBindingObserver(e,t):(i&&this.tearDownBindingObserver(e),this.store.set(e,t))}}setupBindingObserver(e,t){let i=new eQ(t,e,this);return this.bindingObservers.set(e,i),i}tearDownBindingObserver(e){return!!this.bindingObservers.has(e)&&(this.bindingObservers.get(e).disconnect(),this.bindingObservers.delete(e),!0)}}e1.cssCustomPropertyReflector=new class{startReflection(e,t){e.subscribe(this,t),this.handleChange({token:e,target:t})}stopReflection(e,t){e.unsubscribe(this,t),this.remove(e,t)}handleChange(e){let{token:t,target:i}=e;this.add(t,i)}add(e,t){eX.getOrCreate(t).setProperty(e.cssCustomProperty,this.resolveCSSValue(e1.getOrCreate(t).get(e)))}remove(e,t){eX.getOrCreate(t).removeProperty(e.cssCustomProperty)}resolveCSSValue(e){return e&&"function"==typeof e.createCSS?e.createCSS():e}},X([eu],e1.prototype,"children",void 0);let e2=Object.freeze({create:function(e){return eY.from(e)},notifyConnection:e=>!!e.isConnected&&!!e1.existsFor(e)&&(e1.getOrCreate(e).bind(),!0),notifyDisconnection:e=>!e.isConnected&&!!e1.existsFor(e)&&(e1.getOrCreate(e).unbind(),!0),registerRoot(e=eN){eG.registerRoot(e)},unregisterRoot(e=eN){eG.unregisterRoot(e)}});function e5(e,t,i=18){let o=N(e),s=o.c+t*i;return s<0&&(s=0),B(new F(o.l,s,o.h))}function e3(e,t){return new b(e.r*t.r,e.g*t.g,e.b*t.b,1)}function e4(e,t){return e<.5?h(2*t*e,0,1):h(1-2*(1-t)*(1-e),0,1)}function e6(e,t){return new b(e4(e.r,t.r),e4(e.g,t.g),e4(e.b,t.b),1)}function e9(e,t,i,o){var s,r,a,n,l,h,d,c,u,m;if(isNaN(e)||e<=0)return i;if(e>=1)return o;switch(t){case nG.HSL:return A((s=L(i),r=L(o),isNaN(e)||e<=0?s:e>=1?r:new C(g(e,s.h,r.h),p(e,s.s,r.s),p(e,s.l,r.l))));case nG.HSV:return function(e,t=1){let i=e.s*e.v,o=i*(1-Math.abs(e.h/60%2-1)),s=e.v-i,r=0,a=0,n=0;return e.h<60?(r=i,a=o,n=0):e.h<120?(r=o,a=i,n=0):e.h<180?(r=0,a=i,n=o):e.h<240?(r=0,a=o,n=i):e.h<300?(r=o,a=0,n=i):e.h<360&&(r=i,a=0,n=o),new b(r+s,a+s,n+s,t)}((a=V(i),n=V(o),isNaN(e)||e<=0?a:e>=1?n:new I(g(e,a.h,n.h),p(e,a.s,n.s),p(e,a.v,n.v))));case nG.XYZ:return H((l=P(i),h=P(o),isNaN(e)||e<=0?l:e>=1?h:new S(p(e,l.x,h.x),p(e,l.y,h.y),p(e,l.z,h.z))));case nG.LAB:return M((d=z(i),c=z(o),isNaN(e)||e<=0?d:e>=1?c:new T(p(e,d.l,c.l),p(e,d.a,c.a),p(e,d.b,c.b))));case nG.LCH:return B((u=N(i),m=N(o),isNaN(e)||e<=0?u:e>=1?m:new F(p(e,u.l,m.l),p(e,u.c,m.c),g(e,u.h,m.h))));default:return isNaN(e)||e<=0?i:e>=1?o:new b(p(e,i.r,o.r),p(e,i.g,o.g),p(e,i.b,o.b),p(e,i.a,o.a))}}(nY=nU||(nU={})).ltr="ltr",nY.rtl="rtl",(nQ=n_||(n_={}))[nQ.Burn=0]="Burn",nQ[nQ.Color=1]="Color",nQ[nQ.Darken=2]="Darken",nQ[nQ.Dodge=3]="Dodge",nQ[nQ.Lighten=4]="Lighten",nQ[nQ.Multiply=5]="Multiply",nQ[nQ.Overlay=6]="Overlay",nQ[nQ.Screen=7]="Screen",(nZ=nG||(nG={}))[nZ.RGB=0]="RGB",nZ[nZ.HSL=1]="HSL",nZ[nZ.HSV=2]="HSV",nZ[nZ.XYZ=3]="XYZ",nZ[nZ.LAB=4]="LAB",nZ[nZ.LCH=5]="LCH";class e8{constructor(e){if(null==e||0===e.length)throw Error("The stops argument must be non-empty");this.stops=this.sortColorScaleStops(e)}static createBalancedColorScale(e){if(null==e||0===e.length)throw Error("The colors argument must be non-empty");let t=Array(e.length);for(let i=0;i<e.length;i++)0===i?t[i]={color:e[i],position:0}:i===e.length-1?t[i]={color:e[i],position:1}:t[i]={color:e[i],position:i*(1/(e.length-1))};return new e8(t)}getColor(e,t=nG.RGB){if(1===this.stops.length||e<=0)return this.stops[0].color;if(e>=1)return this.stops[this.stops.length-1].color;let i=0;for(let t=0;t<this.stops.length;t++)this.stops[t].position<=e&&(i=t);let o=i+1;return o>=this.stops.length&&(o=this.stops.length-1),e9((e-this.stops[i].position)*(1/(this.stops[o].position-this.stops[i].position)),t,this.stops[i].color,this.stops[o].color)}trim(e,t,i=nG.RGB){if(e<0||t>1||t<e)throw Error("Invalid bounds");if(e===t)return new e8([{color:this.getColor(e,i),position:0}]);let o=[];for(let i=0;i<this.stops.length;i++)this.stops[i].position>=e&&this.stops[i].position<=t&&o.push(this.stops[i]);if(0===o.length)return new e8([{color:this.getColor(e),position:e},{color:this.getColor(t),position:t}]);o[0].position!==e&&o.unshift({color:this.getColor(e),position:e}),o[o.length-1].position!==t&&o.push({color:this.getColor(t),position:t});let s=t-e,r=Array(o.length);for(let t=0;t<o.length;t++)r[t]={color:o[t].color,position:(o[t].position-e)/s};return new e8(r)}findNextColor(e,t,i=!1,o=nG.RGB,s=.005,r=32){isNaN(e)||e<=0?e=0:e>=1&&(e=1);let a=this.getColor(e,o),n=+!i;if(R(a,this.getColor(n,o))<=t)return n;let l=i?0:e,h=i?e:0,d=n,c=0;for(;c<=r;){d=Math.abs(h-l)/2+l;let e=R(a,this.getColor(d,o));if(Math.abs(e-t)<=s)break;e>t?i?l=d:h=d:i?h=d:l=d,c++}return d}clone(){let e=Array(this.stops.length);for(let t=0;t<e.length;t++)e[t]={color:this.stops[t].color,position:this.stops[t].position};return new e8(e)}sortColorScaleStops(e){return e.sort((e,t)=>{let i=e.position,o=t.position;return i<o?-1:+(i>o)})}}class e7{constructor(e){this.config=Object.assign({},e7.defaultPaletteConfig,e),this.palette=[],this.updatePaletteColors()}updatePaletteGenerationValues(e){let t=!1;for(let i in e)this.config[i]&&(this.config[i].equalValue?this.config[i].equalValue(e[i])||(this.config[i]=e[i],t=!0):e[i]!==this.config[i]&&(this.config[i]=e[i],t=!0));return t&&this.updatePaletteColors(),t}updatePaletteColors(){let e=this.generatePaletteColorScale();for(let t=0;t<this.config.steps;t++)this.palette[t]=e.getColor(t/(this.config.steps-1),this.config.interpolationMode)}generatePaletteColorScale(){let e=L(this.config.baseColor),t=new e8([{position:0,color:this.config.scaleColorLight},{position:.5,color:this.config.baseColor},{position:1,color:this.config.scaleColorDark}]).trim(this.config.clipLight,1-this.config.clipDark),i=t.getColor(0),o=t.getColor(1),s=i,r=o;if(e.s>=this.config.saturationAdjustmentCutoff&&(s=e5(s,this.config.saturationLight),r=e5(r,this.config.saturationDark)),0!==this.config.multiplyLight){let e=e3(this.config.baseColor,s);s=e9(this.config.multiplyLight,this.config.interpolationMode,s,e)}if(0!==this.config.multiplyDark){let e=e3(this.config.baseColor,r);r=e9(this.config.multiplyDark,this.config.interpolationMode,r,e)}if(0!==this.config.overlayLight){let e=e6(this.config.baseColor,s);s=e9(this.config.overlayLight,this.config.interpolationMode,s,e)}if(0!==this.config.overlayDark){let e=e6(this.config.baseColor,r);r=e9(this.config.overlayDark,this.config.interpolationMode,r,e)}return this.config.baseScalePosition?new e8(this.config.baseScalePosition<=0?[{position:0,color:this.config.baseColor},{position:1,color:r.clamp()}]:this.config.baseScalePosition>=1?[{position:0,color:s.clamp()},{position:1,color:this.config.baseColor}]:[{position:0,color:s.clamp()},{position:this.config.baseScalePosition,color:this.config.baseColor},{position:1,color:r.clamp()}]):new e8([{position:0,color:s.clamp()},{position:.5,color:this.config.baseColor},{position:1,color:r.clamp()}])}}e7.defaultPaletteConfig={baseColor:w("#808080"),steps:11,interpolationMode:nG.RGB,scaleColorLight:new b(1,1,1,1),scaleColorDark:new b(0,0,0,1),clipLight:.185,clipDark:.16,saturationAdjustmentCutoff:.05,saturationLight:.35,saturationDark:1.25,overlayLight:0,overlayDark:.25,multiplyLight:0,multiplyDark:0,baseScalePosition:.5},e7.greyscalePaletteConfig={baseColor:w("#808080"),steps:11,interpolationMode:nG.RGB,scaleColorLight:new b(1,1,1,1),scaleColorDark:new b(0,0,0,1),clipLight:0,clipDark:0,saturationAdjustmentCutoff:0,saturationLight:0,saturationDark:0,overlayLight:0,overlayDark:0,multiplyLight:0,multiplyDark:0,baseScalePosition:.5},e7.defaultPaletteConfig.scaleColorLight,e7.defaultPaletteConfig.scaleColorDark;class te{constructor(e){this.palette=[],this.config=Object.assign({},te.defaultPaletteConfig,e),this.regenPalettes()}regenPalettes(){let e=this.config.steps;(isNaN(e)||e<3)&&(e=3);let t=new b(.14,.14,.14,1),i=new e7(Object.assign(Object.assign({},e7.greyscalePaletteConfig),{baseColor:t,baseScalePosition:.9148936170212766,steps:e})).palette,o=O(this.config.baseColor),s=L(this.config.baseColor).l,r=this.matchRelativeLuminanceIndex((o+s)/2,i)/(e-1),a=this.matchRelativeLuminanceIndex(.14,i)/(e-1),n=L(this.config.baseColor),l=A(C.fromObject({h:n.h,s:n.s,l:.14})),h=A(C.fromObject({h:n.h,s:n.s,l:.06})),d=[,,,,,];d[0]={position:0,color:new b(1,1,1,1)},d[1]={position:r,color:this.config.baseColor},d[2]={position:a,color:l},d[3]={position:.99,color:h},d[4]={position:1,color:new b(0,0,0,1)};let c=new e8(d);this.palette=Array(e);for(let t=0;t<e;t++){let i=c.getColor(t/(e-1),nG.RGB);this.palette[t]=i}}matchRelativeLuminanceIndex(e,t){let i=Number.MAX_VALUE,o=0,s=0,r=t.length;for(;s<r;s++){let r=Math.abs(O(t[s])-e);r<i&&(i=r,o=s)}return o}}function tt(e){return K(e)?-1:1}te.defaultPaletteConfig={baseColor:w("#808080"),steps:94};let ti=Object.freeze({create:function(e,t,i){return"number"==typeof e?ti.from(q.create(e,t,i)):ti.from(e)},from:function(e){return!function(e){let t={r:0,g:0,b:0,toColorString:()=>"",contrast:()=>0,relativeLuminance:0};for(let i in t)if(typeof t[i]!=typeof e[i])return!1;return!0}(e)?to.from(q.create(e.r,e.g,e.b)):to.from(e)}});class to{constructor(e,t){this.closestIndexCache=new Map,this.source=e,this.swatches=t,this.reversedSwatches=Object.freeze([...this.swatches].reverse()),this.lastIndex=this.swatches.length-1}colorContrast(e,t,i,o){void 0===i&&(i=this.closestIndexOf(e));let s=this.swatches,r=this.lastIndex,a=i;return void 0===o&&(o=tt(e)),-1===o&&(s=this.reversedSwatches,a=r-a),function e(t,i,o=0,s=t.length-1){if(s===o)return t[o];let r=Math.floor((s-o)/2)+o;return i(t[r])?e(t,i,o,r):e(t,i,r+1,s)}(s,i=>j(e,i)>=t,a,r)}get(e){return this.swatches[e]||this.swatches[h(e,0,this.lastIndex)]}closestIndexOf(e){if(this.closestIndexCache.has(e.relativeLuminance))return this.closestIndexCache.get(e.relativeLuminance);let t=this.swatches.indexOf(e);if(-1!==t)return this.closestIndexCache.set(e.relativeLuminance,t),t;let i=this.swatches.reduce((t,i)=>Math.abs(i.relativeLuminance-e.relativeLuminance)<Math.abs(t.relativeLuminance-e.relativeLuminance)?i:t);return t=this.swatches.indexOf(i),this.closestIndexCache.set(e.relativeLuminance,t),t}static from(e){return new to(e,Object.freeze(new te({baseColor:b.fromObject(e)}).palette.map(e=>{let t=w(e.toStringHexRGB());return q.create(t.r,t.g,t.b)})))}}let ts=q.create(1,1,1),tr=q.create(0,0,0),ta=q.from(w("#808080")),tn=q.from(w("#DA1A5F")),tl=q.from(w("#D32F2F"));function th(e,t,i,o,s,r){return Math.max(e.closestIndexOf(_(t))+i,o,s,r)}let{create:td}=e2;function tc(e){return e2.create({name:e,cssCustomPropertyName:null})}let tu=td("body-font").withDefault('aktiv-grotesk, "Segoe UI", Arial, Helvetica, sans-serif'),tp=td("base-height-multiplier").withDefault(10),tg=td("base-horizontal-spacing-multiplier").withDefault(3),tm=td("base-layer-luminance").withDefault(G.DarkMode),tb=td("control-corner-radius").withDefault(4),tf=td("density").withDefault(0),tv=td("design-unit").withDefault(4),ty=td("element-scale").withDefault(0),tx=td("direction").withDefault(nU.ltr),t$=td("disabled-opacity").withDefault(.4),tw=td("stroke-width").withDefault(1),tk=td("focus-stroke-width").withDefault(2),tC=td("type-ramp-base-font-size").withDefault("14px"),tI=td("type-ramp-base-line-height").withDefault("20px"),tT=td("type-ramp-minus-1-font-size").withDefault("12px"),tF=td("type-ramp-minus-1-line-height").withDefault("16px"),tS=td("type-ramp-minus-2-font-size").withDefault("10px"),tO=td("type-ramp-minus-2-line-height").withDefault("16px"),tD=td("type-ramp-plus-1-font-size").withDefault("16px"),tE=td("type-ramp-plus-1-line-height").withDefault("24px"),tR=td("type-ramp-plus-2-font-size").withDefault("20px"),tL=td("type-ramp-plus-2-line-height").withDefault("28px"),tA=td("type-ramp-plus-3-font-size").withDefault("28px"),tV=td("type-ramp-plus-3-line-height").withDefault("36px"),tP=td("type-ramp-plus-4-font-size").withDefault("34px"),tH=td("type-ramp-plus-4-line-height").withDefault("44px"),tz=td("type-ramp-plus-5-font-size").withDefault("46px"),tM=td("type-ramp-plus-5-line-height").withDefault("56px"),tN=td("type-ramp-plus-6-font-size").withDefault("60px"),tB=td("type-ramp-plus-6-line-height").withDefault("72px"),tj=tc("accent-fill-rest-delta").withDefault(0),tq=tc("accent-fill-hover-delta").withDefault(4),tU=tc("accent-fill-active-delta").withDefault(-5),t_=tc("accent-fill-focus-delta").withDefault(0),tG=tc("accent-foreground-rest-delta").withDefault(0),tW=tc("accent-foreground-hover-delta").withDefault(6),tK=tc("accent-foreground-active-delta").withDefault(-4),tX=tc("accent-foreground-focus-delta").withDefault(0),tY=tc("neutral-fill-rest-delta").withDefault(7),tQ=tc("neutral-fill-hover-delta").withDefault(10),tZ=tc("neutral-fill-active-delta").withDefault(5),tJ=tc("neutral-fill-focus-delta").withDefault(0),t0=tc("neutral-fill-input-rest-delta").withDefault(0),t1=tc("neutral-fill-input-hover-delta").withDefault(0),t2=tc("neutral-fill-input-active-delta").withDefault(0),t5=tc("neutral-fill-input-focus-delta").withDefault(0),t3=tc("neutral-fill-stealth-rest-delta").withDefault(0),t4=tc("neutral-fill-stealth-hover-delta").withDefault(5),t6=tc("neutral-fill-stealth-active-delta").withDefault(3),t9=tc("neutral-fill-stealth-focus-delta").withDefault(0),t8=tc("neutral-fill-strong-rest-delta").withDefault(0),t7=tc("neutral-fill-strong-hover-delta").withDefault(8),ie=tc("neutral-fill-strong-active-delta").withDefault(-5),it=tc("neutral-fill-strong-focus-delta").withDefault(0),ii=tc("neutral-fill-layer-rest-delta").withDefault(3),io=tc("neutral-stroke-rest-delta").withDefault(25),is=tc("neutral-stroke-hover-delta").withDefault(40),ir=tc("neutral-stroke-active-delta").withDefault(16),ia=tc("neutral-stroke-focus-delta").withDefault(25),il=tc("neutral-stroke-divider-rest-delta").withDefault(8),ih=td("neutral-color").withDefault(ta),id=tc("neutral-palette").withDefault(e=>ti.from(ih.getValueFor(e))),ic=td("accent-color").withDefault(tn),iu=tc("accent-palette").withDefault(e=>ti.from(ic.getValueFor(e))),ip=tc("neutral-layer-card-container-recipe").withDefault({evaluate:e=>{var t,i,o;return t=id.getValueFor(e),i=tm.getValueFor(e),o=ii.getValueFor(e),t.get(t.closestIndexOf(_(i))+o)}}),ig=td("neutral-layer-card-container").withDefault(e=>ip.getValueFor(e).evaluate(e)),im=tc("neutral-layer-floating-recipe").withDefault({evaluate:e=>{var t,i,o;let s;return t=id.getValueFor(e),i=tm.getValueFor(e),o=ii.getValueFor(e),s=t.closestIndexOf(_(i))-o,t.get(s-o)}}),ib=td("neutral-layer-floating").withDefault(e=>im.getValueFor(e).evaluate(e)),iv=tc("neutral-layer-1-recipe").withDefault({evaluate:e=>{var t,i;return t=id.getValueFor(e),i=tm.getValueFor(e),t.get(t.closestIndexOf(_(i)))}}),iy=td("neutral-layer-1").withDefault(e=>iv.getValueFor(e).evaluate(e)),ix=tc("neutral-layer-2-recipe").withDefault({evaluate:e=>{var t,i,o,s,r,a;return t=id.getValueFor(e),i=tm.getValueFor(e),o=ii.getValueFor(e),s=tY.getValueFor(e),r=tQ.getValueFor(e),a=tZ.getValueFor(e),t.get(th(t,i,o,s,r,a))}}),i$=td("neutral-layer-2").withDefault(e=>ix.getValueFor(e).evaluate(e)),iw=tc("neutral-layer-3-recipe").withDefault({evaluate:e=>{var t,i,o,s,r,a;return t=id.getValueFor(e),i=tm.getValueFor(e),o=ii.getValueFor(e),s=tY.getValueFor(e),r=tQ.getValueFor(e),a=tZ.getValueFor(e),t.get(th(t,i,o,s,r,a)+o)}}),ik=td("neutral-layer-3").withDefault(e=>iw.getValueFor(e).evaluate(e)),iC=tc("neutral-layer-4-recipe").withDefault({evaluate:e=>{var t,i,o,s,r,a;return t=id.getValueFor(e),i=tm.getValueFor(e),o=ii.getValueFor(e),s=tY.getValueFor(e),r=tQ.getValueFor(e),a=tZ.getValueFor(e),t.get(th(t,i,o,s,r,a)+2*o)}}),iI=td("neutral-layer-4").withDefault(e=>iC.getValueFor(e).evaluate(e)),iT=td("fill-color").withDefault(e=>iy.getValueFor(e));(nJ=nW||(nW={}))[nJ.normal=4.5]="normal",nJ[nJ.large=7]="large";let iF=td({name:"accent-fill-recipe",cssCustomPropertyName:null}).withDefault({evaluate:(e,t)=>{var i,o,s,r,a,n,l,h,d;let c,u,p,g;return i=iu.getValueFor(e),o=id.getValueFor(e),s=t||iT.getValueFor(e),r=tq.getValueFor(e),a=tU.getValueFor(e),n=t_.getValueFor(e),l=tY.getValueFor(e),h=tQ.getValueFor(e),d=tZ.getValueFor(e),c=i.source,u=o.closestIndexOf(s)>=Math.max(l,h,d)?-1:1,g=(p=i.closestIndexOf(c))+-1*u*r,{rest:i.get(g),hover:i.get(p),active:i.get(g+u*a),focus:i.get(g+u*n)}}}),iS=td("accent-fill-rest").withDefault(e=>iF.getValueFor(e).evaluate(e).rest),iO=td("accent-fill-hover").withDefault(e=>iF.getValueFor(e).evaluate(e).hover),iD=td("accent-fill-active").withDefault(e=>iF.getValueFor(e).evaluate(e).active),iE=td("accent-fill-focus").withDefault(e=>iF.getValueFor(e).evaluate(e).focus),iR=e=>(t,i)=>{var o;return o=i||iS.getValueFor(t),o.contrast(ts)>=e?ts:tr},iL=tc("foreground-on-accent-recipe").withDefault({evaluate:(e,t)=>iR(nW.normal)(e,t)}),iA=td("foreground-on-accent-rest").withDefault(e=>iL.getValueFor(e).evaluate(e,iS.getValueFor(e))),iV=td("foreground-on-accent-hover").withDefault(e=>iL.getValueFor(e).evaluate(e,iO.getValueFor(e))),iP=td("foreground-on-accent-active").withDefault(e=>iL.getValueFor(e).evaluate(e,iD.getValueFor(e))),iH=td("foreground-on-accent-focus").withDefault(e=>iL.getValueFor(e).evaluate(e,iE.getValueFor(e))),iz=tc("foreground-on-accent-large-recipe").withDefault({evaluate:(e,t)=>iR(nW.large)(e,t)}),iM=td("foreground-on-accent-rest-large").withDefault(e=>iz.getValueFor(e).evaluate(e,iS.getValueFor(e))),iN=td("foreground-on-accent-hover-large").withDefault(e=>iz.getValueFor(e).evaluate(e,iO.getValueFor(e))),iB=td("foreground-on-accent-active-large").withDefault(e=>iz.getValueFor(e).evaluate(e,iD.getValueFor(e))),ij=td("foreground-on-accent-focus-large").withDefault(e=>iz.getValueFor(e).evaluate(e,iE.getValueFor(e))),iq=td({name:"accent-foreground-recipe",cssCustomPropertyName:null}).withDefault({evaluate:(e,t)=>{let i;return(i=nW.normal,(e,t)=>{var o,s,r,a,n,l;let h,d,c,u,p,g,m,b,f;return o=iu.getValueFor(e),s=t||iT.getValueFor(e),r=tG.getValueFor(e),a=tW.getValueFor(e),n=tK.getValueFor(e),l=tX.getValueFor(e),c=o.source,u=o.closestIndexOf(c),g=u+(1===(p=tt(s))?Math.min(r,a):Math.max(p*r,p*a)),m=o.colorContrast(s,i,g,p),f=(b=o.closestIndexOf(m))+p*Math.abs(r-a),(1===p?r<a:p*r>p*a)?(h=b,d=f):(h=f,d=b),{rest:o.get(h),hover:o.get(d),active:o.get(h+p*n),focus:o.get(h+p*l)}})(e,t)}}),iU=td("accent-foreground-rest").withDefault(e=>iq.getValueFor(e).evaluate(e).rest),i_=td("accent-foreground-hover").withDefault(e=>iq.getValueFor(e).evaluate(e).hover),iG=td("accent-foreground-active").withDefault(e=>iq.getValueFor(e).evaluate(e).active),iW=td("accent-foreground-focus").withDefault(e=>iq.getValueFor(e).evaluate(e).focus),iK=td({name:"neutral-fill-recipe",cssCustomPropertyName:null}).withDefault({evaluate:(e,t)=>{var i,o,s,r,a,n;let l,h;return i=id.getValueFor(e),o=t||iT.getValueFor(e),s=tY.getValueFor(e),r=tQ.getValueFor(e),a=tZ.getValueFor(e),n=tJ.getValueFor(e),h=(l=i.closestIndexOf(o))>=Math.max(s,r,a,n)?-1:1,{rest:i.get(l+h*s),hover:i.get(l+h*r),active:i.get(l+h*a),focus:i.get(l+h*n)}}}),iX=td("neutral-fill-rest").withDefault(e=>iK.getValueFor(e).evaluate(e).rest),iY=td("neutral-fill-hover").withDefault(e=>iK.getValueFor(e).evaluate(e).hover),iQ=td("neutral-fill-active").withDefault(e=>iK.getValueFor(e).evaluate(e).active),iZ=td("neutral-fill-focus").withDefault(e=>iK.getValueFor(e).evaluate(e).focus),iJ=td({name:"neutral-fill-input-recipe",cssCustomPropertyName:null}).withDefault({evaluate:(e,t)=>{var i,o,s,r,a,n;let l,h;return i=id.getValueFor(e),o=t||iT.getValueFor(e),s=t0.getValueFor(e),r=t1.getValueFor(e),a=t2.getValueFor(e),n=t5.getValueFor(e),l=tt(o),h=i.closestIndexOf(o),{rest:i.get(h-l*s),hover:i.get(h-l*r),active:i.get(h-l*a),focus:i.get(h-l*n)}}}),i0=td("neutral-fill-input-rest").withDefault(e=>iJ.getValueFor(e).evaluate(e).rest),i1=td("neutral-fill-input-hover").withDefault(e=>iJ.getValueFor(e).evaluate(e).hover),i2=td("neutral-fill-input-active").withDefault(e=>iJ.getValueFor(e).evaluate(e).active),i5=td("neutral-fill-input-focus").withDefault(e=>iJ.getValueFor(e).evaluate(e).focus),i3=td({name:"neutral-fill-stealth-recipe",cssCustomPropertyName:null}).withDefault({evaluate:(e,t)=>{var i,o,s,r,a,n,l,h;let d,c,u;return i=id.getValueFor(e),o=t||iT.getValueFor(e),s=t3.getValueFor(e),r=t4.getValueFor(e),a=t6.getValueFor(e),n=t9.getValueFor(e),l=tY.getValueFor(e),h=tQ.getValueFor(e),d=Math.max(s,r,a,n,l,h,tZ.getValueFor(e),tJ.getValueFor(e)),u=(c=i.closestIndexOf(o))>=d?-1:1,{rest:i.get(c+u*s),hover:i.get(c+u*r),active:i.get(c+u*a),focus:i.get(c+u*n)}}}),i4=td("neutral-fill-stealth-rest").withDefault(e=>i3.getValueFor(e).evaluate(e).rest),i6=td("neutral-fill-stealth-hover").withDefault(e=>i3.getValueFor(e).evaluate(e).hover),i9=td("neutral-fill-stealth-active").withDefault(e=>i3.getValueFor(e).evaluate(e).active),i8=td("neutral-fill-stealth-focus").withDefault(e=>i3.getValueFor(e).evaluate(e).focus),i7=td({name:"neutral-fill-strong-recipe",cssCustomPropertyName:null}).withDefault({evaluate:(e,t)=>{var i,o,s,r,a,n;let l,h,d,c,u;return i=id.getValueFor(e),o=t||iT.getValueFor(e),s=t8.getValueFor(e),r=t7.getValueFor(e),a=ie.getValueFor(e),n=it.getValueFor(e),d=tt(o),u=(c=i.closestIndexOf(i.colorContrast(o,4.5)))+d*Math.abs(s-r),(1===d?s<r:d*s>d*r)?(l=c,h=u):(l=u,h=c),{rest:i.get(l),hover:i.get(h),active:i.get(l+d*a),focus:i.get(l+d*n)}}}),oe=td("neutral-fill-strong-rest").withDefault(e=>i7.getValueFor(e).evaluate(e).rest),ot=td("neutral-fill-strong-hover").withDefault(e=>i7.getValueFor(e).evaluate(e).hover),oi=td("neutral-fill-strong-active").withDefault(e=>i7.getValueFor(e).evaluate(e).active),oo=td("neutral-fill-strong-focus").withDefault(e=>i7.getValueFor(e).evaluate(e).focus),os=tc("neutral-fill-layer-recipe").withDefault({evaluate:(e,t)=>{var i,o,s;let r;return i=id.getValueFor(e),o=t||iT.getValueFor(e),s=ii.getValueFor(e),r=i.closestIndexOf(o),i.get(r-(r<s?-1*s:s))}}),or=td("neutral-fill-layer-rest").withDefault(e=>os.getValueFor(e).evaluate(e)),oa=tc("focus-stroke-outer-recipe").withDefault({evaluate:e=>{var t,i;return t=id.getValueFor(e),i=iT.getValueFor(e),t.colorContrast(i,3.5)}}),on=td("focus-stroke-outer").withDefault(e=>oa.getValueFor(e).evaluate(e)),ol=tc("focus-stroke-inner-recipe").withDefault({evaluate:e=>{var t,i,o;return t=iu.getValueFor(e),i=iT.getValueFor(e),o=on.getValueFor(e),t.colorContrast(o,3.5,t.closestIndexOf(t.source),-1*tt(i))}}),oh=td("focus-stroke-inner").withDefault(e=>ol.getValueFor(e).evaluate(e)),od=tc("neutral-foreground-hint-recipe").withDefault({evaluate:e=>{var t,i;return t=id.getValueFor(e),i=iT.getValueFor(e),t.colorContrast(i,4.5)}}),oc=td("neutral-foreground-hint").withDefault(e=>od.getValueFor(e).evaluate(e)),ou=tc("neutral-foreground-recipe").withDefault({evaluate:e=>{var t,i;return t=id.getValueFor(e),i=iT.getValueFor(e),t.colorContrast(i,14)}}),op=td("neutral-foreground-rest").withDefault(e=>ou.getValueFor(e).evaluate(e)),og=td({name:"neutral-stroke-recipe",cssCustomPropertyName:null}).withDefault({evaluate:e=>{var t,i,o,s,r,a;let n,l,h;return t=id.getValueFor(e),i=iT.getValueFor(e),o=io.getValueFor(e),s=is.getValueFor(e),r=ir.getValueFor(e),a=ia.getValueFor(e),n=t.closestIndexOf(i),h=n+(l=tt(i))*o,{rest:t.get(h),hover:t.get(h+l*(s-o)),active:t.get(h+l*(r-o)),focus:t.get(h+l*(a-o))}}}),om=td("neutral-stroke-rest").withDefault(e=>og.getValueFor(e).evaluate(e).rest),ob=td("neutral-stroke-hover").withDefault(e=>og.getValueFor(e).evaluate(e).hover),of=td("neutral-stroke-active").withDefault(e=>og.getValueFor(e).evaluate(e).active),ov=td("neutral-stroke-focus").withDefault(e=>og.getValueFor(e).evaluate(e).focus),oy=tc("neutral-stroke-divider-recipe").withDefault({evaluate:(e,t)=>{var i,o,s;return i=id.getValueFor(e),o=t||iT.getValueFor(e),s=il.getValueFor(e),i.get(i.closestIndexOf(o)+tt(o)*s)}}),ox=td("neutral-stroke-divider-rest").withDefault(e=>oy.getValueFor(e).evaluate(e)),o$=e2.create({name:"height-number",cssCustomPropertyName:null}).withDefault(e=>(tp.getValueFor(e)+tf.getValueFor(e))*tv.getValueFor(e)),ow=td("error-color").withDefault(tl),ok=tc("error-palette").withDefault(e=>ti.from(ow.getValueFor(e))),oC=td({name:"error-fill-recipe",cssCustomPropertyName:null}).withDefault({evaluate:(e,t)=>{var i,o,s,r,a,n,l,h,d;let c,u,p,g;return i=ok.getValueFor(e),o=id.getValueFor(e),s=t||iT.getValueFor(e),r=tq.getValueFor(e),a=tU.getValueFor(e),n=t_.getValueFor(e),l=tY.getValueFor(e),h=tQ.getValueFor(e),d=tZ.getValueFor(e),c=i.source,u=o.closestIndexOf(s)>=Math.max(l,h,d)?-1:1,g=(p=i.closestIndexOf(c))+-1*u*r,{rest:i.get(g),hover:i.get(p),active:i.get(g+u*a),focus:i.get(g+u*n)}}}),oI=td("error-fill-rest").withDefault(e=>oC.getValueFor(e).evaluate(e).rest),oT=td("error-fill-hover").withDefault(e=>oC.getValueFor(e).evaluate(e).hover),oF=td("error-fill-active").withDefault(e=>oC.getValueFor(e).evaluate(e).active),oS=td("error-fill-focus").withDefault(e=>oC.getValueFor(e).evaluate(e).focus),oO=e=>(t,i)=>{var o;return o=i||oI.getValueFor(t),o.contrast(ts)>=e?ts:tr},oD=td({name:"foreground-on-error-recipe",cssCustomPropertyName:null}).withDefault({evaluate:(e,t)=>oO(nW.normal)(e,t)}),oE=td("foreground-on-error-rest").withDefault(e=>oD.getValueFor(e).evaluate(e,oI.getValueFor(e))),oR=td("foreground-on-error-hover").withDefault(e=>oD.getValueFor(e).evaluate(e,oT.getValueFor(e))),oL=td("foreground-on-error-active").withDefault(e=>oD.getValueFor(e).evaluate(e,oF.getValueFor(e))),oA=td("foreground-on-error-focus").withDefault(e=>oD.getValueFor(e).evaluate(e,oS.getValueFor(e))),oV=td({name:"foreground-on-error-large-recipe",cssCustomPropertyName:null}).withDefault({evaluate:(e,t)=>oO(nW.large)(e,t)}),oP=td("foreground-on-error-rest-large").withDefault(e=>oV.getValueFor(e).evaluate(e,oI.getValueFor(e))),oH=td("foreground-on-error-hover-large").withDefault(e=>oV.getValueFor(e).evaluate(e,oT.getValueFor(e))),oz=td("foreground-on-error-active-large").withDefault(e=>oV.getValueFor(e).evaluate(e,oF.getValueFor(e))),oM=td("foreground-on-error-focus-large").withDefault(e=>oV.getValueFor(e).evaluate(e,oS.getValueFor(e))),oN=td({name:"error-foreground-recipe",cssCustomPropertyName:null}).withDefault({evaluate:(e,t)=>{let i;return(i=nW.normal,(e,t)=>{var o,s,r,a,n,l;let h,d,c,u,p,g,m,b,f;return o=ok.getValueFor(e),s=t||iT.getValueFor(e),r=tG.getValueFor(e),a=tW.getValueFor(e),n=tK.getValueFor(e),l=tX.getValueFor(e),c=o.source,u=o.closestIndexOf(c),g=u+(1==(p=K(s)?-1:1)?Math.min(r,a):Math.max(p*r,p*a)),m=o.colorContrast(s,i,g,p),f=(b=o.closestIndexOf(m))+p*Math.abs(r-a),(1===p?r<a:p*r>p*a)?(h=b,d=f):(h=f,d=b),{rest:o.get(h),hover:o.get(d),active:o.get(h+p*n),focus:o.get(h+p*l)}})(e,t)}}),oB=td("error-foreground-rest").withDefault(e=>oN.getValueFor(e).evaluate(e).rest),oj=td("error-foreground-hover").withDefault(e=>oN.getValueFor(e).evaluate(e).hover),oq=td("error-foreground-active").withDefault(e=>oN.getValueFor(e).evaluate(e).active),oU=td("error-foreground-focus").withDefault(e=>oN.getValueFor(e).evaluate(e).focus),o_="--jp-layout-color1",oG=!1;function oW(){let e;oG||(oG=!0,e=()=>{new MutationObserver(()=>{oY()}).observe(document.body,{attributes:!0,attributeFilter:["data-jp-theme-name"],childList:!1,characterData:!1}),oY()},"complete"===document.readyState?e():window.addEventListener("load",e))}let oK=e=>{let t=parseInt(e,10);return isNaN(t)?null:t},oX={"--jp-border-width":{converter:oK,token:tw},"--jp-border-radius":{converter:oK,token:tb},[o_]:{converter:(e,t)=>{let i=k(e);if(!i)return null;{let e=L(i),t=A(C.fromObject({h:e.h,s:e.s,l:.5}));return q.create(t.r,t.g,t.b)}},token:ih},"--jp-brand-color1":{converter:(e,t)=>{let i=k(e);if(!i)return null;{let e=L(i),o=A(C.fromObject({h:e.h,s:e.s,l:e.l+(t?1:-1)*tq.getValueFor(document.body)/94}));return q.create(o.r,o.g,o.b)}},token:ic},"--jp-error-color1":{converter:(e,t)=>{let i=k(e);if(!i)return null;{let e=L(i),o=A(C.fromObject({h:e.h,s:e.s,l:e.l+(t?1:-1)*tq.getValueFor(document.body)/94}));return q.create(o.r,o.g,o.b)}},token:ow},"--jp-ui-font-family":{token:tu},"--jp-ui-font-size1":{token:tC}};function oY(){var e;let t=getComputedStyle(document.body),i=document.body.getAttribute("data-jp-theme-light"),o=!1;if(i)o="false"===i;else{let e=t.getPropertyValue(o_).toString();if(e){let t=k(e);t&&(o=K(q.create(t.r,t.g,t.b)),console.debug(`Theme is ${o?"dark":"light"} based on '${o_}' value: ${e}.`))}}for(let i in tm.setValueFor(document.body,o?G.DarkMode:G.LightMode),oX){let s=oX[i],r=t.getPropertyValue(i).toString();if(document.body&&""!==r){let t=(null!=(e=s.converter)?e:e=>e)(r.trim(),o);null!==t?s.token.setValueFor(document.body,t):console.error(`Fail to parse value '${r}' for '${i}' as FAST design token.`)}}}(n0=nK||(nK={}))[n0.alt=18]="alt",n0[n0.arrowDown=40]="arrowDown",n0[n0.arrowLeft=37]="arrowLeft",n0[n0.arrowRight=39]="arrowRight",n0[n0.arrowUp=38]="arrowUp",n0[n0.back=8]="back",n0[n0.backSlash=220]="backSlash",n0[n0.break=19]="break",n0[n0.capsLock=20]="capsLock",n0[n0.closeBracket=221]="closeBracket",n0[n0.colon=186]="colon",n0[n0.colon2=59]="colon2",n0[n0.comma=188]="comma",n0[n0.ctrl=17]="ctrl",n0[n0.delete=46]="delete",n0[n0.end=35]="end",n0[n0.enter=13]="enter",n0[n0.equals=187]="equals",n0[n0.equals2=61]="equals2",n0[n0.equals3=107]="equals3",n0[n0.escape=27]="escape",n0[n0.forwardSlash=191]="forwardSlash",n0[n0.function1=112]="function1",n0[n0.function10=121]="function10",n0[n0.function11=122]="function11",n0[n0.function12=123]="function12",n0[n0.function2=113]="function2",n0[n0.function3=114]="function3",n0[n0.function4=115]="function4",n0[n0.function5=116]="function5",n0[n0.function6=117]="function6",n0[n0.function7=118]="function7",n0[n0.function8=119]="function8",n0[n0.function9=120]="function9",n0[n0.home=36]="home",n0[n0.insert=45]="insert",n0[n0.menu=93]="menu",n0[n0.minus=189]="minus",n0[n0.minus2=109]="minus2",n0[n0.numLock=144]="numLock",n0[n0.numPad0=96]="numPad0",n0[n0.numPad1=97]="numPad1",n0[n0.numPad2=98]="numPad2",n0[n0.numPad3=99]="numPad3",n0[n0.numPad4=100]="numPad4",n0[n0.numPad5=101]="numPad5",n0[n0.numPad6=102]="numPad6",n0[n0.numPad7=103]="numPad7",n0[n0.numPad8=104]="numPad8",n0[n0.numPad9=105]="numPad9",n0[n0.numPadDivide=111]="numPadDivide",n0[n0.numPadDot=110]="numPadDot",n0[n0.numPadMinus=109]="numPadMinus",n0[n0.numPadMultiply=106]="numPadMultiply",n0[n0.numPadPlus=107]="numPadPlus",n0[n0.openBracket=219]="openBracket",n0[n0.pageDown=34]="pageDown",n0[n0.pageUp=33]="pageUp",n0[n0.period=190]="period",n0[n0.print=44]="print",n0[n0.quote=222]="quote",n0[n0.scrollLock=145]="scrollLock",n0[n0.shift=16]="shift",n0[n0.space=32]="space",n0[n0.tab=9]="tab",n0[n0.tilde=192]="tilde",n0[n0.windowsLeft=91]="windowsLeft",n0[n0.windowsOpera=219]="windowsOpera",n0[n0.windowsRight=92]="windowsRight";let oQ="ArrowDown",oZ="ArrowLeft",oJ="ArrowRight",o0="ArrowUp",o1={ArrowDown:oQ,ArrowLeft:oZ,ArrowRight:oJ,ArrowUp:o0};function o2(e,t,i){return Math.min(Math.max(i,e),t)}function o5(e,t,i=0){return[t,i]=[t,i].sort((e,t)=>e-t),t<=e&&e<i}let o3=new Map;"metadata"in Reflect||(Reflect.metadata=function(e,t){return function(i){Reflect.defineMetadata(e,t,i)}},Reflect.defineMetadata=function(e,t,i){let o=o3.get(i);void 0===o&&o3.set(i,o=new Map),o.set(e,t)},Reflect.getOwnMetadata=function(e,t){let i=o3.get(t);if(void 0!==i)return i.get(e)});class o4{constructor(e,t){this.container=e,this.key=t}instance(e){return this.registerResolver(0,e)}singleton(e){return this.registerResolver(1,e)}transient(e){return this.registerResolver(2,e)}callback(e){return this.registerResolver(3,e)}cachedCallback(e){return this.registerResolver(3,sy(e))}aliasTo(e){return this.registerResolver(5,e)}registerResolver(e,t){let{container:i,key:o}=this;return this.container=this.key=void 0,i.registerResolver(o,new sn(o,e,t))}}function o6(e){let t,i=e.slice(),o=Object.keys(e),s=o.length;for(let r=0;r<s;++r)sF(t=o[r])||(i[t]=e[t]);return i}let o9=Object.freeze({default:Object.freeze({parentLocator:()=>null,responsibleForOwnerRequests:!1,defaultResolver:Object.freeze({none(e){throw Error(`${e.toString()} not registered, did you forget to add @singleton()?`)},singleton:e=>new sn(e,1,e),transient:e=>new sn(e,2,e)}).singleton})}),o8=new Map;function o7(e){return t=>Reflect.getOwnMetadata(e,t)}let se=null,st=Object.freeze({createContainer:e=>new sf(null,Object.assign({},o9.default,e)),findResponsibleContainer(e){let t=e.$$container$$;return t&&t.responsibleForOwnerRequests?t:st.findParentContainer(e)},findParentContainer(e){let t=new CustomEvent(sm,{bubbles:!0,composed:!0,cancelable:!0,detail:{container:void 0}});return e.dispatchEvent(t),t.detail.container||st.getOrCreateDOMContainer()},getOrCreateDOMContainer:(e,t)=>e?e.$$container$$||new sf(e,Object.assign({},o9.default,t,{parentLocator:st.findParentContainer})):se||(se=new sf(null,Object.assign({},o9.default,t,{parentLocator:()=>null}))),getDesignParamtypes:o7("design:paramtypes"),getAnnotationParamtypes:o7("di:paramtypes"),getOrCreateAnnotationParamTypes(e){let t=this.getAnnotationParamtypes(e);return void 0===t&&Reflect.defineMetadata("di:paramtypes",t=[],e),t},getDependencies(e){let t=o8.get(e);if(void 0===t){let i=e.inject;if(void 0===i){let i=st.getDesignParamtypes(e),o=st.getAnnotationParamtypes(e);if(void 0===i)if(void 0===o){let i=Object.getPrototypeOf(e);t="function"==typeof i&&i!==Function.prototype?o6(st.getDependencies(i)):[]}else t=o6(o);else if(void 0===o)t=o6(i);else{let e,s;t=o6(i);let r=o.length;for(let i=0;i<r;++i)void 0!==(e=o[i])&&(t[i]=e);let a=Object.keys(o);r=a.length;for(let e=0;e<r;++e)sF(s=a[e])||(t[s]=o[s])}}else t=o6(i);o8.set(e,t)}return t},defineProperty(e,t,i,o=!1){let s=`$di_${t}`;Reflect.defineProperty(e,t,{get:function(){let e=this[s];if(void 0===e&&(e=(this instanceof HTMLElement?st.findResponsibleContainer(this):st.getOrCreateDOMContainer()).get(i),this[s]=e,o&&this instanceof ez)){let o=this.$fastController,r=()=>{st.findResponsibleContainer(this).get(i)!==this[s]&&(this[s]=e,o.notify(t))};o.subscribe({handleChange:r},"isConnected")}return e}})},createInterface(e,t){let i="function"==typeof e?e:t,o="string"==typeof e?e:e&&"friendlyName"in e&&e.friendlyName||sk,s="string"!=typeof e&&!!e&&"respectConnection"in e&&(e.respectConnection||!1),r=function(e,t,i){if(null==e||new.target!==void 0)throw Error(`No registration for interface: '${r.friendlyName}'`);t?st.defineProperty(e,t,r,s):st.getOrCreateAnnotationParamTypes(e)[i]=r};return r.$isInterface=!0,r.friendlyName=null==o?"(anonymous)":o,null!=i&&(r.register=function(e,t){return i(new o4(e,null!=t?t:r))}),r.toString=function(){return`InterfaceSymbol<${r.friendlyName}>`},r},inject:(...e)=>function(t,i,o){if("number"==typeof o){let i=st.getOrCreateAnnotationParamTypes(t),s=e[0];void 0!==s&&(i[o]=s)}else if(i)st.defineProperty(t,i,e[0]);else{let i,s=o?st.getOrCreateAnnotationParamTypes(o.value):st.getOrCreateAnnotationParamTypes(t);for(let t=0;t<e.length;++t)void 0!==(i=e[t])&&(s[t]=i)}},transient:e=>(e.register=function(t){return sx.transient(e,e).register(t)},e.registerInRequestor=!1,e),singleton:(e,t=ss)=>(e.register=function(t){return sx.singleton(e,e).register(t)},e.registerInRequestor=t.scoped,e)}),si=st.createInterface("Container");function so(e){return function(t){let i=function(e,t,o){st.inject(i)(e,t,o)};return i.$isResolver=!0,i.resolve=function(i,o){return e(t,i,o)},i}}st.inject;let ss={scoped:!1};function sr(e,t,i){st.inject(sr)(e,t,i)}function sa(e,t){return t.getFactory(e).construct(t)}so((e,t,i)=>()=>i.get(e)),so((e,t,i)=>i.has(e,!0)?i.get(e):void 0),sr.$isResolver=!0,sr.resolve=()=>void 0,so((e,t,i)=>{let o=sa(e,t),s=new sn(e,0,o);return i.registerResolver(e,s),o}),so((e,t,i)=>sa(e,t));class sn{constructor(e,t,i){this.key=e,this.strategy=t,this.state=i,this.resolving=!1}get $isResolver(){return!0}register(e){return e.registerResolver(this.key,this)}resolve(e,t){switch(this.strategy){case 0:return this.state;case 1:if(this.resolving)throw Error(`Cyclic dependency found: ${this.state.name}`);return this.resolving=!0,this.state=e.getFactory(this.state).construct(t),this.strategy=0,this.resolving=!1,this.state;case 2:{let i=e.getFactory(this.state);if(null===i)throw Error(`Resolver for ${String(this.key)} returned a null factory`);return i.construct(t)}case 3:return this.state(e,t,this);case 4:return this.state[0].resolve(e,t);case 5:return t.get(this.state);default:throw Error(`Invalid resolver strategy specified: ${this.strategy}.`)}}getFactory(e){var t,i,o;switch(this.strategy){case 1:case 2:return e.getFactory(this.state);case 5:return null!=(o=null==(i=null==(t=e.getResolver(this.state))?void 0:t.getFactory)?void 0:i.call(t,e))?o:null;default:return null}}}function sl(e){return this.get(e)}function sh(e,t){return t(e)}class sd{constructor(e,t){this.Type=e,this.dependencies=t,this.transformers=null}construct(e,t){let i;return(i=void 0===t?new this.Type(...this.dependencies.map(sl,e)):new this.Type(...this.dependencies.map(sl,e),...t),null==this.transformers)?i:this.transformers.reduce(sh,i)}registerTransformer(e){(this.transformers||(this.transformers=[])).push(e)}}let sc={$isResolver:!0,resolve:(e,t)=>t};function su(e){return"function"==typeof e.register}function sp(e){return su(e)&&"boolean"==typeof e.registerInRequestor&&e.registerInRequestor}let sg=new Set(["Array","ArrayBuffer","Boolean","DataView","Date","Error","EvalError","Float32Array","Float64Array","Function","Int8Array","Int16Array","Int32Array","Map","Number","Object","Promise","RangeError","ReferenceError","RegExp","Set","SharedArrayBuffer","String","SyntaxError","TypeError","Uint8Array","Uint8ClampedArray","Uint16Array","Uint32Array","URIError","WeakMap","WeakSet"]),sm="__DI_LOCATE_PARENT__",sb=new Map;class sf{constructor(e,t){this.owner=e,this.config=t,this._parent=void 0,this.registerDepth=0,this.context=null,null!==e&&(e.$$container$$=this),this.resolvers=new Map,this.resolvers.set(si,sc),e instanceof Node&&e.addEventListener(sm,e=>{e.composedPath()[0]!==this.owner&&(e.detail.container=this,e.stopImmediatePropagation())})}get parent(){return void 0===this._parent&&(this._parent=this.config.parentLocator(this.owner)),this._parent}get depth(){return null===this.parent?0:this.parent.depth+1}get responsibleForOwnerRequests(){return this.config.responsibleForOwnerRequests}registerWithContext(e,...t){return this.context=e,this.register(...t),this.context=null,this}register(...e){let t,i,o,s,r;if(100==++this.registerDepth)throw Error("Unable to autoregister dependency");let a=this.context;for(let n=0,l=e.length;n<l;++n)if(sC(t=e[n]))if(su(t))t.register(this,a);else if(void 0!==t.prototype)sx.singleton(t,t).register(this);else for(i=Object.keys(t),s=0,r=i.length;s<r;++s)sC(o=t[i[s]])&&(su(o)?o.register(this,a):this.register(o));return--this.registerDepth,this}registerResolver(e,t){s$(e);let i=this.resolvers,o=i.get(e);return null==o?i.set(e,t):o instanceof sn&&4===o.strategy?o.state.push(t):i.set(e,new sn(e,4,[o,t])),t}registerTransformer(e,t){let i=this.getResolver(e);if(null==i)return!1;if(i.getFactory){let e=i.getFactory(this);return null!=e&&(e.registerTransformer(t),!0)}return!1}getResolver(e,t=!0){let i;if(s$(e),void 0!==e.resolve)return e;let o=this;for(;null!=o;){if(null!=(i=o.resolvers.get(e)))return i;if(null==o.parent){let i=sp(e)?this:o;return t?this.jitRegister(e,i):null}o=o.parent}return null}has(e,t=!1){return!!this.resolvers.has(e)||!!t&&null!=this.parent&&this.parent.has(e,!0)}get(e){let t;if(s$(e),e.$isResolver)return e.resolve(this,this);let i=this;for(;null!=i;){if(null!=(t=i.resolvers.get(e)))return t.resolve(i,this);if(null==i.parent){let o=sp(e)?this:i;return(t=this.jitRegister(e,o)).resolve(i,this)}i=i.parent}throw Error(`Unable to resolve key: ${String(e)}`)}getAll(e,t=!1){let i;s$(e);let o=this;if(t){let t=ee;for(;null!=o;)null!=(i=o.resolvers.get(e))&&(t=t.concat(sw(i,o,this))),o=o.parent;return t}for(;null!=o;){if(null!=(i=o.resolvers.get(e)))return sw(i,o,this);if(null==(o=o.parent))break}return ee}getFactory(e){let t=sb.get(e);if(void 0===t){if(sI(e))throw Error(`${e.name} is a native function and therefore cannot be safely constructed by DI. If this is intentional, please use a callback or cachedCallback resolver.`);sb.set(e,t=new sd(e,st.getDependencies(e)))}return t}registerFactory(e,t){sb.set(e,t)}createChild(e){return new sf(null,Object.assign({},this.config,e,{parentLocator:()=>this}))}jitRegister(e,t){if("function"!=typeof e)throw Error(`Attempted to jitRegister something that is not a constructor: '${e}'. Did you forget to register this dependency?`);if(sg.has(e.name))throw Error(`Attempted to jitRegister an intrinsic type: ${e.name}. Did you forget to add @inject(Key)`);if(su(e)){let i=e.register(t);if(!(i instanceof Object)||null==i.resolve){let i=t.resolvers.get(e);if(void 0!=i)return i;throw Error("A valid resolver was not returned from the static register method")}return i}if(e.$isInterface)throw Error(`Attempted to jitRegister an interface: ${e.friendlyName}`);{let i=this.config.defaultResolver(e,t);return t.resolvers.set(e,i),i}}}let sv=new WeakMap;function sy(e){return function(t,i,o){if(sv.has(o))return sv.get(o);let s=e(t,i,o);return sv.set(o,s),s}}let sx=Object.freeze({instance:(e,t)=>new sn(e,0,t),singleton:(e,t)=>new sn(e,1,t),transient:(e,t)=>new sn(e,2,t),callback:(e,t)=>new sn(e,3,t),cachedCallback:(e,t)=>new sn(e,3,sy(t)),aliasTo:(e,t)=>new sn(t,5,e)});function s$(e){if(null==e)throw Error("key/value cannot be null or undefined. Are you trying to inject/register something that doesn't exist with DI?")}function sw(e,t,i){if(e instanceof sn&&4===e.strategy){let o=e.state,s=o.length,r=Array(s);for(;s--;)r[s]=o[s].resolve(t,i);return r}return[e.resolve(t,i)]}let sk="(anonymous)";function sC(e){return"object"==typeof e&&null!==e||"function"==typeof e}let sI=(r=new WeakMap,a=!1,n="",l=0,function(e){return void 0===(a=r.get(e))&&(a=(l=(n=e.toString()).length)>=29&&l<=100&&125===n.charCodeAt(l-1)&&32>=n.charCodeAt(l-2)&&93===n.charCodeAt(l-3)&&101===n.charCodeAt(l-4)&&100===n.charCodeAt(l-5)&&111===n.charCodeAt(l-6)&&99===n.charCodeAt(l-7)&&32===n.charCodeAt(l-8)&&101===n.charCodeAt(l-9)&&118===n.charCodeAt(l-10)&&105===n.charCodeAt(l-11)&&116===n.charCodeAt(l-12)&&97===n.charCodeAt(l-13)&&110===n.charCodeAt(l-14)&&88===n.charCodeAt(l-15),r.set(e,a)),a}),sT={};function sF(e){switch(typeof e){case"number":return e>=0&&(0|e)===e;case"string":{let t=sT[e];if(void 0!==t)return t;let i=e.length;if(0===i)return sT[e]=!1;let o=0;for(let t=0;t<i;++t)if(o=e.charCodeAt(t),0===t&&48===o&&i>1||o<48||o>57)return sT[e]=!1;return sT[e]=!0}default:return!1}}function sS(e){return`${e.toLowerCase()}:presentation`}let sO=new Map,sD=Object.freeze({define(e,t,i){let o=sS(e);void 0===sO.get(o)?sO.set(o,t):sO.set(o,!1),i.register(sx.instance(o,t))},forTag(e,t){let i=sS(e),o=sO.get(i);return!1===o?st.findResponsibleContainer(t).get(i):o||null}});class sE{constructor(e,t){this.template=e||null,this.styles=void 0===t?null:Array.isArray(t)?eb.create(t):t instanceof eb?t:eb.create([t])}applyTo(e){let t=e.$fastController;null===t.template&&(t.template=this.template),null===t.styles&&(t.styles=this.styles)}}class sR extends ez{constructor(){super(...arguments),this._presentation=void 0}get $presentation(){return void 0===this._presentation&&(this._presentation=sD.forTag(this.tagName,this)),this._presentation}templateChanged(){void 0!==this.template&&(this.$fastController.template=this.template)}stylesChanged(){void 0!==this.styles&&(this.$fastController.styles=this.styles)}connectedCallback(){null!==this.$presentation&&this.$presentation.applyTo(this),super.connectedCallback()}static compose(e){return (t={})=>new sA(this===sR?class extends sR{}:this,e,t)}}function sL(e,t,i){return"function"==typeof e?e(t,i):e}X([eu],sR.prototype,"template",void 0),X([eu],sR.prototype,"styles",void 0);class sA{constructor(e,t,i){this.type=e,this.elementDefinition=t,this.overrideDefinition=i,this.definition=Object.assign(Object.assign({},this.elementDefinition),this.overrideDefinition)}register(e,t){let i=this.definition,o=this.overrideDefinition,s=i.prefix||t.elementPrefix,r=`${s}-${i.baseName}`;t.tryDefineElement({name:r,type:this.type,baseClass:this.elementDefinition.baseClass,callback:e=>{let t=new sE(sL(i.template,e,i),sL(i.styles,e,i));e.definePresentation(t);let s=sL(i.shadowOptions,e,i);e.shadowRootMode&&(s?o.shadowOptions||(s.mode=e.shadowRootMode):null!==s&&(s={mode:e.shadowRootMode})),e.defineElement({elementOptions:sL(i.elementOptions,e,i),shadowOptions:s,attributes:sL(i.attributes,e,i)})}})}}class sV{constructor(){this.targetIndex=0}}class sP extends sV{constructor(){super(...arguments),this.createPlaceholder=el.createInterpolationPlaceholder}}class sH extends sV{constructor(e,t,i){super(),this.name=e,this.behavior=t,this.options=i}createPlaceholder(e){return el.createCustomAttributePlaceholder(this.name,e)}createBehavior(e){return new this.behavior(e,this.options)}}function sz(e,t){this.source=e,this.context=t,null===this.bindingObserver&&(this.bindingObserver=ec.binding(this.binding,this,this.isBindingVolatile)),this.updateTarget(this.bindingObserver.observe(e,t))}function sM(e,t){this.source=e,this.context=t,this.target.addEventListener(this.targetName,this)}function sN(){this.bindingObserver.disconnect(),this.source=null,this.context=null}function sB(){this.bindingObserver.disconnect(),this.source=null,this.context=null;let e=this.target.$fastView;void 0!==e&&e.isComposed&&(e.unbind(),e.needsBindOnly=!0)}function sj(){this.target.removeEventListener(this.targetName,this),this.source=null,this.context=null}function sq(e){el.setAttribute(this.target,this.targetName,e)}function sU(e){el.setBooleanAttribute(this.target,this.targetName,e)}function s_(e){if(null==e&&(e=""),e.create){this.target.textContent="";let t=this.target.$fastView;void 0===t?t=e.create():this.target.$fastTemplate!==e&&(t.isComposed&&(t.remove(),t.unbind()),t=e.create()),t.isComposed?t.needsBindOnly&&(t.needsBindOnly=!1,t.bind(this.source,this.context)):(t.isComposed=!0,t.bind(this.source,this.context),t.insertBefore(this.target),this.target.$fastView=t,this.target.$fastTemplate=e)}else{let t=this.target.$fastView;void 0!==t&&t.isComposed&&(t.isComposed=!1,t.remove(),t.needsBindOnly?t.needsBindOnly=!1:t.unbind()),this.target.textContent=e}}function sG(e){this.target[this.targetName]=e}function sW(e){let t=this.classVersions||Object.create(null),i=this.target,o=this.version||0;if(null!=e&&e.length){let s=e.split(/\s+/);for(let e=0,r=s.length;e<r;++e){let r=s[e];""!==r&&(t[r]=o,i.classList.add(r))}}if(this.classVersions=t,this.version=o+1,0!==o)for(let e in o-=1,t)t[e]===o&&i.classList.remove(e)}class sK extends sP{constructor(e){super(),this.binding=e,this.bind=sz,this.unbind=sN,this.updateTarget=sq,this.isBindingVolatile=ec.isVolatileBinding(this.binding)}get targetName(){return this.originalTargetName}set targetName(e){if(this.originalTargetName=e,void 0!==e)switch(e[0]){case":":if(this.cleanedTargetName=e.substr(1),this.updateTarget=sG,"innerHTML"===this.cleanedTargetName){let e=this.binding;this.binding=(t,i)=>el.createHTML(e(t,i))}break;case"?":this.cleanedTargetName=e.substr(1),this.updateTarget=sU;break;case"@":this.cleanedTargetName=e.substr(1),this.bind=sM,this.unbind=sj;break;default:this.cleanedTargetName=e,"class"===e&&(this.updateTarget=sW)}}targetAtContent(){this.updateTarget=s_,this.unbind=sB}createBehavior(e){return new sX(e,this.binding,this.isBindingVolatile,this.bind,this.unbind,this.updateTarget,this.cleanedTargetName)}}class sX{constructor(e,t,i,o,s,r,a){this.source=null,this.context=null,this.bindingObserver=null,this.target=e,this.binding=t,this.isBindingVolatile=i,this.bind=o,this.unbind=s,this.updateTarget=r,this.targetName=a}handleChange(){this.updateTarget(this.bindingObserver.observe(this.source,this.context))}handleEvent(e){eg.setEvent(e);let t=this.binding(this.source,this.context);eg.setEvent(null),!0!==t&&e.preventDefault()}}let sY=null;class sQ{addFactory(e){e.targetIndex=this.targetIndex,this.behaviorFactories.push(e)}captureContentBinding(e){e.targetAtContent(),this.addFactory(e)}reset(){this.behaviorFactories=[],this.targetIndex=-1}release(){sY=this}static borrow(e){let t=sY||new sQ;return t.directives=e,t.reset(),sY=null,t}}let sZ=en.length;function sJ(e,t){let i=t.split(ea);if(1===i.length)return null;let o=[];for(let t=0,s=i.length;t<s;++t){let s,r=i[t],a=r.indexOf(en);if(-1===a)s=r;else{let t=parseInt(r.substring(0,a));o.push(e.directives[t]),s=r.substring(a+sZ)}""!==s&&o.push(s)}return o}function s0(e,t,i=!1){let o=t.attributes;for(let s=0,r=o.length;s<r;++s){let a=o[s],n=a.value,l=sJ(e,n),h=null;null===l?i&&((h=new sK(()=>n)).targetName=a.name):h=function(e){let t;if(1===e.length)return e[0];let i=e.length,o=e.map(e=>"string"==typeof e?()=>e:(t=e.targetName||t,e.binding)),s=new sK((e,t)=>{let s="";for(let r=0;r<i;++r)s+=o[r](e,t);return s});return s.targetName=t,s}(l),null!==h&&(t.removeAttributeNode(a),s--,r--,e.addFactory(h))}}let s1=document.createRange();class s2{constructor(e,t){this.fragment=e,this.behaviors=t,this.source=null,this.context=null,this.firstChild=e.firstChild,this.lastChild=e.lastChild}appendTo(e){e.appendChild(this.fragment)}insertBefore(e){if(this.fragment.hasChildNodes())e.parentNode.insertBefore(this.fragment,e);else{let t,i=this.lastChild;if(e.previousSibling===i)return;let o=e.parentNode,s=this.firstChild;for(;s!==i;)t=s.nextSibling,o.insertBefore(s,e),s=t;o.insertBefore(i,e)}}remove(){let e,t=this.fragment,i=this.lastChild,o=this.firstChild;for(;o!==i;)e=o.nextSibling,t.appendChild(o),o=e;t.appendChild(i)}dispose(){let e,t=this.firstChild.parentNode,i=this.lastChild,o=this.firstChild;for(;o!==i;)e=o.nextSibling,t.removeChild(o),o=e;t.removeChild(i);let s=this.behaviors,r=this.source;for(let e=0,t=s.length;e<t;++e)s[e].unbind(r)}bind(e,t){let i=this.behaviors;if(this.source!==e)if(null!==this.source){let o=this.source;this.source=e,this.context=t;for(let s=0,r=i.length;s<r;++s){let r=i[s];r.unbind(o),r.bind(e,t)}}else{this.source=e,this.context=t;for(let o=0,s=i.length;o<s;++o)i[o].bind(e,t)}}unbind(){if(null===this.source)return;let e=this.behaviors,t=this.source;for(let i=0,o=e.length;i<o;++i)e[i].unbind(t);this.source=null}static disposeContiguousBatch(e){if(0!==e.length){s1.setStartBefore(e[0].firstChild),s1.setEndAfter(e[e.length-1].lastChild),s1.deleteContents();for(let t=0,i=e.length;t<i;++t){let i=e[t],o=i.behaviors,s=i.source;for(let e=0,t=o.length;e<t;++e)o[e].unbind(s)}}}}class s5{constructor(e,t){this.behaviorCount=0,this.hasHostBehaviors=!1,this.fragment=null,this.targetOffset=0,this.viewBehaviorFactories=null,this.hostBehaviorFactories=null,this.html=e,this.directives=t}create(e){if(null===this.fragment){let e,t=this.html;if("string"==typeof t){(e=document.createElement("template")).innerHTML=el.createHTML(t);let i=e.content.firstElementChild;null!==i&&"TEMPLATE"===i.tagName&&(e=i)}else e=t;let i=function(e,t){let i,o=e.content;document.adoptNode(o);let s=sQ.borrow(t);s0(s,e,!0);let r=s.behaviorFactories;s.reset();let a=el.createTemplateWalker(o);for(;i=a.nextNode();)switch(s.targetIndex++,i.nodeType){case 1:s0(s,i);break;case 3:!function(e,t,i){let o=sJ(e,t.textContent);if(null!==o){let s=t;for(let r=0,a=o.length;r<a;++r){let a=o[r],n=0===r?t:s.parentNode.insertBefore(document.createTextNode(""),s.nextSibling);"string"==typeof a?n.textContent=a:(n.textContent=" ",e.captureContentBinding(a)),s=n,e.targetIndex++,n!==t&&i.nextNode()}e.targetIndex--}}(s,i,a);break;case 8:el.isMarker(i)&&s.addFactory(t[el.extractDirectiveIndexFromMarker(i)])}let n=0;(el.isMarker(o.firstChild)||1===o.childNodes.length&&t.length)&&(o.insertBefore(document.createComment(""),o.firstChild),n=-1);let l=s.behaviorFactories;return s.release(),{fragment:o,viewBehaviorFactories:l,hostBehaviorFactories:r,targetOffset:n}}(e,this.directives);this.fragment=i.fragment,this.viewBehaviorFactories=i.viewBehaviorFactories,this.hostBehaviorFactories=i.hostBehaviorFactories,this.targetOffset=i.targetOffset,this.behaviorCount=this.viewBehaviorFactories.length+this.hostBehaviorFactories.length,this.hasHostBehaviors=this.hostBehaviorFactories.length>0}let t=this.fragment.cloneNode(!0),i=this.viewBehaviorFactories,o=Array(this.behaviorCount),s=el.createTemplateWalker(t),r=0,a=this.targetOffset,n=s.nextNode();for(let e=i.length;r<e;++r){let e=i[r],t=e.targetIndex;for(;null!==n;)if(a===t){o[r]=e.createBehavior(n);break}else n=s.nextNode(),a++}if(this.hasHostBehaviors){let t=this.hostBehaviorFactories;for(let i=0,s=t.length;i<s;++i,++r)o[r]=t[i].createBehavior(e)}return new s2(t,o)}render(e,t,i){"string"==typeof t&&(t=document.getElementById(t)),void 0===i&&(i=t);let o=this.create(i);return o.bind(e,em),o.appendTo(t),o}}let s3=/([ \x09\x0a\x0c\x0d])([^\0-\x1F\x7F-\x9F "'>=/]+)([ \x09\x0a\x0c\x0d]*=[ \x09\x0a\x0c\x0d]*(?:[^ \x09\x0a\x0c\x0d"'`<>=]*|"[^"]*|'[^']*))$/;function s4(e,...t){let i=[],o="";for(let s=0,r=e.length-1;s<r;++s){let r=e[s],a=t[s];if(o+=r,a instanceof s5){let e=a;a=()=>e}if("function"==typeof a&&(a=new sK(a)),a instanceof sP){let e=s3.exec(r);null!==e&&(a.targetName=e[2])}a instanceof sV?(o+=a.createPlaceholder(i.length),i.push(a)):o+=a}return new s5(o+=e[e.length-1],i)}class s6{constructor(e,t){this.target=e,this.propertyName=t}bind(e){e[this.propertyName]=this.target}unbind(){}}function s9(e){return new sH("fast-ref",s6,e)}class s8{handleStartContentChange(){this.startContainer.classList.toggle("start",this.start.assignedNodes().length>0)}handleEndContentChange(){this.endContainer.classList.toggle("end",this.end.assignedNodes().length>0)}}let s7=(e,t)=>s4`
    <span
        part="end"
        ${s9("endContainer")}
        class=${e=>t.end?"end":void 0}
    >
        <slot name="end" ${s9("end")} @slotchange="${e=>e.handleEndContentChange()}">
            ${t.end||""}
        </slot>
    </span>
`,re=(e,t)=>s4`
    <span
        part="start"
        ${s9("startContainer")}
        class="${e=>t.start?"start":void 0}"
    >
        <slot
            name="start"
            ${s9("start")}
            @slotchange="${e=>e.handleStartContentChange()}"
        >
            ${t.start||""}
        </slot>
    </span>
`;function rt(e,...t){let i=eC.locate(e);t.forEach(t=>{Object.getOwnPropertyNames(t.prototype).forEach(i=>{"constructor"!==i&&Object.defineProperty(e.prototype,i,Object.getOwnPropertyDescriptor(t.prototype,i))}),eC.locate(t).forEach(e=>i.push(e))})}s4`
    <span part="end" ${s9("endContainer")}>
        <slot
            name="end"
            ${s9("end")}
            @slotchange="${e=>e.handleEndContentChange()}"
        ></slot>
    </span>
`,s4`
    <span part="start" ${s9("startContainer")}>
        <slot
            name="start"
            ${s9("start")}
            @slotchange="${e=>e.handleStartContentChange()}"
        ></slot>
    </span>
`;class ri extends sR{constructor(){super(...arguments),this.headinglevel=2,this.expanded=!1,this.clickHandler=e=>{this.expanded=!this.expanded,this.change()},this.change=()=>{this.$emit("change")}}}X([eS({attribute:"heading-level",mode:"fromView",converter:eT})],ri.prototype,"headinglevel",void 0),X([eS({mode:"boolean"})],ri.prototype,"expanded",void 0),X([eS],ri.prototype,"id",void 0),rt(ri,s8);class ro extends sR{constructor(){super(...arguments),this.expandmode="multi",this.activeItemIndex=0,this.change=()=>{this.$emit("change",this.activeid)},this.setItems=()=>{var e;0!==this.accordionItems.length&&(this.accordionIds=this.getItemIds(),this.accordionItems.forEach((e,t)=>{e instanceof ri&&(e.addEventListener("change",this.activeItemChange),this.isSingleExpandMode()&&(this.activeItemIndex!==t?e.expanded=!1:e.expanded=!0));let i=this.accordionIds[t];e.setAttribute("id","string"!=typeof i?`accordion-${t+1}`:i),this.activeid=this.accordionIds[this.activeItemIndex],e.addEventListener("keydown",this.handleItemKeyDown),e.addEventListener("focus",this.handleItemFocus)}),this.isSingleExpandMode()&&(null!=(e=this.findExpandedItem())?e:this.accordionItems[0]).setAttribute("aria-disabled","true"))},this.removeItemListeners=e=>{e.forEach((e,t)=>{e.removeEventListener("change",this.activeItemChange),e.removeEventListener("keydown",this.handleItemKeyDown),e.removeEventListener("focus",this.handleItemFocus)})},this.activeItemChange=e=>{if(e.defaultPrevented||e.target!==e.currentTarget)return;e.preventDefault();let t=e.target;this.activeid=t.getAttribute("id"),this.isSingleExpandMode()&&(this.resetItems(),t.expanded=!0,t.setAttribute("aria-disabled","true"),this.accordionItems.forEach(e=>{e.hasAttribute("disabled")||e.id===this.activeid||e.removeAttribute("aria-disabled")})),this.activeItemIndex=Array.from(this.accordionItems).indexOf(t),this.change()},this.handleItemKeyDown=e=>{if(e.target===e.currentTarget)switch(this.accordionIds=this.getItemIds(),e.key){case o0:e.preventDefault(),this.adjust(-1);break;case oQ:e.preventDefault(),this.adjust(1);break;case"Home":this.activeItemIndex=0,this.focusItem();break;case"End":this.activeItemIndex=this.accordionItems.length-1,this.focusItem()}},this.handleItemFocus=e=>{if(e.target===e.currentTarget){let t=e.target,i=this.activeItemIndex=Array.from(this.accordionItems).indexOf(t);this.activeItemIndex!==i&&-1!==i&&(this.activeItemIndex=i,this.activeid=this.accordionIds[this.activeItemIndex])}}}accordionItemsChanged(e,t){this.$fastController.isConnected&&(this.removeItemListeners(e),this.setItems())}findExpandedItem(){for(let e=0;e<this.accordionItems.length;e++)if("true"===this.accordionItems[e].getAttribute("expanded"))return this.accordionItems[e];return null}resetItems(){this.accordionItems.forEach((e,t)=>{e.expanded=!1})}getItemIds(){return this.accordionItems.map(e=>e.getAttribute("id"))}isSingleExpandMode(){return"single"===this.expandmode}adjust(e){var t,i;this.activeItemIndex=(t=this.accordionItems.length-1,(i=this.activeItemIndex+e)<0?t:i>t?0:i),this.focusItem()}focusItem(){let e=this.accordionItems[this.activeItemIndex];e instanceof ri&&e.expandbutton.focus()}}function rs(e){return e?function(t,i,o){return 1===t.nodeType&&t.matches(e)}:function(e,t,i){return 1===e.nodeType}}X([eS({attribute:"expand-mode"})],ro.prototype,"expandmode",void 0),X([eu],ro.prototype,"accordionItems",void 0);class rr{constructor(e,t){this.target=e,this.options=t,this.source=null}bind(e){let t=this.options.property;this.shouldUpdate=ec.getAccessors(e).some(e=>e.name===t),this.source=e,this.updateTarget(this.computeNodes()),this.shouldUpdate&&this.observe()}unbind(){this.updateTarget(ee),this.source=null,this.shouldUpdate&&this.disconnect()}handleEvent(){this.updateTarget(this.computeNodes())}computeNodes(){let e=this.getNodes();return void 0!==this.options.filter&&(e=e.filter(this.options.filter)),e}updateTarget(e){this.source[this.options.property]=e}}class ra extends rr{constructor(e,t){super(e,t)}observe(){this.target.addEventListener("slotchange",this)}disconnect(){this.target.removeEventListener("slotchange",this)}getNodes(){return this.target.assignedNodes(this.options)}}function rn(e){return"string"==typeof e&&(e={property:e}),new sH("fast-slotted",ra,e)}function rl(e,t){let i=[],o="",s=[];for(let r=0,a=e.length-1;r<a;++r){o+=e[r];let a=t[r];if(a instanceof Y){let e=a.createBehavior();a=a.createCSS(),e&&s.push(e)}a instanceof eb||a instanceof CSSStyleSheet?(""!==o.trim()&&(i.push(o),o=""),i.push(a)):o+=a}return""!==(o+=e[e.length-1]).trim()&&i.push(o),{styles:i,behaviors:s}}function rh(e,...t){let{styles:i,behaviors:o}=rl(e,t),s=eb.create(i);return o.length&&s.withBehaviors(...o),s}class rd extends Y{constructor(e,t){super(),this.behaviors=t,this.css="";const i=e.reduce((e,t)=>("string"==typeof t?this.css+=t:e.push(t),e),[]);i.length&&(this.styles=eb.create(i))}createBehavior(){return this}createCSS(){return this.css}bind(e){this.styles&&e.$fastController.addStyles(this.styles),this.behaviors.length&&e.$fastController.addBehaviors(this.behaviors)}unbind(e){this.styles&&e.$fastController.removeStyles(this.styles),this.behaviors.length&&e.$fastController.removeBehaviors(this.behaviors)}}function rc(e,...t){let{styles:i,behaviors:o}=rl(e,t);return new rd(i,o)}function ru(e){return`:host([hidden]){display:none}:host{display:${e}}`}let rp=(e,t)=>rh`
  ${ru("flex")} :host {
    box-sizing: border-box;
    flex-direction: column;
    font-family: ${tu};
    font-size: ${tT};
    line-height: ${tF};
    color: ${op};
    border-top: calc(${tw} * 1px) solid ${ox};
  }
`;function rg(...e){return e.every(e=>e instanceof HTMLElement)}let rm=!function(){let e;if("boolean"==typeof o)return o;if(!("u">typeof window&&window.document&&window.document.createElement))return o=!1;let t=document.createElement("style"),i=(e=document.querySelector('meta[property="csp-nonce"]'))?e.getAttribute("content"):null;null!==i&&t.setAttribute("nonce",i),document.head.appendChild(t);try{t.sheet.insertRule("foo:focus-visible {color:inherit}",0),o=!0}catch(e){o=!1}finally{document.head.removeChild(t)}return o}()?"focus":"focus-visible";class rb{constructor(e){this.listenerCache=new WeakMap,this.query=e}bind(e){let{query:t}=this,i=this.constructListener(e);i.bind(t)(),t.addListener(i),this.listenerCache.set(e,i)}unbind(e){let t=this.listenerCache.get(e);t&&(this.query.removeListener(t),this.listenerCache.delete(e))}}class rf extends rb{constructor(e,t){super(e),this.styles=t}static with(e){return t=>new rf(e,t)}constructListener(e){let t=!1,i=this.styles;return function(){let{matches:o}=this;o&&!t?(e.$fastController.addStyles(i),t=o):!o&&t&&(e.$fastController.removeStyles(i),t=o)}}unbind(e){super.unbind(e),e.$fastController.removeStyles(this.styles)}}let rv=rf.with(window.matchMedia("(forced-colors)"));rf.with(window.matchMedia("(prefers-color-scheme: dark)")),rf.with(window.matchMedia("(prefers-color-scheme: light)")),(n1=nX||(nX={})).Canvas="Canvas",n1.CanvasText="CanvasText",n1.LinkText="LinkText",n1.VisitedText="VisitedText",n1.ActiveText="ActiveText",n1.ButtonFace="ButtonFace",n1.ButtonText="ButtonText",n1.Field="Field",n1.FieldText="FieldText",n1.Highlight="Highlight",n1.HighlightText="HighlightText",n1.GrayText="GrayText";let ry=rc`(${tp} + ${tf} + ${ty}) * ${tv}`,rx=(e,t)=>rh`
    ${ru("flex")} :host {
      box-sizing: border-box;
      font-family: ${tu};
      flex-direction: column;
      font-size: ${tT};
      line-height: ${tF};
      border-bottom: calc(${tw} * 1px) solid
        ${ox};
    }

    .region {
      display: none;
      padding: calc((6 + (${tv} * 2 * ${tf})) * 1px);
    }

    div.heading {
      display: grid;
      position: relative;
      grid-template-columns: calc(${ry} * 1px) auto 1fr auto;
      color: ${op};
    }

    .button {
      appearance: none;
      border: none;
      background: none;
      grid-column: 3;
      outline: none;
      padding: 0 calc((6 + (${tv} * 2 * ${tf})) * 1px);
      text-align: left;
      height: calc(${ry} * 1px);
      color: currentcolor;
      cursor: pointer;
      font-family: inherit;
    }

    .button:hover {
      color: currentcolor;
    }

    .button:active {
      color: currentcolor;
    }

    .button::before {
      content: '';
      position: absolute;
      top: 0;
      left: 0;
      right: 0;
      bottom: 0;
      cursor: pointer;
    }

    /* prettier-ignore */
    .button:${rm}::before {
      outline: none;
      border: calc(${tk} * 1px) solid ${iE};
      border-radius: calc(${tb} * 1px);
    }

    :host([expanded]) .region {
      display: block;
    }

    .icon {
      display: flex;
      align-items: center;
      justify-content: center;
      grid-column: 1;
      grid-row: 1;
      pointer-events: none;
      position: relative;
    }

    slot[name='expanded-icon'],
    slot[name='collapsed-icon'] {
      fill: currentcolor;
    }

    slot[name='collapsed-icon'] {
      display: flex;
    }

    :host([expanded]) slot[name='collapsed-icon'] {
      display: none;
    }

    slot[name='expanded-icon'] {
      display: none;
    }

    :host([expanded]) slot[name='expanded-icon'] {
      display: flex;
    }

    .start {
      display: flex;
      align-items: center;
      padding-inline-start: calc(${tv} * 1px);
      justify-content: center;
      grid-column: 2;
      position: relative;
    }

    .end {
      display: flex;
      align-items: center;
      justify-content: center;
      grid-column: 4;
      position: relative;
    }
  `.withBehaviors(rv(rh`
      /* prettier-ignore */
      .button:${rm}::before {
          border-color: ${nX.Highlight};
        }
      :host slot[name='collapsed-icon'],
      :host([expanded]) slot[name='expanded-icon'] {
        fill: ${nX.ButtonText};
      }
    `));class r$ extends ri{}let rw=r$.compose({baseName:"accordion-item",baseClass:ri,template:(e,t)=>s4`
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
                ${s9("expandbutton")}
                aria-expanded="${e=>e.expanded}"
                aria-controls="${e=>e.id}-panel"
                id="${e=>e.id}"
                @click="${(e,t)=>e.clickHandler(t.event)}"
            >
                <span class="heading-content" part="heading-content">
                    <slot name="heading"></slot>
                </span>
            </button>
            ${re(e,t)}
            ${s7(e,t)}
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
`,styles:rx,collapsedIcon:`
        <svg
            width="20"
            height="20"
            viewBox="0 0 16 16"
            xmlns="http://www.w3.org/2000/svg"
        >
            <path
                fill-rule="evenodd"
                clip-rule="evenodd"
                d="M5.00001 12.3263C5.00124 12.5147 5.05566 12.699 5.15699 12.8578C5.25831 13.0167 5.40243 13.1437 5.57273 13.2242C5.74304 13.3047 5.9326 13.3354 6.11959 13.3128C6.30659 13.2902 6.4834 13.2152 6.62967 13.0965L10.8988 8.83532C11.0739 8.69473 11.2153 8.51658 11.3124 8.31402C11.4096 8.11146 11.46 7.88966 11.46 7.66499C11.46 7.44033 11.4096 7.21853 11.3124 7.01597C11.2153 6.81341 11.0739 6.63526 10.8988 6.49467L6.62967 2.22347C6.48274 2.10422 6.30501 2.02912 6.11712 2.00691C5.92923 1.9847 5.73889 2.01628 5.56823 2.09799C5.39757 2.17969 5.25358 2.30817 5.153 2.46849C5.05241 2.62882 4.99936 2.8144 5.00001 3.00369V12.3263Z"
            />
        </svg>
    `,expandedIcon:`
        <svg
            width="20"
            height="20"
            viewBox="0 0 16 16"
            xmlns="http://www.w3.org/2000/svg"
        >
            <path
                fill-rule="evenodd"
                clip-rule="evenodd"
                transform="rotate(90,8,8)"
          d="M5.00001 12.3263C5.00124 12.5147 5.05566 12.699 5.15699 12.8578C5.25831 13.0167 5.40243 13.1437 5.57273 13.2242C5.74304 13.3047 5.9326 13.3354 6.11959 13.3128C6.30659 13.2902 6.4834 13.2152 6.62967 13.0965L10.8988 8.83532C11.0739 8.69473 11.2153 8.51658 11.3124 8.31402C11.4096 8.11146 11.46 7.88966 11.46 7.66499C11.46 7.44033 11.4096 7.21853 11.3124 7.01597C11.2153 6.81341 11.0739 6.63526 10.8988 6.49467L6.62967 2.22347C6.48274 2.10422 6.30501 2.02912 6.11712 2.00691C5.92923 1.9847 5.73889 2.01628 5.56823 2.09799C5.39757 2.17969 5.25358 2.30817 5.153 2.46849C5.05241 2.62882 4.99936 2.8144 5.00001 3.00369V12.3263Z"
            />
        </svg>
    `});class rk extends ro{}let rC=rk.compose({baseName:"accordion",baseClass:ro,template:(e,t)=>s4`
    <template>
        <slot ${rn({property:"accordionItems",filter:rs()})}></slot>
        <slot name="item" part="item" ${rn("accordionItems")}></slot>
    </template>
`,styles:rp});function rI(e,t,i,o){var s,r=arguments.length,a=r<3?t:null===o?o=Object.getOwnPropertyDescriptor(t,i):o;if("object"==typeof Reflect&&"function"==typeof Reflect.decorate)a=Reflect.decorate(e,t,i,o);else for(var n=e.length-1;n>=0;n--)(s=e[n])&&(a=(r<3?s(a):r>3?s(t,i,a):s(t,i))||a);return r>3&&a&&Object.defineProperty(t,i,a),a}class rT{}X([eS({attribute:"aria-atomic"})],rT.prototype,"ariaAtomic",void 0),X([eS({attribute:"aria-busy"})],rT.prototype,"ariaBusy",void 0),X([eS({attribute:"aria-controls"})],rT.prototype,"ariaControls",void 0),X([eS({attribute:"aria-current"})],rT.prototype,"ariaCurrent",void 0),X([eS({attribute:"aria-describedby"})],rT.prototype,"ariaDescribedby",void 0),X([eS({attribute:"aria-details"})],rT.prototype,"ariaDetails",void 0),X([eS({attribute:"aria-disabled"})],rT.prototype,"ariaDisabled",void 0),X([eS({attribute:"aria-errormessage"})],rT.prototype,"ariaErrormessage",void 0),X([eS({attribute:"aria-flowto"})],rT.prototype,"ariaFlowto",void 0),X([eS({attribute:"aria-haspopup"})],rT.prototype,"ariaHaspopup",void 0),X([eS({attribute:"aria-hidden"})],rT.prototype,"ariaHidden",void 0),X([eS({attribute:"aria-invalid"})],rT.prototype,"ariaInvalid",void 0),X([eS({attribute:"aria-keyshortcuts"})],rT.prototype,"ariaKeyshortcuts",void 0),X([eS({attribute:"aria-label"})],rT.prototype,"ariaLabel",void 0),X([eS({attribute:"aria-labelledby"})],rT.prototype,"ariaLabelledby",void 0),X([eS({attribute:"aria-live"})],rT.prototype,"ariaLive",void 0),X([eS({attribute:"aria-owns"})],rT.prototype,"ariaOwns",void 0),X([eS({attribute:"aria-relevant"})],rT.prototype,"ariaRelevant",void 0),X([eS({attribute:"aria-roledescription"})],rT.prototype,"ariaRoledescription",void 0);class rF extends sR{constructor(){super(...arguments),this.handleUnsupportedDelegatesFocus=()=>{var e;window.ShadowRoot&&!window.ShadowRoot.prototype.hasOwnProperty("delegatesFocus")&&(null==(e=this.$fastController.definition.shadowOptions)?void 0:e.delegatesFocus)&&(this.focus=()=>{var e;null==(e=this.control)||e.focus()})}}connectedCallback(){super.connectedCallback(),this.handleUnsupportedDelegatesFocus()}}X([eS],rF.prototype,"download",void 0),X([eS],rF.prototype,"href",void 0),X([eS],rF.prototype,"hreflang",void 0),X([eS],rF.prototype,"ping",void 0),X([eS],rF.prototype,"referrerpolicy",void 0),X([eS],rF.prototype,"rel",void 0),X([eS],rF.prototype,"target",void 0),X([eS],rF.prototype,"type",void 0),X([eu],rF.prototype,"defaultSlottedContent",void 0);class rS{}X([eS({attribute:"aria-expanded"})],rS.prototype,"ariaExpanded",void 0),rt(rS,rT),rt(rF,s8,rS);let rO=(e,t)=>s4`
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
        ${s9("control")}
    >
        ${re(e,t)}
        <span class="content" part="content">
            <slot ${rn("defaultSlottedContent")}></slot>
        </span>
        ${s7(e,t)}
    </a>
`,rD=rh`
  ${ru("inline-flex")} :host {
    font-family: ${tu};
    outline: none;
    font-size: ${tC};
    line-height: ${tI};
    height: calc(${ry} * 1px);
    min-width: calc(${ry} * 1px);
    background-color: ${iX};
    color: ${op};
    border-radius: calc(${tb} * 1px);
    fill: currentcolor;
    cursor: pointer;
    margin: calc((${tk} + 2) * 1px);
  }

  .control {
    background: transparent;
    height: inherit;
    flex-grow: 1;
    box-sizing: border-box;
    display: inline-flex;
    justify-content: center;
    align-items: center;
    padding: 0
      max(
        1px,
        calc((10 + (${tv} * 2 * (${tf} + ${ty})))) * 1px
      );
    white-space: nowrap;
    outline: none;
    text-decoration: none;
    border: calc(${tw} * 1px) solid transparent;
    color: inherit;
    border-radius: inherit;
    fill: inherit;
    cursor: inherit;
    font-family: inherit;
    font-size: inherit;
    line-height: inherit;
  }

  :host(:hover) {
    background-color: ${iY};
  }

  :host(:active) {
    background-color: ${iQ};
  }

  :host([aria-pressed='true']) {
    box-shadow: inset 0px 0px 2px 2px ${oi};
  }

  :host([minimal]),
  :host([scale='xsmall']) {
    --element-scale: -4;
  }

  :host([scale='small']) {
    --element-scale: -2;
  }

  :host([scale='medium']) {
    --element-scale: 0;
  }

  :host([scale='large']) {
    --element-scale: 2;
  }

  :host([scale='xlarge']) {
    --element-scale: 4;
  }

  /* prettier-ignore */
  .control:${rm} {
      outline: calc(${tk} * 1px) solid ${oo};
      outline-offset: 2px;
      -moz-outline-radius: 0px;
    }

  .control::-moz-focus-inner {
    border: 0;
  }

  .start,
  .end {
    display: flex;
  }

  .control.icon-only {
    padding: 0;
    line-height: 0;
  }

  ::slotted(svg) {
    ${""} width: 16px;
    height: 16px;
    pointer-events: none;
  }

  .start {
    margin-inline-end: 11px;
  }

  .end {
    margin-inline-start: 11px;
  }
`.withBehaviors(rv(rh`
    :host .control {
      background-color: ${nX.ButtonFace};
      border-color: ${nX.ButtonText};
      color: ${nX.ButtonText};
      fill: currentColor;
    }

    :host(:hover) .control {
      forced-color-adjust: none;
      background-color: ${nX.Highlight};
      color: ${nX.HighlightText};
    }

    /* prettier-ignore */
    .control:${rm} {
          forced-color-adjust: none;
          background-color: ${nX.Highlight};
          outline-color: ${nX.ButtonText};
          color: ${nX.HighlightText};
        }

    .control:hover,
    :host([appearance='outline']) .control:hover {
      border-color: ${nX.ButtonText};
    }

    :host([href]) .control {
      border-color: ${nX.LinkText};
      color: ${nX.LinkText};
    }

    :host([href]) .control:hover,
        :host([href]) .control:${rm} {
      forced-color-adjust: none;
      background: ${nX.ButtonFace};
      outline-color: ${nX.LinkText};
      color: ${nX.LinkText};
      fill: currentColor;
    }
  `)),rE=rh`
  :host([appearance='accent']) {
    background: ${iS};
    color: ${iA};
  }

  :host([appearance='accent']:hover) {
    background: ${iO};
    color: ${iV};
  }

  :host([appearance='accent'][aria-pressed='true']) {
    box-shadow: inset 0px 0px 2px 2px ${iG};
  }

  :host([appearance='accent']:active) .control:active {
    background: ${iD};
    color: ${iP};
  }

  :host([appearance="accent"]) .control:${rm} {
    outline-color: ${iE};
  }
`.withBehaviors(rv(rh`
    :host([appearance='accent']) .control {
      forced-color-adjust: none;
      background: ${nX.Highlight};
      color: ${nX.HighlightText};
    }

    :host([appearance='accent']) .control:hover,
    :host([appearance='accent']:active) .control:active {
      background: ${nX.HighlightText};
      border-color: ${nX.Highlight};
      color: ${nX.Highlight};
    }

    :host([appearance="accent"]) .control:${rm} {
      outline-color: ${nX.Highlight};
    }

    :host([appearance='accent'][href]) .control {
      background: ${nX.LinkText};
      color: ${nX.HighlightText};
    }

    :host([appearance='accent'][href]) .control:hover {
      background: ${nX.ButtonFace};
      border-color: ${nX.LinkText};
      box-shadow: none;
      color: ${nX.LinkText};
      fill: currentColor;
    }

    :host([appearance="accent"][href]) .control:${rm} {
      outline-color: ${nX.HighlightText};
    }
  `)),rR=rh`
  :host([appearance='error']) {
    background: ${oI};
    color: ${iA};
  }

  :host([appearance='error']:hover) {
    background: ${oT};
    color: ${iV};
  }

  :host([appearance='error'][aria-pressed='true']) {
    box-shadow: inset 0px 0px 2px 2px ${oq};
  }

  :host([appearance='error']:active) .control:active {
    background: ${oF};
    color: ${iP};
  }

  :host([appearance="error"]) .control:${rm} {
    outline-color: ${oS};
  }
`.withBehaviors(rv(rh`
    :host([appearance='error']) .control {
      forced-color-adjust: none;
      background: ${nX.Highlight};
      color: ${nX.HighlightText};
    }

    :host([appearance='error']) .control:hover,
    :host([appearance='error']:active) .control:active {
      background: ${nX.HighlightText};
      border-color: ${nX.Highlight};
      color: ${nX.Highlight};
    }

    :host([appearance="error"]) .control:${rm} {
      outline-color: ${nX.Highlight};
    }

    :host([appearance='error'][href]) .control {
      background: ${nX.LinkText};
      color: ${nX.HighlightText};
    }

    :host([appearance='error'][href]) .control:hover {
      background: ${nX.ButtonFace};
      border-color: ${nX.LinkText};
      box-shadow: none;
      color: ${nX.LinkText};
      fill: currentColor;
    }

    :host([appearance="error"][href]) .control:${rm} {
      outline-color: ${nX.HighlightText};
    }
  `)),rL=rh`
  :host([appearance='hypertext']) {
    font-size: inherit;
    line-height: inherit;
    height: auto;
    min-width: 0;
    background: transparent;
  }

  :host([appearance='hypertext']) .control {
    display: inline;
    padding: 0;
    border: none;
    box-shadow: none;
    border-radius: 0;
    line-height: 1;
  }

  :host a.control:not(:link) {
    background-color: transparent;
    cursor: default;
  }
  :host([appearance='hypertext']) .control:link,
  :host([appearance='hypertext']) .control:visited {
    background: transparent;
    color: ${iU};
    border-bottom: calc(${tw} * 1px) solid ${iU};
  }

  :host([appearance='hypertext']:hover),
  :host([appearance='hypertext']) .control:hover {
    background: transparent;
    border-bottom-color: ${i_};
  }

  :host([appearance='hypertext']:active),
  :host([appearance='hypertext']) .control:active {
    background: transparent;
    border-bottom-color: ${iG};
  }

  :host([appearance="hypertext"]) .control:${rm} {
    outline-color: transparent;
    border-bottom: calc(${tk} * 1px) solid ${on};
    margin-bottom: calc(calc(${tw} - ${tk}) * 1px);
  }
`.withBehaviors(rv(rh`
    :host([appearance='hypertext']:hover) {
      background-color: ${nX.ButtonFace};
      color: ${nX.ButtonText};
    }
    :host([appearance="hypertext"][href]) .control:hover,
        :host([appearance="hypertext"][href]) .control:active,
        :host([appearance="hypertext"][href]) .control:${rm} {
      color: ${nX.LinkText};
      border-bottom-color: ${nX.LinkText};
      box-shadow: none;
    }
  `)),rA=rh`
  :host([appearance='lightweight']) {
    background: transparent;
    color: ${iU};
  }

  :host([appearance='lightweight']) .control {
    padding: 0;
    height: initial;
    border: none;
    box-shadow: none;
    border-radius: 0;
  }

  :host([appearance='lightweight']:hover) {
    background: transparent;
    color: ${i_};
  }

  :host([appearance='lightweight']:active) {
    background: transparent;
    color: ${iG};
  }

  :host([appearance='lightweight']) .content {
    position: relative;
  }

  :host([appearance='lightweight']) .content::before {
    content: '';
    display: block;
    height: calc(${tw} * 1px);
    position: absolute;
    top: calc(1em + 4px);
    width: 100%;
  }

  :host([appearance='lightweight']:hover) .content::before {
    background: ${i_};
  }

  :host([appearance='lightweight']:active) .content::before {
    background: ${iG};
  }

  :host([appearance="lightweight"]) .control:${rm} {
    outline-color: transparent;
  }

  :host([appearance="lightweight"]) .control:${rm} .content::before {
    background: ${op};
    height: calc(${tk} * 1px);
  }
`.withBehaviors(rv(rh`
    :host([appearance="lightweight"]) .control:hover,
        :host([appearance="lightweight"]) .control:${rm} {
      forced-color-adjust: none;
      background: ${nX.ButtonFace};
      color: ${nX.Highlight};
    }
    :host([appearance="lightweight"]) .control:hover .content::before,
        :host([appearance="lightweight"]) .control:${rm} .content::before {
      background: ${nX.Highlight};
    }

    :host([appearance="lightweight"][href]) .control:hover,
        :host([appearance="lightweight"][href]) .control:${rm} {
      background: ${nX.ButtonFace};
      box-shadow: none;
      color: ${nX.LinkText};
    }

    :host([appearance="lightweight"][href]) .control:hover .content::before,
        :host([appearance="lightweight"][href]) .control:${rm} .content::before {
      background: ${nX.LinkText};
    }
  `)),rV=rh`
  :host([appearance='outline']) {
    background: transparent;
    border-color: ${iS};
  }

  :host([appearance='outline']:hover) {
    border-color: ${iO};
  }

  :host([appearance='outline']:active) {
    border-color: ${iD};
  }

  :host([appearance='outline']) .control {
    border-color: inherit;
  }

  :host([appearance="outline"]) .control:${rm} {
    outline-color: ${iE};
  }
`.withBehaviors(rv(rh`
    :host([appearance='outline']) .control {
      border-color: ${nX.ButtonText};
    }
    :host([appearance="outline"]) .control:${rm} {
      forced-color-adjust: none;
      background-color: ${nX.Highlight};
      outline-color: ${nX.ButtonText};
      color: ${nX.HighlightText};
      fill: currentColor;
    }
    :host([appearance='outline'][href]) .control {
      background: ${nX.ButtonFace};
      border-color: ${nX.LinkText};
      color: ${nX.LinkText};
      fill: currentColor;
    }
    :host([appearance="outline"][href]) .control:hover,
        :host([appearance="outline"][href]) .control:${rm} {
      forced-color-adjust: none;
      outline-color: ${nX.LinkText};
    }
  `)),rP=rh`
  :host([appearance='stealth']),
  :host([appearance='stealth'][disabled]:active),
  :host([appearance='stealth'][disabled]:hover) {
    background: transparent;
  }

  :host([appearance='stealth']:hover) {
    background: ${i6};
  }

  :host([appearance='stealth']:active) {
    background: ${i9};
  }

  :host([appearance='stealth']) .control:${rm} {
    outline-color: ${iE};
  }

  /* Make the focus outline displayed within the button if
     it is in a start or end slot; e.g. in a tree item
     This will make the focus outline bounded within the container.
   */
  :host([appearance='stealth'][slot="end"]) .control:${rm},
  :host([appearance='stealth'][slot="start"]) .control:${rm} {
    outline-offset: -2px;
  }
`.withBehaviors(rv(rh`
    :host([appearance='stealth']),
    :host([appearance='stealth']) .control {
      forced-color-adjust: none;
      background: ${nX.ButtonFace};
      border-color: transparent;
      color: ${nX.ButtonText};
      fill: currentColor;
    }

    :host([appearance='stealth']:hover) .control {
      background: ${nX.Highlight};
      border-color: ${nX.Highlight};
      color: ${nX.HighlightText};
      fill: currentColor;
    }

    :host([appearance="stealth"]:${rm}) .control {
      outline-color: ${nX.Highlight};
      color: ${nX.HighlightText};
      fill: currentColor;
    }

    :host([appearance='stealth'][href]) .control {
      color: ${nX.LinkText};
    }

    :host([appearance="stealth"][href]:hover) .control,
        :host([appearance="stealth"][href]:${rm}) .control {
      background: ${nX.LinkText};
      border-color: ${nX.LinkText};
      color: ${nX.HighlightText};
      fill: currentColor;
    }

    :host([appearance="stealth"][href]:${rm}) .control {
      forced-color-adjust: none;
      box-shadow: 0 0 0 1px ${nX.LinkText};
    }
  `));class rH{constructor(e,t,i){this.propertyName=e,this.value=t,this.styles=i}bind(e){ec.getNotifier(e).subscribe(this,this.propertyName),this.handleChange(e,this.propertyName)}unbind(e){ec.getNotifier(e).unsubscribe(this,this.propertyName),e.$fastController.removeStyles(this.styles)}handleChange(e,t){e[t]===this.value?e.$fastController.addStyles(this.styles):e.$fastController.removeStyles(this.styles)}}function rz(e,t){return new rH("appearance",e,t)}let rM=(e,t)=>rh`
    ${rD}
  `.withBehaviors(rz("accent",rE),rz("hypertext",rL),rz("lightweight",rA),rz("outline",rV),rz("stealth",rP));class rN extends rF{appearanceChanged(e,t){this.$fastController.isConnected&&(this.classList.remove(e),this.classList.add(t))}connectedCallback(){super.connectedCallback(),this.appearance||(this.appearance="neutral")}defaultSlottedContentChanged(e,t){let i=this.defaultSlottedContent.filter(e=>e.nodeType===Node.ELEMENT_NODE);1===i.length&&i[0]instanceof SVGElement?this.control.classList.add("icon-only"):this.control.classList.remove("icon-only")}}rI([eS],rN.prototype,"appearance",void 0);let rB=rN.compose({baseName:"anchor",baseClass:rF,template:rO,styles:rM,shadowOptions:{delegatesFocus:!0}}),rj=e=>{let t=e.closest("[dir]");return null!==t&&"rtl"===t.dir?nU.rtl:nU.ltr};class rq extends sR{constructor(){super(...arguments),this.anchor="",this.viewport="",this.horizontalPositioningMode="uncontrolled",this.horizontalDefaultPosition="unset",this.horizontalViewportLock=!1,this.horizontalInset=!1,this.horizontalScaling="content",this.verticalPositioningMode="uncontrolled",this.verticalDefaultPosition="unset",this.verticalViewportLock=!1,this.verticalInset=!1,this.verticalScaling="content",this.fixedPlacement=!1,this.autoUpdateMode="anchor",this.anchorElement=null,this.viewportElement=null,this.initialLayoutComplete=!1,this.resizeDetector=null,this.baseHorizontalOffset=0,this.baseVerticalOffset=0,this.pendingPositioningUpdate=!1,this.pendingReset=!1,this.currentDirection=nU.ltr,this.regionVisible=!1,this.forceUpdate=!1,this.updateThreshold=.5,this.update=()=>{this.pendingPositioningUpdate||this.requestPositionUpdates()},this.startObservers=()=>{this.stopObservers(),null!==this.anchorElement&&(this.requestPositionUpdates(),null!==this.resizeDetector&&(this.resizeDetector.observe(this.anchorElement),this.resizeDetector.observe(this)))},this.requestPositionUpdates=()=>{null===this.anchorElement||this.pendingPositioningUpdate||(rq.intersectionService.requestPosition(this,this.handleIntersection),rq.intersectionService.requestPosition(this.anchorElement,this.handleIntersection),null!==this.viewportElement&&rq.intersectionService.requestPosition(this.viewportElement,this.handleIntersection),this.pendingPositioningUpdate=!0)},this.stopObservers=()=>{this.pendingPositioningUpdate&&(this.pendingPositioningUpdate=!1,rq.intersectionService.cancelRequestPosition(this,this.handleIntersection),null!==this.anchorElement&&rq.intersectionService.cancelRequestPosition(this.anchorElement,this.handleIntersection),null!==this.viewportElement&&rq.intersectionService.cancelRequestPosition(this.viewportElement,this.handleIntersection)),null!==this.resizeDetector&&this.resizeDetector.disconnect()},this.getViewport=()=>"string"!=typeof this.viewport||""===this.viewport?document.documentElement:document.getElementById(this.viewport),this.getAnchor=()=>document.getElementById(this.anchor),this.handleIntersection=e=>{!this.pendingPositioningUpdate||(this.pendingPositioningUpdate=!1,this.applyIntersectionEntries(e)&&this.updateLayout())},this.applyIntersectionEntries=e=>{let t=e.find(e=>e.target===this),i=e.find(e=>e.target===this.anchorElement),o=e.find(e=>e.target===this.viewportElement);return void 0!==t&&void 0!==o&&void 0!==i&&!!(!this.regionVisible||this.forceUpdate||void 0===this.regionRect||void 0===this.anchorRect||void 0===this.viewportRect||this.isRectDifferent(this.anchorRect,i.boundingClientRect)||this.isRectDifferent(this.viewportRect,o.boundingClientRect)||this.isRectDifferent(this.regionRect,t.boundingClientRect))&&(this.regionRect=t.boundingClientRect,this.anchorRect=i.boundingClientRect,this.viewportElement===document.documentElement?this.viewportRect=new DOMRectReadOnly(o.boundingClientRect.x+document.documentElement.scrollLeft,o.boundingClientRect.y+document.documentElement.scrollTop,o.boundingClientRect.width,o.boundingClientRect.height):this.viewportRect=o.boundingClientRect,this.updateRegionOffset(),this.forceUpdate=!1,!0)},this.updateRegionOffset=()=>{this.anchorRect&&this.regionRect&&(this.baseHorizontalOffset=this.baseHorizontalOffset+(this.anchorRect.left-this.regionRect.left)+(this.translateX-this.baseHorizontalOffset),this.baseVerticalOffset=this.baseVerticalOffset+(this.anchorRect.top-this.regionRect.top)+(this.translateY-this.baseVerticalOffset))},this.isRectDifferent=(e,t)=>!!(Math.abs(e.top-t.top)>this.updateThreshold||Math.abs(e.right-t.right)>this.updateThreshold||Math.abs(e.bottom-t.bottom)>this.updateThreshold||Math.abs(e.left-t.left)>this.updateThreshold),this.handleResize=e=>{this.update()},this.reset=()=>{this.pendingReset&&(this.pendingReset=!1,null===this.anchorElement&&(this.anchorElement=this.getAnchor()),null===this.viewportElement&&(this.viewportElement=this.getViewport()),this.currentDirection=rj(this),this.startObservers())},this.updateLayout=()=>{let e,t;if("uncontrolled"!==this.horizontalPositioningMode){let e=this.getPositioningOptions(this.horizontalInset);if("center"===this.horizontalDefaultPosition)t="center";else if("unset"!==this.horizontalDefaultPosition){let e=this.horizontalDefaultPosition;if("start"===e||"end"===e){let t=rj(this);if(t!==this.currentDirection){this.currentDirection=t,this.initialize();return}e=this.currentDirection===nU.ltr?"start"===e?"left":"right":"start"===e?"right":"left"}switch(e){case"left":t=this.horizontalInset?"insetStart":"start";break;case"right":t=this.horizontalInset?"insetEnd":"end"}}let i=void 0!==this.horizontalThreshold?this.horizontalThreshold:void 0!==this.regionRect?this.regionRect.width:0,o=void 0!==this.anchorRect?this.anchorRect.left:0,s=void 0!==this.anchorRect?this.anchorRect.right:0,r=void 0!==this.anchorRect?this.anchorRect.width:0,a=void 0!==this.viewportRect?this.viewportRect.left:0,n=void 0!==this.viewportRect?this.viewportRect.right:0;(void 0===t||"locktodefault"!==this.horizontalPositioningMode&&this.getAvailableSpace(t,o,s,r,a,n)<i)&&(t=this.getAvailableSpace(e[0],o,s,r,a,n)>this.getAvailableSpace(e[1],o,s,r,a,n)?e[0]:e[1])}if("uncontrolled"!==this.verticalPositioningMode){let t=this.getPositioningOptions(this.verticalInset);if("center"===this.verticalDefaultPosition)e="center";else if("unset"!==this.verticalDefaultPosition)switch(this.verticalDefaultPosition){case"top":e=this.verticalInset?"insetStart":"start";break;case"bottom":e=this.verticalInset?"insetEnd":"end"}let i=void 0!==this.verticalThreshold?this.verticalThreshold:void 0!==this.regionRect?this.regionRect.height:0,o=void 0!==this.anchorRect?this.anchorRect.top:0,s=void 0!==this.anchorRect?this.anchorRect.bottom:0,r=void 0!==this.anchorRect?this.anchorRect.height:0,a=void 0!==this.viewportRect?this.viewportRect.top:0,n=void 0!==this.viewportRect?this.viewportRect.bottom:0;(void 0===e||"locktodefault"!==this.verticalPositioningMode&&this.getAvailableSpace(e,o,s,r,a,n)<i)&&(e=this.getAvailableSpace(t[0],o,s,r,a,n)>this.getAvailableSpace(t[1],o,s,r,a,n)?t[0]:t[1])}let i=this.getNextRegionDimension(t,e),o=this.horizontalPosition!==t||this.verticalPosition!==e;if(this.setHorizontalPosition(t,i),this.setVerticalPosition(e,i),this.updateRegionStyle(),!this.initialLayoutComplete){this.initialLayoutComplete=!0,this.requestPositionUpdates();return}this.regionVisible||(this.regionVisible=!0,this.style.removeProperty("pointer-events"),this.style.removeProperty("opacity"),this.classList.toggle("loaded",!0),this.$emit("loaded",this,{bubbles:!1})),this.updatePositionClasses(),o&&this.$emit("positionchange",this,{bubbles:!1})},this.updateRegionStyle=()=>{this.style.width=this.regionWidth,this.style.height=this.regionHeight,this.style.transform=`translate(${this.translateX}px, ${this.translateY}px)`},this.updatePositionClasses=()=>{this.classList.toggle("top","start"===this.verticalPosition),this.classList.toggle("bottom","end"===this.verticalPosition),this.classList.toggle("inset-top","insetStart"===this.verticalPosition),this.classList.toggle("inset-bottom","insetEnd"===this.verticalPosition),this.classList.toggle("vertical-center","center"===this.verticalPosition),this.classList.toggle("left","start"===this.horizontalPosition),this.classList.toggle("right","end"===this.horizontalPosition),this.classList.toggle("inset-left","insetStart"===this.horizontalPosition),this.classList.toggle("inset-right","insetEnd"===this.horizontalPosition),this.classList.toggle("horizontal-center","center"===this.horizontalPosition)},this.setHorizontalPosition=(e,t)=>{if(void 0===e||void 0===this.regionRect||void 0===this.anchorRect||void 0===this.viewportRect)return;let i=0;switch(this.horizontalScaling){case"anchor":case"fill":i=this.horizontalViewportLock?this.viewportRect.width:t.width,this.regionWidth=`${i}px`;break;case"content":i=this.regionRect.width,this.regionWidth="unset"}let o=0;switch(e){case"start":this.translateX=this.baseHorizontalOffset-i,this.horizontalViewportLock&&this.anchorRect.left>this.viewportRect.right&&(this.translateX=this.translateX-(this.anchorRect.left-this.viewportRect.right));break;case"insetStart":this.translateX=this.baseHorizontalOffset-i+this.anchorRect.width,this.horizontalViewportLock&&this.anchorRect.right>this.viewportRect.right&&(this.translateX=this.translateX-(this.anchorRect.right-this.viewportRect.right));break;case"insetEnd":this.translateX=this.baseHorizontalOffset,this.horizontalViewportLock&&this.anchorRect.left<this.viewportRect.left&&(this.translateX=this.translateX-(this.anchorRect.left-this.viewportRect.left));break;case"end":this.translateX=this.baseHorizontalOffset+this.anchorRect.width,this.horizontalViewportLock&&this.anchorRect.right<this.viewportRect.left&&(this.translateX=this.translateX-(this.anchorRect.right-this.viewportRect.left));break;case"center":if(o=(this.anchorRect.width-i)/2,this.translateX=this.baseHorizontalOffset+o,this.horizontalViewportLock){let e=this.anchorRect.left+o,t=this.anchorRect.right-o;e<this.viewportRect.left&&!(t>this.viewportRect.right)?this.translateX=this.translateX-(e-this.viewportRect.left):t>this.viewportRect.right&&!(e<this.viewportRect.left)&&(this.translateX=this.translateX-(t-this.viewportRect.right))}}this.horizontalPosition=e},this.setVerticalPosition=(e,t)=>{if(void 0===e||void 0===this.regionRect||void 0===this.anchorRect||void 0===this.viewportRect)return;let i=0;switch(this.verticalScaling){case"anchor":case"fill":i=this.verticalViewportLock?this.viewportRect.height:t.height,this.regionHeight=`${i}px`;break;case"content":i=this.regionRect.height,this.regionHeight="unset"}let o=0;switch(e){case"start":this.translateY=this.baseVerticalOffset-i,this.verticalViewportLock&&this.anchorRect.top>this.viewportRect.bottom&&(this.translateY=this.translateY-(this.anchorRect.top-this.viewportRect.bottom));break;case"insetStart":this.translateY=this.baseVerticalOffset-i+this.anchorRect.height,this.verticalViewportLock&&this.anchorRect.bottom>this.viewportRect.bottom&&(this.translateY=this.translateY-(this.anchorRect.bottom-this.viewportRect.bottom));break;case"insetEnd":this.translateY=this.baseVerticalOffset,this.verticalViewportLock&&this.anchorRect.top<this.viewportRect.top&&(this.translateY=this.translateY-(this.anchorRect.top-this.viewportRect.top));break;case"end":this.translateY=this.baseVerticalOffset+this.anchorRect.height,this.verticalViewportLock&&this.anchorRect.bottom<this.viewportRect.top&&(this.translateY=this.translateY-(this.anchorRect.bottom-this.viewportRect.top));break;case"center":if(o=(this.anchorRect.height-i)/2,this.translateY=this.baseVerticalOffset+o,this.verticalViewportLock){let e=this.anchorRect.top+o,t=this.anchorRect.bottom-o;e<this.viewportRect.top&&!(t>this.viewportRect.bottom)?this.translateY=this.translateY-(e-this.viewportRect.top):t>this.viewportRect.bottom&&!(e<this.viewportRect.top)&&(this.translateY=this.translateY-(t-this.viewportRect.bottom))}}this.verticalPosition=e},this.getPositioningOptions=e=>e?["insetStart","insetEnd"]:["start","end"],this.getAvailableSpace=(e,t,i,o,s,r)=>{let a=t-s,n=r-(t+o);switch(e){case"start":return a;case"insetStart":return a+o;case"insetEnd":return n+o;case"end":return n;case"center":return 2*Math.min(a,n)+o}},this.getNextRegionDimension=(e,t)=>{let i={height:void 0!==this.regionRect?this.regionRect.height:0,width:void 0!==this.regionRect?this.regionRect.width:0};return void 0!==e&&"fill"===this.horizontalScaling?i.width=this.getAvailableSpace(e,void 0!==this.anchorRect?this.anchorRect.left:0,void 0!==this.anchorRect?this.anchorRect.right:0,void 0!==this.anchorRect?this.anchorRect.width:0,void 0!==this.viewportRect?this.viewportRect.left:0,void 0!==this.viewportRect?this.viewportRect.right:0):"anchor"===this.horizontalScaling&&(i.width=void 0!==this.anchorRect?this.anchorRect.width:0),void 0!==t&&"fill"===this.verticalScaling?i.height=this.getAvailableSpace(t,void 0!==this.anchorRect?this.anchorRect.top:0,void 0!==this.anchorRect?this.anchorRect.bottom:0,void 0!==this.anchorRect?this.anchorRect.height:0,void 0!==this.viewportRect?this.viewportRect.top:0,void 0!==this.viewportRect?this.viewportRect.bottom:0):"anchor"===this.verticalScaling&&(i.height=void 0!==this.anchorRect?this.anchorRect.height:0),i},this.startAutoUpdateEventListeners=()=>{window.addEventListener("resize",this.update,{passive:!0}),window.addEventListener("scroll",this.update,{passive:!0,capture:!0}),null!==this.resizeDetector&&null!==this.viewportElement&&this.resizeDetector.observe(this.viewportElement)},this.stopAutoUpdateEventListeners=()=>{window.removeEventListener("resize",this.update),window.removeEventListener("scroll",this.update),null!==this.resizeDetector&&null!==this.viewportElement&&this.resizeDetector.unobserve(this.viewportElement)}}anchorChanged(){this.initialLayoutComplete&&(this.anchorElement=this.getAnchor())}viewportChanged(){this.initialLayoutComplete&&(this.viewportElement=this.getViewport())}horizontalPositioningModeChanged(){this.requestReset()}horizontalDefaultPositionChanged(){this.updateForAttributeChange()}horizontalViewportLockChanged(){this.updateForAttributeChange()}horizontalInsetChanged(){this.updateForAttributeChange()}horizontalThresholdChanged(){this.updateForAttributeChange()}horizontalScalingChanged(){this.updateForAttributeChange()}verticalPositioningModeChanged(){this.requestReset()}verticalDefaultPositionChanged(){this.updateForAttributeChange()}verticalViewportLockChanged(){this.updateForAttributeChange()}verticalInsetChanged(){this.updateForAttributeChange()}verticalThresholdChanged(){this.updateForAttributeChange()}verticalScalingChanged(){this.updateForAttributeChange()}fixedPlacementChanged(){this.$fastController.isConnected&&this.initialLayoutComplete&&this.initialize()}autoUpdateModeChanged(e,t){this.$fastController.isConnected&&this.initialLayoutComplete&&("auto"===e&&this.stopAutoUpdateEventListeners(),"auto"===t&&this.startAutoUpdateEventListeners())}anchorElementChanged(){this.requestReset()}viewportElementChanged(){this.$fastController.isConnected&&this.initialLayoutComplete&&this.initialize()}connectedCallback(){super.connectedCallback(),"auto"===this.autoUpdateMode&&this.startAutoUpdateEventListeners(),this.initialize()}disconnectedCallback(){super.disconnectedCallback(),"auto"===this.autoUpdateMode&&this.stopAutoUpdateEventListeners(),this.stopObservers(),this.disconnectResizeDetector()}adoptedCallback(){this.initialize()}disconnectResizeDetector(){null!==this.resizeDetector&&(this.resizeDetector.disconnect(),this.resizeDetector=null)}initializeResizeDetector(){this.disconnectResizeDetector(),this.resizeDetector=new window.ResizeObserver(this.handleResize)}updateForAttributeChange(){this.$fastController.isConnected&&this.initialLayoutComplete&&(this.forceUpdate=!0,this.update())}initialize(){this.initializeResizeDetector(),null===this.anchorElement&&(this.anchorElement=this.getAnchor()),this.requestReset()}requestReset(){this.$fastController.isConnected&&!1===this.pendingReset&&(this.setInitialState(),el.queueUpdate(()=>this.reset()),this.pendingReset=!0)}setInitialState(){this.initialLayoutComplete=!1,this.regionVisible=!1,this.translateX=0,this.translateY=0,this.baseHorizontalOffset=0,this.baseVerticalOffset=0,this.viewportRect=void 0,this.regionRect=void 0,this.anchorRect=void 0,this.verticalPosition=void 0,this.horizontalPosition=void 0,this.style.opacity="0",this.style.pointerEvents="none",this.forceUpdate=!1,this.style.position=this.fixedPlacement?"fixed":"absolute",this.updatePositionClasses(),this.updateRegionStyle()}}rq.intersectionService=new class{constructor(){this.intersectionDetector=null,this.observedElements=new Map,this.requestPosition=(e,t)=>{var i;if(null!==this.intersectionDetector){if(this.observedElements.has(e)){null==(i=this.observedElements.get(e))||i.push(t);return}this.observedElements.set(e,[t]),this.intersectionDetector.observe(e)}},this.cancelRequestPosition=(e,t)=>{let i=this.observedElements.get(e);if(void 0!==i){let e=i.indexOf(t);-1!==e&&i.splice(e,1)}},this.initializeIntersectionDetector=()=>{Q.IntersectionObserver&&(this.intersectionDetector=new IntersectionObserver(this.handleIntersection,{root:null,rootMargin:"0px",threshold:[0,1]}))},this.handleIntersection=e=>{if(null===this.intersectionDetector)return;let t=[],i=[];e.forEach(e=>{var o;null==(o=this.intersectionDetector)||o.unobserve(e.target);let s=this.observedElements.get(e.target);void 0!==s&&(s.forEach(o=>{let s=t.indexOf(o);-1===s&&(s=t.length,t.push(o),i.push([])),i[s].push(e)}),this.observedElements.delete(e.target))}),t.forEach((e,t)=>{e(i[t])})},this.initializeIntersectionDetector()}},X([eS],rq.prototype,"anchor",void 0),X([eS],rq.prototype,"viewport",void 0),X([eS({attribute:"horizontal-positioning-mode"})],rq.prototype,"horizontalPositioningMode",void 0),X([eS({attribute:"horizontal-default-position"})],rq.prototype,"horizontalDefaultPosition",void 0),X([eS({attribute:"horizontal-viewport-lock",mode:"boolean"})],rq.prototype,"horizontalViewportLock",void 0),X([eS({attribute:"horizontal-inset",mode:"boolean"})],rq.prototype,"horizontalInset",void 0),X([eS({attribute:"horizontal-threshold"})],rq.prototype,"horizontalThreshold",void 0),X([eS({attribute:"horizontal-scaling"})],rq.prototype,"horizontalScaling",void 0),X([eS({attribute:"vertical-positioning-mode"})],rq.prototype,"verticalPositioningMode",void 0),X([eS({attribute:"vertical-default-position"})],rq.prototype,"verticalDefaultPosition",void 0),X([eS({attribute:"vertical-viewport-lock",mode:"boolean"})],rq.prototype,"verticalViewportLock",void 0),X([eS({attribute:"vertical-inset",mode:"boolean"})],rq.prototype,"verticalInset",void 0),X([eS({attribute:"vertical-threshold"})],rq.prototype,"verticalThreshold",void 0),X([eS({attribute:"vertical-scaling"})],rq.prototype,"verticalScaling",void 0),X([eS({attribute:"fixed-placement",mode:"boolean"})],rq.prototype,"fixedPlacement",void 0),X([eS({attribute:"auto-update-mode"})],rq.prototype,"autoUpdateMode",void 0),X([eu],rq.prototype,"anchorElement",void 0),X([eu],rq.prototype,"viewportElement",void 0),X([eu],rq.prototype,"initialLayoutComplete",void 0);let rU=()=>null;function r_(e){return void 0===e?rU:"function"==typeof e?e:()=>e}function rG(e,t,i){let o="function"==typeof e?e:()=>e,s=r_(t),r=r_(i);return(e,t)=>o(e,t)?s(e,t):r(e,t)}let rW=(e,t)=>rh`
  :host {
    contain: layout;
    display: block;
  }
`;class rK extends rq{}let rX=rK.compose({baseName:"anchored-region",baseClass:rq,template:(e,t)=>s4`
    <template class="${e=>e.initialLayoutComplete?"loaded":""}">
        ${rG(e=>e.initialLayoutComplete,s4`
                <slot></slot>
            `)}
    </template>
`,styles:rW});class rY extends sR{connectedCallback(){super.connectedCallback(),this.shape||(this.shape="circle")}}X([eS],rY.prototype,"fill",void 0),X([eS],rY.prototype,"color",void 0),X([eS],rY.prototype,"link",void 0),X([eS],rY.prototype,"shape",void 0);class rQ extends sR{constructor(){super(...arguments),this.generateBadgeStyle=()=>{if(!this.fill&&!this.color)return;let e=`background-color: var(--badge-fill-${this.fill});`,t=`color: var(--badge-color-${this.color});`;return this.fill&&!this.color?e:this.color&&!this.fill?t:`${t} ${e}`}}}X([eS({attribute:"fill"})],rQ.prototype,"fill",void 0),X([eS({attribute:"color"})],rQ.prototype,"color",void 0),X([eS({mode:"boolean"})],rQ.prototype,"circular",void 0);class rZ{constructor(e,t){this.cache=new WeakMap,this.ltr=e,this.rtl=t}bind(e){this.attach(e)}unbind(e){let t=this.cache.get(e);t&&tx.unsubscribe(t)}attach(e){let t=this.cache.get(e)||new rJ(this.ltr,this.rtl,e),i=tx.getValueFor(e);tx.subscribe(t),t.attach(i),this.cache.set(e,t)}}class rJ{constructor(e,t,i){this.ltr=e,this.rtl=t,this.source=i,this.attached=null}handleChange({target:e,token:t}){this.attach(t.getValueFor(e))}attach(e){this.attached!==this[e]&&(null!==this.attached&&this.source.$fastController.removeStyles(this.attached),this.attached=this[e],null!==this.attached&&this.source.$fastController.addStyles(this.attached))}}let r0=(e,t)=>rh`
    ${ru("flex")} :host {
      position: relative;
      height: var(--avatar-size, var(--avatar-size-default));
      width: var(--avatar-size, var(--avatar-size-default));
      --avatar-size-default: calc(
        (
            (${tp} + ${tf}) * ${tv} +
              ((${tv} * 8) - 40)
          ) * 1px
      );
      --avatar-text-size: ${tC};
      --avatar-text-ratio: ${tv};
    }

    .link {
      text-decoration: none;
      color: ${op};
      display: flex;
      flex-direction: row;
      justify-content: center;
      align-items: center;
      min-width: 100%;
    }

    .square {
      border-radius: calc(${tb} * 1px);
      min-width: 100%;
      overflow: hidden;
    }

    .circle {
      border-radius: 100%;
      min-width: 100%;
      overflow: hidden;
    }

    .backplate {
      position: relative;
      display: flex;
      background-color: ${iS};
    }

    .media,
    ::slotted(img) {
      max-width: 100%;
      position: absolute;
      display: block;
    }

    .content {
      font-size: calc(
        (
            var(--avatar-text-size) +
              var(--avatar-size, var(--avatar-size-default))
          ) / var(--avatar-text-ratio)
      );
      line-height: var(--avatar-size, var(--avatar-size-default));
      display: block;
      min-height: var(--avatar-size, var(--avatar-size-default));
    }

    ::slotted(${e.tagFor(rQ)}) {
      position: absolute;
      display: block;
    }
  `.withBehaviors(new rZ(rh`
  ::slotted(${e.tagFor(rQ)}) {
    right: 0;
  }
`,rh`
  ::slotted(${e.tagFor(rQ)}) {
    left: 0;
  }
`));class r1 extends rY{}rI([eS({attribute:"src"})],r1.prototype,"imgSrc",void 0),rI([eS],r1.prototype,"alt",void 0);let r2=s4`
  ${rG(e=>e.imgSrc,s4`
      <img
        src="${e=>e.imgSrc}"
        alt="${e=>e.alt}"
        slot="media"
        class="media"
        part="media"
      />
    `)}
`,r5=r1.compose({baseName:"avatar",baseClass:rY,template:(e,t)=>s4`
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
`,styles:r0,media:r2,shadowOptions:{delegatesFocus:!0}}),r3=(e,t)=>rh`
  ${ru("inline-block")} :host {
    box-sizing: border-box;
    font-family: ${tu};
    font-size: ${tT};
    line-height: ${tF};
  }

  .control {
    border-radius: calc(${tb} * 1px);
    padding: calc(((${tv} * 0.5) - ${tw}) * 1px)
      calc((${tv} - ${tw}) * 1px);
    color: ${op};
    font-weight: 600;
    border: calc(${tw} * 1px) solid transparent;
    background-color: ${iX};
  }

  .control[style] {
    font-weight: 400;
  }

  :host([circular]) .control {
    border-radius: 100px;
    padding: 0 calc(${tv} * 1px);
    height: calc((${ry} - (${tv} * 3)) * 1px);
    min-width: calc((${ry} - (${tv} * 3)) * 1px);
    display: flex;
    align-items: center;
    justify-content: center;
    box-sizing: border-box;
  }
`;class r4 extends rQ{}let r6=r4.compose({baseName:"badge",baseClass:rQ,template:(e,t)=>s4`
    <template class="${e=>e.circular?"circular":""}">
        <div class="control" part="control" style="${e=>e.generateBadgeStyle()}">
            <slot></slot>
        </div>
    </template>
`,styles:r3});class r9 extends rF{constructor(){super(...arguments),this.separator=!0}}X([eu],r9.prototype,"separator",void 0),rt(r9,s8,rS);class r8 extends sR{slottedBreadcrumbItemsChanged(){if(this.$fastController.isConnected){if(void 0===this.slottedBreadcrumbItems||0===this.slottedBreadcrumbItems.length)return;let e=this.slottedBreadcrumbItems[this.slottedBreadcrumbItems.length-1];this.slottedBreadcrumbItems.forEach(t=>{let i=t===e;this.setItemSeparator(t,i),this.setAriaCurrent(t,i)})}}setItemSeparator(e,t){e instanceof r9&&(e.separator=!t)}findChildWithHref(e){var t,i;return e.childElementCount>0?e.querySelector("a[href]"):(null==(t=e.shadowRoot)?void 0:t.childElementCount)?null==(i=e.shadowRoot)?void 0:i.querySelector("a[href]"):null}setAriaCurrent(e,t){let i=this.findChildWithHref(e);null===i&&e.hasAttribute("href")&&e instanceof r9?t?e.setAttribute("aria-current","page"):e.removeAttribute("aria-current"):null!==i&&(t?i.setAttribute("aria-current","page"):i.removeAttribute("aria-current"))}}X([eu],r8.prototype,"slottedBreadcrumbItems",void 0);let r7=(e,t)=>rh`
  ${ru("inline-block")} :host {
    box-sizing: border-box;
    font-family: ${tu};
    font-size: ${tC};
    line-height: ${tI};
  }

  .list {
    display: flex;
    flex-wrap: wrap;
  }
`;class ae extends r8{}let at=ae.compose({baseName:"breadcrumb",baseClass:r8,template:(e,t)=>s4`
    <template role="navigation">
        <div role="list" class="list" part="list">
            <slot
                ${rn({property:"slottedBreadcrumbItems",filter:rs()})}
            ></slot>
        </div>
    </template>
`,styles:r7}),ai=(e,t)=>rh`
    ${ru("inline-flex")} :host {
        background: transparent;
        box-sizing: border-box;
        font-family: ${tu};
        font-size: ${tC};
        fill: currentColor;
        line-height: ${tI};
        min-width: calc(${ry} * 1px);
        outline: none;
        color: ${op}
    }

    .listitem {
        display: flex;
        align-items: center;
        width: max-content;
    }

    .separator {
        margin: 0 6px;
        display: flex;
    }

    .control {
        align-items: center;
        box-sizing: border-box;
        color: ${iU};
        cursor: pointer;
        display: flex;
        fill: inherit;
        outline: none;
        text-decoration: none;
        white-space: nowrap;
    }

    .control:hover {
        color: ${i_};
    }

    .control:active {
        color: ${iG};
    }

    .control .content {
        position: relative;
    }

    .control .content::before {
        content: "";
        display: block;
        height: calc(${tw} * 1px);
        left: 0;
        position: absolute;
        right: 0;
        top: calc(1em + 4px);
        width: 100%;
    }

    .control:hover .content::before {
        background: ${i_};
    }

    .control:active .content::before {
        background: ${iG};
    }

    .control:${rm} .content::before {
        background: ${iW};
        height: calc(${tk} * 1px);
    }

    .control:not([href]) {
        color: ${op};
        cursor: default;
    }

    .control:not([href]) .content::before {
        background: none;
    }

    .start,
    .end {
        display: flex;
    }

    ::slotted(svg) {
        /* TODO: adaptive typography https://github.com/microsoft/fast/issues/2432 */
        width: 16px;
        height: 16px;
    }

    .start {
        margin-inline-end: 6px;
    }

    .end {
        margin-inline-start: 6px;
    }
`.withBehaviors(rv(rh`
      .control:hover .content::before,
                .control:${rm} .content::before {
        background: ${nX.LinkText};
      }
      .start,
      .end {
        fill: ${nX.ButtonText};
      }
    `));class ao extends r9{}let as=ao.compose({baseName:"breadcrumb-item",baseClass:r9,template:(e,t)=>s4`
    <div role="listitem" class="listitem" part="listitem">
        ${rG(e=>e.href&&e.href.length>0,s4`
                ${rO(e,t)}
            `)}
        ${rG(e=>!e.href,s4`
                ${re(e,t)}
                <slot></slot>
                ${s7(e,t)}
            `)}
        ${rG(e=>e.separator,s4`
                <span class="separator" part="separator" aria-hidden="true">
                    <slot name="separator">${t.separator||""}</slot>
                </span>
            `)}
    </div>
`,styles:ai,separator:"/",shadowOptions:{delegatesFocus:!0}}),ar="form-associated-proxy",aa="ElementInternals",an=aa in window&&"setFormValue"in window[aa].prototype,al=new WeakMap;function ah(e){let t=class extends e{constructor(...e){super(...e),this.dirtyValue=!1,this.disabled=!1,this.proxyEventsToBlock=["change","click"],this.proxyInitialized=!1,this.required=!1,this.initialValue=this.initialValue||"",this.elementInternals||(this.formResetCallback=this.formResetCallback.bind(this))}static get formAssociated(){return an}get validity(){return this.elementInternals?this.elementInternals.validity:this.proxy.validity}get form(){return this.elementInternals?this.elementInternals.form:this.proxy.form}get validationMessage(){return this.elementInternals?this.elementInternals.validationMessage:this.proxy.validationMessage}get willValidate(){return this.elementInternals?this.elementInternals.willValidate:this.proxy.willValidate}get labels(){if(this.elementInternals)return Object.freeze(Array.from(this.elementInternals.labels));if(!(this.proxy instanceof HTMLElement)||!this.proxy.ownerDocument||!this.id)return ee;{let e=this.proxy.labels,t=Array.from(this.proxy.getRootNode().querySelectorAll(`[for='${this.id}']`));return Object.freeze(e?t.concat(Array.from(e)):t)}}valueChanged(e,t){this.dirtyValue=!0,this.proxy instanceof HTMLElement&&(this.proxy.value=this.value),this.currentValue=this.value,this.setFormValue(this.value),this.validate()}currentValueChanged(){this.value=this.currentValue}initialValueChanged(e,t){this.dirtyValue||(this.value=this.initialValue,this.dirtyValue=!1)}disabledChanged(e,t){this.proxy instanceof HTMLElement&&(this.proxy.disabled=this.disabled),el.queueUpdate(()=>this.classList.toggle("disabled",this.disabled))}nameChanged(e,t){this.proxy instanceof HTMLElement&&(this.proxy.name=this.name)}requiredChanged(e,t){this.proxy instanceof HTMLElement&&(this.proxy.required=this.required),el.queueUpdate(()=>this.classList.toggle("required",this.required)),this.validate()}get elementInternals(){if(!an)return null;let e=al.get(this);return e||(e=this.attachInternals(),al.set(this,e)),e}connectedCallback(){super.connectedCallback(),this.addEventListener("keypress",this._keypressHandler),this.value||(this.value=this.initialValue,this.dirtyValue=!1),!this.elementInternals&&(this.attachProxy(),this.form&&this.form.addEventListener("reset",this.formResetCallback))}disconnectedCallback(){super.disconnectedCallback(),this.proxyEventsToBlock.forEach(e=>this.proxy.removeEventListener(e,this.stopPropagation)),!this.elementInternals&&this.form&&this.form.removeEventListener("reset",this.formResetCallback)}checkValidity(){return this.elementInternals?this.elementInternals.checkValidity():this.proxy.checkValidity()}reportValidity(){return this.elementInternals?this.elementInternals.reportValidity():this.proxy.reportValidity()}setValidity(e,t,i){this.elementInternals?this.elementInternals.setValidity(e,t,i):"string"==typeof t&&this.proxy.setCustomValidity(t)}formDisabledCallback(e){this.disabled=e}formResetCallback(){this.value=this.initialValue,this.dirtyValue=!1}attachProxy(){var e;this.proxyInitialized||(this.proxyInitialized=!0,this.proxy.style.display="none",this.proxyEventsToBlock.forEach(e=>this.proxy.addEventListener(e,this.stopPropagation)),this.proxy.disabled=this.disabled,this.proxy.required=this.required,"string"==typeof this.name&&(this.proxy.name=this.name),"string"==typeof this.value&&(this.proxy.value=this.value),this.proxy.setAttribute("slot",ar),this.proxySlot=document.createElement("slot"),this.proxySlot.setAttribute("name",ar)),null==(e=this.shadowRoot)||e.appendChild(this.proxySlot),this.appendChild(this.proxy)}detachProxy(){var e;this.removeChild(this.proxy),null==(e=this.shadowRoot)||e.removeChild(this.proxySlot)}validate(e){this.proxy instanceof HTMLElement&&this.setValidity(this.proxy.validity,this.proxy.validationMessage,e)}setFormValue(e,t){this.elementInternals&&this.elementInternals.setFormValue(e,t||e)}_keypressHandler(e){if("Enter"===e.key&&this.form instanceof HTMLFormElement){let e=this.form.querySelector("[type=submit]");null==e||e.click()}}stopPropagation(e){e.stopPropagation()}};return eS({mode:"boolean"})(t.prototype,"disabled"),eS({mode:"fromView",attribute:"value"})(t.prototype,"initialValue"),eS({attribute:"current-value"})(t.prototype,"currentValue"),eS(t.prototype,"name"),eS({mode:"boolean"})(t.prototype,"required"),eu(t.prototype,"value"),t}function ad(e){class t extends ah(e){}class i extends t{constructor(...e){super(e),this.dirtyChecked=!1,this.checkedAttribute=!1,this.checked=!1,this.dirtyChecked=!1}checkedAttributeChanged(){this.defaultChecked=this.checkedAttribute}defaultCheckedChanged(){this.dirtyChecked||(this.checked=this.defaultChecked,this.dirtyChecked=!1)}checkedChanged(e,t){this.dirtyChecked||(this.dirtyChecked=!0),this.currentChecked=this.checked,this.updateForm(),this.proxy instanceof HTMLInputElement&&(this.proxy.checked=this.checked),void 0!==e&&this.$emit("change"),this.validate()}currentCheckedChanged(e,t){this.checked=this.currentChecked}updateForm(){let e=this.checked?this.value:null;this.setFormValue(e,e)}connectedCallback(){super.connectedCallback(),this.updateForm()}formResetCallback(){super.formResetCallback(),this.checked=!!this.checkedAttribute,this.dirtyChecked=!1}}return eS({attribute:"checked",mode:"boolean"})(i.prototype,"checkedAttribute"),eS({attribute:"current-checked",converter:eI})(i.prototype,"currentChecked"),eu(i.prototype,"defaultChecked"),eu(i.prototype,"checked"),i}class ac extends sR{}class au extends ah(ac){constructor(){super(...arguments),this.proxy=document.createElement("input")}}class ap extends au{constructor(){super(...arguments),this.handleClick=e=>{var t;this.disabled&&(null==(t=this.defaultSlottedContent)?void 0:t.length)<=1&&e.stopPropagation()},this.handleSubmission=()=>{if(!this.form)return;let e=this.proxy.isConnected;e||this.attachProxy(),"function"==typeof this.form.requestSubmit?this.form.requestSubmit(this.proxy):this.proxy.click(),e||this.detachProxy()},this.handleFormReset=()=>{var e;null==(e=this.form)||e.reset()},this.handleUnsupportedDelegatesFocus=()=>{var e;window.ShadowRoot&&!window.ShadowRoot.prototype.hasOwnProperty("delegatesFocus")&&(null==(e=this.$fastController.definition.shadowOptions)?void 0:e.delegatesFocus)&&(this.focus=()=>{this.control.focus()})}}formactionChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.formAction=this.formaction)}formenctypeChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.formEnctype=this.formenctype)}formmethodChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.formMethod=this.formmethod)}formnovalidateChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.formNoValidate=this.formnovalidate)}formtargetChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.formTarget=this.formtarget)}typeChanged(e,t){this.proxy instanceof HTMLInputElement&&(this.proxy.type=this.type),"submit"===t&&this.addEventListener("click",this.handleSubmission),"submit"===e&&this.removeEventListener("click",this.handleSubmission),"reset"===t&&this.addEventListener("click",this.handleFormReset),"reset"===e&&this.removeEventListener("click",this.handleFormReset)}validate(){super.validate(this.control)}connectedCallback(){var e;super.connectedCallback(),this.proxy.setAttribute("type",this.type),this.handleUnsupportedDelegatesFocus();let t=Array.from(null==(e=this.control)?void 0:e.children);t&&t.forEach(e=>{e.addEventListener("click",this.handleClick)})}disconnectedCallback(){var e;super.disconnectedCallback();let t=Array.from(null==(e=this.control)?void 0:e.children);t&&t.forEach(e=>{e.removeEventListener("click",this.handleClick)})}}X([eS({mode:"boolean"})],ap.prototype,"autofocus",void 0),X([eS({attribute:"form"})],ap.prototype,"formId",void 0),X([eS],ap.prototype,"formaction",void 0),X([eS],ap.prototype,"formenctype",void 0),X([eS],ap.prototype,"formmethod",void 0),X([eS({mode:"boolean"})],ap.prototype,"formnovalidate",void 0),X([eS],ap.prototype,"formtarget",void 0),X([eS],ap.prototype,"type",void 0),X([eu],ap.prototype,"defaultSlottedContent",void 0);class ag{}X([eS({attribute:"aria-expanded"})],ag.prototype,"ariaExpanded",void 0),X([eS({attribute:"aria-pressed"})],ag.prototype,"ariaPressed",void 0),rt(ag,rT),rt(ap,s8,ag);let am="not-allowed",ab=(e,t)=>rh`
    :host([disabled]),
    :host([disabled]:hover),
    :host([disabled]:active) {
      opacity: ${t$};
      background-color: ${iX};
      cursor: ${am};
    }

    ${rD}
  `.withBehaviors(rv(rh`
      :host([disabled]),
      :host([disabled]) .control,
      :host([disabled]:hover),
      :host([disabled]:active) {
        forced-color-adjust: none;
        background-color: ${nX.ButtonFace};
        outline-color: ${nX.GrayText};
        color: ${nX.GrayText};
        cursor: ${am};
        opacity: 1;
      }
    `),rz("accent",rh`
        :host([appearance='accent'][disabled]),
        :host([appearance='accent'][disabled]:hover),
        :host([appearance='accent'][disabled]:active) {
          background: ${iS};
        }

        ${rE}
      `.withBehaviors(rv(rh`
          :host([appearance='accent'][disabled]) .control,
          :host([appearance='accent'][disabled]) .control:hover {
            background: ${nX.ButtonFace};
            border-color: ${nX.GrayText};
            color: ${nX.GrayText};
          }
        `))),rz("error",rh`
        :host([appearance='error'][disabled]),
        :host([appearance='error'][disabled]:hover),
        :host([appearance='error'][disabled]:active) {
          background: ${oI};
        }

        ${rR}
      `.withBehaviors(rv(rh`
          :host([appearance='error'][disabled]) .control,
          :host([appearance='error'][disabled]) .control:hover {
            background: ${nX.ButtonFace};
            border-color: ${nX.GrayText};
            color: ${nX.GrayText};
          }
        `))),rz("lightweight",rh`
        :host([appearance='lightweight'][disabled]:hover),
        :host([appearance='lightweight'][disabled]:active) {
          background-color: transparent;
          color: ${iU};
        }

        :host([appearance='lightweight'][disabled]) .content::before,
        :host([appearance='lightweight'][disabled]:hover) .content::before,
        :host([appearance='lightweight'][disabled]:active) .content::before {
          background: transparent;
        }

        ${rA}
      `.withBehaviors(rv(rh`
          :host([appearance='lightweight'].disabled) .control {
            forced-color-adjust: none;
            color: ${nX.GrayText};
          }

          :host([appearance='lightweight'].disabled)
            .control:hover
            .content::before {
            background: none;
          }
        `))),rz("outline",rh`
        :host([appearance='outline'][disabled]),
        :host([appearance='outline'][disabled]:hover),
        :host([appearance='outline'][disabled]:active) {
          background: transparent;
          border-color: ${iS};
        }

        ${rV}
      `.withBehaviors(rv(rh`
          :host([appearance='outline'][disabled]) .control {
            border-color: ${nX.GrayText};
          }
        `))),rz("stealth",rh`
        ${rP}
      `.withBehaviors(rv(rh`
          :host([appearance='stealth'][disabled]) {
            background: ${nX.ButtonFace};
          }

          :host([appearance='stealth'][disabled]) .control {
            background: ${nX.ButtonFace};
            border-color: transparent;
            color: ${nX.GrayText};
          }
        `))));class af extends ap{constructor(){super(...arguments),this.appearance="neutral"}defaultSlottedContentChanged(e,t){let i=this.defaultSlottedContent.filter(e=>e.nodeType===Node.ELEMENT_NODE);1===i.length&&(i[0]instanceof SVGElement||i[0].classList.contains("fa")||i[0].classList.contains("fas"))?this.control.classList.add("icon-only"):this.control.classList.remove("icon-only")}}rI([eS],af.prototype,"appearance",void 0),rI([eS({attribute:"minimal",mode:"boolean"})],af.prototype,"minimal",void 0),rI([eS],af.prototype,"scale",void 0);let av=af.compose({baseName:"button",baseClass:ap,template:(e,t)=>s4`
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
        ${s9("control")}
    >
        ${re(e,t)}
        <span class="content" part="content">
            <slot ${rn("defaultSlottedContent")}></slot>
        </span>
        ${s7(e,t)}
    </button>
`,styles:ab,shadowOptions:{delegatesFocus:!0}});class ay extends sR{}let ax="box-shadow: 0 0 calc((var(--elevation) * 0.225px) + 2px) rgba(0, 0, 0, calc(.11 * (2 - var(--background-luminance, 1)))), 0 calc(var(--elevation) * 0.4px) calc((var(--elevation) * 0.9px)) rgba(0, 0, 0, calc(.13 * (2 - var(--background-luminance, 1))));",a$=(e,t)=>rh`
    ${ru("block")} :host {
      --elevation: 4;
      display: block;
      contain: content;
      height: var(--card-height, 100%);
      width: var(--card-width, 100%);
      box-sizing: border-box;
      background: ${iT};
      border-radius: calc(${tb} * 1px);
      ${ax}
    }
  `.withBehaviors(rv(rh`
      :host {
        forced-color-adjust: none;
        background: ${nX.Canvas};
        box-shadow: 0 0 0 1px ${nX.CanvasText};
      }
    `));class aw extends ay{connectedCallback(){super.connectedCallback();let e=eM(this);e&&iT.setValueFor(this,t=>os.getValueFor(t).evaluate(t,iT.getValueFor(e)))}}let ak=aw.compose({baseName:"card",baseClass:ay,template:(e,t)=>s4`
    <slot></slot>
`,styles:a$});class aC extends sR{}class aI extends ad(aC){constructor(){super(...arguments),this.proxy=document.createElement("input")}}class aT extends aI{constructor(){super(),this.initialValue="on",this.indeterminate=!1,this.keypressHandler=e=>{this.readOnly||" "===e.key&&(this.indeterminate&&(this.indeterminate=!1),this.checked=!this.checked)},this.clickHandler=e=>{this.disabled||this.readOnly||(this.indeterminate&&(this.indeterminate=!1),this.checked=!this.checked)},this.proxy.setAttribute("type","checkbox")}readOnlyChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.readOnly=this.readOnly)}}X([eS({attribute:"readonly",mode:"boolean"})],aT.prototype,"readOnly",void 0),X([eu],aT.prototype,"defaultSlottedNodes",void 0),X([eu],aT.prototype,"indeterminate",void 0);let aF=(e,t)=>rh`
    ${ru("inline-flex")} :host {
      align-items: center;
      outline: none;
      margin: calc(${tv} * 1px) 0;
      /* Chromium likes to select label text or the default slot when the checkbox is
            clicked. Maybe there is a better solution here? */
      user-select: none;
    }

    .control {
      position: relative;
      width: calc((${ry} / 2 + ${tv}) * 1px);
      height: calc((${ry} / 2 + ${tv}) * 1px);
      box-sizing: border-box;
      border-radius: calc(${tb} * 1px);
      border: calc(${tw} * 1px) solid ${om};
      background: ${i0};
      outline: none;
      cursor: pointer;
    }

    :host([aria-invalid='true']) .control {
      border-color: ${oI};
    }

    .label {
      font-family: ${tu};
      color: ${op};
      /* Need to discuss with Brian how HorizontalSpacingNumber can work.
            https://github.com/microsoft/fast/issues/2766 */
      padding-inline-start: calc(${tv} * 2px + 2px);
      margin-inline-end: calc(${tv} * 2px + 2px);
      cursor: pointer;
      font-size: ${tC};
      line-height: ${tI};
    }

    .label__hidden {
      display: none;
      visibility: hidden;
    }

    .checked-indicator {
      width: 100%;
      height: 100%;
      display: block;
      fill: ${iA};
      opacity: 0;
      pointer-events: none;
    }

    .indeterminate-indicator {
      border-radius: calc(${tb} * 1px);
      background: ${iA};
      position: absolute;
      top: 50%;
      left: 50%;
      width: 50%;
      height: 50%;
      transform: translate(-50%, -50%);
      opacity: 0;
    }

    :host(:not([disabled])) .control:hover {
      background: ${i1};
      border-color: ${ob};
    }

    :host(:not([disabled])) .control:active {
      background: ${i2};
      border-color: ${of};
    }

    :host([aria-invalid='true']:not([disabled])) .control:hover {
      border-color: ${oT};
    }

    :host([aria-invalid='true']:not([disabled])) .control:active {
      border-color: ${oF};
    }

    :host(:${rm}) .control {
      outline: calc(${tk} * 1px) solid ${iE};
      outline-offset: 2px;
    }

    :host([aria-invalid='true']:${rm}) .control {
      outline-color: ${oS};
    }

    :host([aria-checked='true']) .control {
      background: ${iS};
      border: calc(${tw} * 1px) solid ${iS};
    }

    :host([aria-checked='true']:not([disabled])) .control:hover {
      background: ${iO};
      border: calc(${tw} * 1px) solid ${iO};
    }

    :host([aria-invalid='true'][aria-checked='true']) .control {
      background-color: ${oI};
      border-color: ${oI};
    }

    :host([aria-invalid='true'][aria-checked='true']:not([disabled]))
      .control:hover {
      background-color: ${oT};
      border-color: ${oT};
    }

    :host([aria-checked='true']:not([disabled]))
      .control:hover
      .checked-indicator {
      fill: ${iV};
    }

    :host([aria-checked='true']:not([disabled]))
      .control:hover
      .indeterminate-indicator {
      background: ${iV};
    }

    :host([aria-checked='true']:not([disabled])) .control:active {
      background: ${iD};
      border: calc(${tw} * 1px) solid ${iD};
    }

    :host([aria-invalid='true'][aria-checked='true']:not([disabled]))
      .control:active {
      background-color: ${oF};
      border-color: ${oF};
    }

    :host([aria-checked='true']:not([disabled]))
      .control:active
      .checked-indicator {
      fill: ${iP};
    }

    :host([aria-checked='true']:not([disabled]))
      .control:active
      .indeterminate-indicator {
      background: ${iP};
    }

    :host([aria-checked="true"]:${rm}:not([disabled])) .control {
      outline: calc(${tk} * 1px) solid ${iE};
      outline-offset: 2px;
    }

    :host([aria-invalid='true'][aria-checked="true"]:${rm}:not([disabled])) .control {
      outline-color: ${oS};
    }

    :host([disabled]) .label,
    :host([readonly]) .label,
    :host([readonly]) .control,
    :host([disabled]) .control {
      cursor: ${am};
    }

    :host([aria-checked='true']:not(.indeterminate)) .checked-indicator,
    :host(.indeterminate) .indeterminate-indicator {
      opacity: 1;
    }

    :host([disabled]) {
      opacity: ${t$};
    }
  `.withBehaviors(rv(rh`
      .control {
        forced-color-adjust: none;
        border-color: ${nX.FieldText};
        background: ${nX.Field};
      }
      :host([aria-invalid='true']) .control {
        border-style: dashed;
      }
      .checked-indicator {
        fill: ${nX.FieldText};
      }
      .indeterminate-indicator {
        background: ${nX.FieldText};
      }
      :host(:not([disabled])) .control:hover,
      .control:active {
        border-color: ${nX.Highlight};
        background: ${nX.Field};
      }
      :host(:${rm}) .control {
        outline: calc(${tk} * 1px) solid ${nX.FieldText};
        outline-offset: 2px;
      }
      :host([aria-checked="true"]:${rm}:not([disabled])) .control {
        outline: calc(${tk} * 1px) solid ${nX.FieldText};
        outline-offset: 2px;
      }
      :host([aria-checked='true']) .control {
        background: ${nX.Highlight};
        border-color: ${nX.Highlight};
      }
      :host([aria-checked='true']:not([disabled])) .control:hover,
      .control:active {
        border-color: ${nX.Highlight};
        background: ${nX.HighlightText};
      }
      :host([aria-checked='true']) .checked-indicator {
        fill: ${nX.HighlightText};
      }
      :host([aria-checked='true']:not([disabled]))
        .control:hover
        .checked-indicator {
        fill: ${nX.Highlight};
      }
      :host([aria-checked='true']) .indeterminate-indicator {
        background: ${nX.HighlightText};
      }
      :host([aria-checked='true']) .control:hover .indeterminate-indicator {
        background: ${nX.Highlight};
      }
      :host([disabled]) {
        opacity: 1;
      }
      :host([disabled]) .control {
        forced-color-adjust: none;
        border-color: ${nX.GrayText};
        background: ${nX.Field};
      }
      :host([disabled]) .indeterminate-indicator,
      :host([aria-checked='true'][disabled])
        .control:hover
        .indeterminate-indicator {
        forced-color-adjust: none;
        background: ${nX.GrayText};
      }
      :host([disabled]) .checked-indicator,
      :host([aria-checked='true'][disabled]) .control:hover .checked-indicator {
        forced-color-adjust: none;
        fill: ${nX.GrayText};
      }
    `)),aS=(e,t)=>s4`
  <template
    role="checkbox"
    aria-checked="${e=>e.checked}"
    aria-required="${e=>e.required}"
    aria-disabled="${e=>e.disabled}"
    aria-readonly="${e=>e.readOnly}"
    tabindex="${e=>e.disabled?null:0}"
    @keypress="${(e,t)=>e.keypressHandler(t.event)}"
    @click="${(e,t)=>e.clickHandler(t.event)}"
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
      <slot ${rn("defaultSlottedNodes")}></slot>
    </label>
  </template>
`;class aO extends aT{indeterminateChanged(e,t){this.indeterminate?this.classList.add("indeterminate"):this.classList.remove("indeterminate")}}let aD=aO.compose({baseName:"checkbox",baseClass:aT,template:aS,styles:aF,checkedIndicator:`
        <svg
            part="checked-indicator"
            class="checked-indicator"
            viewBox="0 0 20 20"
            xmlns="http://www.w3.org/2000/svg"
        >
            <path
                fill-rule="evenodd"
                clip-rule="evenodd"
                d="M8.143 12.6697L15.235 4.5L16.8 5.90363L8.23812 15.7667L3.80005 11.2556L5.27591 9.7555L8.143 12.6697Z"
            />
        </svg>
    `,indeterminateIndicator:`
        <div part="indeterminate-indicator" class="indeterminate-indicator"></div>
    `}),aE=0;function aR(e=""){return`${e}${aE++}`}function aL(e){return rg(e)&&("option"===e.getAttribute("role")||e instanceof HTMLOptionElement)}class aA extends sR{constructor(e,t,i,o){super(),this.defaultSelected=!1,this.dirtySelected=!1,this.selected=this.defaultSelected,this.dirtyValue=!1,e&&(this.textContent=e),t&&(this.initialValue=t),i&&(this.defaultSelected=i),o&&(this.selected=o),this.proxy=new Option(`${this.textContent}`,this.initialValue,this.defaultSelected,this.selected),this.proxy.disabled=this.disabled}checkedChanged(e,t){if("boolean"==typeof t){this.ariaChecked=t?"true":"false";return}this.ariaChecked=null}contentChanged(e,t){this.proxy instanceof HTMLOptionElement&&(this.proxy.textContent=this.textContent),this.$emit("contentchange",null,{bubbles:!0})}defaultSelectedChanged(){!this.dirtySelected&&(this.selected=this.defaultSelected,this.proxy instanceof HTMLOptionElement&&(this.proxy.selected=this.defaultSelected))}disabledChanged(e,t){this.ariaDisabled=this.disabled?"true":"false",this.proxy instanceof HTMLOptionElement&&(this.proxy.disabled=this.disabled)}selectedAttributeChanged(){this.defaultSelected=this.selectedAttribute,this.proxy instanceof HTMLOptionElement&&(this.proxy.defaultSelected=this.defaultSelected)}selectedChanged(){this.ariaSelected=this.selected?"true":"false",this.dirtySelected||(this.dirtySelected=!0),this.proxy instanceof HTMLOptionElement&&(this.proxy.selected=this.selected)}initialValueChanged(e,t){this.dirtyValue||(this.value=this.initialValue,this.dirtyValue=!1)}get label(){var e;return null!=(e=this.value)?e:this.text}get text(){var e,t;return null!=(t=null==(e=this.textContent)?void 0:e.replace(/\s+/g," ").trim())?t:""}set value(e){let t=`${null!=e?e:""}`;this._value=t,this.dirtyValue=!0,this.proxy instanceof HTMLOptionElement&&(this.proxy.value=t),ec.notify(this,"value")}get value(){var e;return ec.track(this,"value"),null!=(e=this._value)?e:this.text}get form(){return this.proxy?this.proxy.form:null}}X([eu],aA.prototype,"checked",void 0),X([eu],aA.prototype,"content",void 0),X([eu],aA.prototype,"defaultSelected",void 0),X([eS({mode:"boolean"})],aA.prototype,"disabled",void 0),X([eS({attribute:"selected",mode:"boolean"})],aA.prototype,"selectedAttribute",void 0),X([eu],aA.prototype,"selected",void 0),X([eS({attribute:"value",mode:"fromView"})],aA.prototype,"initialValue",void 0);class aV{}X([eu],aV.prototype,"ariaChecked",void 0),X([eu],aV.prototype,"ariaPosInSet",void 0),X([eu],aV.prototype,"ariaSelected",void 0),X([eu],aV.prototype,"ariaSetSize",void 0),rt(aV,rT),rt(aA,s8,aV);class aP extends sR{constructor(){super(...arguments),this._options=[],this.selectedIndex=-1,this.selectedOptions=[],this.shouldSkipFocus=!1,this.typeaheadBuffer="",this.typeaheadExpired=!0,this.typeaheadTimeout=-1}get firstSelectedOption(){var e;return null!=(e=this.selectedOptions[0])?e:null}get hasSelectableOptions(){return this.options.length>0&&!this.options.every(e=>e.disabled)}get length(){var e,t;return null!=(t=null==(e=this.options)?void 0:e.length)?t:0}get options(){return ec.track(this,"options"),this._options}set options(e){this._options=e,ec.notify(this,"options")}get typeAheadExpired(){return this.typeaheadExpired}set typeAheadExpired(e){this.typeaheadExpired=e}clickHandler(e){let t=e.target.closest("option,[role=option]");if(t&&!t.disabled)return this.selectedIndex=this.options.indexOf(t),!0}focusAndScrollOptionIntoView(e=this.firstSelectedOption){this.contains(document.activeElement)&&null!==e&&(e.focus(),requestAnimationFrame(()=>{e.scrollIntoView({block:"nearest"})}))}focusinHandler(e){this.shouldSkipFocus||e.target!==e.currentTarget||(this.setSelectedOptions(),this.focusAndScrollOptionIntoView()),this.shouldSkipFocus=!1}getTypeaheadMatches(){let e=this.typeaheadBuffer.replace(/[.*+\-?^${}()|[\]\\]/g,"\\$&"),t=RegExp(`^${e}`,"gi");return this.options.filter(e=>e.text.trim().match(t))}getSelectableIndex(e=this.selectedIndex,t){let i=e>t?-1:+(e<t),o=e+i,s=null;switch(i){case -1:s=this.options.reduceRight((e,t,i)=>e||t.disabled||!(i<o)?e:t,s);break;case 1:s=this.options.reduce((e,t,i)=>e||t.disabled||!(i>o)?e:t,s)}return this.options.indexOf(s)}handleChange(e,t){"selected"===t&&(aP.slottedOptionFilter(e)&&(this.selectedIndex=this.options.indexOf(e)),this.setSelectedOptions())}handleTypeAhead(e){this.typeaheadTimeout&&window.clearTimeout(this.typeaheadTimeout),this.typeaheadTimeout=window.setTimeout(()=>this.typeaheadExpired=!0,aP.TYPE_AHEAD_TIMEOUT_MS),e.length>1||(this.typeaheadBuffer=`${this.typeaheadExpired?"":this.typeaheadBuffer}${e}`)}keydownHandler(e){if(this.disabled)return!0;this.shouldSkipFocus=!1;let t=e.key;switch(t){case"Home":e.shiftKey||(e.preventDefault(),this.selectFirstOption());break;case oQ:e.shiftKey||(e.preventDefault(),this.selectNextOption());break;case o0:e.shiftKey||(e.preventDefault(),this.selectPreviousOption());break;case"End":e.preventDefault(),this.selectLastOption();break;case"Tab":return this.focusAndScrollOptionIntoView(),!0;case"Enter":case"Escape":return!0;case" ":if(this.typeaheadExpired)return!0;default:return 1===t.length&&this.handleTypeAhead(`${t}`),!0}}mousedownHandler(e){return this.shouldSkipFocus=!this.contains(document.activeElement),!0}multipleChanged(e,t){this.ariaMultiSelectable=t?"true":null}selectedIndexChanged(e,t){var i;if(!this.hasSelectableOptions){this.selectedIndex=-1;return}if((null==(i=this.options[this.selectedIndex])?void 0:i.disabled)&&"number"==typeof e){let i=this.getSelectableIndex(e,t),o=i>-1?i:e;this.selectedIndex=o,t===o&&this.selectedIndexChanged(t,o);return}this.setSelectedOptions()}selectedOptionsChanged(e,t){var i;let o=t.filter(aP.slottedOptionFilter);null==(i=this.options)||i.forEach(e=>{let t=ec.getNotifier(e);t.unsubscribe(this,"selected"),e.selected=o.includes(e),t.subscribe(this,"selected")})}selectFirstOption(){var e,t;this.disabled||(this.selectedIndex=null!=(t=null==(e=this.options)?void 0:e.findIndex(e=>!e.disabled))?t:-1)}selectLastOption(){this.disabled||(this.selectedIndex=function(e,t){let i=e.length;for(;i--;)if(t(e[i],i,e))return i;return -1}(this.options,e=>!e.disabled))}selectNextOption(){!this.disabled&&this.selectedIndex<this.options.length-1&&(this.selectedIndex+=1)}selectPreviousOption(){!this.disabled&&this.selectedIndex>0&&(this.selectedIndex=this.selectedIndex-1)}setDefaultSelectedOption(){var e,t;this.selectedIndex=null!=(t=null==(e=this.options)?void 0:e.findIndex(e=>e.defaultSelected))?t:-1}setSelectedOptions(){var e,t,i;(null==(e=this.options)?void 0:e.length)&&(this.selectedOptions=[this.options[this.selectedIndex]],this.ariaActiveDescendant=null!=(i=null==(t=this.firstSelectedOption)?void 0:t.id)?i:"",this.focusAndScrollOptionIntoView())}slottedOptionsChanged(e,t){this.options=t.reduce((e,t)=>(aL(t)&&e.push(t),e),[]);let i=`${this.options.length}`;this.options.forEach((e,t)=>{e.id||(e.id=aR("option-")),e.ariaPosInSet=`${t+1}`,e.ariaSetSize=i}),this.$fastController.isConnected&&(this.setSelectedOptions(),this.setDefaultSelectedOption())}typeaheadBufferChanged(e,t){if(this.$fastController.isConnected){let e=this.getTypeaheadMatches();if(e.length){let t=this.options.indexOf(e[0]);t>-1&&(this.selectedIndex=t)}this.typeaheadExpired=!1}}}aP.slottedOptionFilter=e=>aL(e)&&!e.hidden,aP.TYPE_AHEAD_TIMEOUT_MS=1e3,X([eS({mode:"boolean"})],aP.prototype,"disabled",void 0),X([eu],aP.prototype,"selectedIndex",void 0),X([eu],aP.prototype,"selectedOptions",void 0),X([eu],aP.prototype,"slottedOptions",void 0),X([eu],aP.prototype,"typeaheadBuffer",void 0);class aH{}X([eu],aH.prototype,"ariaActiveDescendant",void 0),X([eu],aH.prototype,"ariaDisabled",void 0),X([eu],aH.prototype,"ariaExpanded",void 0),X([eu],aH.prototype,"ariaMultiSelectable",void 0),rt(aH,rT),rt(aP,aH);let az="above",aM="below";class aN extends aP{}class aB extends ah(aN){constructor(){super(...arguments),this.proxy=document.createElement("input")}}class aj extends aB{constructor(){super(...arguments),this._value="",this.filteredOptions=[],this.filter="",this.forcedPosition=!1,this.listboxId=aR("listbox-"),this.maxHeight=0,this.open=!1}formResetCallback(){super.formResetCallback(),this.setDefaultSelectedOption(),this.updateValue()}validate(){super.validate(this.control)}get isAutocompleteInline(){return"inline"===this.autocomplete||this.isAutocompleteBoth}get isAutocompleteList(){return"list"===this.autocomplete||this.isAutocompleteBoth}get isAutocompleteBoth(){return"both"===this.autocomplete}openChanged(){if(this.open){this.ariaControls=this.listboxId,this.ariaExpanded="true",this.setPositioning(),this.focusAndScrollOptionIntoView(),el.queueUpdate(()=>this.focus());return}this.ariaControls="",this.ariaExpanded="false"}get options(){return ec.track(this,"options"),this.filteredOptions.length?this.filteredOptions:this._options}set options(e){this._options=e,ec.notify(this,"options")}placeholderChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.placeholder=this.placeholder)}positionChanged(e,t){this.positionAttribute=t,this.setPositioning()}get value(){return ec.track(this,"value"),this._value}set value(e){var t,i,o;let s=`${this._value}`;if(this.$fastController.isConnected&&this.options){let s=this.options.findIndex(t=>t.text.toLowerCase()===e.toLowerCase()),r=null==(t=this.options[this.selectedIndex])?void 0:t.text,a=null==(i=this.options[s])?void 0:i.text;this.selectedIndex=r!==a?s:this.selectedIndex,e=(null==(o=this.firstSelectedOption)?void 0:o.text)||e}s!==e&&(this._value=e,super.valueChanged(s,e),ec.notify(this,"value"))}clickHandler(e){if(!this.disabled){if(this.open){let t=e.target.closest("option,[role=option]");if(!t||t.disabled)return;this.selectedOptions=[t],this.control.value=t.text,this.clearSelectionRange(),this.updateValue(!0)}return this.open=!this.open,this.open&&this.control.focus(),!0}}connectedCallback(){super.connectedCallback(),this.forcedPosition=!!this.positionAttribute,this.value&&(this.initialValue=this.value)}disabledChanged(e,t){super.disabledChanged&&super.disabledChanged(e,t),this.ariaDisabled=this.disabled?"true":"false"}filterOptions(){this.autocomplete&&"none"!==this.autocomplete||(this.filter="");let e=this.filter.toLowerCase();this.filteredOptions=this._options.filter(e=>e.text.toLowerCase().startsWith(this.filter.toLowerCase())),this.isAutocompleteList&&(this.filteredOptions.length||e||(this.filteredOptions=this._options),this._options.forEach(e=>{e.hidden=!this.filteredOptions.includes(e)}))}focusAndScrollOptionIntoView(){this.contains(document.activeElement)&&(this.control.focus(),this.firstSelectedOption&&requestAnimationFrame(()=>{var e;null==(e=this.firstSelectedOption)||e.scrollIntoView({block:"nearest"})}))}focusoutHandler(e){if(this.syncValue(),!this.open)return!0;let t=e.relatedTarget;this.isSameNode(t)?this.focus():this.options&&this.options.includes(t)||(this.open=!1)}inputHandler(e){if(this.filter=this.control.value,this.filterOptions(),this.isAutocompleteInline||(this.selectedIndex=this.options.map(e=>e.text).indexOf(this.control.value)),e.inputType.includes("deleteContent")||!this.filter.length)return!0;this.isAutocompleteList&&!this.open&&(this.open=!0),this.isAutocompleteInline&&(this.filteredOptions.length?(this.selectedOptions=[this.filteredOptions[0]],this.selectedIndex=this.options.indexOf(this.firstSelectedOption),this.setInlineSelection()):this.selectedIndex=-1)}keydownHandler(e){let t=e.key;if(e.ctrlKey||e.shiftKey)return!0;switch(t){case"Enter":this.syncValue(),this.isAutocompleteInline&&(this.filter=this.value),this.open=!1,this.clearSelectionRange();break;case"Escape":if(this.isAutocompleteInline||(this.selectedIndex=-1),this.open){this.open=!1;break}this.value="",this.control.value="",this.filter="",this.filterOptions();break;case"Tab":if(this.setInputToSelection(),!this.open)return!0;e.preventDefault(),this.open=!1;break;case"ArrowUp":case"ArrowDown":if(this.filterOptions(),!this.open){this.open=!0;break}this.filteredOptions.length>0&&super.keydownHandler(e),this.isAutocompleteInline&&this.setInlineSelection();break;default:return!0}}keyupHandler(e){switch(e.key){case"ArrowLeft":case"ArrowRight":case"Backspace":case"Delete":case"Home":case"End":this.filter=this.control.value,this.selectedIndex=-1,this.filterOptions()}}selectedIndexChanged(e,t){if(this.$fastController.isConnected){if((t=o2(-1,this.options.length-1,t))!==this.selectedIndex){this.selectedIndex=t;return}super.selectedIndexChanged(e,t)}}selectPreviousOption(){!this.disabled&&this.selectedIndex>=0&&(this.selectedIndex=this.selectedIndex-1)}setDefaultSelectedOption(){if(this.$fastController.isConnected&&this.options){let e=this.options.findIndex(e=>null!==e.getAttribute("selected")||e.selected);this.selectedIndex=e,!this.dirtyValue&&this.firstSelectedOption&&(this.value=this.firstSelectedOption.text),this.setSelectedOptions()}}setInputToSelection(){this.firstSelectedOption&&(this.control.value=this.firstSelectedOption.text,this.control.focus())}setInlineSelection(){this.firstSelectedOption&&(this.setInputToSelection(),this.control.setSelectionRange(this.filter.length,this.control.value.length,"backward"))}syncValue(){var e;let t=this.selectedIndex>-1?null==(e=this.firstSelectedOption)?void 0:e.text:this.control.value;this.updateValue(this.value!==t)}setPositioning(){let e=this.getBoundingClientRect(),t=window.innerHeight-e.bottom;this.position=this.forcedPosition?this.positionAttribute:e.top>t?az:aM,this.positionAttribute=this.forcedPosition?this.positionAttribute:this.position,this.maxHeight=this.position===az?~~e.top:~~t}selectedOptionsChanged(e,t){this.$fastController.isConnected&&this._options.forEach(e=>{e.selected=t.includes(e)})}slottedOptionsChanged(e,t){super.slottedOptionsChanged(e,t),this.updateValue()}updateValue(e){var t;this.$fastController.isConnected&&(this.value=(null==(t=this.firstSelectedOption)?void 0:t.text)||this.control.value,this.control.value=this.value),e&&this.$emit("change")}clearSelectionRange(){let e=this.control.value.length;this.control.setSelectionRange(e,e)}}X([eS({attribute:"autocomplete",mode:"fromView"})],aj.prototype,"autocomplete",void 0),X([eu],aj.prototype,"maxHeight",void 0),X([eS({attribute:"open",mode:"boolean"})],aj.prototype,"open",void 0),X([eS],aj.prototype,"placeholder",void 0),X([eS({attribute:"position"})],aj.prototype,"positionAttribute",void 0),X([eu],aj.prototype,"position",void 0);class aq{}X([eu],aq.prototype,"ariaAutoComplete",void 0),X([eu],aq.prototype,"ariaControls",void 0),rt(aq,aH),rt(aj,s8,aq);class aU extends aP{constructor(){super(...arguments),this.activeIndex=-1,this.rangeStartIndex=-1}get activeOption(){return this.options[this.activeIndex]}get checkedOptions(){var e;return null==(e=this.options)?void 0:e.filter(e=>e.checked)}get firstSelectedOptionIndex(){return this.options.indexOf(this.firstSelectedOption)}activeIndexChanged(e,t){var i,o;this.ariaActiveDescendant=null!=(o=null==(i=this.options[t])?void 0:i.id)?o:"",this.focusAndScrollOptionIntoView()}checkActiveIndex(){if(!this.multiple)return;let e=this.activeOption;e&&(e.checked=!0)}checkFirstOption(e=!1){e?(-1===this.rangeStartIndex&&(this.rangeStartIndex=this.activeIndex+1),this.options.forEach((e,t)=>{e.checked=o5(t,this.rangeStartIndex)})):this.uncheckAllOptions(),this.activeIndex=0,this.checkActiveIndex()}checkLastOption(e=!1){e?(-1===this.rangeStartIndex&&(this.rangeStartIndex=this.activeIndex),this.options.forEach((e,t)=>{e.checked=o5(t,this.rangeStartIndex,this.options.length)})):this.uncheckAllOptions(),this.activeIndex=this.options.length-1,this.checkActiveIndex()}connectedCallback(){super.connectedCallback(),this.addEventListener("focusout",this.focusoutHandler)}disconnectedCallback(){this.removeEventListener("focusout",this.focusoutHandler),super.disconnectedCallback()}checkNextOption(e=!1){e?(-1===this.rangeStartIndex&&(this.rangeStartIndex=this.activeIndex),this.options.forEach((e,t)=>{e.checked=o5(t,this.rangeStartIndex,this.activeIndex+1)})):this.uncheckAllOptions(),this.activeIndex+=+(this.activeIndex<this.options.length-1),this.checkActiveIndex()}checkPreviousOption(e=!1){e?(-1===this.rangeStartIndex&&(this.rangeStartIndex=this.activeIndex),1===this.checkedOptions.length&&(this.rangeStartIndex+=1),this.options.forEach((e,t)=>{e.checked=o5(t,this.activeIndex,this.rangeStartIndex)})):this.uncheckAllOptions(),this.activeIndex-=this.activeIndex>0,this.checkActiveIndex()}clickHandler(e){var t;if(!this.multiple)return super.clickHandler(e);let i=null==(t=e.target)?void 0:t.closest("[role=option]");if(i&&!i.disabled)return this.uncheckAllOptions(),this.activeIndex=this.options.indexOf(i),this.checkActiveIndex(),this.toggleSelectedForAllCheckedOptions(),!0}focusAndScrollOptionIntoView(){super.focusAndScrollOptionIntoView(this.activeOption)}focusinHandler(e){if(!this.multiple)return super.focusinHandler(e);this.shouldSkipFocus||e.target!==e.currentTarget||(this.uncheckAllOptions(),-1===this.activeIndex&&(this.activeIndex=-1!==this.firstSelectedOptionIndex?this.firstSelectedOptionIndex:0),this.checkActiveIndex(),this.setSelectedOptions(),this.focusAndScrollOptionIntoView()),this.shouldSkipFocus=!1}focusoutHandler(e){this.multiple&&this.uncheckAllOptions()}keydownHandler(e){if(!this.multiple)return super.keydownHandler(e);if(this.disabled)return!0;let{key:t,shiftKey:i}=e;switch(this.shouldSkipFocus=!1,t){case"Home":return void this.checkFirstOption(i);case oQ:return void this.checkNextOption(i);case o0:return void this.checkPreviousOption(i);case"End":return void this.checkLastOption(i);case"Tab":return this.focusAndScrollOptionIntoView(),!0;case"Escape":return this.uncheckAllOptions(),this.checkActiveIndex(),!0;case" ":if(e.preventDefault(),this.typeAheadExpired)return void this.toggleSelectedForAllCheckedOptions();default:return 1===t.length&&this.handleTypeAhead(`${t}`),!0}}mousedownHandler(e){if(e.offsetX>=0&&e.offsetX<=this.scrollWidth)return super.mousedownHandler(e)}multipleChanged(e,t){var i;this.ariaMultiSelectable=t?"true":null,null==(i=this.options)||i.forEach(e=>{e.checked=!t&&void 0}),this.setSelectedOptions()}setSelectedOptions(){this.multiple?this.$fastController.isConnected&&this.options&&(this.selectedOptions=this.options.filter(e=>e.selected),this.focusAndScrollOptionIntoView()):super.setSelectedOptions()}sizeChanged(e,t){var i;let o=Math.max(0,parseInt(null!=(i=null==t?void 0:t.toFixed())?i:"",10));o!==t&&el.queueUpdate(()=>{this.size=o})}toggleSelectedForAllCheckedOptions(){let e=this.checkedOptions.filter(e=>!e.disabled),t=!e.every(e=>e.selected);e.forEach(e=>e.selected=t),this.selectedIndex=this.options.indexOf(e[e.length-1]),this.setSelectedOptions()}typeaheadBufferChanged(e,t){if(!this.multiple)return void super.typeaheadBufferChanged(e,t);if(this.$fastController.isConnected){let e=this.getTypeaheadMatches(),t=this.options.indexOf(e[0]);t>-1&&(this.activeIndex=t,this.uncheckAllOptions(),this.checkActiveIndex()),this.typeAheadExpired=!1}}uncheckAllOptions(e=!1){this.options.forEach(e=>e.checked=!this.multiple&&void 0),e||(this.rangeStartIndex=-1)}}X([eu],aU.prototype,"activeIndex",void 0),X([eS({mode:"boolean"})],aU.prototype,"multiple",void 0),X([eS({converter:eT})],aU.prototype,"size",void 0);class a_ extends aU{}class aG extends ah(a_){constructor(){super(...arguments),this.proxy=document.createElement("select")}}class aW extends aG{constructor(){super(...arguments),this.open=!1,this.forcedPosition=!1,this.listboxId=aR("listbox-"),this.maxHeight=0}openChanged(e,t){if(this.collapsible){if(this.open){this.ariaControls=this.listboxId,this.ariaExpanded="true",this.setPositioning(),this.focusAndScrollOptionIntoView(),this.indexWhenOpened=this.selectedIndex,el.queueUpdate(()=>this.focus());return}this.ariaControls="",this.ariaExpanded="false"}}get collapsible(){return!(this.multiple||"number"==typeof this.size)}get value(){return ec.track(this,"value"),this._value}set value(e){var t,i,o,s,r,a,n;let l=`${this._value}`;if(null==(t=this._options)?void 0:t.length){let t=this._options.findIndex(t=>t.value===e),l=null!=(o=null==(i=this._options[this.selectedIndex])?void 0:i.value)?o:null,h=null!=(r=null==(s=this._options[t])?void 0:s.value)?r:null;(-1===t||l!==h)&&(e="",this.selectedIndex=t),e=null!=(n=null==(a=this.firstSelectedOption)?void 0:a.value)?n:e}l!==e&&(this._value=e,super.valueChanged(l,e),ec.notify(this,"value"),this.updateDisplayValue())}updateValue(e){var t,i;this.$fastController.isConnected&&(this.value=null!=(i=null==(t=this.firstSelectedOption)?void 0:t.value)?i:""),e&&(this.$emit("input"),this.$emit("change",this,{bubbles:!0,composed:void 0}))}selectedIndexChanged(e,t){super.selectedIndexChanged(e,t),this.updateValue()}positionChanged(e,t){this.positionAttribute=t,this.setPositioning()}setPositioning(){let e=this.getBoundingClientRect(),t=window.innerHeight-e.bottom;this.position=this.forcedPosition?this.positionAttribute:e.top>t?az:aM,this.positionAttribute=this.forcedPosition?this.positionAttribute:this.position,this.maxHeight=this.position===az?~~e.top:~~t}get displayValue(){var e,t;return ec.track(this,"displayValue"),null!=(t=null==(e=this.firstSelectedOption)?void 0:e.text)?t:""}disabledChanged(e,t){super.disabledChanged&&super.disabledChanged(e,t),this.ariaDisabled=this.disabled?"true":"false"}formResetCallback(){this.setProxyOptions(),super.setDefaultSelectedOption(),-1===this.selectedIndex&&(this.selectedIndex=0)}clickHandler(e){if(!this.disabled){if(this.open){let t=e.target.closest("option,[role=option]");if(t&&t.disabled)return}return super.clickHandler(e),this.open=this.collapsible&&!this.open,this.open||this.indexWhenOpened===this.selectedIndex||this.updateValue(!0),!0}}focusoutHandler(e){var t;if(super.focusoutHandler(e),!this.open)return!0;let i=e.relatedTarget;this.isSameNode(i)?this.focus():(null==(t=this.options)?void 0:t.includes(i))||(this.open=!1,this.indexWhenOpened!==this.selectedIndex&&this.updateValue(!0))}handleChange(e,t){super.handleChange(e,t),"value"===t&&this.updateValue()}slottedOptionsChanged(e,t){this.options.forEach(e=>{ec.getNotifier(e).unsubscribe(this,"value")}),super.slottedOptionsChanged(e,t),this.options.forEach(e=>{ec.getNotifier(e).subscribe(this,"value")}),this.setProxyOptions(),this.updateValue()}mousedownHandler(e){var t;return e.offsetX>=0&&e.offsetX<=(null==(t=this.listbox)?void 0:t.scrollWidth)?super.mousedownHandler(e):this.collapsible}multipleChanged(e,t){super.multipleChanged(e,t),this.proxy&&(this.proxy.multiple=t)}selectedOptionsChanged(e,t){var i;super.selectedOptionsChanged(e,t),null==(i=this.options)||i.forEach((e,t)=>{var i;let o=null==(i=this.proxy)?void 0:i.options.item(t);o&&(o.selected=e.selected)})}setDefaultSelectedOption(){var e;let t=null!=(e=this.options)?e:Array.from(this.children).filter(aP.slottedOptionFilter),i=null==t?void 0:t.findIndex(e=>e.hasAttribute("selected")||e.selected||e.value===this.value);if(-1!==i){this.selectedIndex=i;return}this.selectedIndex=0}setProxyOptions(){this.proxy instanceof HTMLSelectElement&&this.options&&(this.proxy.options.length=0,this.options.forEach(e=>{let t=e.proxy||(e instanceof HTMLOptionElement?e.cloneNode():null);t&&this.proxy.options.add(t)}))}keydownHandler(e){super.keydownHandler(e);let t=e.key||e.key.charCodeAt(0);switch(t){case" ":e.preventDefault(),this.collapsible&&this.typeAheadExpired&&(this.open=!this.open);break;case"Home":case"End":e.preventDefault();break;case"Enter":e.preventDefault(),this.open=!this.open;break;case"Escape":this.collapsible&&this.open&&(e.preventDefault(),this.open=!1);break;case"Tab":return this.collapsible&&this.open&&(e.preventDefault(),this.open=!1),!0}return this.open||this.indexWhenOpened===this.selectedIndex||(this.updateValue(!0),this.indexWhenOpened=this.selectedIndex),t!==oQ&&t!==o0}connectedCallback(){super.connectedCallback(),this.forcedPosition=!!this.positionAttribute,this.addEventListener("contentchange",this.updateDisplayValue)}disconnectedCallback(){this.removeEventListener("contentchange",this.updateDisplayValue),super.disconnectedCallback()}sizeChanged(e,t){super.sizeChanged(e,t),this.proxy&&(this.proxy.size=t)}updateDisplayValue(){this.collapsible&&ec.notify(this,"displayValue")}}X([eS({attribute:"open",mode:"boolean"})],aW.prototype,"open",void 0),X([function(e,t,i){return Object.assign({},i,{get:function(){return ec.trackVolatile(),i.get.apply(this)}})}],aW.prototype,"collapsible",null),X([eu],aW.prototype,"control",void 0),X([eS({attribute:"position"})],aW.prototype,"positionAttribute",void 0),X([eu],aW.prototype,"position",void 0),X([eu],aW.prototype,"maxHeight",void 0);class aK{}X([eu],aK.prototype,"ariaControls",void 0),rt(aK,aH),rt(aW,s8,aK);let aX=(e,t)=>{let i=e.tagFor(aA),o=e.name===e.tagFor(aU)?"":".listbox";return rh`
        ${!o?ru("inline-flex"):""}

        :host ${o} {
            background: ${iT};
            border: calc(${tw} * 1px) solid ${om};
            border-radius: calc(${tb} * 1px);
            box-sizing: border-box;
            flex-direction: column;
            padding: calc(${tv} * 1px) 0;
        }

        ${!o?rh`
:host(:${rm}:not([disabled])) {
                outline: none;
            }

            :host(:focus-within:not([disabled])) {
                border-color: ${on};
                box-shadow: 0 0 0
                    calc((${tk} - ${tw}) * 1px)
                    ${on} inset;
            }

            :host([disabled]) ::slotted(*) {
                cursor: ${am};
                opacity: ${t$};
                pointer-events: none;
            }
        `:""}

        ${o||":host([size])"} {
            max-height: calc(
                (var(--size) * ${ry} + (${tv} * ${tw} * 2)) * 1px
            );
            overflow-y: auto;
        }

        :host([size="0"]) ${o} {
            max-height: none;
        }
    `.withBehaviors(rv(rh`
                :host(:not([multiple]):${rm}) ::slotted(${i}[aria-selected="true"]),
                :host([multiple]:${rm}) ::slotted(${i}[aria-checked="true"]) {
                    border-color: ${nX.ButtonText};
                    box-shadow: 0 0 0 calc(${tk} * 1px) inset ${nX.HighlightText};
                }

                :host(:not([multiple]):${rm}) ::slotted(${i}[aria-selected="true"]) {
                    background: ${nX.Highlight};
                    color: ${nX.HighlightText};
                    fill: currentcolor;
                }

                ::slotted(${i}[aria-selected="true"]:not([aria-checked="true"])) {
                    background: ${nX.Highlight};
                    border-color: ${nX.HighlightText};
                    color: ${nX.HighlightText};
                }
            `))},aY=(e,t)=>{let i=e.name===e.tagFor(aW);return rh`
  ${ru("inline-flex")}
  
  :host {
    --elevation: 14;
    background: ${i0};
    border-radius: calc(${tb} * 1px);
    border: calc(${tw} * 1px) solid ${oe};
    box-sizing: border-box;
    color: ${op};
    font-family: ${tu};
    height: calc(${ry} * 1px);
    position: relative;
    user-select: none;
    min-width: 250px;
    outline: none;
    vertical-align: top;
  }

  :host([aria-invalid='true']) {
    border-color: ${oI};
  }
  
  :host(:not([autowidth])) {
    min-width: 250px;
  }
  
  ${i?rh`
  :host(:not([aria-haspopup])) {
    --elevation: 0;
    border: 0;
    height: auto;
    min-width: 0;
  }
  `:""}
  
  ${aX(e,t)}
  
  :host .listbox {
    ${ax}
    border: none;
    display: flex;
    left: 0;
    position: absolute;
    width: 100%;
    z-index: 1;
  }
  
  .control + .listbox {
    --stroke-size: calc(${tv} * ${tw} * 2);
    max-height: calc(
      (var(--listbox-max-height) * ${ry} + var(--stroke-size)) * 1px
      );
  }
  
  ${i?rh`
  :host(:not([aria-haspopup])) .listbox {
    left: auto;
    position: static;
    z-index: auto;
  }
  `:""}
  
  :host(:not([autowidth])) .listbox {
    width: 100%;
  }
  
  :host([autowidth]) ::slotted([role='option']),
  :host([autowidth]) ::slotted(option) {
    padding: 0 calc(1em + ${tv} * 1.25px + 1px);
  }
  
  .listbox[hidden] {
    display: none;
  }
  
  .control {
    align-items: center;
    box-sizing: border-box;
    cursor: pointer;
    display: flex;
    font-size: ${tC};
    font-family: inherit;
    line-height: ${tI};
    min-height: 100%;
    padding: 0 calc(${tv} * 2.25px);
    width: 100%;
  }

  :host([minimal]),
  :host([scale='xsmall']) {
    --element-scale: -4;
  }

  :host([scale='small']) {
    --element-scale: -2;
  }

  :host([scale='medium']) {
    --element-scale: 0;
  }

  :host([scale='large']) {
    --element-scale: 2;
  }

  :host([scale='xlarge']) {
    --element-scale: 4;
  }
  
  :host(:not([disabled]):hover) {
    background: ${i1};
    border-color: ${ot};
  }
  
  :host([aria-invalid='true']:not([disabled]):hover) {
    border-color: ${oT};
  }
  
  :host(:${rm}) {
    border-color: ${iE};
    box-shadow: 0 0 0 calc((${tk} - ${tw}) * 1px)
    ${iE};
  }
  
  :host([aria-invalid='true']:${rm}) {
    border-color: ${oS};
    box-shadow: 0 0 0 calc((${tk} - ${tw}) * 1px)
    ${oS};
  }
  
  :host(:not([size]):not([multiple]):not([open]):${rm}),
  :host([multiple]:${rm}),
  :host([size]:${rm}) {
    box-shadow: 0 0 0 calc((${tk} - ${tw}) * 1px)
    ${iE};
  }
  
  :host([aria-invalid='true']:not([size]):not([multiple]):not([open]):${rm}),
  :host([aria-invalid='true'][multiple]:${rm}),
  :host([aria-invalid='true'][size]:${rm}) {
    box-shadow: 0 0 0 calc((${tk} - ${tw}) * 1px)
    ${oS};
  }
  
  :host(:not([multiple]):not([size]):${rm}) ::slotted(${e.tagFor(aA)}[aria-selected="true"]:not([disabled])) {
    box-shadow: 0 0 0 calc(${tk} * 1px) inset ${iE};
    border-color: ${iE};
    background: ${iE};
    color: ${iH};
  }
    
  :host([disabled]) {
    cursor: ${am};
    opacity: ${t$};
  }
  
  :host([disabled]) .control {
    cursor: ${am};
    user-select: none;
  }
  
  :host([disabled]:hover) {
    background: ${i4};
    color: ${op};
    fill: currentcolor;
  }
  
  :host(:not([disabled])) .control:active {
    background: ${i2};
    border-color: ${iD};
    border-radius: calc(${tb} * 1px);
  }
  
  :host([open][position="above"]) .listbox {
    border-bottom-left-radius: 0;
    border-bottom-right-radius: 0;
    border-bottom: 0;
    bottom: calc(${ry} * 1px);
  }
  
  :host([open][position="below"]) .listbox {
    border-top-left-radius: 0;
    border-top-right-radius: 0;
    border-top: 0;
    top: calc(${ry} * 1px);
  }
  
  .selected-value {
    flex: 1 1 auto;
    font-family: inherit;
    overflow: hidden;
    text-align: start;
    text-overflow: ellipsis;
    white-space: nowrap;
  }
  
  .indicator {
    flex: 0 0 auto;
    margin-inline-start: 1em;
  }
  
  slot[name="listbox"] {
    display: none;
    width: 100%;
  }
  
  :host([open]) slot[name="listbox"] {
    display: flex;
    position: absolute;
    ${ax}
  }
  
  .end {
    margin-inline-start: auto;
  }
  
  .start,
  .end,
  .indicator,
  .select-indicator,
  ::slotted(svg) {
    /* TODO: adaptive typography https://github.com/microsoft/fast/issues/2432 */
    fill: currentcolor;
    height: 1em;
    min-height: calc(${tv} * 4px);
    min-width: calc(${tv} * 4px);
    width: 1em;
  }
  
  ::slotted([role="option"]),
  ::slotted(option) {
    flex: 0 0 auto;
  }
  `.withBehaviors(rv(rh`
      :host(:not([disabled]):hover),
      :host(:not([disabled]):active) {
        border-color: ${nX.Highlight};
      }

      :host([aria-invalid='true']) {
        border-style: dashed;
      }
      
      :host(:not([disabled]):${rm}) {
        background-color: ${nX.ButtonFace};
        box-shadow: 0 0 0 calc(${tk} * 1px) ${nX.Highlight};
        color: ${nX.ButtonText};
        fill: currentcolor;
        forced-color-adjust: none;
      }
      
      :host(:not([disabled]):${rm}) .listbox {
        background: ${nX.ButtonFace};
      }
      
      :host([disabled]) {
        border-color: ${nX.GrayText};
        background-color: ${nX.ButtonFace};
        color: ${nX.GrayText};
        fill: currentcolor;
        opacity: 1;
        forced-color-adjust: none;
      }
      
      :host([disabled]:hover) {
        background: ${nX.ButtonFace};
      }
      
      :host([disabled]) .control {
        color: ${nX.GrayText};
        border-color: ${nX.GrayText};
      }
      
      :host([disabled]) .control .select-indicator {
        fill: ${nX.GrayText};
      }
      
      :host(:${rm}) ::slotted([aria-selected="true"][role="option"]),
      :host(:${rm}) ::slotted(option[aria-selected="true"]),
      :host(:${rm}) ::slotted([aria-selected="true"][role="option"]:not([disabled])) {
        background: ${nX.Highlight};
        border-color: ${nX.ButtonText};
        box-shadow: 0 0 0 calc((${tk} - ${tw}) * 1px)
        ${nX.HighlightText};
        color: ${nX.HighlightText};
        fill: currentcolor;
      }
      
      .start,
      .end,
      .indicator,
      .select-indicator,
      ::slotted(svg) {
        color: ${nX.ButtonText};
        fill: currentcolor;
      }
      `))},aQ=(e,t)=>rh`
  ${aY(e,t)}

  :host(:empty) .listbox {
    display: none;
  }

  :host([disabled]) *,
  :host([disabled]) {
    cursor: ${am};
    user-select: none;
  }

  :host(:focus-within:not([disabled])) {
    border-color: ${iE};
    box-shadow: 0 0 0 calc((${tk} - ${tw}) * 1px)
      ${iE};
  }

  :host([aria-invalid='true']:focus-within:not([disabled])) {
    border-color: ${oS};
    box-shadow: 0 0 0 calc((${tk} - ${tw}) * 1px)
      ${oS};
  }

  .selected-value {
    -webkit-appearance: none;
    background: transparent;
    border: none;
    color: inherit;
    font-size: ${tC};
    line-height: ${tI};
    height: calc(100% - (${tw} * 1px));
    margin: auto 0;
    width: 100%;
  }

  .selected-value:hover,
    .selected-value:${rm},
    .selected-value:disabled,
    .selected-value:active {
    outline: none;
  }
`;class aZ extends aj{connectedCallback(){super.connectedCallback(),this.setAutoWidth()}slottedOptionsChanged(e,t){super.slottedOptionsChanged(e,t),this.setAutoWidth()}autoWidthChanged(e,t){t?this.setAutoWidth():this.style.removeProperty("width")}setAutoWidth(){if(!this.autoWidth||!this.isConnected)return;let e=this.listbox.getBoundingClientRect().width;0===e&&this.listbox.hidden&&(Object.assign(this.listbox.style,{visibility:"hidden"}),this.listbox.removeAttribute("hidden"),e=this.listbox.getBoundingClientRect().width,this.listbox.setAttribute("hidden",""),this.listbox.style.removeProperty("visibility")),e>0&&Object.assign(this.style,{width:`${e}px`})}maxHeightChanged(e,t){this.updateComputedStylesheet()}updateComputedStylesheet(){this.computedStylesheet&&this.$fastController.removeStyles(this.computedStylesheet);let e=Math.floor(this.maxHeight/o$.getValueFor(this)).toString();this.computedStylesheet=rh`
      :host {
        --listbox-max-height: ${e};
      }
    `,this.$fastController.addStyles(this.computedStylesheet)}}rI([eS({attribute:"autowidth",mode:"boolean"})],aZ.prototype,"autoWidth",void 0),rI([eS({attribute:"minimal",mode:"boolean"})],aZ.prototype,"minimal",void 0),rI([eS],aZ.prototype,"scale",void 0);let aJ=aZ.compose({baseName:"combobox",baseClass:aj,template:(e,t)=>s4`
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
            ${re(e,t)}
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
                    ${s9("control")}
                />
                <div class="indicator" part="indicator" aria-hidden="true">
                    <slot name="indicator">
                        ${t.indicator||""}
                    </slot>
                </div>
            </slot>
            ${s7(e,t)}
        </div>
        <div
            class="listbox"
            id="${e=>e.listboxId}"
            part="listbox"
            role="listbox"
            ?disabled="${e=>e.disabled}"
            ?hidden="${e=>!e.open}"
            ${s9("listbox")}
        >
            <slot
                ${rn({filter:aP.slottedOptionFilter,flatten:!0,property:"slottedOptions"})}
            ></slot>
        </div>
    </template>
`,styles:aQ,shadowOptions:{delegatesFocus:!0},indicator:`
        <svg
            class="select-indicator"
            part="select-indicator"
            viewBox="0 0 12 7"
            xmlns="http://www.w3.org/2000/svg"
        >
            <path
                d="M11.85.65c.2.2.2.5 0 .7L6.4 6.84a.55.55 0 01-.78 0L.14 1.35a.5.5 0 11.71-.7L6 5.8 11.15.65c.2-.2.5-.2.7 0z"
            />
        </svg>
    `}),a0="focusin",a1="focusout",a2="keydown",a5="sticky",a3="default",a4="columnheader",a6="default",a9=s4`
    <template>
        ${e=>null===e.rowData||null===e.columnDefinition||null===e.columnDefinition.columnDataKey?null:e.rowData[e.columnDefinition.columnDataKey]}
    </template>
`,a8=s4`
    <template>
        ${e=>null===e.columnDefinition?null:void 0===e.columnDefinition.title?e.columnDefinition.columnDataKey:e.columnDefinition.title}
    </template>
`;class a7 extends sR{constructor(){super(...arguments),this.cellType=a3,this.rowData=null,this.columnDefinition=null,this.isActiveCell=!1,this.customCellView=null,this.updateCellStyle=()=>{this.style.gridColumn=this.gridColumn}}cellTypeChanged(){this.$fastController.isConnected&&this.updateCellView()}gridColumnChanged(){this.$fastController.isConnected&&this.updateCellStyle()}columnDefinitionChanged(e,t){this.$fastController.isConnected&&this.updateCellView()}connectedCallback(){var e;super.connectedCallback(),this.addEventListener(a0,this.handleFocusin),this.addEventListener(a1,this.handleFocusout),this.addEventListener(a2,this.handleKeydown),this.style.gridColumn=`${(null==(e=this.columnDefinition)?void 0:e.gridColumn)===void 0?0:this.columnDefinition.gridColumn}`,this.updateCellView(),this.updateCellStyle()}disconnectedCallback(){super.disconnectedCallback(),this.removeEventListener(a0,this.handleFocusin),this.removeEventListener(a1,this.handleFocusout),this.removeEventListener(a2,this.handleKeydown),this.disconnectCellView()}handleFocusin(e){if(!this.isActiveCell){if(this.isActiveCell=!0,this.cellType===a4){if(null!==this.columnDefinition&&!0!==this.columnDefinition.headerCellInternalFocusQueue&&"function"==typeof this.columnDefinition.headerCellFocusTargetCallback){let e=this.columnDefinition.headerCellFocusTargetCallback(this);null!==e&&e.focus()}}else if(null!==this.columnDefinition&&!0!==this.columnDefinition.cellInternalFocusQueue&&"function"==typeof this.columnDefinition.cellFocusTargetCallback){let e=this.columnDefinition.cellFocusTargetCallback(this);null!==e&&e.focus()}this.$emit("cell-focused",this)}}handleFocusout(e){this===document.activeElement||this.contains(document.activeElement)||(this.isActiveCell=!1)}handleKeydown(e){if(!e.defaultPrevented&&null!==this.columnDefinition&&(this.cellType!==a3||!0===this.columnDefinition.cellInternalFocusQueue)&&(this.cellType!==a4||!0===this.columnDefinition.headerCellInternalFocusQueue))switch(e.key){case"Enter":case"F2":if(this.contains(document.activeElement)&&document.activeElement!==this)return;if(this.cellType===a4){if(void 0!==this.columnDefinition.headerCellFocusTargetCallback){let t=this.columnDefinition.headerCellFocusTargetCallback(this);null!==t&&t.focus(),e.preventDefault()}}else if(void 0!==this.columnDefinition.cellFocusTargetCallback){let t=this.columnDefinition.cellFocusTargetCallback(this);null!==t&&t.focus(),e.preventDefault()}break;case"Escape":this.contains(document.activeElement)&&document.activeElement!==this&&(this.focus(),e.preventDefault())}}updateCellView(){if(this.disconnectCellView(),null!==this.columnDefinition)switch(this.cellType){case a4:void 0!==this.columnDefinition.headerCellTemplate?this.customCellView=this.columnDefinition.headerCellTemplate.render(this,this):this.customCellView=a8.render(this,this);break;case void 0:case"rowheader":case a3:void 0!==this.columnDefinition.cellTemplate?this.customCellView=this.columnDefinition.cellTemplate.render(this,this):this.customCellView=a9.render(this,this)}}disconnectCellView(){null!==this.customCellView&&(this.customCellView.dispose(),this.customCellView=null)}}function ne(e,t,i){return{index:e,removed:t,addedCount:i}}function nt(e,t,i,o,s,r){let a,n=0,l=0,h=Math.min(i-t,r-s);if(0===t&&0===s&&(n=function(e,t,i){for(let o=0;o<i;++o)if(e[o]!==t[o])return o;return i}(e,o,h)),i===e.length&&r===o.length&&(l=function(e,t,i){let o=e.length,s=t.length,r=0;for(;r<i&&e[--o]===t[--s];)r++;return r}(e,o,h-n)),t+=n,s+=n,i-=l,r-=l,i-t==0&&r-s==0)return ee;if(t===i){let e=ne(t,[],0);for(;s<r;)e.removed.push(o[s++]);return[e]}if(s===r)return[ne(t,[],i-t)];let d=function(e){let t=e.length-1,i=e[0].length-1,o=e[t][i],s=[];for(;t>0||i>0;){let r;if(0===t){s.push(2),i--;continue}if(0===i){s.push(3),t--;continue}let a=e[t-1][i-1],n=e[t-1][i],l=e[t][i-1];(r=n<l?n<a?n:a:l<a?l:a)===a?(a===o?s.push(0):(s.push(1),o=a),t--,i--):r===n?(s.push(3),t--,o=n):(s.push(2),i--,o=l)}return s.reverse(),s}(function(e,t,i,o,s,r){let a,n,l=r-s+1,h=i-t+1,d=Array(l);for(let e=0;e<l;++e)d[e]=Array(h),d[e][0]=e;for(let e=0;e<h;++e)d[0][e]=e;for(let i=1;i<l;++i)for(let r=1;r<h;++r)e[t+r-1]===o[s+i-1]?d[i][r]=d[i-1][r-1]:(a=d[i-1][r]+1,n=d[i][r-1]+1,d[i][r]=a<n?a:n);return d}(e,t,i,o,s,r)),c=[],u=t,p=s;for(let e=0;e<d.length;++e)switch(d[e]){case 0:void 0!==a&&(c.push(a),a=void 0),u++,p++;break;case 1:void 0===a&&(a=ne(u,[],0)),a.addedCount++,u++,a.removed.push(o[p]),p++;break;case 2:void 0===a&&(a=ne(u,[],0)),a.addedCount++,u++;break;case 3:void 0===a&&(a=ne(u,[],0)),a.removed.push(o[p]),p++}return void 0!==a&&c.push(a),c}X([eS({attribute:"cell-type"})],a7.prototype,"cellType",void 0),X([eS({attribute:"grid-column"})],a7.prototype,"gridColumn",void 0),X([eu],a7.prototype,"rowData",void 0),X([eu],a7.prototype,"columnDefinition",void 0);let ni=Array.prototype.push,no=!1;function ns(e,t){let i=e.index,o=t.length;return i>o?i=o-e.addedCount:i<0&&(i=o+e.removed.length+i-e.addedCount),i<0&&(i=0),e.index=i,e}class nr extends eh{constructor(e){super(e),this.oldCollection=void 0,this.splices=void 0,this.needsQueue=!0,this.call=this.flush,Reflect.defineProperty(e,"$fastController",{value:this,enumerable:!1})}subscribe(e){this.flush(),super.subscribe(e)}addSplice(e){void 0===this.splices?this.splices=[e]:this.splices.push(e),this.needsQueue&&(this.needsQueue=!1,el.queueUpdate(this))}reset(e){this.oldCollection=e,this.needsQueue&&(this.needsQueue=!1,el.queueUpdate(this))}flush(){let e=this.splices,t=this.oldCollection;if(void 0===e&&void 0===t)return;this.needsQueue=!0,this.splices=void 0,this.oldCollection=void 0;let i=void 0===t?function(e,t){let i=[],o=function(e){let t=[];for(let i=0,o=e.length;i<o;i++){let o=e[i];!function(e,t,i,o){let s=ne(t,i,o),r=!1,a=0;for(let t=0;t<e.length;t++){var n,l,h,d;let i=e[t];if(i.index+=a,r)continue;let o=(n=s.index,l=s.index+s.removed.length,h=i.index,d=i.index+i.addedCount,l<h||d<n?-1:l===h||d===n?0:n<h?l<d?l-h:d-h:d<l?d-n:l-n);if(o>=0){e.splice(t,1),t--,a-=i.addedCount-i.removed.length,s.addedCount+=i.addedCount-o;let n=s.removed.length+i.removed.length-o;if(s.addedCount||n){let e=i.removed;if(s.index<i.index){let t=s.removed.slice(0,i.index-s.index);ni.apply(t,e),e=t}if(s.index+s.removed.length>i.index+i.addedCount){let t=s.removed.slice(i.index+i.addedCount-s.index);ni.apply(e,t)}s.removed=e,i.index<s.index&&(s.index=i.index)}else r=!0}else if(s.index<i.index){r=!0,e.splice(t,0,s),t++;let o=s.addedCount-s.removed.length;i.index+=o,a+=o}}r||e.push(s)}(t,o.index,o.removed,o.addedCount)}return t}(t);for(let t=0,s=o.length;t<s;++t){let s=o[t];if(1===s.addedCount&&1===s.removed.length){s.removed[0]!==e[s.index]&&i.push(s);continue}i=i.concat(nt(e,s.index,s.index+s.addedCount,s.removed,0,s.removed.length))}return i}(this.source,e):nt(this.source,0,this.source.length,t,0,t.length);this.notify(i)}}function na(e,t,i,o){e.bind(t[i],o)}function nn(e,t,i,o){let s=Object.create(o);s.index=i,s.length=t.length,e.bind(t[i],s)}Object.freeze({positioning:!1,recycle:!0});class nl{constructor(e,t,i,o,s,r){this.location=e,this.itemsBinding=t,this.templateBinding=o,this.options=r,this.source=null,this.views=[],this.items=null,this.itemsObserver=null,this.originalContext=void 0,this.childContext=void 0,this.bindView=na,this.itemsBindingObserver=ec.binding(t,this,i),this.templateBindingObserver=ec.binding(o,this,s),r.positioning&&(this.bindView=nn)}bind(e,t){this.source=e,this.originalContext=t,this.childContext=Object.create(t),this.childContext.parent=e,this.childContext.parentContext=this.originalContext,this.items=this.itemsBindingObserver.observe(e,this.originalContext),this.template=this.templateBindingObserver.observe(e,this.originalContext),this.observeItems(!0),this.refreshAllViews()}unbind(){this.source=null,this.items=null,null!==this.itemsObserver&&this.itemsObserver.unsubscribe(this),this.unbindAllViews(),this.itemsBindingObserver.disconnect(),this.templateBindingObserver.disconnect()}handleChange(e,t){e===this.itemsBinding?(this.items=this.itemsBindingObserver.observe(this.source,this.originalContext),this.observeItems(),this.refreshAllViews()):e===this.templateBinding?(this.template=this.templateBindingObserver.observe(this.source,this.originalContext),this.refreshAllViews(!0)):this.updateViews(t)}observeItems(e=!1){if(!this.items){this.items=ee;return}let t=this.itemsObserver,i=this.itemsObserver=ec.getNotifier(this.items),o=t!==i;o&&null!==t&&t.unsubscribe(this),(o||e)&&i.subscribe(this)}updateViews(e){let t=this.childContext,i=this.views,o=this.bindView,s=this.items,r=this.template,a=this.options.recycle,n=[],l=0,h=0;for(let d=0,c=e.length;d<c;++d){let c=e[d],u=c.removed,p=0,g=c.index,m=g+c.addedCount,b=i.splice(c.index,u.length),f=h=n.length+b.length;for(;g<m;++g){let e,d=i[g],c=d?d.firstChild:this.location;a&&h>0?(p<=f&&b.length>0?(e=b[p],p++):(e=n[l],l++),h--):e=r.create(),i.splice(g,0,e),o(e,s,g,t),e.insertBefore(c)}b[p]&&n.push(...b.slice(p))}for(let e=l,t=n.length;e<t;++e)n[e].dispose();if(this.options.positioning)for(let e=0,t=i.length;e<t;++e){let o=i[e].context;o.length=t,o.index=e}}refreshAllViews(e=!1){let t=this.items,i=this.childContext,o=this.template,s=this.location,r=this.bindView,a=t.length,n=this.views,l=n.length;if((0===a||e||!this.options.recycle)&&(s2.disposeContiguousBatch(n),l=0),0===l){this.views=n=Array(a);for(let e=0;e<a;++e){let a=o.create();r(a,t,e,i),n[e]=a,a.insertBefore(s)}}else{let e=0;for(;e<a;++e)if(e<l)r(n[e],t,e,i);else{let a=o.create();r(a,t,e,i),n.push(a),a.insertBefore(s)}let h=n.splice(e,l-e);for(e=0,a=h.length;e<a;++e)h[e].dispose()}}unbindAllViews(){let e=this.views;for(let t=0,i=e.length;t<i;++t)e[t].unbind()}}class nh extends sV{constructor(e,t,i){super(),this.itemsBinding=e,this.templateBinding=t,this.options=i,this.createPlaceholder=el.createBlockPlaceholder,function(){if(no)return;no=!0,ec.setArrayObserverFactory(e=>new nr(e));let e=Array.prototype;if(e.$fastPatch)return;Reflect.defineProperty(e,"$fastPatch",{value:1,enumerable:!1});let t=e.pop,i=e.push,o=e.reverse,s=e.shift,r=e.sort,a=e.splice,n=e.unshift;e.pop=function(){let e=this.length>0,i=t.apply(this,arguments),o=this.$fastController;return void 0!==o&&e&&o.addSplice(ne(this.length,[i],0)),i},e.push=function(){let e=i.apply(this,arguments),t=this.$fastController;return void 0!==t&&t.addSplice(ns(ne(this.length-arguments.length,[],arguments.length),this)),e},e.reverse=function(){let e,t=this.$fastController;void 0!==t&&(t.flush(),e=this.slice());let i=o.apply(this,arguments);return void 0!==t&&t.reset(e),i},e.shift=function(){let e=this.length>0,t=s.apply(this,arguments),i=this.$fastController;return void 0!==i&&e&&i.addSplice(ne(0,[t],0)),t},e.sort=function(){let e,t=this.$fastController;void 0!==t&&(t.flush(),e=this.slice());let i=r.apply(this,arguments);return void 0!==t&&t.reset(e),i},e.splice=function(){let e=a.apply(this,arguments),t=this.$fastController;return void 0!==t&&t.addSplice(ns(ne(+arguments[0],e,arguments.length>2?arguments.length-2:0),this)),e},e.unshift=function(){let e=n.apply(this,arguments),t=this.$fastController;return void 0!==t&&t.addSplice(ns(ne(0,[],arguments.length),this)),e}}(),this.isItemsBindingVolatile=ec.isVolatileBinding(e),this.isTemplateBindingVolatile=ec.isVolatileBinding(t)}createBehavior(e){return new nl(e,this.itemsBinding,this.isItemsBindingVolatile,this.templateBinding,this.isTemplateBindingVolatile,this.options)}}class nd extends sR{constructor(){super(...arguments),this.rowType=a6,this.rowData=null,this.columnDefinitions=null,this.isActiveRow=!1,this.cellsRepeatBehavior=null,this.cellsPlaceholder=null,this.focusColumnIndex=0,this.refocusOnLoad=!1,this.updateRowStyle=()=>{this.style.gridTemplateColumns=this.gridTemplateColumns}}gridTemplateColumnsChanged(){this.$fastController.isConnected&&this.updateRowStyle()}rowTypeChanged(){this.$fastController.isConnected&&this.updateItemTemplate()}rowDataChanged(){if(null!==this.rowData&&this.isActiveRow){this.refocusOnLoad=!0;return}}cellItemTemplateChanged(){this.updateItemTemplate()}headerCellItemTemplateChanged(){this.updateItemTemplate()}connectedCallback(){super.connectedCallback(),null===this.cellsRepeatBehavior&&(this.cellsPlaceholder=document.createComment(""),this.appendChild(this.cellsPlaceholder),this.updateItemTemplate(),this.cellsRepeatBehavior=new nh(e=>e.columnDefinitions,e=>e.activeCellItemTemplate,{positioning:!0}).createBehavior(this.cellsPlaceholder),this.$fastController.addBehaviors([this.cellsRepeatBehavior])),this.addEventListener("cell-focused",this.handleCellFocus),this.addEventListener(a1,this.handleFocusout),this.addEventListener(a2,this.handleKeydown),this.updateRowStyle(),this.refocusOnLoad&&(this.refocusOnLoad=!1,this.cellElements.length>this.focusColumnIndex&&this.cellElements[this.focusColumnIndex].focus())}disconnectedCallback(){super.disconnectedCallback(),this.removeEventListener("cell-focused",this.handleCellFocus),this.removeEventListener(a1,this.handleFocusout),this.removeEventListener(a2,this.handleKeydown)}handleFocusout(e){this.contains(e.target)||(this.isActiveRow=!1,this.focusColumnIndex=0)}handleCellFocus(e){this.isActiveRow=!0,this.focusColumnIndex=this.cellElements.indexOf(e.target),this.$emit("row-focused",this)}handleKeydown(e){if(e.defaultPrevented)return;let t=0;switch(e.key){case oZ:t=Math.max(0,this.focusColumnIndex-1),this.cellElements[t].focus(),e.preventDefault();break;case oJ:t=Math.min(this.cellElements.length-1,this.focusColumnIndex+1),this.cellElements[t].focus(),e.preventDefault();break;case"Home":e.ctrlKey||(this.cellElements[0].focus(),e.preventDefault());break;case"End":e.ctrlKey||(this.cellElements[this.cellElements.length-1].focus(),e.preventDefault())}}updateItemTemplate(){this.activeCellItemTemplate=this.rowType===a6&&void 0!==this.cellItemTemplate?this.cellItemTemplate:this.rowType===a6&&void 0===this.cellItemTemplate?this.defaultCellItemTemplate:void 0!==this.headerCellItemTemplate?this.headerCellItemTemplate:this.defaultHeaderCellItemTemplate}}X([eS({attribute:"grid-template-columns"})],nd.prototype,"gridTemplateColumns",void 0),X([eS({attribute:"row-type"})],nd.prototype,"rowType",void 0),X([eu],nd.prototype,"rowData",void 0),X([eu],nd.prototype,"columnDefinitions",void 0),X([eu],nd.prototype,"cellItemTemplate",void 0),X([eu],nd.prototype,"headerCellItemTemplate",void 0),X([eu],nd.prototype,"rowIndex",void 0),X([eu],nd.prototype,"isActiveRow",void 0),X([eu],nd.prototype,"activeCellItemTemplate",void 0),X([eu],nd.prototype,"defaultCellItemTemplate",void 0),X([eu],nd.prototype,"defaultHeaderCellItemTemplate",void 0),X([eu],nd.prototype,"cellElements",void 0);class nc extends rr{constructor(e,t){super(e,t),this.observer=null,t.childList=!0}observe(){null===this.observer&&(this.observer=new MutationObserver(this.handleEvent.bind(this))),this.observer.observe(this.target,this.options)}disconnect(){this.observer.disconnect()}getNodes(){return"subtree"in this.options?Array.from(this.target.querySelectorAll(this.options.selector)):Array.from(this.target.childNodes)}}function nu(e){return"string"==typeof e&&(e={property:e}),new sH("fast-children",nc,e)}class np extends sR{constructor(){super(),this.noTabbing=!1,this.generateHeader="default",this.rowsData=[],this.columnDefinitions=null,this.focusRowIndex=0,this.focusColumnIndex=0,this.rowsPlaceholder=null,this.generatedHeader=null,this.isUpdatingFocus=!1,this.pendingFocusUpdate=!1,this.rowindexUpdateQueued=!1,this.columnDefinitionsStale=!0,this.generatedGridTemplateColumns="",this.focusOnCell=(e,t,i)=>{if(0===this.rowElements.length){this.focusRowIndex=0,this.focusColumnIndex=0;return}let o=Math.max(0,Math.min(this.rowElements.length-1,e)),s=this.rowElements[o].querySelectorAll('[role="cell"], [role="gridcell"], [role="columnheader"], [role="rowheader"]'),r=Math.max(0,Math.min(s.length-1,t)),a=s[r];i&&this.scrollHeight!==this.clientHeight&&(o<this.focusRowIndex&&this.scrollTop>0||o>this.focusRowIndex&&this.scrollTop<this.scrollHeight-this.clientHeight)&&a.scrollIntoView({block:"center",inline:"center"}),a.focus()},this.onChildListChange=(e,t)=>{e&&e.length&&(e.forEach(e=>{e.addedNodes.forEach(e=>{1===e.nodeType&&"row"===e.getAttribute("role")&&(e.columnDefinitions=this.columnDefinitions)})}),this.queueRowIndexUpdate())},this.queueRowIndexUpdate=()=>{this.rowindexUpdateQueued||(this.rowindexUpdateQueued=!0,el.queueUpdate(this.updateRowIndexes))},this.updateRowIndexes=()=>{let e=this.gridTemplateColumns;if(void 0===e){if(""===this.generatedGridTemplateColumns&&this.rowElements.length>0){let e=this.rowElements[0];this.generatedGridTemplateColumns=Array(e.cellElements.length).fill("1fr").join(" ")}e=this.generatedGridTemplateColumns}this.rowElements.forEach((t,i)=>{t.rowIndex=i,t.gridTemplateColumns=e,this.columnDefinitionsStale&&(t.columnDefinitions=this.columnDefinitions)}),this.rowindexUpdateQueued=!1,this.columnDefinitionsStale=!1}}static generateTemplateColumns(e){let t="";return e.forEach(e=>{t=`${t}${""===t?"":" "}1fr`}),t}noTabbingChanged(){this.$fastController.isConnected&&(this.noTabbing?this.setAttribute("tabIndex","-1"):this.setAttribute("tabIndex",this.contains(document.activeElement)||this===document.activeElement?"-1":"0"))}generateHeaderChanged(){this.$fastController.isConnected&&this.toggleGeneratedHeader()}gridTemplateColumnsChanged(){this.$fastController.isConnected&&this.updateRowIndexes()}rowsDataChanged(){null===this.columnDefinitions&&this.rowsData.length>0&&(this.columnDefinitions=np.generateColumns(this.rowsData[0])),this.$fastController.isConnected&&this.toggleGeneratedHeader()}columnDefinitionsChanged(){if(null===this.columnDefinitions){this.generatedGridTemplateColumns="";return}this.generatedGridTemplateColumns=np.generateTemplateColumns(this.columnDefinitions),this.$fastController.isConnected&&(this.columnDefinitionsStale=!0,this.queueRowIndexUpdate())}headerCellItemTemplateChanged(){this.$fastController.isConnected&&null!==this.generatedHeader&&(this.generatedHeader.headerCellItemTemplate=this.headerCellItemTemplate)}focusRowIndexChanged(){this.$fastController.isConnected&&this.queueFocusUpdate()}focusColumnIndexChanged(){this.$fastController.isConnected&&this.queueFocusUpdate()}connectedCallback(){super.connectedCallback(),void 0===this.rowItemTemplate&&(this.rowItemTemplate=this.defaultRowItemTemplate),this.rowsPlaceholder=document.createComment(""),this.appendChild(this.rowsPlaceholder),this.toggleGeneratedHeader(),this.rowsRepeatBehavior=new nh(e=>e.rowsData,e=>e.rowItemTemplate,{positioning:!0}).createBehavior(this.rowsPlaceholder),this.$fastController.addBehaviors([this.rowsRepeatBehavior]),this.addEventListener("row-focused",this.handleRowFocus),this.addEventListener("focus",this.handleFocus),this.addEventListener(a2,this.handleKeydown),this.addEventListener(a1,this.handleFocusOut),this.observer=new MutationObserver(this.onChildListChange),this.observer.observe(this,{childList:!0}),this.noTabbing&&this.setAttribute("tabindex","-1"),el.queueUpdate(this.queueRowIndexUpdate)}disconnectedCallback(){super.disconnectedCallback(),this.removeEventListener("row-focused",this.handleRowFocus),this.removeEventListener("focus",this.handleFocus),this.removeEventListener(a2,this.handleKeydown),this.removeEventListener(a1,this.handleFocusOut),this.observer.disconnect(),this.rowsPlaceholder=null,this.generatedHeader=null}handleRowFocus(e){this.isUpdatingFocus=!0;let t=e.target;this.focusRowIndex=this.rowElements.indexOf(t),this.focusColumnIndex=t.focusColumnIndex,this.setAttribute("tabIndex","-1"),this.isUpdatingFocus=!1}handleFocus(e){this.focusOnCell(this.focusRowIndex,this.focusColumnIndex,!0)}handleFocusOut(e){null!==e.relatedTarget&&this.contains(e.relatedTarget)||this.setAttribute("tabIndex",this.noTabbing?"-1":"0")}handleKeydown(e){let t;if(e.defaultPrevented)return;let i=this.rowElements.length-1,o=this.offsetHeight+this.scrollTop,s=this.rowElements[i];switch(e.key){case o0:e.preventDefault(),this.focusOnCell(this.focusRowIndex-1,this.focusColumnIndex,!0);break;case oQ:e.preventDefault(),this.focusOnCell(this.focusRowIndex+1,this.focusColumnIndex,!0);break;case"PageUp":if(e.preventDefault(),0===this.rowElements.length){this.focusOnCell(0,0,!1);break}if(0===this.focusRowIndex)return void this.focusOnCell(0,this.focusColumnIndex,!1);for(t=this.focusRowIndex-1;t>=0;t--){let e=this.rowElements[t];if(e.offsetTop<this.scrollTop){this.scrollTop=e.offsetTop+e.clientHeight-this.clientHeight;break}}this.focusOnCell(t,this.focusColumnIndex,!1);break;case"PageDown":if(e.preventDefault(),0===this.rowElements.length){this.focusOnCell(0,0,!1);break}if(this.focusRowIndex>=i||s.offsetTop+s.offsetHeight<=o)return void this.focusOnCell(i,this.focusColumnIndex,!1);for(t=this.focusRowIndex+1;t<=i;t++){let e=this.rowElements[t];if(e.offsetTop+e.offsetHeight>o){let t=0;this.generateHeader===a5&&null!==this.generatedHeader&&(t=this.generatedHeader.clientHeight),this.scrollTop=e.offsetTop-t;break}}this.focusOnCell(t,this.focusColumnIndex,!1);break;case"Home":e.ctrlKey&&(e.preventDefault(),this.focusOnCell(0,0,!0));break;case"End":e.ctrlKey&&null!==this.columnDefinitions&&(e.preventDefault(),this.focusOnCell(this.rowElements.length-1,this.columnDefinitions.length-1,!0))}}queueFocusUpdate(){this.isUpdatingFocus&&(this.contains(document.activeElement)||this===document.activeElement)||!1===this.pendingFocusUpdate&&(this.pendingFocusUpdate=!0,el.queueUpdate(()=>this.updateFocus()))}updateFocus(){this.pendingFocusUpdate=!1,this.focusOnCell(this.focusRowIndex,this.focusColumnIndex,!0)}toggleGeneratedHeader(){if(null!==this.generatedHeader&&(this.removeChild(this.generatedHeader),this.generatedHeader=null),"none"!==this.generateHeader&&this.rowsData.length>0){let e=document.createElement(this.rowElementTag);this.generatedHeader=e,this.generatedHeader.columnDefinitions=this.columnDefinitions,this.generatedHeader.gridTemplateColumns=this.gridTemplateColumns,this.generatedHeader.rowType=this.generateHeader===a5?"sticky-header":"header",(null!==this.firstChild||null!==this.rowsPlaceholder)&&this.insertBefore(e,null!==this.firstChild?this.firstChild:this.rowsPlaceholder);return}}}np.generateColumns=e=>Object.getOwnPropertyNames(e).map((e,t)=>({columnDataKey:e,gridColumn:`${t}`})),X([eS({attribute:"no-tabbing",mode:"boolean"})],np.prototype,"noTabbing",void 0),X([eS({attribute:"generate-header"})],np.prototype,"generateHeader",void 0),X([eS({attribute:"grid-template-columns"})],np.prototype,"gridTemplateColumns",void 0),X([eu],np.prototype,"rowsData",void 0),X([eu],np.prototype,"columnDefinitions",void 0),X([eu],np.prototype,"rowItemTemplate",void 0),X([eu],np.prototype,"cellItemTemplate",void 0),X([eu],np.prototype,"headerCellItemTemplate",void 0),X([eu],np.prototype,"focusRowIndex",void 0),X([eu],np.prototype,"focusColumnIndex",void 0),X([eu],np.prototype,"defaultRowItemTemplate",void 0),X([eu],np.prototype,"rowElementTag",void 0),X([eu],np.prototype,"rowElements",void 0);let ng=(e,t)=>rh`
  :host {
    display: flex;
    position: relative;
    flex-direction: column;
  }
`,nm=(e,t)=>rh`
  :host {
    display: grid;
    padding: 1px 0;
    box-sizing: border-box;
    width: 100%;
    border-bottom: calc(${tw} * 1px) solid ${ox};
  }

  :host(.header) {
  }

  :host(.sticky-header) {
    background: ${iX};
    position: sticky;
    top: 0;
  }
`,nb=(e,t)=>rh`
    :host {
      padding: calc(${tv} * 1px) calc(${tv} * 3px);
      color: ${op};
      box-sizing: border-box;
      font-family: ${tu};
      font-size: ${tC};
      line-height: ${tI};
      font-weight: 400;
      border: transparent calc(${tk} * 1px) solid;
      overflow: hidden;
      white-space: nowrap;
      border-radius: calc(${tb} * 1px);
    }

    :host(.column-header) {
      font-weight: 600;
    }

    :host(:${rm}) {
      outline: calc(${tk} * 1px) solid ${iE};
      color: ${op};
    }
  `.withBehaviors(rv(rh`
      :host {
        forced-color-adjust: none;
        border-color: transparent;
        background: ${nX.Field};
        color: ${nX.FieldText};
      }

      :host(:${rm}) {
        border-color: ${nX.FieldText};
        box-shadow: 0 0 0 2px inset ${nX.Field};
        color: ${nX.FieldText};
      }
    `));class nf extends a7{}let nv=nf.compose({baseName:"data-grid-cell",baseClass:a7,template:(e,t)=>s4`
        <template
            tabindex="-1"
            role="${e=>e.cellType&&"default"!==e.cellType?e.cellType:"gridcell"}"
            class="
            ${e=>"columnheader"===e.cellType?"column-header":"rowheader"===e.cellType?"row-header":""}
            "
        >
            <slot></slot>
        </template>
    `,styles:nb});class ny extends nd{}let nx=ny.compose({baseName:"data-grid-row",baseClass:nd,template:(e,t)=>{let i,o,s=(i=e.tagFor(a7),s4`
    <${i}
        cell-type="${e=>e.isRowHeader?"rowheader":void 0}"
        grid-column="${(e,t)=>t.index+1}"
        :rowData="${(e,t)=>t.parent.rowData}"
        :columnDefinition="${e=>e}"
    ></${i}>
`),r=(o=e.tagFor(a7),s4`
    <${o}
        cell-type="columnheader"
        grid-column="${(e,t)=>t.index+1}"
        :columnDefinition="${e=>e}"
    ></${o}>
`);return s4`
        <template
            role="row"
            class="${e=>"default"!==e.rowType?e.rowType:""}"
            :defaultCellItemTemplate="${s}"
            :defaultHeaderCellItemTemplate="${r}"
            ${nu({property:"cellElements",filter:rs('[role="cell"],[role="gridcell"],[role="columnheader"],[role="rowheader"]')})}
        >
            <slot ${rn("slottedCellElements")}></slot>
        </template>
    `},styles:nm});class n$ extends np{}let nw=n$.compose({baseName:"data-grid",baseClass:np,template:(e,t)=>{let i,o=(i=e.tagFor(nd),s4`
    <${i}
        :rowData="${e=>e}"
        :cellItemTemplate="${(e,t)=>t.parent.cellItemTemplate}"
        :headerCellItemTemplate="${(e,t)=>t.parent.headerCellItemTemplate}"
    ></${i}>
`),s=e.tagFor(nd);return s4`
        <template
            role="grid"
            tabindex="0"
            :rowElementTag="${()=>s}"
            :defaultRowItemTemplate="${o}"
            ${nu({property:"rowElements",filter:rs("[role=row]")})}
        >
            <slot></slot>
        </template>
    `},styles:ng});class nk extends sR{}class nC extends ah(nk){constructor(){super(...arguments),this.proxy=document.createElement("input")}}class nI extends nC{constructor(){super(...arguments),this.type="text"}readOnlyChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.readOnly=this.readOnly,this.validate())}autofocusChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.autofocus=this.autofocus,this.validate())}placeholderChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.placeholder=this.placeholder)}typeChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.type=this.type,this.validate())}listChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.setAttribute("list",this.list),this.validate())}maxlengthChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.maxLength=this.maxlength,this.validate())}minlengthChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.minLength=this.minlength,this.validate())}patternChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.pattern=this.pattern,this.validate())}sizeChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.size=this.size)}spellcheckChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.spellcheck=this.spellcheck)}connectedCallback(){super.connectedCallback(),this.proxy.setAttribute("type",this.type),this.validate(),this.autofocus&&el.queueUpdate(()=>{this.focus()})}select(){this.control.select(),this.$emit("select")}handleTextInput(){this.value=this.control.value}handleChange(){this.$emit("change")}validate(){super.validate(this.control)}}X([eS({attribute:"readonly",mode:"boolean"})],nI.prototype,"readOnly",void 0),X([eS({mode:"boolean"})],nI.prototype,"autofocus",void 0),X([eS],nI.prototype,"placeholder",void 0),X([eS],nI.prototype,"type",void 0),X([eS],nI.prototype,"list",void 0),X([eS({converter:eT})],nI.prototype,"maxlength",void 0),X([eS({converter:eT})],nI.prototype,"minlength",void 0),X([eS],nI.prototype,"pattern",void 0),X([eS({converter:eT})],nI.prototype,"size",void 0),X([eS({mode:"boolean"})],nI.prototype,"spellcheck",void 0),X([eu],nI.prototype,"defaultSlottedNodes",void 0);class nT{}rt(nT,rT),rt(nI,s8,nT);class nF extends sR{}class nS extends ah(nF){constructor(){super(...arguments),this.proxy=document.createElement("input")}}let nO={toView(e){if(null==e)return null;let t=new Date(e);return"Invalid Date"===t.toString()?null:`${t.getFullYear().toString().padStart(4,"0")}-${(t.getMonth()+1).toString().padStart(2,"0")}-${t.getDate().toString().padStart(2,"0")}`},fromView(e){if(null==e)return null;let t=new Date(e);return"Invalid Date"===t.toString()?null:t}},nD="Invalid Date";class nE extends nS{constructor(){super(...arguments),this.step=1,this.isUserInput=!1}readOnlyChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.readOnly=this.readOnly,this.validate())}autofocusChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.autofocus=this.autofocus,this.validate())}listChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.setAttribute("list",this.list),this.validate())}maxChanged(e,t){var i;this.max=t<(null!=(i=this.min)?i:t)?this.min:t,this.value=this.getValidValue(this.value)}minChanged(e,t){var i;this.min=t>(null!=(i=this.max)?i:t)?this.max:t,this.value=this.getValidValue(this.value)}get valueAsNumber(){return new Date(super.value).valueOf()}set valueAsNumber(e){this.value=new Date(e).toString()}get valueAsDate(){return new Date(super.value)}set valueAsDate(e){this.value=e.toString()}valueChanged(e,t){this.value=this.getValidValue(t),t===this.value&&(this.control&&!this.isUserInput&&(this.control.value=this.value),super.valueChanged(e,this.value),void 0===e||this.isUserInput||this.$emit("change"),this.isUserInput=!1)}getValidValue(e){var t,i;let o=new Date(e);return o.toString()===nD?o="":(o=(o=o>(null!=(t=this.max)?t:o)?this.max:o)<(null!=(i=this.min)?i:o)?this.min:o,o=`${o.getFullYear().toString().padStart(4,"0")}-${(o.getMonth()+1).toString().padStart(2,"0")}-${o.getDate().toString().padStart(2,"0")}`),o}stepUp(){let e=864e5*this.step,t=new Date(this.value);this.value=new Date(t.toString()!==nD?t.valueOf()+e:0).toString()}stepDown(){let e=864e5*this.step,t=new Date(this.value);this.value=new Date(t.toString()!==nD?Math.max(t.valueOf()-e,0):0).toString()}connectedCallback(){super.connectedCallback(),this.validate(),this.control.value=this.value,this.autofocus&&el.queueUpdate(()=>{this.focus()}),this.appearance||(this.appearance="outline")}handleTextInput(){this.isUserInput=!0,this.value=this.control.value}handleChange(){this.$emit("change")}handleKeyDown(e){switch(e.key){case o0:return this.stepUp(),!1;case oQ:return this.stepDown(),!1}return!0}handleBlur(){this.control.value=this.value}}rI([eS],nE.prototype,"appearance",void 0),rI([eS({attribute:"readonly",mode:"boolean"})],nE.prototype,"readOnly",void 0),rI([eS({mode:"boolean"})],nE.prototype,"autofocus",void 0),rI([eS],nE.prototype,"list",void 0),rI([eS({converter:eT})],nE.prototype,"step",void 0),rI([eS({converter:nO})],nE.prototype,"max",void 0),rI([eS({converter:nO})],nE.prototype,"min",void 0),rI([eu],nE.prototype,"defaultSlottedNodes",void 0),rt(nE,s8,nT);let nR=rh`
  ${ru("inline-block")} :host {
    font-family: ${tu};
    outline: none;
    user-select: none;
    /* Ensure to display focus highlight */
    margin: calc((${tk} - ${tw}) * 1px);
  }

  .root {
    box-sizing: border-box;
    position: relative;
    display: flex;
    flex-direction: row;
    color: ${op};
    background: ${i0};
    border-radius: calc(${tb} * 1px);
    border: calc(${tw} * 1px) solid ${oe};
    height: calc(${ry} * 1px);
  }

  :host([aria-invalid='true']) .root {
    border-color: ${oI};
  }

  .control {
    -webkit-appearance: none;
    font: inherit;
    background: transparent;
    border: 0;
    color: inherit;
    height: calc(100% - 4px);
    width: 100%;
    margin-top: auto;
    margin-bottom: auto;
    border: none;
    padding: 0 calc(${tv} * 2px + 1px);
    font-size: ${tC};
    line-height: ${tI};
  }

  .control:placeholder-shown {
    text-overflow: ellipsis;
  }

  .control:hover,
  .control:${rm},
  .control:disabled,
  .control:active {
    outline: none;
  }

  .label {
    display: block;
    color: ${op};
    cursor: pointer;
    font-size: ${tC};
    line-height: ${tI};
    margin-bottom: 4px;
  }

  .label__hidden {
    display: none;
    visibility: hidden;
  }

  .start,
  .end {
    margin: auto;
    fill: currentcolor;
  }

  ::slotted(svg) {
    /* TODO: adaptive typography https://github.com/microsoft/fast/issues/2432 */
    width: 16px;
    height: 16px;
  }

  .start {
    margin-inline-start: 11px;
  }

  .end {
    margin-inline-end: 11px;
  }

  :host(:hover:not([disabled])) .root {
    background: ${i1};
    border-color: ${ot};
  }

  :host([aria-invalid='true']:hover:not([disabled])) .root {
    border-color: ${oT};
  }

  :host(:active:not([disabled])) .root {
    background: ${i1};
    border-color: ${oi};
  }

  :host([aria-invalid='true']:active:not([disabled])) .root {
    border-color: ${oF};
  }

  :host(:focus-within:not([disabled])) .root {
    border-color: ${iE};
    box-shadow: 0 0 0 calc((${tk} - ${tw}) * 1px)
      ${iE};
  }

  :host([aria-invalid='true']:focus-within:not([disabled])) .root {
    border-color: ${oS};
    box-shadow: 0 0 0 calc((${tk} - ${tw}) * 1px)
      ${oS};
  }

  :host([appearance='filled']) .root {
    background: ${iX};
  }

  :host([appearance='filled']:hover:not([disabled])) .root {
    background: ${iY};
  }

  :host([disabled]) .label,
  :host([readonly]) .label,
  :host([readonly]) .control,
  :host([disabled]) .control {
    cursor: ${am};
  }

  :host([disabled]) {
    opacity: ${t$};
  }

  :host([disabled]) .control {
    border-color: ${om};
  }
`.withBehaviors(rv(rh`
    .root,
    :host([appearance='filled']) .root {
      forced-color-adjust: none;
      background: ${nX.Field};
      border-color: ${nX.FieldText};
    }
    :host([aria-invalid='true']) .root {
      border-style: dashed;
    }
    :host(:hover:not([disabled])) .root,
    :host([appearance='filled']:hover:not([disabled])) .root,
    :host([appearance='filled']:hover) .root {
      background: ${nX.Field};
      border-color: ${nX.Highlight};
    }
    .start,
    .end {
      fill: currentcolor;
    }
    :host([disabled]) {
      opacity: 1;
    }
    :host([disabled]) .root,
    :host([appearance='filled']:hover[disabled]) .root {
      border-color: ${nX.GrayText};
      background: ${nX.Field};
    }
    :host(:focus-within:enabled) .root {
      border-color: ${nX.Highlight};
      box-shadow: 0 0 0 calc((${tk} - ${tw}) * 1px)
        ${nX.Highlight};
    }
    input::placeholder {
      color: ${nX.GrayText};
    }
  `)),nL=(e,t)=>rh`
  ${nR}
`;function nA(e,t,i){return e.nodeType!==Node.TEXT_NODE||"string"==typeof e.nodeValue&&!!e.nodeValue.trim().length}let nV=(e,t)=>s4`
  <template class="${e=>e.readOnly?"readonly":""}">
    <label
      part="label"
      for="control"
      class="${e=>e.defaultSlottedNodes&&e.defaultSlottedNodes.length?"label":"label label__hidden"}"
    >
      <slot
        ${rn({property:"defaultSlottedNodes",filter:nA})}
      ></slot>
    </label>
    <div class="root" part="root">
      ${re(e,t)}
      <input
        class="control"
        part="control"
        id="control"
        @input="${e=>e.handleTextInput()}"
        @change="${e=>e.handleChange()}"
        ?autofocus="${e=>e.autofocus}"
        ?disabled="${e=>e.disabled}"
        list="${e=>e.list}"
        ?readonly="${e=>e.readOnly}"
        ?required="${e=>e.required}"
        :value="${e=>e.value}"
        type="date"
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
        ${s9("control")}
      />
      ${s7(e,t)}
    </div>
  </template>
`,nP=nE.compose({baseName:"date-field",styles:nL,template:nV,shadowOptions:{delegatesFocus:!0}}),nH={toView:e=>null==e?null:null==e?void 0:e.toColorString(),fromView(e){if(null==e)return null;let t=w(e);return t?q.create(t.r,t.g,t.b):null}},nz=rh`
  :host {
    background-color: ${iT};
    color: ${op};
  }
`.withBehaviors(rv(rh`
    :host {
      background-color: ${nX.ButtonFace};
      box-shadow: 0 0 0 1px ${nX.CanvasText};
      color: ${nX.ButtonText};
    }
  `));function nM(e){return(t,i)=>{t[i+"Changed"]=function(t,i){null!=i?e.setValueFor(this,i):e.deleteValueFor(this)}}}class nN extends sR{constructor(){super(),this.noPaint=!1;const e={handleChange:this.noPaintChanged.bind(this)};ec.getNotifier(this).subscribe(e,"fillColor"),ec.getNotifier(this).subscribe(e,"baseLayerLuminance")}noPaintChanged(){!this.noPaint&&(void 0!==this.fillColor||this.baseLayerLuminance)?this.$fastController.addStyles(nz):this.$fastController.removeStyles(nz)}}rI([eS({attribute:"no-paint",mode:"boolean"})],nN.prototype,"noPaint",void 0),rI([eS({attribute:"fill-color",converter:nH}),nM(iT)],nN.prototype,"fillColor",void 0),rI([eS({attribute:"accent-color",converter:nH,mode:"fromView"}),nM(ic)],nN.prototype,"accentColor",void 0),rI([eS({attribute:"neutral-color",converter:nH,mode:"fromView"}),nM(ih)],nN.prototype,"neutralColor",void 0),rI([eS({attribute:"error-color",converter:nH,mode:"fromView"}),nM(ow)],nN.prototype,"errorColor",void 0),rI([eS({converter:eT}),nM(tf)],nN.prototype,"density",void 0),rI([eS({attribute:"design-unit",converter:eT}),nM(tv)],nN.prototype,"designUnit",void 0),rI([eS({attribute:"direction"}),nM(tx)],nN.prototype,"direction",void 0),rI([eS({attribute:"base-height-multiplier",converter:eT}),nM(tp)],nN.prototype,"baseHeightMultiplier",void 0),rI([eS({attribute:"base-horizontal-spacing-multiplier",converter:eT}),nM(tg)],nN.prototype,"baseHorizontalSpacingMultiplier",void 0),rI([eS({attribute:"control-corner-radius",converter:eT}),nM(tb)],nN.prototype,"controlCornerRadius",void 0),rI([eS({attribute:"stroke-width",converter:eT}),nM(tw)],nN.prototype,"strokeWidth",void 0),rI([eS({attribute:"focus-stroke-width",converter:eT}),nM(tk)],nN.prototype,"focusStrokeWidth",void 0),rI([eS({attribute:"disabled-opacity",converter:eT}),nM(t$)],nN.prototype,"disabledOpacity",void 0),rI([eS({attribute:"type-ramp-minus-2-font-size"}),nM(tS)],nN.prototype,"typeRampMinus2FontSize",void 0),rI([eS({attribute:"type-ramp-minus-2-line-height"}),nM(tO)],nN.prototype,"typeRampMinus2LineHeight",void 0),rI([eS({attribute:"type-ramp-minus-1-font-size"}),nM(tT)],nN.prototype,"typeRampMinus1FontSize",void 0),rI([eS({attribute:"type-ramp-minus-1-line-height"}),nM(tF)],nN.prototype,"typeRampMinus1LineHeight",void 0),rI([eS({attribute:"type-ramp-base-font-size"}),nM(tC)],nN.prototype,"typeRampBaseFontSize",void 0),rI([eS({attribute:"type-ramp-base-line-height"}),nM(tI)],nN.prototype,"typeRampBaseLineHeight",void 0),rI([eS({attribute:"type-ramp-plus-1-font-size"}),nM(tD)],nN.prototype,"typeRampPlus1FontSize",void 0),rI([eS({attribute:"type-ramp-plus-1-line-height"}),nM(tE)],nN.prototype,"typeRampPlus1LineHeight",void 0),rI([eS({attribute:"type-ramp-plus-2-font-size"}),nM(tR)],nN.prototype,"typeRampPlus2FontSize",void 0),rI([eS({attribute:"type-ramp-plus-2-line-height"}),nM(tL)],nN.prototype,"typeRampPlus2LineHeight",void 0),rI([eS({attribute:"type-ramp-plus-3-font-size"}),nM(tA)],nN.prototype,"typeRampPlus3FontSize",void 0),rI([eS({attribute:"type-ramp-plus-3-line-height"}),nM(tV)],nN.prototype,"typeRampPlus3LineHeight",void 0),rI([eS({attribute:"type-ramp-plus-4-font-size"}),nM(tP)],nN.prototype,"typeRampPlus4FontSize",void 0),rI([eS({attribute:"type-ramp-plus-4-line-height"}),nM(tH)],nN.prototype,"typeRampPlus4LineHeight",void 0),rI([eS({attribute:"type-ramp-plus-5-font-size"}),nM(tz)],nN.prototype,"typeRampPlus5FontSize",void 0),rI([eS({attribute:"type-ramp-plus-5-line-height"}),nM(tM)],nN.prototype,"typeRampPlus5LineHeight",void 0),rI([eS({attribute:"type-ramp-plus-6-font-size"}),nM(tN)],nN.prototype,"typeRampPlus6FontSize",void 0),rI([eS({attribute:"type-ramp-plus-6-line-height"}),nM(tB)],nN.prototype,"typeRampPlus6LineHeight",void 0),rI([eS({attribute:"accent-fill-rest-delta",converter:eT}),nM(tj)],nN.prototype,"accentFillRestDelta",void 0),rI([eS({attribute:"accent-fill-hover-delta",converter:eT}),nM(tq)],nN.prototype,"accentFillHoverDelta",void 0),rI([eS({attribute:"accent-fill-active-delta",converter:eT}),nM(tU)],nN.prototype,"accentFillActiveDelta",void 0),rI([eS({attribute:"accent-fill-focus-delta",converter:eT}),nM(t_)],nN.prototype,"accentFillFocusDelta",void 0),rI([eS({attribute:"accent-foreground-rest-delta",converter:eT}),nM(tG)],nN.prototype,"accentForegroundRestDelta",void 0),rI([eS({attribute:"accent-foreground-hover-delta",converter:eT}),nM(tW)],nN.prototype,"accentForegroundHoverDelta",void 0),rI([eS({attribute:"accent-foreground-active-delta",converter:eT}),nM(tK)],nN.prototype,"accentForegroundActiveDelta",void 0),rI([eS({attribute:"accent-foreground-focus-delta",converter:eT}),nM(tX)],nN.prototype,"accentForegroundFocusDelta",void 0),rI([eS({attribute:"neutral-fill-rest-delta",converter:eT}),nM(tY)],nN.prototype,"neutralFillRestDelta",void 0),rI([eS({attribute:"neutral-fill-hover-delta",converter:eT}),nM(tQ)],nN.prototype,"neutralFillHoverDelta",void 0),rI([eS({attribute:"neutral-fill-active-delta",converter:eT}),nM(tZ)],nN.prototype,"neutralFillActiveDelta",void 0),rI([eS({attribute:"neutral-fill-focus-delta",converter:eT}),nM(tJ)],nN.prototype,"neutralFillFocusDelta",void 0),rI([eS({attribute:"neutral-fill-input-rest-delta",converter:eT}),nM(t0)],nN.prototype,"neutralFillInputRestDelta",void 0),rI([eS({attribute:"neutral-fill-input-hover-delta",converter:eT}),nM(t1)],nN.prototype,"neutralFillInputHoverDelta",void 0),rI([eS({attribute:"neutral-fill-input-active-delta",converter:eT}),nM(t2)],nN.prototype,"neutralFillInputActiveDelta",void 0),rI([eS({attribute:"neutral-fill-input-focus-delta",converter:eT}),nM(t5)],nN.prototype,"neutralFillInputFocusDelta",void 0),rI([eS({attribute:"neutral-fill-stealth-rest-delta",converter:eT}),nM(t3)],nN.prototype,"neutralFillStealthRestDelta",void 0),rI([eS({attribute:"neutral-fill-stealth-hover-delta",converter:eT}),nM(t4)],nN.prototype,"neutralFillStealthHoverDelta",void 0),rI([eS({attribute:"neutral-fill-stealth-active-delta",converter:eT}),nM(t6)],nN.prototype,"neutralFillStealthActiveDelta",void 0),rI([eS({attribute:"neutral-fill-stealth-focus-delta",converter:eT}),nM(t9)],nN.prototype,"neutralFillStealthFocusDelta",void 0),rI([eS({attribute:"neutral-fill-strong-hover-delta",converter:eT}),nM(t7)],nN.prototype,"neutralFillStrongHoverDelta",void 0),rI([eS({attribute:"neutral-fill-strong-active-delta",converter:eT}),nM(ie)],nN.prototype,"neutralFillStrongActiveDelta",void 0),rI([eS({attribute:"neutral-fill-strong-focus-delta",converter:eT}),nM(it)],nN.prototype,"neutralFillStrongFocusDelta",void 0),rI([eS({attribute:"base-layer-luminance",converter:eT}),nM(tm)],nN.prototype,"baseLayerLuminance",void 0),rI([eS({attribute:"neutral-fill-layer-rest-delta",converter:eT}),nM(ii)],nN.prototype,"neutralFillLayerRestDelta",void 0),rI([eS({attribute:"neutral-stroke-divider-rest-delta",converter:eT}),nM(il)],nN.prototype,"neutralStrokeDividerRestDelta",void 0),rI([eS({attribute:"neutral-stroke-rest-delta",converter:eT}),nM(io)],nN.prototype,"neutralStrokeRestDelta",void 0),rI([eS({attribute:"neutral-stroke-hover-delta",converter:eT}),nM(is)],nN.prototype,"neutralStrokeHoverDelta",void 0),rI([eS({attribute:"neutral-stroke-active-delta",converter:eT}),nM(ir)],nN.prototype,"neutralStrokeActiveDelta",void 0),rI([eS({attribute:"neutral-stroke-focus-delta",converter:eT}),nM(ia)],nN.prototype,"neutralStrokeFocusDelta",void 0);let nB=(e,t)=>s4` <slot></slot> `,nj=(e,t)=>rh`
  ${ru("block")}
`,nq=nN.compose({baseName:"design-system-provider",template:nB,styles:nj});var nU,n_,nG,nW,nK,nX,nY,nQ,nZ,nJ,n0,n1,n2=["input","select","textarea","a[href]","button","[tabindex]:not(slot)","audio[controls]","video[controls]",'[contenteditable]:not([contenteditable="false"])',"details>summary:first-of-type","details"],n5=n2.join(","),n3="u"<typeof Element,n4=n3?function(){}:Element.prototype.matches||Element.prototype.msMatchesSelector||Element.prototype.webkitMatchesSelector,n6=!n3&&Element.prototype.getRootNode?function(e){return e.getRootNode()}:function(e){return e.ownerDocument},n9=function(e){return"INPUT"===e.tagName},n8=function(e,t){for(var i=0;i<e.length;i++)if(e[i].checked&&e[i].form===t)return e[i]},n7=function(e){if(!e.name)return!0;var t,i=e.form||n6(e),o=function(e){return i.querySelectorAll('input[type="radio"][name="'+e+'"]')};if("u">typeof window&&void 0!==window.CSS&&"function"==typeof window.CSS.escape)t=o(window.CSS.escape(e.name));else try{t=o(e.name)}catch(e){return console.error("Looks like you have a radio button with a name attribute containing invalid CSS selector characters and need the CSS.escape polyfill: %s",e.message),!1}var s=n8(t,e.form);return!s||s===e},le=function(e){return n9(e)&&"radio"===e.type&&!n7(e)},lt=function(e){var t=e.getBoundingClientRect(),i=t.width,o=t.height;return 0===i&&0===o},li=function(e,t){var i=t.displayCheck,o=t.getShadowRoot;if("hidden"===getComputedStyle(e).visibility)return!0;var s=n4.call(e,"details>summary:first-of-type")?e.parentElement:e;if(n4.call(s,"details:not([open]) *"))return!0;var r=n6(e).host,a=(null==r?void 0:r.ownerDocument.contains(r))||e.ownerDocument.contains(e);if(i&&"full"!==i){if("non-zero-area"===i)return lt(e)}else{if("function"==typeof o){for(var n=e;e;){var l=e.parentElement,h=n6(e);if(l&&!l.shadowRoot&&!0===o(l))return lt(e);e=e.assignedSlot?e.assignedSlot:l||h===e.ownerDocument?l:h.host}e=n}if(a)return!e.getClientRects().length}return!1},lo=function(e){if(/^(INPUT|BUTTON|SELECT|TEXTAREA)$/.test(e.tagName))for(var t=e.parentElement;t;){if("FIELDSET"===t.tagName&&t.disabled){for(var i=0;i<t.children.length;i++){var o=t.children.item(i);if("LEGEND"===o.tagName)return!!n4.call(t,"fieldset[disabled] *")||!o.contains(e)}return!0}t=t.parentElement}return!1},ls=function(e,t){return!(t.disabled||n9(t)&&"hidden"===t.type||li(t,e)||"DETAILS"===t.tagName&&Array.prototype.slice.apply(t.children).some(function(e){return"SUMMARY"===e.tagName})||lo(t))},lr=function(e,t){return!le(t)&&!(0>(t.tabIndex<0&&(/^(AUDIO|VIDEO|DETAILS)$/.test(t.tagName)||t.isContentEditable)&&isNaN(parseInt(t.getAttribute("tabindex"),10))?0:t.tabIndex))&&!!ls(e,t)},la=function(e,t){if(t=t||{},!e)throw Error("No node provided");return!1!==n4.call(e,n5)&&lr(t,e)},ln=n2.concat("iframe").join(","),ll=function(e,t){if(t=t||{},!e)throw Error("No node provided");return!1!==n4.call(e,ln)&&ls(t,e)};class lh extends sR{constructor(){super(...arguments),this.modal=!0,this.hidden=!1,this.trapFocus=!0,this.trapFocusChanged=()=>{this.$fastController.isConnected&&this.updateTrapFocus()},this.isTrappingFocus=!1,this.handleDocumentKeydown=e=>{if(!e.defaultPrevented&&!this.hidden)switch(e.key){case"Escape":this.dismiss(),e.preventDefault();break;case"Tab":this.handleTabKeyDown(e)}},this.handleDocumentFocus=e=>{!e.defaultPrevented&&this.shouldForceFocus(e.target)&&(this.focusFirstElement(),e.preventDefault())},this.handleTabKeyDown=e=>{if(!this.trapFocus||this.hidden)return;let t=this.getTabQueueBounds();if(0!==t.length){if(1===t.length){t[0].focus(),e.preventDefault();return}e.shiftKey&&e.target===t[0]?(t[t.length-1].focus(),e.preventDefault()):e.shiftKey||e.target!==t[t.length-1]||(t[0].focus(),e.preventDefault())}},this.getTabQueueBounds=()=>lh.reduceTabbableItems([],this),this.focusFirstElement=()=>{let e=this.getTabQueueBounds();e.length>0?e[0].focus():this.dialog instanceof HTMLElement&&this.dialog.focus()},this.shouldForceFocus=e=>this.isTrappingFocus&&!this.contains(e),this.shouldTrapFocus=()=>this.trapFocus&&!this.hidden,this.updateTrapFocus=e=>{let t=void 0===e?this.shouldTrapFocus():e;t&&!this.isTrappingFocus?(this.isTrappingFocus=!0,document.addEventListener("focusin",this.handleDocumentFocus),el.queueUpdate(()=>{this.shouldForceFocus(document.activeElement)&&this.focusFirstElement()})):!t&&this.isTrappingFocus&&(this.isTrappingFocus=!1,document.removeEventListener("focusin",this.handleDocumentFocus))}}dismiss(){this.$emit("dismiss"),this.$emit("cancel")}show(){this.hidden=!1}hide(){this.hidden=!0,this.$emit("close")}connectedCallback(){super.connectedCallback(),document.addEventListener("keydown",this.handleDocumentKeydown),this.notifier=ec.getNotifier(this),this.notifier.subscribe(this,"hidden"),this.updateTrapFocus()}disconnectedCallback(){super.disconnectedCallback(),document.removeEventListener("keydown",this.handleDocumentKeydown),this.updateTrapFocus(!1),this.notifier.unsubscribe(this,"hidden")}handleChange(e,t){"hidden"===t&&this.updateTrapFocus()}static reduceTabbableItems(e,t){return"-1"===t.getAttribute("tabindex")?e:la(t)||lh.isFocusableFastElement(t)&&lh.hasTabbableShadow(t)?(e.push(t),e):t.childElementCount?e.concat(Array.from(t.children).reduce(lh.reduceTabbableItems,[])):e}static isFocusableFastElement(e){var t,i;return!!(null==(i=null==(t=e.$fastController)?void 0:t.definition.shadowOptions)?void 0:i.delegatesFocus)}static hasTabbableShadow(e){var t,i;return Array.from(null!=(i=null==(t=e.shadowRoot)?void 0:t.querySelectorAll("*"))?i:[]).some(e=>la(e))}}X([eS({mode:"boolean"})],lh.prototype,"modal",void 0),X([eS({mode:"boolean"})],lh.prototype,"hidden",void 0),X([eS({attribute:"trap-focus",mode:"boolean"})],lh.prototype,"trapFocus",void 0),X([eS({attribute:"aria-describedby"})],lh.prototype,"ariaDescribedby",void 0),X([eS({attribute:"aria-labelledby"})],lh.prototype,"ariaLabelledby",void 0),X([eS({attribute:"aria-label"})],lh.prototype,"ariaLabel",void 0);let ld=(e,t)=>rh`
  :host([hidden]) {
    display: none;
  }

  :host {
    --elevation: 14;
    --dialog-height: 480px;
    --dialog-width: 640px;
    display: block;
  }

  .overlay {
    position: fixed;
    top: 0;
    left: 0;
    right: 0;
    bottom: 0;
    background: rgba(0, 0, 0, 0.3);
    touch-action: none;
  }

  .positioning-region {
    display: flex;
    justify-content: center;
    position: fixed;
    top: 0;
    bottom: 0;
    left: 0;
    right: 0;
    overflow: auto;
  }

  .control {
    ${ax}
    margin-top: auto;
    margin-bottom: auto;
    width: var(--dialog-width);
    height: var(--dialog-height);
    background-color: ${iT};
    z-index: 1;
    border-radius: calc(${tb} * 1px);
    border: calc(${tw} * 1px) solid transparent;
  }
`;class lc extends lh{}let lu=lc.compose({baseName:"dialog",baseClass:lh,template:(e,t)=>s4`
    <div class="positioning-region" part="positioning-region">
        ${rG(e=>e.modal,s4`
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
            ${s9("dialog")}
        >
            <slot></slot>
        </div>
    </div>
`,styles:ld});class lp extends sR{connectedCallback(){super.connectedCallback(),this.setup()}disconnectedCallback(){super.disconnectedCallback(),this.details.removeEventListener("toggle",this.onToggle)}show(){this.details.open=!0}hide(){this.details.open=!1}toggle(){this.details.open=!this.details.open}setup(){this.onToggle=this.onToggle.bind(this),this.details.addEventListener("toggle",this.onToggle),this.expanded&&this.show()}onToggle(){this.expanded=this.details.open,this.$emit("toggle")}}X([eS({mode:"boolean"})],lp.prototype,"expanded",void 0),X([eS],lp.prototype,"title",void 0);let lg=(e,t)=>rh`
  .disclosure {
    transition: height 0.35s;
  }

  .disclosure .invoker::-webkit-details-marker {
    display: none;
  }

  .disclosure .invoker {
    list-style-type: none;
  }

  :host([appearance='accent']) .invoker {
    background: ${iS};
    color: ${iA};
    font-family: ${tu};
    font-size: ${tC};
    border-radius: calc(${tb} * 1px);
    outline: none;
    cursor: pointer;
    margin: 16px 0;
    padding: 12px;
    max-width: max-content;
  }

  :host([appearance='accent']) .invoker:active {
    background: ${iD};
    color: ${iP};
  }

  :host([appearance='accent']) .invoker:hover {
    background: ${iO};
    color: ${iV};
  }

  :host([appearance='lightweight']) .invoker {
    background: transparent;
    color: ${iU};
    border-bottom: calc(${tw} * 1px) solid ${iU};
    cursor: pointer;
    width: max-content;
    margin: 16px 0;
  }

  :host([appearance='lightweight']) .invoker:active {
    border-bottom-color: ${iG};
  }

  :host([appearance='lightweight']) .invoker:hover {
    border-bottom-color: ${i_};
  }

  .disclosure[open] .invoker ~ * {
    animation: fadeIn 0.5s ease-in-out;
  }

  @keyframes fadeIn {
    0% {
      opacity: 0;
    }
    100% {
      opacity: 1;
    }
  }
`;class lm extends lp{constructor(){super(...arguments),this.height=0,this.totalHeight=0}connectedCallback(){super.connectedCallback(),this.appearance||(this.appearance="accent")}appearanceChanged(e,t){e!==t&&(this.classList.add(t),this.classList.remove(e))}onToggle(){super.onToggle(),this.details.style.setProperty("height",`${this.disclosureHeight}px`)}setup(){super.setup();let e=()=>this.details.getBoundingClientRect().height;this.show(),this.totalHeight=e(),this.hide(),this.height=e(),this.expanded&&this.show()}get disclosureHeight(){return this.expanded?this.totalHeight:this.height}}rI([eS],lm.prototype,"appearance",void 0);let lb=lm.compose({baseName:"disclosure",baseClass:lp,template:(e,t)=>s4`
    <details class="disclosure" ${s9("details")}>
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
`,styles:lg}),lf="horizontal",lv="vertical";class ly extends sR{constructor(){super(...arguments),this.role="separator",this.orientation=lf}}X([eS],ly.prototype,"role",void 0),X([eS],ly.prototype,"orientation",void 0);let lx=(e,t)=>rh`
  ${ru("block")} :host {
    box-sizing: content-box;
    height: 0;
    margin: calc(${tv} * 1px) 0;
    border-top: calc(${tw} * 1px) solid ${ox};
    border-left: none;
  }

  :host([orientation='vertical']) {
    height: 100%;
    margin: 0 calc(${tv} * 1px);
    border-top: none;
    border-left: calc(${tw} * 1px) solid ${ox};
  }
`;class l$ extends ly{}let lw=l$.compose({baseName:"divider",baseClass:ly,template:(e,t)=>s4`
    <template role="${e=>e.role}" aria-orientation="${e=>e.orientation}"></template>
`,styles:lx});class lk extends aU{sizeChanged(e,t){super.sizeChanged(e,t),this.updateComputedStylesheet()}updateComputedStylesheet(){this.computedStylesheet&&this.$fastController.removeStyles(this.computedStylesheet);let e=`${this.size}`;this.computedStylesheet=rh`
      :host {
        --size: ${e};
      }
    `,this.$fastController.addStyles(this.computedStylesheet)}}let lC=lk.compose({baseName:"listbox",baseClass:aU,template:(e,t)=>s4`
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
            ${rn({filter:aU.slottedOptionFilter,flatten:!0,property:"slottedOptions"})}
        ></slot>
    </template>
`,styles:aX}),lI="menuitem",lT="menuitemcheckbox",lF="menuitemradio";class lS extends sR{constructor(){super(...arguments),this.role=lI,this.hasSubmenu=!1,this.currentDirection=nU.ltr,this.focusSubmenuOnLoad=!1,this.handleMenuItemKeyDown=e=>{if(e.defaultPrevented)return!1;switch(e.key){case"Enter":case" ":return this.invoke(),!1;case oJ:return this.expandAndFocus(),!1;case oZ:if(this.expanded)return this.expanded=!1,this.focus(),!1}return!0},this.handleMenuItemClick=e=>!e.defaultPrevented&&!this.disabled&&(this.invoke(),!1),this.submenuLoaded=()=>{this.focusSubmenuOnLoad&&(this.focusSubmenuOnLoad=!1,this.hasSubmenu&&(this.submenu.focus(),this.setAttribute("tabindex","-1")))},this.handleMouseOver=e=>!this.disabled&&!!this.hasSubmenu&&!this.expanded&&(this.expanded=!0,!1),this.handleMouseOut=e=>!(!this.expanded||this.contains(document.activeElement))&&(this.expanded=!1,!1),this.expandAndFocus=()=>{this.hasSubmenu&&(this.focusSubmenuOnLoad=!0,this.expanded=!0)},this.invoke=()=>{if(!this.disabled)switch(this.role){case lT:this.checked=!this.checked;break;case lI:this.updateSubmenu(),this.hasSubmenu?this.expandAndFocus():this.$emit("change");break;case lF:this.checked||(this.checked=!0)}},this.updateSubmenu=()=>{this.submenu=this.domChildren().find(e=>"menu"===e.getAttribute("role")),this.hasSubmenu=void 0!==this.submenu}}expandedChanged(e){this.$fastController.isConnected&&void 0!==this.submenu&&(!1===this.expanded?this.submenu.collapseExpandedItem():this.currentDirection=rj(this),this.$emit("expanded-change",this,{bubbles:!1}))}checkedChanged(e,t){this.$fastController.isConnected&&this.$emit("change")}connectedCallback(){super.connectedCallback(),el.queueUpdate(()=>{this.updateSubmenu()}),this.startColumnCount||(this.startColumnCount=1),this.observer=new MutationObserver(this.updateSubmenu)}disconnectedCallback(){super.disconnectedCallback(),this.submenu=void 0,void 0!==this.observer&&(this.observer.disconnect(),this.observer=void 0)}domChildren(){return Array.from(this.children).filter(e=>!e.hasAttribute("hidden"))}}X([eS({mode:"boolean"})],lS.prototype,"disabled",void 0),X([eS({mode:"boolean"})],lS.prototype,"expanded",void 0),X([eu],lS.prototype,"startColumnCount",void 0),X([eS],lS.prototype,"role",void 0),X([eS({mode:"boolean"})],lS.prototype,"checked",void 0),X([eu],lS.prototype,"submenuRegion",void 0),X([eu],lS.prototype,"hasSubmenu",void 0),X([eu],lS.prototype,"currentDirection",void 0),X([eu],lS.prototype,"submenu",void 0),rt(lS,s8);class lO extends sR{constructor(){super(...arguments),this.expandedItem=null,this.focusIndex=-1,this.isNestedMenu=()=>null!==this.parentElement&&rg(this.parentElement)&&"menuitem"===this.parentElement.getAttribute("role"),this.handleFocusOut=e=>{if(!this.contains(e.relatedTarget)&&void 0!==this.menuItems){this.collapseExpandedItem();let e=this.menuItems.findIndex(this.isFocusableElement);this.menuItems[this.focusIndex].setAttribute("tabindex","-1"),this.menuItems[e].setAttribute("tabindex","0"),this.focusIndex=e}},this.handleItemFocus=e=>{let t=e.target;void 0!==this.menuItems&&t!==this.menuItems[this.focusIndex]&&(this.menuItems[this.focusIndex].setAttribute("tabindex","-1"),this.focusIndex=this.menuItems.indexOf(t),t.setAttribute("tabindex","0"))},this.handleExpandedChanged=e=>{if(e.defaultPrevented||null===e.target||void 0===this.menuItems||0>this.menuItems.indexOf(e.target))return;e.preventDefault();let t=e.target;if(null!==this.expandedItem&&t===this.expandedItem&&!1===t.expanded){this.expandedItem=null;return}t.expanded&&(null!==this.expandedItem&&this.expandedItem!==t&&(this.expandedItem.expanded=!1),this.menuItems[this.focusIndex].setAttribute("tabindex","-1"),this.expandedItem=t,this.focusIndex=this.menuItems.indexOf(t),t.setAttribute("tabindex","0"))},this.removeItemListeners=()=>{void 0!==this.menuItems&&this.menuItems.forEach(e=>{e.removeEventListener("expanded-change",this.handleExpandedChanged),e.removeEventListener("focus",this.handleItemFocus)})},this.setItems=()=>{let e=this.domChildren();this.removeItemListeners(),this.menuItems=e;let t=this.menuItems.filter(this.isMenuItemElement);t.length&&(this.focusIndex=0);let i=t.reduce((e,t)=>{var i;let o,s,r=(o=(i=t).getAttribute("role"),s=i.querySelector("[slot=start]"),o!==lI&&null===s||o===lI&&null!==s?1:2*(o!==lI&&null!==s));return e>r?e:r},0);t.forEach((e,t)=>{e.setAttribute("tabindex",0===t?"0":"-1"),e.addEventListener("expanded-change",this.handleExpandedChanged),e.addEventListener("focus",this.handleItemFocus),(e instanceof lS||"startColumnCount"in e)&&(e.startColumnCount=i)})},this.changeHandler=e=>{if(void 0===this.menuItems)return;let t=e.target,i=this.menuItems.indexOf(t);if(-1!==i&&"menuitemradio"===t.role&&!0===t.checked){for(let e=i-1;e>=0;--e){let t=this.menuItems[e],i=t.getAttribute("role");if(i===lF&&(t.checked=!1),"separator"===i)break}let e=this.menuItems.length-1;for(let t=i+1;t<=e;++t){let e=this.menuItems[t],i=e.getAttribute("role");if(i===lF&&(e.checked=!1),"separator"===i)break}}},this.isMenuItemElement=e=>rg(e)&&lO.focusableElementRoles.hasOwnProperty(e.getAttribute("role")),this.isFocusableElement=e=>this.isMenuItemElement(e)}itemsChanged(e,t){this.$fastController.isConnected&&void 0!==this.menuItems&&this.setItems()}connectedCallback(){super.connectedCallback(),el.queueUpdate(()=>{this.setItems()}),this.addEventListener("change",this.changeHandler)}disconnectedCallback(){super.disconnectedCallback(),this.removeItemListeners(),this.menuItems=void 0,this.removeEventListener("change",this.changeHandler)}focus(){this.setFocus(0,1)}collapseExpandedItem(){null!==this.expandedItem&&(this.expandedItem.expanded=!1,this.expandedItem=null)}handleMenuKeyDown(e){if(!e.defaultPrevented&&void 0!==this.menuItems)switch(e.key){case oQ:this.setFocus(this.focusIndex+1,1);return;case o0:this.setFocus(this.focusIndex-1,-1);return;case"End":this.setFocus(this.menuItems.length-1,-1);return;case"Home":this.setFocus(0,1);return;default:return!0}}domChildren(){return Array.from(this.children).filter(e=>!e.hasAttribute("hidden"))}setFocus(e,t){if(void 0!==this.menuItems)for(;e>=0&&e<this.menuItems.length;){let i=this.menuItems[e];if(this.isFocusableElement(i)){this.focusIndex>-1&&this.menuItems.length>=this.focusIndex-1&&this.menuItems[this.focusIndex].setAttribute("tabindex","-1"),this.focusIndex=e,i.setAttribute("tabindex","0"),i.focus();break}e+=t}}}lO.focusableElementRoles={[lI]:"menuitem",[lT]:"menuitemcheckbox",[lF]:"menuitemradio"},X([eu],lO.prototype,"items",void 0);let lD=(e,t)=>rh`
    ${ru("block")} :host {
      --elevation: 11;
      background: ${iT};
      border: calc(${tw} * 1px) solid transparent;
      ${ax}
      margin: 0;
      border-radius: calc(${tb} * 1px);
      padding: calc(${tv} * 1px) 0;
      max-width: 368px;
      min-width: 64px;
    }

    :host([slot='submenu']) {
      width: max-content;
      margin: 0 calc(${tv} * 1px);
    }

    ::slotted(hr) {
      box-sizing: content-box;
      height: 0;
      margin: 0;
      border: none;
      border-top: calc(${tw} * 1px) solid ${ox};
    }
  `.withBehaviors(rv(rh`
      :host {
        background: ${nX.Canvas};
        border-color: ${nX.CanvasText};
      }
    `));class lE extends lO{connectedCallback(){super.connectedCallback(),iT.setValueFor(this,ib)}}let lR=lE.compose({baseName:"menu",baseClass:lO,template:(e,t)=>s4`
    <template
        slot="${e=>e.slot?e.slot:e.isNestedMenu()?"submenu":void 0}"
        role="menu"
        @keydown="${(e,t)=>e.handleMenuKeyDown(t.event)}"
        @focusout="${(e,t)=>e.handleFocusOut(t.event)}"
    >
        <slot ${rn("items")}></slot>
    </template>
`,styles:lD}),lL=(e,t)=>rh`
    ${ru("grid")} :host {
      contain: layout;
      overflow: visible;
      font-family: ${tu};
      outline: none;
      box-sizing: border-box;
      height: calc(${ry} * 1px);
      grid-template-columns: minmax(42px, auto) 1fr minmax(42px, auto);
      grid-template-rows: auto;
      justify-items: center;
      align-items: center;
      padding: 0;
      margin: 0 calc(${tv} * 1px);
      white-space: nowrap;
      background: ${i4};
      color: ${op};
      fill: currentcolor;
      cursor: pointer;
      font-size: ${tC};
      line-height: ${tI};
      border-radius: calc(${tb} * 1px);
      border: calc(${tk} * 1px) solid transparent;
    }

    :host(:hover) {
      position: relative;
      z-index: 1;
    }

    :host(.indent-0) {
      grid-template-columns: auto 1fr minmax(42px, auto);
    }
    :host(.indent-0) .content {
      grid-column: 1;
      grid-row: 1;
      margin-inline-start: 10px;
    }
    :host(.indent-0) .expand-collapse-glyph-container {
      grid-column: 5;
      grid-row: 1;
    }
    :host(.indent-2) {
      grid-template-columns:
        minmax(42px, auto) minmax(42px, auto) 1fr minmax(42px, auto)
        minmax(42px, auto);
    }
    :host(.indent-2) .content {
      grid-column: 3;
      grid-row: 1;
      margin-inline-start: 10px;
    }
    :host(.indent-2) .expand-collapse-glyph-container {
      grid-column: 5;
      grid-row: 1;
    }
    :host(.indent-2) .start {
      grid-column: 2;
    }
    :host(.indent-2) .end {
      grid-column: 4;
    }

    :host(:${rm}) {
      border-color: ${iE};
      background: ${i8};
      color: ${op};
    }

    :host(:hover) {
      background: ${i6};
      color: ${op};
    }

    :host(:active) {
      background: ${i9};
    }

    :host([aria-checked='true']),
    :host(.expanded) {
      background: ${iX};
      color: ${op};
    }

    :host([disabled]) {
      cursor: ${am};
      opacity: ${t$};
    }

    :host([disabled]:hover) {
      color: ${op};
      fill: currentcolor;
      background: ${i4};
    }

    :host([disabled]:hover) .start,
    :host([disabled]:hover) .end,
    :host([disabled]:hover)::slotted(svg) {
      fill: ${op};
    }

    .expand-collapse-glyph {
      /* TODO: adaptive typography https://github.com/microsoft/fast/issues/2432 */
      width: calc((16 + ${tf}) * 1px);
      height: calc((16 + ${tf}) * 1px);
      fill: currentcolor;
    }

    .content {
      grid-column-start: 2;
      justify-self: start;
      overflow: hidden;
      text-overflow: ellipsis;
    }

    .start,
    .end {
      display: flex;
      justify-content: center;
    }

    ::slotted(svg) {
      /* TODO: adaptive typography https://github.com/microsoft/fast/issues/2432 */
      width: 16px;
      height: 16px;

      /* Something like that would do if the typography is adaptive
      font-size: inherit;
      width: ${tD};
      height: ${tD};
      */
    }

    :host(:hover) .start,
    :host(:hover) .end,
    :host(:hover)::slotted(svg),
    :host(:active) .start,
    :host(:active) .end,
    :host(:active)::slotted(svg) {
      fill: ${op};
    }

    :host(.indent-0[aria-haspopup='menu']) {
      display: grid;
      grid-template-columns: minmax(42px, auto) auto 1fr minmax(42px, auto) minmax(
          42px,
          auto
        );
      align-items: center;
      min-height: 32px;
    }

    :host(.indent-1[aria-haspopup='menu']),
    :host(.indent-1[role='menuitemcheckbox']),
    :host(.indent-1[role='menuitemradio']) {
      display: grid;
      grid-template-columns: minmax(42px, auto) auto 1fr minmax(42px, auto) minmax(
          42px,
          auto
        );
      align-items: center;
      min-height: 32px;
    }

    :host(.indent-2:not([aria-haspopup='menu'])) .end {
      grid-column: 5;
    }

    :host .input-container,
    :host .expand-collapse-glyph-container {
      display: none;
    }

    :host([aria-haspopup='menu']) .expand-collapse-glyph-container,
    :host([role='menuitemcheckbox']) .input-container,
    :host([role='menuitemradio']) .input-container {
      display: grid;
      margin-inline-end: 10px;
    }

    :host([aria-haspopup='menu']) .content,
    :host([role='menuitemcheckbox']) .content,
    :host([role='menuitemradio']) .content {
      grid-column-start: 3;
    }

    :host([aria-haspopup='menu'].indent-0) .content {
      grid-column-start: 1;
    }

    :host([aria-haspopup='menu']) .end,
    :host([role='menuitemcheckbox']) .end,
    :host([role='menuitemradio']) .end {
      grid-column-start: 4;
    }

    :host .expand-collapse,
    :host .checkbox,
    :host .radio {
      display: flex;
      align-items: center;
      justify-content: center;
      position: relative;
      width: 20px;
      height: 20px;
      box-sizing: border-box;
      outline: none;
      margin-inline-start: 10px;
    }

    :host .checkbox,
    :host .radio {
      border: calc(${tw} * 1px) solid ${op};
    }

    :host([aria-checked='true']) .checkbox,
    :host([aria-checked='true']) .radio {
      background: ${iS};
      border-color: ${iS};
    }

    :host .checkbox {
      border-radius: calc(${tb} * 1px);
    }

    :host .radio {
      border-radius: 999px;
    }

    :host .checkbox-indicator,
    :host .radio-indicator,
    :host .expand-collapse-indicator,
    ::slotted([slot='checkbox-indicator']),
    ::slotted([slot='radio-indicator']),
    ::slotted([slot='expand-collapse-indicator']) {
      display: none;
    }

    ::slotted([slot='end']:not(svg)) {
      margin-inline-end: 10px;
      color: ${oc};
    }

    :host([aria-checked='true']) .checkbox-indicator,
    :host([aria-checked='true']) ::slotted([slot='checkbox-indicator']) {
      width: 100%;
      height: 100%;
      display: block;
      fill: ${iA};
      pointer-events: none;
    }

    :host([aria-checked='true']) .radio-indicator {
      position: absolute;
      top: 4px;
      left: 4px;
      right: 4px;
      bottom: 4px;
      border-radius: 999px;
      display: block;
      background: ${iA};
      pointer-events: none;
    }

    :host([aria-checked='true']) ::slotted([slot='radio-indicator']) {
      display: block;
      pointer-events: none;
    }
  `.withBehaviors(rv(rh`
      :host {
        border-color: transparent;
        color: ${nX.ButtonText};
        forced-color-adjust: none;
      }

      :host(:hover) {
        background: ${nX.Highlight};
        color: ${nX.HighlightText};
      }

      :host(:hover) .start,
      :host(:hover) .end,
      :host(:hover)::slotted(svg),
      :host(:active) .start,
      :host(:active) .end,
      :host(:active)::slotted(svg) {
        fill: ${nX.HighlightText};
      }

      :host(.expanded) {
        background: ${nX.Highlight};
        border-color: ${nX.Highlight};
        color: ${nX.HighlightText};
      }

      :host(:${rm}) {
        background: ${nX.Highlight};
        border-color: ${nX.ButtonText};
        box-shadow: 0 0 0 calc(${tk} * 1px) inset
          ${nX.HighlightText};
        color: ${nX.HighlightText};
        fill: currentcolor;
      }

      :host([disabled]),
      :host([disabled]:hover),
      :host([disabled]:hover) .start,
      :host([disabled]:hover) .end,
      :host([disabled]:hover)::slotted(svg) {
        background: ${nX.Canvas};
        color: ${nX.GrayText};
        fill: currentcolor;
        opacity: 1;
      }

      :host .expanded-toggle,
      :host .checkbox,
      :host .radio {
        border-color: ${nX.ButtonText};
        background: ${nX.HighlightText};
      }

      :host([checked='true']) .checkbox,
      :host([checked='true']) .radio {
        background: ${nX.HighlightText};
        border-color: ${nX.HighlightText};
      }

      :host(:hover) .expanded-toggle,
            :host(:hover) .checkbox,
            :host(:hover) .radio,
            :host(:${rm}) .expanded-toggle,
            :host(:${rm}) .checkbox,
            :host(:${rm}) .radio,
            :host([checked="true"]:hover) .checkbox,
            :host([checked="true"]:hover) .radio,
            :host([checked="true"]:${rm}) .checkbox,
            :host([checked="true"]:${rm}) .radio {
        border-color: ${nX.HighlightText};
      }

      :host([aria-checked='true']) {
        background: ${nX.Highlight};
        color: ${nX.HighlightText};
      }

      :host([aria-checked='true']) .checkbox-indicator,
      :host([aria-checked='true']) ::slotted([slot='checkbox-indicator']),
      :host([aria-checked='true']) ::slotted([slot='radio-indicator']) {
        fill: ${nX.Highlight};
      }

      :host([aria-checked='true']) .radio-indicator {
        background: ${nX.Highlight};
      }

      ::slotted([slot='end']:not(svg)) {
        color: ${nX.ButtonText};
      }

      :host(:hover) ::slotted([slot="end"]:not(svg)),
            :host(:${rm}) ::slotted([slot="end"]:not(svg)) {
        color: ${nX.HighlightText};
      }
    `),new rZ(rh`
        .expand-collapse-glyph {
          transform: rotate(0deg);
        }
      `,rh`
        .expand-collapse-glyph {
          transform: rotate(180deg);
        }
      `));class lA extends lS{}let lV=lA.compose({baseName:"menu-item",baseClass:lS,template:(e,t)=>s4`
    <template
        role="${e=>e.role}"
        aria-haspopup="${e=>e.hasSubmenu?"menu":void 0}"
        aria-checked="${e=>e.role!==lI?e.checked:void 0}"
        aria-disabled="${e=>e.disabled}"
        aria-expanded="${e=>e.expanded}"
        @keydown="${(e,t)=>e.handleMenuItemKeyDown(t.event)}"
        @click="${(e,t)=>e.handleMenuItemClick(t.event)}"
        @mouseover="${(e,t)=>e.handleMouseOver(t.event)}"
        @mouseout="${(e,t)=>e.handleMouseOut(t.event)}"
        class="${e=>e.disabled?"disabled":""} ${e=>e.expanded?"expanded":""} ${e=>`indent-${e.startColumnCount}`}"
    >
            ${rG(e=>e.role===lT,s4`
                    <div part="input-container" class="input-container">
                        <span part="checkbox" class="checkbox">
                            <slot name="checkbox-indicator">
                                ${t.checkboxIndicator||""}
                            </slot>
                        </span>
                    </div>
                `)}
            ${rG(e=>e.role===lF,s4`
                    <div part="input-container" class="input-container">
                        <span part="radio" class="radio">
                            <slot name="radio-indicator">
                                ${t.radioIndicator||""}
                            </slot>
                        </span>
                    </div>
                `)}
        </div>
        ${re(e,t)}
        <span class="content" part="content">
            <slot></slot>
        </span>
        ${s7(e,t)}
        ${rG(e=>e.hasSubmenu,s4`
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
        ${rG(e=>e.expanded,s4`
                <${e.tagFor(rq)}
                    :anchorElement="${e=>e}"
                    vertical-positioning-mode="dynamic"
                    vertical-default-position="bottom"
                    vertical-inset="true"
                    horizontal-positioning-mode="dynamic"
                    horizontal-default-position="end"
                    class="submenu-region"
                    dir="${e=>e.currentDirection}"
                    @loaded="${e=>e.submenuLoaded()}"
                    ${s9("submenuRegion")}
                    part="submenu-region"
                >
                    <slot name="submenu"></slot>
                </${e.tagFor(rq)}>
            `)}
    </template>
`,styles:lL,checkboxIndicator:`
        <svg
            part="checkbox-indicator"
            class="checkbox-indicator"
            viewBox="0 0 20 20"
            xmlns="http://www.w3.org/2000/svg"
        >
            <path
                fill-rule="evenodd"
                clip-rule="evenodd"
                d="M8.143 12.6697L15.235 4.5L16.8 5.90363L8.23812 15.7667L3.80005 11.2556L5.27591 9.7555L8.143 12.6697Z"
            />
        </svg>
    `,expandCollapseGlyph:`
        <svg
            viewBox="0 0 16 16"
            xmlns="http://www.w3.org/2000/svg"
            class="expand-collapse-glyph"
            part="expand-collapse-glyph"
        >
            <path
                d="M5.00001 12.3263C5.00124 12.5147 5.05566 12.699 5.15699 12.8578C5.25831 13.0167 5.40243 13.1437 5.57273 13.2242C5.74304 13.3047 5.9326 13.3354 6.11959 13.3128C6.30659 13.2902 6.4834 13.2152 6.62967 13.0965L10.8988 8.83532C11.0739 8.69473 11.2153 8.51658 11.3124 8.31402C11.4096 8.11146 11.46 7.88966 11.46 7.66499C11.46 7.44033 11.4096 7.21853 11.3124 7.01597C11.2153 6.81341 11.0739 6.63526 10.8988 6.49467L6.62967 2.22347C6.48274 2.10422 6.30501 2.02912 6.11712 2.00691C5.92923 1.9847 5.73889 2.01628 5.56823 2.09799C5.39757 2.17969 5.25358 2.30817 5.153 2.46849C5.05241 2.62882 4.99936 2.8144 5.00001 3.00369V12.3263Z"
            />
        </svg>
    `,radioIndicator:`
        <span part="radio-indicator" class="radio-indicator"></span>
    `});class lP extends sR{}class lH extends ah(lP){constructor(){super(...arguments),this.proxy=document.createElement("input")}}class lz extends lH{constructor(){super(...arguments),this.hideStep=!1,this.step=1,this.isUserInput=!1}maxChanged(e,t){var i;this.max=Math.max(t,null!=(i=this.min)?i:t);let o=Math.min(this.min,this.max);void 0!==this.min&&this.min!==o&&(this.min=o),this.value=this.getValidValue(this.value)}minChanged(e,t){var i;this.min=Math.min(t,null!=(i=this.max)?i:t);let o=Math.max(this.min,this.max);void 0!==this.max&&this.max!==o&&(this.max=o),this.value=this.getValidValue(this.value)}get valueAsNumber(){return parseFloat(super.value)}set valueAsNumber(e){this.value=e.toString()}valueChanged(e,t){this.value=this.getValidValue(t),t===this.value&&(this.control&&!this.isUserInput&&(this.control.value=this.value),super.valueChanged(e,this.value),void 0===e||this.isUserInput||(this.$emit("input"),this.$emit("change")),this.isUserInput=!1)}validate(){super.validate(this.control)}getValidValue(e){var t,i;let o=parseFloat(parseFloat(e).toPrecision(12));return isNaN(o)?"":Math.max(o=Math.min(o,null!=(t=this.max)?t:o),null!=(i=this.min)?i:o).toString()}stepUp(){let e=parseFloat(this.value),t=isNaN(e)?this.min>0?this.min:this.max<0?this.max:this.min?0:this.step:e+this.step;this.value=t.toString()}stepDown(){let e=parseFloat(this.value),t=isNaN(e)?this.min>0?this.min:this.max<0?this.max:this.min?0:0-this.step:e-this.step;this.value=t.toString()}connectedCallback(){super.connectedCallback(),this.proxy.setAttribute("type","number"),this.validate(),this.control.value=this.value,this.autofocus&&el.queueUpdate(()=>{this.focus()})}select(){this.control.select(),this.$emit("select")}handleTextInput(){this.control.value=this.control.value.replace(/[^0-9\-+e.]/g,""),this.isUserInput=!0,this.value=this.control.value}handleChange(){this.$emit("change")}handleKeyDown(e){switch(e.key){case o0:return this.stepUp(),!1;case oQ:return this.stepDown(),!1}return!0}handleBlur(){this.control.value=this.value}}X([eS({attribute:"readonly",mode:"boolean"})],lz.prototype,"readOnly",void 0),X([eS({mode:"boolean"})],lz.prototype,"autofocus",void 0),X([eS({attribute:"hide-step",mode:"boolean"})],lz.prototype,"hideStep",void 0),X([eS],lz.prototype,"placeholder",void 0),X([eS],lz.prototype,"list",void 0),X([eS({converter:eT})],lz.prototype,"maxlength",void 0),X([eS({converter:eT})],lz.prototype,"minlength",void 0),X([eS({converter:eT})],lz.prototype,"size",void 0),X([eS({converter:eT})],lz.prototype,"step",void 0),X([eS({converter:eT})],lz.prototype,"max",void 0),X([eS({converter:eT})],lz.prototype,"min",void 0),X([eu],lz.prototype,"defaultSlottedNodes",void 0),rt(lz,s8,nT);let lM=(e,t)=>rh`
  ${nR}

  .controls {
    opacity: 0;
    margin: auto 0;
  }

  .step-up-glyph,
  .step-down-glyph {
    display: block;
    padding: calc(
        (${tv} + 0.5 * ${tf} + 0.5 * ${ty}) * 1px
      )
      calc((2 + 2 * ${tv} + ${tf} + ${ty}) * 1px);
    cursor: pointer;
  }

  .step-up-glyph:before,
  .step-down-glyph:before {
    content: '';
    display: block;
    border: solid transparent
      calc((2 + ${tv} + 0.5 * ${tf} + 0.5 * ${ty}) * 1px);
  }

  .step-up-glyph:hover:before,
  .step-down-glyph:hover:before {
    background-color: ${iY};
  }

  .step-up-glyph:active:before,
  .step-down-glyph:active:before {
    background-color: ${iQ};
  }

  .step-up-glyph:before {
    border-bottom-color: ${op};
  }

  .step-down-glyph:before {
    border-top-color: ${op};
  }

  :host(:hover:not([disabled])) .controls,
  :host(:focus-within:not([disabled])) .controls {
    opacity: 1;
  }
`;class lN extends lz{constructor(){super(...arguments),this.appearance="outline"}}rI([eS],lN.prototype,"appearance",void 0);let lB=lN.compose({baseName:"number-field",baseClass:lz,styles:lM,template:(e,t)=>s4`
    <template class="${e=>e.readOnly?"readonly":""}">
        <label
            part="label"
            for="control"
            class="${e=>e.defaultSlottedNodes&&e.defaultSlottedNodes.length?"label":"label label__hidden"}"
        >
            <slot ${rn("defaultSlottedNodes")}></slot>
        </label>
        <div class="root" part="root">
            ${re(e,t)}
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
                ${s9("control")}
            />
            ${rG(e=>!e.hideStep&&!e.readOnly&&!e.disabled,s4`
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
            ${s7(e,t)}
        </div>
    </template>
`,shadowOptions:{delegatesFocus:!0},stepDownGlyph:`
        <span class="step-down-glyph" part="step-down-glyph"></span>
    `,stepUpGlyph:`
        <span class="step-up-glyph" part="step-up-glyph"></span>
    `}),lj=(e,t)=>rh`
    ${ru("inline-flex")} :host {
      align-items: center;
      font-family: ${tu};
      border-radius: calc(${tb} * 1px);
      border: calc(${tk} * 1px) solid transparent;
      box-sizing: border-box;
      background: ${i4};
      color: ${op};
      cursor: pointer;
      flex: 0 0 auto;
      fill: currentcolor;
      font-size: ${tC};
      height: calc(${ry} * 1px);
      line-height: ${tI};
      margin: 0 calc((${tv} - ${tk}) * 1px);
      outline: none;
      overflow: hidden;
      padding: 0 1ch;
      user-select: none;
      white-space: nowrap;
    }

    :host(:not([disabled]):not([aria-selected='true']):hover) {
      background: ${i6};
    }

    :host(:not([disabled]):not([aria-selected='true']):active) {
      background: ${i9};
    }

    :host([aria-selected='true']) {
      background: ${iS};
      color: ${iA};
    }

    :host(:not([disabled])[aria-selected='true']:hover) {
      background: ${iO};
      color: ${iV};
    }

    :host(:not([disabled])[aria-selected='true']:active) {
      background: ${iD};
      color: ${iP};
    }

    :host([disabled]) {
      cursor: ${am};
      opacity: ${t$};
    }

    .content {
      grid-column-start: 2;
      justify-self: start;
      overflow: hidden;
      text-overflow: ellipsis;
    }

    .start,
    .end,
    ::slotted(svg) {
      display: flex;
    }

    ::slotted(svg) {
      /* TODO: adaptive typography https://github.com/microsoft/fast/issues/2432 */
      height: calc(${tv} * 4px);
      width: calc(${tv} * 4px);
    }

    ::slotted([slot='end']) {
      margin-inline-start: 1ch;
    }

    ::slotted([slot='start']) {
      margin-inline-end: 1ch;
    }

    :host([aria-checked='true'][aria-selected='false']) {
      border-color: ${on};
    }

    :host([aria-checked='true'][aria-selected='true']) {
      border-color: ${on};
      box-shadow: 0 0 0 calc(${tk} * 2 * 1px) inset
        ${oh};
    }
  `.withBehaviors(rv(rh`
      :host {
        border-color: transparent;
        forced-color-adjust: none;
        color: ${nX.ButtonText};
        fill: currentcolor;
      }

      :host(:not([aria-selected='true']):hover),
      :host([aria-selected='true']) {
        background: ${nX.Highlight};
        color: ${nX.HighlightText};
      }

      :host([disabled]),
      :host([disabled][aria-selected='false']:hover) {
        background: ${nX.Canvas};
        color: ${nX.GrayText};
        fill: currentcolor;
        opacity: 1;
      }

      :host([aria-checked='true'][aria-selected='false']) {
        background: ${nX.ButtonFace};
        color: ${nX.ButtonText};
        border-color: ${nX.ButtonText};
      }

      :host([aria-checked='true'][aria-selected='true']),
      :host([aria-checked='true'][aria-selected='true']:hover) {
        background: ${nX.Highlight};
        color: ${nX.HighlightText};
        border-color: ${nX.ButtonText};
      }
    `));class lq extends aA{}let lU=lq.compose({baseName:"option",baseClass:aA,template:(e,t)=>s4`
    <template
        aria-checked="${e=>e.ariaChecked}"
        aria-disabled="${e=>e.ariaDisabled}"
        aria-posinset="${e=>e.ariaPosInSet}"
        aria-selected="${e=>e.ariaSelected}"
        aria-setsize="${e=>e.ariaSetSize}"
        class="${e=>[e.checked&&"checked",e.selected&&"selected",e.disabled&&"disabled"].filter(Boolean).join(" ")}"
        role="option"
    >
        ${re(e,t)}
        <span class="content" part="content">
            <slot ${rn("content")}></slot>
        </span>
        ${s7(e,t)}
    </template>
`,styles:lj});class l_ extends sR{constructor(){super(...arguments),this.percentComplete=0}valueChanged(){this.$fastController.isConnected&&this.updatePercentComplete()}minChanged(){this.$fastController.isConnected&&this.updatePercentComplete()}maxChanged(){this.$fastController.isConnected&&this.updatePercentComplete()}connectedCallback(){super.connectedCallback(),this.updatePercentComplete()}updatePercentComplete(){let e="number"==typeof this.min?this.min:0,t="number"==typeof this.max?this.max:100,i="number"==typeof this.value?this.value:0,o=t-e;this.percentComplete=0===o?0:Math.fround((i-e)/o*100)}}X([eS({converter:eT})],l_.prototype,"value",void 0),X([eS({converter:eT})],l_.prototype,"min",void 0),X([eS({converter:eT})],l_.prototype,"max",void 0),X([eS({mode:"boolean"})],l_.prototype,"paused",void 0),X([eu],l_.prototype,"percentComplete",void 0);let lG=(e,t)=>rh`
    ${ru("flex")} :host {
      align-items: center;
      outline: none;
      height: calc(${tv} * 1px);
      margin: calc(${tv} * 1px) 0;
    }

    .progress {
      background-color: ${iX};
      border-radius: calc(${tv} * 1px);
      width: 100%;
      height: 100%;
      display: flex;
      align-items: center;
      position: relative;
    }

    .determinate {
      background-color: ${iU};
      border-radius: calc(${tv} * 1px);
      height: 100%;
      transition: all 0.2s ease-in-out;
      display: flex;
    }

    .indeterminate {
      height: 100%;
      border-radius: calc(${tv} * 1px);
      display: flex;
      width: 100%;
      position: relative;
      overflow: hidden;
    }

    .indeterminate-indicator-1 {
      position: absolute;
      opacity: 0;
      height: 100%;
      background-color: ${iU};
      border-radius: calc(${tv} * 1px);
      animation-timing-function: cubic-bezier(0.4, 0, 0.6, 1);
      width: 40%;
      animation: indeterminate-1 2s infinite;
    }

    .indeterminate-indicator-2 {
      position: absolute;
      opacity: 0;
      height: 100%;
      background-color: ${iU};
      border-radius: calc(${tv} * 1px);
      animation-timing-function: cubic-bezier(0.4, 0, 0.6, 1);
      width: 60%;
      animation: indeterminate-2 2s infinite;
    }

    :host([paused]) .indeterminate-indicator-1,
    :host([paused]) .indeterminate-indicator-2 {
      animation-play-state: paused;
      background-color: ${iX};
    }

    :host([paused]) .determinate {
      background-color: ${oc};
    }

    @keyframes indeterminate-1 {
      0% {
        opacity: 1;
        transform: translateX(-100%);
      }
      70% {
        opacity: 1;
        transform: translateX(300%);
      }
      70.01% {
        opacity: 0;
      }
      100% {
        opacity: 0;
        transform: translateX(300%);
      }
    }

    @keyframes indeterminate-2 {
      0% {
        opacity: 0;
        transform: translateX(-150%);
      }
      29.99% {
        opacity: 0;
      }
      30% {
        opacity: 1;
        transform: translateX(-150%);
      }
      100% {
        transform: translateX(166.66%);
        opacity: 1;
      }
    }
  `.withBehaviors(rv(rh`
      .progress {
        forced-color-adjust: none;
        background-color: ${nX.Field};
        box-shadow: 0 0 0 1px inset ${nX.FieldText};
      }
      .determinate,
      .indeterminate-indicator-1,
      .indeterminate-indicator-2 {
        forced-color-adjust: none;
        background-color: ${nX.FieldText};
      }
      :host([paused]) .determinate,
      :host([paused]) .indeterminate-indicator-1,
      :host([paused]) .indeterminate-indicator-2 {
        background-color: ${nX.GrayText};
      }
    `));class lW extends l_{}let lK=lW.compose({baseName:"progress",baseClass:l_,template:(e,t)=>s4`
    <template
        role="progressbar"
        aria-valuenow="${e=>e.value}"
        aria-valuemin="${e=>e.min}"
        aria-valuemax="${e=>e.max}"
        class="${e=>e.paused?"paused":""}"
    >
        ${rG(e=>"number"==typeof e.value,s4`
                <div class="progress" part="progress" slot="determinate">
                    <div
                        class="determinate"
                        part="determinate"
                        style="width: ${e=>e.percentComplete}%"
                    ></div>
                </div>
            `,s4`
                <div class="progress" part="progress" slot="indeterminate">
                    <slot class="indeterminate" name="indeterminate">
                        ${t.indeterminateIndicator1||""}
                        ${t.indeterminateIndicator2||""}
                    </slot>
                </div>
            `)}
    </template>
`,styles:lG,indeterminateIndicator1:`
        <span class="indeterminate-indicator-1" part="indeterminate-indicator-1"></span>
    `,indeterminateIndicator2:`
        <span class="indeterminate-indicator-2" part="indeterminate-indicator-2"></span>
    `}),lX=(e,t)=>rh`
    ${ru("flex")} :host {
      align-items: center;
      outline: none;
      height: calc(${ry} * 1px);
      width: calc(${ry} * 1px);
      margin: calc(${ry} * 1px) 0;
    }

    .progress {
      height: 100%;
      width: 100%;
    }

    .background {
      stroke: ${iX};
      fill: none;
      stroke-width: 2px;
    }

    .determinate {
      stroke: ${iU};
      fill: none;
      stroke-width: 2px;
      stroke-linecap: round;
      transform-origin: 50% 50%;
      transform: rotate(-90deg);
      transition: all 0.2s ease-in-out;
    }

    .indeterminate-indicator-1 {
      stroke: ${iU};
      fill: none;
      stroke-width: 2px;
      stroke-linecap: round;
      transform-origin: 50% 50%;
      transform: rotate(-90deg);
      transition: all 0.2s ease-in-out;
      animation: spin-infinite 2s linear infinite;
    }

    :host([paused]) .indeterminate-indicator-1 {
      animation-play-state: paused;
      stroke: ${iX};
    }

    :host([paused]) .determinate {
      stroke: ${oc};
    }

    @keyframes spin-infinite {
      0% {
        stroke-dasharray: 0.01px 43.97px;
        transform: rotate(0deg);
      }
      50% {
        stroke-dasharray: 21.99px 21.99px;
        transform: rotate(450deg);
      }
      100% {
        stroke-dasharray: 0.01px 43.97px;
        transform: rotate(1080deg);
      }
    }
  `.withBehaviors(rv(rh`
      .indeterminate-indicator-1,
      .determinate {
        stroke: ${nX.FieldText};
      }
      .background {
        stroke: ${nX.Field};
      }
      :host([paused]) .indeterminate-indicator-1 {
        stroke: ${nX.Field};
      }
      :host([paused]) .determinate {
        stroke: ${nX.GrayText};
      }
    `));class lY extends l_{}let lQ=lY.compose({baseName:"progress-ring",baseClass:l_,template:(e,t)=>s4`
    <template
        role="progressbar"
        aria-valuenow="${e=>e.value}"
        aria-valuemin="${e=>e.min}"
        aria-valuemax="${e=>e.max}"
        class="${e=>e.paused?"paused":""}"
    >
        ${rG(e=>"number"==typeof e.value,s4`
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
            `,s4`
                <slot name="indeterminate" slot="indeterminate">
                    ${t.indeterminateIndicator||""}
                </slot>
            `)}
    </template>
`,styles:lX,indeterminateIndicator:`
        <svg class="progress" part="progress" viewBox="0 0 16 16">
            <circle
                class="background"
                part="background"
                cx="8px"
                cy="8px"
                r="7px"
            ></circle>
            <circle
                class="indeterminate-indicator-1"
                part="indeterminate-indicator-1"
                cx="8px"
                cy="8px"
                r="7px"
            ></circle>
        </svg>
    `});class lZ extends sR{}class lJ extends ad(lZ){constructor(){super(...arguments),this.proxy=document.createElement("input")}}class l0 extends lJ{constructor(){super(),this.initialValue="on",this.keypressHandler=e=>{if(" "===e.key){this.checked||this.readOnly||(this.checked=!0);return}return!0},this.proxy.setAttribute("type","radio")}readOnlyChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.readOnly=this.readOnly)}defaultCheckedChanged(){var e;!this.$fastController.isConnected||this.dirtyChecked||this.isInsideRadioGroup()||(this.checked=null!=(e=this.defaultChecked)&&e,this.dirtyChecked=!1)}connectedCallback(){var e,t;super.connectedCallback(),this.validate(),(null==(e=this.parentElement)?void 0:e.getAttribute("role"))==="radiogroup"||null!==this.getAttribute("tabindex")||this.disabled||this.setAttribute("tabindex","0"),!this.checkedAttribute||this.dirtyChecked||this.isInsideRadioGroup()||(this.checked=null!=(t=this.defaultChecked)&&t,this.dirtyChecked=!1)}isInsideRadioGroup(){return null!==this.closest("[role=radiogroup]")}clickHandler(e){this.disabled||this.readOnly||this.checked||(this.checked=!0)}}X([eS({attribute:"readonly",mode:"boolean"})],l0.prototype,"readOnly",void 0),X([eu],l0.prototype,"name",void 0),X([eu],l0.prototype,"defaultSlottedNodes",void 0);let l1=(e,t)=>rh`
    ${ru("inline-flex")} :host {
      --input-size: calc((${ry} / 2) + ${tv});
      align-items: center;
      outline: none;
      margin: calc(${tv} * 1px) 0;
      /* Chromium likes to select label text or the default slot when
                the radio is clicked. Maybe there is a better solution here? */
      user-select: none;
      position: relative;
      flex-direction: row;
      transition: all 0.2s ease-in-out;
    }

    .control {
      position: relative;
      width: calc((${ry} / 2 + ${tv}) * 1px);
      height: calc((${ry} / 2 + ${tv}) * 1px);
      box-sizing: border-box;
      border-radius: 999px;
      border: calc(${tw} * 1px) solid ${om};
      background: ${i0};
      outline: none;
      cursor: pointer;
    }

    :host([aria-invalid='true']) .control {
      border-color: ${oI};
    }

    .label {
      font-family: ${tu};
      color: ${op};
      padding-inline-start: calc(${tv} * 2px + 2px);
      margin-inline-end: calc(${tv} * 2px + 2px);
      cursor: pointer;
      font-size: ${tC};
      line-height: ${tI};
    }

    .label__hidden {
      display: none;
      visibility: hidden;
    }

    .control,
    .checked-indicator {
      flex-shrink: 0;
    }

    .checked-indicator {
      position: absolute;
      top: 5px;
      left: 5px;
      right: 5px;
      bottom: 5px;
      border-radius: 999px;
      display: inline-block;
      background: ${iA};
      fill: ${iA};
      opacity: 0;
      pointer-events: none;
    }

    :host(:not([disabled])) .control:hover {
      background: ${i1};
      border-color: ${ob};
    }

    :host([aria-invalid='true']:not([disabled])) .control:hover {
      border-color: ${oT};
    }

    :host(:not([disabled])) .control:active {
      background: ${i2};
      border-color: ${of};
    }

    :host([aria-invalid='true']:not([disabled])) .control:active {
      border-color: ${oF};
    }

    :host(:${rm}) .control {
      outline: solid calc(${tk} * 1px) ${iE};
    }

    :host([aria-invalid='true']:${rm}) .control {
      outline-color: ${oS};
    }

    :host([aria-checked='true']) .control {
      background: ${iS};
      border: calc(${tw} * 1px) solid ${iS};
    }

    :host([aria-invalid='true'][aria-checked='true']) .control {
      background-color: ${oI};
      border-color: ${oI};
    }

    :host([aria-checked='true']:not([disabled])) .control:hover {
      background: ${iO};
      border: calc(${tw} * 1px) solid ${iO};
    }

    :host([aria-invalid='true'][aria-checked='true']:not([disabled]))
      .control:hover {
      background-color: ${oT};
      border-color: ${oT};
    }

    :host([aria-checked='true']:not([disabled]))
      .control:hover
      .checked-indicator {
      background: ${iV};
      fill: ${iV};
    }

    :host([aria-checked='true']:not([disabled])) .control:active {
      background: ${iD};
      border: calc(${tw} * 1px) solid ${iD};
    }

    :host([aria-invalid='true'][aria-checked='true']:not([disabled]))
      .control:active {
      background-color: ${oF};
      border-color: ${oF};
    }

    :host([aria-checked='true']:not([disabled]))
      .control:active
      .checked-indicator {
      background: ${iP};
      fill: ${iP};
    }

    :host([aria-checked="true"]:${rm}:not([disabled])) .control {
      outline-offset: 2px;
      outline: solid calc(${tk} * 1px) ${iE};
    }

    :host([aria-invalid='true'][aria-checked="true"]:${rm}:not([disabled])) .control {
      outline-color: ${oS};
    }

    :host([disabled]) .label,
    :host([readonly]) .label,
    :host([readonly]) .control,
    :host([disabled]) .control {
      cursor: ${am};
    }

    :host([aria-checked='true']) .checked-indicator {
      opacity: 1;
    }

    :host([disabled]) {
      opacity: ${t$};
    }
  `.withBehaviors(rv(rh`
      .control,
      :host([aria-checked='true']:not([disabled])) .control {
        forced-color-adjust: none;
        border-color: ${nX.FieldText};
        background: ${nX.Field};
      }
      :host([aria-invalid='true']) {
        border-style: dashed;
      }
      :host(:not([disabled])) .control:hover {
        border-color: ${nX.Highlight};
        background: ${nX.Field};
      }
      :host([aria-checked='true']:not([disabled])) .control:hover,
      :host([aria-checked='true']:not([disabled])) .control:active {
        border-color: ${nX.Highlight};
        background: ${nX.Highlight};
      }
      :host([aria-checked='true']) .checked-indicator {
        background: ${nX.Highlight};
        fill: ${nX.Highlight};
      }
      :host([aria-checked='true']:not([disabled]))
        .control:hover
        .checked-indicator,
      :host([aria-checked='true']:not([disabled]))
        .control:active
        .checked-indicator {
        background: ${nX.HighlightText};
        fill: ${nX.HighlightText};
      }
      :host(:${rm}) .control {
        border-color: ${nX.Highlight};
        outline-offset: 2px;
        outline: solid calc(${tk} * 1px) ${nX.FieldText};
      }
      :host([aria-checked="true"]:${rm}:not([disabled])) .control {
        border-color: ${nX.Highlight};
        outline: solid calc(${tk} * 1px) ${nX.FieldText};
      }
      :host([disabled]) {
        forced-color-adjust: none;
        opacity: 1;
      }
      :host([disabled]) .label {
        color: ${nX.GrayText};
      }
      :host([disabled]) .control,
      :host([aria-checked='true'][disabled]) .control:hover,
      .control:active {
        background: ${nX.Field};
        border-color: ${nX.GrayText};
      }
      :host([disabled]) .checked-indicator,
      :host([aria-checked='true'][disabled]) .control:hover .checked-indicator {
        fill: ${nX.GrayText};
        background: ${nX.GrayText};
      }
    `)),l2=(e,t)=>s4`
  <template
    role="radio"
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
      <slot ${rn("defaultSlottedNodes")}></slot>
    </label>
  </template>
`;class l5 extends l0{}let l3=l5.compose({baseName:"radio",baseClass:l0,template:l2,styles:l1,checkedIndicator:`
        <div part="checked-indicator" class="checked-indicator"></div>
    `});class l4 extends sR{constructor(){super(...arguments),this.orientation=lf,this.radioChangeHandler=e=>{let t=e.target;t.checked&&(this.slottedRadioButtons.forEach(e=>{e!==t&&(e.checked=!1,this.isInsideFoundationToolbar||e.setAttribute("tabindex","-1"))}),this.selectedRadio=t,this.value=t.value,t.setAttribute("tabindex","0"),this.focusedRadio=t),e.stopPropagation()},this.moveToRadioByIndex=(e,t)=>{let i=e[t];this.isInsideToolbar||(i.setAttribute("tabindex","0"),i.readOnly?this.slottedRadioButtons.forEach(e=>{e!==i&&e.setAttribute("tabindex","-1")}):(i.checked=!0,this.selectedRadio=i)),this.focusedRadio=i,i.focus()},this.moveRightOffGroup=()=>{var e;null==(e=this.nextElementSibling)||e.focus()},this.moveLeftOffGroup=()=>{var e;null==(e=this.previousElementSibling)||e.focus()},this.focusOutHandler=e=>{let t=this.slottedRadioButtons,i=e.target,o=null!==i?t.indexOf(i):0,s=this.focusedRadio?t.indexOf(this.focusedRadio):-1;return(0===s&&o===s||s===t.length-1&&s===o)&&(this.selectedRadio?(this.focusedRadio=this.selectedRadio,this.isInsideFoundationToolbar||(this.selectedRadio.setAttribute("tabindex","0"),t.forEach(e=>{e!==this.selectedRadio&&e.setAttribute("tabindex","-1")}))):(this.focusedRadio=t[0],this.focusedRadio.setAttribute("tabindex","0"),t.forEach(e=>{e!==this.focusedRadio&&e.setAttribute("tabindex","-1")}))),!0},this.clickHandler=e=>{let t=e.target;if(t){let e=this.slottedRadioButtons;t.checked||0===e.indexOf(t)?(t.setAttribute("tabindex","0"),this.selectedRadio=t):(t.setAttribute("tabindex","-1"),this.selectedRadio=null),this.focusedRadio=t}e.preventDefault()},this.shouldMoveOffGroupToTheRight=(e,t,i)=>e===t.length&&this.isInsideToolbar&&i===oJ,this.shouldMoveOffGroupToTheLeft=(e,t)=>(this.focusedRadio?e.indexOf(this.focusedRadio)-1:0)<0&&this.isInsideToolbar&&t===oZ,this.checkFocusedRadio=()=>{null===this.focusedRadio||this.focusedRadio.readOnly||this.focusedRadio.checked||(this.focusedRadio.checked=!0,this.focusedRadio.setAttribute("tabindex","0"),this.focusedRadio.focus(),this.selectedRadio=this.focusedRadio)},this.moveRight=e=>{let t=this.slottedRadioButtons,i=0;if(i=this.focusedRadio?t.indexOf(this.focusedRadio)+1:1,this.shouldMoveOffGroupToTheRight(i,t,e.key))return void this.moveRightOffGroup();for(i===t.length&&(i=0);i<t.length&&t.length>1;)if(t[i].disabled)if(this.focusedRadio&&i===t.indexOf(this.focusedRadio))break;else if(i+1>=t.length)if(this.isInsideToolbar)break;else i=0;else i+=1;else{this.moveToRadioByIndex(t,i);break}},this.moveLeft=e=>{let t=this.slottedRadioButtons,i=0;if(i=(i=this.focusedRadio?t.indexOf(this.focusedRadio)-1:0)<0?t.length-1:i,this.shouldMoveOffGroupToTheLeft(t,e.key))return void this.moveLeftOffGroup();for(;i>=0&&t.length>1;)if(t[i].disabled)if(this.focusedRadio&&i===t.indexOf(this.focusedRadio))break;else i-1<0?i=t.length-1:i-=1;else{this.moveToRadioByIndex(t,i);break}},this.keydownHandler=e=>{let t=e.key;if(t in o1&&this.isInsideFoundationToolbar)return!0;switch(t){case"Enter":this.checkFocusedRadio();break;case oJ:case oQ:this.direction===nU.ltr?this.moveRight(e):this.moveLeft(e);break;case oZ:case o0:this.direction===nU.ltr?this.moveLeft(e):this.moveRight(e);break;default:return!0}}}readOnlyChanged(){void 0!==this.slottedRadioButtons&&this.slottedRadioButtons.forEach(e=>{this.readOnly?e.readOnly=!0:e.readOnly=!1})}disabledChanged(){void 0!==this.slottedRadioButtons&&this.slottedRadioButtons.forEach(e=>{this.disabled?e.disabled=!0:e.disabled=!1})}nameChanged(){this.slottedRadioButtons&&this.slottedRadioButtons.forEach(e=>{e.setAttribute("name",this.name)})}valueChanged(){this.slottedRadioButtons&&this.slottedRadioButtons.forEach(e=>{e.value===this.value&&(e.checked=!0,this.selectedRadio=e)}),this.$emit("change")}slottedRadioButtonsChanged(e,t){this.slottedRadioButtons&&this.slottedRadioButtons.length>0&&this.setupRadioButtons()}get parentToolbar(){return this.closest('[role="toolbar"]')}get isInsideToolbar(){var e;return null!=(e=this.parentToolbar)&&e}get isInsideFoundationToolbar(){var e;return!!(null==(e=this.parentToolbar)?void 0:e.$fastController)}connectedCallback(){super.connectedCallback(),this.direction=rj(this),this.setupRadioButtons()}disconnectedCallback(){this.slottedRadioButtons.forEach(e=>{e.removeEventListener("change",this.radioChangeHandler)})}setupRadioButtons(){let e=this.slottedRadioButtons.filter(e=>e.hasAttribute("checked")),t=e?e.length:0;t>1&&(e[t-1].checked=!0);let i=!1;if(this.slottedRadioButtons.forEach(e=>{void 0!==this.name&&e.setAttribute("name",this.name),this.disabled&&(e.disabled=!0),this.readOnly&&(e.readOnly=!0),this.value&&this.value===e.value?(this.selectedRadio=e,this.focusedRadio=e,e.checked=!0,e.setAttribute("tabindex","0"),i=!0):(this.isInsideFoundationToolbar||e.setAttribute("tabindex","-1"),e.checked=!1),e.addEventListener("change",this.radioChangeHandler)}),void 0===this.value&&this.slottedRadioButtons.length>0){let e=this.slottedRadioButtons.filter(e=>e.hasAttribute("checked")),t=null!==e?e.length:0;if(t>0&&!i){let i=e[t-1];i.checked=!0,this.focusedRadio=i,i.setAttribute("tabindex","0")}else this.slottedRadioButtons[0].setAttribute("tabindex","0"),this.focusedRadio=this.slottedRadioButtons[0]}}}X([eS({attribute:"readonly",mode:"boolean"})],l4.prototype,"readOnly",void 0),X([eS({attribute:"disabled",mode:"boolean"})],l4.prototype,"disabled",void 0),X([eS],l4.prototype,"name",void 0),X([eS],l4.prototype,"value",void 0),X([eS],l4.prototype,"orientation",void 0),X([eu],l4.prototype,"childItems",void 0),X([eu],l4.prototype,"slottedRadioButtons",void 0);let l6=(e,t)=>rh`
  ${ru("flex")} :host {
    align-items: flex-start;
    margin: calc(${tv} * 1px) 0;
    flex-direction: column;
  }
  .positioning-region {
    display: flex;
    flex-wrap: wrap;
  }
  :host([orientation='vertical']) .positioning-region {
    flex-direction: column;
  }
  :host([orientation='horizontal']) .positioning-region {
    flex-direction: row;
  }
`;class l9 extends l4{constructor(){super(),ec.getNotifier(this).subscribe({handleChange(e,t){"slottedRadioButtons"===t&&e.ariaInvalidChanged()}},"slottedRadioButtons")}ariaInvalidChanged(){this.slottedRadioButtons&&this.slottedRadioButtons.forEach(e=>{var t;e.setAttribute("aria-invalid",null!=(t=this.getAttribute("aria-invalid"))?t:"false")})}}let l8=l9.compose({baseName:"radio-group",baseClass:l4,template:(e,t)=>s4`
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
            class="positioning-region ${e=>e.orientation===lf?"horizontal":"vertical"}"
            part="positioning-region"
        >
            <slot
                ${rn({property:"slottedRadioButtons",filter:rs("[role=radio]")})}
            ></slot>
        </div>
    </template>
`,styles:l6});class l7 extends sR{}class he extends ah(l7){constructor(){super(...arguments),this.proxy=document.createElement("input")}}class ht extends he{readOnlyChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.readOnly=this.readOnly,this.validate())}autofocusChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.autofocus=this.autofocus,this.validate())}placeholderChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.placeholder=this.placeholder)}listChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.setAttribute("list",this.list),this.validate())}maxlengthChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.maxLength=this.maxlength,this.validate())}minlengthChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.minLength=this.minlength,this.validate())}patternChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.pattern=this.pattern,this.validate())}sizeChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.size=this.size)}spellcheckChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.spellcheck=this.spellcheck)}connectedCallback(){super.connectedCallback(),this.validate(),this.autofocus&&el.queueUpdate(()=>{this.focus()})}validate(){super.validate(this.control)}handleTextInput(){this.value=this.control.value}handleClearInput(){this.value="",this.control.focus(),this.handleChange()}handleChange(){this.$emit("change")}}X([eS({attribute:"readonly",mode:"boolean"})],ht.prototype,"readOnly",void 0),X([eS({mode:"boolean"})],ht.prototype,"autofocus",void 0),X([eS],ht.prototype,"placeholder",void 0),X([eS],ht.prototype,"list",void 0),X([eS({converter:eT})],ht.prototype,"maxlength",void 0),X([eS({converter:eT})],ht.prototype,"minlength",void 0),X([eS],ht.prototype,"pattern",void 0),X([eS({converter:eT})],ht.prototype,"size",void 0),X([eS({mode:"boolean"})],ht.prototype,"spellcheck",void 0),X([eu],ht.prototype,"defaultSlottedNodes",void 0);class hi{}rt(hi,rT),rt(ht,s8,hi);let ho=e2.create("clear-button-hover").withDefault(e=>{let t=i3.getValueFor(e),i=iK.getValueFor(e);return t.evaluate(e,i.evaluate(e).hover).hover}),hs=e2.create("clear-button-active").withDefault(e=>{let t=i3.getValueFor(e),i=iK.getValueFor(e);return t.evaluate(e,i.evaluate(e).hover).active}),hr=(e,t)=>rh`
  ${nR}

  .control::-webkit-search-cancel-button {
    -webkit-appearance: none;
  }

  .control:hover,
    .control:${rm},
    .control:disabled,
    .control:active {
    outline: none;
  }

  .clear-button {
    height: calc(100% - 2px);
    opacity: 0;
    margin: 1px;
    background: transparent;
    color: ${op};
    fill: currentcolor;
    border: none;
    border-radius: calc(${tb} * 1px);
    min-width: calc(${ry} * 1px);
    font-size: ${tC};
    line-height: ${tI};
    outline: none;
    font-family: ${tu};
    padding: 0 calc((10 + (${tv} * 2 * ${tf})) * 1px);
  }

  .clear-button:hover {
    background: ${i6};
  }

  .clear-button:active {
    background: ${i9};
  }

  :host([appearance='filled']) .clear-button:hover {
    background: ${ho};
  }

  :host([appearance='filled']) .clear-button:active {
    background: ${hs};
  }

  .input-wrapper {
    display: flex;
    position: relative;
    width: 100%;
  }

  .start,
  .end {
    display: flex;
    margin: 1px;
    fill: currentcolor;
  }

  ::slotted([slot='end']) {
    height: 100%;
  }

  .end {
    margin-inline-end: 1px;
    height: calc(100% - 2px);
  }

  ::slotted(svg) {
    /* TODO: adaptive typography https://github.com/microsoft/fast/issues/2432 */
    width: 16px;
    height: 16px;
    margin-inline-end: 11px;
    margin-inline-start: 11px;
    margin-top: auto;
    margin-bottom: auto;
  }

  .clear-button__hidden {
    opacity: 0;
  }

  :host(:hover:not([disabled], [readOnly])) .clear-button,
  :host(:active:not([disabled], [readOnly])) .clear-button,
  :host(:focus-within:not([disabled], [readOnly])) .clear-button {
    opacity: 1;
  }

  :host(:hover:not([disabled], [readOnly])) .clear-button__hidden,
  :host(:active:not([disabled], [readOnly])) .clear-button__hidden,
  :host(:focus-within:not([disabled], [readOnly])) .clear-button__hidden {
    opacity: 0;
  }
`;class ha extends ht{constructor(){super(...arguments),this.appearance="outline"}}rI([eS],ha.prototype,"appearance",void 0);let hn=ha.compose({baseName:"search",baseClass:ht,template:(e,t)=>s4`
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
                ${rn({property:"defaultSlottedNodes",filter:nA})}
            ></slot>
        </label>
        <div class="root" part="root" ${s9("root")}>
            ${re(e,t)}
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
                    ${s9("control")}
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
            ${s7(e,t)}
        </div>
    </template>
`,styles:hr,shadowOptions:{delegatesFocus:!0}});class hl extends aW{constructor(){super(...arguments),this.listboxScrollWidth=""}autoWidthChanged(e,t){t?this.setAutoWidth():this.style.removeProperty("width")}setAutoWidth(){if(!this.autoWidth||!this.isConnected)return;let e=this.listbox.getBoundingClientRect().width;0===e&&this.listbox.hidden&&(Object.assign(this.listbox.style,{visibility:"hidden"}),this.listbox.removeAttribute("hidden"),e=this.listbox.getBoundingClientRect().width,this.listbox.setAttribute("hidden",""),this.listbox.style.removeProperty("visibility")),e>0&&Object.assign(this.style,{width:`${e}px`})}connectedCallback(){super.connectedCallback(),this.setAutoWidth(),this.listbox&&iT.setValueFor(this.listbox,ib)}slottedOptionsChanged(e,t){super.slottedOptionsChanged(e,t),this.setAutoWidth()}get listboxMaxHeight(){return Math.floor(this.maxHeight/o$.getValueFor(this)).toString()}listboxScrollWidthChanged(){this.updateComputedStylesheet()}get selectSize(){var e;return`${null!=(e=this.size)?e:4*!!this.multiple}`}multipleChanged(e,t){super.multipleChanged(e,t),this.updateComputedStylesheet()}maxHeightChanged(e,t){this.collapsible&&this.updateComputedStylesheet()}setPositioning(){super.setPositioning(),this.updateComputedStylesheet()}sizeChanged(e,t){(super.sizeChanged(e,t),this.updateComputedStylesheet(),this.collapsible)?requestAnimationFrame(()=>{this.listbox.style.setProperty("display","flex"),this.listbox.style.setProperty("overflow","visible"),this.listbox.style.setProperty("visibility","hidden"),this.listbox.style.setProperty("width","auto"),this.listbox.hidden=!1,this.listboxScrollWidth=`${this.listbox.scrollWidth}`,this.listbox.hidden=!0,this.listbox.style.removeProperty("display"),this.listbox.style.removeProperty("overflow"),this.listbox.style.removeProperty("visibility"),this.listbox.style.removeProperty("width")}):this.listboxScrollWidth=""}updateComputedStylesheet(){this.computedStylesheet&&this.$fastController.removeStyles(this.computedStylesheet),this.computedStylesheet=rh`
      :host {
        --listbox-max-height: ${this.listboxMaxHeight};
        --listbox-scroll-width: ${this.listboxScrollWidth};
        --size: ${this.selectSize};
      }
    `,this.$fastController.addStyles(this.computedStylesheet)}}rI([eS({attribute:"autowidth",mode:"boolean"})],hl.prototype,"autoWidth",void 0),rI([eS({attribute:"minimal",mode:"boolean"})],hl.prototype,"minimal",void 0),rI([eS],hl.prototype,"scale",void 0),rI([eu],hl.prototype,"listboxScrollWidth",void 0);let hh=hl.compose({baseName:"select",baseClass:aW,template:(e,t)=>s4`
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
        ${rG(e=>e.collapsible,s4`
                <div
                    class="control"
                    part="control"
                    ?disabled="${e=>e.disabled}"
                    ${s9("control")}
                >
                    ${re(e,t)}
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
                    ${s7(e,t)}
                </div>
            `)}
        <div
            class="listbox"
            id="${e=>e.listboxId}"
            part="listbox"
            role="listbox"
            ?disabled="${e=>e.disabled}"
            ?hidden="${e=>!!e.collapsible&&!e.open}"
            ${s9("listbox")}
        >
            <slot
                ${rn({filter:aP.slottedOptionFilter,flatten:!0,property:"slottedOptions"})}
            ></slot>
        </div>
    </template>
`,styles:aY,indicator:`
        <svg
            class="select-indicator"
            part="select-indicator"
            viewBox="0 0 12 7"
            xmlns="http://www.w3.org/2000/svg"
        >
            <path
                d="M11.85.65c.2.2.2.5 0 .7L6.4 6.84a.55.55 0 01-.78 0L.14 1.35a.5.5 0 11.71-.7L6 5.8 11.15.65c.2-.2.5-.2.7 0z"
            />
        </svg>
    `});class hd extends sR{constructor(){super(...arguments),this.shape="rect"}}X([eS],hd.prototype,"fill",void 0),X([eS],hd.prototype,"shape",void 0),X([eS],hd.prototype,"pattern",void 0),X([eS({mode:"boolean"})],hd.prototype,"shimmer",void 0);let hc=(e,t)=>rh`
    ${ru("block")} :host {
      --skeleton-fill-default: #e1dfdd;
      overflow: hidden;
      width: 100%;
      position: relative;
      background-color: var(--skeleton-fill, var(--skeleton-fill-default));
      --skeleton-animation-gradient-default: linear-gradient(
        270deg,
        var(--skeleton-fill, var(--skeleton-fill-default)) 0%,
        #f3f2f1 51.13%,
        var(--skeleton-fill, var(--skeleton-fill-default)) 100%
      );
      --skeleton-animation-timing-default: ease-in-out;
    }

    :host([shape='rect']) {
      border-radius: calc(${tb} * 1px);
    }

    :host([shape='circle']) {
      border-radius: 100%;
      overflow: hidden;
    }

    object {
      position: absolute;
      width: 100%;
      height: auto;
      z-index: 2;
    }

    object img {
      width: 100%;
      height: auto;
    }

    ${ru("block")} span.shimmer {
      position: absolute;
      width: 100%;
      height: 100%;
      background-image: var(
        --skeleton-animation-gradient,
        var(--skeleton-animation-gradient-default)
      );
      background-size: 0px 0px / 90% 100%;
      background-repeat: no-repeat;
      background-color: var(--skeleton-animation-fill, ${iX});
      animation: shimmer 2s infinite;
      animation-timing-function: var(
        --skeleton-animation-timing,
        var(--skeleton-timing-default)
      );
      animation-direction: normal;
      z-index: 1;
    }

    ::slotted(svg) {
      z-index: 2;
    }

    ::slotted(.pattern) {
      width: 100%;
      height: 100%;
    }

    @keyframes shimmer {
      0% {
        transform: translateX(-100%);
      }
      100% {
        transform: translateX(100%);
      }
    }
  `.withBehaviors(rv(rh`
      :host {
        forced-color-adjust: none;
        background-color: ${nX.ButtonFace};
        box-shadow: 0 0 0 1px ${nX.ButtonText};
      }

      ${ru("block")} span.shimmer {
        display: none;
      }
    `));class hu extends hd{}let hp=hu.compose({baseName:"skeleton",baseClass:hd,template:(e,t)=>s4`
    <template
        class="${e=>"circle"===e.shape?"circle":"rect"}"
        pattern="${e=>e.pattern}"
        ?shimmer="${e=>e.shimmer}"
    >
        ${rG(e=>!0===e.shimmer,s4`
                <span class="shimmer"></span>
            `)}
        <object type="image/svg+xml" data="${e=>e.pattern}" role="presentation">
            <img class="pattern" src="${e=>e.pattern}" />
        </object>
        <slot></slot>
    </template>
`,styles:hc});function hg(e,t,i,o){let s=o2(0,1,(e-t)/(i-t));return o===nU.rtl&&(s=1-s),s}class hm extends sR{}class hb extends ah(hm){constructor(){super(...arguments),this.proxy=document.createElement("input")}}class hf extends hb{constructor(){super(...arguments),this.direction=nU.ltr,this.isDragging=!1,this.trackWidth=0,this.trackMinWidth=0,this.trackHeight=0,this.trackLeft=0,this.trackMinHeight=0,this.valueTextFormatter=()=>null,this.min=0,this.max=10,this.step=1,this.orientation=lf,this.mode="single-value",this.keypressHandler=e=>{if(!this.readOnly){if("Home"===e.key)e.preventDefault(),this.value=`${this.min}`;else if("End"===e.key)e.preventDefault(),this.value=`${this.max}`;else if(!e.shiftKey)switch(e.key){case oJ:case o0:e.preventDefault(),this.increment();break;case oZ:case oQ:e.preventDefault(),this.decrement()}}},this.setupTrackConstraints=()=>{let e=this.track.getBoundingClientRect();this.trackWidth=this.track.clientWidth,this.trackMinWidth=this.track.clientLeft,this.trackHeight=e.bottom,this.trackMinHeight=e.top,this.trackLeft=this.getBoundingClientRect().left,0===this.trackWidth&&(this.trackWidth=1)},this.setupListeners=(e=!1)=>{let t=`${e?"remove":"add"}EventListener`;this[t]("keydown",this.keypressHandler),this[t]("mousedown",this.handleMouseDown),this.thumb[t]("mousedown",this.handleThumbMouseDown,{passive:!0}),this.thumb[t]("touchstart",this.handleThumbMouseDown,{passive:!0}),e&&(this.handleMouseDown(null),this.handleThumbMouseDown(null))},this.initialValue="",this.handleThumbMouseDown=e=>{if(e){if(this.readOnly||this.disabled||e.defaultPrevented)return;e.target.focus()}let t=`${null!==e?"add":"remove"}EventListener`;window[t]("mouseup",this.handleWindowMouseUp),window[t]("mousemove",this.handleMouseMove,{passive:!0}),window[t]("touchmove",this.handleMouseMove,{passive:!0}),window[t]("touchend",this.handleWindowMouseUp),this.isDragging=null!==e},this.handleMouseMove=e=>{if(this.readOnly||this.disabled||e.defaultPrevented)return;let t=window.TouchEvent&&e instanceof TouchEvent?e.touches[0]:e,i=this.orientation===lf?t.pageX-document.documentElement.scrollLeft-this.trackLeft:t.pageY-document.documentElement.scrollTop;this.value=`${this.calculateNewValue(i)}`},this.calculateNewValue=e=>{let t=hg(e,this.orientation===lf?this.trackMinWidth:this.trackMinHeight,this.orientation===lf?this.trackWidth:this.trackHeight,this.direction),i=(this.max-this.min)*t+this.min;return this.convertToConstrainedValue(i)},this.handleWindowMouseUp=e=>{this.stopDragging()},this.stopDragging=()=>{this.isDragging=!1,this.handleMouseDown(null),this.handleThumbMouseDown(null)},this.handleMouseDown=e=>{let t=`${null!==e?"add":"remove"}EventListener`;if((null===e||!this.disabled&&!this.readOnly)&&(window[t]("mouseup",this.handleWindowMouseUp),window.document[t]("mouseleave",this.handleWindowMouseUp),window[t]("mousemove",this.handleMouseMove),e)){e.preventDefault(),this.setupTrackConstraints(),e.target.focus();let t=this.orientation===lf?e.pageX-document.documentElement.scrollLeft-this.trackLeft:e.pageY-document.documentElement.scrollTop;this.value=`${this.calculateNewValue(t)}`}},this.convertToConstrainedValue=e=>{isNaN(e)&&(e=this.min);let t=e-this.min,i=Math.round(t/this.step),o=t-i*(this.stepMultiplier*this.step)/this.stepMultiplier;return(t=o>=Number(this.step)/2?t-o+Number(this.step):t-o)+this.min}}readOnlyChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.readOnly=this.readOnly)}get valueAsNumber(){return parseFloat(super.value)}set valueAsNumber(e){this.value=e.toString()}valueChanged(e,t){super.valueChanged(e,t),this.$fastController.isConnected&&this.setThumbPositionForOrientation(this.direction),this.$emit("change")}minChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.min=`${this.min}`),this.validate()}maxChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.max=`${this.max}`),this.validate()}stepChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.step=`${this.step}`),this.updateStepMultiplier(),this.validate()}orientationChanged(){this.$fastController.isConnected&&this.setThumbPositionForOrientation(this.direction)}connectedCallback(){super.connectedCallback(),this.proxy.setAttribute("type","range"),this.direction=rj(this),this.updateStepMultiplier(),this.setupTrackConstraints(),this.setupListeners(),this.setupDefaultValue(),this.setThumbPositionForOrientation(this.direction)}disconnectedCallback(){this.setupListeners(!0)}increment(){let e=this.direction!==nU.rtl&&this.orientation!==lv?Number(this.value)+Number(this.step):Number(this.value)-Number(this.step),t=this.convertToConstrainedValue(e),i=t<Number(this.max)?`${t}`:`${this.max}`;this.value=i}decrement(){let e=this.direction!==nU.rtl&&this.orientation!==lv?Number(this.value)-Number(this.step):Number(this.value)+Number(this.step),t=this.convertToConstrainedValue(e),i=t>Number(this.min)?`${t}`:`${this.min}`;this.value=i}setThumbPositionForOrientation(e){let t=(1-hg(Number(this.value),Number(this.min),Number(this.max),e))*100;this.orientation===lf?this.position=this.isDragging?`right: ${t}%; transition: none;`:`right: ${t}%; transition: all 0.2s ease;`:this.position=this.isDragging?`bottom: ${t}%; transition: none;`:`bottom: ${t}%; transition: all 0.2s ease;`}updateStepMultiplier(){let e=this.step+"",t=this.step%1?e.length-e.indexOf(".")-1:0;this.stepMultiplier=Math.pow(10,t)}get midpoint(){return`${this.convertToConstrainedValue((this.max+this.min)/2)}`}setupDefaultValue(){if("string"==typeof this.value)if(0===this.value.length)this.initialValue=this.midpoint;else{let e=parseFloat(this.value);!Number.isNaN(e)&&(e<this.min||e>this.max)&&(this.value=this.midpoint)}}}X([eS({attribute:"readonly",mode:"boolean"})],hf.prototype,"readOnly",void 0),X([eu],hf.prototype,"direction",void 0),X([eu],hf.prototype,"isDragging",void 0),X([eu],hf.prototype,"position",void 0),X([eu],hf.prototype,"trackWidth",void 0),X([eu],hf.prototype,"trackMinWidth",void 0),X([eu],hf.prototype,"trackHeight",void 0),X([eu],hf.prototype,"trackLeft",void 0),X([eu],hf.prototype,"trackMinHeight",void 0),X([eu],hf.prototype,"valueTextFormatter",void 0),X([eS({converter:eT})],hf.prototype,"min",void 0),X([eS({converter:eT})],hf.prototype,"max",void 0),X([eS({converter:eT})],hf.prototype,"step",void 0),X([eS],hf.prototype,"orientation",void 0),X([eS],hf.prototype,"mode",void 0);let hv=rh`
  .track-start {
    left: 0;
  }
`,hy=rh`
  .track-start {
    right: 0;
  }
`,hx=(e,t)=>rh`
    :host([hidden]) {
      display: none;
    }

    ${ru("inline-grid")} :host {
      --thumb-size: calc(${ry} * 0.5 - ${tv});
      --thumb-translate: calc(
        var(--thumb-size) * -0.5 + var(--track-width) / 2
      );
      --track-overhang: calc((${tv} / 2) * -1);
      --track-width: ${tv};
      --jp-slider-height: calc(var(--thumb-size) * 10);
      align-items: center;
      width: 100%;
      margin: calc(${tv} * 1px) 0;
      user-select: none;
      box-sizing: border-box;
      border-radius: calc(${tb} * 1px);
      outline: none;
      cursor: pointer;
    }
    :host([orientation='horizontal']) .positioning-region {
      position: relative;
      margin: 0 8px;
      display: grid;
      grid-template-rows: calc(var(--thumb-size) * 1px) 1fr;
    }
    :host([orientation='vertical']) .positioning-region {
      position: relative;
      margin: 0 8px;
      display: grid;
      height: 100%;
      grid-template-columns: calc(var(--thumb-size) * 1px) 1fr;
    }

    :host(:${rm}) .thumb-cursor {
      box-shadow:
        0 0 0 2px ${iT},
        0 0 0 calc((2 + ${tk}) * 1px) ${iE};
    }

    :host([aria-invalid='true']:${rm}) .thumb-cursor {
      box-shadow:
        0 0 0 2px ${iT},
        0 0 0 calc((2 + ${tk}) * 1px) ${oS};
    }

    .thumb-container {
      position: absolute;
      height: calc(var(--thumb-size) * 1px);
      width: calc(var(--thumb-size) * 1px);
      transition: all 0.2s ease;
      color: ${op};
      fill: currentcolor;
    }
    .thumb-cursor {
      border: none;
      width: calc(var(--thumb-size) * 1px);
      height: calc(var(--thumb-size) * 1px);
      background: ${op};
      border-radius: calc(${tb} * 1px);
    }
    .thumb-cursor:hover {
      background: ${op};
      border-color: ${ob};
    }
    .thumb-cursor:active {
      background: ${op};
    }
    .track-start {
      background: ${iU};
      position: absolute;
      height: 100%;
      left: 0;
      border-radius: calc(${tb} * 1px);
    }
    :host([aria-invalid='true']) .track-start {
      background-color: ${oI};
    }
    :host([orientation='horizontal']) .thumb-container {
      transform: translateX(calc(var(--thumb-size) * 0.5px))
        translateY(calc(var(--thumb-translate) * 1px));
    }
    :host([orientation='vertical']) .thumb-container {
      transform: translateX(calc(var(--thumb-translate) * 1px))
        translateY(calc(var(--thumb-size) * 0.5px));
    }
    :host([orientation='horizontal']) {
      min-width: calc(var(--thumb-size) * 1px);
    }
    :host([orientation='horizontal']) .track {
      right: calc(var(--track-overhang) * 1px);
      left: calc(var(--track-overhang) * 1px);
      align-self: start;
      height: calc(var(--track-width) * 1px);
    }
    :host([orientation='vertical']) .track {
      top: calc(var(--track-overhang) * 1px);
      bottom: calc(var(--track-overhang) * 1px);
      width: calc(var(--track-width) * 1px);
      height: 100%;
    }
    .track {
      background: ${om};
      position: absolute;
      border-radius: calc(${tb} * 1px);
    }
    :host([orientation='vertical']) {
      height: calc(var(--fast-slider-height) * 1px);
      min-height: calc(var(--thumb-size) * 1px);
      min-width: calc(${tv} * 20px);
    }
    :host([orientation='vertical']) .track-start {
      height: auto;
      width: 100%;
      top: 0;
    }
    :host([disabled]),
    :host([readonly]) {
      cursor: ${am};
    }
    :host([disabled]) {
      opacity: ${t$};
    }
  `.withBehaviors(new rZ(hv,hy),rv(rh`
      .thumb-cursor {
        forced-color-adjust: none;
        border-color: ${nX.FieldText};
        background: ${nX.FieldText};
      }
      .thumb-cursor:hover,
      .thumb-cursor:active {
        background: ${nX.Highlight};
      }
      .track {
        forced-color-adjust: none;
        background: ${nX.FieldText};
      }
      :host(:${rm}) .thumb-cursor {
        border-color: ${nX.Highlight};
      }
      :host([disabled]) {
        opacity: 1;
      }
      :host([disabled]) .track,
      :host([disabled]) .thumb-cursor {
        forced-color-adjust: none;
        background: ${nX.GrayText};
      }

      :host(:${rm}) .thumb-cursor {
        background: ${nX.Highlight};
        border-color: ${nX.Highlight};
        box-shadow:
          0 0 0 2px ${nX.Field},
          0 0 0 4px ${nX.FieldText};
      }
    `));class h$ extends hf{}let hw=h$.compose({baseName:"slider",baseClass:hf,template:(e,t)=>s4`
    <template
        role="slider"
        class="${e=>e.readOnly?"readonly":""}
        ${e=>e.orientation||lf}"
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
            <div ${s9("track")} part="track-container" class="track">
                <slot name="track"></slot>
                <div part="track-start" class="track-start" style="${e=>e.position}">
                    <slot name="track-start"></slot>
                </div>
            </div>
            <slot></slot>
            <div
                ${s9("thumb")}
                part="thumb-container"
                class="thumb-container"
                style="${e=>e.position}"
            >
                <slot name="thumb">${t.thumb||""}</slot>
            </div>
        </div>
    </template>
`,styles:hx,thumb:`
        <div class="thumb-cursor"></div>
    `}),hk={min:0,max:0,direction:nU.ltr,orientation:lf,disabled:!1};class hC extends sR{constructor(){super(...arguments),this.hideMark=!1,this.sliderDirection=nU.ltr,this.getSliderConfiguration=()=>{if(this.isSliderConfig(this.parentNode)){let{min:e,max:t,direction:i,orientation:o,disabled:s}=this.parentNode;void 0!==s&&(this.disabled=s),this.sliderDirection=i||nU.ltr,this.sliderOrientation=o||lf,this.sliderMaxPosition=t,this.sliderMinPosition=e}else this.sliderDirection=hk.direction||nU.ltr,this.sliderOrientation=hk.orientation||lf,this.sliderMaxPosition=hk.max,this.sliderMinPosition=hk.min},this.positionAsStyle=()=>{let e=this.sliderDirection?this.sliderDirection:nU.ltr,t=hg(Number(this.position),Number(this.sliderMinPosition),Number(this.sliderMaxPosition)),i=Math.round((1-t)*100),o=Math.round(100*t);return(Number.isNaN(o)&&Number.isNaN(i)&&(i=50,o=50),this.sliderOrientation===lf)?e===nU.rtl?`right: ${o}%; left: ${i}%;`:`left: ${o}%; right: ${i}%;`:`top: ${o}%; bottom: ${i}%;`}}positionChanged(){this.positionStyle=this.positionAsStyle()}sliderOrientationChanged(){}connectedCallback(){super.connectedCallback(),this.getSliderConfiguration(),this.positionStyle=this.positionAsStyle(),this.notifier=ec.getNotifier(this.parentNode),this.notifier.subscribe(this,"orientation"),this.notifier.subscribe(this,"direction"),this.notifier.subscribe(this,"max"),this.notifier.subscribe(this,"min")}disconnectedCallback(){super.disconnectedCallback(),this.notifier.unsubscribe(this,"orientation"),this.notifier.unsubscribe(this,"direction"),this.notifier.unsubscribe(this,"max"),this.notifier.unsubscribe(this,"min")}handleChange(e,t){switch(t){case"direction":this.sliderDirection=e.direction;break;case"orientation":this.sliderOrientation=e.orientation;break;case"max":this.sliderMaxPosition=e.max;break;case"min":this.sliderMinPosition=e.min}this.positionStyle=this.positionAsStyle()}isSliderConfig(e){return void 0!==e.max&&void 0!==e.min}}X([eu],hC.prototype,"positionStyle",void 0),X([eS],hC.prototype,"position",void 0),X([eS({attribute:"hide-mark",mode:"boolean"})],hC.prototype,"hideMark",void 0),X([eS({attribute:"disabled",mode:"boolean"})],hC.prototype,"disabled",void 0),X([eu],hC.prototype,"sliderOrientation",void 0),X([eu],hC.prototype,"sliderMinPosition",void 0),X([eu],hC.prototype,"sliderMaxPosition",void 0),X([eu],hC.prototype,"sliderDirection",void 0);let hI=rh`
  :host {
    align-self: start;
    grid-row: 2;
    margin-top: -2px;
    height: calc((${ry} / 2 + ${tv}) * 1px);
    width: auto;
  }
  .container {
    grid-template-rows: auto auto;
    grid-template-columns: 0;
  }
  .label {
    margin: 2px 0;
  }
`,hT=rh`
  :host {
    justify-self: start;
    grid-column: 2;
    margin-left: 2px;
    height: auto;
    width: calc((${ry} / 2 + ${tv}) * 1px);
  }
  .container {
    grid-template-columns: auto auto;
    grid-template-rows: 0;
    min-width: calc(var(--thumb-size) * 1px);
    height: calc(var(--thumb-size) * 1px);
  }
  .mark {
    transform: rotate(90deg);
    align-self: center;
  }
  .label {
    margin-left: calc((${tv} / 2) * 3px);
    align-self: center;
  }
`,hF=(e,t)=>rh`
    ${ru("block")} :host {
      font-family: ${tu};
      color: ${op};
      fill: currentcolor;
    }
    .root {
      position: absolute;
      display: grid;
    }
    .container {
      display: grid;
      justify-self: center;
    }
    .label {
      justify-self: center;
      align-self: center;
      white-space: nowrap;
      max-width: 30px;
    }
    .mark {
      width: calc((${tv} / 4) * 1px);
      height: calc(${ry} * 0.25 * 1px);
      background: ${om};
      justify-self: center;
    }
    :host(.disabled) {
      opacity: ${t$};
    }
  `.withBehaviors(rv(rh`
      .mark {
        forced-color-adjust: none;
        background: ${nX.FieldText};
      }
      :host(.disabled) {
        forced-color-adjust: none;
        opacity: 1;
      }
      :host(.disabled) .label {
        color: ${nX.GrayText};
      }
      :host(.disabled) .mark {
        background: ${nX.GrayText};
      }
    `));class hS extends hC{sliderOrientationChanged(){this.sliderOrientation===lf?(this.$fastController.addStyles(hI),this.$fastController.removeStyles(hT)):(this.$fastController.addStyles(hT),this.$fastController.removeStyles(hI))}}let hO=hS.compose({baseName:"slider-label",baseClass:hC,template:(e,t)=>s4`
    <template
        aria-disabled="${e=>e.disabled}"
        class="${e=>e.sliderOrientation||lf}
            ${e=>e.disabled?"disabled":""}"
    >
        <div ${s9("root")} part="root" class="root" style="${e=>e.positionStyle}">
            <div class="container">
                ${rG(e=>!e.hideMark,s4`
                        <div class="mark"></div>
                    `)}
                <div class="label">
                    <slot></slot>
                </div>
            </div>
        </div>
    </template>
`,styles:hF});class hD extends sR{}class hE extends ad(hD){constructor(){super(...arguments),this.proxy=document.createElement("input")}}class hR extends hE{constructor(){super(),this.initialValue="on",this.keypressHandler=e=>{if(!this.readOnly)switch(e.key){case"Enter":case" ":this.checked=!this.checked}},this.clickHandler=e=>{this.disabled||this.readOnly||(this.checked=!this.checked)},this.proxy.setAttribute("type","checkbox")}readOnlyChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.readOnly=this.readOnly),this.readOnly?this.classList.add("readonly"):this.classList.remove("readonly")}checkedChanged(e,t){super.checkedChanged(e,t),this.checked?this.classList.add("checked"):this.classList.remove("checked")}}X([eS({attribute:"readonly",mode:"boolean"})],hR.prototype,"readOnly",void 0),X([eu],hR.prototype,"defaultSlottedNodes",void 0);let hL=(e,t)=>rh`
    :host([hidden]) {
      display: none;
    }

    ${ru("inline-flex")} :host {
      align-items: center;
      outline: none;
      font-family: ${tu};
      margin: calc(${tv} * 1px) 0;
      ${""} user-select: none;
    }

    :host([disabled]) {
      opacity: ${t$};
    }

    :host([disabled]) .label,
    :host([readonly]) .label,
    :host([readonly]) .switch,
    :host([disabled]) .switch {
      cursor: ${am};
    }

    .switch {
      position: relative;
      outline: none;
      box-sizing: border-box;
      width: calc(${ry} * 1px);
      height: calc((${ry} / 2 + ${tv}) * 1px);
      background: ${i0};
      border-radius: calc(${tb} * 1px);
      border: calc(${tw} * 1px) solid ${om};
    }

    :host([aria-invalid='true']) .switch {
      border-color: ${oI};
    }

    .switch:hover {
      background: ${i1};
      border-color: ${ob};
      cursor: pointer;
    }

    :host([disabled]) .switch:hover,
    :host([readonly]) .switch:hover {
      background: ${i1};
      border-color: ${ob};
      cursor: ${am};
    }

    :host([aria-invalid='true'][disabled]) .switch:hover,
    :host([aria-invalid='true'][readonly]) .switch:hover {
      border-color: ${oT};
    }

    :host(:not([disabled])) .switch:active {
      background: ${i2};
      border-color: ${of};
    }

    :host([aria-invalid='true']:not([disabled])) .switch:active {
      border-color: ${oF};
    }

    :host(:${rm}) .switch {
      outline-offset: 2px;
      outline: solid calc(${tk} * 1px) ${iE};
    }

    :host([aria-invalid='true']:${rm}) .switch {
      outline-color: ${oS};
    }

    .checked-indicator {
      position: absolute;
      top: 5px;
      bottom: 5px;
      background: ${op};
      border-radius: calc(${tb} * 1px);
      transition: all 0.2s ease-in-out;
    }

    .status-message {
      color: ${op};
      cursor: pointer;
      font-size: ${tC};
      line-height: ${tI};
    }

    :host([disabled]) .status-message,
    :host([readonly]) .status-message {
      cursor: ${am};
    }

    .label {
      color: ${op};
      margin-inline-end: calc(${tv} * 2px + 2px);
      font-size: ${tC};
      line-height: ${tI};
      cursor: pointer;
    }

    .label__hidden {
      display: none;
      visibility: hidden;
    }

    ::slotted([slot='checked-message']),
    ::slotted([slot='unchecked-message']) {
      margin-inline-start: calc(${tv} * 2px + 2px);
    }

    :host([aria-checked='true']) .checked-indicator {
      background: ${iA};
    }

    :host([aria-checked='true']) .switch {
      background: ${iS};
      border-color: ${iS};
    }

    :host([aria-checked='true']:not([disabled])) .switch:hover {
      background: ${iO};
      border-color: ${iO};
    }

    :host([aria-invalid='true'][aria-checked='true']) .switch {
      background-color: ${oI};
      border-color: ${oI};
    }

    :host([aria-invalid='true'][aria-checked='true']:not([disabled]))
      .switch:hover {
      background-color: ${oT};
      border-color: ${oT};
    }

    :host([aria-checked='true']:not([disabled]))
      .switch:hover
      .checked-indicator {
      background: ${iV};
    }

    :host([aria-checked='true']:not([disabled])) .switch:active {
      background: ${iD};
      border-color: ${iD};
    }

    :host([aria-invalid='true'][aria-checked='true']:not([disabled]))
      .switch:active {
      background-color: ${oF};
      border-color: ${oF};
    }

    :host([aria-checked='true']:not([disabled]))
      .switch:active
      .checked-indicator {
      background: ${iP};
    }

    :host([aria-checked="true"]:${rm}:not([disabled])) .switch {
      outline: solid calc(${tk} * 1px) ${iE};
    }

    :host([aria-invalid='true'][aria-checked="true"]:${rm}:not([disabled])) .switch {
      outline-color: ${oS};
    }

    .unchecked-message {
      display: block;
    }

    .checked-message {
      display: none;
    }

    :host([aria-checked='true']) .unchecked-message {
      display: none;
    }

    :host([aria-checked='true']) .checked-message {
      display: block;
    }
  `.withBehaviors(rv(rh`
      .checked-indicator,
      :host(:not([disabled])) .switch:active .checked-indicator {
        forced-color-adjust: none;
        background: ${nX.FieldText};
      }
      .switch {
        forced-color-adjust: none;
        background: ${nX.Field};
        border-color: ${nX.FieldText};
      }
      :host([aria-invalid='true']) .switch {
        border-style: dashed;
      }
      :host(:not([disabled])) .switch:hover {
        background: ${nX.HighlightText};
        border-color: ${nX.Highlight};
      }
      :host([aria-checked='true']) .switch {
        background: ${nX.Highlight};
        border-color: ${nX.Highlight};
      }
      :host([aria-checked='true']:not([disabled])) .switch:hover,
      :host(:not([disabled])) .switch:active {
        background: ${nX.HighlightText};
        border-color: ${nX.Highlight};
      }
      :host([aria-checked='true']) .checked-indicator {
        background: ${nX.HighlightText};
      }
      :host([aria-checked='true']:not([disabled]))
        .switch:hover
        .checked-indicator {
        background: ${nX.Highlight};
      }
      :host([disabled]) {
        opacity: 1;
      }
      :host(:${rm}) .switch {
        border-color: ${nX.Highlight};
        outline-offset: 2px;
        outline: solid calc(${tk} * 1px) ${nX.FieldText};
      }
      :host([aria-checked="true"]:${rm}:not([disabled])) .switch {
        outline: solid calc(${tk} * 1px) ${nX.FieldText};
      }
      :host([disabled]) .checked-indicator {
        background: ${nX.GrayText};
      }
      :host([disabled]) .switch {
        background: ${nX.Field};
        border-color: ${nX.GrayText};
      }
    `),new rZ(rh`
        .checked-indicator {
          left: 5px;
          right: calc(((${ry} / 2) + 1) * 1px);
        }

        :host([aria-checked='true']) .checked-indicator {
          left: calc(((${ry} / 2) + 1) * 1px);
          right: 5px;
        }
      `,rh`
        .checked-indicator {
          right: 5px;
          left: calc(((${ry} / 2) + 1) * 1px);
        }

        :host([aria-checked='true']) .checked-indicator {
          right: calc(((${ry} / 2) + 1) * 1px);
          left: 5px;
        }
      `));class hA extends hR{}let hV=hA.compose({baseName:"switch",baseClass:hR,template:(e,t)=>s4`
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
            <slot ${rn("defaultSlottedNodes")}></slot>
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
`,styles:hL,switch:`
        <span class="checked-indicator" part="checked-indicator"></span>
    `});class hP extends sR{}let hH=(e,t)=>rh`
  ${ru("block")} :host {
    box-sizing: border-box;
    font-size: ${tC};
    line-height: ${tI};
    padding: 0 calc((6 + (${tv} * 2 * ${tf})) * 1px);
  }
`;class hz extends hP{}let hM=hz.compose({baseName:"tab-panel",baseClass:hP,template:(e,t)=>s4`
    <template slot="tabpanel" role="tabpanel">
        <slot></slot>
    </template>
`,styles:hH});class hN extends sR{}X([eS({mode:"boolean"})],hN.prototype,"disabled",void 0);let hB=(e,t)=>rh`
    ${ru("inline-flex")} :host {
      box-sizing: border-box;
      font-family: ${tu};
      font-size: ${tC};
      line-height: ${tI};
      height: calc(${ry} * 1px);
      padding: calc(${tv} * 5px) calc(${tv} * 4px);
      color: ${oc};
      fill: currentcolor;
      border-radius: 0 0 calc(${tb} * 1px)
        calc(${tb} * 1px);
      border: calc(${tw} * 1px) solid transparent;
      align-items: center;
      justify-content: center;
      grid-row: 2;
      cursor: pointer;
    }

    :host(:hover) {
      color: ${op};
      fill: currentcolor;
    }

    :host(:active) {
      color: ${op};
      fill: currentcolor;
    }

    :host([disabled]) {
      cursor: ${am};
      opacity: ${t$};
    }

    :host([disabled]:hover) {
      color: ${oc};
      background: ${i4};
    }

    :host([aria-selected='true']) {
      background: ${iX};
      color: ${op};
      fill: currentcolor;
    }

    :host([aria-selected='true']:hover) {
      background: ${iY};
      color: ${op};
      fill: currentcolor;
    }

    :host([aria-selected='true']:active) {
      background: ${iQ};
      color: ${op};
      fill: currentcolor;
    }

    :host(:${rm}) {
      outline: none;
      border-color: ${iE};
      box-shadow: 0 0 0 calc((${tk} - ${tw}) * 1px)
        ${iE};
    }

    :host(:focus) {
      outline: none;
    }

    :host(.vertical) {
      justify-content: end;
      grid-column: 2;
      border-bottom-left-radius: 0;
      border-top-right-radius: calc(${tb} * 1px);
    }

    :host(.vertical[aria-selected='true']) {
      z-index: 2;
    }

    :host(.vertical:hover) {
      color: ${op};
    }

    :host(.vertical:active) {
      color: ${op};
    }

    :host(.vertical:hover[aria-selected='true']) {
    }
  `.withBehaviors(rv(rh`
      :host {
        forced-color-adjust: none;
        border-color: transparent;
        color: ${nX.ButtonText};
        fill: currentcolor;
      }
      :host(:hover),
      :host(.vertical:hover),
      :host([aria-selected='true']:hover) {
        background: ${nX.Highlight};
        color: ${nX.HighlightText};
        fill: currentcolor;
      }
      :host([aria-selected='true']) {
        background: ${nX.HighlightText};
        color: ${nX.Highlight};
        fill: currentcolor;
      }
      :host(:${rm}) {
        border-color: ${nX.ButtonText};
        box-shadow: none;
      }
      :host([disabled]),
      :host([disabled]:hover) {
        opacity: 1;
        color: ${nX.GrayText};
        background: ${nX.ButtonFace};
      }
    `));class hj extends hN{}let hq=hj.compose({baseName:"tab",baseClass:hN,template:(e,t)=>s4`
    <template slot="tab" role="tab" aria-disabled="${e=>e.disabled}">
        <slot></slot>
    </template>
`,styles:hB}),hU="horizontal";class h_ extends sR{constructor(){super(...arguments),this.orientation=hU,this.activeindicator=!0,this.showActiveIndicator=!0,this.prevActiveTabIndex=0,this.activeTabIndex=0,this.ticking=!1,this.change=()=>{this.$emit("change",this.activetab)},this.isDisabledElement=e=>"true"===e.getAttribute("aria-disabled"),this.isHiddenElement=e=>e.hasAttribute("hidden"),this.isFocusableElement=e=>!this.isDisabledElement(e)&&!this.isHiddenElement(e),this.setTabs=()=>{let e="gridColumn",t="gridRow",i=this.isHorizontal()?e:t;this.activeTabIndex=this.getActiveIndex(),this.showActiveIndicator=!1,this.tabs.forEach((o,s)=>{if("tab"===o.slot){let e=this.activeTabIndex===s&&this.isFocusableElement(o);this.activeindicator&&this.isFocusableElement(o)&&(this.showActiveIndicator=!0);let t=this.tabIds[s],i=this.tabpanelIds[s];o.setAttribute("id",t),o.setAttribute("aria-selected",e?"true":"false"),o.setAttribute("aria-controls",i),o.addEventListener("click",this.handleTabClick),o.addEventListener("keydown",this.handleTabKeyDown),o.setAttribute("tabindex",e?"0":"-1"),e&&(this.activetab=o,this.activeid=t)}o.style[e]="",o.style[t]="",o.style[i]=`${s+1}`,this.isHorizontal()?o.classList.remove("vertical"):o.classList.add("vertical")})},this.setTabPanels=()=>{this.tabpanels.forEach((e,t)=>{let i=this.tabIds[t],o=this.tabpanelIds[t];e.setAttribute("id",o),e.setAttribute("aria-labelledby",i),this.activeTabIndex!==t?e.setAttribute("hidden",""):e.removeAttribute("hidden")})},this.handleTabClick=e=>{let t=e.currentTarget;1===t.nodeType&&this.isFocusableElement(t)&&(this.prevActiveTabIndex=this.activeTabIndex,this.activeTabIndex=this.tabs.indexOf(t),this.setComponent())},this.handleTabKeyDown=e=>{if(this.isHorizontal())switch(e.key){case oZ:e.preventDefault(),this.adjustBackward(e);break;case oJ:e.preventDefault(),this.adjustForward(e)}else switch(e.key){case o0:e.preventDefault(),this.adjustBackward(e);break;case oQ:e.preventDefault(),this.adjustForward(e)}switch(e.key){case"Home":e.preventDefault(),this.adjust(-this.activeTabIndex);break;case"End":e.preventDefault(),this.adjust(this.tabs.length-this.activeTabIndex-1)}},this.adjustForward=e=>{let t=this.tabs,i=0;for((i=this.activetab?t.indexOf(this.activetab)+1:1)===t.length&&(i=0);i<t.length&&t.length>1;)if(this.isFocusableElement(t[i])){this.moveToTabByIndex(t,i);break}else if(this.activetab&&i===t.indexOf(this.activetab))break;else i+1>=t.length?i=0:i+=1},this.adjustBackward=e=>{let t=this.tabs,i=0;for(i=(i=this.activetab?t.indexOf(this.activetab)-1:0)<0?t.length-1:i;i>=0&&t.length>1;)if(this.isFocusableElement(t[i])){this.moveToTabByIndex(t,i);break}else i-1<0?i=t.length-1:i-=1},this.moveToTabByIndex=(e,t)=>{let i=e[t];this.activetab=i,this.prevActiveTabIndex=this.activeTabIndex,this.activeTabIndex=t,i.focus(),this.setComponent()}}orientationChanged(){this.$fastController.isConnected&&(this.setTabs(),this.setTabPanels(),this.handleActiveIndicatorPosition())}activeidChanged(e,t){this.$fastController.isConnected&&this.tabs.length<=this.tabpanels.length&&(this.prevActiveTabIndex=this.tabs.findIndex(t=>t.id===e),this.setTabs(),this.setTabPanels(),this.handleActiveIndicatorPosition())}tabsChanged(){this.$fastController.isConnected&&this.tabs.length<=this.tabpanels.length&&(this.tabIds=this.getTabIds(),this.tabpanelIds=this.getTabPanelIds(),this.setTabs(),this.setTabPanels(),this.handleActiveIndicatorPosition())}tabpanelsChanged(){this.$fastController.isConnected&&this.tabpanels.length<=this.tabs.length&&(this.tabIds=this.getTabIds(),this.tabpanelIds=this.getTabPanelIds(),this.setTabs(),this.setTabPanels(),this.handleActiveIndicatorPosition())}getActiveIndex(){return void 0!==this.activeid?-1===this.tabIds.indexOf(this.activeid)?0:this.tabIds.indexOf(this.activeid):0}getTabIds(){return this.tabs.map(e=>{var t;return null!=(t=e.getAttribute("id"))?t:`tab-${aR()}`})}getTabPanelIds(){return this.tabpanels.map(e=>{var t;return null!=(t=e.getAttribute("id"))?t:`panel-${aR()}`})}setComponent(){this.activeTabIndex!==this.prevActiveTabIndex&&(this.activeid=this.tabIds[this.activeTabIndex],this.focusTab(),this.change())}isHorizontal(){return this.orientation===hU}handleActiveIndicatorPosition(){this.showActiveIndicator&&this.activeindicator&&this.activeTabIndex!==this.prevActiveTabIndex&&(this.ticking?this.ticking=!1:(this.ticking=!0,this.animateActiveIndicator()))}animateActiveIndicator(){this.ticking=!0;let e=this.isHorizontal()?"gridColumn":"gridRow",t=this.isHorizontal()?"translateX":"translateY",i=this.isHorizontal()?"offsetLeft":"offsetTop",o=this.activeIndicatorRef[i];this.activeIndicatorRef.style[e]=`${this.activeTabIndex+1}`;let s=this.activeIndicatorRef[i];this.activeIndicatorRef.style[e]=`${this.prevActiveTabIndex+1}`,this.activeIndicatorRef.style.transform=`${t}(${s-o}px)`,this.activeIndicatorRef.classList.add("activeIndicatorTransition"),this.activeIndicatorRef.addEventListener("transitionend",()=>{this.ticking=!1,this.activeIndicatorRef.style[e]=`${this.activeTabIndex+1}`,this.activeIndicatorRef.style.transform=`${t}(0px)`,this.activeIndicatorRef.classList.remove("activeIndicatorTransition")})}adjust(e){let t=this.tabs.filter(e=>this.isFocusableElement(e)),i=t.indexOf(this.activetab),o=o2(0,t.length-1,i+e),s=this.tabs.indexOf(t[o]);s>-1&&this.moveToTabByIndex(this.tabs,s)}focusTab(){this.tabs[this.activeTabIndex].focus()}connectedCallback(){super.connectedCallback(),this.tabIds=this.getTabIds(),this.tabpanelIds=this.getTabPanelIds(),this.activeTabIndex=this.getActiveIndex()}}X([eS],h_.prototype,"orientation",void 0),X([eS],h_.prototype,"activeid",void 0),X([eu],h_.prototype,"tabs",void 0),X([eu],h_.prototype,"tabpanels",void 0),X([eS({mode:"boolean"})],h_.prototype,"activeindicator",void 0),X([eu],h_.prototype,"activeIndicatorRef",void 0),X([eu],h_.prototype,"showActiveIndicator",void 0),rt(h_,s8);let hG=(e,t)=>rh`
    ${ru("grid")} :host {
      box-sizing: border-box;
      font-family: ${tu};
      font-size: ${tC};
      line-height: ${tI};
      color: ${op};
      grid-template-columns: auto 1fr auto;
      grid-template-rows: auto 1fr;
    }

    .tablist {
      display: grid;
      grid-template-rows: auto auto;
      grid-template-columns: auto;
      position: relative;
      width: max-content;
      align-self: end;
      padding: calc(${tv} * 4px) calc(${tv} * 4px) 0;
      box-sizing: border-box;
    }

    .start,
    .end {
      align-self: center;
    }

    .activeIndicator {
      grid-row: 1;
      grid-column: 1;
      width: 100%;
      height: 4px;
      justify-self: center;
      background: ${iS};
      margin-top: 0;
      border-radius: calc(${tb} * 1px)
        calc(${tb} * 1px) 0 0;
    }

    .activeIndicatorTransition {
      transition: transform 0.01s ease-in-out;
    }

    .tabpanel {
      grid-row: 2;
      grid-column-start: 1;
      grid-column-end: 4;
      position: relative;
    }

    :host([orientation='vertical']) {
      grid-template-rows: auto 1fr auto;
      grid-template-columns: auto 1fr;
    }

    :host([orientation='vertical']) .tablist {
      grid-row-start: 2;
      grid-row-end: 2;
      display: grid;
      grid-template-rows: auto;
      grid-template-columns: auto 1fr;
      position: relative;
      width: max-content;
      justify-self: end;
      align-self: flex-start;
      width: 100%;
      padding: 0 calc(${tv} * 4px)
        calc((${ry} - ${tv}) * 1px) 0;
    }

    :host([orientation='vertical']) .tabpanel {
      grid-column: 2;
      grid-row-start: 1;
      grid-row-end: 4;
    }

    :host([orientation='vertical']) .end {
      grid-row: 3;
    }

    :host([orientation='vertical']) .activeIndicator {
      grid-column: 1;
      grid-row: 1;
      width: 4px;
      height: 100%;
      margin-inline-end: 0px;
      align-self: center;
      background: ${iS};
      border-radius: calc(${tb} * 1px) 0 0
        calc(${tb} * 1px);
    }

    :host([orientation='vertical']) .activeIndicatorTransition {
      transition: transform 0.01s ease-in-out;
    }
  `.withBehaviors(rv(rh`
      .activeIndicator,
      :host([orientation='vertical']) .activeIndicator {
        forced-color-adjust: none;
        background: ${nX.Highlight};
      }
    `));class hW extends h_{}let hK=hW.compose({baseName:"tabs",baseClass:h_,template:(e,t)=>s4`
    <template class="${e=>e.orientation}">
        ${re(e,t)}
        <div class="tablist" part="tablist" role="tablist">
            <slot class="tab" name="tab" part="tab" ${rn("tabs")}></slot>

            ${rG(e=>e.showActiveIndicator,s4`
                    <div
                        ${s9("activeIndicatorRef")}
                        class="activeIndicator"
                        part="activeIndicator"
                    ></div>
                `)}
        </div>
        ${s7(e,t)}
        <div class="tabpanel" part="tabpanel">
            <slot name="tabpanel" ${rn("tabpanels")}></slot>
        </div>
    </template>
`,styles:hG});class hX extends sR{}class hY extends ah(hX){constructor(){super(...arguments),this.proxy=document.createElement("textarea")}}let hQ="none";class hZ extends hY{constructor(){super(...arguments),this.resize=hQ,this.cols=20,this.handleTextInput=()=>{this.value=this.control.value}}readOnlyChanged(){this.proxy instanceof HTMLTextAreaElement&&(this.proxy.readOnly=this.readOnly)}autofocusChanged(){this.proxy instanceof HTMLTextAreaElement&&(this.proxy.autofocus=this.autofocus)}listChanged(){this.proxy instanceof HTMLTextAreaElement&&this.proxy.setAttribute("list",this.list)}maxlengthChanged(){this.proxy instanceof HTMLTextAreaElement&&(this.proxy.maxLength=this.maxlength)}minlengthChanged(){this.proxy instanceof HTMLTextAreaElement&&(this.proxy.minLength=this.minlength)}spellcheckChanged(){this.proxy instanceof HTMLTextAreaElement&&(this.proxy.spellcheck=this.spellcheck)}select(){this.control.select(),this.$emit("select")}handleChange(){this.$emit("change")}validate(){super.validate(this.control)}}X([eS({mode:"boolean"})],hZ.prototype,"readOnly",void 0),X([eS],hZ.prototype,"resize",void 0),X([eS({mode:"boolean"})],hZ.prototype,"autofocus",void 0),X([eS({attribute:"form"})],hZ.prototype,"formId",void 0),X([eS],hZ.prototype,"list",void 0),X([eS({converter:eT})],hZ.prototype,"maxlength",void 0),X([eS({converter:eT})],hZ.prototype,"minlength",void 0),X([eS],hZ.prototype,"name",void 0),X([eS],hZ.prototype,"placeholder",void 0),X([eS({converter:eT,mode:"fromView"})],hZ.prototype,"cols",void 0),X([eS({converter:eT,mode:"fromView"})],hZ.prototype,"rows",void 0),X([eS({mode:"boolean"})],hZ.prototype,"spellcheck",void 0),X([eu],hZ.prototype,"defaultSlottedNodes",void 0),rt(hZ,nT);let hJ=(e,t)=>rh`
    ${ru("inline-block")} :host {
      font-family: ${tu};
      outline: none;
      user-select: none;
    }

    .control {
      box-sizing: border-box;
      position: relative;
      color: ${op};
      background: ${i0};
      border-radius: calc(${tb} * 1px);
      border: calc(${tw} * 1px) solid ${oe};
      height: calc(${ry} * 2px);
      font: inherit;
      font-size: ${tC};
      line-height: ${tI};
      padding: calc(${tv} * 2px + 1px);
      width: 100%;
      resize: none;
    }

    :host([aria-invalid='true']) .control {
      border-color: ${oI};
    }

    .control:hover:enabled {
      background: ${i1};
      border-color: ${ot};
    }

    :host([aria-invalid='true']) .control:hover:enabled {
      border-color: ${oT};
    }

    .control:active:enabled {
      background: ${i2};
      border-color: ${oi};
    }

    :host([aria-invalid='true']) .control:active:enabled {
      border-color: ${oF};
    }

    .control:hover,
    .control:${rm},
    .control:disabled,
    .control:active {
      outline: none;
    }

    :host(:focus-within) .control {
      border-color: ${iE};
      box-shadow: 0 0 0 calc((${tk} - ${tw}) * 1px)
        ${iE};
    }

    :host([aria-invalid='true']:focus-within) .control {
      border-color: ${oS};
      box-shadow: 0 0 0 calc((${tk} - ${tw}) * 1px)
        ${oS};
    }

    :host([appearance='filled']) .control {
      background: ${iX};
    }

    :host([appearance='filled']:hover:not([disabled])) .control {
      background: ${iY};
    }

    :host([resize='both']) .control {
      resize: both;
    }

    :host([resize='horizontal']) .control {
      resize: horizontal;
    }

    :host([resize='vertical']) .control {
      resize: vertical;
    }

    .label {
      display: block;
      color: ${op};
      cursor: pointer;
      font-size: ${tC};
      line-height: ${tI};
      margin-bottom: 4px;
    }

    .label__hidden {
      display: none;
      visibility: hidden;
    }

    :host([disabled]) .label,
    :host([readonly]) .label,
    :host([readonly]) .control,
    :host([disabled]) .control {
      cursor: ${am};
    }
    :host([disabled]) {
      opacity: ${t$};
    }
    :host([disabled]) .control {
      border-color: ${om};
    }

    :host([cols]) {
      width: initial;
    }

    :host([rows]) .control {
      height: initial;
    }
  `.withBehaviors(rv(rh`
      :host([disabled]) {
        opacity: 1;
      }

      :host([aria-invalid='true']) .control {
        border-style: dashed;
      }
    `));class h0 extends hZ{constructor(){super(...arguments),this.appearance="outline"}}rI([eS],h0.prototype,"appearance",void 0);let h1=h0.compose({baseName:"text-area",baseClass:hZ,template:(e,t)=>s4`
    <template
        class="
            ${e=>e.readOnly?"readonly":""}
            ${e=>e.resize!==hQ?`resize-${e.resize}`:""}"
    >
        <label
            part="label"
            for="control"
            class="${e=>e.defaultSlottedNodes&&e.defaultSlottedNodes.length?"label":"label label__hidden"}"
        >
            <slot ${rn("defaultSlottedNodes")}></slot>
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
            ${s9("control")}
        ></textarea>
    </template>
`,styles:hJ,shadowOptions:{delegatesFocus:!0}}),h2=(e,t)=>rh`
  ${nR}

  .start,
    .end {
    display: flex;
  }
`;class h5 extends nI{constructor(){super(...arguments),this.appearance="outline"}}rI([eS],h5.prototype,"appearance",void 0);let h3=h5.compose({baseName:"text-field",baseClass:nI,template:(e,t)=>s4`
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
                ${rn({property:"defaultSlottedNodes",filter:nA})}
            ></slot>
        </label>
        <div class="root" part="root">
            ${re(e,t)}
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
                ${s9("control")}
            />
            ${s7(e,t)}
        </div>
    </template>
`,styles:h2,shadowOptions:{delegatesFocus:!0}}),h4=(e,t)=>rh`
    ${ru("inline-flex")} :host {
      --toolbar-item-gap: calc(
        (var(--design-unit) + calc(var(--density) + 2)) * 1px
      );
      background-color: ${iT};
      border-radius: calc(${tb} * 1px);
      fill: currentcolor;
      padding: var(--toolbar-item-gap);
    }

    :host(${rm}) {
      outline: calc(${tw} * 1px) solid ${iE};
    }

    .positioning-region {
      align-items: flex-start;
      display: inline-flex;
      flex-flow: row wrap;
      justify-content: flex-start;
      width: 100%;
      height: 100%;
    }

    :host([orientation='vertical']) .positioning-region {
      flex-direction: column;
    }

    ::slotted(:not([slot])) {
      flex: 0 0 auto;
      margin: 0 var(--toolbar-item-gap);
    }

    :host([orientation='vertical']) ::slotted(:not([slot])) {
      margin: var(--toolbar-item-gap) 0;
    }

    .start,
    .end {
      display: flex;
      margin: auto;
      margin-inline: 0;
    }

    ::slotted(svg) {
      /* TODO: adaptive typography https://github.com/microsoft/fast/issues/2432 */
      width: 16px;
      height: 16px;
    }
  `.withBehaviors(rv(rh`
      :host(:${rm}) {
        box-shadow: 0 0 0 calc(${tk} * 1px)
          ${nX.Highlight};
        color: ${nX.ButtonText};
        forced-color-adjust: none;
      }
    `)),h6=Object.freeze({[o1.ArrowUp]:{[lv]:-1},[o1.ArrowDown]:{[lv]:1},[o1.ArrowLeft]:{[lf]:{[nU.ltr]:-1,[nU.rtl]:1}},[o1.ArrowRight]:{[lf]:{[nU.ltr]:1,[nU.rtl]:-1}}});class h9 extends sR{constructor(){super(...arguments),this._activeIndex=0,this.direction=nU.ltr,this.orientation=lf}get activeIndex(){return ec.track(this,"activeIndex"),this._activeIndex}set activeIndex(e){this.$fastController.isConnected&&(this._activeIndex=o2(0,this.focusableElements.length-1,e),ec.notify(this,"activeIndex"))}slottedItemsChanged(){this.$fastController.isConnected&&this.reduceFocusableElements()}mouseDownHandler(e){var t;let i=null==(t=this.focusableElements)?void 0:t.findIndex(t=>t.contains(e.target));return i>-1&&this.activeIndex!==i&&this.setFocusedElement(i),!0}childItemsChanged(e,t){this.$fastController.isConnected&&this.reduceFocusableElements()}connectedCallback(){super.connectedCallback(),this.direction=rj(this)}focusinHandler(e){let t=e.relatedTarget;!t||this.contains(t)||this.setFocusedElement()}getDirectionalIncrementer(e){var t,i,o,s,r;return null!=(r=null!=(o=null==(i=null==(t=h6[e])?void 0:t[this.orientation])?void 0:i[this.direction])?o:null==(s=h6[e])?void 0:s[this.orientation])?r:0}keydownHandler(e){let t=e.key;if(!(t in o1)||e.defaultPrevented||e.shiftKey)return!0;let i=this.getDirectionalIncrementer(t);if(!i)return!e.target.closest("[role=radiogroup]");let o=this.activeIndex+i;return this.focusableElements[o]&&e.preventDefault(),this.setFocusedElement(o),!0}get allSlottedItems(){return[...this.start.assignedElements(),...this.slottedItems,...this.end.assignedElements()]}reduceFocusableElements(){var e;let t=null==(e=this.focusableElements)?void 0:e[this.activeIndex];this.focusableElements=this.allSlottedItems.reduce(h9.reduceFocusableItems,[]);let i=this.focusableElements.indexOf(t);this.activeIndex=Math.max(0,i),this.setFocusableElements()}setFocusedElement(e=this.activeIndex){this.activeIndex=e,this.setFocusableElements(),this.focusableElements[this.activeIndex]&&this.contains(document.activeElement)&&this.focusableElements[this.activeIndex].focus()}static reduceFocusableItems(e,t){var i,o,s,r;let a="radio"===t.getAttribute("role"),n=null==(o=null==(i=t.$fastController)?void 0:i.definition.shadowOptions)?void 0:o.delegatesFocus,l=Array.from(null!=(r=null==(s=t.shadowRoot)?void 0:s.querySelectorAll("*"))?r:[]).some(e=>ll(e));return!t.hasAttribute("disabled")&&!t.hasAttribute("hidden")&&(ll(t)||a||n||l)?(e.push(t),e):t.childElementCount?e.concat(Array.from(t.children).reduce(h9.reduceFocusableItems,[])):e}setFocusableElements(){this.$fastController.isConnected&&this.focusableElements.length>0&&this.focusableElements.forEach((e,t)=>{e.tabIndex=this.activeIndex===t?0:-1})}}rI([eu],h9.prototype,"direction",void 0),rI([eS],h9.prototype,"orientation",void 0),rI([eu],h9.prototype,"slottedItems",void 0),rI([eu],h9.prototype,"slottedLabel",void 0),rI([eu],h9.prototype,"childItems",void 0);class h8{}rI([eS({attribute:"aria-labelledby"})],h8.prototype,"ariaLabelledby",void 0),rI([eS({attribute:"aria-label"})],h8.prototype,"ariaLabel",void 0),rt(h8,rT),rt(h9,s8,h8);class h7 extends h9{connectedCallback(){super.connectedCallback();let e=eM(this);e&&iT.setValueFor(this,t=>os.getValueFor(t).evaluate(t,iT.getValueFor(e)))}}let de=h7.compose({baseName:"toolbar",baseClass:h9,template:(e,t)=>s4`
    <template
        aria-label="${e=>e.ariaLabel}"
        aria-labelledby="${e=>e.ariaLabelledby}"
        aria-orientation="${e=>e.orientation}"
        orientation="${e=>e.orientation}"
        role="toolbar"
        @mousedown="${(e,t)=>e.mouseDownHandler(t.event)}"
        @focusin="${(e,t)=>e.focusinHandler(t.event)}"
        @keydown="${(e,t)=>e.keydownHandler(t.event)}"
        ${nu({property:"childItems",attributeFilter:["disabled","hidden"],filter:rs(),subtree:!0})}
    >
        <slot name="label"></slot>
        <div class="positioning-region" part="positioning-region">
            ${re(e,t)}
            <slot
                ${rn({filter:rs(),property:"slottedItems"})}
            ></slot>
            ${s7(e,t)}
        </div>
    </template>
`,styles:h4,shadowOptions:{delegatesFocus:!0}});class dt extends sR{constructor(){super(...arguments),this.anchor="",this.delay=300,this.autoUpdateMode="anchor",this.anchorElement=null,this.viewportElement=null,this.verticalPositioningMode="dynamic",this.horizontalPositioningMode="dynamic",this.horizontalInset="false",this.verticalInset="false",this.horizontalScaling="content",this.verticalScaling="content",this.verticalDefaultPosition=void 0,this.horizontalDefaultPosition=void 0,this.tooltipVisible=!1,this.currentDirection=nU.ltr,this.showDelayTimer=null,this.hideDelayTimer=null,this.isAnchorHoveredFocused=!1,this.isRegionHovered=!1,this.handlePositionChange=e=>{this.classList.toggle("top","start"===this.region.verticalPosition),this.classList.toggle("bottom","end"===this.region.verticalPosition),this.classList.toggle("inset-top","insetStart"===this.region.verticalPosition),this.classList.toggle("inset-bottom","insetEnd"===this.region.verticalPosition),this.classList.toggle("center-vertical","center"===this.region.verticalPosition),this.classList.toggle("left","start"===this.region.horizontalPosition),this.classList.toggle("right","end"===this.region.horizontalPosition),this.classList.toggle("inset-left","insetStart"===this.region.horizontalPosition),this.classList.toggle("inset-right","insetEnd"===this.region.horizontalPosition),this.classList.toggle("center-horizontal","center"===this.region.horizontalPosition)},this.handleRegionMouseOver=e=>{this.isRegionHovered=!0},this.handleRegionMouseOut=e=>{this.isRegionHovered=!1,this.startHideDelayTimer()},this.handleAnchorMouseOver=e=>{if(this.tooltipVisible){this.isAnchorHoveredFocused=!0;return}this.startShowDelayTimer()},this.handleAnchorMouseOut=e=>{this.isAnchorHoveredFocused=!1,this.clearShowDelayTimer(),this.startHideDelayTimer()},this.handleAnchorFocusIn=e=>{this.startShowDelayTimer()},this.handleAnchorFocusOut=e=>{this.isAnchorHoveredFocused=!1,this.clearShowDelayTimer(),this.startHideDelayTimer()},this.startHideDelayTimer=()=>{this.clearHideDelayTimer(),this.tooltipVisible&&(this.hideDelayTimer=window.setTimeout(()=>{this.updateTooltipVisibility()},60))},this.clearHideDelayTimer=()=>{null!==this.hideDelayTimer&&(clearTimeout(this.hideDelayTimer),this.hideDelayTimer=null)},this.startShowDelayTimer=()=>{if(!this.isAnchorHoveredFocused){if(this.delay>1){null===this.showDelayTimer&&(this.showDelayTimer=window.setTimeout(()=>{this.startHover()},this.delay));return}this.startHover()}},this.startHover=()=>{this.isAnchorHoveredFocused=!0,this.updateTooltipVisibility()},this.clearShowDelayTimer=()=>{null!==this.showDelayTimer&&(clearTimeout(this.showDelayTimer),this.showDelayTimer=null)},this.getAnchor=()=>{let e=this.getRootNode();return e instanceof ShadowRoot?e.getElementById(this.anchor):document.getElementById(this.anchor)},this.handleDocumentKeydown=e=>{!e.defaultPrevented&&this.tooltipVisible&&"Escape"===e.key&&(this.isAnchorHoveredFocused=!1,this.updateTooltipVisibility(),this.$emit("dismiss"))},this.updateTooltipVisibility=()=>{if(!1===this.visible)this.hideTooltip();else{if(!0===this.visible||this.isAnchorHoveredFocused||this.isRegionHovered)return void this.showTooltip();this.hideTooltip()}},this.showTooltip=()=>{this.tooltipVisible||(this.currentDirection=rj(this),this.tooltipVisible=!0,document.addEventListener("keydown",this.handleDocumentKeydown),el.queueUpdate(this.setRegionProps))},this.hideTooltip=()=>{this.tooltipVisible&&(this.clearHideDelayTimer(),null!==this.region&&void 0!==this.region&&(this.region.removeEventListener("positionchange",this.handlePositionChange),this.region.viewportElement=null,this.region.anchorElement=null,this.region.removeEventListener("mouseover",this.handleRegionMouseOver),this.region.removeEventListener("mouseout",this.handleRegionMouseOut)),document.removeEventListener("keydown",this.handleDocumentKeydown),this.tooltipVisible=!1)},this.setRegionProps=()=>{this.tooltipVisible&&(this.region.viewportElement=this.viewportElement,this.region.anchorElement=this.anchorElement,this.region.addEventListener("positionchange",this.handlePositionChange),this.region.addEventListener("mouseover",this.handleRegionMouseOver,{passive:!0}),this.region.addEventListener("mouseout",this.handleRegionMouseOut,{passive:!0}))}}visibleChanged(){this.$fastController.isConnected&&(this.updateTooltipVisibility(),this.updateLayout())}anchorChanged(){this.$fastController.isConnected&&(this.anchorElement=this.getAnchor())}positionChanged(){this.$fastController.isConnected&&this.updateLayout()}anchorElementChanged(e){if(this.$fastController.isConnected){if(null!=e&&(e.removeEventListener("mouseover",this.handleAnchorMouseOver),e.removeEventListener("mouseout",this.handleAnchorMouseOut),e.removeEventListener("focusin",this.handleAnchorFocusIn),e.removeEventListener("focusout",this.handleAnchorFocusOut)),null!==this.anchorElement&&void 0!==this.anchorElement){this.anchorElement.addEventListener("mouseover",this.handleAnchorMouseOver,{passive:!0}),this.anchorElement.addEventListener("mouseout",this.handleAnchorMouseOut,{passive:!0}),this.anchorElement.addEventListener("focusin",this.handleAnchorFocusIn,{passive:!0}),this.anchorElement.addEventListener("focusout",this.handleAnchorFocusOut,{passive:!0});let e=this.anchorElement.id;null!==this.anchorElement.parentElement&&this.anchorElement.parentElement.querySelectorAll(":hover").forEach(t=>{t.id===e&&this.startShowDelayTimer()})}null!==this.region&&void 0!==this.region&&this.tooltipVisible&&(this.region.anchorElement=this.anchorElement),this.updateLayout()}}viewportElementChanged(){null!==this.region&&void 0!==this.region&&(this.region.viewportElement=this.viewportElement),this.updateLayout()}connectedCallback(){super.connectedCallback(),this.anchorElement=this.getAnchor(),this.updateTooltipVisibility()}disconnectedCallback(){this.hideTooltip(),this.clearShowDelayTimer(),this.clearHideDelayTimer(),super.disconnectedCallback()}updateLayout(){switch(this.verticalPositioningMode="locktodefault",this.horizontalPositioningMode="locktodefault",this.position){case"top":case"bottom":this.verticalDefaultPosition=this.position,this.horizontalDefaultPosition="center";break;case"right":case"left":case"start":case"end":this.verticalDefaultPosition="center",this.horizontalDefaultPosition=this.position;break;case"top-left":this.verticalDefaultPosition="top",this.horizontalDefaultPosition="left";break;case"top-right":this.verticalDefaultPosition="top",this.horizontalDefaultPosition="right";break;case"bottom-left":this.verticalDefaultPosition="bottom",this.horizontalDefaultPosition="left";break;case"bottom-right":this.verticalDefaultPosition="bottom",this.horizontalDefaultPosition="right";break;case"top-start":this.verticalDefaultPosition="top",this.horizontalDefaultPosition="start";break;case"top-end":this.verticalDefaultPosition="top",this.horizontalDefaultPosition="end";break;case"bottom-start":this.verticalDefaultPosition="bottom",this.horizontalDefaultPosition="start";break;case"bottom-end":this.verticalDefaultPosition="bottom",this.horizontalDefaultPosition="end";break;default:this.verticalPositioningMode="dynamic",this.horizontalPositioningMode="dynamic",this.verticalDefaultPosition=void 0,this.horizontalDefaultPosition="center"}}}X([eS({mode:"boolean"})],dt.prototype,"visible",void 0),X([eS],dt.prototype,"anchor",void 0),X([eS],dt.prototype,"delay",void 0),X([eS],dt.prototype,"position",void 0),X([eS({attribute:"auto-update-mode"})],dt.prototype,"autoUpdateMode",void 0),X([eS({attribute:"horizontal-viewport-lock"})],dt.prototype,"horizontalViewportLock",void 0),X([eS({attribute:"vertical-viewport-lock"})],dt.prototype,"verticalViewportLock",void 0),X([eu],dt.prototype,"anchorElement",void 0),X([eu],dt.prototype,"viewportElement",void 0),X([eu],dt.prototype,"verticalPositioningMode",void 0),X([eu],dt.prototype,"horizontalPositioningMode",void 0),X([eu],dt.prototype,"horizontalInset",void 0),X([eu],dt.prototype,"verticalInset",void 0),X([eu],dt.prototype,"horizontalScaling",void 0),X([eu],dt.prototype,"verticalScaling",void 0),X([eu],dt.prototype,"verticalDefaultPosition",void 0),X([eu],dt.prototype,"horizontalDefaultPosition",void 0),X([eu],dt.prototype,"tooltipVisible",void 0),X([eu],dt.prototype,"currentDirection",void 0);let di=(e,t)=>{let i=e.tagFor(rq);return rh`
    :host {
      contain: size;
      overflow: visible;
      height: 0;
      width: 0;
    }

    .tooltip {
      box-sizing: border-box;
      border-radius: calc(${tb} * 1px);
      border: calc(${tw} * 1px) solid ${on};
      box-shadow: 0 0 0 1px ${on} inset;
      background: ${iX};
      color: ${op};
      padding: 4px;
      height: fit-content;
      width: fit-content;
      font-family: ${tu};
      font-size: ${tC};
      line-height: ${tI};
      white-space: nowrap;
      /* TODO: a mechanism to manage z-index across components
                    https://github.com/microsoft/fast/issues/3813 */
      z-index: 10000;
    }

    ${i} {
      display: flex;
      justify-content: center;
      align-items: center;
      overflow: visible;
      flex-direction: row;
    }

    ${i}.right,
    ${i}.left {
      flex-direction: column;
    }

    ${i}.top .tooltip {
      margin-bottom: 4px;
    }

    ${i}.bottom .tooltip {
      margin-top: 4px;
    }

    ${i}.left .tooltip {
      margin-right: 4px;
    }

    ${i}.right .tooltip {
      margin-left: 4px;
    }

    ${i}.top.left .tooltip,
            ${i}.top.right .tooltip {
      margin-bottom: 0px;
    }

    ${i}.bottom.left .tooltip,
            ${i}.bottom.right .tooltip {
      margin-top: 0px;
    }

    ${i}.top.left .tooltip,
            ${i}.bottom.left .tooltip {
      margin-right: 0px;
    }

    ${i}.top.right .tooltip,
            ${i}.bottom.right .tooltip {
      margin-left: 0px;
    }
  `.withBehaviors(rv(rh`
      :host([disabled]) {
        opacity: 1;
      }
    `))};class ds extends dt{}let dr=ds.compose({baseName:"tooltip",baseClass:dt,template:(e,t)=>s4`
        ${rG(e=>e.tooltipVisible,s4`
            <${e.tagFor(rq)}
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
                ${s9("region")}
            >
                <div class="tooltip" part="tooltip" role="tooltip">
                    <slot></slot>
                </div>
            </${e.tagFor(rq)}>
        `)}
    `,styles:di});function da(e){return rg(e)&&"treeitem"===e.getAttribute("role")}class dn extends sR{constructor(){super(...arguments),this.expanded=!1,this.focusable=!1,this.isNestedItem=()=>da(this.parentElement),this.handleExpandCollapseButtonClick=e=>{this.disabled||e.defaultPrevented||(this.expanded=!this.expanded)},this.handleFocus=e=>{this.setAttribute("tabindex","0")},this.handleBlur=e=>{this.setAttribute("tabindex","-1")}}expandedChanged(){this.$fastController.isConnected&&this.$emit("expanded-change",this)}selectedChanged(){this.$fastController.isConnected&&this.$emit("selected-change",this)}itemsChanged(e,t){this.$fastController.isConnected&&this.items.forEach(e=>{da(e)&&(e.nested=!0)})}static focusItem(e){e.focusable=!0,e.focus()}childItemLength(){let e=this.childItems.filter(e=>da(e));return e?e.length:0}}X([eS({mode:"boolean"})],dn.prototype,"expanded",void 0),X([eS({mode:"boolean"})],dn.prototype,"selected",void 0),X([eS({mode:"boolean"})],dn.prototype,"disabled",void 0),X([eu],dn.prototype,"focusable",void 0),X([eu],dn.prototype,"childItems",void 0),X([eu],dn.prototype,"items",void 0),X([eu],dn.prototype,"nested",void 0),X([eu],dn.prototype,"renderCollapsedChildren",void 0),rt(dn,s8);let dl=rc`(((${tp} + ${tf}) * 0.5 + 2) * ${tv})`,dh=rh`
  .expand-collapse-glyph {
    transform: rotate(0deg);
  }
  :host(.nested) .expand-collapse-button {
    left: var(
      --expand-collapse-button-nested-width,
      calc(
        (
            ${dl} +
              ((${tp} + ${tf}) * 1.25)
          ) * -1px
      )
    );
  }
  :host([selected])::after {
    left: calc(${tk} * 1px);
  }
  :host([expanded]) > .positioning-region .expand-collapse-glyph {
    transform: rotate(90deg);
  }
`,dd=rh`
  .expand-collapse-glyph {
    transform: rotate(180deg);
  }
  :host(.nested) .expand-collapse-button {
    right: var(
      --expand-collapse-button-nested-width,
      calc(
        (
            ${dl} +
              ((${tp} + ${tf}) * 1.25)
          ) * -1px
      )
    );
  }
  :host([selected])::after {
    right: calc(${tk} * 1px);
  }
  :host([expanded]) > .positioning-region .expand-collapse-glyph {
    transform: rotate(90deg);
  }
`,dc=e2.create("tree-item-expand-collapse-hover").withDefault(e=>{let t=i3.getValueFor(e);return t.evaluate(e,t.evaluate(e).hover).hover}),du=e2.create("tree-item-expand-collapse-selected-hover").withDefault(e=>{let t=iK.getValueFor(e);return i3.getValueFor(e).evaluate(e,t.evaluate(e).rest).hover}),dp=(e,t)=>rh`
    /**
     * This animation exists because when tree item children are conditionally loaded
     * there is a visual bug where the DOM exists but styles have not yet been applied (essentially FOUC).
     * This subtle animation provides a ever so slight timing adjustment for loading that solves the issue.
     */
    @keyframes treeItemLoading {
      0% {
        opacity: 0;
      }
      100% {
        opacity: 1;
      }
    }

    ${ru("block")} :host {
      contain: content;
      position: relative;
      outline: none;
      color: ${op};
      background: ${i4};
      cursor: pointer;
      font-family: ${tu};
      --tree-item-nested-width: 0;
    }

    :host(:focus) > .positioning-region {
      outline: none;
    }

    :host(:focus) .content-region {
      outline: none;
    }

    :host(:${rm}) .positioning-region {
      border-color: ${iE};
      box-shadow: 0 0 0 calc((${tk} - ${tw}) * 1px)
        ${iE} inset;
      color: ${op};
    }

    .positioning-region {
      display: flex;
      position: relative;
      box-sizing: border-box;
      background: ${i4};
      border: transparent calc(${tw} * 1px) solid;
      border-radius: calc(${tb} * 1px);
      height: calc((${ry} + 1) * 1px);
    }

    .positioning-region::before {
      content: '';
      display: block;
      width: var(--tree-item-nested-width);
      flex-shrink: 0;
    }

    :host(:not([disabled])) .positioning-region:hover {
      background: ${i6};
    }

    :host(:not([disabled])) .positioning-region:active {
      background: ${i9};
    }

    .content-region {
      display: inline-flex;
      align-items: center;
      white-space: nowrap;
      width: 100%;
      min-width: 0;
      height: calc(${ry} * 1px);
      margin-inline-start: calc(${tv} * 2px + 8px);
      font-size: ${tC};
      line-height: ${tI};
      font-weight: 400;
    }

    .items {
      /* TODO: adaptive typography https://github.com/microsoft/fast/issues/2432 */
      font-size: calc(1em + (${tv} + 16) * 1px);
    }

    .expand-collapse-button {
      background: none;
      border: none;
      outline: none;
      /* TODO: adaptive typography https://github.com/microsoft/fast/issues/2432 */
      width: calc(${dl} * 1px);
      height: calc(${dl} * 1px);
      padding: 0;
      display: flex;
      justify-content: center;
      align-items: center;
      cursor: pointer;
      margin-left: 6px;
      margin-right: 6px;
    }

    .expand-collapse-glyph {
      /* TODO: adaptive typography https://github.com/microsoft/fast/issues/2432 */
      width: calc((16 + ${tf}) * 1px);
      height: calc((16 + ${tf}) * 1px);
      transition: transform 0.1s linear;

      pointer-events: none;
      fill: currentcolor;
    }

    .start,
    .end {
      display: flex;
      fill: currentcolor;
    }

    ::slotted(svg) {
      /* TODO: adaptive typography https://github.com/microsoft/fast/issues/2432 */
      width: 16px;
      height: 16px;

      /* Something like that would do if the typography is adaptive
      font-size: inherit;
      width: ${tD};
      height: ${tD};
      */
    }

    .start {
      /* TODO: horizontalSpacing https://github.com/microsoft/fast/issues/2766 */
      margin-inline-end: calc(${tv} * 2px + 2px);
    }

    .end {
      /* TODO: horizontalSpacing https://github.com/microsoft/fast/issues/2766 */
      margin-inline-start: calc(${tv} * 2px + 2px);
    }

    :host([expanded]) > .items {
      animation: treeItemLoading ease-in 10ms;
      animation-iteration-count: 1;
      animation-fill-mode: forwards;
    }

    :host([disabled]) .content-region {
      opacity: ${t$};
      cursor: ${am};
    }

    :host(.nested) .content-region {
      position: relative;
      /* Add left margin to collapse button size */
      margin-inline-start: calc(
        (
            ${dl} +
              ((${tp} + ${tf}) * 1.25)
          ) * 1px
      );
    }

    :host(.nested) .expand-collapse-button {
      position: absolute;
    }

    :host(.nested:not([disabled])) .expand-collapse-button:hover {
      background: ${dc};
    }

    :host([selected]) .positioning-region {
      background: ${iX};
    }

    :host([selected]:not([disabled])) .positioning-region:hover {
      background: ${iY};
    }

    :host([selected]:not([disabled])) .positioning-region:active {
      background: ${iQ};
    }

    :host([selected]:not([disabled])) .expand-collapse-button:hover {
      background: ${du};
    }

    :host([selected])::after {
      /* The background needs to be calculated based on the selected background state
         for this control. We currently have no way of changing that, so setting to
         accent-foreground-rest for the time being */
      background: ${iU};
      border-radius: calc(${tb} * 1px);
      content: '';
      display: block;
      position: absolute;
      top: calc((${ry} / 4) * 1px);
      width: 3px;
      height: calc((${ry} / 2) * 1px);
    }

    ::slotted(${e.tagFor(dn)}) {
      --tree-item-nested-width: 1em;
      --expand-collapse-button-nested-width: calc(
        (
            ${dl} +
              ((${tp} + ${tf}) * 1.25)
          ) * -1px
      );
    }
  `.withBehaviors(new rZ(dh,dd),rv(rh`
      :host {
        forced-color-adjust: none;
        border-color: transparent;
        background: ${nX.Field};
        color: ${nX.FieldText};
      }
      :host .content-region .expand-collapse-glyph {
        fill: ${nX.FieldText};
      }
      :host .positioning-region:hover,
      :host([selected]) .positioning-region {
        background: ${nX.Highlight};
      }
      :host .positioning-region:hover .content-region,
      :host([selected]) .positioning-region .content-region {
        color: ${nX.HighlightText};
      }
      :host .positioning-region:hover .content-region .expand-collapse-glyph,
      :host .positioning-region:hover .content-region .start,
      :host .positioning-region:hover .content-region .end,
      :host([selected]) .content-region .expand-collapse-glyph,
      :host([selected]) .content-region .start,
      :host([selected]) .content-region .end {
        fill: ${nX.HighlightText};
      }
      :host([selected])::after {
        background: ${nX.Field};
      }
      :host(:${rm}) .positioning-region {
        border-color: ${nX.FieldText};
        box-shadow: 0 0 0 2px inset ${nX.Field};
        color: ${nX.FieldText};
      }
      :host([disabled]) .content-region,
      :host([disabled]) .positioning-region:hover .content-region {
        opacity: 1;
        color: ${nX.GrayText};
      }
      :host([disabled]) .content-region .expand-collapse-glyph,
      :host([disabled]) .content-region .start,
      :host([disabled]) .content-region .end,
      :host([disabled])
        .positioning-region:hover
        .content-region
        .expand-collapse-glyph,
      :host([disabled]) .positioning-region:hover .content-region .start,
      :host([disabled]) .positioning-region:hover .content-region .end {
        fill: ${nX.GrayText};
      }
      :host([disabled]) .positioning-region:hover {
        background: ${nX.Field};
      }
      .expand-collapse-glyph,
      .start,
      .end {
        fill: ${nX.FieldText};
      }
      :host(.nested) .expand-collapse-button:hover {
        background: ${nX.Field};
      }
      :host(.nested) .expand-collapse-button:hover .expand-collapse-glyph {
        fill: ${nX.FieldText};
      }
    `));class dg extends dn{}let dm=dg.compose({baseName:"tree-item",baseClass:dn,template:(e,t)=>s4`
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
        ${nu({property:"childItems",filter:rs()})}
    >
        <div class="positioning-region" part="positioning-region">
            <div class="content-region" part="content-region">
                ${rG(e=>e.childItems&&e.childItemLength()>0,s4`
                        <div
                            aria-hidden="true"
                            class="expand-collapse-button"
                            part="expand-collapse-button"
                            @click="${(e,t)=>e.handleExpandCollapseButtonClick(t.event)}"
                            ${s9("expandCollapseButton")}
                        >
                            <slot name="expand-collapse-glyph">
                                ${t.expandCollapseGlyph||""}
                            </slot>
                        </div>
                    `)}
                ${re(e,t)}
                <slot></slot>
                ${s7(e,t)}
            </div>
        </div>
        ${rG(e=>e.childItems&&e.childItemLength()>0&&(e.expanded||e.renderCollapsedChildren),s4`
                <div role="group" class="items" part="items">
                    <slot name="item" ${rn("items")}></slot>
                </div>
            `)}
    </template>
`,styles:dp,expandCollapseGlyph:`
        <svg
            viewBox="0 0 16 16"
            xmlns="http://www.w3.org/2000/svg"
            class="expand-collapse-glyph"
        >
            <path
                d="M5.00001 12.3263C5.00124 12.5147 5.05566 12.699 5.15699 12.8578C5.25831 13.0167 5.40243 13.1437 5.57273 13.2242C5.74304 13.3047 5.9326 13.3354 6.11959 13.3128C6.30659 13.2902 6.4834 13.2152 6.62967 13.0965L10.8988 8.83532C11.0739 8.69473 11.2153 8.51658 11.3124 8.31402C11.4096 8.11146 11.46 7.88966 11.46 7.66499C11.46 7.44033 11.4096 7.21853 11.3124 7.01597C11.2153 6.81341 11.0739 6.63526 10.8988 6.49467L6.62967 2.22347C6.48274 2.10422 6.30501 2.02912 6.11712 2.00691C5.92923 1.9847 5.73889 2.01628 5.56823 2.09799C5.39757 2.17969 5.25358 2.30817 5.153 2.46849C5.05241 2.62882 4.99936 2.8144 5.00001 3.00369V12.3263Z"
            />
        </svg>
    `});class db extends sR{constructor(){super(...arguments),this.currentFocused=null,this.handleFocus=e=>{if(!(this.slottedTreeItems.length<1)){if(e.target===this){null===this.currentFocused&&(this.currentFocused=this.getValidFocusableItem()),null!==this.currentFocused&&dn.focusItem(this.currentFocused);return}this.contains(e.target)&&(this.setAttribute("tabindex","-1"),this.currentFocused=e.target)}},this.handleBlur=e=>{e.target instanceof HTMLElement&&(null===e.relatedTarget||!this.contains(e.relatedTarget))&&this.setAttribute("tabindex","0")},this.handleKeyDown=e=>{if(e.defaultPrevented)return;if(this.slottedTreeItems.length<1)return!0;let t=this.getVisibleNodes();switch(e.key){case"Home":t.length&&dn.focusItem(t[0]);return;case"End":t.length&&dn.focusItem(t[t.length-1]);return;case oZ:if(e.target&&this.isFocusableElement(e.target)){let t=e.target;t instanceof dn&&t.childItemLength()>0&&t.expanded?t.expanded=!1:t instanceof dn&&t.parentElement instanceof dn&&dn.focusItem(t.parentElement)}return!1;case oJ:if(e.target&&this.isFocusableElement(e.target)){let t=e.target;t instanceof dn&&t.childItemLength()>0&&!t.expanded?t.expanded=!0:t instanceof dn&&t.childItemLength()>0&&this.focusNextNode(1,e.target)}return;case oQ:e.target&&this.isFocusableElement(e.target)&&this.focusNextNode(1,e.target);return;case o0:e.target&&this.isFocusableElement(e.target)&&this.focusNextNode(-1,e.target);return;case"Enter":this.handleClick(e);return}return!0},this.handleSelectedChange=e=>{if(e.defaultPrevented)return;if(!(e.target instanceof Element)||!da(e.target))return!0;let t=e.target;t.selected?(this.currentSelected&&this.currentSelected!==t&&(this.currentSelected.selected=!1),this.currentSelected=t):t.selected||this.currentSelected!==t||(this.currentSelected=null)},this.setItems=()=>{let e=this.treeView.querySelector("[aria-selected='true']");this.currentSelected=e,null!==this.currentFocused&&this.contains(this.currentFocused)||(this.currentFocused=this.getValidFocusableItem()),this.nested=this.checkForNestedItems(),this.getVisibleNodes().forEach(e=>{da(e)&&(e.nested=this.nested)})},this.isFocusableElement=e=>da(e),this.isSelectedElement=e=>e.selected}slottedTreeItemsChanged(){this.$fastController.isConnected&&this.setItems()}connectedCallback(){super.connectedCallback(),this.setAttribute("tabindex","0"),el.queueUpdate(()=>{this.setItems()})}handleClick(e){if(e.defaultPrevented)return;if(!(e.target instanceof Element)||!da(e.target))return!0;let t=e.target;t.disabled||(t.selected=!t.selected)}focusNextNode(e,t){let i=this.getVisibleNodes();if(!i)return;let o=i[i.indexOf(t)+e];rg(o)&&dn.focusItem(o)}getValidFocusableItem(){let e=this.getVisibleNodes(),t=e.findIndex(this.isSelectedElement);return(-1===t&&(t=e.findIndex(this.isFocusableElement)),-1!==t)?e[t]:null}checkForNestedItems(){return this.slottedTreeItems.some(e=>da(e)&&e.querySelector("[role='treeitem']"))}getVisibleNodes(){return function(e,t){if(e&&t&&rg(e))return Array.from(e.querySelectorAll(t)).filter(e=>null!==e.offsetParent)}(this,"[role='treeitem']")||[]}}X([eS({attribute:"render-collapsed-nodes"})],db.prototype,"renderCollapsedNodes",void 0),X([eu],db.prototype,"currentSelected",void 0),X([eu],db.prototype,"slottedTreeItems",void 0);let df=(e,t)=>rh`
  ${ru("flex")} :host {
    flex-direction: column;
    align-items: stretch;
    min-width: fit-content;
    font-size: 0;
  }

  :host:focus-visible {
    outline: none;
  }
`;class dv extends db{handleClick(e){if(e.defaultPrevented)return;if(!(e.target instanceof Element))return!0;let t=e.target;for(;t&&!da(t);)(t=t.parentElement)===this&&(t=null);t&&!t.disabled&&(t.selected=!0)}}let dy=dv.compose({baseName:"tree-view",baseClass:db,template:(e,t)=>s4`
    <template
        role="tree"
        ${s9("treeView")}
        @keydown="${(e,t)=>e.handleKeyDown(t.event)}"
        @focusin="${(e,t)=>e.handleFocus(t.event)}"
        @focusout="${(e,t)=>e.handleBlur(t.event)}"
        @click="${(e,t)=>e.handleClick(t.event)}"
        @selected-change="${(e,t)=>e.handleSelectedChange(t.event)}"
    >
        <slot ${rn("slottedTreeItems")}></slot>
    </template>
`,styles:df}),dx={horizontalDefaultPosition:"center",horizontalPositioningMode:"locktodefault",horizontalInset:!1,horizontalScaling:"anchor"},d$=Object.assign(Object.assign({},dx),{verticalDefaultPosition:"top",verticalPositioningMode:"locktodefault",verticalInset:!1,verticalScaling:"content"}),dw=Object.assign(Object.assign({},dx),{verticalDefaultPosition:"bottom",verticalPositioningMode:"locktodefault",verticalInset:!1,verticalScaling:"content"}),dk=Object.assign(Object.assign({},dx),{verticalPositioningMode:"dynamic",verticalInset:!1,verticalScaling:"content"}),dC=Object.assign(Object.assign({},d$),{verticalScaling:"fill"}),dI=Object.assign(Object.assign({},dw),{verticalScaling:"fill"}),dT=Object.assign(Object.assign({},dk),{verticalScaling:"fill"}),dF=s4`
    <template>
        ${e=>e.value}
    </template>
`;class dS extends sR{contentsTemplateChanged(){this.$fastController.isConnected&&this.updateView()}connectedCallback(){super.connectedCallback(),this.updateView()}disconnectedCallback(){super.disconnectedCallback(),this.disconnectView()}handleClick(e){return!e.defaultPrevented&&(this.handleInvoked(),!1)}handleInvoked(){this.$emit("pickeroptioninvoked")}updateView(){var e,t;this.disconnectView(),this.customView=null!=(t=null==(e=this.contentsTemplate)?void 0:e.render(this,this))?t:dF.render(this,this)}disconnectView(){var e;null==(e=this.customView)||e.dispose(),this.customView=void 0}}X([eS({attribute:"value"})],dS.prototype,"value",void 0),X([eu],dS.prototype,"contentsTemplate",void 0);let dO=s4`
    <template>
        ${e=>e.value}
    </template>
`;class dD extends sR{contentsTemplateChanged(){this.$fastController.isConnected&&this.updateView()}connectedCallback(){super.connectedCallback(),this.updateView()}disconnectedCallback(){this.disconnectView(),super.disconnectedCallback()}handleKeyDown(e){return!e.defaultPrevented&&("Enter"!==e.key||(this.handleInvoke(),!1))}handleClick(e){return e.defaultPrevented||this.handleInvoke(),!1}handleInvoke(){this.$emit("pickeriteminvoked")}updateView(){var e,t;this.disconnectView(),this.customView=null!=(t=null==(e=this.contentsTemplate)?void 0:e.render(this,this))?t:dO.render(this,this)}disconnectView(){var e;null==(e=this.customView)||e.dispose(),this.customView=void 0}}X([eS({attribute:"value"})],dD.prototype,"value",void 0),X([eu],dD.prototype,"contentsTemplate",void 0);class dE extends sR{}class dR extends ah(dE){constructor(){super(...arguments),this.proxy=document.createElement("input")}}let dL=s4`
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
        ${s9("inputElement")}
    ></input>
`;class dA extends dR{constructor(){super(...arguments),this.selection="",this.filterSelected=!0,this.filterQuery=!0,this.noSuggestionsText="No suggestions available",this.suggestionsAvailableText="Suggestions available",this.loadingText="Loading suggestions",this.menuPlacement="bottom-fill",this.showLoading=!1,this.optionsList=[],this.filteredOptionsList=[],this.flyoutOpen=!1,this.menuFocusIndex=-1,this.showNoOptions=!1,this.selectedItems=[],this.inputElementView=null,this.handleTextInput=e=>{this.query=this.inputElement.value},this.handleInputClick=e=>{e.preventDefault(),this.toggleFlyout(!0)},this.setRegionProps=()=>{if(this.flyoutOpen){if(null===this.region||void 0===this.region)return void el.queueUpdate(this.setRegionProps);this.region.anchorElement=this.inputElement}},this.configLookup={top:d$,bottom:dw,tallest:dk,"top-fill":dC,"bottom-fill":dI,"tallest-fill":dT}}selectionChanged(){this.$fastController.isConnected&&(this.handleSelectionChange(),this.proxy instanceof HTMLInputElement&&(this.proxy.value=this.selection,this.validate()))}optionsChanged(){this.optionsList=this.options.split(",").map(e=>e.trim()).filter(e=>""!==e)}menuPlacementChanged(){this.$fastController.isConnected&&this.updateMenuConfig()}showLoadingChanged(){this.$fastController.isConnected&&el.queueUpdate(()=>{this.setFocusedOption(0)})}listItemTemplateChanged(){this.updateListItemTemplate()}defaultListItemTemplateChanged(){this.updateListItemTemplate()}menuOptionTemplateChanged(){this.updateOptionTemplate()}defaultMenuOptionTemplateChanged(){this.updateOptionTemplate()}optionsListChanged(){this.updateFilteredOptions()}queryChanged(){this.$fastController.isConnected&&(this.inputElement.value!==this.query&&(this.inputElement.value=this.query),this.updateFilteredOptions(),this.$emit("querychange",{bubbles:!1}))}filteredOptionsListChanged(){this.$fastController.isConnected&&(this.showNoOptions=0===this.filteredOptionsList.length&&0===this.menuElement.querySelectorAll('[role="listitem"]').length,this.setFocusedOption(this.showNoOptions?-1:0))}flyoutOpenChanged(){this.flyoutOpen?(el.queueUpdate(this.setRegionProps),this.$emit("menuopening",{bubbles:!1})):this.$emit("menuclosing",{bubbles:!1})}showNoOptionsChanged(){this.$fastController.isConnected&&el.queueUpdate(()=>{this.setFocusedOption(0)})}connectedCallback(){super.connectedCallback(),this.listElement=document.createElement(this.selectedListTag),this.appendChild(this.listElement),this.itemsPlaceholderElement=document.createComment(""),this.listElement.append(this.itemsPlaceholderElement),this.inputElementView=dL.render(this,this.listElement);let e=this.menuTag.toUpperCase();this.menuElement=Array.from(this.children).find(t=>t.tagName===e),void 0===this.menuElement&&(this.menuElement=document.createElement(this.menuTag),this.appendChild(this.menuElement)),""===this.menuElement.id&&(this.menuElement.id=aR("listbox-")),this.menuId=this.menuElement.id,this.optionsPlaceholder=document.createComment(""),this.menuElement.append(this.optionsPlaceholder),this.updateMenuConfig(),el.queueUpdate(()=>this.initialize())}disconnectedCallback(){super.disconnectedCallback(),this.toggleFlyout(!1),this.inputElement.removeEventListener("input",this.handleTextInput),this.inputElement.removeEventListener("click",this.handleInputClick),null!==this.inputElementView&&(this.inputElementView.dispose(),this.inputElementView=null)}focus(){this.inputElement.focus()}initialize(){this.updateListItemTemplate(),this.updateOptionTemplate(),this.itemsRepeatBehavior=new nh(e=>e.selectedItems,e=>e.activeListItemTemplate,{positioning:!0}).createBehavior(this.itemsPlaceholderElement),this.inputElement.addEventListener("input",this.handleTextInput),this.inputElement.addEventListener("click",this.handleInputClick),this.$fastController.addBehaviors([this.itemsRepeatBehavior]),this.menuElement.suggestionsAvailableText=this.suggestionsAvailableText,this.menuElement.addEventListener("optionsupdated",this.handleMenuOptionsUpdated),this.optionsRepeatBehavior=new nh(e=>e.filteredOptionsList,e=>e.activeMenuOptionTemplate,{positioning:!0}).createBehavior(this.optionsPlaceholder),this.$fastController.addBehaviors([this.optionsRepeatBehavior]),this.handleSelectionChange()}toggleFlyout(e){if(this.flyoutOpen!==e){if(e&&document.activeElement===this.inputElement){this.flyoutOpen=e,el.queueUpdate(()=>{void 0!==this.menuElement?this.setFocusedOption(0):this.disableMenu()});return}this.flyoutOpen=!1,this.disableMenu()}}handleMenuOptionsUpdated(e){e.preventDefault(),this.flyoutOpen&&this.setFocusedOption(0)}handleKeyDown(e){if(e.defaultPrevented)return!1;switch(e.key){case oQ:if(this.flyoutOpen){let e=this.flyoutOpen?Math.min(this.menuFocusIndex+1,this.menuElement.optionElements.length-1):0;this.setFocusedOption(e)}else this.toggleFlyout(!0);return!1;case o0:if(this.flyoutOpen){let e=this.flyoutOpen?Math.max(this.menuFocusIndex-1,0):0;this.setFocusedOption(e)}else this.toggleFlyout(!0);return!1;case"Escape":return this.toggleFlyout(!1),!1;case"Enter":return -1!==this.menuFocusIndex&&this.menuElement.optionElements.length>this.menuFocusIndex&&this.menuElement.optionElements[this.menuFocusIndex].click(),!1;case oJ:if(document.activeElement!==this.inputElement)return this.incrementFocusedItem(1),!1;return!0;case oZ:if(0===this.inputElement.selectionStart)return this.incrementFocusedItem(-1),!1;return!0;case"Delete":case"Backspace":{if(null===document.activeElement)return!0;if(document.activeElement===this.inputElement){if(0===this.inputElement.selectionStart)return this.selection=this.selectedItems.slice(0,this.selectedItems.length-1).toString(),this.toggleFlyout(!1),!1;return!0}let e=Array.from(this.listElement.children),t=e.indexOf(document.activeElement);if(t>-1)return this.selection=this.selectedItems.splice(t,1).toString(),el.queueUpdate(()=>{e[Math.min(e.length,t)].focus()}),!1;return!0}}return this.toggleFlyout(!0),!0}handleFocusIn(e){return!1}handleFocusOut(e){return void 0!==this.menuElement&&this.menuElement.contains(e.relatedTarget)||this.toggleFlyout(!1),!1}handleSelectionChange(){this.selectedItems.toString()!==this.selection&&(this.selectedItems=""===this.selection?[]:this.selection.split(","),this.updateFilteredOptions(),el.queueUpdate(()=>{this.checkMaxItems()}),this.$emit("selectionchange",{bubbles:!1}))}handleRegionLoaded(e){el.queueUpdate(()=>{this.setFocusedOption(0),this.$emit("menuloaded",{bubbles:!1})})}checkMaxItems(){if(void 0!==this.inputElement)if(void 0!==this.maxSelected&&this.selectedItems.length>=this.maxSelected){if(document.activeElement===this.inputElement){let e=Array.from(this.listElement.querySelectorAll("[role='listitem']"));e[e.length-1].focus()}this.inputElement.hidden=!0}else this.inputElement.hidden=!1}handleItemInvoke(e){if(e.defaultPrevented)return!1;if(e.target instanceof dD){let t=Array.from(this.listElement.querySelectorAll("[role='listitem']")).indexOf(e.target);if(-1!==t){let e=this.selectedItems.slice();e.splice(t,1),this.selection=e.toString(),el.queueUpdate(()=>this.incrementFocusedItem(0))}return!1}return!0}handleOptionInvoke(e){return!e.defaultPrevented&&(!(e.target instanceof dS)||(void 0!==e.target.value&&(this.selection=`${this.selection}${""===this.selection?"":","}${e.target.value}`),this.inputElement.value="",this.query="",this.inputElement.focus(),this.toggleFlyout(!1),!1))}incrementFocusedItem(e){if(0===this.selectedItems.length)return void this.inputElement.focus();let t=Array.from(this.listElement.querySelectorAll("[role='listitem']"));if(null!==document.activeElement){let i=t.indexOf(document.activeElement);-1===i&&(i=t.length);let o=Math.min(t.length,Math.max(0,i+e));o===t.length?void 0!==this.maxSelected&&this.selectedItems.length>=this.maxSelected?t[o-1].focus():this.inputElement.focus():t[o].focus()}}disableMenu(){var e,t,i;this.menuFocusIndex=-1,this.menuFocusOptionId=void 0,null==(e=this.inputElement)||e.removeAttribute("aria-activedescendant"),null==(t=this.inputElement)||t.removeAttribute("aria-owns"),null==(i=this.inputElement)||i.removeAttribute("aria-expanded")}setFocusedOption(e){if(!this.flyoutOpen||-1===e||this.showNoOptions||this.showLoading)return void this.disableMenu();if(0===this.menuElement.optionElements.length)return;this.menuElement.optionElements.forEach(e=>{e.setAttribute("aria-selected","false")}),this.menuFocusIndex=e,this.menuFocusIndex>this.menuElement.optionElements.length-1&&(this.menuFocusIndex=this.menuElement.optionElements.length-1),this.menuFocusOptionId=this.menuElement.optionElements[this.menuFocusIndex].id,this.inputElement.setAttribute("aria-owns",this.menuId),this.inputElement.setAttribute("aria-expanded","true"),this.inputElement.setAttribute("aria-activedescendant",this.menuFocusOptionId);let t=this.menuElement.optionElements[this.menuFocusIndex];t.setAttribute("aria-selected","true"),this.menuElement.scrollTo(0,t.offsetTop)}updateListItemTemplate(){var e;this.activeListItemTemplate=null!=(e=this.listItemTemplate)?e:this.defaultListItemTemplate}updateOptionTemplate(){var e;this.activeMenuOptionTemplate=null!=(e=this.menuOptionTemplate)?e:this.defaultMenuOptionTemplate}updateFilteredOptions(){this.filteredOptionsList=this.optionsList.slice(0),this.filterSelected&&(this.filteredOptionsList=this.filteredOptionsList.filter(e=>-1===this.selectedItems.indexOf(e))),this.filterQuery&&""!==this.query&&void 0!==this.query&&(this.filteredOptionsList=this.filteredOptionsList.filter(e=>-1!==e.indexOf(this.query)))}updateMenuConfig(){let e=this.configLookup[this.menuPlacement];null===e&&(e=dI),this.menuConfig=Object.assign(Object.assign({},e),{autoUpdateMode:"auto",fixedPlacement:!0,horizontalViewportLock:!1,verticalViewportLock:!1})}}X([eS({attribute:"selection"})],dA.prototype,"selection",void 0),X([eS({attribute:"options"})],dA.prototype,"options",void 0),X([eS({attribute:"filter-selected",mode:"boolean"})],dA.prototype,"filterSelected",void 0),X([eS({attribute:"filter-query",mode:"boolean"})],dA.prototype,"filterQuery",void 0),X([eS({attribute:"max-selected"})],dA.prototype,"maxSelected",void 0),X([eS({attribute:"no-suggestions-text"})],dA.prototype,"noSuggestionsText",void 0),X([eS({attribute:"suggestions-available-text"})],dA.prototype,"suggestionsAvailableText",void 0),X([eS({attribute:"loading-text"})],dA.prototype,"loadingText",void 0),X([eS({attribute:"label"})],dA.prototype,"label",void 0),X([eS({attribute:"labelledby"})],dA.prototype,"labelledBy",void 0),X([eS({attribute:"placeholder"})],dA.prototype,"placeholder",void 0),X([eS({attribute:"menu-placement"})],dA.prototype,"menuPlacement",void 0),X([eu],dA.prototype,"showLoading",void 0),X([eu],dA.prototype,"listItemTemplate",void 0),X([eu],dA.prototype,"defaultListItemTemplate",void 0),X([eu],dA.prototype,"activeListItemTemplate",void 0),X([eu],dA.prototype,"menuOptionTemplate",void 0),X([eu],dA.prototype,"defaultMenuOptionTemplate",void 0),X([eu],dA.prototype,"activeMenuOptionTemplate",void 0),X([eu],dA.prototype,"listItemContentsTemplate",void 0),X([eu],dA.prototype,"menuOptionContentsTemplate",void 0),X([eu],dA.prototype,"optionsList",void 0),X([eu],dA.prototype,"query",void 0),X([eu],dA.prototype,"filteredOptionsList",void 0),X([eu],dA.prototype,"flyoutOpen",void 0),X([eu],dA.prototype,"menuId",void 0),X([eu],dA.prototype,"selectedListTag",void 0),X([eu],dA.prototype,"menuTag",void 0),X([eu],dA.prototype,"menuFocusIndex",void 0),X([eu],dA.prototype,"menuFocusOptionId",void 0),X([eu],dA.prototype,"showNoOptions",void 0),X([eu],dA.prototype,"menuConfig",void 0),X([eu],dA.prototype,"selectedItems",void 0);class dV extends sR{constructor(){super(...arguments),this.optionElements=[]}menuElementsChanged(){this.updateOptions()}headerElementsChanged(){this.updateOptions()}footerElementsChanged(){this.updateOptions()}updateOptions(){this.optionElements.splice(0,this.optionElements.length),this.addSlottedListItems(this.headerElements),this.addSlottedListItems(this.menuElements),this.addSlottedListItems(this.footerElements),this.$emit("optionsupdated",{bubbles:!1})}addSlottedListItems(e){void 0!==e&&e.forEach(e=>{1===e.nodeType&&"listitem"===e.getAttribute("role")&&(e.id=e.id||aR("option-"),this.optionElements.push(e))})}}X([eu],dV.prototype,"menuElements",void 0),X([eu],dV.prototype,"headerElements",void 0),X([eu],dV.prototype,"footerElements",void 0),X([eu],dV.prototype,"suggestionsAvailableText",void 0);class dP extends sR{}let dH=(e,t)=>rh`
  .region {
    z-index: 1000;
    overflow: hidden;
    display: flex;
    font-family: ${tu};
    font-size: ${tC};
  }

  .loaded {
    opacity: 1;
    pointer-events: none;
  }

  .loading-display,
  .no-options-display {
    background: ${iT};
    width: 100%;
    min-height: calc(${ry} * 1px);
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-items: center;
    padding: calc(${tv} * 1px);
  }

  .loading-progress {
    width: 42px;
    height: 42px;
  }

  .bottom {
    flex-direction: column;
  }

  .top {
    flex-direction: column-reverse;
  }
`,dz=(e,t)=>rh`
    :host {
      background: ${iT};
      --elevation: 11;
      /* TODO: a mechanism to manage z-index across components
            https://github.com/microsoft/fast/issues/3813 */
      z-index: 1000;
      display: flex;
      width: 100%;
      max-height: 100%;
      min-height: 58px;
      box-sizing: border-box;
      flex-direction: column;
      overflow-y: auto;
      overflow-x: hidden;
      pointer-events: auto;
      border-radius: calc(${tb} * 1px);
      padding: calc(${tv} * 1px) 0;
      border: calc(${tw} * 1px) solid transparent;
      ${ax}
    }

    .suggestions-available-alert {
      height: 0;
      opacity: 0;
      overflow: hidden;
    }
  `.withBehaviors(rv(rh`
      :host {
        background: ${nX.Canvas};
        border-color: ${nX.CanvasText};
      }
    `)),dM=(e,t)=>rh`
    :host {
      display: flex;
      align-items: center;
      justify-items: center;
      font-family: ${tu};
      border-radius: calc(${tb} * 1px);
      border: calc(${tk} * 1px) solid transparent;
      box-sizing: border-box;
      background: ${i4};
      color: ${op};
      cursor: pointer;
      fill: currentcolor;
      font-size: ${tC};
      min-height: calc(${ry} * 1px);
      line-height: ${tI};
      margin: 0 calc(${tv} * 1px);
      outline: none;
      overflow: hidden;
      padding: 0 calc(${tv} * 2.25px);
      user-select: none;
      white-space: nowrap;
    }

    :host(:${rm}[role="listitem"]) {
      border-color: ${on};
      background: ${i8};
    }

    :host(:hover) {
      background: ${i6};
    }

    :host(:active) {
      background: ${i9};
    }

    :host([aria-selected='true']) {
      background: ${iS};
      color: ${iA};
    }

    :host([aria-selected='true']:hover) {
      background: ${iO};
      color: ${iV};
    }

    :host([aria-selected='true']:active) {
      background: ${iD};
      color: ${iP};
    }
  `.withBehaviors(rv(rh`
      :host {
        border-color: transparent;
        forced-color-adjust: none;
        color: ${nX.ButtonText};
        fill: currentcolor;
      }

      :host(:not([aria-selected='true']):hover),
      :host([aria-selected='true']) {
        background: ${nX.Highlight};
        color: ${nX.HighlightText};
      }

      :host([disabled]),
      :host([disabled]:not([aria-selected='true']):hover) {
        background: ${nX.Canvas};
        color: ${nX.GrayText};
        fill: currentcolor;
        opacity: 1;
      }
    `)),dN=(e,t)=>rh`
    :host {
      display: flex;
      align-items: center;
      justify-items: center;
      font-family: ${tu};
      border-radius: calc(${tb} * 1px);
      border: calc(${tk} * 1px) solid transparent;
      box-sizing: border-box;
      background: ${i4};
      color: ${op};
      cursor: pointer;
      fill: currentcolor;
      font-size: ${tC};
      height: calc(${ry} * 1px);
      line-height: ${tI};
      outline: none;
      overflow: hidden;
      padding: 0 calc(${tv} * 2.25px);
      user-select: none;
      white-space: nowrap;
    }

    :host(:hover) {
      background: ${i6};
    }

    :host(:active) {
      background: ${i9};
    }

    :host(:${rm}) {
      background: ${i8};
      border-color: ${on};
    }

    :host([aria-selected='true']) {
      background: ${iS};
      color: ${iP};
    }
  `.withBehaviors(rv(rh`
      :host {
        border-color: transparent;
        forced-color-adjust: none;
        color: ${nX.ButtonText};
        fill: currentcolor;
      }

      :host(:not([aria-selected='true']):hover),
      :host([aria-selected='true']) {
        background: ${nX.Highlight};
        color: ${nX.HighlightText};
      }

      :host([disabled]),
      :host([disabled]:not([aria-selected='true']):hover) {
        background: ${nX.Canvas};
        color: ${nX.GrayText};
        fill: currentcolor;
        opacity: 1;
      }
    `));class dB extends dA{}let dj=dB.compose({baseName:"draft-picker",baseClass:dA,template:(e,t)=>{let i,o,s=e.tagFor(rq),r=e.tagFor(dV),a=e.tagFor(dP),n=e.tagFor(dP),l=(i=e.tagFor(dD),s4`
    <${i}
        value="${e=>e}"
        :contentsTemplate="${(e,t)=>t.parent.listItemContentsTemplate}"
    >
    </${i}>
    `),h=(o=e.tagFor(dS),s4`
    <${o}
        value="${e=>e}"
        :contentsTemplate="${(e,t)=>t.parent.menuOptionContentsTemplate}"
    >
    </${o}>
    `);return s4`
        <template
            :selectedListTag="${()=>a}"
            :menuTag="${()=>r}"
            :defaultListItemTemplate="${l}"
            :defaultMenuOptionTemplate="${h}"
            @focusin="${(e,t)=>e.handleFocusIn(t.event)}"
            @focusout="${(e,t)=>e.handleFocusOut(t.event)}"
            @keydown="${(e,t)=>e.handleKeyDown(t.event)}"
            @pickeriteminvoked="${(e,t)=>e.handleItemInvoke(t.event)}"
            @pickeroptioninvoked="${(e,t)=>e.handleOptionInvoke(t.event)}"
        >
            <slot name="list-region"></slot>

            ${rG(e=>e.flyoutOpen,s4`
                <${s}
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
                    ${s9("region")}
                >
                    ${rG(e=>!e.showNoOptions&&!e.showLoading,s4`
                            <slot name="menu-region"></slot>
                        `)}
                    ${rG(e=>e.showNoOptions&&!e.showLoading,s4`
                            <div class="no-options-display" part="no-options-display">
                                <slot name="no-options-region">
                                    ${e=>e.noSuggestionsText}
                                </slot>
                            </div>
                        `)}
                    ${rG(e=>e.showLoading,s4`
                            <div class="loading-display" part="loading-display">
                                <slot name="loading-region">
                                    <${n}
                                        part="loading-progress"
                                        class="loading-progress
                                        slot="loading-region"
                                    ></${n}>
                                        ${e=>e.loadingText}
                                </slot>
                            </div>
                        `)}
                </${s}>
            `)}
        </template>
    `},styles:dH,shadowOptions:{}});class dq extends dV{connectedCallback(){iT.setValueFor(this,ib),super.connectedCallback()}}let dU=dq.compose({baseName:"draft-picker-menu",baseClass:dV,template:(e,t)=>s4`
        <template role="list" slot="menu-region">
            <div class="options-display" part="options-display">
                <div class="header-region" part="header-region">
                    <slot name="header-region" ${rn("headerElements")}></slot>
                </div>

                <slot ${rn("menuElements")}></slot>
                <div class="footer-region" part="footer-region">
                    <slot name="footer-region" ${rn("footerElements")}></slot>
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
    `,styles:dz});class d_ extends dS{}let dG=d_.compose({baseName:"draft-picker-menu-option",baseClass:dS,template:(e,t)=>s4`
        <template
            role="listitem"
            tabindex="-1"
            @click="${(e,t)=>e.handleClick(t.event)}"
        >
            <slot></slot>
        </template>
    `,styles:dM});class dW extends dP{}let dK=dW.compose({baseName:"draft-picker-list",baseClass:dP,template:(e,t)=>s4`
        <template slot="list-region" role="list" class="picker-list">
            <slot></slot>
            <slot name="input-region"></slot>
        </template>
    `,styles:(e,t)=>rh`
        :host {
            display: flex;
            flex-direction: row;
            column-gap: calc(${tv} * 1px);
            row-gap: calc(${tv} * 1px);
            flex-wrap: wrap;
        }

        ::slotted([role="combobox"]) {
            min-width: 260px;
            width: auto;
            box-sizing: border-box;
            color: ${op};
            background: ${i0};
            border-radius: calc(${tb} * 1px);
            border: calc(${tw} * 1px) solid ${iS};
            height: calc(${ry} * 1px);
            font-family: ${tu};
            outline: none;
            user-select: none;
            font-size: ${tC};
            line-height: ${tI};
            padding: 0 calc(${tv} * 2px + 1px);
        }

        ::slotted([role="combobox"]:active) { {
            background: ${i1};
            border-color: ${iD};
        }

        ::slotted([role="combobox"]:focus-within) {
            border-color: ${on};
            box-shadow: 0 0 0 1px ${on} inset;
        }
    `.withBehaviors(rv(rh`
      ::slotted([role='combobox']:active) {
        background: ${nX.Field};
        border-color: ${nX.Highlight};
      }
      ::slotted([role='combobox']:focus-within) {
        border-color: ${nX.Highlight};
        box-shadow: 0 0 0 1px ${nX.Highlight} inset;
      }
      ::slotted(input:placeholder) {
        color: ${nX.GrayText};
      }
    `))});class dX extends dD{}let dY=dX.compose({baseName:"draft-picker-list-item",baseClass:dD,template:(e,t)=>s4`
        <template
            role="listitem"
            tabindex="0"
            @click="${(e,t)=>e.handleClick(t.event)}"
            @keydown="${(e,t)=>e.handleKeyDown(t.event)}"
        >
            <slot></slot>
        </template>
    `,styles:dN}),dQ={jpAccordion:rC,jpAccordionItem:rw,jpAnchor:rB,jpAnchoredRegion:rX,jpAvatar:r5,jpBadge:r6,jpBreadcrumb:at,jpBreadcrumbItem:as,jpButton:av,jpCard:ak,jpCheckbox:aD,jpCombobox:aJ,jpDataGrid:nw,jpDataGridCell:nv,jpDataGridRow:nx,jpDateField:nP,jpDesignSystemProvider:nq,jpDialog:lu,jpDisclosure:lb,jpDivider:lw,jpListbox:lC,jpMenu:lR,jpMenuItem:lV,jpNumberField:lB,jpOption:lU,jpPicker:dj,jpPickerList:dK,jpPickerListItem:dY,jpPickerMenu:dU,jpPickerMenuOption:dG,jpProgress:lK,jpProgressRing:lQ,jpRadio:l3,jpRadioGroup:l8,jpSearch:hn,jpSelect:hh,jpSkeleton:hp,jpSlider:hw,jpSliderLabel:hO,jpSwitch:hV,jpTab:hq,jpTabPanel:hM,jpTabs:hK,jpTextArea:h1,jpTextField:h3,jpToolbar:de,jpTooltip:dr,jpTreeItem:dm,jpTreeView:dy,register(e,...t){if(e)for(let i in this)"register"!==i&&this[i]().register(e,...t)}},dZ=Object.freeze({definitionCallbackOnly:null,ignoreDuplicate:Symbol()}),dJ=new Map,d0=new Map,d1=null,d2=st.createInterface(e=>e.cachedCallback(e=>(null===d1&&(d1=new d3(null,e)),d1))),d5=Object.freeze({tagFor:e=>d0.get(e),responsibleFor(e){let t=e.$$designSystem$$;return t||st.findResponsibleContainer(e).get(d2)},getOrCreate(e){if(!e)return null===d1&&(d1=st.getOrCreateDOMContainer().get(d2)),d1;let t=e.$$designSystem$$;if(t)return t;let i=st.getOrCreateDOMContainer(e);if(i.has(d2,!1))return i.get(d2);{let t=new d3(e,i);return i.register(sx.instance(d2,t)),t}}});class d3{constructor(e,t){this.owner=e,this.container=t,this.designTokensInitialized=!1,this.prefix="fast",this.shadowRootMode=void 0,this.disambiguate=()=>dZ.definitionCallbackOnly,null!==e&&(e.$$designSystem$$=this)}withPrefix(e){return this.prefix=e,this}withShadowRootMode(e){return this.shadowRootMode=e,this}withElementDisambiguation(e){return this.disambiguate=e,this}withDesignTokenRoot(e){return this.designTokenRoot=e,this}register(...e){let t=this.container,i=[],o=this.disambiguate,s=this.shadowRootMode,r={elementPrefix:this.prefix,tryDefineElement(e,r,a){let n="string"==typeof e?{name:e,type:r,callback:a}:e,{name:l,callback:h,baseClass:d}=n,{type:c}=n,u=l,p=dJ.get(u),g=!0;for(;p;){let e=o(u,c,p);switch(e){case dZ.ignoreDuplicate:return;case dZ.definitionCallbackOnly:g=!1,p=void 0;break;default:u=e,p=dJ.get(u)}}g&&((d0.has(c)||c===sR)&&(c=class extends c{}),dJ.set(u,c),d0.set(c,u),d&&d0.set(d,u)),i.push(new d4(t,u,c,s,h,g))}};for(let o of(this.designTokensInitialized||(this.designTokensInitialized=!0,null!==this.designTokenRoot&&e2.registerRoot(this.designTokenRoot)),t.registerWithContext(r,...e),i))o.callback(o),o.willDefine&&null!==o.definition&&o.definition.define();return this}}class d4{constructor(e,t,i,o,s,r){this.container=e,this.name=t,this.type=i,this.shadowRootMode=o,this.callback=s,this.willDefine=r,this.definition=null}definePresentation(e){sD.define(this.name,e,this.container)}defineElement(e){this.definition=new eR(this.type,Object.assign(Object.assign({},e),{name:this.name}))}tagFor(e){return d5.tagFor(e)}}function d6(e){return d5.getOrCreate(e).withPrefix("jp")}}}]);