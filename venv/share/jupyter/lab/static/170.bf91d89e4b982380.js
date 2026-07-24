"use strict";(self.rspackChunk_jupyterlab_application_top=self.rspackChunk_jupyterlab_application_top||[]).push([[170],{99313(e,t,o){function r(e,t,o){return isNaN(e)||e<=t?t:e>=o?o:e}function a(e,t,o){return isNaN(e)||e<=t?0:e>=o?1:e/(o-t)}function i(e,t,o){return isNaN(e)?t:t+e*(o-t)}function l(e){return Math.PI/180*e}function s(e,t,o){return isNaN(e)||e<=0?t:e>=1?o:t+e*(o-t)}function n(e,t,o){if(e<=0)return t%360;if(e>=1)return o%360;let r=(t-o+360)%360;return r<=(o-t+360)%360?(t-r*e+360)%360:(t+r*e+360)%360}function c(e,t){let o=Math.pow(10,t);return Math.round(e*o)/o}o.r(t),o.d(t,{Accordion:()=>o9,AccordionItem:()=>o3,Anchor:()=>rc,AnchoredRegion:()=>ru,Avatar:()=>rm,Badge:()=>ry,Breadcrumb:()=>rF,BreadcrumbItem:()=>rD,Button:()=>rz,Card:()=>rH,Checkbox:()=>rR,Combobox:()=>rE,DataGrid:()=>rK,DataGridCell:()=>rX,DataGridRow:()=>rY,DateField:()=>r3,DelegatesARIAToolbar:()=>iC,DesignSystemProvider:()=>ar,Dialog:()=>an,DirectionalStyleSheetBehavior:()=>rg,Disclosure:()=>ah,Divider:()=>ag,FoundationToolbar:()=>iF,Listbox:()=>af,Menu:()=>a$,MenuItem:()=>ak,NumberField:()=>aC,Option:()=>aS,PaletteRGB:()=>en,Picker:()=>iW,PickerList:()=>iK,PickerListItem:()=>i0,PickerMenu:()=>iX,PickerMenuOption:()=>iY,Progress:()=>aj,ProgressRing:()=>aH,Radio:()=>aR,RadioGroup:()=>aM,Search:()=>aW,Select:()=>aX,Skeleton:()=>aJ,Slider:()=>a2,SliderLabel:()=>a8,StandardLuminance:()=>A,SwatchRGB:()=>I,Switch:()=>it,Tab:()=>is,TabPanel:()=>ia,Tabs:()=>ih,TextArea:()=>ig,TextField:()=>iv,Toolbar:()=>iV,Tooltip:()=>iT,TreeItem:()=>iI,TreeView:()=>iA,accentColor:()=>tf,accentFillActive:()=>tH,accentFillActiveDelta:()=>eY,accentFillFocus:()=>tN,accentFillFocusDelta:()=>eJ,accentFillHover:()=>tL,accentFillHoverDelta:()=>eZ,accentFillRecipe:()=>tj,accentFillRest:()=>tB,accentFillRestDelta:()=>eX,accentForegroundActive:()=>tY,accentForegroundActiveDelta:()=>e0,accentForegroundFocus:()=>tJ,accentForegroundFocusDelta:()=>e1,accentForegroundHover:()=>tZ,accentForegroundHoverDelta:()=>eQ,accentForegroundRecipe:()=>tU,accentForegroundRest:()=>tX,accentForegroundRestDelta:()=>eK,accentPalette:()=>tm,accordionItemStyles:()=>o6,accordionStyles:()=>o2,addJupyterLabThemeChangeListener:()=>oJ,allComponents:()=>i2,anchorStyles:()=>rn,anchoredRegionStyles:()=>rh,applyJupyterTheme:()=>o0,avatarStyles:()=>rf,badgeStyles:()=>rx,baseHeightMultiplier:()=>e$,baseHorizontalSpacingMultiplier:()=>ex,baseLayerLuminance:()=>ey,bodyFont:()=>ev,breadcrumbItemStyles:()=>rV,breadcrumbStyles:()=>rw,buttonStyles:()=>rT,cardStyles:()=>rL,checkboxStyles:()=>rO,checkboxTemplate:()=>rI,comboboxStyles:()=>rG,controlCornerRadius:()=>ek,dataGridCellStyles:()=>rU,dataGridRowStyles:()=>rW,dataGridStyles:()=>rq,dateFieldStyles:()=>r9,dateFieldTemplate:()=>r8,density:()=>ew,designSystemProviderStyles:()=>ai,designSystemProviderTemplate:()=>aa,designUnit:()=>eF,dialogStyles:()=>as,direction:()=>eV,disabledOpacity:()=>eD,disclosureStyles:()=>ad,dividerStyles:()=>ap,elementScale:()=>eC,errorColor:()=>oV,errorFillActive:()=>oj,errorFillFocus:()=>oB,errorFillHover:()=>oz,errorFillRecipe:()=>oS,errorFillRest:()=>oT,errorForegroundActive:()=>oU,errorForegroundFocus:()=>oX,errorForegroundHover:()=>oW,errorForegroundRecipe:()=>o_,errorForegroundRest:()=>oq,errorPalette:()=>oD,fillColor:()=>tz,focusStrokeInner:()=>op,focusStrokeInnerRecipe:()=>ou,focusStrokeOuter:()=>oh,focusStrokeOuterRecipe:()=>od,focusStrokeWidth:()=>eT,foregroundOnAccentActive:()=>tA,foregroundOnAccentActiveLarge:()=>tq,foregroundOnAccentFocus:()=>tM,foregroundOnAccentFocusLarge:()=>tW,foregroundOnAccentHover:()=>tP,foregroundOnAccentHoverLarge:()=>t_,foregroundOnAccentLargeRecipe:()=>tG,foregroundOnAccentRecipe:()=>tI,foregroundOnAccentRest:()=>tR,foregroundOnAccentRestLarge:()=>tE,foregroundOnErrorActive:()=>oI,foregroundOnErrorActiveLarge:()=>oG,foregroundOnErrorFocus:()=>oR,foregroundOnErrorFocusLarge:()=>oE,foregroundOnErrorHover:()=>oO,foregroundOnErrorHoverLarge:()=>oM,foregroundOnErrorLargeRecipe:()=>oP,foregroundOnErrorRecipe:()=>oH,foregroundOnErrorRest:()=>oN,foregroundOnErrorRestLarge:()=>oA,heightNumberAsToken:()=>oC,horizontalSliderLabelStyles:()=>a3,imgTemplate:()=>rv,isDark:()=>G,jpAccordion:()=>o8,jpAccordionItem:()=>o4,jpAnchor:()=>rd,jpAnchoredRegion:()=>rp,jpAvatar:()=>r$,jpBadge:()=>rk,jpBreadcrumb:()=>rC,jpBreadcrumbItem:()=>rS,jpButton:()=>rj,jpCard:()=>rN,jpCheckbox:()=>rP,jpCombobox:()=>r_,jpDataGrid:()=>rQ,jpDataGridCell:()=>rZ,jpDataGridRow:()=>rJ,jpDateField:()=>r7,jpDesignSystemProvider:()=>al,jpDialog:()=>ac,jpDisclosure:()=>au,jpDivider:()=>ab,jpListbox:()=>am,jpMenu:()=>ax,jpMenuItem:()=>aw,jpNumberField:()=>aV,jpOption:()=>aT,jpPicker:()=>iU,jpPickerList:()=>iQ,jpPickerListItem:()=>i1,jpPickerMenu:()=>iZ,jpPickerMenuOption:()=>iJ,jpProgress:()=>aB,jpProgressRing:()=>aN,jpRadio:()=>aP,jpRadioGroup:()=>aG,jpSearch:()=>aU,jpSelect:()=>aZ,jpSkeleton:()=>aK,jpSlider:()=>a5,jpSliderLabel:()=>a7,jpSwitch:()=>io,jpTab:()=>ic,jpTabPanel:()=>ii,jpTabs:()=>iu,jpTextArea:()=>ib,jpTextField:()=>i$,jpToolbar:()=>iD,jpTooltip:()=>iz,jpTreeItem:()=>iR,jpTreeView:()=>iM,listboxStyles:()=>rA,menuItemStyles:()=>ay,menuStyles:()=>av,neutralColor:()=>tg,neutralFillActive:()=>t1,neutralFillActiveDelta:()=>e6,neutralFillFocus:()=>t2,neutralFillFocusDelta:()=>e3,neutralFillHover:()=>t0,neutralFillHoverDelta:()=>e5,neutralFillInputActive:()=>t4,neutralFillInputActiveDelta:()=>e8,neutralFillInputFocus:()=>t9,neutralFillInputFocusDelta:()=>e7,neutralFillInputHover:()=>t3,neutralFillInputHoverDelta:()=>e9,neutralFillInputRecipe:()=>t5,neutralFillInputRest:()=>t6,neutralFillInputRestDelta:()=>e4,neutralFillLayerRecipe:()=>on,neutralFillLayerRest:()=>oc,neutralFillLayerRestDelta:()=>tn,neutralFillRecipe:()=>tK,neutralFillRest:()=>tQ,neutralFillRestDelta:()=>e2,neutralFillStealthActive:()=>ot,neutralFillStealthActiveDelta:()=>to,neutralFillStealthFocus:()=>oo,neutralFillStealthFocusDelta:()=>tr,neutralFillStealthHover:()=>oe,neutralFillStealthHoverDelta:()=>tt,neutralFillStealthRecipe:()=>t8,neutralFillStealthRest:()=>t7,neutralFillStealthRestDelta:()=>te,neutralFillStrongActive:()=>ol,neutralFillStrongActiveDelta:()=>tl,neutralFillStrongFocus:()=>os,neutralFillStrongFocusDelta:()=>ts,neutralFillStrongHover:()=>oi,neutralFillStrongHoverDelta:()=>ti,neutralFillStrongRecipe:()=>or,neutralFillStrongRest:()=>oa,neutralFillStrongRestDelta:()=>ta,neutralForegroundHint:()=>ob,neutralForegroundHintRecipe:()=>og,neutralForegroundRecipe:()=>of,neutralForegroundRest:()=>om,neutralLayer1:()=>tw,neutralLayer1Recipe:()=>tk,neutralLayer2:()=>tC,neutralLayer2Recipe:()=>tF,neutralLayer3:()=>tD,neutralLayer3Recipe:()=>tV,neutralLayer4:()=>tT,neutralLayer4Recipe:()=>tS,neutralLayerCardContainer:()=>t$,neutralLayerCardContainerRecipe:()=>tv,neutralLayerFloating:()=>ty,neutralLayerFloatingRecipe:()=>tx,neutralPalette:()=>tb,neutralStrokeActive:()=>oy,neutralStrokeActiveDelta:()=>th,neutralStrokeDividerRecipe:()=>ow,neutralStrokeDividerRest:()=>oF,neutralStrokeDividerRestDelta:()=>tp,neutralStrokeFocus:()=>ok,neutralStrokeFocusDelta:()=>tu,neutralStrokeHover:()=>ox,neutralStrokeHoverDelta:()=>td,neutralStrokeRecipe:()=>ov,neutralStrokeRest:()=>o$,neutralStrokeRestDelta:()=>tc,numberFieldStyles:()=>aF,optionStyles:()=>aD,pickerListItemStyles:()=>iq,pickerMenuOptionStyles:()=>i_,pickerMenuStyles:()=>iE,pickerStyles:()=>iG,progressRingStyles:()=>aL,progressStyles:()=>az,provideJupyterDesignSystem:()=>i5,radioGroupStyles:()=>aA,radioStyles:()=>aO,radioTemplate:()=>aI,searchStyles:()=>aq,selectStyles:()=>rM,skeletonStyles:()=>aY,sliderLabelStyles:()=>a9,sliderStyles:()=>a1,strokeWidth:()=>eS,switchStyles:()=>ie,tabPanelStyles:()=>ir,tabStyles:()=>il,tabsStyles:()=>id,textAreaStyles:()=>ip,textFieldStyles:()=>im,toolbarStyles:()=>ik,tooltipStyles:()=>iS,treeItemStyles:()=>iO,treeViewStyles:()=>iP,typeRampBaseFontSize:()=>ez,typeRampBaseLineHeight:()=>ej,typeRampMinus1FontSize:()=>eB,typeRampMinus1LineHeight:()=>eL,typeRampMinus2FontSize:()=>eH,typeRampMinus2LineHeight:()=>eN,typeRampPlus1FontSize:()=>eO,typeRampPlus1LineHeight:()=>eI,typeRampPlus2FontSize:()=>eR,typeRampPlus2LineHeight:()=>eP,typeRampPlus3FontSize:()=>eA,typeRampPlus3LineHeight:()=>eM,typeRampPlus4FontSize:()=>eG,typeRampPlus4LineHeight:()=>eE,typeRampPlus5FontSize:()=>e_,typeRampPlus5LineHeight:()=>eq,typeRampPlus6FontSize:()=>eW,typeRampPlus6LineHeight:()=>eU,verticalSliderLabelStyles:()=>a4});class d{constructor(e,t,o,r){this.r=e,this.g=t,this.b=o,this.a="number"!=typeof r||isNaN(r)?1:r}static fromObject(e){return!e||isNaN(e.r)||isNaN(e.g)||isNaN(e.b)?null:new d(e.r,e.g,e.b,e.a)}equalValue(e){return this.r===e.r&&this.g===e.g&&this.b===e.b&&this.a===e.a}toStringHexRGB(){return"#"+[this.r,this.g,this.b].map(this.formatHexValue).join("")}toStringHexRGBA(){return this.toStringHexRGB()+this.formatHexValue(this.a)}toStringHexARGB(){return"#"+[this.a,this.r,this.g,this.b].map(this.formatHexValue).join("")}toStringWebRGB(){return`rgb(${Math.round(i(this.r,0,255))},${Math.round(i(this.g,0,255))},${Math.round(i(this.b,0,255))})`}toStringWebRGBA(){return`rgba(${Math.round(i(this.r,0,255))},${Math.round(i(this.g,0,255))},${Math.round(i(this.b,0,255))},${r(this.a,0,1)})`}roundToPrecision(e){return new d(c(this.r,e),c(this.g,e),c(this.b,e),c(this.a,e))}clamp(){return new d(r(this.r,0,1),r(this.g,0,1),r(this.b,0,1),r(this.a,0,1))}toObject(){return{r:this.r,g:this.g,b:this.b,a:this.a}}formatHexValue(e){let t;return 1===(t=Math.round(r(i(e,0,255),0,255)).toString(16)).length?"0"+t:t}}let h={aliceblue:{r:.941176,g:.972549,b:1},antiquewhite:{r:.980392,g:.921569,b:.843137},aqua:{r:0,g:1,b:1},aquamarine:{r:.498039,g:1,b:.831373},azure:{r:.941176,g:1,b:1},beige:{r:.960784,g:.960784,b:.862745},bisque:{r:1,g:.894118,b:.768627},black:{r:0,g:0,b:0},blanchedalmond:{r:1,g:.921569,b:.803922},blue:{r:0,g:0,b:1},blueviolet:{r:.541176,g:.168627,b:.886275},brown:{r:.647059,g:.164706,b:.164706},burlywood:{r:.870588,g:.721569,b:.529412},cadetblue:{r:.372549,g:.619608,b:.627451},chartreuse:{r:.498039,g:1,b:0},chocolate:{r:.823529,g:.411765,b:.117647},coral:{r:1,g:.498039,b:.313725},cornflowerblue:{r:.392157,g:.584314,b:.929412},cornsilk:{r:1,g:.972549,b:.862745},crimson:{r:.862745,g:.078431,b:.235294},cyan:{r:0,g:1,b:1},darkblue:{r:0,g:0,b:.545098},darkcyan:{r:0,g:.545098,b:.545098},darkgoldenrod:{r:.721569,g:.52549,b:.043137},darkgray:{r:.662745,g:.662745,b:.662745},darkgreen:{r:0,g:.392157,b:0},darkgrey:{r:.662745,g:.662745,b:.662745},darkkhaki:{r:.741176,g:.717647,b:.419608},darkmagenta:{r:.545098,g:0,b:.545098},darkolivegreen:{r:.333333,g:.419608,b:.184314},darkorange:{r:1,g:.54902,b:0},darkorchid:{r:.6,g:.196078,b:.8},darkred:{r:.545098,g:0,b:0},darksalmon:{r:.913725,g:.588235,b:.478431},darkseagreen:{r:.560784,g:.737255,b:.560784},darkslateblue:{r:.282353,g:.239216,b:.545098},darkslategray:{r:.184314,g:.309804,b:.309804},darkslategrey:{r:.184314,g:.309804,b:.309804},darkturquoise:{r:0,g:.807843,b:.819608},darkviolet:{r:.580392,g:0,b:.827451},deeppink:{r:1,g:.078431,b:.576471},deepskyblue:{r:0,g:.74902,b:1},dimgray:{r:.411765,g:.411765,b:.411765},dimgrey:{r:.411765,g:.411765,b:.411765},dodgerblue:{r:.117647,g:.564706,b:1},firebrick:{r:.698039,g:.133333,b:.133333},floralwhite:{r:1,g:.980392,b:.941176},forestgreen:{r:.133333,g:.545098,b:.133333},fuchsia:{r:1,g:0,b:1},gainsboro:{r:.862745,g:.862745,b:.862745},ghostwhite:{r:.972549,g:.972549,b:1},gold:{r:1,g:.843137,b:0},goldenrod:{r:.854902,g:.647059,b:.12549},gray:{r:.501961,g:.501961,b:.501961},green:{r:0,g:.501961,b:0},greenyellow:{r:.678431,g:1,b:.184314},grey:{r:.501961,g:.501961,b:.501961},honeydew:{r:.941176,g:1,b:.941176},hotpink:{r:1,g:.411765,b:.705882},indianred:{r:.803922,g:.360784,b:.360784},indigo:{r:.294118,g:0,b:.509804},ivory:{r:1,g:1,b:.941176},khaki:{r:.941176,g:.901961,b:.54902},lavender:{r:.901961,g:.901961,b:.980392},lavenderblush:{r:1,g:.941176,b:.960784},lawngreen:{r:.486275,g:.988235,b:0},lemonchiffon:{r:1,g:.980392,b:.803922},lightblue:{r:.678431,g:.847059,b:.901961},lightcoral:{r:.941176,g:.501961,b:.501961},lightcyan:{r:.878431,g:1,b:1},lightgoldenrodyellow:{r:.980392,g:.980392,b:.823529},lightgray:{r:.827451,g:.827451,b:.827451},lightgreen:{r:.564706,g:.933333,b:.564706},lightgrey:{r:.827451,g:.827451,b:.827451},lightpink:{r:1,g:.713725,b:.756863},lightsalmon:{r:1,g:.627451,b:.478431},lightseagreen:{r:.12549,g:.698039,b:.666667},lightskyblue:{r:.529412,g:.807843,b:.980392},lightslategray:{r:.466667,g:.533333,b:.6},lightslategrey:{r:.466667,g:.533333,b:.6},lightsteelblue:{r:.690196,g:.768627,b:.870588},lightyellow:{r:1,g:1,b:.878431},lime:{r:0,g:1,b:0},limegreen:{r:.196078,g:.803922,b:.196078},linen:{r:.980392,g:.941176,b:.901961},magenta:{r:1,g:0,b:1},maroon:{r:.501961,g:0,b:0},mediumaquamarine:{r:.4,g:.803922,b:.666667},mediumblue:{r:0,g:0,b:.803922},mediumorchid:{r:.729412,g:.333333,b:.827451},mediumpurple:{r:.576471,g:.439216,b:.858824},mediumseagreen:{r:.235294,g:.701961,b:.443137},mediumslateblue:{r:.482353,g:.407843,b:.933333},mediumspringgreen:{r:0,g:.980392,b:.603922},mediumturquoise:{r:.282353,g:.819608,b:.8},mediumvioletred:{r:.780392,g:.082353,b:.521569},midnightblue:{r:.098039,g:.098039,b:.439216},mintcream:{r:.960784,g:1,b:.980392},mistyrose:{r:1,g:.894118,b:.882353},moccasin:{r:1,g:.894118,b:.709804},navajowhite:{r:1,g:.870588,b:.678431},navy:{r:0,g:0,b:.501961},oldlace:{r:.992157,g:.960784,b:.901961},olive:{r:.501961,g:.501961,b:0},olivedrab:{r:.419608,g:.556863,b:.137255},orange:{r:1,g:.647059,b:0},orangered:{r:1,g:.270588,b:0},orchid:{r:.854902,g:.439216,b:.839216},palegoldenrod:{r:.933333,g:.909804,b:.666667},palegreen:{r:.596078,g:.984314,b:.596078},paleturquoise:{r:.686275,g:.933333,b:.933333},palevioletred:{r:.858824,g:.439216,b:.576471},papayawhip:{r:1,g:.937255,b:.835294},peachpuff:{r:1,g:.854902,b:.72549},peru:{r:.803922,g:.521569,b:.247059},pink:{r:1,g:.752941,b:.796078},plum:{r:.866667,g:.627451,b:.866667},powderblue:{r:.690196,g:.878431,b:.901961},purple:{r:.501961,g:0,b:.501961},red:{r:1,g:0,b:0},rosybrown:{r:.737255,g:.560784,b:.560784},royalblue:{r:.254902,g:.411765,b:.882353},saddlebrown:{r:.545098,g:.270588,b:.07451},salmon:{r:.980392,g:.501961,b:.447059},sandybrown:{r:.956863,g:.643137,b:.376471},seagreen:{r:.180392,g:.545098,b:.341176},seashell:{r:1,g:.960784,b:.933333},sienna:{r:.627451,g:.321569,b:.176471},silver:{r:.752941,g:.752941,b:.752941},skyblue:{r:.529412,g:.807843,b:.921569},slateblue:{r:.415686,g:.352941,b:.803922},slategray:{r:.439216,g:.501961,b:.564706},slategrey:{r:.439216,g:.501961,b:.564706},snow:{r:1,g:.980392,b:.980392},springgreen:{r:0,g:1,b:.498039},steelblue:{r:.27451,g:.509804,b:.705882},tan:{r:.823529,g:.705882,b:.54902},teal:{r:0,g:.501961,b:.501961},thistle:{r:.847059,g:.74902,b:.847059},tomato:{r:1,g:.388235,b:.278431},transparent:{r:0,g:0,b:0,a:0},turquoise:{r:.25098,g:.878431,b:.815686},violet:{r:.933333,g:.509804,b:.933333},wheat:{r:.960784,g:.870588,b:.701961},white:{r:1,g:1,b:1},whitesmoke:{r:.960784,g:.960784,b:.960784},yellow:{r:1,g:1,b:0},yellowgreen:{r:.603922,g:.803922,b:.196078}},u=/^rgb\(\s*((?:(?:25[0-5]|2[0-4]\d|1\d\d|\d{1,2})\s*,\s*){2}(?:25[0-5]|2[0-4]\d|1\d\d|\d{1,2})\s*)\)$/i,p=/^rgba\(\s*((?:(?:25[0-5]|2[0-4]\d|1\d\d|\d{1,2})\s*,\s*){3}(?:0|1|0?\.\d*)\s*)\)$/i,g=/^#((?:[0-9a-f]{6}|[0-9a-f]{3}))$/i,b=/^#((?:[0-9a-f]{8}|[0-9a-f]{4}))$/i;function f(e){let t=g.exec(e);if(null===t)return null;let o=t[1];if(3===o.length){let e=o.charAt(0),t=o.charAt(1),r=o.charAt(2);o=e.concat(e,t,t,r,r)}let r=parseInt(o,16);return isNaN(r)?null:new d(a((0xff0000&r)>>>16,0,255),a((65280&r)>>>8,0,255),a(255&r,0,255),1)}function m(e){let t,o=e.toLowerCase();return g.test(o)?f(o):b.test(o)?function(e){let t=b.exec(e);if(null===t)return null;let o=t[1];if(4===o.length){let e=o.charAt(0),t=o.charAt(1),r=o.charAt(2),a=o.charAt(3);o=e.concat(e,t,t,r,r,a,a)}let r=parseInt(o,16);return isNaN(r)?null:new d(a((0xff0000&r)>>>16,0,255),a((65280&r)>>>8,0,255),a(255&r,0,255),a((0xff000000&r)>>>24,0,255))}(o):u.test(o)?function(e){let t=u.exec(e);if(null===t)return null;let o=t[1].split(",");return new d(a(Number(o[0]),0,255),a(Number(o[1]),0,255),a(Number(o[2]),0,255),1)}(o):p.test(o)?function(e){let t=p.exec(e);if(null===t)return null;let o=t[1].split(",");return 4===o.length?new d(a(Number(o[0]),0,255),a(Number(o[1]),0,255),a(Number(o[2]),0,255),Number(o[3])):null}(o):h.hasOwnProperty(o)&&(t=h[o.toLowerCase()])?new d(t.r,t.g,t.b,t.hasOwnProperty("a")?t.a:void 0):null}class v{constructor(e,t,o){this.h=e,this.s=t,this.l=o}static fromObject(e){return!e||isNaN(e.h)||isNaN(e.s)||isNaN(e.l)?null:new v(e.h,e.s,e.l)}equalValue(e){return this.h===e.h&&this.s===e.s&&this.l===e.l}roundToPrecision(e){return new v(c(this.h,e),c(this.s,e),c(this.l,e))}toObject(){return{h:this.h,s:this.s,l:this.l}}}class ${constructor(e,t,o){this.h=e,this.s=t,this.v=o}static fromObject(e){return!e||isNaN(e.h)||isNaN(e.s)||isNaN(e.v)?null:new $(e.h,e.s,e.v)}equalValue(e){return this.h===e.h&&this.s===e.s&&this.v===e.v}roundToPrecision(e){return new $(c(this.h,e),c(this.s,e),c(this.v,e))}toObject(){return{h:this.h,s:this.s,v:this.v}}}class x{constructor(e,t,o){this.l=e,this.a=t,this.b=o}static fromObject(e){return!e||isNaN(e.l)||isNaN(e.a)||isNaN(e.b)?null:new x(e.l,e.a,e.b)}equalValue(e){return this.l===e.l&&this.a===e.a&&this.b===e.b}roundToPrecision(e){return new x(c(this.l,e),c(this.a,e),c(this.b,e))}toObject(){return{l:this.l,a:this.a,b:this.b}}}x.epsilon=216/24389,x.kappa=24389/27;class y{constructor(e,t,o){this.l=e,this.c=t,this.h=o}static fromObject(e){return!e||isNaN(e.l)||isNaN(e.c)||isNaN(e.h)?null:new y(e.l,e.c,e.h)}equalValue(e){return this.l===e.l&&this.c===e.c&&this.h===e.h}roundToPrecision(e){return new y(c(this.l,e),c(this.c,e),c(this.h,e))}toObject(){return{l:this.l,c:this.c,h:this.h}}}class k{constructor(e,t,o){this.x=e,this.y=t,this.z=o}static fromObject(e){return!e||isNaN(e.x)||isNaN(e.y)||isNaN(e.z)?null:new k(e.x,e.y,e.z)}equalValue(e){return this.x===e.x&&this.y===e.y&&this.z===e.z}roundToPrecision(e){return new k(c(this.x,e),c(this.y,e),c(this.z,e))}toObject(){return{x:this.x,y:this.y,z:this.z}}}function w(e){return .2126*e.r+.7152*e.g+.0722*e.b}function F(e){function t(e){return e<=.03928?e/12.92:Math.pow((e+.055)/1.055,2.4)}return w(new d(t(e.r),t(e.g),t(e.b),1))}k.whitePoint=new k(.95047,1,1.08883);let C=(e,t)=>(e+.05)/(t+.05);function V(e,t){let o=F(e),r=F(t);return o>r?C(o,r):C(r,o)}function D(e){let t=Math.max(e.r,e.g,e.b),o=Math.min(e.r,e.g,e.b),r=t-o,a=0;0!==r&&(a=t===e.r?60*((e.g-e.b)/r%6):t===e.g?60*((e.b-e.r)/r+2):60*((e.r-e.g)/r+4)),a<0&&(a+=360);let i=(t+o)/2,l=0;return 0!==r&&(l=r/(1-Math.abs(2*i-1))),new v(a,l,i)}function S(e,t=1){let o=(1-Math.abs(2*e.l-1))*e.s,r=o*(1-Math.abs(e.h/60%2-1)),a=e.l-o/2,i=0,l=0,s=0;return e.h<60?(i=o,l=r,s=0):e.h<120?(i=r,l=o,s=0):e.h<180?(i=0,l=o,s=r):e.h<240?(i=0,l=r,s=o):e.h<300?(i=r,l=0,s=o):e.h<360&&(i=o,l=0,s=r),new d(i+a,l+a,s+a,t)}function T(e){let t=Math.max(e.r,e.g,e.b),o=t-Math.min(e.r,e.g,e.b),r=0;0!==o&&(r=t===e.r?60*((e.g-e.b)/o%6):t===e.g?60*((e.b-e.r)/o+2):60*((e.r-e.g)/o+4)),r<0&&(r+=360);let a=0;return 0!==t&&(a=o/t),new $(r,a,t)}function z(e){function t(e){return e<=.04045?e/12.92:Math.pow((e+.055)/1.055,2.4)}let o=t(e.r),r=t(e.g),a=t(e.b);return new k(.4124564*o+.3575761*r+.1804375*a,.2126729*o+.7151522*r+.072175*a,.0193339*o+.119192*r+.9503041*a)}function j(e,t=1){function o(e){return e<=.0031308?12.92*e:1.055*Math.pow(e,1/2.4)-.055}return new d(o(3.2404542*e.x-1.5371385*e.y-.4985314*e.z),o(-.969266*e.x+1.8760108*e.y+.041556*e.z),o(.0556434*e.x-.2040259*e.y+1.0572252*e.z),t)}function B(e){return function(e){function t(e){return e>x.epsilon?Math.pow(e,1/3):(x.kappa*e+16)/116}let o=t(e.x/k.whitePoint.x),r=t(e.y/k.whitePoint.y);return new x(116*r-16,500*(o-r),200*(r-t(e.z/k.whitePoint.z)))}(z(e))}function L(e,t=1){let o,r,a,i,l,s,n,c,d;return j((r=(o=(e.l+16)/116)+e.a/500,a=o-e.b/200,i=Math.pow(r,3),l=Math.pow(o,3),s=Math.pow(a,3),n=0,n=i>x.epsilon?i:(116*r-16)/x.kappa,c=0,c=e.l>x.epsilon*x.kappa?l:e.l/x.kappa,d=0,d=s>x.epsilon?s:(116*a-16)/x.kappa,n=k.whitePoint.x*n,c=k.whitePoint.y*c,d=k.whitePoint.z*d,new k(n,c,d)),t)}function H(e){var t=B(e);let o=0;(Math.abs(t.b)>.001||Math.abs(t.a)>.001)&&(o=180/Math.PI*Math.atan2(t.b,t.a)),o<0&&(o+=360);let r=Math.sqrt(t.a*t.a+t.b*t.b);return new y(t.l,r,o)}function N(e,t=1){let o,r;return L((o=0,r=0,0!==e.h&&(o=Math.cos(l(e.h))*e.c,r=Math.sin(l(e.h))*e.c),new x(e.l,o,r)),t)}function O(e,t){let o=e.relativeLuminance>t.relativeLuminance?e:t,r=e.relativeLuminance>t.relativeLuminance?t:e;return(o.relativeLuminance+.05)/(r.relativeLuminance+.05)}let I=Object.freeze({create:(e,t,o)=>new R(e,t,o),from:e=>new R(e.r,e.g,e.b)});class R extends d{constructor(e,t,o){super(e,t,o,1),this.toColorString=this.toStringHexRGB,this.contrast=O.bind(null,this),this.createCSS=this.toColorString,this.relativeLuminance=F(this)}static fromObject(e){return new R(e.r,e.g,e.b)}}function P(e){return I.create(e,e,e)}let A={LightMode:1,DarkMode:.23},M=(-.1+Math.sqrt(.21))/2;function G(e){return e.relativeLuminance<=M}var E,_,q,W,U,X,Z,Y,J=o(45226),K=o(30086);function Q(e,t,o=18){let r=H(e),a=r.c+t*o;return a<0&&(a=0),N(new y(r.l,a,r.h))}function ee(e,t){return new d(e.r*t.r,e.g*t.g,e.b*t.b,1)}function et(e,t){return e<.5?r(2*t*e,0,1):r(1-2*(1-t)*(1-e),0,1)}function eo(e,t){return new d(et(e.r,t.r),et(e.g,t.g),et(e.b,t.b),1)}function er(e,t,o,r){var a,i,l,c,h,u,p,g,b,f;if(isNaN(e)||e<=0)return o;if(e>=1)return r;switch(t){case X.HSL:return S((a=D(o),i=D(r),isNaN(e)||e<=0?a:e>=1?i:new v(n(e,a.h,i.h),s(e,a.s,i.s),s(e,a.l,i.l))));case X.HSV:return function(e,t=1){let o=e.s*e.v,r=o*(1-Math.abs(e.h/60%2-1)),a=e.v-o,i=0,l=0,s=0;return e.h<60?(i=o,l=r,s=0):e.h<120?(i=r,l=o,s=0):e.h<180?(i=0,l=o,s=r):e.h<240?(i=0,l=r,s=o):e.h<300?(i=r,l=0,s=o):e.h<360&&(i=o,l=0,s=r),new d(i+a,l+a,s+a,t)}((l=T(o),c=T(r),isNaN(e)||e<=0?l:e>=1?c:new $(n(e,l.h,c.h),s(e,l.s,c.s),s(e,l.v,c.v))));case X.XYZ:return j((h=z(o),u=z(r),isNaN(e)||e<=0?h:e>=1?u:new k(s(e,h.x,u.x),s(e,h.y,u.y),s(e,h.z,u.z))));case X.LAB:return L((p=B(o),g=B(r),isNaN(e)||e<=0?p:e>=1?g:new x(s(e,p.l,g.l),s(e,p.a,g.a),s(e,p.b,g.b))));case X.LCH:return N((b=H(o),f=H(r),isNaN(e)||e<=0?b:e>=1?f:new y(s(e,b.l,f.l),s(e,b.c,f.c),n(e,b.h,f.h))));default:return isNaN(e)||e<=0?o:e>=1?r:new d(s(e,o.r,r.r),s(e,o.g,r.g),s(e,o.b,r.b),s(e,o.a,r.a))}}(E=U||(U={}))[E.Burn=0]="Burn",E[E.Color=1]="Color",E[E.Darken=2]="Darken",E[E.Dodge=3]="Dodge",E[E.Lighten=4]="Lighten",E[E.Multiply=5]="Multiply",E[E.Overlay=6]="Overlay",E[E.Screen=7]="Screen",(_=X||(X={}))[_.RGB=0]="RGB",_[_.HSL=1]="HSL",_[_.HSV=2]="HSV",_[_.XYZ=3]="XYZ",_[_.LAB=4]="LAB",_[_.LCH=5]="LCH";class ea{constructor(e){if(null==e||0===e.length)throw Error("The stops argument must be non-empty");this.stops=this.sortColorScaleStops(e)}static createBalancedColorScale(e){if(null==e||0===e.length)throw Error("The colors argument must be non-empty");let t=Array(e.length);for(let o=0;o<e.length;o++)0===o?t[o]={color:e[o],position:0}:o===e.length-1?t[o]={color:e[o],position:1}:t[o]={color:e[o],position:o*(1/(e.length-1))};return new ea(t)}getColor(e,t=X.RGB){if(1===this.stops.length||e<=0)return this.stops[0].color;if(e>=1)return this.stops[this.stops.length-1].color;let o=0;for(let t=0;t<this.stops.length;t++)this.stops[t].position<=e&&(o=t);let r=o+1;return r>=this.stops.length&&(r=this.stops.length-1),er((e-this.stops[o].position)*(1/(this.stops[r].position-this.stops[o].position)),t,this.stops[o].color,this.stops[r].color)}trim(e,t,o=X.RGB){if(e<0||t>1||t<e)throw Error("Invalid bounds");if(e===t)return new ea([{color:this.getColor(e,o),position:0}]);let r=[];for(let o=0;o<this.stops.length;o++)this.stops[o].position>=e&&this.stops[o].position<=t&&r.push(this.stops[o]);if(0===r.length)return new ea([{color:this.getColor(e),position:e},{color:this.getColor(t),position:t}]);r[0].position!==e&&r.unshift({color:this.getColor(e),position:e}),r[r.length-1].position!==t&&r.push({color:this.getColor(t),position:t});let a=t-e,i=Array(r.length);for(let t=0;t<r.length;t++)i[t]={color:r[t].color,position:(r[t].position-e)/a};return new ea(i)}findNextColor(e,t,o=!1,r=X.RGB,a=.005,i=32){isNaN(e)||e<=0?e=0:e>=1&&(e=1);let l=this.getColor(e,r),s=+!o;if(V(l,this.getColor(s,r))<=t)return s;let n=o?0:e,c=o?e:0,d=s,h=0;for(;h<=i;){d=Math.abs(c-n)/2+n;let e=V(l,this.getColor(d,r));if(Math.abs(e-t)<=a)break;e>t?o?n=d:c=d:o?c=d:n=d,h++}return d}clone(){let e=Array(this.stops.length);for(let t=0;t<e.length;t++)e[t]={color:this.stops[t].color,position:this.stops[t].position};return new ea(e)}sortColorScaleStops(e){return e.sort((e,t)=>{let o=e.position,r=t.position;return o<r?-1:+(o>r)})}}class ei{constructor(e){this.config=Object.assign({},ei.defaultPaletteConfig,e),this.palette=[],this.updatePaletteColors()}updatePaletteGenerationValues(e){let t=!1;for(let o in e)this.config[o]&&(this.config[o].equalValue?this.config[o].equalValue(e[o])||(this.config[o]=e[o],t=!0):e[o]!==this.config[o]&&(this.config[o]=e[o],t=!0));return t&&this.updatePaletteColors(),t}updatePaletteColors(){let e=this.generatePaletteColorScale();for(let t=0;t<this.config.steps;t++)this.palette[t]=e.getColor(t/(this.config.steps-1),this.config.interpolationMode)}generatePaletteColorScale(){let e=D(this.config.baseColor),t=new ea([{position:0,color:this.config.scaleColorLight},{position:.5,color:this.config.baseColor},{position:1,color:this.config.scaleColorDark}]).trim(this.config.clipLight,1-this.config.clipDark),o=t.getColor(0),r=t.getColor(1),a=o,i=r;if(e.s>=this.config.saturationAdjustmentCutoff&&(a=Q(a,this.config.saturationLight),i=Q(i,this.config.saturationDark)),0!==this.config.multiplyLight){let e=ee(this.config.baseColor,a);a=er(this.config.multiplyLight,this.config.interpolationMode,a,e)}if(0!==this.config.multiplyDark){let e=ee(this.config.baseColor,i);i=er(this.config.multiplyDark,this.config.interpolationMode,i,e)}if(0!==this.config.overlayLight){let e=eo(this.config.baseColor,a);a=er(this.config.overlayLight,this.config.interpolationMode,a,e)}if(0!==this.config.overlayDark){let e=eo(this.config.baseColor,i);i=er(this.config.overlayDark,this.config.interpolationMode,i,e)}return this.config.baseScalePosition?new ea(this.config.baseScalePosition<=0?[{position:0,color:this.config.baseColor},{position:1,color:i.clamp()}]:this.config.baseScalePosition>=1?[{position:0,color:a.clamp()},{position:1,color:this.config.baseColor}]:[{position:0,color:a.clamp()},{position:this.config.baseScalePosition,color:this.config.baseColor},{position:1,color:i.clamp()}]):new ea([{position:0,color:a.clamp()},{position:.5,color:this.config.baseColor},{position:1,color:i.clamp()}])}}ei.defaultPaletteConfig={baseColor:f("#808080"),steps:11,interpolationMode:X.RGB,scaleColorLight:new d(1,1,1,1),scaleColorDark:new d(0,0,0,1),clipLight:.185,clipDark:.16,saturationAdjustmentCutoff:.05,saturationLight:.35,saturationDark:1.25,overlayLight:0,overlayDark:.25,multiplyLight:0,multiplyDark:0,baseScalePosition:.5},ei.greyscalePaletteConfig={baseColor:f("#808080"),steps:11,interpolationMode:X.RGB,scaleColorLight:new d(1,1,1,1),scaleColorDark:new d(0,0,0,1),clipLight:0,clipDark:0,saturationAdjustmentCutoff:0,saturationLight:0,saturationDark:0,overlayLight:0,overlayDark:0,multiplyLight:0,multiplyDark:0,baseScalePosition:.5},ei.defaultPaletteConfig.scaleColorLight,ei.defaultPaletteConfig.scaleColorDark;class el{constructor(e){this.palette=[],this.config=Object.assign({},el.defaultPaletteConfig,e),this.regenPalettes()}regenPalettes(){let e=this.config.steps;(isNaN(e)||e<3)&&(e=3);let t=new d(.14,.14,.14,1),o=new ei(Object.assign(Object.assign({},ei.greyscalePaletteConfig),{baseColor:t,baseScalePosition:.9148936170212766,steps:e})).palette,r=w(this.config.baseColor),a=D(this.config.baseColor).l,i=this.matchRelativeLuminanceIndex((r+a)/2,o)/(e-1),l=this.matchRelativeLuminanceIndex(.14,o)/(e-1),s=D(this.config.baseColor),n=S(v.fromObject({h:s.h,s:s.s,l:.14})),c=S(v.fromObject({h:s.h,s:s.s,l:.06})),h=[,,,,,];h[0]={position:0,color:new d(1,1,1,1)},h[1]={position:i,color:this.config.baseColor},h[2]={position:l,color:n},h[3]={position:.99,color:c},h[4]={position:1,color:new d(0,0,0,1)};let u=new ea(h);this.palette=Array(e);for(let t=0;t<e;t++){let o=u.getColor(t/(e-1),X.RGB);this.palette[t]=o}}matchRelativeLuminanceIndex(e,t){let o=Number.MAX_VALUE,r=0,a=0,i=t.length;for(;a<i;a++){let i=Math.abs(w(t[a])-e);i<o&&(o=i,r=a)}return r}}function es(e){return G(e)?-1:1}el.defaultPaletteConfig={baseColor:f("#808080"),steps:94};let en=Object.freeze({create:function(e,t,o){return"number"==typeof e?en.from(I.create(e,t,o)):en.from(e)},from:function(e){return!function(e){let t={r:0,g:0,b:0,toColorString:()=>"",contrast:()=>0,relativeLuminance:0};for(let o in t)if(typeof t[o]!=typeof e[o])return!1;return!0}(e)?ec.from(I.create(e.r,e.g,e.b)):ec.from(e)}});class ec{constructor(e,t){this.closestIndexCache=new Map,this.source=e,this.swatches=t,this.reversedSwatches=Object.freeze([...this.swatches].reverse()),this.lastIndex=this.swatches.length-1}colorContrast(e,t,o,r){void 0===o&&(o=this.closestIndexOf(e));let a=this.swatches,i=this.lastIndex,l=o;return void 0===r&&(r=es(e)),-1===r&&(a=this.reversedSwatches,l=i-l),function e(t,o,r=0,a=t.length-1){if(a===r)return t[r];let i=Math.floor((a-r)/2)+r;return o(t[i])?e(t,o,r,i):e(t,o,i+1,a)}(a,o=>O(e,o)>=t,l,i)}get(e){return this.swatches[e]||this.swatches[r(e,0,this.lastIndex)]}closestIndexOf(e){if(this.closestIndexCache.has(e.relativeLuminance))return this.closestIndexCache.get(e.relativeLuminance);let t=this.swatches.indexOf(e);if(-1!==t)return this.closestIndexCache.set(e.relativeLuminance,t),t;let o=this.swatches.reduce((t,o)=>Math.abs(o.relativeLuminance-e.relativeLuminance)<Math.abs(t.relativeLuminance-e.relativeLuminance)?o:t);return t=this.swatches.indexOf(o),this.closestIndexCache.set(e.relativeLuminance,t),t}static from(e){return new ec(e,Object.freeze(new el({baseColor:d.fromObject(e)}).palette.map(e=>{let t=f(e.toStringHexRGB());return I.create(t.r,t.g,t.b)})))}}let ed=I.create(1,1,1),eh=I.create(0,0,0),eu=I.from(f("#808080")),ep=I.from(f("#DA1A5F")),eg=I.from(f("#D32F2F"));function eb(e,t,o,r,a,i){return Math.max(e.closestIndexOf(P(t))+o,r,a,i)}let{create:ef}=J.DesignToken;function em(e){return J.DesignToken.create({name:e,cssCustomPropertyName:null})}let ev=ef("body-font").withDefault('aktiv-grotesk, "Segoe UI", Arial, Helvetica, sans-serif'),e$=ef("base-height-multiplier").withDefault(10),ex=ef("base-horizontal-spacing-multiplier").withDefault(3),ey=ef("base-layer-luminance").withDefault(A.DarkMode),ek=ef("control-corner-radius").withDefault(4),ew=ef("density").withDefault(0),eF=ef("design-unit").withDefault(4),eC=ef("element-scale").withDefault(0),eV=ef("direction").withDefault(K.O.ltr),eD=ef("disabled-opacity").withDefault(.4),eS=ef("stroke-width").withDefault(1),eT=ef("focus-stroke-width").withDefault(2),ez=ef("type-ramp-base-font-size").withDefault("14px"),ej=ef("type-ramp-base-line-height").withDefault("20px"),eB=ef("type-ramp-minus-1-font-size").withDefault("12px"),eL=ef("type-ramp-minus-1-line-height").withDefault("16px"),eH=ef("type-ramp-minus-2-font-size").withDefault("10px"),eN=ef("type-ramp-minus-2-line-height").withDefault("16px"),eO=ef("type-ramp-plus-1-font-size").withDefault("16px"),eI=ef("type-ramp-plus-1-line-height").withDefault("24px"),eR=ef("type-ramp-plus-2-font-size").withDefault("20px"),eP=ef("type-ramp-plus-2-line-height").withDefault("28px"),eA=ef("type-ramp-plus-3-font-size").withDefault("28px"),eM=ef("type-ramp-plus-3-line-height").withDefault("36px"),eG=ef("type-ramp-plus-4-font-size").withDefault("34px"),eE=ef("type-ramp-plus-4-line-height").withDefault("44px"),e_=ef("type-ramp-plus-5-font-size").withDefault("46px"),eq=ef("type-ramp-plus-5-line-height").withDefault("56px"),eW=ef("type-ramp-plus-6-font-size").withDefault("60px"),eU=ef("type-ramp-plus-6-line-height").withDefault("72px"),eX=em("accent-fill-rest-delta").withDefault(0),eZ=em("accent-fill-hover-delta").withDefault(4),eY=em("accent-fill-active-delta").withDefault(-5),eJ=em("accent-fill-focus-delta").withDefault(0),eK=em("accent-foreground-rest-delta").withDefault(0),eQ=em("accent-foreground-hover-delta").withDefault(6),e0=em("accent-foreground-active-delta").withDefault(-4),e1=em("accent-foreground-focus-delta").withDefault(0),e2=em("neutral-fill-rest-delta").withDefault(7),e5=em("neutral-fill-hover-delta").withDefault(10),e6=em("neutral-fill-active-delta").withDefault(5),e3=em("neutral-fill-focus-delta").withDefault(0),e4=em("neutral-fill-input-rest-delta").withDefault(0),e9=em("neutral-fill-input-hover-delta").withDefault(0),e8=em("neutral-fill-input-active-delta").withDefault(0),e7=em("neutral-fill-input-focus-delta").withDefault(0),te=em("neutral-fill-stealth-rest-delta").withDefault(0),tt=em("neutral-fill-stealth-hover-delta").withDefault(5),to=em("neutral-fill-stealth-active-delta").withDefault(3),tr=em("neutral-fill-stealth-focus-delta").withDefault(0),ta=em("neutral-fill-strong-rest-delta").withDefault(0),ti=em("neutral-fill-strong-hover-delta").withDefault(8),tl=em("neutral-fill-strong-active-delta").withDefault(-5),ts=em("neutral-fill-strong-focus-delta").withDefault(0),tn=em("neutral-fill-layer-rest-delta").withDefault(3),tc=em("neutral-stroke-rest-delta").withDefault(25),td=em("neutral-stroke-hover-delta").withDefault(40),th=em("neutral-stroke-active-delta").withDefault(16),tu=em("neutral-stroke-focus-delta").withDefault(25),tp=em("neutral-stroke-divider-rest-delta").withDefault(8),tg=ef("neutral-color").withDefault(eu),tb=em("neutral-palette").withDefault(e=>en.from(tg.getValueFor(e))),tf=ef("accent-color").withDefault(ep),tm=em("accent-palette").withDefault(e=>en.from(tf.getValueFor(e))),tv=em("neutral-layer-card-container-recipe").withDefault({evaluate:e=>{var t,o,r;return t=tb.getValueFor(e),o=ey.getValueFor(e),r=tn.getValueFor(e),t.get(t.closestIndexOf(P(o))+r)}}),t$=ef("neutral-layer-card-container").withDefault(e=>tv.getValueFor(e).evaluate(e)),tx=em("neutral-layer-floating-recipe").withDefault({evaluate:e=>{var t,o,r;let a;return t=tb.getValueFor(e),o=ey.getValueFor(e),r=tn.getValueFor(e),a=t.closestIndexOf(P(o))-r,t.get(a-r)}}),ty=ef("neutral-layer-floating").withDefault(e=>tx.getValueFor(e).evaluate(e)),tk=em("neutral-layer-1-recipe").withDefault({evaluate:e=>{var t,o;return t=tb.getValueFor(e),o=ey.getValueFor(e),t.get(t.closestIndexOf(P(o)))}}),tw=ef("neutral-layer-1").withDefault(e=>tk.getValueFor(e).evaluate(e)),tF=em("neutral-layer-2-recipe").withDefault({evaluate:e=>{var t,o,r,a,i,l;return t=tb.getValueFor(e),o=ey.getValueFor(e),r=tn.getValueFor(e),a=e2.getValueFor(e),i=e5.getValueFor(e),l=e6.getValueFor(e),t.get(eb(t,o,r,a,i,l))}}),tC=ef("neutral-layer-2").withDefault(e=>tF.getValueFor(e).evaluate(e)),tV=em("neutral-layer-3-recipe").withDefault({evaluate:e=>{var t,o,r,a,i,l;return t=tb.getValueFor(e),o=ey.getValueFor(e),r=tn.getValueFor(e),a=e2.getValueFor(e),i=e5.getValueFor(e),l=e6.getValueFor(e),t.get(eb(t,o,r,a,i,l)+r)}}),tD=ef("neutral-layer-3").withDefault(e=>tV.getValueFor(e).evaluate(e)),tS=em("neutral-layer-4-recipe").withDefault({evaluate:e=>{var t,o,r,a,i,l;return t=tb.getValueFor(e),o=ey.getValueFor(e),r=tn.getValueFor(e),a=e2.getValueFor(e),i=e5.getValueFor(e),l=e6.getValueFor(e),t.get(eb(t,o,r,a,i,l)+2*r)}}),tT=ef("neutral-layer-4").withDefault(e=>tS.getValueFor(e).evaluate(e)),tz=ef("fill-color").withDefault(e=>tw.getValueFor(e));(q=Z||(Z={}))[q.normal=4.5]="normal",q[q.large=7]="large";let tj=ef({name:"accent-fill-recipe",cssCustomPropertyName:null}).withDefault({evaluate:(e,t)=>{var o,r,a,i,l,s,n,c,d;let h,u,p,g;return o=tm.getValueFor(e),r=tb.getValueFor(e),a=t||tz.getValueFor(e),i=eZ.getValueFor(e),l=eY.getValueFor(e),s=eJ.getValueFor(e),n=e2.getValueFor(e),c=e5.getValueFor(e),d=e6.getValueFor(e),h=o.source,u=r.closestIndexOf(a)>=Math.max(n,c,d)?-1:1,g=(p=o.closestIndexOf(h))+-1*u*i,{rest:o.get(g),hover:o.get(p),active:o.get(g+u*l),focus:o.get(g+u*s)}}}),tB=ef("accent-fill-rest").withDefault(e=>tj.getValueFor(e).evaluate(e).rest),tL=ef("accent-fill-hover").withDefault(e=>tj.getValueFor(e).evaluate(e).hover),tH=ef("accent-fill-active").withDefault(e=>tj.getValueFor(e).evaluate(e).active),tN=ef("accent-fill-focus").withDefault(e=>tj.getValueFor(e).evaluate(e).focus),tO=e=>(t,o)=>{var r;return r=o||tB.getValueFor(t),r.contrast(ed)>=e?ed:eh},tI=em("foreground-on-accent-recipe").withDefault({evaluate:(e,t)=>tO(Z.normal)(e,t)}),tR=ef("foreground-on-accent-rest").withDefault(e=>tI.getValueFor(e).evaluate(e,tB.getValueFor(e))),tP=ef("foreground-on-accent-hover").withDefault(e=>tI.getValueFor(e).evaluate(e,tL.getValueFor(e))),tA=ef("foreground-on-accent-active").withDefault(e=>tI.getValueFor(e).evaluate(e,tH.getValueFor(e))),tM=ef("foreground-on-accent-focus").withDefault(e=>tI.getValueFor(e).evaluate(e,tN.getValueFor(e))),tG=em("foreground-on-accent-large-recipe").withDefault({evaluate:(e,t)=>tO(Z.large)(e,t)}),tE=ef("foreground-on-accent-rest-large").withDefault(e=>tG.getValueFor(e).evaluate(e,tB.getValueFor(e))),t_=ef("foreground-on-accent-hover-large").withDefault(e=>tG.getValueFor(e).evaluate(e,tL.getValueFor(e))),tq=ef("foreground-on-accent-active-large").withDefault(e=>tG.getValueFor(e).evaluate(e,tH.getValueFor(e))),tW=ef("foreground-on-accent-focus-large").withDefault(e=>tG.getValueFor(e).evaluate(e,tN.getValueFor(e))),tU=ef({name:"accent-foreground-recipe",cssCustomPropertyName:null}).withDefault({evaluate:(e,t)=>{let o;return(o=Z.normal,(e,t)=>{var r,a,i,l,s,n;let c,d,h,u,p,g,b,f,m;return r=tm.getValueFor(e),a=t||tz.getValueFor(e),i=eK.getValueFor(e),l=eQ.getValueFor(e),s=e0.getValueFor(e),n=e1.getValueFor(e),h=r.source,u=r.closestIndexOf(h),g=u+(1===(p=es(a))?Math.min(i,l):Math.max(p*i,p*l)),b=r.colorContrast(a,o,g,p),m=(f=r.closestIndexOf(b))+p*Math.abs(i-l),(1===p?i<l:p*i>p*l)?(c=f,d=m):(c=m,d=f),{rest:r.get(c),hover:r.get(d),active:r.get(c+p*s),focus:r.get(c+p*n)}})(e,t)}}),tX=ef("accent-foreground-rest").withDefault(e=>tU.getValueFor(e).evaluate(e).rest),tZ=ef("accent-foreground-hover").withDefault(e=>tU.getValueFor(e).evaluate(e).hover),tY=ef("accent-foreground-active").withDefault(e=>tU.getValueFor(e).evaluate(e).active),tJ=ef("accent-foreground-focus").withDefault(e=>tU.getValueFor(e).evaluate(e).focus),tK=ef({name:"neutral-fill-recipe",cssCustomPropertyName:null}).withDefault({evaluate:(e,t)=>{var o,r,a,i,l,s;let n,c;return o=tb.getValueFor(e),r=t||tz.getValueFor(e),a=e2.getValueFor(e),i=e5.getValueFor(e),l=e6.getValueFor(e),s=e3.getValueFor(e),c=(n=o.closestIndexOf(r))>=Math.max(a,i,l,s)?-1:1,{rest:o.get(n+c*a),hover:o.get(n+c*i),active:o.get(n+c*l),focus:o.get(n+c*s)}}}),tQ=ef("neutral-fill-rest").withDefault(e=>tK.getValueFor(e).evaluate(e).rest),t0=ef("neutral-fill-hover").withDefault(e=>tK.getValueFor(e).evaluate(e).hover),t1=ef("neutral-fill-active").withDefault(e=>tK.getValueFor(e).evaluate(e).active),t2=ef("neutral-fill-focus").withDefault(e=>tK.getValueFor(e).evaluate(e).focus),t5=ef({name:"neutral-fill-input-recipe",cssCustomPropertyName:null}).withDefault({evaluate:(e,t)=>{var o,r,a,i,l,s;let n,c;return o=tb.getValueFor(e),r=t||tz.getValueFor(e),a=e4.getValueFor(e),i=e9.getValueFor(e),l=e8.getValueFor(e),s=e7.getValueFor(e),n=es(r),c=o.closestIndexOf(r),{rest:o.get(c-n*a),hover:o.get(c-n*i),active:o.get(c-n*l),focus:o.get(c-n*s)}}}),t6=ef("neutral-fill-input-rest").withDefault(e=>t5.getValueFor(e).evaluate(e).rest),t3=ef("neutral-fill-input-hover").withDefault(e=>t5.getValueFor(e).evaluate(e).hover),t4=ef("neutral-fill-input-active").withDefault(e=>t5.getValueFor(e).evaluate(e).active),t9=ef("neutral-fill-input-focus").withDefault(e=>t5.getValueFor(e).evaluate(e).focus),t8=ef({name:"neutral-fill-stealth-recipe",cssCustomPropertyName:null}).withDefault({evaluate:(e,t)=>{var o,r,a,i,l,s,n,c;let d,h,u;return o=tb.getValueFor(e),r=t||tz.getValueFor(e),a=te.getValueFor(e),i=tt.getValueFor(e),l=to.getValueFor(e),s=tr.getValueFor(e),n=e2.getValueFor(e),c=e5.getValueFor(e),d=Math.max(a,i,l,s,n,c,e6.getValueFor(e),e3.getValueFor(e)),u=(h=o.closestIndexOf(r))>=d?-1:1,{rest:o.get(h+u*a),hover:o.get(h+u*i),active:o.get(h+u*l),focus:o.get(h+u*s)}}}),t7=ef("neutral-fill-stealth-rest").withDefault(e=>t8.getValueFor(e).evaluate(e).rest),oe=ef("neutral-fill-stealth-hover").withDefault(e=>t8.getValueFor(e).evaluate(e).hover),ot=ef("neutral-fill-stealth-active").withDefault(e=>t8.getValueFor(e).evaluate(e).active),oo=ef("neutral-fill-stealth-focus").withDefault(e=>t8.getValueFor(e).evaluate(e).focus),or=ef({name:"neutral-fill-strong-recipe",cssCustomPropertyName:null}).withDefault({evaluate:(e,t)=>{var o,r,a,i,l,s;let n,c,d,h,u;return o=tb.getValueFor(e),r=t||tz.getValueFor(e),a=ta.getValueFor(e),i=ti.getValueFor(e),l=tl.getValueFor(e),s=ts.getValueFor(e),d=es(r),u=(h=o.closestIndexOf(o.colorContrast(r,4.5)))+d*Math.abs(a-i),(1===d?a<i:d*a>d*i)?(n=h,c=u):(n=u,c=h),{rest:o.get(n),hover:o.get(c),active:o.get(n+d*l),focus:o.get(n+d*s)}}}),oa=ef("neutral-fill-strong-rest").withDefault(e=>or.getValueFor(e).evaluate(e).rest),oi=ef("neutral-fill-strong-hover").withDefault(e=>or.getValueFor(e).evaluate(e).hover),ol=ef("neutral-fill-strong-active").withDefault(e=>or.getValueFor(e).evaluate(e).active),os=ef("neutral-fill-strong-focus").withDefault(e=>or.getValueFor(e).evaluate(e).focus),on=em("neutral-fill-layer-recipe").withDefault({evaluate:(e,t)=>{var o,r,a;let i;return o=tb.getValueFor(e),r=t||tz.getValueFor(e),a=tn.getValueFor(e),i=o.closestIndexOf(r),o.get(i-(i<a?-1*a:a))}}),oc=ef("neutral-fill-layer-rest").withDefault(e=>on.getValueFor(e).evaluate(e)),od=em("focus-stroke-outer-recipe").withDefault({evaluate:e=>{var t,o;return t=tb.getValueFor(e),o=tz.getValueFor(e),t.colorContrast(o,3.5)}}),oh=ef("focus-stroke-outer").withDefault(e=>od.getValueFor(e).evaluate(e)),ou=em("focus-stroke-inner-recipe").withDefault({evaluate:e=>{var t,o,r;return t=tm.getValueFor(e),o=tz.getValueFor(e),r=oh.getValueFor(e),t.colorContrast(r,3.5,t.closestIndexOf(t.source),-1*es(o))}}),op=ef("focus-stroke-inner").withDefault(e=>ou.getValueFor(e).evaluate(e)),og=em("neutral-foreground-hint-recipe").withDefault({evaluate:e=>{var t,o;return t=tb.getValueFor(e),o=tz.getValueFor(e),t.colorContrast(o,4.5)}}),ob=ef("neutral-foreground-hint").withDefault(e=>og.getValueFor(e).evaluate(e)),of=em("neutral-foreground-recipe").withDefault({evaluate:e=>{var t,o;return t=tb.getValueFor(e),o=tz.getValueFor(e),t.colorContrast(o,14)}}),om=ef("neutral-foreground-rest").withDefault(e=>of.getValueFor(e).evaluate(e)),ov=ef({name:"neutral-stroke-recipe",cssCustomPropertyName:null}).withDefault({evaluate:e=>{var t,o,r,a,i,l;let s,n,c;return t=tb.getValueFor(e),o=tz.getValueFor(e),r=tc.getValueFor(e),a=td.getValueFor(e),i=th.getValueFor(e),l=tu.getValueFor(e),s=t.closestIndexOf(o),c=s+(n=es(o))*r,{rest:t.get(c),hover:t.get(c+n*(a-r)),active:t.get(c+n*(i-r)),focus:t.get(c+n*(l-r))}}}),o$=ef("neutral-stroke-rest").withDefault(e=>ov.getValueFor(e).evaluate(e).rest),ox=ef("neutral-stroke-hover").withDefault(e=>ov.getValueFor(e).evaluate(e).hover),oy=ef("neutral-stroke-active").withDefault(e=>ov.getValueFor(e).evaluate(e).active),ok=ef("neutral-stroke-focus").withDefault(e=>ov.getValueFor(e).evaluate(e).focus),ow=em("neutral-stroke-divider-recipe").withDefault({evaluate:(e,t)=>{var o,r,a;return o=tb.getValueFor(e),r=t||tz.getValueFor(e),a=tp.getValueFor(e),o.get(o.closestIndexOf(r)+es(r)*a)}}),oF=ef("neutral-stroke-divider-rest").withDefault(e=>ow.getValueFor(e).evaluate(e)),oC=J.DesignToken.create({name:"height-number",cssCustomPropertyName:null}).withDefault(e=>(e$.getValueFor(e)+ew.getValueFor(e))*eF.getValueFor(e)),oV=ef("error-color").withDefault(eg),oD=em("error-palette").withDefault(e=>en.from(oV.getValueFor(e))),oS=ef({name:"error-fill-recipe",cssCustomPropertyName:null}).withDefault({evaluate:(e,t)=>{var o,r,a,i,l,s,n,c,d;let h,u,p,g;return o=oD.getValueFor(e),r=tb.getValueFor(e),a=t||tz.getValueFor(e),i=eZ.getValueFor(e),l=eY.getValueFor(e),s=eJ.getValueFor(e),n=e2.getValueFor(e),c=e5.getValueFor(e),d=e6.getValueFor(e),h=o.source,u=r.closestIndexOf(a)>=Math.max(n,c,d)?-1:1,g=(p=o.closestIndexOf(h))+-1*u*i,{rest:o.get(g),hover:o.get(p),active:o.get(g+u*l),focus:o.get(g+u*s)}}}),oT=ef("error-fill-rest").withDefault(e=>oS.getValueFor(e).evaluate(e).rest),oz=ef("error-fill-hover").withDefault(e=>oS.getValueFor(e).evaluate(e).hover),oj=ef("error-fill-active").withDefault(e=>oS.getValueFor(e).evaluate(e).active),oB=ef("error-fill-focus").withDefault(e=>oS.getValueFor(e).evaluate(e).focus),oL=e=>(t,o)=>{var r;return r=o||oT.getValueFor(t),r.contrast(ed)>=e?ed:eh},oH=ef({name:"foreground-on-error-recipe",cssCustomPropertyName:null}).withDefault({evaluate:(e,t)=>oL(Z.normal)(e,t)}),oN=ef("foreground-on-error-rest").withDefault(e=>oH.getValueFor(e).evaluate(e,oT.getValueFor(e))),oO=ef("foreground-on-error-hover").withDefault(e=>oH.getValueFor(e).evaluate(e,oz.getValueFor(e))),oI=ef("foreground-on-error-active").withDefault(e=>oH.getValueFor(e).evaluate(e,oj.getValueFor(e))),oR=ef("foreground-on-error-focus").withDefault(e=>oH.getValueFor(e).evaluate(e,oB.getValueFor(e))),oP=ef({name:"foreground-on-error-large-recipe",cssCustomPropertyName:null}).withDefault({evaluate:(e,t)=>oL(Z.large)(e,t)}),oA=ef("foreground-on-error-rest-large").withDefault(e=>oP.getValueFor(e).evaluate(e,oT.getValueFor(e))),oM=ef("foreground-on-error-hover-large").withDefault(e=>oP.getValueFor(e).evaluate(e,oz.getValueFor(e))),oG=ef("foreground-on-error-active-large").withDefault(e=>oP.getValueFor(e).evaluate(e,oj.getValueFor(e))),oE=ef("foreground-on-error-focus-large").withDefault(e=>oP.getValueFor(e).evaluate(e,oB.getValueFor(e))),o_=ef({name:"error-foreground-recipe",cssCustomPropertyName:null}).withDefault({evaluate:(e,t)=>{let o;return(o=Z.normal,(e,t)=>{var r,a,i,l,s,n;let c,d,h,u,p,g,b,f,m;return r=oD.getValueFor(e),a=t||tz.getValueFor(e),i=eK.getValueFor(e),l=eQ.getValueFor(e),s=e0.getValueFor(e),n=e1.getValueFor(e),h=r.source,u=r.closestIndexOf(h),g=u+(1==(p=G(a)?-1:1)?Math.min(i,l):Math.max(p*i,p*l)),b=r.colorContrast(a,o,g,p),m=(f=r.closestIndexOf(b))+p*Math.abs(i-l),(1===p?i<l:p*i>p*l)?(c=f,d=m):(c=m,d=f),{rest:r.get(c),hover:r.get(d),active:r.get(c+p*s),focus:r.get(c+p*n)}})(e,t)}}),oq=ef("error-foreground-rest").withDefault(e=>o_.getValueFor(e).evaluate(e).rest),oW=ef("error-foreground-hover").withDefault(e=>o_.getValueFor(e).evaluate(e).hover),oU=ef("error-foreground-active").withDefault(e=>o_.getValueFor(e).evaluate(e).active),oX=ef("error-foreground-focus").withDefault(e=>o_.getValueFor(e).evaluate(e).focus),oZ="--jp-layout-color1",oY=!1;function oJ(){let e;oY||(oY=!0,e=()=>{new MutationObserver(()=>{o0()}).observe(document.body,{attributes:!0,attributeFilter:["data-jp-theme-name"],childList:!1,characterData:!1}),o0()},"complete"===document.readyState?e():window.addEventListener("load",e))}let oK=e=>{let t=parseInt(e,10);return isNaN(t)?null:t},oQ={"--jp-border-width":{converter:oK,token:eS},"--jp-border-radius":{converter:oK,token:ek},[oZ]:{converter:(e,t)=>{let o=m(e);if(!o)return null;{let e=D(o),t=S(v.fromObject({h:e.h,s:e.s,l:.5}));return I.create(t.r,t.g,t.b)}},token:tg},"--jp-brand-color1":{converter:(e,t)=>{let o=m(e);if(!o)return null;{let e=D(o),r=S(v.fromObject({h:e.h,s:e.s,l:e.l+(t?1:-1)*eZ.getValueFor(document.body)/94}));return I.create(r.r,r.g,r.b)}},token:tf},"--jp-error-color1":{converter:(e,t)=>{let o=m(e);if(!o)return null;{let e=D(o),r=S(v.fromObject({h:e.h,s:e.s,l:e.l+(t?1:-1)*eZ.getValueFor(document.body)/94}));return I.create(r.r,r.g,r.b)}},token:oV},"--jp-ui-font-family":{token:ev},"--jp-ui-font-size1":{token:ez}};function o0(){var e;let t=getComputedStyle(document.body),o=document.body.getAttribute("data-jp-theme-light"),r=!1;if(o)r="false"===o;else{let e=t.getPropertyValue(oZ).toString();if(e){let t=m(e);t&&(r=G(I.create(t.r,t.g,t.b)),console.debug(`Theme is ${r?"dark":"light"} based on '${oZ}' value: ${e}.`))}}for(let o in ey.setValueFor(document.body,r?A.DarkMode:A.LightMode),oQ){let a=oQ[o],i=t.getPropertyValue(o).toString();if(document.body&&""!==i){let t=(null!=(e=a.converter)?e:e=>e)(i.trim(),r);null!==t?a.token.setValueFor(document.body,t):console.error(`Fail to parse value '${i}' for '${o}' as FAST design token.`)}}}var o1=o(23453);let o2=(e,t)=>(0,o1.css)`
  ${(0,J.display)("flex")} :host {
    box-sizing: border-box;
    flex-direction: column;
    font-family: ${ev};
    font-size: ${eB};
    line-height: ${eL};
    color: ${om};
    border-top: calc(${eS} * 1px) solid ${oF};
  }
`;(W=Y||(Y={})).Canvas="Canvas",W.CanvasText="CanvasText",W.LinkText="LinkText",W.VisitedText="VisitedText",W.ActiveText="ActiveText",W.ButtonFace="ButtonFace",W.ButtonText="ButtonText",W.Field="Field",W.FieldText="FieldText",W.Highlight="Highlight",W.HighlightText="HighlightText",W.GrayText="GrayText";let o5=(0,o1.cssPartial)`(${e$} + ${ew} + ${eC}) * ${eF}`,o6=(e,t)=>(0,o1.css)`
    ${(0,J.display)("flex")} :host {
      box-sizing: border-box;
      font-family: ${ev};
      flex-direction: column;
      font-size: ${eB};
      line-height: ${eL};
      border-bottom: calc(${eS} * 1px) solid
        ${oF};
    }

    .region {
      display: none;
      padding: calc((6 + (${eF} * 2 * ${ew})) * 1px);
    }

    div.heading {
      display: grid;
      position: relative;
      grid-template-columns: calc(${o5} * 1px) auto 1fr auto;
      color: ${om};
    }

    .button {
      appearance: none;
      border: none;
      background: none;
      grid-column: 3;
      outline: none;
      padding: 0 calc((6 + (${eF} * 2 * ${ew})) * 1px);
      text-align: left;
      height: calc(${o5} * 1px);
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
    .button:${J.focusVisible}::before {
      outline: none;
      border: calc(${eT} * 1px) solid ${tN};
      border-radius: calc(${ek} * 1px);
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
      padding-inline-start: calc(${eF} * 1px);
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
  `.withBehaviors((0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
      /* prettier-ignore */
      .button:${J.focusVisible}::before {
          border-color: ${Y.Highlight};
        }
      :host slot[name='collapsed-icon'],
      :host([expanded]) slot[name='expanded-icon'] {
        fill: ${Y.ButtonText};
      }
    `));class o3 extends J.AccordionItem{}let o4=o3.compose({baseName:"accordion-item",baseClass:J.AccordionItem,template:J.accordionItemTemplate,styles:o6,collapsedIcon:`
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
    `});class o9 extends J.Accordion{}let o8=o9.compose({baseName:"accordion",baseClass:J.Accordion,template:J.accordionTemplate,styles:o2});function o7(e,t,o,r){var a,i=arguments.length,l=i<3?t:null===r?r=Object.getOwnPropertyDescriptor(t,o):r;if("object"==typeof Reflect&&"function"==typeof Reflect.decorate)l=Reflect.decorate(e,t,o,r);else for(var s=e.length-1;s>=0;s--)(a=e[s])&&(l=(i<3?a(l):i>3?a(t,o,l):a(t,o))||l);return i>3&&l&&Object.defineProperty(t,o,l),l}let re=(0,o1.css)`
  ${(0,J.display)("inline-flex")} :host {
    font-family: ${ev};
    outline: none;
    font-size: ${ez};
    line-height: ${ej};
    height: calc(${o5} * 1px);
    min-width: calc(${o5} * 1px);
    background-color: ${tQ};
    color: ${om};
    border-radius: calc(${ek} * 1px);
    fill: currentcolor;
    cursor: pointer;
    margin: calc((${eT} + 2) * 1px);
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
        calc((10 + (${eF} * 2 * (${ew} + ${eC})))) * 1px
      );
    white-space: nowrap;
    outline: none;
    text-decoration: none;
    border: calc(${eS} * 1px) solid transparent;
    color: inherit;
    border-radius: inherit;
    fill: inherit;
    cursor: inherit;
    font-family: inherit;
    font-size: inherit;
    line-height: inherit;
  }

  :host(:hover) {
    background-color: ${t0};
  }

  :host(:active) {
    background-color: ${t1};
  }

  :host([aria-pressed='true']) {
    box-shadow: inset 0px 0px 2px 2px ${ol};
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
  .control:${J.focusVisible} {
      outline: calc(${eT} * 1px) solid ${os};
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
`.withBehaviors((0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
    :host .control {
      background-color: ${Y.ButtonFace};
      border-color: ${Y.ButtonText};
      color: ${Y.ButtonText};
      fill: currentColor;
    }

    :host(:hover) .control {
      forced-color-adjust: none;
      background-color: ${Y.Highlight};
      color: ${Y.HighlightText};
    }

    /* prettier-ignore */
    .control:${J.focusVisible} {
          forced-color-adjust: none;
          background-color: ${Y.Highlight};
          outline-color: ${Y.ButtonText};
          color: ${Y.HighlightText};
        }

    .control:hover,
    :host([appearance='outline']) .control:hover {
      border-color: ${Y.ButtonText};
    }

    :host([href]) .control {
      border-color: ${Y.LinkText};
      color: ${Y.LinkText};
    }

    :host([href]) .control:hover,
        :host([href]) .control:${J.focusVisible} {
      forced-color-adjust: none;
      background: ${Y.ButtonFace};
      outline-color: ${Y.LinkText};
      color: ${Y.LinkText};
      fill: currentColor;
    }
  `)),rt=(0,o1.css)`
  :host([appearance='accent']) {
    background: ${tB};
    color: ${tR};
  }

  :host([appearance='accent']:hover) {
    background: ${tL};
    color: ${tP};
  }

  :host([appearance='accent'][aria-pressed='true']) {
    box-shadow: inset 0px 0px 2px 2px ${tY};
  }

  :host([appearance='accent']:active) .control:active {
    background: ${tH};
    color: ${tA};
  }

  :host([appearance="accent"]) .control:${J.focusVisible} {
    outline-color: ${tN};
  }
`.withBehaviors((0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
    :host([appearance='accent']) .control {
      forced-color-adjust: none;
      background: ${Y.Highlight};
      color: ${Y.HighlightText};
    }

    :host([appearance='accent']) .control:hover,
    :host([appearance='accent']:active) .control:active {
      background: ${Y.HighlightText};
      border-color: ${Y.Highlight};
      color: ${Y.Highlight};
    }

    :host([appearance="accent"]) .control:${J.focusVisible} {
      outline-color: ${Y.Highlight};
    }

    :host([appearance='accent'][href]) .control {
      background: ${Y.LinkText};
      color: ${Y.HighlightText};
    }

    :host([appearance='accent'][href]) .control:hover {
      background: ${Y.ButtonFace};
      border-color: ${Y.LinkText};
      box-shadow: none;
      color: ${Y.LinkText};
      fill: currentColor;
    }

    :host([appearance="accent"][href]) .control:${J.focusVisible} {
      outline-color: ${Y.HighlightText};
    }
  `)),ro=(0,o1.css)`
  :host([appearance='error']) {
    background: ${oT};
    color: ${tR};
  }

  :host([appearance='error']:hover) {
    background: ${oz};
    color: ${tP};
  }

  :host([appearance='error'][aria-pressed='true']) {
    box-shadow: inset 0px 0px 2px 2px ${oU};
  }

  :host([appearance='error']:active) .control:active {
    background: ${oj};
    color: ${tA};
  }

  :host([appearance="error"]) .control:${J.focusVisible} {
    outline-color: ${oB};
  }
`.withBehaviors((0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
    :host([appearance='error']) .control {
      forced-color-adjust: none;
      background: ${Y.Highlight};
      color: ${Y.HighlightText};
    }

    :host([appearance='error']) .control:hover,
    :host([appearance='error']:active) .control:active {
      background: ${Y.HighlightText};
      border-color: ${Y.Highlight};
      color: ${Y.Highlight};
    }

    :host([appearance="error"]) .control:${J.focusVisible} {
      outline-color: ${Y.Highlight};
    }

    :host([appearance='error'][href]) .control {
      background: ${Y.LinkText};
      color: ${Y.HighlightText};
    }

    :host([appearance='error'][href]) .control:hover {
      background: ${Y.ButtonFace};
      border-color: ${Y.LinkText};
      box-shadow: none;
      color: ${Y.LinkText};
      fill: currentColor;
    }

    :host([appearance="error"][href]) .control:${J.focusVisible} {
      outline-color: ${Y.HighlightText};
    }
  `)),rr=(0,o1.css)`
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
    color: ${tX};
    border-bottom: calc(${eS} * 1px) solid ${tX};
  }

  :host([appearance='hypertext']:hover),
  :host([appearance='hypertext']) .control:hover {
    background: transparent;
    border-bottom-color: ${tZ};
  }

  :host([appearance='hypertext']:active),
  :host([appearance='hypertext']) .control:active {
    background: transparent;
    border-bottom-color: ${tY};
  }

  :host([appearance="hypertext"]) .control:${J.focusVisible} {
    outline-color: transparent;
    border-bottom: calc(${eT} * 1px) solid ${oh};
    margin-bottom: calc(calc(${eS} - ${eT}) * 1px);
  }
`.withBehaviors((0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
    :host([appearance='hypertext']:hover) {
      background-color: ${Y.ButtonFace};
      color: ${Y.ButtonText};
    }
    :host([appearance="hypertext"][href]) .control:hover,
        :host([appearance="hypertext"][href]) .control:active,
        :host([appearance="hypertext"][href]) .control:${J.focusVisible} {
      color: ${Y.LinkText};
      border-bottom-color: ${Y.LinkText};
      box-shadow: none;
    }
  `)),ra=(0,o1.css)`
  :host([appearance='lightweight']) {
    background: transparent;
    color: ${tX};
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
    color: ${tZ};
  }

  :host([appearance='lightweight']:active) {
    background: transparent;
    color: ${tY};
  }

  :host([appearance='lightweight']) .content {
    position: relative;
  }

  :host([appearance='lightweight']) .content::before {
    content: '';
    display: block;
    height: calc(${eS} * 1px);
    position: absolute;
    top: calc(1em + 4px);
    width: 100%;
  }

  :host([appearance='lightweight']:hover) .content::before {
    background: ${tZ};
  }

  :host([appearance='lightweight']:active) .content::before {
    background: ${tY};
  }

  :host([appearance="lightweight"]) .control:${J.focusVisible} {
    outline-color: transparent;
  }

  :host([appearance="lightweight"]) .control:${J.focusVisible} .content::before {
    background: ${om};
    height: calc(${eT} * 1px);
  }
`.withBehaviors((0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
    :host([appearance="lightweight"]) .control:hover,
        :host([appearance="lightweight"]) .control:${J.focusVisible} {
      forced-color-adjust: none;
      background: ${Y.ButtonFace};
      color: ${Y.Highlight};
    }
    :host([appearance="lightweight"]) .control:hover .content::before,
        :host([appearance="lightweight"]) .control:${J.focusVisible} .content::before {
      background: ${Y.Highlight};
    }

    :host([appearance="lightweight"][href]) .control:hover,
        :host([appearance="lightweight"][href]) .control:${J.focusVisible} {
      background: ${Y.ButtonFace};
      box-shadow: none;
      color: ${Y.LinkText};
    }

    :host([appearance="lightweight"][href]) .control:hover .content::before,
        :host([appearance="lightweight"][href]) .control:${J.focusVisible} .content::before {
      background: ${Y.LinkText};
    }
  `)),ri=(0,o1.css)`
  :host([appearance='outline']) {
    background: transparent;
    border-color: ${tB};
  }

  :host([appearance='outline']:hover) {
    border-color: ${tL};
  }

  :host([appearance='outline']:active) {
    border-color: ${tH};
  }

  :host([appearance='outline']) .control {
    border-color: inherit;
  }

  :host([appearance="outline"]) .control:${J.focusVisible} {
    outline-color: ${tN};
  }
`.withBehaviors((0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
    :host([appearance='outline']) .control {
      border-color: ${Y.ButtonText};
    }
    :host([appearance="outline"]) .control:${J.focusVisible} {
      forced-color-adjust: none;
      background-color: ${Y.Highlight};
      outline-color: ${Y.ButtonText};
      color: ${Y.HighlightText};
      fill: currentColor;
    }
    :host([appearance='outline'][href]) .control {
      background: ${Y.ButtonFace};
      border-color: ${Y.LinkText};
      color: ${Y.LinkText};
      fill: currentColor;
    }
    :host([appearance="outline"][href]) .control:hover,
        :host([appearance="outline"][href]) .control:${J.focusVisible} {
      forced-color-adjust: none;
      outline-color: ${Y.LinkText};
    }
  `)),rl=(0,o1.css)`
  :host([appearance='stealth']),
  :host([appearance='stealth'][disabled]:active),
  :host([appearance='stealth'][disabled]:hover) {
    background: transparent;
  }

  :host([appearance='stealth']:hover) {
    background: ${oe};
  }

  :host([appearance='stealth']:active) {
    background: ${ot};
  }

  :host([appearance='stealth']) .control:${J.focusVisible} {
    outline-color: ${tN};
  }

  /* Make the focus outline displayed within the button if
     it is in a start or end slot; e.g. in a tree item
     This will make the focus outline bounded within the container.
   */
  :host([appearance='stealth'][slot="end"]) .control:${J.focusVisible},
  :host([appearance='stealth'][slot="start"]) .control:${J.focusVisible} {
    outline-offset: -2px;
  }
`.withBehaviors((0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
    :host([appearance='stealth']),
    :host([appearance='stealth']) .control {
      forced-color-adjust: none;
      background: ${Y.ButtonFace};
      border-color: transparent;
      color: ${Y.ButtonText};
      fill: currentColor;
    }

    :host([appearance='stealth']:hover) .control {
      background: ${Y.Highlight};
      border-color: ${Y.Highlight};
      color: ${Y.HighlightText};
      fill: currentColor;
    }

    :host([appearance="stealth"]:${J.focusVisible}) .control {
      outline-color: ${Y.Highlight};
      color: ${Y.HighlightText};
      fill: currentColor;
    }

    :host([appearance='stealth'][href]) .control {
      color: ${Y.LinkText};
    }

    :host([appearance="stealth"][href]:hover) .control,
        :host([appearance="stealth"][href]:${J.focusVisible}) .control {
      background: ${Y.LinkText};
      border-color: ${Y.LinkText};
      color: ${Y.HighlightText};
      fill: currentColor;
    }

    :host([appearance="stealth"][href]:${J.focusVisible}) .control {
      forced-color-adjust: none;
      box-shadow: 0 0 0 1px ${Y.LinkText};
    }
  `));function rs(e,t){return new J.PropertyStyleSheetBehavior("appearance",e,t)}let rn=(e,t)=>(0,o1.css)`
    ${re}
  `.withBehaviors(rs("accent",rt),rs("hypertext",rr),rs("lightweight",ra),rs("outline",ri),rs("stealth",rl));class rc extends J.Anchor{appearanceChanged(e,t){this.$fastController.isConnected&&(this.classList.remove(e),this.classList.add(t))}connectedCallback(){super.connectedCallback(),this.appearance||(this.appearance="neutral")}defaultSlottedContentChanged(e,t){let o=this.defaultSlottedContent.filter(e=>e.nodeType===Node.ELEMENT_NODE);1===o.length&&o[0]instanceof SVGElement?this.control.classList.add("icon-only"):this.control.classList.remove("icon-only")}}o7([o1.attr],rc.prototype,"appearance",void 0);let rd=rc.compose({baseName:"anchor",baseClass:J.Anchor,template:J.anchorTemplate,styles:rn,shadowOptions:{delegatesFocus:!0}}),rh=(e,t)=>(0,o1.css)`
  :host {
    contain: layout;
    display: block;
  }
`;class ru extends J.AnchoredRegion{}let rp=ru.compose({baseName:"anchored-region",baseClass:J.AnchoredRegion,template:J.anchoredRegionTemplate,styles:rh});class rg{constructor(e,t){this.cache=new WeakMap,this.ltr=e,this.rtl=t}bind(e){this.attach(e)}unbind(e){let t=this.cache.get(e);t&&eV.unsubscribe(t)}attach(e){let t=this.cache.get(e)||new rb(this.ltr,this.rtl,e),o=eV.getValueFor(e);eV.subscribe(t),t.attach(o),this.cache.set(e,t)}}class rb{constructor(e,t,o){this.ltr=e,this.rtl=t,this.source=o,this.attached=null}handleChange({target:e,token:t}){this.attach(t.getValueFor(e))}attach(e){this.attached!==this[e]&&(null!==this.attached&&this.source.$fastController.removeStyles(this.attached),this.attached=this[e],null!==this.attached&&this.source.$fastController.addStyles(this.attached))}}let rf=(e,t)=>(0,o1.css)`
    ${(0,J.display)("flex")} :host {
      position: relative;
      height: var(--avatar-size, var(--avatar-size-default));
      width: var(--avatar-size, var(--avatar-size-default));
      --avatar-size-default: calc(
        (
            (${e$} + ${ew}) * ${eF} +
              ((${eF} * 8) - 40)
          ) * 1px
      );
      --avatar-text-size: ${ez};
      --avatar-text-ratio: ${eF};
    }

    .link {
      text-decoration: none;
      color: ${om};
      display: flex;
      flex-direction: row;
      justify-content: center;
      align-items: center;
      min-width: 100%;
    }

    .square {
      border-radius: calc(${ek} * 1px);
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
      background-color: ${tB};
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

    ::slotted(${e.tagFor(J.Badge)}) {
      position: absolute;
      display: block;
    }
  `.withBehaviors(new rg((0,o1.css)`
  ::slotted(${e.tagFor(J.Badge)}) {
    right: 0;
  }
`,(0,o1.css)`
  ::slotted(${e.tagFor(J.Badge)}) {
    left: 0;
  }
`));class rm extends J.Avatar{}o7([(0,o1.attr)({attribute:"src"})],rm.prototype,"imgSrc",void 0),o7([o1.attr],rm.prototype,"alt",void 0);let rv=(0,o1.html)`
  ${(0,o1.when)(e=>e.imgSrc,(0,o1.html)`
      <img
        src="${e=>e.imgSrc}"
        alt="${e=>e.alt}"
        slot="media"
        class="media"
        part="media"
      />
    `)}
`,r$=rm.compose({baseName:"avatar",baseClass:J.Avatar,template:J.avatarTemplate,styles:rf,media:rv,shadowOptions:{delegatesFocus:!0}}),rx=(e,t)=>(0,o1.css)`
  ${(0,J.display)("inline-block")} :host {
    box-sizing: border-box;
    font-family: ${ev};
    font-size: ${eB};
    line-height: ${eL};
  }

  .control {
    border-radius: calc(${ek} * 1px);
    padding: calc(((${eF} * 0.5) - ${eS}) * 1px)
      calc((${eF} - ${eS}) * 1px);
    color: ${om};
    font-weight: 600;
    border: calc(${eS} * 1px) solid transparent;
    background-color: ${tQ};
  }

  .control[style] {
    font-weight: 400;
  }

  :host([circular]) .control {
    border-radius: 100px;
    padding: 0 calc(${eF} * 1px);
    height: calc((${o5} - (${eF} * 3)) * 1px);
    min-width: calc((${o5} - (${eF} * 3)) * 1px);
    display: flex;
    align-items: center;
    justify-content: center;
    box-sizing: border-box;
  }
`;class ry extends J.Badge{}let rk=ry.compose({baseName:"badge",baseClass:J.Badge,template:J.badgeTemplate,styles:rx}),rw=(e,t)=>(0,o1.css)`
  ${(0,J.display)("inline-block")} :host {
    box-sizing: border-box;
    font-family: ${ev};
    font-size: ${ez};
    line-height: ${ej};
  }

  .list {
    display: flex;
    flex-wrap: wrap;
  }
`;class rF extends J.Breadcrumb{}let rC=rF.compose({baseName:"breadcrumb",baseClass:J.Breadcrumb,template:J.breadcrumbTemplate,styles:rw}),rV=(e,t)=>(0,o1.css)`
    ${(0,J.display)("inline-flex")} :host {
        background: transparent;
        box-sizing: border-box;
        font-family: ${ev};
        font-size: ${ez};
        fill: currentColor;
        line-height: ${ej};
        min-width: calc(${o5} * 1px);
        outline: none;
        color: ${om}
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
        color: ${tX};
        cursor: pointer;
        display: flex;
        fill: inherit;
        outline: none;
        text-decoration: none;
        white-space: nowrap;
    }

    .control:hover {
        color: ${tZ};
    }

    .control:active {
        color: ${tY};
    }

    .control .content {
        position: relative;
    }

    .control .content::before {
        content: "";
        display: block;
        height: calc(${eS} * 1px);
        left: 0;
        position: absolute;
        right: 0;
        top: calc(1em + 4px);
        width: 100%;
    }

    .control:hover .content::before {
        background: ${tZ};
    }

    .control:active .content::before {
        background: ${tY};
    }

    .control:${J.focusVisible} .content::before {
        background: ${tJ};
        height: calc(${eT} * 1px);
    }

    .control:not([href]) {
        color: ${om};
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
`.withBehaviors((0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
      .control:hover .content::before,
                .control:${J.focusVisible} .content::before {
        background: ${Y.LinkText};
      }
      .start,
      .end {
        fill: ${Y.ButtonText};
      }
    `));class rD extends J.BreadcrumbItem{}let rS=rD.compose({baseName:"breadcrumb-item",baseClass:J.BreadcrumbItem,template:J.breadcrumbItemTemplate,styles:rV,separator:"/",shadowOptions:{delegatesFocus:!0}}),rT=(e,t)=>(0,o1.css)`
    :host([disabled]),
    :host([disabled]:hover),
    :host([disabled]:active) {
      opacity: ${eD};
      background-color: ${tQ};
      cursor: ${J.disabledCursor};
    }

    ${re}
  `.withBehaviors((0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
      :host([disabled]),
      :host([disabled]) .control,
      :host([disabled]:hover),
      :host([disabled]:active) {
        forced-color-adjust: none;
        background-color: ${Y.ButtonFace};
        outline-color: ${Y.GrayText};
        color: ${Y.GrayText};
        cursor: ${J.disabledCursor};
        opacity: 1;
      }
    `),rs("accent",(0,o1.css)`
        :host([appearance='accent'][disabled]),
        :host([appearance='accent'][disabled]:hover),
        :host([appearance='accent'][disabled]:active) {
          background: ${tB};
        }

        ${rt}
      `.withBehaviors((0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
          :host([appearance='accent'][disabled]) .control,
          :host([appearance='accent'][disabled]) .control:hover {
            background: ${Y.ButtonFace};
            border-color: ${Y.GrayText};
            color: ${Y.GrayText};
          }
        `))),rs("error",(0,o1.css)`
        :host([appearance='error'][disabled]),
        :host([appearance='error'][disabled]:hover),
        :host([appearance='error'][disabled]:active) {
          background: ${oT};
        }

        ${ro}
      `.withBehaviors((0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
          :host([appearance='error'][disabled]) .control,
          :host([appearance='error'][disabled]) .control:hover {
            background: ${Y.ButtonFace};
            border-color: ${Y.GrayText};
            color: ${Y.GrayText};
          }
        `))),rs("lightweight",(0,o1.css)`
        :host([appearance='lightweight'][disabled]:hover),
        :host([appearance='lightweight'][disabled]:active) {
          background-color: transparent;
          color: ${tX};
        }

        :host([appearance='lightweight'][disabled]) .content::before,
        :host([appearance='lightweight'][disabled]:hover) .content::before,
        :host([appearance='lightweight'][disabled]:active) .content::before {
          background: transparent;
        }

        ${ra}
      `.withBehaviors((0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
          :host([appearance='lightweight'].disabled) .control {
            forced-color-adjust: none;
            color: ${Y.GrayText};
          }

          :host([appearance='lightweight'].disabled)
            .control:hover
            .content::before {
            background: none;
          }
        `))),rs("outline",(0,o1.css)`
        :host([appearance='outline'][disabled]),
        :host([appearance='outline'][disabled]:hover),
        :host([appearance='outline'][disabled]:active) {
          background: transparent;
          border-color: ${tB};
        }

        ${ri}
      `.withBehaviors((0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
          :host([appearance='outline'][disabled]) .control {
            border-color: ${Y.GrayText};
          }
        `))),rs("stealth",(0,o1.css)`
        ${rl}
      `.withBehaviors((0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
          :host([appearance='stealth'][disabled]) {
            background: ${Y.ButtonFace};
          }

          :host([appearance='stealth'][disabled]) .control {
            background: ${Y.ButtonFace};
            border-color: transparent;
            color: ${Y.GrayText};
          }
        `))));class rz extends J.Button{constructor(){super(...arguments),this.appearance="neutral"}defaultSlottedContentChanged(e,t){let o=this.defaultSlottedContent.filter(e=>e.nodeType===Node.ELEMENT_NODE);1===o.length&&(o[0]instanceof SVGElement||o[0].classList.contains("fa")||o[0].classList.contains("fas"))?this.control.classList.add("icon-only"):this.control.classList.remove("icon-only")}}o7([o1.attr],rz.prototype,"appearance",void 0),o7([(0,o1.attr)({attribute:"minimal",mode:"boolean"})],rz.prototype,"minimal",void 0),o7([o1.attr],rz.prototype,"scale",void 0);let rj=rz.compose({baseName:"button",baseClass:J.Button,template:J.buttonTemplate,styles:rT,shadowOptions:{delegatesFocus:!0}}),rB="box-shadow: 0 0 calc((var(--elevation) * 0.225px) + 2px) rgba(0, 0, 0, calc(.11 * (2 - var(--background-luminance, 1)))), 0 calc(var(--elevation) * 0.4px) calc((var(--elevation) * 0.9px)) rgba(0, 0, 0, calc(.13 * (2 - var(--background-luminance, 1))));",rL=(e,t)=>(0,o1.css)`
    ${(0,J.display)("block")} :host {
      --elevation: 4;
      display: block;
      contain: content;
      height: var(--card-height, 100%);
      width: var(--card-width, 100%);
      box-sizing: border-box;
      background: ${tz};
      border-radius: calc(${ek} * 1px);
      ${rB}
    }
  `.withBehaviors((0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
      :host {
        forced-color-adjust: none;
        background: ${Y.Canvas};
        box-shadow: 0 0 0 1px ${Y.CanvasText};
      }
    `));class rH extends J.Card{connectedCallback(){super.connectedCallback();let e=(0,J.composedParent)(this);e&&tz.setValueFor(this,t=>on.getValueFor(t).evaluate(t,tz.getValueFor(e)))}}let rN=rH.compose({baseName:"card",baseClass:J.Card,template:J.cardTemplate,styles:rL}),rO=(e,t)=>(0,o1.css)`
    ${(0,J.display)("inline-flex")} :host {
      align-items: center;
      outline: none;
      margin: calc(${eF} * 1px) 0;
      /* Chromium likes to select label text or the default slot when the checkbox is
            clicked. Maybe there is a better solution here? */
      user-select: none;
    }

    .control {
      position: relative;
      width: calc((${o5} / 2 + ${eF}) * 1px);
      height: calc((${o5} / 2 + ${eF}) * 1px);
      box-sizing: border-box;
      border-radius: calc(${ek} * 1px);
      border: calc(${eS} * 1px) solid ${o$};
      background: ${t6};
      outline: none;
      cursor: pointer;
    }

    :host([aria-invalid='true']) .control {
      border-color: ${oT};
    }

    .label {
      font-family: ${ev};
      color: ${om};
      /* Need to discuss with Brian how HorizontalSpacingNumber can work.
            https://github.com/microsoft/fast/issues/2766 */
      padding-inline-start: calc(${eF} * 2px + 2px);
      margin-inline-end: calc(${eF} * 2px + 2px);
      cursor: pointer;
      font-size: ${ez};
      line-height: ${ej};
    }

    .label__hidden {
      display: none;
      visibility: hidden;
    }

    .checked-indicator {
      width: 100%;
      height: 100%;
      display: block;
      fill: ${tR};
      opacity: 0;
      pointer-events: none;
    }

    .indeterminate-indicator {
      border-radius: calc(${ek} * 1px);
      background: ${tR};
      position: absolute;
      top: 50%;
      left: 50%;
      width: 50%;
      height: 50%;
      transform: translate(-50%, -50%);
      opacity: 0;
    }

    :host(:not([disabled])) .control:hover {
      background: ${t3};
      border-color: ${ox};
    }

    :host(:not([disabled])) .control:active {
      background: ${t4};
      border-color: ${oy};
    }

    :host([aria-invalid='true']:not([disabled])) .control:hover {
      border-color: ${oz};
    }

    :host([aria-invalid='true']:not([disabled])) .control:active {
      border-color: ${oj};
    }

    :host(:${J.focusVisible}) .control {
      outline: calc(${eT} * 1px) solid ${tN};
      outline-offset: 2px;
    }

    :host([aria-invalid='true']:${J.focusVisible}) .control {
      outline-color: ${oB};
    }

    :host([aria-checked='true']) .control {
      background: ${tB};
      border: calc(${eS} * 1px) solid ${tB};
    }

    :host([aria-checked='true']:not([disabled])) .control:hover {
      background: ${tL};
      border: calc(${eS} * 1px) solid ${tL};
    }

    :host([aria-invalid='true'][aria-checked='true']) .control {
      background-color: ${oT};
      border-color: ${oT};
    }

    :host([aria-invalid='true'][aria-checked='true']:not([disabled]))
      .control:hover {
      background-color: ${oz};
      border-color: ${oz};
    }

    :host([aria-checked='true']:not([disabled]))
      .control:hover
      .checked-indicator {
      fill: ${tP};
    }

    :host([aria-checked='true']:not([disabled]))
      .control:hover
      .indeterminate-indicator {
      background: ${tP};
    }

    :host([aria-checked='true']:not([disabled])) .control:active {
      background: ${tH};
      border: calc(${eS} * 1px) solid ${tH};
    }

    :host([aria-invalid='true'][aria-checked='true']:not([disabled]))
      .control:active {
      background-color: ${oj};
      border-color: ${oj};
    }

    :host([aria-checked='true']:not([disabled]))
      .control:active
      .checked-indicator {
      fill: ${tA};
    }

    :host([aria-checked='true']:not([disabled]))
      .control:active
      .indeterminate-indicator {
      background: ${tA};
    }

    :host([aria-checked="true"]:${J.focusVisible}:not([disabled])) .control {
      outline: calc(${eT} * 1px) solid ${tN};
      outline-offset: 2px;
    }

    :host([aria-invalid='true'][aria-checked="true"]:${J.focusVisible}:not([disabled])) .control {
      outline-color: ${oB};
    }

    :host([disabled]) .label,
    :host([readonly]) .label,
    :host([readonly]) .control,
    :host([disabled]) .control {
      cursor: ${J.disabledCursor};
    }

    :host([aria-checked='true']:not(.indeterminate)) .checked-indicator,
    :host(.indeterminate) .indeterminate-indicator {
      opacity: 1;
    }

    :host([disabled]) {
      opacity: ${eD};
    }
  `.withBehaviors((0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
      .control {
        forced-color-adjust: none;
        border-color: ${Y.FieldText};
        background: ${Y.Field};
      }
      :host([aria-invalid='true']) .control {
        border-style: dashed;
      }
      .checked-indicator {
        fill: ${Y.FieldText};
      }
      .indeterminate-indicator {
        background: ${Y.FieldText};
      }
      :host(:not([disabled])) .control:hover,
      .control:active {
        border-color: ${Y.Highlight};
        background: ${Y.Field};
      }
      :host(:${J.focusVisible}) .control {
        outline: calc(${eT} * 1px) solid ${Y.FieldText};
        outline-offset: 2px;
      }
      :host([aria-checked="true"]:${J.focusVisible}:not([disabled])) .control {
        outline: calc(${eT} * 1px) solid ${Y.FieldText};
        outline-offset: 2px;
      }
      :host([aria-checked='true']) .control {
        background: ${Y.Highlight};
        border-color: ${Y.Highlight};
      }
      :host([aria-checked='true']:not([disabled])) .control:hover,
      .control:active {
        border-color: ${Y.Highlight};
        background: ${Y.HighlightText};
      }
      :host([aria-checked='true']) .checked-indicator {
        fill: ${Y.HighlightText};
      }
      :host([aria-checked='true']:not([disabled]))
        .control:hover
        .checked-indicator {
        fill: ${Y.Highlight};
      }
      :host([aria-checked='true']) .indeterminate-indicator {
        background: ${Y.HighlightText};
      }
      :host([aria-checked='true']) .control:hover .indeterminate-indicator {
        background: ${Y.Highlight};
      }
      :host([disabled]) {
        opacity: 1;
      }
      :host([disabled]) .control {
        forced-color-adjust: none;
        border-color: ${Y.GrayText};
        background: ${Y.Field};
      }
      :host([disabled]) .indeterminate-indicator,
      :host([aria-checked='true'][disabled])
        .control:hover
        .indeterminate-indicator {
        forced-color-adjust: none;
        background: ${Y.GrayText};
      }
      :host([disabled]) .checked-indicator,
      :host([aria-checked='true'][disabled]) .control:hover .checked-indicator {
        forced-color-adjust: none;
        fill: ${Y.GrayText};
      }
    `)),rI=(e,t)=>(0,o1.html)`
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
      <slot ${(0,o1.slotted)("defaultSlottedNodes")}></slot>
    </label>
  </template>
`;class rR extends J.Checkbox{indeterminateChanged(e,t){this.indeterminate?this.classList.add("indeterminate"):this.classList.remove("indeterminate")}}let rP=rR.compose({baseName:"checkbox",baseClass:J.Checkbox,template:rI,styles:rO,checkedIndicator:`
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
    `}),rA=(e,t)=>{let o=e.tagFor(J.ListboxOption),r=e.name===e.tagFor(J.ListboxElement)?"":".listbox";return(0,o1.css)`
        ${!r?(0,J.display)("inline-flex"):""}

        :host ${r} {
            background: ${tz};
            border: calc(${eS} * 1px) solid ${o$};
            border-radius: calc(${ek} * 1px);
            box-sizing: border-box;
            flex-direction: column;
            padding: calc(${eF} * 1px) 0;
        }

        ${!r?(0,o1.css)`
:host(:${J.focusVisible}:not([disabled])) {
                outline: none;
            }

            :host(:focus-within:not([disabled])) {
                border-color: ${oh};
                box-shadow: 0 0 0
                    calc((${eT} - ${eS}) * 1px)
                    ${oh} inset;
            }

            :host([disabled]) ::slotted(*) {
                cursor: ${J.disabledCursor};
                opacity: ${eD};
                pointer-events: none;
            }
        `:""}

        ${r||":host([size])"} {
            max-height: calc(
                (var(--size) * ${o5} + (${eF} * ${eS} * 2)) * 1px
            );
            overflow-y: auto;
        }

        :host([size="0"]) ${r} {
            max-height: none;
        }
    `.withBehaviors((0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
                :host(:not([multiple]):${J.focusVisible}) ::slotted(${o}[aria-selected="true"]),
                :host([multiple]:${J.focusVisible}) ::slotted(${o}[aria-checked="true"]) {
                    border-color: ${Y.ButtonText};
                    box-shadow: 0 0 0 calc(${eT} * 1px) inset ${Y.HighlightText};
                }

                :host(:not([multiple]):${J.focusVisible}) ::slotted(${o}[aria-selected="true"]) {
                    background: ${Y.Highlight};
                    color: ${Y.HighlightText};
                    fill: currentcolor;
                }

                ::slotted(${o}[aria-selected="true"]:not([aria-checked="true"])) {
                    background: ${Y.Highlight};
                    border-color: ${Y.HighlightText};
                    color: ${Y.HighlightText};
                }
            `))},rM=(e,t)=>{let o=e.name===e.tagFor(J.Select);return(0,o1.css)`
  ${(0,J.display)("inline-flex")}
  
  :host {
    --elevation: 14;
    background: ${t6};
    border-radius: calc(${ek} * 1px);
    border: calc(${eS} * 1px) solid ${oa};
    box-sizing: border-box;
    color: ${om};
    font-family: ${ev};
    height: calc(${o5} * 1px);
    position: relative;
    user-select: none;
    min-width: 250px;
    outline: none;
    vertical-align: top;
  }

  :host([aria-invalid='true']) {
    border-color: ${oT};
  }
  
  :host(:not([autowidth])) {
    min-width: 250px;
  }
  
  ${o?(0,o1.css)`
  :host(:not([aria-haspopup])) {
    --elevation: 0;
    border: 0;
    height: auto;
    min-width: 0;
  }
  `:""}
  
  ${rA(e,t)}
  
  :host .listbox {
    ${rB}
    border: none;
    display: flex;
    left: 0;
    position: absolute;
    width: 100%;
    z-index: 1;
  }
  
  .control + .listbox {
    --stroke-size: calc(${eF} * ${eS} * 2);
    max-height: calc(
      (var(--listbox-max-height) * ${o5} + var(--stroke-size)) * 1px
      );
  }
  
  ${o?(0,o1.css)`
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
    padding: 0 calc(1em + ${eF} * 1.25px + 1px);
  }
  
  .listbox[hidden] {
    display: none;
  }
  
  .control {
    align-items: center;
    box-sizing: border-box;
    cursor: pointer;
    display: flex;
    font-size: ${ez};
    font-family: inherit;
    line-height: ${ej};
    min-height: 100%;
    padding: 0 calc(${eF} * 2.25px);
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
    background: ${t3};
    border-color: ${oi};
  }
  
  :host([aria-invalid='true']:not([disabled]):hover) {
    border-color: ${oz};
  }
  
  :host(:${J.focusVisible}) {
    border-color: ${tN};
    box-shadow: 0 0 0 calc((${eT} - ${eS}) * 1px)
    ${tN};
  }
  
  :host([aria-invalid='true']:${J.focusVisible}) {
    border-color: ${oB};
    box-shadow: 0 0 0 calc((${eT} - ${eS}) * 1px)
    ${oB};
  }
  
  :host(:not([size]):not([multiple]):not([open]):${J.focusVisible}),
  :host([multiple]:${J.focusVisible}),
  :host([size]:${J.focusVisible}) {
    box-shadow: 0 0 0 calc((${eT} - ${eS}) * 1px)
    ${tN};
  }
  
  :host([aria-invalid='true']:not([size]):not([multiple]):not([open]):${J.focusVisible}),
  :host([aria-invalid='true'][multiple]:${J.focusVisible}),
  :host([aria-invalid='true'][size]:${J.focusVisible}) {
    box-shadow: 0 0 0 calc((${eT} - ${eS}) * 1px)
    ${oB};
  }
  
  :host(:not([multiple]):not([size]):${J.focusVisible}) ::slotted(${e.tagFor(J.ListboxOption)}[aria-selected="true"]:not([disabled])) {
    box-shadow: 0 0 0 calc(${eT} * 1px) inset ${tN};
    border-color: ${tN};
    background: ${tN};
    color: ${tM};
  }
    
  :host([disabled]) {
    cursor: ${J.disabledCursor};
    opacity: ${eD};
  }
  
  :host([disabled]) .control {
    cursor: ${J.disabledCursor};
    user-select: none;
  }
  
  :host([disabled]:hover) {
    background: ${t7};
    color: ${om};
    fill: currentcolor;
  }
  
  :host(:not([disabled])) .control:active {
    background: ${t4};
    border-color: ${tH};
    border-radius: calc(${ek} * 1px);
  }
  
  :host([open][position="above"]) .listbox {
    border-bottom-left-radius: 0;
    border-bottom-right-radius: 0;
    border-bottom: 0;
    bottom: calc(${o5} * 1px);
  }
  
  :host([open][position="below"]) .listbox {
    border-top-left-radius: 0;
    border-top-right-radius: 0;
    border-top: 0;
    top: calc(${o5} * 1px);
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
    ${rB}
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
    min-height: calc(${eF} * 4px);
    min-width: calc(${eF} * 4px);
    width: 1em;
  }
  
  ::slotted([role="option"]),
  ::slotted(option) {
    flex: 0 0 auto;
  }
  `.withBehaviors((0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
      :host(:not([disabled]):hover),
      :host(:not([disabled]):active) {
        border-color: ${Y.Highlight};
      }

      :host([aria-invalid='true']) {
        border-style: dashed;
      }
      
      :host(:not([disabled]):${J.focusVisible}) {
        background-color: ${Y.ButtonFace};
        box-shadow: 0 0 0 calc(${eT} * 1px) ${Y.Highlight};
        color: ${Y.ButtonText};
        fill: currentcolor;
        forced-color-adjust: none;
      }
      
      :host(:not([disabled]):${J.focusVisible}) .listbox {
        background: ${Y.ButtonFace};
      }
      
      :host([disabled]) {
        border-color: ${Y.GrayText};
        background-color: ${Y.ButtonFace};
        color: ${Y.GrayText};
        fill: currentcolor;
        opacity: 1;
        forced-color-adjust: none;
      }
      
      :host([disabled]:hover) {
        background: ${Y.ButtonFace};
      }
      
      :host([disabled]) .control {
        color: ${Y.GrayText};
        border-color: ${Y.GrayText};
      }
      
      :host([disabled]) .control .select-indicator {
        fill: ${Y.GrayText};
      }
      
      :host(:${J.focusVisible}) ::slotted([aria-selected="true"][role="option"]),
      :host(:${J.focusVisible}) ::slotted(option[aria-selected="true"]),
      :host(:${J.focusVisible}) ::slotted([aria-selected="true"][role="option"]:not([disabled])) {
        background: ${Y.Highlight};
        border-color: ${Y.ButtonText};
        box-shadow: 0 0 0 calc((${eT} - ${eS}) * 1px)
        ${Y.HighlightText};
        color: ${Y.HighlightText};
        fill: currentcolor;
      }
      
      .start,
      .end,
      .indicator,
      .select-indicator,
      ::slotted(svg) {
        color: ${Y.ButtonText};
        fill: currentcolor;
      }
      `))},rG=(e,t)=>(0,o1.css)`
  ${rM(e,t)}

  :host(:empty) .listbox {
    display: none;
  }

  :host([disabled]) *,
  :host([disabled]) {
    cursor: ${J.disabledCursor};
    user-select: none;
  }

  :host(:focus-within:not([disabled])) {
    border-color: ${tN};
    box-shadow: 0 0 0 calc((${eT} - ${eS}) * 1px)
      ${tN};
  }

  :host([aria-invalid='true']:focus-within:not([disabled])) {
    border-color: ${oB};
    box-shadow: 0 0 0 calc((${eT} - ${eS}) * 1px)
      ${oB};
  }

  .selected-value {
    -webkit-appearance: none;
    background: transparent;
    border: none;
    color: inherit;
    font-size: ${ez};
    line-height: ${ej};
    height: calc(100% - (${eS} * 1px));
    margin: auto 0;
    width: 100%;
  }

  .selected-value:hover,
    .selected-value:${J.focusVisible},
    .selected-value:disabled,
    .selected-value:active {
    outline: none;
  }
`;class rE extends J.Combobox{connectedCallback(){super.connectedCallback(),this.setAutoWidth()}slottedOptionsChanged(e,t){super.slottedOptionsChanged(e,t),this.setAutoWidth()}autoWidthChanged(e,t){t?this.setAutoWidth():this.style.removeProperty("width")}setAutoWidth(){if(!this.autoWidth||!this.isConnected)return;let e=this.listbox.getBoundingClientRect().width;0===e&&this.listbox.hidden&&(Object.assign(this.listbox.style,{visibility:"hidden"}),this.listbox.removeAttribute("hidden"),e=this.listbox.getBoundingClientRect().width,this.listbox.setAttribute("hidden",""),this.listbox.style.removeProperty("visibility")),e>0&&Object.assign(this.style,{width:`${e}px`})}maxHeightChanged(e,t){this.updateComputedStylesheet()}updateComputedStylesheet(){this.computedStylesheet&&this.$fastController.removeStyles(this.computedStylesheet);let e=Math.floor(this.maxHeight/oC.getValueFor(this)).toString();this.computedStylesheet=(0,o1.css)`
      :host {
        --listbox-max-height: ${e};
      }
    `,this.$fastController.addStyles(this.computedStylesheet)}}o7([(0,o1.attr)({attribute:"autowidth",mode:"boolean"})],rE.prototype,"autoWidth",void 0),o7([(0,o1.attr)({attribute:"minimal",mode:"boolean"})],rE.prototype,"minimal",void 0),o7([o1.attr],rE.prototype,"scale",void 0);let r_=rE.compose({baseName:"combobox",baseClass:J.Combobox,template:J.comboboxTemplate,styles:rG,shadowOptions:{delegatesFocus:!0},indicator:`
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
    `}),rq=(e,t)=>(0,o1.css)`
  :host {
    display: flex;
    position: relative;
    flex-direction: column;
  }
`,rW=(e,t)=>(0,o1.css)`
  :host {
    display: grid;
    padding: 1px 0;
    box-sizing: border-box;
    width: 100%;
    border-bottom: calc(${eS} * 1px) solid ${oF};
  }

  :host(.header) {
  }

  :host(.sticky-header) {
    background: ${tQ};
    position: sticky;
    top: 0;
  }
`,rU=(e,t)=>(0,o1.css)`
    :host {
      padding: calc(${eF} * 1px) calc(${eF} * 3px);
      color: ${om};
      box-sizing: border-box;
      font-family: ${ev};
      font-size: ${ez};
      line-height: ${ej};
      font-weight: 400;
      border: transparent calc(${eT} * 1px) solid;
      overflow: hidden;
      white-space: nowrap;
      border-radius: calc(${ek} * 1px);
    }

    :host(.column-header) {
      font-weight: 600;
    }

    :host(:${J.focusVisible}) {
      outline: calc(${eT} * 1px) solid ${tN};
      color: ${om};
    }
  `.withBehaviors((0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
      :host {
        forced-color-adjust: none;
        border-color: transparent;
        background: ${Y.Field};
        color: ${Y.FieldText};
      }

      :host(:${J.focusVisible}) {
        border-color: ${Y.FieldText};
        box-shadow: 0 0 0 2px inset ${Y.Field};
        color: ${Y.FieldText};
      }
    `));class rX extends J.DataGridCell{}let rZ=rX.compose({baseName:"data-grid-cell",baseClass:J.DataGridCell,template:J.dataGridCellTemplate,styles:rU});class rY extends J.DataGridRow{}let rJ=rY.compose({baseName:"data-grid-row",baseClass:J.DataGridRow,template:J.dataGridRowTemplate,styles:rW});class rK extends J.DataGrid{}let rQ=rK.compose({baseName:"data-grid",baseClass:J.DataGrid,template:J.dataGridTemplate,styles:rq});var r0=o(74291);class r1 extends J.FoundationElement{}class r2 extends(0,J.FormAssociated)(r1){constructor(){super(...arguments),this.proxy=document.createElement("input")}}let r5={toView(e){if(null==e)return null;let t=new Date(e);return"Invalid Date"===t.toString()?null:`${t.getFullYear().toString().padStart(4,"0")}-${(t.getMonth()+1).toString().padStart(2,"0")}-${t.getDate().toString().padStart(2,"0")}`},fromView(e){if(null==e)return null;let t=new Date(e);return"Invalid Date"===t.toString()?null:t}},r6="Invalid Date";class r3 extends r2{constructor(){super(...arguments),this.step=1,this.isUserInput=!1}readOnlyChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.readOnly=this.readOnly,this.validate())}autofocusChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.autofocus=this.autofocus,this.validate())}listChanged(){this.proxy instanceof HTMLInputElement&&(this.proxy.setAttribute("list",this.list),this.validate())}maxChanged(e,t){var o;this.max=t<(null!=(o=this.min)?o:t)?this.min:t,this.value=this.getValidValue(this.value)}minChanged(e,t){var o;this.min=t>(null!=(o=this.max)?o:t)?this.max:t,this.value=this.getValidValue(this.value)}get valueAsNumber(){return new Date(super.value).valueOf()}set valueAsNumber(e){this.value=new Date(e).toString()}get valueAsDate(){return new Date(super.value)}set valueAsDate(e){this.value=e.toString()}valueChanged(e,t){this.value=this.getValidValue(t),t===this.value&&(this.control&&!this.isUserInput&&(this.control.value=this.value),super.valueChanged(e,this.value),void 0===e||this.isUserInput||this.$emit("change"),this.isUserInput=!1)}getValidValue(e){var t,o;let r=new Date(e);return r.toString()===r6?r="":(r=(r=r>(null!=(t=this.max)?t:r)?this.max:r)<(null!=(o=this.min)?o:r)?this.min:r,r=`${r.getFullYear().toString().padStart(4,"0")}-${(r.getMonth()+1).toString().padStart(2,"0")}-${r.getDate().toString().padStart(2,"0")}`),r}stepUp(){let e=864e5*this.step,t=new Date(this.value);this.value=new Date(t.toString()!==r6?t.valueOf()+e:0).toString()}stepDown(){let e=864e5*this.step,t=new Date(this.value);this.value=new Date(t.toString()!==r6?Math.max(t.valueOf()-e,0):0).toString()}connectedCallback(){super.connectedCallback(),this.validate(),this.control.value=this.value,this.autofocus&&o1.DOM.queueUpdate(()=>{this.focus()}),this.appearance||(this.appearance="outline")}handleTextInput(){this.isUserInput=!0,this.value=this.control.value}handleChange(){this.$emit("change")}handleKeyDown(e){switch(e.key){case r0.I5:return this.stepUp(),!1;case r0.HX:return this.stepDown(),!1}return!0}handleBlur(){this.control.value=this.value}}o7([o1.attr],r3.prototype,"appearance",void 0),o7([(0,o1.attr)({attribute:"readonly",mode:"boolean"})],r3.prototype,"readOnly",void 0),o7([(0,o1.attr)({mode:"boolean"})],r3.prototype,"autofocus",void 0),o7([o1.attr],r3.prototype,"list",void 0),o7([(0,o1.attr)({converter:o1.nullableNumberConverter})],r3.prototype,"step",void 0),o7([(0,o1.attr)({converter:r5})],r3.prototype,"max",void 0),o7([(0,o1.attr)({converter:r5})],r3.prototype,"min",void 0),o7([o1.observable],r3.prototype,"defaultSlottedNodes",void 0),(0,J.applyMixins)(r3,J.StartEnd,J.DelegatesARIATextbox);let r4=(0,o1.css)`
  ${(0,J.display)("inline-block")} :host {
    font-family: ${ev};
    outline: none;
    user-select: none;
    /* Ensure to display focus highlight */
    margin: calc((${eT} - ${eS}) * 1px);
  }

  .root {
    box-sizing: border-box;
    position: relative;
    display: flex;
    flex-direction: row;
    color: ${om};
    background: ${t6};
    border-radius: calc(${ek} * 1px);
    border: calc(${eS} * 1px) solid ${oa};
    height: calc(${o5} * 1px);
  }

  :host([aria-invalid='true']) .root {
    border-color: ${oT};
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
    padding: 0 calc(${eF} * 2px + 1px);
    font-size: ${ez};
    line-height: ${ej};
  }

  .control:placeholder-shown {
    text-overflow: ellipsis;
  }

  .control:hover,
  .control:${J.focusVisible},
  .control:disabled,
  .control:active {
    outline: none;
  }

  .label {
    display: block;
    color: ${om};
    cursor: pointer;
    font-size: ${ez};
    line-height: ${ej};
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
    background: ${t3};
    border-color: ${oi};
  }

  :host([aria-invalid='true']:hover:not([disabled])) .root {
    border-color: ${oz};
  }

  :host(:active:not([disabled])) .root {
    background: ${t3};
    border-color: ${ol};
  }

  :host([aria-invalid='true']:active:not([disabled])) .root {
    border-color: ${oj};
  }

  :host(:focus-within:not([disabled])) .root {
    border-color: ${tN};
    box-shadow: 0 0 0 calc((${eT} - ${eS}) * 1px)
      ${tN};
  }

  :host([aria-invalid='true']:focus-within:not([disabled])) .root {
    border-color: ${oB};
    box-shadow: 0 0 0 calc((${eT} - ${eS}) * 1px)
      ${oB};
  }

  :host([appearance='filled']) .root {
    background: ${tQ};
  }

  :host([appearance='filled']:hover:not([disabled])) .root {
    background: ${t0};
  }

  :host([disabled]) .label,
  :host([readonly]) .label,
  :host([readonly]) .control,
  :host([disabled]) .control {
    cursor: ${J.disabledCursor};
  }

  :host([disabled]) {
    opacity: ${eD};
  }

  :host([disabled]) .control {
    border-color: ${o$};
  }
`.withBehaviors((0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
    .root,
    :host([appearance='filled']) .root {
      forced-color-adjust: none;
      background: ${Y.Field};
      border-color: ${Y.FieldText};
    }
    :host([aria-invalid='true']) .root {
      border-style: dashed;
    }
    :host(:hover:not([disabled])) .root,
    :host([appearance='filled']:hover:not([disabled])) .root,
    :host([appearance='filled']:hover) .root {
      background: ${Y.Field};
      border-color: ${Y.Highlight};
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
      border-color: ${Y.GrayText};
      background: ${Y.Field};
    }
    :host(:focus-within:enabled) .root {
      border-color: ${Y.Highlight};
      box-shadow: 0 0 0 calc((${eT} - ${eS}) * 1px)
        ${Y.Highlight};
    }
    input::placeholder {
      color: ${Y.GrayText};
    }
  `)),r9=(e,t)=>(0,o1.css)`
  ${r4}
`,r8=(e,t)=>(0,o1.html)`
  <template class="${e=>e.readOnly?"readonly":""}">
    <label
      part="label"
      for="control"
      class="${e=>e.defaultSlottedNodes&&e.defaultSlottedNodes.length?"label":"label label__hidden"}"
    >
      <slot
        ${(0,o1.slotted)({property:"defaultSlottedNodes",filter:J.whitespaceFilter})}
      ></slot>
    </label>
    <div class="root" part="root">
      ${(0,J.startSlotTemplate)(e,t)}
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
        ${(0,o1.ref)("control")}
      />
      ${(0,J.endSlotTemplate)(e,t)}
    </div>
  </template>
`,r7=r3.compose({baseName:"date-field",styles:r9,template:r8,shadowOptions:{delegatesFocus:!0}}),ae={toView:e=>null==e?null:null==e?void 0:e.toColorString(),fromView(e){if(null==e)return null;let t=f(e);return t?I.create(t.r,t.g,t.b):null}},at=(0,o1.css)`
  :host {
    background-color: ${tz};
    color: ${om};
  }
`.withBehaviors((0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
    :host {
      background-color: ${Y.ButtonFace};
      box-shadow: 0 0 0 1px ${Y.CanvasText};
      color: ${Y.ButtonText};
    }
  `));function ao(e){return(t,o)=>{t[o+"Changed"]=function(t,o){null!=o?e.setValueFor(this,o):e.deleteValueFor(this)}}}class ar extends J.FoundationElement{constructor(){super(),this.noPaint=!1;const e={handleChange:this.noPaintChanged.bind(this)};o1.Observable.getNotifier(this).subscribe(e,"fillColor"),o1.Observable.getNotifier(this).subscribe(e,"baseLayerLuminance")}noPaintChanged(){!this.noPaint&&(void 0!==this.fillColor||this.baseLayerLuminance)?this.$fastController.addStyles(at):this.$fastController.removeStyles(at)}}o7([(0,o1.attr)({attribute:"no-paint",mode:"boolean"})],ar.prototype,"noPaint",void 0),o7([(0,o1.attr)({attribute:"fill-color",converter:ae}),ao(tz)],ar.prototype,"fillColor",void 0),o7([(0,o1.attr)({attribute:"accent-color",converter:ae,mode:"fromView"}),ao(tf)],ar.prototype,"accentColor",void 0),o7([(0,o1.attr)({attribute:"neutral-color",converter:ae,mode:"fromView"}),ao(tg)],ar.prototype,"neutralColor",void 0),o7([(0,o1.attr)({attribute:"error-color",converter:ae,mode:"fromView"}),ao(oV)],ar.prototype,"errorColor",void 0),o7([(0,o1.attr)({converter:o1.nullableNumberConverter}),ao(ew)],ar.prototype,"density",void 0),o7([(0,o1.attr)({attribute:"design-unit",converter:o1.nullableNumberConverter}),ao(eF)],ar.prototype,"designUnit",void 0),o7([(0,o1.attr)({attribute:"direction"}),ao(eV)],ar.prototype,"direction",void 0),o7([(0,o1.attr)({attribute:"base-height-multiplier",converter:o1.nullableNumberConverter}),ao(e$)],ar.prototype,"baseHeightMultiplier",void 0),o7([(0,o1.attr)({attribute:"base-horizontal-spacing-multiplier",converter:o1.nullableNumberConverter}),ao(ex)],ar.prototype,"baseHorizontalSpacingMultiplier",void 0),o7([(0,o1.attr)({attribute:"control-corner-radius",converter:o1.nullableNumberConverter}),ao(ek)],ar.prototype,"controlCornerRadius",void 0),o7([(0,o1.attr)({attribute:"stroke-width",converter:o1.nullableNumberConverter}),ao(eS)],ar.prototype,"strokeWidth",void 0),o7([(0,o1.attr)({attribute:"focus-stroke-width",converter:o1.nullableNumberConverter}),ao(eT)],ar.prototype,"focusStrokeWidth",void 0),o7([(0,o1.attr)({attribute:"disabled-opacity",converter:o1.nullableNumberConverter}),ao(eD)],ar.prototype,"disabledOpacity",void 0),o7([(0,o1.attr)({attribute:"type-ramp-minus-2-font-size"}),ao(eH)],ar.prototype,"typeRampMinus2FontSize",void 0),o7([(0,o1.attr)({attribute:"type-ramp-minus-2-line-height"}),ao(eN)],ar.prototype,"typeRampMinus2LineHeight",void 0),o7([(0,o1.attr)({attribute:"type-ramp-minus-1-font-size"}),ao(eB)],ar.prototype,"typeRampMinus1FontSize",void 0),o7([(0,o1.attr)({attribute:"type-ramp-minus-1-line-height"}),ao(eL)],ar.prototype,"typeRampMinus1LineHeight",void 0),o7([(0,o1.attr)({attribute:"type-ramp-base-font-size"}),ao(ez)],ar.prototype,"typeRampBaseFontSize",void 0),o7([(0,o1.attr)({attribute:"type-ramp-base-line-height"}),ao(ej)],ar.prototype,"typeRampBaseLineHeight",void 0),o7([(0,o1.attr)({attribute:"type-ramp-plus-1-font-size"}),ao(eO)],ar.prototype,"typeRampPlus1FontSize",void 0),o7([(0,o1.attr)({attribute:"type-ramp-plus-1-line-height"}),ao(eI)],ar.prototype,"typeRampPlus1LineHeight",void 0),o7([(0,o1.attr)({attribute:"type-ramp-plus-2-font-size"}),ao(eR)],ar.prototype,"typeRampPlus2FontSize",void 0),o7([(0,o1.attr)({attribute:"type-ramp-plus-2-line-height"}),ao(eP)],ar.prototype,"typeRampPlus2LineHeight",void 0),o7([(0,o1.attr)({attribute:"type-ramp-plus-3-font-size"}),ao(eA)],ar.prototype,"typeRampPlus3FontSize",void 0),o7([(0,o1.attr)({attribute:"type-ramp-plus-3-line-height"}),ao(eM)],ar.prototype,"typeRampPlus3LineHeight",void 0),o7([(0,o1.attr)({attribute:"type-ramp-plus-4-font-size"}),ao(eG)],ar.prototype,"typeRampPlus4FontSize",void 0),o7([(0,o1.attr)({attribute:"type-ramp-plus-4-line-height"}),ao(eE)],ar.prototype,"typeRampPlus4LineHeight",void 0),o7([(0,o1.attr)({attribute:"type-ramp-plus-5-font-size"}),ao(e_)],ar.prototype,"typeRampPlus5FontSize",void 0),o7([(0,o1.attr)({attribute:"type-ramp-plus-5-line-height"}),ao(eq)],ar.prototype,"typeRampPlus5LineHeight",void 0),o7([(0,o1.attr)({attribute:"type-ramp-plus-6-font-size"}),ao(eW)],ar.prototype,"typeRampPlus6FontSize",void 0),o7([(0,o1.attr)({attribute:"type-ramp-plus-6-line-height"}),ao(eU)],ar.prototype,"typeRampPlus6LineHeight",void 0),o7([(0,o1.attr)({attribute:"accent-fill-rest-delta",converter:o1.nullableNumberConverter}),ao(eX)],ar.prototype,"accentFillRestDelta",void 0),o7([(0,o1.attr)({attribute:"accent-fill-hover-delta",converter:o1.nullableNumberConverter}),ao(eZ)],ar.prototype,"accentFillHoverDelta",void 0),o7([(0,o1.attr)({attribute:"accent-fill-active-delta",converter:o1.nullableNumberConverter}),ao(eY)],ar.prototype,"accentFillActiveDelta",void 0),o7([(0,o1.attr)({attribute:"accent-fill-focus-delta",converter:o1.nullableNumberConverter}),ao(eJ)],ar.prototype,"accentFillFocusDelta",void 0),o7([(0,o1.attr)({attribute:"accent-foreground-rest-delta",converter:o1.nullableNumberConverter}),ao(eK)],ar.prototype,"accentForegroundRestDelta",void 0),o7([(0,o1.attr)({attribute:"accent-foreground-hover-delta",converter:o1.nullableNumberConverter}),ao(eQ)],ar.prototype,"accentForegroundHoverDelta",void 0),o7([(0,o1.attr)({attribute:"accent-foreground-active-delta",converter:o1.nullableNumberConverter}),ao(e0)],ar.prototype,"accentForegroundActiveDelta",void 0),o7([(0,o1.attr)({attribute:"accent-foreground-focus-delta",converter:o1.nullableNumberConverter}),ao(e1)],ar.prototype,"accentForegroundFocusDelta",void 0),o7([(0,o1.attr)({attribute:"neutral-fill-rest-delta",converter:o1.nullableNumberConverter}),ao(e2)],ar.prototype,"neutralFillRestDelta",void 0),o7([(0,o1.attr)({attribute:"neutral-fill-hover-delta",converter:o1.nullableNumberConverter}),ao(e5)],ar.prototype,"neutralFillHoverDelta",void 0),o7([(0,o1.attr)({attribute:"neutral-fill-active-delta",converter:o1.nullableNumberConverter}),ao(e6)],ar.prototype,"neutralFillActiveDelta",void 0),o7([(0,o1.attr)({attribute:"neutral-fill-focus-delta",converter:o1.nullableNumberConverter}),ao(e3)],ar.prototype,"neutralFillFocusDelta",void 0),o7([(0,o1.attr)({attribute:"neutral-fill-input-rest-delta",converter:o1.nullableNumberConverter}),ao(e4)],ar.prototype,"neutralFillInputRestDelta",void 0),o7([(0,o1.attr)({attribute:"neutral-fill-input-hover-delta",converter:o1.nullableNumberConverter}),ao(e9)],ar.prototype,"neutralFillInputHoverDelta",void 0),o7([(0,o1.attr)({attribute:"neutral-fill-input-active-delta",converter:o1.nullableNumberConverter}),ao(e8)],ar.prototype,"neutralFillInputActiveDelta",void 0),o7([(0,o1.attr)({attribute:"neutral-fill-input-focus-delta",converter:o1.nullableNumberConverter}),ao(e7)],ar.prototype,"neutralFillInputFocusDelta",void 0),o7([(0,o1.attr)({attribute:"neutral-fill-stealth-rest-delta",converter:o1.nullableNumberConverter}),ao(te)],ar.prototype,"neutralFillStealthRestDelta",void 0),o7([(0,o1.attr)({attribute:"neutral-fill-stealth-hover-delta",converter:o1.nullableNumberConverter}),ao(tt)],ar.prototype,"neutralFillStealthHoverDelta",void 0),o7([(0,o1.attr)({attribute:"neutral-fill-stealth-active-delta",converter:o1.nullableNumberConverter}),ao(to)],ar.prototype,"neutralFillStealthActiveDelta",void 0),o7([(0,o1.attr)({attribute:"neutral-fill-stealth-focus-delta",converter:o1.nullableNumberConverter}),ao(tr)],ar.prototype,"neutralFillStealthFocusDelta",void 0),o7([(0,o1.attr)({attribute:"neutral-fill-strong-hover-delta",converter:o1.nullableNumberConverter}),ao(ti)],ar.prototype,"neutralFillStrongHoverDelta",void 0),o7([(0,o1.attr)({attribute:"neutral-fill-strong-active-delta",converter:o1.nullableNumberConverter}),ao(tl)],ar.prototype,"neutralFillStrongActiveDelta",void 0),o7([(0,o1.attr)({attribute:"neutral-fill-strong-focus-delta",converter:o1.nullableNumberConverter}),ao(ts)],ar.prototype,"neutralFillStrongFocusDelta",void 0),o7([(0,o1.attr)({attribute:"base-layer-luminance",converter:o1.nullableNumberConverter}),ao(ey)],ar.prototype,"baseLayerLuminance",void 0),o7([(0,o1.attr)({attribute:"neutral-fill-layer-rest-delta",converter:o1.nullableNumberConverter}),ao(tn)],ar.prototype,"neutralFillLayerRestDelta",void 0),o7([(0,o1.attr)({attribute:"neutral-stroke-divider-rest-delta",converter:o1.nullableNumberConverter}),ao(tp)],ar.prototype,"neutralStrokeDividerRestDelta",void 0),o7([(0,o1.attr)({attribute:"neutral-stroke-rest-delta",converter:o1.nullableNumberConverter}),ao(tc)],ar.prototype,"neutralStrokeRestDelta",void 0),o7([(0,o1.attr)({attribute:"neutral-stroke-hover-delta",converter:o1.nullableNumberConverter}),ao(td)],ar.prototype,"neutralStrokeHoverDelta",void 0),o7([(0,o1.attr)({attribute:"neutral-stroke-active-delta",converter:o1.nullableNumberConverter}),ao(th)],ar.prototype,"neutralStrokeActiveDelta",void 0),o7([(0,o1.attr)({attribute:"neutral-stroke-focus-delta",converter:o1.nullableNumberConverter}),ao(tu)],ar.prototype,"neutralStrokeFocusDelta",void 0);let aa=(e,t)=>(0,o1.html)` <slot></slot> `,ai=(e,t)=>(0,o1.css)`
  ${(0,J.display)("block")}
`,al=ar.compose({baseName:"design-system-provider",template:aa,styles:ai}),as=(e,t)=>(0,o1.css)`
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
    ${rB}
    margin-top: auto;
    margin-bottom: auto;
    width: var(--dialog-width);
    height: var(--dialog-height);
    background-color: ${tz};
    z-index: 1;
    border-radius: calc(${ek} * 1px);
    border: calc(${eS} * 1px) solid transparent;
  }
`;class an extends J.Dialog{}let ac=an.compose({baseName:"dialog",baseClass:J.Dialog,template:J.dialogTemplate,styles:as}),ad=(e,t)=>(0,o1.css)`
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
    background: ${tB};
    color: ${tR};
    font-family: ${ev};
    font-size: ${ez};
    border-radius: calc(${ek} * 1px);
    outline: none;
    cursor: pointer;
    margin: 16px 0;
    padding: 12px;
    max-width: max-content;
  }

  :host([appearance='accent']) .invoker:active {
    background: ${tH};
    color: ${tA};
  }

  :host([appearance='accent']) .invoker:hover {
    background: ${tL};
    color: ${tP};
  }

  :host([appearance='lightweight']) .invoker {
    background: transparent;
    color: ${tX};
    border-bottom: calc(${eS} * 1px) solid ${tX};
    cursor: pointer;
    width: max-content;
    margin: 16px 0;
  }

  :host([appearance='lightweight']) .invoker:active {
    border-bottom-color: ${tY};
  }

  :host([appearance='lightweight']) .invoker:hover {
    border-bottom-color: ${tZ};
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
`;class ah extends J.Disclosure{constructor(){super(...arguments),this.height=0,this.totalHeight=0}connectedCallback(){super.connectedCallback(),this.appearance||(this.appearance="accent")}appearanceChanged(e,t){e!==t&&(this.classList.add(t),this.classList.remove(e))}onToggle(){super.onToggle(),this.details.style.setProperty("height",`${this.disclosureHeight}px`)}setup(){super.setup();let e=()=>this.details.getBoundingClientRect().height;this.show(),this.totalHeight=e(),this.hide(),this.height=e(),this.expanded&&this.show()}get disclosureHeight(){return this.expanded?this.totalHeight:this.height}}o7([o1.attr],ah.prototype,"appearance",void 0);let au=ah.compose({baseName:"disclosure",baseClass:J.Disclosure,template:J.disclosureTemplate,styles:ad}),ap=(e,t)=>(0,o1.css)`
  ${(0,J.display)("block")} :host {
    box-sizing: content-box;
    height: 0;
    margin: calc(${eF} * 1px) 0;
    border-top: calc(${eS} * 1px) solid ${oF};
    border-left: none;
  }

  :host([orientation='vertical']) {
    height: 100%;
    margin: 0 calc(${eF} * 1px);
    border-top: none;
    border-left: calc(${eS} * 1px) solid ${oF};
  }
`;class ag extends J.Divider{}let ab=ag.compose({baseName:"divider",baseClass:J.Divider,template:J.dividerTemplate,styles:ap});class af extends J.ListboxElement{sizeChanged(e,t){super.sizeChanged(e,t),this.updateComputedStylesheet()}updateComputedStylesheet(){this.computedStylesheet&&this.$fastController.removeStyles(this.computedStylesheet);let e=`${this.size}`;this.computedStylesheet=(0,o1.css)`
      :host {
        --size: ${e};
      }
    `,this.$fastController.addStyles(this.computedStylesheet)}}let am=af.compose({baseName:"listbox",baseClass:J.ListboxElement,template:J.listboxTemplate,styles:rA}),av=(e,t)=>(0,o1.css)`
    ${(0,J.display)("block")} :host {
      --elevation: 11;
      background: ${tz};
      border: calc(${eS} * 1px) solid transparent;
      ${rB}
      margin: 0;
      border-radius: calc(${ek} * 1px);
      padding: calc(${eF} * 1px) 0;
      max-width: 368px;
      min-width: 64px;
    }

    :host([slot='submenu']) {
      width: max-content;
      margin: 0 calc(${eF} * 1px);
    }

    ::slotted(hr) {
      box-sizing: content-box;
      height: 0;
      margin: 0;
      border: none;
      border-top: calc(${eS} * 1px) solid ${oF};
    }
  `.withBehaviors((0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
      :host {
        background: ${Y.Canvas};
        border-color: ${Y.CanvasText};
      }
    `));class a$ extends J.Menu{connectedCallback(){super.connectedCallback(),tz.setValueFor(this,ty)}}let ax=a$.compose({baseName:"menu",baseClass:J.Menu,template:J.menuTemplate,styles:av}),ay=(e,t)=>(0,o1.css)`
    ${(0,J.display)("grid")} :host {
      contain: layout;
      overflow: visible;
      font-family: ${ev};
      outline: none;
      box-sizing: border-box;
      height: calc(${o5} * 1px);
      grid-template-columns: minmax(42px, auto) 1fr minmax(42px, auto);
      grid-template-rows: auto;
      justify-items: center;
      align-items: center;
      padding: 0;
      margin: 0 calc(${eF} * 1px);
      white-space: nowrap;
      background: ${t7};
      color: ${om};
      fill: currentcolor;
      cursor: pointer;
      font-size: ${ez};
      line-height: ${ej};
      border-radius: calc(${ek} * 1px);
      border: calc(${eT} * 1px) solid transparent;
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

    :host(:${J.focusVisible}) {
      border-color: ${tN};
      background: ${oo};
      color: ${om};
    }

    :host(:hover) {
      background: ${oe};
      color: ${om};
    }

    :host(:active) {
      background: ${ot};
    }

    :host([aria-checked='true']),
    :host(.expanded) {
      background: ${tQ};
      color: ${om};
    }

    :host([disabled]) {
      cursor: ${J.disabledCursor};
      opacity: ${eD};
    }

    :host([disabled]:hover) {
      color: ${om};
      fill: currentcolor;
      background: ${t7};
    }

    :host([disabled]:hover) .start,
    :host([disabled]:hover) .end,
    :host([disabled]:hover)::slotted(svg) {
      fill: ${om};
    }

    .expand-collapse-glyph {
      /* TODO: adaptive typography https://github.com/microsoft/fast/issues/2432 */
      width: calc((16 + ${ew}) * 1px);
      height: calc((16 + ${ew}) * 1px);
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
      width: ${eO};
      height: ${eO};
      */
    }

    :host(:hover) .start,
    :host(:hover) .end,
    :host(:hover)::slotted(svg),
    :host(:active) .start,
    :host(:active) .end,
    :host(:active)::slotted(svg) {
      fill: ${om};
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
      border: calc(${eS} * 1px) solid ${om};
    }

    :host([aria-checked='true']) .checkbox,
    :host([aria-checked='true']) .radio {
      background: ${tB};
      border-color: ${tB};
    }

    :host .checkbox {
      border-radius: calc(${ek} * 1px);
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
      color: ${ob};
    }

    :host([aria-checked='true']) .checkbox-indicator,
    :host([aria-checked='true']) ::slotted([slot='checkbox-indicator']) {
      width: 100%;
      height: 100%;
      display: block;
      fill: ${tR};
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
      background: ${tR};
      pointer-events: none;
    }

    :host([aria-checked='true']) ::slotted([slot='radio-indicator']) {
      display: block;
      pointer-events: none;
    }
  `.withBehaviors((0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
      :host {
        border-color: transparent;
        color: ${Y.ButtonText};
        forced-color-adjust: none;
      }

      :host(:hover) {
        background: ${Y.Highlight};
        color: ${Y.HighlightText};
      }

      :host(:hover) .start,
      :host(:hover) .end,
      :host(:hover)::slotted(svg),
      :host(:active) .start,
      :host(:active) .end,
      :host(:active)::slotted(svg) {
        fill: ${Y.HighlightText};
      }

      :host(.expanded) {
        background: ${Y.Highlight};
        border-color: ${Y.Highlight};
        color: ${Y.HighlightText};
      }

      :host(:${J.focusVisible}) {
        background: ${Y.Highlight};
        border-color: ${Y.ButtonText};
        box-shadow: 0 0 0 calc(${eT} * 1px) inset
          ${Y.HighlightText};
        color: ${Y.HighlightText};
        fill: currentcolor;
      }

      :host([disabled]),
      :host([disabled]:hover),
      :host([disabled]:hover) .start,
      :host([disabled]:hover) .end,
      :host([disabled]:hover)::slotted(svg) {
        background: ${Y.Canvas};
        color: ${Y.GrayText};
        fill: currentcolor;
        opacity: 1;
      }

      :host .expanded-toggle,
      :host .checkbox,
      :host .radio {
        border-color: ${Y.ButtonText};
        background: ${Y.HighlightText};
      }

      :host([checked='true']) .checkbox,
      :host([checked='true']) .radio {
        background: ${Y.HighlightText};
        border-color: ${Y.HighlightText};
      }

      :host(:hover) .expanded-toggle,
            :host(:hover) .checkbox,
            :host(:hover) .radio,
            :host(:${J.focusVisible}) .expanded-toggle,
            :host(:${J.focusVisible}) .checkbox,
            :host(:${J.focusVisible}) .radio,
            :host([checked="true"]:hover) .checkbox,
            :host([checked="true"]:hover) .radio,
            :host([checked="true"]:${J.focusVisible}) .checkbox,
            :host([checked="true"]:${J.focusVisible}) .radio {
        border-color: ${Y.HighlightText};
      }

      :host([aria-checked='true']) {
        background: ${Y.Highlight};
        color: ${Y.HighlightText};
      }

      :host([aria-checked='true']) .checkbox-indicator,
      :host([aria-checked='true']) ::slotted([slot='checkbox-indicator']),
      :host([aria-checked='true']) ::slotted([slot='radio-indicator']) {
        fill: ${Y.Highlight};
      }

      :host([aria-checked='true']) .radio-indicator {
        background: ${Y.Highlight};
      }

      ::slotted([slot='end']:not(svg)) {
        color: ${Y.ButtonText};
      }

      :host(:hover) ::slotted([slot="end"]:not(svg)),
            :host(:${J.focusVisible}) ::slotted([slot="end"]:not(svg)) {
        color: ${Y.HighlightText};
      }
    `),new rg((0,o1.css)`
        .expand-collapse-glyph {
          transform: rotate(0deg);
        }
      `,(0,o1.css)`
        .expand-collapse-glyph {
          transform: rotate(180deg);
        }
      `));class ak extends J.MenuItem{}let aw=ak.compose({baseName:"menu-item",baseClass:J.MenuItem,template:J.menuItemTemplate,styles:ay,checkboxIndicator:`
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
    `}),aF=(e,t)=>(0,o1.css)`
  ${r4}

  .controls {
    opacity: 0;
  }

  .step-up-glyph,
  .step-down-glyph {
    display: block;
    padding: 4px 10px;
    cursor: pointer;
  }

  .step-up-glyph:before,
  .step-down-glyph:before {
    content: '';
    display: block;
    border: solid transparent 6px;
  }

  .step-up-glyph:before {
    border-bottom-color: ${om};
  }

  .step-down-glyph:before {
    border-top-color: ${om};
  }

  :host(:hover:not([disabled])) .controls,
  :host(:focus-within:not([disabled])) .controls {
    opacity: 1;
  }
`;class aC extends J.NumberField{constructor(){super(...arguments),this.appearance="outline"}}o7([o1.attr],aC.prototype,"appearance",void 0);let aV=aC.compose({baseName:"number-field",baseClass:J.NumberField,styles:aF,template:J.numberFieldTemplate,shadowOptions:{delegatesFocus:!0},stepDownGlyph:`
        <span class="step-down-glyph" part="step-down-glyph"></span>
    `,stepUpGlyph:`
        <span class="step-up-glyph" part="step-up-glyph"></span>
    `}),aD=(e,t)=>(0,o1.css)`
    ${(0,J.display)("inline-flex")} :host {
      align-items: center;
      font-family: ${ev};
      border-radius: calc(${ek} * 1px);
      border: calc(${eT} * 1px) solid transparent;
      box-sizing: border-box;
      background: ${t7};
      color: ${om};
      cursor: pointer;
      flex: 0 0 auto;
      fill: currentcolor;
      font-size: ${ez};
      height: calc(${o5} * 1px);
      line-height: ${ej};
      margin: 0 calc((${eF} - ${eT}) * 1px);
      outline: none;
      overflow: hidden;
      padding: 0 1ch;
      user-select: none;
      white-space: nowrap;
    }

    :host(:not([disabled]):not([aria-selected='true']):hover) {
      background: ${oe};
    }

    :host(:not([disabled]):not([aria-selected='true']):active) {
      background: ${ot};
    }

    :host([aria-selected='true']) {
      background: ${tB};
      color: ${tR};
    }

    :host(:not([disabled])[aria-selected='true']:hover) {
      background: ${tL};
      color: ${tP};
    }

    :host(:not([disabled])[aria-selected='true']:active) {
      background: ${tH};
      color: ${tA};
    }

    :host([disabled]) {
      cursor: ${J.disabledCursor};
      opacity: ${eD};
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
      height: calc(${eF} * 4px);
      width: calc(${eF} * 4px);
    }

    ::slotted([slot='end']) {
      margin-inline-start: 1ch;
    }

    ::slotted([slot='start']) {
      margin-inline-end: 1ch;
    }

    :host([aria-checked='true'][aria-selected='false']) {
      border-color: ${oh};
    }

    :host([aria-checked='true'][aria-selected='true']) {
      border-color: ${oh};
      box-shadow: 0 0 0 calc(${eT} * 2 * 1px) inset
        ${op};
    }
  `.withBehaviors((0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
      :host {
        border-color: transparent;
        forced-color-adjust: none;
        color: ${Y.ButtonText};
        fill: currentcolor;
      }

      :host(:not([aria-selected='true']):hover),
      :host([aria-selected='true']) {
        background: ${Y.Highlight};
        color: ${Y.HighlightText};
      }

      :host([disabled]),
      :host([disabled][aria-selected='false']:hover) {
        background: ${Y.Canvas};
        color: ${Y.GrayText};
        fill: currentcolor;
        opacity: 1;
      }

      :host([aria-checked='true'][aria-selected='false']) {
        background: ${Y.ButtonFace};
        color: ${Y.ButtonText};
        border-color: ${Y.ButtonText};
      }

      :host([aria-checked='true'][aria-selected='true']),
      :host([aria-checked='true'][aria-selected='true']:hover) {
        background: ${Y.Highlight};
        color: ${Y.HighlightText};
        border-color: ${Y.ButtonText};
      }
    `));class aS extends J.ListboxOption{}let aT=aS.compose({baseName:"option",baseClass:J.ListboxOption,template:J.listboxOptionTemplate,styles:aD}),az=(e,t)=>(0,o1.css)`
    ${(0,J.display)("flex")} :host {
      align-items: center;
      outline: none;
      height: calc(${eF} * 1px);
      margin: calc(${eF} * 1px) 0;
    }

    .progress {
      background-color: ${tQ};
      border-radius: calc(${eF} * 1px);
      width: 100%;
      height: 100%;
      display: flex;
      align-items: center;
      position: relative;
    }

    .determinate {
      background-color: ${tX};
      border-radius: calc(${eF} * 1px);
      height: 100%;
      transition: all 0.2s ease-in-out;
      display: flex;
    }

    .indeterminate {
      height: 100%;
      border-radius: calc(${eF} * 1px);
      display: flex;
      width: 100%;
      position: relative;
      overflow: hidden;
    }

    .indeterminate-indicator-1 {
      position: absolute;
      opacity: 0;
      height: 100%;
      background-color: ${tX};
      border-radius: calc(${eF} * 1px);
      animation-timing-function: cubic-bezier(0.4, 0, 0.6, 1);
      width: 40%;
      animation: indeterminate-1 2s infinite;
    }

    .indeterminate-indicator-2 {
      position: absolute;
      opacity: 0;
      height: 100%;
      background-color: ${tX};
      border-radius: calc(${eF} * 1px);
      animation-timing-function: cubic-bezier(0.4, 0, 0.6, 1);
      width: 60%;
      animation: indeterminate-2 2s infinite;
    }

    :host([paused]) .indeterminate-indicator-1,
    :host([paused]) .indeterminate-indicator-2 {
      animation-play-state: paused;
      background-color: ${tQ};
    }

    :host([paused]) .determinate {
      background-color: ${ob};
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
  `.withBehaviors((0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
      .progress {
        forced-color-adjust: none;
        background-color: ${Y.Field};
        box-shadow: 0 0 0 1px inset ${Y.FieldText};
      }
      .determinate,
      .indeterminate-indicator-1,
      .indeterminate-indicator-2 {
        forced-color-adjust: none;
        background-color: ${Y.FieldText};
      }
      :host([paused]) .determinate,
      :host([paused]) .indeterminate-indicator-1,
      :host([paused]) .indeterminate-indicator-2 {
        background-color: ${Y.GrayText};
      }
    `));class aj extends J.BaseProgress{}let aB=aj.compose({baseName:"progress",baseClass:J.BaseProgress,template:J.progressTemplate,styles:az,indeterminateIndicator1:`
        <span class="indeterminate-indicator-1" part="indeterminate-indicator-1"></span>
    `,indeterminateIndicator2:`
        <span class="indeterminate-indicator-2" part="indeterminate-indicator-2"></span>
    `}),aL=(e,t)=>(0,o1.css)`
    ${(0,J.display)("flex")} :host {
      align-items: center;
      outline: none;
      height: calc(${o5} * 1px);
      width: calc(${o5} * 1px);
      margin: calc(${o5} * 1px) 0;
    }

    .progress {
      height: 100%;
      width: 100%;
    }

    .background {
      stroke: ${tQ};
      fill: none;
      stroke-width: 2px;
    }

    .determinate {
      stroke: ${tX};
      fill: none;
      stroke-width: 2px;
      stroke-linecap: round;
      transform-origin: 50% 50%;
      transform: rotate(-90deg);
      transition: all 0.2s ease-in-out;
    }

    .indeterminate-indicator-1 {
      stroke: ${tX};
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
      stroke: ${tQ};
    }

    :host([paused]) .determinate {
      stroke: ${ob};
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
  `.withBehaviors((0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
      .indeterminate-indicator-1,
      .determinate {
        stroke: ${Y.FieldText};
      }
      .background {
        stroke: ${Y.Field};
      }
      :host([paused]) .indeterminate-indicator-1 {
        stroke: ${Y.Field};
      }
      :host([paused]) .determinate {
        stroke: ${Y.GrayText};
      }
    `));class aH extends J.BaseProgress{}let aN=aH.compose({baseName:"progress-ring",baseClass:J.BaseProgress,template:J.progressRingTemplate,styles:aL,indeterminateIndicator:`
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
    `}),aO=(e,t)=>(0,o1.css)`
    ${(0,J.display)("inline-flex")} :host {
      --input-size: calc((${o5} / 2) + ${eF});
      align-items: center;
      outline: none;
      margin: calc(${eF} * 1px) 0;
      /* Chromium likes to select label text or the default slot when
                the radio is clicked. Maybe there is a better solution here? */
      user-select: none;
      position: relative;
      flex-direction: row;
      transition: all 0.2s ease-in-out;
    }

    .control {
      position: relative;
      width: calc((${o5} / 2 + ${eF}) * 1px);
      height: calc((${o5} / 2 + ${eF}) * 1px);
      box-sizing: border-box;
      border-radius: 999px;
      border: calc(${eS} * 1px) solid ${o$};
      background: ${t6};
      outline: none;
      cursor: pointer;
    }

    :host([aria-invalid='true']) .control {
      border-color: ${oT};
    }

    .label {
      font-family: ${ev};
      color: ${om};
      padding-inline-start: calc(${eF} * 2px + 2px);
      margin-inline-end: calc(${eF} * 2px + 2px);
      cursor: pointer;
      font-size: ${ez};
      line-height: ${ej};
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
      background: ${tR};
      fill: ${tR};
      opacity: 0;
      pointer-events: none;
    }

    :host(:not([disabled])) .control:hover {
      background: ${t3};
      border-color: ${ox};
    }

    :host([aria-invalid='true']:not([disabled])) .control:hover {
      border-color: ${oz};
    }

    :host(:not([disabled])) .control:active {
      background: ${t4};
      border-color: ${oy};
    }

    :host([aria-invalid='true']:not([disabled])) .control:active {
      border-color: ${oj};
    }

    :host(:${J.focusVisible}) .control {
      outline: solid calc(${eT} * 1px) ${tN};
    }

    :host([aria-invalid='true']:${J.focusVisible}) .control {
      outline-color: ${oB};
    }

    :host([aria-checked='true']) .control {
      background: ${tB};
      border: calc(${eS} * 1px) solid ${tB};
    }

    :host([aria-invalid='true'][aria-checked='true']) .control {
      background-color: ${oT};
      border-color: ${oT};
    }

    :host([aria-checked='true']:not([disabled])) .control:hover {
      background: ${tL};
      border: calc(${eS} * 1px) solid ${tL};
    }

    :host([aria-invalid='true'][aria-checked='true']:not([disabled]))
      .control:hover {
      background-color: ${oz};
      border-color: ${oz};
    }

    :host([aria-checked='true']:not([disabled]))
      .control:hover
      .checked-indicator {
      background: ${tP};
      fill: ${tP};
    }

    :host([aria-checked='true']:not([disabled])) .control:active {
      background: ${tH};
      border: calc(${eS} * 1px) solid ${tH};
    }

    :host([aria-invalid='true'][aria-checked='true']:not([disabled]))
      .control:active {
      background-color: ${oj};
      border-color: ${oj};
    }

    :host([aria-checked='true']:not([disabled]))
      .control:active
      .checked-indicator {
      background: ${tA};
      fill: ${tA};
    }

    :host([aria-checked="true"]:${J.focusVisible}:not([disabled])) .control {
      outline-offset: 2px;
      outline: solid calc(${eT} * 1px) ${tN};
    }

    :host([aria-invalid='true'][aria-checked="true"]:${J.focusVisible}:not([disabled])) .control {
      outline-color: ${oB};
    }

    :host([disabled]) .label,
    :host([readonly]) .label,
    :host([readonly]) .control,
    :host([disabled]) .control {
      cursor: ${J.disabledCursor};
    }

    :host([aria-checked='true']) .checked-indicator {
      opacity: 1;
    }

    :host([disabled]) {
      opacity: ${eD};
    }
  `.withBehaviors((0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
      .control,
      :host([aria-checked='true']:not([disabled])) .control {
        forced-color-adjust: none;
        border-color: ${Y.FieldText};
        background: ${Y.Field};
      }
      :host([aria-invalid='true']) {
        border-style: dashed;
      }
      :host(:not([disabled])) .control:hover {
        border-color: ${Y.Highlight};
        background: ${Y.Field};
      }
      :host([aria-checked='true']:not([disabled])) .control:hover,
      :host([aria-checked='true']:not([disabled])) .control:active {
        border-color: ${Y.Highlight};
        background: ${Y.Highlight};
      }
      :host([aria-checked='true']) .checked-indicator {
        background: ${Y.Highlight};
        fill: ${Y.Highlight};
      }
      :host([aria-checked='true']:not([disabled]))
        .control:hover
        .checked-indicator,
      :host([aria-checked='true']:not([disabled]))
        .control:active
        .checked-indicator {
        background: ${Y.HighlightText};
        fill: ${Y.HighlightText};
      }
      :host(:${J.focusVisible}) .control {
        border-color: ${Y.Highlight};
        outline-offset: 2px;
        outline: solid calc(${eT} * 1px) ${Y.FieldText};
      }
      :host([aria-checked="true"]:${J.focusVisible}:not([disabled])) .control {
        border-color: ${Y.Highlight};
        outline: solid calc(${eT} * 1px) ${Y.FieldText};
      }
      :host([disabled]) {
        forced-color-adjust: none;
        opacity: 1;
      }
      :host([disabled]) .label {
        color: ${Y.GrayText};
      }
      :host([disabled]) .control,
      :host([aria-checked='true'][disabled]) .control:hover,
      .control:active {
        background: ${Y.Field};
        border-color: ${Y.GrayText};
      }
      :host([disabled]) .checked-indicator,
      :host([aria-checked='true'][disabled]) .control:hover .checked-indicator {
        fill: ${Y.GrayText};
        background: ${Y.GrayText};
      }
    `)),aI=(e,t)=>(0,o1.html)`
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
      <slot ${(0,o1.slotted)("defaultSlottedNodes")}></slot>
    </label>
  </template>
`;class aR extends J.Radio{}let aP=aR.compose({baseName:"radio",baseClass:J.Radio,template:aI,styles:aO,checkedIndicator:`
        <div part="checked-indicator" class="checked-indicator"></div>
    `}),aA=(e,t)=>(0,o1.css)`
  ${(0,J.display)("flex")} :host {
    align-items: flex-start;
    margin: calc(${eF} * 1px) 0;
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
`;class aM extends J.RadioGroup{constructor(){super(),o1.Observable.getNotifier(this).subscribe({handleChange(e,t){"slottedRadioButtons"===t&&e.ariaInvalidChanged()}},"slottedRadioButtons")}ariaInvalidChanged(){this.slottedRadioButtons&&this.slottedRadioButtons.forEach(e=>{var t;e.setAttribute("aria-invalid",null!=(t=this.getAttribute("aria-invalid"))?t:"false")})}}let aG=aM.compose({baseName:"radio-group",baseClass:J.RadioGroup,template:J.radioGroupTemplate,styles:aA}),aE=J.DesignToken.create("clear-button-hover").withDefault(e=>{let t=t8.getValueFor(e),o=tK.getValueFor(e);return t.evaluate(e,o.evaluate(e).hover).hover}),a_=J.DesignToken.create("clear-button-active").withDefault(e=>{let t=t8.getValueFor(e),o=tK.getValueFor(e);return t.evaluate(e,o.evaluate(e).hover).active}),aq=(e,t)=>(0,o1.css)`
  ${r4}

  .control::-webkit-search-cancel-button {
    -webkit-appearance: none;
  }

  .control:hover,
    .control:${J.focusVisible},
    .control:disabled,
    .control:active {
    outline: none;
  }

  .clear-button {
    height: calc(100% - 2px);
    opacity: 0;
    margin: 1px;
    background: transparent;
    color: ${om};
    fill: currentcolor;
    border: none;
    border-radius: calc(${ek} * 1px);
    min-width: calc(${o5} * 1px);
    font-size: ${ez};
    line-height: ${ej};
    outline: none;
    font-family: ${ev};
    padding: 0 calc((10 + (${eF} * 2 * ${ew})) * 1px);
  }

  .clear-button:hover {
    background: ${oe};
  }

  .clear-button:active {
    background: ${ot};
  }

  :host([appearance='filled']) .clear-button:hover {
    background: ${aE};
  }

  :host([appearance='filled']) .clear-button:active {
    background: ${a_};
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
`;class aW extends J.Search{constructor(){super(...arguments),this.appearance="outline"}}o7([o1.attr],aW.prototype,"appearance",void 0);let aU=aW.compose({baseName:"search",baseClass:J.Search,template:J.searchTemplate,styles:aq,shadowOptions:{delegatesFocus:!0}});class aX extends J.Select{constructor(){super(...arguments),this.listboxScrollWidth=""}autoWidthChanged(e,t){t?this.setAutoWidth():this.style.removeProperty("width")}setAutoWidth(){if(!this.autoWidth||!this.isConnected)return;let e=this.listbox.getBoundingClientRect().width;0===e&&this.listbox.hidden&&(Object.assign(this.listbox.style,{visibility:"hidden"}),this.listbox.removeAttribute("hidden"),e=this.listbox.getBoundingClientRect().width,this.listbox.setAttribute("hidden",""),this.listbox.style.removeProperty("visibility")),e>0&&Object.assign(this.style,{width:`${e}px`})}connectedCallback(){super.connectedCallback(),this.setAutoWidth(),this.listbox&&tz.setValueFor(this.listbox,ty)}slottedOptionsChanged(e,t){super.slottedOptionsChanged(e,t),this.setAutoWidth()}get listboxMaxHeight(){return Math.floor(this.maxHeight/oC.getValueFor(this)).toString()}listboxScrollWidthChanged(){this.updateComputedStylesheet()}get selectSize(){var e;return`${null!=(e=this.size)?e:4*!!this.multiple}`}multipleChanged(e,t){super.multipleChanged(e,t),this.updateComputedStylesheet()}maxHeightChanged(e,t){this.collapsible&&this.updateComputedStylesheet()}setPositioning(){super.setPositioning(),this.updateComputedStylesheet()}sizeChanged(e,t){(super.sizeChanged(e,t),this.updateComputedStylesheet(),this.collapsible)?requestAnimationFrame(()=>{this.listbox.style.setProperty("display","flex"),this.listbox.style.setProperty("overflow","visible"),this.listbox.style.setProperty("visibility","hidden"),this.listbox.style.setProperty("width","auto"),this.listbox.hidden=!1,this.listboxScrollWidth=`${this.listbox.scrollWidth}`,this.listbox.hidden=!0,this.listbox.style.removeProperty("display"),this.listbox.style.removeProperty("overflow"),this.listbox.style.removeProperty("visibility"),this.listbox.style.removeProperty("width")}):this.listboxScrollWidth=""}updateComputedStylesheet(){this.computedStylesheet&&this.$fastController.removeStyles(this.computedStylesheet),this.computedStylesheet=(0,o1.css)`
      :host {
        --listbox-max-height: ${this.listboxMaxHeight};
        --listbox-scroll-width: ${this.listboxScrollWidth};
        --size: ${this.selectSize};
      }
    `,this.$fastController.addStyles(this.computedStylesheet)}}o7([(0,o1.attr)({attribute:"autowidth",mode:"boolean"})],aX.prototype,"autoWidth",void 0),o7([(0,o1.attr)({attribute:"minimal",mode:"boolean"})],aX.prototype,"minimal",void 0),o7([o1.attr],aX.prototype,"scale",void 0),o7([o1.observable],aX.prototype,"listboxScrollWidth",void 0);let aZ=aX.compose({baseName:"select",baseClass:J.Select,template:J.selectTemplate,styles:rM,indicator:`
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
    `}),aY=(e,t)=>(0,o1.css)`
    ${(0,J.display)("block")} :host {
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
      border-radius: calc(${ek} * 1px);
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

    ${(0,J.display)("block")} span.shimmer {
      position: absolute;
      width: 100%;
      height: 100%;
      background-image: var(
        --skeleton-animation-gradient,
        var(--skeleton-animation-gradient-default)
      );
      background-size: 0px 0px / 90% 100%;
      background-repeat: no-repeat;
      background-color: var(--skeleton-animation-fill, ${tQ});
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
  `.withBehaviors((0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
      :host {
        forced-color-adjust: none;
        background-color: ${Y.ButtonFace};
        box-shadow: 0 0 0 1px ${Y.ButtonText};
      }

      ${(0,J.display)("block")} span.shimmer {
        display: none;
      }
    `));class aJ extends J.Skeleton{}let aK=aJ.compose({baseName:"skeleton",baseClass:J.Skeleton,template:J.skeletonTemplate,styles:aY}),aQ=(0,o1.css)`
  .track-start {
    left: 0;
  }
`,a0=(0,o1.css)`
  .track-start {
    right: 0;
  }
`,a1=(e,t)=>(0,o1.css)`
    :host([hidden]) {
      display: none;
    }

    ${(0,J.display)("inline-grid")} :host {
      --thumb-size: calc(${o5} * 0.5 - ${eF});
      --thumb-translate: calc(
        var(--thumb-size) * -0.5 + var(--track-width) / 2
      );
      --track-overhang: calc((${eF} / 2) * -1);
      --track-width: ${eF};
      --jp-slider-height: calc(var(--thumb-size) * 10);
      align-items: center;
      width: 100%;
      margin: calc(${eF} * 1px) 0;
      user-select: none;
      box-sizing: border-box;
      border-radius: calc(${ek} * 1px);
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

    :host(:${J.focusVisible}) .thumb-cursor {
      box-shadow:
        0 0 0 2px ${tz},
        0 0 0 calc((2 + ${eT}) * 1px) ${tN};
    }

    :host([aria-invalid='true']:${J.focusVisible}) .thumb-cursor {
      box-shadow:
        0 0 0 2px ${tz},
        0 0 0 calc((2 + ${eT}) * 1px) ${oB};
    }

    .thumb-container {
      position: absolute;
      height: calc(var(--thumb-size) * 1px);
      width: calc(var(--thumb-size) * 1px);
      transition: all 0.2s ease;
      color: ${om};
      fill: currentcolor;
    }
    .thumb-cursor {
      border: none;
      width: calc(var(--thumb-size) * 1px);
      height: calc(var(--thumb-size) * 1px);
      background: ${om};
      border-radius: calc(${ek} * 1px);
    }
    .thumb-cursor:hover {
      background: ${om};
      border-color: ${ox};
    }
    .thumb-cursor:active {
      background: ${om};
    }
    .track-start {
      background: ${tX};
      position: absolute;
      height: 100%;
      left: 0;
      border-radius: calc(${ek} * 1px);
    }
    :host([aria-invalid='true']) .track-start {
      background-color: ${oT};
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
      background: ${o$};
      position: absolute;
      border-radius: calc(${ek} * 1px);
    }
    :host([orientation='vertical']) {
      height: calc(var(--fast-slider-height) * 1px);
      min-height: calc(var(--thumb-size) * 1px);
      min-width: calc(${eF} * 20px);
    }
    :host([orientation='vertical']) .track-start {
      height: auto;
      width: 100%;
      top: 0;
    }
    :host([disabled]),
    :host([readonly]) {
      cursor: ${J.disabledCursor};
    }
    :host([disabled]) {
      opacity: ${eD};
    }
  `.withBehaviors(new rg(aQ,a0),(0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
      .thumb-cursor {
        forced-color-adjust: none;
        border-color: ${Y.FieldText};
        background: ${Y.FieldText};
      }
      .thumb-cursor:hover,
      .thumb-cursor:active {
        background: ${Y.Highlight};
      }
      .track {
        forced-color-adjust: none;
        background: ${Y.FieldText};
      }
      :host(:${J.focusVisible}) .thumb-cursor {
        border-color: ${Y.Highlight};
      }
      :host([disabled]) {
        opacity: 1;
      }
      :host([disabled]) .track,
      :host([disabled]) .thumb-cursor {
        forced-color-adjust: none;
        background: ${Y.GrayText};
      }

      :host(:${J.focusVisible}) .thumb-cursor {
        background: ${Y.Highlight};
        border-color: ${Y.Highlight};
        box-shadow:
          0 0 0 2px ${Y.Field},
          0 0 0 4px ${Y.FieldText};
      }
    `));class a2 extends J.Slider{}let a5=a2.compose({baseName:"slider",baseClass:J.Slider,template:J.sliderTemplate,styles:a1,thumb:`
        <div class="thumb-cursor"></div>
    `});var a6=o(67002);let a3=(0,o1.css)`
  :host {
    align-self: start;
    grid-row: 2;
    margin-top: -2px;
    height: calc((${o5} / 2 + ${eF}) * 1px);
    width: auto;
  }
  .container {
    grid-template-rows: auto auto;
    grid-template-columns: 0;
  }
  .label {
    margin: 2px 0;
  }
`,a4=(0,o1.css)`
  :host {
    justify-self: start;
    grid-column: 2;
    margin-left: 2px;
    height: auto;
    width: calc((${o5} / 2 + ${eF}) * 1px);
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
    margin-left: calc((${eF} / 2) * 3px);
    align-self: center;
  }
`,a9=(e,t)=>(0,o1.css)`
    ${(0,J.display)("block")} :host {
      font-family: ${ev};
      color: ${om};
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
      width: calc((${eF} / 4) * 1px);
      height: calc(${o5} * 0.25 * 1px);
      background: ${o$};
      justify-self: center;
    }
    :host(.disabled) {
      opacity: ${eD};
    }
  `.withBehaviors((0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
      .mark {
        forced-color-adjust: none;
        background: ${Y.FieldText};
      }
      :host(.disabled) {
        forced-color-adjust: none;
        opacity: 1;
      }
      :host(.disabled) .label {
        color: ${Y.GrayText};
      }
      :host(.disabled) .mark {
        background: ${Y.GrayText};
      }
    `));class a8 extends J.SliderLabel{sliderOrientationChanged(){this.sliderOrientation===a6.t.horizontal?(this.$fastController.addStyles(a3),this.$fastController.removeStyles(a4)):(this.$fastController.addStyles(a4),this.$fastController.removeStyles(a3))}}let a7=a8.compose({baseName:"slider-label",baseClass:J.SliderLabel,template:J.sliderLabelTemplate,styles:a9}),ie=(e,t)=>(0,o1.css)`
    :host([hidden]) {
      display: none;
    }

    ${(0,J.display)("inline-flex")} :host {
      align-items: center;
      outline: none;
      font-family: ${ev};
      margin: calc(${eF} * 1px) 0;
      ${""} user-select: none;
    }

    :host([disabled]) {
      opacity: ${eD};
    }

    :host([disabled]) .label,
    :host([readonly]) .label,
    :host([readonly]) .switch,
    :host([disabled]) .switch {
      cursor: ${J.disabledCursor};
    }

    .switch {
      position: relative;
      outline: none;
      box-sizing: border-box;
      width: calc(${o5} * 1px);
      height: calc((${o5} / 2 + ${eF}) * 1px);
      background: ${t6};
      border-radius: calc(${ek} * 1px);
      border: calc(${eS} * 1px) solid ${o$};
    }

    :host([aria-invalid='true']) .switch {
      border-color: ${oT};
    }

    .switch:hover {
      background: ${t3};
      border-color: ${ox};
      cursor: pointer;
    }

    :host([disabled]) .switch:hover,
    :host([readonly]) .switch:hover {
      background: ${t3};
      border-color: ${ox};
      cursor: ${J.disabledCursor};
    }

    :host([aria-invalid='true'][disabled]) .switch:hover,
    :host([aria-invalid='true'][readonly]) .switch:hover {
      border-color: ${oz};
    }

    :host(:not([disabled])) .switch:active {
      background: ${t4};
      border-color: ${oy};
    }

    :host([aria-invalid='true']:not([disabled])) .switch:active {
      border-color: ${oj};
    }

    :host(:${J.focusVisible}) .switch {
      outline-offset: 2px;
      outline: solid calc(${eT} * 1px) ${tN};
    }

    :host([aria-invalid='true']:${J.focusVisible}) .switch {
      outline-color: ${oB};
    }

    .checked-indicator {
      position: absolute;
      top: 5px;
      bottom: 5px;
      background: ${om};
      border-radius: calc(${ek} * 1px);
      transition: all 0.2s ease-in-out;
    }

    .status-message {
      color: ${om};
      cursor: pointer;
      font-size: ${ez};
      line-height: ${ej};
    }

    :host([disabled]) .status-message,
    :host([readonly]) .status-message {
      cursor: ${J.disabledCursor};
    }

    .label {
      color: ${om};
      margin-inline-end: calc(${eF} * 2px + 2px);
      font-size: ${ez};
      line-height: ${ej};
      cursor: pointer;
    }

    .label__hidden {
      display: none;
      visibility: hidden;
    }

    ::slotted([slot='checked-message']),
    ::slotted([slot='unchecked-message']) {
      margin-inline-start: calc(${eF} * 2px + 2px);
    }

    :host([aria-checked='true']) .checked-indicator {
      background: ${tR};
    }

    :host([aria-checked='true']) .switch {
      background: ${tB};
      border-color: ${tB};
    }

    :host([aria-checked='true']:not([disabled])) .switch:hover {
      background: ${tL};
      border-color: ${tL};
    }

    :host([aria-invalid='true'][aria-checked='true']) .switch {
      background-color: ${oT};
      border-color: ${oT};
    }

    :host([aria-invalid='true'][aria-checked='true']:not([disabled]))
      .switch:hover {
      background-color: ${oz};
      border-color: ${oz};
    }

    :host([aria-checked='true']:not([disabled]))
      .switch:hover
      .checked-indicator {
      background: ${tP};
    }

    :host([aria-checked='true']:not([disabled])) .switch:active {
      background: ${tH};
      border-color: ${tH};
    }

    :host([aria-invalid='true'][aria-checked='true']:not([disabled]))
      .switch:active {
      background-color: ${oj};
      border-color: ${oj};
    }

    :host([aria-checked='true']:not([disabled]))
      .switch:active
      .checked-indicator {
      background: ${tA};
    }

    :host([aria-checked="true"]:${J.focusVisible}:not([disabled])) .switch {
      outline: solid calc(${eT} * 1px) ${tN};
    }

    :host([aria-invalid='true'][aria-checked="true"]:${J.focusVisible}:not([disabled])) .switch {
      outline-color: ${oB};
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
  `.withBehaviors((0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
      .checked-indicator,
      :host(:not([disabled])) .switch:active .checked-indicator {
        forced-color-adjust: none;
        background: ${Y.FieldText};
      }
      .switch {
        forced-color-adjust: none;
        background: ${Y.Field};
        border-color: ${Y.FieldText};
      }
      :host([aria-invalid='true']) .switch {
        border-style: dashed;
      }
      :host(:not([disabled])) .switch:hover {
        background: ${Y.HighlightText};
        border-color: ${Y.Highlight};
      }
      :host([aria-checked='true']) .switch {
        background: ${Y.Highlight};
        border-color: ${Y.Highlight};
      }
      :host([aria-checked='true']:not([disabled])) .switch:hover,
      :host(:not([disabled])) .switch:active {
        background: ${Y.HighlightText};
        border-color: ${Y.Highlight};
      }
      :host([aria-checked='true']) .checked-indicator {
        background: ${Y.HighlightText};
      }
      :host([aria-checked='true']:not([disabled]))
        .switch:hover
        .checked-indicator {
        background: ${Y.Highlight};
      }
      :host([disabled]) {
        opacity: 1;
      }
      :host(:${J.focusVisible}) .switch {
        border-color: ${Y.Highlight};
        outline-offset: 2px;
        outline: solid calc(${eT} * 1px) ${Y.FieldText};
      }
      :host([aria-checked="true"]:${J.focusVisible}:not([disabled])) .switch {
        outline: solid calc(${eT} * 1px) ${Y.FieldText};
      }
      :host([disabled]) .checked-indicator {
        background: ${Y.GrayText};
      }
      :host([disabled]) .switch {
        background: ${Y.Field};
        border-color: ${Y.GrayText};
      }
    `),new rg((0,o1.css)`
        .checked-indicator {
          left: 5px;
          right: calc(((${o5} / 2) + 1) * 1px);
        }

        :host([aria-checked='true']) .checked-indicator {
          left: calc(((${o5} / 2) + 1) * 1px);
          right: 5px;
        }
      `,(0,o1.css)`
        .checked-indicator {
          right: 5px;
          left: calc(((${o5} / 2) + 1) * 1px);
        }

        :host([aria-checked='true']) .checked-indicator {
          right: calc(((${o5} / 2) + 1) * 1px);
          left: 5px;
        }
      `));class it extends J.Switch{}let io=it.compose({baseName:"switch",baseClass:J.Switch,template:J.switchTemplate,styles:ie,switch:`
        <span class="checked-indicator" part="checked-indicator"></span>
    `}),ir=(e,t)=>(0,o1.css)`
  ${(0,J.display)("block")} :host {
    box-sizing: border-box;
    font-size: ${ez};
    line-height: ${ej};
    padding: 0 calc((6 + (${eF} * 2 * ${ew})) * 1px);
  }
`;class ia extends J.TabPanel{}let ii=ia.compose({baseName:"tab-panel",baseClass:J.TabPanel,template:J.tabPanelTemplate,styles:ir}),il=(e,t)=>(0,o1.css)`
    ${(0,J.display)("inline-flex")} :host {
      box-sizing: border-box;
      font-family: ${ev};
      font-size: ${ez};
      line-height: ${ej};
      height: calc(${o5} * 1px);
      padding: calc(${eF} * 5px) calc(${eF} * 4px);
      color: ${ob};
      fill: currentcolor;
      border-radius: 0 0 calc(${ek} * 1px)
        calc(${ek} * 1px);
      border: calc(${eS} * 1px) solid transparent;
      align-items: center;
      justify-content: center;
      grid-row: 2;
      cursor: pointer;
    }

    :host(:hover) {
      color: ${om};
      fill: currentcolor;
    }

    :host(:active) {
      color: ${om};
      fill: currentcolor;
    }

    :host([disabled]) {
      cursor: ${J.disabledCursor};
      opacity: ${eD};
    }

    :host([disabled]:hover) {
      color: ${ob};
      background: ${t7};
    }

    :host([aria-selected='true']) {
      background: ${tQ};
      color: ${om};
      fill: currentcolor;
    }

    :host([aria-selected='true']:hover) {
      background: ${t0};
      color: ${om};
      fill: currentcolor;
    }

    :host([aria-selected='true']:active) {
      background: ${t1};
      color: ${om};
      fill: currentcolor;
    }

    :host(:${J.focusVisible}) {
      outline: none;
      border-color: ${tN};
      box-shadow: 0 0 0 calc((${eT} - ${eS}) * 1px)
        ${tN};
    }

    :host(:focus) {
      outline: none;
    }

    :host(.vertical) {
      justify-content: end;
      grid-column: 2;
      border-bottom-left-radius: 0;
      border-top-right-radius: calc(${ek} * 1px);
    }

    :host(.vertical[aria-selected='true']) {
      z-index: 2;
    }

    :host(.vertical:hover) {
      color: ${om};
    }

    :host(.vertical:active) {
      color: ${om};
    }

    :host(.vertical:hover[aria-selected='true']) {
    }
  `.withBehaviors((0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
      :host {
        forced-color-adjust: none;
        border-color: transparent;
        color: ${Y.ButtonText};
        fill: currentcolor;
      }
      :host(:hover),
      :host(.vertical:hover),
      :host([aria-selected='true']:hover) {
        background: ${Y.Highlight};
        color: ${Y.HighlightText};
        fill: currentcolor;
      }
      :host([aria-selected='true']) {
        background: ${Y.HighlightText};
        color: ${Y.Highlight};
        fill: currentcolor;
      }
      :host(:${J.focusVisible}) {
        border-color: ${Y.ButtonText};
        box-shadow: none;
      }
      :host([disabled]),
      :host([disabled]:hover) {
        opacity: 1;
        color: ${Y.GrayText};
        background: ${Y.ButtonFace};
      }
    `));class is extends J.Tab{}let ic=is.compose({baseName:"tab",baseClass:J.Tab,template:J.tabTemplate,styles:il}),id=(e,t)=>(0,o1.css)`
    ${(0,J.display)("grid")} :host {
      box-sizing: border-box;
      font-family: ${ev};
      font-size: ${ez};
      line-height: ${ej};
      color: ${om};
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
      padding: calc(${eF} * 4px) calc(${eF} * 4px) 0;
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
      background: ${tB};
      margin-top: 0;
      border-radius: calc(${ek} * 1px)
        calc(${ek} * 1px) 0 0;
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
      padding: 0 calc(${eF} * 4px)
        calc((${o5} - ${eF}) * 1px) 0;
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
      background: ${tB};
      border-radius: calc(${ek} * 1px) 0 0
        calc(${ek} * 1px);
    }

    :host([orientation='vertical']) .activeIndicatorTransition {
      transition: transform 0.01s ease-in-out;
    }
  `.withBehaviors((0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
      .activeIndicator,
      :host([orientation='vertical']) .activeIndicator {
        forced-color-adjust: none;
        background: ${Y.Highlight};
      }
    `));class ih extends J.Tabs{}let iu=ih.compose({baseName:"tabs",baseClass:J.Tabs,template:J.tabsTemplate,styles:id}),ip=(e,t)=>(0,o1.css)`
    ${(0,J.display)("inline-block")} :host {
      font-family: ${ev};
      outline: none;
      user-select: none;
    }

    .control {
      box-sizing: border-box;
      position: relative;
      color: ${om};
      background: ${t6};
      border-radius: calc(${ek} * 1px);
      border: calc(${eS} * 1px) solid ${oa};
      height: calc(${o5} * 2px);
      font: inherit;
      font-size: ${ez};
      line-height: ${ej};
      padding: calc(${eF} * 2px + 1px);
      width: 100%;
      resize: none;
    }

    :host([aria-invalid='true']) .control {
      border-color: ${oT};
    }

    .control:hover:enabled {
      background: ${t3};
      border-color: ${oi};
    }

    :host([aria-invalid='true']) .control:hover:enabled {
      border-color: ${oz};
    }

    .control:active:enabled {
      background: ${t4};
      border-color: ${ol};
    }

    :host([aria-invalid='true']) .control:active:enabled {
      border-color: ${oj};
    }

    .control:hover,
    .control:${J.focusVisible},
    .control:disabled,
    .control:active {
      outline: none;
    }

    :host(:focus-within) .control {
      border-color: ${tN};
      box-shadow: 0 0 0 calc((${eT} - ${eS}) * 1px)
        ${tN};
    }

    :host([aria-invalid='true']:focus-within) .control {
      border-color: ${oB};
      box-shadow: 0 0 0 calc((${eT} - ${eS}) * 1px)
        ${oB};
    }

    :host([appearance='filled']) .control {
      background: ${tQ};
    }

    :host([appearance='filled']:hover:not([disabled])) .control {
      background: ${t0};
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
      color: ${om};
      cursor: pointer;
      font-size: ${ez};
      line-height: ${ej};
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
      cursor: ${J.disabledCursor};
    }
    :host([disabled]) {
      opacity: ${eD};
    }
    :host([disabled]) .control {
      border-color: ${o$};
    }

    :host([cols]) {
      width: initial;
    }

    :host([rows]) .control {
      height: initial;
    }
  `.withBehaviors((0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
      :host([disabled]) {
        opacity: 1;
      }

      :host([aria-invalid='true']) .control {
        border-style: dashed;
      }
    `));class ig extends J.TextArea{constructor(){super(...arguments),this.appearance="outline"}}o7([o1.attr],ig.prototype,"appearance",void 0);let ib=ig.compose({baseName:"text-area",baseClass:J.TextArea,template:J.textAreaTemplate,styles:ip,shadowOptions:{delegatesFocus:!0}}),im=(e,t)=>(0,o1.css)`
  ${r4}

  .start,
    .end {
    display: flex;
  }
`;class iv extends J.TextField{constructor(){super(...arguments),this.appearance="outline"}}o7([o1.attr],iv.prototype,"appearance",void 0);let i$=iv.compose({baseName:"text-field",baseClass:J.TextField,template:J.textFieldTemplate,styles:im,shadowOptions:{delegatesFocus:!0}});var ix=o(83021),iy=o(49054);let ik=(e,t)=>(0,o1.css)`
    ${(0,J.display)("inline-flex")} :host {
      --toolbar-item-gap: calc(
        (var(--design-unit) + calc(var(--density) + 2)) * 1px
      );
      background-color: ${tz};
      border-radius: calc(${ek} * 1px);
      fill: currentcolor;
      padding: var(--toolbar-item-gap);
    }

    :host(${J.focusVisible}) {
      outline: calc(${eS} * 1px) solid ${tN};
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
  `.withBehaviors((0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
      :host(:${J.focusVisible}) {
        box-shadow: 0 0 0 calc(${eT} * 1px)
          ${Y.Highlight};
        color: ${Y.ButtonText};
        forced-color-adjust: none;
      }
    `)),iw=Object.freeze({[r0.Is.ArrowUp]:{[a6.t.vertical]:-1},[r0.Is.ArrowDown]:{[a6.t.vertical]:1},[r0.Is.ArrowLeft]:{[a6.t.horizontal]:{[K.O.ltr]:-1,[K.O.rtl]:1}},[r0.Is.ArrowRight]:{[a6.t.horizontal]:{[K.O.ltr]:1,[K.O.rtl]:-1}}});class iF extends J.FoundationElement{constructor(){super(...arguments),this._activeIndex=0,this.direction=K.O.ltr,this.orientation=a6.t.horizontal}get activeIndex(){return o1.Observable.track(this,"activeIndex"),this._activeIndex}set activeIndex(e){this.$fastController.isConnected&&(this._activeIndex=(0,ix.AB)(0,this.focusableElements.length-1,e),o1.Observable.notify(this,"activeIndex"))}slottedItemsChanged(){this.$fastController.isConnected&&this.reduceFocusableElements()}mouseDownHandler(e){var t;let o=null==(t=this.focusableElements)?void 0:t.findIndex(t=>t.contains(e.target));return o>-1&&this.activeIndex!==o&&this.setFocusedElement(o),!0}childItemsChanged(e,t){this.$fastController.isConnected&&this.reduceFocusableElements()}connectedCallback(){super.connectedCallback(),this.direction=(0,J.getDirection)(this)}focusinHandler(e){let t=e.relatedTarget;!t||this.contains(t)||this.setFocusedElement()}getDirectionalIncrementer(e){var t,o,r,a,i;return null!=(i=null!=(r=null==(o=null==(t=iw[e])?void 0:t[this.orientation])?void 0:o[this.direction])?r:null==(a=iw[e])?void 0:a[this.orientation])?i:0}keydownHandler(e){let t=e.key;if(!(t in r0.Is)||e.defaultPrevented||e.shiftKey)return!0;let o=this.getDirectionalIncrementer(t);if(!o)return!e.target.closest("[role=radiogroup]");let r=this.activeIndex+o;return this.focusableElements[r]&&e.preventDefault(),this.setFocusedElement(r),!0}get allSlottedItems(){return[...this.start.assignedElements(),...this.slottedItems,...this.end.assignedElements()]}reduceFocusableElements(){var e;let t=null==(e=this.focusableElements)?void 0:e[this.activeIndex];this.focusableElements=this.allSlottedItems.reduce(iF.reduceFocusableItems,[]);let o=this.focusableElements.indexOf(t);this.activeIndex=Math.max(0,o),this.setFocusableElements()}setFocusedElement(e=this.activeIndex){this.activeIndex=e,this.setFocusableElements(),this.focusableElements[this.activeIndex]&&this.contains(document.activeElement)&&this.focusableElements[this.activeIndex].focus()}static reduceFocusableItems(e,t){var o,r,a,i;let l="radio"===t.getAttribute("role"),s=null==(r=null==(o=t.$fastController)?void 0:o.definition.shadowOptions)?void 0:r.delegatesFocus,n=Array.from(null!=(i=null==(a=t.shadowRoot)?void 0:a.querySelectorAll("*"))?i:[]).some(e=>(0,iy.tp)(e));return!t.hasAttribute("disabled")&&!t.hasAttribute("hidden")&&((0,iy.tp)(t)||l||s||n)?(e.push(t),e):t.childElementCount?e.concat(Array.from(t.children).reduce(iF.reduceFocusableItems,[])):e}setFocusableElements(){this.$fastController.isConnected&&this.focusableElements.length>0&&this.focusableElements.forEach((e,t)=>{e.tabIndex=this.activeIndex===t?0:-1})}}o7([o1.observable],iF.prototype,"direction",void 0),o7([o1.attr],iF.prototype,"orientation",void 0),o7([o1.observable],iF.prototype,"slottedItems",void 0),o7([o1.observable],iF.prototype,"slottedLabel",void 0),o7([o1.observable],iF.prototype,"childItems",void 0);class iC{}o7([(0,o1.attr)({attribute:"aria-labelledby"})],iC.prototype,"ariaLabelledby",void 0),o7([(0,o1.attr)({attribute:"aria-label"})],iC.prototype,"ariaLabel",void 0),(0,J.applyMixins)(iC,J.ARIAGlobalStatesAndProperties),(0,J.applyMixins)(iF,J.StartEnd,iC);class iV extends iF{connectedCallback(){super.connectedCallback();let e=(0,J.composedParent)(this);e&&tz.setValueFor(this,t=>on.getValueFor(t).evaluate(t,tz.getValueFor(e)))}}let iD=iV.compose({baseName:"toolbar",baseClass:iF,template:J.toolbarTemplate,styles:ik,shadowOptions:{delegatesFocus:!0}}),iS=(e,t)=>{let o=e.tagFor(J.AnchoredRegion);return(0,o1.css)`
    :host {
      contain: size;
      overflow: visible;
      height: 0;
      width: 0;
    }

    .tooltip {
      box-sizing: border-box;
      border-radius: calc(${ek} * 1px);
      border: calc(${eS} * 1px) solid ${oh};
      box-shadow: 0 0 0 1px ${oh} inset;
      background: ${tQ};
      color: ${om};
      padding: 4px;
      height: fit-content;
      width: fit-content;
      font-family: ${ev};
      font-size: ${ez};
      line-height: ${ej};
      white-space: nowrap;
      /* TODO: a mechanism to manage z-index across components
                    https://github.com/microsoft/fast/issues/3813 */
      z-index: 10000;
    }

    ${o} {
      display: flex;
      justify-content: center;
      align-items: center;
      overflow: visible;
      flex-direction: row;
    }

    ${o}.right,
    ${o}.left {
      flex-direction: column;
    }

    ${o}.top .tooltip {
      margin-bottom: 4px;
    }

    ${o}.bottom .tooltip {
      margin-top: 4px;
    }

    ${o}.left .tooltip {
      margin-right: 4px;
    }

    ${o}.right .tooltip {
      margin-left: 4px;
    }

    ${o}.top.left .tooltip,
            ${o}.top.right .tooltip {
      margin-bottom: 0px;
    }

    ${o}.bottom.left .tooltip,
            ${o}.bottom.right .tooltip {
      margin-top: 0px;
    }

    ${o}.top.left .tooltip,
            ${o}.bottom.left .tooltip {
      margin-right: 0px;
    }

    ${o}.top.right .tooltip,
            ${o}.bottom.right .tooltip {
      margin-left: 0px;
    }
  `.withBehaviors((0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
      :host([disabled]) {
        opacity: 1;
      }
    `))};class iT extends J.Tooltip{}let iz=iT.compose({baseName:"tooltip",baseClass:J.Tooltip,template:J.tooltipTemplate,styles:iS}),ij=(0,o1.cssPartial)`(((${e$} + ${ew}) * 0.5 + 2) * ${eF})`,iB=(0,o1.css)`
  .expand-collapse-glyph {
    transform: rotate(0deg);
  }
  :host(.nested) .expand-collapse-button {
    left: var(
      --expand-collapse-button-nested-width,
      calc(
        (
            ${ij} +
              ((${e$} + ${ew}) * 1.25)
          ) * -1px
      )
    );
  }
  :host([selected])::after {
    left: calc(${eT} * 1px);
  }
  :host([expanded]) > .positioning-region .expand-collapse-glyph {
    transform: rotate(90deg);
  }
`,iL=(0,o1.css)`
  .expand-collapse-glyph {
    transform: rotate(180deg);
  }
  :host(.nested) .expand-collapse-button {
    right: var(
      --expand-collapse-button-nested-width,
      calc(
        (
            ${ij} +
              ((${e$} + ${ew}) * 1.25)
          ) * -1px
      )
    );
  }
  :host([selected])::after {
    right: calc(${eT} * 1px);
  }
  :host([expanded]) > .positioning-region .expand-collapse-glyph {
    transform: rotate(90deg);
  }
`,iH=J.DesignToken.create("tree-item-expand-collapse-hover").withDefault(e=>{let t=t8.getValueFor(e);return t.evaluate(e,t.evaluate(e).hover).hover}),iN=J.DesignToken.create("tree-item-expand-collapse-selected-hover").withDefault(e=>{let t=tK.getValueFor(e);return t8.getValueFor(e).evaluate(e,t.evaluate(e).rest).hover}),iO=(e,t)=>(0,o1.css)`
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

    ${(0,J.display)("block")} :host {
      contain: content;
      position: relative;
      outline: none;
      color: ${om};
      background: ${t7};
      cursor: pointer;
      font-family: ${ev};
      --tree-item-nested-width: 0;
    }

    :host(:focus) > .positioning-region {
      outline: none;
    }

    :host(:focus) .content-region {
      outline: none;
    }

    :host(:${J.focusVisible}) .positioning-region {
      border-color: ${tN};
      box-shadow: 0 0 0 calc((${eT} - ${eS}) * 1px)
        ${tN} inset;
      color: ${om};
    }

    .positioning-region {
      display: flex;
      position: relative;
      box-sizing: border-box;
      background: ${t7};
      border: transparent calc(${eS} * 1px) solid;
      border-radius: calc(${ek} * 1px);
      height: calc((${o5} + 1) * 1px);
    }

    .positioning-region::before {
      content: '';
      display: block;
      width: var(--tree-item-nested-width);
      flex-shrink: 0;
    }

    :host(:not([disabled])) .positioning-region:hover {
      background: ${oe};
    }

    :host(:not([disabled])) .positioning-region:active {
      background: ${ot};
    }

    .content-region {
      display: inline-flex;
      align-items: center;
      white-space: nowrap;
      width: 100%;
      min-width: 0;
      height: calc(${o5} * 1px);
      margin-inline-start: calc(${eF} * 2px + 8px);
      font-size: ${ez};
      line-height: ${ej};
      font-weight: 400;
    }

    .items {
      /* TODO: adaptive typography https://github.com/microsoft/fast/issues/2432 */
      font-size: calc(1em + (${eF} + 16) * 1px);
    }

    .expand-collapse-button {
      background: none;
      border: none;
      outline: none;
      /* TODO: adaptive typography https://github.com/microsoft/fast/issues/2432 */
      width: calc(${ij} * 1px);
      height: calc(${ij} * 1px);
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
      width: calc((16 + ${ew}) * 1px);
      height: calc((16 + ${ew}) * 1px);
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
      width: ${eO};
      height: ${eO};
      */
    }

    .start {
      /* TODO: horizontalSpacing https://github.com/microsoft/fast/issues/2766 */
      margin-inline-end: calc(${eF} * 2px + 2px);
    }

    .end {
      /* TODO: horizontalSpacing https://github.com/microsoft/fast/issues/2766 */
      margin-inline-start: calc(${eF} * 2px + 2px);
    }

    :host([expanded]) > .items {
      animation: treeItemLoading ease-in 10ms;
      animation-iteration-count: 1;
      animation-fill-mode: forwards;
    }

    :host([disabled]) .content-region {
      opacity: ${eD};
      cursor: ${J.disabledCursor};
    }

    :host(.nested) .content-region {
      position: relative;
      /* Add left margin to collapse button size */
      margin-inline-start: calc(
        (
            ${ij} +
              ((${e$} + ${ew}) * 1.25)
          ) * 1px
      );
    }

    :host(.nested) .expand-collapse-button {
      position: absolute;
    }

    :host(.nested:not([disabled])) .expand-collapse-button:hover {
      background: ${iH};
    }

    :host([selected]) .positioning-region {
      background: ${tQ};
    }

    :host([selected]:not([disabled])) .positioning-region:hover {
      background: ${t0};
    }

    :host([selected]:not([disabled])) .positioning-region:active {
      background: ${t1};
    }

    :host([selected]:not([disabled])) .expand-collapse-button:hover {
      background: ${iN};
    }

    :host([selected])::after {
      /* The background needs to be calculated based on the selected background state
         for this control. We currently have no way of changing that, so setting to
         accent-foreground-rest for the time being */
      background: ${tX};
      border-radius: calc(${ek} * 1px);
      content: '';
      display: block;
      position: absolute;
      top: calc((${o5} / 4) * 1px);
      width: 3px;
      height: calc((${o5} / 2) * 1px);
    }

    ::slotted(${e.tagFor(J.TreeItem)}) {
      --tree-item-nested-width: 1em;
      --expand-collapse-button-nested-width: calc(
        (
            ${ij} +
              ((${e$} + ${ew}) * 1.25)
          ) * -1px
      );
    }
  `.withBehaviors(new rg(iB,iL),(0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
      :host {
        forced-color-adjust: none;
        border-color: transparent;
        background: ${Y.Field};
        color: ${Y.FieldText};
      }
      :host .content-region .expand-collapse-glyph {
        fill: ${Y.FieldText};
      }
      :host .positioning-region:hover,
      :host([selected]) .positioning-region {
        background: ${Y.Highlight};
      }
      :host .positioning-region:hover .content-region,
      :host([selected]) .positioning-region .content-region {
        color: ${Y.HighlightText};
      }
      :host .positioning-region:hover .content-region .expand-collapse-glyph,
      :host .positioning-region:hover .content-region .start,
      :host .positioning-region:hover .content-region .end,
      :host([selected]) .content-region .expand-collapse-glyph,
      :host([selected]) .content-region .start,
      :host([selected]) .content-region .end {
        fill: ${Y.HighlightText};
      }
      :host([selected])::after {
        background: ${Y.Field};
      }
      :host(:${J.focusVisible}) .positioning-region {
        border-color: ${Y.FieldText};
        box-shadow: 0 0 0 2px inset ${Y.Field};
        color: ${Y.FieldText};
      }
      :host([disabled]) .content-region,
      :host([disabled]) .positioning-region:hover .content-region {
        opacity: 1;
        color: ${Y.GrayText};
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
        fill: ${Y.GrayText};
      }
      :host([disabled]) .positioning-region:hover {
        background: ${Y.Field};
      }
      .expand-collapse-glyph,
      .start,
      .end {
        fill: ${Y.FieldText};
      }
      :host(.nested) .expand-collapse-button:hover {
        background: ${Y.Field};
      }
      :host(.nested) .expand-collapse-button:hover .expand-collapse-glyph {
        fill: ${Y.FieldText};
      }
    `));class iI extends J.TreeItem{}let iR=iI.compose({baseName:"tree-item",baseClass:J.TreeItem,template:J.treeItemTemplate,styles:iO,expandCollapseGlyph:`
        <svg
            viewBox="0 0 16 16"
            xmlns="http://www.w3.org/2000/svg"
            class="expand-collapse-glyph"
        >
            <path
                d="M5.00001 12.3263C5.00124 12.5147 5.05566 12.699 5.15699 12.8578C5.25831 13.0167 5.40243 13.1437 5.57273 13.2242C5.74304 13.3047 5.9326 13.3354 6.11959 13.3128C6.30659 13.2902 6.4834 13.2152 6.62967 13.0965L10.8988 8.83532C11.0739 8.69473 11.2153 8.51658 11.3124 8.31402C11.4096 8.11146 11.46 7.88966 11.46 7.66499C11.46 7.44033 11.4096 7.21853 11.3124 7.01597C11.2153 6.81341 11.0739 6.63526 10.8988 6.49467L6.62967 2.22347C6.48274 2.10422 6.30501 2.02912 6.11712 2.00691C5.92923 1.9847 5.73889 2.01628 5.56823 2.09799C5.39757 2.17969 5.25358 2.30817 5.153 2.46849C5.05241 2.62882 4.99936 2.8144 5.00001 3.00369V12.3263Z"
            />
        </svg>
    `}),iP=(e,t)=>(0,o1.css)`
  ${(0,J.display)("flex")} :host {
    flex-direction: column;
    align-items: stretch;
    min-width: fit-content;
    font-size: 0;
  }

  :host:focus-visible {
    outline: none;
  }
`;class iA extends J.TreeView{handleClick(e){if(e.defaultPrevented)return;if(!(e.target instanceof Element))return!0;let t=e.target;for(;t&&!(0,J.isTreeItemElement)(t);)(t=t.parentElement)===this&&(t=null);t&&!t.disabled&&(t.selected=!0)}}let iM=iA.compose({baseName:"tree-view",baseClass:J.TreeView,template:J.treeViewTemplate,styles:iP}),iG=(e,t)=>(0,o1.css)`
  .region {
    z-index: 1000;
    overflow: hidden;
    display: flex;
    font-family: ${ev};
    font-size: ${ez};
  }

  .loaded {
    opacity: 1;
    pointer-events: none;
  }

  .loading-display,
  .no-options-display {
    background: ${tz};
    width: 100%;
    min-height: calc(${o5} * 1px);
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-items: center;
    padding: calc(${eF} * 1px);
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
`,iE=(e,t)=>(0,o1.css)`
    :host {
      background: ${tz};
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
      border-radius: calc(${ek} * 1px);
      padding: calc(${eF} * 1px) 0;
      border: calc(${eS} * 1px) solid transparent;
      ${rB}
    }

    .suggestions-available-alert {
      height: 0;
      opacity: 0;
      overflow: hidden;
    }
  `.withBehaviors((0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
      :host {
        background: ${Y.Canvas};
        border-color: ${Y.CanvasText};
      }
    `)),i_=(e,t)=>(0,o1.css)`
    :host {
      display: flex;
      align-items: center;
      justify-items: center;
      font-family: ${ev};
      border-radius: calc(${ek} * 1px);
      border: calc(${eT} * 1px) solid transparent;
      box-sizing: border-box;
      background: ${t7};
      color: ${om};
      cursor: pointer;
      fill: currentcolor;
      font-size: ${ez};
      min-height: calc(${o5} * 1px);
      line-height: ${ej};
      margin: 0 calc(${eF} * 1px);
      outline: none;
      overflow: hidden;
      padding: 0 calc(${eF} * 2.25px);
      user-select: none;
      white-space: nowrap;
    }

    :host(:${J.focusVisible}[role="listitem"]) {
      border-color: ${oh};
      background: ${oo};
    }

    :host(:hover) {
      background: ${oe};
    }

    :host(:active) {
      background: ${ot};
    }

    :host([aria-selected='true']) {
      background: ${tB};
      color: ${tR};
    }

    :host([aria-selected='true']:hover) {
      background: ${tL};
      color: ${tP};
    }

    :host([aria-selected='true']:active) {
      background: ${tH};
      color: ${tA};
    }
  `.withBehaviors((0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
      :host {
        border-color: transparent;
        forced-color-adjust: none;
        color: ${Y.ButtonText};
        fill: currentcolor;
      }

      :host(:not([aria-selected='true']):hover),
      :host([aria-selected='true']) {
        background: ${Y.Highlight};
        color: ${Y.HighlightText};
      }

      :host([disabled]),
      :host([disabled]:not([aria-selected='true']):hover) {
        background: ${Y.Canvas};
        color: ${Y.GrayText};
        fill: currentcolor;
        opacity: 1;
      }
    `)),iq=(e,t)=>(0,o1.css)`
    :host {
      display: flex;
      align-items: center;
      justify-items: center;
      font-family: ${ev};
      border-radius: calc(${ek} * 1px);
      border: calc(${eT} * 1px) solid transparent;
      box-sizing: border-box;
      background: ${t7};
      color: ${om};
      cursor: pointer;
      fill: currentcolor;
      font-size: ${ez};
      height: calc(${o5} * 1px);
      line-height: ${ej};
      outline: none;
      overflow: hidden;
      padding: 0 calc(${eF} * 2.25px);
      user-select: none;
      white-space: nowrap;
    }

    :host(:hover) {
      background: ${oe};
    }

    :host(:active) {
      background: ${ot};
    }

    :host(:${J.focusVisible}) {
      background: ${oo};
      border-color: ${oh};
    }

    :host([aria-selected='true']) {
      background: ${tB};
      color: ${tA};
    }
  `.withBehaviors((0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
      :host {
        border-color: transparent;
        forced-color-adjust: none;
        color: ${Y.ButtonText};
        fill: currentcolor;
      }

      :host(:not([aria-selected='true']):hover),
      :host([aria-selected='true']) {
        background: ${Y.Highlight};
        color: ${Y.HighlightText};
      }

      :host([disabled]),
      :host([disabled]:not([aria-selected='true']):hover) {
        background: ${Y.Canvas};
        color: ${Y.GrayText};
        fill: currentcolor;
        opacity: 1;
      }
    `));class iW extends J.Picker{}let iU=iW.compose({baseName:"draft-picker",baseClass:J.Picker,template:J.pickerTemplate,styles:iG,shadowOptions:{}});class iX extends J.PickerMenu{connectedCallback(){tz.setValueFor(this,ty),super.connectedCallback()}}let iZ=iX.compose({baseName:"draft-picker-menu",baseClass:J.PickerMenu,template:J.pickerMenuTemplate,styles:iE});class iY extends J.PickerMenuOption{}let iJ=iY.compose({baseName:"draft-picker-menu-option",baseClass:J.PickerMenuOption,template:J.pickerMenuOptionTemplate,styles:i_});class iK extends J.PickerList{}let iQ=iK.compose({baseName:"draft-picker-list",baseClass:J.PickerList,template:J.pickerListTemplate,styles:(e,t)=>(0,o1.css)`
        :host {
            display: flex;
            flex-direction: row;
            column-gap: calc(${eF} * 1px);
            row-gap: calc(${eF} * 1px);
            flex-wrap: wrap;
        }

        ::slotted([role="combobox"]) {
            min-width: 260px;
            width: auto;
            box-sizing: border-box;
            color: ${om};
            background: ${t6};
            border-radius: calc(${ek} * 1px);
            border: calc(${eS} * 1px) solid ${tB};
            height: calc(${o5} * 1px);
            font-family: ${ev};
            outline: none;
            user-select: none;
            font-size: ${ez};
            line-height: ${ej};
            padding: 0 calc(${eF} * 2px + 1px);
        }

        ::slotted([role="combobox"]:active) { {
            background: ${t3};
            border-color: ${tH};
        }

        ::slotted([role="combobox"]:focus-within) {
            border-color: ${oh};
            box-shadow: 0 0 0 1px ${oh} inset;
        }
    `.withBehaviors((0,J.forcedColorsStylesheetBehavior)((0,o1.css)`
      ::slotted([role='combobox']:active) {
        background: ${Y.Field};
        border-color: ${Y.Highlight};
      }
      ::slotted([role='combobox']:focus-within) {
        border-color: ${Y.Highlight};
        box-shadow: 0 0 0 1px ${Y.Highlight} inset;
      }
      ::slotted(input:placeholder) {
        color: ${Y.GrayText};
      }
    `))});class i0 extends J.PickerListItem{}let i1=i0.compose({baseName:"draft-picker-list-item",baseClass:J.PickerListItem,template:J.pickerListItemTemplate,styles:iq}),i2={jpAccordion:o8,jpAccordionItem:o4,jpAnchor:rd,jpAnchoredRegion:rp,jpAvatar:r$,jpBadge:rk,jpBreadcrumb:rC,jpBreadcrumbItem:rS,jpButton:rj,jpCard:rN,jpCheckbox:rP,jpCombobox:r_,jpDataGrid:rQ,jpDataGridCell:rZ,jpDataGridRow:rJ,jpDateField:r7,jpDesignSystemProvider:al,jpDialog:ac,jpDisclosure:au,jpDivider:ab,jpListbox:am,jpMenu:ax,jpMenuItem:aw,jpNumberField:aV,jpOption:aT,jpPicker:iU,jpPickerList:iQ,jpPickerListItem:i1,jpPickerMenu:iZ,jpPickerMenuOption:iJ,jpProgress:aB,jpProgressRing:aN,jpRadio:aP,jpRadioGroup:aG,jpSearch:aU,jpSelect:aZ,jpSkeleton:aK,jpSlider:a5,jpSliderLabel:a7,jpSwitch:io,jpTab:ic,jpTabPanel:ii,jpTabs:iu,jpTextArea:ib,jpTextField:i$,jpToolbar:iD,jpTooltip:iz,jpTreeItem:iR,jpTreeView:iM,register(e,...t){if(e)for(let o in this)"register"!==o&&this[o]().register(e,...t)}};function i5(e){return J.DesignSystem.getOrCreate(e).withPrefix("jp")}}}]);