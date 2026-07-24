"use strict";(self.rspackChunk_JUPYTERLAB_CORE_OUTPUT=self.rspackChunk_JUPYTERLAB_CORE_OUTPUT||[]).push([[3295],{23850(e,t,a){function i(e,t){e.accDescr&&t.setAccDescription?.(e.accDescr),e.accTitle&&t.setAccTitle?.(e.accTitle),e.title&&t.setDiagramTitle?.(e.title)}a.d(t,{S:()=>i}),(0,a(87724).K2)(i,"populateCommonDb")},6834(e,t,a){a.d(t,{diagram:()=>v});var i=a(74571),l=a(23850),r=a(2961),s=a(74725),n=a(87724),o=a(11150),c=a(91864),p=s.UI.pie,d={sections:new Map,showData:!1,config:p},u=d.sections,g=d.showData,h=structuredClone(p),f=(0,n.K2)(()=>structuredClone(h),"getConfig"),m=(0,n.K2)(()=>{u=new Map,g=d.showData,(0,s.IU)()},"clear"),x=(0,n.K2)(({label:e,value:t})=>{if(t<0)throw Error(`"${e}" has invalid value: ${t}. Negative values are not allowed in pie charts. All slice values must be >= 0.`);u.has(e)||(u.set(e,t),n.Rm.debug(`added new section: ${e}, with value: ${t}`))},"addSection"),T=(0,n.K2)(()=>u,"getSections"),w=(0,n.K2)(e=>{g=e},"setShowData"),S=(0,n.K2)(()=>g,"getShowData"),$={getConfig:f,clear:m,setDiagramTitle:s.ke,getDiagramTitle:s.ab,setAccTitle:s.SV,getAccTitle:s.iN,setAccDescription:s.EI,getAccDescription:s.m7,addSection:x,getSections:T,setShowData:w,getShowData:S},C=(0,n.K2)((e,t)=>{(0,l.S)(e,t),t.setShowData(e.showData),e.sections.map(t.addSection)},"populateDb"),D={parse:(0,n.K2)(async e=>{let t=await (0,o.qg)("pie",e);n.Rm.debug(t),C(t,$)},"parse")},y=(0,n.K2)(e=>`
  .pieCircle{
    stroke: ${e.pieStrokeColor};
    stroke-width : ${e.pieStrokeWidth};
    opacity : ${e.pieOpacity};
  }
  .pieOuterCircle{
    stroke: ${e.pieOuterStrokeColor};
    stroke-width: ${e.pieOuterStrokeWidth};
    fill: none;
  }
  .pieTitleText {
    text-anchor: middle;
    font-size: ${e.pieTitleTextSize};
    fill: ${e.pieTitleTextColor};
    font-family: ${e.fontFamily};
  }
  .slice {
    font-family: ${e.fontFamily};
    fill: ${e.pieSectionTextColor};
    font-size:${e.pieSectionTextSize};
    // fill: white;
  }
  .legend text {
    fill: ${e.pieLegendTextColor};
    font-family: ${e.fontFamily};
    font-size: ${e.pieLegendTextSize};
  }
`,"getStyles"),k=(0,n.K2)(e=>{let t=[...e.values()].reduce((e,t)=>e+t,0),a=[...e.entries()].map(([e,t])=>({label:e,value:t})).filter(e=>e.value/t*100>=1);return(0,c.rLf)().value(e=>e.value).sort(null)(a)},"createPieArcs"),v={parser:D,db:$,renderer:{draw:(0,n.K2)((e,t,a,l)=>{n.Rm.debug("rendering pie chart\n"+e);let o=l.db,p=(0,s.D7)(),d=(0,r.$t)(o.getConfig(),p.pie),u=(0,i.D)(t),g=u.append("g");g.attr("transform","translate(225,225)");let{themeVariables:h}=p,[f]=(0,r.I5)(h.pieOuterStrokeWidth);f??=2;let m=d.textPosition,x=(0,c.JLW)().innerRadius(0).outerRadius(185),T=(0,c.JLW)().innerRadius(185*m).outerRadius(185*m);g.append("circle").attr("cx",0).attr("cy",0).attr("r",185+f/2).attr("class","pieOuterCircle");let w=o.getSections(),S=k(w),$=[h.pie1,h.pie2,h.pie3,h.pie4,h.pie5,h.pie6,h.pie7,h.pie8,h.pie9,h.pie10,h.pie11,h.pie12],C=0;w.forEach(e=>{C+=e});let D=S.filter(e=>"0"!==(e.data.value/C*100).toFixed(0)),y=(0,c.UMr)($).domain([...w.keys()]);g.selectAll("mySlices").data(D).enter().append("path").attr("d",x).attr("fill",e=>y(e.data.label)).attr("class","pieCircle"),g.selectAll("mySlices").data(D).enter().append("text").text(e=>(e.data.value/C*100).toFixed(0)+"%").attr("transform",e=>"translate("+T.centroid(e)+")").style("text-anchor","middle").attr("class","slice");let v=g.append("text").text(o.getDiagramTitle()).attr("x",0).attr("y",-200).attr("class","pieTitleText"),A=[...w.entries()].map(([e,t])=>({label:e,value:t})),b=g.selectAll(".legend").data(A).enter().append("g").attr("class","legend").attr("transform",(e,t)=>"translate(216,"+(22*t-22*A.length/2)+")");b.append("rect").attr("width",18).attr("height",18).style("fill",e=>y(e.label)).style("stroke",e=>y(e.label)),b.append("text").attr("x",22).attr("y",14).text(e=>o.getShowData()?`${e.label} [${e.value}]`:e.label);let R=Math.max(...b.selectAll("text").nodes().map(e=>e?.getBoundingClientRect().width??0)),K=v.node()?.getBoundingClientRect().width??0,O=Math.min(0,225-K/2),U=Math.max(512+R,225+K/2)-O;u.attr("viewBox",`${O} 0 ${U} 450`),(0,s.a$)(u,450,U,d.useMaxWidth)},"draw")},styles:y}}}]);