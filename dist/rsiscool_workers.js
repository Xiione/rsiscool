
var rsiscool = (() => {
  var _scriptName = import.meta.url;
  
  return (
function(moduleArg = {}) {
  var moduleRtn;

var k=moduleArg,aa,ba,ca=new Promise((a,b)=>{aa=a;ba=b}),da=Object.assign({},k),ea="./this.program",p="",fa,ha;p=self.location.href;_scriptName&&(p=_scriptName);p.startsWith("blob:")?p="":p=p.substr(0,p.replace(/[?#].*/,"").lastIndexOf("/")+1);ha=a=>{var b=new XMLHttpRequest;b.open("GET",a,!1);b.responseType="arraybuffer";b.send(null);return new Uint8Array(b.response)};fa=async a=>{a=await fetch(a,{credentials:"same-origin"});if(a.ok)return a.arrayBuffer();throw Error(a.status+" : "+a.url);};
var ja=k.print||console.log.bind(console),q=k.printErr||console.error.bind(console);Object.assign(k,da);da=null;k.thisProgram&&(ea=k.thisProgram);var ka=k.wasmBinary,la,ma=!1,na,r,v,B,oa,C,D,pa,qa,ra=[],sa=[],ta=[];function ua(){var a=k.preRun.shift();ra.unshift(a)}var E=0,F=null;function G(a){k.onAbort?.(a);a="Aborted("+a+")";q(a);ma=!0;a=new WebAssembly.RuntimeError(a+". Build with -sASSERTIONS for more info.");ba(a);throw a;}var va=a=>a.startsWith("data:application/octet-stream;base64,"),wa;
async function xa(a){if(!ka)try{var b=await fa(a);return new Uint8Array(b)}catch{}if(a==wa&&ka)a=new Uint8Array(ka);else if(ha)a=ha(a);else throw"both async and sync fetching of the wasm failed";return a}async function ya(a,b){try{var c=await xa(a);return await WebAssembly.instantiate(c,b)}catch(d){q(`failed to asynchronously prepare wasm: ${d}`),G(d)}}
async function za(a){var b=wa;if(!ka&&"function"==typeof WebAssembly.instantiateStreaming&&!va(b)&&"function"==typeof fetch)try{var c=fetch(b,{credentials:"same-origin"});return await WebAssembly.instantiateStreaming(c,a)}catch(d){q(`wasm streaming compile failed: ${d}`),q("falling back to ArrayBuffer instantiation")}return ya(b,a)}var H,Aa;class Ba{name="ExitStatus";constructor(a){this.message=`Program terminated with exit(${a})`;this.status=a}}
var Ca=a=>{for(;0<a.length;)a.shift()(k)},Da=k.noExitRuntime||!0;class Ea{constructor(a){this.Z=a-24}}var Fa=0,Ga=0,Ha={},Ia=a=>{for(;a.length;){var b=a.pop();a.pop()(b)}};function I(a){return this.fromWireType(D[a>>2])}
var J={},K={},Ja={},L,N=(a,b,c)=>{function d(h){h=c(h);if(h.length!==a.length)throw new L("Mismatched type converter count");for(var m=0;m<a.length;++m)M(a[m],h[m])}a.forEach(h=>Ja[h]=b);var e=Array(b.length),f=[],g=0;b.forEach((h,m)=>{K.hasOwnProperty(h)?e[m]=K[h]:(f.push(h),J.hasOwnProperty(h)||(J[h]=[]),J[h].push(()=>{e[m]=K[h];++g;g===f.length&&d(e)}))});0===f.length&&d(e)},Ka,P=a=>{for(var b="";v[a];)b+=Ka[v[a++]];return b},Q,La=a=>{throw new Q(a);};
function Ma(a,b,c={}){var d=b.name;if(!a)throw new Q(`type "${d}" must have a positive integer typeid pointer`);if(K.hasOwnProperty(a)){if(c.sb)return;throw new Q(`Cannot register type '${d}' twice`);}K[a]=b;delete Ja[a];J.hasOwnProperty(a)&&(b=J[a],delete J[a],b.forEach(e=>e()))}function M(a,b,c={}){return Ma(a,b,c)}
var Na=a=>{throw new Q(a.V.aa.$.name+" instance already deleted");},Oa=!1,Pa=()=>{},Qa=(a,b,c)=>{if(b===c)return a;if(void 0===c.fa)return null;a=Qa(a,b,c.fa);return null===a?null:c.lb(a)},Ra={},Sa={},Ta=(a,b)=>{if(void 0===b)throw new Q("ptr should not be undefined");for(;a.fa;)b=a.Ea(b),a=a.fa;return Sa[b]},Va=(a,b)=>{if(!b.aa||!b.Z)throw new L("makeClassHandle requires ptr and ptrType");if(!!b.ga!==!!b.da)throw new L("Both smartPtrType and smartPtr must be specified");b.count={value:1};return Ua(Object.create(a,
{V:{value:b,writable:!0}}))},Ua=a=>{if("undefined"===typeof FinalizationRegistry)return Ua=b=>b,a;Oa=new FinalizationRegistry(b=>{b=b.V;--b.count.value;0===b.count.value&&(b.da?b.ga.na(b.da):b.aa.$.na(b.Z))});Ua=b=>{var c=b.V;c.da&&Oa.register(b,{V:c},b);return b};Pa=b=>{Oa.unregister(b)};return Ua(a)},Wa=[];function Xa(){}
var Ya=(a,b)=>Object.defineProperty(b,"name",{value:a}),Za=(a,b,c)=>{if(void 0===a[b].ca){var d=a[b];a[b]=function(...e){if(!a[b].ca.hasOwnProperty(e.length))throw new Q(`Function '${c}' called with an invalid number of arguments (${e.length}) - expects one of (${a[b].ca})!`);return a[b].ca[e.length].apply(this,e)};a[b].ca=[];a[b].ca[d.ya]=d}},$a=(a,b,c)=>{if(k.hasOwnProperty(a)){if(void 0===c||void 0!==k[a].ca&&void 0!==k[a].ca[c])throw new Q(`Cannot register public name '${a}' twice`);Za(k,a,a);
if(k[a].ca.hasOwnProperty(c))throw new Q(`Cannot register multiple overloads of a function with the same number of arguments (${c})!`);k[a].ca[c]=b}else k[a]=b,k[a].ya=c},ab=a=>{a=a.replace(/[^a-zA-Z0-9_]/g,"$");var b=a.charCodeAt(0);return 48<=b&&57>=b?`_${a}`:a};function cb(a,b,c,d,e,f,g,h){this.name=a;this.constructor=b;this.va=c;this.na=d;this.fa=e;this.nb=f;this.Ea=g;this.lb=h;this.wb=[]}
var db=(a,b,c)=>{for(;b!==c;){if(!b.Ea)throw new Q(`Expected null or instance of ${c.name}, got an instance of ${b.name}`);a=b.Ea(a);b=b.fa}return a};function eb(a,b){if(null===b){if(this.Pa)throw new Q(`null is not a valid ${this.name}`);return 0}if(!b.V)throw new Q(`Cannot pass "${fb(b)}" as a ${this.name}`);if(!b.V.Z)throw new Q(`Cannot pass deleted object as a pointer of type ${this.name}`);return db(b.V.Z,b.V.aa.$,this.$)}
function gb(a,b){if(null===b){if(this.Pa)throw new Q(`null is not a valid ${this.name}`);if(this.Ia){var c=this.Ra();null!==a&&a.push(this.na,c);return c}return 0}if(!b||!b.V)throw new Q(`Cannot pass "${fb(b)}" as a ${this.name}`);if(!b.V.Z)throw new Q(`Cannot pass deleted object as a pointer of type ${this.name}`);if(!this.Ha&&b.V.aa.Ha)throw new Q(`Cannot convert argument of type ${b.V.ga?b.V.ga.name:b.V.aa.name} to parameter type ${this.name}`);c=db(b.V.Z,b.V.aa.$,this.$);if(this.Ia){if(void 0===
b.V.da)throw new Q("Passing raw pointer to smart pointer is illegal");switch(this.Bb){case 0:if(b.V.ga===this)c=b.V.da;else throw new Q(`Cannot convert argument of type ${b.V.ga?b.V.ga.name:b.V.aa.name} to parameter type ${this.name}`);break;case 1:c=b.V.da;break;case 2:if(b.V.ga===this)c=b.V.da;else{var d=b.clone();c=this.xb(c,hb(()=>d["delete"]()));null!==a&&a.push(this.na,c)}break;default:throw new Q("Unsupporting sharing policy");}}return c}
function ib(a,b){if(null===b){if(this.Pa)throw new Q(`null is not a valid ${this.name}`);return 0}if(!b.V)throw new Q(`Cannot pass "${fb(b)}" as a ${this.name}`);if(!b.V.Z)throw new Q(`Cannot pass deleted object as a pointer of type ${this.name}`);if(b.V.aa.Ha)throw new Q(`Cannot convert argument of type ${b.V.aa.name} to parameter type ${this.name}`);return db(b.V.Z,b.V.aa.$,this.$)}
function jb(a,b,c,d,e,f,g,h,m,l,n){this.name=a;this.$=b;this.Pa=c;this.Ha=d;this.Ia=e;this.vb=f;this.Bb=g;this.fb=h;this.Ra=m;this.xb=l;this.na=n;e||void 0!==b.fa?this.toWireType=gb:(this.toWireType=d?eb:ib,this.ja=null)}
var kb=(a,b,c)=>{if(!k.hasOwnProperty(a))throw new L("Replacing nonexistent public symbol");void 0!==k[a].ca&&void 0!==c?k[a].ca[c]=b:(k[a]=b,k[a].ya=c)},lb=[],mb,nb=a=>{var b=lb[a];b||(a>=lb.length&&(lb.length=a+1),lb[a]=b=mb.get(a));return b},ob=(a,b,c=[])=>{a.includes("j")?(a=a.replace(/p/g,"i"),b=(0,k["dynCall_"+a])(b,...c)):b=nb(b)(...c);return b},pb=(a,b)=>(...c)=>ob(a,b,c),R=(a,b)=>{a=P(a);var c=a.includes("j")?pb(a,b):nb(b);if("function"!=typeof c)throw new Q(`unknown function pointer with signature ${a}: ${b}`);
return c},qb,sb=a=>{a=rb(a);var b=P(a);S(a);return b},tb=(a,b)=>{function c(f){e[f]||K[f]||(Ja[f]?Ja[f].forEach(c):(d.push(f),e[f]=!0))}var d=[],e={};b.forEach(c);throw new qb(`${a}: `+d.map(sb).join([", "]));},ub=(a,b)=>{for(var c=[],d=0;d<a;d++)c.push(D[b+4*d>>2]);return c};function vb(a){for(var b=1;b<a.length;++b)if(null!==a[b]&&void 0===a[b].ja)return!0;return!1}
function wb(a){var b=Function;if(!(b instanceof Function))throw new TypeError(`new_ called with constructor type ${typeof b} which is not a function`);var c=Ya(b.name||"unknownFunctionName",function(){});c.prototype=b.prototype;c=new c;a=b.apply(c,a);return a instanceof Object?a:c}
function xb(a,b,c,d,e,f){var g=b.length;if(2>g)throw new Q("argTypes array size mismatch! Must at least get return value and 'this' types!");var h=null!==b[1]&&null!==c,m=vb(b);c="void"!==b[0].name;d=[a,La,d,e,Ia,b[0],b[1]];for(e=0;e<g-2;++e)d.push(b[e+2]);if(!m)for(e=h?1:2;e<b.length;++e)null!==b[e].ja&&d.push(b[e].ja);m=vb(b);e=b.length-2;var l=[],n=["fn"];h&&n.push("thisWired");for(g=0;g<e;++g)l.push(`arg${g}`),n.push(`arg${g}Wired`);l=l.join(",");n=n.join(",");l=`return function (${l}) {\n`;m&&
(l+="var destructors = [];\n");var t=m?"destructors":"null",u="humanName throwBindingError invoker fn runDestructors retType classParam".split(" ");h&&(l+=`var thisWired = classParam['toWireType'](${t}, this);\n`);for(g=0;g<e;++g)l+=`var arg${g}Wired = argType${g}['toWireType'](${t}, arg${g});\n`,u.push(`argType${g}`);l+=(c||f?"var rv = ":"")+`invoker(${n});\n`;if(m)l+="runDestructors(destructors);\n";else for(g=h?1:2;g<b.length;++g)f=1===g?"thisWired":"arg"+(g-2)+"Wired",null!==b[g].ja&&(l+=`${f}_dtor(${f});\n`,
u.push(`${f}_dtor`));c&&(l+="var ret = retType['fromWireType'](rv);\nreturn ret;\n");let [x,w]=[u,l+"}\n"];x.push(w);b=wb(x)(...d);return Ya(a,b)}
var yb=a=>{a=a.trim();const b=a.indexOf("(");return-1!==b?a.substr(0,b):a},zb=[],T=[],Ab=a=>{9<a&&0===--T[a+1]&&(T[a]=void 0,zb.push(a))},Bb=a=>{if(!a)throw new Q("Cannot use deleted val. handle = "+a);return T[a]},hb=a=>{switch(a){case void 0:return 2;case null:return 4;case !0:return 6;case !1:return 8;default:const b=zb.pop()||T.length;T[b]=a;T[b+1]=1;return b}},Cb={name:"emscripten::val",fromWireType:a=>{var b=Bb(a);Ab(a);return b},toWireType:(a,b)=>hb(b),pa:8,readValueFromPointer:I,ja:null},
fb=a=>{if(null===a)return"null";var b=typeof a;return"object"===b||"array"===b||"function"===b?a.toString():""+a},Db=(a,b)=>{switch(b){case 4:return function(c){return this.fromWireType(pa[c>>2])};case 8:return function(c){return this.fromWireType(qa[c>>3])};default:throw new TypeError(`invalid float width (${b}): ${a}`);}},Eb=(a,b,c)=>{switch(b){case 1:return c?d=>r[d]:d=>v[d];case 2:return c?d=>B[d>>1]:d=>oa[d>>1];case 4:return c?d=>C[d>>2]:d=>D[d>>2];default:throw new TypeError(`invalid integer width (${b}): ${a}`);
}},Fb=Object.assign({optional:!0},Cb),Gb=(a,b,c)=>{var d=v;if(0<c){c=b+c-1;for(var e=0;e<a.length;++e){var f=a.charCodeAt(e);if(55296<=f&&57343>=f){var g=a.charCodeAt(++e);f=65536+((f&1023)<<10)|g&1023}if(127>=f){if(b>=c)break;d[b++]=f}else{if(2047>=f){if(b+1>=c)break;d[b++]=192|f>>6}else{if(65535>=f){if(b+2>=c)break;d[b++]=224|f>>12}else{if(b+3>=c)break;d[b++]=240|f>>18;d[b++]=128|f>>12&63}d[b++]=128|f>>6&63}d[b++]=128|f&63}}d[b]=0}},Hb="undefined"!=typeof TextDecoder?new TextDecoder:void 0,Ib=(a,
b=0,c=NaN)=>{var d=b+c;for(c=b;a[c]&&!(c>=d);)++c;if(16<c-b&&a.buffer&&Hb)return Hb.decode(a.subarray(b,c));for(d="";b<c;){var e=a[b++];if(e&128){var f=a[b++]&63;if(192==(e&224))d+=String.fromCharCode((e&31)<<6|f);else{var g=a[b++]&63;e=224==(e&240)?(e&15)<<12|f<<6|g:(e&7)<<18|f<<12|g<<6|a[b++]&63;65536>e?d+=String.fromCharCode(e):(e-=65536,d+=String.fromCharCode(55296|e>>10,56320|e&1023))}}else d+=String.fromCharCode(e)}return d},Jb="undefined"!=typeof TextDecoder?new TextDecoder("utf-16le"):void 0,
Kb=(a,b)=>{var c=a>>1;for(var d=c+b/2;!(c>=d)&&oa[c];)++c;c<<=1;if(32<c-a&&Jb)return Jb.decode(v.subarray(a,c));c="";for(d=0;!(d>=b/2);++d){var e=B[a+2*d>>1];if(0==e)break;c+=String.fromCharCode(e)}return c},Lb=(a,b,c)=>{c??=2147483647;if(2>c)return 0;c-=2;var d=b;c=c<2*a.length?c/2:a.length;for(var e=0;e<c;++e)B[b>>1]=a.charCodeAt(e),b+=2;B[b>>1]=0;return b-d},Mb=a=>2*a.length,Nb=(a,b)=>{for(var c=0,d="";!(c>=b/4);){var e=C[a+4*c>>2];if(0==e)break;++c;65536<=e?(e-=65536,d+=String.fromCharCode(55296|
e>>10,56320|e&1023)):d+=String.fromCharCode(e)}return d},Ob=(a,b,c)=>{c??=2147483647;if(4>c)return 0;var d=b;c=d+c-4;for(var e=0;e<a.length;++e){var f=a.charCodeAt(e);if(55296<=f&&57343>=f){var g=a.charCodeAt(++e);f=65536+((f&1023)<<10)|g&1023}C[b>>2]=f;b+=4;if(b+4>c)break}C[b>>2]=0;return b-d},Pb=a=>{for(var b=0,c=0;c<a.length;++c){var d=a.charCodeAt(c);55296<=d&&57343>=d&&++c;b+=4}return b},Qb=0,Rb=(a,b)=>{var c=K[a];if(void 0===c)throw a=`${b} has unknown type ${sb(a)}`,new Q(a);return c},Sb={},
Tb=a=>{if(!(a instanceof Ba||"unwind"==a))throw a;},Ub=a=>{na=a;Da||0<Qb||(k.onExit?.(a),ma=!0);throw new Ba(a);},Vb=a=>{if(!ma)try{if(a(),!(Da||0<Qb))try{na=a=na,Ub(a)}catch(b){Tb(b)}}catch(b){Tb(b)}},Wb={},Yb=()=>{if(!Xb){var a={USER:"web_user",LOGNAME:"web_user",PATH:"/",PWD:"/",HOME:"/home/web_user",LANG:("object"==typeof navigator&&navigator.languages&&navigator.languages[0]||"C").replace("-","_")+".UTF-8",_:ea||"./this.program"},b;for(b in Wb)void 0===Wb[b]?delete a[b]:a[b]=Wb[b];var c=[];for(b in a)c.push(`${b}=${a[b]}`);
Xb=c}return Xb},Xb,Zb=(a,b)=>{for(var c=0,d=a.length-1;0<=d;d--){var e=a[d];"."===e?a.splice(d,1):".."===e?(a.splice(d,1),c++):c&&(a.splice(d,1),c--)}if(b)for(;c;c--)a.unshift("..");return a},$b=a=>{var b="/"===a.charAt(0),c="/"===a.substr(-1);(a=Zb(a.split("/").filter(d=>!!d),!b).join("/"))||b||(a=".");a&&c&&(a+="/");return(b?"/":"")+a},ac=a=>{var b=/^(\/?|)([\s\S]*?)((?:\.{1,2}|[^\/]+?|)(\.[^.\/]*|))(?:[\/]*)$/.exec(a).slice(1);a=b[0];b=b[1];if(!a&&!b)return".";b&&=b.substr(0,b.length-1);return a+
b},bc=a=>{if("/"===a)return"/";a=$b(a);a=a.replace(/\/$/,"");var b=a.lastIndexOf("/");return-1===b?a:a.substr(b+1)},cc=()=>{if("object"==typeof crypto&&"function"==typeof crypto.getRandomValues)return a=>crypto.getRandomValues(a);G("initRandomDevice")},ec=a=>(ec=cc())(a),fc=(...a)=>{for(var b="",c=!1,d=a.length-1;-1<=d&&!c;d--){c=0<=d?a[d]:"/";if("string"!=typeof c)throw new TypeError("Arguments to path.resolve must be strings");if(!c)return"";b=c+"/"+b;c="/"===c.charAt(0)}b=Zb(b.split("/").filter(e=>
!!e),!c).join("/");return(c?"/":"")+b||"."},gc=[],hc=[];function ic(a,b){hc[a]={input:[],ea:[],xa:b};jc(a,kc)}
var kc={open(a){var b=hc[a.node.La];if(!b)throw new U(43);a.ha=b;a.seekable=!1},close(a){a.ha.xa.Ga(a.ha)},Ga(a){a.ha.xa.Ga(a.ha)},read(a,b,c,d){if(!a.ha||!a.ha.xa.$a)throw new U(60);for(var e=0,f=0;f<d;f++){try{var g=a.ha.xa.$a(a.ha)}catch(h){throw new U(29);}if(void 0===g&&0===e)throw new U(6);if(null===g||void 0===g)break;e++;b[c+f]=g}e&&(a.node.sa=Date.now());return e},write(a,b,c,d){if(!a.ha||!a.ha.xa.Qa)throw new U(60);try{for(var e=0;e<d;e++)a.ha.xa.Qa(a.ha,b[c+e])}catch(f){throw new U(29);
}d&&(a.node.ka=a.node.ia=Date.now());return e}},lc={$a(){return gc.length?gc.shift():null},Qa(a,b){null===b||10===b?(ja(Ib(a.ea)),a.ea=[]):0!=b&&a.ea.push(b)},Ga(a){a.ea&&0<a.ea.length&&(ja(Ib(a.ea)),a.ea=[])},Ob(){return{Hb:25856,Jb:5,Gb:191,Ib:35387,Fb:[3,28,127,21,4,0,1,0,17,19,26,0,18,15,23,22,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]}},Pb(){return 0},Qb(){return[24,80]}},mc={Qa(a,b){null===b||10===b?(q(Ib(a.ea)),a.ea=[]):0!=b&&a.ea.push(b)},Ga(a){a.ea&&0<a.ea.length&&(q(Ib(a.ea)),a.ea=[])}};
function nc(a,b){var c=a.X?a.X.length:0;c>=b||(b=Math.max(b,c*(1048576>c?2:1.125)>>>0),0!=c&&(b=Math.max(b,256)),c=a.X,a.X=new Uint8Array(b),0<a.ba&&a.X.set(c.subarray(0,a.ba),0))}
var V={ma:null,qa(){return V.createNode(null,"/",16895,0)},createNode(a,b,c,d){if(24576===(c&61440)||4096===(c&61440))throw new U(63);V.ma||(V.ma={dir:{node:{ra:V.W.ra,oa:V.W.oa,Ba:V.W.Ba,Ja:V.W.Ja,gb:V.W.gb,ib:V.W.ib,hb:V.W.hb,Sa:V.W.Sa,Ma:V.W.Ma},stream:{la:V.Y.la}},file:{node:{ra:V.W.ra,oa:V.W.oa},stream:{la:V.Y.la,read:V.Y.read,write:V.Y.write,Ua:V.Y.Ua,bb:V.Y.bb,eb:V.Y.eb}},link:{node:{ra:V.W.ra,oa:V.W.oa,Da:V.W.Da},stream:{}},Va:{node:{ra:V.W.ra,oa:V.W.oa},stream:oc}});c=pc(a,b,c,d);16384===
(c.mode&61440)?(c.W=V.ma.dir.node,c.Y=V.ma.dir.stream,c.X={}):32768===(c.mode&61440)?(c.W=V.ma.file.node,c.Y=V.ma.file.stream,c.ba=0,c.X=null):40960===(c.mode&61440)?(c.W=V.ma.link.node,c.Y=V.ma.link.stream):8192===(c.mode&61440)&&(c.W=V.ma.Va.node,c.Y=V.ma.Va.stream);c.sa=c.ka=c.ia=Date.now();a&&(a.X[b]=c,a.sa=a.ka=a.ia=c.sa);return c},Lb(a){return a.X?a.X.subarray?a.X.subarray(0,a.ba):new Uint8Array(a.X):new Uint8Array(0)},W:{ra(a){var b={};b.Kb=8192===(a.mode&61440)?a.id:1;b.Nb=a.id;b.mode=a.mode;
b.Sb=1;b.uid=0;b.Mb=0;b.La=a.La;16384===(a.mode&61440)?b.size=4096:32768===(a.mode&61440)?b.size=a.ba:40960===(a.mode&61440)?b.size=a.link.length:b.size=0;b.sa=new Date(a.sa);b.ka=new Date(a.ka);b.ia=new Date(a.ia);b.jb=4096;b.Eb=Math.ceil(b.size/b.jb);return b},oa(a,b){for(var c of["mode","atime","mtime","ctime"])b[c]&&(a[c]=b[c]);void 0!==b.size&&(b=b.size,a.ba!=b&&(0==b?(a.X=null,a.ba=0):(c=a.X,a.X=new Uint8Array(b),c&&a.X.set(c.subarray(0,Math.min(b,a.ba))),a.ba=b)))},Ba(){throw V.Xa;},Ja(a,b,
c,d){return V.createNode(a,b,c,d)},gb(a,b,c){try{var d=qc(b,c)}catch(f){}if(d){if(16384===(a.mode&61440))for(var e in d.X)throw new U(55);e=rc(d.parent.id,d.name);if(W[e]===d)W[e]=d.wa;else for(e=W[e];e;){if(e.wa===d){e.wa=d.wa;break}e=e.wa}}delete a.parent.X[a.name];b.X[c]=a;a.name=c;b.ia=b.ka=a.parent.ia=a.parent.ka=Date.now()},ib(a,b){delete a.X[b];a.ia=a.ka=Date.now()},hb(a,b){var c=qc(a,b),d;for(d in c.X)throw new U(55);delete a.X[b];a.ia=a.ka=Date.now()},Sa(a){return[".","..",...Object.keys(a.X)]},
Ma(a,b,c){a=V.createNode(a,b,41471,0);a.link=c;return a},Da(a){if(40960!==(a.mode&61440))throw new U(28);return a.link}},Y:{read(a,b,c,d,e){var f=a.node.X;if(e>=a.node.ba)return 0;a=Math.min(a.node.ba-e,d);if(8<a&&f.subarray)b.set(f.subarray(e,e+a),c);else for(d=0;d<a;d++)b[c+d]=f[e+d];return a},write(a,b,c,d,e,f){if(!d)return 0;a=a.node;a.ka=a.ia=Date.now();if(b.subarray&&(!a.X||a.X.subarray)){if(f)return a.X=b.subarray(c,c+d),a.ba=d;if(0===a.ba&&0===e)return a.X=b.slice(c,c+d),a.ba=d;if(e+d<=a.ba)return a.X.set(b.subarray(c,
c+d),e),d}nc(a,e+d);if(a.X.subarray&&b.subarray)a.X.set(b.subarray(c,c+d),e);else for(f=0;f<d;f++)a.X[e+f]=b[c+f];a.ba=Math.max(a.ba,e+d);return d},la(a,b,c){1===c?b+=a.position:2===c&&32768===(a.node.mode&61440)&&(b+=a.node.ba);if(0>b)throw new U(28);return b},Ua(a,b,c){nc(a.node,b+c);a.node.ba=Math.max(a.node.ba,b+c)},bb(a,b,c,d,e){if(32768!==(a.node.mode&61440))throw new U(43);a=a.node.X;if(e&2||!a||a.buffer!==r.buffer){d=!0;G();e=void 0;if(!e)throw new U(48);if(a){if(0<c||c+b<a.length)a.subarray?
a=a.subarray(c,c+b):a=Array.prototype.slice.call(a,c,c+b);r.set(a,e)}}else d=!1,e=a.byteOffset;return{Z:e,Db:d}},eb(a,b,c,d){V.Y.write(a,b,0,d,c,!1);return 0}}},sc=(a,b)=>{var c=0;a&&(c|=365);b&&(c|=146);return c},tc=null,uc={},vc=[],wc=1,W=null,xc=!1,yc=!0,U=class{name="ErrnoError";constructor(a){this.Aa=a}},zc={},Ac=class{Fa={};node=null;get flags(){return this.Fa.flags}set flags(a){this.Fa.flags=a}get position(){return this.Fa.position}set position(a){this.Fa.position=a}},Bc=class{W={};Y={};Ka=null;constructor(a,
b,c,d){a||=this;this.parent=a;this.qa=a.qa;this.id=wc++;this.name=b;this.mode=c;this.La=d;this.sa=this.ka=this.ia=Date.now()}get read(){return 365===(this.mode&365)}set read(a){a?this.mode|=365:this.mode&=-366}get write(){return 146===(this.mode&146)}set write(a){a?this.mode|=146:this.mode&=-147}};
function Cc(a,b={}){if(!a)return{path:"",node:null};b.Na??(b.Na=!0);"/"===a.charAt(0)||(a="//"+a);var c=0;a:for(;40>c;c++){a=a.split("/").filter(h=>!!h&&"."!==h);for(var d=tc,e="/",f=0;f<a.length;f++){var g=f===a.length-1;if(g&&b.parent)break;if(".."===a[f])e=ac(e),d=d.parent;else{e=$b(e+"/"+a[f]);try{d=qc(d,a[f])}catch(h){if(44===h?.Aa&&g&&b.ub)return{path:e};throw h;}!d.Ka||g&&!b.Na||(d=d.Ka.root);if(40960===(d.mode&61440)&&(!g||b.Za)){if(!d.W.Da)throw new U(52);d=d.W.Da(d);"/"===d.charAt(0)||(d=
ac(e)+"/"+d);a=d+"/"+a.slice(f+1).join("/");continue a}}}return{path:e,node:d}}throw new U(32);}function rc(a,b){for(var c=0,d=0;d<b.length;d++)c=(c<<5)-c+b.charCodeAt(d)|0;return(a+c>>>0)%W.length}function qc(a,b){var c=16384===(a.mode&61440)?(c=Dc(a,"x"))?c:a.W.Ba?0:2:54;if(c)throw new U(c);for(c=W[rc(a.id,b)];c;c=c.wa){var d=c.name;if(c.parent.id===a.id&&d===b)return c}return a.W.Ba(a,b)}function pc(a,b,c,d){a=new Bc(a,b,c,d);b=rc(a.parent.id,a.name);a.wa=W[b];return W[b]=a}
function Ec(a){var b=["r","w","rw"][a&3];a&512&&(b+="w");return b}function Dc(a,b){if(yc)return 0;if(!b.includes("r")||a.mode&292){if(b.includes("w")&&!(a.mode&146)||b.includes("x")&&!(a.mode&73))return 2}else return 2;return 0}function Fc(a,b){if(16384!==(a.mode&61440))return 54;try{return qc(a,b),20}catch(c){}return Dc(a,"wx")}function Gc(a){a=vc[a];if(!a)throw new U(8);return a}var oc={open(a){a.Y=uc[a.node.La].Y;a.Y.open?.(a)},la(){throw new U(70);}};function jc(a,b){uc[a]={Y:b}}
function Hc(a,b){var c="/"===b;if(c&&tc)throw new U(10);if(!c&&b){var d=Cc(b,{Na:!1});b=d.path;d=d.node;if(d.Ka)throw new U(10);if(16384!==(d.mode&61440))throw new U(54);}b={type:a,Tb:{},cb:b,tb:[]};a=a.qa(b);a.qa=b;b.root=a;c?tc=a:d&&(d.Ka=b,d.qa&&d.qa.tb.push(b))}function Ic(a,b,c){var d=Cc(a,{parent:!0}).node;a=bc(a);if(!a||"."===a||".."===a)throw new U(28);var e=Fc(d,a);if(e)throw new U(e);if(!d.W.Ja)throw new U(63);return d.W.Ja(d,a,b,c)}function X(a){return Ic(a,16895,0)}
function Jc(a,b,c){"undefined"==typeof c&&(c=b,b=438);Ic(a,b|8192,c)}function Kc(a,b){if(!fc(a))throw new U(44);var c=Cc(b,{parent:!0}).node;if(!c)throw new U(44);b=bc(b);var d=Fc(c,b);if(d)throw new U(d);if(!c.W.Ma)throw new U(63);c.W.Ma(c,b,a)}
function Lc(a,b){if(""===a)throw new U(44);if("string"==typeof b){var c={r:0,"r+":2,w:577,"w+":578,a:1089,"a+":1090}[b];if("undefined"==typeof c)throw Error(`Unknown file open mode: ${b}`);b=c}var d=b&64?33206:0;"object"==typeof a?c=a:(a=Cc(a,{Za:!(b&131072),ub:!0}),c=a.node,a=a.path);var e=!1;if(b&64)if(c){if(b&128)throw new U(20);}else c=Ic(a,d,0),e=!0;if(!c)throw new U(44);8192===(c.mode&61440)&&(b&=-513);if(b&65536&&16384!==(c.mode&61440))throw new U(54);if(!e&&(d=c?40960===(c.mode&61440)?32:
16384===(c.mode&61440)&&("r"!==Ec(b)||b&512)?31:Dc(c,Ec(b)):44))throw new U(d);if(b&512&&!e){d=c;d="string"==typeof d?Cc(d,{Za:!0}).node:d;if(!d.W.oa)throw new U(63);if(16384===(d.mode&61440))throw new U(31);if(32768!==(d.mode&61440))throw new U(28);if(e=Dc(d,"w"))throw new U(e);d.W.oa(d,{size:0,timestamp:Date.now()})}b&=-131713;a:for(d=c;;){if(d===d.parent){d=d.qa.cb;var f=f?"/"!==d[d.length-1]?`${d}/${f}`:d+f:d;break a}f=f?`${d.name}/${f}`:d.name;d=d.parent}f={node:c,path:f,flags:b,seekable:!0,
position:0,Y:c.Y,Cb:[],error:!1};c=-1;f=Object.assign(new Ac,f);if(-1==c)a:{for(c=0;4096>=c;c++)if(!vc[c])break a;throw new U(33);}f.ua=c;vc[c]=f;f.Y.open&&f.Y.open(f);!k.logReadFiles||b&1||a in zc||(zc[a]=1)}function Mc(a,b,c){if(null===a.ua)throw new U(8);if(!a.seekable||!a.Y.la)throw new U(70);if(0!=c&&1!=c&&2!=c)throw new U(28);a.position=a.Y.la(a,b,c);a.Cb=[]}
function Y(a,b,c){a=$b("/dev/"+a);var d=sc(!!b,!!c);Y.ab??(Y.ab=64);var e=Y.ab++<<8|0;jc(e,{open(f){f.seekable=!1},close(){c?.buffer?.length&&c(10)},read(f,g,h,m){for(var l=0,n=0;n<m;n++){try{var t=b()}catch(u){throw new U(29);}if(void 0===t&&0===l)throw new U(6);if(null===t||void 0===t)break;l++;g[h+n]=t}l&&(f.node.sa=Date.now());return l},write(f,g,h,m){for(var l=0;l<m;l++)try{c(g[h+l])}catch(n){throw new U(29);}m&&(f.node.ka=f.node.ia=Date.now());return l}});Jc(a,d,e)}var Nc={};
L=k.InternalError=class extends Error{constructor(a){super(a);this.name="InternalError"}};for(var Oc=Array(256),Pc=0;256>Pc;++Pc)Oc[Pc]=String.fromCharCode(Pc);Ka=Oc;Q=k.BindingError=class extends Error{constructor(a){super(a);this.name="BindingError"}};
Object.assign(Xa.prototype,{isAliasOf:function(a){if(!(this instanceof Xa&&a instanceof Xa))return!1;var b=this.V.aa.$,c=this.V.Z;a.V=a.V;var d=a.V.aa.$;for(a=a.V.Z;b.fa;)c=b.Ea(c),b=b.fa;for(;d.fa;)a=d.Ea(a),d=d.fa;return b===d&&c===a},clone:function(){this.V.Z||Na(this);if(this.V.Ca)return this.V.count.value+=1,this;var a=Ua,b=Object,c=b.create,d=Object.getPrototypeOf(this),e=this.V;a=a(c.call(b,d,{V:{value:{count:e.count,za:e.za,Ca:e.Ca,Z:e.Z,aa:e.aa,da:e.da,ga:e.ga}}}));a.V.count.value+=1;a.V.za=
!1;return a},["delete"](){this.V.Z||Na(this);if(this.V.za&&!this.V.Ca)throw new Q("Object already scheduled for deletion");Pa(this);var a=this.V;--a.count.value;0===a.count.value&&(a.da?a.ga.na(a.da):a.aa.$.na(a.Z));this.V.Ca||(this.V.da=void 0,this.V.Z=void 0)},isDeleted:function(){return!this.V.Z},deleteLater:function(){this.V.Z||Na(this);if(this.V.za&&!this.V.Ca)throw new Q("Object already scheduled for deletion");Wa.push(this);this.V.za=!0;return this}});
Object.assign(jb.prototype,{ob(a){this.fb&&(a=this.fb(a));return a},Wa(a){this.na?.(a)},pa:8,readValueFromPointer:I,fromWireType:function(a){function b(){return this.Ia?Va(this.$.va,{aa:this.vb,Z:c,ga:this,da:a}):Va(this.$.va,{aa:this,Z:a})}var c=this.ob(a);if(!c)return this.Wa(a),null;var d=Ta(this.$,c);if(void 0!==d){if(0===d.V.count.value)return d.V.Z=c,d.V.da=a,d.clone();d=d.clone();this.Wa(a);return d}d=this.$.nb(c);d=Ra[d];if(!d)return b.call(this);d=this.Ha?d.kb:d.pointerType;var e=Qa(c,this.$,
d.$);return null===e?b.call(this):this.Ia?Va(d.$.va,{aa:d,Z:e,ga:this,da:a}):Va(d.$.va,{aa:d,Z:e})}});qb=k.UnboundTypeError=((a,b)=>{var c=Ya(b,function(d){this.name=b;this.message=d;d=Error(d).stack;void 0!==d&&(this.stack=this.toString()+"\n"+d.replace(/^Error(:[^\n]*)?\n/,""))});c.prototype=Object.create(a.prototype);c.prototype.constructor=c;c.prototype.toString=function(){return void 0===this.message?this.name:`${this.name}: ${this.message}`};return c})(Error,"UnboundTypeError");
T.push(0,1,void 0,1,null,1,!0,1,!1,1);k.count_emval_handles=()=>T.length/2-5-zb.length;W=Array(4096);Hc(V,"/");X("/tmp");X("/home");X("/home/web_user");(function(){X("/dev");jc(259,{read:()=>0,write:(d,e,f,g)=>g,la:()=>0});Jc("/dev/null",259);ic(1280,lc);ic(1536,mc);Jc("/dev/tty",1280);Jc("/dev/tty1",1536);var a=new Uint8Array(1024),b=0,c=()=>{0===b&&(b=ec(a).byteLength);return a[--b]};Y("random",c);Y("urandom",c);X("/dev/shm");X("/dev/shm/tmp")})();
(function(){X("/proc");var a=X("/proc/self");X("/proc/self/fd");Hc({qa(){var b=pc(a,"fd",16895,73);b.Y={la:V.Y.la};b.W={Ba(c,d){c=+d;var e=Gc(c);c={parent:null,qa:{cb:"fake"},W:{Da:()=>e.path},id:c+1};return c.parent=c},Sa(){return Array.from(vc.entries()).filter(([,c])=>c).map(([c])=>c.toString())}};return b}},"/proc/self/fd")})();V.Xa=new U(44);V.Xa.stack="<generic error, no stack>";
var Sc={i:(a,b,c)=>{var d=new Ea(a);D[d.Z+16>>2]=0;D[d.Z+4>>2]=b;D[d.Z+8>>2]=c;Fa=a;Ga++;throw Fa;},A:()=>G(""),B:a=>{var b=Ha[a];delete Ha[a];var c=b.Ra,d=b.na,e=b.Ya,f=e.map(g=>g.rb).concat(e.map(g=>g.zb));N([a],f,g=>{var h={};e.forEach((m,l)=>{var n=g[l],t=m.pb,u=m.qb,x=g[l+e.length],w=m.yb,y=m.Ab;h[m.mb]={read:z=>n.fromWireType(t(u,z)),write:(z,O)=>{var A=[];w(y,z,x.toWireType(A,O));Ia(A)}}});return[{name:b.name,fromWireType:m=>{var l={},n;for(n in h)l[n]=h[n].read(m);d(m);return l},toWireType:(m,
l)=>{for(var n in h)if(!(n in l))throw new TypeError(`Missing field: "${n}"`);var t=c();for(n in h)h[n].write(t,l[n]);null!==m&&m.push(d,t);return t},pa:8,readValueFromPointer:I,ja:d}]})},o:()=>{},E:(a,b,c,d)=>{b=P(b);M(a,{name:b,fromWireType:function(e){return!!e},toWireType:function(e,f){return f?c:d},pa:8,readValueFromPointer:function(e){return this.fromWireType(v[e])},ja:null})},m:(a,b,c,d,e,f,g,h,m,l,n,t,u)=>{n=P(n);f=R(e,f);h&&=R(g,h);l&&=R(m,l);u=R(t,u);var x=ab(n);$a(x,function(){tb(`Cannot construct ${n} due to unbound types`,
[d])});N([a,b,c],d?[d]:[],w=>{w=w[0];if(d){var y=w.$;var z=y.va}else z=Xa.prototype;w=Ya(n,function(...bb){if(Object.getPrototypeOf(this)!==O)throw new Q("Use 'new' to construct "+n);if(void 0===A.ta)throw new Q(n+" has no accessible constructor");var dc=A.ta[bb.length];if(void 0===dc)throw new Q(`Tried to invoke ctor of ${n} with invalid number of parameters (${bb.length}) - expected (${Object.keys(A.ta).toString()}) parameters instead!`);return dc.apply(this,bb)});var O=Object.create(z,{constructor:{value:w}});
w.prototype=O;var A=new cb(n,w,O,u,y,f,h,l);if(A.fa){var ia;(ia=A.fa).Ta??(ia.Ta=[]);A.fa.Ta.push(A)}y=new jb(n,A,!0,!1,!1);ia=new jb(n+"*",A,!1,!1,!1);z=new jb(n+" const*",A,!1,!0,!1);Ra[a]={pointerType:ia,kb:z};kb(x,w);return[y,ia,z]})},l:(a,b,c,d,e,f)=>{var g=ub(b,c);e=R(d,e);N([],[a],h=>{h=h[0];var m=`constructor ${h.name}`;void 0===h.$.ta&&(h.$.ta=[]);if(void 0!==h.$.ta[b-1])throw new Q(`Cannot register multiple constructors with identical number of parameters (${b-1}) for class '${h.name}'! Overload resolution is currently only performed using the parameter count, not actual type info!`);
h.$.ta[b-1]=()=>{tb(`Cannot construct ${h.name} due to unbound types`,g)};N([],g,l=>{l.splice(1,0,null);h.$.ta[b-1]=xb(m,l,null,e,f);return[]});return[]})},c:(a,b,c,d,e,f,g,h,m)=>{var l=ub(c,d);b=P(b);b=yb(b);f=R(e,f);N([],[a],n=>{function t(){tb(`Cannot call ${u} due to unbound types`,l)}n=n[0];var u=`${n.name}.${b}`;b.startsWith("@@")&&(b=Symbol[b.substring(2)]);h&&n.$.wb.push(b);var x=n.$.va,w=x[b];void 0===w||void 0===w.ca&&w.className!==n.name&&w.ya===c-2?(t.ya=c-2,t.className=n.name,x[b]=t):
(Za(x,b,u),x[b].ca[c-2]=t);N([],l,y=>{y=xb(u,y,n,f,g,m);void 0===x[b].ca?(y.ya=c-2,x[b]=y):x[b].ca[c-2]=y;return[]});return[]})},C:a=>M(a,Cb),h:(a,b,c)=>{b=P(b);M(a,{name:b,fromWireType:d=>d,toWireType:(d,e)=>e,pa:8,readValueFromPointer:Db(b,c),ja:null})},f:(a,b,c,d,e,f,g)=>{var h=ub(b,c);a=P(a);a=yb(a);e=R(d,e);$a(a,function(){tb(`Cannot call ${a} due to unbound types`,h)},b-1);N([],h,m=>{kb(a,xb(a,[m[0],null].concat(m.slice(1)),null,e,f,g),b-1);return[]})},b:(a,b,c,d,e)=>{b=P(b);-1===e&&(e=4294967295);
e=h=>h;if(0===d){var f=32-8*c;e=h=>h<<f>>>f}var g=b.includes("unsigned")?function(h,m){return m>>>0}:function(h,m){return m};M(a,{name:b,fromWireType:e,toWireType:g,pa:8,readValueFromPointer:Eb(b,c,0!==d),ja:null})},a:(a,b,c)=>{function d(f){return new e(r.buffer,D[f+4>>2],D[f>>2])}var e=[Int8Array,Uint8Array,Int16Array,Uint16Array,Int32Array,Uint32Array,Float32Array,Float64Array][b];c=P(c);M(a,{name:c,fromWireType:d,pa:8,readValueFromPointer:d},{sb:!0})},g:a=>{M(a,Fb)},D:(a,b)=>{b=P(b);M(a,{name:b,
fromWireType:function(c){for(var d=D[c>>2],e=c+4,f,g=e,h=0;h<=d;++h){var m=e+h;if(h==d||0==v[m])g=g?Ib(v,g,m-g):"",void 0===f?f=g:(f+=String.fromCharCode(0),f+=g),g=m+1}S(c);return f},toWireType:function(c,d){d instanceof ArrayBuffer&&(d=new Uint8Array(d));var e,f="string"==typeof d;if(!(f||d instanceof Uint8Array||d instanceof Uint8ClampedArray||d instanceof Int8Array))throw new Q("Cannot pass non-string to std::string");if(f)for(var g=e=0;g<d.length;++g){var h=d.charCodeAt(g);127>=h?e++:2047>=h?
e+=2:55296<=h&&57343>=h?(e+=4,++g):e+=3}else e=d.length;g=Qc(4+e+1);h=g+4;D[g>>2]=e;if(f)Gb(d,h,e+1);else if(f)for(f=0;f<e;++f){var m=d.charCodeAt(f);if(255<m)throw S(h),new Q("String has UTF-16 code units that do not fit in 8 bits");v[h+f]=m}else for(f=0;f<e;++f)v[h+f]=d[f];null!==c&&c.push(S,g);return g},pa:8,readValueFromPointer:I,ja(c){S(c)}})},d:(a,b,c)=>{c=P(c);if(2===b){var d=Kb;var e=Lb;var f=Mb;var g=h=>oa[h>>1]}else 4===b&&(d=Nb,e=Ob,f=Pb,g=h=>D[h>>2]);M(a,{name:c,fromWireType:h=>{for(var m=
D[h>>2],l,n=h+4,t=0;t<=m;++t){var u=h+4+t*b;if(t==m||0==g(u))n=d(n,u-n),void 0===l?l=n:(l+=String.fromCharCode(0),l+=n),n=u+b}S(h);return l},toWireType:(h,m)=>{if("string"!=typeof m)throw new Q(`Cannot pass non-string to C++ string type ${c}`);var l=f(m),n=Qc(4+l+b);D[n>>2]=l/b;e(m,n+4,l+b);null!==h&&h.push(S,n);return n},pa:8,readValueFromPointer:I,ja(h){S(h)}})},I:(a,b,c,d,e,f)=>{Ha[a]={name:P(b),Ra:R(c,d),na:R(e,f),Ya:[]}},j:(a,b,c,d,e,f,g,h,m,l)=>{Ha[a].Ya.push({mb:P(b),rb:c,pb:R(d,e),qb:f,zb:g,
yb:R(h,m),Ab:l})},F:(a,b)=>{b=P(b);M(a,{Rb:!0,name:b,pa:0,fromWireType:()=>{},toWireType:()=>{}})},z:(a,b,c)=>v.copyWithin(a,b,b+c),q:()=>{Da=!1;Qb=0},k:(a,b,c)=>{a=Bb(a);b=Rb(b,"emval::as");var d=[];a=b.toWireType(d,a);d.length&&(D[c>>2]=hb(d));return a},G:Ab,H:a=>{var b=Bb(a);Ia(b);Ab(a)},e:(a,b)=>{a=Rb(a,"_emval_take_value");a=a.readValueFromPointer(b);return hb(a)},r:(a,b)=>{Sb[a]&&(clearTimeout(Sb[a].id),delete Sb[a]);if(!b)return 0;var c=setTimeout(()=>{delete Sb[a];Vb(()=>Rc(a,performance.now()))},
b);Sb[a]={id:c,Ub:b};return 0},s:(a,b,c,d)=>{var e=(new Date).getFullYear(),f=(new Date(e,0,1)).getTimezoneOffset();e=(new Date(e,6,1)).getTimezoneOffset();D[a>>2]=60*Math.max(f,e);C[b>>2]=Number(f!=e);b=g=>{var h=Math.abs(g);return`UTC${0<=g?"-":"+"}${String(Math.floor(h/60)).padStart(2,"0")}${String(h%60).padStart(2,"0")}`};a=b(f);b=b(e);e<f?(Gb(a,c,17),Gb(b,d,17)):(Gb(a,d,17),Gb(b,c,17))},y:()=>{G("OOM")},t:(a,b)=>{var c=0;Yb().forEach((d,e)=>{var f=b+c;e=D[a+4*e>>2]=f;for(f=0;f<d.length;++f)r[e++]=
d.charCodeAt(f);r[e]=0;c+=d.length+1});return 0},u:(a,b)=>{var c=Yb();D[a>>2]=c.length;var d=0;c.forEach(e=>d+=e.length+1);D[b>>2]=d;return 0},v:function(a){try{var b=Gc(a);if(null===b.ua)throw new U(8);b.Oa&&(b.Oa=null);try{b.Y.close&&b.Y.close(b)}catch(c){throw c;}finally{vc[b.ua]=null}b.ua=null;return 0}catch(c){if("undefined"==typeof Nc||"ErrnoError"!==c.name)throw c;return c.Aa}},x:function(a,b,c,d){try{a:{var e=Gc(a);a=b;for(var f,g=b=0;g<c;g++){var h=D[a>>2],m=D[a+4>>2];a+=8;var l=e,n=f,t=
r;if(0>m||0>n)throw new U(28);if(null===l.ua)throw new U(8);if(1===(l.flags&2097155))throw new U(8);if(16384===(l.node.mode&61440))throw new U(31);if(!l.Y.read)throw new U(28);var u="undefined"!=typeof n;if(!u)n=l.position;else if(!l.seekable)throw new U(70);var x=l.Y.read(l,t,h,m,n);u||(l.position+=x);var w=x;if(0>w){var y=-1;break a}b+=w;if(w<m)break;"undefined"!=typeof f&&(f+=w)}y=b}D[d>>2]=y;return 0}catch(z){if("undefined"==typeof Nc||"ErrnoError"!==z.name)throw z;return z.Aa}},n:function(a,
b,c,d,e){b=c+2097152>>>0<4194305-!!b?(b>>>0)+4294967296*c:NaN;try{if(isNaN(b))return 61;var f=Gc(a);Mc(f,b,d);Aa=[f.position>>>0,(H=f.position,1<=+Math.abs(H)?0<H?+Math.floor(H/4294967296)>>>0:~~+Math.ceil((H-+(~~H>>>0))/4294967296)>>>0:0)];C[e>>2]=Aa[0];C[e+4>>2]=Aa[1];f.Oa&&0===b&&0===d&&(f.Oa=null);return 0}catch(g){if("undefined"==typeof Nc||"ErrnoError"!==g.name)throw g;return g.Aa}},w:function(a,b,c,d){try{a:{var e=Gc(a);a=b;for(var f,g=b=0;g<c;g++){var h=D[a>>2],m=D[a+4>>2];a+=8;var l=e,n=
h,t=m,u=f,x=r;if(0>t||0>u)throw new U(28);if(null===l.ua)throw new U(8);if(0===(l.flags&2097155))throw new U(8);if(16384===(l.node.mode&61440))throw new U(31);if(!l.Y.write)throw new U(28);l.seekable&&l.flags&1024&&Mc(l,0,2);var w="undefined"!=typeof u;if(!w)u=l.position;else if(!l.seekable)throw new U(70);var y=l.Y.write(l,x,n,t,u,void 0);w||(l.position+=y);var z=y;if(0>z){var O=-1;break a}b+=z;if(z<m)break;"undefined"!=typeof f&&(f+=z)}O=b}D[d>>2]=O;return 0}catch(A){if("undefined"==typeof Nc||
"ErrnoError"!==A.name)throw A;return A.Aa}},p:Ub},Z;
(async function(){function a(d){Z=d.exports;la=Z.J;d=la.buffer;k.HEAP8=r=new Int8Array(d);k.HEAP16=B=new Int16Array(d);k.HEAPU8=v=new Uint8Array(d);k.HEAPU16=oa=new Uint16Array(d);k.HEAP32=C=new Int32Array(d);k.HEAPU32=D=new Uint32Array(d);k.HEAPF32=pa=new Float32Array(d);k.HEAPF64=qa=new Float64Array(d);mb=Z.M;sa.unshift(Z.K);E--;k.monitorRunDependencies?.(E);0==E&&F&&(d=F,F=null,d());return Z}E++;k.monitorRunDependencies?.(E);var b={a:Sc};if(k.instantiateWasm)try{return k.instantiateWasm(b,a)}catch(d){q(`Module.instantiateWasm callback failed with error: ${d}`),
ba(d)}wa??=k.locateFile?va("rsiscool_workers.wasm")?"rsiscool_workers.wasm":k.locateFile?k.locateFile("rsiscool_workers.wasm",p):p+"rsiscool_workers.wasm":(new URL("rsiscool_workers.wasm",import.meta.url)).href;try{var c=await za(b);a(c.instance);return c}catch(d){ba(d)}})();var rb=a=>(rb=Z.L)(a),Rc=(a,b)=>(Rc=Z.N)(a,b),Qc=a=>(Qc=Z.O)(a),S=a=>(S=Z.P)(a);k.dynCall_viijii=(a,b,c,d,e,f,g)=>(k.dynCall_viijii=Z.Q)(a,b,c,d,e,f,g);k.dynCall_jiji=(a,b,c,d,e)=>(k.dynCall_jiji=Z.R)(a,b,c,d,e);
k.dynCall_iiiiij=(a,b,c,d,e,f,g)=>(k.dynCall_iiiiij=Z.S)(a,b,c,d,e,f,g);k.dynCall_iiiiijj=(a,b,c,d,e,f,g,h,m)=>(k.dynCall_iiiiijj=Z.T)(a,b,c,d,e,f,g,h,m);k.dynCall_iiiiiijj=(a,b,c,d,e,f,g,h,m,l)=>(k.dynCall_iiiiiijj=Z.U)(a,b,c,d,e,f,g,h,m,l);var Tc;F=function Uc(){Tc||Vc();Tc||(F=Uc)};
function Vc(){function a(){if(!Tc&&(Tc=!0,k.calledRun=!0,!ma)){if(!k.noFSInit&&!xc){var b,c;xc=!0;d??=k.stdin;b??=k.stdout;c??=k.stderr;d?Y("stdin",d):Kc("/dev/tty","/dev/stdin");b?Y("stdout",null,b):Kc("/dev/tty","/dev/stdout");c?Y("stderr",null,c):Kc("/dev/tty1","/dev/stderr");Lc("/dev/stdin",0);Lc("/dev/stdout",1);Lc("/dev/stderr",1)}yc=!1;Ca(sa);aa(k);k.onRuntimeInitialized?.();if(k.postRun)for("function"==typeof k.postRun&&(k.postRun=[k.postRun]);k.postRun.length;){var d=k.postRun.shift();ta.unshift(d)}Ca(ta)}}
if(!(0<E)){if(k.preRun)for("function"==typeof k.preRun&&(k.preRun=[k.preRun]);k.preRun.length;)ua();Ca(ra);0<E||(k.setStatus?(k.setStatus("Running..."),setTimeout(()=>{setTimeout(()=>k.setStatus(""),1);a()},1)):a())}}if(k.preInit)for("function"==typeof k.preInit&&(k.preInit=[k.preInit]);0<k.preInit.length;)k.preInit.pop()();Vc();moduleRtn=ca;


  return moduleRtn;
}
);
})();
export default rsiscool;
