/**
 * nmr-auto-assignment - Automatic assignment for Nuclear Magnetic Resonance spectra for small molecules 
 * @version v0.2.2
 * @link https://github.com/cheminfo-js/nmr-auto-assignment
 * @license MIT
 */
(function webpackUniversalModuleDefinition(root, factory) {
	if(typeof exports === 'object' && typeof module === 'object')
		module.exports = factory();
	else if(typeof define === 'function' && define.amd)
		define([], factory);
	else if(typeof exports === 'object')
		exports["nmrAutoAssignment"] = factory();
	else
		root["nmrAutoAssignment"] = factory();
})(typeof self !== 'undefined' ? self : this, function() {
return /******/ (function(modules) { // webpackBootstrap
/******/ 	// The module cache
/******/ 	var installedModules = {};
/******/
/******/ 	// The require function
/******/ 	function __webpack_require__(moduleId) {
/******/
/******/ 		// Check if module is in cache
/******/ 		if(installedModules[moduleId]) {
/******/ 			return installedModules[moduleId].exports;
/******/ 		}
/******/ 		// Create a new module (and put it into the cache)
/******/ 		var module = installedModules[moduleId] = {
/******/ 			i: moduleId,
/******/ 			l: false,
/******/ 			exports: {}
/******/ 		};
/******/
/******/ 		// Execute the module function
/******/ 		modules[moduleId].call(module.exports, module, module.exports, __webpack_require__);
/******/
/******/ 		// Flag the module as loaded
/******/ 		module.l = true;
/******/
/******/ 		// Return the exports of the module
/******/ 		return module.exports;
/******/ 	}
/******/
/******/
/******/ 	// expose the modules object (__webpack_modules__)
/******/ 	__webpack_require__.m = modules;
/******/
/******/ 	// expose the module cache
/******/ 	__webpack_require__.c = installedModules;
/******/
/******/ 	// define getter function for harmony exports
/******/ 	__webpack_require__.d = function(exports, name, getter) {
/******/ 		if(!__webpack_require__.o(exports, name)) {
/******/ 			Object.defineProperty(exports, name, { enumerable: true, get: getter });
/******/ 		}
/******/ 	};
/******/
/******/ 	// define __esModule on exports
/******/ 	__webpack_require__.r = function(exports) {
/******/ 		if(typeof Symbol !== 'undefined' && Symbol.toStringTag) {
/******/ 			Object.defineProperty(exports, Symbol.toStringTag, { value: 'Module' });
/******/ 		}
/******/ 		Object.defineProperty(exports, '__esModule', { value: true });
/******/ 	};
/******/
/******/ 	// create a fake namespace object
/******/ 	// mode & 1: value is a module id, require it
/******/ 	// mode & 2: merge all properties of value into the ns
/******/ 	// mode & 4: return value when already ns object
/******/ 	// mode & 8|1: behave like require
/******/ 	__webpack_require__.t = function(value, mode) {
/******/ 		if(mode & 1) value = __webpack_require__(value);
/******/ 		if(mode & 8) return value;
/******/ 		if((mode & 4) && typeof value === 'object' && value && value.__esModule) return value;
/******/ 		var ns = Object.create(null);
/******/ 		__webpack_require__.r(ns);
/******/ 		Object.defineProperty(ns, 'default', { enumerable: true, value: value });
/******/ 		if(mode & 2 && typeof value != 'string') for(var key in value) __webpack_require__.d(ns, key, function(key) { return value[key]; }.bind(null, key));
/******/ 		return ns;
/******/ 	};
/******/
/******/ 	// getDefaultExport function for compatibility with non-harmony modules
/******/ 	__webpack_require__.n = function(module) {
/******/ 		var getter = module && module.__esModule ?
/******/ 			function getDefault() { return module['default']; } :
/******/ 			function getModuleExports() { return module; };
/******/ 		__webpack_require__.d(getter, 'a', getter);
/******/ 		return getter;
/******/ 	};
/******/
/******/ 	// Object.prototype.hasOwnProperty.call
/******/ 	__webpack_require__.o = function(object, property) { return Object.prototype.hasOwnProperty.call(object, property); };
/******/
/******/ 	// __webpack_public_path__
/******/ 	__webpack_require__.p = "";
/******/
/******/
/******/ 	// Load entry module and return exports
/******/ 	return __webpack_require__(__webpack_require__.s = 13);
/******/ })
/************************************************************************/
/******/ ([
/* 0 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


__webpack_require__(36);

var abstractMatrix = __webpack_require__(8);

var util = __webpack_require__(3);

class Matrix extends abstractMatrix(Array) {
  constructor(nRows, nColumns) {
    var i;

    if (arguments.length === 1 && typeof nRows === 'number') {
      return new Array(nRows);
    }

    if (Matrix.isMatrix(nRows)) {
      return nRows.clone();
    } else if (Number.isInteger(nRows) && nRows > 0) {
      // Create an empty matrix
      super(nRows);

      if (Number.isInteger(nColumns) && nColumns > 0) {
        for (i = 0; i < nRows; i++) {
          this[i] = new Array(nColumns);
        }
      } else {
        throw new TypeError('nColumns must be a positive integer');
      }
    } else if (Array.isArray(nRows)) {
      // Copy the values from the 2D array
      const matrix = nRows;
      nRows = matrix.length;
      nColumns = matrix[0].length;

      if (typeof nColumns !== 'number' || nColumns === 0) {
        throw new TypeError('Data must be a 2D array with at least one element');
      }

      super(nRows);

      for (i = 0; i < nRows; i++) {
        if (matrix[i].length !== nColumns) {
          throw new RangeError('Inconsistent array dimensions');
        }

        this[i] = [].concat(matrix[i]);
      }
    } else {
      throw new TypeError('First argument must be a positive number or an array');
    }

    this.rows = nRows;
    this.columns = nColumns;
    return this;
  }

  set(rowIndex, columnIndex, value) {
    this[rowIndex][columnIndex] = value;
    return this;
  }

  get(rowIndex, columnIndex) {
    return this[rowIndex][columnIndex];
  }
  /**
   * Creates an exact and independent copy of the matrix
   * @return {Matrix}
   */


  clone() {
    var newMatrix = new this.constructor[Symbol.species](this.rows, this.columns);

    for (var row = 0; row < this.rows; row++) {
      for (var column = 0; column < this.columns; column++) {
        newMatrix.set(row, column, this.get(row, column));
      }
    }

    return newMatrix;
  }
  /**
   * Removes a row from the given index
   * @param {number} index - Row index
   * @return {Matrix} this
   */


  removeRow(index) {
    util.checkRowIndex(this, index);

    if (this.rows === 1) {
      throw new RangeError('A matrix cannot have less than one row');
    }

    this.splice(index, 1);
    this.rows -= 1;
    return this;
  }
  /**
   * Adds a row at the given index
   * @param {number} [index = this.rows] - Row index
   * @param {Array|Matrix} array - Array or vector
   * @return {Matrix} this
   */


  addRow(index, array) {
    if (array === undefined) {
      array = index;
      index = this.rows;
    }

    util.checkRowIndex(this, index, true);
    array = util.checkRowVector(this, array, true);
    this.splice(index, 0, array);
    this.rows += 1;
    return this;
  }
  /**
   * Removes a column from the given index
   * @param {number} index - Column index
   * @return {Matrix} this
   */


  removeColumn(index) {
    util.checkColumnIndex(this, index);

    if (this.columns === 1) {
      throw new RangeError('A matrix cannot have less than one column');
    }

    for (var i = 0; i < this.rows; i++) {
      this[i].splice(index, 1);
    }

    this.columns -= 1;
    return this;
  }
  /**
   * Adds a column at the given index
   * @param {number} [index = this.columns] - Column index
   * @param {Array|Matrix} array - Array or vector
   * @return {Matrix} this
   */


  addColumn(index, array) {
    if (typeof array === 'undefined') {
      array = index;
      index = this.columns;
    }

    util.checkColumnIndex(this, index, true);
    array = util.checkColumnVector(this, array);

    for (var i = 0; i < this.rows; i++) {
      this[i].splice(index, 0, array[i]);
    }

    this.columns += 1;
    return this;
  }

}

exports.Matrix = Matrix;
Matrix.abstractMatrix = abstractMatrix;

/***/ }),
/* 1 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


var abstractMatrix = __webpack_require__(8);

var Matrix = __webpack_require__(0);

class BaseView extends abstractMatrix() {
  constructor(matrix, rows, columns) {
    super();
    this.matrix = matrix;
    this.rows = rows;
    this.columns = columns;
  }

  static get [Symbol.species]() {
    return Matrix.Matrix;
  }

}

module.exports = BaseView;

/***/ }),
/* 2 */
/***/ (function(module, exports, __webpack_require__) {

/* WEBPACK VAR INJECTION */(function(global) {/**
 * openchemlib - Manipulate molecules
 * @version v4.4.0
 * @date 2016-09-01T11:25:09.510Z
 * @link https://github.com/cheminfo/openchemlib-js
 * @license BSD-3-Clause
*/(function(root){'use strict';function getExports($wnd){var $doc=$wnd.document;var $gwt={};var navigator={userAgent:'webkit'};function noop(){}var __gwtModuleFunction=noop;__gwtModuleFunction.__moduleStartupDone=noop;var $sendStats=noop;var $moduleName,$moduleBase;// Start GWT code 
function IG(){}function FG(){}function Fd(){}function nc(){}function oB(){}function vB(){}function tf(){}function wf(){}function yf(){}function Bf(){}function ag(){}function pg(){}function mh(){}function _m(){}function _q(){}function Xn(){}function Eo(){}function No(){}function Fz(){}function cA(){}function dA(){}function nA(){}function eB(){}function AB(){}function TB(){}function hH(){}function sH(){}function yH(){}function yO(){}function EN(){}function Tg(){this.g=1;}function QI(a){this.a=a;}function Lz(a){this.a=a;}function FH(a){this.a=a;}function aJ(a){this.a=a;}function CL(a){this.a=a;}function IL(a){this.a=a;}function ML(a){this.a=a;}function RL(a){this.a=a;}function mM(a){this.a=a;}function sM(a){this.a=a;}function hM(a){this.b=a;}function om(a){this.b=a;}function UK(a){this.c=a;}function dN(a){this.c=a;}function EP(a){this.a=a;}function wh(){th(this);}function pJ(){GA.call(this);}function pk(a,b){a.M=b;}function kk(a,b){a.G=b;}function mk(a,b){a.P=b;}function uo(a,b){a.i=b;}function Zi(a,b){a.Q&=~b;}function uj(a,b){a.A[b]=-1;}function vj(a,b){a.F[b]=128;}function hk(a,b){a.C[b]|=GR;}function qk(a,b){a.s[b]|=qR;}function wN(a,b){a.sort(b);}function hz(a){return ho(a);}function jz(a){return io(a);}function kz(a){return jo(a);}function Cz(){Bz();new $q();}function bz(){this.a=new Vn();}function gz(){this.a=new go();}function oq(){this.a=new LM();}function sA(){this.a=new LM();}function _O(){this.a=new LM();}function Uo(){this.c=new LM();}function IH(){GA.call(this);}function MI(){GA.call(this);}function OI(){GA.call(this);}function eK(){GA.call(this);}function HN(){GA.call(this);}function KH(){IH.call(this);}function MH(){MH=FG;LH=false;}function bB(){bB=FG;aB=new eB();}function yB(){yB=FG;xB=new AB();}function KA(){KA=FG;JA=new nc();}function cK(){cK=FG;bK=new yH();}function DN(){DN=FG;CN=new EN();}function kN(a){nN(a,a.length);}function mN(a){pN(a,a.length);}function jB(a){iB();hB.pb(a);}function UP(a){OP(a);return a;}function We(a,b){return a.k[b];}function af(a,b){return a.W[b];}function Pf(a,b){return a.e[b];}function ji(a,b){return a.q[b];}function ti(a,b){return a.v[b];}function vi(a,b){return a.w[b];}function Ai(a,b){return a.A[b];}function Oi(a,b){return a.D[b];}function Pi(a,b){return a.F[b];}function Qk(a,b){return a.c[b];}function Tk(a,b){return a.k[b];}function bl(a,b){return a.g[b];}function jn(a,b){return a.a[b];}function kn(a,b){return a.b[b];}function on(a,b){return a.d[b];}function rn(a,b){return a.e[b];}function kJ(a){return a<0?-a:a;}function lz(){return $n(),Yn;}function WN(a){a.b=null;a.c=0;}function Fj(a,b){a.o=b;a.Q=0;}function Gj(a,b){a.p=b;a.Q=0;}function HA(a){FA.call(this,a);}function qJ(a){IA.call(this,a);}function rJ(a){HA.call(this,a);}function oJ(a){HA.call(this,a);}function JH(a){HA.call(this,a);}function NI(a){HA.call(this,a);}function fK(a){HA.call(this,a);}function aK(a){JH.call(this,a);}function nK(a){mK.call(this,a);}function tJ(a){NI.call(this,a);}function FP(a){EP.call(this,a);}function SJ(){FH.call(this,'');}function ZJ(){FH.call(this,'');}function $J(){FH.call(this,'');}function yG(){wG==null&&(wG=[]);}function TA(){TA=FG;!!(iB(),hB);}function aQ(){aQ=FG;ZP={};_P={};}function bI(a){aI(a);return a.k;}function _N(a){return!!a&&a.b;}function Xe(a){return Ye(a,a._);}function VB(a,b){return jI(a,b);}function mz(a,b){return ko(a,b);}function nz(a,b){return lo(a,b);}function lJ(a,b){return a>b?a:b;}function mJ(a,b){return a<b?a:b;}function yi(a,b){return a.H[b].b;}function xi(a,b){return a.H[b].a;}function zi(a,b){return a.H[b].c;}function ui(a,b){return a.s[b]&3;}function Ni(a,b){return a.C[b]&3;}function cj(a,b){return a.u[b]<0;}function gp(a,b){a.a=b;a.Q&=-136;}function Xl(a,b){a.Q|=252&(4|b);}function SO(a){this.a=new iO(a);}function AH(a){this.b=a;this.a=0;}function qf(a){rf.call(this,a,0);}function qu(){ru.call(this,null);}function hO(){iO.call(this,null);}function WG(a){TG();this.c=oR|a;}function er(){dr();this.a=new hq();}function Nv(a){ou();return Ck(a);}function gx(a){ou();return Dk(a);}function hx(a){ou();return Ek(a);}function UJ(a,b){a.a+=b;return a;}function BP(a,b,c){a.splice(b,c);}function Jj(a,b,c){a.q[b]=c;a.Q=0;}function jk(a,b,c){a.F[b]=c;a.Q=0;}function th(a){a.a=0;a.b=0;a.c=1;}function lN(a){oN(a,a.length,-1);}function tl(a){Xo(a,3);return a.n;}function ln(a,b){return BM(a.g,b);}function mn(a,b){return BM(a.i,b);}function Yo(a,b){return af(a.b,b);}function Zo(a,b){return We(a.b,b);}function fp(a,b){return _e(a.b,b);}function lG(a,b){return gG(a,b)<0;}function Jn(a){return Fn(a,a.b)>0;}function BH(a,b){return yJ(a.a,b);}function kK(a,b){return GB(a.a,b);}function lK(a,b){return GB(a.a,b);}function bM(a,b){return!!XN(a,b);}function pf(a,b){return a<b?a:a-b;}function vh(a,b){return b*a.c+a.b;}function uh(a,b){return b*a.c+a.a;}function RO(a,b){return bM(a.a,b);}function ZO(a,b){return yM(a.a,b);}function $O(a,b){return KM(a.a,b);}function si(a,b){return kJ(a.u[b]);}function pi(a,b){return Ch[a.A[b]];}function wi(a,b){return a.s[b]&48;}function ki(a,b){return a.s[b]&448;}function Qi(a,b){return a.F[b]&103;}function nj(a,b){return a.v[b]==0;}function _i(a,b){return a.A[b]==-1;}function jj(a,b){return Dk(a.A[b]);}function kj(a,b){return Ek(a.A[b]);}function UB(a){return a[4]||a[1];}function oO(a){return a.b=SK(a.a);}function CC(a){return a.l|a.m<<22;}function iG(a,b){return gG(a,b)==0;}function nG(a,b){return gG(a,b)!=0;}function pm(a,b){this.a=a;this.b=b;}function lq(a,b){this.a=a;this.b=b;}function mq(a,b){this.a=a;this.b=b;}function iH(a,b){this.a=a;this.b=b;}function FO(a,b){this.a=a;this.b=b;}function tm(a,b){this.b=a;this.a=b;}function _L(a,b){this.c=a;this.d=b;}function xH(a){this.b=a;this.a=-2;}function _A(){QA!=0&&(QA=0);SA=-1;}function QC(a){return typeof a===lQ;}function RC(a){return typeof a===mQ;}function UC(a){return typeof a===nQ;}function DJ(a,b){return OP(a),a===b;}function jq(a,b){return BM(a.a,b).b;}function iq(a,b){return BM(a.a,b).a;}function Ei(a,b,c){return a.B[b][c];}function al(a,b,c){return a.f[b][c];}function cl(a,b,c){return a.i[b][c];}function dl(a,b,c){return a.j[b][c];}function gj(a,b){return a.F[b]==128;}function _I(a,b){return bJ(a.a,b.a);}function cM(a,b){return mL(XN(a,b));}function kG(a){return typeof a===mQ;}function Hz(a){Am();Em.call(this,a);}function JP(){EP.call(this,'UTF-8');}function $A(a){$wnd.clearTimeout(a);}function FJ(a,b){return a.indexOf(b);}function VC(a){return a==null?null:a;}function VJ(a,b){a.a+=''+b;return a;}function RJ(a,b){a.a+=''+b;return a;}function WJ(a,b){a.a+=''+b;return a;}function nv(a,b){ou();return Bk(a,b);}function zN(a,b){uN(a,0,a.length,b);}function yN(a){uN(a,0,a.length,null);}function RK(a){return a.a<a.c.size();}function mL(a){return!a?null:a.Gb();}function yI(a){return DJ(mQ,typeof a);}function HJ(a){return DJ(nQ,typeof a);}function tA(a){a.g=ZB(gF,iQ,43,0,0,1);}function JO(){FO.call(this,'Head',1);}function OO(){FO.call(this,'Tail',3);}function Tj(a,b,c){a.v[b]=c;a.Q&=3;}function Zj(a,b,c){a.H[b].a=c;a.Q&=3;}function $j(a,b,c){a.H[b].b=c;a.Q&=3;}function _j(a,b,c){a.H[b].c=c;a.Q&=3;}function Rj(a,b,c,d){a.u[b]=d?-c:c;}function zP(a,b,c){a.splice(b,0,c);}function tH(a,b,c,d){rH(this,a,b,c,d);}function Hm(a){Im.call(this,a,new $J());}function LO(){FO.call(this,'Range',2);}function go(){$n();this.g=new Vn();ao();}function eC(a){return fC(a.l,a.m,a.h);}function Bi(a){return Ci(a,a.o,a.p,Fh);}function PJ(a){return QJ(a,0,a.length);}function JN(a){return a!=null?tc(a):0;}function RI(a,b){return a<b?-1:a>b?1:0;}function aj(a,b){return(a.s[b]&4)!=0;}function hj(a,b){return(a.C[b]&4)!=0;}function dj(a,b){return(a.C[b]&qR)!=0;}function fj(a,b){return(a.C[b]&cR)!=0;}function ij(a,b){return(a.C[b]&GR)!=0;}function bj(a,b){return(a.s[b]&FR)!=0;}function lj(a,b){return(a.s[b]&cR)!=0;}function Wi(a,b){return(a.s[b]&qR)!=0;}function $i(a,b){return(a.s[b]&ER)!=0;}function Dl(a,b){return(a.s[b]&iR)!=0;}function Kl(a,b){return(a.s[b]&BR)!=0;}function Nl(a,b){return(a.s[b]&8)!=0;}function El(a,b){return(a.s[b]&xQ)!=0;}function Pl(a,b){return(a.s[b]&yQ)!=0;}function ej(a,b){return(a.D[b]&$Q)!=0;}function Ll(a,b){return(a.C[b]&64)!=0;}function Hi(a,b){return(a.C[b]&48)>>4;}function dG(a){return a.backingJsObject;}function yJ(a,b){return a.charCodeAt(b);}function nq(a,b,c){yM(a.a,new mq(b,c));}function QN(a,b,c){a.a=b^1502;a.b=c^jU;}function kh(a,b,c,d){a.a=b;a.b=c;a.c=d;}function bk(a,b,c,d){a.B[b][c]=d;a.Q=0;}function Pe(a){a.q=new ZJ();a.p=6;a.r=0;}function ai(a,b){a.F[b]=128;Uh(a);a.Q=0;}function VP(a,b){return a==b?0:a<b?-1:1;}function oi(a,b){return(a.s[b]&zR)>>19;}function Ji(a,b){return(a.C[b]&BR)>>10;}function Gi(a,b){return(a.D[b]&960)>>6;}function qj(a,b){return(a.s[b]&512)!=0;}function Fl(a,b){return(a.C[b]&256)!=0;}function Hl(a,b){return(a.C[b]&512)!=0;}function Ol(a,b){return(a.C[b]&128)!=0;}function Ti(a,b){return Ui(a,b)+Si(a,b);}function zo(a,b){Ic();yo.call(this,a,b);}function _J(a){FH.call(this,(OP(a),a));}function lp(a,b){Gh();zk.call(this,a,b);}function Wo(a,b){Yh(a,b);!!a.b&&(b.Q=0);}function hp(a){Xo(a,15);!!a.b&&of(a.b);}function Xm(a,b,c){a.c=b;return Zm(a,c);}function fC(a,b,c){return{l:a,m:b,h:c};}function oC(a){return a.l+a.m*XT+a.h*YT;}function YP(a){return a.$H||(a.$H=++XP);}function IJ(a,b){return a.lastIndexOf(b);}function GJ(a,b,c){return a.indexOf(b,c);}function PC(a,b){return a!=null&&NC(a,b);}function Df(a,b){a.d[a.c]=eG(a.d[a.c],b);}function Gf(a){this.d=ZB(aD,gR,6,a,14,1);}function LM(){this.a=ZB(eF,fR,1,0,5,1);}function vO(a){this.a=a;hM.call(this,a);}function qO(a){rO.call(this,a,(EO(),AO));}function BB(){BB=FG;zB((yB(),yB(),xB));}function GA(){tA(this);uA(this);this.nb();}function aI(a){if(a.k!=null){return;}nI(a);}function KP(a){if(!a){throw dG(new MI());}}function SP(a){if(!a){throw dG(new OI());}}function MP(a){if(!a){throw dG(new HN());}}function RH(a){if(!a){throw dG(new MI());}}function mB(a){iB();return parseInt(a)||-1;}function NH(a){MH();return DJ(lQ,typeof a);}function OH(a,b){MH();return a==b?0:a?1:-1;}function Kc(a,b){return a==null?b:a+','+b;}function fA(a,b){return{type:b,value:a};}function BJ(a,b){return OP(a),a+(OP(b),b);}function JM(a){return xP(a.a,a.a.length);}function XH(a){return a>=56320&&a<=57343;}function Ym(a){return a.length==0?0:xI(a);}function LJ(a,b,c){return a.substr(b,c-b);}function ii(a,b){return(a.s[b]&98304)>>15;}function _e(a,b){return a.d==null?-1:a.d[b];}function Ff(a,b){a.a=b;a.c=0;a.b=63;mN(a.d);}function dk(a,b,c){a.C[b]&=-49;a.C[b]|=c<<4;}function Kj(a,b,c){a.s[b]&=-449;a.s[b]|=c;}function WP(b,c,d){try{b[c]=d;}catch(a){}}function ar(a,b,c,d,e){return bq(a,b,c,d,e);}function JJ(a,b,c){return a.lastIndexOf(b,c);}function TC(a,b){return a&&b&&a instanceof b;}function hi(a,b){return((a.s[b]&yR)>>>28)-1;}function lh(a,b){kh(a,b.a,b.b,b.c);return a;}function VL(a,b){var c;c=a.d;a.d=b;return c;}function Xz(a){!a.a&&(a.a=new Tg());return a.a;}function Zz(a){!a.d&&(a.d=new _m());return a.d;}function $z(a){!a.e&&(a.e=new go());return a.e;}function _z(a){!a.f&&(a.f=new Eo());return a.f;}function aA(a){!a.g&&(a.g=new No());return a.g;}function bA(a){!a.i&&(a.i=new dA());return a.i;}function Vy(a,b){this.a=new Vp(new AH(a),b);}function eH(a){this.a='Helvetica';this.b=a;}function pO(a){TK(a.a);cO(a.c,a.b);a.b=null;}function zB(a){!a.a&&(a.a=new TB());return a.a;}function UA(a,b,c){return a.apply(b,c);var d;}function DM(a,b){return EM(a,b,a.a.length-1);}function mv(a,b,c,d){ou();return Ak(a,b,c,d);}function Oc(a,b,c){qo(a,b-a.O/2,c-a.O/2,a.O);}function Ed(a,b,c){this.b=a;this.c=b;this.a=c;}function RB(a,b){BB();QB.call(this,a,b,true);}function HH(){HA.call(this,'divide by zero');}function JG(){$wnd.setTimeout(eQ(LG));KG();}function kp(){Gh();this.K=this.L=256;Yi(this);}function $I(){$I=FG;ZI=ZB(ZE,eU,31,256,0,1);}function jJ(){jJ=FG;iJ=ZB(_E,eU,45,256,0,1);}function dr(){dr=FG;cr=(!Wz&&(Wz=new cA()),Wz);}function fI(a){var b;b=eI(a);rI(a,b);return b;}function Tp(a,b){a.c==null&&Pp(a,b);return a.c;}function yM(a,b){a.a[a.a.length]=b;return true;}function KJ(a,b){return a.substr(b,a.length-b);}function qi(a,b){return a.t==null?null:a.t[b];}function mi(a,b){return a.r==null?null:a.r[b];}function sj(a,b){return a.F[b]==17||a.F[b]==9;}function Rk(a,b){return a.c[b]-a.g[b]+ml(a,b);}function jl(a,b){return Ui(a,b)+Si(a,b)-ol(a,b);}function kH(a,b){return nH(a,b.c,b.d,b.b,b.a);}function zJ(a,b){return VP((OP(a),a),(OP(b),b));}function pG(a,b){return hG(xC(kG(a)?rG(a):a,b));}function qG(a,b){return hG(yC(kG(a)?rG(a):a,b));}function nJ(a){return a==0||isNaN(a)?a:a<0?-1:1;}function cp(a){Xo(a,7);return!a.b?null:Ze(a.b);}function dp(a){Xo(a,7);return!a.b?null:Xe(a.b);}function lg(a,b){return b==1?a.a+a.f++:a.i+a.g++;}function rh(a,b){b.a=b.a*a.c+a.a;b.b=b.b*a.c+a.b;}function rH(a,b,c,d,e){a.c=b;a.d=c;a.b=d;a.a=e;}function Wj(a,b,c){a.s[b]&=-49;a.s[b]|=c;a.Q&=3;}function Ij(a,b,c){a.s[b]&=-98305;a.s[b]|=c<<15;}function XG(a,b,c){TG();YG.call(this,a,b,c,255);}function xM(a,b,c){QP(b,a.a.length);zP(a.a,b,c);}function BM(a,b){NP(b,a.a.length);return a.a[b];}function AP(a,b,c){yP(c,0,a,b,c.length,false);}function vo(a,b){if(a.j!=b){a.j=b;a.e=new eH(b);}}function eP(a,b){if(a<0||a>=b){throw dG(new KH());}}function tG(a){if(kG(a)){return a|0;}return CC(a);}function xN(c){c.sort(function(a,b){return a-b;});}function Sp(a){a.c==null&&Pp(a,10240);return a.c;}function hI(a){var b;b=eI(a);b.j=a;b.e=1;return b;}function pN(a,b){var c;for(c=0;c<b;++c){a[c]=0;}}function nN(a,b){var c;for(c=0;c<b;++c){a[c]=-1;}}function oN(a,b,c){var d;for(d=0;d<b;++d){a[d]=c;}}function Xj(a,b,c){c?a.s[b]|=512:a.s[b]&=-513;}function XB(a,b,c,d,e,f){return YB(a,b,c,d,e,0,f);}function uG(a){if(kG(a)){return''+a;}return DC(a);}function PP(a,b){if(a==null){throw dG(new rJ(b));}}function Xk(a,b){return!!a.n&&b<a.d?jn(a.n,b):0;}function _k(a,b){return!!a.n&&b<a.e?kn(a.n,b):0;}function bJ(a,b){return gG(a,b)<0?-1:gG(a,b)>0?1:0;}function QO(a,b){return aO(a.a,b,(MH(),LH))==null;}function _B(a){return Array.isArray(a)&&a.Lb===IG;}function OC(a){return!Array.isArray(a)&&a.Lb===IG;}function iO(a){this.b=null;this.a=(DN(),!a?CN:a);}function FA(a){tA(this);this.f=a;uA(this);this.nb();}function ru(a){ou();!a&&(a=new lp(32,32));this.a=a;}function fB(a,b){!a&&(a=[]);a[a.length]=b;return a;}function xo(a,b){RJ(a.c,'\t');RJ(a.c,b);RJ(a.c,kQ);}function Sj(a,b,c){c?a.s[b]|=cR:a.s[b]&=-262145;}function fk(a,b,c){c?a.C[b]|=cR:a.C[b]&=-262145;}function ck(a,b,c){c?a.C[b]|=qR:a.C[b]&=-131073;}function Aj(a){var b;for(b=0;b<a.o;b++)a.s[b]&=-513;}function yj(a){var b;for(b=0;b<a.o;b++)a.s[b]&=-449;}function qN(a,b){var c;for(c=0;c<b;++c){a[c]=false;}}function bO(a,b){var c;c=new yO();dO(a,b,c);return c.d;}function xP(a,b){var c;c=a.slice(0,b);return bC(c,a);}function UH(a){var b;b=a-10;return(b<0?48+a:97+b)&AQ;}function ff(a){var b;b=0;while(a>0){a>>=1;++b;}return b;}function OP(a){if(a==null){throw dG(new pJ());}return a;}function TJ(a,b){a.a+=String.fromCharCode(b);return a;}function DH(a,b,c){CH(a,b,b+1,String.fromCharCode(c));}function Yj(a,b,c){a.s[b]&=-134217729;c&&(a.s[b]|=FR);}function Of(a,b){return a.f[b]&&(a.o[b]==1||a.o[b]==2);}function rC(a,b){return fC(a.l&b.l,a.m&b.m,a.h&b.h);}function wC(a,b){return fC(a.l|b.l,a.m|b.m,a.h|b.h);}function EC(a,b){return fC(a.l^b.l,a.m^b.m,a.h^b.h);}function zH(a){return a.a==a.b.length?-1:yJ(a.b,a.a++);}function CH(a,b,c,d){a.a=LJ(a.a,0,b)+(''+d)+KJ(a.a,c);}function hm(a,b,c){a.c=6;a.d=c;a.a=b;a.e=b[a.d]-64<<11;}function Bj(a){var b;for(b=0;b<a.p;b++)a.C[b]&=-393217;}function zj(a){var b;for(b=0;b<a.o;b++)a.s[b]&=-262145;}function _o(a){var b;b=new lp(a.o,a.p);Xh(a,b);return b;}function Ye(a,b){if(a.e==null){Ue(a);Qe(a,b);}return a.e;}function Rp(a,b){if(a.b==null)return null;return a.b[b];}function Se(a,b){if(!b){Ne(a,1,1);Ne(a,15,4);}return true;}function Oe(a){a.r<<=a.p;TJ(a.q,a.r+64&AQ);return a.q.a;}function uJ(a,b,c){this.a=jQ;this.d=a;this.b=b;this.c=c;}function zh(a,b,c,d){this.b=a;this.a=b;this.c=c;this.d=d;}function Qd(a,b){this.d=a;this.c=b;Xo(this.d,1);Gd(this);}function zk(a,b){this.K=1>a?1:a;this.L=1>b?1:b;Yi(this);}function gI(a,b){var c;c=eI(a);rI(a,c);c.e=b?8:0;return c;}function xf(a,b){if(a.c!=b.c)return a.c>b.c?1:-1;return 0;}function XJ(a,b,c){a.a=LJ(a.a,0,b)+''+KJ(a.a,c);return a;}function Uj(a,b,c,d){a.s[b]&=-8;a.s[b]|=c;d&&(a.s[b]|=4);}function qh(a,b){b.c*=a.c;b.a=b.a*a.c+a.a;b.b=b.b*a.c+a.b;}function YK(a,b){UK.call(this,a);QP(b,a.size());this.a=b;}function pp(){Uo.call(this);this.b=new LM();this.a=new LM();}function iB(){iB=FG;var a,b;b=!nB();a=new vB();hB=b?new oB():a;}function Bz(){Bz=FG;Az=(!Wz&&(Wz=new cA()),Wz);tz=(Vq(),Eq);}function dQ(){if($P==256){ZP=_P;_P={};$P=0;}++$P;}function IN(a,b){return VC(a)===VC(b)||a!=null&&pc(a,b);}function xA(a,b){a.backingJsObject=b;b!=null&&WP(b,pQ,a);}function yA(a,b){var c;c=bI(a.Jb);return b==null?c:c+': '+b;}function hN(a,b){LP(b);return jN(a,ZB(_C,DQ,6,b,15,1),b);}function CJ(a){var b;return PJ(GP(a,0,(b=a.length,DP(),b)));}function NJ(a){return String.fromCharCode.apply(null,a);}function AJ(a,b){return zJ(a.toLowerCase(),b.toLowerCase());}function gL(a,b){return b===a?'(this Map)':b==null?rQ:vc(b);}function fG(a,b){return hG(rC(kG(a)?rG(a):a,kG(b)?rG(b):b));}function oG(a,b){return hG(wC(kG(a)?rG(a):a,kG(b)?rG(b):b));}function vG(a,b){return hG(EC(kG(a)?rG(a):a,kG(b)?rG(b):b));}function vI(a){return DJ(mQ,typeof a)||a instanceof Number;}function mc(a){return bI(rc(a))+'@'+(tc(a)>>>0).toString(16);}function WC(a){return Math.max(Math.min(a,oQ),-2147483648)|0;}function ZA(a){TA();$wnd.setTimeout(function(){throw a;},0);}function Uc(a){var b,c;b=(TG(),PG);c=new WG(a);return kA(c,b);}function pu(a){var b;b=Xz(nu);b.j=new SN(0);Ig(b,a.a);$l(a.a);}function md(a){var b;b=a.a;a.a=a.b;a.b=b;b=a.c;a.c=a.d;a.d=b;}function td(a,b){var c;c=a.b;a.b=b.b;b.b=c;c=a.d;a.d=b.d;b.d=c;}function jI(a,b){var c=a.a=a.a||[];return c[b]||(c[b]=a.sb(b));}function HM(a,b,c){var d;RP(b,c,a.a.length);d=c-b;BP(a.a,b,d);}function nf(a){var b;for(b=0;b<a.L.d;b++){Yj(a.L,b,a.H[b]);}}function YJ(a,b,c){a.a=LJ(a.a,0,b)+(''+c)+KJ(a.a,b);return a;}function Hh(a,b,c,d){var e;e=Ih(a,6);kh(a.H[e],b,c,d);return e;}function Vj(a,b,c,d){d?a.w[b]|=c:a.w[b]&=~c;a.Q=0;a.I=true;}function ik(a,b,c,d){d?a.D[b]|=c:a.D[b]&=~c;a.Q=0;a.I=true;}function Lj(a,b,c){c?a.s[b]|=ER:a.s[b]&=-67108865;a.Q&=3;}function tj(a,b,c){return(a.F[b]==17||a.F[b]==9)&&a.B[0][b]==c;}function rj(a,b){return(a.s[a.B[0][b]]&a.s[a.B[1][b]]&512)!=0;}function Fi(a,b){return((a.D[b]&960)>>6)+((a.D[b]&15360)>>10);}function IO(){EO();return aC(VB(UF,1),eU,42,0,[AO,BO,CO,DO]);}function SN(a){NN();QN(this,tG(fG(qG(a,24),EQ)),tG(fG(a,EQ)));}function MA(a){KA();IA.call(this,a);this.a='';this.b=a;this.a='';}function MG(a,b){tA(this);this.e=b;this.f=a;uA(this);this.nb();}function eo(a,b,c){a.b=null;a.a=b;c==null?a.c=_n(a,b):a.c=c;}function fo(a,b,c){a.e=null;a.d=b;c==null?a.f=_n(a,b):a.f=c;}function Di(a,b,c){return Ak(a.H[b].a,a.H[b].b,a.H[c].a,a.H[c].b);}function gk(a,b,c,d){a.C[b]&=-16777224;a.C[b]|=c;d&&(a.C[b]|=4);}function fO(a,b){var c;c=1-b;a.a[c]=gO(a.a[c],c);return gO(a,b);}function gA(a,b){var c;c=a-b;c>=MQ?c-=LQ:c<tR&&(c+=LQ);return c;}function LG(){var a;a=OG();if(!DJ(_T,a)){throw dG(new NG(a));}}function sG(a){var b;if(kG(a)){b=a;return b==-0.?0:b;}return BC(a);}function GG(a){function b(){};b.prototype=a||{};return new b();}function AA(b){if(!('stack'in b)){try{throw b;}catch(a){}}return b;}function VH(a){return null!=String.fromCharCode(a).match(/\d/);}function WH(a){return null!=String.fromCharCode(a).match(/[A-Z]/i);}function cN(a){MP(a.a<a.c.a.length);a.b=a.a++;return a.c.a[a.b];}function SK(a){MP(a.a<a.c.size());return a.c.getAtIndex(a.b=a.a++);}function TK(a){SP(a.b!=-1);a.c.removeAtIndex(a.b);a.a=a.b;a.b=-1;}function $K(a,b,c){RP(b,c,a.size());this.c=a;this.a=b;this.b=c-b;}function sh(a,b){b.c=b.c*a.c+a.a;b.d=b.d*a.c+a.b;b.b*=a.c;b.a*=a.c;}function Ej(a,b){var c;for(c=0;c<a.o;c++){a.H[c].a*=b;a.H[c].b*=b;}}function Vk(a,b){var c;c=a.s[b]&BR;return c==0?0:c==OQ?2:c==YQ?3:4;}function XA(a,b,c){var d;d=VA();try{return UA(a,b,c);}finally{YA(d);}}function Km(a,b){a.c=null;return Zm(a,new xH(new AH(b)))?a.c:null;}function Uk(a,b){Xo(a,3);return a.k[b]==2&&a.g[b]==2?Rl(a,b):Tl(a,b);}function li(a,b){return a.r==null?null:a.r[b]==null?null:CJ(a.r[b]);}function lI(a){if(a.xb()){return null;}var b=a.j;var c=BG[b];return c;}function Pn(a,b){if(b.o==0){a.A=null;return;}a.A=b;a.F=false;Xo(a.A,1);}function sm(a){this.a=ZB(_C,DQ,6,a,15,1);this.b=ZB(_C,DQ,6,a,15,1);}function kf(a,b){var c;c=ZB(_C,DQ,6,b,15,1);dK(a,c,a.length);return c;}function Lm(a,b){var c;c=!a.a?null:cM(a.a,new QI(b));return!c?b-1:c.a;}function Mm(a,b){var c;c=!a.b?null:cM(a.b,new QI(b));return!c?b-1:c.a;}function qA(a,b){var c;c=rA(a,b);if(c<0){c=-(c+1);xM(a.a,c,b);}return c;}function vH(a){var b;if(a.a!=-2){b=a.a;a.a=-2;}else{b=zH(a.b);}return b;}function Og(a){var b,c;for(c=0;c<a.f.a.length;c++){b=BM(a.f,c);bh(b);}}function Jp(a){Ip();var b,c;Xo(a,1);c=a.d;for(b=0;b<c;b++){Kp(a,b);}}function dC(a){var b,c,d;b=a&OR;c=a>>22&OR;d=a<0?VT:0;return fC(b,c,d);}function aM(a,b){var c,d;c=b.Fb();d=XN(a,c);return!!d&&IN(d.d,b.Gb());}function EM(a,b,c){for(;c>=0;--c){if(IN(b,a.a[c])){return c;}}return-1;}function IM(a,b,c){var d;d=(NP(b,a.a.length),a.a[b]);a.a[b]=c;return d;}function fh(a,b,c){var d;for(d=0;d<a.a.length;d++){a.c[d]+=b;a.d[d]+=c;}}function wP(a,b,c,d){Array.prototype.splice.apply(a,[b,c].concat(d));}function LP(a){if(a<0){throw dG(new oJ('Negative array size: '+a));}}function mK(a){this.a=(BB(),new RB(a,['USD','US$',2,'US$','$']));}function YA(a){a&&dB((bB(),aB));--QA;if(a){if(SA!=-1){$A(SA);SA=-1;}}}function ve(a){var b,c;b=le(a);do{c=b;he(a);b=le(a);}while(c!=b);return b;}function Ag(a,b){var c;c=a-b;while(c<tR)c+=LQ;while(c>MQ)c-=LQ;return c;}function FM(a,b){var c;c=(NP(b,a.a.length),a.a[b]);BP(a.a,b,1);return c;}function TN(a,b){!a.a?a.a=new _J(a.d):WJ(a.a,a.b);VJ(a.a,b);return a;}function uA(a){if(a.j){a.backingJsObject!==uQ&&a.nb();a.g=null;}return a;}function Ml(a,b){return a.A[b]==1&&a.v[b]==0&&(a.r==null||a.r[b]==null);}function ZH(a){return String.fromCharCode(a).toLowerCase().charCodeAt(0);}function PA(){if(Date.now){return Date.now();}return new Date().getTime();}function CP(){if(Date.now){return Date.now();}return new Date().getTime();}function WA(b){TA();return function(){return XA(b,this,arguments);var a;};}function Ze(a){if(a.D==null){Ue(a);ef(a);gf(a,1);gf(a,2);df(a);}return a.D;}function mp(a){Gh();zk.call(this,!a?256:a.K,!a?256:a.L);!!a&&Xh(a,this);}function xO(a,b){_L.call(this,a,b);this.a=ZB(PF,fR,60,2,0,1);this.b=true;}function Vp(a,b){Np();this.c=b;this.g=new xH(a);this.f=new $J();this.a=new $J();}function UN(a,b){this.b=', ';this.d=a;this.e=b;this.c=this.d+(''+this.e);}function cO(a,b){var c;c=new yO();c.c=true;c.d=b.Gb();return dO(a,b.Fb(),c);}function jg(a,b){var c,d;d=a.j.k[b];c=a.j.j[b];return d==0?a.b:d==1?c:a.a+c;}function mH(a,b,c){var d,e;d=a.c;e=a.d;return b>=d&&c>=e&&b<d+a.b&&c<e+a.a;}function jN(a,b,c){var d,e;e=a.length;d=c<e?c:e;yP(a,0,b,0,d,true);return b;}function ZB(a,b,c,d,e,f){var g;g=$B(e,d);e!=10&&aC(VB(a,f),b,c,e,g);return g;}function uN(a,b,c,d){var e;d=(DN(),!d?CN:d);e=a.slice(b,c);vN(e,a,b,c,-b,d);}function Cj(a,b){var c;for(c=0;c<a.o;c++)kJ(a.u[c])==(b<0?-b:b)&&(a.u[c]=0);}function tp(a){var b;Xo(a,15);for(b=0;b<a.o;b++){(a.s[b]&3)!=0&&Oj(a,b,1,0);}}function GM(a,b){var c;c=CM(a,b,0);if(c==-1){return false;}FM(a,c);return true;}function Bk(a,b){Gh();var c;c=a-b;while(c<tR)c+=LQ;while(c>MQ)c-=LQ;return c;}function ev(a){ou();var b;b=new qu();Xm(Zz(nu),b.a,new xH(new AH(a)));return b;}function bC(a,b){WB(b)!=10&&aC(rc(b),b.Kb,b.__elementTypeId$,WB(b),a);return a;}function _g(a,b){var c;for(c=0;c<a.a.length;c++)if(b==a.a[c])return c;return-1;}function CM(a,b,c){for(;c<a.a.length;++c){if(IN(b,a.a[c])){return c;}}return-1;}function ol(a,b){var c,d;a.gb(1);d=0;for(c=0;c<a.c[b];c++)d+=a.j[b][c];return d;}function iN(a,b){var c,d;LP(b);return c=(d=a.slice(0,b),bC(d,a)),c.length=b,c;}function ap(a){var b,c;b=ZB(_C,DQ,6,a.o,15,1);c=il(a,b,false);return bp(a,b,c);}function to(a,b){a.d='rgb('+(b.c>>16&255)+','+(b.c>>8&255)+','+(b.c&255)+')';}function QP(a,b){if(a<0||a>b){throw dG(new JH('Index: '+a+', Size: '+b));}}function NP(a,b){if(a<0||a>=b){throw dG(new JH('Index: '+a+', Size: '+b));}}function TP(a,b,c){if(a<0||b>c||b<a){throw dG(new aK(BQ+a+CQ+b+', length: '+c));}}function PB(a,b){var c;if(a.d>a.b+a.i&&BH(b,a.b+a.i)>=53){c=a.b+a.i-1;OB(a,b,c);}}function qe(a,b){var c;c=bl(a.L,b);while(c<Qk(a.L,b)&&dl(a.L,b,c)==0)++c;return c;}function ni(a,b){return(a.s[b]&zR)>>19!=1&&(a.s[b]&zR)>>19!=2?-1:(a.s[b]&AR)>>21;}function Ii(a,b){return(a.C[b]&BR)>>10!=1&&(a.C[b]&BR)>>10!=2?-1:(a.C[b]&CR)>>12;}function WB(a){return a.__elementTypeCategory$==null?10:a.__elementTypeCategory$;}function EO(){EO=FG;AO=new FO('All',0);BO=new JO();CO=new LO();DO=new OO();}function DP(){DP=FG;new JP();new FP('ISO-LATIN-1');new FP('ISO-8859-1');}function Vn(){this.b=8;this.v=new LM();this.H=new SO(new nA());this.c=new SO(new nA());}function Wn(){this.b=1;this.v=new LM();this.H=new SO(new nA());this.c=new SO(new nA());}function YG(a,b,c,d){this.c=(d&255)<<24|(a&255)<<16|(b&255)<<8|c&255;aH(a,b,c,d);}function Fm(a,b){var c,d;d=kK(a.a,b);for(c=d.length;c<10;c++)TJ(a.b,32);WJ(a.b,d);}function el(a,b){var c,d;c=0;for(d=a.g[b];d<a.c[b];d++)a.j[b][d]!=0&&++c;return c;}function hg(a,b){var c,d;c=0;for(d=0;d<a.b;d++)a.e[d][b]&&a.c[d]==-3&&++c;return c;}function vk(a,b,c){var d;for(d=0;d<a.o;d++){a.H[d].a+=b;a.H[d].b+=c;}a.R+=b;a.S+=c;}function Kn(a,b){var c;for(c=0;c<b.length;c++)if(b[c]==a)return true;return false;}function Ck(a){Gh();var b;for(b=1;b<Ch.length;b++)if(EJ(a,Ch[b]))return b;return 0;}function Bm(a){var b,c;c=a.a;for(b=0;b<a.b.length;b++)c+=a.b[b]*xm[a.c[b]];return c;}function Dm(a){var b,c;c=a.d;for(b=0;b<a.b.length;b++)c+=a.b[b]*zm[a.c[b]];return c;}function tq(a){var b,c,d;b=wq(a);d=0;for(c=0;c<qq.length;c++)d+=b[c]*qq[c];return d;}function $l(a){var b,c;Xo(a,3);for(b=0;b<a.d;b++)Yl(a,b);for(c=0;c<a.e;c++)Zl(a,c);}function cB(a){var b,c;if(a.a){c=null;do{b=a.a;a.a=null;c=gB(b,c);}while(a.a);a.a=c;}}function dB(a){var b,c;if(a.b){c=null;do{b=a.b;a.b=null;c=gB(b,c);}while(a.b);a.b=c;}}function jm(a,b){var c;return b==null||b.length==0?null:lm(a,IP((c=b,DP(),c)),null);}function un(a,b,c){var d;d=BM(a.i,b).length;while(c>=d)c-=d;while(c<0)c+=d;return c;}function Bg(a,b,c){var d,e;d=0;for(e=0;e<a.e[c];e++){ah(b,al(a.i,c,e))&&++d;}return d;}function rI(a,b){var c;if(!a){return;}b.j=a;var d=lI(b);if(!d){BG[a]=[b];return;}d.Jb=b;}function ah(a,b){var c;for(c=0;c<a.a.length;c++)if(b==a.a[c])return true;return false;}function xG(){yG();var a=wG;for(var b=0;b<arguments.length;b++){a.push(arguments[b]);}}function Yk(a,b){if(b){Xo(a,1);return Ci(a,a.d,a.e,Fh);}else{return Ci(a,a.o,a.p,Fh);}}function BC(a){if(sC(a,(JC(),IC))<0){return-oC(vC(a));}return a.l+a.m*XT+a.h*YT;}function Qn(a,b){if(!a.F){Tn(a,b);a.F=true;}if(!a.n){Rn(a,b);An(a);zn(a);a.n=true;}}function Ce(a,b,c){if(a.a==null){a.a=ZB(XC,pR,6,a.L.d,15,1);kN(a.a);}a.a[b]=c<<24>>24;}function Pj(a,b,c){a.t==null&&(a.t=ZB(_C,mR,7,a.K,0,2));xN(c);a.t[b]=c;a.Q=0;a.I=true;}function cm(a,b,c){var d;d=wk(a,b,c);if(d&&c==26){Xo(a,3);d=d&(a.C[b]&128)==0;}return d;}function eI(a){var b;b=new cI();b.k='Class$'+(a?'S'+a:''+b.g);b.b=b.k;b.i=b.k;return b;}function cG(a){var b;if(PC(a,14)){return a;}b=a&&a[pQ];if(!b){b=new MA(a);jB(b);}return b;}function Sk(a){var b,c;Xo(a,3);b=0;for(c=0;c<a.n.g.a.length;c++)on(a.n,c)&&++b;return b;}function WI(a){var b,c;if(a==0){return 32;}else{c=0;for(b=1;(b&a)==0;b<<=1){++c;}return c;}}function lH(a,b,c,d,e){var f;if(d<b){f=b;b=d;d=f;}if(e<c){f=c;c=e;e=f;}rH(a,b,c,d-b,e-c);}function $k(a,b,c){var d;for(d=0;d<a.c[b];d++)if(a.f[b][d]==c)return a.i[b][d];return-1;}function Bl(a,b){var c;for(c=0;c<a.g[b];c++)if(a.q[a.f[b][c]]<0)return true;return false;}function Cl(a,b){var c;for(c=0;c<a.g[b];c++)if(a.q[a.f[b][c]]>0)return true;return false;}function Ri(a,b){var c;c=a.A[b]<Dh.length?Dh[a.A[b]]:null;return c==null?6:c[c.length-1];}function di(a){a.o=0;a.p=0;a.I=false;a.J=false;a.G=0;a.t=null;a.r=null;a.M=null;a.Q=0;}function Yh(a,b){b.I=a.I;b.J=a.J;b.P=a.P;b.G=a.G;b.M=a.M==null?null:UP(a.M);b.Q=a.Q&12;}function Ad(a,b){a.w=-5;a.d='rgb('+(b.c>>16&255)+','+(b.c>>8&255)+','+(b.c&255)+')';}function JC(){JC=FG;FC=fC(OR,OR,524287);GC=fC(0,0,tQ);HC=dC(1);dC(2);IC=dC(0);}function SC(a){return a!=null&&(typeof a===fQ||typeof a==='function')&&!(a.Lb===IG);}function AG(a,b){typeof window===fQ&&typeof window['$gwt']===fQ&&(window['$gwt'][a]=b);}function EG(a,b){for(var c in b){b[c]['configurable']=true;}Object.defineProperties(a,b);}function xq(a,b){var c;for(c=0;c<a.g[b];c++)if(ji(a,a.f[b][c])<0)return true;return false;}function Hd(a,b){var c;for(c=0;c<bl(a.d,b);c++)if(a.c[cl(a.d,b,c)])return true;return false;}function an(a){var b,c;c=0;for(b=0;b<a.a.o;b++)(Ai(a.a,b)==7||Ai(a.a,b)==8)&&++c;return c;}function gO(a,b){var c,d;c=1-b;d=a.a[c];a.a[c]=d.a[b];d.a[b]=a;a.b=true;d.b=false;return d;}function im(a,b){var c,d,e,f;d=b/2|0;e=a>=d;e&&(a-=d);f=b/32|0;c=f*a/(d-a);return e?-c:c;}function rG(a){var b,c,d,e;e=a;d=0;if(e<0){e+=YT;d=VT;}c=WC(e/XT);b=WC(e-c*XT);return fC(b,c,d);}function lB(a){var b=/function(?:\s+([\w$]+))?\s*\(/;var c=b.exec(a);return c&&c[1]||gQ;}function _h(a,b){var c;if(b.length==0)return null;for(c=0;c<b.length;c++)uj(a,b[c]);return ci(a);}function hG(a){var b;b=a.h;if(b==0){return a.l+a.m*XT;}if(b==VT){return a.l+a.m*XT-YT;}return a;}function Qf(a){var b,c;c=true;for(b=0;b<a.i.d;b++){if(a.o[b]!=0&&!a.e[b]){c=false;break;}}return c;}function ep(a){var b,c;Xo(a,15);c=0;for(b=0;b<a.d;b++)(a.s[b]&3)!=0&&(a.s[b]&4)==0&&++c;return c;}function AN(a){var b,c,d;d=0;for(c=a.yb();c.Bb();){b=c.Cb();d=d+(b!=null?tc(b):0);d=d|0;}return d;}function pn(a,b,c){var d,e;e=BM(a.g,b);for(d=0;d<e.length;d++)if(c==e[d])return true;return false;}function qn(a,b,c){var d,e;e=BM(a.i,b);for(d=0;d<e.length;d++)if(c==e[d])return true;return false;}function qm(a,b,c,d){var e,f;this.a=rm(a,b,c,d);e=c-a;f=d-b;this.b=$wnd.Math.sqrt(e*e+f*f);}function fm(a,b){return $wnd.Math.pow(10,$wnd.Math.log(2000)*$wnd.Math.LOG10E*a/(b-1)-1);}function jG(a){if(ZT<a&&a<YT){return a<0?$wnd.Math.ceil(a):$wnd.Math.floor(a);}return hG(tC(a));}function qd(a,b){var c;if(b>0)return(a[b]+a[b-1])/2;c=MQ+(a[0]+a[a.length-1])/2;return c>MQ?c-LQ:c;}function qC(a,b){var c,d,e;c=a.l+b.l;d=a.m+b.m+(c>>22);e=a.h+b.h+(d>>22);return fC(c&OR,d&OR,e&VT);}function AC(a,b){var c,d,e;c=a.l-b.l;d=a.m-b.m+(c>>22);e=a.h-b.h+(d>>22);return fC(c&OR,d&OR,e&VT);}function km(a,b,c){var d,e;return b==null?null:lm(a,IP((e=b,DP(),e)),c==null?null:IP((d=c,d)));}function rO(a,b){var c;this.c=a;c=new LM();YN(a,c,b,a.b,null,false,null,false);this.a=new YK(c,0);}function ig(a,b){var c;for(c=0;c<a.b;c++)if(a.e[c][b]&&a.c[c]==-3)return c<a.a?1:c<a.b?2:0;return-1;}function kg(a,b,c){var d;for(d=0;d<a.j.g.length;d++)if(a.e[b][d]&&a.e[c][d])return true;return false;}function rk(a,b){var c,d;d=0;for(c=0;c<a.p;c++)(a.B[0][c]==b||a.B[1][c]==b)&&(d+=Mi(a,c));return d;}function BN(a){var b,c,d;d=1;for(c=a.yb();c.Bb();){b=c.Cb();d=31*d+(b!=null?tc(b):0);d=d|0;}return d;}function Nd(a,b){var c,d;--a.a;for(d=0;d<bl(a.d,b);d++){c=cl(a.d,b,d);if(a.c[c]){a.c[c]=false;--a.b;}}}function Md(a,b){var c;for(c=0;c<a.d.e;c++){if(a.c[c]&&sn(b,c)){Nd(a,Ei(a.d,0,c));Nd(a,Ei(a.d,1,c));}}}function fd(a){var b,c;for(c=new dN(a.N);c.a<c.c.a.length;){b=cN(c);zd(a,b.a);Oc(a,b.b,b.c);}zd(a,a.J);}function EH(a){var b;b=a.a.length;0<b?a.a=a.a.substr(0,0):0>b&&(a.a+=PJ(ZB(YC,gR,6,-b,15,1)));}function JB(a,b){var c,d,e;e=a.a.length;for(d=0;d<e;++d){c=yJ(a.a,d);c>=48&&c<=57&&DH(a,d,c-48+b&AQ);}}function oj(a,b){var c;c=a.A[b];return c==1||c>=5&&c<=9||c>=14&&c<=17||c>=33&&c<=35||c>=52&&c<=53;}function vC(a){var b,c,d;b=~a.l+1&OR;c=~a.m+(b==0?1:0)&OR;d=~a.h+(b==0&&c==0?1:0)&VT;return fC(b,c,d);}function wq(a){var b,c;c=ZB(_C,DQ,6,qq.length+2,15,1);Xo(a,3);for(b=0;b<a.d;b++)++c[vq(a,b)];return c;}function bn(a){var b,c;c=0;for(b=0;b<a.a.o;b++)(Ai(a.a,b)==7||Ai(a.a,b)==8)&&Rk(a.a,b)>0&&++c;return c;}function aO(a,b,c){var d,e;d=new xO(b,c);e=new yO();a.b=$N(a,a.b,d,e);e.b||++a.c;a.b.b=false;return e.d;}function Yz(a,b){if(b){!a.c&&(a.c=new om(true));return a.c;}else{!a.b&&(a.b=new om(false));return a.b;}}function QH(a){if(DJ(typeof a,nQ)){return true;}return a!=null&&a.$implements__java_lang_CharSequence;}function aC(a,b,c,d,e){e.Jb=a;e.Kb=b;e.Lb=IG;e.__elementTypeId$=c;e.__elementTypeCategory$=d;return e;}function kd(a,b,c,d,e){yM(a.T,new tH(b-a.O,c-a.O,2*a.O,2*a.O));e&&yM(a.N,new Ed(b,c,$c(a,d)?-3:a.o[d]));}function jh(a,b){return $wnd.Math.sqrt((b.a-a.a)*(b.a-a.a)+(b.b-a.b)*(b.b-a.b)+(b.c-a.c)*(b.c-a.c));}function ok(a,b){a.B[0]=hN(a.B[0],b);a.B[1]=hN(a.B[1],b);a.F=hN(a.F,b);a.C=hN(a.C,b);a.D=hN(a.D,b);a.L=b;}function lC(a){var b,c,d;b=~a.l+1&OR;c=~a.m+(b==0?1:0)&OR;d=~a.h+(b==0&&c==0?1:0)&VT;a.l=b;a.m=c;a.h=d;}function mC(a){var b,c;c=VI(a.h);if(c==32){b=VI(a.m);return b==32?VI(a.l)+32:b+20-10;}else{return c-12;}}function AM(a,b){var c,d;c=b.toArray();d=c.length;if(d==0){return false;}AP(a.a,a.a.length,c);return true;}function gg(a,b){var c;for(c=0;c<a.b;c++)if(a.e[c][b]&&a.c[c]==-3)return c<a.a?c:c<a.b?c-a.a:-1;return-1;}function zI(a,b){if(a<b){return-1;}if(a>b){return 1;}if(a==b){return 0;}return isNaN(a)?isNaN(b)?0:1:-1;}function DA(a){var b;if(a!=null){b=a[pQ];if(b){return b;}}return TC(a,$wnd.TypeError)?new qJ(a):new IA(a);}function wA(a){var b,c,d,e;for(b=(a.g==null&&(a.g=(iB(),e=hB.qb(a),kB(e))),a.g),c=0,d=b.length;c<d;++c);}function rN(a){var b,c,d,e;e=1;for(c=0,d=a.length;c<d;++c){b=a[c];e=31*e+(b!=null?tc(b):0);e=e|0;}return e;}function Nf(a,b,c,d){var e,f;for(f=0;f<bl(a.i,c);f++){e=al(a.i,c,f);if(!d[e]&&Rf(a,b,e))return e;}return-1;}function pK(a,b){var c,d;OP(b);for(d=b.yb();d.Bb();){c=d.Cb();if(!a.contains(c)){return false;}}return true;}function XN(a,b){var c,d,e;e=a.b;while(e){c=a.a.eb(b,e.c);if(c==0){return e;}d=c<0?0:1;e=e.a[d];}return null;}function iC(a,b,c,d,e){var f;f=yC(a,b);c&&lC(f);if(e){a=kC(a,b);d?cC=vC(a):cC=fC(a.l,a.m,a.h);}return f;}function sd(a,b,c,d){var e;if((a.B&1)!=0)return false;e=BM(a.T,d);return b>e.c&&b<e.c+e.b&&c>e.d&&c<e.d+e.a;}function IB(a,b,c,d){var e;if(d>0){for(e=d;e<a.b;e+=d+1){YJ(b,a.b-e,String.fromCharCode(c));++a.b;++a.d;}}}function Ne(a,b,c){while(c!=0){if(a.p==0){TJ(a.q,a.r+64&AQ);a.p=6;a.r=0;}a.r<<=1;a.r|=b&1;b>>=1;--c;--a.p;}}function fe(a,b){if(Ah(a)==-1||Ah(b)==-1)return 3;if(((Ah(a)|Ah(b))&1)!=0)return 3;return Ah(a)==Ah(b)?1:2;}function Pm(a){if(a.indexOf('ATOMS=(')!=-1)return LR;if(a.indexOf('BONDS=(')!=-1)return'BONDS';return null;}function op(a,b){var c,d;c=So(a,b);if(c==-1)return-1;d=a.b.a.length;yM(a.b,b);xM(a.a,c,new QI(d));return d;}function Lp(a){Ip();var b,c;Xo(a,1);c=0;for(b=0;b<a.o;b++){a.A[b]==1?++c:c+=a.c[b]-a.g[b]+ml(a,b);}return c;}function lo(a,b){$n();var c,d,e;e=0;c=0;for(d=0;d<a.length;d++){e+=ho(a[d]&b[d]);c+=ho(a[d]|b[d]);}return e/c;}function wm(a,b){vm();var c,d;d=b-a;for(c=0;c<um[a].length;c++)if(um[a][c].b==d)return um[a][c].a;return NaN;}function Wf(a,b,c){var d,e;for(e=0;e<a.g[b].length;e++){d=a.g[b][e];if(a.k[d]==2){a.k[d]=1;a.j[d]=c<<24>>24;}}}function ei(a){var b,c;c=false;for(b=0;b<a.o;b++){if((a.s[b]&512)!=0){a.A[b]=-1;c=true;}}return c&&ci(a)!=null;}function tk(a){var b,c;a.J=false;for(b=0;b<a.o;b++)a.s[b]&=-133693441;for(c=0;c<a.p;c++)a.F[c]&=-25;a.Q&=-253;}function Ik(a,b){var c,d;for(d=0;d<a.c[b];d++){c=a.i[b][d];(a.F[c]==17||a.F[c]==9)&&a.B[0][c]==b&&(a.F[c]=1);}}function gG(a,b){var c;if(kG(a)&&kG(b)){c=a-b;if(!isNaN(c)){return c;}}return sC(kG(a)?rG(a):a,kG(b)?rG(b):b);}function PH(a,b){MH();return UC(a)?zJ(a,b):RC(a)?zI((OP(a),a),(OP(b),b)):QC(a)?OH((OP(a),a),(OP(b),b)):a.fb(b);}function rc(a){return UC(a)?kF:RC(a)?SE:QC(a)?QE:OC(a)?a.Jb:_B(a)?a.Jb:a.Jb||Array.isArray(a)&&VB(nE,1)||nE;}function NG(a){MG.call(this,aU+a+bU+cU==null?rQ:vc(aU+a+bU+cU),PC(aU+a+bU+cU,14)?aU+a+bU+cU:null);}function IA(a){tA(this);uA(this);this.backingJsObject=a;a!=null&&WP(a,pQ,this);this.f=a==null?rQ:vc(a);}function cI(){this.g=_H++;this.k=null;this.i=null;this.f=null;this.d=null;this.b=null;this.j=null;this.a=null;}function RN(){NN();var a,b,c;c=MN++ +CP();a=WC($wnd.Math.floor(c*kU))&EQ;b=WC(c-a*GR);this.a=a^1502;this.b=b^jU;}function aq(){Zp();var a,b,c,d;if(!Yp){if(!Yp){Yp=new sA();for(b=Wp,c=0,d=b.length;c<d;++c){a=hJ(b[c]);qA(Yp,a);}}}}function CB(a,b){var c,d;b.a+='E';if(a.e<0){a.e=-a.e;b.a+='-';}c=''+a.e;for(d=c.length;d<a.k;++d){b.a+='0';}b.a+=c;}function ad(a,b){var c;if(bl(a.G,b)!=2)return false;for(c=0;c<2;c++)if(dl(a.G,b,c)!=2)return false;return true;}function eA(a){var b,c,d;d=a.a.a.length;b=new Array(d);for(c=0;c<d;c++){b[c]=fA(BM(a.a,c).a,BM(a.a,c).b);}return b;}function Om(a,b){var c;for(c=b;c<a.length;c++){if(a.charCodeAt(c)==32||a.charCodeAt(c)==9){return c;}}return-1;}function YI(a){var b,c;if(a>-129&&a<128){b=a+128;c=($I(),ZI)[b];!c&&(c=ZI[b]=new QI(a));return c;}return new QI(a);}function qK(a,b){var c,d,e;OP(b);c=false;for(d=a.yb();d.Bb();){e=d.Cb();if(b.contains(e)){d.Db();c=true;}}return c;}function xp(a,b){if(a.A[b]!=6)return false;if(a.q[b]!=0)return false;if(ml(a,b)+a.g[b]!=4)return false;return true;}function TH(a){if(a>=48&&a<58){return a-48;}if(a>=97&&a<97){return a-97+10;}if(a>=65&&a<65){return a-65+10;}return-1;}function Mi(a,b){switch(a.F[b]&103){case 1:case 64:return 1;case 2:return 2;case 4:return 3;default:return 0;}}function eG(a,b){var c;if(kG(a)&&kG(b)){c=a+b;if(ZT<c&&c<YT){return c;}}return hG(qC(kG(a)?rG(a):a,kG(b)?rG(b):b));}function mG(a,b){var c;if(kG(a)&&kG(b)){c=a*b;if(ZT<c&&c<YT){return c;}}return hG(uC(kG(a)?rG(a):a,kG(b)?rG(b):b));}function zM(a,b,c){var d,e;QP(b,a.a.length);d=c.toArray();e=d.length;if(e==0){return false;}AP(a.a,b,d);return true;}function Sh(a,b,c){if(c){if(a.q[b]>3)return false;++a.q[b];}else{if(a.q[b]<-3)return false;--a.q[b];}a.Q=0;return true;}function up(a,b,c){var d,e;d=XB(kF,[iQ,vR],[28,2],6,[a.d,b],2);Xo(a,3);for(e=0;e<a.d;e++){d[e]=vp(a,e,b,c);}return d;}function Qo(a,b,c){var d,e,f;e=false;for(d=1;d<c;d++){for(f=0;f<d;f++){a[f]>a[d]&&(e=!e);b[f]>b[d]&&(e=!e);}}return e;}function Vi(a,b){var c,d,e,f;f=3;for(d=0;d<2;d++){c=a.B[d][b];e=Mi(a,b)+(Ui(a,c)+Si(a,c))-ol(a,c);f>e&&(f=e);}return f;}function sq(a,b){var c,d;c=wq(a);for(d=0;d<qq.length;d++)c[d]!=0&&nq(b,''+c[d]+' * '+qq[d]+'   AtomType: '+pq[d],2);}function Cm(a){var b,c;b=new SJ();for(c=0;c<a.b.length;c++){RJ(b,(Gh(),Ch)[a.c[c]]);a.b[c]>1&&RJ(b,''+a.b[c]);}return b.a;}function _c(a){var b;a.p=ZB(aG,HQ,6,a.G.o,16,1);for(b=0;b<a.G.p;b++){a.p[Ei(a.G,0,b)]=true;a.p[Ei(a.G,1,b)]=true;}}function le(a){var b,c;b=0;yN(a.b);for(c=0;c<a.b.length;c++){(c==0||Ef(a.b[c],a.b[c-1])!=0)&&++b;a.c[a.b[c].a]=b;}return b;}function tN(a,b,c,d,e,f,g,h){var i;i=c;while(f<g){i>=d||b<c&&h.eb(a[b],a[i])<=0?e[f++]=a[b++]:e[f++]=a[i++];}}function sN(a,b,c,d){var e,f,g;for(e=b+1;e<c;++e){for(f=e;f>b&&d.eb(a[f-1],a[f])>0;--f){g=a[f];a[f]=a[f-1];a[f-1]=g;}}}function qo(a,b,c,d){var e;e='<circle cx="'+WC(b)+_Q+'cy="'+WC(c)+_Q+'r="'+WC(d)+_Q+'fill="'+a.d+'" />';xo(a,e);}function Oo(a,b,c,d){if(a.b)return;if(a.g==4||a.g==3&&a.c!=-1){a.b=true;return;}a.i[a.g]=d;a.f[a.g]=b;a.j[a.g]=c;++a.g;}function ih(a,b){if(a.a!=b.a)return a.a<b.a?-1:1;if(a.b!=b.b)return a.b<b.b?-1:1;if(a.c!=b.c)return a.c<b.c?-1:1;return 0;}function $c(a,b){var c;if(Qk(a.G,b)==0)return false;for(c=0;c<Qk(a.G,b);c++)if(!fj(a.G,cl(a.G,b,c)))return false;return true;}function $f(a,b){var c,d;d=ZB(_C,DQ,6,a==null?1:a.length+1,15,1);for(c=0;c<d.length-1;c++)d[c]=a[c];d[d.length-1]=b;return d;}function wn(a){var b;b=hN(a.w,a.w.length);if(a.j!=0){xN(b);if(!RO(a.H,b)){QO(a.H,b);yM(a.v,hN(a.w,a.w.length));}return;}return;}function lf(a){var b,c;if(a.R!=null)for(b=0;b<a.L.d;b++)Ij(a.L,b,a.R[b]);if(a.f!=null)for(c=0;c<a.L.e;c++)dk(a.L,c,a.f[c]);}function hq(){if(!eq){try{dq=new kq(cq);eq=true;}catch(a){a=cG(a);if(PC(a,12)){cK();}else throw dG(a);}}}function RP(a,b,c){if(a<0||b>c){throw dG(new JH(BQ+a+CQ+b+', size: '+c));}if(a>b){throw dG(new NI(BQ+a+' > toIndex: '+b));}}function ak(a,b,c){if(c>=0&&c<=190){if(c==151||c==152){a.A[b]=1;a.v[b]=c-149;}else{a.A[b]=c;a.v[b]=0;}a.s[b]&=268435455;a.Q=0;}}function mj(a,b){var c;c=a.A[b];return c>=3&&c<=4||c>=11&&c<=13||c>=19&&c<=31||c>=37&&c<=51||c>=55&&c<=84||c>=87&&c<=103;}function wg(a,b){var c,d;d=0;uN(a,0,a.length,null);for(c=0;c<a.length;c++){(c==0||Ef(a[c],a[c-1])!=0)&&++d;b[a[c].a]=d;}return d;}function QJ(a,b,c){var d,e,f,g;f=b+c;TP(b,f,a.length);g='';for(e=b;e<f;){d=e+10000<f?e+10000:f;g+=NJ(a.slice(e,d));e=d;}return g;}function rm(a,b,c,d){var e,f,g;f=c-a;g=d-b;if(g!=0){e=$wnd.Math.atan(f/g);g<0&&(f<0?e-=MQ:e+=MQ);}else e=f>0?NQ:ZQ;return e;}function Nn(a,b){var c,d,e,f;e=0;f=0;while(e<a.length&&f<b.length){c=a[e];d=b[f];if(c==d)return true;c<d?++e:++f;}return false;}function Ki(a,b){var c,d,e,f;c=a.B[0][b];d=a.B[1][b];e=a.H[d].a-a.H[c].a;f=a.H[d].b-a.H[c].b;return $wnd.Math.sqrt(e*e+f*f);}function Ui(a,b){var c,d;c=((a.s[b]&yR)>>>28)-1;c==-1&&(c=(d=a.A[b]<Dh.length?Dh[a.A[b]]:null,d==null?6:d[d.length-1]));return c;}function cQ(a){aQ();var b,c,d;c=':'+a;d=_P[c];if(!(d===undefined)){return d;}d=ZP[c];b=d===undefined?bQ(a):d;dQ();_P[c]=b;return b;}function CG(){BG={};!Array.isArray&&(Array.isArray=function(a){return Object.prototype.toString.call(a)==='[object Array]';});}function hC(a,b){if(a.h==tQ&&a.m==0&&a.l==0){b&&(cC=fC(0,0,0));return eC((JC(),HC));}b&&(cC=fC(a.l,a.m,a.h));return fC(0,0,0);}function hJ(a){var b,c;if(gG(a,-129)>0&&gG(a,128)<0){b=tG(a)+128;c=(jJ(),iJ)[b];!c&&(c=iJ[b]=new aJ(a));return c;}return new aJ(a);}function oK(a,b,c){var d,e;for(e=a.yb();e.Bb();){d=e.Cb();if(VC(b)===VC(d)||b!=null&&pc(b,d)){c&&e.Db();return true;}}return false;}function hf(a,b,c){var d,e;a.N=b;for(d=0;d<a.L.d;d++){a.c[d]=c[d];a.W[d]=0;a.$[d]=false;}for(e=0;e<a.L.e;e++){a.k[e]=0;a.o[e]=false;}}function lk(a,b){var c,d;a.I=b;if(!b){a.t=null;for(c=0;c<a.o;c++)a.w[c]=0;for(d=0;d<a.p;d++){a.D[d]=0;a.F[d]==64&&(a.F[d]=1);}}a.Q=0;}function sk(a){var b,c,d;c=false;d=false;for(b=0;b<a.o;b++){if(a.v[b]!=0){a.v[b]=0;c=true;a.A[b]==1&&(d=true);}}d&&(a.Q=0);return c;}function gm(a,b){var c,d;c=b;d=0;while(b!=0){if(a.c==0){a.e=a.a[++a.d]-64<<11;a.c=6;}d|=(zQ&a.e)>>16-c+b;a.e<<=1;--b;--a.c;}return d;}function Sc(a,b){var c,d;for(d=0;d<a.T.a.length;d++)a.t=qH(a.t,BM(a.T,d));Tc(a,b);c=0.1*b;a.t.c-=c;a.t.d-=c;a.t.b+=2*c;a.t.a+=2*c;}function nH(a,b,c,d,e){var f,g;if(a.b<=0||a.a<=0||d<=0||e<=0){return false;}f=a.c;g=a.d;return b>=f&&c>=g&&b+d<=f+a.b&&c+e<=g+a.a;}function ql(a,b,c,d){var e,f;Xo(a,1);for(e=0;e<d;e++){for(f=0;f<a.g[b[e]];f++){if(a.f[b[e]][f]===b[e+1]){c[e]=a.i[b[e]][f];break;}}}}function Nm(a,b){var c;if(b==-1){return-1;}for(c=b+1;c<a.length;c++){if(a.charCodeAt(c)!=32&&a.charCodeAt(c)!=9){return c;}}return-1;}function SH(a,b,c){var d,e;d=yJ(a,b++);if(d>=55296&&d<=56319&&b<c&&XH(e=a.charCodeAt(b))){return zQ+((d&1023)<<10)+(e&1023);}return d;}function LI(a){var b;b=wI(a);if(b>3.4028234663852886E38){return Infinity;}else if(b<-3.4028234663852886E38){return-Infinity;}return b;}function tI(a){var b;b=typeof a;if(DJ(b,lQ)||DJ(b,mQ)||DJ(b,nQ)){return true;}return a!=null&&a.$implements__java_lang_Comparable;}function Np(){Np=FG;Mp=aC(VB(kF,1),vR,2,6,['Actelion No','ID','IDNUMBER','COMPOUND_ID','NAME','COMPND']);}function iA(){iA=FG;hA=aC(VB($C,1),gR,6,15,[0.29899999499320984,0.5870000123977661,0.11400000005960464]);}function Re(a,b){var c,d,e,f;c=b/2|0;e=a<0;a=$wnd.Math.abs(a);f=b/32|0;d=mJ(c-1,tG(jG($wnd.Math.round(a*c/(a+f)))));return e?c+d:d;}function Ak(a,b,c,d){Gh();var e,f,g;f=c-a;g=d-b;if(g!=0){e=$wnd.Math.atan(f/g);g<0&&(f<0?e-=MQ:e+=MQ);}else e=f>0?NQ:ZQ;return e;}function Li(a,b,c){var d;for(d=0;d<a.p;d++)if(a.B[0][d]==b&&a.B[1][d]==c||a.B[0][d]==c&&a.B[1][d]==b)if(a.F[d]!=128)return d;return-1;}function Qg(a){var b,c;for(b=0;b<a.b;b++){if(a.e[b]==0){c=new hh(a,a.i,1);a.a[b]=true;c.a[0]=b;c.c[0]=0;c.d[0]=0;c.q[0]=0;yM(a.f,c);}}}function Kf(a,b){var c,d;for(d=0;d<a.g[b].length;d++){c=a.g[b][d];if(a.f[c]&&(a.o[c]==1||a.o[c]==2)&&a.k[c]==0)return true;}return false;}function VA(){var a;if(QA!=0){a=PA();if(a-RA>2000){RA=a;SA=$wnd.setTimeout(_A,10);}}if(QA++==0){cB((bB(),aB));return true;}return false;}function Dk(a){Gh();switch(a){case 7:case 8:case 9:case 15:case 16:case 17:case 33:case 34:case 35:case 53:return true;}return false;}function Co(a){switch(a){case 5:case 6:case 7:case 8:case 9:case 15:case 16:case 17:case 36:case 53:return true;default:return false;}}function Zf(a,b,c,d,e,f,g,h,i,j){this.i=a;this.a=b;this.f=c;this.o=d;this.c=e;this.k=f;this.j=g;this.p=h;this.d=i;this.n=j;Mf(this);}function cn(b){var c;try{return $p((new aq(),b.a));}catch(a){a=cG(a);if(PC(a,12)){c=a;vA(c,(cK(),bK),'');return-999;}else throw dG(a);}}function DB(a,b,c){if(a.d==0){b.a=b.a.substr(0,0)+'0'+KJ(b.a,0);++a.b;++a.d;}if(a.b<a.d||a.c){YJ(b,a.b,String.fromCharCode(c));++a.d;}}function nB(){if(Error.stackTraceLimit>0){$wnd.Error.stackTraceLimit=Error.stackTraceLimit=64;return true;}return'stack'in new Error();}function qI(a,b){var c=0;while(!b[c]||b[c]==''){c++;}var d=b[c++];for(;c<b.length;c++){if(!b[c]||b[c]==''){continue;}d+=a+b[c];}return d;}function ko(a,b){$n();var c,d,e,f;f=0;d=0;e=0;for(c=0;c<a.length;c++){f+=ho(a[c]&b[c]);d+=ho(a[c]);e+=ho(b[c]);}return f/$wnd.Math.sqrt(d*e);}function EJ(a,b){OP(a);if(b==null){return false;}if(DJ(a,b)){return true;}return a.length==b.length&&DJ(a.toLowerCase(),b.toLowerCase());}function ao(){var a,b;if(Zn==null){b=new om(false);Zn=ZB(PD,iQ,24,Yn.length,0,1);for(a=0;a<Yn.length;a++){Zn[a]=jm(b,Yn[a]);Xo(Zn[a],1);}}}function Gz(a,b){b=b||{};var c=(typeof b.maxSphereSize===TT?5:b.maxSphereSize)|0;var d=(typeof b.type===TT?0:b.type)|0;return wp(a,c,d);}function Lc(a){var b;b=a.K.c*Bi(a.G);a.R=b*0.06;a.M=b*0.15;a.L=b*0.38;a.P=b*0.47;a.Q=WC(b*a.F*0.6+0.5);a.O=b*0.12;a.S=b*0.4;a.v=b*0.5+0.5;}function Pd(a){var b;this.d=a;Xo(a,1);this.c=ZB(aG,HQ,6,a.e,16,1);for(b=0;b<a.e;b++){if(a.F[b]==64){this.c[b]=true;a.F[b]=1;a.Q=0;}}Gd(this);}function vl(a,b){var c,d,e;e=ZB(_C,DQ,6,a.c[b],15,1);for(d=0;d<a.c[b];d++)e[d]=(a.f[b][d]<<16)+d;xN(e);for(c=0;c<a.c[b];c++)e[c]&=AQ;return e;}function Rd(a,b){var c,d;c=0;for(d=0;d<a.g[b];d++)a.j[b][d]==2&&(Ai(a,a.f[b][d])==7||Ai(a,a.f[b][d])==8||Ai(a,a.f[b][d])==16)&&++c;return c;}function Ln(a,b){var c,d,e;e=0;for(d=0;d<a.length;d++){c=a[d];while(b[e]<c){++e;if(e==b.length)return false;}if(b[e]>c)return false;}return true;}function Lf(a,b,c){var d,e,f,g,h;e=0;g=0;for(h=0;h<a.g[b].length;h++){d=a.g[b][h];if(a.k[d]==c){f=1<<a.j[d];if((g&f)==0){g|=f;++e;}}}return e;}function ON(a,b){var c,d;KP(b>0);if((b&-b)==b){return WC(b*PN(a)*4.6566128730773926E-10);}do{c=PN(a);d=c%b;}while(c-d+(b-1)<0);return WC(d);}function Nj(a,b,c){c!=null&&c.length==0&&(c=null);if(c==null){a.r!=null&&(a.r[b]=null);}else{a.r==null&&(a.r=ZB(XC,xR,13,a.K,0,2));a.r[b]=c;}}function fv(a,b){var e,f;ou();b=b||{};var c=!b.noCoordinates;var d=!b.noStereo;return e=new qu(),Io(aA(nu),e.a,IP((f=a,DP(),f)),d),c&&pu(e),e;}function of(a){var b,c;for(b=0;b<a.L.d;b++)!$i(a.L,b)&&a.W[b]==3&&Lj(a.L,b,true);for(c=0;c<a.L.e;c++){a.k[c]==3&&Mi(a.L,c)==2&&jk(a.L,c,26);}}function Be(a){var b,c;for(b=0;b<a.L.d;b++)(!a.H[b]||a.W[b]==3)&&(a.U[b]=0);for(c=0;c<a.L.e;c++)(Pi(a.L,c)!=1||a.k[c]==0||a.k[c]==3)&&(a.j[c]=0);}function co(a){var b,c;if(a.I){for(b=0;b<a.o;b++){if((a.w[b]&IQ)!=0){a=new mp(a);for(c=b;c<a.o;c++)(a.w[c]&IQ)!=0&&(a.A[c]=-1);ci(a);}}}return a;}function _f(a,b){var c;if(a.length!=b.length)return a.length<b.length?-1:1;for(c=0;c<a.length;c++)if(a[c]!==b[c])return a[c]<b[c]?-1:1;return 0;}function KM(a,b){var c,d,e;e=a.a.length;b.length<e&&(b=(d=new Array(e),bC(d,b)));for(c=0;c<e;++c){b[c]=a.a[c];}b.length>e&&(b[e]=null);return b;}function sg(a,b,c){var d,e,f;f=b.length;d=new hh(a,a.i,f);d.c[0]=0;d.d[0]=0;for(e=0;e<f;e++){d.q[e]=128-f;d.a[e]=b[e];}f<8?zg(d):yg(a,d,c);yM(a.f,d);}function If(a,b,c,d){var e,f;this.a=ZB(_C,DQ,6,b,15,1);this.b=ZB(_C,DQ,6,d,15,1);for(e=0;e<b;e++)this.a[e]=a[e];for(f=0;f<d;f++)this.b[f]=c[f];}function hh(a,b,c){this.r=a;this.p=b;this.a=ZB(_C,DQ,6,c,15,1);this.q=ZB(_C,DQ,6,c,15,1);this.c=ZB(ZC,GQ,6,c,15,1);this.d=ZB(ZC,GQ,6,c,15,1);}function tc(a){return UC(a)?cQ(a):RC(a)?WC((OP(a),a)):QC(a)?(OP(a),a)?1231:1237:OC(a)?a.cb():_B(a)?YP(a):!!a&&!!a.hashCode?a.hashCode():YP(a);}function pc(a,b){return UC(a)?DJ(a,b):RC(a)?(OP(a),a===b):QC(a)?(OP(a),a===b):OC(a)?a.ab(b):_B(a)?a===b:!!a&&!!a.equals?a.equals(b):VC(a)===VC(b);}function vc(a){return UC(a)?(OP(a),a):RC(a)?''+(OP(a),a):QC(a)?''+(OP(a),a):OC(a)?a.db():_B(a)?mc(a):a.toString?a.toString():'[JavaScriptObject]';}function YH(a,b,c){RH(a>=0&&a<=1114111);if(a>=zQ){b[c++]=55296+(a-zQ>>10&1023)&AQ;b[c]=56320+(a-zQ&1023)&AQ;return 2;}else{b[c]=a&AQ;return 1;}}function Um(a,b){var c,d,e;if(!a.c){if(DJ(b.substr(0,6),'COUNTS')){c=Nm(b,Om(b,7));d=xI(LJ(b,7,Om(b,7)));e=xI(LJ(b,c,Om(b,c)));a.c=new lp(d,e);}}}function NC(a,b){if(UC(a)){return!!MC[b];}else if(a.Kb){return!!a.Kb[b];}else if(RC(a)){return!!LC[b];}else if(QC(a)){return!!KC[b];}return false;}function Jo(a,b){var c,d,e;Pi(a.b,b)==1&&jk(a.b,b,2);for(d=0;d<2;d++){c=Ei(a.b,d,b);Sj(a.b,c,false);for(e=0;e<bl(a.b,c);e++)a.a[cl(a.b,c,e)]=false;}}function YN(a,b,c,d,e,f,g,h){var i,j;if(!d){return;}i=d.a[0];!!i&&YN(a,b,c,i,e,f,g,h);ZN(a,c,d.c,e,f,g,h)&&b.add(d);j=d.a[1];!!j&&YN(a,b,c,j,e,f,g,h);}function Ef(a,b){var c;for(c=0;c<a.c;c++)if(nG(a.d[c],b.d[c]))return lG(a.d[c],b.d[c])?-1:1;return iG(a.d[a.c],b.d[a.c])?0:lG(a.d[a.c],b.d[a.c])?-1:1;}function Nk(a,b){var c;if(a.g[b]==3&&(a.s[b]&xQ)!=0&&(!!a.n&&b<a.d?jn(a.n,b):0)>=6)for(c=0;c<a.g[b];c++)if(Gl(a,a.i[b][c]))return a.i[b][c];return-1;}function dK(a,b,c){cK();var d,e;PP(a,'src');PP(b,'dest');rc(a);rc(b);e=a.length;d=b.length;if(c<0||c>e||c>d){throw dG(new IH());}c>0&&yP(a,0,b,0,c,true);}function eL(a,b){var c,d,e;for(d=new qO(new vO(a).b);RK(d.a);){c=d.b=SK(d.a);e=c.Gb();if(VC(b)===VC(e)||b!=null&&pc(b,e)){return true;}}return false;}function QB(a,b,c){if(!b){throw dG(new NI('Unknown currency code'));}this.s=a;this.a=b;LB(this,this.s);if(!c&&this.g){this.n=this.a[2]&7;this.i=this.n;}}function MJ(a){var b,c,d;c=a.length;d=0;while(d<c&&a.charCodeAt(d)<=32){++d;}b=c;while(b>d&&a.charCodeAt(b-1)<=32){--b;}return d>0||b<c?a.substr(d,b-d):a;}function Vf(a,b,c,d){var e,f;for(f=0;f<a.g[b].length;f++){e=a.g[b][f];if(a.f[e]&&(a.o[e]==1||a.o[e]==2)&&a.k[e]==0){a.k[e]=d<<24>>24;a.j[e]=c<<24>>24;}}}function Ek(a){Gh();if(a==1||a==6)return false;if(Dk(a))return false;if(a==2||a==10||a==18||a==36||a==54)return false;if(a>103)return false;return true;}function lm(a,b,c){var d,e,f,g,h;if(b==null)return null;hm(a,b,0);d=gm(a,4);g=gm(a,4);d>8&&(d=g);e=gm(a,d);f=gm(a,g);h=new lp(e,f);nm(a,h,b,c);return h;}function kC(a,b){var c,d,e;if(b<=22){c=a.l&(1<<b)-1;d=e=0;}else if(b<=44){c=a.l;d=a.m&(1<<b-22)-1;e=0;}else{c=a.l;d=a.m;e=a.h&(1<<b-44)-1;}return fC(c,d,e);}function vA(a,b,c){var d,e,f,g,h;wA(a);for(e=(a.i==null&&(a.i=ZB(lF,iQ,14,0,0,1)),a.i),f=0,g=e.length;f<g;++f){d=e[f];vA(d,b,'\t'+c);}h=a.e;!!h&&vA(h,b,c);}function $B(a,b){var c=new Array(b);var d;switch(a){case 14:case 15:d=0;break;case 16:d=false;break;default:return c;}for(var e=0;e<b;++e){c[e]=d;}return c;}function od(a){var b;b=new sH();if(a.a<=a.b){b.c=a.a;b.b=a.b-a.a;}else{b.c=a.b;b.b=a.a-a.b;}if(a.c<=a.d){b.d=a.c;b.a=a.d-a.c;}else{b.d=a.d;b.a=a.c-a.d;}return b;}function pC(a,b){var c,d,e;e=a.h-b.h;if(e<0){return false;}c=a.l-b.l;d=a.m-b.m+(c>>22);e+=d>>22;if(e<0){return false;}a.l=c&OR;a.m=d&OR;a.h=e&VT;return true;}function NN(){NN=FG;var a,b,c,d;KN=ZB(ZC,GQ,6,25,15,1);LN=ZB(ZC,GQ,6,33,15,1);d=dU;for(b=32;b>=0;b--){LN[b]=d;d*=0.5;}c=1;for(a=24;a>=0;a--){KN[a]=c;c*=0.5;}}function Up(a){var b;if(a.e)return a.e;a.e=Km(new _m(),(b=a.f.a,b));!!a.e&&(a.e.M==null||a.e.M.length==0)&&pk(a.e,a.d==-1||a.b==null?null:a.b[a.d]);return a.e;}function ZN(a,b,c,d,e,f,g){var h,i;if(b.Hb()&&(i=a.a.eb(c,d),i<0||!e&&i==0)){return false;}if(b.Ib()&&(h=a.a.eb(c,f),h>0||!g&&h==0)){return false;}return true;}function tn(a,b,c){var d,e,f;f=b.length;for(e=0;e<f;e++)(a.a[b[e]]==0||a.a[b[e]]>f)&&(a.a[b[e]]=f);for(d=0;d<f;d++)(a.b[c[d]]==0||a.b[c[d]]>f)&&(a.b[c[d]]=f);}function nd(a,b,c){var d;d=b==0?LQ+a[0]-a[a.length-1]:a[b]-a[b-1];c>-2.0943951023931953&&c<dR?d-=2*$wnd.Math.cos(c+eR):d-=0.5*$wnd.Math.cos(c+eR);return d;}function yP(a,b,c,d,e,f){var g,h,i;if(VC(a)===VC(c)){a=a.slice(b,b+e);b=0;}for(h=b,i=b+e;h<i;){g=h+10000<i?h+10000:i;e=g-h;wP(c,d,f?e:0,a.slice(h,g));h=g;d+=e;}}function ro(a,b,c,d){var e;e='<circle id="'+(a.g!=null?a.g:'mol'+mo)+':Atom:'+b+_Q+QR+'cx="'+WC(c)+_Q+'cy="'+WC(d)+_Q+'r="'+8+_Q+'fill-opacity="0"/>';yM(a.a,e);}function Gd(a){var b,c,d;a.e=0;d=ZB(aG,HQ,6,a.d.d,16,1);for(b=0;b<a.d.e;b++){if(a.c[b]){++a.b;for(c=0;c<2;c++){if(!d[Ei(a.d,c,b)]){d[Ei(a.d,c,b)]=true;++a.a;}}}}}function Xd(a,b){var c,d,e;for(d=0;d<a.g[b];d++){if(a.j[b][d]!=1){c=a.f[b][d];for(e=0;e<a.g[c];e++)if(a.j[c][e]==1&&Rd(a,a.f[c][e])!=0)return true;}}return false;}function Xh(a,b){var c,d;b.t=null;b.r=null;b.I=a.I;b.o=0;for(c=0;c<a.o;c++)Vh(a,b,c,0,0);b.p=0;for(d=0;d<a.p;d++)Wh(a,b,d,0,0,null,false);Yh(a,b);!!a.b&&(b.Q=0);}function cd(a,b,c,d){var e;if(b==0){c<0?d.a=a.M:d.a=-a.M;d.b=0;return;}e=$wnd.Math.atan(c/b);b<0&&(e+=MQ);d.a=-(a.M*$wnd.Math.sin(e));d.b=a.M*$wnd.Math.cos(e);}function zg(a){var b,c;b=MQ-MQ*(a.a.length-2)/a.a.length;for(c=1;c<a.a.length;c++){a.c[c]=a.c[c-1]+$wnd.Math.sin(b*(c-1));a.d[c]=a.d[c-1]+$wnd.Math.cos(b*(c-1));}}function sC(a,b){var c,d,e,f,g,h,i,j;i=a.h>>19;j=b.h>>19;if(i!=j){return j-i;}e=a.h;h=b.h;if(e!=h){return e-h;}d=a.m;g=b.m;if(d!=g){return d-g;}c=a.l;f=b.l;return c-f;}function dH(e,a){var b=cH;if(!b){b=$doc.createElement('canvas');cH=b;}var c=''+e.b+'px '+e.a;var d=b.getContext('2d');d.font=c;var a=d.measureText(a);return a.width;}function cf(a,b){var c,d,e;if(Tk(a.L,b)<2)return false;if(bl(a.L,b)==2)return true;c=0;for(e=0;e<bl(a.L,b);e++){d=cl(a.L,b,e);Fl(a.L,d)&&(c+=Mi(a.L,d)-1);}return c>1;}function Jd(a,b){var c,d,e,f;if(Pi(a.d,b)==1){jk(a.d,b,2);a.e+=2;}for(e=0;e<2;e++){c=Ei(a.d,e,b);for(f=0;f<bl(a.d,c);f++){d=cl(a.d,c,f);if(a.c[d]){a.c[d]=false;--a.b;}}}}function Od(a){var b,c,d,e,f;for(c=0;c<a.d.e;c++){if(Mi(a.d,c)==2){for(e=0;e<2;e++){b=Ei(a.d,e,c);for(f=0;f<bl(a.d,b);f++){d=cl(a.d,b,f);if(a.c[d]){Nd(a,b);break;}}}}}}function zn(a){var b,c;a.a=null;for(b=0;b<a.d.e;b++){if(ej(a.d,b)){!a.a&&(a.a=new LM());c=new Xn();c.a=Ei(a.d,0,b);c.b=Ei(a.d,1,b);c.d=Gi(a.d,b);c.c=Fi(a.d,b);yM(a.a,c);}}}function am(a,b){var c;if((a.s[b]&yR)!=0)return true;if(a.A[b]==1)return false;return c=a.A[b],c==1||c>=5&&c<=9||c>=14&&c<=17||c>=33&&c<=35||c>=52&&c<=53||a.A[b]==13;}function EB(a,b){var c,d;c=a.b+a.n;if(a.d<c){while(a.d<c){b.a+='0';++a.d;}}else{d=a.b+a.i;d>a.d&&(d=a.d);while(d>c&&yJ(b.a,d-1)==48){--d;}if(d<a.d){XJ(b,d,a.d);a.d=d;}}}function OB(a,b,c){var d,e;d=true;while(d&&c>=0){e=yJ(b.a,c);if(e==57){DH(b,c--,48);}else{DH(b,c,e+1&AQ);d=false;}}if(d){b.a=b.a.substr(0,0)+'1'+KJ(b.a,0);++a.b;++a.d;}}function zG(b,c,d,e){yG();var f=wG;$moduleName=c;$moduleBase=d;bG=e;function g(){for(var a=0;a<f.length;a++){f[a]();}}if(b){try{eQ(g)();}catch(a){b(c,a);}}else{eQ(g)();}}function VG(a,b,c,d){TG();YG.call(this,WC(a*255+0.5),WC(b*255+0.5),WC(c*255+0.5),WC(d*255+0.5));this.b=ZB($C,gR,6,3,15,1);this.b[0]=a;this.b[1]=b;this.b[2]=c;this.a=d;}function kB(a){var b,c,d,e;b='jB';c='CA';e=mJ(a.length,5);for(d=e-1;d>=0;d--){if(DJ(a[d].d,b)||DJ(a[d].d,c)){a.length>=d+1&&(a.splice(0,d+1),undefined);break;}}return a;}function wI(a){uI==null&&(uI=/^\s*[+-]?(NaN|Infinity|((\d+\.?\d*)|(\.\d+))([eE][+-]?\d+)?[dDfF]?)\s*$/);if(!uI.test(a)){throw dG(new tJ(qQ+a+'"'));}return parseFloat(a);}function Lh(a,b,c,d,e,f,g){var h;h=fi(a,b,c);if(h==-1){a.o>=a.K&&nk(a,a.K*2);h=Ih(a,d);kh(a.H[h],b,c,0);a.v[h]=e;Hj(a,h,f);Wj(a,h,g);return true;}return Rh(a,h,d,e,f,g);}function Cd(a,b){var c,d;if(a.G.o==0)return null;Bd(a);c=a.K.c*Bi(a.G);d=new xh(a.t,b,c);if(d.c==1&&d.a==0&&d.b==0){d=null;}else{qh(d,a.K);sh(d,a.t);}yd(a,b,c,zQ);return d;}function Rl(a,b){var c,d,e,f,g,h,i,j;i=-1;d=0;for(g=0;g<2;g++){c=a.f[b][g];for(h=0;h<a.c[c];h++){e=a.f[c][h];if(e!=b){f=a.i[c][h];j=xl(a,f,e);if(d<j){d=j;i=f;}}}}return i;}function YB(a,b,c,d,e,f,g){var h,i,j,k,l;k=e[f];j=f==g-1;h=j?d:0;l=$B(h,k);d!=10&&aC(VB(a,g-f),b[f],c[f],h,l);if(!j){++f;for(i=0;i<k;++i){l[i]=YB(a,b,c,d,e,f,g);}}return l;}function Yq(b){var c,d,e;e=new xH(new AH(b));d=new Uo();while(true){try{c=wH(e);if(c!=null)So(d,c);else break;}catch(a){a=cG(a);if(PC(a,52)){break;}else throw dG(a);}}return d;}function nl(a){var b,c,d;Xo(a,1);d=0;for(b=0;b<a.o;b++){c=a.v[b]!=0?a.v[b]:Eh[a.A[b]];d+=c+ml(a,b)*Eh[1];a.A[b]>=171&&a.A[b]<=190&&a.c[b]>2&&(d-=(a.c[b]-2)*Eh[1]);}return d;}function Pc(a,b,c,d){if(fj(a.G,$k(a.G,c,d))){zd(a,-3);no(a,b);zd(a,a.J);}else if(a.o[c]!==a.o[d]){Mc(a,b,c,d);}else if(a.o[c]!=0){zd(a,a.o[c]);no(a,b);zd(a,a.J);}else{no(a,b);}}function xk(a,b,c,d){var e,f,g;for(e=0;e<a.o;e++){if(!d||(a.s[e]&512)!=0){g=a.O[e]*b;f=a.N[e]-c;a.H[e].a=a.R+g*$wnd.Math.sin(f);a.H[e].b=a.S+g*$wnd.Math.cos(f);}}d&&(a.Q&=3);}function xl(a,b,c){if(Mi(a,b)!=1)return 0;return 16-a.c[c]+(a.A[c]==1?xQ:0)+((a.F[b]&24)==0||a.B[0][b]!=c?YQ:0)+((a.s[c]&3)==0?OQ:0)+((a.C[b]&64)!=0?0:512)+(a.A[c]!=6?256:0);}function rK(a){var b,c,d;d=new UN('[',']');for(c=a.yb();c.Bb();){b=c.Cb();TN(d,b===a?'(this Collection)':b==null?rQ:vc(b));}return!d.a?d.c:d.e.length==0?d.a.a:d.a.a+(''+d.e);}function fi(a,b,c){var d,e,f,g,h,i,j,k;g=-1;e=Ci(a,a.o,a.p,Fh);h=DR;i=e*e/12;for(d=0;d<a.o;d++){j=a.H[d].a;k=a.H[d].b;f=(b-j)*(b-j)+(c-k)*(c-k);if(f<i&&f<h){h=f;g=d;}}return g;}function ho(a){$n();a=(a&1431655765)+(a>>>1&1431655765);a=(a&858993459)+(a>>>2&858993459);a=(a&117901063)+(a>>>4&117901063);a=(a&983055)+(a>>>8&983055);return(a&31)+(a>>>16);}function sp(a){var b,c,d,e;rp(a);d=a.o;c=ZB(kF,vR,2,d,6,1);for(b=0;b<d;b++){e=_o(a);Xo(e,15);Mj(e,b,(Gh(),Ch)[e.A[b]]+'*');ak(e,b,Ck('X'));tp(e);c[b]=Ze(new rf(e,8));}return c;}function lA(a,b){iA();return new XG(WC((a.c>>16&255)+FQ*((b.c>>16&255)-(a.c>>16&255))),WC((a.c>>8&255)+FQ*((b.c>>8&255)-(a.c>>8&255))),WC((a.c&255)+FQ*((b.c&255)-(a.c&255))));}function no(a,b){var c,d,e,f,g;d=WC(b.a);e=WC(b.b);f=WC(b.c);g=WC(b.d);c='<line x1="'+d+_Q+'y1="'+f+_Q+'x2="'+e+_Q+'y2="'+g+_Q+'style="stroke:'+a.d+';'+aR+WC(a.i)+'"/>';xo(a,c);}function $p(b){var c,d,e,f;d=0;Ql(b);Xo(b,3);for(c=0;c<b.d;c++){try{f=(e=rA(Yp,hJ(_d(b,c,6241))),e<0?-1:e);f!=-1&&(d+=Xp[f]);}catch(a){a=cG(a);if(!PC(a,12))throw dG(a);}}return d;}function Wc(a,b){var c,d,e;e=-1;d=-1;if((a.B&128)!=0)return-1;if(bj(a.G,b)){e=oi(a.G,b);d=ni(a.G,b);}c=Nk(a.G,b);if(c!=-1){e=Ji(a.G,c);d=Ii(a.G,c);}e!=-1&&e!=0&&(e|=d<<8);return e;}function Ie(a,b,c,d,e,f,g,h){var i,j;for(j=1;j<h;j++){for(i=g[j];i<g[j+1];i++)c[i]=f[e[i]]+(c[d[i]]<<8);Je(a,b,c,d,e,g,j);if(c[1]!==c[2])return true;j>1&&Ge(c,d,g,j);}return false;}function Gm(a,b){var c,d,e;if(b<0||b>999){WJ(a.b,'  ?');return;}c=false;for(d=0;d<3;d++){e=b/100|0;if(e==0){d==2||c?TJ(a.b,48):TJ(a.b,32);}else{TJ(a.b,48+e&AQ);c=true;}b=10*(b%100);}}function qH(a,b){var c,d,e,f,g;c=new sH();d=$wnd.Math.min(a.c,b.c);e=$wnd.Math.min(a.d,b.d);f=$wnd.Math.max(a.c+a.b,b.c+b.b);g=$wnd.Math.max(a.d+a.a,b.d+b.a);lH(c,d,e,f,g);return c;}function OJ(a){var b,c;if(a>=zQ){b=55296+(a-zQ>>10&1023)&AQ;c=56320+(a-zQ&1023)&AQ;return String.fromCharCode(b)+(''+String.fromCharCode(c));}else{return String.fromCharCode(a&AQ);}}function Sl(a,b){var c,d,e,f,g,h,i,j;i=-1;d=0;for(g=0;g<2;g++){c=a.B[g][b];for(h=0;h<a.c[c];h++){e=a.f[c][h];if(e!=a.B[1-g][b]){f=a.i[c][h];j=xl(a,f,e);if(d<j){d=j;i=f;}}}}return i;}function Zq(b){var c,d,e;e=new xH(new AH(b));c=new LM();while(true){try{d=wH(e);if(d!=null)c.a[c.a.length]=d;else break;}catch(a){a=cG(a);if(PC(a,52)){break;}else throw dG(a);}}return c;}function Vc(a,b){var c,d;if((a.B&128)!=0)return a.o[b];c=Mk(a.G,b);c!=-1&&(b=c);d=Wc(a,b);if(d==-1)return a.o[b];switch(d&255){case 1:return 384;case 2:return 64;default:return 448;}}function fL(a,b,c){var d,e,f;for(e=new qO(new vO(a).b);RK(e.a);){d=e.b=SK(e.a);f=d.Fb();if(VC(b)===VC(f)||b!=null&&pc(b,f)){if(c){d=new _L(d.Fb(),d.Gb());pO(e);}return d;}}return null;}function Sd(a,b){var c,d,e,f,g,h;if(a.q[b]==0){return false;}h=true;c=a.q[b];f=a.g[b];g=0;for(d=0;d<f;d++){e=a.f[b][d];g+=a.q[e];}(c<0?-c:c)<=(g<0?-g:g)&&nJ(c)!=nJ(g)&&(h=false);return h;}function _n(a,b){var c,d;if(!b)return null;d=ZB(_C,DQ,6,(Yn.length+31)/32|0,15,1);b=co(b);Pn(a.g,b);for(c=0;c<Yn.length;c++){On(a.g,Zn[c]);Fn(a.g,4)>0&&(d[c/32|0]|=1<<31-c%32);}return d;}function ld(a,b,c,d,e,f){var g,h,i,j,k;if(f){h=(g=(i=dH(a.e,d),new tH(0,0,i,0)).b,g);j=h/2+(a.j/8|0);k=a.j/2|0;(d=='+'||d=='-')&&(k=k*2/3);yM(a.T,new tH(b-j,c-k,2*j,2*k));}e&&po(a,d,b,c);}function Jc(a){var b,c;if((a.B&32)!=0)return;c=$o(a.G);if(c!=null){if(a.u.a==0&&a.u.b==0){b=a.K.c*Bi(a.G);Bd(a);Sc(a,b);yd(a,null,b,0);}vo(a,WC(a.v));zd(a,448);po(a,c,a.u.a,a.u.b+FQ*a.v);}}function Ke(a){var b,c;c=XB(_C,[iQ,mR],[23,7],0,[2,32],2);for(b=0;b<a.L.d;b++){a.H[b]&&(a.U[b]==1?c[0][a.T[b]]=$f(c[0][a.T[b]],b):a.U[b]==2&&(c[1][a.T[b]]=$f(c[0][a.T[b]],b)));}return c;}function eO(a,b,c,d){var e,f;f=b;e=f.c==null||a.a.eb(c.c,f.c)>0?1:0;while(f.a[e]!=c){f=f.a[e];e=a.a.eb(c.c,f.c)>0?1:0;}f.a[e]=d;d.b=c.b;d.a[0]=c.a[0];d.a[1]=c.a[1];c.a[0]=null;c.a[1]=null;}function nn(a,b){var c,d,e,f,g;f=b.length;g=ZB(_C,DQ,6,f,15,1);for(d=0;d<f;d++){c=d==f-1?b[0]:b[d+1];for(e=0;e<bl(a.f,b[d]);e++){if(al(a.f,b[d],e)==c){g[d]=cl(a.f,b[d],e);break;}}}return g;}function mA(a,b){var c;if(a==null)return b==null?0:1;if(b==null)return-1;for(c=0;c<a.length;c++){if(b.length==c)return 1;if(a[c]!==b[c])return a[c]<b[c]?-1:1;}return b.length>a.length?-1:0;}function wk(a,b,c){var d,e;e=c&103;d=Vi(a,b);switch(e){case 1:case 64:return d>=1;case 2:return d>=2;case 4:return d>=3;case 32:return mj(a,a.B[0][b])^mj(a,a.B[1][b]);default:return false;}}function ri(a,b){var c,d,e;if(a.t==null||a.t[b]==null)return(a.w[b]&1)!=0?'':Ch[a.A[b]];e='';for(d=0;d<a.t[b].length;d++){d>0&&(e=(OP(e),e+(OP(','),',')));c=a.t[b][d];e=BJ(e,Ch[c]);}return e;}function Bn(a,b,c){var d,e,f;if(a.a){for(e=new dN(a.a);e.a<e.c.a.length;){d=cN(e);if((a.u[d.a]||a.u[d.b])==c){f=sl(a.A,a.w[d.a],a.w[d.b],d.c+1,b)-1;if(f<d.d||f>d.c)return false;}}}return true;}function Un(a,b,c,d,e,f){var g,h;g=al(a.d,a.o[b],d);if(g!=a.q[b]){h=cl(a.d,a.o[b],d);if(!f[h]&&!ej(a.d,h)){a.o[++c]=g;a.q[c]=a.o[b];a.r[c]=h;f[h]=true;e[g]?a.p[c]=true:e[g]=true;}}return c;}function xC(a,b){var c,d,e;b&=63;if(b<22){c=a.l<<b;d=a.m<<b|a.l>>22-b;e=a.h<<b|a.m>>22-b;}else if(b<44){c=0;d=a.l<<b-22;e=a.m<<b-22|a.l>>44-b;}else{c=0;d=0;e=a.l<<b-44;}return fC(c&OR,d&OR,e&VT);}function zC(a,b){var c,d,e,f;b&=63;c=a.h&VT;if(b<22){f=c>>>b;e=a.m>>b|c<<22-b;d=a.l>>b|a.m<<22-b;}else if(b<44){f=0;e=c>>>b-22;d=a.m>>b-22|a.h<<44-b;}else{f=0;e=0;d=c>>>b-44;}return fC(d&OR,e&OR,f&VT);}function fJ(a){var b,c;if(gG(sQ,a)<=0&&gG(a,oQ)<=0){return tG(a).toString(16);}b=ZB(YC,gR,6,17,15,1);c=17;do{b[--c]=UH(tG(a)&15);a=hG(zC(kG(a)?rG(a):a,4));}while(gG(a,0)!=0);return QJ(b,c,17-c);}function dv(a,b){ou();var c;typeof b===TT&&(b=true);if(typeof b===lQ){c=new ru(jm(Yz(nu,false),a));b===true&&c.inventCoordinates();}else typeof b===nQ&&(c=new ru(km(Yz(nu,false),a,b)));return c;}function gB(b,c){var d,e,f,g;for(e=0,f=b.length;e<f;e++){g=b[e];try{g[1]?g[0].Mb()&&(c=fB(c,g)):g[0].Mb();}catch(a){a=cG(a);if(PC(a,14)){d=a;TA();ZA(PC(d,75)?d.ob():d);}else throw dG(a);}}return c;}function HG(a,b){var c=$wnd;if(a===''){return c;}var d=a.split('.');!(d[0]in c)&&c.execScript&&c.execScript('var '+d[0]);for(var e;d.length&&(e=d.shift());){c=c[e]=c[e]||!d.length&&b||{};}return c;}function wH(a){var b,c,d;c=vH(a);if(c==-1)return null;d=new ZJ();b=false;while(!b){if(c==10){b=true;}else if(c==13){b=true;c=vH(a);c!=10&&(a.a=c);}if(!b){if(c==-1){break;}TJ(d,c&AQ);c=vH(a);}}return d.a;}function xh(a,b,c){var d,e,f,g;th(this);e=b.b/a.b;g=b.a/a.a;f=0;f==0&&(f=24);d=f/c;this.c=$wnd.Math.min(d,$wnd.Math.min(e,g));this.a=b.c+b.b/2-this.c*(a.c+a.b/2);this.b=b.d+b.a/2-this.c*(a.d+a.a/2);}function Ih(a,b){a.o>=a.K&&nk(a,a.K*2);a.A[a.o]=0;ak(a,a.o,b);a.q[a.o]=0;a.s[a.o]=0;a.w[a.o]=0;a.u[a.o]=0;kh(a.H[a.o],0,0,0);a.t!=null&&(a.t[a.o]=null);a.r!=null&&(a.r[a.o]=null);a.Q=0;return a.o++;}function Kh(a,b){var c,d,e,f,g;a.I=a.I|b.I;d=ZB(_C,DQ,6,b.o,15,1);f=Dj(a,1);g=Dj(a,2);for(c=0;c<b.o;c++){d[c]=Vh(b,a,c,f,g);}for(e=0;e<b.p;e++){Wh(b,a,e,f,g,d,false);}a.J=a.J&&b.J;a.G=0;a.Q=0;return d;}function bo(a){var b;for(b=0;b<a.f.length;b++)if((a.c[b]&~a.f[b])!=0)return false;!a.d&&(a.d=lm(new om(false),a.e,null));!a.a&&(a.a=lm(new om(false),a.b,null));Pn(a.g,a.d);On(a.g,a.a);return Jn(a.g);}function Xg(a,b,c){var d,e,f;for(f=0;f<a.a.length;f++){e=$wnd.Math.sqrt((a.c[f]-b)*(a.c[f]-b)+(a.d[f]-c)*(a.d[f]-c));d=0-rm(b,c,a.c[f],a.d[f]);a.c[f]=b+e*$wnd.Math.sin(d);a.d[f]=c+e*$wnd.Math.cos(d);}}function eh(a,b,c,d){var e,f,g;for(g=0;g<a.a.length;g++){f=$wnd.Math.sqrt((a.c[g]-b)*(a.c[g]-b)+(a.d[g]-c)*(a.d[g]-c));e=rm(b,c,a.c[g],a.d[g])+d;a.c[g]=b+f*$wnd.Math.sin(e);a.d[g]=c+f*$wnd.Math.cos(e);}}function jo(a){$n();var b,c,d,e;if(a.length==0||(a.length&7)!=0)return null;d=ZB(_C,DQ,6,a.length/8|0,15,1);for(c=0;c<a.length;c++){e=c/8|0;b=a.charCodeAt(c)-48;b>16&&(b-=7);d[e]<<=4;d[e]+=b;}return d;}function LA(a){var b;if(a.c==null){b=VC(a.b)===VC(JA)?null:a.b;a.d=b==null?rQ:SC(b)?b==null?null:b.name:UC(b)?'String':bI(rc(b));a.a=a.a+': '+(SC(b)?b==null?null:b.message:b+'');a.c='('+a.d+') '+a.a;}}function yk(a,b,c){var d,e,f;a.R=b;a.S=c;a.N=ZB(ZC,GQ,6,a.o,15,1);a.O=ZB(ZC,GQ,6,a.o,15,1);for(d=0;d<a.o;d++){e=b-a.H[d].a;f=c-a.H[d].b;a.O[d]=$wnd.Math.sqrt(e*e+f*f);a.N[d]=Ak(b,c,a.H[d].a,a.H[d].b);}}function Bq(b){var c,d,e,f;e=-0.5299999713897705;for(c=0;c<b.d;c++){f=-1;try{f=_d(b,c,2144);}catch(a){a=cG(a);if(!PC(a,12))throw dG(a);}for(d=0;d<zq.length;d++){if(iG(yq[d],f)){e+=zq[d];break;}}}return e;}function ll(a,b,c){var d,e,f,g;e=ol(a,b)-Si(a,b);c&&(e-=a.c[b]-a.g[b]);g=a.A[b]<Dh.length?Dh[a.A[b]]:null;f=g==null?6:g[0];if(e<=f)return-1;if(g!=null)for(d=1;f<e&&d<g.length;d++)f=g[d];return f>e?f:e;}function UG(a){var b;b=ZB($C,gR,6,4,15,1);if(a.b==null){b[0]=(a.c>>16&255)/255;b[1]=(a.c>>8&255)/255;b[2]=(a.c&255)/255;b[3]=(a.c>>24&255)/255;}else{b[0]=a.b[0];b[1]=a.b[1];b[2]=a.b[2];b[3]=a.a;}return b;}function PN(a){var b,c,d,e,f,g;e=a.a*jU+a.b*1502;g=a.b*jU+11;b=$wnd.Math.floor(g*kU);e+=b;g-=b*GR;e%=GR;a.a=e;a.b=g;d=a.a*128;f=$wnd.Math.floor(a.b*LN[31]);c=d+f;c>=2147483648&&(c-=4294967296);return c;}function Mc(a,b,c,d){var e,f;e=new Fd();f=new Fd();e.a=b.a;e.c=b.c;e.b=(b.a+b.b)/2;e.d=(b.c+b.d)/2;f.a=e.b;f.c=e.d;f.b=b.b;f.d=b.d;if(vd(a,e)){zd(a,a.o[c]);no(a,e);}if(vd(a,f)){zd(a,a.o[d]);no(a,f);}zd(a,a.J);}function $h(a,b,c){var d,e;d=fi(a,b,c);if(d!=-1){(a.s[d]&512)!=0?ei(a):Zh(a,d);a.Q=0;return true;}e=gi(a,b,c);if(e!=-1){(a.s[a.B[0][e]]&a.s[a.B[1][e]]&512)!=0?ei(a):bi(a,e);a.Q=0;return true;}return false;}function oo(a,b,c,d){var e,f;f=new _J('<polygon points="');for(e=0;e<d;e++){UJ(f,WC(b[e]));f.a+=',';UJ(f,WC(c[e]));f.a+=' ';}WJ(f,'" style="fill:'+a.d+';'+'stroke:'+a.d+';'+'stroke-width:1"/>');xo(a,f.a);}function Mk(a,b){var c,d,e,f,g;c=-1;if(a.k[b]==1){for(f=0;f<a.g[b];f++){if(a.j[b][f]==2){d=a.f[b][f];if(a.g[d]==2&&a.k[d]==2){for(g=0;g<2;g++){e=a.f[d][g];if(e!=b&&a.k[e]==1){c=d;break;}}}break;}}}return c;}function oe(a,b){var c,d,e,f;e=false;for(d=0;d<a.L.e;d++)if(ee(a,d,false)){a.o[d]=a.G;b&&De(a,d);e=true;}f=false;for(c=0;c<a.L.d;c++)if(ie(a,c,false)){a.$[c]=a.G;b&&Ee(a,c);f=true;}f&&(a.G=!a.G);return e||f;}function vN(a,b,c,d,e,f){var g,h,i,j;g=d-c;if(g<7){sN(b,c,d,f);return;}i=c+e;h=d+e;j=i+(h-i>>1);vN(b,a,i,j,-e,f);vN(b,a,j,h,-e,f);if(f.eb(a[j-1],a[j])<=0){while(c<d){b[c++]=a[i++];}return;}tN(a,i,j,h,b,c,d,f);}function bi(a,b){var c,d,e;for(d=0;d<2;d++){c=0;for(e=0;e<a.p;e++){if(e==b)continue;(a.B[0][e]===a.B[d][b]||a.B[1][e]===a.B[d][b])&&++c;}if(c==0){Cj(a,a.u[a.B[d][b]]);a.A[a.B[d][b]]=-1;}}a.F[b]=128;Uh(a);a.Q=0;}function Vq(){Vq=FG;Dq=aC(VB(kF,1),vR,2,6,['mutagenic','tumorigenic','irritant','reproductive effective']);Eq=aC(VB(kF,1),vR,2,6,['Mutagenicity','Tumorigenicity','Irritating effects','Reproductive effects']);}function bQ(a){var b,c,d,e;b=0;d=a.length;e=d-4;c=0;while(c<e){b=a.charCodeAt(c+3)+31*(a.charCodeAt(c+2)+31*(a.charCodeAt(c+1)+31*(a.charCodeAt(c)+31*b)));b=b|0;c+=4;}while(c<d){b=b*31+yJ(a,c++);}b=b|0;return b;}function Qh(a,b,c){var d,e,f,g,h;e=ZB(_C,DQ,6,b.o,15,1);g=Dj(a,1);h=Dj(a,2);for(d=0;d<b.o;d++){b.A[d]!=0?e[d]=Vh(b,a,d,g,h):e[d]=c;}for(f=0;f<b.p;f++){Wh(b,a,f,g,h,e,false);}a.J=a.J&&b.J;a.G=0;a.Q=0;return e;}function Ul(a){var b,c,d,e;Xo(a,7);a.o=a.d;a.p=a.e;for(c=0;c<a.d;c++){if(a.c[c]!==a.g[c]){b=ll(a,c,false);a.c[c]=a.g[c];if(b!=-1){e=ll(a,c,true);if(b!=e){d=((a.s[c]&yR)>>>28)-1;(d==-1||d<b)&&Hj(a,c,b);}}}}$l(a);}function To(a,b){var c,d,e,f;f=a.c.a.length;if(f==0)return-1;e=1;while(2*e<=f)e<<=1;d=e;--e;while(d!=0){d>>=1;if(e>=f){e-=d;continue;}c=zJ(b,BM(a.c,e));if(c==0)return e;if(d==0)break;c<0?e-=d:e+=d;}return-1;}function pj(a){var b;for(b=0;b<a.o;b++){switch(a.A[b]){case 1:case 5:case 6:case 7:case 8:case 9:case 14:case 15:case 16:case 17:case 33:case 34:case 35:case 52:case 53:continue;default:return false;}}return true;}function po(a,b,c,d){var e,f,g,h;g=(e=(h=dH(a.e,b),new tH(0,0,h,0)).b,e);f='<text x="'+WC(c-g/2)+_Q+'y="'+WC(d+(a.j/3|0))+_Q+'font-family=" '+a.e.a+_Q+'font-size="'+a.e.b+_Q+'fill="'+a.d+'">'+b+'<\/text>';xo(a,f);}function pd(a,b){var c,d,e,f,g,h,i;c=ZB(ZC,GQ,6,Qk(a.G,b),15,1);for(e=0;e<Qk(a.G,b);e++)c[e]=Di(a.G,b,al(a.G,b,e));xN(c);f=qd(c,0);g=nd(c,0,f);for(d=1;d<c.length;d++){h=qd(c,d);i=nd(c,d,h);if(g<i){g=i;f=h;}}return f;}function Hn(a,b){var c;c=0;if((a.C[b]&512)!=0||a.F[b]==64)c|=8;else switch(Mi(a,b)){case 1:c|=1;break;case 2:c|=2;break;case 3:c|=4;}(a.C[b]&64)!=0?c|=32:a.I||(c|=16);(a.C[b]&256)!=0?c|=cR:a.I||(c|=tQ);return c;}function Wd(a,b,c){var d,e,f;d=false;for(f=0;f<a.g[b];f++){if(!Fl(a,a.i[b][f])&&a.j[b][f]==1){e=a.f[b][f];if((a.s[e]&xQ)==0&&(a.A[e]==6&&Rd(a,e)==1||a.A[e]==16&&Rd(a,e)==2)){if(d||!c)return true;d=true;}}}return false;}function Xc(a){var b,c,d,e;uo(a,2*a.L);e=new Fd();for(d=0;d<a.G.p;d++){b=Ei(a.G,0,d);c=Ei(a.G,1,d);if(dj(a.G,d)){e.a=uh(a.K,xi(a.G,b));e.c=vh(a.K,yi(a.G,b));e.b=uh(a.K,xi(a.G,c));e.d=vh(a.K,yi(a.G,c));zd(a,-2);no(a,e);}}}function so(a,b,c,d,e,f,g){var h;h='<line id="'+(a.g!=null?a.g:'mol'+mo)+':Bond:'+b+'-'+c+_Q+QR+'x1="'+WC(d)+_Q+'y1="'+WC(e)+_Q+'x2="'+WC(f)+_Q+'y2="'+WC(g)+_Q+'stroke-width="'+8+_Q+'stroke-opacity="0"'+'/>';yM(a.b,h);}function Hj(a,b,c){var d;if(c>=-1&&c<=14){a.s[b]&=268435455;c!=(d=a.A[b]<Dh.length?Dh[a.A[b]]:null,d==null?6:d[d.length-1])&&(a.s[b]|=1+c<<28);if(a.A[b]==6){if(c==-1||c==0||c==2||c==4){a.s[b]&=-49;c==2&&(a.s[b]|=16);}}}}function mm(a,b,c){var d,e,f,g;if(c==null||c.length==0){nm(a,b,null,null);return;}d=FJ(c,OJ(32));d>0&&d<c.length-1?nm(a,b,IP((f=c.substr(0,d),DP(),f)),IP((g=c.substr(d+1,c.length-(d+1)),g))):nm(a,b,IP((e=c,DP(),e)),null);}function og(a,b){var c,d;if(a==null)return b==null?0:1;if(b==null)return-1;c=mJ(a.length,b.length);for(d=0;d<c;d++)if((a[d]&rR)!=(b[d]&rR))return(a[d]&rR)<(b[d]&rR)?-1:1;return a.length==b.length?0:a.length<b.length?-1:1;}function He(a,b){var c,d,e,f,g,h,i;g=Qk(a.L,b);h=ZB(_C,DQ,6,g,15,1);for(e=0;e<g;e++)h[e]=al(a.L,b,e);for(d=g;d>1;d--){c=false;for(f=1;f<d;f++){if(Fe(a,b,h[f-1],h[f])){c=true;i=h[f-1];h[f-1]=h[f];h[f]=i;}}if(!c)break;}return h;}function Vm(a,b){var c,d,e,f,g,h;f=a.indexOf(b+'=(')+b.length+2;g=GJ(a,OJ(41),f);e=Om(a,f);c=xI(a.substr(f,e-f));h=ZB(_C,DQ,6,c,15,1);for(d=0;d<c;d++){f=Nm(a,e);e=Om(a,f);(e==-1||e>g)&&(e=g);h[d]=xI(a.substr(f,e-f));}return h;}function io(a){$n();var b,c,d,e,f,g;if(a==null)return null;b=ZB(XC,pR,6,a.length*8,15,1);for(d=0;d<a.length;d++){g=a[d];for(e=7;e>=0;e--){c=g&15;c>9&&(c+=7);b[d*8+e]=48+c<<24>>24;g>>=4;}}return PJ(GP(b,0,(f=b.length,DP(),f)));}function VI(a){var b,c,d;if(a<0){return 0;}else if(a==0){return 32;}else{d=-(a>>16);b=d>>16&16;c=16-b;a=a>>b;d=a-256;b=d>>16&8;c+=b;a<<=b;d=a-xQ;b=d>>16&4;c+=b;a<<=b;d=a-yQ;b=d>>16&2;c+=b;a<<=b;d=a>>14;b=d&~(d>>1);return c+2-b;}}function nk(a,b){var c,d;a.A=hN(a.A,b);a.q=hN(a.q,b);a.u=hN(a.u,b);d=a.H.length;a.H=iN(a.H,b);for(c=d;c<b;c++)a.H[c]=new mh();a.v=hN(a.v,b);a.s=hN(a.s,b);a.w=hN(a.w,b);a.t!=null&&(a.t=iN(a.t,b));a.r!=null&&(a.r=iN(a.r,b));a.K=b;}function Ug(a,b){var c,d,e,f,g;g=0;c=0;for(d=0;d<b;d++){g+=a[d].b*$wnd.Math.sin(a[d].a);c+=a[d].b*$wnd.Math.cos(a[d].a);}if(c==0)f=g>0?NQ:ZQ;else{f=$wnd.Math.atan(g/c);c<0&&(f+=MQ);}e=$wnd.Math.sqrt(g*g+c*c)/b;return new pm(f,e);}function Mj(a,b,c){var d,e;if(c!=null){if(c.length==0)c=null;else{d=Ck(c);if(d!=0&&DJ(c,Ch[d])||DJ(c,'?')){ak(a,b,d);c=null;}}}if(c==null){a.r!=null&&(a.r[b]=null);}else{a.r==null&&(a.r=ZB(XC,xR,13,a.K,0,2));a.r[b]=IP((e=c,DP(),e));}}function sf(a,b){var c,d,e,f;if(a.d!=b.d)return a.d>b.d?1:-1;e=a.a.length;f=b.a.length;c=e<f?e:f;for(d=0;d<c;d++){--e;--f;if(a.a[e]!==b.a[f])return a.a[e]>b.a[f]?1:-1;}if(e!=f)return e>f?1:-1;if(a.b!=b.b)return a.b>b.b?1:-1;return 0;}function Mh(a,b,c,d){var e;for(e=0;e<a.p;e++){if(a.B[0][e]==b&&a.B[1][e]==c||a.B[0][e]==c&&a.B[1][e]==b){Th(a,e,d);a.Q=0;return e;}}a.p>=a.L&&ok(a,a.L*2);a.B[0][a.p]=b;a.B[1][a.p]=c;a.F[a.p]=d;a.C[a.p]=0;a.D[a.p]=0;a.Q=0;return a.p++;}function ci(a){var b,c,d;d=false;for(b=0;b<a.o;b++){if(a.A[b]==-1){d=true;Cj(a,a.u[b]);}}for(c=0;c<a.p;c++){if(a.F[c]==128){d=true;}else if(a.A[a.B[0][c]]==-1||a.A[a.B[1][c]]==-1){a.F[c]=128;d=true;}}if(d){a.Q=0;return Uh(a);}return null;}function Qj(a,b,c,d){var e;if(c==null){a.t!=null&&(a.t[b]=null);return;}if(c.length==1&&!d){e=c[0];a.A[b]!=e&&Rh(a,b,e,0,-1,0);a.t!=null&&(a.t[b]=null);return;}a.t==null&&(a.t=ZB(_C,mR,7,a.K,0,2));a.t[b]=c;d&&(a.w[b]|=1);a.Q=0;a.I=true;}function Nh(a,b,c,d,e){var f,g,h;while(a.o+d>a.K)nk(a,a.K*2);while(a.p+d>a.L)ok(a,a.L*2);f=fi(a,b,c);if(f!=-1)return Oh(a,f,d,e);g=gi(a,b,c);if(g!=-1)return Ph(a,g,d,e);f=Hh(a,b,c,0);h=MQ*(d-2)/d;wj(a,f,d,f,e,0,MQ-h);a.Q=0;return true;}function Ro(a,b,c,d,e,f){this.k=a;if(d!=0&&d!=1){this.b=true;}else{this.a=b;this.c=c;this.d=d;this.e=f;this.g=0;this.i=ZB(aG,HQ,6,4,16,1);this.f=ZB(_C,DQ,6,4,15,1);this.j=ZB(_C,DQ,6,4,15,1);if(c!=-1&&d==1){Oo(this,oQ,e,true);this.d=0;}}}function bf(a,b){var c,d,e,f,g,h,i;i=tl(a.L);for(c=0;c<i.g.a.length;c++){if(i.d[c]&&pn(i,c,b)){for(e=BM(i.g,c),f=0,g=e.length;f<g;++f){d=e[f];if(d!=b)for(h=0;h<bl(a.L,d);h++)if(Gl(a.L,cl(a.L,d,h)))return true;}return false;}}return false;}function gq(a){var b,c,d,e,f,g;if(!eq)return NT;f=0;d=0;g=new Wn();c=new kp();for(e=0;e<dq.a.a.length;e++){mm(new om(false),c,iq(dq,e));Pn(g,a);On(g,c);if(Fn(g,g.b)>0){f+=jq(dq,e);++d;}}b=d==0?-1:f/$wnd.Math.sqrt(d);return b+'\t'+d+'\t'+a.d;}function yC(a,b){var c,d,e,f,g;b&=63;c=a.h;d=(c&tQ)!=0;d&&(c|=-1048576);if(b<22){g=c>>b;f=a.m>>b|c<<22-b;e=a.l>>b|a.m<<22-b;}else if(b<44){g=d?VT:0;f=c>>b-22;e=a.m>>b-22|c<<44-b;}else{g=d?VT:0;f=d?OR:0;e=c>>b-44;}return fC(e&OR,f&OR,g&VT);}function eg(a,b,c){var d,e,f,g,h;h=false;g=1;b[c]=1;d=true;while(d){d=false;for(e=0;e<a.b;e++){if(b[e]==g){for(f=0;f<a.b;f++){if(b[f]==0&&kg(a,e,f)){if(a.c[f]==-2){b[f]=g+1;d=true;}else if(a.c[f]!==a.c[c]){b[f]=g+1;h=true;}}}}}++g;}return h;}function wp(a,b,c){var d,e,f,g,h;g=jm(new om(true),a);e=-1;for(;0<g.o;f++){d=g.r==null?null:g.r[0]==null?null:CJ(g.r[0]);d!=null&&(h='*'.length,DJ(d.substr(d.length-h,h),'*'));e=0;break;}if(e>=0){return vp(g,e,b,c);}return ZB(kF,vR,2,0,6,1);}function kq(a){var b,c,d,e,f;f=new xH(new AH(a));this.a=new LM();while(true){e=wH(f);if(e==null)break;d=FJ(e,OJ(9));if(d==-1)throw dG(new FA('line without TAB'));b=e.substr(0,d);c=UP(wI(e.substr(d+1,e.length-(d+1))));yM(this.a,new lq(b,c));}}function Wl(a,b,c,d){var e,f,g,h,i;g=ZB(_C,DQ,6,a.c[b],15,1);i=gl(a,b,c,d,g);if(i==3)return false;f=(a.s[b]&3)==i?17:9;for(h=0;h<a.c[b];h++){if((g[h]&1)==1){e=a.i[b][c[h]];a.F[e]=f;if(a.B[0][e]!=b){a.B[1][e]=a.B[0][e];a.B[0][e]=b;}}}return true;}function rA(a,b){var c,d,e,f;f=a.a.a.length;if(f==0){return-1;}e=1;while(2*e<=f)e<<=1;d=e;--e;while(d!=0){d>>=1;if(e>=f){e-=d;continue;}c=_I(b,BM(a.a,e));if(c==0)return e;if(d==0)break;c<0?e-=d:e+=d;}e<f&&_I(b,BM(a.a,e))>0&&++e;return-(e+1);}function Pp(b,c){var d,e,f,g;g=0;e=new pp();while(g<c){try{f=wH(b.g);}catch(a){a=cG(a);if(PC(a,52)){break;}else throw dG(a);}if(f==null){break;}DJ(f.substr(0,4),NR)&&++g;if(DJ(f.substr(0,1),'>')){d=Qp(f);d!=null&&op(e,d);}}b.c=KM(e.b,ZB(kF,vR,2,0,6,1));}function IP(a){var b,c,d,e,f,g,h;g=a.length;b=0;for(f=0;f<g;){d=SH(a,f,a.length);f+=d>=zQ?2:1;d<128?++b:d<YQ?b+=2:d<zQ?b+=3:d<UQ?b+=4:d<ER&&(b+=5);}c=ZB(XC,pR,6,b,15,1);h=0;for(e=0;e<g;){d=SH(a,e,a.length);e+=d>=zQ?2:1;h+=HP(c,h,d);}return c;}function nI(a){if(a.wb()){var b=a.c;b.xb()?a.k='['+b.j:!b.wb()?a.k='[L'+b.ub()+';':a.k='['+b.ub();a.b=b.tb()+'[]';a.i=b.vb()+'[]';return;}var c=a.f;var d=a.d;d=d.split('/');a.k=qI('.',[c,qI('$',d)]);a.b=qI('.',[c,qI('.',d)]);a.i=d[d.length-1];}function Ud(a,b){var c,d,e,f,g,h;c=false;if(a.A[b]!=8)return false;if(a.g[b]!=1)return false;g=a.f[b][0];if(a.A[g]==15){h=a.g[g];for(d=0;d<h;d++){e=a.f[g][d];if(e==b){continue;}if(a.A[e]!=8){continue;}f=$k(a,g,e);if(a.F[f]==2){c=true;break;}}}return c;}function aH(a,b,c,d){var e,f;f=false;e='';if(d<0||d>255){f=true;e=' Alpha';}if(a<0||a>255){f=true;e=e+' Red';}if(b<0||b>255){f=true;e=e+' Green';}if(c<0||c>255){f=true;e=e+' Blue';}if(f){throw dG(new NI('Color parameter outside of expected range:'+e));}}function sn(a,b){var c,d,e,f,g,h;for(g=0;g<2;g++){c=Ei(a.f,g,b);if(Ai(a.f,c)==7&&bl(a.f,c)==2){d=Ei(a.f,1-g,b);for(h=0;h<bl(a.f,d);h++){e=al(a.f,d,h);f=cl(a.f,d,h);if((Ai(a.f,e)==8||Ai(a.f,e)==16)&&Mi(a.f,f)==2&&bl(a.f,e)==1)return true;}}}return false;}function DG(a,b,c){var d=BG,h;var e=d[a];var f=e instanceof Array?e[0]:null;if(e&&!f){_=e;}else{_=(h=b&&b.prototype,!h&&(h=BG[b]),GG(h));_.Kb=c;_.constructor=_;!b&&(_.Lb=IG);d[a]=_;}for(var g=3;g<arguments.length;++g){arguments[g].prototype=_;}f&&(_.Jb=f);}function Ld(a){var b,c,d,e,f,g,h;do{h=false;for(c=0;c<a.d.e;c++){if(a.c[c]){f=false;for(e=0;e<2;e++){d=Ei(a.d,e,c);b=false;for(g=0;g<bl(a.d,d);g++){if(c!=cl(a.d,d,g)&&a.c[cl(a.d,d,g)]){b=true;break;}}if(!b){f=true;break;}}if(f){h=true;Jd(a,c);}}}}while(h);}function Ko(a){var b,c,d,e,f,g,h;do{h=false;for(c=0;c<a.b.e;c++){if(a.a[c]){f=false;for(e=0;e<2;e++){b=false;d=Ei(a.b,e,c);for(g=0;g<bl(a.b,d);g++){if(c!=cl(a.b,d,g)&&a.a[cl(a.b,d,g)]){b=true;break;}}if(!b){f=true;break;}}if(f){h=true;Jo(a,c);}}}}while(h);}function zl(a,b,c){var d,e,f,g,h,i;Xo(a,1);f=ZB(_C,DQ,6,a.d,15,1);i=ZB(aG,HQ,6,a.d,16,1);f[0]=b;f[1]=c;i[b]=true;i[c]=true;e=1;g=1;while(e<=g){for(h=0;h<a.g[f[e]];h++){d=a.f[f[e]][h];if(d==b){if(e!=1)return-1;}if(!i[d]){i[d]=true;f[++g]=d;}}++e;}return g;}function rl(a,b,c){var d,e,f,g,h,i;if(b==c)return 0;Xo(a,1);g=ZB(_C,DQ,6,a.d,15,1);f=ZB(_C,DQ,6,a.d,15,1);f[0]=b;g[b]=1;e=0;h=0;while(e<=h){for(i=0;i<a.g[f[e]];i++){d=a.f[f[e]][i];if(d==c)return g[f[e]];if(g[d]==0){f[++h]=d;g[d]=g[f[e]]+1;}}++e;}return-1;}function ae(a,b,c,d,e,f,g){var h,i,j,k;j=0;for(i=0;i<a.L.d;i++)(vi(a.L,a.t[i])&e)!=0&&++j;if(j==0)return false;if(b>15){Se(a,c);b-=16;}Ne(a,1,1);Ne(a,b,4);Ne(a,j,d);for(h=0;h<a.L.d;h++){k=vi(a.L,a.t[h])&e;if(k!=0){Ne(a,h,d);f!=1&&Ne(a,k>>g,f);}}return true;}function be(a,b,c,d,e,f,g){var h,i,j,k;j=0;for(i=0;i<a.L.e;i++)(Oi(a.L,a.u[i])&e)!=0&&++j;if(j==0)return false;if(b>15){Se(a,c);b-=16;}Ne(a,1,1);Ne(a,b,4);Ne(a,j,d);for(h=0;h<a.L.e;h++){k=Oi(a.L,a.u[h])&e;if(k!=0){Ne(a,h,d);f!=1&&Ne(a,k>>g,f);}}return true;}function LB(a,b){var c,d;d=0;c=new ZJ();d+=KB(a,b,0,c,false);a.t=c.a;d+=MB(a,b,d,false);d+=KB(a,b,d,c,false);a.u=c.a;if(d<b.length&&b.charCodeAt(d)==59){++d;d+=KB(a,b,d,c,true);a.q=c.a;d+=MB(a,b,d,true);d+=KB(a,b,d,c,true);a.r=c.a;}else{a.q='-'+a.t;a.r=a.u;}}function Jh(a,b,c,d){var e;if(b==c)return-1;for(e=0;e<a.p;e++){if(a.B[0][e]==b&&a.B[1][e]==c||a.B[0][e]==c&&a.B[1][e]==b){a.F[e]<d&&(a.F[e]=d);return e;}}a.p>=a.L&&ok(a,a.L*2);a.B[0][a.p]=b;a.B[1][a.p]=c;a.F[a.p]=d;a.C[a.p]=0;a.D[a.p]=0;a.Q=0;return a.p++;}function wl(a,b){var c,d;Xo(a,1);if(a.g[b]==2&&a.j[b][0]==2&&a.j[b][1]==2){for(c=0;c<2;c++)for(d=0;d<a.c[a.f[b][c]];d++)if(tj(a,a.i[a.f[b][c]][d],a.f[b][c]))return a.i[a.f[b][c]][d];}else{for(c=0;c<a.c[b];c++)if(tj(a,a.i[b][c],b))return a.i[b][c];}return-1;}function Tf(a,b){var c,d,e,f,g,h;if(!a.b)return false;e=false;for(f=a.b.a.length-1;f>=0;f--){d=false;g=BM(a.b,f);g.a==2?d=Sf(a,g.b,g.c,g.d,b):g.a==1&&(d=Xf(a,g.b,b));if(d){GM(a.b,g);for(h=0;h<a.g[g.b].length;h++){c=a.g[g.b][h];a.n[c]=false;}e=true;}}return e;}function Hg(a,b,c,d){var e,f,g,h,i,j;g=ZB(_C,DQ,6,a.b,15,1);h=ZB(_C,DQ,6,a.b,15,1);g[0]=c;g[1]=b;h[c]=1;h[b]=2;f=1;i=1;while(f<=i){for(j=0;j<a.e[g[f]];j++){e=al(a.i,g[f],j);if(e==d)return 1+h[g[f]];if(h[e]==0&&Kl(a.i,e)){g[++i]=e;h[e]=h[g[f]]+1;}}++f;}return 0;}function Si(a,b){var c,d;if(a.A[b]>=171&&a.A[b]<=190)return 0;d=0;(a.s[b]&48)==32&&(d-=1);((a.s[b]&48)==16||(a.s[b]&48)==48)&&(d-=2);c=a.q[b];if(c==0&&a.I){(a.w[b]&PQ)==RQ&&(c=-1);(a.w[b]&PQ)==QQ&&(c=1);}a.A[b]==6?d-=c<0?-c:c:Dk(a.A[b])?d+=c:d-=c;return d;}function FB(a,b){var c,d;d=0;while(d<a.d-1&&yJ(b.a,d)==48){++d;}if(d>0){b.a=b.a.substr(0,0)+''+KJ(b.a,d);a.d-=d;a.e-=d;}if(a.j>a.o&&a.j>0){a.e+=a.b-1;c=a.e%a.j;c<0&&(c+=a.j);a.b=c+1;a.e-=c;}else{a.e+=a.b-a.o;a.b=a.o;}if(a.d==1&&b.a.charCodeAt(0)==48){a.e=0;a.b=a.o;}}function So(a,b){var c,d,e,f;f=a.c.a.length;if(f==0){xM(a.c,0,b);return 0;}e=1;while(2*e<=f)e<<=1;d=e;--e;while(d!=0){d>>=1;if(e>=f){e-=d;continue;}c=zJ(b,BM(a.c,e));if(c==0)return-1;if(d==0)break;c<0?e-=d:e+=d;}e<f&&zJ(b,BM(a.c,e))>0&&++e;xM(a.c,e,b);return e;}function ZG(a,b,c,d){TG();var e,f,g,h,i,j,k,l;g=a>b?a:b;c>g&&(g=c);h=a<b?a:b;c<h&&(h=c);f=g/255;g!=0?l=(g-h)/g:l=0;if(l==0)j=0;else{k=(g-a)/(g-h);i=(g-b)/(g-h);e=(g-c)/(g-h);a==g?j=e-i:b==g?j=2+k-e:j=4+i-k;j=j/6;j<0&&(j=j+1);}d[0]=j;d[1]=l;d[2]=f;return d;}function Ip(){Ip=FG;yp=$wnd.Math.cos(eR);Dp=$wnd.Math.sin(eR);Ap=$wnd.Math.cos(UR);Fp=$wnd.Math.sin(UR);Cp=$wnd.Math.cos(dR);Hp=$wnd.Math.sin(dR);zp=$wnd.Math.cos(sR);Ep=$wnd.Math.sin(sR);Bp=$wnd.Math.cos(IR);Gp=$wnd.Math.sin(IR);$wnd.Math.cos(VR);$wnd.Math.sin(VR);}function TG(){TG=FG;SG=new XG(255,255,255);PG=SG;new XG(192,192,192);RG=new XG(128,128,128);new XG(64,64,64);QG=new XG(0,0,0);new XG(255,0,0);new XG(255,175,175);new XG(255,200,0);new XG(255,255,0);new XG(0,255,0);new XG(255,0,255);new XG(0,255,255);new XG(0,0,255);}function nC(a){var b,c,d;c=a.l;if((c&c-1)!=0){return-1;}d=a.m;if((d&d-1)!=0){return-1;}b=a.h;if((b&b-1)!=0){return-1;}if(b==0&&d==0&&c==0){return-1;}if(b==0&&d==0&&c!=0){return WI(c);}if(b==0&&d!=0&&c==0){return WI(d)+22;}if(b!=0&&d==0&&c==0){return WI(b)+44;}return-1;}function yn(a,b,c){var d,e,f,g;if((a.D[b]&~a.g[c])!=0)return false;g=(Oi(a.d,c)&SQ)>>14;if(g!=0){if(a.A.I&&g==(Oi(a.A,c)&SQ)>>14)return true;d=false;f=tl(a.A);for(e=0;e<f.g.a.length;e++){if(BM(f.i,e).length==g){if(qn(f,e,b)){d=true;break;}}}if(!d)return false;}return true;}function vg(a){var b,c,d,e,f,g,h;c=ZB(jD,nR,67,a.b,0,1);for(b=0;b<a.b;b++){c[b]=new Gf(2);Ff(c[b],b);}h=ZB(_C,DQ,6,a.b,15,1);for(d=0;d<a.i.e;d++){e=Ni(a.i,d);if(e==1||e==2){Df(c[Ei(a.i,0,d)],e);Df(c[Ei(a.i,1,d)],e);}}f=wg(c,h);do{g=f;ug(a,c,h);f=wg(c,h);}while(g!=f);return h;}function tC(a){var b,c,d,e,f;if(isNaN(a)){return JC(),IC;}if(a<-9223372036854775808){return JC(),GC;}if(a>=9223372036854775807){return JC(),FC;}e=false;if(a<0){e=true;a=-a;}d=0;if(a>=YT){d=WC(a/YT);a-=d*YT;}c=0;if(a>=XT){c=WC(a/XT);a-=c*XT;}b=WC(a);f=fC(b,c,d);e&&lC(f);return f;}function bq(a,b,c,d,e){var f,g,h,i,j,k;f=1/(1+$wnd.Math.exp(a-5));k=1-1/(1+$wnd.Math.exp(b+5));j=1/(1+$wnd.Math.exp(0.012*c-6));g=1-1/(1+$wnd.Math.exp(d));h=(0.5+f/2)*(0.5+k/2)*(0.5+j/2)*(0.5+g/2);for(i=0;e!=null&&i<e.length;i++){e[i]==2?h*=0.8:e[i]==3&&(h*=0.6);}return h;}function Me(a,b,c,d,e,f){var g,h,i;g=c==-1?(xi(a.L,b)-xi(a.L,a.t[0]))/8:xi(a.L,b)-xi(a.L,c);h=c==-1?(yi(a.L,b)-yi(a.L,a.t[0]))/8:yi(a.L,b)-yi(a.L,c);Ne(a,WC((d+g)/e),f);Ne(a,WC((d+h)/e),f);if(a._){i=c==-1?(zi(a.L,b)-zi(a.L,a.t[0]))/8:zi(a.L,b)-zi(a.L,c);Ne(a,WC((d+i)/e),f);}}function Zk(a){var b,c,d,e,f,g,h,i;Xo(a,1);h=ZB($C,gR,6,a.d,15,1);d=ZB(_C,DQ,6,a.d,15,1);for(i=0;i<a.d;i++){d[0]=i;e=ZB(_C,DQ,6,a.d,15,1);e[i]=1;c=0;f=0;while(c<=f){for(g=0;g<a.g[d[c]];g++){b=a.f[d[c]][g];if(e[b]==0){e[b]=e[d[c]]+1;d[++f]=b;h[i]+=e[b]-1;}}++c;}h[i]/=f;}return h;}function Go(a){var b,c,d,e;for(b=0;b<a.b.d;b++){if(Ai(a.b,b)==7&&ji(a.b,b)==0&&ol(a.b,b)>3&&Tk(a.b,b)>0){for(e=0;e<bl(a.b,b);e++){c=al(a.b,b,e);d=cl(a.b,b,e);if(Mi(a.b,d)>1&&jj(a.b,c)){Pi(a.b,d)==4?jk(a.b,d,2):jk(a.b,d,1);Jj(a.b,b,ji(a.b,b)+1);Jj(a.b,c,ji(a.b,c)-1);break;}}}}}function NB(a,b){var c,d,e;if(a.b>a.d){while(a.d<a.b){b.a+='0';++a.d;}}if(!a.v){if(a.b<a.o){d=new ZJ();while(a.b<a.o){d.a+='0';++a.b;++a.d;}YJ(b,0,d.a);}else if(a.b>a.o){e=a.b-a.o;for(c=0;c<e;++c){if(yJ(b.a,c)!=48){e=c;break;}}if(e>0){b.a=b.a.substr(0,0)+''+KJ(b.a,e);a.d-=e;a.b-=e;}}}}function Wg(a){var b,c;if(a.k)return;a.n=a.c[0];a.i=a.c[0];a.o=a.d[0];a.j=a.d[0];for(b=0;b<a.a.length;b++){c=vi(a.p,a.a[b])!=0?0.6:Ai(a.p,a.a[b])!=6?0.25:0;a.n>a.c[b]-c&&(a.n=a.c[b]-c);a.i<a.c[b]+c&&(a.i=a.c[b]+c);a.o>a.d[b]-c&&(a.o=a.d[b]-c);a.j<a.d[b]+c&&(a.j=a.d[b]+c);}a.k=true;}function hl(a,b){var c,d,e,f,g,h,i,j,k;Xo(a,1);k=ZB(aG,HQ,6,a.o,16,1);h=ZB(_C,DQ,6,a.o,15,1);h[0]=b;k[b]=true;e=0;i=0;g=1;while(e<=i){for(j=0;j<a.c[h[e]];j++){d=a.f[h[e]][j];if(!k[d]){h[++i]=d;k[d]=true;++g;}}++e;}f=ZB(_C,DQ,6,g,15,1);g=0;for(c=0;c<a.o;c++)k[c]&&(f[g++]=c);return f;}function Cf(a,b,c){if(a.b==0){++a.c;a.b=63;}if(a.b==63){a.d[a.c]=oG(a.d[a.c],c);a.b-=b;}else{if(a.b>=b){a.d[a.c]=pG(a.d[a.c],b);a.d[a.c]=oG(a.d[a.c],c);a.b-=b;}else{a.d[a.c]=pG(a.d[a.c],a.b);a.d[a.c]=oG(a.d[a.c],qG(c,b-a.b));b-=a.b;++a.c;a.b=63-b;a.d[a.c]=oG(a.d[a.c],fG(c,(1<<b)-1));}}}function ug(a,b,c){var d,e,f,g,h,i,j,k;e=ZB(_C,DQ,6,16,15,1);for(d=0;d<a.b;d++){for(g=0;g<a.e[d];g++){k=c[al(a.i,d,g)];for(h=0;h<g;h++)if(k<e[h])break;for(i=g;i>h;i--)e[i]=e[i-1];e[h]=k;}j=mJ(6,a.e[d]);Ff(b[d],d);Cf(b[d],16,c[d]);Cf(b[d],(6-j)*17,0);for(f=0;f<j;f++)Cf(b[d],17,e[f]);}}function sl(a,b,c,d,e){var f,g,h,i,j,k;if(b==c)return 0;Xo(a,1);i=ZB(_C,DQ,6,a.d,15,1);h=ZB(_C,DQ,6,a.d,15,1);h[0]=b;i[b]=1;g=0;j=0;while(g<=j&&i[h[g]]<=d){for(k=0;k<a.g[h[g]];k++){f=a.f[h[g]][k];if(f==c)return i[h[g]];if(i[f]==0&&(e==null||!e[f])){h[++j]=f;i[f]=i[h[g]]+1;}}++g;}return-1;}function dg(a,b){var c,d,e,f,g,h;for(e=0;e<a.b;e++){if(a.e[e][b]&&a.c[e]!=-3){for(d=0;d<=a.j.g.length;d++){if(d!=b&&a.e[e][d]){a.e[e][b]=false;h=e<a.a?e:e<a.b?e-a.a:-1;g=lg(a,e<a.a?1:e<a.b?2:0);for(f=0;f<a.j.g[b].length;f++){c=a.j.g[b][f];Of(a.j,c)&&a.j.j[c]==h&&(a.j.j[c]=g<<24>>24);}}}}}}function uq(a){var b;b=new oq();yM(b.a,new mq('The polar surface area prediction is based on an atom-type based',2));yM(b.a,new mq('increment system, published by P. Ertl, B. Rohde, P. Selzer',2));yM(b.a,new mq('in J. Med. Chem. 2000, 43, 3714-3717',2));yM(b.a,new mq(MT,2));sq(a,b);return b;}function jf(a,b){var c,d,e;c=Ei(a.L,0,b);if(c>=a.L.d)return false;if(a.W[c]==1||a.W[c]==2)return true;if(a.W[c]==3)return false;d=Nk(a.L,c);if(d!=-1)return a.k[d]==1||a.k[d]==2;for(e=0;e<bl(a.L,c);e++){if(dl(a.L,c,e)==2){if(a.W[al(a.L,c,e)]==1||a.W[al(a.L,c,e)]==2)return true;}}return false;}function SB(a,b){var c,d,e,f,g;g=a.a.length;WJ(a,b.toPrecision(20));f=0;e=GJ(a.a,'e',g);e<0&&(e=GJ(a.a,'E',g));if(e>=0){d=e+1;d<a.a.length&&yJ(a.a,d)==43&&++d;d<a.a.length&&(f=xI(KJ(a.a,d)));XJ(a,e,a.a.length);}c=GJ(a.a,'.',g);if(c>=0){a.a=LJ(a.a,0,c)+''+KJ(a.a,c+1);f-=a.a.length-c;}return f;}function $N(a,b,c,d){var e,f;if(!b){return c;}else{e=a.a.eb(c.c,b.c);if(e==0){d.d=VL(b,c.d);d.b=true;return b;}f=e<0?0:1;b.a[f]=$N(a,b.a[f],c,d);if(_N(b.a[f])){if(_N(b.a[1-f])){b.b=true;b.a[0].b=false;b.a[1].b=false;}else{_N(b.a[f].a[f])?b=gO(b,1-f):_N(b.a[f].a[1-f])&&(b=fO(b,1-f));}}}return b;}function Tc(a,b){var c,d,e,f,g,h,i;e=ZB(aG,HQ,6,a.G.o,16,1);for(d=0;d<a.G.p;d++){if(dj(a.G,d)){e[Ei(a.G,0,d)]=true;e[Ei(a.G,1,d)]=true;}}g=new sH();for(c=0;c<a.G.o;c++){f=(vi(a.G,c)&IQ)!=0?b*0.47:e[c]?b*0.38:0;if(f!=0){h=uh(a.K,xi(a.G,c));i=vh(a.K,yi(a.G,c));rH(g,h-f,i-f,f*2,f*2);a.t=qH(a.t,g);}}}function Zg(a){var b,c,d,e,f,g,h,i;a.f=0;c=new LM();for(e=1;e<a.a.length;e++){for(f=0;f<e;f++){h=$wnd.Math.abs(a.c[e]-a.c[f]);i=$wnd.Math.abs(a.d[e]-a.d[f]);d=$wnd.Math.sqrt(h*h+i*i);if(d<0.8){b=ZB(_C,DQ,6,2,15,1);b[0]=a.a[e];b[1]=a.a[f];c.a[c.a.length]=b;}g=1-$wnd.Math.min(d,1);a.f+=g*g;}}return c;}function bh(a){var b,c,d,e,f;b=0;for(d=0;d<a.a.length;d++)for(f=0;f<a.r.e[a.a[d]];f++)al(a.p,a.a[d],f)>a.a[d]&&++b;a.e=ZB(_C,DQ,6,b,15,1);a.b=ZB(_C,DQ,6,a.r.b,15,1);b=0;for(c=0;c<a.a.length;c++){a.b[a.a[c]]=c;for(e=0;e<a.r.e[a.a[c]];e++){if(al(a.p,a.a[c],e)>a.a[c]){a.e[b]=cl(a.p,a.a[c],e);++b;}}}}function il(a,b,c){var d,e,f,g,h,i,j,k;Xo(a,1);for(e=0;e<a.o;e++)b[e]=-1;h=0;for(d=0;d<a.o;d++){if(b[d]==-1&&(!c||(a.s[d]&cR)!=0)){b[d]=h;i=ZB(_C,DQ,6,a.o,15,1);i[0]=d;g=0;j=0;while(g<=j){for(k=0;k<a.c[i[g]];k++){f=a.f[i[g]][k];if(b[f]==-1&&(!c||(a.s[f]&cR)!=0)){i[++j]=f;b[f]=h;}}++g;}++h;}}return h;}function Ok(a,b,c,d,e){var f,g,h,i,j,k;Xo(a,3);if((a.s[b]&BR)==0||c&&(a.s[b]&xQ)==0)return;i=ZB(_C,DQ,6,a.d,15,1);i[0]=b;d[b]=true;h=0;j=0;while(h<=j){for(k=0;k<a.g[i[h]];k++){g=a.i[i[h]][k];if(!e[g]&&(a.C[g]&64)!=0&&(!c||(a.C[g]&256)!=0)){e[g]=true;f=a.f[i[h]][k];if(!d[f]){d[f]=true;i[++j]=f;}}}++h;}}function kl(a){var b,c,d,e,f,g;g=ZB(_C,DQ,6,a.o,15,1);e=ZB(aG,HQ,6,a.o,16,1);for(c=0;c<a.p;c++)for(d=0;d<2;d++)Ml(a,a.B[d][c])&&!Ml(a,a.B[1-d][c])&&(e[a.B[d][c]]=true);f=a.o-1;while(f>=0&&e[f]){g[f]=f;--f;}for(b=0;b<=f;b++){if(e[b]){g[b]=f;g[f]=b;--f;while(f>=0&&e[f]){g[f]=f;--f;}}else{g[b]=b;}}return g;}function Jl(a,b){var c,d,e,f,g,h;if(Mi(a,b)!=1)return false;for(f=0;f<2;f++){c=a.B[f][b];h=a.B[1-f][b];while(a.k[c]==2&&a.g[c]==2&&a.A[c]<10){for(g=0;g<2;g++){d=a.f[c][g];if(d!=h){if(a.g[d]==1)return true;e=a.i[c][g];if(Mi(a,e)==1&&e<b)return true;h=c;c=d;break;}}}if(a.g[c]==1)return true;}return false;}function Cn(a,b){var c,d,e,f,g,h,i,j;for(e=0;e<a.d.e;e++){if((Oi(a.d,e)&qR)!=0){f=Ni(a.d,e);if(f==0)continue;c=Ei(a.d,0,e);d=Ei(a.d,1,e);if((a.u[c]||a.u[d])==b){g=a.w[c];h=a.w[d];i=$k(a.A,g,h);j=Ni(a.A,i);if(j==0)continue;if(f==3)continue;if(j==3)continue;if(In(a,e,i)==(f==j))return false;}}}return true;}function Ci(a,b,c,d){var e,f,g,h,i,j,k;for(i=0;i<c;i++)(a.D[i]&$Q)!=0&&--c;if(c==0){if(b<2)return d;k=0;j=0;for(e=1;e<b;e++){for(f=0;f<e;f++){k+=jh(a.H[e],a.H[f]);++j;}}return $wnd.Math.min(d,$wnd.Math.sqrt(b)*k/(2*j));}g=0;for(h=0;h<c;h++){(a.D[h]&$Q)==0&&(g+=jh(a.H[a.B[1][h]],a.H[a.B[0][h]]));}return g/c;}function gl(a,b,c,d,e){var f,g,h,i,j,k,l;e==null&&(e=ZB(_C,DQ,6,a.c[b],15,1));if(!fl(a,b,c,d,e))return 3;h=-1;for(i=0;i<a.c[b];i++){if((e[i]&1)==1){f=a.F[a.i[b][c[i]]];if(h!=-1&&h!=f)return 3;h=f;}}j=kJ(e[0]-e[1])==2?1:0;g=e[j]-e[j+1];l=(g<0?-g:g)==3^e[j]<e[j+1];k=a.c[b]==3||(e[3]&1)==1;return l^k^h==9?1:2;}function Rm(a){var b,c,d,e,f,g,h,i;h=null;c=a.indexOf('[');d=a.indexOf(']',c);if(c>=0&&d>0){b=ZB(_C,DQ,6,16,15,1);i=a.substr(c+1,d-(c+1));e=0;g=true;while(g&&e<16){c=i.indexOf(',');if(c==-1){f=i;g=false;}else{f=i.substr(0,c);i=i.substr(c+1,i.length-(c+1));}b[e++]=Ck(f);}h=ZB(_C,DQ,6,e,15,1);dK(b,h,e);}return h;}function ze(a){var b,c,d,e;a.G=true;a.R=ZB(XC,pR,6,a.L.d,15,1);a.f=ZB(XC,pR,6,a.L.e,15,1);e=oe(a,true);while(a.N<a.L.d&&e){for(b=0;b<a.L.d;b++){Ff(a.b[b],b);Cf(a.b[b],20,a.c[b]<<4|a.W[b]<<2);}for(c=0;c<a.L.e;c++){Df(a.b[Ei(a.L,0,c)],a.k[c]);Df(a.b[Ei(a.L,1,c)],a.k[c]);}d=ve(a);if(a.N==d)break;a.N=d;e=oe(a,true);}}function Lo(a,b){var c;if(Ai(a.b,b)==16&&ji(a.b,b)<=0||Ai(a.b,b)==6&&ji(a.b,b)!=0||!lj(a.b,b))return false;c=li(a.b,b)==null?0:mi(a.b,b)[0];if(jl(a.b,b)-c<1)return false;if(Ai(a.b,b)!=5&&Ai(a.b,b)!=6&&Ai(a.b,b)!=7&&Ai(a.b,b)!=8&&Ai(a.b,b)!=15&&Ai(a.b,b)!=16&&Ai(a.b,b)!=33&&Ai(a.b,b)!=34)return false;return true;}function yd(a,b,c,d){var e;e=c/2;switch(d&bR){case bR:if(b){a.u.a=b.c+b.b/2;a.u.b=b.d+b.a-e;break;}case 0:a.u.a=a.t.c+a.t.b/2;a.u.b=a.t.d+a.t.a+e;!!b&&a.u.b>b.d+b.a-e&&(a.u.b=b.d+b.a-e);break;case tQ:if(b){a.u.a=b.c+b.b/2;a.u.b=b.d+e;break;}case cR:a.u.a=a.t.c+a.t.b/2;a.u.b=a.t.d-e;!!b&&a.u.b<b.d+e&&(a.u.b=b.d+e);}}function xj(a,b,c){var d,e;if(Mi(a,b)!=1)return false;if((a.s[c]&3)!=0)return true;for(e=0;e<a.p;e++)if(e!=b&&a.F[e]==2&&(a.B[0][e]==c&&(a.s[a.B[1][e]]&3)!=0||a.B[1][e]==c&&(a.s[a.B[0][e]]&3)!=0))return true;for(d=0;d<a.p;d++)if(d!=b&&a.F[d]==1&&(a.B[0][d]==c||a.B[1][d]==c)&&(a.C[d]&3)!=0)return true;return false;}function $q(){Vq();if(!Pq){try{Rq=ZB(OD,iQ,58,4,0,1);Oq=ZB(IF,iQ,17,4,0,1);Qq=ZB(IF,iQ,17,4,0,1);Oq[0]=Zq(Iq);Oq[1]=Zq(Sq);Oq[2]=Zq(Fq);Oq[3]=Zq(Lq);Qq[0]=Zq(Jq);Qq[1]=Zq(Tq);Qq[2]=Zq(Gq);Qq[3]=Zq(Mq);Rq[0]=Yq(Kq);Rq[1]=Yq(Uq);Rq[2]=Yq(Hq);Rq[3]=Yq(Nq);Pq=true;}catch(a){a=cG(a);if(PC(a,12)){cK();}else throw dG(a);}}}function Wq(a,b){var c,d,e,f;if(!Pq)return 0;if(Rq[b].hb(Ze(new qf(a)))!=-1)return 3;f=new Wn();c=new kp();for(e=0;e<Oq[b].a.length;e++){mm(new om(false),c,BM(Oq[b],e));Pn(f,a);On(f,c);if(Fn(f,f.b)>0)return 3;}for(d=0;d<Qq[b].a.length;d++){mm(new om(false),c,BM(Qq[b],d));Pn(f,a);On(f,c);if(Fn(f,f.b)>0)return 2;}return 1;}function xI(a){var b,c,d,e,f;if(a==null){throw dG(new tJ(rQ));}d=a.length;e=d>0&&(a.charCodeAt(0)==45||a.charCodeAt(0)==43)?1:0;for(b=e;b<d;b++){if(TH(a.charCodeAt(b))==-1){throw dG(new tJ(qQ+a+'"'));}}f=parseInt(a,10);c=f<sQ;if(isNaN(f)){throw dG(new tJ(qQ+a+'"'));}else if(c||f>oQ){throw dG(new tJ(qQ+a+'"'));}return f;}function Ge(a,b,c,d){var e,f,g,h,i,j,k,l,m;l=c[d];g=c[d+1]-l;m=ZB(hD,fR,89,g,0,1);for(i=0;i<g;i++){m[i]=new Bf();m[i].c=a[i+l];m[i].b=b[i+l];m[i].a=i+l;}e=new yf();for(k=d;k>1;k--){for(j=0;j<g;j++){m[j].c+=a[m[j].b]<<16;m[j].b=b[m[j].b];}uN(m,0,m.length,e);f=1;for(h=0;h<g;h++){a[m[h].a]=f;h!=g-1&&xf(m[h],m[h+1])!=0&&++f;}}}function Eg(a,b,c,d){var e,f,g,h,i;f=new hh(a,a.i,b.a.length+c.a.length-d);e=0;for(h=0;h<b.a.length;h++){f.a[e]=b.a[h];f.q[e]=b.q[h];f.c[e]=b.c[h];f.d[e++]=b.d[h];}for(g=0;g<c.a.length;g++){i=_g(b,c.a[g]);if(i==-1){f.a[e]=c.a[g];f.q[e]=c.q[g];f.c[e]=c.c[g];f.d[e++]=c.d[g];}else{f.q[i]<c.q[g]&&(f.q[i]=c.q[g]);}}return f;}function GB(a,b){var c,d,e,f;if(isNaN(b)){return'NaN';}d=b<0||b==0&&1/b<0;d&&(b=-b);c=new ZJ();if(!isNaN(b)&&!isFinite(b)){WJ(c,d?a.q:a.t);c.a+='\u221E';WJ(c,d?a.r:a.u);return c.a;}b*=a.p;f=SB(c,b);e=c.a.length+f+a.i+3;if(e>0&&e<c.a.length&&yJ(c.a,e)==57){OB(a,c,e-1);f+=c.a.length-e;XJ(c,e,c.a.length);}HB(a,d,c,f);return c.a;}function Tn(a,b){var c,d,e,f;Xo(a.A,a.G);e=a.A.d;a.C=ZB(_C,DQ,6,e,15,1);a.B=ZB(_C,DQ,6,e,15,1);for(c=0;c<e;c++){a.B[c]=(Gn(a.A,c)|vi(a.A,c))&PR^PR;a.C[c]=Ai(a.A,c);(b&1)!=0&&(a.C[c]+=ji(a.A,c)+16<<8);(b&2)!=0&&(a.C[c]+=ti(a.A,c)<<16);}f=a.A.e;a.D=ZB(_C,DQ,6,f,15,1);for(d=0;d<f;d++)a.D[d]=(Hn(a.A,d)|Oi(a.A,d))&802815^786480;}function HB(a,b,c,d){var e,f,g,h,i;if(a.g){f='.'.charCodeAt(0);g=','.charCodeAt(0);}else{f='.'.charCodeAt(0);g=','.charCodeAt(0);}a.e=0;a.d=c.a.length;a.b=a.d+d;h=a.v;e=a.f;a.b>OQ&&(h=true);h&&FB(a,c);NB(a,c);PB(a,c);IB(a,c,g,e);EB(a,c);DB(a,c,f);h&&CB(a,c);i='0'.charCodeAt(0);i!=48&&JB(c,i);YJ(c,0,b?a.q:a.t);WJ(c,b?a.r:a.u);}function Bd(a){var b,c,d,e,f;e=uh(a.K,xi(a.G,0));c=uh(a.K,xi(a.G,0));f=vh(a.K,yi(a.G,0));d=vh(a.K,yi(a.G,0));for(b=0;b<a.G.o;b++){e>uh(a.K,xi(a.G,b))&&(e=uh(a.K,xi(a.G,b)));c<uh(a.K,xi(a.G,b))&&(c=uh(a.K,xi(a.G,b)));f>vh(a.K,yi(a.G,b))&&(f=vh(a.K,yi(a.G,b)));d<vh(a.K,yi(a.G,b))&&(d=vh(a.K,yi(a.G,b)));}a.t=new tH(e,f,c-e,d-f);}function DC(a){var b,c,d,e,f;if(a.l==0&&a.m==0&&a.h==0){return'0';}if(a.h==tQ&&a.m==0&&a.l==0){return'-9223372036854775808';}if(a.h>>19!=0){return'-'+DC(vC(a));}c=a;d='';while(!(c.l==0&&c.m==0&&c.h==0)){e=dC(1000000000);c=gC(c,e,true);b=''+CC(cC);if(!(c.l==0&&c.m==0&&c.h==0)){f=9-b.length;for(;f>0;f--){b='0'+b;}}d=b+d;}return d;}function DI(){DI=FG;CI=aC(VB(ZC,1),GQ,6,15,[1.3407807929942597E154,1.157920892373162E77,3.4028236692093846E38,1.8446744073709552E19,4294967296,zQ,256,16,4,2]);BI=aC(VB(ZC,1),GQ,6,15,[7.458340731200207E-155,8.636168555094445E-78,2.9387358770557188E-39,5.421010862427522E-20,2.3283064365386963E-10,dU,0.00390625,0.0625,0.25,0.5]);}function Zh(a,b){var c,d,e,f;for(c=0;c<a.p;c++){for(e=0;e<2;e++){if(a.B[e][c]==b){a.F[c]=128;d=0;for(f=0;f<a.p;f++){if(f==c)continue;(a.B[0][f]===a.B[1-e][c]||a.B[1][f]===a.B[1-e][c])&&++d;}if(d==0){Cj(a,a.u[a.B[1-e][c]]);a.A[a.B[1-e][c]]=-1;}}}}Cj(a,a.u[b]);a.A[b]=-1;a.t!=null&&(a.t[b]=null);a.r!=null&&(a.r[b]=null);Uh(a);a.Q=0;}function Xo(a,b){var c,d,e,f;Lk(a,b);if((b&~a.Q)==0)return;a.a&&(b|=128);for(c=0;c<a.o;c++)a.s[c]&=-134447112;for(d=0;d<a.e;d++)a.C[d]&=-64;e=0;f=0;if((b&16)!=0){e=16;f=1;}else if((b&32)!=0){e=32;f=3;}else if((b&64)!=0){e=64;f=5;}if((b&128)!=0){e|=128;f|=32;}a.b=new rf(a,f);mf(a.b);nf(a.b);lf(a.b);jp(a)&&(a.b=new rf(a,f));a.Q|=12|e;}function rp(a){var b,c,d,e,f;for(d=0;d<a.o;d++){f=_o(a);Xo(f,15);Mj(f,d,(Gh(),Ch)[f.A[d]]+'*');ak(f,d,Ck('X'));if(ep(f)>0){for(c=0;c<f.d;c++){if((f.s[c]&FR)!=0&&wl(f,c)==-1){e=(Xo(f,3),f.k[c]==2&&f.g[c]==2?Rl(f,c):Tl(f,c));if(e!=-1){a.F[e]=17;a.Q=0;if(a.B[1][e]==c){b=a.B[0][e];a.B[0][e]=c;a.Q=0;a.B[1][e]=b;a.Q=0;}Oj(a,c,1,0);}}}}}}function $g(a,b){var c,d,e,f;e=9999;for(c=0;c<a.a.length;c++){f=vi(a.p,a.a[c])!=0?0.6:Ai(a.p,a.a[c])!=6?0.25:0;d=0;switch(b){case 0:d=a.i-0.5*(a.i+a.o+a.c[c]-a.d[c]);break;case 1:d=a.i-0.5*(a.i-a.j+a.c[c]+a.d[c]);break;case 2:d=0.5*(a.n+a.j+a.c[c]-a.d[c])-a.n;break;case 3:d=0.5*(a.n-a.o+a.c[c]+a.d[c])-a.n;}e>d-f&&(e=d-f);}return e;}function Zc(a){var b,c,d,e;if(a.G.I){zd(a,320);if((a.B&8)!=0)for(b=0;b<a.G.d;b++)(vi(a.G,b)&-536870913)!=0&&qo(a,uh(a.K,xi(a.G,b))-a.S/2,vh(a.K,yi(a.G,b))-a.S/2,a.S);for(e=0;e<a.G.e;e++){if(Oi(a.G,e)!=0){c=Ei(a.G,0,e);d=Ei(a.G,1,e);qo(a,(uh(a.K,xi(a.G,c))+uh(a.K,xi(a.G,d))-a.S)/2,(vh(a.K,yi(a.G,c))+vh(a.K,yi(a.G,d))-a.S)/2,a.S);}}}}function Wm(a){var b,c;if(a.indexOf('[')>=0){b=a.indexOf(' NOT[');c=a.indexOf(']',b);if(b>=0&&c>0){return-(c+1);}else{b=a.indexOf(' [');c=a.indexOf(']',b);if(b>=0&&c>0){return c+1;}}b=a.indexOf(" 'NOT[");c=a.indexOf("]'",b);if(b>=0&&c>0){return-(c+2);}else{b=a.indexOf(" '[");c=a.indexOf("]'",b);if(b>=0&&c>0){return c+2;}}cK();}return 0;}function $e(a,b,c,d){var e,f,g;e=c==-1?$wnd.Math.abs(xi(a.L,b)-xi(a.L,a.t[0]))/8:$wnd.Math.abs(xi(a.L,b)-xi(a.L,c));d<e&&(d=e);f=c==-1?$wnd.Math.abs(yi(a.L,b)-yi(a.L,a.t[0]))/8:$wnd.Math.abs(yi(a.L,b)-yi(a.L,c));d<f&&(d=f);if(a._){g=c==-1?$wnd.Math.abs(zi(a.L,b)-zi(a.L,a.t[0]))/8:$wnd.Math.abs(zi(a.L,b)-zi(a.L,c));d<g&&(d=g);}return d;}function Mn(a,b){var c,d,e,f,g,h,i,j;g=false;if(Tk(a.d,b)==0){for(f=1;f<bl(a.d,b);f++){for(h=0;h<f;h++){d=al(a.d,b,f);e=al(a.d,b,h);a.w[d]>a.w[e]^d>e&&(g=!g);}}}else{for(f=0;f<bl(a.d,b);f++){c=al(a.d,b,f);j=0;i=ZB(_C,DQ,6,3,15,1);for(h=0;h<bl(a.d,c);h++){i[j]=al(a.d,c,h);i[j]!=b&&++j;}j==2&&a.w[i[0]]>a.w[i[1]]^i[0]>i[1]&&(g=!g);}}return g;}function Jf(a,b,c){var d,e,f,g,h,i,j,k,l;if(b==null)return;h=0;for(e=0;e<a.i.d;e++)b[e]&&++h;l=ZB(_C,DQ,6,h,15,1);h=0;for(d=0;d<a.i.d;d++)b[d]&&(l[h++]=d);j=false;for(g=new dN(c);g.a<g.c.a.length;){f=cN(g);if(f.length==l.length){i=false;for(k=0;k<f.length;k++){if(f[k]!==l[k]){i=true;break;}}if(!i){j=true;break;}}}j||(c.a[c.a.length]=l,true);}function Rh(a,b,c,d,e,f){if((c==1||c==151||c==152)&&ol(a,b)>1)return false;a.w[b]&=-2;a.t!=null&&(a.t[b]=null);a.r!=null&&(a.r[b]=null);if(c==a.A[b]&&d==a.v[b]&&e==((a.s[b]&yR)>>>28)-1&&f==(a.s[b]&48))return false;if(c==151||c==152){d=c-149;c=1;}a.s[b]&=960;a.A[b]=c;a.v[b]=d;a.q[b]=0;a.w[b]=0;Hj(a,b,e);Wj(a,b,f);Cj(a,a.u[b]);a.Q=0;return true;}function ul(a){var b,c,d,e,f,g,h,i,j;j=0;Xo(a,3);for(d=0;d<a.e;d++){if(Mi(a,d)==1&&(a.C[d]&64)==0){h=true;for(g=0;g<2;g++){b=a.B[g][d];if(a.g[b]==1){h=false;break;}if(a.A[b]==7&&(a.s[b]&xQ)==0){c=a.B[1-g][d];for(i=0;i<a.g[c];i++){e=a.f[c][i];f=a.i[c][i];if(f!=d&&Mi(a,f)>1&&(a.s[e]&xQ)==0&&Dk(a.A[e])){h=false;break;}}}}h&&!Jl(a,d)&&++j;}}return j;}function Yc(a){var b,c,d,e,f,g;if(a.G.I){g=a.P;zd(a,-7);for(b=0;b<a.G.d;b++)(vi(a.G,b)&IQ)!=0&&qo(a,uh(a.K,xi(a.G,b))-g,vh(a.K,yi(a.G,b))-g,2*g);uo(a,2*a.P);f=new Fd();for(e=0;e<a.G.p;e++){c=Ei(a.G,0,e);d=Ei(a.G,1,e);if((vi(a.G,c)&vi(a.G,d)&IQ)!=0){f.a=uh(a.K,xi(a.G,c));f.c=vh(a.K,yi(a.G,c));f.b=uh(a.K,xi(a.G,d));f.d=vh(a.K,yi(a.G,d));no(a,f);}}}}function he(a){var b,c,d,e,f,g,h,i,j,k;d=ZB(_C,DQ,6,16,15,1);for(b=0;b<a.L.d;b++){j=qe(a,b);for(f=0;f<j;f++){k=2*a.c[al(a.L,b,f)];c=cl(a.L,b,f);Mi(a.L,c)==2&&(Fl(a.L,c)||++k);for(h=0;h<f;h++)if(k<d[h])break;for(i=f;i>h;i--)d[i]=d[i-1];d[h]=k;}Ff(a.b[b],b);Cf(a.b[b],16,a.c[b]);for(g=j;g<a.I;g++)Cf(a.b[b],17,0);for(e=0;e<j;e++)Cf(a.b[b],17,d[e]);}}function Qc(a,b,c,d){var e,f,g,h,i,j,k,l,m,n,o,p;l=b.b-b.a;o=b.d-b.c;i=$wnd.Math.sqrt(l*l+o*o);j=2*tG(jG($wnd.Math.round(i/(4*a.R))));m=l/(j-1);p=o/(j-1);if(fj(a.G,$k(a.G,c,d))){e=-3;f=-3;}else{e=a.o[c];f=a.o[d];}k=b.a-a.R/2;n=b.c-a.R/2;zd(a,e);for(h=0;h<(j/2|0);h++){qo(a,k,n,a.R);k+=m;n+=p;}zd(a,f);for(g=0;g<(j/2|0);g++){qo(a,k,n,a.R);k+=m;n+=p;}zd(a,a.J);}function rd(a,b,c,d){var e,f,g,h;h=new Fd();if(b.a==b.b&&b.c==b.d)return;h.a=b.a;h.c=b.c;h.b=b.b;h.d=b.d;g=od(h);for(e=0;e<a.T.a.length;e++){f=BM(a.T,e);if(f.c>g.c+g.b||f.d>g.d+g.a||g.c>f.c+f.b||g.d>f.d+f.a)continue;if(sd(a,h.a,h.c,e)){if(sd(a,h.b,h.d,e))return;wd(a,h,0,e);rd(a,h,c,d);return;}if(sd(a,h.b,h.d,e)){wd(a,h,1,e);rd(a,h,c,d);return;}}Rc(a,h,c,d);}function De(b,c){var d,e,f,g,h,i;if((b.k[c]==1||b.k[c]==2)&&!Ol(b.L,c)){h=false;try{for(g=0;g<2;g++){d=Ei(b.L,g,c);if(bl(b.L,d)==3){e=ZB(_C,DQ,6,2,15,1);f=0;for(i=0;i<bl(b.L,d);i++)cl(b.L,d,i)!=c&&(e[f++]=al(b.L,d,i));b.c[e[0]]>b.c[e[1]]^Fe(b,d,e[0],e[1])&&(h=!h);}}}catch(a){a=cG(a);if(PC(a,12)){b.f[c]=3;return;}else throw dG(a);}b.k[c]==1^h?b.f[c]=1:b.f[c]=2;}}function Qp(a){var b,c,d,e;if(a.length==0||a.charCodeAt(0)!=62)return null;d=1;e=0;b=0;while(d<a.length){if(a.charCodeAt(d)==60){if(e!=0)return null;e=d;}else if(a.charCodeAt(d)==62){if(b!=0)return null;b=d;}++d;}if(e!=0&&e<b)return a.substr(e+1,b-(e+1));d=a.indexOf('DT',1);if(d==-1)return null;c=d+2;while(VH(a.charCodeAt(c)))++c;return c==d+2?null:a.substr(d,c-d);}function Nc(a,b,c,d){var e,f,g,h,i;h=(b.b-b.a)/10;i=(b.d-b.c)/10;e=new Fd();if(fj(a.G,$k(a.G,c,d))){f=-3;g=-3;}else{f=a.o[c];g=a.o[d];}zd(a,f);e.a=b.a;e.c=b.c;e.b=b.a+h*2;e.d=b.c+i*2;no(a,e);e.a=b.a+h*4;e.c=b.c+i*4;e.b=b.a+h*5;e.d=b.c+i*5;no(a,e);zd(a,g);e.a=b.a+h*5;e.c=b.c+i*5;e.b=b.a+h*6;e.d=b.c+i*6;no(a,e);e.a=b.a+h*8;e.c=b.c+i*8;e.b=b.b;e.d=b.d;no(a,e);zd(a,a.J);}function _l(b){var c,d,e,f,g,h,i,j,k,l;h=ZB(_C,DQ,6,b.o,15,1);g=il(b,h,false);if(g<=1)return null;i=ZB(_C,DQ,6,g,15,1);for(d=0;d<b.d;d++)++i[h[d]];k=0;l=i[0];for(j=1;j<g;j++){if(l<i[j]){l=i[j];k=j;}}for(c=0;c<b.o;c++)h[c]!=k&&(b.A[c]=-1);for(f=0;f<b.p;f++)h[b.B[0][f]]!=k&&(b.F[f]=128);e=Uh(b);b.Q=0;try{Hk(b,true);}catch(a){a=cG(a);if(!PC(a,12))throw dG(a);}return e;}function $o(a){var b;Xo(a,15);b=a.G&AQ;switch(a.G&rR){case zQ:return null;case qR:return b==1?'meso':''+b+' meso diastereomers';case 0:return'unknown chirality';case 196608:return'racemate';case cR:return'this enantiomer';case 327680:return'this or other enantiomer';case VQ:return'two epimers';default:return b==1?'one stereo isomer':''+b+' stereo isomers';}}function yo(a,b){var c;this.t=new sH();this.G=a;this.B=0;this.F=1;this.K=new wh();this.T=new LM();this.N=new LM();this.q=ZB(aG,HQ,6,this.G.o,16,1);this.u=new hH();this.J=0;this.w=-1;c=(TG(),PG);this.r=lA(c,yc);this.s=kA(zc,c);this.C=xc;this.D=Hc;this.i=1;this.j=10;this.k=400;this.f=400;this.d='black';this.b=new LM();this.a=new LM();this.c=new SJ();this.e=new eH(12);this.g=b;++mo;}function Th(a,b,c){var d,e,f,g,h;f=false;g=a.F[b];if(c==127){f=Xi(a,b);}else if(cm(a,b,c)){if(c==17||c==9){d=xj(a,b,a.B[0][b]);e=xj(a,b,a.B[1][b]);if(c==g){if(d==e||e){h=a.B[0][b];a.B[0][b]=a.B[1][b];a.B[1][b]=h;f=true;}}else{if(!d&&e){h=a.B[0][b];a.B[0][b]=a.B[1][b];a.B[1][b]=h;}a.F[b]=c;f=true;}}else{a.F[b]=c;f=true;}}if(f){a.Q=(g&103)==(c&103)?a.Q&3:0;a.D[b]=0;}return f;}function Oh(a,b,c,d){var e,f,g,h,i,j;if(d&&ol(a,b)>1||!d&&ol(a,b)>2)return false;f=0;e=ZB(ZC,GQ,6,4,15,1);for(h=0;h<a.p;h++){for(i=0;i<2;i++){if(a.B[i][h]==b){if(f==2){f=3;break;}e[f++]=Di(a,b,a.B[1-i][h]);}}if(f==3)break;}if(f==3)return false;j=f==1?e[0]+MQ:$wnd.Math.abs(e[0]-e[1])>MQ?(e[0]+e[1])/2:(e[0]+e[1])/2+MQ;g=MQ*(c-2)/c;wj(a,b,c,b,d,j-g/2,MQ-g);a.Q=0;return true;}function Kd(a){var b,c,d,e,f,g,h,i;for(c=0;c<a.d.e;c++){if(a.c[c]){for(e=0;e<2;e++){h=Ei(a.d,e,c);b=false;for(g=0;g<bl(a.d,h);g++){if(c!=cl(a.d,h,g)&&a.c[cl(a.d,h,g)]){b=true;break;}}if(!b){i=c;d=Ei(a.d,1-e,c);while(i!=-1){a.c[i]=false;--a.b;jk(a.d,i,64);i=-1;h=d;for(f=0;f<bl(a.d,h);f++){if(a.c[cl(a.d,h,f)]){if(i==-1){i=cl(a.d,h,f);d=al(a.d,h,f);}else{i=-1;break;}}}}break;}}}}}function Bh(a,b,c,d){var e,f,g,h;this.e=a;this.g=c;this.a=d;g=-1;for(h=0;h<Qk(this.e,this.a);h++){e=al(this.e,this.a,h);f=cl(this.e,this.a,h);if(e==this.g){Pi(this.e,f)==26&&(this.j=-1);continue;}if(tj(this.e,f,this.a)){this.i&&(a.s[d]|=qR);this.i=true;}if(g==b[e]){this.d=e;this.f=true;this.c=Ll(this.e,f);continue;}else if(g<b[e]){g=b[e];this.d=this.b;this.b=e;}else{this.d=e;}}}function Fk(a,b,c,d){var e,f,g,h,i,j,k,l,m;Xo(b,1);d==null&&(d=ZB(_C,DQ,6,b.o,15,1));h=Dj(a,1);i=Dj(a,2);m=ZB(aG,HQ,6,b.o,16,1);j=ZB(_C,DQ,6,b.o,15,1);j[0]=c;m[c]=true;d[c]=Vh(b,a,c,h,i);g=0;k=0;while(g<=k){for(l=0;l<b.c[j[g]];l++){f=b.f[j[g]][l];if(!m[f]){j[++k]=f;m[f]=true;d[f]=Vh(b,a,f,h,i);}}++g;}for(e=0;e<b.p;e++)m[b.B[0][e]]&&Wh(b,a,e,h,i,d,false);Dj(a,1);Dj(a,2);a.Q=0;}function bm(a){var b,c,d,e,f,g,h,i;f=Ci(a,a.o,a.p,Fh);g=f*f/16;for(d=1;d<a.o;d++){for(e=0;e<d;e++){h=a.H[e].a-a.H[d].a;i=a.H[e].b-a.H[d].b;if(h*h+i*i<g)throw dG(new FA('The distance between two atoms is too close.'));}}Xo(a,1);b=0;for(c=0;c<a.d;c++){if(ol(a,c)>Ui(a,c)+Si(a,c))throw dG(new FA('atom valence exceeded'));b+=a.q[c];}if(b!=0)throw dG(new FA('unbalanced atom charge'));}function Xf(a,b,c){var d,e,f,g,h,i,j,k;f=a.g[b];e=1;for(i=0;i<f.length;i++){d=f[i];if(a.f[d]&&a.k[d]==2){e=2;break;}}g=ZB(_C,mR,7,32,0,2);for(j=0;j<f.length;j++){d=f[j];a.f[d]&&a.k[d]==e&&(g[a.j[d]]=$f(g[a.j[d]],(c[d]<<16)+d));}for(k=0;k<32;k++)g[k]!=null&&xN(g[k]);zN(g,new pg());if(og(g[0],g[1])==0)return false;for(h=0;h<g[0].length;h++){d=g[0][h]&AQ;a.k[d]=0;a.j[d]=-1;}return true;}function Jg(a,b,c,d){var e,f,g,h,i,j,k,l,m;e=ZB(_C,DQ,6,d,15,1);f=0;for(g=0;g<b.a.length;g++)for(h=0;h<c.a.length;h++)b.a[g]===c.a[h]&&(e[f++]=b.a[g]);d==1?yM(a.f,(i=_g(b,e[0]),j=_g(c,e[0]),fh(c,b.c[i]-c.c[j],b.d[i]-c.d[j]),k=Sg(a,b,e[0]),l=Sg(a,c,e[0]),m=0,Bg(a,b,e[0])==1&&Bg(a,c,e[0])==1&&(m=dR),eh(c,c.c[j],c.d[j],k-l+m+MQ),Eg(a,b,c,1))):yM(a.f,Cg(a,b,c,e,d));GM(a.f,b);GM(a.f,c);}function jC(a,b,c,d,e,f){var g,h,i,j,k,l,m;j=mC(b)-mC(a);g=xC(b,j);i=fC(0,0,0);while(j>=0){h=pC(a,g);if(h){j<22?(i.l|=1<<j,undefined):j<44?(i.m|=1<<j-22,undefined):(i.h|=1<<j-44,undefined);if(a.l==0&&a.m==0&&a.h==0){break;}}k=g.m;l=g.h;m=g.l;g.h=l>>>1;g.m=k>>>1|(l&1)<<21;g.l=m>>>1|(k&1)<<21;--j;}c&&lC(i);if(f){if(d){cC=vC(a);e&&(cC=AC(cC,(JC(),HC)));}else{cC=fC(a.l,a.m,a.h);}}return i;}function Lg(a){var b,c,d,e,f,g,h;while(true){f=null;for(b=0;b<a.b;b++){h=0;for(e=0;e<a.e[b];e++)a.c[cl(a.i,b,e)]||++h;if(h==1){g=Dg(a,b);(!f||g.a.length>f.a.length)&&(f=g);}}if(!f)break;c=new hh(a,a.i,f.a.length);for(d=0;d<f.a.length;d++){a.a[f.a[d]]=true;d<f.a.length-1&&(a.c[f.b[d]]=true);c.a[d]=f.a[d];c.c[d]=$wnd.Math.cos(eR)*d;c.d[d]=(d&1)==1?0:0.5;c.q[d]=128+f.a.length;}yM(a.f,c);}}function On(a,b){var c,d,e,f;if(b.o==0||!b.I){a.d=null;return;}a.d=b;a.n=false;Xo(a.d,1);a.G=3;for(d=0;d<a.d.d;d++)(vi(a.d,d)&iR)!=0&&(a.G=7);for(f=0;f<a.d.e;f++)(Oi(a.d,f)&qR)!=0&&(a.G=7);a.F&&a.G!=3&&Xo(a.A,a.G);a.j=0;a.u=ZB(aG,HQ,6,a.d.d,16,1);for(c=0;c<a.d.d;c++){a.u[c]=(vi(a.d,c)&IQ)!=0;a.u[c]&&++a.j;}a.k=0;if(a.j!=0)for(e=0;e<a.d.e;e++)(a.u[Ei(a.d,0,e)]||a.u[Ei(a.d,1,e)])&&++a.k;}function te(a){var b,c,d,e,f,g,h,i,j,k;f=0;for(c=0;c<a.L.d;c++)a.U[c]!=0&&++f;if(f==0)return;k=ZB(_C,DQ,6,f,15,1);f=0;for(d=0;d<a.L.d;d++){if(a.U[d]!=0){k[f]=a.U[d]<<29|a.T[d]<<24|a.c[d]<<12|d;++f;}}xN(k);g=0;j=0;h=k[0]&oR;while(true){++j;if(j==k.length||h!=(k[j]&oR)){e=ZB(_C,DQ,6,j-g,15,1);for(i=g;i<j;i++){b=k[i]&4095;e[i-g]=b;a.Y[b]=true;}yM(a.Z,e);if(j==k.length)break;h=k[j]&oR;g=j;}}}function Sn(a,b,c){var d,e,f,g;f=b.d;a.e=ZB(_C,DQ,6,b.d,15,1);a.f=ZB(_C,DQ,6,b.d,15,1);for(d=0;d<f;d++){a.e[d]=(Gn(b,d)|b.w[d])&PR^PR;a.f[d]=b.A[d];(c&1)!=0&&(a.f[d]+=b.q[d]+16<<8);(c&2)!=0&&(a.f[d]+=b.v[d]<<16);}g=b.e;a.g=ZB(_C,DQ,6,b.e,15,1);for(e=0;e<g;e++){a.g[e]=(Hn(b,e)|b.D[e])&786495^786480;(c&4)!=0?(a.g[e]&2)!=0&&(a.g[e]|=8):(c&8)!=0&&(a.g[e]&2)!=0&&(b.C[e]&256)!=0&&(a.g[e]|=8);}}function Dd(a,b){var c,d,e,f;if(a.G.o==0)return null;e=(a.k=WC(b.b),a.f=WC(b.a),Cd(a,b));Xo(a.G,(a.B&256)!=0?31:(a.B&512)!=0?47:(a.B&OQ)!=0?79:15);_c(a);a.N.a=ZB(eF,fR,1,0,5,1);a.T.a=ZB(eF,fR,1,0,5,1);Lc(a);vo(a,a.Q);for(d=0;d<a.G.o;d++)gd(a,d,false);c=a.K.c*Bi(a.G);Sc(a,c);yd(a,b,c,zQ);if(kH(b,a.t))return e;f=new xh(a.t,b,c);qh(f,a.K);sh(f,a.t);rh(f,a.u);if(!e)return f;qh(f,e);return e;}function fl(a,b,c,d,e){var f,g,h;if(a.k[b]!=0||(a.s[b]&xQ)!=0||a.g[b]<3||a.c[b]>4)return false;h=ZB(aG,HQ,6,4,16,1);for(g=0;g<a.c[b];g++){f=3.9269908169872414-d[g];if($wnd.Math.abs(IR-f%NQ)>0.0872664675116539)return false;e[g]=3&WC(f/NQ);if(h[e[g]])return false;h[e[g]]=true;if((e[g]&1)==0){if(a.F[a.i[b][c[g]]]!=1)return false;}else{if(!tj(a,a.i[b][c[g]],b))return false;}}return h[0]&&h[2];}function Gl(a,b){var c,d,e,f,g,h;if(a.F[b]!=1||(a.C[b]&256)!=0||(a.C[b]&64)!=0&&(!!a.n&&b<a.e?kn(a.n,b):0)<7)return false;c=a.B[0][b];if((a.s[c]&xQ)==0||(!!a.n&&c<a.d?jn(a.n,c):0)<6)return false;d=a.B[1][b];if((a.s[d]&xQ)==0||(!!a.n&&d<a.d?jn(a.n,d):0)<6)return false;h=0;for(g=0;g<a.g[c];g++){e=a.f[c][g];e!=d&&a.g[e]>2&&++h;}for(f=0;f<a.g[d];f++){e=a.f[d][f];e!=c&&a.g[e]>2&&++h;}return h>2;}function Yi(a){var b;a.Q=0;a.A=ZB(_C,DQ,6,a.K,15,1);a.q=ZB(_C,DQ,6,a.K,15,1);a.u=ZB(_C,DQ,6,a.K,15,1);a.H=ZB(sD,{196:1,4:1,10:1,5:1,18:1,8:1},65,a.K,0,1);for(b=0;b<a.K;b++)a.H[b]=new mh();a.v=ZB(_C,DQ,6,a.K,15,1);a.s=ZB(_C,DQ,6,a.K,15,1);a.w=ZB(_C,DQ,6,a.K,15,1);a.t=null;a.r=null;a.B=XB(_C,[mR,DQ],[7,6],15,[2,a.L],2);a.F=ZB(_C,DQ,6,a.L,15,1);a.C=ZB(_C,DQ,6,a.L,15,1);a.D=ZB(_C,DQ,6,a.L,15,1);}function pl(a,b,c,d,e,f){var g,h,i,j,k,l,m,n,o;if(c==d){b[0]=c;return 0;}Xo(a,1);j=ZB(_C,DQ,6,a.d,15,1);i=ZB(_C,DQ,6,a.d,15,1);o=ZB(_C,DQ,6,a.d,15,1);i[0]=c;j[c]=1;h=0;k=0;while(h<=k&&j[i[h]]<=e){n=i[h];for(l=0;l<a.g[n];l++){if(f==null||!f[a.i[n][l]]){g=a.f[n][l];if(g==d){m=j[n];b[m]=g;b[--m]=n;while(m>0){b[m-1]=o[b[m]];--m;}return j[n];}if(j[g]==0){i[++k]=g;j[g]=j[n]+1;o[g]=n;}}}++h;}return-1;}function ml(a,b){var c,d,e,f,g,h;if(a.I&&(a.w[b]&YQ)==0)return 0;if(!am(a,b))return 0;Xo(a,1);g=0;for(e=0;e<a.c[b];e++)g+=a.j[b][e];if(a.I){c=1;for(d=0;d<a.g[b];d++)a.F[a.i[b][d]]==64&&++c;g+=c>>1;}g-=Si(a,b);f=((a.s[b]&yR)>>>28)-1;if(f==-1){if(a.A[b]>=171&&a.A[b]<=190){f=2;}else{h=a.A[b]<Dh.length?Dh[a.A[b]]:null;if(h==null){f=6;}else{f=h[0];for(d=1;f<g&&d<h.length;d++)f=h[d];}}}return 0>f-g?0:f-g;}function yl(a,b,c,d,e,f){var g,h,i,j,k;Xo(a,1);if(e){di(e);e.I=false;}i=ZB(_C,DQ,6,a.o,15,1);d==null?d=ZB(aG,HQ,6,a.o,16,1):qN(d,d.length);i[0]=b;i[1]=c;d[b]=true;d[c]=true;h=1;j=1;while(h<=j){for(k=0;k<a.c[i[h]];k++){g=a.f[i[h]][k];if(g==b){if(h!=1)return-1;}if(!d[g]){d[g]=true;i[++j]=g;}}++h;}if(e){f==null&&(f=ZB(_C,DQ,6,d.length,15,1));Jk(a,e,d,false,f);Rh(e,f[b],0,0,-1,0);}d[b]=false;return j;}function dn(a,b){var c,d,e,f,g,h,i,j,k,l,m,n,o,p;m=b.length;j=a.f.K;k=0;for(e=0;e<m;e++){if(j>b[e]){j=b[e];k=e;}}p=ZB(_C,DQ,6,m,15,1);i=k>0?k-1:m-1;l=k<m-1?k+1:0;g=b[i]<b[l];for(f=0;f<m;f++){p[f]=b[k];g?--k<0&&(k=m-1):++k==m&&(k=0);}for(d=0;d<a.g.a.length;d++){o=BM(a.g,d);if(o.length!=m)continue;c=true;for(h=0;h<m;h++){if(o[h]!==p[h]){c=false;break;}}if(c)return;}yM(a.g,p);n=nn(a,p);yM(a.i,n);tn(a,p,n);}function Bo(a,b){var c,d,e,f;f=false;a.b=b;Xo(a.b,7);c=a.b.d;d=a.b.e;a.j=ZB(aG,HQ,6,d,16,1);for(e=0;e<d;++e)a.j[e]=false;a.g=ZB(aG,HQ,6,c,16,1);a.c=ZB(_C,DQ,6,c,15,1);for(e=0;e<c;++e){a.g[e]=false;a.c[e]=-1;}a.e=ZB(kF,vR,2,3*c,6,1);a.i=0;a.d=0;a.a=0;while(!f){for(e=0;e<c;++e){if(!a.g[e]){a.a>0&&(a.e[a.i++]='.');Do(a,e,-1);++a.a;break;}}e==c&&(f=true);}a.f='';for(e=0;e<a.i;++e)a.f+=''+a.e[e];return a.f;}function Je(a,b,c,d,e,f,g){var h,i,j,k,l,m,n,o,p,q,r;for(l=g;l>1;l--){p=f[l]-f[l-1];r=ZB(fD,fR,88,p,0,1);h=f[l];for(o=0;o<p;o++){q=f[l-1]+o;m=h;while(m<f[l+1]&&d[m]==q)++m;r[o]=new wf();r[o].c=q;r[o].d=c[q];r[o].b=b[q]?0:Rk(a.L,e[q]);r[o].a=ZB(_C,DQ,6,m-h,15,1);for(k=h;k<m;k++)r[o].a[k-h]=c[k];xN(r[o].a);h=m;}i=new tf();uN(r,0,r.length,i);j=1;for(n=0;n<p;n++){c[r[n].c]=j;n!=p-1&&sf(r[n],r[n+1])!=0&&++j;}}}function Jm(a,b,c,d,e,f){var g,h,i,j;j=1;h=false;switch(e){case 1:j=17;break;case 3:j=26;break;case 4:j=17;h=true;break;case 6:j=9;break;default:switch(d){case 1:j=1;break;case 2:j=2;break;case 3:j=4;break;case 4:j=64;}}g=Jh(a.c,b,c,j);i=0;h&&Oj(a.c,b,1,-1);if(d>4){switch(d){case 5:i|=3;break;case 6:i|=9;break;case 7:i|=10;break;case 8:i|=15;}}f==1&&(i|=32);f==2&&(i|=16);i!=0&&ik(a.c,g,i,true);return g;}function jp(a){var b,c,d,e,f,g;g=false;for(c=0;c<a.d;c++)((a.s[c]&FR)==0||(a.s[c]&3)==3)&&(a.s[c]&=HR);if(a.J){if((a.G&rR)!=qR){f=ZB(aG,HQ,6,a.d,16,1);for(d=0;d<a.d;d++)(a.s[d]&FR)!=0&&(a.s[d]&3)!=3&&(a.s[d]&zR)>>19==1&&(f[d]=true);for(e=0;e<a.d;e++){if((a.s[e]&FR)!=0&&(a.s[e]&3)!=3){Oj(a,e,1,0);g=true;}}for(b=0;b<a.d;b++){if(f[b]){Uj(a,b,1,false);Oj(a,b,1,-1);g=true;}}}a.J=false;}Dj(a,1);Dj(a,2);return g;}function wd(a,b,c,d){var e,f,g,h,i,j,k,l,m,n,o;if(c==0){l=b.a;n=b.c;m=b.b;o=b.d;}else{l=b.b;n=b.d;m=b.a;o=b.c;}k=BM(a.T,d);i=m>l?k.c+k.b:k.c;j=o>n?k.d+k.a:k.d;e=m-l;f=o-n;if($wnd.Math.abs(e)>$wnd.Math.abs(f)){if(n==o){g=i;h=n;}else{g=l+e*(j-n)/f;if(m>l==i>g){h=j;}else{g=i;h=n+f*(i-l)/e;}}}else{if(l==m){g=l;h=j;}else{h=n+f*(i-l)/e;if(o>n==j>h){g=i;}else{g=l+e*(j-n)/f;h=j;}}}if(c==0){b.a=g;b.c=h;}else{b.b=g;b.d=h;}}function Fg(a,b,c){var d,e,f,g,h,i,j,k,l,m,n,o;h=ZB(_C,DQ,6,a.b,15,1);i=ZB(_C,DQ,6,a.b,15,1);j=ZB(_C,DQ,6,a.b,15,1);k=ZB(_C,DQ,6,a.b,15,1);h[0]=c;j[c]=1;k[0]=-1;g=0;l=0;while(g<=l){for(m=0;m<a.e[h[g]];m++){e=al(a.i,h[g],m);o=cl(a.i,h[g],m);if(e==b){f=j[h[g]];d=ZB(_C,DQ,6,f,15,1);d[0]=o;for(n=1;n<f;n++){d[n]=i[g];g=k[g];}return d;}if(j[e]==0){h[++l]=e;i[l]=o;j[e]=j[h[g]]+1;k[l]=g;}}if(g==l)return null;++g;}return null;}function ye(a){var b,c,d,e,f,g;a.G=true;f=oe(a,false);while(a.N<a.L.d&&f){for(b=0;b<a.L.d;b++){Ff(a.b[b],b);Cf(a.b[b],16,a.c[b]);g=a.W[b]<<7;if((a.W[b]==1||a.W[b]==2)&&a.U[b]!=0){g|=a.U[b]<<5;g|=a.T[b];}Cf(a.b[b],18,g<<9);}for(c=0;c<a.L.e;c++){d=a.k[c]<<7;if((a.k[c]==1||a.k[c]==2)&&Pi(a.L,c)==1&&a.j[c]!=0){d|=a.j[c]<<5;d|=a.i[c];}Df(a.b[Ei(a.L,0,c)],d);Df(a.b[Ei(a.L,1,c)],d);}e=ve(a);if(a.N==e)break;a.N=e;f=oe(a,false);}}function hn(a,b,c){var d,e,f,g,h,i,j,k,l,m,n,o;e=Ei(a.f,0,b);f=Ei(a.f,1,b);i=ZB(_C,DQ,6,a.f.d,15,1);j=ZB(_C,DQ,6,a.f.d,15,1);k=ZB(_C,DQ,6,a.f.d,15,1);i[0]=e;i[1]=f;j[e]=1;j[f]=2;k[e]=-1;k[f]=e;h=1;l=1;while(h<=l){for(m=0;m<bl(a.f,i[h]);m++){g=al(a.f,i[h],m);if(h>1&&g==e){o=ZB(_C,DQ,6,j[i[h]],15,1);d=i[h];for(n=0;n<o.length;n++){o[n]=d;d=k[d];}return o;}if(j[g]==0&&!c[g]){i[++l]=g;j[g]=j[i[h]]+1;k[g]=i[h];}}++h;}return null;}function en(a,b,c){var d,e,f,g,h,i,j,k,l,m,n,o;e=Ei(a.f,0,b);f=Ei(a.f,1,b);i=ZB(_C,DQ,6,a.f.d,15,1);j=ZB(_C,DQ,6,a.f.d,15,1);k=ZB(_C,DQ,6,a.f.d,15,1);i[0]=e;i[1]=f;j[e]=1;j[f]=2;k[e]=-1;k[f]=e;h=1;l=1;while(h<=l){if(j[i[h]]>7)return;for(m=0;m<bl(a.f,i[h]);m++){g=al(a.f,i[h],m);if(h>1&&g==e){o=ZB(_C,DQ,6,j[i[h]],15,1);d=i[h];for(n=0;n<o.length;n++){o[n]=d;d=k[d];}dn(a,o);continue;}if(j[g]==0&&!c[g]){i[++l]=g;j[g]=j[i[h]]+1;k[g]=i[h];}}++h;}}function Ng(a,b,c){var d,e,f,g,h,i,j,k,l,m;for(f=0;f<a.d;f++){d=Ei(a.i,0,f);e=Ei(a.i,1,f);if(Ll(a.i,f)||Mi(a.i,f)!=1||a.e[d]==1||a.e[e]==1)continue;if((a.g&2)!=0&&lj(a.i,d)&&lj(a.i,e))continue;l=false;for(j=0;j<2;j++){g=Ei(a.i,j,f);if(a.e[g]>2){m=true;i=-1;for(k=0;k<a.e[g];k++){h=al(a.i,g,k);h!=Ei(a.i,1-j,f)&&(i==-1?i=c[h]:i!=c[h]&&(m=false));}if(m){l=true;break;}}}l||((a.g&4)!=0&&lj(a.i,d)&&lj(a.i,e)?b[f]=1:b[f]=2);}}function mg(a,b){var c,d,e,f,g,h,i,j,k,l;k=oQ;i=-1;l=-1;j=-1;for(d=0;d<a.j.i.d;d++){if(Of(a.j,d)&&a.j.k[d]!=0){for(h=0;h<b.length;h++){e=b[h]<a.a?b[h]:b[h]<a.b?b[h]-a.a:-1;f=b[h]<a.a?1:b[h]<a.b?2:0;if(a.j.k[d]==f&&a.j.j[d]==e){if(k>a.j.a[d]+(f==1?zQ:0)){k=a.j.a[d]+(f==1?zQ:0);i=e;l=f;j=b[h];}}}}}for(c=0;c<a.j.i.d;c++){if(Of(a.j,c)&&a.j.k[c]==l&&a.j.j[c]==i){a.j.k[c]=0;a.j.j[c]=-1;}}for(g=0;g<a.j.g.length;g++)a.e[j][g]=false;}function gi(a,b,c){var d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s;m=-1;o=Ci(a,a.o,a.p,Fh);n=DR;for(d=0;d<a.p;d++){p=a.H[a.B[0][d]].a;r=a.H[a.B[0][d]].b;q=a.H[a.B[1][d]].a;s=a.H[a.B[1][d]].b;k=q-p;l=s-r;e=$wnd.Math.sqrt(k*k+l*l);f=(p+q)/2;g=(r+s)/2;k=b-f;l=c-g;if($wnd.Math.sqrt(k*k+l*l)>e/2)continue;if(q==p)j=$wnd.Math.abs(p-b);else{h=(s-r)/(p-q);i=-h*p-r;j=$wnd.Math.abs((h*b+c+i)/$wnd.Math.sqrt(h*h+1));}if(j<o&&j<n){n=j;m=d;}}return m;}function Wk(a,b,c){var d,e,f,g,h,i,j,k,l,m,n,o;Xo(a,3);f=ZB(aG,HQ,6,a.e,16,1);l=ZB(aG,HQ,6,a.e,16,1);o=ZB(_C,DQ,6,a.d,15,1);g=0;for(h=1;h<a.g[b];h++){d=a.i[b][h];if((a.C[d]&64)!=0){for(j=0;j<h;j++){e=a.i[b][j];if((a.C[e]&64)!=0){l[d]=true;l[e]=true;n=pl(a,o,a.f[b][h],a.f[b][j],c-2,l);l[d]=false;l[e]=false;if(n!=-1){i=false;m=ZB(_C,DQ,6,n,15,1);ql(a,o,m,n);for(k=0;k<n;k++){if(!f[m[k]]){f[m[k]]=true;i=true;}}i&&++g;}}}}}return g;}function Xi(a,b){var c,d,e;d=Vi(a,b);c=mj(a,a.B[0][b])||mj(a,a.B[1][b]);e=c?32:1;if(a.F[b]==4){a.F[b]=e;a.Q=0;return true;}if(a.F[b]==2){a.F[b]=26;a.Q&=3;if((a.C[b]&128)==0)return true;}if(a.F[b]==26){d==3?a.F[b]=4:a.F[b]=e;a.Q=0;return true;}if((24&a.F[b])!=0){a.F[b]=1;a.Q&=3;return true;}if(!c&&d<2)return false;if(a.F[b]==1){a.F[b]=2;a.Q=0;return true;}if(d<1)return false;if(a.F[b]==32){a.F[b]=1;a.Q=0;return true;}return false;}function ue(a){var b,c,d,e,f,g,h,i,j,k,l;e=false;for(f=0;f<a.Z.a.length;f++){d=BM(a.Z,f);b=true;l=-1;g=false;for(j=0;j<d.length;j++){c=d[j];if(a.W[c]==0){b=false;break;}if(a.W[c]!=3){h=true;for(k=0;k<d.length;k++){if(k!=j&&a.c[c]===a.c[d[k]]){h=false;break;}}if(h&&l<a.c[c]){l=a.c[c];g=a.W[c]==1;}}}if(b&&l!=-1){for(i=0;i<d.length;i++){c=d[i];g&&(a.W[c]==1?a.W[c]=2:a.W[c]==2&&(a.W[c]=1));a.Y[c]=false;}GM(a.Z,d);e=true;--f;}}return e;}function Sf(a,b,c,d,e){var f,g,h,i,j,k;i=null;f=null;for(k=0;k<a.g[b].length;k++){g=a.g[b][k];a.f[g]&&(a.o[g]==1||a.o[g]==2)&&(a.k[g]==0?f=$f(f,(e[g]<<16)+g):a.k[g]==d&&a.j[g]==c&&(i=$f(i,(e[g]<<16)+g)));}h=og(i,f);if(h==0)return false;if(h<0){for(j=0;j<a.g[b].length;j++){g=a.g[b][j];if(a.f[g]&&(a.o[g]==1||a.o[g]==2)){if(a.k[g]==0){a.k[g]=d<<24>>24;a.j[g]=c<<24>>24;}else if(a.k[g]==d&&a.j[g]==c){a.k[g]=0;a.j[g]=-1;}}}}return true;}function Dg(a,b){var c,d,e,f,g,h,i,j,k,l,m,n;e=ZB(_C,DQ,6,a.b,15,1);f=ZB(_C,DQ,6,a.b,15,1);g=ZB(_C,DQ,6,a.b,15,1);h=ZB(_C,DQ,6,a.b,15,1);e[0]=b;g[b]=1;h[0]=-1;d=0;i=0;while(d<=i){if(d==0||!a.a[e[d]]){for(j=0;j<a.e[e[d]];j++){c=al(a.i,e[d],j);m=cl(a.i,e[d],j);if(g[c]==0&&!a.c[m]){e[++i]=c;f[i]=m;g[c]=g[e[d]]+1;h[i]=d;}}}if(d==i){n=new sm(g[e[d]]);k=d;for(l=0;l<n.a.length;l++){n.a[l]=e[k];n.b[l]=f[k];k=h[k];}return n;}++d;}return null;}function OG(){var a=navigator.userAgent.toLowerCase();var b=$doc.documentMode;if(function(){return a.indexOf('webkit')!=-1;}())return _T;if(function(){return a.indexOf('msie')!=-1&&b>=10&&b<11;}())return'ie10';if(function(){return a.indexOf('msie')!=-1&&b>=9&&b<11;}())return'ie9';if(function(){return a.indexOf('msie')!=-1&&b>=8&&b<11;}())return'ie8';if(function(){return a.indexOf('gecko')!=-1||b>=11;}())return'gecko1_8';return'unknown';}function Pk(a){var b,c,d,e,f,g,h,i;a.n=new vn(a,7);c=ZB(_C,DQ,6,a.d,15,1);for(d=0;d<a.e;d++){if(kn(a.n,d)!=0){a.C[d]|=64;++c[a.B[0][d]];++c[a.B[1][d]];}}for(b=0;b<a.d;b++){c[b]==2?a.s[b]|=OQ:c[b]==3?a.s[b]|=YQ:c[b]>3&&(a.s[b]|=BR);}for(i=0;i<a.n.g.a.length;i++){f=ln(a.n,i);h=mn(a.n,i);g=f.length;for(e=0;e<g;e++){a.s[f[e]]|=8;a.C[h[e]]|=128;if(on(a.n,i)){a.s[f[e]]|=xQ;a.C[h[e]]|=256;}rn(a.n,i)&&(a.C[h[e]]|=512);a.F[h[e]]==26&&(a.F[h[e]]=2);}}}function uk(a,b,c){var d,e,f,g,h;g=a.A[b];a.A[b]=a.A[c];a.A[c]=g;g=a.q[b];a.q[b]=a.q[c];a.q[c]=g;g=a.v[b];a.v[b]=a.v[c];a.v[c]=g;g=a.s[b];a.s[b]=a.s[c];a.s[c]=g;g=a.w[b];a.w[b]=a.w[c];a.w[c]=g;g=a.u[b];a.u[b]=a.u[c];a.u[c]=g;f=a.H[b];a.H[b]=a.H[c];a.H[c]=f;if(a.t!=null){h=a.t[b];a.t[b]=a.t[c];a.t[c]=h;}if(a.r!=null){h=a.r[b];a.r[b]=a.r[c];a.r[c]=h;}for(d=0;d<a.p;d++){for(e=0;e<2;e++){a.B[e][d]==b?a.B[e][d]=c:a.B[e][d]==c&&(a.B[e][d]=b);}}a.Q=0;}function In(a,b,c){var d,e,f,g,h,i,j,k,l,m;h=false;for(g=0;g<2;g++){d=Ei(a.d,g,b);k=a.w[d];if(bl(a.d,d)==2){if(bl(a.A,k)==2)continue;e=-1;for(j=0;j<2;j++)cl(a.d,d,j)!=b&&(e=al(a.d,d,j));m=0;l=ZB(_C,DQ,6,2,15,1);for(i=0;i<3;i++)cl(a.A,k,i)!=c&&(l[m++]=al(a.A,k,i));a.w[e]!==l[0]&&(h=!h);}else if(bl(a.d,d)==3&&bl(a.A,k)==3){e=ZB(_C,DQ,6,2,15,1);f=0;for(i=0;i<3;i++)cl(a.d,d,i)!=b&&(e[f++]=al(a.d,d,i));a.w[e[0]]>a.w[e[1]]^e[0]>e[1]&&(h=!h);}}return h;}function vp(a,b,c,d){var e,f,g,h,i,j,k,l,m,n,o,p;i=new lp(a.d,a.e);k=new _O();n=0;m=0;g=ZB(aG,HQ,6,a.d,16,1);f=ZB(_C,DQ,6,a.d,15,1);for(p=0;p<c&&m<a.d;p++){if(m==0){f[0]=b;g[b]=true;m=1;}else{o=m;for(j=n;j<m;j++){e=f[j];for(l=0;l<a.g[e];l++){h=a.f[e][l];if(!g[h]){switch(d){case 0:g[h]=true;f[o++]=h;break;case 1:if(!(xp(a,e)&&xp(a,h))){g[h]=true;f[o++]=h;}}}}}n=m;m=o;}Jk(a,i,g,true,null);ZO(k,Ze(new rf(i,8)));}return $O(k,ZB(kF,vR,2,k.a.a.length,6,1));}function Jk(a,b,c,d,e){var f,g,h,i,j;d&&Xo(a,3);b.t=null;a.I&&lk(b,true);i=c.length;e==null&&(e=ZB(_C,DQ,6,i,15,1));b.o=0;for(f=0;f<i;f++)e[f]=c[f]?Vh(a,b,f,0,0):-1;b.p=0;for(j=0;j<a.p;j++){g=a.B[0][j];h=a.B[1][j];if(g<i&&h<i){if(c[g]&&c[h])Wh(a,b,j,0,0,e,d);else if(a.q[g]!=0&&a.q[h]!=0&&a.q[g]<0^a.q[h]<0){c[g]&&(a.q[g]+=a.q[g]<0?1:-1);c[h]&&(a.q[h]+=a.q[h]<0?1:-1);}}}Yh(a,b);!!a.b&&(b.Q=0);b.Q=0;Dj(b,1);Dj(b,2);b.o!=i&&lk(b,true);d&&Id(new Pd(b));}function Ah(a){var b,c,d,e,f,g;if(a.j!=0)return a.j;if(a.i&&Ai(a.e,a.a)!=15&&Ai(a.e,a.a)!=16){for(g=0;g<Qk(a.e,a.a);g++){f=cl(a.e,a.a,g);if(tj(a.e,f,a.a)){al(a.e,a.a,g)==a.b?a.j=Pi(a.e,f)==17?3:1:a.j=Pi(a.e,f)==17?1:3;return a.j;}}}b=Di(a.e,a.a,a.g);d=Di(a.e,a.a,a.b);d<b&&(d+=LQ);if(Qk(a.e,a.a)==2){c=d-b;if(c>3.0915926535897933&&c<3.191592653589793){a.j=-1;return a.j;}a.j=c<MQ?4:2;return a.j;}else{e=Di(a.e,a.a,a.d);e<b&&(e+=LQ);a.j=e<d?2:4;return a.j;}}function me(a){var b,c,d,e,f,g,h,i,j,k,l,m;if(a.s)return;a.s=new LM();k=0;l=ZB(_C,DQ,6,a.L.d,15,1);g=ZB(_C,DQ,6,a.L.d,15,1);i=ZB(_C,DQ,6,a.L.e,15,1);for(b=0;b<a.L.d;b++){if(l[b]==0&&(Kl(a.L,b)||Tk(a.L,b)==1)){g[0]=b;h=1;j=0;l[b]=++k;c=ZB(aG,HQ,6,a.L.e,16,1);for(f=0;f<h;f++){for(m=0;m<qe(a,g[f]);m++){e=cl(a.L,g[f],m);if(Ll(a.L,e)||Mi(a.L,e)==2||Gl(a.L,e)){d=al(a.L,g[f],m);if(!c[e]){i[j++]=e;c[e]=true;}if(l[d]==0){g[h++]=d;l[d]=k;}}}}yM(a.s,new If(g,h,i,j));}}}function Rc(a,b,c,d){var e,f,g,h,i,j,k,l;k=(b.c-b.d)/9;l=(b.b-b.a)/9;g=ZB(ZC,GQ,6,3,15,1);h=ZB(ZC,GQ,6,3,15,1);i=ZB(ZC,GQ,6,4,15,1);j=ZB(ZC,GQ,6,4,15,1);g[0]=b.a;h[0]=b.c;i[2]=b.b+k;j[2]=b.d+l;i[3]=b.b-k;j[3]=b.d-l;g[1]=(g[0]+i[2])/2;h[1]=(h[0]+j[2])/2;g[2]=(g[0]+i[3])/2;h[2]=(h[0]+j[3])/2;i[0]=g[2];j[0]=h[2];i[1]=g[1];j[1]=h[1];if(fj(a.G,$k(a.G,c,d))){e=-3;f=-3;}else{e=a.o[c];f=Vc(a,c);e==ki(a.G,c)&&(e=f);}zd(a,e);oo(a,g,h,3);zd(a,f);oo(a,i,j,4);zd(a,a.J);}function gh(a,b){var c,d;this.r=a;this.p=b.p;this.a=ZB(_C,DQ,6,b.a.length,15,1);this.q=ZB(_C,DQ,6,b.a.length,15,1);this.c=ZB(ZC,GQ,6,b.a.length,15,1);this.d=ZB(ZC,GQ,6,b.a.length,15,1);for(d=0;d<b.a.length;d++){this.a[d]=b.a[d];this.q[d]=b.q[d];this.c[d]=b.c[d];this.d[d]=b.d[d];}if(b.e!=null){this.e=ZB(_C,DQ,6,b.e.length,15,1);for(c=0;c<b.e.length;c++)this.e[c]=b.e[c];}if(b.b!=null){this.b=ZB(_C,DQ,6,b.b.length,15,1);for(c=0;c<b.b.length;c++)this.b[c]=b.b[c];}}function Tl(a,b){var c,d,e,f,g,h,i,j,k,l,m,n,o,p,q;c=ZB(ZC,GQ,6,a.c[b],15,1);for(l=0;l<a.c[b];l++)c[l]=Di(a,b,a.f[b][l]);for(m=1;m<a.c[b];m++){for(n=0;n<m;n++){d=$wnd.Math.abs(Bk(c[m],c[n]));if(d<0.08||d>JR){e=0;f=0;for(o=0;o<a.c[b];o++){if(o!=m&&o!=n){e+=$wnd.Math.abs(gA(c[m],c[o]));f+=$wnd.Math.abs(gA(c[n],c[o]));}}h=e<f?a.i[b][m]:a.i[b][n];if(Mi(a,h)==1)return h;}}}p=-1;g=0;for(k=0;k<a.c[b];k++){i=a.f[b][k];j=a.i[b][k];q=xl(a,j,i);if(g<q){g=q;p=j;}}return p;}function Fo(a,b){var c,d,e,f,g,h,i,j,k,l,m,n,o;j=ZB(_C,DQ,6,a.b.d,15,1);h=ZB(_C,DQ,6,a.b.d,15,1);i=ZB(_C,DQ,6,a.b.d,15,1);k=ZB(_C,DQ,6,a.b.d,15,1);c=Ei(a.b,0,b);d=Ei(a.b,1,b);h[0]=c;h[1]=d;i[0]=-1;i[1]=b;j[c]=1;j[d]=2;k[c]=-1;k[d]=c;g=1;l=1;while(g<=l&&j[h[g]]<15){o=h[g];for(m=0;m<bl(a.b,o);m++){e=al(a.b,o,m);if(e!=k[o]){f=cl(a.b,o,m);if(e==c){i[0]=f;for(n=0;n<=l;n++)a.a[i[m]]=true;return;}if(lj(a.b,e)&&j[e]==0){++l;h[l]=e;i[l]=f;j[e]=j[o]+1;k[e]=o;}}}++g;}return;}function vd(a,b){var c,d,e,f,g,h;if(b.a==b.b&&b.c==b.d){for(e=0;e<a.T.a.length;e++){g=BM(a.T,e);if(mH(g,b.a,b.c))return false;}return true;}h=od(b);c=false;if(b.a>b.b){md(b);c=true;}for(d=0;d<a.T.a.length;d++){g=BM(a.T,d);if(g.c>h.c+h.b||g.d>h.d+h.a||h.c>g.c+g.b||h.d>g.d+g.a)continue;if(sd(a,b.a,b.c,d)){if(sd(a,b.b,b.d,d)){c&&md(b);return false;}wd(a,b,0,d);f=vd(a,b);c&&md(b);return f;}if(sd(a,b.b,b.d,d)){wd(a,b,1,d);f=vd(a,b);c&&md(b);return f;}}c&&md(b);return true;}function Wh(a,b,c,d,e,f,g){var h,i,j,k,l;i=b.p;i>=b.L&&ok(b,b.L*2);k=(a.C[c]&BR)>>10;j=-1;k==1&&(d==-1?j=Dj(b,1):j=mJ(32,d+((a.C[c]&BR)>>10!=1&&(a.C[c]&BR)>>10!=2?-1:(a.C[c]&CR)>>12)));k==2&&(e==-1?j=Dj(b,2):j=mJ(32,e+((a.C[c]&BR)>>10!=1&&(a.C[c]&BR)>>10!=2?-1:(a.C[c]&CR)>>12)));for(l=0;l<2;l++)b.B[l][i]=f==null?a.B[l][c]:f[a.B[l][c]];h=g&&(a.C[c]&512)!=0?64:a.F[c];b.F[i]=h;b.C[i]=a.C[c];b.D[i]=b.I?a.D[c]:0;if(j!=-1){b.C[i]&=-126977;b.C[i]|=j<<12;}++b.p;b.Q=0;return i;}function re(a,b){var c,d,e,f,g,h,i,j,k,l,m;f=XB(_C,[mR,DQ],[7,6],15,[2,32],2);for(g=0;g<2;g++){c=ZB(_C,mR,7,32,0,2);m=0;for(e=0;e<32;e++){if(b[g][e]!=null){k=b[g][e].length;c[e]=ZB(_C,DQ,6,k,15,1);for(h=0;h<k;h++)c[e][h]=a.c[b[g][e][h]];xN(c[e]);++m;}}for(l=m;l>0;l--){j=0;i=null;for(d=0;d<32;d++){if(c[d]!=null){if(i==null||i.length<c[d].length){i=c[d];j=d;}else if(i.length==c[d].length){for(h=i.length-1;h>=0;h--){if(i[h]<c[d][h]){i=c[d];j=d;break;}}}}}f[g][j]=l;c[j]=null;}}return f;}function Ae(a){var b,c,d,e,f,g;a.G=true;d=Ke(a);!!a.J&&Tf(a.J,a.c)&&(d=Ke(a));oe(a,false)&&ue(a);g=true;while(a.N<a.L.d&&g){e=re(a,d);for(b=0;b<a.L.d;b++){Ff(a.b[b],b);Cf(a.b[b],16,a.c[b]);Cf(a.b[b],20,0);!a.V[b]&&a.U[b]!=0&&Df(a.b[b],(a.U[b]<<18)+(e[a.U[b]==1?0:1][a.T[b]]<<8));Df(a.b[b],a.W[b]<<4);}for(c=0;c<a.L.e;c++){Df(a.b[Ei(a.L,0,c)],a.k[c]);Df(a.b[Ei(a.L,1,c)],a.k[c]);}f=ve(a);if(a.N==f)break;a.N=f;g=false;if(!!a.J&&Tf(a.J,a.c)){g=true;d=Ke(a);}if(oe(a,false)){g=true;ue(a);}}}function Td(a,b){var c,d,e,f,g,h,i,j;c=false;if(a.A[b]!=8)return false;if(a.g[b]!=1)return false;if(a.j[b][0]!=1)return false;g=a.f[b][0];if(a.A[g]==6){h=a.g[g];for(d=0;d<h;d++){e=a.f[g][d];if(e==b){continue;}if(a.A[e]!=8){continue;}f=$k(a,g,e);if(a.F[f]==2){c=true;break;}}}else if(a.A[g]==7){a.q[g]==1&&(c=true);}else if(a.A[g]==16){i=a.g[g];j=0;for(d=0;d<i;d++){e=a.f[g][d];if(e==b){continue;}if(a.A[e]!=8){continue;}f=$k(a,g,e);a.F[f]==2&&++j;}j==2&&(c=true);}else Ud(a,b)&&(c=true);return c;}function Oj(a,b,c,d){var e,f,g;if(c==0){a.s[b]&=HR;a.s[b]|=0;}else{if(d>=32)return;if(d==-1){g=-1;for(f=0;f<a.o;f++)f!=b&&c==(a.s[f]&zR)>>19&&g<((a.s[f]&zR)>>19!=1&&(a.s[f]&zR)>>19!=2?-1:(a.s[f]&AR)>>21)&&(g=(a.s[f]&zR)>>19!=1&&(a.s[f]&zR)>>19!=2?-1:(a.s[f]&AR)>>21);for(e=0;e<a.p;e++)c==(a.C[e]&BR)>>10&&g<((a.C[e]&BR)>>10!=1&&(a.C[e]&BR)>>10!=2?-1:(a.C[e]&CR)>>12)&&(g=(a.C[e]&BR)>>10!=1&&(a.C[e]&BR)>>10!=2?-1:(a.C[e]&CR)>>12);d=g+1;if(d>=32)return;}a.s[b]&=HR;a.s[b]|=c<<19|d<<21;}a.Q&=3;}function ek(a,b,c,d){var e,f,g;if(c==0){a.C[b]&=-130049;a.C[b]|=0;}else{if(d>=32)return;if(d==-1){g=-1;for(f=0;f<a.o;f++)c==(a.s[f]&zR)>>19&&g<((a.s[f]&zR)>>19!=1&&(a.s[f]&zR)>>19!=2?-1:(a.s[f]&AR)>>21)&&(g=(a.s[f]&zR)>>19!=1&&(a.s[f]&zR)>>19!=2?-1:(a.s[f]&AR)>>21);for(e=0;e<a.p;e++)e!=b&&c==(a.C[e]&BR)>>10&&g<((a.C[e]&BR)>>10!=1&&(a.C[e]&BR)>>10!=2?-1:(a.C[e]&CR)>>12)&&(g=(a.C[e]&BR)>>10!=1&&(a.C[e]&BR)>>10!=2?-1:(a.C[e]&CR)>>12);d=g+1;if(d>=32)return;}a.C[b]&=-130049;a.C[b]|=c<<10|d<<12;}a.Q&=3;}function tg(a){var b,c,d,e,f,g,h,i,j;while(a.f.a.length>1){g=ZB(ZC,GQ,6,2,15,1);f=ZB(qD,fR,29,2,0,1);b=BM(a.f,0);c=BM(a.f,1);h=(Wg(b),b.i-b.n+1+(Wg(b),b.j-b.o+1));i=(Wg(c),c.i-c.n+1+(Wg(c),c.j-c.o+1));if(h>i){f[0]=b;g[0]=h;f[1]=c;g[1]=i;}else{f[0]=c;g[0]=i;f[1]=b;g[1]=h;}for(e=2;e<a.f.a.length;e++){d=BM(a.f,e);j=(Wg(d),d.i-d.n+1+(Wg(d),d.j-d.o+1));if(g[0]<j){f[1]=f[0];f[0]=d;g[1]=g[0];g[0]=j;}else if(g[1]<j){f[1]=d;g[1]=j;}}Vg(f[0],f[1]);yM(a.f,Eg(a,f[0],f[1],0));GM(a.f,f[0]);GM(a.f,f[1]);}}function Gg(a,b){var c,d,e,f,g,h,i,j,k,l,m,n,o;c=Ei(a.i,0,b);d=Ei(a.i,1,b);g=ZB(_C,DQ,6,a.b,15,1);h=ZB(_C,DQ,6,a.b,15,1);i=ZB(_C,DQ,6,a.b,15,1);j=ZB(_C,DQ,6,a.b,15,1);g[0]=c;g[1]=d;h[1]=b;i[c]=1;i[d]=2;j[0]=-1;j[1]=0;f=1;k=1;while(f<=k){for(l=0;l<a.e[g[f]];l++){e=al(a.i,g[f],l);if(f>1&&e==c){o=new sm(i[g[f]]);h[0]=cl(a.i,g[f],l);m=f;for(n=0;n<o.a.length;n++){o.a[n]=g[m];o.b[n]=h[m];m=j[m];}return o;}if(i[e]==0&&Kl(a.i,e)){g[++k]=e;h[k]=cl(a.i,g[f],l);i[e]=i[g[f]]+1;j[k]=f;}}++f;}return null;}function fg(a,b,c,d){var e,f,g,h,i,j,k,l,m,n,o;for(i=d+1;i<a.j.g.length;i++){if(i!=d&&a.e[b][i]&&a.e[c][i]){g=ZB(_C,DQ,6,2,15,1);g[0]=c;g[1]=b;return g;}}o=ZB(_C,DQ,6,a.b,15,1);k=ZB(_C,DQ,6,a.b,15,1);j=ZB(_C,DQ,6,a.b,15,1);f=0;l=0;j[0]=b;k[b]=1;while(f<=l){for(m=0;m<a.d[j[f]].length;m++){e=a.d[j[f]][m];if(e==c){if(f==0)continue;h=k[j[f]]+1;g=ZB(_C,DQ,6,h,15,1);g[0]=e;g[1]=j[f];for(n=2;n<h;n++)g[n]=o[g[n-1]];return g;}if(k[e]==0&&a.c[e]!=-3){k[e]=k[j[f]]+1;j[++l]=e;o[e]=j[f];}}++f;}return null;}function wo(a){var b,c,d,e,f,g;f='<svg id="'+(a.g!=null?a.g:'mol'+mo)+_Q+'xmlns="http://www.w3.org/2000/svg" version="1.1" '+'width="'+a.k+'px" '+'height="'+a.f+'px" '+'viewBox="0 0 '+a.k+' '+a.f+'">\n';g='<style> #'+(a.g!=null?a.g:'mol'+mo)+' {pointer-events:none; } '+' #'+(a.g!=null?a.g:'mol'+mo)+' .event '+' { pointer-events:all;} '+' <\/style>\n';f+='\t';f+=g;for(e=new dN(a.b);e.a<e.c.a.length;){d=cN(e);xo(a,d);}for(c=new dN(a.a);c.a<c.c.a.length;){b=cN(c);xo(a,b);}return f+a.c.a+'<\/svg>';}function bp(a,b,c){var d,e,f,g,h,i,j,k,l,m,n,o,p,q;p=ZB(PD,iQ,24,c,0,1);g=ZB(_C,DQ,6,c,15,1);j=ZB(_C,DQ,6,c,15,1);f=ZB(_C,DQ,6,a.o,15,1);for(e=0;e<a.o;e++)b[e]!=-1&&(f[e]=g[b[e]]++);for(i=0;i<a.p;i++){n=b[a.B[0][i]];o=b[a.B[1][i]];n==o&&n!=-1&&++j[n];}for(q=0;q<c;q++){p[q]=new lp(g[q],j[q]);Wo(a,p[q]);}for(d=0;d<a.o;d++)b[d]!=-1&&Vh(a,p[b[d]],d,0,0);for(h=0;h<a.p;h++){n=b[a.B[0][h]];o=b[a.B[1][h]];n==o&&n!=-1&&Wh(a,p[n],h,0,0,f,false);}for(l=0,m=p.length;l<m;++l){k=p[l];Dj(k,1);Dj(k,2);}return p;}function jA(a,b){var c,d,e,f,g,h,i,j,k;c=UG(a);h=!a?1:(hA[0]*(a.c>>16&255)+hA[1]*(a.c>>8&255)+hA[2]*(a.c&255))/255;if(h==0)return new VG(h,h,h,c[3]);d=b/(!a?1:(hA[0]*(a.c>>16&255)+hA[1]*(a.c>>8&255)+hA[2]*(a.c&255))/255);k=0;j=0;for(f=0;f<3;f++){c[f]*=d;if(c[f]<1){j+=hA[f];}else{k+=(c[f]-1)*hA[f];c[f]=1;}}if(k!=0){i=0;for(g=0;g<3;g++){if(c[g]<1){c[g]+=k/j;if(c[g]>1){i+=(c[g]-1)*hA[g];c[g]=1;}}}if(i!=0){for(e=0;e<3;e++){if(c[e]<1){c[e]+=i/hA[e];c[e]>1&&(c[e]=1);}}}}return new VG(c[0],c[1],c[2],c[3]);}function ud(a,b){var c,d,e,f,g,h,i,j,k,l,m,n,o;k=ZB(aG,HQ,6,16,16,1);l=ZB(aG,HQ,6,16,16,1);c=ZB(ZC,GQ,6,16,15,1);f=ZB(ZC,GQ,6,2,15,1);d=0;for(j=0;j<2;j++){e=Ei(a.G,j,b);for(m=0;m<bl(a.G,e);m++){h=cl(a.G,e,m);if(h==b)continue;if(d==4)return 0;k[d]=Fl(a.G,h);l[d]=Ll(a.G,h);c[d++]=Di(a.G,e,al(a.G,e,m));}}f[0]=Di(a.G,Ei(a.G,0,b),Ei(a.G,1,b));if(f[0]<0){f[1]=f[0]+MQ;g=false;}else{f[1]=f[0];f[0]=f[1]-MQ;g=true;}n=0;for(i=0;i<d;i++){k[i]?o=20:l[i]?o=17:o=16;c[i]>f[0]&&c[i]<f[1]?n-=o:n+=o;}return g?-n:n;}function Rn(a,b){var c,d,e,f,g,h,i,j,k,l;f=null;i=null;g=null;Xo(a.d,a.G);a.i=ZB(_C,DQ,6,a.d.d,15,1);for(d=0;d<a.d.d;d++)a.i[d]=bl(a.d,d);if(a.j!=0){j=new lp(a.d.o,a.d.p);l=ZB(aG,HQ,6,a.d.o,16,1);for(e=0;e<a.d.o;e++)l[e]=!a.u[e];Jk(a.d,j,l,true,null);Xo(j,a.G);Sn(a,j,b);f=a.e;i=a.g;g=a.f;k=0;for(c=0;c<a.d.d;c++)a.u[c]||(a.i[c]=bl(j,k++));}Sn(a,a.d,b);if(a.j!=0){k=0;for(c=0;c<a.d.o;c++){if(!a.u[c]){a.e[c]=f[k];a.f[c]=g[k++];}}k=0;for(h=0;h<a.d.p;h++){!a.u[Ei(a.d,0,h)]&&!a.u[Ei(a.d,1,h)]&&(a.g[h]=i[k++]);}}}function fq(a,b){var c,d,e,f,g,h,i,j,k;c=new oq();if(!eq){yM(c.a,new mq(NT,2));return-999;}yM(c.a,new mq('Found sub-structure fragments and their contributions:',2));yM(c.a,new mq('(yellow atoms carry at least one more substituent)',2));j=0;i=0;f=0;k=new Wn();e=new kp();for(g=0;g<dq.a.a.length;g++){mm(new om(false),e,iq(dq,g));Pn(k,b);On(k,e);if(Fn(k,k.b)>0){h=jq(dq,g);if(h<-1)j+=h;else{i+=h;++f;}nq(c,iq(dq,g),1);yM(c.a,new mq(''+h,3));}}if(f==0)return-1;d=j+i/$wnd.Math.sqrt(f);d=d+0.0625*(f-40);a.a=c;return d;}function Ee(b,c){var d,e,f,g,h,i,j;if(b.W[c]==1||b.W[c]==2){i=false;if(Tk(b.L,c)==2){try{for(h=0;h<2;h++){d=al(b.L,c,h);if(bl(b.L,d)==3){f=ZB(_C,DQ,6,2,15,1);g=0;for(j=0;j<bl(b.L,d);j++)dl(b.L,d,j)==1&&(f[g++]=al(b.L,d,j));b.c[f[0]]>b.c[f[1]]^Fe(b,d,f[0],f[1])&&(i=!i);}}}catch(a){a=cG(a);if(PC(a,12)){b.R[c]=3;return;}else throw dG(a);}}else{try{e=He(b,c);}catch(a){a=cG(a);if(PC(a,12)){b.R[c]=3;return;}else throw dG(a);}for(h=1;h<e.length;h++)for(j=0;j<h;j++)b.c[e[h]]<b.c[e[j]]&&(i=!i);}b.W[c]==1^i?b.R[c]=1:b.R[c]=2;}}function Ve(a){var b,c,d,e,f,g,h,i,j,k,l,m,n;h=null;n=tl(a.L);for(k=0;k<n.g.a.length;k++){if(n.e[k]){e=0;l=BM(n.g,k);for(c=0,d=l.length;c<d;++c){b=l[c];cf(a,b)&&++e;}if(e!=0){m=BM(n.i,k);h==null&&(h=ZB(aG,HQ,6,a.L.e,16,1));if(e==l.length){i=-1;j=oQ;for(f=0;f<l.length;f++){if(j>a.t[m[f]]){j=a.t[m[f]];i=f;}}while(e>0){h[m[i]]=true;i=pf(i+2,l.length);e-=2;}}else{g=0;while(cf(a,l[g]))++g;while(!cf(a,l[g]))g=pf(g+1,l.length);while(e>0){h[m[g]]=true;g=pf(g+2,l.length);e-=2;while(!cf(a,l[g]))g=pf(g+1,l.length);}}}}}return h;}function Uf(a){var b,c,d,e,f,g,h,i;if(a.g!=null){g=new ng(a);a.b=new LM();for(e=0;e<a.g.length;e++){d=hg(g,e);if(d==0){dg(g,e);h=Lf(a,e,2);b=Lf(a,e,1);c=Kf(a,e);if(h==1&&b==1&&!c){Wf(a,e,g.a+g.f++);yM(a.b,new zh(e,1,-1,-1));}if(h>0){if(c){Vf(a,e,g.i+g.g++,2);++h;}yM(a.b,new zh(e,1,-1,-1));}else if(b>0){c&&Vf(a,e,g.a+g.f++,1);yM(a.b,new zh(e,1,-1,-1));}else if(c){Vf(a,e,g.a+g.f++,1);yM(a.b,new zh(e,1,-1,-1));}}else if(d==1){if(Kf(a,e)){f=gg(g,e);i=ig(g,e);yM(a.b,new zh(e,2,f,i));}else{dg(g,e);yM(a.b,new zh(e,1,-1,-1));}}}}}function Kk(a,b,c,d,e){var f,g,h,i,j;d&&Xo(a,3);b.t=null;a.I&&lk(b,true);e==null&&(e=ZB(_C,DQ,6,a.o,15,1));b.o=0;for(f=0;f<a.o;f++){e[f]=-1;for(j=0;j<a.g[f];j++){if(c[a.i[f][j]]){e[f]=Vh(a,b,f,0,0);break;}}}b.p=0;for(i=0;i<a.p;i++)if(c[i]){Wh(a,b,i,0,0,e,d);}else{g=a.B[0][i];h=a.B[1][i];if(e[g]==-1^e[h]==-1){if(a.q[g]!=0&&a.q[h]!=0&&a.q[g]<0^a.q[h]<0){e[g]!=-1&&(b.q[e[g]]+=a.q[g]<0?1:-1);e[h]!=-1&&(b.q[e[h]]+=a.q[h]<0?1:-1);}}}Yh(a,b);!!a.b&&(b.Q=0);b.Q=0;Dj(b,1);Dj(b,2);b.o!=a.o&&lk(b,true);d&&Id(new Pd(b));return e;}function Cq(b){var c,d,e,f,g,h,i;e=new oq();yM(e.a,new mq('Solubility values are estimated applying an atom-type based increment system.',2));yM(e.a,new mq(LT,2));yM(e.a,new mq(MT,2));yM(e.a,new mq('Base value = -0.530',2));d=ZB(_C,DQ,6,zq.length,15,1);if(b){for(c=0;c<b.d;c++){i=-1;try{i=_d(b,c,2144);}catch(a){a=cG(a);if(!PC(a,12))throw dG(a);}for(h=0;h<zq.length;h++){if(iG(yq[h],i)){++d[h];break;}}}}f=new mK('#0.000');for(g=0;g<zq.length;g++)d[g]!=0&&nq(e,''+d[g]+' * '+lK(f,zq[g])+'   AtomType: 0x'+fJ(yq[g]),2);return e;}function rf(a,b){var c;if(a.o>AQ)throw dG(new NI('Cannot canonize a molecule having more than 65535 atoms'));if(a.p>AQ)throw dG(new NI('Cannot canonize a molecule having more than 65535 bonds'));this.L=a;this.K=b;Xo(this.L,3);ne(this);for(c=0;c<this.L.d;c++){if(zi(this.L,c)!=0){this._=true;break;}}this.W=ZB(XC,pR,6,this.L.d,15,1);this.X=ZB(aG,HQ,6,this.L.d,16,1);this.$=ZB(aG,HQ,6,this.L.d,16,1);this.k=ZB(XC,pR,6,this.L.e,15,1);this.o=ZB(aG,HQ,6,this.L.e,16,1);this.n=ZB(aG,HQ,6,this.L.e,16,1);se(this);xe(this);we(this);}function zd(a,b){if(b==-10){a.w=-999;b=a.J;}if(b==a.w)return;if(a.w==-8&&b!=-9)return;b==-8&&(a.I=a.w);b==-9&&(b=a.I);a.w=b;switch(b){case 0:to(a,(TG(),QG));break;case-6:to(a,a.A);break;case-4:to(a,a.H);break;case-2:to(a,a.r);break;case-3:to(a,a.s);break;case-7:to(a,a.C);break;case-8:to(a,a.D);break;case 64:to(a,Ac);break;case 128:to(a,Gc);break;case 256:to(a,Ec);break;case 192:to(a,Dc);break;case 320:to(a,Fc);break;case 384:to(a,Bc);break;case 448:to(a,Cc);break;case 1:to(a,(TG(),RG));break;default:to(a,(TG(),QG));}}function dm(a){var b,c,d,e,f,g,h;if(!a.I)return false;for(c=0;c<a.o;c++)ol(a,c)>=Ui(a,c)+Si(a,c)&&(a.w[c]&=-6145);e=false;for(b=0;b<a.d;b++){f=el(a,b);if(!a.P&&f>0){if((a.w[b]&YQ)==0){if(Ui(a,b)+Si(a,b)-ol(a,b)==0)a.w[b]|=YQ;else{h=0;(a.w[b]&128)==128&&++h;(a.w[b]&1920)==384&&++h;a.w[b]&=-1921;Ui(a,b)+Si(a,b)-ol(a,b)<=h?a.w[b]|=YQ:h==0?a.w[b]|=128:a.w[b]|=384;}}for(g=a.g[b];g<a.c[b];g++){d=a.i[b][g];if(a.F[d]==1){a.A[a.f[b][g]]=-1;a.F[d]=128;e=true;}}}(a.w[b]&2)!=0&&(a.w[b]&=-9);a.q[b]!=0&&(a.s[b]&=-234881025);}e&&Uh(a);return e;}function ip(a){var b,c,d,e,f,g;bm(a);Xo(a,15);for(d=0;d<a.d;d++){if(((a.s[d]&zR)>>19==1||(a.s[d]&zR)>>19==2)&&((a.s[d]&FR)==0||(a.s[d]&3)==3))throw dG(new FA(SR));if((a.s[d]&qR)!=0)throw dG(new FA('Over- or under-specified stereofeature or more than one racemic type bond'));if(((a.s[d]&3)==1||(a.s[d]&3)==2)&&a.k[d]==0){b=ZB(ZC,GQ,6,a.g[d],15,1);for(f=0;f<a.g[d];f++)b[f]=Di(a,d,a.f[d][f]);for(e=1;e<a.g[d];e++)if(!tj(a,a.i[d][e],d))for(g=0;g<e;g++)if(!tj(a,a.i[d][g],d)){c=$wnd.Math.abs(Bk(b[e],b[g]));if(c<0.08||c>JR)throw dG(new FA(TR));}}}}function kA(a,b){iA();var c,d,e,f,g,h,i,j,k,l,m;c=!b?1:(hA[0]*(b.c>>16&255)+hA[1]*(b.c>>8&255)+hA[2]*(b.c&255))/255;g=!a?1:(hA[0]*(a.c>>16&255)+hA[1]*(a.c>>8&255)+hA[2]*(a.c&255))/255;e=$wnd.Math.abs(c-g);if(e>FQ)return a;h=ZB($C,gR,6,3,15,1);ZG(b.c>>16&255,b.c>>8&255,b.c&255,h);i=ZB($C,gR,6,3,15,1);ZG(a.c>>16&255,a.c>>8&255,a.c&255,i);j=$wnd.Math.abs(i[0]-h[0]);j>0.5&&(j=1-j);m=1-$wnd.Math.max(i[1],h[1]);d=$wnd.Math.abs(g+c-1);k=$wnd.Math.cos(MQ*j*3);l=FQ*$wnd.Math.max(m,$wnd.Math.max(d,k));if(e>l)return a;f=g>c?g+l>1:g-l>0;return jA(a,f?c-l:c+l);}function Ig(a,b){var c,d,e,f;!a.j&&(a.j=new RN());a.i=b;Xo(a.i,3);if((a.g&1)==0){a.b=a.i.o;a.d=a.i.p;a.e=ZB(_C,DQ,6,a.b,15,1);for(c=0;c<a.b;c++)a.e[c]=Qk(a.i,c);}else{a.b=a.i.d;a.d=a.i.e;a.e=ZB(_C,DQ,6,a.b,15,1);for(c=0;c<a.b;c++)a.e[c]=bl(a.i,c);}a.f=new LM();a.a=ZB(aG,HQ,6,a.b,16,1);a.c=ZB(aG,HQ,6,a.d,16,1);(a.g&6)!=0&&Mg(a);Pg(a);Kg(a);Lg(a);Kg(a);Og(a);xg(a);Rg(a);Qg(a);tg(a);for(e=0;e<a.f.a.length;e++){d=BM(a.f,e);for(f=0;f<d.a.length;f++){Zj(a.i,d.a[f],d.c[f]);$j(a.i,d.a[f],d.d[f]);_j(a.i,d.a[f],0);}}if((a.g&1)!=0){if(a.b!=0){Fj(a.i,a.b);Gj(a.i,a.d);}}}function uC(a,b){var c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,A,B,C,D,F,G;c=a.l&8191;d=a.l>>13|(a.m&15)<<9;e=a.m>>4&8191;f=a.m>>17|(a.h&255)<<5;g=(a.h&1048320)>>8;h=b.l&8191;i=b.l>>13|(b.m&15)<<9;j=b.m>>4&8191;k=b.m>>17|(b.h&255)<<5;l=(b.h&1048320)>>8;B=c*h;C=d*h;D=e*h;F=f*h;G=g*h;if(i!=0){C+=c*i;D+=d*i;F+=e*i;G+=f*i;}if(j!=0){D+=c*j;F+=d*j;G+=e*j;}if(k!=0){F+=c*k;G+=d*k;}l!=0&&(G+=c*l);n=B&OR;o=(C&511)<<13;m=n+o;q=B>>22;r=C>>9;s=(D&262143)<<4;t=(F&31)<<17;p=q+r+s+t;v=D>>18;w=F>>5;A=(G&4095)<<8;u=v+w+A;p+=m>>22;m&=OR;u+=p>>22;p&=OR;u&=VT;return fC(m,p,u);}function Te(a){var b,c,d;for(b=0;b<a.L.d;b++){$i(a.L,b)^a.W[b]==3&&qk(a.L,b);(oi(a.L,b)==1||oi(a.L,b)==2)&&(!a.H[b]||a.W[b]==3)&&qk(a.L,b);}for(d=0;d<a.L.p;d++)sj(a.L,d)&&!jf(a,d)&&qk(a.L,Ei(a.L,0,d));for(c=0;c<a.L.e;c++){if(Mi(a.L,c)==2){if(ij(a.L,c)&&(a.k[c]==1||a.k[c]==2)){a.k[c]=3;jk(a.L,c,26);}if(a.k[c]==3&&!a.n[c]){if(Pi(a.L,c)!=26){qk(a.L,Ei(a.L,0,c));qk(a.L,Ei(a.L,1,c));}}}if(Pi(a.L,c)==1&&a.k[c]==3){qk(a.L,Ei(a.L,0,c));qk(a.L,Ei(a.L,1,c));}if((Ji(a.L,c)==1||Ji(a.L,c)==2)&&(Pi(a.L,c)!=1||a.k[c]!=1&&a.k[c]!=2)){qk(a.L,Ei(a.L,0,c));qk(a.L,Ei(a.L,1,c));}}}function gC(a,b,c){var d,e,f,g,h,i;if(b.l==0&&b.m==0&&b.h==0){throw dG(new HH());}if(a.l==0&&a.m==0&&a.h==0){c&&(cC=fC(0,0,0));return fC(0,0,0);}if(b.h==tQ&&b.m==0&&b.l==0){return hC(a,c);}i=false;if(b.h>>19!=0){b=vC(b);i=true;}g=nC(b);f=false;e=false;d=false;if(a.h==tQ&&a.m==0&&a.l==0){e=true;f=true;if(g==-1){a=eC((JC(),FC));d=true;i=!i;}else{h=yC(a,g);i&&lC(h);c&&(cC=fC(0,0,0));return h;}}else if(a.h>>19!=0){f=true;a=vC(a);d=true;i=!i;}if(g!=-1){return iC(a,g,i,f,c);}if(sC(a,b)<0){c&&(f?cC=vC(a):cC=fC(a.l,a.m,a.h));return fC(0,0,0);}return jC(d?a:fC(a.l,a.m,a.h),b,i,f,e,c);}function AI(a){var b,c,d,e,f,g;if(isNaN(a)){return{l:0,m:0,h:524160};}g=false;if(a==0){return 1/a==-Infinity?{l:0,m:0,h:tQ}:0;}if(a<0){g=true;a=-a;}if(!isNaN(a)&&!isFinite(a)){return g?{l:0,m:0,h:1048320}:{l:0,m:0,h:524032};}c=0;if(a<1){b=512;for(d=0;d<10;++d,b>>=1){if(a<(DI(),BI)[d]&&c-b>=-1023){a*=CI[d];c-=b;}}if(a<1&&c-1>=-1023){a*=2;--c;}}else if(a>=2){b=512;for(d=0;d<10;++d,b>>=1){if(a>=(DI(),CI)[d]){a*=BI[d];c+=b;}}}c>-1023?a-=1:a*=0.5;e=jG(a*1048576);a-=sG(e)*9.5367431640625E-7;f=jG(a*4503599627370496);e=oG(e,c+1023<<20);g&&(e=oG(e,2147483648));return oG(pG(e,32),f);}function bd(a,b,c,d,e){var f,g,h,i,j,k,l,m,n,o;m=false;e.a=0;e.b=0;d>0?f=JQ:f=KQ;o=Di(a.G,b,c);for(k=0;k<bl(a.G,b);k++){g=cl(a.G,b,k);h=o;Ei(a.G,0,g)==b?l=Ei(a.G,1,g):l=Ei(a.G,0,g);if(l==c)continue;n=Di(a.G,b,l);o<n&&(h+=LQ);i=h-n;if(d>0){i<MQ&&(m=true);i>JQ&&(i=JQ);i<0.523598776&&(i=0.523598776);if(i<=f){f=i;j=a.M*$wnd.Math.tan(i-NQ)/2;e.a=-(j*$wnd.Math.sin(h));e.b=-(j*$wnd.Math.cos(h));}}else{i>=MQ&&(m=true);i<KQ&&(i=KQ);i>5.759586531&&(i=5.759586531);if(i>=f){f=i;j=a.M*$wnd.Math.tan(4.712388981-i)/2;e.a=-(j*$wnd.Math.sin(h));e.b=-(j*$wnd.Math.cos(h));}}}return m;}function Mo(a){var b,c,d,e,f,g,h,i,j,k,l,m,n;Xo(a.b,3);l=false;m=ZB(_C,DQ,6,2,15,1);n=ZB(_C,DQ,6,2,15,1);k=ZB(_C,DQ,6,2,15,1);for(d=0;d<a.b.e;d++){if(!Ol(a.b,d)&&Pi(a.b,d)==2){for(g=0;g<2;g++){m[g]=-1;k[g]=-1;b=Ei(a.b,g,d);for(j=0;j<bl(a.b,b);j++){e=cl(a.b,b,j);if(e!=d){if(Pi(a.b,e)==17||Pi(a.b,e)==9){m[g]=al(a.b,b,j);n[g]=e;}else{k[g]=al(a.b,b,j);}}}if(m[g]==-1)break;}if(m[0]!=-1&&m[1]!=-1){i=Pi(a.b,n[0])!=Pi(a.b,n[1]);h=false;for(f=0;f<2;f++){k[f]!=-1&&k[f]<m[f]&&(h=!h);}gk(a.b,d,i^h?2:1,false);l=true;}}}for(c=0;c<a.b.e;c++)(Pi(a.b,c)==17||Pi(a.b,c)==9)&&jk(a.b,c,1);return l;}function ke(a,b,c){var d,e,f,g,h,i;d=ZB(_C,DQ,6,4,15,1);for(h=0;h<Qk(a.L,b);h++)d[h]=al(a.L,b,c[h]);Qk(a.L,b)==3&&(d[3]=b);e=XB(ZC,[iQ,GQ],[15,6],15,[3,3],2);for(g=0;g<3;g++){e[g][0]=xi(a.L,d[g+1])-xi(a.L,d[0]);e[g][1]=yi(a.L,d[g+1])-yi(a.L,d[0]);e[g][2]=zi(a.L,d[g+1])-zi(a.L,d[0]);}i=ZB(ZC,GQ,6,3,15,1);i[0]=e[0][1]*e[1][2]-e[0][2]*e[1][1];i[1]=e[0][2]*e[1][0]-e[0][0]*e[1][2];i[2]=e[0][0]*e[1][1]-e[0][1]*e[1][0];f=(e[2][0]*i[0]+e[2][1]*i[1]+e[2][2]*i[2])/($wnd.Math.sqrt(e[2][0]*e[2][0]+e[2][1]*e[2][1]+e[2][2]*e[2][2])*$wnd.Math.sqrt(i[0]*i[0]+i[1]*i[1]+i[2]*i[2]));return f>0?1:2;}function gn(a,b,c,d){var e,f,g,h,i,j,k,l,m,n,o,p;e=ZB(_C,mR,7,a.g.a.length,0,2);for(i=0;i<a.g.a.length;i++){e[i]=ZB(_C,DQ,6,BM(a.g,i).length,15,1);for(j=0;j<BM(a.g,i).length;j++)e[i][j]=-1;}o=ZB(_C,DQ,6,a.f.e,15,1);for(m=0;m<a.i.a.length;m++){n=BM(a.i,m);if(n.length>=5&&n.length<=7){for(h=0;h<n.length;h++){g=n[h];if(bl(a.f,Ei(a.f,0,g))==3&&bl(a.f,Ei(a.f,1,g))==3){if(o[g]>0){e[o[g]>>>16][o[g]&32767]=m;e[m][h]=o[g]>>>16;}else{o[g]=(m<<16)+32768+h;}}}}}f=ZB(aG,HQ,6,a.g.a.length,16,1);p=0;k=-1;while(p>k){k=p;for(l=0;l<a.g.a.length;l++){if(!f[l]){if(fn(a,l,e,f,b,c,d)){f[l]=true;++p;}}}}}function HP(a,b,c){if(c<128){a[b]=(c&127)<<24>>24;return 1;}else if(c<YQ){a[b++]=(c>>6&31|192)<<24>>24;a[b]=(c&63|128)<<24>>24;return 2;}else if(c<zQ){a[b++]=(c>>12&15|224)<<24>>24;a[b++]=(c>>6&63|128)<<24>>24;a[b]=(c&63|128)<<24>>24;return 3;}else if(c<UQ){a[b++]=(c>>18&7|240)<<24>>24;a[b++]=(c>>12&63|128)<<24>>24;a[b++]=(c>>6&63|128)<<24>>24;a[b]=(c&63|128)<<24>>24;return 4;}else if(c<ER){a[b++]=(c>>24&3|248)<<24>>24;a[b++]=(c>>18&63|128)<<24>>24;a[b++]=(c>>12&63|128)<<24>>24;a[b++]=(c>>6&63|128)<<24>>24;a[b]=(c&63|128)<<24>>24;return 5;}throw dG(new NI('Character out of range: '+c));}function Rf(a,b,c){var d,e,f,g,h;if(b==c)return false;if(a.a[b]!==a.a[c])return false;if(a.o[b]!=0){if(a.o[b]==3||a.o[c]==3)return false;if(a.p[b]^a.o[b]!==a.o[c])return false;if(a.k[b]!==a.k[c]||a.j[b]!==a.j[c])return false;}d=$k(a.i,b,c);if(d!=-1){if(Mi(a.i,d)==1&&a.c[d]!=0)return false;if(Mi(a.i,d)==2&&a.c[d]==2)return false;}if(Tk(a.i,b)==1&&!El(a.i,b)){e=-1;for(h=0;h<bl(a.i,b);h++){if(al(a.i,b,h)!=c&&dl(a.i,b,h)==2){e=cl(a.i,b,h);break;}}f=-1;for(g=0;g<bl(a.i,c);g++){if(al(a.i,c,g)!=b&&dl(a.i,c,g)==2){f=cl(a.i,c,g);break;}}if(e!=-1&&a.c[e]!=0&&a.d[e]^a.c[e]===a.c[f])return false;}return true;}function Ph(a,b,c,d){var e,f,g,h,i,j,k,l,m,n,o;i=ZB(_C,DQ,6,2,15,1);h=ZB(ZC,GQ,6,2,15,1);i[0]=a.B[0][b];i[1]=a.B[1][b];if(ol(a,i[0])>3)return false;if(ol(a,i[1])>3)return false;f=0;e=ZB(ZC,GQ,6,4,15,1);for(l=0;l<a.p;l++){if(l==b)continue;for(m=0;m<2;m++){for(n=0;n<2;n++){if(a.B[m][l]===i[n]){if(f==4){f=5;break;}e[f++]=Di(a,i[n],a.B[1-m][l]);}}if(f==5)break;}if(f==5)break;}if(f==5)return false;h[0]=Di(a,i[0],i[1]);if(h[0]<0){h[1]=h[0]+MQ;g=0;}else{h[1]=h[0];h[0]=h[1]-MQ;g=1;}o=0;for(k=0;k<f;k++){e[k]>h[0]&&e[k]<h[1]?--o:++o;}g=o>0?1-g:g;j=MQ*(c-2)/c;wj(a,i[g],c-1,i[1-g],d,h[o>0?0:1]+MQ-j,MQ-j);a.Q=0;return true;}function Po(a,b){var c,d,e,f,g,h,i,j,k,l,m,n,o,p;if(a.b)return 3;a.c!=-1&&(a.c=b[a.c]);for(g=0;g<a.g;g++)a.f[g]!=oQ&&(a.f[g]=b[a.f[g]]);if(a.c==-1&&a.d==0){n=oQ;m=-1;for(h=0;h<a.g;h++){if(n>a.j[h]){n=a.j[h];m=h;}}a.c=a.f[m];for(i=m+1;i<a.g;i++){a.f[i-1]=a.f[i];a.j[i-1]=a.j[i];a.i[i-1]=a.i[i];}--a.g;}p=(a.c==-1?0:1)+a.d+a.g;if(p>4||p<3)return 3;c=a.c==-1&&a.d==1||a.c!=-1&&Ml(a.k.b,a.c);e=-1;for(j=0;j<a.g;j++){if(a.i[j]){if(e!=-1||c)return 3;e=j;}}l=false;if(e!=-1)for(k=0;k<a.g;k++)!a.i[k]&&a.f[e]<a.f[k]&&(l=!l);d=false;if(a.c!=-1&&!c)for(f=0;f<a.g;f++)a.c<a.f[f]&&(d=!d);o=a.e^Qo(a.f,a.j,a.g)^d^l?2:1;return o;}function gf(a,b){var c,d,e,f,g,h,i,j,k,l,m,n,o;i=ZB(_C,DQ,6,32,15,1);j=0;for(d=0;d<a.L.d;d++){if((a.S[d]==1||a.S[d]==2)&&a.U[d]==b){h=a.T[d];if(i[h]<a.c[d]){i[h]==0&&++j;i[h]=a.c[d];}}}for(f=0;f<a.L.e;f++){if((a.g[f]==1||a.g[f]==2)&&a.j[f]==b&&Pi(a.L,f)==1){h=a.i[f];o=lJ(a.c[Ei(a.L,0,f)],a.c[Ei(a.L,1,f)]);if(i[h]<o){i[h]==0&&++j;i[h]=o;}}}g=ZB(XC,pR,6,32,15,1);for(k=0;k<j;k++){m=-1;n=0;for(l=0;l<32;l++){if(n<i[l]){n=i[l];m=l;}}i[m]=0;g[m]=k<<24>>24;}for(c=0;c<a.L.d;c++)(a.S[c]==1||a.S[c]==2)&&a.U[c]==b&&(a.T[c]=g[a.T[c]]);for(e=0;e<a.L.e;e++)(a.g[e]==1||a.g[e]==2)&&a.j[e]==b&&Pi(a.L,e)==1&&(a.i[e]=g[a.i[e]]);}function Al(a){var b,c,d,e,f,g,h,i,j,k,l,m,n;j=ZB(aG,HQ,6,a.o,16,1);for(f=0;f<a.p;f++)for(h=0;h<2;h++)Ml(a,a.B[h][f])&&!Ml(a,a.B[1-h][f])&&(j[a.B[h][f]]=true);k=a.o;do--k;while(k>=0&&j[k]);for(b=0;b<k;b++){if(j[b]){uk(a,b,k);m=j[b];j[b]=j[k];j[k]=m;do--k;while(j[k]);}}a.d=k+1;if(a.o==a.d){a.e=a.p;return;}i=ZB(aG,HQ,6,a.p,16,1);for(g=0;g<a.p;g++){c=a.B[0][g];d=a.B[1][g];(j[c]||j[d])&&(i[g]=true);}l=a.p;do--l;while(l>=0&&i[l]);for(e=0;e<l;e++){if(i[e]){n=a.B[0][e];a.B[0][e]=a.B[0][l];a.B[0][l]=n;n=a.B[1][e];a.B[1][e]=a.B[1][l];a.B[1][l]=n;n=a.F[e];a.F[e]=a.F[l];a.F[l]=n;i[e]=false;do--l;while(i[l]);}}a.e=l+1;}function we(a){var b,c,d,e,f,g,h;a.P=ZB(aG,HQ,6,a.L.d,16,1);a.O=ZB(aG,HQ,6,a.L.e,16,1);if((a.K&6)!=0){for(b=0;b<a.L.d;b++){Ff(a.b[b],b);Cf(a.b[b],28,a.c[b]<<12);}}if(a.N<a.L.d){f=0;for(b=0;b<a.L.d;b++)ie(a,b,true)&&++f;for(e=0;e<a.L.e;e++)ee(a,e,true)&&++f;}(a.K&6)!=0&&(a.N=ve(a));if((a.K&1)!=0){a.d=ZB(_C,DQ,6,a.L.d,15,1);for(b=0;b<a.L.d;b++)a.d[b]=a.c[b];}while(a.N<a.L.d){for(c=0;c<a.L.d;c++){Ff(a.b[c],c);Cf(a.b[c],17,2*a.c[c]);}h=ZB(_C,DQ,6,a.N+1,15,1);for(d=0;d<a.L.d;d++)++h[a.c[d]];g=1;while(h[g]==1)++g;for(b=0;b<a.L.d;b++){if(a.c[b]==g){Df(a.b[b],1);break;}}a.N=ve(a);ue(a);!!a.J&&Tf(a.J,a.c);}ue(a);pe(a);Te(a);}function Mf(a){var b,c,d,e,f,g,h,i,j,k,l,m;k=new LM();for(l=0;l<a.i.d;l++){if(Tk(a.i,l)<2||bl(a.i,l)>2){for(g=1;g<bl(a.i,l);g++){b=al(a.i,l,g);for(j=0;j<g;j++){c=al(a.i,l,j);Rf(a,b,c)&&Jf(a,Yf(a,b,c),k);}}}}for(m=0;m<a.i.e;m++){if(a.c[m]!=0){if(Mi(a.i,m)!=2||a.c[m]!=2)continue;}b=Ei(a.i,0,m);c=Ei(a.i,1,m);Rf(a,b,c)&&Jf(a,Yf(a,b,c),k);}for(h=k.a.length-1;h>=0;h--){d=(NP(h,k.a.length),k.a[h]);e=false;for(j=0;j<d.length;j++){if(a.f[d[j]]){e=true;break;}}e||GM(k,d);}a.g=KM(k,XB(_C,[mR,DQ],[7,6],15,[0,0],2));zN(a.g,new ag());a.e=ZB(aG,HQ,6,a.i.d,16,1);for(f=0;f<a.g.length;f++)for(i=0;i<a.g[f].length;i++)a.e[a.g[f][i]]=true;}function de(a,b,c){var d,e,f,g,h,i,j;if(!Gl(a.L,b))return false;d=Ei(a.L,0,b);e=Ei(a.L,1,b);g=new Bh(a.L,a.c,d,e);if(g.f&&!c)return false;h=new Bh(a.L,a.c,e,d);if(h.f&&!c)return false;if(g.f&&h.f)return false;if(c){g.f&&(a.O[b]=bf(a,e));h.f&&(a.O[b]=bf(a,d));}i=Ah(g);j=Ah(h);if(i==-1||j==-1||(i+j&1)==0){c||(a.k[b]=3);return true;}f=0;switch(i+j){case 3:case 7:f=1;break;case 5:f=2;}if(c){if(a.Q&&(a.K&2)!=0||!a.Q&&(a.K&4)!=0){if(g.f){if(f==2){Df(a.b[g.b],4);Df(a.b[g.d],1);}else{Df(a.b[g.b],1);Df(a.b[g.d],4);}}if(h.f){if(f==2){Df(a.b[h.b],4);Df(a.b[h.d],1);}else{Df(a.b[h.b],1);Df(a.b[h.d],4);}}}}else{a.k[b]=f;}return true;}function Sg(a,b,c){var d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t;e=ZB(ZC,GQ,6,a.e[c]+1,15,1);g=ZB(_C,DQ,6,a.e[c]+1,15,1);h=ZB(_C,DQ,6,a.e[c]+1,15,1);q=_g(b,c);f=0;for(j=0;j<a.e[c];j++){g[f]=al(a.i,c,j);h[f]=cl(a.i,c,j);l=_g(b,g[f]);l!=-1&&(e[f++]=rm(b.c[q],b.d[q],b.c[l],b.d[l]));}if(f==1)return e[0]+MQ;for(k=f-1;k>0;k--){for(m=0;m<k;m++){if(e[m]>e[m+1]){r=e[m];e[m]=e[m+1];e[m+1]=r;s=g[m];g[m]=g[m+1];g[m+1]=s;t=h[m];h[m]=h[m+1];h[m+1]=t;}}}e[f]=e[0]+LQ;g[f]=g[0];h[f]=h[0];n=-100;o=0;for(i=0;i<f;i++){d=e[i+1]-e[i];if(f>2&&Ll(a.i,h[i])&&Ll(a.i,h[i+1])){p=Hg(a,g[i],c,g[i+1]);p!=0&&(d-=100-p);}if(n<d){n=d;o=i;}}return(e[o]+e[o+1])/2;}function rB(a,b){var c,d,e,f,g,h,i,j,k;if(b.length==0){return a.rb(jQ,gQ,-1,-1);}k=MJ(b);DJ(k.substr(0,3),'at ')&&(k=k.substr(3,k.length-3));k=k.replace(/\[.*?\]/g,'');g=k.indexOf('(');if(g==-1){g=k.indexOf('@');if(g==-1){j=k;k='';}else{j=MJ(k.substr(g+1,k.length-(g+1)));k=MJ(k.substr(0,g));}}else{c=k.indexOf(')',g);j=k.substr(g+1,c-(g+1));k=MJ(k.substr(0,g));}g=FJ(k,OJ(46));g!=-1&&(k=k.substr(g+1,k.length-(g+1)));(k.length==0||DJ(k,'Anonymous function'))&&(k=gQ);h=IJ(j,OJ(58));e=JJ(j,OJ(58),h-1);i=-1;d=-1;f=jQ;if(h!=-1&&e!=-1){f=j.substr(0,e);i=mB(j.substr(e+1,h-(e+1)));d=mB(j.substr(h+1,j.length-(h+1)));}return a.rb(f,k,i,d);}function GP(a,b,c){var d,e,f,g,h,i,j,k;f=0;for(j=0;j<c;){++f;e=a[b+j];if((e&192)==128){throw dG(new NI(mU));}else if((e&128)==0){++j;}else if((e&224)==192){j+=2;}else if((e&240)==224){j+=3;}else if((e&248)==240){j+=4;}else{throw dG(new NI(mU));}if(j>c){throw dG(new JH(mU));}}g=ZB(YC,gR,6,f,15,1);k=0;h=0;for(i=0;i<c;){e=a[b+i++];if((e&128)==0){h=1;e&=127;}else if((e&224)==192){h=2;e&=31;}else if((e&240)==224){h=3;e&=15;}else if((e&248)==240){h=4;e&=7;}else if((e&252)==248){h=5;e&=3;}while(--h>0){d=a[b+i++];if((d&192)!=128){throw dG(new NI('Invalid UTF8 sequence at '+(b+i-1)+', byte='+(d>>>0).toString(16)));}e=e<<6|d&63;}k+=YH(e,g,k);}return g;}function Kg(a){var b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s;while(true){s=0;n=0;q=null;r=null;for(g=1;g<a.f.a.length;g++){d=BM(a.f,g);for(h=0;h<g;h++){e=BM(a.f,h);b=0;c=0;o=0;p=0;for(k=0;k<d.a.length;k++){for(m=0;m<e.a.length;m++){if(d.a[k]===e.a[m]){++c;b=d.a[k];o<d.q[k]&&(o=d.q[k]);p<e.q[m]&&(p=e.q[m]);}}}if(c>0){f=c==1&&Bg(a,d,b)==1&&Bg(a,e,b)==1?0:1;o>p?i=(f<<24)+(o<<16)+(p<<8)+c:i=(f<<24)+(p<<16)+(o<<8)+c;if(s<i){s=i;n=c;o=0;p=0;for(l=0;l<d.a.length;l++)o<d.q[l]&&(o=d.q[l]);for(j=0;j<e.a.length;j++)p<e.q[j]&&(p=e.q[j]);if(o>p){q=d;r=e;}else{q=e;r=d;}}}}}if(s==0)break;n==q.a.length?GM(a.f,q):n==r.a.length?GM(a.f,r):Jg(a,q,r,n);}}function Uh(a){var b,c,d,e,f,g,h,i;for(g=0;g<a.p;g++){if(a.F[g]==128){c=a.B[0][g];d=a.B[1][g];if(a.A[c]==-1^a.A[d]==-1){if(a.q[c]!=0&&a.q[d]!=0){if(a.q[c]<0^a.q[d]<0){if(a.q[c]<0){++a.q[c];--a.q[d];}else{--a.q[c];++a.q[d];}}}}}}i=ZB(_C,DQ,6,a.o,15,1);e=0;for(b=0;b<a.o;b++){if(a.A[b]==-1){i[b]=-1;continue;}if(e<b){a.A[e]=a.A[b];a.q[e]=a.q[b];a.v[e]=a.v[b];a.s[e]=a.s[b];a.w[e]=a.w[b];a.u[e]=a.u[b];lh(a.H[e],a.H[b]);a.t!=null&&(a.t[e]=a.t[b]);a.r!=null&&(a.r[e]=a.r[b]);}i[b]=e;++e;}a.o=e;h=0;for(f=0;f<a.p;f++){if(a.F[f]==128)continue;a.F[h]=a.F[f];a.C[h]=a.C[f];a.D[h]=a.D[f];a.B[0][h]=i[a.B[0][f]];a.B[1][h]=i[a.B[1][f]];++h;}a.p=h;return i;}function Sm(a,b){var c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t;p=(i=Om(b,1),i==-1?b.length:i);f=xI(b.substr(0,p));o=Nm(b,p);p=(j=Om(b,o+1),j==-1?b.length:j);g=xI(b.substr(o,p-o));o=Nm(b,p);p=(k=Om(b,o+1),k==-1?b.length:k);c=Lm(a,xI(b.substr(o,p-o)));o=Nm(b,p);p=(l=Om(b,o+1),l==-1?b.length:l);d=Lm(a,xI(b.substr(o,p-o)));r=0;s=0;while((o=Nm(b,p))!=-1){p=(h=Om(b,o+1),h==-1?b.length:h);q=b.substr(o,p-o);n=FJ(q,OJ(61));m=q.substr(0,n);t=xI(q.substr(n+1,q.length-(n+1)));if(DJ(m,'CFG')){switch(t){case 1:r=1;break;case 2:r=g==2?3:4;break;case 3:r=6;}}else DJ(m,'TOPO')?s=t:undefined;}e=Jm(a,c,d,g,r,s);e+1!=f&&(!a.b&&(a.b=new hO()),aO(a.b,new QI(f),new QI(e)));}function Mg(a){var b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r;h=0;f=0;for(g=0;g<a.d;g++){if(lj(a.i,Ei(a.i,0,g))&&lj(a.i,Ei(a.i,1,g))){a.c[g]=true;f+=Ki(a.i,g);++h;}}if(h==0||f==0)return;f/=h;for(c=0;c<a.b;c++){lj(a.i,c)&&(a.e[c]==0?Sj(a.i,c,false):a.a[c]=true);}p=ZB(_C,DQ,6,a.b,15,1);i=il(a.i,p,true);o=ZB(_C,DQ,6,i,15,1);for(d=0;d<a.b;d++)p[d]!=-1&&++o[p[d]];n=ZB(qD,fR,29,i,0,1);for(k=0;k<i;k++)n[k]=new hh(a,a.i,o[k]);e=ZB(_C,DQ,6,i,15,1);for(b=0;b<a.b;b++){l=p[b];if(l!=-1){n[l].q[e[l]]=256;n[l].a[e[l]]=b;n[l].c[e[l]]=xi(a.i,b)/f;n[l].d[e[l]]=yi(a.i,b)/f;++e[l];}}q=-1;r=0;for(m=0;m<i;m++){if(r<o[m]){r=o[m];q=m;}}yM(a.f,n[q]);for(j=0;j<i;j++)j!=q&&yM(a.f,n[j]);}function xe(a){var b,c,d,e,f,g;g=a.N;f=ZB(_C,DQ,6,a.L.d,15,1);for(c=0;c<a.L.d;c++)f[c]=a.c[c];if(!a.L.I){ze(a);hf(a,g,f);}a.U=ZB(XC,pR,6,a.L.d,15,1);a.T=ZB(XC,pR,6,a.L.d,15,1);for(d=0;d<a.L.d;d++){a.U[d]=oi(a.L,d)<<24>>24;a.T[d]=ni(a.L,d)<<24>>24;}a.j=ZB(XC,pR,6,a.L.e,15,1);a.i=ZB(XC,pR,6,a.L.e,15,1);for(e=0;e<a.L.e;e++){a.j[e]=Ji(a.L,e)<<24>>24;a.i[e]=Ii(a.L,e)<<24>>24;}ye(a);a.Q=false;a.H=ZB(aG,HQ,6,a.L.d,16,1);for(b=0;b<a.L.d;b++){if(a.W[b]!=0){a.H[b]=true;a.Q=true;}}Be(a);a.J=null;a.V=ZB(aG,HQ,6,a.L.d,16,1);if(a.Q){a.J=new Zf(a.L,f,a.H,a.W,a.k,a.U,a.T,a.$,a.o,a.V);Uf(a.J);}a.Y=ZB(aG,HQ,6,a.L.d,16,1);a.Z=new LM();te(a);hf(a,g,f);Ae(a);!!a.J&&(a.F=Qf(a.J));Le(a);}function xg(a){var b,c,d,e,f,g,h,i,j,k,l,m;for(i=0;i<a.f.a.length;i++){h=BM(a.f,i);for(j=0;j<h.e.length;j++){d=h.e[j];if(Mi(a.i,d)==2){!Ol(a.i,d)&&Ni(a.i,d)==0&&hk(a.i,d);if(!Ll(a.i,d)&&bl(a.i,Ei(a.i,0,d))>1&&bl(a.i,Ei(a.i,1,d))>1&&(Ni(a.i,d)==1||Ni(a.i,d)==2)){m=ZB(_C,DQ,6,2,15,1);e=ZB(_C,DQ,6,2,15,1);for(k=0;k<2;k++){m[k]=a.i.K;e[k]=Ei(a.i,k,d);for(l=0;l<a.e[e[k]];l++){f=al(a.i,e[k],l);f!=Ei(a.i,1-k,d)&&m[k]>f&&(m[k]=f);}}g=rm(h.c[h.b[e[0]]],h.d[h.b[e[0]]],h.c[h.b[e[1]]],h.d[h.b[e[1]]]);b=rm(h.c[h.b[m[0]]],h.d[h.b[m[0]]],h.c[h.b[e[0]]],h.d[h.b[e[0]]]);c=rm(h.c[h.b[e[1]]],h.d[h.b[e[1]]],h.c[h.b[m[1]]],h.d[h.b[m[1]]]);Ag(g,b)<0^Ag(g,c)<0^Ni(a.i,d)==2&&Yg(h,d);}}}}}function Vg(a,b){var c,d,e,f,g,h,i,j,k,l,m;h=0;g=0;for(e=0;e<4;e++){f=$g(a,e)+$g(b,e>=2?e-2:e+2);if(h<f){h=f;g=e;}}k=(Wg(a),a.j-a.o+1+(Wg(b),b.j-b.o+1));l=(Wg(a),0.75*(a.i-a.n+1+(Wg(b),b.i-b.n+1)));i=$wnd.Math.max((Wg(a),a.j-a.o+1),(Wg(b),b.j-b.o+1));j=0.75*$wnd.Math.max((Wg(a),a.i-a.n+1),(Wg(b),b.i-b.n+1));d=$wnd.Math.sqrt((k-h)*(k-h)+(l-0.75*h)*(l-0.75*h));m=$wnd.Math.max(j,k);c=$wnd.Math.max(i,l);if(d<m&&d<c){switch(g){case 0:fh(b,a.i-b.n-h+1,a.o-b.j+h-1);break;case 1:fh(b,a.i-b.n-h+1,a.j-b.o-h+1);break;case 2:fh(b,a.n-b.i+h-1,a.j-b.o-h+1);break;case 3:fh(b,a.n-b.i+h-1,a.o-b.j+h-1);}}else c<m?fh(b,a.i-b.n+1,(a.j+a.o-b.j-b.o)/2):fh(b,(a.i+a.n-b.i-b.n)/2,a.j-b.o+1);}function xd(a){var b,c,d;if(a.G.o==0)return;Xo(a.G,(a.B&256)!=0?31:(a.B&512)!=0?47:(a.B&OQ)!=0?79:15);Lc(a);c=false;a.o=ZB(_C,DQ,6,a.G.o,15,1);for(b=0;b<a.G.o;b++){a.o[b]=ki(a.G,b);a.o[b]!=0&&(c=true);qj(a.G,b)&&(a.o[b]=128);Wi(a.G,b)&&(a.B&xQ)==0&&(a.o[b]=256);}zd(a,-10);Yc(a);Xc(a);Zc(a);Jc(a);vo(a,a.Q);uo(a,a.R);zd(a,a.J);_c(a);a.N.a=ZB(eF,fR,1,0,5,1);a.T.a=ZB(eF,fR,1,0,5,1);for(d=0;d<a.G.o;d++){if($c(a,d)){zd(a,-3);gd(a,d,true);zd(a,a.J);}else if(a.o[d]!=0){zd(a,a.o[d]);gd(a,d,true);zd(a,a.J);}else if(!c&&Ai(a.G,d)!=1&&Ai(a.G,d)!=6&&(a.B&YQ)==0&&qi(a.G,d)==null&&Ai(a.G,d)<wc.length){Ad(a,Uc(wc[Ai(a.G,d)]));gd(a,d,true);zd(a,a.J);}else{gd(a,d,true);}}fd(a);jd(a);ed(a);}function dO(a,b,c){var d,e,f,g,h,i,j,k,l,m,n;if(!a.b){return false;}g=null;m=null;i=new xO(null,null);e=1;i.a[1]=a.b;l=i;while(l.a[e]){j=e;h=m;m=l;l=l.a[e];d=a.a.eb(b,l.c);e=d<0?0:1;d==0&&(!c.c||IN(l.d,c.d))&&(g=l);if(!(!!l&&l.b)&&!_N(l.a[e])){if(_N(l.a[1-e])){m=m.a[j]=gO(l,e);}else if(!_N(l.a[1-e])){n=m.a[1-j];if(n){if(!_N(n.a[1-j])&&!_N(n.a[j])){m.b=false;n.b=true;l.b=true;}else{f=h.a[1]==m?1:0;_N(n.a[j])?h.a[f]=fO(m,j):_N(n.a[1-j])&&(h.a[f]=gO(m,j));l.b=h.a[f].b=true;h.a[f].a[0].b=false;h.a[f].a[1].b=false;}}}}}if(g){c.b=true;c.d=g.d;if(l!=g){k=new xO(l.c,l.d);eO(a,i,g,k);m==g&&(m=k);}m.a[m.a[1]==l?1:0]=l.a[!l.a[0]?1:0];--a.c;}a.b=i.a[1];!!a.b&&(a.b.b=false);return c.b;}function KB(a,b,c,d,e){var f,g,h,i;XJ(d,0,d.a.length);g=false;h=b.length;for(i=c;i<h;++i){f=b.charCodeAt(i);if(f==39){if(i+1<h&&b.charCodeAt(i+1)==39){++i;d.a+="'";}else{g=!g;}continue;}if(g){d.a+=String.fromCharCode(f);}else{switch(f){case 35:case 48:case 44:case 46:case 59:return i-c;case 164:a.g=true;if(i+1<h&&b.charCodeAt(i+1)==164){++i;if(i<h-2&&b.charCodeAt(i+1)==164&&b.charCodeAt(i+2)==164){i+=2;WJ(d,UB(a.a));}else{WJ(d,a.a[0]);}}else{WJ(d,a.a[1]);}break;case 37:if(!e){if(a.p!=1){throw dG(new NI(WT+b+'"'));}a.p=100;}d.a+='%';break;case 8240:if(!e){if(a.p!=1){throw dG(new NI(WT+b+'"'));}a.p=1000;}d.a+='\u2030';break;case 45:d.a+='-';break;default:d.a+=String.fromCharCode(f);}}}return h-c;}function An(a){var b,c,d,e,f,g,h,i,j;Xo(a.d,a.G);h=a.d.e+12;a.o=ZB(_C,DQ,6,h,15,1);a.q=ZB(_C,DQ,6,h,15,1);a.r=ZB(_C,DQ,6,h,15,1);a.p=ZB(aG,HQ,6,h+1,16,1);f=ZB(aG,HQ,6,a.d.d,16,1);g=ZB(aG,HQ,6,a.d.e,16,1);e=0;for(c=0;c<a.d.d;c++){if(!a.u[c]&&!f[c]){a.o[e]=c;a.r[e]=-1;a.q[e]=-1;i=e;while(e<=i){for(j=0;j<bl(a.d,a.o[e]);j++){d=al(a.d,a.o[e],j);a.u[d]||(i=Un(a,e,i,j,f,g));}while(a.p[++e]);}}}a.s=e;if(a.j!=0){i=e-1;e=0;while(e<=i){for(j=0;j<bl(a.d,a.o[e]);j++){d=al(a.d,a.o[e],j);(a.u[d]||a.u[a.o[e]])&&(i=Un(a,e,i,j,f,g));}while(a.p[++e]);}for(b=0;b<a.d.d;b++){if(a.u[b]&&!f[b]){a.o[e]=b;a.r[e]=-1;a.q[e]=-1;i=e;while(e<=i){for(j=0;j<bl(a.d,a.o[e]);j++){i=Un(a,e,i,j,f,g);}while(a.p[++e]);}}}}a.t=e;}function Xq(a,b){var c,d,e,f,g,h,i;i=new oq();if(!Pq){yM(i.a,new mq('Toxicity predictor not properly initialized.',2));return i;}g=Ze(new qf(a));if(Rq[b].hb(g)!=-1){nq(i,'This molecule is known to be '+Dq[b]+':',2);yM(i.a,new mq(g,1));return i;}h=new Wn();c=false;d=new kp();for(f=0;f<Oq[b].a.length;f++){mm(new om(false),d,BM(Oq[b],f));Pn(h,a);On(h,d);if(Fn(h,h.b)>0){c||nq(i,'High-risk fragments indicating '+Eq[b]+':',2);c=true;nq(i,BM(Oq[b],f),1);}}c=false;for(e=0;e<Qq[b].a.length;e++){mm(new om(false),d,BM(Qq[b],e));Pn(h,a);On(h,d);if(Fn(h,h.b)>0){c||nq(i,'Medium-risk fragments indicating '+Eq[b]+':',2);c=true;nq(i,BM(Qq[b],e),1);}}i.a.a.length==0&&nq(i,'No indication for '+Eq[b]+' found.',2);return i;}function Dj(a,b){var c,d,e,f,g,h,i,j,k;if(b==0)return 0;h=null;for(d=0;d<a.o;d++){if((a.s[d]&zR)>>19==b){h==null&&(h=ZB(aG,HQ,6,32,16,1));h[(a.s[d]&zR)>>19!=1&&(a.s[d]&zR)>>19!=2?-1:(a.s[d]&AR)>>21]=true;}}for(f=0;f<a.p;f++){if((a.C[f]&BR)>>10==b){h==null&&(h=ZB(aG,HQ,6,32,16,1));h[(a.C[f]&BR)>>10!=1&&(a.C[f]&BR)>>10!=2?-1:(a.C[f]&CR)>>12]=true;}}k=0;if(h!=null){j=ZB(_C,DQ,6,32,15,1);for(i=0;i<32;i++)h[i]&&(j[i]=k++);for(c=0;c<a.o;c++){if((a.s[c]&zR)>>19==b){g=j[(a.s[c]&zR)>>19!=1&&(a.s[c]&zR)>>19!=2?-1:(a.s[c]&AR)>>21];a.s[c]&=-65011713;a.s[c]|=g<<21;}}for(e=0;e<a.p;e++){if((a.C[e]&BR)>>10==b){g=j[(a.C[e]&BR)>>10!=1&&(a.C[e]&BR)>>10!=2?-1:(a.C[e]&CR)>>12];a.C[e]&=-126977;a.C[e]|=g<<12;}}}return k;}function Dn(a,b,c){var d,e,f,g,h;for(g=a.s;g<a.t;g++)c[g]=-1;f=a.s;while(true){++c[f];h=a.q[f]==-1?a.A.d:bl(a.A,a.w[a.q[f]]);if(c[f]==h){c[f]=-1;if(f==a.s)break;--f;if(!a.p[f]){b[a.w[a.o[f]]]=false;a.w[a.o[f]]=-1;}continue;}if(a.q[f]==-1){if(!b[c[f]]){if(xn(a,c[f],a.o[f])){a.w[a.o[f]]=c[f];b[c[f]]=true;++f;}}}else if(a.p[f]){e=al(a.A,a.w[a.q[f]],c[f]);e==a.w[a.o[f]]&&yn(a,cl(a.A,a.w[a.q[f]],c[f]),a.r[f])&&++f;}else{e=al(a.A,a.w[a.q[f]],c[f]);if(!b[e]){if(xn(a,e,a.o[f])&&yn(a,cl(a.A,a.w[a.q[f]],c[f]),a.r[f])){b[e]=true;a.w[a.o[f]]=e;++f;}}}if(f==a.t){if(En(a,true)&&Cn(a,true)&&Bn(a,b,true)){for(d=0;d<a.d.d;d++){if(a.u[d]){b[a.w[d]]=false;a.w[d]=-1;}}return true;}--f;if(!a.p[f]){b[a.w[a.o[f]]]=false;a.w[a.o[f]]=-1;}}}return false;}function Ql(a){var b,c,d,e,f,g,h,i,j,k;Xo(a,1);i=false;for(c=0;c<a.d;c++){if(a.A[c]==7&&a.q[c]==0){k=ol(a,c);if(k==4){for(j=0;j<a.g[c];j++){g=a.f[c][j];if(a.j[c][j]==1&&a.A[g]==8&&a.g[g]==1&&a.q[g]==0){i=true;++a.q[c];--a.q[g];break;}}}else if(k==5){for(j=0;j<a.g[c];j++){g=a.f[c][j];h=a.i[c][j];if(a.j[c][j]==2&&a.A[g]==8){i=true;++a.q[c];--a.q[g];a.F[h]=1;break;}if(a.j[c][j]==3&&a.A[g]==7){i=true;++a.q[c];--a.q[g];a.F[h]=2;break;}}}}}f=false;for(e=0;e<a.e;e++){for(j=0;j<2;j++){if(jj(a,a.B[j][e])){b=a.B[1-j][e];d=a.A[b];if(d==3||d==11||d==12||d==19||d==20||d==37||d==38||d==55||d==56){if(Mi(a,e)==1){++a.q[b];--a.q[a.B[j][e]];a.F[e]=128;f=true;}else if(a.F[e]==32){a.F[e]=128;f=true;}}break;}}}if(f){Uh(a);i=true;}i&&(a.Q=0);return i;}function ee(a,b,c){var d,e,f,g,h;if(a.k[b]!=0)return false;if(Mi(a.L,b)==1)return de(a,b,c);if(Mi(a.L,b)!=2)return false;if(Fl(a.L,b))return false;e=Ei(a.L,0,b);f=Ei(a.L,1,b);if(bl(a.L,e)==1||bl(a.L,f)==1)return false;if(bl(a.L,e)>3||bl(a.L,f)>3)return false;if(Tk(a.L,e)==2||Tk(a.L,f)==2)return false;g=new Bh(a.L,a.c,f,e);if(g.f&&!c)return false;h=new Bh(a.L,a.c,e,f);if(h.f&&!c)return false;if(g.f&&h.f)return false;if(c){g.f&&g.c&&(a.O[b]=true);h.f&&h.c&&(a.O[b]=true);}d=ij(a.L,b)?3:a._?ge(a,g,h):fe(g,h);if(c){if((a.K&2)!=0){if(g.f){if(d==1){Df(a.b[g.b],4);Df(a.b[g.d],1);}else if(d==2){Df(a.b[g.b],1);Df(a.b[g.d],4);}}if(h.f){if(d==1){Df(a.b[h.b],4);Df(a.b[h.d],1);}else if(d==2){Df(a.b[h.b],1);Df(a.b[h.d],4);}}}}else{a.k[b]=d;}return true;}function jd(a){var b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r;o=false;for(d=0;d<a.G.e;d++){j=null;if(ej(a.G,d)){l=Gi(a.G,d);k=Fi(a.G,d);j=l==k?'['+l+']':'['+l+':'+k+']';}else(Oi(a.G,d)&bR)!=0?j=(Oi(a.G,d)&bR)==cR?'a':(Oi(a.G,d)&48)==32?'r!a':'!a':(Oi(a.G,d)&48)!=0&&(j=(Oi(a.G,d)&48)==32?'r':'!r');n=(Oi(a.G,d)&SQ)>>14;n!=0&&(j=(j==null?'':j)+n);if(j!=null){b=Ei(a.G,0,d);c=Ei(a.G,1,d);if(!o){vo(a,(a.Q*2+1)/3|0);o=true;}q=(uh(a.K,xi(a.G,b))+uh(a.K,xi(a.G,c)))/2;r=(vh(a.K,yi(a.G,b))+vh(a.K,yi(a.G,c)))/2;f=uh(a.K,xi(a.G,c))-uh(a.K,xi(a.G,b));g=vh(a.K,yi(a.G,c))-vh(a.K,yi(a.G,b));e=$wnd.Math.sqrt(f*f+g*g);i=(m=(p=dH(a.e,j),new tH(0,0,p,0)).b,0.6*m);h=0.55*a.j;e!=0&&(f>0?ld(a,q+i*g/e,r-h*f/e,j,true,true):ld(a,q-i*g/e,r+h*f/e,j,true,true));}}o&&vo(a,a.Q);}function Op(b){var c,d,e,f,g,h,i,j,k;if(!b.g)return false;EH(b.f);EH(b.a);b.e=null;k=false;d=-1;b.b=b.c==null?null:ZB(kF,vR,2,b.c.length,6,1);b.d=-1;do{try{j=wH(b.g);if(j==null){EH(b.f);return false;}}catch(a){a=cG(a);if(PC(a,52)){EH(b.f);return false;}else throw dG(a);}if(k){WJ(b.a,j);TJ(b.a,10);}else{if(DJ(j.substr(0,1),'>')){k=true;WJ(b.f,MR);TJ(b.f,10);WJ(b.a,j);TJ(b.a,10);}else{WJ(b.f,j);TJ(b.f,10);DJ(j.substr(0,6),MR)&&(k=true);continue;}}if(b.c!=null){if(j.length==0){d=-1;}else if(d==-1){e=Qp(j);if(e!=null){d=-1;for(c=0;c<b.c.length;c++){if(DJ(e,b.c[c])){d=c;break;}}if(b.d==-1){for(g=Mp,h=0,i=g.length;h<i;++h){f=g[h];if(DJ(e,f)){b.d=d;break;}}}}}else{b.b[d]==null?b.b[d]=j:b.b[d]=BJ(BJ(b.b[d],kQ),j);}}}while(!DJ(j.substr(0,4),NR));return true;}function $m(a,b){var c,d,e,f,g;!!a.a&&WN(a.a);!!a.b&&WN(a.b);e=0;d=wH(b);while(d!=null&&DJ(d.substr(0,7),'M  V30 ')){d=MJ(d.substr(7,d.length-7));while(g='-'.length,DJ(d.substr(d.length-g,g),'-')){c=wH(b);if(!DJ(c.substr(0,7),'M  V30 ')){return false;}d=MJ(BJ(LJ(d,0,d.length-1),c.substr(7,c.length-7)));}if(DJ(d.substr(0,5),'BEGIN')){f=MJ(d.substr(6,d.length-6));if(DJ(f.substr(0,4),'CTAB')){e=1;}else if(DJ(f.substr(0,4),'ATOM')){e=2;}else if(DJ(f.substr(0,4),'BOND')){e=3;}else if(DJ(f.substr(0,10),'COLLECTION')){e=4;}else{return false;}}else if(DJ(d.substr(0,3),'END')){e=0;}else if(e==1){Um(a,d);}else if(e==2){Qm(a,d);}else if(e==3){Sm(a,d);}else if(e==4){Tm(a,d);}else{return false;}d=wH(b);}while(d!=null&&!(DJ(d.substr(0,6),MR)||DJ(d,NR))){d=wH(b);}return true;}function Tm(a,b){var c,d,e,f,g,h;h=Pm(b);if(h!=null){g=Vm(b,h);if(DJ(b.substr(0,13),'MDLV30/STEABS')){if(DJ(h,LR))for(f=0;f<g.length;f++)Oj(a.c,Lm(a,g[f]),0,-1);else for(e=0;e<g.length;e++)ek(a.c,Mm(a,g[e]),0,-1);}else if(DJ(b.substr(0,13),'MDLV30/STERAC')){d=xI(LJ(b,13,Om(b,13)));if(DJ(h,LR))for(f=0;f<g.length;f++)Oj(a.c,Lm(a,g[f]),1,d-1);else for(e=0;e<g.length;e++)ek(a.c,Mm(a,g[e]),1,d-1);}else if(DJ(b.substr(0,13),'MDLV30/STEREL')){d=xI(LJ(b,13,Om(b,13)));if(DJ(h,LR))for(f=0;f<g.length;f++)Oj(a.c,Lm(a,g[f]),2,d-1);else for(e=0;e<g.length;e++)ek(a.c,Mm(a,g[e]),2,d-1);}else if(DJ(b.substr(0,13),'MDLV30/HILITE')){if(DJ(h,LR)){for(e=0;e<g.length;e++)Kj(a.c,Lm(a,g[e]),448);}else{for(e=0;e<g.length;e++){c=Mm(a,g[e]);Kj(a.c,Ei(a.c,0,c),448);Kj(a.c,Ei(a.c,1,c),448);}}}}}function ce(a,b,c){var d,e,f,g,h,i,j;if(Ai(a.L,b)!=6&&Ai(a.L,b)!=7)return false;e=al(a.L,b,0);f=al(a.L,b,1);if(Tk(a.L,e)!=1||Tk(a.L,f)!=1)return false;if(bl(a.L,e)==1||bl(a.L,f)==1)return false;if(Qk(a.L,e)>3||Qk(a.L,f)>3)return false;g=new Bh(a.L,a.c,b,e);if(g.f&&!c)return false;h=new Bh(a.L,a.c,b,f);if(h.f&&!c)return false;if(g.f&&h.f)return false;if(c){g.f&&g.c&&(a.P[b]=true);h.f&&h.c&&(a.P[b]=true);}i=Ah(g);j=Ah(h);if(i==-1||j==-1||(i+j&1)==0){c||(a.W[b]=3);return true;}d=0;switch(i+j){case 3:case 7:d=2;break;case 5:d=1;}if(c){if(a.Q&&(a.K&2)!=0||!a.Q&&(a.K&4)!=0){if(g.f){if(d==1){Df(a.b[g.b],64);Df(a.b[g.d],16);}else{Df(a.b[g.b],16);Df(a.b[g.d],64);}}if(h.f){if(d==2){Df(a.b[h.b],64);Df(a.b[h.d],16);}else{Df(a.b[h.b],16);Df(a.b[h.d],64);}}}}else{a.W[b]=d;}return true;}function En(a,b){var c,d,e,f,g,h,i,j,k,l,m,n,o,p;g=0;for(i=0;i<a.d.d;i++){if(a.u[i]==b&&(vi(a.d,i)&iR)!=0){m=a.w[i];l=ui(a.d,i);o=ui(a.A,m);if(l==0)continue;if(o==0)continue;if(l==3)continue;if(o==3)continue;if(oi(a.d,i)==1){++g;continue;}if(oi(a.A,m)==1)return false;if(oi(a.d,i)==2){++g;continue;}if(oi(a.A,m)==2)return false;if(Mn(a,i)==(l==o))return false;}}if(g!=0){e=ZB(_C,DQ,6,g,15,1);f=0;for(j=0;j<a.d.d;j++){if(a.u[j]==b&&(vi(a.d,j)&iR)!=0){l=ui(a.d,j);l!=0&&l!=3&&(e[f++]=ni(a.d,j)<<24|oi(a.d,j)<<22|j);}}xN(e);f=0;while(f<e.length){k=e[f]&OR;n=a.w[k];c=e[f]&-4194304;d=Mn(a,k)^ui(a.d,k)==ui(a.A,n);for(++f;f<e.length&&(e[f]&-4194304)==c;f++){h=e[f]&OR;m=a.w[h];if(oi(a.A,m)!=oi(a.A,n)||ni(a.A,m)!=ni(a.A,n))return false;p=Mn(a,h)^ui(a.d,h)==ui(a.A,m);if(p!=d)return false;}}}return true;}function _p(b){var c,d,e,f,g,h,i,j,k,l,m,n,o;f=new oq();yM(f.a,new mq('cLogP Values are estimated applying an atom-type based increment system.',2));yM(f.a,new mq(LT,2));yM(f.a,new mq(MT,2));Ql(b);Xo(b,3);if(b){i=0;e=new hO();j=new mK('#0.000');for(c=0;c<b.d;c++){try{d=_d(b,c,6241);o=cM(e,new aJ(d));!o?aO(e,new aJ(d),new QI(1)):aO(e,new aJ(d),new QI(o.a+1));}catch(a){a=cG(a);if(PC(a,12)){++i;}else throw dG(a);}}i!=0&&yM(f.a,new mq('Warning: '+i+' atom type(s) could not be determined.',2));for(n=(h=new qO(new vO(new mM(e).a).b),new sM(h));RK(n.a.a);){m=(g=oO(n.a),g.Fb());l=rA(Yp,m);(l<0?-1:l)!=-1?nq(f,mL(XN(e,m))+' * '+lK(j,Xp[(k=rA(Yp,m),k<0?-1:k)])+' AtomType: 0x'+fJ(m.a),2):nq(f,'Warning: For atom type 0x'+fJ(m.a)+' ('+mL(XN(e,m))+' times found) is no increment available.',2);}}return f;}function mf(a){var b,c,d,e,f,g,h,i,j,k,l;for(b=0;b<a.L.d;b++){if(a.W[b]==1||a.W[b]==2){i=false;if(Tk(a.L,b)!=0&&bl(a.L,b)==2&&dl(a.L,b,0)==2&&dl(a.L,b,1)==2){for(h=0;h<bl(a.L,b);h++){e=al(a.L,b,h);l=0;k=ZB(_C,DQ,6,3,15,1);for(j=0;j<bl(a.L,e);j++){k[l]=al(a.L,e,j);k[l]!=b&&++l;}l==2&&a.c[k[0]]>a.c[k[1]]^k[0]<k[1]&&(i=!i);}}else{for(h=1;h<bl(a.L,b);h++){for(j=0;j<h;j++){f=al(a.L,b,h);g=al(a.L,b,j);a.c[f]>a.c[g]&&(i=!i);f<g&&(i=!i);}}}Uj(a.L,b,a.W[b]==1^i?1:2,a.X[b]);}else{Uj(a.L,b,a.W[b],a.X[b]);}}for(c=0;c<a.L.e;c++){if(a.k[c]==1||a.k[c]==2){i=false;for(h=0;h<2;h++){d=Ei(a.L,h,c);if(bl(a.L,d)==3){k=ZB(_C,DQ,6,2,15,1);l=0;for(j=0;j<3;j++)al(a.L,d,j)!=Ei(a.L,1-h,c)&&(k[l++]=al(a.L,d,j));a.c[k[0]]>a.c[k[1]]&&(i=!i);k[0]<k[1]&&(i=!i);}}gk(a.L,c,a.k[c]==1^i?1:2,a.n[c]);}else{gk(a.L,c,a.k[c],a.n[c]);}}}function Rg(a){var b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u;b=vg(a);f=ZB(XC,pR,6,a.d,15,1);Ng(a,f,b);for(e=0;e<a.d;e++)f[e]==2&&(Kl(a.i,Ei(a.i,0,e))||Kl(a.i,Ei(a.i,1,e)))&&(f[e]=3);for(n=0;n<a.f.a.length;n++){l=BM(a.f,n);i=Zg(l);r=l.f;q=new gh(a,l);p=-1;for(m=0;m<224&&i.a.length!=0;m++){j=ON(a.j,i.a.length);h=(NP(j,i.a.length),i.a[j]);g=Fg(a,h[0],h[1]);c=ZB(_C,DQ,6,g.length,15,1);d=0;if(m<32){for(o=1;o<g.length-1;o++)f[g[o]]==3&&(c[d++]=g[o]);}else if(m<96){for(o=1;o<g.length-1;o++)f[g[o]]>=2&&(c[d++]=g[o]);}else{for(o=1;o<g.length-1;o++)f[g[o]]>=1&&(c[d++]=g[o]);}if(d!=0){t=c[0];if(d>1){do{t=c[ON(a.j,d)];}while(t==p);}if(t!=p){p=t;Yg(l,t);i=Zg(l);if(r>l.f){r=l.f;q=new gh(a,l);}}}}IM(a.f,n,q);l=q;k=1;do{s=9999;for(o=0;o<l.a.length;o++){u=b[l.a[o]];u==k?dh(l,o):u>k&&u<s&&(s=u);}k=s;}while(s!=9999);}}function Zl(a,b){var c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,A;if((a.C[b]&3)==0||(a.C[b]&3)==3||!Gl(a,b))return;v=-1;t=-1;u=-1;s=-1;e=0;for(l=0;l<2;l++){d=a.B[l][b];for(o=0;o<a.c[d];o++){h=a.i[d][o];if(h!=b&&Mi(a,h)==1){g=a.f[d][o];w=xl(a,h,g);if(e<w){e=w;t=g;v=h;u=d;s=a.B[1-l][b];}}}}if(t==-1)return;for(m=0;m<2;m++){for(o=0;o<a.c[a.B[m][b]];o++){h=a.i[a.B[m][b]][o];h!=b&&Mi(a,h)==1&&(a.F[h]=1);}}if(a.B[1][v]!=t){a.B[0][v]=a.B[1][v];a.B[1][v]=t;}i=oQ;for(n=0;n<a.g[u];n++){g=a.f[u][n];a.i[u][n]!=b&&i>g&&(i=g);}q=ZB(_C,DQ,6,2,15,1);r=0;for(k=0;k<a.g[s];k++)a.i[s][k]!=b&&(q[r++]=a.f[s][k]);f=Ak(a.H[u].a,a.H[u].b,a.H[s].a,a.H[s].b);if(r==2){if(q[0]>q[1]){A=q[0];q[0]=q[1];q[1]=A;}j=Bk(f,Di(a,s,q[0]));p=Bk(f,Di(a,s,q[1]));c=j-p;}else{c=Bk(f,Di(a,s,q[0]));}c<0^(a.C[b]&3)==2^i==t?a.F[v]=17:a.F[v]=9;}function Il(a,b){var c,d,e,f,g,h,i,j,k,l,m;if(a.A[b]!=7)return false;if((a.s[b]&xQ)!=0||a.k[b]!=0||(a.w[b]&XQ)!=0)return true;if(a.q[b]==1)return false;f=0;for(h=0;h<a.g[b];h++){if(a.j[b][h]==1){c=a.A[a.f[b][h]];(c==8||c==9||c==17)&&++f;}}if(f==0){for(g=0;g<a.g[b];g++){d=a.f[b][g];if(a.k[d]!=0){if((a.s[d]&xQ)!=0){if((!!a.n&&d<a.d?jn(a.n,d):0)>=5){m=0;for(k=0;k<a.g[d];k++){l=a.f[d][k];l!=b&&a.g[l]>=3&&++m;}if(m==2||m==1&&a.g[b]==3)continue;}return true;}for(j=0;j<a.g[d];j++){if((a.j[d][j]==2||Fl(a,a.i[d][j]))&&Pl(a,a.f[d][j]))return true;}}}}if(f<2){for(g=0;g<a.g[b];g++){d=a.f[b][g];i=false;e=false;for(j=0;j<a.g[d];j++){if(a.f[d][j]!=b){a.j[d][j]!=1&&(a.A[a.f[d][j]]==7||a.A[a.f[d][j]]==8||a.A[a.f[d][j]]==16)&&(i=true);a.j[d][j]==1&&a.A[a.f[d][j]]==7&&(e=true);}}if(i&&(!e||f==0))return true;}}return false;}function Lk(a,b){var c,d,e,f,g,h,i,j,k,l,m;if((b&~a.Q)==0)return;if((a.Q&1)==0){Al(a);Gk(a);a.Q|=1;if(dm(a)){Al(a);Gk(a);}}if((b&~a.Q)==0)return;if((a.Q&2)==0){for(d=0;d<a.d;d++)a.s[d]&=-31753;for(g=0;g<a.e;g++)a.C[g]&=-961;Pk(a);for(f=0;f<a.e;f++){if(a.F[f]==64){a.s[a.B[0][f]]|=xQ;a.s[a.B[1][f]]|=xQ;a.C[f]|=256;a.C[f]|=512;}}for(e=0;e<a.d;e++){for(l=0;l<a.g[e];l++){j=a.i[e][l];if((a.C[j]&256)!=0)continue;i=a.f[e][l];for(m=0;m<a.g[i];m++){if(a.i[i][m]==j)continue;a.j[i][m]>1&&(a.A[a.f[i][m]]==6?a.s[e]|=iR:!Fl(a,a.i[i][m])&&jj(a,a.f[i][m])&&(a.s[e]|=yQ));}}}while(true){k=false;for(c=0;c<a.d;c++){if(a.k[c]>0&&(20480&a.s[c])==yQ){for(l=0;l<a.g[c];l++){if(a.j[c][l]>1){i=a.f[c][l];j=a.i[c][l];for(m=0;m<a.g[i];m++){if(a.i[i][m]!=j){h=a.f[i][m];if((a.s[h]&yQ)==0){a.s[h]|=yQ;k=true;}}}}}}}if(!k)break;}a.Q|=2;}}function Vh(a,b,c,d,e){var f,g,h,i;f=b.o;f>=b.K&&nk(b,b.K*2);h=(a.s[c]&zR)>>19;g=-1;h==1?d==-1?g=Dj(b,1):g=mJ(32,d+((a.s[c]&zR)>>19!=1&&(a.s[c]&zR)>>19!=2?-1:(a.s[c]&AR)>>21)):h==2&&(e==-1?g=Dj(b,2):g=mJ(32,e+((a.s[c]&zR)>>19!=1&&(a.s[c]&zR)>>19!=2?-1:(a.s[c]&AR)>>21)));b.A[f]=a.A[c];b.q[f]=a.q[c];b.v[f]=a.v[c];b.s[f]=a.s[c];b.w[f]=b.I?a.w[c]:0;lh(b.H[f],a.H[c]);b.u[f]=a.u[c];b.t!=null&&(b.t[f]=null);if(a.t!=null&&a.t[c]!=null&&b.I){b.t==null&&(b.t=ZB(_C,mR,7,b.A.length,0,2));b.t[f]=ZB(_C,DQ,6,a.t[c].length,15,1);for(i=0;i<a.t[c].length;i++)b.t[f][i]=a.t[c][i];}b.r!=null&&(b.r[f]=null);if(a.r!=null&&a.r[c]!=null){b.r==null&&(b.r=ZB(XC,xR,13,b.A.length,0,2));b.r[f]=ZB(XC,pR,6,a.r[c].length,15,1);for(i=0;i<a.r[c].length;i++)b.r[f][i]=a.r[c][i];}if(g!=-1){b.s[f]&=-65011713;b.s[f]|=g<<21;}++b.o;b.Q=0;return f;}function ge(a,b,c){var d,e,f,g,h,i,j;f=ZB(ZC,GQ,6,3,15,1);f[0]=xi(a.L,c.a)-xi(a.L,b.a);f[1]=yi(a.L,c.a)-yi(a.L,b.a);f[2]=zi(a.L,c.a)-zi(a.L,b.a);i=ZB(ZC,GQ,6,3,15,1);i[0]=xi(a.L,b.b)-xi(a.L,b.a);i[1]=yi(a.L,b.b)-yi(a.L,b.a);i[2]=zi(a.L,b.b)-zi(a.L,b.a);j=ZB(ZC,GQ,6,3,15,1);j[0]=xi(a.L,c.b)-xi(a.L,c.a);j[1]=yi(a.L,c.b)-yi(a.L,c.a);j[2]=zi(a.L,c.b)-zi(a.L,c.a);g=ZB(ZC,GQ,6,3,15,1);g[0]=f[1]*i[2]-f[2]*i[1];g[1]=f[2]*i[0]-f[0]*i[2];g[2]=f[0]*i[1]-f[1]*i[0];h=ZB(ZC,GQ,6,3,15,1);h[0]=f[1]*g[2]-f[2]*g[1];h[1]=f[2]*g[0]-f[0]*g[2];h[2]=f[0]*g[1]-f[1]*g[0];d=(i[0]*h[0]+i[1]*h[1]+i[2]*h[2])/($wnd.Math.sqrt(i[0]*i[0]+i[1]*i[1]+i[2]*i[2])*$wnd.Math.sqrt(h[0]*h[0]+h[1]*h[1]+h[2]*h[2]));e=(j[0]*h[0]+j[1]*h[1]+j[2]*h[2])/($wnd.Math.sqrt(j[0]*j[0]+j[1]*j[1]+j[2]*j[2])*$wnd.Math.sqrt(h[0]*h[0]+h[1]*h[1]+h[2]*h[2]));return d<0^e<0?1:2;}function ef(a){var b,c,d,e,f,g,h,i,j,k,l;a.S=ZB(XC,pR,6,a.L.d,15,1);for(b=0;b<a.L.d;b++){if(a.W[b]==1||a.W[b]==2){i=false;if(bl(a.L,b)==2&&dl(a.L,b,0)==2&&dl(a.L,b,1)==2){for(h=0;h<bl(a.L,b);h++){e=al(a.L,b,h);l=0;k=ZB(_C,DQ,6,3,15,1);for(j=0;j<bl(a.L,e);j++){k[l]=al(a.L,e,j);k[l]!=b&&++l;}l==2&&a.c[k[0]]>a.c[k[1]]^a.B[k[0]]<a.B[k[1]]&&(i=!i);}}else{for(h=1;h<bl(a.L,b);h++){for(j=0;j<h;j++){f=al(a.L,b,h);g=al(a.L,b,j);a.c[f]>a.c[g]&&(i=!i);a.B[f]<a.B[g]&&(i=!i);}}}a.S[b]=a.W[b]==1^i?1:2;}else{a.S[b]=a.W[b];}}a.g=ZB(XC,pR,6,a.L.e,15,1);for(c=0;c<a.L.e;c++){if(a.k[c]==1||a.k[c]==2){i=false;for(h=0;h<2;h++){d=Ei(a.L,h,c);if(bl(a.L,d)==3){k=ZB(_C,DQ,6,2,15,1);l=0;for(j=0;j<3;j++)al(a.L,d,j)!=Ei(a.L,1-h,c)&&(k[l++]=al(a.L,d,j));a.c[k[0]]>a.c[k[1]]&&(i=!i);a.B[k[0]]<a.B[k[1]]&&(i=!i);}}a.g[c]=a.k[c]==1^i?1:2;}else{a.g[c]=a.k[c];}}}function Id(a){var b,c,d,e,f,g,h,i,j,k,l,m,n,o;if(a.b==0)return true;o=new vn(a.d,1);a.d.I&&Kd(a);for(c=0;c<a.d.d;c++){if(Hd(a,c)){o.a[c]==7&&(Ai(a.d,c)==5&&ji(a.d,c)==0||Ai(a.d,c)==6&&ji(a.d,c)==1)&&Nd(a,c);o.a[c]==5&&(Ai(a.d,c)==6&&ji(a.d,c)==-1||Ai(a.d,c)==7&&ji(a.d,c)==0&&Qk(a.d,c)==3||Ai(a.d,c)==8&&ji(a.d,c)==0&&bl(a.d,c)==2||Ai(a.d,c)==16&&ji(a.d,c)==0&&bl(a.d,c)==2)&&Nd(a,c);}}Md(a,o);Od(a);Ld(a);while(a.b!=0){g=false;for(e=0;e<a.d.e;e++){if(a.c[e]){b=0;for(k=0;k<2;k++){f=Ei(a.d,k,e);for(l=0;l<bl(a.d,f);l++)a.c[cl(a.d,f,l)]&&++b;}if(b==4){Jd(a,e);Ld(a);g=true;break;}}}if(!g){for(m=0;m<o.g.a.length;m++){if(BM(o.i,m).length==6){j=true;n=BM(o.i,m);for(i=0;i<6;i++){if(!a.c[n[i]]){j=false;break;}}if(j){for(h=0;h<6;h+=2)Jd(a,n[h]);g=true;break;}}}}if(!g){for(d=0;d<a.d.e;d++){if(a.c[d]){Jd(a,d);Ld(a);break;}}}}return a.a==a.e;}function Hk(a,b){var c,d,e,f,g,h,i,j,k,l,m,n,o;Xo(a,1);for(h=0;h<a.p;h++){if(Mi(a,h)<3){if(a.q[a.B[0][h]]>0&&a.q[a.B[1][h]]<0){e=a.B[0][h];f=a.B[1][h];}else if(a.q[a.B[0][h]]<0&&a.q[a.B[1][h]]>0){e=a.B[1][h];f=a.B[0][h];}else continue;if(a.A[e]<9)if(ol(a,e)>3)continue;a.q[e]-=1;a.q[f]+=1;Mi(a,h)==1?a.F[h]=2:a.F[h]=4;a.Q=0;}}n=0;m=0;k=0;for(d=0;d<a.o;d++){n+=a.q[d];if(a.q[d]<0){++m;Dk(a.A[d])&&(k+=a.q[d]);}}if(!b&&n!=0)throw dG(new FA("molecule's overall charges are not balanced"));Xo(a,1);o=0;for(g=0;g<a.o;g++){if(a.q[g]>0){if(!Bl(a,g)&&Dk(a.A[g])){i=mJ(ml(a,g),a.q[g]);if(i!=0&&k>=i){o-=i;k+=i;a.q[g]-=i;a.Q&=1;}}}}if(o<0){l=ZB(_C,DQ,6,m,15,1);m=0;for(e=0;e<a.o;e++){a.q[e]<0&&(Cl(a,e)||(l[m++]=(a.A[e]<<16)+e));}xN(l);for(j=l.length-1;n<0&&j>=l.length-m;j--){c=l[j]&AQ;if(Dk(a.A[c])){i=mJ(-o,-a.q[c]);o+=i;a.q[c]+=i;a.Q&=1;}}}return n;}function Gn(a,b){var c,d,e,f,g,h,i,j;i=0;if(a.I){(a.s[b]&xQ)!=0&&(i|=2);j=(d=a.s[b]&BR,d==0?0:d==OQ?2:d==YQ?3:4);if(j!=0){i|=8;j>2&&(i|=16);j>3&&(i|=32);}c=a.q[b];c<0?i|=RQ:c>0&&(i|=QQ);f=a.g[b];switch(f){case 0:break;case 1:i|=qR;break;case 2:i|=VQ;break;case 3:i|=917504;break;default:i|=1966080;}}else{(a.s[b]&xQ)!=0?i|=2:i|=4;j=(d=a.s[b]&BR,d==0?0:d==OQ?2:d==YQ?3:4);j==0?i|=112:j==2?i|=104:j==3?i|=88:i|=56;c=a.q[b];c==0?i|=167772160:c<0?i|=RQ:c>0&&(i|=QQ);e=a.c[b]-a.g[b]+ml(a,b);switch(e){case 0:i|=1792;break;case 1:i|=1664;break;case 2:i|=1408;break;default:i|=896;}f=a.g[b];switch(f){case 0:i|=3932160;break;case 1:i|=3801088;break;case 2:i|=3538944;break;case 3:i|=3014656;break;default:i|=1966080;}h=a.k[b];switch(h){case 0:i|=98304;break;case 1:i|=81920;break;default:i|=49152;}}g=a.k[b];g>0&&(i|=yQ);g>1&&(i|=32768);return i;}function Vl(a,b){var c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,A,B;if(a.g[b]!=2||a.j[b][0]!=2||a.j[b][1]!=2||a.g[a.f[b][0]]<2||a.g[a.f[b][1]]<2||a.k[a.f[b][0]]!=1||a.k[a.f[b][1]]!=1){Uj(a,b,0,false);return;}w=-1;v=-1;u=-1;r=-1;f=0;for(l=0;l<2;l++){d=a.f[b][l];for(p=0;p<a.c[d];p++){g=a.f[d][p];if(g!=b){h=a.i[d][p];A=xl(a,h,g);if(f<A){f=A;v=g;w=h;u=d;r=a.f[b][1-l];}}}}if(v==-1)return;for(m=0;m<2;m++)for(o=0;o<a.c[a.f[b][m]];o++)a.f[a.f[b][m]][o]!=b&&(a.F[a.i[a.f[b][m]][o]]=1);if(a.B[1][w]!=v){a.B[0][w]=a.B[1][w];a.B[1][w]=v;}i=oQ;for(n=0;n<a.g[u];n++){g=a.f[u][n];g!=b&&i>g&&(i=g);}s=ZB(_C,DQ,6,2,15,1);t=0;for(k=0;k<a.g[r];k++){g=a.f[r][k];g!=b&&(s[t++]=g);}c=Ak(a.H[b].a,a.H[b].b,a.H[r].a,a.H[r].b);if(t==2){if(s[0]>s[1]){B=s[0];s[0]=s[1];s[1]=B;}j=Bk(c,Di(a,r,s[0]));q=Bk(c,Di(a,r,s[1]));e=j-q;}else{e=Bk(c,Di(a,r,s[0]));}e<0^(a.s[b]&3)==1^i==v?a.F[w]=17:a.F[w]=9;}function Le(a){var b,c,d,e,f,g,h,i,j,k,l,m;f=0;k=0;g=0;h=0;i=0;j=0;l=0;m=false;b=ZB(aG,HQ,6,32,16,1);for(c=0;c<a.L.d;c++){if(a.W[c]!=0){++f;if(a.W[c]==3){++k;}else{if(a.U[c]==0){++g;!!a.J&&Pf(a.J,c)&&++h;}else if(a.U[c]==2){a.T[c]==0&&++j;}else if(a.U[c]==1){e=a.T[c];if(!b[e]){++l;b[e]=true;}a.T[c]==0&&++i;!!a.J&&Pf(a.J,c)&&(m=true);}}}}for(d=0;d<a.L.e;d++){if(a.k[d]!=0&&Pi(a.L,d)==1){++f;if(a.k[d]==3){++k;}else{if(a.j[d]==0){++g;!!a.J&&Pf(a.J,Ei(a.L,0,d))&&Pf(a.J,Ei(a.L,1,d))&&++h;}else if(a.j[d]==2){a.i[d]==0&&++j;}else if(a.j[d]==1){e=a.i[d];if(!b[e]){++l;b[e]=true;}a.i[d]==0&&++i;!!a.J&&Pf(a.J,Ei(a.L,0,d))&&Pf(a.J,Ei(a.L,1,d))&&(m=true);}}}}if(f==0){a.L.G=zQ;return;}if(k!=0){a.L.G=0;return;}if(a.F){kk(a.L,qR+(1<<l));return;}i+h==f&&!m?(a.L.G=196608,undefined):g==f?(a.L.G=cR,undefined):j==f?(a.L.G=327680,undefined):g==f-1&&i==1?(a.L.G=VQ,undefined):kk(a.L,458752+(1<<l));}function Qe(a,b){var c,d,e,f,g,h,i,j,k,l,m,n,o,p;if(a.L.d==0){a.e='';return;}k=false;if(a._&&a.L.o>a.L.d&&!a.L.I){k=true;for(h=0;h<a.L.d;h++){if(ml(a.L,h)!=0){k=false;break;}}}p=a._?16:8;Pe(a);TJ(a.q,k?35:33);Ne(a,a._?1:0,1);Ne(a,b?1:0,1);Ne(a,p/2|0,4);n=0;for(i=1;i<a.L.d;i++)n=$e(a,a.t[i],a.w[i]==-1?-1:a.t[a.w[i]],n);if(k){for(h=0;h<a.L.d;h++){c=a.t[h];for(m=qe(a,c);m<Qk(a.L,c);m++)n=$e(a,al(a.L,c,m),c,n);}}if(n==0){a.e='';return;}f=1<<p;l=n/(f/2-1);o=n+l/2;for(j=1;j<a.L.d;j++)Me(a,a.t[j],a.w[j]==-1?-1:a.t[a.w[j]],o,l,p);if(k){for(g=0;g<a.L.d;g++){c=a.t[g];for(m=qe(a,c);m<Qk(a.L,c);m++)Me(a,al(a.L,c,m),c,o,l,p);}}if(b){e=a._?1.5:(Gh(),Gh(),Fh);d=Ci(a.L,a.L.d,a.L.e,e);Ne(a,mJ(f-1,lJ(0,WC(0.5+$wnd.Math.log(d/0.1)*$wnd.Math.LOG10E/($wnd.Math.log(2000)*$wnd.Math.LOG10E)*(f-1)))),p);Ne(a,Re(xi(a.L,a.t[0])/d,f),p);Ne(a,Re(yi(a.L,a.t[0])/d,f),p);a._&&Ne(a,Re(zi(a.L,a.t[0]),f),p);}a.e=Oe(a);}function Gk(a){var b,c,d,e,f,g,h,i,j,k,l;a.g=ZB(_C,DQ,6,a.o,15,1);a.c=ZB(_C,DQ,6,a.o,15,1);a.f=ZB(_C,mR,7,a.o,0,2);a.i=ZB(_C,mR,7,a.o,0,2);a.j=ZB(_C,mR,7,a.o,0,2);a.k=ZB(_C,DQ,6,a.d,15,1);i=ZB(_C,DQ,6,a.o,15,1);for(f=0;f<a.p;f++){++i[a.B[0][f]];++i[a.B[1][f]];}for(c=0;c<a.o;c++){a.f[c]=ZB(_C,DQ,6,i[c],15,1);a.i[c]=ZB(_C,DQ,6,i[c],15,1);a.j[c]=ZB(_C,DQ,6,i[c],15,1);}k=false;for(g=0;g<a.e;g++){l=Mi(a,g);if(l==0){k=true;continue;}for(j=0;j<2;j++){d=a.B[j][g];a.j[d][a.c[d]]=l;a.f[d][a.c[d]]=a.B[1-j][g];a.i[d][a.c[d]]=g;++a.c[d];++a.g[d];d<a.d&&(l>1?a.k[d]+=l+l-2:a.F[g]==64&&(a.k[d]=2));}}if(k){for(h=0;h<a.e;h++){l=Mi(a,h);if(l==0){for(j=0;j<2;j++){d=a.B[j][h];a.j[d][a.c[d]]=0;a.f[d][a.c[d]]=a.B[1-j][h];a.i[d][a.c[d]]=h;++a.c[d];}}}}for(e=a.e;e<a.p;e++){for(j=0;j<2;j++){d=a.B[j][e];a.j[d][a.c[d]]=1;a.f[d][a.c[d]]=a.B[1-j][e];a.i[d][a.c[d]]=e;++a.c[d];}}for(b=0;b<a.d;b++)a.k[b]=a.k[b]/2|0;}function Yg(a,b){var c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v;a.g==null&&(a.g=ZB(_C,mR,7,a.r.d,0,2));if(a.g[b]==null){m=ZB(_C,DQ,6,a.a.length,15,1);s=ZB(aG,HQ,6,a.r.b,16,1);d=Ei(a.p,0,b);e=Ei(a.p,1,b);m[0]=d;s[d]=true;j=0;n=0;while(j<=n){for(p=0;p<a.r.e[m[j]];p++){f=al(a.p,m[j],p);if(!s[f]&&f!=e){m[++n]=f;s[f]=true;}}if(j==n)break;++j;}l=n+1>(a.a.length/2|0);if((a.r.g&6)!=0){h=false;g=false;for(p=0;p<a.a.length;p++){lj(a.p,a.a[p])&&(s[a.a[p]]?h=true:g=true);}h!=g&&(l=h);}i=2;a.g[b]=ZB(_C,DQ,6,l?a.a.length-n:n+2,15,1);for(q=0;q<a.a.length;q++){a.a[q]==d?a.g[b][l?0:1]=q:a.a[q]==e?a.g[b][l?1:0]=q:l^s[a.a[q]]&&(a.g[b][i++]=q);}}u=a.c[a.g[b][0]];v=a.d[a.g[b][0]];t=rm(u,v,a.c[a.g[b][1]],a.d[a.g[b][1]]);for(o=2;o<a.g[b].length;o++){r=a.g[b][o];k=$wnd.Math.sqrt((a.c[r]-u)*(a.c[r]-u)+(a.d[r]-v)*(a.d[r]-v));c=2*t-rm(u,v,a.c[r],a.d[r]);a.c[r]=u+k*$wnd.Math.sin(c);a.d[r]=v+k*$wnd.Math.cos(c);}}function Do(a,b,c){var d,e,f,g,h,i,j,k,l,m,n,o,p;d=true;i=0;p=0;m=a.i;a.c[b]=m;h=Ai(a.b,b);g=pi(a.b,b);e=ji(a.b,b);f=ti(a.b,b);k=bl(a.b,b);e==0&&f==0&&Co(h)&&(d=false);a.e[m]='';if(c!=-1){switch(Mi(a.b,c)){case 0:a.e[m]+='~';break;case 2:a.e[m]+='=';break;case 3:a.e[m]+='#';}}d&&(a.e[m]+='[');f!=0&&(a.e[m]+=''+f);a.e[m]+=''+g;if(d){if(0<(o=ml(a.b,b))){a.e[m]+='H';1<o&&(a.e[m]+=o);}}if(e!=0){e>0?a.e[m]+='+':a.e[m]+='-';(e<0?-e:e)>1&&(a.e[m]+=''+(e<0?-e:e));}d&&(a.e[m]+=']');c!=-1&&(a.j[c]=true);a.g[b]=true;++a.i;for(n=0;n<k;++n)a.j[cl(a.b,b,n)]||++i;for(n=0;n<k;++n){j=al(a.b,b,n);l=cl(a.b,b,n);if(a.j[l]){++p;continue;}if(a.g[j]){++a.d;a.j[l]=true;switch(Mi(a.b,l)){case 0:a.e[a.c[j]]+='~';a.e[m]+='~';break;case 2:a.e[a.c[j]]+='=';a.e[m]+='=';break;case 3:a.e[a.c[j]]+='#';a.e[m]+='3';}if(a.d>9){a.e[a.c[j]]+='%';a.e[m]+='%';}a.e[a.c[j]]+=''+a.d;a.e[m]+=''+a.d;continue;}n-p<i-1&&(a.e[a.i++]='(');Do(a,j,l);n-p<i-1&&(a.e[a.i++]=')');}}function ed(a){var b,c,d,e,f,g,h,i,j,k,l;a.n=ZB(BE,iQ,40,a.G.o,0,1);for(h=0;h<a.G.p;h++)(Pi(a.G,h)==2||Pi(a.G,h)==26||Pi(a.G,h)==64)&&hd(a,h);for(i=0;i<a.G.p;i++)Pi(a.G,i)!=2&&Pi(a.G,i)!=26&&Pi(a.G,i)!=64&&hd(a,i);if((a.B&64)==0){for(g=0;g<a.G.p;g++){if(Hi(a.G,g)!=0){switch(Hi(a.G,g)){case 1:d=Mi(a.G,g)==2?'E':hj(a.G,g)?'p':'P';break;case 2:d=Mi(a.G,g)==2?'Z':hj(a.G,g)?'m':'M';break;default:d='?';}vo(a,(a.Q*2+1)/3|0);zd(a,fj(a.G,g)?-3:448);b=Ei(a.G,0,g);c=Ei(a.G,1,g);k=(uh(a.K,xi(a.G,b))+uh(a.K,xi(a.G,c)))/2;l=(vh(a.K,yi(a.G,b))+vh(a.K,yi(a.G,c)))/2;e=(uh(a.K,xi(a.G,b))-uh(a.K,xi(a.G,c)))/3;f=(vh(a.K,yi(a.G,b))-vh(a.K,yi(a.G,c)))/3;ld(a,k+f,l-e,d,true,true);zd(a,a.J);vo(a,a.Q);}}}if((a.B&4)!=0){vo(a,(a.Q*2+1)/3|0);zd(a,384);for(g=0;g<a.G.p;g++){b=Ei(a.G,0,g);c=Ei(a.G,1,g);j=Hl(a.G,g)?'d':Fl(a.G,g)?'a':'';k=(uh(a.K,xi(a.G,b))+uh(a.K,xi(a.G,c)))/2;l=(vh(a.K,yi(a.G,b))+vh(a.K,yi(a.G,c)))/2;ld(a,k,l,j+(''+g),true,true);}zd(a,a.J);vo(a,a.Q);}}function xn(a,b,c){var d,e,f,g,h,i,j,k,l,m;i=bl(a.A,b);e=a.i[c];if(e>i)return false;k=vi(a.A,b);g=vi(a.d,c);f=qi(a.d,c);j=qi(a.A,b);if((g&1)!=0){if(f!=null){if((k&1)!=0){if(j==null)return false;if(!Ln(f,j))return false;}else{if(j!=null){if(Nn(j,f))return false;}else{if(Kn(Ai(a.A,b),f))return false;}}}}else{if((k&1)!=0)return false;if(f!=null){if(j!=null){if(!Ln(j,f))return false;}else{if(!Kn(Ai(a.A,b),f))return false;}}else{if(j!=null)return false;if(a.C[b]!==a.f[c])return false;}}if((k|g)!=0){if((g&YQ)!=0){if(a.A.I&&(k&YQ)==0)return false;else if(e!=i)return false;}if((g&xQ)!=0){if(e>=i&&(k&xQ)==0)return false;}}if((a.B[b]&~a.e[c])!=0)return false;if(ji(a.d,c)!=0&&ji(a.d,c)!=ji(a.A,b))return false;if(ti(a.d,c)!=0&&ti(a.d,c)!=ti(a.A,b))return false;m=(vi(a.d,c)&WQ)>>22;if(m!=0){if(a.A.I&&m==(vi(a.A,c)&WQ)>>22)return true;d=false;l=tl(a.A);for(h=0;h<l.g.a.length;h++){if(BM(l.i,h).length==m){if(pn(l,h,b)){d=true;break;}}}if(!d)return false;}return true;}function KG(){var a,b,c;b=$doc.compatMode;a=aC(VB(kF,1),vR,2,6,[$T]);for(c=0;c<a.length;c++){if(DJ(a[c],b)){return;}}a.length==1&&DJ($T,a[0])&&DJ('BackCompat',b)?"GWT no longer supports Quirks Mode (document.compatMode=' BackCompat').<br>Make sure your application's host HTML page has a Standards Mode (document.compatMode=' CSS1Compat') doctype,<br>e.g. by using &lt;!doctype html&gt; at the start of your application's HTML page.<br><br>To continue using this unsupported rendering mode and risk layout problems, suppress this message by adding<br>the following line to your*.gwt.xml module file:<br>&nbsp;&nbsp;&lt;extend-configuration-property name=\"document.compatMode\" value=\""+b+'"/&gt;':"Your *.gwt.xml module configuration prohibits the use of the current document rendering mode (document.compatMode=' "+b+"').<br>Modify your application's host HTML page doctype, or update your custom "+"'document.compatMode' configuration property settings.";}function fn(a,b,c,d,e,f,g){var h,i,j,k,l,m,n,o,p,q,r,s,t;q=BM(a.g,b);r=BM(a.i,b);s=r.length;j=0;i=0;t=false;for(o=0;o<s;o++){j<<=1;i<<=1;if(Mi(a.f,r[o])>1||Pi(a.f,r[o])==64){j|=1;}else{h=c[b][o];if(h!=-1){if(d[h]){if(e[h]){j|=1;f[h]||(i|=1);}}else{t=true;}}}}n=false;switch(s){case 5:k=aC(VB(_C,1),DQ,6,15,[10,5,18,9,20]);n=true;for(p=0;p<5;p++){if((j&k[p])==k[p]){switch(Ai(a.f,q[p])){case 6:if(ji(a.f,q[p])==-1){e[b]=true;g[b]=p;(i&k[p])==0&&(n=false);}break;case 7:if(ji(a.f,q[p])<=0){e[b]=true;g[b]=p;}break;case 8:e[b]=true;g[b]=p;break;case 16:if(bl(a.f,q[p])==2){e[b]=true;g[b]=p;}}}}break;case 6:n=true;if((j&21)==21){e[b]=true;(i&21)==0&&(n=false);}if((j&42)==42){e[b]=true;(i&42)==0&&(n=false);}break;case 7:l=aC(VB(_C,1),DQ,6,15,[42,21,74,37,82,41,84]);n=true;for(m=0;m<7;m++){if((j&l[m])==l[m]){if(Ai(a.f,q[m])==6&&ji(a.f,q[m])==1||Ai(a.f,q[m])==5&&ji(a.f,q[m])==0){e[b]=true;g[b]=m;(i&l[m])==0&&(n=false);}}}}e[b]&&!n&&(f[b]=true);if(e[b])return true;return!t;}function MB(a,b,c,d){var e,f,g,h,i,j,k,l,m,n,o,p;f=-1;g=0;p=0;h=0;j=-1;k=b.length;n=c;l=true;for(;n<k&&l;++n){e=b.charCodeAt(n);switch(e){case 35:p>0?++h:++g;j>=0&&f<0&&++j;break;case 48:if(h>0){throw dG(new NI("Unexpected '0' in pattern \""+b+'"'));}++p;j>=0&&f<0&&++j;break;case 44:j=0;break;case 46:if(f>=0){throw dG(new NI('Multiple decimal separators in pattern "'+b+'"'));}f=g+p+h;break;case 69:if(!d){if(a.v){throw dG(new NI('Multiple exponential symbols in pattern "'+b+'"'));}a.v=true;a.k=0;}while(n+1<k&&b.charCodeAt(n+1)==48){++n;d||++a.k;}if(!d&&g+p<1||a.k<1){throw dG(new NI('Malformed exponential pattern "'+b+'"'));}l=false;break;default:--n;l=false;}}if(p==0&&g>0&&f>=0){m=f;f==0&&++m;h=g-m;g=m-1;p=1;}if(f<0&&h>0||f>=0&&(f<g||f>g+p)||j==0){throw dG(new NI('Malformed pattern "'+b+'"'));}if(d){return n-c;}o=g+p+h;a.i=f>=0?o-f:0;if(f>=0){a.n=g+p-f;a.n<0&&(a.n=0);}i=f>=0?f:o;a.o=i-g;if(a.v){a.j=g+a.o;a.i==0&&a.o==0&&(a.o=1);}a.f=j>0?j:0;a.c=f==0||f==o;return n-c;}function dh(a,b){var c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v;o=a.c[b];s=a.d[b];f=ZB(yD,fR,20,4,0,1);k=0;for(l=0;l<a.e.length;l++){if(k>=4)break;if(b==a.b[Ei(a.p,0,a.e[l])]||b==a.b[Ei(a.p,1,a.e[l])])continue;p=a.c[a.b[Ei(a.p,0,a.e[l])]];t=a.d[a.b[Ei(a.p,0,a.e[l])]];q=a.c[a.b[Ei(a.p,1,a.e[l])]];u=a.d[a.b[Ei(a.p,1,a.e[l])]];h=$wnd.Math.sqrt((p-o)*(p-o)+(t-s)*(t-s));i=$wnd.Math.sqrt((q-o)*(q-o)+(u-s)*(u-s));e=$wnd.Math.sqrt((q-p)*(q-p)+(u-t)*(u-t));if(h<e&&i<e){if(p==q){g=$wnd.Math.abs(o-p);g<0.5&&(f[k++]=new pm(rm(p,s,o,s),(0.5-g)/2));}else if(t==u){g=$wnd.Math.abs(s-t);g<0.5&&(f[k++]=new pm(rm(o,t,o,s),(0.5-g)/2));}else{m=(u-t)/(q-p);n=-1/m;c=t-m*p;d=s-n*o;r=(d-c)/(m-n);v=m*r+c;g=$wnd.Math.sqrt((r-o)*(r-o)+(v-s)*(v-s));g<0.5&&(f[k++]=new pm(rm(r,v,o,s),(0.5-g)/2));}continue;}if(h<0.5){f[k++]=new pm(rm(p,t,o,s),(0.5-h)/2);continue;}if(i<0.5){f[k++]=new pm(rm(q,u,o,s),(0.5-i)/2);continue;}}if(k>0){j=Ug(f,k);a.c[b]+=j.b*$wnd.Math.sin(j.a);a.d[b]+=j.b*$wnd.Math.cos(j.a);}}function je(a,b,c){var d,e,f,g,h,i,j,k,l,m;m=aC(VB(_C,2),mR,7,0,[aC(VB(_C,1),DQ,6,15,[2,1,2,1]),aC(VB(_C,1),DQ,6,15,[1,2,2,1]),aC(VB(_C,1),DQ,6,15,[1,1,2,2]),aC(VB(_C,1),DQ,6,15,[2,1,1,2]),aC(VB(_C,1),DQ,6,15,[2,2,1,1]),aC(VB(_C,1),DQ,6,15,[1,2,1,2])]);d=ZB(ZC,GQ,6,Qk(a.L,b),15,1);for(g=0;g<Qk(a.L,b);g++)d[g]=Di(a.L,al(a.L,b,c[g]),b);j=gl(a.L,b,c,d,null)<<24>>24;if(j!=3)return j;k=0;l=0;for(h=0;h<Qk(a.L,b);h++){e=cl(a.L,b,c[h]);if(Ei(a.L,0,e)==b){if(Pi(a.L,e)==9){l!=0&&qk(a.L,b);k=h;l=1;}if(Pi(a.L,e)==17){l!=0&&qk(a.L,b);k=h;l=2;}}}if(l==0)return 3;for(f=1;f<Qk(a.L,b);f++)d[f]<d[0]&&(d[f]+=LQ);if(Qk(a.L,b)==3){switch(k){case 0:(d[1]<d[2]&&d[2]-d[1]<MQ||d[1]>d[2]&&d[1]-d[2]>MQ)&&(l=3-l);break;case 1:d[2]-d[0]>MQ&&(l=3-l);break;case 2:d[1]-d[0]<MQ&&(l=3-l);}return l==1?2:1;}i=0;d[1]<=d[2]&&d[2]<=d[3]?i=0:d[1]<=d[3]&&d[3]<=d[2]?i=1:d[2]<=d[1]&&d[1]<=d[3]?i=2:d[2]<=d[3]&&d[3]<=d[1]?i=3:d[3]<=d[1]&&d[1]<=d[2]?i=4:d[3]<=d[2]&&d[2]<=d[1]&&(i=5);return m[i][k]==l?2:1;}function wj(a,b,c,d,e,f,g){var h,i,j,k,l,m,n,o,p,q,r,s,t,u,v;if(b==d){n=Ci(a,a.o,a.p,Fh);}else{u=a.H[b].a-a.H[d].a;v=a.H[b].b-a.H[d].b;n=$wnd.Math.sqrt(u*u+v*v);}h=b;o=rk(a,b)!=3;for(t=1;t<c;t++){q=a.H[h].a+n*$wnd.Math.sin(f);r=a.H[h].b+n*$wnd.Math.cos(f);s=-1;for(p=0;p<a.o;p++){if($wnd.Math.abs(q-a.H[p].a)<4&&$wnd.Math.abs(r-a.H[p].b)<4){s=p;break;}}if(s==-1){s=Hh(a,q,r,0);a.H[s].a=q;a.H[s].b=r;a.H[s].c=0;}m=Li(a,h,s);if(m==-1){m=Jh(a,h,s,(j=a.A[h],j>=3&&j<=4||j>=11&&j<=13||j>=19&&j<=31||j>=37&&j<=51||j>=55&&j<=84||j>=87&&j<=103||(k=a.A[s],k>=3&&k<=4||k>=11&&k<=13||k>=19&&k<=31||k>=37&&k<=51||k>=55&&k<=84||k>=87&&k<=103)?32:1));if(e){o&&rk(a,a.B[0][m])<4&&rk(a,a.B[1][m])<3&&(a.F[m]=2);o=!o;}}h=s;f+=g;}m=Li(a,h,d);m==-1&&(m=Jh(a,h,d,(l=a.A[h],l>=3&&l<=4||l>=11&&l<=13||l>=19&&l<=31||l>=37&&l<=51||l>=55&&l<=84||l>=87&&l<=103||(i=a.A[d],i>=3&&i<=4||i>=11&&i<=13||i>=19&&i<=31||i>=37&&i<=51||i>=55&&i<=84||i>=87&&i<=103)?32:1)));e&&o&&rk(a,a.B[0][m])<4&&rk(a,a.B[1][m])<4&&(a.F[m]=2);}function Ue(a){var b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v;if(a.L.d==0)return;if(a.A)return;a.C=0;v=0;for(c=1;c<a.L.d;c++)a.c[c]>a.c[v]&&(v=c);d=ZB(aG,HQ,6,a.L.d,16,1);g=ZB(aG,HQ,6,a.L.e,16,1);a.B=ZB(_C,DQ,6,a.L.d,15,1);a.t=ZB(_C,DQ,6,a.L.d,15,1);a.w=ZB(_C,DQ,6,a.L.d,15,1);a.u=ZB(_C,DQ,6,a.L.e,15,1);a.t[0]=v;a.B[v]=0;d[v]=true;e=1;i=0;j=1;k=0;while(i<a.L.d){if(i<j){while(true){o=0;p=0;m=-1;for(q=0;q<qe(a,a.t[i]);q++){h=al(a.L,a.t[i],q);if(!d[h]&&a.c[h]>m){o=h;p=cl(a.L,a.t[i],q);m=a.c[h];}}if(m==-1)break;a.B[o]=j;a.w[j]=i;a.t[j++]=o;a.u[k++]=p;d[o]=true;g[p]=true;}++i;}else{n=0;m=-1;for(b=0;b<a.L.d;b++){if(!d[b]&&a.c[b]>m){n=b;m=a.c[b];}}++e;a.B[n]=j;a.w[j]=-1;a.t[j++]=n;d[n]=true;}}a.v=ZB(_C,DQ,6,2*(a.L.e-k),15,1);while(true){s=a.L.K;t=a.L.K;u=-1;for(f=0;f<a.L.e;f++){if(!g[f]){if(a.B[Ei(a.L,0,f)]<a.B[Ei(a.L,1,f)]){r=a.B[Ei(a.L,0,f)];l=a.B[Ei(a.L,1,f)];}else{r=a.B[Ei(a.L,1,f)];l=a.B[Ei(a.L,0,f)];}if(r<s||r==s&&l<t){s=r;t=l;u=f;}}}if(u==-1)break;g[u]=true;a.u[k++]=u;a.v[2*a.C]=s;a.v[2*a.C+1]=t;++a.C;}a.A=true;}function vn(a,b){var c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r;this.f=a;this.g=new LM();this.i=new LM();this.a=ZB(_C,DQ,6,this.f.d,15,1);this.b=ZB(_C,DQ,6,this.f.e,15,1);this.f.gb(1);m=ZB(aG,HQ,6,this.f.d,16,1);n=ZB(aG,HQ,6,this.f.e,16,1);do{g=false;for(c=0;c<this.f.d;c++){if(!m[c]){p=0;for(l=0;l<bl(this.f,c);l++)m[al(this.f,c,l)]||++p;if(p<2){m[c]=true;for(k=0;k<bl(this.f,c);k++)n[cl(this.f,c,k)]=true;g=true;}}}}while(g);r=0;while(r<this.f.d&&m[r])++r;if(r==this.f.d)return;i=ZB(_C,DQ,6,this.f.d,15,1);i[0]=r;h=ZB(_C,DQ,6,this.f.d,15,1);h[r]=1;f=0;j=0;o=1;while(f<=j){for(k=0;k<bl(this.f,i[f]);k++){e=al(this.f,i[f],k);if(h[e]!=0){en(this,cl(this.f,i[f],k),m);continue;}if(!m[e]){h[e]=o;i[++j]=e;}}++f;if(f>j){for(c=0;c<this.f.d;c++){if(h[c]==0&&!m[c]){h[c]=++o;i[++j]=c;break;}}}}if((b&4)!=0){this.d=ZB(aG,HQ,6,this.g.a.length,16,1);this.e=ZB(aG,HQ,6,this.g.a.length,16,1);this.c=ZB(_C,DQ,6,this.g.a.length,15,1);gn(this,this.d,this.e,this.c);}if((b&2)!=0){for(d=0;d<this.f.e;d++){if(!n[d]){q=hn(this,d,m);q!=null&&tn(this,q,nn(this,q));}}}}function Fn(a,b){var c,d,e,f,g,h,i,j,k,l;a.v=new LM();WN(a.H.a);WN(a.c.a);if(!a.A||!a.d)return 0;if(a.d.d-a.j>a.A.d||a.d.e-a.k>a.A.e)return 0;if(a.d.d-a.j==0)return 0;Qn(a,b);c=ZB(aG,HQ,6,a.A.d,16,1);a.w=ZB(_C,DQ,6,a.d.d,15,1);lN(a.w);g=ZB(_C,DQ,6,a.t,15,1);oN(g,g.length,-1);e=0;while(true){++g[e];j=a.q[e]==-1?a.A.d:bl(a.A,a.w[a.q[e]]);if(g[e]==j){g[e]=-1;if(e==0)break;--e;a.p[e]||(c[a.w[a.o[e]]]=false);continue;}if(a.q[e]==-1){if(!c[g[e]]){if(xn(a,g[e],a.o[e])){a.w[a.o[e]]=g[e];c[g[e]]=true;++e;}}}else if(a.p[e]){d=al(a.A,a.w[a.q[e]],g[e]);d==a.w[a.o[e]]&&yn(a,cl(a.A,a.w[a.q[e]],g[e]),a.r[e])&&++e;}else{d=al(a.A,a.w[a.q[e]],g[e]);if(!c[d]){if(xn(a,d,a.o[e])&&yn(a,cl(a.A,a.w[a.q[e]],g[e]),a.r[e])){c[d]=true;a.w[a.o[e]]=d;++e;}}}if(e==a.s){if(En(a,false)&&Cn(a,false)&&Bn(a,c,false)){if(a.j==0)return 1;h=false;if(a.j!=0){k=hN(a.w,a.w.length);xN(k);if(RO(a.c,k)){h=true;}else if(Dn(a,c,g)){QO(a.c,k);l=ZB(_C,DQ,6,k.length,15,1);for(f=a.v.a.length-1;f>=0;f--){i=BM(a.v,f);dK(i,l,l.length);xN(l);mA(l,k)==0&&FM(a.v,f);}h=true;}}h||wn(a);}--e;a.p[e]||(c[a.w[a.o[e]]]=false);}}return a.v.a.length;}function ie(a,b,c){var d,e,f,g,h,i,j,k,l,m,n,o;if(a.W[b]!=0)return false;if(Ai(a.L,b)!=6&&Ai(a.L,b)!=7&&Ai(a.L,b)!=14&&Ai(a.L,b)!=15&&Ai(a.L,b)!=16)return false;if(Tk(a.L,b)!=0){if(bl(a.L,b)==2&&dl(a.L,b,0)==2&&dl(a.L,b,1)==2)return ce(a,b,c);if(Ai(a.L,b)!=15&&Ai(a.L,b)!=16)return false;}if(bl(a.L,b)<3||Qk(a.L,b)>4)return false;if(Ai(a.L,b)==7&&!a.M[b])return false;n=ZB(_C,DQ,6,4,15,1);o=ZB(_C,DQ,6,4,15,1);j=ZB(aG,HQ,6,4,16,1);for(h=0;h<Qk(a.L,b);h++){f=-1;e=0;for(i=0;i<Qk(a.L,b);i++){if(!j[i]){if(f<a.c[al(a.L,b,i)]){f=a.c[al(a.L,b,i)];e=i;}}}n[h]=e;o[h]=f;j[e]=true;}if(Qk(a.L,b)==4&&o[0]===o[1]&&o[2]===o[3])return false;if(Qk(a.L,b)==4&&(o[0]===o[2]||o[1]===o[3]))return false;if(Qk(a.L,b)==3&&o[0]===o[2])return false;k=0;l=0;m=false;for(g=1;g<Qk(a.L,b);g++){if(o[g-1]===o[g]){if(!c||o[g]==0)return false;k=al(a.L,b,n[g-1]);l=al(a.L,b,n[g]);Ll(a.L,cl(a.L,b,n[g]))&&(a.P[b]=true);m=true;}}if(c&&!m)return false;d=a._?ke(a,b,n):je(a,b,n);if(c){if(a.Q&&(a.K&2)!=0||!a.Q&&(a.K&4)!=0){if(d==1){Df(a.b[k],OQ);Df(a.b[l],256);}else if(d==2){Df(a.b[k],256);Df(a.b[l],OQ);}}}else{a.W[b]=d;}return true;}function $d(){$d=FG;Yd=aC(VB(_F,1),gR,6,15,[-1,-1,-1,0,0,1,2,3,4,5,-1,0,0,0,6,7,8,9,-1,0,0,10,10,10,10,10,10,10,10,10,10,1,11,11,12,13,-1,0,0,10,10,10,10,10,10,10,10,10,10,0,0,0,11,14,-1,0,0,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,10,10,10,10,10,10,10,10,1,1,1,1,-1,-1,-1,-1,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1]);Zd=aC(VB(_F,1),gR,6,15,[-1,-1,-1,0,0,0,2,5,5,5,-1,0,0,0,0,9,9,9,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,9,-1,0,0,0,0,0,0,10,0,0,0,0,0,0,0,0,0,9,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1]);}function Ho(a){var b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w;Xo(a.b,1);a.a=ZB(aG,HQ,6,a.b.e,16,1);for(g=0;g<a.b.e;g++){if(Pi(a.b,g)==64){jk(a.b,g,1);a.a[g]=true;}}w=new vn(a.b,3);m=ZB(aG,HQ,6,w.g.a.length,16,1);for(s=0;s<w.g.a.length;s++){u=BM(w.g,s);m[s]=true;for(l=0;l<u.length;l++){if(!lj(a.b,u[l])){m[s]=false;break;}}if(m[s]){v=BM(w.i,s);for(k=0;k<v.length;k++)a.a[v[k]]=true;}}for(h=0;h<a.b.e;h++){!a.a[h]&&w.b[h]!=0&&lj(a.b,Ei(a.b,0,h))&&lj(a.b,Ei(a.b,1,h))&&Fo(a,h);}Xo(a.b,3);for(t=0;t<w.g.a.length;t++){if(m[t]){u=BM(w.g,t);for(k=0;k<u.length;k++){if(!Lo(a,u[k])){Sj(a.b,u[k],false);for(o=0;o<bl(a.b,u[k]);o++)a.a[cl(a.b,u[k],o)]=false;}}}}Ko(a);for(r=0;r<w.g.a.length;r++){if(m[r]&&BM(w.i,r).length==6){v=BM(w.i,r);n=true;for(e=0,f=v.length;e<f;++e){i=v[e];if(!a.a[i]){n=false;break;}}if(n){Jo(a,v[0]);Jo(a,v[2]);Jo(a,v[4]);Ko(a);}}}for(q=5;q>=4;q--){do{p=false;for(i=0;i<a.b.e;i++){if(a.a[i]){b=0;for(k=0;k<2;k++){j=Ei(a.b,k,i);for(o=0;o<bl(a.b,j);o++)a.a[cl(a.b,j,o)]&&++b;}if(b==q){Jo(a,i);Ko(a);p=true;break;}}}}while(p);}for(d=0;d<a.b.e;d++)if(a.a[d])throw dG(new FA(RR));for(c=0;c<a.b.d;c++)if(lj(a.b,c))throw dG(new FA(RR));}function dd(a,b,c,d){var e,f,g,h,i,j,k,l,m,n,o,p;e=new Fd();i=new Fd();k=new hH();j=new hH();g=Ei(a.G,0,c);h=Ei(a.G,1,c);if(d){m=b.a;b.a=b.b;b.b=m;m=b.c;b.c=b.d;b.d=m;n=g;g=h;h=n;}if(!vd(a,b))return;if(Ll(a.G,c)){e.a=b.a;e.c=b.c;e.b=b.b;e.d=b.d;l=d?-ud(a,c):ud(a,c);l==0&&(l=1);cd(a,b.b-b.a,b.d-b.c,k);if(l>0){i.a=b.a+k.a;i.c=b.c+k.b;i.b=b.b+k.a;i.d=b.d+k.b;if(bd(a,g,h,1,j)||bl(a.G,g)>1){i.a+=j.a+k.b;i.c+=j.b-k.a;}}else{i.a=b.a-k.a;i.c=b.c-k.b;i.b=b.b-k.a;i.d=b.d-k.b;if(bd(a,g,h,-1,j)||bl(a.G,g)>1){i.a+=j.a+k.b;i.c+=j.b-k.a;}}Pi(a.G,c)==26&&td(e,i);vd(a,e)&&Pc(a,e,g,h);Pi(a.G,c)==64?vd(a,i)&&Nc(a,i,g,h):vd(a,i)&&Pc(a,i,g,h);}else{cd(a,b.b-b.a,b.d-b.c,k);o=k.a/2;p=k.b/2;f=false;e.a=b.a+o;e.c=b.c+p;e.b=b.b+o;e.d=b.d+p;if(bl(a.G,g)>1){if(bd(a,g,h,1,j)){e.a+=j.a;e.c+=j.b;if(bl(a.G,g)==2){if(j.a!=0||j.b!=0){e.a+=k.b;e.c-=k.a;}}}else{a.n[g]=new iH(e.a,e.c);}}i.a=b.a-o;i.c=b.c-p;i.b=b.b-o;i.d=b.d-p;if(bl(a.G,g)>1){if(bd(a,g,h,0,j)){i.a+=j.a;i.c+=j.b;if(bl(a.G,g)==2){if(j.a!=0||j.b!=0){i.a+=k.b;i.c-=k.a;}}}else{a.n[g]=new iH(i.a,i.c);f=true;}}Pi(a.G,c)==26&&td(e,i);if(Pi(a.G,c)==64){if(f){Nc(a,e,g,h);Pc(a,i,g,h);}else{Pc(a,e,g,h);Nc(a,i,g,h);}}else{Pc(a,e,g,h);Pc(a,i,g,h);}}}function Cg(a,b,c,d,e){var f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,A,B,C,D,F,G,H,I,J,K,L;t=ZB(_C,DQ,6,e,15,1);u=ZB(_C,DQ,6,e,15,1);for(p=0;p<e;p++){t[p]=_g(b,d[p]);u[p]=_g(c,d[p]);}F=0;H=0;G=0;I=0;for(q=0;q<e;q++){F+=b.c[t[q]];H+=b.d[t[q]];G+=c.c[u[q]];I+=c.d[u[q]];}F/=e;H/=e;G/=e;I/=e;fh(c,F-G,H-I);j=ZB(yD,fR,20,e,0,1);l=ZB(yD,fR,20,e,0,1);f=ZB(yD,fR,20,e,0,1);g=ZB(yD,fR,20,e,0,1);for(r=0;r<e;r++){j[r]=new qm(F,H,b.c[t[r]],b.d[t[r]]);l[r]=new qm(F,H,c.c[u[r]],c.d[u[r]]);f[r]=new pm(j[r].a-l[r].a,j[r].b*l[r].b);g[r]=new pm(j[r].a+l[r].a,j[r].b*l[r].b);}w=Ug(f,e);A=Ug(g,e);K=0;L=0;for(s=0;s<e;s++){for(v=0;v<a.e[d[s]];v++){h=al(a.i,d[s],v);ah(b,h)&&!ah(c,h)&&++K;!ah(b,h)&&ah(c,h)&&++L;}}k=ZB(yD,fR,20,K,0,1);m=ZB(yD,fR,20,L,0,1);n=ZB(yD,fR,20,L,0,1);K=0;L=0;for(o=0;o<e;o++){for(v=0;v<a.e[d[o]];v++){h=al(a.i,d[o],v);if(ah(b,h)&&!ah(c,h)){i=_g(b,h);k[K]=new qm(b.c[t[o]],b.d[t[o]],b.c[i],b.d[i]);++K;}if(!ah(b,h)&&ah(c,h)){i=_g(c,h);J=new qm(c.c[u[o]],c.d[u[o]],c.c[i],c.d[i]);m[L]=new pm(w.a+J.a,J.b);n[L]=new pm(A.a-J.a,J.b);++L;}}}B=Ug(k,K);C=Ug(m,L);D=Ug(n,L);if($wnd.Math.abs(Ag(B.a,C.a))>$wnd.Math.abs(Ag(B.a,D.a))){eh(c,F,H,w.a);}else{Xg(c,F,H);eh(c,F,H,A.a);}return Eg(a,b,c,e);}function _d(a,b,c){$d();var d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t;Xo(a,3);p=ZB(aD,gR,6,a.g[b],14,1);for(l=0;l<a.g[b];l++){g=0;i=0;if((c&32)!=0){h=a.j[b][l];gG(h,3)<0&&Fl(a,a.i[b][l])&&(h=0);i=eG(i,h);}f=a.f[b][l];if((c&128)!=0){if(Yd[a.A[f]]==-1)throw dG(new FA(hR+a.A[f]));g=eG(g,Yd[a.A[f]]);}else if((c&64)!=0){if(Zd[a.A[f]]==-1)throw dG(new FA(hR+a.A[f]));g=eG(g,Zd[a.A[f]]);}if((c&256)!=0){r=a.g[f]-1;r>3&&(r=3);(c&512)==0&&r>1&&(r=1);g=eG(g,r<<4);}(c&OQ)!=0&&(a.s[f]&8)!=0&&(g=eG(g,64));(c&YQ)!=0&&(a.s[f]&xQ)!=0&&(g=eG(g,128));t=eG(g,pG(i,8));n=0;while(gG(t,p[n])<0)++n;for(o=l;o>n;o--)p[o]=p[o-1];p[n]=t;}q=a.g[b]<4?a.g[b]:4;e=0;for(m=0;m<q;m++){e=pG(e,10);e=eG(e,p[m]);}e=pG(e,10);if(Yd[a.A[b]]==-1)throw dG(new FA(hR+a.A[b]));e=oG(e,Yd[a.A[b]]);if((c&2)!=0){s=!!a.n&&b<a.d?jn(a.n,b):0;s>9&&(s=9);s>2&&(s-=2);e=oG(e,s<<4);}else(c&1)!=0&&(a.s[b]&8)!=0&&(e=oG(e,64));(c&4)!=0&&(a.s[b]&xQ)!=0&&(e=eG(e,128));(c&8)!=0&&(a.s[b]&iR)!=0&&(e=eG(e,256));(c&16)!=0&&(a.s[b]&yQ)!=0&&(e=eG(e,512));if(nG(fG(e,jR),0)){j=new HA(kR);vA(j,(cK(),bK),'');}if(nG(fG(e,lR),0)){j=new HA(kR);vA(j,(cK(),bK),'');}if((c&xQ)!=0){Sd(a,b)&&(e=eG(e,jR));d=false;if(Vd(a,b)){for(k=0;k<a.d;k++){if(Td(a,k)){d=true;break;}}}d&&(e=eG(e,lR));}return e;}function ne(a){var b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v;a.M=ZB(aG,HQ,6,a.L.d,16,1);for(b=0;b<a.L.d;b++){if(Ai(a.L,b)==7){if(bl(a.L,b)==4){a.M[b]=true;continue;}if(bl(a.L,b)==3){if(ji(a.L,b)==1){a.M[b]=true;continue;}if(Il(a.L,b))continue;if((a.K&32)!=0){a.M[b]=true;continue;}if(Vk(a.L,b)!=3)continue;v=Xk(a.L,b);if(v>7)continue;t=tl(a.L);u=0;while(u<t.g.a.length){if(BM(t.i,u).length==v&&pn(t,u,b))break;++u;}i=-1;j=-1;for(l=0;l<3;l++){h=cl(a.L,b,l);if(!qn(t,u,h)){i=al(a.L,b,l);j=h;break;}}n=ZB(aG,HQ,6,a.L.e,16,1);n[j]=true;o=ZB(_C,DQ,6,11,15,1);p=pl(a.L,o,i,b,10,n);if(p==-1)continue;d=1;while(!pn(t,u,o[d]))++d;c=p-d;e=o[d];if(v==6&&c==2&&d==3){if(Vk(a.L,o[1])>=3){m=false;s=BM(t.g,u);for(k=0;k<6;k++){if(b==s[k]){r=un(t,u,e==s[un(t,u,k+2)]?k-2:k+2);q=s[r];Vk(a.L,q)>=3&&sl(a.L,o[1],q,2,null)==2&&(m=true);break;}}if(m){a.M[b]=true;continue;}}}f=Tk(a.L,e)==1||El(a.L,e)||Il(a.L,e);g=!f&&Ai(a.L,e)==7&&ji(a.L,e)!=1;if(c==1){!f&&!g&&v<=4&&d<=3&&(a.M[b]=true);continue;}switch(v){case 4:!f&&!g&&d<=4&&(a.M[b]=true);break;case 5:g?d<=3&&(a.M[b]=true):f||d<=4&&(a.M[b]=true);break;case 6:c==2?f?d<=4&&(a.M[b]=true):g||d<=3&&(a.M[b]=true):c==3&&(f?d<=6&&(a.M[b]=true):d<=4&&(a.M[b]=true));break;case 7:c==3&&d<=3&&(a.M[b]=true);}}}}}function Yl(a,b){var c,d,e,f,g,h,i,j,k,l,m,n,o,p;if((a.s[b]&3)==0||(a.s[b]&3)==3)return;if(a.k[b]==2&&a.g[b]==2){Vl(a,b);return;}if(a.g[b]<3||a.g[b]>4){Uj(a,b,0,false);return;}o=vl(a,b);c=ZB(ZC,GQ,6,a.c[b],15,1);for(g=0;g<a.c[b];g++)c[g]=Di(a,a.f[b][o[g]],b);for(h=0;h<a.c[b];h++)a.B[0][a.i[b][h]]==b&&Mi(a,a.i[b][h])==1&&(a.F[a.i[b][h]]=1);if(Wl(a,b,o,c))return;m=-1;for(i=0;i<a.c[b];i++){e=a.i[b][i];if((a.F[e]==17||a.F[e]==9)&&a.B[0][e]==b){a.F[a.i[b][i]]=1;m==-1?m=e:m=-2;}}m<0&&(m=Tl(a,b));if(a.B[0][m]!=b){a.B[1][m]=a.B[0][m];a.B[0][m]=b;}n=-1;for(j=0;j<a.c[b];j++){if(m==a.i[b][o[j]]){n=j;break;}}p=aC(VB(_C,2),mR,7,0,[aC(VB(_C,1),DQ,6,15,[2,1,2,1]),aC(VB(_C,1),DQ,6,15,[1,2,2,1]),aC(VB(_C,1),DQ,6,15,[1,1,2,2]),aC(VB(_C,1),DQ,6,15,[2,1,1,2]),aC(VB(_C,1),DQ,6,15,[2,2,1,1]),aC(VB(_C,1),DQ,6,15,[1,2,1,2])]);for(f=1;f<a.c[b];f++)c[f]<c[0]&&(c[f]+=LQ);if(a.c[b]==3){k=false;switch(n){case 0:k=c[1]<c[2]&&c[2]-c[1]<MQ||c[1]>c[2]&&c[1]-c[2]>MQ;break;case 1:k=c[2]-c[0]>MQ;break;case 2:k=c[1]-c[0]<MQ;}d=(a.s[b]&3)==1^k?17:9;}else{l=0;c[1]<=c[2]&&c[2]<=c[3]?l=0:c[1]<=c[3]&&c[3]<=c[2]?l=1:c[2]<=c[1]&&c[1]<=c[3]?l=2:c[2]<=c[3]&&c[3]<=c[1]?l=3:c[3]<=c[1]&&c[1]<=c[2]?l=4:c[3]<=c[2]&&c[2]<=c[1]&&(l=5);d=(a.s[b]&3)==1^p[l][n]==1?9:17;}a.F[m]=d;}function ng(a){var b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q;this.j=a;for(d=0;d<a.i.d;d++){a.f[d]&&(a.o[d]==1||a.o[d]==2)&&(a.k[d]==1?this.a<=a.j[d]&&(this.a=1+a.j[d]):a.k[d]==2&&this.i<=a.j[d]&&(this.i=1+a.j[d]));}this.b=this.a+this.i;this.e=XB(aG,[iQ,HQ],[11,6],16,[this.b+1,a.g.length+1],2);for(e=0;e<a.i.d;e++)a.f[e]&&(a.o[e]==1||a.o[e]==2)&&!a.e[e]&&(this.e[jg(this,e)][a.g.length]=true);for(i=0;i<a.g.length;i++){for(q=0;q<a.g[i].length;q++){c=a.g[i][q];a.f[c]&&(a.o[c]==1||a.o[c]==2)&&(this.e[jg(this,c)][i]=true);}}this.d=ZB(_C,mR,7,this.b,0,2);for(j=0;j<a.g.length;j++){for(n=1;n<this.b;n++){if(this.e[n][j]){for(o=0;o<n;o++){if(this.e[o][j]){this.d[n]=$f(this.d[n],o);this.d[o]=$f(this.d[o],n);}}}}}this.c=ZB(_C,DQ,6,this.b+1,15,1);for(m=0;m<this.b;m++){this.e[m][a.g.length]?this.c[m]=-1:this.c[m]=-2;}for(k=0;k<a.g.length;k++){if(this.e[this.b][k]){for(l=0;l<this.b;l++){this.e[l][k]&&this.c[l]!=k&&(this.c[l]==-2?this.c[l]=k:this.c[l]=-3);}}}for(b=0;b<this.b;b++){if(this.c[b]>=-1){f=ZB(_C,DQ,6,this.b,15,1);if(eg(this,f,b)){for(l=0;l<this.b;l++){f[l]!=0&&(this.c[l]=-3);}}}}for(h=0;h<a.g.length-1;h++){for(n=1;n<this.b;n++){if(this.e[n][h]&&this.c[n]!=-3){for(o=0;o<n;o++){if(this.e[o][h]&&this.c[o]!=-3){g=fg(this,n,o,h);if(g!=null){for(p=0;p<g.length;p++)this.c[g[p]]=-3;mg(this,g);break;}}}}}}}function pe(a){var b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,A,B;s=ZB(aG,HQ,6,a.L.d,16,1);t=ZB(aG,HQ,6,a.L.e,16,1);b=0;v=false;for(d=0;d<a.L.d;d++){if(a.P[d]){if(!a.X[d]){if(ie(a,d,false)){a.X[d]=true;s[d]=true;++b;}}}}for(f=0;f<a.L.e;f++){if(a.O[f]){if(!a.n[f]){if(ee(a,f,false)){a.n[f]=true;t[f]=true;++b;}}}}if(b==1){for(c=0;c<a.L.d;c++){if(s[c]){a.W[c]=0;break;}}for(e=0;e<a.L.e;e++){if(t[e]){a.k[e]=0;break;}}}else if(b>1){me(a);for(h=new dN(a.s);h.a<h.c.a.length;){g=cN(h);u=0;w=0;k=0;j=0;l=-1;i=-1;for(o=0;o<g.a.length;o++){if(s[g.a[o]]){++u;if(a.W[g.a[o]]==1||a.W[g.a[o]]==2){++w;v=true;if(l<a.c[g.a[o]]){l=a.c[g.a[o]];k=g.a[o];}}}}for(p=0;p<g.b.length;p++){if(t[g.b[p]]){++u;A=a.c[Ei(a.L,0,g.b[p])];B=a.c[Ei(a.L,1,g.b[p])];m=A>B?(A<<16)+B:(B<<16)+A;if(a.k[g.b[p]]==1||a.k[g.b[p]]==2){++w;v=true;if(i<m){i=m;j=g.b[p];}}}}if(u==0)continue;if(u==1){for(q=0;q<g.a.length;q++)s[g.a[q]]&&(a.W[g.a[q]]=0);for(n=0;n<g.b.length;n++)t[g.b[n]]&&(a.k[g.b[n]]=0);}else{if(w==1){for(q=0;q<g.a.length;q++)s[g.a[q]]&&(a.W[g.a[q]]=3);for(n=0;n<g.b.length;n++)t[g.b[n]]&&(a.k[g.b[n]]=3);}else{r=false;l!=-1?a.W[k]==2&&(r=true):a.k[j]==2&&(r=true);if(r){for(q=0;q<g.a.length;q++){if(s[g.a[q]]){switch(a.W[g.a[q]]){case 1:a.W[g.a[q]]=2;break;case 2:a.W[g.a[q]]=1;}}}for(n=0;n<g.b.length;n++){if(t[g.b[n]]){switch(a.k[g.b[n]]){case 1:a.k[g.b[n]]=2;break;case 2:a.k[g.b[n]]=1;}}}}}}}}return v;}function Yf(a,b,c){var d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v;i=ZB(_C,DQ,6,a.i.d,15,1);s=ZB(_C,DQ,6,a.i.d,15,1);o=ZB(aG,HQ,6,a.i.d,16,1);p=ZB(aG,HQ,6,a.i.d,16,1);j=ZB(aG,HQ,6,a.i.d,16,1);i[0]=b;s[b]=c;s[c]=-2;o[b]=true;o[c]=true;f=0;k=0;while(f<=k){g=i[f];if(s[g]==g){for(l=0;l<bl(a.i,g);l++){d=al(a.i,g,l);if(!o[d]){if(dl(a.i,g,l)==2&&Ai(a.i,d)<10){i[++k]=d;s[d]=d;j[d]=j[g]||Tk(a.i,d)==2;p[d]=j[g]&&!p[g];o[d]=true;}else if(j[g]&&p[g]){t=Nf(a,d,s[g],o);if(t==-1)return null;i[++k]=d;s[d]=t;s[t]=-2;j[d]=false;o[d]=true;o[t]=true;}else if(Ll(a.i,cl(a.i,g,l))){i[++k]=d;s[d]=d;j[d]=false;o[d]=true;if((Ai(a.i,d)==6&&Tk(a.i,d)==0||Ai(a.i,d)==7&&ji(a.i,d)==1||Ai(a.i,d)==14||Ai(a.i,d)==15&&bl(a.i,d)>2||Ai(a.i,d)==16&&bl(a.i,d)>2)&&bl(a.i,d)>2){h=false;for(q=1;q<bl(a.i,d);q++){u=al(a.i,d,q);if(!o[u]){for(r=0;r<q;r++){v=al(a.i,d,r);if(!o[v]){if(Rf(a,u,v)){i[++k]=u;s[u]=v;s[v]=-2;j[u]=false;o[u]=true;o[v]=true;h=true;}}}}}if(!h)return null;}}}}}else{e=ZB(aG,HQ,6,bl(a.i,g),16,1);for(m=0;m<bl(a.i,g);m++){d=al(a.i,g,m);if(o[d]){e[m]=s[d]==d;}else{for(q=0;q<bl(a.i,d);q++){if(al(a.i,d,q)==s[g]){e[m]=true;break;}}}}for(n=0;n<bl(a.i,g);n++){if(e[n]){d=al(a.i,g,n);if(o[d]){if($k(a.i,d,s[g])==-1)return null;}else{i[++k]=d;s[d]=d;p[d]=false;j[d]=true;o[d]=true;}}}for(l=0;l<bl(a.i,g);l++){if(!e[l]){d=al(a.i,g,l);if(!o[d]){t=Nf(a,d,s[g],o);if(t==-1)return null;i[++k]=d;s[d]=t;s[t]=-2;j[d]=false;o[d]=true;o[t]=true;}}}}++f;}return o;}function Ic(){Ic=FG;wc=aC(VB(_C,1),DQ,6,15,[0,EQ,14286847,13402367,12779264,16758197,9474192,3166456,16715021,9494608,11789301,11230450,9109248,12560038,15780000,16744448,16777008,2093087,8442339,9388244,4062976,15132390,12567239,10921643,9083335,10255047,14706227,15765664,5296208,13140019,8224944,12750735,6721423,12419299,16752896,10889513,6076625,7351984,65280,9764863,9756896,7586505,5551541,3907230,2396047,687500,27013,12632256,16767375,10909043,6717568,10380213,13924864,9699476,4366000,5707663,51456,7394559,16777159,14286791,13107143,10747847,9437127,6422471,4587463,3211207,2097095,65436,58997,54354,48952,43812,5096191,5089023,2200790,2522539,2516630,1528967,13684960,16765219,12105936,10900557,5724513,10375093,11230208,7688005,4358806,4325478,32000,7384058,47871,41471,36863,33023,27647,5528818,7888099,9064419,10565332,11739092,11739066,11734438,12389767,13041766,13369433,13697103,14221381,14680120,15073326,15400998,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,13158600,1334015,56540,15075850,15132160,56540,15075850,15461355,8553170,1016335,1016335,1334015,15132160,3289770,14456450,16422400,16422400,11819700,3289770,1016335]);zc=new XG(255,128,0);yc=new XG(92,160,255);Hc=new XG(160,0,64);xc=new XG(255,160,255);Ac=new XG(32,96,255);Gc=new XG(255,0,0);Dc=new XG(0,255,0);Ec=new XG(192,0,255);Fc=new XG(255,160,0);Bc=new XG(0,128,0);Cc=new XG(160,0,0);}function rq(){rq=FG;pq=aC(VB(kF,1),vR,2,6,['[N](-*)(-*)-*','[N](-*)=*','[N]#*','[N](-*)(=*)=* as in nitro','[N](=*)#* middle atom of azide','[N]1(-*)-*-*-1 3-membered ring','[NH](-*)-*','[NH]1-*-*-1 3-membered ring','[NH]=*','[NH2]-*','[N+](-*)(-*)(-*)-*','[N+](-*)(-*)=*','[N+](-*)#* N in isocyano','[NH+](-*)(-*)-*','[NH+](-*)=*','[NH2+](-*)-*','[NH2+]=*','[NH3+]-*','[n](:*):*','[n](:*)(:*):*','[n](-*)(:*):*','[n](=*)(:*):* as in pyridine-N-oxid','[nH](:*):*','[n+](:*)(:*):*','[n+](-*)(:*):*','[nH+](:*):*','[O](-*)-*','[O]1-*-*-1 3-membered ring','[O]=*','[OH]-*','[O-]-*','[o](:*):*','[S](-*)-*','[S]=*','[S](-*)(-*)=*','[S](-*)(-*)(=*)=*','[SH]-*','[s](:*):*','[s](=*)(:*):*','[P](-*)(-*)-*','[P](-*)=*','[P](-*)(-*)(-*)=*','[PH](-*)(-*)=*']);qq=aC(VB($C,1),gR,6,15,[3.240000009536743,12.359999656677246,23.790000915527344,11.680000305175781,13.600000381469727,OT,12.029999732971191,21.940000534057617,23.850000381469727,26.020000457763672,0,OT,4.360000133514404,4.440000057220459,13.970000267028809,16.610000610351562,25.59000015258789,27.639999389648438,12.890000343322754,4.409999847412109,4.929999828338623,8.390000343322754,15.789999961853027,4.099999904632568,3.880000114440918,14.140000343322754,9.229999542236328,12.529999732971191,17.06999969482422,20.229999542236328,23.059999465942383,13.140000343322754,25.299999237060547,32.09000015258789,19.209999084472656,8.380000114440918,38.79999923706055,28.239999771118164,21.700000762939453,13.59000015258789,34.13999938964844,9.8100004196167,23.469999313354492]);}function ou(){ou=FG;nu=(!Wz&&(Wz=new cA()),Wz);Ar=aC(VB(kF,1),vR,2,6,['?','H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr',wR,wR,wR,wR,wR,wR,wR,wR,wR,wR,wR,wR,wR,wR,wR,wR,wR,wR,wR,wR,wR,wR,wR,wR,wR,'R4','R5','R6','R7','R8','R9','R10','R11','R12','R13','R14','R15','R16','R1','R2','R3','A','A1','A2','A3',wR,wR,'D','T','X','R','H2','H+','Nnn','HYD','Pol',wR,wR,wR,wR,wR,wR,wR,wR,wR,wR,wR,'Ala','Arg','Asn','Asp','Cys','Gln','Glu','Gly','His','Ile','Leu','Lys','Met','Phe','Pro','Ser','Thr','Trp','Tyr','Val']);mu=aC(VB(_F,1),gR,6,15,[0,1,4,7,9,11,12,14,16,19,20,23,24,27,28,31,32,35,40,39,40,45,48,51,52,55,56,59,58,63,64,69,74,75,80,79,84,85,88,89,90,93,98,0,102,103,106,107,114,115,120,121,130,127,132,133,138,139,140,141,142,0,152,153,158,159,164,165,166,169,174,175,180,181,184,187,192,193,195,197,202,205,208,209,0,0,0,0,0,0,232,0,238,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,71,156,114,115,103,128,129,57,137,113,113,128,131,147,97,87,101,186,163,99]);kr=aC(VB(kF,1),vR,2,6,[SR,ST,TR]);}function Qm(a,b){var c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,A,B,C,D,F;q=(h=Om(b,1),h==-1?b.length:h);d=xI(b.substr(0,q));p=Nm(b,q);q=(i=Om(b,p+1),i==-1?b.length:i);s=b.substr(p,q-p);A=null;e=false;r=Wm(b);if(r!=0){A=Rm(b);r<0&&(e=true);q=r<0?-r:r;}p=Nm(b,q);q=(j=Om(b,p+1),j==-1?b.length:j);C=LI(b.substr(p,q-p));p=Nm(b,q);q=(k=Om(b,p+1),k==-1?b.length:k);D=LI(b.substr(p,q-p));p=Nm(b,q);q=(l=Om(b,p+1),l==-1?b.length:l);F=LI(b.substr(p,q-p));p=Nm(b,q);q=(m=Om(b,p+1),m==-1?b.length:m);u=xI(b.substr(p,q-p));c=Hh(a.c,C,-D,-F);c+1!=d&&(!a.a&&(a.a=new hO()),aO(a.a,new QI(d),new QI(c)));A!=null&&Qj(a.c,c,A,e);u!=0&&Rj(a.c,c,u,false);if(DJ(s,'A')){Vj(a.c,c,1,true);}else if(DJ(s,'Q')){t=ZB(_C,DQ,6,1,15,1);t[0]=6;Qj(a.c,c,t,true);}else{ak(a.c,c,Ck(s));}while((p=Nm(b,q))!=-1){q=(g=Om(b,p+1),g==-1?b.length:g);v=b.substr(p,q-p);o=FJ(v,OJ(61));n=v.substr(0,o);B=xI(v.substr(o+1,v.length-(o+1)));if(DJ(n,'CHG')){Jj(a.c,c,B);}else if(DJ(n,'RAD')){switch(B){case 1:Wj(a.c,c,16);break;case 2:Wj(a.c,c,32);break;case 3:Wj(a.c,c,48);}}else if(DJ(n,'CFG'));else if(DJ(n,'MASS')){Tj(a.c,c,B);}else if(DJ(n,'VAL')){Hj(a.c,c,B==-1?0:B==0?-1:B);}else if(DJ(n,'HCOUNT')){switch(B){case 0:break;case-1:Vj(a.c,c,1792,true);break;case 1:Vj(a.c,c,128,true);break;case 2:Vj(a.c,c,384,true);break;default:Vj(a.c,c,896,true);}}else if(DJ(n,'SUBST')){if(B==-1){Vj(a.c,c,YQ,true);}else if(B>0){w=0;for(f=0;f<a.c.p;f++){(Ei(a.c,0,f)==c||Ei(a.c,1,f)==c)&&++w;}B>w&&Vj(a.c,c,xQ,true);}}else if(DJ(n,'RBCNT')){switch(B){case 3:case-1:Vj(a.c,c,112,true);break;case 1:Vj(a.c,c,8,true);break;case 2:Vj(a.c,c,104,true);break;case 4:Vj(a.c,c,56,true);}}}}function yg(a,b,c){var d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t;j=aC(VB(_C,2),mR,7,0,[aC(VB(_C,1),DQ,6,15,[627]),null,aC(VB(_C,1),DQ,6,15,[2457]),null,aC(VB(_C,1),DQ,6,15,[2451,8643,2519]),null,aC(VB(_C,1),DQ,6,15,[34377,-2147448999]),null,aC(VB(_C,1),DQ,6,15,[37449,137313,95703,34371,37815,54891,132867,-2147309741,54857,55129,-2147449005,-2147449065]),null,aC(VB(_C,1),DQ,6,15,[530697,531819,899169,137289,694617,-2146951863,-2146952797,-2146939175,-2146929547,-2146929564,-2146625111,-2146931799,-2146940503,-2146931935]),null,aC(VB(_C,1),DQ,6,15,[542985,137283,2122017,530691,2206773,-2144711351,219209,2840841,137555,-2146871031,-2147264167,613705,-2145360543,-2146625271,694611,2454837,-2145356703,-2147345133,-2146928951,-2146931805,-2144641719,-2146951869,-2146625237,-2146624183,2841963,1074905,-2146625117,2799955,-2144723645,138583,859225,-2145264843,-2145216253,-2146624149,-2144700727,-2146928917,-2143905527,-2144045771,-2146789097,2288547,544407,2104323,-2146911977,-2144479405,3633737,-2146870089,-2146952169]),null,aC(VB(_C,1),DQ,6,15,[8487297,2172633,2116611,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8829813])]);q=b.a.length-10;if(b.a.length>=10&&b.a.length<=24&&j[q]!=null){o=1<<b.a.length;f=0;h=0;for(l=0;l<b.a.length;l++){if(Mi(a.i,c[l])==2){g=Ni(a.i,c[l]);g==1&&(f+=o);g==2&&(h+=o);}f>>>=1;h>>>=1;}for(t=0;t<j[q].length;t++){n=(sQ&j[q][t])==0;i=oQ&j[q][t];for(m=false;!m;m=!m){if(m){if(n)break;p=0;for(d=1;d!=o;d<<=1){p<<=1;(i&d)!=0&&(p|=1);}i=p;}for(r=0;r<b.a.length;r++){if((i&f)==0&&(~i&h)==0){e=0;s=true;for(k=1;k<b.a.length;k++){b.c[k]=b.c[k-1]+$wnd.Math.sin(e);b.d[k]=b.d[k-1]+$wnd.Math.cos(e);(i&1)==0&&(s=!s);e+=s?dR:sR;i>>>=1;}return;}(i&1)!=0&&(i|=o);i>>>=1;}}}}zg(b);}function Em(a){var b,c,d,e,f,g,h,i,j,k,l;Xo(a,1);e=ZB(_C,DQ,6,191,15,1);for(c=0;c<a.o;c++){switch(a.A[c]){case 171:e[1]+=5;e[6]+=3;e[7]+=1;e[8]+=1;break;case 172:e[1]+=12;e[6]+=6;e[7]+=4;e[8]+=1;break;case 173:e[1]+=6;e[6]+=4;e[7]+=2;e[8]+=2;break;case 174:e[1]+=5;e[6]+=4;e[7]+=1;e[8]+=3;break;case 175:e[1]+=5;e[6]+=3;e[7]+=1;e[8]+=1;e[16]+=1;break;case 176:e[1]+=8;e[6]+=5;e[7]+=2;e[8]+=2;break;case 177:e[1]+=7;e[6]+=5;e[7]+=1;e[8]+=3;break;case 178:e[1]+=3;e[6]+=2;e[7]+=1;e[8]+=1;break;case 179:e[1]+=7;e[6]+=6;e[7]+=3;e[8]+=1;break;case 181:case 180:e[1]+=11;e[6]+=6;e[7]+=1;e[8]+=1;break;case 182:e[1]+=12;e[6]+=6;e[7]+=2;e[8]+=1;break;case 183:e[1]+=9;e[6]+=5;e[7]+=1;e[8]+=1;e[16]+=1;break;case 184:e[1]+=9;e[6]+=9;e[7]+=1;e[8]+=1;break;case 185:e[1]+=7;e[6]+=5;e[7]+=1;e[8]+=1;break;case 186:e[1]+=5;e[6]+=3;e[7]+=1;e[8]+=2;break;case 187:e[1]+=7;e[6]+=4;e[7]+=1;e[8]+=2;break;case 188:e[1]+=10;e[6]+=11;e[7]+=2;e[8]+=1;break;case 189:e[1]+=9;e[6]+=9;e[7]+=1;e[8]+=2;break;case 190:e[1]+=9;e[6]+=5;e[7]+=1;e[8]+=1;break;case 1:switch(a.v[c]){case 0:case 1:++e[1];break;case 2:++e[151];break;case 3:++e[152];}break;default:++e[a.A[c]];}}for(d=0;d<a.o;d++)a.A[d]>=171&&a.A[d]<=190?e[1]+=2-ol(a,d):e[1]+=ml(a,d);h=0;for(j=1;j<=190;j++)e[j]!=0&&++h;this.b=ZB(_C,DQ,6,h,15,1);this.c=ZB(_C,DQ,6,h,15,1);h=0;for(i=0;i<ym.length;i++){if(e[ym[i]]!=0){this.b[h]=e[ym[i]];this.c[h]=ym[i];++h;e[ym[i]]=0;}}while(true){l='zzz';k=-1;for(g=1;g<=190;g++)if(e[g]>0&&zJ(l,(Gh(),Ch)[g])>0){l=(Gh(),Ch)[g];k=g;}if(k==-1)break;this.b[h]=e[k];this.c[h]=k;++h;e[k]=0;}this.a=0;this.d=0;for(b=0;b<a.d;b++){if(a.A[b]!=1&&a.v[b]!=0){g=a.A[b];f=a.v[b];this.a+=wm(g,f)-xm[g];this.d+=wm(g,f)-zm[g];}}}function vq(a,b){var c;switch(a.A[b]){case 7:if((a.s[b]&xQ)!=0){if(a.q[b]==0){if(a.c[b]-a.g[b]+ml(a,b)==0){if(a.g[b]==2)return 18;else{for(c=0;c<a.g[b];c++)if(!Fl(a,a.i[b][c]))return 20;return 19;}}else return 22;}else if(a.q[b]==1){if(a.c[b]-a.g[b]+ml(a,b)==0){for(c=0;c<a.g[b];c++)if(!Fl(a,a.i[b][c]))return ji(a,a.f[b][c])<0?21:24;return 23;}else return 25;}}else{if(a.q[b]==0){switch(a.c[b]-a.g[b]+ml(a,b)){case 0:switch(a.k[b]){case 0:return(!!a.n&&b<a.d?jn(a.n,b):0)==3?5:0;case 1:return 1;case 2:return 2;}break;case 1:switch(a.k[b]){case 0:return(!!a.n&&b<a.d?jn(a.n,b):0)==3?7:6;case 1:return 8;}break;case 2:return 9;}}else if(a.q[b]==1){switch(a.c[b]-a.g[b]+ml(a,b)){case 0:switch(a.k[b]){case 0:return 10;case 1:return xq(a,b)?3:11;case 2:return a.j[b][0]==2?xq(a,b)?4:qq.length+1:12;}break;case 1:switch(a.k[b]){case 0:return 13;case 1:return 14;}break;case 2:return a.k[b]==0?15:16;case 3:return 17;}}}return qq.length+1;case 8:if((a.s[b]&xQ)!=0){if(a.q[b]==0)return 31;}else{if(a.q[b]==0){if(a.k[b]>0)return 28;if(a.g[b]==1)return 29;if((!!a.n&&b<a.d?jn(a.n,b):0)==3)return 27;return 26;}else if(a.q[b]==-1){if(a.g[b]==1&&ji(a,a.f[b][0])>0)return 28;return 30;}}return qq.length+1;case 15:if(a.q[b]==0){if(a.c[b]-a.g[b]+ml(a,b)==0){if(a.g[b]==3&&a.k[b]==0)return 39;if(a.g[b]==2&&a.k[b]==1)return 40;if(a.g[b]==4&&a.k[b]==1)return 41;}else if(a.c[b]-a.g[b]+ml(a,b)==1){if(a.g[b]==3&&a.k[b]==1)return 42;}}return qq.length+1;case 16:if(a.q[b]==0){if((a.s[b]&xQ)!=0){return a.g[b]==2?37:38;}else{if(a.c[b]-a.g[b]+ml(a,b)==0){if(a.g[b]==2&&a.k[b]==0)return 32;if(a.g[b]==1&&a.k[b]==1)return 33;if(a.g[b]==3&&a.k[b]==1)return 34;if(a.g[b]==4&&a.k[b]==2)return 35;}else if(a.c[b]-a.g[b]+ml(a,b)==1){if(a.g[b]==1)return 36;}}}return qq.length+1;}return qq.length;}function Fe(a,b,c,d){var e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,A,B,C,D,F,G,H,I,J,K,L,M;if(Ai(a.L,c)!=Ai(a.L,d))return Ai(a.L,c)>Ai(a.L,d);if(ti(a.L,c)!=ti(a.L,d)){H=nj(a.L,c)?(Gh(),Eh)[Ai(a.L,c)]:ti(a.L,c);I=nj(a.L,d)?(Gh(),Eh)[Ai(a.L,d)]:ti(a.L,d);return H>I;}w=a.L.d;s=ZB(_C,DQ,6,w,15,1);u=ZB(_C,DQ,6,w,15,1);v=ZB(_C,DQ,6,w,15,1);t=ZB(aG,HQ,6,w,16,1);i=ZB(aG,HQ,6,a.L.o,16,1);s[0]=b;s[1]=c;s[2]=d;u[0]=-1;u[1]=0;u[2]=0;i[b]=true;i[c]=true;i[d]=true;m=1;A=2;G=ZB(_C,DQ,6,64,15,1);G[1]=1;G[2]=3;o=2;while(m<=A){while(m<G[o]){n=s[m];if(!t[m]){p=0;q=0;for(C=0;C<bl(a.L,n);C++){k=al(a.L,n,C);if(A+dl(a.L,n,C)+1>=w){w+=a.L.d;s=kf(s,w);u=kf(u,w);v=kf(v,w);t=(M=ZB(aG,HQ,6,w,16,1),dK(t,M,t.length),M);}if(Hl(a.L,cl(a.L,n,C))){++p;q+=Ai(a.L,k);}else{for(F=1;F<dl(a.L,n,C);F++){++A;s[A]=k;u[A]=m;t[A]=true;}}K=u[m];if(k==s[K])continue;h=false;if(i[k]){J=u[K];while(J!=-1){if(k==s[J]){h=true;break;}J=u[J];}}if(h){++A;s[A]=k;u[A]=m;t[A]=true;}else{++A;s[A]=k;u[A]=m;i[k]=true;}}if(p!=0){++A;v[A]=(q<<2)/p|0;u[A]=m;t[A]=true;}}++m;if(m==10000){throw dG(new FA('Emergency break in while loop.'));}}G.length==o+1&&(G=kf(G,G.length+64));G[o+1]=A+1;for(B=G[o];B<G[o+1];B++){v[B]==0&&(v[B]=(Ai(a.L,s[B])==151?1:Ai(a.L,s[B])==152?1:Ai(a.L,s[B]))<<2);v[B]+=v[u[B]]<<16;}Je(a,t,v,u,s,G,o);if(v[1]!==v[2])return v[1]>v[2];o>1&&Ge(v,u,G,o);++o;}l=ZB(_C,DQ,6,a.L.d,15,1);D=false;for(f=0;f<a.L.d;f++){if(i[f]&&!nj(a.L,f)){D=true;break;}}if(D){for(g=0;g<a.L.d;g++)l[g]=nj(a.L,g)?(Gh(),Eh)[Ai(a.L,g)]:ti(a.L,g);if(Ie(a,t,v,u,s,l,G,o))return v[1]>v[2];}oN(l,l.length,0);r=false;for(j=0;j<a.L.e;j++){if(i[Ei(a.L,0,j)]||i[Ei(a.L,1,j)]){if(a.f[j]==1){l[Ei(a.L,0,j)]=1;l[Ei(a.L,1,j)]=1;r=true;}else if(a.f[j]==2){l[Ei(a.L,0,j)]=2;l[Ei(a.L,1,j)]=2;r=true;}}}if(r&&Ie(a,t,v,u,s,l,G,o))return v[1]>v[2];oN(l,l.length,0);L=false;for(e=0;e<a.L.d;e++){if(i[e]){if(a.R[e]==2){l[e]=1;L=true;}else if(a.R[e]==1){l[e]=2;L=true;}}}if(L&&Ie(a,t,v,u,s,l,G,o))return v[1]>v[2];throw dG(new FA('no distinction applying CIP rules'));}function Vd(a,b){var c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,A;if(a.A[b]!=7)return false;if(a.g[b]+a.k[b]>3)return false;if((a.s[b]&xQ)!=0){if(a.k[b]!=1)return false;if(Wk(a,b,7)!=1)return false;u=(Xo(a,3),a.n);for(s=0;s<u.g.a.length;s++){if(pn(u,s,b)){if(BM(u.i,s).length==5||BM(u.i,s).length==6){v=BM(u.g,s);q=-1;for(i=0;i<v.length;i++){if(v[i]==b){q=i;break;}}e=0;r=null;p=null;if(v.length==5){r=ZB(_C,DQ,6,2,15,1);r[0]=v[q-1<0?q+4:q-1];r[1]=v[q-4<0?q+1:q-4];p=ZB(_C,DQ,6,2,15,1);p[0]=v[q-2<0?q+3:q-2];p[1]=v[q-3<0?q+2:q-3];}if(v.length==6){r=ZB(_C,DQ,6,3,15,1);r[0]=v[q-1<0?q+5:q-1];r[1]=v[q-3<0?q+3:q-3];r[2]=v[q-5<0?q+1:q-5];p=ZB(_C,DQ,6,2,15,1);p[0]=v[q-2<0?q+4:q-2];p[1]=v[q-4<0?q+2:q-4];}for(j=0;j<v.length;j++)b!=v[j]&&a.A[v[j]]==7&&a.k[v[j]]==1&&--e;for(k=0;k<r.length;k++){f=-1;g=-1;for(o=0;o<a.g[r[k]];o++){if(!Fl(a,a.i[r[k]][o])){f=a.f[r[k]][o];g=a.i[r[k]][o];break;}}if(f!=-1){if(a.A[f]==7&&a.k[f]==0&&a.g[f]+a.k[f]<=3&&!Wd(a,f,false)){++e;continue;}if(a.A[f]==8&&a.g[f]==1){e+=2;continue;}if((a.C[g]&256)!=0){for(w=0;w<u.g.a.length;w++){if(u.d[w]&&pn(u,w,f)){t=BM(u.g,w);for(n=0;n<t.length;n++){if(a.A[t[n]]==7&&a.k[t[n]]==1){--e;break;}}break;}}}}}for(l=0;l<p.length;l++){f=-1;for(n=0;n<a.g[p[l]];n++)Fl(a,a.i[p[l]][n])||(f=a.f[p[l]][n]);if(a.A[p[l]]==7){a.k[p[l]]==0&&(f==-1||Rd(a,f)==0)&&++e;continue;}if(f!=-1&&Rd(a,f)!=0){--e;continue;}}return e>0;}break;}}return false;}if(a.k[b]>1)return false;if(a.k[b]==1){m=-1;A=0;for(i=0;i<a.g[b];i++){d=a.f[b][i];if(a.j[b][i]==2){if(a.A[d]!=6)return false;m=d;continue;}if(a.A[d]==8)return false;if(a.A[d]==7){--A;Wd(a,d,false)&&--A;continue;}(a.s[d]&xQ)!=0&&--A;}if(m==-1)return false;c=0;for(j=0;j<a.g[m];j++){if(a.j[m][j]==1){d=a.f[m][j];if(Rd(a,d)!=0)return false;(a.s[d]&xQ)!=0&&++c;a.A[d]==7&&!Wd(a,d,true)&&++A;(a.A[d]==8||a.A[d]==16)&&--A;}}c==2&&--A;return A>=0;}for(h=0;h<a.g[b];h++){d=a.f[b][h];if((a.s[d]&xQ)!=0)return false;if(a.A[d]!=6)return false;if(Rd(a,d)!=0)return false;if(a.k[d]!=0&&Xd(a,d))return false;}return true;}function se(a){var b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w;l=false;if(a.L.I){for(j=0;j<a.L.e;j++){if(Oi(a.L,j)!=0){l=true;break;}}}a.I=2;for(c=0;c<a.L.d;c++)a.I=lJ(a.I,qe(a,c));i=lJ(2,l?(78+a.I*36)/63|0:(78+a.I*21)/63|0);a.c=ZB(_C,DQ,6,a.L.o,15,1);a.b=ZB(jD,nR,67,a.L.d,0,1);for(d=0;d<a.L.d;d++)a.b[d]=new Gf(i);h=false;for(e=0;e<a.L.d;e++){Ff(a.b[e],e);(vi(a.L,e)&1)!=0||qi(a.L,e)!=null?Cf(a.b[e],8,6):Cf(a.b[e],8,Ai(a.L,e));Cf(a.b[e],8,ti(a.L,e));Cf(a.b[e],2,Tk(a.L,e));Cf(a.b[e],3,qe(a,e));(vi(a.L,e)&1)!=0?Cf(a.b[e],4,8):Cf(a.b[e],4,8+ji(a.L,e));Cf(a.b[e],5,mJ(31,Xk(a.L,e)));Cf(a.b[e],4,(t=hi(a.L,e),u=ll(a.L,e,false),v=ll(a.L,e,true),w=-1,u!=v?t!=-1&&t>u?w=t<<24>>24:w=u<<24>>24:t!=-1?(t>v||t<v&&t>=ol(a.L,e))&&(w=t<<24>>24):!am(a.L,e)&&ml(a.L,e)!=0&&(w=ol(a.L,e)-Si(a.L,e)),Ce(a,e,w),w)+1);Cf(a.b[e],2,wi(a.L,e)>>4);if(a.L.I){Cf(a.b[e],30,vi(a.L,e));qi(a.L,e)!=null&&(h=true);}}a.N=ve(a);if(a.N<a.L.d){for(b=0;b<a.L.d;b++){Ff(a.b[b],b);Cf(a.b[b],16,a.c[b]);m=ZB(_C,DQ,6,bl(a.L,b),15,1);for(o=0;o<bl(a.L,b);o++){m[o]=a.c[al(a.L,b,o)]<<5;m[o]|=mJ(31,_k(a.L,cl(a.L,b,o)));}xN(m);for(p=a.I;p>m.length;p--)Cf(a.b[b],21,0);for(n=m.length-1;n>=0;n--)Cf(a.b[b],21,m[n]);}a.N=ve(a);}if(h&&a.N<a.L.d){for(b=0;b<a.L.d;b++){Ff(a.b[b],b);Cf(a.b[b],16,a.c[b]);g=qi(a.L,b);r=g==null?0:mJ(12,g.length);for(o=12;o>r;o--)Cf(a.b[b],8,0);for(n=r-1;n>=0;n--)Cf(a.b[b],8,g[n]);}a.N=ve(a);}if(l&&a.N<a.L.d){for(b=0;b<a.L.d;b++){Ff(a.b[b],b);Cf(a.b[b],16,a.c[b]);k=ZB(aD,gR,6,qe(a,b),14,1);for(o=0;o<qe(a,b);o++){k[o]=a.c[al(a.L,b,o)];k[o]=pG(k[o],20);k[o]=oG(k[o],Oi(a.L,cl(a.L,b,o)));}wN(k,gG);for(p=a.I;p>k.length;p--)Cf(a.b[b],36,0);for(n=k.length-1;n>=0;n--)Cf(a.b[b],36,k[n]);}a.N=ve(a);}if((a.K&8)!=0&&a.N<a.L.d){q=new Uo();for(f=0;f<a.L.d;f++)li(a.L,f)!=null&&So(q,li(a.L,f));for(b=0;b<a.L.d;b++){s=li(a.L,b)==null?0:1+To(q,li(a.L,b));Ff(a.b[b],b);Cf(a.b[b],16,a.c[b]);Cf(a.b[b],16,s);}a.N=ve(a);}if((a.K&16)!=0&&a.N<a.L.d){for(b=0;b<a.L.d;b++){Ff(a.b[b],b);Cf(a.b[b],16,a.c[b]);Cf(a.b[b],1,qj(a.L,b)?1:0);}a.N=ve(a);}}function hd(a,b){var c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u;n=new Fd();c=new Fd();f=new Fd();l=new hH();k=new hH();d=Ei(a.G,0,b);e=Ei(a.G,1,b);((vi(a.G,d)|vi(a.G,e))&IQ)!=0;so(a,d,e,uh(a.K,xi(a.G,d)),vh(a.K,yi(a.G,d)),uh(a.K,xi(a.G,e)),vh(a.K,yi(a.G,e)));!qj(a.G,d)&&!qj(a.G,e)&&((vi(a.G,d)|vi(a.G,e))&IQ)!=0&&zd(a,-8);if(!a.n[d]){n.a=uh(a.K,xi(a.G,d));n.c=vh(a.K,yi(a.G,d));}else{n.a=a.n[d].a;n.c=a.n[d].b;}if(!a.n[e]){n.b=uh(a.K,xi(a.G,e));n.d=vh(a.K,yi(a.G,e));}else{n.b=a.n[e].a;n.d=a.n[e].b;}if((Oi(a.G,b)&$Q)!=0){vd(a,n)&&(q=WC(n.a),r=WC(n.b),s=WC(n.c),t=WC(n.d),u='<line stroke-dasharray="3, 3" x1="'+q+_Q+'y1="'+s+_Q+'x2="'+r+_Q+'y2="'+t+_Q+'stroke="'+a.d+_Q+aR+WC(a.i)+'"/>',xo(a,u),undefined);zd(a,-9);return;}g=Pi(a.G,b)==64?0:Pi(a.G,b)==32?1:Mi(a.G,b);switch(g){case 1:switch(Pi(a.G,b)){case 1:vd(a,n)&&Pc(a,n,d,e);break;case 17:rd(a,n,d,e);break;case 9:o=n.b-n.a;p=n.d-n.c;if(fj(a.G,$k(a.G,d,e))){h=-3;i=-3;}else{h=a.o[d];i=Vc(a,d);h==ki(a.G,d)&&(h=i);}for(j=2;j<17;j+=2){c.a=n.a+j*o/17-j*p/128;c.c=n.c+j*p/17+j*o/128;c.b=n.a+j*o/17+j*p/128;c.d=n.c+j*p/17-j*o/128;if(vd(a,c)){zd(a,j<9?h:i);no(a,c);zd(a,a.J);}}break;case 32:vd(a,n)&&Qc(a,n,d,e);}break;case 0:case 2:if((a.q[d]||Tk(a.G,d)==2)&&(a.q[e]||Tk(a.G,e)==2)&&!Ll(a.G,b)&&g==2){if(!vd(a,n))break;cd(a,n.b-n.a,n.d-n.c,l);o=l.a/2;p=l.b/2;c.a=n.a+o;c.c=n.c+p;c.b=n.b+o;c.d=n.d+p;f.a=n.a-o;f.c=n.c-p;f.b=n.b-o;f.d=n.d-p;Pi(a.G,b)==26&&td(c,f);Pc(a,c,d,e);Pc(a,f,d,e);}else if((a.q[e]||Tk(a.G,e)==2)&&g==2){dd(a,n,b,false);}else if((a.q[d]||Tk(a.G,d)==2)&&g==2){dd(a,n,b,true);}else{m=ud(a,b);m==0&&(m=1);c.a=n.a;c.c=n.c;c.b=n.b;c.d=n.d;cd(a,n.b-n.a,n.d-n.c,l);if(m>0){f.a=n.a+l.a;f.c=n.c+l.b;f.b=n.b+l.a;f.d=n.d+l.b;if(bd(a,d,e,1,k)||bl(a.G,d)>1){f.a+=k.a+l.b;f.c+=k.b-l.a;}if(bd(a,e,d,-1,k)||bl(a.G,e)>1){f.b+=k.a-l.b;f.d+=k.b+l.a;}}else{f.a=n.a-l.a;f.c=n.c-l.b;f.b=n.b-l.a;f.d=n.d-l.b;if(bd(a,d,e,-1,k)||bl(a.G,d)>1){f.a+=k.a+l.b;f.c+=k.b-l.a;}if(bd(a,e,d,1,k)||bl(a.G,e)>1){f.b+=k.a-l.b;f.d+=k.b+l.a;}}Pi(a.G,b)==26&&td(c,f);vd(a,c)&&Pc(a,c,d,e);g==2?vd(a,f)&&Pc(a,f,d,e):vd(a,f)&&Nc(a,f,d,e);}break;case 3:if(vd(a,n)){Pc(a,n,d,e);cd(a,n.b-n.a,n.d-n.c,l);c.a=n.a+l.a;c.c=n.c+l.b;c.b=n.b+l.a;c.d=n.d+l.b;Pc(a,c,d,e);c.a=n.a-l.a;c.c=n.c-l.b;c.b=n.b-l.a;c.d=n.d-l.b;Pc(a,c,d,e);}}a.w==-8&&zd(a,-9);}function Am(){Am=FG;zm=aC(VB(ZC,1),GQ,6,15,[0,1.00794,4.0026,6.941,9.0122,10.811,12.011,14.007,15.999,18.998,20.18,22.99,24.305,26.982,28.086,30.974,32.066,35.453,39.948,39.098,40.078,44.956,47.867,50.942,51.996,54.938,55.845,58.933,58.693,63.546,65.39,69.723,72.61,74.922,78.96,79.904,83.8,85.468,87.62,88.906,91.224,92.906,95.94,98.906,101.07,102.91,106.42,107.87,112.41,114.82,118.71,121.76,127.6,126.9,131.29,132.91,137.33,138.91,140.12,140.91,144.24,146.92,150.36,151.96,157.25,158.93,162.5,164.93,167.26,168.93,173.04,174.97,178.49,180.95,183.84,186.21,190.23,192.22,195.08,196.97,200.59,204.38,207.2,208.98,209.98,209.99,222.02,223.02,226.03,227.03,232.04,231.04,238.03,237.05,239.05,241.06,244.06,249.08,252.08,252.08,257.1,258.1,259.1,262.11,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.0141,3.016,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]);xm=aC(VB(ZC,1),GQ,6,15,[0,1.007825,4.0026,7.016003,9.012182,11.009305,12,14.003074,15.994915,18.998403,19.992435,22.989767,23.985042,26.98153,27.976927,30.973762,31.97207,34.968852,39.962384,38.963707,39.962591,44.95591,47.947947,50.943962,51.940509,54.938047,55.934939,58.933198,57.935346,62.939598,63.929145,68.92558,73.921177,74.921594,79.91652,78.918336,83.911507,84.911794,87.905619,88.905849,89.904703,92.906377,97.905406,89.92381,101.904348,102.9055,105.903478,106.905092,113.903357,114.90388,119.9022,120.903821,129.906229,126.904473,131.904144,132.905429,137.905232,138.906346,139.905433,140.907647,141.907719,135.92398,151.919729,152.921225,157.924099,158.925342,163.929171,164.930319,165.93029,168.934212,173.938859,174.94077,179.946545,180.947992,183.950928,186.955744,191.961467,192.962917,194.964766,196.966543,201.970617,204.974401,207.976627,208.980374,193.98818,195.99573,199.9957,201.00411,206.0038,210.00923,232.038054,216.01896,238.050784,229.03623,232.041169,237.05005,238.05302,242.06194,240.06228,243.06947,243.07446,248.08275,251.08887,253.09515,257.10295,257.10777,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.014,3.01605,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]);ym=aC(VB(_C,1),DQ,6,15,[6,1,7,8]);}function Kp(a,b){Ip();var c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,A,B,C,D,F,G,H,I,J,K,L,M,N,O,P;H=ml(a,b);O=a.H[b].a;P=a.H[b].b;switch(H){case 1:{J=a.g[b];if(J==0){n=DR;B=-1;for(f=0;f<a.o;f++){if(f==b)continue;l=O-a.H[f].a;m=P-a.H[f].b;j=$wnd.Math.sqrt(l*l+m*m);if(n>j){n=j;B=f;}}p=O-a.H[B].a;s=P-a.H[B].b;}else{p=O-xi(a,a.f[b][0]);s=P-yi(a,a.f[b][0]);}if(J==1){C=Hh(a,O+Bp*p+Gp*s,P-Gp*p+Bp*s,0);}else if(J==2){p=O-0.5*(xi(a,a.f[b][0])+xi(a,a.f[b][1]));s=P-0.5*(yi(a,a.f[b][0])+yi(a,a.f[b][1]));C=Hh(a,O+p,P+s,0);}else if(J==3){L=a.f[b][0];for(w=1;w<3;w++){i=a.i[b][w];(a.F[i]==9||a.F[i]==17)&&(L=a.f[b][w]);}c=$wnd.Math.abs(Bk(Ak(a.H[b].a,a.H[b].b,xi(a,a.f[b][0]),yi(a,a.f[b][0])),Ak(a.H[b].a,a.H[b].b,xi(a,a.f[b][1]),yi(a,a.f[b][1]))));d=$wnd.Math.abs(Bk(Ak(a.H[b].a,a.H[b].b,xi(a,a.f[b][0]),yi(a,a.f[b][0])),Ak(a.H[b].a,a.H[b].b,xi(a,a.f[b][2]),yi(a,a.f[b][2]))));e=$wnd.Math.abs(Bk(Ak(a.H[b].a,a.H[b].b,xi(a,a.f[b][1]),yi(a,a.f[b][1])),Ak(a.H[b].a,a.H[b].b,xi(a,a.f[b][2]),yi(a,a.f[b][2]))));K=true;if(c>d&&c>e){if(d+e<MQ){K=false;p=O-0.5*(xi(a,a.f[b][0])+xi(a,a.f[b][1]));s=P-0.5*(yi(a,a.f[b][0])+yi(a,a.f[b][1]));}}else if(d>c&&d>e){if(c+e<MQ){K=false;p=O-0.5*(xi(a,a.f[b][0])+xi(a,a.f[b][2]));s=P-0.5*(yi(a,a.f[b][0])+yi(a,a.f[b][2]));}}else{if(c+d<MQ){K=false;p=O-0.5*(xi(a,a.f[b][1])+xi(a,a.f[b][2]));s=P-0.5*(yi(a,a.f[b][1])+yi(a,a.f[b][2]));}}if(K){M=a.f[b][0];o=DR;for(v=0;v<3;v++){g=a.f[b][v];if(g!=L){k=$wnd.Math.pow(a.H[b].a-a.H[g].a,2)+$wnd.Math.pow(a.H[b].b-a.H[g].b,2);if(k<o){M=g;o=k;cK();}}}C=Hh(a,(a.H[L].a+a.H[M].a)/2,(a.H[L].b+a.H[M].b)/2,0);}else{C=Hh(a,O+p,P+s,0);}}else{C=Hh(a,O+p,P+s,0);}ak(a,C,1);Jh(a,b,C,1);}break;case 2:I=a.g[b];if(I==1){p=O-xi(a,a.f[b][0]);s=P-yi(a,a.f[b][0]);C=Hh(a,O+(Cp*p-Hp*s)*0.7,P+(Hp*p+Cp*s)*0.7,0);ak(a,C,1);Jh(a,b,C,1);C=Hh(a,O+(zp*p-Ep*s)*0.7,P+(Ep*p+zp*s)*0.7,0);ak(a,C,1);Jh(a,b,C,1);}else if(I==2){q=O-xi(a,a.f[b][0]);t=P-yi(a,a.f[b][0]);r=O-xi(a,a.f[b][1]);u=P-yi(a,a.f[b][1]);F=$wnd.Math.sqrt(q*q+t*t)*0.7;G=$wnd.Math.sqrt(r*r+u*u)*0.7;p=q+r;s=t+u;D=$wnd.Math.sqrt(p*p+s*s);h=(F+G)/2;p=p/D*h;s=s/D*h;N=wl(a,b);C=Hh(a,O+yp*p-Dp*s,P+Dp*p+yp*s,0);ak(a,C,1);N>-1?Jh(a,b,C,1):Jh(a,b,C,17);C=Hh(a,O+Ap*p-Fp*s,P+Fp*p+Ap*s,0);ak(a,C,1);Jh(a,b,C,1);}else{for(A=0;A<2;A++){C=Hh(a,O,P,0);ak(a,C,1);Jh(a,b,C,1);}}break;case 3:{p=(O-xi(a,a.f[b][0]))*0.7;s=(P-yi(a,a.f[b][0]))*0.7;C=Hh(a,O+p,P+s,0);ak(a,C,1);Jh(a,b,C,1);C=Hh(a,O-s,P+p,0);ak(a,C,1);Jh(a,b,C,1);C=Hh(a,O+s,P-p,0);ak(a,C,1);Jh(a,b,C,1);}break;default:{for(A=0;A<H;A++){C=Hh(a,O,P,0);ak(a,C,1);Jh(a,b,C,1);}break;}}}function Gh(){Gh=FG;Ch=aC(VB(kF,1),vR,2,6,['?','H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr',wR,wR,wR,wR,wR,wR,wR,wR,wR,wR,wR,wR,wR,wR,wR,wR,wR,wR,wR,wR,wR,wR,wR,wR,wR,'R4','R5','R6','R7','R8','R9','R10','R11','R12','R13','R14','R15','R16','R1','R2','R3','A','A1','A2','A3',wR,wR,'D','T','X','R','H2','H+','Nnn','HYD','Pol',wR,wR,wR,wR,wR,wR,wR,wR,wR,wR,wR,'Ala','Arg','Asn','Asp','Cys','Gln','Glu','Gly','His','Ile','Leu','Lys','Met','Phe','Pro','Ser','Thr','Trp','Tyr','Val']);Eh=aC(VB(_F,1),gR,6,15,[0,1,4,7,9,11,12,14,16,19,20,23,24,27,28,31,32,35,40,39,40,45,48,51,52,55,56,59,58,63,64,69,74,75,80,79,84,85,88,89,90,93,98,0,102,103,106,107,114,115,120,121,130,127,132,133,138,139,140,141,142,0,152,153,158,159,164,165,166,169,174,175,180,181,184,187,192,193,195,197,202,205,208,209,0,0,0,0,0,0,232,0,238,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,71,156,114,115,103,128,129,57,137,113,113,128,131,147,97,87,101,186,163,99]);Dh=aC(VB(XC,2),xR,13,0,[null,aC(VB(XC,1),pR,6,15,[1]),aC(VB(XC,1),pR,6,15,[0]),aC(VB(XC,1),pR,6,15,[1]),aC(VB(XC,1),pR,6,15,[2]),aC(VB(XC,1),pR,6,15,[3]),aC(VB(XC,1),pR,6,15,[4]),aC(VB(XC,1),pR,6,15,[3]),aC(VB(XC,1),pR,6,15,[2]),aC(VB(XC,1),pR,6,15,[1]),aC(VB(XC,1),pR,6,15,[0]),aC(VB(XC,1),pR,6,15,[1]),aC(VB(XC,1),pR,6,15,[2]),aC(VB(XC,1),pR,6,15,[3]),aC(VB(XC,1),pR,6,15,[4]),aC(VB(XC,1),pR,6,15,[3,5]),aC(VB(XC,1),pR,6,15,[2,4,6]),aC(VB(XC,1),pR,6,15,[1,3,5,7]),aC(VB(XC,1),pR,6,15,[0]),aC(VB(XC,1),pR,6,15,[1]),aC(VB(XC,1),pR,6,15,[2]),null,null,null,null,null,null,null,null,null,null,aC(VB(XC,1),pR,6,15,[2,3]),aC(VB(XC,1),pR,6,15,[2,4]),aC(VB(XC,1),pR,6,15,[3,5]),aC(VB(XC,1),pR,6,15,[2,4,6]),aC(VB(XC,1),pR,6,15,[1,3,5,7]),aC(VB(XC,1),pR,6,15,[0,2]),aC(VB(XC,1),pR,6,15,[1,2,3,4]),aC(VB(XC,1),pR,6,15,[2]),null,null,null,null,null,null,null,null,null,null,aC(VB(XC,1),pR,6,15,[1,2,3]),aC(VB(XC,1),pR,6,15,[2,4]),aC(VB(XC,1),pR,6,15,[3,5]),aC(VB(XC,1),pR,6,15,[2,4,6]),aC(VB(XC,1),pR,6,15,[1,3,5,7]),aC(VB(XC,1),pR,6,15,[0,2,4,6]),aC(VB(XC,1),pR,6,15,[1]),aC(VB(XC,1),pR,6,15,[2])]);}function Pg(a){var b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,A,B,C,D,F,G,H,I,J,K,L,M;for(d=0;d<a.b;d++){if(a.e[d]>4){m=new hh(a,a.i,1+a.e[d]);m.c[a.e[d]]=0;m.d[a.e[d]]=0;m.q[a.e[d]]=32;m.a[a.e[d]]=d;a.a[d]=true;for(o=0;o<a.e[d];o++){j=al(a.i,d,o);m.c[o]=$wnd.Math.sin(dR*o-uR);m.d[o]=$wnd.Math.cos(dR*o-uR);m.q[o]=32;m.a[o]=j;a.a[j]=true;a.c[cl(a.i,d,o)]=true;}yM(a.f,m);}}I=tl(a.i);for(H=0;H<I.g.a.length;H++){J=BM(I.i,H).length;F=BM(I.g,H);K=false;if((a.g&6)!=0){K=true;for(o=0;o<J;o++){if(!lj(a.i,F[o])){K=false;break;}}}if(!K){r=false;for(p=0;p<J;p++){if(Xk(a.i,F[p])==J){r=true;break;}}if(r){G=BM(I.i,H);sg(a,F,G);for(o=0;o<J;o++){a.a[F[o]]=true;a.c[G[o]]=true;}}}}for(h=0;h<a.d;h++){if(Ll(a.i,h)&&!a.c[h]){M=Gg(a,h);F=M.a;G=M.b;sg(a,F,G);for(o=0;o<M.a.length;o++){a.a[F[o]]=true;a.c[G[o]]=true;}}}for(i=0;i<a.d;i++){if(!a.c[i]&&Mi(a.i,i)==3){e=Ei(a.i,0,i);f=Ei(a.i,1,i);w=a.e[e]+a.e[f];if(w>2){m=new hh(a,a.i,w);k=0;for(p=0;p<a.e[e];p++){j=al(a.i,e,p);if(j!=f){m.a[k++]=j;a.a[j]=true;a.c[cl(a.i,e,p)]=true;}}m.a[k++]=e;m.a[k++]=f;for(q=0;q<a.e[f];q++){j=al(a.i,f,q);if(j!=e){m.a[k++]=j;a.a[j]=true;a.c[cl(a.i,f,q)]=true;}}for(o=0;o<w;o++){m.c[o]=o;m.d[o]=0;m.q[o]=1;}a.a[e]=true;a.a[f]=true;a.c[i]=true;yM(a.f,m);}}}for(g=0;g<a.d;g++){if(!a.c[g]&&Mi(a.i,g)==2){b=ZB(_C,DQ,6,a.b,15,1);for(o=0;o<2;o++){b[0]=Ei(a.i,o,g);b[1]=Ei(a.i,1-o,g);if(Tk(a.i,b[0])==1&&Tk(a.i,b[1])==2&&a.e[b[1]]==2){a.a[b[0]]=true;a.a[b[1]]=true;a.c[g]=true;v=1;do{A=al(a.i,b[v],0)==b[v-1]?1:0;b[v+1]=al(a.i,b[v],A);if(Tk(a.i,b[v+1])==2&&a.e[b[v+1]]>2)break;a.a[b[v+1]]=true;a.c[cl(a.i,b[v],A)]=true;++v;}while(Tk(a.i,b[v])==2&&a.e[b[v]]==2);w=a.e[b[0]]+a.e[b[v]]+v-1;m=new hh(a,a.i,w);for(t=0;t<=v;t++){m.c[t]=t;m.d[t]=0;m.q[t]=64;m.a[t]=b[t];}l=v+1;n=false;for(u=0;u<a.e[b[0]];u++){j=al(a.i,b[0],u);if(j!=b[1]){m.c[l]=-0.5;m.d[l]=n?$wnd.Math.sin(dR):-$wnd.Math.sin(dR);m.q[l]=64;m.a[l]=j;++l;n=true;}}n=false;for(s=0;s<a.e[b[v]];s++){j=al(a.i,b[v],s);if(j!=b[v-1]){m.c[l]=v+0.5;m.d[l]=n?-$wnd.Math.sin(dR):$wnd.Math.sin(dR);m.q[l]=64;m.a[l]=j;++l;n=true;}}yM(a.f,m);}}}}for(c=0;c<a.b;c++){if(a.e[c]==4){B=ZB(_C,DQ,6,4,15,1);C=ZB(_C,DQ,6,4,15,1);D=0;for(p=0;p<4;p++){B[D]=al(a.i,c,p);C[D]=cl(a.i,c,p);a.e[B[D]]==1&&!a.c[C[D]]&&++D;}if(D==2){m=new hh(a,a.i,3);for(o=0;o<2;o++){a.a[B[o]]=true;a.c[C[o]]=true;m.a[o]=B[o];m.q[o]=32;}m.c[0]=-0.5;m.d[0]=0.866;m.c[1]=0.5;m.d[1]=0.866;m.c[2]=0;m.d[2]=0;m.q[2]=32;m.a[2]=c;yM(a.f,m);}if(D==3){for(q=0;q<2;q++){if(Mi(a.i,C[q])==1){L=B[q];B[q]=B[2];B[2]=L;L=C[q];C[q]=C[2];C[2]=L;}}m=new hh(a,a.i,4);for(o=0;o<3;o++){a.a[B[o]]=true;a.c[C[o]]=true;m.a[o]=B[o];m.q[o]=32;}m.c[0]=-1;m.d[0]=0;m.c[1]=1;m.d[1]=0;m.c[2]=0;m.d[2]=1;m.c[3]=0;m.d[3]=0;m.q[3]=32;m.a[3]=c;yM(a.f,m);}}}}function Zm(b,c){var d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,A,B,C,D,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X;try{if(b.c){di(b.c);lk(b.c,false);}D=wH(c);if(null==D){return false;}if(null==wH(c)){return false;}if(null==wH(c)){return false;}if(null==(w=wH(c))){return false;}try{F=xI(MJ(w.substr(0,3)));G=xI(MJ(w.substr(3,3)));H=Ym(MJ(w.substr(6,3)));n=Ym(MJ(w.substr(12,3)));T=w.length>=39&&DJ(w.substr(34,5),'V3000')?3:2;}catch(a){a=cG(a);if(PC(a,12)){return false;}else throw dG(a);}if(T==3){K=$m(b,c);pk(b.c,D);return K;}!b.c&&(b.c=new lp(F,G));pk(b.c,D);n==0&&(b.c.J=true);if(0==F){while(w!=null&&!(DJ(w,MR)||DJ(w,NR)||DJ(w.substr(1,w.length-1),'$'))){w=wH(c);}return true;}for(r=0;r<F;r++){if(null==(w=wH(c))){return false;}V=LI(MJ(w.substr(0,10)));W=LI(MJ(w.substr(10,10)));X=LI(MJ(w.substr(20,10)));e=Hh(b.c,V,-W,-X);v=MJ(w.substr(31,3));h=Ck(v);ak(b.c,e,h);DJ(v,'A')&&Vj(b.c,e,1,true);C=Ym(MJ(w.substr(34,2)));C!=0&&Tj(b.c,e,(Gh(),Eh)[h]+C);m=Ym(MJ(w.substr(36,3)));m!=0&&Jj(b.c,e,4-m);A=w.length<63?0:Ym(MJ(w.substr(60,3)));Rj(b.c,e,A,false);p=w.length<45?0:Ym(MJ(w.substr(42,3)));switch(p){case 0:break;case 1:Vj(b.c,e,768,true);break;case 2:Vj(b.c,e,128,true);break;case 3:Vj(b.c,e,384,true);break;default:Vj(b.c,e,896,true);}w.length>=48&&w.charCodeAt(47)==49&&Vj(b.c,e,iR,true);S=w.length<51?0:Ym(MJ(w.substr(48,3)));switch(S){case 0:break;case 15:Hj(b.c,e,0);break;default:Hj(b.c,e,S);}}for(s=0;s<G;s++){if(null==(w=wH(c))){return false;}f=xI(MJ(w.substr(0,3)))-1;g=xI(MJ(w.substr(3,3)))-1;k=xI(MJ(w.substr(6,3)));M=w.length<12?0:Ym(MJ(w.substr(9,3)));Q=w.length<18?0:Ym(MJ(w.substr(15,3)));Jm(b,f,g,k,M,Q);}for(q=0;q<H;q++){if(null==wH(c)){return false;}}if(null==(w=wH(c))){n==0&&Xo(b.c,7);return true;}while(w!=null&&!(DJ(w,MR)||DJ(w,NR))){if(DJ(w.substr(0,6),'M  CHG')){t=xI(MJ(w.substr(6,3)));if(t>0){d=10;U=14;for(u=1;u<=t;++u,d+=8,U+=8){e=xI(MJ(w.substr(d,d+3-d)))-1;l=xI(MJ(w.substr(U,U+3-U)));Jj(b.c,e,l);}}}if(DJ(w.substr(0,6),'M  ISO')){t=xI(MJ(w.substr(6,3)));if(t>0){d=10;U=14;for(u=1;u<=t;++u,d+=8,U+=8){e=xI(MJ(w.substr(d,d+3-d)))-1;B=xI(MJ(w.substr(U,U+3-U)));Tj(b.c,e,B);}}}if(DJ(w.substr(0,6),'M  RAD')){t=xI(MJ(w.substr(6,3)));if(t>0){d=10;U=14;for(u=1;u<=t;++u,d+=8,U+=8){e=xI(MJ(w.substr(d,d+3-d)))-1;J=xI(MJ(w.substr(U,U+3-U)));switch(J){case 1:Wj(b.c,e,16);break;case 2:Wj(b.c,e,32);break;case 3:Wj(b.c,e,48);}}}}if(DJ(w.substr(0,6),'M  RBD')){t=xI(MJ(w.substr(6,3)));if(t>0){d=10;U=14;for(u=1;u<=t;++u,d+=8,U+=8){e=xI(MJ(w.substr(d,d+3-d)))-1;L=xI(MJ(w.substr(U,U+3-U)));switch(L){case 3:case-1:Vj(b.c,e,112,true);break;case 1:Vj(b.c,e,8,true);break;case 2:Vj(b.c,e,104,true);break;case 4:Vj(b.c,e,56,true);}}}}if(DJ(w.substr(0,6),'M  ALS')){e=xI(MJ(w.substr(7,3)))-1;if(e>=0){I=xI(MJ(w.substr(10,3)));i=w.charCodeAt(14)==84;R=ZB(_C,DQ,6,I,15,1);d=16;for(u=0;u<I;++u,d+=4){P=MJ(w.substr(d,d+4-d));R[u]=Ck(P);}Qj(b.c,e,R,i);}}if(DJ(w.substr(0,6),'M  SUB')){t=xI(MJ(w.substr(6,3)));if(t>0){d=10;U=14;for(u=1;u<=t;++u,d+=8,U+=8){e=xI(MJ(w.substr(d,d+3-d)))-1;N=xI(MJ(w.substr(U,U+3-U)));if(N==-2){Vj(b.c,e,YQ,true);}else if(N>0){O=0;for(j=0;j<b.c.p;j++){(Ei(b.c,0,j)==e||Ei(b.c,1,j)==e)&&++O;}N>O&&Vj(b.c,e,xQ,true);}}}}w=wH(c);}}catch(a){a=cG(a);if(PC(a,12)){o=a;vA(o,(cK(),bK),'');return false;}else throw dG(a);}Xo(b.c,7);return true;}function Io(a,b,c,d){var e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,A,B,C,D,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X;a.b=b;di(a.b);K=null;i=ZB(_C,DQ,6,64,15,1);i[0]=-1;N=ZB(_C,DQ,6,64,15,1);O=ZB(_C,DQ,6,64,15,1);for(u=0;u<64;u++)N[u]=-1;M=0;g=0;R=false;L=false;P=false;k=0;Q=c.length;j=1;while(c[M]<=32)++M;while(M<Q){S=c[M++]&AQ;if(WH(S)||S==42){h=0;q=-1;w=false;J=false;v=false;if(R){if(S==82&&VH(c[M]&AQ)){D=VH(c[M+1]&AQ)?2:1;h=Ck(PJ(GP(c,M-1,(B=1+D,DP(),B))));M+=D;}else{A=ZH(c[M]&AQ)==(c[M]&AQ)&&WH(c[M]&AQ)?2:1;h=Ck(PJ(GP(c,M-1,(DP(),A))));M+=A-1;q=0;}if(c[M]==64){++M;if(c[M]==64){v=true;++M;}J=true;}if(c[M]==72){++M;q=1;if(VH(c[M]&AQ)){q=c[M]-48;++M;}}}else if(S==42){h=6;w=true;}else{switch(String.fromCharCode(S).toUpperCase().charCodeAt(0)){case 66:if(M<Q&&c[M]==114){h=35;++M;}else h=5;break;case 67:if(M<Q&&c[M]==108){h=17;++M;}else h=6;break;case 70:h=9;break;case 73:h=53;break;case 78:h=7;break;case 79:h=8;break;case 80:h=15;break;case 83:h=16;}}if(h==0)throw dG(new FA('SmilesParser: unknown element label found'));f=Ih(a.b,h);if(w){P=true;Vj(a.b,f,1,true);}else{Sj(a.b,f,ZH(S)==S&&WH(S));}if(q!=-1&&h!=1){l=ZB(XC,pR,6,1,15,1);l[0]=q<<24>>24;Nj(a.b,f,l);}r=i[k];i[k]!=-1&&j!=128&&Jh(a.b,f,i[k],j);j=1;i[k]=f;if(g!=0){Tj(a.b,f,g);g=0;}if(d){H=!K?null:cM(K,YI(r));!!H&&Oo(H,f,M,h==1);if(J){!K&&(K=new hO());aO(K,YI(f),new Ro(a,f,r,q,M,v));}}continue;}if(S==46){j=128;continue;}if(S==61){j=2;continue;}if(S==35){j=4;continue;}if(VH(S)){F=S-48;if(R){while(M<Q&&VH(c[M]&AQ)){F=10*F+c[M]-48;++M;}g=F;}else{if(L&&M<Q&&VH(c[M]&AQ)){F=10*F+c[M]-48;++M;}L=false;if(F>=64)throw dG(new FA('SmilesParser: ringClosureAtom number out of range'));if(N[F]==-1){N[F]=i[k];O[F]=M-1;}else{if(N[F]===i[k])throw dG(new FA('SmilesParser: ring closure to same atom'));if(d&&!!K){H=cM(K,YI(N[F]));!!H&&Oo(H,i[k],O[F],false);H=cM(K,YI(i[k]));!!H&&Oo(H,N[F],M-1,false);}Jh(a.b,i[k],N[F],j);N[F]=-1;}j=1;}continue;}if(S==43){if(!R)throw dG(new FA("SmilesParser: '+' found outside brackets"));m=1;while(c[M]==43){++m;++M;}if(m==1&&VH(c[M]&AQ)){m=c[M]-48;++M;}Jj(a.b,i[k],m);continue;}if(S==45){if(!R)continue;m=-1;while(c[M]==45){--m;++M;}if(m==-1&&VH(c[M]&AQ)){m=48-c[M];++M;}Jj(a.b,i[k],m);continue;}if(S==40){if(i[k]==-1)throw dG(new FA('Smiles with leading parenthesis are not supported'));i[k+1]=i[k];++k;continue;}if(S==41){--k;continue;}if(S==91){if(R)throw dG(new FA('SmilesParser: nested square brackets found'));R=true;continue;}if(S==93){if(!R)throw dG(new FA('SmilesParser: closing bracket without opening one'));R=false;continue;}if(S==37){L=true;continue;}if(S==58){if(!R){j=64;continue;}C=0;while(VH(c[M]&AQ)){C=10*C+c[M]-48;++M;}Rj(a.b,i[k],C,false);continue;}if(S==47){d&&(j=17);continue;}if(S==92){d&&(j=9);continue;}throw dG(new FA("SmilesParser: unexpected character found: '"+String.fromCharCode(S)+"'"));}if(j!=1)throw dG(new FA('SmilesParser: dangling open bond'));for(t=0;t<64;t++)if(N[t]!=-1)throw dG(new FA('SmilesParser: dangling ring closure'));s=kl(a.b);a.b.P=true;Xo(a.b,1);for(e=0;e<a.b.o;e++){if((b.r==null?null:b.r[e]==null?null:CJ(b.r[e]))!=null){if(!lj(a.b,e)){p=mi(a.b,e)[0];if(Ai(a.b,e)<(Gh(),Dh).length&&Dh[Ai(a.b,e)]!=null){n=false;T=ol(a.b,e)-Si(a.b,e);for(V=Dh[Ai(a.b,e)],W=0,X=V.length;W<X;++W){U=V[W];if(T<=U){n=true;U!=T+p&&Hj(a.b,e,T+p);break;}}n||Hj(a.b,e,T+p);}}}}Go(a);Ho(a);a.b.r=null;a.b.P=false;d&&Mo(a)&&Xl(a.b,0);if(d){Ig(new Tg(),a.b);if(K){for(I=(G=new qO(new vO(new ML(K).a).b),new RL(G));RK(I.a.a);){H=(o=oO(I.a),o.Gb());Uj(a.b,H.a,Po(H,s),false);}Xl(a.b,0);}$l(a.b);hp(a.b);}P&&lk(a.b,true);}function df(a){var b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,A,B,C,D,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W;Pe(a);Ne(a,9,4);S=lJ(ff(a.L.d),ff(a.L.e));Ne(a,S,4);if(S==0){Ne(a,a.L.I?1:0,1);Ne(a,0,1);a.D=Oe(a);return;}T=V=U=F=0;for(f=0;f<a.L.d;f++){if((vi(a.L,f)&1)==0){switch(Ai(a.L,f)){case 6:break;case 7:++T;break;case 8:++V;break;default:++U;}ji(a.L,f)!=0&&++F;}}Ne(a,a.L.d,S);Ne(a,a.L.e,S);Ne(a,T,S);Ne(a,V,S);Ne(a,U,S);Ne(a,F,S);for(g=0;g<a.L.d;g++)Ai(a.L,a.t[g])==7&&(vi(a.L,a.t[g])&1)==0&&Ne(a,g,S);for(l=0;l<a.L.d;l++)Ai(a.L,a.t[l])==8&&(vi(a.L,a.t[l])&1)==0&&Ne(a,l,S);for(m=0;m<a.L.d;m++)if(Ai(a.L,a.t[m])!=6&&Ai(a.L,a.t[m])!=7&&Ai(a.L,a.t[m])!=8&&(vi(a.L,a.t[m])&1)==0){Ne(a,m,S);Ne(a,Ai(a.L,a.t[m]),8);}for(n=0;n<a.L.d;n++)if(ji(a.L,a.t[n])!=0&&(vi(a.L,a.t[n])&1)==0){Ne(a,n,S);Ne(a,8+ji(a.L,a.t[n]),4);}R=0;u=0;for(o=1;o<a.L.d;o++){if(a.w[o]==-1){J=0;}else{J=1+a.w[o]-u;u=a.w[o];}R<J&&(R=J);}I=ff(R);Ne(a,I,4);u=0;for(p=1;p<a.L.d;p++){if(a.w[p]==-1){J=0;}else{J=1+a.w[p]-u;u=a.w[p];}Ne(a,J,I);}for(L=0;L<2*a.C;L++)Ne(a,a.v[L],S);for(w=0;w<a.L.e;w++){D=(Oi(a.L,w)&$Q)!=0?1:Hl(a.L,a.u[w])?0:Mi(a.L,a.u[w]);Ne(a,D,2);}c=0;for(q=0;q<a.L.d;q++)a.S[a.t[q]]!=0&&a.S[a.t[q]]!=3&&++c;Ne(a,c,S);for(r=0;r<a.L.d;r++)if(a.S[a.t[r]]!=0&&a.S[a.t[r]]!=3){Ne(a,r,S);if(a.U[a.t[r]]==0){Ne(a,a.S[a.t[r]],3);}else{W=a.S[a.t[r]]==1?a.U[a.t[r]]==1?4:6:a.U[a.t[r]]==1?5:7;Ne(a,W,3);Ne(a,a.T[a.t[r]],3);}}b=0;for(A=0;A<a.L.e;A++)a.g[a.u[A]]!=0&&a.g[a.u[A]]!=3&&(!Ol(a.L,a.u[A])||Pi(a.L,a.u[A])==1)&&++b;Ne(a,b,S);for(B=0;B<a.L.e;B++)if(a.g[a.u[B]]!=0&&a.g[a.u[B]]!=3&&(!Ol(a.L,a.u[B])||Pi(a.L,a.u[B])==1)){Ne(a,B,S);if(Pi(a.L,a.u[B])==1){if(a.j[a.u[B]]==0){Ne(a,a.g[a.u[B]],3);}else{W=a.g[a.u[B]]==1?a.j[a.u[B]]==1?4:6:a.j[a.u[B]]==1?5:7;Ne(a,W,3);Ne(a,a.i[a.u[B]],3);}}else{Ne(a,a.g[a.u[B]],2);}}Ne(a,a.L.I?1:0,1);G=0;for(s=0;s<a.L.d;s++)ti(a.L,a.t[s])!=0&&++G;if(G!=0){Ne(a,1,1);Ne(a,1,4);Ne(a,G,S);for(h=0;h<a.L.d;h++){if(ti(a.L,a.t[h])!=0){Ne(a,h,S);Ne(a,ti(a.L,a.t[h]),8);}}}N=false;if(a.L.I){ae(a,0,false,S,YQ,1,-1);be(a,2,false,S,64,1,-1);ae(a,3,false,S,xQ,1,-1);ae(a,4,false,S,120,4,3);ae(a,5,false,S,6,2,1);ae(a,6,false,S,1,1,-1);ae(a,7,false,S,1920,4,7);G=0;for(h=0;h<a.L.d;h++)qi(a.L,a.t[h])!=null&&++G;if(G>0){Ne(a,1,1);Ne(a,8,4);Ne(a,G,S);for(i=0;i<a.L.d;i++){t=qi(a.L,a.t[i]);if(t!=null){Ne(a,i,S);Ne(a,t.length,4);for(K=0;K<t.length;K++)Ne(a,t[K],8);}}}be(a,9,false,S,48,2,4);be(a,10,false,S,15,4,0);ae(a,11,false,S,iR,1,-1);be(a,12,false,S,$Q,8,6);ae(a,13,false,S,SQ,3,14);ae(a,14,false,S,TQ,5,17);N=N|ae(a,16,false,S,WQ,3,22);}G=0;for(j=0;j<a.L.d;j++)a.a!=null&&a.a[a.t[j]]!=-1&&++G;if(G!=0){N=Se(a,N);Ne(a,1,1);Ne(a,1,4);Ne(a,G,S);for(h=0;h<a.L.d;h++){if(a.a!=null&&a.a[a.t[h]]!=-1){Ne(a,h,S);Ne(a,a.a[a.t[h]],4);}}}if((a.K&8)!=0){G=0;Q=0;for(h=0;h<a.L.d;h++){O=li(a.L,a.t[h]);if(O!=null){++G;Q=lJ(Q,O.length);}}if(G!=0){N=Se(a,N);P=ff(Q);Ne(a,1,1);Ne(a,2,4);Ne(a,G,S);Ne(a,P,4);for(i=0;i<a.L.d;i++){H=li(a.L,a.t[i]);if(H!=null){Ne(a,i,S);Ne(a,H.length,P);for(K=0;K<H.length;K++)Ne(a,H.charCodeAt(K),7);}}}}if(a.L.I){N=N|ae(a,19,N,S,PQ,3,25);N=N|be(a,20,N,S,SQ,3,14);}G=0;for(k=0;k<a.L.d;k++)wi(a.L,a.t[k])!=0&&++G;if(G!=0){N=Se(a,N);Ne(a,1,1);Ne(a,5,4);Ne(a,G,S);for(e=0;e<a.L.d;e++){if(wi(a.L,a.t[e])!=0){Ne(a,e,S);Ne(a,wi(a.L,a.t[e])>>4,2);}}}if(a.L.I){N=N|ae(a,22,N,S,XQ,1,-1);N=N|be(a,23,N,S,qR,1,-1);N=N|be(a,24,N,S,bR,2,18);}if((a.K&16)!=0){for(e=0;e<a.L.d;e++){if(qj(a.L,a.t[e])){N=Se(a,N);Ne(a,1,1);Ne(a,9,4);for(d=0;d<a.L.d;d++)Ne(a,qj(a.L,a.t[d])?1:0,1);break;}}}M=Ve(a);if(M!=null){G=0;for(C=0;C<a.L.e;C++)M[a.u[C]]&&++G;N=Se(a,N);Ne(a,1,1);Ne(a,10,4);Ne(a,G,S);for(v=0;v<a.L.e;v++)M[a.u[v]]&&Ne(a,v,S);}a.L.I&&(N=N|ae(a,27,N,S,IQ,1,-1));Ne(a,0,1);a.D=Oe(a);}function Im(a,b){var c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,A,B,C,D,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T;this.a=new nK('0.0000');Xo(a,7);H=true;for(d=0;d<a.d;d++){if((a.s[d]&3)!=0&&(a.s[d]&3)!=3&&(a.s[d]&zR)>>19!=1){H=false;break;}}J=-1;if(H){A=ZB(_C,DQ,6,32,15,1);for(e=0;e<a.d;e++){if((a.s[e]&3)!=0&&(a.s[e]&3)!=3&&(a.s[e]&zR)>>19==1){C=(a.s[e]&zR)>>19!=1&&(a.s[e]&zR)>>19!=2?-1:(a.s[e]&AR)>>21;++A[C];0<A[C]&&(J=C);break;}}}this.b=b;L=a.M!=null?a.M:'';WJ(this.b,L+kQ);WJ(this.b,'Actelion Java MolfileCreator 1.0\n\n');Gm(this,a.o);Gm(this,a.p);WJ(this.b,'  0  0');Gm(this,H?0:1);WJ(this.b,'  0  0  0  0  0999 V2000\n');D=a.o==1;for(g=1;g<a.o;g++){if(a.H[g].a!=a.H[0].a||a.H[g].b!=a.H[0].b||a.H[g].c!=a.H[0].c){D=true;break;}}B=1;if(D){p=Ci(a,a.o,a.p,(Gh(),Fh));if(p!=0){(p<1||p>3)&&(B=1.5/p);}else{K=DR;for(e=1;e<a.o;e++){for(f=0;f<e;f++){u=a.H[f].a-a.H[e].a;v=a.H[f].b-a.H[e].b;w=a.H[f].c-a.H[e].c;t=u*u+v*v+w*w;K>t&&(K=t);}}B=3/K;}}for(h=0;h<a.o;h++){if(D){Fm(this,B*a.H[h].a);Fm(this,B*-a.H[h].b);Fm(this,B*-a.H[h].c);}else{WJ(this.b,'    0.0000    0.0000    0.0000');}if((a.t==null?null:a.t[h])!=null)WJ(this.b,' L  ');else if((a.w[h]&1)!=0)WJ(this.b,' A  ');else{n=(Gh(),Ch)[a.A[h]];WJ(this.b,' '+n);n.length==1?WJ(this.b,'  '):n.length==2&&WJ(this.b,' ');}WJ(this.b,' 0  0  0');F=1920&a.w[h];F==0?WJ(this.b,'  0'):F==384?WJ(this.b,'  3'):F==128?WJ(this.b,'  2'):F==1792?WJ(this.b,'  1'):F==1664&&WJ(this.b,'  2');WJ(this.b,(a.w[h]&iR)!=0?'  1':'  0');T=((a.s[h]&yR)>>>28)-1;T==-1?WJ(this.b,'  0'):T==0?WJ(this.b,' 15'):Gm(this,T);WJ(this.b,'  0  0  0');Gm(this,kJ(a.u[h]));WJ(this.b,'  0  0\n');}for(q=0;q<a.p;q++){switch(a.F[q]){case 1:N=1;Q=0;break;case 2:N=2;Q=0;break;case 4:N=3;Q=0;break;case 9:N=1;Q=6;break;case 17:N=1;Q=1;break;case 26:N=2;Q=3;break;case 64:N=4;Q=0;break;default:N=1;Q=0;}H&&(Q==1||Q==6)&&ni(a,a.B[0][q])!=J&&(Q=0);r=a.D[q]&15;r!=0&&(r==8?N=4:r==3?N=5:r==9?N=6:r==10?N=7:N=8);P=a.D[q]&48;S=P==0?0:P==32?1:2;Gm(this,1+a.B[0][q]);Gm(this,1+a.B[1][q]);Gm(this,N);Gm(this,Q);WJ(this.b,'  0');Gm(this,S);WJ(this.b,'  0\n');}M=0;for(i=0;i<a.o;i++)a.q[i]!=0&&++M;if(M!=0){WJ(this.b,'M  CHG');Gm(this,M);for(e=0;e<a.o;e++){if(a.q[e]!=0){WJ(this.b,' ');Gm(this,e+1);s=a.q[e];if(s<0){WJ(this.b,'  -');s=-s;}else WJ(this.b,'   ');TJ(this.b,48+s&AQ);}}WJ(this.b,kQ);}M=0;for(j=0;j<a.o;j++)a.v[j]==0||++M;if(M!=0){WJ(this.b,'M  ISO');Gm(this,M);for(e=0;e<a.o;e++){if(a.v[e]!=0){WJ(this.b,' ');Gm(this,e+1);WJ(this.b,' ');Gm(this,a.v[e]);}}WJ(this.b,kQ);}M=0;for(k=0;k<a.o;k++)(a.s[k]&48)!=0&&++M;if(M!=0){WJ(this.b,'M  RAD');Gm(this,M);for(c=0;c<a.o;c++){if((a.s[c]&48)!=0){WJ(this.b,' ');Gm(this,c+1);switch(a.s[c]&48){case 16:WJ(this.b,'   1');break;case 32:WJ(this.b,'   2');break;case 48:WJ(this.b,'   3');}}}WJ(this.b,kQ);}if(a.I){M=0;for(e=0;e<a.o;e++)(a.w[e]&120)!=0&&++M;if(M!=0){WJ(this.b,'M  RBD');Gm(this,M);for(f=0;f<a.o;f++){O=a.w[f]&120;if(O!=0){WJ(this.b,' ');Gm(this,f+1);switch(O){case 112:WJ(this.b,'  -1');break;case 8:WJ(this.b,'   1');break;case 104:WJ(this.b,'   2');break;case 88:WJ(this.b,'   3');break;case 56:WJ(this.b,'   4');}}}WJ(this.b,kQ);}for(l=0;l<a.o;l++){o=a.t==null?null:a.t[l];if(o!=null){WJ(this.b,'M  ALS ');Gm(this,l+1);Gm(this,o.length);WJ(this.b,(a.w[l]&1)!=0?' T ':' F ');for(G=0;G<o.length;G++){I=(Gh(),Ch)[o[G]];switch(I.length){case 1:WJ(this.b,I+'   ');break;case 2:WJ(this.b,I+'  ');break;case 3:WJ(this.b,I+' ');break;default:WJ(this.b,'   ?');}}WJ(this.b,kQ);}}M=0;for(m=0;m<a.o;m++)(a.w[m]&6144)!=0&&++M;if(M!=0){WJ(this.b,'M  SUB');Gm(this,M);for(c=0;c<a.o;c++){R=a.w[c]&6144;if(R!=0){WJ(this.b,' ');Gm(this,c+1);(R&xQ)!=0?WJ(this.b,'   '+(a.c[c]+1)):WJ(this.b,'  -2');}}WJ(this.b,kQ);}}WJ(this.b,'M  END\n');}function gd(a,b,c){var d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,A,B,C,D,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X;c&&ro(a,b,(pi(a.G,b),uh(a.K,xi(a.G,b))),vh(a.K,yi(a.G,b)));J=null;if(ji(a.G,b)!=0){R=kJ(ji(a.G,b))==1?'':''+kJ(ji(a.G,b));J=ji(a.G,b)<0?R+'-':R+'+';}B=null;K=vi(a.G,b);if(K!=0){(K&2)!=0&&(B='a');(K&4)!=0&&(B=B==null?'!a':B+','+'!a');(K&xQ)!=0&&(B=B==null?'s':B+','+'s');if((K&1920)!=0){u=K&1920;u==1792?B=B==null?'h0':B+','+'h0':u==1664?B=B==null?'h1':B+','+'h1':u==1408?B=B==null?'h2':B+','+'h2':u==128?B=B==null?'h>0':B+','+'h>0':u==384?B=B==null?'h>1':B+','+'h>1':u==OQ?B=B==null?'h<3':B+','+'h<3':u==1536&&(B=B==null?'h<2':B+','+'h<2');}if((K&PQ)!=0){i=K&PQ;i==167772160?B=B==null?'c0':B+','+'c0':i==QQ?B=B==null?'c+':B+','+'c+':i==RQ&&(B=B==null?'c-':B+','+'c-');}if((K&SQ)!=0){I=K&SQ;I==98304?B=B==null?'pi0':B+','+'pi0':I==81920?B=B==null?'pi1':B+','+'pi1':I==49152?B=B==null?'pi2':B+','+'pi2':I==yQ&&(B=B==null?'pi>0':B+','+'pi>0');}if((K&TQ)!=0){H=K&TQ;H==3801088?B=B==null?'n1':B+','+'n1':H==3538944?B=B==null?'n2':B+','+'n2':H==3014656?B=B==null?'n3':B+','+'n3':H==3145728?B=B==null?'n<3':B+','+'n<3':H==UQ?B=B==null?'n<4':B+','+'n<4':H==VQ?B=B==null?'n>1':B+','+'n>1':H==917504?B=B==null?'n>2':B+','+'n>2':H==1966080&&(B=B==null?'n>3':B+','+'n>3');}if((K&120)!=0){N=K&120;N==112?B=B==null?'c':B+','+'c':N==8?B=B==null?'r':B+','+'r':N==104?B=B==null?'rb2':B+','+'rb2':N==88?B=B==null?'rb3':B+','+'rb3':N==56&&(B=B==null?'rb4':B+','+'rb4');}(K&WQ)!=0&&(B=B==null?'rs'+((K&WQ)>>22):B+','+('rs'+((K&WQ)>>22)));(K&XQ)!=0&&(B=B==null?'sp2':B+','+'sp2');}ti(a.G,b)!=0&&(B=Kc(B,''+ti(a.G,b)));Q=0;if(wi(a.G,b)!=0){switch(wi(a.G,b)){case 16:J=J==null?'|':J+','+'|';break;case 32:Q=1;break;case 48:Q=2;}}l=null;if((a.B&64)==0){if($i(a.G,b))l='?';else if(ii(a.G,b)!=0){if(bl(a.G,b)==2){switch(ii(a.G,b)){case 2:l=aj(a.G,b)?'p':'P';break;case 1:l=aj(a.G,b)?'m':'M';break;default:l='*';}}else{switch(ii(a.G,b)){case 1:l=aj(a.G,b)?'r':'R';break;case 2:l=aj(a.G,b)?'s':'S';break;default:l='*';}}}}(a.B&1792)!=0&&(l=Kc(l,''+fp(a.G,b)));F=null;(a.B&16)!=0&&si(a.G,b)!=0&&(F=''+si(a.G,b));p=null;if(wl(a.G,b)!=-1){o=Wc(a,b);o!=-1&&(p=o==0?'abs':((o&255)==1?'&':'or')+(1+(o>>8)));}v=0;(Ai(a.G,b)!=6||!a.p[b]||(vi(a.G,b)&YQ)!=0||wi(a.G,b)!=0)&&(v=ml(a.G,b));f=li(a.G,b);if(f!=null){v=0;}else if(qi(a.G,b)!=null){e=(vi(a.G,b)&1)!=0?'[!':'[';f=e+ri(a.G,b)+']';f.length>5&&(f=e+qi(a.G,b).length+']');(vi(a.G,b)&YQ)!=0&&(v=-1);}else if((vi(a.G,b)&1)!=0){f='?';(vi(a.G,b)&YQ)!=0&&(v=-1);}else(Ai(a.G,b)!=6||J!=null||B!=null||v>0||!a.p[b])&&(f=pi(a.G,b));D=0;!qj(a.G,b)&(vi(a.G,b)&IQ)!=0&&zd(a,-8);if(f!=null){D=(L=(S=dH(a.e,f),new tH(0,0,S,0)).b,L);ld(a,uh(a.K,xi(a.G,b)),vh(a.K,yi(a.G,b)),f,c,true);a.q[b]=true;}else ad(a,b)&&kd(a,uh(a.K,xi(a.G,b)),vh(a.K,yi(a.G,b)),b,c);if(J!=null){vo(a,(a.Q*2+1)/3|0);U=uh(a.K,xi(a.G,b))+((D+(L=(S=dH(a.e,J),new tH(0,0,S,0)).b,L))/2+1);W=vh(a.K,yi(a.G,b))-((a.j*4-4)/8|0);ld(a,U,W,J,c,true);vo(a,a.Q);}(a.B&2)!=0&&(B=''+b);if(B!=null){vo(a,(a.Q*2+1)/3|0);U=uh(a.K,xi(a.G,b))-(D+(L=(S=dH(a.e,B),new tH(0,0,S,0)).b,L))/2;W=vh(a.K,yi(a.G,b))-((a.j*4-4)/8|0);ld(a,U,W,B,c,true);vo(a,a.Q);}if(l!=null){vo(a,(a.Q*2+1)/3|0);U=uh(a.K,xi(a.G,b))-(D+(L=(S=dH(a.e,l),new tH(0,0,S,0)).b,L))/2;W=vh(a.K,yi(a.G,b))+((a.j*4+4)/8|0);P=a.w;zd(a,448);ld(a,U,W,l,c,false);zd(a,P);vo(a,a.Q);}if(F!=null){vo(a,(a.Q*2+1)/3|0);U=uh(a.K,xi(a.G,b))+((D+(L=(S=dH(a.e,F),new tH(0,0,S,0)).b,L))/2+1);W=vh(a.K,yi(a.G,b))+((a.j*4+4)/8|0);P=a.w;zd(a,cj(a.G,b)?384:448);ld(a,U,W,F,c,true);zd(a,P);vo(a,a.Q);}if(p!=null){d=pd(a,b);vo(a,(a.Q*2+1)/3|0);U=uh(a.K,xi(a.G,b))+0.7*a.j*$wnd.Math.sin(d);W=vh(a.K,yi(a.G,b))+0.7*a.j*$wnd.Math.cos(d);P=a.w;zd(a,Vc(a,b));ld(a,U,W,p,c,false);zd(a,P);vo(a,a.Q);}if(v==0&&Q==0){a.w==-8&&zd(a,-9);return;}s=ZB(ZC,GQ,6,4,15,1);for(A=0;A<Qk(a.G,b);A++){h=cl(a.G,b,A);for(C=0;C<2;C++){if(Ei(a.G,C,h)==b){O=Di(a.G,Ei(a.G,C,h),Ei(a.G,1-C,h));if(O<ZQ){s[0]-=O+NQ;s[3]+=O+MQ;}else if(O<0){s[2]+=O+NQ;s[3]-=O;}else if(O<NQ){s[1]+=O;s[2]+=NQ-O;}else{s[0]+=O-NQ;s[1]+=MQ-O;}}}}bl(a.G,b)==0?jj(a.G,b)?s[3]-=0.2:s[1]-=0.2:s[1]-=0.1;(J!=null||F!=null)&&(s[1]+=10);(B!=null||l!=null)&&(s[3]+=10);q='';if(v!=0){t=(M=(T=dH(a.e,'H'),new tH(0,0,T,0)).b,M);r=0;if(v==-1){q='n';vo(a,(a.Q*2+1)/3|0);r=(L=(S=dH(a.e,'n'),new tH(0,0,S,0)).b,L);}else if(v>1){q=''+v;vo(a,(a.Q*2+1)/3|0);r=(L=(S=dH(a.e,q),new tH(0,0,S,0)).b,L);}if(s[1]<0.6||s[3]<0.6){k=vh(a.K,yi(a.G,b));if(s[1]<=s[3]){s[1]+=10;j=uh(a.K,xi(a.G,b))+(D+t)/2;}else{s[3]+=10;j=uh(a.K,xi(a.G,b))-(D+t)/2-r;}}else{j=uh(a.K,xi(a.G,b));if(s[0]<s[2]){s[0]+=10;k=vh(a.K,yi(a.G,b))-a.j;}else{s[2]+=10;k=vh(a.K,yi(a.G,b))+a.j;}}if(r>0){U=j+(t+r)/2;W=k+((a.j*4+4)/8|0);ld(a,U,W,q,c,true);vo(a,a.Q);}ld(a,j,k,'H',c,true);}g=0;if(Q!=0){G=50;m=0;for(w=0;w<4;w++){n=w>1?w-2:w+2;if(s[w]<G){g=w;G=s[w];m=s[n];}else if(s[w]==G){if(s[n]>m){g=w;m=s[n];}}}switch(g){case 0:j=uh(a.K,xi(a.G,b));k=vh(a.K,yi(a.G,b))-a.O-D/2;break;case 1:j=uh(a.K,xi(a.G,b))+a.O+D/2;k=vh(a.K,yi(a.G,b));break;case 2:j=uh(a.K,xi(a.G,b));k=vh(a.K,yi(a.G,b))+a.O+D/2;break;default:j=uh(a.K,xi(a.G,b))-a.O-D/2;k=vh(a.K,yi(a.G,b));}if(Q==1){yM(a.T,new tH(j-a.O,k-a.O,2*a.O,2*a.O));c&&yM(a.N,new Ed(j,k,$c(a,b)?-3:a.o[b]));}else{switch(g){case 2:case 0:V=2*a.O;X=0;j-=a.O;break;case 1:V=0;X=2*a.O;k-=a.O;break;default:V=0;X=2*a.O;k-=a.O;}yM(a.T,new tH(j-a.O,k-a.O,2*a.O,2*a.O));c&&yM(a.N,new Ed(j,k,$c(a,b)?-3:a.o[b]));yM(a.T,new tH(j+V-a.O,k+X-a.O,2*a.O,2*a.O));c&&yM(a.N,new Ed(j+V,k+X,$c(a,b)?-3:a.o[b]));}}a.w==-8&&zd(a,-9);}function nm(b,c,d,e){var f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,A,B,C,D,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,$,ab,bb,cb,db,eb,fb,gb,hb,ib,jb,kb,lb,mb,nb,ob,pb,qb,rb,sb,tb,ub,vb,wb,xb,yb,zb,Ab,Bb,Cb,Db,Eb,Fb,Gb,Hb,Ib,Jb,Kb,Lb,Mb,Nb,Ob,Pb,Qb,Rb,Sb,Tb,Ub,Vb,Wb,Xb,Yb,Zb,$b,_b,ac,bc,cc,dc,ec,fc,gc,hc,ic,jc,kc,lc;fc=8;b.f=c;di(b.f);if(d==null||d.length==0)return;e!=null&&e.length==0&&(e=null);hm(b,d,0);h=gm(b,4);w=gm(b,4);if(h>8){fc=h;h=w;}if(h==0){lk(b.f,gm(b,1)==1);return;}i=gm(b,h);j=gm(b,w);Wb=gm(b,h);$b=gm(b,h);Zb=gm(b,h);K=gm(b,h);for(n=0;n<i;n++)Ih(b.f,6);for(fb=0;fb<Wb;fb++)ak(b.f,gm(b,h),7);for(gb=0;gb<$b;gb++)ak(b.f,gm(b,h),8);for(rb=0;rb<Zb;rb++)ak(b.f,gm(b,h),gm(b,8));for(Cb=0;Cb<K;Cb++)Jj(b.f,gm(b,h),gm(b,4)-8);L=1+j-i;S=gm(b,4);v=0;Zj(b.f,0,0);$j(b.f,0,0);_j(b.f,0,0);T=e!=null&&e[0]>=39;ec=0;hc=0;jc=0;lc=0;O=false;P=false;if(T){if(e.length>2*i-2&&e[2*i-2]==39||e.length>3*i-3&&e[3*i-3]==39){P=true;O=e.length==3*i-3+9;Nb=O?3*i-3:2*i-2;u=86*(e[Nb+1]-40)+e[Nb+2]-40;ec=$wnd.Math.pow(10,u/2000-1);Nb+=2;gc=86*(e[Nb+1]-40)+e[Nb+2]-40;hc=$wnd.Math.pow(10,gc/1500-1);Nb+=2;ic=86*(e[Nb+1]-40)+e[Nb+2]-40;jc=$wnd.Math.pow(10,ic/1500-1);if(O){Nb+=2;kc=86*(e[Nb+1]-40)+e[Nb+2]-40;lc=$wnd.Math.pow(10,kc/1500-1);}}else{O=e.length==3*i-3;}}if(b.b&&O){e=null;T=false;}for(Hb=1;Hb<i;Hb++){U=gm(b,S);if(U==0){if(T){Zj(b.f,Hb,b.f.H[0].a+8*(e[Hb*2-2]-83));$j(b.f,Hb,b.f.H[0].b+8*(e[Hb*2-1]-83));O&&_j(b.f,Hb,b.f.H[0].c+8*(e[2*i-3+Hb]-83));}++L;continue;}v+=U-1;if(T){Zj(b.f,Hb,xi(b.f,v)+e[Hb*2-2]-83);$j(b.f,Hb,yi(b.f,v)+e[Hb*2-1]-83);O&&_j(b.f,Hb,zi(b.f,v)+(e[2*i-3+Hb]-83));}Jh(b.f,v,Hb,1);}for(Ib=0;Ib<L;Ib++)Jh(b.f,gm(b,h),gm(b,h),1);Pb=ZB(aG,HQ,6,j,16,1);for(F=0;F<j;F++){G=gm(b,2);switch(G){case 0:mj(b.f,Ei(b.f,0,F))||mj(b.f,Ei(b.f,1,F))?jk(b.f,F,32):Pb[F]=true;break;case 2:jk(b.f,F,2);break;case 3:jk(b.f,F,4);}}g=gm(b,h);for(Jb=0;Jb<g;Jb++){m=gm(b,h);if(fc==8){_b=gm(b,2);if(_b==3){Oj(b.f,m,1,0);Uj(b.f,m,1,false);}else{Uj(b.f,m,_b,false);}}else{_b=gm(b,3);switch(_b){case 4:Uj(b.f,m,1,false);Oj(b.f,m,1,gm(b,3));break;case 5:Uj(b.f,m,2,false);Oj(b.f,m,1,gm(b,3));break;case 6:Uj(b.f,m,1,false);Oj(b.f,m,2,gm(b,3));break;case 7:Uj(b.f,m,2,false);Oj(b.f,m,2,gm(b,3));break;default:Uj(b.f,m,_b,false);}}}fc==8&&gm(b,1)==0&&(b.f.J=true);f=gm(b,w);for(Kb=0;Kb<f;Kb++){B=gm(b,w);if(Pi(b.f,B)==1){_b=gm(b,3);switch(_b){case 4:gk(b.f,B,1,false);ek(b.f,B,1,gm(b,3));break;case 5:gk(b.f,B,2,false);ek(b.f,B,1,gm(b,3));break;case 6:gk(b.f,B,1,false);ek(b.f,B,2,gm(b,3));break;case 7:gk(b.f,B,2,false);ek(b.f,B,2,gm(b,3));break;default:gk(b.f,B,_b,false);}}else{gk(b.f,B,gm(b,2),false);}}lk(b.f,gm(b,1)==1);l=null;Yb=0;while(gm(b,1)==1){R=Yb+gm(b,4);switch(R){case 0:Xb=gm(b,h);for(Lb=0;Lb<Xb;Lb++){m=gm(b,h);Vj(b.f,m,YQ,true);}break;case 1:Xb=gm(b,h);for(Mb=0;Mb<Xb;Mb++){m=gm(b,h);Ub=gm(b,8);Tj(b.f,m,Ub);}break;case 2:Xb=gm(b,w);for(hb=0;hb<Xb;hb++){B=gm(b,w);jk(b.f,B,64);}break;case 3:Xb=gm(b,h);for(ib=0;ib<Xb;ib++){m=gm(b,h);Vj(b.f,m,xQ,true);}break;case 4:Xb=gm(b,h);for(jb=0;jb<Xb;jb++){m=gm(b,h);dc=gm(b,4)<<3;Vj(b.f,m,dc,true);}break;case 5:Xb=gm(b,h);for(kb=0;kb<Xb;kb++){m=gm(b,h);k=gm(b,2)<<1;Vj(b.f,m,k,true);}break;case 6:Xb=gm(b,h);for(lb=0;lb<Xb;lb++){m=gm(b,h);Vj(b.f,m,1,true);}break;case 7:Xb=gm(b,h);for(mb=0;mb<Xb;mb++){m=gm(b,h);cb=gm(b,4)<<7;Vj(b.f,m,cb,true);}break;case 8:Xb=gm(b,h);for(nb=0;nb<Xb;nb++){m=gm(b,h);s=gm(b,4);q=ZB(_C,DQ,6,s,15,1);for(Qb=0;Qb<s;Qb++){r=gm(b,8);q[Qb]=r;}Pj(b.f,m,q);}break;case 9:Xb=gm(b,w);for(ob=0;ob<Xb;ob++){B=gm(b,w);dc=gm(b,2)<<4;ik(b.f,B,dc,true);}break;case 10:Xb=gm(b,w);for(pb=0;pb<Xb;pb++){B=gm(b,w);H=gm(b,4);ik(b.f,B,H,true);}break;case 11:Xb=gm(b,h);for(qb=0;qb<Xb;qb++){m=gm(b,h);Vj(b.f,m,iR,true);}break;case 12:Xb=gm(b,w);for(sb=0;sb<Xb;sb++){B=gm(b,w);I=gm(b,8)<<6;ik(b.f,B,I,true);}break;case 13:Xb=gm(b,h);for(tb=0;tb<Xb;tb++){m=gm(b,h);ac=gm(b,3)<<14;Vj(b.f,m,ac,true);}break;case 14:Xb=gm(b,h);for(ub=0;ub<Xb;ub++){m=gm(b,h);Vb=gm(b,5)<<17;Vj(b.f,m,Vb,true);}break;case 15:Yb=16;break;case 16:Xb=gm(b,h);for(vb=0;vb<Xb;vb++){m=gm(b,h);cc=gm(b,3)<<22;Vj(b.f,m,cc,true);}break;case 17:Xb=gm(b,h);for(wb=0;wb<Xb;wb++){m=gm(b,h);Hj(b.f,m,gm(b,4));}break;case 18:Xb=gm(b,h);Tb=gm(b,4);for(xb=0;xb<Xb;xb++){m=gm(b,h);Q=gm(b,Tb);Rb=ZB(XC,pR,6,Q,15,1);for(Qb=0;Qb<Q;Qb++)Rb[Qb]=gm(b,7)<<24>>24;Mj(b.f,m,PJ(GP(Rb,0,(Sb=Rb.length,DP(),Sb))));}break;case 19:Xb=gm(b,h);for(yb=0;yb<Xb;yb++){m=gm(b,h);J=gm(b,3)<<25;Vj(b.f,m,J,true);}break;case 20:Xb=gm(b,w);for(zb=0;zb<Xb;zb++){B=gm(b,w);cc=gm(b,3)<<14;ik(b.f,B,cc,true);}break;case 21:Xb=gm(b,h);for(Ab=0;Ab<Xb;Ab++){m=gm(b,h);Wj(b.f,m,gm(b,2)<<4);}break;case 22:Xb=gm(b,h);for(Bb=0;Bb<Xb;Bb++){m=gm(b,h);Vj(b.f,m,XQ,true);}break;case 23:Xb=gm(b,w);for(Db=0;Db<Xb;Db++){B=gm(b,w);ik(b.f,B,qR,true);}break;case 24:Xb=gm(b,w);for(Eb=0;Eb<Xb;Eb++){B=gm(b,w);k=gm(b,2)<<18;ik(b.f,B,k,true);}break;case 25:for(Fb=0;Fb<i;Fb++)gm(b,1)==1&&Xj(b.f,Fb,true);break;case 26:Xb=gm(b,w);l=ZB(_C,DQ,6,Xb,15,1);for(Gb=0;Gb<Xb;Gb++)l[Gb]=gm(b,w);break;case 27:Xb=gm(b,h);for(eb=0;eb<Xb;eb++){m=gm(b,h);Vj(b.f,m,IQ,true);}}}Id(new Qd(b.f,Pb));if(l!=null)for(C=0,D=l.length;C<D;++C){B=l[C];jk(b.f,B,Pi(b.f,B)==2?4:2);}M=0;if(e==null&&d.length>b.d+1&&(d[b.d+1]==32||d[b.d+1]==9)){e=d;M=b.d+2;}if(e!=null){try{if(e[M]==33||e[M]==35){hm(b,e,M+1);O=gm(b,1)==1;P=gm(b,1)==1;bc=2*gm(b,4);A=1<<bc;B=0;for(o=1;o<i;o++){if(B<j&&Ei(b.f,1,B)==o){ab=Ei(b.f,0,B++);$=1;}else{ab=0;$=8;}Zj(b.f,o,xi(b.f,ab)+$*(gm(b,bc)-(A/2|0)));$j(b.f,o,yi(b.f,ab)+$*(gm(b,bc)-(A/2|0)));O&&_j(b.f,o,zi(b.f,ab)+$*(gm(b,bc)-(A/2|0)));}if(e[M]==35){db=0;bb=ZB(_C,DQ,6,i,15,1);for(p=0;p<i;p++)db+=bb[p]=ml(b.f,p);for(m=0;m<i;m++){for(eb=0;eb<bb[m];eb++){cb=Ih(b.f,1);Jh(b.f,m,cb,1);Zj(b.f,cb,xi(b.f,m)+(gm(b,bc)-(A/2|0)));$j(b.f,cb,yi(b.f,m)+(gm(b,bc)-(A/2|0)));O&&_j(b.f,cb,zi(b.f,m)+(gm(b,bc)-(A/2|0)));}}i+=db;j+=db;}if(P){ec=fm(gm(b,bc),A);hc=ec*im(gm(b,bc),A);jc=ec*im(gm(b,bc),A);O&&(lc=ec*im(gm(b,bc),A));$=ec/Yk(b.f,true);for(m=0;m<i;m++){Zj(b.f,m,hc+$*xi(b.f,m));$j(b.f,m,jc+$*yi(b.f,m));O&&_j(b.f,m,lc+$*zi(b.f,m));}}else{$=1.5/Yk(b.f,true);for(m=0;m<i;m++){Zj(b.f,m,$*xi(b.f,m));$j(b.f,m,$*yi(b.f,m));O&&_j(b.f,m,$*zi(b.f,m));}}}else{O&&!P&&ec==0&&(ec=1.5);if(ec!=0&&b.f.p!=0){t=0;for(B=0;B<b.f.p;B++){V=xi(b.f,Ei(b.f,0,B))-xi(b.f,Ei(b.f,1,B));W=yi(b.f,Ei(b.f,0,B))-yi(b.f,Ei(b.f,1,B));X=O?zi(b.f,Ei(b.f,0,B))-zi(b.f,Ei(b.f,1,B)):0;t+=$wnd.Math.sqrt(V*V+W*W+X*X);}t/=b.f.p;Z=ec/t;for(m=0;m<b.f.o;m++){Zj(b.f,m,xi(b.f,m)*Z+hc);$j(b.f,m,yi(b.f,m)*Z+jc);O&&_j(b.f,m,zi(b.f,m)*Z+lc);}}}}catch(a){a=cG(a);if(PC(a,12)){Y=a;cK();'Faulty id-coordinates:'+Y.mb();e=null;O=false;}else throw dG(a);}}N=e!=null&&!O;if(N||b.b){Xo(b.f,3);for(B=0;B<b.f.e;B++)Mi(b.f,B)==2&&!Ol(b.f,B)&&Ni(b.f,B)==0&&hk(b.f,B);}if(!N&&b.b){Ob=new Tg();Ob.j=new SN(78187493520);Ig(Ob,b.f);N=true;}if(N){$l(b.f);hp(b.f);}else O||Xl(b.f,0);}function Aq(){Aq=FG;zq=aC(VB($C,1),gR,6,15,[-0.1899999976158142,1.2699999809265137,-0.7009999752044678,2.690999984741211,-0.22699999809265137,0.029999999329447746,0.10599999874830246,-0.47600001096725464,-0.44699999690055847,-0.19099999964237213,-0.3330000042915344,0.0860000029206276,0.24699999392032623,-0.06199999898672104,0.01600000075995922,0.3869999945163727,0.23499999940395355,-0.4320000112056732,-0.902999997138977,0.38999998569488525,0.5809999704360962,4.52400016784668,-0.6349999904632568,0.7919999957084656,0.5920000076293945,0.9639999866485596,HT,-0.6850000023841858,-0.3149999976158142,-0.4129999876022339,-0.5950000286102295,0.2199999988079071,-0.2800000011920929,0.7699999809265137,-0.05000000074505806,1.0870000123977661,0.19200000166893005,0.19599999487400055,-0.5199999809265137,0.5419999957084656,0.3630000054836273,PT,2.384000062942505,1.75,-1.6660000085830688,-1.065999984741211,1.3270000219345093,0.8029999732971191,-1.5049999952316284,-2.5369999408721924,QT,0.14900000393390656,0.5210000276565552,2.9049999713897705,-0.25200000405311584,-1.4320000410079956,-2.253999948501587,0.4399999976158142,-0.27000001072883606,-0.13300000131130219,-0.26899999380111694,0.2669999897480011,0.5720000267028809,-0.5680000185966492,0.17399999499320984,-0.1850000023841858,-0.23499999940395355,QT,PT,-0.34200000762939453,-0.3479999899864197,-0.43700000643730164,-0.8040000200271606,-0.41200000047683716,-0.2150000035762787,-0.625,-0.8309999704360962,0.4970000088214874,-0.4309999942779541,-1.3309999704360962,0.5070000290870667,-0.6320000290870667,-0.5989999771118164,0.8600000143051147,0.3610000014305115,0.40299999713897705,0.004999999888241291,1.1460000276565552,0.9359999895095825,-0.30000001192092896,0.20900000631809235,-0.5830000042915344,-0.024000000208616257,-0.009999999776482582,1.6469999551773071,0.843999981880188,0.125,0.1420000046491623,-0.17100000381469727,0.44200000166893005,0.08799999952316284,3.065999984741211,1.6519999504089355,-0.15600000321865082,-0.3529999852180481,-0.164000004529953,-0.4410000145435333,-0.4970000088214874,-1.059999942779541,0.6110000014305115,0.4860000014305115,0.11500000208616257,-0.22499999403953552,-0.15399999916553497,-0.03099999949336052,0.8619999885559082,-0.03500000014901161,-0.5960000157356262,-1.5740000009536743,-1.093000054359436,1.1610000133514404,KT,-0.44999998807907104,-0.5559999942779541,-0.621999979019165,2.121999979019165,-1.4019999504089355,2.072999954223633,-3.131999969482422,-2.119999885559082,0.34700000286102295,-1.2649999856948853,-1.3170000314712524,2.500999927520752,-2.2260000705718994,0.9129999876022339,-2.9570000171661377,0.29100000858306885,-0.7250000238418579,-1.4249999523162842,RT,-0.017999999225139618,-0.8489999771118164,-2.259000062942505,-3.4760000705718994,-0.296999990940094,-1.659999966621399,0.023000000044703484,0.0729999989271164,0.2540000081062317,0.5540000200271606,0.5950000286102295,JT,-1.25,1.3940000534057617,-2.7269999980926514,0.08299999684095383,-1.281999945640564,-0.4059999883174896,-0.6370000243186951,-0.17399999499320984,-0.10100000351667404,-0.5429999828338623,-2.4059998989105225,-3.2920000553131104,-0.6809999942779541,-1.2580000162124634,1.0700000524520874,-3.0959999561309814,-0.2280000001192093,0.718999981880188,0.1379999965429306,1.3020000457763672,0.859000027179718,1.3589999675750732,0.6589999794960022,-0.9399999976158142,0.8999999761581421,0.3190000057220459,-2.571000099182129,1.9329999685287476,0.11900000274181366,2.1080000400543213,0.11299999803304672,3.3359999656677246,0.7540000081062317,-0.4650000035762787,-0.05299999937415123,-0.19300000369548798,1.850000023841858,-1.2610000371932983,-0.656000018119812,-0.7300000190734863,-0.9380000233650208,1.1089999675750732,0.972000002861023,1.652999997138977,2.6019999980926514,1.628000020980835,-0.3970000147819519,0.12800000607967377,1.1540000438690186,0.24199999868869781,-0.5289999842643738,-0.27799999713897705,-0.8019999861717224,0.9120000004768372,-1.38100004196167,0.46299999952316284,1.0740000009536743,-0.628000020980835,-0.9620000123977661,0.7289999723434448,1.065999984741211,1.0670000314712524,-0.3109999895095825,0.03099999949336052,1.3079999685287476,0.07699999958276749,-0.4790000021457672,RT,-1.8320000171661377,-1.4989999532699585,-2.115999937057495,-2.2070000171661377,-0.15299999713897705,0.14100000262260437,2.134999990463257,0.23399999737739563,0.460999995470047,0.6700000166893005,-0.3610000014305115,-1.0390000343322754,-0.4830000102519989,0.13699999451637268,-0.7680000066757202,-0.5109999775886536,3.4240000247955322,-0.8550000190734863,-0.5849999785423279,-1.5670000314712524,0.6570000052452087,1.1150000095367432,1.9759999513626099,1.7860000133514404,-0.035999998450279236,-1.0499999523162842,2.5390000343322754,2.234999895095825,2.2899999618530273,3.121000051498413,3.931999921798706,2.75,3.3429999351501465,1.840000033378601,0.3889999985694885,1.121999979019165,1.6299999952316284,1.3350000381469727,0.3659999966621399,-0.5569999814033508,1.0449999570846558,0.4320000112056732,0.20399999618530273,0.8820000290870667,0.4659999907016754,-0.4580000042915344,0.04399999976158142,1.0329999923706055,-1.0800000429153442,0.40400001406669617]);yq=aC(VB(aD,1),gR,6,14,[262146,262148,262153,262157,264194,264195,264196,264197,264200,264201,264205,264206,267266,267267,267268,267273,267277,271362,271363,271364,271365,271368,271369,395266,395267,395268,395269,395272,395273,395277,395278,398338,526338,526339,526340,526344,529412,533508,533512,788482,788483,136448002,136448003,136448004,136448008,139593730,139593731,139593732,139593736,139596802,139596803,139596804,143788034,143788035,143791106,268697604,270794754,270794756,270796802,270796803,270796804,270796808,270796812,273940482,273942530,273942531,273942532,273942536,273945602,273945608,273945612,278136834,278136835,278136836,278136840,278136844,278139906,278139907,278139908,278144002,278144003,278144004,278144008,405014530,405014531,405014532,405014536,405017602,405017603,405017604,405021698,405021699,405021700,405021704,405145602,405145603,405145604,405145608,408158210,408160258,408163330,408167426,408291330,539232258,539232259,539235330,539235331,539239426,539239427,539363330,539363331,542377986,542377987,542381058,542381059,542509058,542509059,542509070,546837506,807667714,807798786,810813442,810816514,810820610,139722885122,139722885123,142944110594,142944110595,142947256322,142947259394,147239077890,147242223618,277161838594,277161838595,277164984322,277164984323,277164987394,277164987395,277169178626,277169181698,WR,XR,280383064066,280386209794,280386212866,280390404098,280390407170,YR,ZR,$R,_R,aS,284678031362,284681177090,284681177091,284681180162,284685371394,284685374466,bS,cS,dS,eS,fS,gS,284819727363,414600792066,414603937794,414603937795,414603940866,414603940867,hS,iS,jS,kS,lS,mS,nS,414742483970,414742488066,oS,pS,qS,rS,414869361667,sS,tS,uS,vS,wS,552177240067,xS,yS,552181437442,zS,AS,BS,CS,DS,ES,552311457794,FS,GS,HS,555398465539,IS,555398468611,JS,555398468616,KS,555402659848,LS,555402667010,MS,NS,OS,555529540615,PS,QS,RS,SS,555667032078,TS,US,VS,559697634306,WS,XS,{l:2361347,m:1376832,h:16},YS,ZS,$S,_S,aT,bT,cT,dT,eT,{l:1315842,m:2427201,h:16},fT,gT,hT,iT,jT,kT,lT,{l:2361346,m:592192,h:24},mT,nT,oT,{l:1315842,m:623169,h:24},{l:2361346,m:623200,h:24},pT,{l:2368514,m:623200,h:24},{l:2361351,m:1376832,h:32},qT,rT,sT,tT,uT,vT,{l:1319943,m:1378626,h:32},wT,xT,yT,zT,{l:1312776,m:1443138,h:32},{l:1315848,m:1443138,h:32},AT,BT,{l:2368520,m:1443168,h:32},CT]);}function $n(){$n=FG;Yn=aC(VB(kF,1),vR,2,6,['QM@HzAmdqjF@','RF@Q``','qC`@ISTAlQE`','`J@H','QM@HzAmdqbF@','qC`@ISTAlQEhqPp@','sJP@DiZhAmQEb','RF@QPvR@','QM@HzA@','qC`@ISTAlQEhpPp@','qC`@Qz`MbHl','sJP@DiZhAmQEcFZF@','RFPDXH','qC`@IVtAlQE`','QM@HvAmdqfF@','sGP@DiVj`FsDVM@','`L@H','sJP@DizhAmQEcFBF@','sJP@DjvhAmQEb','sFp@DiTt@@AlqEcP','sGP@LdbMU@MfHlZ','QMHAIhD','QM@HzAy@','sJP@DkVhAmQEb','sNp@DiUjj@[\\QXu`','sJP@DiZhAmQEcFBF@','sGP@DjVj`FsDVM@','RFPDTH','RG@DXOH@','sGP@Divj`FsDVMcAC@','sGP@Dj}j`FsDVM@','qC`@Qz`MbHmFRF@','sNp@LdbJjj@[\\QXu`','QMHAIhGe@','QM@HzAyd`','QM`AIhD','qC`@ISTA@','sGP@DkUj`FsDVM@','qC`@IVtAlQEhqPp@','sNp@DiUjj@[\\QXuqea`@','KAx@@IRjuUPAlHPfES\\','QM`BN`P','sJP@DjZhAmQEcFJF@','Hid@@DjU^nBBH@FtaBXUMp`','sNp@Diujj@[\\QXuq`a`@','sJP@DjvhAmQEcFZF@','sJP@DjZhAmQEcFFF@','sOp@DjWkB@@FwDVM\\YhX@','sNp@Dj}Zj@[\\QXu`','sNp@DiWjj@[\\QXuq`a`@','sOp@DjWkB@@D','KAx@@ITouUPAlHPfES\\','KAx@@YIDTjjh@vDHSBin@','sNp@DkUZj@[\\QXu`','RFPDXOH@','QM`BN`^L`','qC`@ISTAy@','sGP@LdbMU@MfHl[FVF@','qCb@AIZ`H','KAx@@IRjuUPAlHPfES]FFa`@','KAx@@ITnuUPAlHPfES\\','HiD@@DiUVjj`AmHPfES\\H','sNp@DjUjj@[\\QXu`','sJP@DkVhAmQEcFJF@','sGP@DjVj`FsDVMcCC@','qC`@Qz`MbHmFBF@','sJP@DkfhAmQEb','qC`@IVtAlQEhsPp@','sGP@Djuj`FsDVM@','sGP@Dj}j`FsDVMcMC@','sJP@DiZhA@','KAx@@ISjuUPAlHPfES]F@a`@','sJP@DjZhAmQEcFRF@','KAx@@IRnuUPAlHPfES]F@a`@','HiD@@DjWvjj`AmHPfES\\H','QMHAIhGd@','sNp@DiUjj@[\\QXuq`a`@','KAx@@IVjmUPAlHPfES\\','sGP@DjVj`FsDVMcMC@','QM`AIhGe@','HiD@@LdbJRjjh@[RDIaTwB','qCp@AIZ`H','sGP@LdbMU@MfHl[FFF@','QMDARVA@','sNp@LdbJjj@[\\QXuqba`@','sNp@LdbJjj@[\\QXuqca`@','sGP@Dkej`FsDVM@','qCb@AIZ`OI@','HaD@@DjUZxHH@AlHPfES]FLa`@','sGP@DkYj`FsDVM@','qCb@AIV`H','sNp@LdbJjj@[\\QXuqea`@','sGP@DkUj`FsDVMcEC@','sFp@DiTt@@Axa@','Hmt@@DjU_ZxHHj@AmhPfES\\Lj','QM`BN`^P','qCb@AIZ`OH`','sFp@DiTt@@AxaP','sGP@Djuj`FsDVMcEC@','sGP@Djuj`FsDVMcIC@','sGP@DkUj`FsDVMcKC@','sJP@DkfhAmQEcFRF@','sGP@DjVj`FsDVMcIC@','HaD@@DjUZxHH@AlHPfES]FFa`@','qC`@IRtDVqDV@','sNp@Dj}Zj@[\\QXuqfa`@','KAx@@ITnuUPAlHPfES]FFa`@','HiD@@DkUUjj`AmHPfES\\H','sJQ@@dkU@H','qC`@Qz`H','KAx@@IUkmUPAlHPfES\\','KAx@@ITouUPAlHPfES]FJa`@','sJP@H~j@[TQX`','sGP@DjZj`FsDVM@','sJP@DkVhAmQEcFFF@','sJX@@eKU@H','sJP@DizhAy@','QMHAIhGbP','KAx@@ITouUPAlHPfES]FNa`@','HaD@@DjUZxHD@AlHPfES\\','HaD@@DjUZxHH@A@','sNp@LdbJjj@[\\QXuqaa`@','Hed@@LdbRQUUUP@vTHSBinFP','KAx@@ITouUPAlHPfES]FLa`@','sNp@DkUZj@[\\QXuqba`@','KAx@@ITjuUPAlHPfES]FNa`@','KAx@@YIDTjjh@vDHSBincGPp@','HaD@@DjYvxH`@AlHPfES]FLa`@','RF@QP`','qCb@AIj`H','sNp@DjUjj@[\\QXuqaa`@','sNp@DkVZj@[\\QXu`','KAx@@YIDUJjh@vDHSBin@','sGP@DkYj`FsDVMcIC@','sGP@DjVj`FsDVMcAC@','sGP@DiVj`D','sJP@DkVhAmQEcFZF@','sNp@LdbLjj@[\\QXu`','QM@HvAmdqbF@','HaD@@DjWjXHB@AlHPfES\\','sNp@DjwZj@[\\QXuqba`@','sNp@LdbJjj@[\\QXuqda`@','sFp@DiTt@@Axa`','HiD@@Djuujj`AmHPfES\\H','sNp@DkUZj@[\\QXuqca`@','sJP@DiZhAy@','KAx@@YIDTjjh@vDHSBincCPp@','KAx@@IWNmUPAlHPfES\\','KAx@@IVkMUPAlHPfES\\','sJQ@@dju@H','qCb@AIZ`OH@','qC`@ISTAxa@','sNp@DjyZj@[\\QXu`','Hid@@DjUfaBB`@FtaBXUMp`','HiD@@DiUVjj`AmHPfES\\LXBF@','KAx@@IUjmUPAlHPfES\\','HiD@@DjWvjj`AmHPfES\\LXjF@','sJP@DjVhAmQEb','qCb@AIV`OH`','HiD@@LdbJRjjh@[RDIaTwCFDa`@','KAx@@YIDTjjh@vDHSBinc@Pp@','sNp@DjUjj@[\\QXuqda`@','qC`@Qz`OED','sJP@DkfhAmQEcFZF@','KAx@@YIDbjjh@vDHSBincDPp@','sGP@Djyj`FsDVMcMC@','KAx@@IVrmUPAlHPfES\\','qCp@AIZ`OI@','sJX@@dkU@H','sJQ@@dkU@OH`','sNp@Di]ZjBBvxbqk@','Hkl@@DjU_Uk``bj`@[VDIaTwCJzX','sGP@DjZj`FsDVMcEC@','Hid@@DjU^nBBH@FtaBXUMpqcHX@','sNp@DkeZj@[\\QXu`','sNp@DjYjj@[\\QXuqca`@','sGQ@@djuT@`','HiD@@LdbJTjjh@[RDIaTwB','sOp@DjWkB@@Gd`','HeT@@LdbbRKBDQD@CYPaLJfxY@','qCr@XIKTA@','HiD@@DjW^jj`AmHPfES\\LXJF@','HeT@@DjU]k``b`@[JDIaTwCH','sGP@Djuj`FsDVMcCC@','`IH`B','sOp@DjWkB@@GdX','sJQ@@eKU@H','KAx@@YIDUJjh@vDHSBincBPp@','sJX@@eKU@OH@','KAx@@YIDTjjh@vDHSBincAPp@','sOq@@drm\\@@@`','KAx@@IUkMUPAlHPfES\\','qCp@AIj`H','Hed@@DjUUjjj@FraBXUMpr','sGX@@eJuT@`','sGP@DkUj`FsDVMcCC@','HiD@@Dj}Ujj`AmHPfES\\LXrF@','KAx@@ITouUPAlHPfES]FHa`@','Hed@@DjWujjj@FraBXUMpsFIa`@','sGP@DiUj``mfHlZ','sFp@DiTvjhAlqEcP','Hid@@DjU^nBBH@FtaBXUMpq`XX@','sJP@DkVdAmQEb','qCp@AIZ`OH`','QMhDRVA@','qC`@ISJAlQE`','qCp@BOTAyhl','sJX@@eOU@ODB','sFp@DiTt@@AyaB','sGP@DkUj`FsDVMcMC@','Hid@@DjYUaBH`@FtaBXUMpqcHX@','qC`@Qz`OH@','HiD@@DjUVjj`AmHPfES\\LXZF@','sJP@H~j@[TQXqda`@','sJX@@eKU@OI@','sNp@Djejj@[\\QXu`','sJQ@@dsU@H','sJQ@@dkU@OI`','KAx@@YIMDVjh@vDHSBin@','Hid@@DjU^nBBD@FtaBXUMp`','sNp@DkgZj@[\\QXuqca`@','qC`@IRtDVqDVcEC@','Hed@@LdbRQeUUP@vTHSBinFP','sNp@DiUjj@P','qC`@IRtDT','sNp@DkYZj@[\\QXuqca`@','KAx@@IUkmUPAlHPfES]FDa`@','KAx@@IVjmUPAlHPfES]FNa`@','sOx@@drm\\@@@`','KAx@@ITjuUPAlHPfES]FBa`@','QMDARVAyH','sJP`@dfvhA@','HeT@@DjU_k``b`@[JDIaTwCLXfF@','KAx@@IToUUPAlHPfES]FJa`@','sGP@DkYj`FsDVMcEC@','qCb@AIZ`ODH','`I@`B','KAx@@IUzmUPAlHPfES]FFa`@','sNp@DkfZj@[\\QXu`','KAx@@ITnuUPAlHPfES]F@a`@','HiD@@LddURjjh@[RDIaTwB','sNp@Dj~Zj@[\\QXuqfa`@','Hed@@Dj{uZjj@FraBXUMpr','KAx@@ITsUUPAlHPfES\\','Hid@@LdbRQk``b@AmHPfES\\LXrF@','sOp@DjWkB@@GdH','sJQ@@dkU@OH@','Hid@@DjU^nBBH@FtaBXUMpqahX@','sGP@DiYj``mfHlZ','KAx@@IToUUPAlHPfES]FLa`@','qCp@AJZ`ODH','Hmt@@DjU]ZxHHj@AmhPfES\\Lj','sGP@DkUjPFsDVM@','qC`@IVtA@','Hed@@LdbJReUUP@vTHSBinFP','sNp@DjuZj@[\\QXuqea`@','KAx@@IUkmUPAlHPfES]FNa`@','HiD@@DkVUjj`AmHPfES\\H','Hed@@DkUeZjj@FraBXUMpr','sNp@DkVZj@[\\QXuqea`@','sJP@DiVhHKZbKFLLL@','HiD@@Djuyjj`AmHPfES\\H','sNp@DjUjj@[\\QXuq`a`@','HeT@@DjYUXPbH`@[JDIaTwCH','HiD@@DjwUjj`AmHPfES\\LXRF@','sNq@@djmUPB','KAx@@YIEEZjh@vDHSBincCPp@','sGP@Di^V`dmfHlZ','Hid@@DjYUaBHP@FtaBXUMp`','sNp@DjYjj@[\\QXuqba`@','sGP@Dkej`FsDVMcKC@','HeT@@DjU^k``b`@[JDIaTwCH','qC`@Qv`MbHmFBF@','sGQ@@djmT@`','qCr@XIKTAyH','qC`@IVtAlQEhpPp@','Hid@@LdbbQxXF@@AmHPfES\\LXjF@','sGP@DkYj`FsDVMcCC@','KAx@@IVsMUPAlHPfES\\','qCp@AIj`ODl','HiD@@DkeUjj`AmHPfES\\H','HeT@@DjU[kjjjh@ZLDXSSYPaLJfxY@','sJP@DkVdAmQEcFRF@','HiD@@LdbJTjjh@[RDIaTwCFDa`@','HiD@@DkYyjj`AmHPfES\\H','sJP@DjZhAyH','KAx@@IVkMUPAlHPfES]FDa`@','sJX@@dkU@OI@','Hed@@LdbRQUUUP@vTHSBinFXpLL@','Hed@@DjuUZjj@FraBXUMpr','sGP@Djfj`FsDVMcKC@','sNp@DkVZj@[\\QXuqba`@','sNp@DjyZj@[\\QXuqfa`@','qCb@AIj`OH@','sNp@DjUZj@[\\QXu`','KAx@@IWOMUPAlHPfES\\','Hid@@DjU^nBBH@D','Hed@@DjuvZjj@FraBXUMpr','sJP@DiVhHKZbKFLtL@','Hmt@@DjU_Zzjjj`AhpQaLmmBDpj[aeXplL@','sNp@DjuZj@[\\QXuqca`@','sJP@DkfhAmQEcFJF@','sNp@LdbJZj@[\\QXu`','HeT@@DjU_k``b`@[JDIaTwCLXFF@','KAx@@IVlmUPAlHPfES]FNa`@','HeT@@LdbbRKBDQD@CYPaLJfxYcEPp@','Hid@@DjUZnBBH@FtaBXUMpqcHX@','qCa@CIKTA@','HiD@@Dj~]jj`AmHPfES\\LXFF@','sKP@Di\\Zj@[TQX`','sGP@Djfj`FsDVMcEC@','HiD@@DkgYjj`AmHPfES\\H','sNp@DjuZj@[\\QXuqaa`@','KAx@@YIMDVjh@vDHSBincDPp@','sJP@DjVhHKZbKFLTL@','Hid@@LdbRQk``b@AmHPfES\\LXZF@','HiD@@Dj}Ujj`AmHPfES\\LXzF@','HeT@@DjU_k``bP@[JDIaTwCH','sNp@DkUZi@[\\QXu`','HiD@@DjYfjj`AmHPfES\\H','sGP@DjZj`FsDVMcAC@','Hmt@@DjU_jxHHj@AmhPfES\\Lj','Hid@@LdbRQk``R@AmHPfES\\H','KAx@@YIDUJjh@vDHSBincDPp@','qCr@XIKTAyD','sOq@@drm\\@@@|`@','Hed@@DjW^jjj@FraBXUMpsFBa`@','HeT@@DjY]zXFB@@[JDIaTwCH','Hkl@@DjU_Vk``bj`@[VDIaTwCJzX','Hid@@DjY}nBHH@FtaBXUMpqcHX@','sGX@@eKuT@|d@','sGP@Dj^Y`FsDVM@','HcL@@DjU_ZnBBJh@FqaBXUMprn`','sJP@DkVdAmQEcFJF@','sOq@@drm\\@@@|b@','sNp@DjyZj@[\\QXuqaa`@','HaD@@DjUZxHH@AyD@','qC`@Qv`H','Hmt@@DjU_Zzjjj`AhpQaLmmBDpj[aeXqdL@','sGP@Dkej`FsDVMcMC@','Hed@@DjUUjjj@FraBXUMpsFHa`@','HeT@@LdbbRkBDQD@CYPaLJfxY@','KAx@@IU{MUPAlHPfES]FLa`@','RG@DTH','sJY@DDeVhA@','KAx@@YIDUJjh@vDHSBinc@Pp@','sJX@@dkU@OI`','sJQ@@dju@OI`','HeT@@LdbbRKBDQD@CYPaLJfxYcFPp@','sFp@DiTvjhAlqEcXpPp@','HaD@@DjUZxHH@AyG@','sNx@@eJ}UPB','sNp@LddUjj@[\\QXuqca`@','HaDH@@RVU[j@@@D','sNp@DkgZi@[\\QXu`','sGY@LDeVj`D','sNp@LdbJfZBZvxbqk@','sJP`@dfvhAyL','sGX@AddQjhAxe`','Hmt@@DjU_ZxHHj@AmhPfES\\LkFIa`@','qCh@CIKTA@','sNp@LdbLjj@[\\QXuq`a`@','sOq@@drm\\@@@|a@','KAx@@IUzmUPAlHPfES]FJa`@','sNx@AddQUUPB','sGP@Di]jP`mfHlZ','sJP`@TeZhA@','KAx@@IRjmUPHKXPaLJfx','HeT@@LdbRTM\\DDT@CYPaLJfxY@','HaF@@@Rfu[j@@@D','Hid@@DjYUaBH`@FtaBXUMpqchX@','KAx@@IUjmTpAlHPfES\\','Hid@@DjU^nBBD@FtaBXUMpqcHX@','sGP@DiUj``mfHl[FFF@','KAx@@IUvmUPAlHPfES]FLa`@','Hed@@LdbQTUUUP@vTHSBinFXqDL@','sJP@DkVhA@','sOx@@drm\\@@@|b@','KAx@@IUkMUPAlHPfES]FDa`@','HeT@@LdbRQU\\DDT@CYPaLJfxY@','HiD@@Dj}Yjj`AmHPfES\\LXrF@','HiD@@Dj{ujj`AmHPfES\\LXFF@','KAx@@IWNmUPAlHPfES]FFa`@','KAx@@IRkMUPHKXPaLJfx','sJP@DjYdAmQEcFZF@','sJY@LDeZhAyL','HaDH@@RVU[f@@@D','sJP`@deVhAyB','HaD@@DjWjZjj`AlHPfES\\','sGP@DkYj`FsDVMcMC@','sNp@DkgZj@[\\QXuqea`@','sJQ@@dlu@H','HeT@@DjU]k``b`@[JDIaTwCLXrF@','sJX@@dkU@OH`','RFDDQFCr`','sJP@DiYXIKZbKFLLL@','KAx@@YIHjjjh@vDHSBincGPp@','Hk\\@@DjU^ukmLHH@@@AmXPfES\\Lki`','sGQ@@djmT@|b@','Hid@@DjUfaBB`@FtaBXUMpqahX@','sNx@@eRmUPB','Hmt@@LdbRVak``ah@FvaBXUMprh','qCr@XIJtA@','KAx@@IWMmUPAlHPfES]FNa`@','HeT@@DjYYZPbJ@@[JDIaTwCH','sNp@DkfZj@[\\QXuqea`@','Hid@@DjU^nBAHAEVtaBXUMp`','Hmt@@DjYU^Vjjj`AhtISRmmBDpj[aeP','sGP@DkejPFsDVM@','sNx@@eJmUPB','qCb@AIf`H','HcL@@DjU_VnBBJh@FqaBXUMprnqcXX@','Hid@@DjUZnBBH@FtaBXUMpqahX@','sNp@LdbQZjBBvxbqkcGC@','sOx@@drm\\@@@|c@','sJP@H~j@^R@','KAx@@YIDcFjhDElHPfES\\','Hid@@DjUZnBAH@FtaBXUMp`','sNp@LddUji@[\\QXu`','sGP@DjfjPFsDVM@','HeT@@DjYUXPbD`@[JDIaTwCH','KAx@@IUoMUPAlHPfES]FDa`@','sFp@DiTt@@AyaD','Hed@@DjuuZjj@FraBXUMpsFIa`@','HeT@@DjUghP`h`@[JDIaTwCLXfF@','sOp@DjWkjj`FwDVM\\YhX@','sGP@Djfj`FsDVMcIC@','KAx@@IRkmUPHKXPaLJfzL]C@','sNx@@djmUPB','QM`AIdD','sOp@DjWkB@@Gbe@','sNp@DjyZj@[\\QXuqca`@','QM@HuAmd`','sNp@LddUjj@[\\QXuqea`@','HaD@@DkeVyjj`AhrXUMuaBDpj[hpDL@','qCb@AIZPH','HiD@@LdbJTjjh@[RDIaTwCF@a`@','Hmt@@DjU_ZxHHi@AmhPfES\\Lj','HaDH@@RYWih@H@D','HiD@@LdbJTjjh@[RDIaTwCFHa`@','sGX@@djuT@|a@','sNp@DkfZj@[\\QXuqaa`@','Hid@@DjU^nBBH@GdL','KAx@@IVkMUPAlHPfES]FJa`@','qCr@XIKTAy@','HmT@@Dj{uVjjh@[ZDIaTwCJqaXX@','Hmt@@DjYWVFjjj`AhpQe\\mmBDpj[aeP','Hif@@@RUe^Fh@@@P','HaDH@@Rfu[j@@@GdH','KAx@@IVsMUPAlHPfES]FDa`@','sKP@Di\\Zj@[TQXq`a`@','sJX@@eMU@OH@','HeT@@DjU^k``b`@[JDIaTwCLXFF@','Hmt@@LdbbRJXPbHh@FvaBXUMprh','sJP@DjvhAmQEcFBF@','Hmt@@LdbbRNXZjjj@FcAFUrvtHSBinFUcBpp@','sJP`@dfvhAyD','sGP@Di^V`dmfHl[FVF@','KAx@@IVsmUPAlHPfES]FBa`@','sOq@@drm\\@@@|PP','sJY@BDeZhA@','HeT@@LdbRbmBDED@CYPaLJfxY@','Hed@@Djy[Zjj@FraBXUMpr','HeT@@DjU]k``b`@[JDIaTwCLXFF@','Hid@@DjUfaBB`@D','qCa@CIJtA@','QMPARVA@','Hid@@DjUfaBB`@FtaBXUMpqcHX@','sJY@BDfZhA@','HeT@@DjUghP`hP@[JDIaTwCH','Hed@@Dj{uZjj@FraBXUMpsFIa`@','Hmt@@LdbbRUXZjjj@FcAFUrvtHSBinFUcFPp@','sNp`@dfuZj@P','sJQ@@dmU@OH@','sJX@@dmU@H','HeT@@DjU]k``b`@[JDIaTwCLXZF@','HiD@@LdfbJZjh@[RDIaTwCFAa`@','sOx@@drm\\@@@|a@','HeT@@LdbbQgCUUU@CQhRfz[JDIaTwCH','Hmt@@DjU]Zzjjj`AhpQaLmmBDpj[aeXplL@','sOp@DjWkjj`FwDVM\\XHX@','HcL@@LdbbRNSBDQEP@McBDpj[ae]cFpp@','HiD@@Dj}Yji`AmHPfES\\H','HaDH@@RYe[hB@@D','Hid@@DjU^njjj@FtaBXUMpq`XX@','HeT@@DkYeFVjjh@ZMaUpsYPaLJfxY@','QMPARZA@','sOq@@drm\\@@@|QX','HaD@@DjYvxH`@A@','HcL@@LdbbRNcBDQEP@McBDpj[ae]@','QMhDRZA@','RG@DXLHmP','QM`BN`XQYd','RG@DTLHmP','QMHAIXFEVd','QMDARVAaH','RFPDXLHmP','RF@Q`vRbdLEC@','RF@QpvR@','QO@HyjAmd`','`II@B','`II@CFspqJp','`II@CF[@hM@prB`','`H@[T[|B`XN@PdM@p|@bHrBcDk@','RG@DXMj}F@','QM`BN`[L~b@','RG@DTMj}D@','QMHAIXFt~j@','QMDARVA}L@','RFPDXMj}D@','sKP@Di\\YZ@[TQXqaa`@','RG@DXMH']);}function Zp(){Zp=FG;Wp=aC(VB(aD,1),gR,6,14,[262146,262148,264194,264195,264196,264197,264200,264201,264205,264206,267266,267267,267268,267273,267278,271362,271363,271364,271365,271374,395266,395267,395268,395269,395272,395273,395277,395278,398338,398340,398345,524292,526338,526339,526340,526344,529412,533508,533512,788482,788483,136448066,136448067,136448068,136448072,139593794,139593795,139593796,139593800,139596866,139596867,139596868,139596872,143788098,143788099,143791170,143791171,143791172,270794754,270794760,270794824,270796802,270796803,270796804,270796808,270796812,270796866,270796867,270796868,270796872,273942530,273942531,273942532,273942536,273942594,273942595,273942596,273942600,273945602,273945608,273945666,278134788,278136834,278136835,278136836,278136840,278136898,278136899,278136900,278136904,278139906,278139907,278139908,278139970,278144002,278144008,278144010,278144066,405014530,405014531,405014532,405014536,405014540,405014594,405014595,405014596,405014600,405014604,405017602,405017603,405017604,405017666,405017672,405021698,405021699,405021700,405021704,405021762,405021763,405021764,405021768,405145602,405145603,405145604,405145608,405145666,405145667,405145668,405145672,408160258,408160260,408160322,408163330,408291330,408291331,539232258,539232259,539232322,539232323,539235330,539235331,539235394,539235395,539235396,539239426,539239427,539239490,539239491,539363330,539363331,539363394,539363395,539366402,542377986,542377987,542378050,542378051,542378052,542378056,542381058,542381059,542381122,542381123,542509058,542509059,542509123,542512130,546837506,807667714,807798786,810813442,810816514,810820610,810944514,139722885186,139722885187,142944110658,142944110659,142947256386,142947256392,142947259458,147239077954,147242223682,147246417986,277161838658,277161838659,277164984386,277164984387,277164987458,277164987459,277169178690,277169181762,277296185411,WR,XR,277296187458,277296187459,280383064130,280383064131,280386209858,280386209859,280386212930,280390404162,280390407234,YR,ZR,280517412930,280517412931,$R,_R,280520558658,aS,280520561730,284678031426,284681177154,284681180226,284685371458,284685374530,bS,cS,284812380226,284812380227,dS,284815525891,284815525954,eS,fS,gS,414600792130,414600792131,414603937858,414603937859,414603940930,414603940931,414608132162,414608135234,hS,iS,414735140930,414735140931,jS,kS,414738286658,414738286659,lS,414738289730,mS,nS,414742480962,414742480963,414742483971,414742484034,oS,pS,qS,414869358658,414869358659,rS,414869361730,sS,417822017602,417825163330,417825166402,tS,417956366402,uS,417959512130,417963706434,418090584066,418090715138,vS,552174094402,552174094403,wS,552177240130,552177240131,xS,552177243202,yS,552181434434,552181437506,zS,552181441602,AS,552308312130,BS,552308315202,CS,DS,552308319298,ES,552308443202,552311588930,FS,GS,555395319874,555395319875,555395319880,HS,555398465602,IS,JS,555398468674,555398468680,KS,555402659906,LS,555402662978,555402667074,MS,NS,555529537602,555529537603,555529537608,OS,555529540616,555529540674,PS,555529544770,QS,555529668616,555529668674,RS,SS,555532814338,559690287170,TS,559693432898,US,559693435970,559697627202,VS,559697630274,559824507906,559824507970,559827653634,{l:2359363,m:590400,h:16},{l:2361345,m:590400,h:16},WS,{l:2361410,m:590400,h:16},XS,{l:2361410,m:1376832,h:16},YS,{l:2361410,m:1377600,h:16},ZS,{l:1312834,m:1377601,h:16},$S,{l:2230339,m:2425376,h:16},_S,{l:2361347,m:2425408,h:16},{l:2361410,m:2425408,h:16},{l:2361411,m:2425408,h:16},aT,bT,cT,dT,{l:2361410,m:2427200,h:16},eT,fT,gT,hT,iT,{l:2361410,m:590400,h:24},jT,{l:2361410,m:591168,h:24},kT,{l:1312834,m:591169,h:24},lT,mT,nT,{l:2361410,m:623168,h:24},oT,{l:1312834,m:623169,h:24},pT,{l:2364482,m:623200,h:24},{l:2361346,m:1376832,h:24},{l:2361410,m:1376832,h:24},{l:2361347,m:592192,h:32},qT,rT,sT,{l:1315911,m:1377601,h:32},tT,uT,vT,{l:1315914,m:1378626,h:32},wT,xT,{l:2361416,m:1443136,h:32},yT,zT,{l:1315912,m:1443137,h:32},AT,{l:2361416,m:1443168,h:32},BT,{l:2364488,m:1443168,h:32},{l:2492424,m:1443168,h:32},{l:2492488,m:1443168,h:32},CT,{l:1315847,m:2426177,h:32},{l:1315847,m:2427201,h:32},{l:1315847,m:2458177,h:32},{l:264195,m:0,h:64},{l:264196,m:0,h:64},{l:264200,m:0,h:64},{l:267268,m:0,h:64},{l:271364,m:0,h:64},{l:395268,m:0,h:64},{l:398340,m:0,h:64},{l:529411,m:0,h:64},{l:2230339,m:32,h:64},{l:2361347,m:64,h:64},{l:1312771,m:65,h:64},{l:1574915,m:129,h:64},{l:1577987,m:129,h:64},{l:2364419,m:192,h:64},{l:2230339,m:66080,h:64},{l:1181763,m:66081,h:64},{l:2361347,m:66112,h:64},{l:2361411,m:66112,h:64},{l:2230339,m:66848,h:64},{l:1181763,m:66849,h:64},{l:1181763,m:98849,h:64},{l:2361347,m:131648,h:64},{l:1312771,m:131649,h:64},{l:1312835,m:131649,h:64},{l:2364419,m:131680,h:64},{l:1312771,m:132417,h:64},{l:1315843,m:132417,h:64},{l:2364419,m:132448,h:64},{l:2361347,m:590400,h:80},{l:2361411,m:590400,h:80},{l:2361347,m:1376832,h:80},{l:2361411,m:1376832,h:80},{l:2361347,m:590400,h:88},{l:264195,m:0,h:128},{l:2230339,m:32,h:128},{l:2361347,m:64,h:128},{l:2361411,m:64,h:128},{l:2368579,m:128,h:128},{l:2361347,m:66112,h:128},{l:2361411,m:66112,h:128}]);Xp=aC(VB($C,1),gR,6,15,[0.6966999769210815,0,0.4885999858379364,-0.47269999980926514,-0.07490000128746033,DT,0.273499995470047,0.5699999928474426,0.7009999752044678,0.9534000158309937,-0.2809000015258789,-0.8259999752044678,-0.1784999966621399,-1.620300054550171,-1.0959999561309814,0.13950000703334808,-0.29750001430511475,-1.2907999753952026,1.0161999464035034,ET,0.5110999941825867,-0.435699999332428,-0.10409999638795853,0.3424000144004822,-0.061500001698732376,0.6035000085830688,0.7226999998092651,0.43459999561309814,-0.3310000002384186,-0.49799999594688416,FT,GT,0.4291999936103821,-0.5824000239372253,-0.1834000051021576,0.1306000053882599,-0.5015000104904175,-0.5257999897003174,0.4244000017642975,-0.16099999845027924,-0.2777999937534332,0.2766000032424927,0.35929998755455017,0.7714999914169312,0.3149999976158142,-0.2651999890804291,-0.09650000184774399,0.420199990272522,0.18709999322891235,-0.3684000074863434,-0.07779999822378159,0.8942999839782715,0.3693999946117401,0.28790000081062317,0.4489000141620636,-0.26010000705718994,0.4771000146865845,0.1923000067472458,0.45969998836517334,0.3384000062942505,0.6632999777793884,0.4544000029563904,0.15970000624656677,0.633899986743927,0.35040000081062317,0.04490000009536743,0.34200000762939453,0.26109999418258667,0.40459999442100525,0.5218999981880188,-0.36320000886917114,-0.4108000099658966,0.30570000410079956,-0.14560000598430634,-0.27129998803138733,-0.5192999839782715,0.45260000228881836,0.5539000034332275,-0.7070000171661377,-0.48809999227523804,-0.4099999964237213,0,0.14790000021457672,0.3447999954223633,0.42980000376701355,0.5579000115394592,-0.1264999955892563,-0.042500000447034836,0.07670000195503235,0.6635000109672546,-0.38119998574256897,-0.8367999792098999,1.0286999940872192,-0.10209999978542328,0.3587000072002411,-0.5945000052452087,0.16920000314712524,-0.121799997985363,0.43810001015663147,0.16949999332427979,0.45249998569488525,0.3352000117301941,0.1582999974489212,0.4036000072956085,-0.04800000041723251,0.5023000240325928,-0.26489999890327454,0.76910001039505,-0.35519999265670776,1.0300999879837036,-0.11410000175237656,-0.5932000279426575,0.17489999532699585,0.13130000233650208,-0.18039999902248383,0.399399995803833,0.22910000383853912,0.31690001487731934,0.35989999771118164,-0.0038999998942017555,-0.2955999970436096,0.49070000648498535,HT,0.219200000166893,0.15649999678134918,0.6934999823570251,0.3617999851703644,0.6735000014305115,0.5777999758720398,-0.5636000037193298,0.5569000244140625,0.30379998683929443,-0.32760000228881836,-0.4659000039100647,IT,0.32829999923706055,0.22390000522136688,0.20430000126361847,0.05900000035762787,-0.48350000381469727,0.6165000200271606,-0.4011000096797943,0.5577999949455261,-0.21639999747276306,-0.017500000074505806,0.29809999465942383,0.10999999940395355,0.27149999141693115,0.4408999979496002,-0.16089999675750732,0.3774999976158142,-0.13459999859333038,-0.6991999745368958,-0.46700000762939453,0.1565999984741211,0.046799998730421066,-0.13210000097751617,1.3686000108718872,0,-0.4115999937057495,1.0185999870300293,-0.3935000002384186,0.5223000049591064,0.2838999927043915,0.5128999948501587,0.1265999972820282,0.010300000198185444,1.5192999839782715,0.2705000042915344,0.4293999969959259,0.012000000104308128,-0.33970001339912415,0.14830000698566437,0.28060001134872437,0.3206000030040741,0.5662000179290771,-0.09870000183582306,-0.10050000250339508,-0.35760000348091125,0.09610000252723694,-0.6401000022888184,0.19210000336170197,-0.15330000221729279,-0.4169999957084656,0.10939999669790268,0.8230999708175659,-0.3783999979496002,0.4032000005245209,-0.6460999846458435,0.8034999966621399,0.2029000073671341,-0.37450000643730164,GT,0.18410000205039978,0.707099974155426,0.12269999831914902,0.7949000000953674,0.03500000014901161,IT,-0.15539999306201935,0.3785000145435333,-0.24050000309944153,0.23589999973773956,0.34630000591278076,-0.4925000071525574,-0.09290000051259995,-0.4352000057697296,-0.2206999957561493,-0.9959999918937683,-0.723800003528595,-0.5468999743461609,-1.2939000129699707,-0.01360000018030405,0.2791000008583069,-0.16529999673366547,-0.12380000203847885,0.4950999915599823,0.289900004863739,0.065700002014637,0.7189000248908997,0.05700000002980232,0.661899983882904,-0.6381000280380249,-0.8072999715805054,0.23549999296665192,0.30480000376701355,-0.019899999722838402,-0.07519999891519547,0.27639999985694885,0.8011000156402588,-0.17440000176429749,0.15809999406337738,-0.384799987077713,0.5993000268936157,0.5267999768257141,-0.04170000180602074,0.37700000405311584,0.6998000144958496,0.593999981880188,0.5911999940872192,-0.5570999979972839,0.023800000548362732,-0.2475000023841858,0.030700000002980232,-0.38749998807907104,-0.7437000274658203,0.5144000053405762,0.00570000009611249,0.765500009059906,0.1720000058412552,-2.5624001026153564,-0.30660000443458557,0.36469998955726624,0.4733000099658966,-0.3400999903678894,-0.14499999582767487,0.7088000178337097,-0.13179999589920044,0.04259999841451645,-0.12030000239610672,-0.36239999532699585,0.5357999801635742,-0.3700999915599823,-0.5648000240325928,-0.1972000002861023,-0.8769000172615051,-0.3675000071525574,-0.2003999948501587,0.13359999656677246,-0.16990000009536743,0.44609999656677246,0.1559000015258789,1.1167999505996704,0.23649999499320984,-0.22059999406337738,0.4480000138282776,-0.40529999136924744,-0.13609999418258667,0.2198999971151352,0.053599998354911804,-0.020999999716877937,0.6984999775886536,0.9642999768257141,0.17270000278949738,-0.03290000185370445,-0.18930000066757202,0.07020000368356705,0.14959999918937683,ET,0.4146000146865845,-0.5027999877929688,0.3831999897956848,0.9545000195503235,-0.41519999504089355,-1.0369000434875488,-0.18299999833106995,0.5882999897003174,-0.29179999232292175,-0.5293999910354614,-0.6541000008583069,-0.19059999287128448,-0.8483999967575073,-0.3456999957561493,0.9541000127792358,-0.7924000024795532,JT,0.07999999821186066,-0.2596000134944916,0.8381999731063843,-0.2667999863624573,-0.11060000211000443,0.03620000183582306,-0.3188999891281128,-0.7278000116348267,-0.08940000087022781,-0.22769999504089355,-0.2393999993801117,-0.2962000072002411,0.7775999903678894,-0.011800000444054604,-0.4357999861240387,0.3749000132083893,-0.6069999933242798,-0.18569999933242798,0.11389999836683273,-0.4415999948978424,-0.37040001153945923,-0.7487000226974487,-0.10790000110864639,-0.29919999837875366,-0.3276999890804291,0.025100000202655792,-0.9187999963760376,0.2944999933242798,-0.22339999675750732,0.3467999994754791,GT,0.2890999913215637,0.2612999975681305,-0.03440000116825104,-0.6004999876022339,-0.6258000135421753,-0.54339998960495,-0.7712000012397766,-0.9057000279426575,-0.16680000722408295,-0.9904999732971191,FT,-0.03720000013709068,-1.1638000011444092,0.12620000541210175,-0.5248000025749207,-0.15379999577999115,-0.36820000410079956,0.3249000012874603,0.06499999761581421,0.051100000739097595,-0.46070000529289246,0.22310000658035278,0.28220000863075256,0.1396999955177307,0.2833000123500824,-0.1225999966263771,-0.459199994802475,-0.3434999883174896,-0.6654000282287598,-0.5055999755859375,-0.863099992275238,0.15360000729560852,-0.4050999879837036,0.08910000324249268,-0.6972000002861023,-0.4699000120162964,-0.6773999929428101,-0.062199998646974564,-0.9300000071525574,0.13369999825954437,-0.49380001425743103,0.39480000734329224,-0.4074999988079071,-0.6410999894142151,-0.009100000374019146,-0.13330000638961792,-0.5192000269889832,-0.16609999537467957,GT,-0.6427000164985657,-0.07069999724626541,0.4805999994277954,0.38280001282691956,0.22290000319480896,0.6159999966621399,-0.08839999884366989,-0.0471000000834465,0.11060000211000443,0.382099986076355,0.09220000356435776,0.08060000091791153,0.33709999918937683,0.188400000333786,0.13809999823570251,-0.23919999599456787,-3.6435000896453857,-2.150899887084961,0.4975000023841858,-0.3619999885559082,-2.5383999347686768,-1.6821000576019287,-0.3650999963283539,DT,-1.618499994277954,-1.065000057220459,0.8374000191688538,0.3684999942779541,0.25769999623298645,0.3790999948978424,-3.233299970626831,-1.794800043106079,-0.6592000126838684,-1.3148000240325928,KT,0.05339999869465828,-1.7552000284194946,-1.8039000034332275,-1.1339999437332153,-0.5652999877929688,-1.2453999519348145,0.9473000168800354,1.4230999946594238,1.011199951171875,-1.9498000144958496,-2.024899959564209,-1.2348999977111816,0.328000009059906,-3.9189000129699707,-2.19950008392334,0.18889999389648438,-1.2314000129699707,-1.802299976348877,-0.2994999885559082,-0.4066999852657318,-0.1316000074148178]);}function vm(){vm=FG;um=aC(VB(AD,2),iQ,9,0,[null,aC(VB(AD,1),KR,3,0,[new tm(0,1.007825032),new tm(1,2.014101778),new tm(2,3.016049268),new tm(3,4.027834627),new tm(4,5.039542911),new tm(5,6.044942608)]),aC(VB(AD,1),KR,3,0,[new tm(1,3.01602931),new tm(2,4.00260325),new tm(3,5.012223628),new tm(4,6.018888072),new tm(5,7.028030527),new tm(6,8.033921838),new tm(7,9.043820323),new tm(8,10.052399713)]),aC(VB(AD,1),KR,3,0,[new tm(1,4.027182329),new tm(2,5.012537796),new tm(3,6.015122281),new tm(4,7.016004049),new tm(5,8.02248667),new tm(6,9.026789122),new tm(7,10.035480884),new tm(8,11.043796166),new tm(9,12.05378)]),aC(VB(AD,1),KR,3,0,[new tm(1,5.04079),new tm(2,6.019725804),new tm(3,7.016929246),new tm(4,8.005305094),new tm(5,9.012182135),new tm(6,10.01353372),new tm(7,11.021657653),new tm(8,12.026920631),new tm(9,13.036133834),new tm(10,14.042815522)]),aC(VB(AD,1),KR,3,0,[new tm(2,7.029917389),new tm(3,8.024606713),new tm(4,9.013328806),new tm(5,10.012937027),new tm(6,11.009305466),new tm(7,12.014352109),new tm(8,13.017780267),new tm(9,14.025404064),new tm(10,15.031097291),new tm(11,16.039808836),new tm(12,17.046931399),new tm(13,18.05617),new tm(14,19.06373)]),aC(VB(AD,1),KR,3,0,[new tm(2,8.037675026),new tm(3,9.031040087),new tm(4,10.01685311),new tm(5,11.011433818),new tm(6,12),new tm(7,13.003354838),new tm(8,14.003241988),new tm(9,15.010599258),new tm(10,16.014701243),new tm(11,17.022583712),new tm(12,18.026757058),new tm(13,19.035248094),new tm(14,20.040322395),new tm(15,21.04934),new tm(16,22.05645)]),aC(VB(AD,1),KR,3,0,[new tm(3,10.042618),new tm(4,11.026796226),new tm(5,12.018613202),new tm(6,13.005738584),new tm(7,14.003074005),new tm(8,15.000108898),new tm(9,16.006101417),new tm(10,17.008449673),new tm(11,18.014081827),new tm(12,19.017026896),new tm(13,20.023367295),new tm(14,21.027087574),new tm(15,22.034440259),new tm(16,23.04051),new tm(17,24.0505)]),aC(VB(AD,1),KR,3,0,[new tm(4,12.034404776),new tm(5,13.0248104),new tm(6,14.008595285),new tm(7,15.003065386),new tm(8,15.994914622),new tm(9,16.999131501),new tm(10,17.999160419),new tm(11,19.00357873),new tm(12,20.00407615),new tm(13,21.008654631),new tm(14,22.009967157),new tm(15,23.015691325),new tm(16,24.020369922),new tm(17,25.02914),new tm(18,26.03775)]),aC(VB(AD,1),KR,3,0,[new tm(5,14.03608),new tm(6,15.018010856),new tm(7,16.01146573),new tm(8,17.002095238),new tm(9,18.000937667),new tm(10,18.998403205),new tm(11,19.999981324),new tm(12,20.999948921),new tm(13,22.00299925),new tm(14,23.003574385),new tm(15,24.008099371),new tm(16,25.012094963),new tm(17,26.019633157),new tm(18,27.026892316),new tm(19,28.03567),new tm(20,29.04326)]),aC(VB(AD,1),KR,3,0,[new tm(6,16.025756907),new tm(7,17.017697565),new tm(8,18.005697066),new tm(9,19.001879839),new tm(10,19.992440176),new tm(11,20.993846744),new tm(12,21.99138551),new tm(13,22.994467337),new tm(14,23.993615074),new tm(15,24.997789899),new tm(16,26.000461498),new tm(17,27.0076152),new tm(18,28.012108072),new tm(19,29.019345902),new tm(20,30.023872),new tm(21,31.03311),new tm(22,32.03991)]),aC(VB(AD,1),KR,3,0,[new tm(7,18.02718),new tm(8,19.01387945),new tm(9,20.00734826),new tm(10,20.997655099),new tm(11,21.994436782),new tm(12,22.989769675),new tm(13,23.990963332),new tm(14,24.989954352),new tm(15,25.992589898),new tm(16,26.994008702),new tm(17,27.99889041),new tm(18,29.002811301),new tm(19,30.009226487),new tm(20,31.013595108),new tm(21,32.019649792),new tm(22,33.027386),new tm(23,34.0349),new tm(24,35.04418)]),aC(VB(AD,1),KR,3,0,[new tm(8,20.018862744),new tm(9,21.011714174),new tm(10,21.999574055),new tm(11,22.99412485),new tm(12,23.985041898),new tm(13,24.985837023),new tm(14,25.98259304),new tm(15,26.984340742),new tm(16,27.983876703),new tm(17,28.988554743),new tm(18,29.990464529),new tm(19,30.996548459),new tm(20,31.999145889),new tm(21,33.005586975),new tm(22,34.00907244),new tm(23,35.018669),new tm(24,36.02245),new tm(25,37.03124)]),aC(VB(AD,1),KR,3,0,[new tm(8,21.02804),new tm(9,22.01952),new tm(10,23.0072649),new tm(11,23.999940911),new tm(12,24.990428555),new tm(13,25.986891659),new tm(14,26.981538441),new tm(15,27.981910184),new tm(16,28.980444848),new tm(17,29.982960304),new tm(18,30.983946023),new tm(19,31.988124379),new tm(20,32.990869587),new tm(21,33.996927255),new tm(22,34.99993765),new tm(23,36.006351501),new tm(24,37.01031),new tm(25,38.0169),new tm(26,39.0219)]),aC(VB(AD,1),KR,3,0,[new tm(8,22.03453),new tm(9,23.02552),new tm(10,24.011545711),new tm(11,25.00410664),new tm(12,25.992329935),new tm(13,26.986704764),new tm(14,27.976926533),new tm(15,28.976494719),new tm(16,29.973770218),new tm(17,30.975363275),new tm(18,31.974148129),new tm(19,32.97800052),new tm(20,33.978575745),new tm(21,34.984584158),new tm(22,35.986687363),new tm(23,36.99299599),new tm(24,37.99598),new tm(25,39.0023),new tm(26,40.0058),new tm(27,41.0127),new tm(28,42.0161)]),aC(VB(AD,1),KR,3,0,[new tm(9,24.03435),new tm(10,25.02026),new tm(11,26.01178),new tm(12,26.999191645),new tm(13,27.99231233),new tm(14,28.981801376),new tm(15,29.978313807),new tm(16,30.973761512),new tm(17,31.973907163),new tm(18,32.971725281),new tm(19,33.973636381),new tm(20,34.973314249),new tm(21,35.978259824),new tm(22,36.979608338),new tm(23,37.98447),new tm(24,38.98642),new tm(25,39.99105),new tm(26,40.9948),new tm(27,42.00009),new tm(28,43.00331),new tm(29,44.00988),new tm(30,45.01514),new tm(31,46.02383)]),aC(VB(AD,1),KR,3,0,[new tm(10,26.02788),new tm(11,27.018795),new tm(12,28.004372661),new tm(13,28.996608805),new tm(14,29.984902954),new tm(15,30.979554421),new tm(16,31.97207069),new tm(17,32.971458497),new tm(18,33.967866831),new tm(19,34.96903214),new tm(20,35.96708088),new tm(21,36.971125716),new tm(22,37.971163443),new tm(23,38.975135275),new tm(24,39.97547),new tm(25,40.98003),new tm(26,41.98149),new tm(27,42.9866),new tm(28,43.98832),new tm(29,44.99482),new tm(30,45.99957),new tm(31,47.00762),new tm(32,48.01299),new tm(33,49.02201)]),aC(VB(AD,1),KR,3,0,[new tm(11,28.02851),new tm(12,29.01411),new tm(13,30.00477),new tm(14,30.992416014),new tm(15,31.985688908),new tm(16,32.977451798),new tm(17,33.973761967),new tm(18,34.968852707),new tm(19,35.968306945),new tm(20,36.9659026),new tm(21,37.96801055),new tm(22,38.968007677),new tm(23,39.970415555),new tm(24,40.970650212),new tm(25,41.973174994),new tm(26,42.974203385),new tm(27,43.978538712),new tm(28,44.9797),new tm(29,45.98412),new tm(30,46.98795),new tm(31,47.99485),new tm(32,48.99989),new tm(33,50.00773),new tm(34,51.01353)]),aC(VB(AD,1),KR,3,0,[new tm(12,30.02156),new tm(13,31.012126),new tm(14,31.99766066),new tm(15,32.989928719),new tm(16,33.980270118),new tm(17,34.975256726),new tm(18,35.967546282),new tm(19,36.966775912),new tm(20,37.962732161),new tm(21,38.964313413),new tm(22,39.962383123),new tm(23,40.964500828),new tm(24,41.963046386),new tm(25,42.965670701),new tm(26,43.965365269),new tm(27,44.968094979),new tm(28,45.968093467),new tm(29,46.972186238),new tm(30,47.97507),new tm(31,48.98218),new tm(32,49.98594),new tm(33,50.99324),new tm(34,51.99817),new tm(35,53.006227)]),aC(VB(AD,1),KR,3,0,[new tm(13,32.02192),new tm(14,33.00726),new tm(15,33.99841),new tm(16,34.988011615),new tm(17,35.981293405),new tm(18,36.973376915),new tm(19,37.969080107),new tm(20,38.963706861),new tm(21,39.963998672),new tm(22,40.961825972),new tm(23,41.962403059),new tm(24,42.960715746),new tm(25,43.961556146),new tm(26,44.960699658),new tm(27,45.961976203),new tm(28,46.961677807),new tm(29,47.965512946),new tm(30,48.967450084),new tm(31,49.972782832),new tm(32,50.97638),new tm(33,51.98261),new tm(34,52.98712),new tm(35,53.99399),new tm(36,54.999388)]),aC(VB(AD,1),KR,3,0,[new tm(14,34.01412),new tm(15,35.004765),new tm(16,35.993087234),new tm(17,36.985871505),new tm(18,37.976318637),new tm(19,38.970717729),new tm(20,39.962591155),new tm(21,40.962278349),new tm(22,41.958618337),new tm(23,42.958766833),new tm(24,43.955481094),new tm(25,44.956185938),new tm(26,45.953692759),new tm(27,46.954546459),new tm(28,47.952533512),new tm(29,48.955673302),new tm(30,49.957518286),new tm(31,50.961474238),new tm(32,51.9651),new tm(33,52.97005),new tm(34,53.97468),new tm(35,54.98055),new tm(36,55.98579),new tm(37,56.992356)]),aC(VB(AD,1),KR,3,0,[new tm(15,36.01492),new tm(16,37.00305),new tm(17,37.9947),new tm(18,38.984790009),new tm(19,39.977964014),new tm(20,40.969251316),new tm(21,41.965516761),new tm(22,42.96115098),new tm(23,43.959403048),new tm(24,44.955910243),new tm(25,45.95517025),new tm(26,46.952408027),new tm(27,47.952234991),new tm(28,48.950024065),new tm(29,49.952187008),new tm(30,50.9536027),new tm(31,51.95665),new tm(32,52.95817),new tm(33,53.963),new tm(34,54.9694),new tm(35,55.97266),new tm(36,56.97704),new tm(37,57.98307),new tm(38,58.988041)]),aC(VB(AD,1),KR,3,0,[new tm(16,38.00977),new tm(17,39.001323),new tm(18,39.990498907),new tm(19,40.983131),new tm(20,41.973031622),new tm(21,42.968523342),new tm(22,43.959690235),new tm(23,44.958124349),new tm(24,45.952629491),new tm(25,46.951763792),new tm(26,47.947947053),new tm(27,48.947870789),new tm(28,49.944792069),new tm(29,50.946616017),new tm(30,51.946898175),new tm(31,52.949731709),new tm(32,53.95087),new tm(33,54.95512),new tm(34,55.95799),new tm(35,56.9643),new tm(36,57.96611),new tm(37,58.97196),new tm(38,59.97564),new tm(39,60.982018)]),aC(VB(AD,1),KR,3,0,[new tm(17,40.01109),new tm(18,40.99974),new tm(19,41.99123),new tm(20,42.98065),new tm(21,43.9744),new tm(22,44.965782286),new tm(23,45.960199491),new tm(24,46.954906918),new tm(25,47.95225448),new tm(26,48.948516914),new tm(27,49.947162792),new tm(28,50.943963675),new tm(29,51.944779658),new tm(30,52.944342517),new tm(31,53.946444381),new tm(32,54.947238194),new tm(33,55.95036),new tm(34,56.95236),new tm(35,57.95665),new tm(36,58.9593),new tm(37,59.9645),new tm(38,60.96741),new tm(39,61.97314),new tm(40,62.97675)]),aC(VB(AD,1),KR,3,0,[new tm(18,42.00643),new tm(19,42.997707),new tm(20,43.98547),new tm(21,44.97916),new tm(22,45.968361649),new tm(23,46.962906512),new tm(24,47.954035861),new tm(25,48.951341135),new tm(26,49.946049607),new tm(27,50.944771767),new tm(28,51.940511904),new tm(29,52.940653781),new tm(30,53.938884921),new tm(31,54.940844164),new tm(32,55.940645238),new tm(33,56.9437538),new tm(34,57.94425),new tm(35,58.94863),new tm(36,59.94973),new tm(37,60.95409),new tm(38,61.9558),new tm(39,62.96186),new tm(40,63.9642),new tm(41,64.97037)]),aC(VB(AD,1),KR,3,0,[new tm(19,44.00687),new tm(20,44.99451),new tm(21,45.98672),new tm(22,46.9761),new tm(23,47.96887),new tm(24,48.959623415),new tm(25,49.95424396),new tm(26,50.948215487),new tm(27,51.945570079),new tm(28,52.941294702),new tm(29,53.940363247),new tm(30,54.938049636),new tm(31,55.938909366),new tm(32,56.938287458),new tm(33,57.939986451),new tm(34,58.940447166),new tm(35,59.943193998),new tm(36,60.94446),new tm(37,61.94797),new tm(38,62.94981),new tm(39,63.95373),new tm(40,64.9561),new tm(41,65.96082),new tm(42,66.96382)]),aC(VB(AD,1),KR,3,0,[new tm(19,45.01456),new tm(20,46.00081),new tm(21,46.99289),new tm(22,47.98056),new tm(23,48.97361),new tm(24,49.962993316),new tm(25,50.956824936),new tm(26,51.948116526),new tm(27,52.945312282),new tm(28,53.939614836),new tm(29,54.938298029),new tm(30,55.934942133),new tm(31,56.935398707),new tm(32,57.933280458),new tm(33,58.934880493),new tm(34,59.934076943),new tm(35,60.936749461),new tm(36,61.936770495),new tm(37,62.940118442),new tm(38,63.94087),new tm(39,64.94494),new tm(40,65.94598),new tm(41,66.95),new tm(42,67.95251),new tm(43,68.9577)]),aC(VB(AD,1),KR,3,0,[new tm(21,48.00176),new tm(22,48.98972),new tm(23,49.98154),new tm(24,50.97072),new tm(25,51.96359),new tm(26,52.954224985),new tm(27,53.948464147),new tm(28,54.942003149),new tm(29,55.939843937),new tm(30,56.936296235),new tm(31,57.935757571),new tm(32,58.933200194),new tm(33,59.933822196),new tm(34,60.932479381),new tm(35,61.934054212),new tm(36,62.933615218),new tm(37,63.935813523),new tm(38,64.936484581),new tm(39,65.939825412),new tm(40,66.94061),new tm(41,67.94436),new tm(42,68.9452),new tm(43,69.94981),new tm(44,70.95173),new tm(45,71.95641)]),aC(VB(AD,1),KR,3,0,[new tm(22,49.99593),new tm(23,50.98772),new tm(24,51.97568),new tm(25,52.96846),new tm(26,53.957910508),new tm(27,54.951336329),new tm(28,55.942136339),new tm(29,56.939800489),new tm(30,57.935347922),new tm(31,58.934351553),new tm(32,59.930790633),new tm(33,60.931060442),new tm(34,61.928348763),new tm(35,62.929672948),new tm(36,63.927969574),new tm(37,64.930088013),new tm(38,65.929115232),new tm(39,66.931569638),new tm(40,67.931844932),new tm(41,68.935181837),new tm(42,69.93614),new tm(43,70.94),new tm(44,71.9413),new tm(45,72.94608),new tm(46,73.94791),new tm(47,74.95297),new tm(48,75.95533),new tm(49,76.96083),new tm(50,77.9638)]),aC(VB(AD,1),KR,3,0,[new tm(23,51.99718),new tm(24,52.98555),new tm(25,53.97671),new tm(26,54.96605),new tm(27,55.95856),new tm(28,56.949215695),new tm(29,57.944540734),new tm(30,58.939504114),new tm(31,59.937368123),new tm(32,60.933462181),new tm(33,61.932587299),new tm(34,62.929601079),new tm(35,63.929767865),new tm(36,64.927793707),new tm(37,65.928873041),new tm(38,66.927750294),new tm(39,67.929637875),new tm(40,68.929425281),new tm(41,69.932409287),new tm(42,70.932619818),new tm(43,71.93552),new tm(44,72.93649),new tm(45,73.9402),new tm(46,74.9417),new tm(47,75.94599),new tm(48,76.94795),new tm(49,77.95281),new tm(50,78.95528),new tm(51,79.96189)]),aC(VB(AD,1),KR,3,0,[new tm(24,53.99295),new tm(25,54.98398),new tm(26,55.97238),new tm(27,56.96491),new tm(28,57.954596465),new tm(29,58.949267074),new tm(30,59.941832031),new tm(31,60.939513907),new tm(32,61.934334132),new tm(33,62.933215563),new tm(34,63.929146578),new tm(35,64.929245079),new tm(36,65.926036763),new tm(37,66.927130859),new tm(38,67.924847566),new tm(39,68.926553538),new tm(40,69.92532487),new tm(41,70.927727195),new tm(42,71.926861122),new tm(43,72.929779469),new tm(44,73.929458261),new tm(45,74.932937379),new tm(46,75.933394207),new tm(47,76.937085857),new tm(48,77.938569576),new tm(49,78.942095175),new tm(50,79.944414722),new tm(51,80.95048),new tm(52,81.95484)]),aC(VB(AD,1),KR,3,0,[new tm(25,55.99491),new tm(26,56.98293),new tm(27,57.97425),new tm(28,58.96337),new tm(29,59.95706),new tm(30,60.94917),new tm(31,61.944179608),new tm(32,62.939141527),new tm(33,63.936838307),new tm(34,64.932739322),new tm(35,65.931592355),new tm(36,66.928204915),new tm(37,67.927983497),new tm(38,68.925580912),new tm(39,69.926027741),new tm(40,70.92470501),new tm(41,71.92636935),new tm(42,72.925169832),new tm(43,73.926940999),new tm(44,74.926500645),new tm(45,75.928928262),new tm(46,76.929281189),new tm(47,77.93165595),new tm(48,78.932916371),new tm(49,79.936588154),new tm(50,80.937752955),new tm(51,81.94316),new tm(52,82.94687),new tm(53,83.95234)]),aC(VB(AD,1),KR,3,0,[new tm(26,57.99101),new tm(27,58.98175),new tm(28,59.97019),new tm(29,60.96379),new tm(30,61.95465),new tm(31,62.94964),new tm(32,63.941572638),new tm(33,64.939440762),new tm(34,65.933846798),new tm(35,66.932738415),new tm(36,67.928097266),new tm(37,68.927972002),new tm(38,69.924250365),new tm(39,70.924953991),new tm(40,71.922076184),new tm(41,72.923459361),new tm(42,73.921178213),new tm(43,74.922859494),new tm(44,75.921402716),new tm(45,76.923548462),new tm(46,77.922852886),new tm(47,78.92540156),new tm(48,79.925444764),new tm(49,80.928821065),new tm(50,81.929550326),new tm(51,82.93451),new tm(52,83.93731),new tm(53,84.94269),new tm(54,85.94627)]),aC(VB(AD,1),KR,3,0,[new tm(27,59.99313),new tm(28,60.98062),new tm(29,61.9732),new tm(30,62.96369),new tm(31,63.957572),new tm(32,64.949484),new tm(33,65.944099147),new tm(34,66.939190417),new tm(35,67.936792976),new tm(36,68.932280154),new tm(37,69.930927811),new tm(38,70.927114724),new tm(39,71.926752647),new tm(40,72.923825288),new tm(41,73.923929076),new tm(42,74.921596417),new tm(43,75.922393933),new tm(44,76.920647703),new tm(45,77.921828577),new tm(46,78.920948498),new tm(47,79.922578162),new tm(48,80.922132884),new tm(49,81.924504668),new tm(50,82.924980625),new tm(51,83.92906),new tm(52,84.93181),new tm(53,85.93623),new tm(54,86.93958),new tm(55,87.94456),new tm(56,88.94923)]),aC(VB(AD,1),KR,3,0,[new tm(31,64.96466),new tm(32,65.95521),new tm(33,66.95009),new tm(34,67.94187),new tm(35,68.939562155),new tm(36,69.933504),new tm(37,70.931868378),new tm(38,71.927112313),new tm(39,72.9267668),new tm(40,73.922476561),new tm(41,74.922523571),new tm(42,75.919214107),new tm(43,76.91991461),new tm(44,77.917309522),new tm(45,78.918499802),new tm(46,79.916521828),new tm(47,80.917992931),new tm(48,81.9167),new tm(49,82.919119072),new tm(50,83.918464523),new tm(51,84.922244678),new tm(52,85.924271165),new tm(53,86.928520749),new tm(54,87.931423982),new tm(55,88.93602),new tm(56,89.93942),new tm(57,90.94537),new tm(58,91.94933)]),aC(VB(AD,1),KR,3,0,[new tm(32,66.96479),new tm(33,67.958248),new tm(34,68.950178),new tm(35,69.944208),new tm(36,70.939246),new tm(37,71.936496876),new tm(38,72.931794889),new tm(39,73.929891152),new tm(40,74.92577641),new tm(41,75.924541974),new tm(42,76.921380123),new tm(43,77.92114613),new tm(44,78.918337647),new tm(45,79.918529952),new tm(46,80.91629106),new tm(47,81.916804666),new tm(48,82.915180219),new tm(49,83.916503685),new tm(50,84.915608027),new tm(51,85.918797162),new tm(52,86.920710713),new tm(53,87.924065908),new tm(54,88.92638726),new tm(55,89.930634988),new tm(56,90.9339653),new tm(57,91.939255258),new tm(58,92.9431),new tm(59,93.94868)]),aC(VB(AD,1),KR,3,0,[new tm(33,68.96532),new tm(34,69.95601),new tm(35,70.95051),new tm(36,71.94190754),new tm(37,72.938931115),new tm(38,73.933258225),new tm(39,74.931033794),new tm(40,75.925948304),new tm(41,76.92466788),new tm(42,77.920386271),new tm(43,78.920082992),new tm(44,79.91637804),new tm(45,80.916592419),new tm(46,81.913484601),new tm(47,82.914135952),new tm(48,83.911506627),new tm(49,84.912526954),new tm(50,85.910610313),new tm(51,86.913354251),new tm(52,87.914446951),new tm(53,88.917632505),new tm(54,89.919523803),new tm(55,90.923442418),new tm(56,91.926152752),new tm(57,92.931265246),new tm(58,93.934362),new tm(59,94.93984),new tm(60,95.94307),new tm(61,96.94856)]),aC(VB(AD,1),KR,3,0,[new tm(34,70.96532),new tm(35,71.95908),new tm(36,72.950366),new tm(37,73.944470376),new tm(38,74.938569199),new tm(39,75.935071448),new tm(40,76.930406599),new tm(41,77.928141485),new tm(42,78.923996719),new tm(43,79.922519322),new tm(44,80.918994165),new tm(45,81.918207691),new tm(46,82.915111951),new tm(47,83.914384676),new tm(48,84.911789341),new tm(49,85.91116708),new tm(50,86.909183465),new tm(51,87.911318556),new tm(52,88.912279939),new tm(53,89.914808941),new tm(54,90.91653416),new tm(55,91.919725442),new tm(56,92.922032765),new tm(57,93.926407326),new tm(58,94.92931926),new tm(59,95.934283962),new tm(60,96.937342863),new tm(61,97.941703557),new tm(62,98.945420616),new tm(63,99.94987),new tm(64,100.953195994),new tm(65,101.95921)]),aC(VB(AD,1),KR,3,0,[new tm(35,72.96597),new tm(36,73.95631),new tm(37,74.94992),new tm(38,75.94161),new tm(39,76.937761511),new tm(40,77.932179362),new tm(41,78.929707076),new tm(42,79.924524588),new tm(43,80.923213095),new tm(44,81.918401258),new tm(45,82.917555029),new tm(46,83.913424778),new tm(47,84.912932689),new tm(48,85.909262351),new tm(49,86.908879316),new tm(50,87.905614339),new tm(51,88.907452906),new tm(52,89.907737596),new tm(53,90.910209845),new tm(54,91.911029895),new tm(55,92.91402241),new tm(56,93.915359856),new tm(57,94.919358213),new tm(58,95.921680473),new tm(59,96.926148757),new tm(60,97.928471177),new tm(61,98.933315038),new tm(62,99.935351729),new tm(63,100.940517434),new tm(64,101.943018795),new tm(65,102.94895),new tm(66,103.95233)]),aC(VB(AD,1),KR,3,0,[new tm(38,76.94962),new tm(39,77.9435),new tm(40,78.937350712),new tm(41,79.931982402),new tm(42,80.929128719),new tm(43,81.926792071),new tm(44,82.922352572),new tm(45,83.920387768),new tm(46,84.916427076),new tm(47,85.914887724),new tm(48,86.910877833),new tm(49,87.909503361),new tm(50,88.905847902),new tm(51,89.907151443),new tm(52,90.907303415),new tm(53,91.908946832),new tm(54,92.909581582),new tm(55,93.911594008),new tm(56,94.912823709),new tm(57,95.915897787),new tm(58,96.918131017),new tm(59,97.922219525),new tm(60,98.924634736),new tm(61,99.927756402),new tm(62,100.930313395),new tm(63,101.933555501),new tm(64,102.93694),new tm(65,103.94145),new tm(66,104.94509),new tm(67,105.95022)]),aC(VB(AD,1),KR,3,0,[new tm(39,78.94916),new tm(40,79.94055),new tm(41,80.936815296),new tm(42,81.931086249),new tm(43,82.92865213),new tm(44,83.92325),new tm(45,84.92146522),new tm(46,85.916472851),new tm(47,86.914816578),new tm(48,87.910226179),new tm(49,88.908888916),new tm(50,89.904703679),new tm(51,90.905644968),new tm(52,91.905040106),new tm(53,92.906475627),new tm(54,93.906315765),new tm(55,94.908042739),new tm(56,95.908275675),new tm(57,96.910950716),new tm(58,97.912746366),new tm(59,98.916511084),new tm(60,99.917761704),new tm(61,100.921139958),new tm(62,101.922981089),new tm(63,102.926597062),new tm(64,103.92878),new tm(65,104.93305),new tm(66,105.93591),new tm(67,106.94086),new tm(68,107.94428)]),aC(VB(AD,1),KR,3,0,[new tm(40,80.94905),new tm(41,81.94313),new tm(42,82.936703713),new tm(43,83.93357),new tm(44,84.927906486),new tm(45,85.925037588),new tm(46,86.920361435),new tm(47,87.91833144),new tm(48,88.913495503),new tm(49,89.911264109),new tm(50,90.906990538),new tm(51,91.907193214),new tm(52,92.906377543),new tm(53,93.907283457),new tm(54,94.906835178),new tm(55,95.908100076),new tm(56,96.908097144),new tm(57,97.91033069),new tm(58,98.911617864),new tm(59,99.914181434),new tm(60,100.915251567),new tm(61,101.918037417),new tm(62,102.919141297),new tm(63,103.922459464),new tm(64,104.923934023),new tm(65,105.92819),new tm(66,106.93031),new tm(67,107.93501),new tm(68,108.93763),new tm(69,109.94268)]),aC(VB(AD,1),KR,3,0,[new tm(41,82.94874),new tm(42,83.94009),new tm(43,84.93659),new tm(44,85.930695167),new tm(45,86.92732683),new tm(46,87.921952728),new tm(47,88.919480562),new tm(48,89.913936161),new tm(49,90.911750754),new tm(50,91.90681048),new tm(51,92.906812213),new tm(52,93.905087578),new tm(53,94.905841487),new tm(54,95.904678904),new tm(55,96.906021033),new tm(56,97.905407846),new tm(57,98.907711598),new tm(58,99.907477149),new tm(59,100.910346543),new tm(60,101.910297162),new tm(61,102.913204596),new tm(62,103.913758387),new tm(63,104.916972087),new tm(64,105.918134284),new tm(65,106.921694724),new tm(66,107.923973837),new tm(67,108.92781),new tm(68,109.92973),new tm(69,110.93451),new tm(70,111.93684),new tm(71,112.94203)]),aC(VB(AD,1),KR,3,0,[new tm(42,84.94894),new tm(43,85.94288),new tm(44,86.93653),new tm(45,87.93283),new tm(46,88.92754288),new tm(47,89.92355583),new tm(48,90.9184282),new tm(49,91.915259655),new tm(50,92.910248473),new tm(51,93.909656309),new tm(52,94.907656454),new tm(53,95.907870803),new tm(54,96.906364843),new tm(55,97.907215692),new tm(56,98.906254554),new tm(57,99.907657594),new tm(58,100.90731438),new tm(59,101.909212938),new tm(60,102.909178805),new tm(61,103.911444898),new tm(62,104.911658043),new tm(63,105.914355408),new tm(64,106.915081691),new tm(65,107.918479973),new tm(66,108.919980998),new tm(67,109.92339),new tm(68,110.92505),new tm(69,111.92924),new tm(70,112.93133),new tm(71,113.93588),new tm(72,114.93828)]),aC(VB(AD,1),KR,3,0,[new tm(43,86.94918),new tm(44,87.94042),new tm(45,88.93611),new tm(46,89.92978),new tm(47,90.926377434),new tm(48,91.92012),new tm(49,92.917051523),new tm(50,93.911359569),new tm(51,94.910412729),new tm(52,95.907597681),new tm(53,96.907554546),new tm(54,97.905287111),new tm(55,98.905939307),new tm(56,99.904219664),new tm(57,100.905582219),new tm(58,101.904349503),new tm(59,102.906323677),new tm(60,103.905430145),new tm(61,104.907750341),new tm(62,105.907326913),new tm(63,106.909907207),new tm(64,107.910192211),new tm(65,108.913201565),new tm(66,109.913966185),new tm(67,110.91756),new tm(68,111.918821673),new tm(69,112.92254),new tm(70,113.923891981),new tm(71,114.92831),new tm(72,115.93016),new tm(73,116.93479),new tm(74,117.93703)]),aC(VB(AD,1),KR,3,0,[new tm(44,88.94938),new tm(45,89.94287),new tm(46,90.93655),new tm(47,91.93198),new tm(48,92.92574),new tm(49,93.921698),new tm(50,94.915898541),new tm(51,95.914518212),new tm(52,96.911336643),new tm(53,97.910716431),new tm(54,98.908132101),new tm(55,99.90811663),new tm(56,100.906163526),new tm(57,101.906842845),new tm(58,102.905504182),new tm(59,103.906655315),new tm(60,104.905692444),new tm(61,105.907284615),new tm(62,106.90675054),new tm(63,107.908730768),new tm(64,108.908735621),new tm(65,109.910949525),new tm(66,110.91166),new tm(67,111.913969253),new tm(68,112.91542),new tm(69,113.91734336),new tm(70,114.920124676),new tm(71,115.922746643),new tm(72,116.92535),new tm(73,117.92943),new tm(74,118.93136),new tm(75,119.93578),new tm(76,120.93808)]),aC(VB(AD,1),KR,3,0,[new tm(45,90.94948),new tm(46,91.94042),new tm(47,92.93591),new tm(48,93.92877),new tm(49,94.92469),new tm(50,95.91822194),new tm(51,96.916478921),new tm(52,97.912720751),new tm(53,98.911767757),new tm(54,99.908504596),new tm(55,100.908289144),new tm(56,101.905607716),new tm(57,102.906087204),new tm(58,103.904034912),new tm(59,104.905084046),new tm(60,105.903483087),new tm(61,106.905128453),new tm(62,107.903894451),new tm(63,108.905953535),new tm(64,109.905152385),new tm(65,110.907643952),new tm(66,111.907313277),new tm(67,112.910151346),new tm(68,113.910365322),new tm(69,114.91368341),new tm(70,115.914158288),new tm(71,116.91784),new tm(72,117.918983915),new tm(73,118.92268),new tm(74,119.92403),new tm(75,120.92818),new tm(76,121.9298),new tm(77,122.93426)]),aC(VB(AD,1),KR,3,0,[new tm(47,93.94278),new tm(48,94.93548),new tm(49,95.93068),new tm(50,96.924),new tm(51,97.921759995),new tm(52,98.917597103),new tm(53,99.916069387),new tm(54,100.912802135),new tm(55,101.911999996),new tm(56,102.908972453),new tm(57,103.908628228),new tm(58,104.906528234),new tm(59,105.906666431),new tm(60,106.90509302),new tm(61,107.905953705),new tm(62,108.904755514),new tm(63,109.90611046),new tm(64,110.905294679),new tm(65,111.907004132),new tm(66,112.906565708),new tm(67,113.908807907),new tm(68,114.908762282),new tm(69,115.911359558),new tm(70,116.911684187),new tm(71,117.914582383),new tm(72,118.915666045),new tm(73,119.918788609),new tm(74,120.919851074),new tm(75,121.92332),new tm(76,122.9249),new tm(77,123.92853),new tm(78,124.93054),new tm(79,125.9345),new tm(80,126.93688)]),aC(VB(AD,1),KR,3,0,[new tm(48,95.93977),new tm(49,96.93494),new tm(50,97.927579),new tm(51,98.92501),new tm(52,99.920230232),new tm(53,100.918681442),new tm(54,101.914777255),new tm(55,102.913418952),new tm(56,103.909848091),new tm(57,104.909467818),new tm(58,105.906458007),new tm(59,106.906614232),new tm(60,107.904183403),new tm(61,108.904985569),new tm(62,109.903005578),new tm(63,110.904181628),new tm(64,111.902757226),new tm(65,112.904400947),new tm(66,113.903358121),new tm(67,114.905430553),new tm(68,115.904755434),new tm(69,116.907218242),new tm(70,117.906914144),new tm(71,118.909922582),new tm(72,119.909851352),new tm(73,120.91298039),new tm(74,121.9135),new tm(75,122.917003675),new tm(76,123.917648302),new tm(77,124.92124717),new tm(78,125.922353996),new tm(79,126.926434822),new tm(80,127.927760617),new tm(81,128.93226),new tm(82,129.93398)]),aC(VB(AD,1),KR,3,0,[new tm(49,97.94224),new tm(50,98.93461),new tm(51,99.931149033),new tm(52,100.92656),new tm(53,101.924707541),new tm(54,102.919913896),new tm(55,103.918338416),new tm(56,104.914673434),new tm(57,105.913461134),new tm(58,106.910292195),new tm(59,107.909719683),new tm(60,108.907154078),new tm(61,109.907168783),new tm(62,110.905110677),new tm(63,111.905533338),new tm(64,112.904061223),new tm(65,113.904916758),new tm(66,114.903878328),new tm(67,115.905259995),new tm(68,116.904515731),new tm(69,117.906354623),new tm(70,118.905846334),new tm(71,119.907961505),new tm(72,120.907848847),new tm(73,121.910277103),new tm(74,122.910438951),new tm(75,123.913175916),new tm(76,124.913601387),new tm(77,125.916464532),new tm(78,126.917344048),new tm(79,127.920170658),new tm(80,128.921657958),new tm(81,129.924854941),new tm(82,130.926767408),new tm(83,131.932919005),new tm(84,132.93834),new tm(85,133.94466)]),aC(VB(AD,1),KR,3,0,[new tm(50,99.938954),new tm(51,100.93606),new tm(52,101.93049),new tm(53,102.92813),new tm(54,103.923185469),new tm(55,104.921390409),new tm(56,105.916880472),new tm(57,106.915666702),new tm(58,107.911965339),new tm(59,108.911286879),new tm(60,109.907852688),new tm(61,110.907735404),new tm(62,111.90482081),new tm(63,112.905173373),new tm(64,113.902781816),new tm(65,114.903345973),new tm(66,115.901744149),new tm(67,116.902953765),new tm(68,117.901606328),new tm(69,118.90330888),new tm(70,119.902196571),new tm(71,120.904236867),new tm(72,121.903440138),new tm(73,122.905721901),new tm(74,123.90527463),new tm(75,124.907784924),new tm(76,125.907653953),new tm(77,126.91035098),new tm(78,127.910534953),new tm(79,128.913439976),new tm(80,129.913852185),new tm(81,130.916919144),new tm(82,131.917744455),new tm(83,132.923814085),new tm(84,133.928463576),new tm(85,134.93473),new tm(86,135.93934),new tm(87,136.94579)]),aC(VB(AD,1),KR,3,0,[new tm(52,102.94012),new tm(53,103.936287),new tm(54,104.931528593),new tm(55,105.928183134),new tm(56,106.92415),new tm(57,107.92216),new tm(58,108.918136092),new tm(59,109.917533911),new tm(60,110.912534147),new tm(61,111.91239464),new tm(62,112.909377941),new tm(63,113.909095876),new tm(64,114.906598812),new tm(65,115.906797235),new tm(66,116.90483959),new tm(67,117.905531885),new tm(68,118.90394646),new tm(69,119.905074315),new tm(70,120.903818044),new tm(71,121.905175415),new tm(72,122.904215696),new tm(73,123.905937525),new tm(74,124.905247804),new tm(75,125.907248153),new tm(76,126.906914564),new tm(77,127.90916733),new tm(78,128.909150092),new tm(79,129.911546459),new tm(80,130.911946487),new tm(81,131.914413247),new tm(82,132.915236466),new tm(83,133.920551554),new tm(84,134.925167962),new tm(85,135.93066),new tm(86,136.93531),new tm(87,137.94096),new tm(88,138.94571)]),aC(VB(AD,1),KR,3,0,[new tm(54,105.937702),new tm(55,106.935036),new tm(56,107.929486838),new tm(57,108.927456483),new tm(58,109.922407164),new tm(59,110.921120589),new tm(60,111.917061617),new tm(61,112.915452551),new tm(62,113.912498025),new tm(63,114.911578627),new tm(64,115.908420253),new tm(65,116.90863418),new tm(66,117.905825187),new tm(67,118.90640811),new tm(68,119.904019891),new tm(69,120.904929815),new tm(70,121.903047064),new tm(71,122.904272951),new tm(72,123.902819466),new tm(73,124.904424718),new tm(74,125.903305543),new tm(75,126.90521729),new tm(76,127.904461383),new tm(77,128.906595593),new tm(78,129.906222753),new tm(79,130.90852188),new tm(80,131.908523782),new tm(81,132.910939068),new tm(82,133.911540546),new tm(83,134.916450782),new tm(84,135.920103155),new tm(85,136.925324769),new tm(86,137.92922),new tm(87,138.93473),new tm(88,139.9387),new tm(89,140.94439),new tm(90,141.9485)]),aC(VB(AD,1),KR,3,0,[new tm(55,107.943291),new tm(56,108.938191658),new tm(57,109.934634181),new tm(58,110.930276),new tm(59,111.92797),new tm(60,112.923644245),new tm(61,113.92185),new tm(62,114.918272),new tm(63,115.916735014),new tm(64,116.913647692),new tm(65,117.91337523),new tm(66,118.910180837),new tm(67,119.910047843),new tm(68,120.907366063),new tm(69,121.907592451),new tm(70,122.905597944),new tm(71,123.906211423),new tm(72,124.90462415),new tm(73,125.905619387),new tm(74,126.90446842),new tm(75,127.905805254),new tm(76,128.904987487),new tm(77,129.906674018),new tm(78,130.906124168),new tm(79,131.907994525),new tm(80,132.907806465),new tm(81,133.909876552),new tm(82,134.91005031),new tm(83,135.914655105),new tm(84,136.917872653),new tm(85,137.922383666),new tm(86,138.926093402),new tm(87,139.93121),new tm(88,140.93483),new tm(89,141.94018),new tm(90,142.94407),new tm(91,143.94961)]),aC(VB(AD,1),KR,3,0,[new tm(56,109.944476),new tm(57,110.941632),new tm(58,111.93566535),new tm(59,112.933382836),new tm(60,113.928145),new tm(61,114.926979032),new tm(62,115.921394197),new tm(63,116.920564355),new tm(64,117.91657092),new tm(65,118.915554295),new tm(66,119.91215199),new tm(67,120.911386497),new tm(68,121.908548396),new tm(69,122.908470748),new tm(70,123.905895774),new tm(71,124.906398236),new tm(72,125.904268868),new tm(73,126.905179581),new tm(74,127.903530436),new tm(75,128.904779458),new tm(76,129.903507903),new tm(77,130.90508192),new tm(78,131.904154457),new tm(79,132.90590566),new tm(80,133.905394504),new tm(81,134.907207499),new tm(82,135.907219526),new tm(83,136.911562939),new tm(84,137.913988549),new tm(85,138.918786859),new tm(86,139.921635665),new tm(87,140.926646282),new tm(88,141.929702981),new tm(89,142.93489),new tm(90,143.93823),new tm(91,144.94367),new tm(92,145.9473),new tm(93,146.95301)]),aC(VB(AD,1),KR,3,0,[new tm(57,111.950331),new tm(58,112.944535512),new tm(59,113.940841319),new tm(60,114.935939),new tm(61,115.932914152),new tm(62,116.928639484),new tm(63,117.926554883),new tm(64,118.922370879),new tm(65,119.920678219),new tm(66,120.917183637),new tm(67,121.916121946),new tm(68,122.912990168),new tm(69,123.912245731),new tm(70,124.909724871),new tm(71,125.909447953),new tm(72,126.9074176),new tm(73,127.907747919),new tm(74,128.906063369),new tm(75,129.906706163),new tm(76,130.905460232),new tm(77,131.906429799),new tm(78,132.90544687),new tm(79,133.906713419),new tm(80,134.905971903),new tm(81,135.907305741),new tm(82,136.907083505),new tm(83,137.911010537),new tm(84,138.913357921),new tm(85,139.917277075),new tm(86,140.920043984),new tm(87,141.924292317),new tm(88,142.927330292),new tm(89,143.932027373),new tm(90,144.935388226),new tm(91,145.940162028),new tm(92,146.943864435),new tm(93,147.948899539),new tm(94,148.95272),new tm(95,149.95797),new tm(96,150.962)]),aC(VB(AD,1),KR,3,0,[new tm(58,113.950941),new tm(59,114.94771),new tm(60,115.94168),new tm(61,116.937700229),new tm(62,117.93344),new tm(63,118.931051927),new tm(64,119.926045941),new tm(65,120.924485908),new tm(66,121.92026),new tm(67,122.91885),new tm(68,123.915088437),new tm(69,124.914620234),new tm(70,125.911244146),new tm(71,126.911121328),new tm(72,127.90830887),new tm(73,128.908673749),new tm(74,129.906310478),new tm(75,130.906930798),new tm(76,131.905056152),new tm(77,132.906002368),new tm(78,133.904503347),new tm(79,134.905682749),new tm(80,135.904570109),new tm(81,136.905821414),new tm(82,137.905241273),new tm(83,138.908835384),new tm(84,139.910599485),new tm(85,140.914406439),new tm(86,141.916448175),new tm(87,142.920617184),new tm(88,143.922940468),new tm(89,144.926923807),new tm(90,145.930106645),new tm(91,146.933992519),new tm(92,147.937682377),new tm(93,148.94246),new tm(94,149.94562),new tm(95,150.9507),new tm(96,151.95416),new tm(97,152.95961)]),aC(VB(AD,1),KR,3,0,[new tm(60,116.95001),new tm(61,117.94657),new tm(62,118.94099),new tm(63,119.93807),new tm(64,120.93301),new tm(65,121.93071),new tm(66,122.92624),new tm(67,123.92453),new tm(68,124.92067),new tm(69,125.91937),new tm(70,126.91616),new tm(71,127.91544794),new tm(72,128.912667334),new tm(73,129.91232),new tm(74,130.910108489),new tm(75,131.910110399),new tm(76,132.908396372),new tm(77,133.908489607),new tm(78,134.906971003),new tm(79,135.907651181),new tm(80,136.906465656),new tm(81,137.907106826),new tm(82,138.90634816),new tm(83,139.909472552),new tm(84,140.910957016),new tm(85,141.914074489),new tm(86,142.916058646),new tm(87,143.919591666),new tm(88,144.92163837),new tm(89,145.925700146),new tm(90,146.927819639),new tm(91,147.932191197),new tm(92,148.93437),new tm(93,149.93857),new tm(94,150.94156),new tm(95,151.94611),new tm(96,152.94945),new tm(97,153.9544),new tm(98,154.95813)]),aC(VB(AD,1),KR,3,0,[new tm(61,118.95276),new tm(62,119.94664),new tm(63,120.94367),new tm(64,121.93801),new tm(65,122.93551),new tm(66,123.93052),new tm(67,124.92854),new tm(68,125.9241),new tm(69,126.92275),new tm(70,127.91887),new tm(71,128.918679183),new tm(72,129.914339361),new tm(73,130.914424137),new tm(74,131.91149),new tm(75,132.91155),new tm(76,133.909026379),new tm(77,134.909145555),new tm(78,135.907143574),new tm(79,136.907777634),new tm(80,137.905985574),new tm(81,138.906646605),new tm(82,139.905434035),new tm(83,140.908271103),new tm(84,141.909239733),new tm(85,142.912381158),new tm(86,143.913642686),new tm(87,144.917227871),new tm(88,145.918689722),new tm(89,146.922510962),new tm(90,147.924394738),new tm(91,148.928289207),new tm(92,149.930226399),new tm(93,150.93404),new tm(94,151.93638),new tm(95,152.94058),new tm(96,153.94332),new tm(97,154.94804),new tm(98,155.95126),new tm(99,156.95634)]),aC(VB(AD,1),KR,3,0,[new tm(62,120.955364),new tm(63,121.95165),new tm(64,122.94596),new tm(65,123.94296),new tm(66,124.93783),new tm(67,125.93531),new tm(68,126.93083),new tm(69,127.9288),new tm(70,128.92486),new tm(71,129.92338),new tm(72,130.920060245),new tm(73,131.91912),new tm(74,132.9162),new tm(75,133.915672),new tm(76,134.91313914),new tm(77,135.912646935),new tm(78,136.910678351),new tm(79,137.910748891),new tm(80,138.908932181),new tm(81,139.909071204),new tm(82,140.907647726),new tm(83,141.910039865),new tm(84,142.910812233),new tm(85,143.913300595),new tm(86,144.914506897),new tm(87,145.917588016),new tm(88,146.918979001),new tm(89,147.922183237),new tm(90,148.923791056),new tm(91,149.926995031),new tm(92,150.928227869),new tm(93,151.9316),new tm(94,152.93365),new tm(95,153.93739),new tm(96,154.93999),new tm(97,155.94412),new tm(98,156.94717),new tm(99,157.95178),new tm(100,158.95523)]),aC(VB(AD,1),KR,3,0,[new tm(66,125.94307),new tm(67,126.9405),new tm(68,127.93539),new tm(69,128.932385),new tm(70,129.92878),new tm(71,130.927102697),new tm(72,131.92312),new tm(73,132.92221),new tm(74,133.918645),new tm(75,134.91824),new tm(76,135.915020542),new tm(77,136.91463973),new tm(78,137.91291745),new tm(79,138.91192415),new tm(80,139.909309824),new tm(81,140.9096048),new tm(82,141.907718643),new tm(83,142.909809626),new tm(84,143.910082629),new tm(85,144.912568847),new tm(86,145.913112139),new tm(87,146.916095794),new tm(88,147.916888516),new tm(89,148.92014419),new tm(90,149.920886563),new tm(91,150.923824739),new tm(92,151.924682428),new tm(93,152.927694534),new tm(94,153.929483295),new tm(95,154.932629551),new tm(96,155.9352),new tm(97,156.93927),new tm(98,157.94187),new tm(99,158.94639),new tm(100,159.94939),new tm(101,160.95433)]),aC(VB(AD,1),KR,3,0,[new tm(67,127.94826),new tm(68,128.94316),new tm(69,129.94045),new tm(70,130.9358),new tm(71,131.93375),new tm(72,132.92972),new tm(73,133.92849),new tm(74,134.924617),new tm(75,135.923447865),new tm(76,136.920713),new tm(77,137.920432261),new tm(78,138.916759814),new tm(79,139.915801649),new tm(80,140.913606636),new tm(81,141.912950738),new tm(82,142.910927571),new tm(83,143.912585768),new tm(84,144.912743879),new tm(85,145.914692165),new tm(86,146.915133898),new tm(87,147.917467786),new tm(88,148.918329195),new tm(89,149.920979477),new tm(90,150.921202693),new tm(91,151.923490557),new tm(92,152.924113189),new tm(93,153.926547019),new tm(94,154.928097047),new tm(95,155.931060357),new tm(96,156.9332),new tm(97,157.93669),new tm(98,158.93913),new tm(99,159.94299),new tm(100,160.94586),new tm(101,161.95029),new tm(102,162.95352)]),aC(VB(AD,1),KR,3,0,[new tm(68,129.94863),new tm(69,130.94589),new tm(70,131.94082),new tm(71,132.93873),new tm(72,133.93402),new tm(73,134.93235),new tm(74,135.9283),new tm(75,136.927046709),new tm(76,137.92354),new tm(77,138.922302),new tm(78,139.918991),new tm(79,140.918468512),new tm(80,141.915193274),new tm(81,142.914623555),new tm(82,143.91199473),new tm(83,144.913405611),new tm(84,145.91303676),new tm(85,146.914893275),new tm(86,147.914817914),new tm(87,148.917179521),new tm(88,149.917271454),new tm(89,150.919928351),new tm(90,151.919728244),new tm(91,152.922093907),new tm(92,153.922205303),new tm(93,154.92463594),new tm(94,155.925526236),new tm(95,156.928354506),new tm(96,157.929987938),new tm(97,158.9332),new tm(98,159.93514),new tm(99,160.93883),new tm(100,161.94122),new tm(101,162.94536),new tm(102,163.94828),new tm(103,164.95298)]),aC(VB(AD,1),KR,3,0,[new tm(69,131.95416),new tm(70,132.9489),new tm(71,133.94632),new tm(72,134.94172),new tm(73,135.9395),new tm(74,136.93521),new tm(75,137.93345),new tm(76,138.92882915),new tm(77,139.928083921),new tm(78,140.924885867),new tm(79,141.923400033),new tm(80,142.920286634),new tm(81,143.918774116),new tm(82,144.916261285),new tm(83,145.917199714),new tm(84,146.916741206),new tm(85,147.918153775),new tm(86,148.917925922),new tm(87,149.919698294),new tm(88,150.919846022),new tm(89,151.921740399),new tm(90,152.921226219),new tm(91,153.922975386),new tm(92,154.922889429),new tm(93,155.924750855),new tm(94,156.925419435),new tm(95,157.927841923),new tm(96,158.9290845),new tm(97,159.931460406),new tm(98,160.93368),new tm(99,161.93704),new tm(100,162.93921),new tm(101,163.94299),new tm(102,164.94572),new tm(103,165.94997),new tm(104,166.95305)]),aC(VB(AD,1),KR,3,0,[new tm(72,135.94707),new tm(73,136.94465),new tm(74,137.93997),new tm(75,138.93808),new tm(76,139.933236934),new tm(77,140.93221),new tm(78,141.927908919),new tm(79,142.926738636),new tm(80,143.923390357),new tm(81,144.921687498),new tm(82,145.918305344),new tm(83,146.919089446),new tm(84,147.918109771),new tm(85,148.919336427),new tm(86,149.918655455),new tm(87,150.920344273),new tm(88,151.919787882),new tm(89,152.921746283),new tm(90,153.920862271),new tm(91,154.922618801),new tm(92,155.922119552),new tm(93,156.923956686),new tm(94,157.924100533),new tm(95,158.926385075),new tm(96,159.927050616),new tm(97,160.929665688),new tm(98,161.930981211),new tm(99,162.93399),new tm(100,163.93586),new tm(101,164.93938),new tm(102,165.9416),new tm(103,166.94557),new tm(104,167.94836),new tm(105,168.95287)]),aC(VB(AD,1),KR,3,0,[new tm(73,137.95287),new tm(74,138.94803),new tm(75,139.945367985),new tm(76,140.94116),new tm(77,141.939073781),new tm(78,142.93475),new tm(79,143.93253),new tm(80,144.92888),new tm(81,145.927180629),new tm(82,146.924037176),new tm(83,147.924298636),new tm(84,148.92324163),new tm(85,149.923654158),new tm(86,150.923098169),new tm(87,151.924071324),new tm(88,152.923430858),new tm(89,153.924686236),new tm(90,154.923500411),new tm(91,155.924743749),new tm(92,156.924021155),new tm(93,157.92541026),new tm(94,158.925343135),new tm(95,159.927164021),new tm(96,160.927566289),new tm(97,161.929484803),new tm(98,162.930643942),new tm(99,163.933347253),new tm(100,164.93488),new tm(101,165.93805),new tm(102,166.94005),new tm(103,167.94364),new tm(104,168.94622),new tm(105,169.95025),new tm(106,170.9533)]),aC(VB(AD,1),KR,3,0,[new tm(74,139.95379),new tm(75,140.95119),new tm(76,141.946695946),new tm(77,142.94383),new tm(78,143.93907),new tm(79,144.936717),new tm(80,145.932720118),new tm(81,146.930878496),new tm(82,147.927177882),new tm(83,148.927333981),new tm(84,149.925579728),new tm(85,150.92617963),new tm(86,151.924713874),new tm(87,152.925760865),new tm(88,153.924422046),new tm(89,154.92574895),new tm(90,155.924278273),new tm(91,156.925461256),new tm(92,157.924404637),new tm(93,158.92573566),new tm(94,159.925193718),new tm(95,160.926929595),new tm(96,161.926794731),new tm(97,162.928727532),new tm(98,163.929171165),new tm(99,164.931699828),new tm(100,165.932803241),new tm(101,166.935649025),new tm(102,167.93723),new tm(103,168.940303648),new tm(104,169.94267),new tm(105,170.94648),new tm(106,171.94911),new tm(107,172.95344)]),aC(VB(AD,1),KR,3,0,[new tm(75,141.95986),new tm(76,142.95469),new tm(77,143.95164),new tm(78,144.94688),new tm(79,145.9441),new tm(80,146.93984),new tm(81,147.937269),new tm(82,148.933789944),new tm(83,149.932760914),new tm(84,150.931680791),new tm(85,151.931740598),new tm(86,152.930194506),new tm(87,153.930596268),new tm(88,154.929079084),new tm(89,155.929001869),new tm(90,156.928188059),new tm(91,157.92894573),new tm(92,158.927708537),new tm(93,159.928725679),new tm(94,160.927851662),new tm(95,161.92909242),new tm(96,162.928730286),new tm(97,163.930230577),new tm(98,164.930319169),new tm(99,165.932281267),new tm(100,166.933126195),new tm(101,167.935496424),new tm(102,168.936868306),new tm(103,169.939614951),new tm(104,170.941461227),new tm(105,171.94482),new tm(106,172.94729),new tm(107,173.95115),new tm(108,174.95405)]),aC(VB(AD,1),KR,3,0,[new tm(76,143.96059),new tm(77,144.95746),new tm(78,145.95212),new tm(79,146.94931),new tm(80,147.94444),new tm(81,148.942780527),new tm(82,149.937171034),new tm(83,150.93746),new tm(84,151.935078452),new tm(85,152.935093125),new tm(86,153.932777294),new tm(87,154.933204273),new tm(88,155.931015001),new tm(89,156.931945517),new tm(90,157.929912),new tm(91,158.930680718),new tm(92,159.929078924),new tm(93,160.930001348),new tm(94,161.928774923),new tm(95,162.930029273),new tm(96,163.929196996),new tm(97,164.9307228),new tm(98,165.93028997),new tm(99,166.932045448),new tm(100,167.932367781),new tm(101,168.934588082),new tm(102,169.935460334),new tm(103,170.938025885),new tm(104,171.939352149),new tm(105,172.9424),new tm(106,173.94434),new tm(107,174.94793),new tm(108,175.95029),new tm(109,176.95437)]),aC(VB(AD,1),KR,3,0,[new tm(77,145.966495),new tm(78,146.961081),new tm(79,147.95755),new tm(80,148.95265),new tm(81,149.94967),new tm(82,150.944842),new tm(83,151.9443),new tm(84,152.942027631),new tm(85,153.940832325),new tm(86,154.939191562),new tm(87,155.939006895),new tm(88,156.936756069),new tm(89,157.936996),new tm(90,158.934808966),new tm(91,159.935090772),new tm(92,160.933398042),new tm(93,161.933970147),new tm(94,162.932647648),new tm(95,163.933450972),new tm(96,164.932432463),new tm(97,165.933553133),new tm(98,166.932848844),new tm(99,167.934170375),new tm(100,168.934211117),new tm(101,169.935797877),new tm(102,170.936425817),new tm(103,171.938396118),new tm(104,172.939600336),new tm(105,173.942164618),new tm(106,174.943832897),new tm(107,175.946991412),new tm(108,176.94904),new tm(109,177.95264),new tm(110,178.95534)]),aC(VB(AD,1),KR,3,0,[new tm(78,147.96676),new tm(79,148.96348),new tm(80,149.95799),new tm(81,150.954657965),new tm(82,151.950167),new tm(83,152.94921),new tm(84,153.945651145),new tm(85,154.945792),new tm(86,155.942847109),new tm(87,156.94265865),new tm(88,157.939857897),new tm(89,158.940153735),new tm(90,159.93756),new tm(91,160.937357719),new tm(92,161.93575),new tm(93,162.936265492),new tm(94,163.93452),new tm(95,164.935397592),new tm(96,165.933879623),new tm(97,166.934946862),new tm(98,167.933894465),new tm(99,168.93518712),new tm(100,169.934758652),new tm(101,170.936322297),new tm(102,171.936377696),new tm(103,172.938206756),new tm(104,173.938858101),new tm(105,174.941272494),new tm(106,175.942568409),new tm(107,176.945257126),new tm(108,177.946643396),new tm(109,178.95017),new tm(110,179.95233),new tm(111,180.95615)]),aC(VB(AD,1),KR,3,0,[new tm(79,149.972668),new tm(80,150.967147),new tm(81,151.96361),new tm(82,152.95869),new tm(83,153.9571),new tm(84,154.953641324),new tm(85,155.952907),new tm(86,156.950101536),new tm(87,157.948577981),new tm(88,158.946615113),new tm(89,159.945383),new tm(90,160.943047504),new tm(91,161.943222),new tm(92,162.941203796),new tm(93,163.941215),new tm(94,164.939605886),new tm(95,165.939762646),new tm(96,166.938307056),new tm(97,167.938698576),new tm(98,168.937648757),new tm(99,169.93847219),new tm(100,170.937909903),new tm(101,171.939082239),new tm(102,172.938926901),new tm(103,173.940333522),new tm(104,174.940767904),new tm(105,175.942682399),new tm(106,176.943754987),new tm(107,177.945951366),new tm(108,178.947324216),new tm(109,179.949879968),new tm(110,180.95197),new tm(111,181.95521),new tm(112,182.95757),new tm(113,183.96117)]),aC(VB(AD,1),KR,3,0,[new tm(82,153.96425),new tm(83,154.96276),new tm(84,155.959247),new tm(85,156.958127),new tm(86,157.95405528),new tm(87,158.954003),new tm(88,159.950713588),new tm(89,160.950330852),new tm(90,161.947202977),new tm(91,162.947057),new tm(92,163.944422),new tm(93,164.94454),new tm(94,165.94225),new tm(95,166.9426),new tm(96,167.94063),new tm(97,168.941158567),new tm(98,169.93965),new tm(99,170.94049),new tm(100,171.93945798),new tm(101,172.94065),new tm(102,173.940040159),new tm(103,174.941502991),new tm(104,175.941401828),new tm(105,176.943220013),new tm(106,177.943697732),new tm(107,178.945815073),new tm(108,179.94654876),new tm(109,180.949099124),new tm(110,181.950552893),new tm(111,182.953531012),new tm(112,183.95544788),new tm(113,184.95878),new tm(114,185.96092)]),aC(VB(AD,1),KR,3,0,[new tm(83,155.971689),new tm(84,156.968145),new tm(85,157.966368),new tm(86,158.96232309),new tm(87,159.961358),new tm(88,160.958372992),new tm(89,161.956556553),new tm(90,162.95431665),new tm(91,163.95357),new tm(92,164.950817),new tm(93,165.95047),new tm(94,166.948639),new tm(95,167.947787),new tm(96,168.94592),new tm(97,169.94609),new tm(98,170.94446),new tm(99,171.944739818),new tm(100,172.94459),new tm(101,173.944167937),new tm(102,174.94365),new tm(103,175.944740551),new tm(104,176.944471766),new tm(105,177.945750349),new tm(106,178.945934113),new tm(107,179.947465655),new tm(108,180.947996346),new tm(109,181.950152414),new tm(110,182.951373188),new tm(111,183.954009331),new tm(112,184.955559086),new tm(113,185.9585501),new tm(114,186.96041),new tm(115,187.96371)]),aC(VB(AD,1),KR,3,0,[new tm(84,157.973939),new tm(85,158.97228),new tm(86,159.968369),new tm(87,160.967089),new tm(88,161.962750303),new tm(89,162.962532),new tm(90,163.95898381),new tm(91,164.958335962),new tm(92,165.955019896),new tm(93,166.954672),new tm(94,167.951863),new tm(95,168.951759),new tm(96,169.948473988),new tm(97,170.94946),new tm(98,171.948228837),new tm(99,172.948884),new tm(100,173.94616),new tm(101,174.94677),new tm(102,175.94559),new tm(103,176.94662),new tm(104,177.945848364),new tm(105,178.947071733),new tm(106,179.946705734),new tm(107,180.948198054),new tm(108,181.948205519),new tm(109,182.950224458),new tm(110,183.950932553),new tm(111,184.953420586),new tm(112,185.954362204),new tm(113,186.957158365),new tm(114,187.958486954),new tm(115,188.96191222),new tm(116,189.963179541)]),aC(VB(AD,1),KR,3,0,[new tm(85,159.981485),new tm(86,160.977661),new tm(87,161.975707),new tm(88,162.971375872),new tm(89,163.970319),new tm(90,164.967050268),new tm(91,165.965211372),new tm(92,166.962564),new tm(93,167.961609),new tm(94,168.95883),new tm(95,169.958163),new tm(96,170.955547),new tm(97,171.955285),new tm(98,172.953062),new tm(99,173.952114),new tm(100,174.951393),new tm(101,175.95157),new tm(102,176.95027),new tm(103,177.950851081),new tm(104,178.949981038),new tm(105,179.95078768),new tm(106,180.950064596),new tm(107,181.951211444),new tm(108,182.950821349),new tm(109,183.952524289),new tm(110,184.952955747),new tm(111,185.954986529),new tm(112,186.955750787),new tm(113,187.958112287),new tm(114,188.959228359),new tm(115,189.961816139),new tm(116,190.963123592),new tm(117,191.96596)]),aC(VB(AD,1),KR,3,0,[new tm(86,161.983819),new tm(87,162.982048),new tm(88,163.977927),new tm(89,164.976475),new tm(90,165.971934911),new tm(91,166.971554),new tm(92,167.967832911),new tm(93,168.967076205),new tm(94,169.963569716),new tm(95,170.96304),new tm(96,171.960078),new tm(97,172.959791),new tm(98,173.956307704),new tm(99,174.95708),new tm(100,175.953757941),new tm(101,176.955045),new tm(102,177.953348225),new tm(103,178.953951),new tm(104,179.952308241),new tm(105,180.953274494),new tm(106,181.952186222),new tm(107,182.95311),new tm(108,183.952490808),new tm(109,184.954043023),new tm(110,185.953838355),new tm(111,186.955747928),new tm(112,187.955835993),new tm(113,188.958144866),new tm(114,189.95844521),new tm(115,190.960927951),new tm(116,191.961479047),new tm(117,192.964148083),new tm(118,193.965179314),new tm(119,194.968123889),new tm(120,195.96962255)]),aC(VB(AD,1),KR,3,0,[new tm(88,164.98758),new tm(89,165.985506),new tm(90,166.980951577),new tm(91,167.979966),new tm(92,168.976390868),new tm(93,169.974441697),new tm(94,170.971779),new tm(95,171.970643),new tm(96,172.967707),new tm(97,173.966804),new tm(98,174.964279),new tm(99,175.963511),new tm(100,176.96117),new tm(101,177.960084944),new tm(102,178.95915),new tm(103,179.958555615),new tm(104,180.957642156),new tm(105,181.958127689),new tm(106,182.956814),new tm(107,183.957388318),new tm(108,184.95659),new tm(109,185.957951104),new tm(110,186.95736083),new tm(111,187.958851962),new tm(112,188.958716473),new tm(113,189.960592299),new tm(114,190.960591191),new tm(115,191.962602198),new tm(116,192.9629237),new tm(117,193.96507561),new tm(118,194.9659768),new tm(119,195.968379906),new tm(120,196.969636496),new tm(121,197.97228),new tm(122,198.973787159)]),aC(VB(AD,1),KR,3,0,[new tm(90,167.988035),new tm(91,168.986421),new tm(92,169.981734918),new tm(93,170.981251),new tm(94,171.977376138),new tm(95,172.976499642),new tm(96,173.972811276),new tm(97,174.972276),new tm(98,175.969),new tm(99,176.968453),new tm(100,177.964894223),new tm(101,178.965475),new tm(102,179.962023729),new tm(103,180.963177),new tm(104,181.961267637),new tm(105,182.961729),new tm(106,183.959851685),new tm(107,184.960753782),new tm(108,185.959432346),new tm(109,186.960697),new tm(110,187.959395697),new tm(111,188.9608319),new tm(112,189.959930073),new tm(113,190.961684653),new tm(114,191.961035158),new tm(115,192.962984504),new tm(116,193.962663581),new tm(117,194.964774449),new tm(118,195.964934884),new tm(119,196.967323401),new tm(120,197.967876009),new tm(121,198.970576213),new tm(122,199.971423885),new tm(123,200.974496467),new tm(124,201.97574)]),aC(VB(AD,1),KR,3,0,[new tm(92,170.991183),new tm(93,171.990109),new tm(94,172.986398138),new tm(95,173.984325861),new tm(96,174.981552),new tm(97,175.980269),new tm(98,176.977215),new tm(99,177.975975),new tm(100,178.973412),new tm(101,179.972396),new tm(102,180.969948),new tm(103,181.968621416),new tm(104,182.96762),new tm(105,183.966776046),new tm(106,184.965806956),new tm(107,185.965997671),new tm(108,186.964562),new tm(109,187.965321662),new tm(110,188.9642243),new tm(111,189.964698757),new tm(112,190.963649239),new tm(113,191.964810107),new tm(114,192.964131745),new tm(115,193.96533889),new tm(116,194.965017928),new tm(117,195.966551315),new tm(118,196.966551609),new tm(119,197.968225244),new tm(120,198.968748016),new tm(121,199.970717886),new tm(122,200.971640839),new tm(123,201.973788431),new tm(124,202.975137256),new tm(125,203.977705),new tm(126,204.97961)]),aC(VB(AD,1),KR,3,0,[new tm(95,174.991411),new tm(96,175.987413248),new tm(97,176.986336874),new tm(98,177.982476325),new tm(99,178.981783),new tm(100,179.978322),new tm(101,180.977806),new tm(102,181.97393546),new tm(103,182.974561),new tm(104,183.970705219),new tm(105,184.971983),new tm(106,185.969460021),new tm(107,186.969785),new tm(108,187.967511693),new tm(109,188.968733187),new tm(110,189.966958568),new tm(111,190.96706311),new tm(112,191.965921572),new tm(113,192.966644169),new tm(114,193.965381832),new tm(115,194.966638981),new tm(116,195.965814846),new tm(117,196.967195333),new tm(118,197.96675183),new tm(119,198.968262489),new tm(120,199.968308726),new tm(121,200.970285275),new tm(122,201.970625604),new tm(123,202.972857096),new tm(124,203.97347564),new tm(125,204.976056104),new tm(126,205.977498672),new tm(127,206.982577025),new tm(128,207.98594)]),aC(VB(AD,1),KR,3,0,[new tm(96,176.996881),new tm(97,177.994637),new tm(98,178.991466),new tm(99,179.990194),new tm(100,180.986904),new tm(101,181.98561),new tm(102,182.982697),new tm(103,183.98176),new tm(104,184.9791),new tm(105,185.977549881),new tm(106,186.97617),new tm(107,187.97592),new tm(108,188.974290451),new tm(109,189.974473379),new tm(110,190.972261952),new tm(111,191.972770785),new tm(112,192.970548),new tm(113,193.971053),new tm(114,194.96965),new tm(115,195.970515),new tm(116,196.9695362),new tm(117,197.970466294),new tm(118,198.969813837),new tm(119,199.970945394),new tm(120,200.97080377),new tm(121,201.972090569),new tm(122,202.972329088),new tm(123,203.973848646),new tm(124,204.97441227),new tm(125,205.976095321),new tm(126,206.977407908),new tm(127,207.982004653),new tm(128,208.985349125),new tm(129,209.990065574)]),aC(VB(AD,1),KR,3,0,[new tm(99,180.996714),new tm(100,181.992676101),new tm(101,182.99193),new tm(102,183.988198),new tm(103,184.98758),new tm(104,185.983485388),new tm(105,186.98403),new tm(106,187.979869108),new tm(107,188.98088),new tm(108,189.978180008),new tm(109,190.9782),new tm(110,191.975719811),new tm(111,192.97608),new tm(112,193.974648056),new tm(113,194.975920279),new tm(114,195.97271),new tm(115,196.97338),new tm(116,197.97198),new tm(117,198.972909384),new tm(118,199.97181556),new tm(119,200.972846589),new tm(120,201.972143786),new tm(121,202.973375491),new tm(122,203.973028761),new tm(123,204.974467112),new tm(124,205.974449002),new tm(125,206.975880605),new tm(126,207.97663585),new tm(127,208.981074801),new tm(128,209.984173129),new tm(129,210.988731474),new tm(130,211.991887495),new tm(131,212.9965),new tm(132,213.999798147)]),aC(VB(AD,1),KR,3,0,[new tm(102,184.997708),new tm(103,185.99648),new tm(104,186.993458),new tm(105,187.992173),new tm(106,188.989505),new tm(107,189.987520007),new tm(108,190.986053),new tm(109,191.985368),new tm(110,192.983662229),new tm(111,193.983430186),new tm(112,194.98112697),new tm(113,195.981236107),new tm(114,196.978934287),new tm(115,197.979024396),new tm(116,198.977576953),new tm(117,199.978141983),new tm(118,200.976970721),new tm(119,201.977674504),new tm(120,202.976868118),new tm(121,203.977805161),new tm(122,204.977374688),new tm(123,205.978482854),new tm(124,206.978455217),new tm(125,207.979726699),new tm(126,208.980383241),new tm(127,209.984104944),new tm(128,210.987258139),new tm(129,211.991271542),new tm(130,212.994374836),new tm(131,213.998698664),new tm(132,215.001832349),new tm(133,216.006199)]),aC(VB(AD,1),KR,3,0,[new tm(106,189.994293888),new tm(107,190.994653),new tm(108,191.99033039),new tm(109,192.991102),new tm(110,193.988284107),new tm(111,194.988045),new tm(112,195.985469432),new tm(113,196.985567),new tm(114,197.984024384),new tm(115,198.985044507),new tm(116,199.981735),new tm(117,200.982209),new tm(118,201.980704),new tm(119,202.981412863),new tm(120,203.980307113),new tm(121,204.981165396),new tm(122,205.980465241),new tm(123,206.981578228),new tm(124,207.981231059),new tm(125,208.982415788),new tm(126,209.982857396),new tm(127,210.986636869),new tm(128,211.988851755),new tm(129,212.992842522),new tm(130,213.995185949),new tm(131,214.999414609),new tm(132,216.001905198),new tm(133,217.006253),new tm(134,218.008965773)]),aC(VB(AD,1),KR,3,0,[new tm(108,193.000188),new tm(109,193.997973),new tm(110,194.996554),new tm(111,195.995702),new tm(112,196.993891293),new tm(113,197.99343368),new tm(114,198.991008569),new tm(115,199.990920883),new tm(116,200.988486908),new tm(117,201.988448629),new tm(118,202.986847216),new tm(119,203.987261559),new tm(120,204.986036352),new tm(121,205.986599242),new tm(122,206.985775861),new tm(123,207.986582508),new tm(124,208.986158678),new tm(125,209.987131308),new tm(126,210.987480806),new tm(127,211.990734657),new tm(128,212.99292115),new tm(129,213.996356412),new tm(130,214.998641245),new tm(131,216.002408839),new tm(132,217.004709619),new tm(133,218.008681458),new tm(134,219.011296478),new tm(135,220.015301),new tm(136,221.01814),new tm(137,222.02233),new tm(138,223.02534)]),aC(VB(AD,1),KR,3,0,[new tm(110,196.001117268),new tm(111,197.001661),new tm(112,197.998779978),new tm(113,198.998309),new tm(114,199.995634148),new tm(115,200.995535),new tm(116,201.993899382),new tm(117,202.994765192),new tm(118,203.991365),new tm(119,204.991668),new tm(120,205.99016),new tm(121,206.990726826),new tm(122,207.989631237),new tm(123,208.990376634),new tm(124,209.989679862),new tm(125,210.99058541),new tm(126,211.990688899),new tm(127,212.993868354),new tm(128,213.995346275),new tm(129,214.998729195),new tm(130,216.000258153),new tm(131,217.003914555),new tm(132,218.005586315),new tm(133,219.009474831),new tm(134,220.011384149),new tm(135,221.015455),new tm(136,222.017570472),new tm(137,223.02179),new tm(138,224.02409),new tm(139,225.02844),new tm(140,226.03089),new tm(141,227.035407),new tm(142,228.038084)]),aC(VB(AD,1),KR,3,0,[new tm(113,200.006499),new tm(114,201.00458692),new tm(115,202.00396885),new tm(116,203.001423829),new tm(117,204.001221209),new tm(118,204.998663961),new tm(119,205.998486886),new tm(120,206.996859385),new tm(121,207.997133849),new tm(122,208.995915421),new tm(123,209.996398327),new tm(124,210.995529332),new tm(125,211.996194988),new tm(126,212.996174845),new tm(127,213.99895474),new tm(128,215.000326029),new tm(129,216.003187873),new tm(130,217.004616452),new tm(131,218.007563326),new tm(132,219.009240843),new tm(133,220.012312978),new tm(134,221.014245654),new tm(135,222.017543957),new tm(136,223.019730712),new tm(137,224.023235513),new tm(138,225.025606914),new tm(139,226.029343423),new tm(140,227.031833167),new tm(141,228.034776087),new tm(142,229.038426),new tm(143,230.04251),new tm(144,231.045407),new tm(145,232.049654)]),aC(VB(AD,1),KR,3,0,[new tm(115,203.00921),new tm(116,204.006434513),new tm(117,205.006187),new tm(118,206.004463814),new tm(119,207.005176607),new tm(120,208.001776),new tm(121,209.001944),new tm(122,210.000446),new tm(123,211.000893996),new tm(124,211.999783492),new tm(125,213.000345847),new tm(126,214.000091141),new tm(127,215.002704195),new tm(128,216.003518402),new tm(129,217.00630601),new tm(130,218.007123948),new tm(131,219.010068787),new tm(132,220.011014669),new tm(133,221.013907762),new tm(134,222.01536182),new tm(135,223.01849714),new tm(136,224.020202004),new tm(137,225.023604463),new tm(138,226.025402555),new tm(139,227.029170677),new tm(140,228.031064101),new tm(141,229.034820309),new tm(142,230.037084774),new tm(143,231.04122),new tm(144,232.043693),new tm(145,233.047995),new tm(146,234.050547)]),aC(VB(AD,1),KR,3,0,[new tm(118,207.012469754),new tm(119,208.012112949),new tm(120,209.009568736),new tm(121,210.009256802),new tm(122,211.007648196),new tm(123,212.007811441),new tm(124,213.006573689),new tm(125,214.006893072),new tm(126,215.006450832),new tm(127,216.008721268),new tm(128,217.009332676),new tm(129,218.011625045),new tm(130,219.012404918),new tm(131,220.014752105),new tm(132,221.015575746),new tm(133,222.017828852),new tm(134,223.01912603),new tm(135,224.021708435),new tm(136,225.023220576),new tm(137,226.026089848),new tm(138,227.027746979),new tm(139,228.031014825),new tm(140,229.032930871),new tm(141,230.036025144),new tm(142,231.038551503),new tm(143,232.042022474),new tm(144,233.04455),new tm(145,234.04842),new tm(146,235.051102),new tm(147,236.055178)]),aC(VB(AD,1),KR,3,0,[new tm(120,210.015711883),new tm(121,211.016306912),new tm(122,212.012916),new tm(123,213.012962),new tm(124,214.011451),new tm(125,215.011726597),new tm(126,216.011050963),new tm(127,217.013066169),new tm(128,218.013267744),new tm(129,219.015521253),new tm(130,220.015733126),new tm(131,221.018171499),new tm(132,222.018454131),new tm(133,223.020795153),new tm(134,224.02145925),new tm(135,225.023941441),new tm(136,226.024890681),new tm(137,227.027698859),new tm(138,228.028731348),new tm(139,229.03175534),new tm(140,230.033126574),new tm(141,231.03629706),new tm(142,232.03805036),new tm(143,233.041576923),new tm(144,234.043595497),new tm(145,235.04750442),new tm(146,236.04971),new tm(147,237.053894),new tm(148,238.056243)]),aC(VB(AD,1),KR,3,0,[new tm(122,213.021183209),new tm(123,214.02073923),new tm(124,215.019097612),new tm(125,216.019109649),new tm(126,217.018288571),new tm(127,218.020007906),new tm(128,219.019880348),new tm(129,220.021876493),new tm(130,221.021863742),new tm(131,222.023726),new tm(132,223.023963748),new tm(133,224.025614854),new tm(134,225.026115172),new tm(135,226.02793275),new tm(136,227.028793151),new tm(137,228.031036942),new tm(138,229.032088601),new tm(139,230.034532562),new tm(140,231.035878898),new tm(141,232.03858172),new tm(142,233.040240235),new tm(143,234.043302325),new tm(144,235.045436759),new tm(145,236.048675176),new tm(146,237.05113943),new tm(147,238.054497046),new tm(148,239.05713),new tm(149,240.06098)]),aC(VB(AD,1),KR,3,0,[new tm(126,218.023487),new tm(127,219.024915423),new tm(128,220.024712),new tm(129,221.026351),new tm(130,222.02607),new tm(131,223.027722956),new tm(132,224.027590139),new tm(133,225.029384369),new tm(134,226.02933975),new tm(135,227.031140069),new tm(136,228.031366357),new tm(137,229.033496137),new tm(138,230.033927392),new tm(139,231.036289158),new tm(140,232.03714628),new tm(141,233.039628196),new tm(142,234.040945606),new tm(143,235.043923062),new tm(144,236.045561897),new tm(145,237.048723955),new tm(146,238.050782583),new tm(147,239.054287777),new tm(148,240.056585734),new tm(149,241.06033),new tm(150,242.062925)]),aC(VB(AD,1),KR,3,0,[new tm(132,225.033899689),new tm(133,226.035129),new tm(134,227.034958261),new tm(135,228.03618),new tm(136,229.036246866),new tm(137,230.037812591),new tm(138,231.038233161),new tm(139,232.040099),new tm(140,233.04073235),new tm(141,234.042888556),new tm(142,235.044055876),new tm(143,236.046559724),new tm(144,237.048167253),new tm(145,238.050940464),new tm(146,239.052931399),new tm(147,240.056168828),new tm(148,241.058246266),new tm(149,242.061635),new tm(150,243.064273),new tm(151,244.06785)]),aC(VB(AD,1),KR,3,0,[new tm(134,228.038727686),new tm(135,229.040138934),new tm(136,230.039645603),new tm(137,231.041258),new tm(138,232.041179445),new tm(139,233.04298757),new tm(140,234.043304681),new tm(141,235.0452815),new tm(142,236.046048088),new tm(143,237.048403774),new tm(144,238.0495534),new tm(145,239.052156519),new tm(146,240.05380746),new tm(147,241.056845291),new tm(148,242.058736847),new tm(149,243.061997013),new tm(150,244.06419765),new tm(151,245.067738657),new tm(152,246.070198429),new tm(153,247.07407)]),aC(VB(AD,1),KR,3,0,[new tm(136,231.04556),new tm(137,232.04659),new tm(138,233.046472),new tm(139,234.047794),new tm(140,235.048029),new tm(141,236.049569),new tm(142,237.049970748),new tm(143,238.051977839),new tm(144,239.053018481),new tm(145,240.055287826),new tm(146,241.056822944),new tm(147,242.059543039),new tm(148,243.061372686),new tm(149,244.064279429),new tm(150,245.066445398),new tm(151,246.069768438),new tm(152,247.072086),new tm(153,248.075745),new tm(154,249.07848)]),aC(VB(AD,1),KR,3,0,[new tm(137,233.0508),new tm(138,234.05024),new tm(139,235.051591),new tm(140,236.051405),new tm(141,237.052891),new tm(142,238.053016298),new tm(143,239.054951),new tm(144,240.055519046),new tm(145,241.057646736),new tm(146,242.058829326),new tm(147,243.061382249),new tm(148,244.062746349),new tm(149,245.065485586),new tm(150,246.067217551),new tm(151,247.070346811),new tm(152,248.072342247),new tm(153,249.075947062),new tm(154,250.078350687),new tm(155,251.082277873),new tm(156,252.08487)]),aC(VB(AD,1),KR,3,0,[new tm(138,235.05658),new tm(139,236.05733),new tm(140,237.057127),new tm(141,238.058266),new tm(142,239.058362),new tm(143,240.059749),new tm(144,241.060223),new tm(145,242.06205),new tm(146,243.06300157),new tm(147,244.065167882),new tm(148,245.066355386),new tm(149,246.068666836),new tm(150,247.070298533),new tm(151,248.07308),new tm(152,249.074979937),new tm(153,250.078310529),new tm(154,251.08075344),new tm(155,252.084303),new tm(156,253.08688),new tm(157,254.0906)]),aC(VB(AD,1),KR,3,0,[new tm(139,237.06207),new tm(140,238.06141),new tm(141,239.062579),new tm(142,240.062295),new tm(143,241.063716),new tm(144,242.063688713),new tm(145,243.065421),new tm(146,244.06599039),new tm(147,245.068039),new tm(148,246.068798807),new tm(149,247.070992043),new tm(150,248.07217808),new tm(151,249.074846818),new tm(152,250.076399951),new tm(153,251.079580056),new tm(154,252.081619582),new tm(155,253.085126791),new tm(156,254.087316198),new tm(157,255.091039),new tm(158,256.09344)]),aC(VB(AD,1),KR,3,0,[new tm(141,240.06892),new tm(142,241.068662),new tm(143,242.069699),new tm(144,243.069631),new tm(145,244.070969),new tm(146,245.071317),new tm(147,246.072965),new tm(148,247.07365),new tm(149,248.075458),new tm(150,249.076405),new tm(151,250.078654),new tm(152,251.079983592),new tm(153,252.082972247),new tm(154,253.084817974),new tm(155,254.088016026),new tm(156,255.090266386),new tm(157,256.093592),new tm(158,257.095979)]),aC(VB(AD,1),KR,3,0,[new tm(142,242.07343),new tm(143,243.07451),new tm(144,244.074077),new tm(145,245.075375),new tm(146,246.075281634),new tm(147,247.076819),new tm(148,248.077184411),new tm(149,249.079024),new tm(150,250.079514759),new tm(151,251.081566467),new tm(152,252.082460071),new tm(153,253.085176259),new tm(154,254.086847795),new tm(155,255.089955466),new tm(156,256.091766522),new tm(157,257.095098635),new tm(158,258.097069),new tm(159,259.100588)]),aC(VB(AD,1),KR,3,0,[new tm(144,245.081017),new tm(145,246.081933),new tm(146,247.081804),new tm(147,248.082909),new tm(148,249.083002),new tm(149,250.084488),new tm(150,251.084919),new tm(151,252.08663),new tm(152,253.08728),new tm(153,254.089725),new tm(154,255.091075196),new tm(155,256.094052757),new tm(156,257.095534643),new tm(157,258.098425321),new tm(158,259.100503),new tm(159,260.103645)]),aC(VB(AD,1),KR,3,0,[new tm(147,249.087823),new tm(148,250.087493),new tm(149,251.08896),new tm(150,252.088965909),new tm(151,253.090649),new tm(152,254.090948746),new tm(153,255.093232449),new tm(154,256.094275879),new tm(155,257.096852778),new tm(156,258.0982),new tm(157,259.101024),new tm(158,260.102636),new tm(159,261.105743),new tm(160,262.10752)]),aC(VB(AD,1),KR,3,0,[new tm(148,251.09436),new tm(149,252.09533),new tm(150,253.095258),new tm(151,254.096587),new tm(152,255.096769),new tm(153,256.098763),new tm(154,257.099606),new tm(155,258.101883),new tm(156,259.10299),new tm(157,260.105572),new tm(158,261.106941),new tm(159,262.109692),new tm(160,263.111394)]),aC(VB(AD,1),KR,3,0,[new tm(149,253.100679),new tm(150,254.100166),new tm(151,255.101492),new tm(152,256.101179573),new tm(153,257.103072),new tm(154,258.103568),new tm(155,259.105628),new tm(156,260.106434),new tm(157,261.108752),new tm(158,262.109918),new tm(159,263.11254),new tm(160,264.113978)]),aC(VB(AD,1),KR,3,0,[new tm(150,255.107398),new tm(151,256.10811),new tm(152,257.107858),new tm(153,258.109438),new tm(154,259.109721),new tm(155,260.111427),new tm(156,261.112106),new tm(157,262.114153),new tm(158,263.115078),new tm(159,264.117473),new tm(160,265.118659)]),aC(VB(AD,1),KR,3,0,[new tm(152,258.113151),new tm(153,259.114652),new tm(154,260.114435447),new tm(155,261.116199),new tm(156,262.116477),new tm(157,263.118313),new tm(158,264.118924),new tm(159,265.121066),new tm(160,266.121928)]),aC(VB(AD,1),KR,3,0,[new tm(153,260.121803),new tm(154,261.1218),new tm(155,262.123009),new tm(156,263.123146),new tm(157,264.12473),new tm(158,265.125198),new tm(159,266.127009),new tm(160,267.12774)]),aC(VB(AD,1),KR,3,0,[new tm(155,263.12871),new tm(156,264.128408258),new tm(157,265.130001),new tm(158,266.130042),new tm(159,267.131774),new tm(160,268.132156),new tm(161,269.134114)]),aC(VB(AD,1),KR,3,0,[new tm(156,265.136567),new tm(157,266.13794),new tm(158,267.137526),new tm(159,268.138816),new tm(160,269.139106),new tm(161,270.140723),new tm(162,271.141229)])]);}var fQ='object',gQ='anonymous',hQ='fnStack',iQ={4:1,10:1,5:1,8:1},jQ='Unknown',kQ='\n',lQ='boolean',mQ='number',nQ='string',oQ=2147483647,pQ='__java$exception',qQ='For input string: "',rQ='null',sQ=-2147483648,tQ=524288,uQ='__noinit__',vQ={4:1,14:1},wQ={4:1,12:1,14:1},xQ=4096,yQ=16384,zQ=65536,AQ=65535,BQ='fromIndex: ',CQ=', toIndex: ',DQ={7:1,4:1,5:1},EQ=16777215,FQ=0.30000001192092896,GQ={15:1,4:1,5:1},HQ={11:1,4:1,5:1},IQ=536870912,JQ=2.617993878,KQ=3.665191429,LQ=6.283185307179586,MQ=3.141592653589793,NQ=1.5707963267948966,OQ=1024,PQ=234881024,QQ=100663296,RQ=201326592,SQ=114688,TQ=4063232,UQ=2097152,VQ=393216,WQ=29360128,XQ=268435456,YQ=2048,ZQ=-1.5707963267948966,$Q=16320,_Q='" ',aR='stroke-width:',bR=786432,cR=262144,dR=1.0471975511965976,eR=0.5235987755982988,fR={4:1,5:1,8:1},gR={4:1,5:1},hR='unsupported atomicNo:',iR=8192,jR={l:0,m:0,h:64},kR='Bit already set!',lR={l:0,m:0,h:128},mR={23:1,4:1,10:1,5:1,8:1},nR={4:1,5:1,18:1,8:1},oR=-16777216,pR={13:1,4:1,5:1},qR=131072,rR=-65536,sR=-1.0471975511965976,tR=-3.141592653589793,uR=2.0943951023931953,vR={4:1,10:1,27:1,5:1,18:1,8:1,28:1},wR='??',xR={99:1,4:1,10:1,5:1,8:1},yR=-268435456,zR=1572864,AR=65011712,BR=3072,CR=126976,DR=1.7976931348623157E308,ER=67108864,FR=134217728,GR=16777216,HR=-66584577,IR=0.7853981633974483,JR=3.061592653589793,KR={9:1,4:1,5:1,8:1},LR='ATOMS',MR='M  END',NR='$$$$',OR=4194303,PR=239060990,QR='class="event" ',RR='Assignment of aromatic double bonds failed',SR='Members of ESR groups must only be stereo centers with known configuration.',TR='Ambiguous configuration at stereo center because of 2 parallel bonds',UR=-0.5235987755982988,VR=-0.7853981633974483,WR=277296187394,XR=277296187395,YR=280517412866,ZR=280517412867,$R=280520558594,_R=280520558595,aS=280520561666,bS=284812380162,cS=284812380163,dS=284815525890,eS=284815528962,fS=284819720194,gS=284819727362,hS=414735140866,iS=414735140867,jS=414738286594,kS=414738286595,lS=414738289666,mS=414742480898,nS=414742480899,oS=414742488067,pS=414869358594,qS=414869358595,rS=414869361666,sS=414869489666,tS=417956366338,uS=417959512066,vS=552174094338,wS=552177240066,xS=552177243138,yS=552181434370,zS=552181441538,AS=552308312066,BS=552308315138,CS=552308319234,DS=552308319240,ES=552308443138,FS=555395319810,GS=555395319816,HS=555398465538,IS=555398468610,JS=555398468615,KS=555402659842,LS=555402662914,MS=555529537538,NS=555529537544,OS=555529540610,PS=555529544706,QS=555529668610,RS=555532683266,SS=555532686338,TS=559693432834,US=559693435906,VS=559697630210,WS={l:2361346,m:590400,h:16},XS={l:2361346,m:1376832,h:16},YS={l:2361346,m:1377600,h:16},ZS={l:1312770,m:1377601,h:16},$S={l:1315842,m:1377601,h:16},_S={l:2361346,m:2425408,h:16},aT={l:2361346,m:2426176,h:16},bT={l:1312770,m:2426177,h:16},cT={l:1315842,m:2426177,h:16},dT={l:2361346,m:2427200,h:16},eT={l:1312770,m:2427201,h:16},fT={l:1312770,m:2427202,h:16},gT={l:1315842,m:2427202,h:16},hT={l:1319938,m:2427202,h:16},iT={l:2361346,m:590400,h:24},jT={l:2361346,m:591168,h:24},kT={l:1312770,m:591169,h:24},lT={l:1315842,m:591169,h:24},mT={l:1319938,m:592194,h:24},nT={l:2361346,m:623168,h:24},oT={l:1312770,m:623169,h:24},pT={l:2364418,m:623200,h:24},qT={l:2361351,m:1377600,h:32},rT={l:1312775,m:1377601,h:32},sT={l:1315847,m:1377601,h:32},tT={l:1312775,m:1378625,h:32},uT={l:1315847,m:1378625,h:32},vT={l:1315847,m:1378626,h:32},wT={l:1315847,m:1409601,h:32},xT={l:2361352,m:1443136,h:32},yT={l:1312776,m:1443137,h:32},zT={l:1315848,m:1443137,h:32},AT={l:2361352,m:1443168,h:32},BT={l:2364424,m:1443168,h:32},CT={l:1312775,m:2426177,h:32},DT=0.6262000203132629,ET=-1.3825000524520874,FT=-1.4915000200271606,GT=0.33169999718666077,HT=0.3540000021457672,IT=0.38179999589920044,JT=-0.6019999980926514,KT=-0.7379999756813049,LT='Atom-types are 64-bit numbers describing atoms and their near surrounding.',MT='Recognized atom types and their contributions are:',NT='Druglikeness predictor not properly initialized.',OT=3.009999990463257,PT=-0.1809999942779541,QT=-0.17000000178813934,RT=-0.2029999941587448,ST='Over- or under-specified stereo feature or more than one racemic type bond',TT='undefined',UT=0.08726646502812703,VT=1048575,WT='Too many percent/per mille characters in pattern "',XT=4194304,YT=17592186044416,ZT=-17592186044416,$T='CSS1Compat',_T='safari',aU='Possible problem with your *.gwt.xml module file.\nThe compile time user.agent value (safari) does not match the runtime user.agent value (',bU=').\n',cU='Expect more errors.',dU=1.52587890625E-5,eU={4:1,10:1,5:1,18:1,8:1},fU={30:1,51:1},gU={37:1,33:1,36:1},hU={37:1,33:1,74:1,36:1,93:1},iU={37:1,33:1,36:1,62:1},jU=15525485,kU=5.9604644775390625E-8,lU={4:1,30:1,53:1,42:1},mU='Invalid UTF8 sequence';var _,BG,wG,bG=-1;CG();DG(1,null,{},nc);_.ab=function oc(a){return this===a;};_.bb=function qc(){return this.Jb;};_.cb=function sc(){return YP(this);};_.db=function uc(){return mc(this);};_.equals=function(a){return this.ab(a);};_.hashCode=function(){return this.cb();};_.toString=function(){return this.db();};var hB;DG(186,1,{});DG(138,186,{},oB);_.pb=function pB(a){var b={},j;var c=[];a[hQ]=c;var d=arguments.callee.caller;while(d){var e=(iB(),d.name||(d.name=lB(d.toString())));c.push(e);var f=':'+e;var g=b[f];if(g){var h,i;for(h=0,i=g.length;h<i;h++){if(g[h]===d){return;}}}(g||(b[f]=[])).push(d);d=d.caller;}};_.qb=function qB(a){var b,c,d,e;d=(iB(),a&&a[hQ]?a[hQ]:[]);c=d.length;e=ZB(gF,iQ,43,c,0,1);for(b=0;b<c;b++){e[b]=new uJ(d[b],null,-1);}return e;};DG(187,186,{});_.pb=function sB(a){};_.rb=function tB(a,b,c,d){return new uJ(b,a+'@'+d,c<0?-1:c);};_.qb=function uB(a){var b,c,d,e,f,g,h;e=(iB(),h=a.backingJsObject,h&&h.stack?h.stack.split(kQ):[]);f=ZB(gF,iQ,43,0,0,1);b=0;d=e.length;if(d==0){return f;}g=rB(this,e[0]);DJ(g.d,gQ)||(f[b++]=g);for(c=1;c<d;c++){f[b++]=rB(this,e[c]);}return f;};DG(139,187,{},vB);_.rb=function wB(a,b,c,d){return new uJ(b,a,-1);};var KC,LC,MC;KC={4:1,135:1,30:1};var LH;DG(96,1,{},cI);_.sb=function dI(a){var b;b=new cI();b.e=4;a>1?b.c=jI(this,a-1):b.c=this;return b;};_.tb=function iI(){aI(this);return this.b;};_.ub=function kI(){return bI(this);};_.vb=function mI(){aI(this);return this.i;};_.wb=function oI(){return(this.e&4)!=0;};_.xb=function pI(){return(this.e&1)!=0;};_.db=function sI(){return((this.e&2)!=0?'interface ':(this.e&1)!=0?'':'class ')+(aI(this),this.k);};_.e=0;_.g=0;var _H=1;DG(70,1,{4:1,70:1});var uI;LC={4:1,30:1,136:1,70:1};DG(14,1,vQ);_.lb=function zA(a){return new $wnd.Error(a);};_.mb=function BA(){return this.f;};_.nb=function CA(){var a,b,c;c=this.f==null?null:this.f.replace(new $wnd.RegExp(kQ,'g'),' ');b=(a=bI(this.Jb),c==null?a:a+': '+c);xA(this,AA(this.lb(b)));jB(this);};_.db=function EA(){return yA(this,this.mb());};_.backingJsObject=uQ;_.j=true;DG(12,14,wQ,FA);DG(31,70,{4:1,30:1,31:1,70:1},QI);_.fb=function SI(a){return RI(this.a,a.a);};_.ab=function TI(a){return PC(a,31)&&a.a==this.a;};_.cb=function UI(){return this.a;};_.db=function XI(){return''+this.a;};_.a=0;DG(26,12,wQ,HA);DG(61,26,wQ,IA);DG(80,61,wQ,pJ,qJ,rJ);_.lb=function sJ(a){return new $wnd.TypeError(a);};MC={4:1,94:1,30:1,2:1};var XP=0;var ZP,$P=0,_P;var eF=fI(1);var nE=fI(0);var uE=fI(186);var rE=fI(138);var tE=fI(187);var sE=fI(139);var QE=fI(135);var RE=fI(96);var dF=fI(70);var SE=fI(136);var lF=fI(14);var VE=fI(12);var ZE=fI(31);var fF=fI(26);var $E=fI(61);var bF=fI(80);var kF=fI(2);DG(140,1,{});_.v=0;_.w=0;_.B=0;_.F=0;_.I=0;_.J=0;_.L=0;_.M=0;_.O=0;_.P=0;_.Q=0;_.R=0;_.S=0;var wc,xc,yc,zc,Ac,Bc,Cc,Dc,Ec,Fc,Gc,Hc;var dD=fI(140);DG(68,1,{68:1},Ed);_.a=0;_.b=0;_.c=0;var bD=fI(68);DG(35,1,{},Fd);_.a=0;_.b=0;_.c=0;_.d=0;var cD=fI(35);DG(84,1,{},Pd,Qd);_.a=0;_.b=0;_.e=0;var eD=fI(84);var Yd,Zd;DG(47,1,{},qf,rf);_.p=0;_.r=0;_.A=false;_.C=0;_.F=false;_.G=false;_.I=0;_.K=0;_.N=0;_.Q=false;_._=false;var pD=fI(47);DG(144,1,{},tf);_.eb=function uf(a,b){return sf(a,b);};_.ab=function vf(a){return this===a;};var gD=fI(144);DG(88,1,{88:1},wf);_.b=0;_.c=0;_.d=0;var fD=fI(88);DG(145,1,{},yf);_.eb=function zf(a,b){return xf(a,b);};_.ab=function Af(a){return this===a;};var iD=fI(145);DG(89,1,{89:1},Bf);_.a=0;_.b=0;_.c=0;var hD=fI(89);DG(67,1,{67:1,30:1},Gf);_.fb=function Hf(a){return Ef(this,a);};_.a=0;_.b=0;_.c=0;var jD=fI(67);DG(114,1,{114:1},If);var kD=fI(114);DG(153,1,{},Zf);var nD=fI(153);DG(156,1,{},ag);_.eb=function bg(a,b){return _f(a,b);};_.ab=function cg(a){return this===a;};var lD=fI(156);DG(154,1,{},ng);_.a=0;_.b=0;_.f=0;_.g=0;_.i=0;var mD=fI(154);DG(155,1,{},pg);_.eb=function qg(a,b){return og(a,b);};_.ab=function rg(a){return this===a;};var oD=fI(155);DG(86,1,{},Tg);_.b=0;_.d=0;_.g=0;var rD=fI(86);DG(29,1,{29:1},gh,hh);_.f=0;_.i=0;_.j=0;_.k=false;_.n=0;_.o=0;var qD=fI(29);DG(65,1,{65:1,4:1,30:1},mh);_.fb=function nh(a){return ih(this,a);};_.ab=function oh(a){var b;if(a==null||!PC(a,65))return false;b=a;return $wnd.Math.abs(b.a-this.a)+$wnd.Math.abs(b.b-this.b)+$wnd.Math.abs(b.c-this.c)<1.0E-6;};_.db=function ph(){var a;a=new mK('0.00');return'['+kK(a,this.a)+', '+kK(a,this.b)+', '+kK(a,this.c)+']';};_.a=0;_.b=0;_.c=0;var sD=fI(65);DG(87,1,{},wh,xh);_.db=function yh(){return'DepictorTransformation Offset: '+this.a+','+this.b+' Scaling: '+this.c;};_.a=0;_.b=0;_.c=0;var tD=fI(87);DG(49,1,{49:1},zh);_.a=0;_.b=0;_.c=0;_.d=0;var uD=fI(49);DG(59,1,{},Bh);_.a=0;_.b=0;_.c=false;_.d=0;_.f=false;_.g=0;_.i=false;_.j=0;var vD=fI(59);DG(56,1,{56:1,4:1});_.o=0;_.p=0;_.G=0;_.I=false;_.J=false;_.K=0;_.L=0;_.P=false;_.Q=0;_.R=0;_.S=0;var Ch,Dh,Eh,Fh=24;var CD=fI(56);DG(64,56,{64:1,56:1,4:1});_.gb=function em(a){Lk(this,a);};_.d=0;_.e=0;var wD=fI(64);DG(34,1,{},om);_.b=false;_.c=0;_.d=0;_.e=0;var xD=fI(34);DG(20,1,{20:1},pm,qm);_.a=0;_.b=0;var yD=fI(20);DG(106,1,{},sm);var zD=fI(106);DG(3,1,{3:1},tm);_.a=0;_.b=0;var AD=fI(3);var um;DG(126,1,{});_.a=0;_.d=0;var xm,ym,zm;var BD=fI(126);DG(82,1,{},Hm);var DD=fI(82);DG(101,1,{},_m);var ED=fI(101);DG(127,1,{});var FD=fI(127);DG(85,1,{},vn);var GD=fI(85);DG(57,1,{},Vn,Wn);_.b=0;_.j=0;_.k=0;_.n=false;_.s=0;_.t=0;_.F=false;_.G=0;var JD=fI(57);DG(103,1,{103:1},Xn);_.a=0;_.b=0;_.c=0;_.d=0;var HD=fI(103);DG(102,1,{},go);var Yn,Zn;var ID=fI(102);DG(98,140,{},zo);_.db=function Ao(){return wo(this);};_.f=0;_.i=0;_.j=0;_.k=0;var mo=0;var KD=fI(98);DG(129,1,{},Eo);_.a=0;_.d=0;_.i=0;var LD=fI(129);DG(128,1,{},No);var ND=fI(128);DG(109,1,{109:1},Ro);_.a=0;_.b=false;_.c=0;_.d=0;_.e=false;_.g=0;var MD=fI(109);DG(58,1,{58:1,4:1},Uo);_.hb=function Vo(a){return To(this,a);};var OD=fI(58);DG(24,64,{64:1,56:1,24:1,4:1},kp,lp,mp);_.gb=function np(a){Xo(this,a);};_.a=false;var PD=fI(24);DG(133,58,{58:1,4:1},pp);_.hb=function qp(a){var b;b=To(this,a);if(b==-1)return-1;return BM(this.a,b).a;};var QD=fI(133);var yp=0,zp=0,Ap=0,Bp=0,Cp=0,Dp=0,Ep=0,Fp=0,Gp=0,Hp=0;DG(172,1,{});var RD=fI(172);DG(132,172,{},Vp);_.d=0;var Mp;var SD=fI(132);DG(100,1,{},aq);var Wp,Xp,Yp;var TD=fI(100);DG(123,1,{},hq);var cq='sNy@LDeVjj@XTKU|TH\t5.89\nKA|P@QvRRJjjj@LFaLJnyC@\t5.86\nHiFH@H{IIEUJjj@FGPfES]rA@\t5.82\nJoB@@BUJssoPvxPTA@@@FFpaUpriv\x7FDsP\t5.61\nHeTH@@RV[TYZ`@@AaUMpssHP`\t5.19\nHmtH@@RVYWeVhH@@FCESCJ|rDx\t5.18\n.1\nsNpOAxRPTai@rEdGHCh\x7Fjj@P\t-10.77',dq,eq=false;var UD=fI(123);DG(131,1,{},kq);var WD=fI(131);DG(110,1,{110:1},lq);_.b=0;var VD=fI(110);DG(22,1,{22:1},mq);_.b=0;var YD=fI(22);DG(66,1,{},oq);var XD=fI(66);var pq,qq;var yq,zq;DG(134,1,{},$q);var Dq,Eq,Fq='eOHBNZ`pge@\ngCa@@dmHFFbwH@\ngFp@DiTujhCBbKWdH\ngJPXHlQxQ{TAaeb\ngCi`HEdfZ@pRp\nfJ@FD\ngJPHAbIJuPFADVyB\ngC`@Die@ppfyD\ngJQ@@dkU@XFKGd@\ngJP@DizhC@qX|`@\ngJQ@@dru@XI[dH\ngGQ@@eMuTA`Xl^R`\ngJQ@@eKU@XYX|d@\ngGQ@@djuTAaQEcrT\neMABHXaIhH\ngCa@@dkPFBbyL\ngGQ@@drmTAaekrD\ngJQ@@eKS@XZK\\a@\ngCh@@dmPFFDwH`\neMACD\\QIhH\ngJU@DPdju`P\ngGX`HDdwMLA@\ngGY@JDivjpH\ngCi@LDek@`\ndeTH@@RYWVf`@j@CC`pjxYyB@\ndid@p@bBbFbDfYoa`b@@LJ@fx^QP\ndeTH@@RYe]aZZj`cJ\ngFq@@drfmU@X[F|b@\ndeTH@@rJJIHdsSUL@P\ndeTL@@QdfygFV``@@pjfxYyB@\ndaDH@@RVU[f@@@LBcB[bUp\ndifD@HADfyWaZjj@H\ngGQ@@dkLtAad[rP\ngGQ@@dlltHkCXwd`\ndifH@JAJ[gxZB@@CBdJf{dB\ndidHPBBHFHRYgVzB@`@phLKayE@\ndeVD`FFPbDfUnkh@a`@`\ndiD@@LdbJTZjh@pzDpjX^Qp\neFDBcA`d\ndidL@@KdiuVDjjj@H\ndaDH@@RYe[hB@@LJpj[nP@\ndiT@@DiYXfifjjh@`\ngF|@ABeKZsU@P\ngNq`@fdvkSHf\\EcqJ`\ndiDL@@PTfU]jZP`f\ngJP`@TeZhCCKGbD@\ngNp`@dfujj@plVMx`p\ngGQ`@jdjmTA`l^JT\ngNp`@dfUZi@pVOEV\ngNq@@djkUPFEVM_I@\ngNq@@djmUPFCbqky@\ndiDD@@QIeuZfhHJpkBiny@`\ndaxD@@QIgUjj@LLpfxe\\\ndaxD@@QImeijBLlBSJ[d\\\ngNp`@deVZj@pNM_I@\ngNp`@dfvZj@pNM_H`\ndaxD@@QIeejj@LBpj[d\\\ngNq@@djmSPf\\GEcWr@\ngJPdE`DRPcyXCrD\ngGQ@@dsmLIkCKGdh\ngNq@@dsmUPFEDVOEV\ndaxH@@RYuji`bgA`SBknIC@\ngGQLJIARAdDfzj`LInJT\ngGQLLIAREdDfvj`LINP`\ndiDHhLBPRPjPzPFPR[[jjj@H\ngNq`@jdjsUPFAqkyH\ngGP`ATiVj`LMEnHt\ndaxH@@RYWjj`CC`aLInyE@\nfHc`A@\ndiF@@@rRHiJjj`CChSBiny@`\nday@@@{IEHkUU@FC@fES\\c`\ngGXDL@aABS\\uPFFAoH`\ngNxLL@aAABDfVZj@`\ngGPB@DHHpQPaIUZdB\ngGPD@DXHRfZf@ppUxlP\ngGX@@dsuTAaQEcrT\ngJY@BDeZlCBSbB@\ngGY@LDeVj`LCD[qC@\ngNx@@eRmUPFEbu_HP\ngNx`HDdvkSPfX\neMhDRZCAKd`\ngGQ`@ZdruTA`qEcP\ndiDDpLH@bOA@aIkUZjh@`\ndaDH@@RVU[f@@@LJcBinQp\ngFp`@dfTujXCAZ|a@\ngOx`FDdrikTtA@\ngJP@DkfhC@bH|f@\ngC`D@DSpRnhB\ngNq`@bdvkUPFFV_IP\ngGPP@cTfyi`ODj\ngJP`@TeZhC@qX|`@\ngCaHL@aIZPLDIrH\ngNp`@dfUZj@pvMyF\ndaxH@@RUVZj`CAdpj[nQ@\ngJQHLHaIVj@`\ngNphMQDIK[UTA@\ngJQ@@dkSBJpHbwHP\ngNqBLIAREdGHIMmUTA@\ndaxL@@RdfuVjh@pILInQp\ngGQ`@jdvmTA`l^IT\ngNp`@|dTQjj@plVMyF\ngJQ@@eOU@XZH|f@\ngGT@ATiVj`LJHm^P`\ngGT@Ade[j`LHnHt\ngCe@H`dkPGbV@\ngGUHLZHaIUjdB\ndazH@LAIV^jj@LFBDInyE@\ngJP@DiVhPVFbwDB\neMDARVCAd\ngNphBqDILsTrA@\nfHghA@\ndeTD@@EIYe^efjZ@H\ndeTD@@QIgeQej@@@LJrfF^Qp\ndcLD@@UIUe]FVX@J@aKCdrfx]yD@\ndmtD@@QImYVUZX@@Hrp{B[ae^P`\ngGP`@deUjPLEcqR`\ndeTH@@RUYTYi`@@aK@XSJ[agd@\ndmtD@@QIee^UZ``@@pXjXYWd\\\ngOp@DjWkjj`LJEc^JL\ndmtHPBBHfHRYeUXXHHh@H\ngOp@DjWkZj`LCD[qY`\ngJQHBJqIVz@`\ngCa`@ldkPD\ndmtL@@RTeYW^Eh@J@CBbinWdP\ndid@@DjU^nBBD@LFaLiaxa\\\ndeTD@@qI[eQej@@@LNSJ[agdP\ndcLD@@QInUuxV`@j@CCdJfzUyF`\ndmwD@HePQInUwaZ@B`@`\ndaE@@@YIeZn`B@@pILinHG@\ndmvL@HAFR[f^FV``H@H\ndeVD`NFPbDfUvih@I`@`\ndmvD@DATfYUQUjjj`CAlJfx^QH\ndeVD`BxPbDfYYZXHF@@`\ndmt@H@bAdIdEdDfUvjZ@Bj@C@`pjyB@\ndcNH@DAIee^eVhHB@CCbine]yB@\ndmvH@JAJUuTjjjjh@pFDpfxYyB`\ndmvH`HX@cIEDdTljjVj`B\ngJQ@@dkU@XDSGdP\ndiDB@@SaRUYVjj@LJrfx^IA@\ngGQ@@dkUTA`Xm^HD\ngGP@H\x7FUPD\ngGP@DkYj`LJEc^R@\ndaxD@@QIgUjfBJlLpfxe\\\ngGQ@@eNuLIkCKWba@\ndctH@@RgYujfjZBLX\ndiDD@@QImiZjh@piLiny@`\ngJP`@dfvdCBGbV@\ngNq`@fdrkThD\ngJT@@defhCCSGd@\ndeTL`HS@BLddlRPrm@@@FEYSCOHx\ndmLD@@QIee\\jeVhHB@CBbXYWd\\\ndet@@Dje^ifzjjj`B\ngOp@DiUMZj`LFHlWrT\ngC`H@DIMTAa`mrP\ngNq`@fdr{UHFBqyJ\ngJPXHlPDQztAxlP\ngJQ@@dmS@XZX|`@\ngJQ@@eMU@XYX|d@\ngNphBpDISkURA`QEnR`\ngGPhLQDIKUU@P\ngJP@DkVhC@bK\\a@\ngJQ@@djsBJpTwDb\ngNq@@drmUPFCDVMyF\ngGP`@TfVj`LJHl^R`\ngNq`@fdrkUHD\ngGY@LDenj`LJHl^R`\ngGPhCQbILmU@P\ngOx@@drm]UTA`plZ~R@\ngGP@DjZj`LCEkq@`\ngNq@@dkMMPb\\CD[qA`\ngJQ@@dju@XJGbV@\ngNq@@ds]UPFCDVMyF\ngNx@@eLsUP`lLb~KT\nfHdXA@\ndidH@@RUe^Ejjh@pyLJfx^Q@\ndidH@@Rf~hRjjh@`\ndmv@`EBHrJJIHin`HFpHy`\ndeVD@HADfyeFV`H@@piJ[iy@`\ngOp@DjWkB@@LBmWqP`\ndeVD@HADfyUxV`@`@piJXYyA@\ndeVD@HADfyWxV`@`@piJ[iyA@\ndmv@@@rQQJEJUjh@@@pdHPfxYWdT\ngNqhHl@cIIJeiaCP\ngNpP@jtfvZf@pfxdp\ngGP`@dfUjpLH^R`\ngNp@DkUzj@pjqkyH\ngJQHHOAJuj@prqy@\ngJQHHOBOZ`H\ngJQ@@dlu@XJGbV@\ngNq@@eJuUPFEbu_HP\ndeU@@@aJyenF`HH@H\ngGY`HEdf]j`LEkrD\ngJX@@dksBIpROI`\ndedD@@QImUVjj`C@Tpj[ad\neO`BNZ`pYy@\ndaDH@@RYe[jfj@LFBLJnyC@\ngJQ@@dsT`XDQ[dH\ndid@@DjU^nBBD@LNaLJfGd\\\ndidL@@RdeVWaZjj@LJpfx^Q@\ndaF@@@Rfu[j@@@LFABinyF@\ndedD@@QImiVji`B\ngGP@DjZj`LCL[qA@\ndiDL@@RdeVyjf`B\ndmtH@@ReY}Jjjjj@H\ndcLL@@Sdf^YV]Z```@LJqae^Qh\ngGT@@dfuj`LBEcrT\ndif@PACDJHRYgvzB@`@`\ndcMH@DTLbbbRRHjuUUT@X]aTwCoHd\ndkmD`LND@HrRPjZIE]VhBB@B\ndidHPBBHFHRYgVzB@`@pHLx^HU@\ndmvH@DCHhhhTiUjjj`CAl[ae^HU@\ndknH@DAIfYuUMjjjj@H\ndcND@DCdf^YV]Z```@LFQae]yF@\ngF|@AbeJf`@@P\ndifD@HADfyWaZ@@@LBRnGbDp\nfHapA@\nf`i@`@@VRYfYU]`eNMyh@`AB@@H\ndk\x7F@@@cLdRfbTQragSfhJ@T@H\ndmLH@@rJJIQEneX@@@@C@`xYxVTr`\ndk]H@BDLbbbbRQZ]NB@P@@B\nfoA@R@HHqH@QddebRbrPeV\\m@D@@@A@\ndaFH`JHHaIf]n``@@`\ndmv@HBBHFPfPVPRYUzih@Jh@H\ndeTHPIBPzPRYeea`Ha@CCLX^QP\ndeVHPIHHchaIf^VFBBH@LBCNGbEH\nek`PJ@@@GNimlbbRfbbebrTRLrThXTlBbrjZVjjjjjfijh@CAICBc`RQSPrqspIHjhX{dFP\ndmvL`EaL@HrRRqIXYV`@`@`\nfewAP@@LtT^QQQRUJQYSQZXpgCNCeNVfjjijjjZj@B\nfoQ`@@@YIEDeDTdqWAF]UUAAE@@XB`cENRUkaFFlx\nfde``@C@heMrklk|dYpXtDDUUTT@D\ndk^H@EAJ[UVVySh@JjH@pyBhYT~If`\nfoA@@@LdbRdVeRiJs``j@@@@LEPPRgAJtZwDR\\\nf`qA@@@YEDeeHhTjL]z@Bj@@@@`\ndcnH@EAIfV^XYv@@B@@H\ndmtL`HS@BLddlRTFUh@H@H\ndeT@H@bBbAbIbDfYu[hHB`@`\ndcLHpJBPRPrPrJPsQIKmTp@@A`WBn|bP\ngGP@Di]ja@xTQF|f@\ndcLB@@Q]R[e[neh@a`@piBj]yG@\ndeU@@@eJYW~F``H@LFCBinxbR\ndeUD@HdDR[eWaZ@B@C@diaxaR\ndeVH@DAIgeQej@@@LJSJX^It`\ndcNL@HAFR[fUqUhHH`@`\ndg^L`EaC@HrRPsIYJCJt@EL@FCTwBm|``\ndaFHpBxHa@bhcHheBSTuL@P\ndg^L@D@mRY^UueVj@Bh@H\ngOp@DjWkB@@LMc^ZI`\ndidD@@yIfVXXBH@C@`[axfT\ndaDH@@RYWih@H@LJSBknP@\ndeTD@@YIfUqehD@@LJJnF^Q`\ngFtHHPDIRnMKPFBEyD\ndcnDPJa@BBBLdabRRS\\nkSP@@A@\ngNxhMD@cIHUEj`LCcWr@\ndclL@@pTign_JWZjB@`@`\ndcnD`HI`BDfYoVnWZfX@@@`\ndidH@@RVY^Ejjh@pyLJfx^Q@\ndedB@@PiRUi]jjT@`\ndmTD@@QIeUyjZZBBh\neMFINbMP|`\ndo}H`BMPbLbbbRfebXXHHfj@B\ndo}H@DhDfUfUWWZhJB`@H\ngGT`EaTf]j`LDkqX`\ndayL@DpFyIeUjj@H\ngKP@H~Jj`LEcqQ@\ndid@@LdbLTifjj`B\ndaDD@@QInXjZjh@`\neMJDBDfPpce@\ndaEH@FxDiebjiV`CCDJnHw@\ndk^H@FAJY}e\\kSie`@`HipyLhYW^HQ`\ndaFH@LAIVUnZjh@pHj[nIB@\ndieH@DHDfvWaZjj@H\ndg~H@LCHeEEDdcJg[UPTA@A@\ndifD@HADfyWaZ@@@LLPjyc\\H\ndg^L@HAER[e[[xV`@jh@H\ndeVH@LAIVYQejjj@LNPj[agdP\ndeVH@BAIV[Qejjj@H\ndmtHPEBHzHrJJISEa`HbP@`\ngGQHLHaI[ihCCBWdH\ndmuD@HXDR[fUEV```@LBRiWfMp`\ndg}@@@mJYeU|]Tz@@@H@B\nfhy`B@N@BLdTTTRRVqirUmNh@`BBh@@pt@cARUhugAyCp\ngKP@LdabjhC@bH|f@\ndmtH`ABHRYW[ih@Jh@LNALJaWb]H\ndg|@@DjWmijXYB@jjb@H\ndmtD@@gHhhhdVEjjj`C@TpfxYWdT\ndeT@@DjU_k``RPHjpFDpj[ayD`\ndeTJ@@qaeJYyzzjjj@LIaLJfxY@\ndev`@@rfeJY{ZxYBBJD@`\ndefH@LAIVYfjj`B\nf`aQB@BFTBHrJJIJZUJLEAADuT@A@\nfHdHA@\ndidHPFBHJHRYf~FBH@@`\ngCe@H`dkPFDwH`\ndcND@DCTfVutYZ`@d@LJJfx^Hb`\ngGT`IPdfuj`LDmqF`\nda{D@Hi`QImejj@H\ndkm@`DDHaIe]YZZ@Bjh@H\ndcL@@LdbRbjeBDEEP@XCBES\\n|SI@\ndmvHPBTH`XaIfUmi``hP@pXLIagd\\\ngGYHCaDIK]M@XHSbM@\ndcND@DCdefV]]Z`b@@H\ndmtD@@QIVYdUZ@b@@`\neFPBca@\ndiDJ`HSJDOCIIJdfjj@LJSBx^Pp\nf`aqR@AF}AFJZAxYIITdhhmkNZjjfV`@pxcANJm{dF@\ndid@p@bFbAbDfUfn`BH@LLAaybYp\ndeVL@HAFR[f\\YZB@@C@diixgB\ndif@`ABHRYevz@``@phLJny@`\ndcNH`BdHaIfUyXXHHZ@B\ndo}H`FMPbLbbRfRaRkh@bfj@B\ndayH@DhDfWVjh@pKB[nHe@\ngNx`BDdvkUPFFu_Db\ndigD@Dq`yIeUifff`B\ndk^D@D@\\bbbbRQImMj`XB@C@jWSxgR\nfgA@@@DjYU_VByHu`@@@@@@H\nfoA`@@@ILkjrmFV]@AL@@@ar`\ndcLD@@SHhmEDcJmPDD@FA\\L{qLX\ndeUD@BdDR[YTYZjj`B\nfoA@b@HHAxIRlrjzkF]U@@@@@A@\ndif@@@RifvFjjh@pzDJfx^Q`\ndaD@`@bDfUZZ@B@C@`qnxbD\ndcLD`BTHaIfUVXXHHf@B\ndidH@@rJJIEn`HH@LJ@jX^Q`\nfbu@`@@NrJJJIPjFKQLDFLADUA@T@@P\ndidH`ACDRYWZZ@B`@pXDpjGd\\\ndg}@@@mIe]e^ftx@H@H@B\ndmvH`FdHaIe[zn`BI`@pXDrf{dJ\ndcNH@EAJYYeGZBHh@B\ndeTL`HRPBLdabRwBl@D@FEES]OH@\ngOp`Adigujj`LCEWrD\ndmtH@@RYeUEV`P`BJlNpjxYWdL\ndeTL`HS@BLddlRPrm@@@FETwCOHp\ndeV@@@rQQQHcMAP@@XBBXUMp|bP\ndmvD@H@de[eYVZ`@@bJ\ndaE@@@yIe^f`@`@pKBknHB@\ndcLH@@RfUW~f``b`@pEBinE]xfR\ndo|H@@rJJIQPrEn`HJjh@H\ndg^D@MADfVU}iUjB@j@B\ngGY@LDenjPLBHcqZ`\nf`qh@@@XIQfRJJKZJEJgG^ejj`@`@H\ndcnH@LAIYe_x^fjjjj@H\ndg^B`LaAl@cIIBheEeikTBPH@XUgSi|Re@\ndg^B`LaEt@cIIKEDhcIkPPLP@P\ndo}B`LfDEpBLddJbbbJNujAbA@B\ndieH`LE`BDiU_Bjij@LBAJ{bPH\ndk~@@@RfYU_JGUN`@@B@@pDB[aeSyF`\ndmtB`HSE@HrRRqIXYV`@`@pKFy^IP`\ndmtL@@jTef_^E``J@C@biexd\\\ndcl@@DjYU_egX@@@@@pxjxYW^XfL\nf`q@`@@HRYyV{TRg^Z`B@@@@H\nf`qA@@@ILk\\joSagPA@PP`@D\ndmt@H@bAdIdEdDfU^jZ@Bj@C@`pfyB@\ndcLHHBBHfPVPvPRYg^fzB@j`@`\ndcn@@@Re]eRi]jj@B@CBXPjxYW^HF@\nfbuQB@BFTBHrJJIJZUIQILDgLDDSU@@@@P\ndcMH@DhDfufU]Zjjj@H\ndmtD@@[HhhhdYUhJ@@CAlJae^Hw@\ndmtD@@qIYyVUZh@@@pYFxYWbPP\ndeT@@LdbRTm\\DDT@FGIeMpsrH@\ndmtHpEBHJPFHRYgW[hHBd@H\ndid@@DiUWajjj@LAaLJfx^P@\nfde@P@@BLEIfYfY\x7FyiWgQZ@HJ@b@@H\ndieH@JxLbTTQkfej`CAFGbPP\ndaD@P@bNbDfUzZ@B@C@`pnxbT\ndmMH@DhDfUmZZU`@@@@H\ngNqdEb@b^FQRHmU@P\nfoAab@GPQ`@QddebbTVLmFlm@@@P@A@\nfH`pA@\ndg\\B@@SSrJISISPbkT@Pt@D\ndeV@@@RiU\\Yjjj`C@XSBkagdL\ndmLD@@IJ]YVDeZj@B@B\ndmLH@@RgYVaAfj@B`@`\ndid@@DjWxjZjZ@H\ndaxD@@QIUijj@LLqnxdT\ndmtD@@QIee^UZ``@@pZfxYWdD\ngC`dE`DSpRZXCsBX\ndeV@@@RfV\\YhH`@CCdJfxYy@@\ndcMH@FxLdTReRQUTkUT@P\ndcNL@HAErJZJIECKPDE@A@\ndefD`BpPbDemgijj@H\ndo~H@DCHeEEEEBmmjBbH`@`\ndo~D@D@|bTTTTRqvvhJHb@B\ndknH@ICHhhhdUFF@bJh@H\ndeT@@LdbRbmBDED@FGPfTwSrB@\nfoAab@NPQ`@QddebbrRTeV\\mA@@@@A@\ndeVD`Aa@BLdabRgRl@P@D\nfoAPB@NJ@DYHhhheEcJqgT@PP@@@P\neMAALhbN``\ndidH@@RYWZZ@B`@pXDpjGd\\\nfde`B@N@BLdTTTRRVqeNRmiu@D@PUP@A@\ndk^L@IANRY[f~]tvjjjj@LApjXUt~It`\nf`iPB@N^@DISLro\x7FSdcZmPA@@W@@D\ndcLD@@QIVYVFV`HJ@CAfxYW^IE@\ndklD@@QIgfUiUj@`h@H\ndk\\D@@QIVZVVfSZjjjh@`\ndk\\D@@sIEDdXdj[Sjj`@`@`\ndeL`@@JfRiUfnXVfjjjh@`\ndif@PBBPFPRYgvzB@`@psJ[dD\ndmvHPBdIAYAIfUua``a`@`\ndieHPJD`bFbDfYoaiZi@H\ngJX`LDdvu@XI[dH\ndcnH@BAJ]fuaEvjh@H@H\nfHcxA@\neFAADdRLD\nfHdxA@\neMJD|Df`pYy@\ngJPhHaxIRuPFBqy@\ngJ]@EbDfVhCAH|f@',Gq='fHe`A@\ngFq`@ldrfmU@XR|a@\ngCl@@ldsPFFBp\ngC`HADIKTAaaMrH\ngNxHLHaIYjj`H\ngNxHLHaIVjj`H\neFDBcAaWH@\neMDARVCBnR\ngGP`@TfYi`LI^S@\ngJQ@@dkT`XDKGd`\ngGQ@@dkUTAaXl[r@\ngJQ@@djsBJprqyH\neMHAIdLF^P\ngJQ@@dls@XKGd`\ngGP`ATiVj`LCEkrD\ngJP`AdejhC@qX|`@\ngGP`ATiVj`LCEcrT\ngJQ@@dsU@XDKGd`\neMJDbDfP`\ngC`DADZHRVhB\ngC`DAb[DRVhB\ngCa@@dkHFBbyL\ngCa@@dkPFBfyD\neMABHYAIhH\ngJPH@DIRuPFABqyH\ngCa@@dkHFBVyH\ngGQ`@jdjmTAal[rP\ngGY@LDeVj`LJHm^P`\ngJXHLHaIVj@pHbOI`\ngOx@@drm\\@@A`plZp\ngGXHLHaIUjhCBbHwdp\ngCh@@doHFDwH`\ngJY@BDeZlCAQ\\`@\nfHgdA@\ngJQ@@dsT`XDKGd`\ngOp`@dfUMZj`LMc^Q`\ndaD@@DjUZxHD@CAhSBinQp\ngOp`@dfUkZj`LMj~P`\ndeTH@@RVYWahBA@CC``j[ayD`\ndid@p@bBbFbDfYoa`b@@H\ndaDH@@RYg[ffj@LBrf{dD\ngFq@@drfmM@X[F|b@\ngGQhHj@cIHTmPFFqoH@\ngJQhHl@cIHUhCBGd@\ndaE@@@yIe^f`@`@piLJny@@\neMhDRUB\ndie@@@aJvUxZ`@@CChPj[ay@@\ndifH@DAIf_Ifjj`CBlJf{dB\ndifD@BADfyWaZ@@@LBQnGdT\ndeVD@FADfygFV``@@pjfxYyB@\ndeVD@AADfyVzV`B@@piBkayD`\ngF|@AbeJfuU@P\ndetD@@QInYvDYZjjh@`\ndedH@@RUUUfjhHRpELJfxYyD@\ngJQ@@dkU@XZX|PP\ngNp@DiWjj@p\\VM_H@\ngNq@@djmUPFEfM_DD\ngJQ@@dkT`XZK\\PH\ngJQ@@dkU@XDQGdp\neFA@HoBJD\ngJP@DjvhCCKGd`\ngGQLLIAREdDfvj`H\ndiDHhABPRPjPZPzPRZyjjj@H\ndiDB@@SaRY]fifBBX\ngNxHF@aJUzZqDxXH|Tp\ndmv@@@Rf~UeZj@@@LEBDpfxYT\ndaxD@@QIUYjj@LBrf{bPP\ngNq@AdTbMUPFEBq_IP\ngGP@DiUjaAXFKGbE@\ndaE@@@YIeZn`B@@piLiny@@\ngCa@@dsPFBV@\ngOp@DjWkfZ`LKEc^Q`\ngOp@DjWkZjPLKEc^Q`\ndaD@`@bDeeVz`@@CA`cBinQp\neMABHPaIhLDnR\ndaDH@@RVU[f@@@LJ`j[nQ`\ngC``AdeZ@pTWI`\ngFp`AdiTvjhCCQWdH\neMbDBDfp`\ndaEH@JXDiWRjjj`CBhSB[d\\\ndaFH@NAIe^f`@`@piLJny@@\ndaF@@@RYe[hB@@LJCBknPp\ndaF@`NBHRYUih@H@LJCB[nP`\ngC`@Die@ptVy@\neFA@HoBLD\ngGQ`@jdvmTAaecrT\ndaxB@@RfRYYZf`B\ngJQ@@eKS@XJKbq@\ngCaHLLQIZ`LDEqS@\ndaxD@@QIeUjj@LBpj{dL\ndaxD@@QIUUjj@LJpj[nQ@\ngNp`@df]Zj@pvMyF\ngJP`@dfzhCCA[ba@\ngGQ@@eKuTA`Uc^R@\ngCaHLHaIZ`LLHnS@\ngGPhMQDIK]U@XTQX|e@\ngNp@DiUZjDC`qEc^Q`\ngNpH@DIRoUTA`qEj~P`\ngNp`AdiWjj@pJM_I@\ngGP`@TeZj`LKEc^P@\ngGQ@@eJuTA`Xm^P`\neFJBhHp^I@\ngCah@mJAIj`H\ngCahHlOBOTAaAsQX\ngJT`H`TeVdB\ngNx@@eJmThFCbqky@\ndif@@@RUe^Fh@@@pDHPj[a@\ngNy`LDtfuZj@pNM_H`\ngNx`LFdjmUPFCDQkyL\ngOx@@drm]UTAaqEcV\ngGT@ATivj`LKEc^P@\ndiDDHJXIAICi@YAIkfjfh@`\ndcLL@@STfVyVUZ`HD@H\ndaDH@@RVU[fZj@LBJf{bQ`\ndmtD@@QIUYVUZh@@@p{B[ae^Q@\ndeTD@@QIgeQejjj@LFpj[ayD`\ndmtD@@QIgeTYZjjh@p{B[ae^QP\ngOp@DjWkZj`LCL[qI`\ngNp`@df^Zj@pvkyB\ndeT@@DjWvifjjh@pFDpjXYyG@\ndidD@@qJY~rjjZ`B\ngGQ@@dmltA`h^KT\ndcL@@DjU_ZnZjij@H\ndcL@@LdbRbjUBDEEP@XSBXUMt{rE@\ndid@@DjWZfZjj@LFcB[ayB@\ndid@p@bBbAbDfYun``H@H\ndcMB@HDDWTfyV{iZ@HX@H\ngC`dEaDPHRZTB\ndmtD`NDHcHhheDVfBAb@CB`rfWdR\ngNxhGD@cIHTefqMP\ndaFH`HX@aJYWJjeh@pHDJnQp\ndax@@DkUfjh@pZDpf{bAP\ngNp`@df^Zj@pV_DZ\ndedD@@qJ[^ZjZ`bf\ndaxD@@QImUjj@LJSBknPp\ndiD@@LdRbJZjhHBpxHpj[nH``\ndaxD@@QInejj@LBRf{dD\nfH`TA@\ngJQ`@bdjt`P\ndaD@`@bDfUjZ@B@CB`SJ{dL\ndkm@`ATHaIe[ujZ@BfhBAh\ndif@`BBHRYgfzB@`@`\ndetD@@eIYe~DYZjjh@`\ndaFH@DAIYUnZjX@pkBinyD@\ngC`HAVIMTAaaMrH\ngNq@@djuUPFCDqkyD\ngJQ@@dju@XZX|b@\ngCa@@dmHFFDwH`\ngJQ@@dsMBRppVyB\ngCahHlOBNtA`anQ@\ngGQDJH`qBSKMHdX\neMIDbKpRYB\ngNx@@eJ}UPFCDVkyB\ngJXhEbLQIZf@`\ndnD@@DiYrbYjj`CA`aLinPP\ngF}@EbDfTuiXB\ngJP`@TfZhCCQ[bA@\ngJPhLQDIKUPD\ndiDB@@RnRYufjf@LDp^PP\ndmTD@@SHheHjfjjh@pXj[agdJ\ngNpTHjpDDHrREQZTB\ngGP@LdbMU@XTQZ|a@\ndiDH@@rJJQUjj`CBlJf{dB\ngJQDDH`qBS]LHj\ngNplJqHJPtadTaeTpGdX\ngJYHLPDIStpblDEqP`\ndiFDpJXPdDdLdLbdJRjfdHI`\ngGX@@eKMTHGCKWba@\neFHBJFE@\nfHfpAa@\nfH`XA@\ngOu@DPdrykURA`l~Q@\nfHgHA@\ngCa@@dmHFBVxa@\ndmtD@@QIgYVUZh@@@p[FxYWdD\ndmtH@@RgfueZj@@@LASJ[ae^Q@\ndid@p@bFbAbDfUfn`BH@H\ndidH@@RYUZZ@B`@phLKayB@\ndmuL@DpIUIfVTfZjjX@`\ndeVD@HADfyVxV`@`@piBkiy@`\ndeVL@HAIR[e_aZ@B@CBdJngdL\ndmvL`NaL@HrRRqQZUV``@@`\ndco@`LK`BLdTTRRITntpTA@Pe@\ngGP@Djuj`LLm^JD\ngKP@Di\\Vj@pHfOH`\ngNp`@dfVZf@pQ_IP\ndaxD`Fx@aJUzjf@LJBDsnPP\ngGQ@@djuTAaQL[rH\ndaxH@@RUUjj`CC`cBinyB@\ndiDH@@RYujjj@LAALJfx^Q@\ngJX@@dku@XIGdp\ngJX`DBdju@XI[ba@\ndazH@LAIUjjj@LFBDinyF@\ngNq`@fdjkUHFBqxiP\ndaDD@@aJyUnh@@@`\ndid@@LdbbQxXF@@CAdrfx^PP\ndid@@Dj{WaZjf@LNaLJf{dB\ndaDH@@RYUifjj@LJpj[nP`\ndaDD@@YIeZn`B@@piLiny@@\ndidD@@iJ[gxZB@@CBdJnGdL\ndaFD@FADfyVyjj`CBdpj[d\\\ndewH@HP`RY[TjFZd@H@H\ndiD@@DiuejjP`GChSB[ax`T\ndmtDPNDHaXaIfVUi``X`@`\ndmt@@DjU_jxHHj@CBXPj[ae^Q`\ndmt@@DjU_ZxHDj@bg@XSBkagbLh\ndkl@@LdbRdSRjP`jJ`@pfDJfxUOdZ\ngOp@DiWMZj`LKEb~HT\ndmtH@@RYfWXXBHh@LF@fFUxe\\\ndax@@LddUeUT@XMBXS]rJ@\ngOq@@eJqmUTA`Xl~Ht\ndaD@P@qBbDfYvzB@@CB`pj[d\\\ndmuL@HDDWHihdh^eh@b@B\ndcNL@HALRYymUujh@@@`\ngNyhMDpDYIBdmTA`\\Z~P@\ngNx@@dlmUPFEbq_DJ\ndknL@CAErJIQIIF]Z``b@B\ndeVD@DCdeeY[fjjh@`\ndid@@DjUfaBB`@LNaLinGdD\ndifL@DCaRY]bijih@`\ndifH@DAIf_Ifjj`C@lJnxcB\ndcMB@HXDeTfyed]ZBA`@H\nfHcdAa@\nfHchA@\neF`BNFE@\ndeTD`HP@cIHXdepk@A@A`ULL|PB@\ndeTL`HS@BLddlRPrm@@@FAXwSqJD\ndcnH@NCHhheEBtkl@D@@@P\nfoA@R@HHqX@QddebRfR`iF\\m@@A@@A@\ndieH@BxDfYUa``P@LBCJ[b\\H\ndeTB@@pYRf[^njjj`CBXSBinFP\ndaE@@@aJyUnX@@@pkBinyD@\nfduA@@@ILsLjm{AJ\\XOhm@@@@A@@A@\nfhy`@@@ISLjm{btjw`t@@@@@@@P\ndmvH@DAIf{VUZh@@@p[J[agdJ\nfnkA`@@N[dTRTtTTlVRbUFJlFNZmKUUUUUUT@F@TXipTeZMYw`iKWbLm@\nfig@P@@NZOHhdihhiXleDbjLUXL\\uZVjjjjjij`@`\ndeTL`HS@BLddlRPrm@@@FAXwCqJd\ndet@@DjYUX^dHbH`CAdJfx^Id`\nffsA`@@LudTTTeRdVTtLIps`ySeijjjZjjj@B\nf`ia@@E@RfuUe]gEF]z`@jjhh@LAHIPTmFmsFG^@\nf`ia@@M@RfuU[UgEZ]z`@jjhh@LAHIPTeFmsFG^@\ndmLH@@RYegXYV@@@@@pXJXYWbYp\ndmtD@@QIee^UZ``@@pXfxYWdT\ndmtL@@QTfyeQehBA@C@jXYxb\\\ndcLB@@RUR[fVQuhHF@@`\ndklH@@rJIJQQNfZjji`B\ndmt@@DjU_ZxHHj@CBXSBine^PH\ndaDH@@RYWih@H@LB@j{bI`\ngNq@@dsKSPFFu_HP\ndeTL@@QdfygFV``@@piJXYyG@\ndmv@@@RfYWEZB``@LIALJfxYyB`\ndmuL@HTDYInYtYZB@`@`\ndeVH@IAJYW~F``H@LJPj[nId`\ngFuHC\\@aJYMif@`\ndifH@HAIYexV@`@C@biny@`\ndmvD@D@dfWeYUj`@@CCdpfxYyB`\nfoA``@H@PdwJ{J|EYsP@P@@@D\ndif@pDBHjHFHrJIQEn`HH@H\ndclL@@{TivY~DeZhHB`@`\ngFp`ATiTvjhCAH|Tp\nfgA`B@N@BDifYWz\\d[Uj@H@B@@LI@Hs`eZM[dI@\ngC`HAxIKTAahmr@\ndg|@`@|DjYmUyO[j@@@@@LM@j[ae]N~Q`\ndeTD@@YIfUqehH@@LFJfxYyF@\ndeTL@@jTef_xVB@`@phj[iyD@\ndeTL@@RdfVUFVjZd@pHfF^Qp\ndeT`@@biRnY\x7FaXHB@B\ndif@PJ@H`HRnY~F``@@pXBinGdP\ndk^@@@RfYU\\]Tz@@@@@LECBinFUOdZ\ndieDPJZD@HHHrRFIYnVfi@H\ndmND@DCdfVUrjUZjZi@LFrfFUyB@\ngOp@DjWkjj`LFEcWrP\ngOx@@eJqmUTA`xlZ~P@\ngOx@@eLm]UTA`xlZ~P@\ndcnH@DAIYegzUujBHH@LNInFUwdP\nday@@@aJVYjjB@h\ndedB@@PYR[UYjjX@pILZ^HU@\ndkLF@@RUttfyenZjif@LFSBiSyG@\ndmTL@@QdfUivijdHJ`\ngGY@LDeej`LBl[rP\ngGY@BDeUj`LBl[rP\ngGY@LDeUj`LBl[rP\ndkmH`NMPbLbbbTNfXXBHf`@`\ngGY@LDeUj`LLc^JX\ndmv@@@rRJIIFUjB`@@pELJfzUy@`\ngOx@@drm]KTA`Pl^Jl\ndaF@@@RYUijVj@H\ndif@@@RYWZZejh@`\nfoAa@@D@RUfV]qZlyhH@@@@CBRJ\\DkQkNxbD`\nfoAa@@D@RYYeUuVLyj@@@@@CARJLEIVcV]rG@\nfj}a@@D@rJIQQIQKQEYSkQrUj@@`@@@@@`\ndie@@@iJYWxYB@@bK@lInybTH\ndifH@NAIe]ih@I@C@dpjx`B\ndcn@@@rRIHqIER{UT@D@FDpaUprn|c@\ndig@@@aDiyWaV@@@LJPj[nPH\ndg~DPFvpbEBLbbRfbRM\\JpAESI@D\ndeVD@JADeUeFVjjh@p[B[agdP\ndmuD@LXDRYueeVjjj@H\ndmtD`BTHaIfUma``bP@`\ndaE@@@{IHdbUUUT@XMBXUMrN@\ndifH@LAIVUxVjj`CAdJfx^Q@\ndeVH@IAJYW~F``H@LFSBinyD`\ndifD@HADfyWaZ@@@LLRayaMp\ndeTD@@QIgeQej@@@LLsayeMH\ndeT@@DjUghP`h`@pYL[agfPU@\ndeVH@LAIUeQejjj@LASBinF^Q@\ndmuH@DHDfvYYUjjj`B\ngOx@@eLvmUTA`xlZ~P@\ngOx@@eRimUTAaXcWrX\nfle``@C@Pdrrj\x7FLlmQRuAAEUTP@P\ndiFH`JpHaI[kijh@`\ndeVL@BAIR[YTYZjj`B\ndg|H`ABHRYW[ficn@BjjH@`\nfoAqB@EZ\\HDPdrnvrtYYt@EURd@FD`HpR`iFmrG@\ndk_H@FdprJJISPkatzjjiZ@H\nfoAab@GPQ`@QddebRfRpiFlm@@@P@A@\ngNy@LDeVZj@pJM_EL\nf`aQ@@DT@drllsNoMTDQT@@P\ngJX@@eKU@XYX|P`\ndid@@LddRL[jjj`CChSBkayC@\ndk]H@DdDefueFUujjjj`B\ndg|@@DjU_eZx{BAH@@BJlARne]N~EFDp\ndmvD@HADfyeQehBA@C@dXYye\\h\ngBX@@eLUTAahmr@\ndmLH@@RYiYKnVjjjh@pFLinFUy@@\ndg_@`DGPbDfUueZZ@Bjj@B\ndeVH`BdHaIfUvFBBD@H\ndaFH@LAIVUnZjh@piBinyD@\ndg^H@LAJUyfUSjhBH`@`\ngGPHAbIKUU@P\ndmtH@@rIQQQWiXBH`@pjfxUyD@\ndcOD@Ds`wHheELUPmMUT`FCTwBn|bP\ndeTH@@rJJIHmtAAP@XTALL|sJx\ngOu@E`drm[SRAalWrT\ndo}H`AMPbDfUo[Vf`@jZh@H\ndg]H`AMPbDfY_[Vz@`ij@B\ndkmH`NMPbLbbRrbaih@Jf`@pxDpjg^PX\nf`q@@@BbHRDRHQHaIXkf`bJbjb@CARB\\EKQkN}rB`\ndcLL@@G\\dTRRbOKPPTP@XRfES\\LkrM@\ndklL@@Ptfym]eVj@BP@pyLinF^Qh\ndg}B@HTDf|bfbbTThfnmA@tI@D\nfoAP@@@HR[ieUuVLyjhJBH@B\nfoAab@KPQ`@QddebRfR`iF\\m@@A@@A@\ndg~H@FCIEDhcLdLg]AP@@@A@\ndmvH@DCHhhhdYUjjj`CClJfxYyB`\nfhiA`@@Hddjrm|jIW`mPAD@@@A@\ndk~@@@RfYU_JGUN`XJJH@pDB[aeSyF`\nf`i``@E@PdwJvvoAJt{sP@TuQP@P\ndcNL@FAMRUUeUujjjh@`\ndidH@@rJJIEn`HH@LBJfGbA`\nfhyA`@@BMdTTTTTTVoMJ|xKPAAPA@@F@aRTZsoAyCP\nflu@P@@BLgHhhhhhhml^ZUytV`BB`HX@@`\ndeTL`HS@|LddlRPru@@@FETwCODS@\ndk^LPLaC@HTHrRPqQYKiWUjfYj@LJXUt~QP\ndcNH@BAIUfYgVjjj`B\ndmt@@DjU^jxHHj@CChPjFUxf\\\ndg\\D@@QInUukaZ@Bj`@`\ndknJ@HAIT|bfbbTUGV``X`@`\ndidH@@Rge^Fh@@@pZDJf{dB\nf`qAA@A@bOQBSJ{\\ktYYt@EP@P@A@\ndmtD`ATHaIe]nf`@jP@pXDpjGb]H\ndcND@NADfUyU]Zj@@@H\ndcn@@@RigVRX]fBBb@CClkae]yB@\ngCa@@duPFADV@\ndidL@@cDkkWajjj@H\ndg|H@@RVYYwySn``@@@@pTHIne]N~PH\ndmtLPNePbABDfUujZ@Be@B\ndaDh@DqnAIeZfZZd@`\ndmM@PBx@c@aJYg\\jeZdHB@B\neFAAx`bLD\ndmv@@@rQQJEJUjh@@@pdHJfxYWdH\ndmLH@@rJJIQEneX@@@@CA`rnF^Hr`\nfoAa@@D@rJJJJHqYQkNZj`bB@@`\ndklH@@RUYffSYjZj`aJ\ndeTH@@rJIJFTt@EP@XLFTpsqDx\ndo|H@@rJJIQRFIn`HJjh@H\ndmvH@DAIgfVUZ`H@@pYFxYWdX\ndeTD`HP@cIHXhdLk@P@A`UMp|PI@\ndcLD@@QIeeUgVhHH@CAbinF^Hn`\ndknD@FADfye_EV``bP@pyJ[aW^II`\ndo~D@FADfye_TUZBBIhBMX\ndg^H@LAIYVUW[jiBB@Hi`\ndkmD@NLJrQSQQITUhHbI@B\ndiW@@@cDi[WBxYjeX@`\ndev@@@rRIIHus]UUUPA`JXUMpsrH@\ndifD@B@TfYun``H@LJCJx^PP\ndo~@@@RV^UviUj`@j`@`\ndaDH`BBHRYg[hH@@LJCJ[nPP\ndif@`ABHrJJIEn`HH@LJCBinPH\ngJXLBIARFdDfjhB\ndcLHPBBHzHrJJKQFLLDDU@A@\ndcLHPBCD{DrJJKQFLLDDU@A@\ndklJ`HSNL@cIIKEDdYuZBBH@H\nfhiAb@B^BBHRYe^unXHshHH@HX@@`\ndg\\L`AWPbDfUv{ZZ@Bij@B\ndcn@@@Re[mRY]jj@B@CBXSBxYW^HC@\nfde`@@@IRmrkNyFZ\\FMUPP@@P@A@\ndk_@`LI`BDigvUrmNfV@B@bj\nfoApA@EZ\\BHeDILrk|kNV]@PMTi@ar`\ndaE@@@yJeVnjjh@pZDpj{dL\ndmvH@LAIUYVUZh@@@pELInFUyD@\ndmL@@DjYUVGi@@@`@LNpjxYWdH\ngJX@@djsBIptQxa`\ndmV@@@RUgVYjf`aJ\ndeTD@@qI[eQfj@@@LNSJ[agbA@\ndazD@J@dfWjjh@pjLInx`H\ndg^D@EAdfYewiuhJBJ@B\ndidD@@QIVUxV`@@CCFx^YAT\ngJX@@eST`XZK\\a@\ndmvD@HALbfbbQFV`HH@H\ndk]H@BxLbbbRaRX]NBA@@@B\ndeu@@@gIHhikWLMUTu@D\ndid@P@bNbDfYYa`H`@LJBfx^Q`\ndmtHPIBHVHRYfUXXBHX@H\ndeUD@FxJRVYmnYjZ`B\nf`q@`@@LrQQJIQKHbL{uAPA@@@A`iANBdkUg^y@p\ndk\\@@DjWmkiadHBjh`B\ndk]H@DXDfYYwz]MjhHB@CClhYWSyC@\ndg~HpFlI@i@YCHhhihdUtz\\DQEJPA@\ndcLL@@S\\bbTrTHru@AH@XDUMODEP\ndg\\`@@SFRYueUNvjjjj@H\ndk\\@`LhDjU^ukmLHH@@@CAXPjxYWSyF@\ndk]H`AdpqDfUmUiev@Bfd`B\ngFp@DjTujhCCKWba@\ndkLD@@SHdiDbeFjff`B\ngGP`Ademj`LBl[rP\nfhyH@@@X\\EJYnUWoEEVMyj@@@b`@B\nfhyH@@@X|EJYnUWwEEVLyj@@@J`@B\ndk^H@FAJY}e\\kSie`@`Hj`\ndieH@FxLdTReJjeZ`B\ndeVD@HADfyeFV`H@@pqNgfTp`\ndeU@@@qJYejxBHh@LJJfF^Qp\ndmuH@DXDfUgjZ@Bj@B\ndmvH@ACHhhdcFz@`j@CA`Jfe^Ig@\ndmvH@AAIe]Zf`@j`@pXBkiWbI`\ndcn@@@RfumVy]d@@@@CClIae]yB@\nfoAP@@@TRfUVu~RlzP```@ABC@\ndknDpItpdDdLdLbdLRTtEZh@a`@`\nf`qa`@H@PqInYWmQJ]yhHH@@@@`\nfbuac@HjSdD`bPqHYEHXdleDeeAR|Fj@B`@f@@H\ndeTD`HP@cIHXhdLk@P@A`UMt|PA@\nfdy`b@LPP@HrRPjJJIIDf|xJu@A@A@@D\nf`qa@@M@rQSIYQIPlxZu@AUUQ@A@\ndig@@@`Tke]nX@H@LBpfGbQP\ndaF@`B@HRYg[hH@@LBpj[bAp\ndid@@Ldbbq[`bB@C@lJaxbL\ndidD`HPOAJvUxVjj`CAlJfx^P@\ndmtB`HSBCprRSFIJUZh@@@pinFUxfD\nfHbXA@\ndaFH@BAIf]n``@@pKBknHC@\nfjsQ@@DB@dsLsKjvldIUhaJgKU@p@@@D@@D\ndnD@@DiYrbYjj`CB`aLkbDp\ndmLD@@eIYfUayVjjZPB\ndaDH@@RVU[f@@@LJcB[nQP\ngFp`AdiTvjhCAF|TP\neMBBHRYCAKd`\ndieH@HPDeYWaZ@@@H\ngOu@E`dsu[UTA`TZ~S@\ngOx`DFdrikTlA`e^S@\ndaD@@DjWZXHB@CBdpf{dP\ndmL@@DjYeVdUBHhb@LNpjXYWd\\\ndiTH@@RVYV{ajjf`B\ndmN@@@RYVuiiV@@@@@pxLkae^P`\ndml@@LdfbTJifzUZjjj`B\nfHcDA@\nfHgPAa@\nfHe@Aa@\neMBBHRYCAGe@\ngC`DAbZHRVhB\nfHbxA@\ngCa@@dtpFBVy@\ngC`HADIMTAa`mrP\ngJQHDHaIjj@`\ngC`HADIKRAaaMrH\ngJP@DivhPNAbqy@\ngJQ@@dkSBJpHVOI@\neMJAhHzB\ngJPhLQxIRuPD',Hq='daD@@DiUVyjjPPGd\ndaD@@DjUZxHD@@\ndaD@@DjUZxHH@@\ndaD@@DjWjXHB@@\ndaD@@DjWzXHB@@\ndaD@D@BHBDBLBBBJBFBNBDiUVzjj`@\ndaD@P@bBbDfYvzB@@@\ndaD@P@bNBDfUzZ@B@@\ndaD@P@bNbDfUzZ@B@@\ndaD@P@qBdDfYvzB@@@\ndaD@P@qFdDfUjz@H@@\ndaD@`@BDeeVz`@@@\ndaD@`@bDfYVz@`@@\ndaDD@@IJVVfijh@@\ndaDD@@QIe\\jZjh@@\ndaDD@@YIeZn`B@@@\ndaDD@@YJZUnjjh@@\ndaDD@@qIYUnZjX@@\ndaDD@@qJYoJjjX@@\ndaDD@@qJZ_Fjjh@@\ndaDD@@yIe^f`@`@@\ndaDD@@yJYVfjjh@@\ndaDD@@yJYfnjjh@@\ndaDH@@RVU[f@@@@\ndaDH@@RVU[j@@@@\ndaDH@@RYVih@H@@\ndaDH@@RYWih@H@@\ndaDH@@RYe[hB@@@\ndaDH@@Rfu[j@@@@\ndaDH`L@HRf][jZj@@\ndaDL@@SDfUrijj`@\ndaE@@@IIf]njjh@@\ndaE@@@YIeZn`B@@@\ndaE@@@yIe^f`@`@@\ndaE@@@yIe^fjjh@@\ndaEH@DpDfYbYjj`@\ndaF@@@RVU[n@@@@\ndaF@@@RYVih@H@@\ndaF@@@RYWih@H@@\ndaF@@@RYe[hB@@@\ndaF@`NBHRYWih@H@@\ndaF@`NBPRYWih@H@@\ndaFD@DCdeeVyjj`@\ndaFH@BAIf]n``@@@\ndaFH@DAIYUnZjh@@\ndaFH@DAIeUnZjh@@\ndaFH@DAIfVfZjX@@\ndaFH@DAIf\\fZjh@@\ndaFH@FAIeZn`B@@@\ndaFH@NAIe^f`@`@@\ndaFH@NAJYfnjjh@@\ndae@@@yJeVn[jjj`@\ndaf@@@RiUkfzjjh@@\ndax@@Djuvjh@@\ndaxB@@QnR[VZY`cD\ndaxB@@QnR[VZY`cH\ndaxB@@QnR[VZi`@\ndaxB@@RnRUUZjP@\ndaxD@@KHhhbtu@@\ndaxD@@KHhhcUT`@\ndaxD@@QIUUjj@@\ndaxD@@QIUYjj@@\ndaxD@@QIeUjVBB`\ndaxD@@QIeUjZBB`\ndaxD@@QIeUjf@@\ndaxD@@QIeUjj@@\ndaxD@@QIeejj@@\ndaxD@@QIgijj@@\ndaxD@@QImUifALj`\ndaxD@@[HiDYUU@@\ndaxD@@iIijjj@@\ndaxD@@iJUfjj@@\ndaxD@@iJUvjj@@\ndaxD@@yIi^jj@@\ndaxDPDxHahaImZYiBL`\ndaxDPFxH`haIf^jj@@\ndaxDPFxLPlQIf^jj@@\ndaxDpBhHa@chaIffeZBDP\ndaxH@@RUUYf`QHh\ndaxH@@RV[jj`@\ndaxH@@RYUZj`@\ndaxH@@ReWii`QJh\ndaxH@@rJJHmUP@@\ndaxH`HALRkUZj`@\ndaxH`HALRkfZj`@\ndaxJ`HSFx@cIIJhmS@@\ndaxL@@SdfUvih@@\ndaxL@@SdfvVjh@@\ndax`@@PjRUnzZ`@\nday@@@[HiDYUU@@\nday@@@kIHbiUU@@\ndayL`BzDp@cIICdmU@@\ndaz@@@RfVjj`@\ndaz@PLBHzHRUgjj`@\ndaz@`JBHRV^jj`@\ndazD@DAdfUvih@@\ndazD@LCdeYzjh@@\ndazD@NADf{Vfl@@\ndazH@DAIeYjZ@@\ndazH@LAIUjjj@@\ndazH`BpJqI[Vfj@@\ndazH`LPHaInVZj@@\ndcL@@DiUUU]ZjejC@^QP\ndcL@@DjYYYiBHhh@@\ndcL@H@qBqEqMqDfYn]a``bh@@\ndcLB@@Q]R[e[neh@a`@@\ndcLB@@Q]R[e]nEh@I`@@\ndcLD@@EJY[WWZjB@@@\ndcLD@@GIEDTiMBUUUS@@\ndcLD@@IIf]z[hHBj@@\ndcLD@@QIUVUWVj`@@@\ndcLD@@QIe]UJfjjj`@\ndcLD@@QIgVUWVj`@@@\ndcLD@@QIge]FVh@I@@\ndcLD@@SHhdhTeSMUUU@@\ndcLD@@YJYYwEZB`b@@\ndcLD@@iJ[g]xZB@f@bX\ndcLD@@uIUfUeVXBB@aH\ndcLDPDtHaXaInUnzY`BI@@\ndcLH@@RYWZZYjjjh@@\ndcLH@@RYeZvz@`j`@@\ndcLH@@RYguRYijjh@@\ndcLL@@G\\dTRRbOKPPTP@@\ndcLL@@QTfVUV]Zjjj@@\ndcLL@@QTfVutYZ`@h@@\ndcLL@@QTfvUtYZ`@h@@\ndcLL`HQPBLddLRRTzmAA@@@\ndcLL`HS@BLddJfRtjmP@P@@\ndcM@@@eIi_Z[jjjj`@\ndcM@@@wIHhdd]JuAPD@@\ndcM@pItIAICICHiCDeDJuPAD@@\ndcMB@hDDWXeNF]yInUnzV`BF@@\ndcMD@DdMRUe]^EX@IP@@\ndcMD@LDMRUe]^Fjjjh@@\ndcMD@LtARem[~fl@bp@@\ndcMDBLHDSdf{YU]Zj@@@@\ndcMH@DdDfY}XfZjjj@@\ndcMLAHtDVISdfyW[aZ@BX@@\ndcNB@BAEuInVVFV`HF@@\ndcNBAHAEvISdfyW[aZ@BX@@\ndcND@AALbbTRbJzuAP@@@\ndcND@MADfVU~UZ``H@@\ndcNDaFePBFe^RfV^Qv`hF@@@\ndcNH@DAIf_UIfjjj`@\ndcNLAHAEbTyInUnzV`BF@@\ndcNL`NaL@HrRRqQ[RjtDA@@@\ndcO@@@JldbRTQU]UUUT@@\ndcO@@@xTjU]Znjjjj@@\ndcOLpHkb[PbAbEbDfYn\x7FijYjf@@\ndcl@@DjUg^aWP`hjH@@\ndcl@@DjYU_egX@@@@@@\ndclD@@QIe[UiiUjP@h@@\ndclD@@QIe[WiiUjP@h@@\ndclD@@SHheDhbJkkUUUT@@\ndclD@@aJVfnKfFjjjf@@\ndclD@@iJYW]rnF``IhBI`\ndclD@@iJYW]rnF``Jh@@\ndcm@@@UJfUWyYvhJBH@@\ndcmH@DHDfUyTjWVi`@`@@\ndcn@@@rJJJJJlJ{@A@@@@\ndcnD`EGPqDfUn^neX@aiH@@\ndcnH@AAJYYwhUvj``H@@\ndcnH@DCHhheBeSKkUUUT@@\ndctB@@PYRYU{Vjij@@\ndctD@@QIeWUZjjh@@\ndctD@@QIgUUZjjh@@\ndctD@@QImUUZjjh@@\ndctD@@QImUVZjjh@@\ndctD@@wHiDThmUUUP@@\ndctH@@RUgVVZijADn`\ndctH@@RV]UVjjj@@\ndctH@@RgUUZjjj@@\ndctH`HALRkUUVjjj@@\ndctH`HALRkUYfjjj@@\ndctL@@Pdf{e]jjj`@\ndctL@@X\\dTQfReUUU@@\ndctL@@jTiV[Vjjj`@\ndcvD@MADfuUUjjjp@\ndcvH@LAIVUUjjjh@@\ndeL@@Di[ernDYZjij@@\ndeL@@DjUYkfEijjjj@@\ndeT@@DjUghP`h`@@\ndeT@@DjWvifjih@@\ndeT@@DjWviifjh@@\ndeT@@DjYUXPbH`@@\ndeT@@Dj[[[ifjd@@\ndeT@@DjeUZZjjh@@\ndeT@@LdbRTm\\DBR@@\ndeT@@LdbRTm\\DBT@@\ndeT@@LdbRbmBDED@@\ndeT@@LdbTRoBuUM@@\ndeT@@LdbTRoBuUT`@\ndeT@@LdbbQwBuSU@@\ndeT@@LdbbQwCUUU@@\ndeT@P@bIbDee][j@B`@@\ndeT@p@bNBIbDfUuih@J`@@\ndeTD@@EJYU^f```@@\ndeTD@@IJVWZZfjj@@\ndeTD@@QIe]Rijjj@@\ndeTD@@QIgeQej@@@@\ndeTD@@QIgeQejjj@@\ndeTD@@QImeQej@@@@\ndeTD@@SHhdhUSMUUP@@\ndeTD@@eJ[WVz`@h@@\ndeTD@@iIYe^e```@@\ndeTD@@qIUeQej@@@@\ndeTD`AdHaIe]jZ@BX@@\ndeTD`DpHaImeQfZ@@@@\ndeTD`NDHaIfVVfBA`@@\ndeTH@@RUYTYY`@@aH\ndeTH@@RUYTYj`@@@\ndeTH@@RYVZfZZj`@\ndeTH@@RYVZfZij`@\ndeTH@@RYVffjjj`@\ndeTH@@RYWZfjjj`@\ndeTH@@RYe\\YZB@@@\ndeTH@@RYm_aZ@B@@\ndeTH@@RYm_ajjj`@\ndeTH@@RgYTYj`@@@\ndeTH@@rJJIHmtAAH@@\ndeTH@@rJQPiCMT@@@@\ndeTHPABHxHRYWZf`@f@@\ndeTH`DBHRZ{TYf`@@@\ndeTH`IBHrJJJJlLADP@@\ndeTHpACDKD[DRYf{i`bH@@\ndeTL@@QdfygFV``@@@\ndeTL`HS@BLddlRPrm@@@@\ndeTL`HS@|LddlRPru@@@@\ndeU@@@eIYVvG`BL@@\ndeU@@@eIYWVz`@h@@\ndeU@@@eIYWVzjjj@@\ndeU@@@eIYWV{`@l@@\ndeU@@@gHeDeBwT@E@@@\ndeU@@@qJYVjXBBh@@\ndeU@`Ld@aJVu~Fl@H@@\ndeUB@DpFFTfUgkfjYX@@\ndeUH@JXDiWUJjjjh@@\ndeV@@@RUYTYy`@@aH\ndeV@@@RVUenh@J@@\ndeV@@@rQQIHtpDCP@@\ndeV@@@rQQQHcMAP@@@\ndeV@@@rQQRItMMUR@@\ndeV@@@rRHqICMT@@@@\ndeV@PABPRPR[e[ij@H@@\ndeV@pABPJPZPRYf{i`bH@@\ndeVD@DAdfygFV``@@@\ndeVH@AAIfuneh@`@@\ndeVH@IAJYW~F``H@@\ndeVH@JAJUuRjjjj@@\ndeVH@LAIUeQfj@@@@\ndeVH@LAIVUVzjjj@@\ndeVH@NCIELeBpt@Q@@@\ndeVHAH@NbTfY_[hBB`@@\ndeVH`IHHaIfUVFBBH@@\ndeVH`NdHaIe]ZZ@BT@@\ndeVL@D@YrJJHjDsUTt@@\ndeVLAHAAbTyInUneh@`@@\ndeWH@DJPRY[TYZ`@@@\nded@@Dj_VfZZ@@\nded@@LdbQRdsSPQD@\nded@P@SHBDjuUZjj@@\ndedB@@PYRYWYjZX@@\ndedB@@PYR[UYjjX@@\ndedBpJxYCDKD[pRY]jjih@@\ndedD@@QIUUVjj`@\ndedD@@QIkWZjj`@\ndedD@@QInUvjj`@\ndedD@@aJVfjjj`@\ndedD@@eIffZjj`@\ndedD@@qJY]zjj`@\ndedH@@RUUUfihDRY@\ndedH@@RUUUfjhHR@\ndedH@@RUgVZjhHb@\ndedH@@RV]Ujjh@@\ndedH@@RYVYjjh@@\ndedH@@RZW^jkh@@\ndedH@@RZWvjjh@@\ndedH@@rJQEJUUU@@\ndedH`HALRkUUjjh@@\ndedH`HALRkVYjjh@@\ndedL@@PTfUvZfj@@\ndedL@@pdie]jjj@@\ndedL`@JdhDi[Ujjj@@\ndedL`HS@|DjUUjjj@@\ndee@@@eIUUZjj`@\ndee@@@{IHh\\iUUT@@\ndee@PJdLQlSHiDdUU]V@@\ndef@@@RiUVjiT@@\ndef@@@Ri]Vjjh@@\ndef@pLBHFHfHRUeVjjh@@\ndefB@LAAeIeUVjj`@\ndefB`HKad@aIf^Zff`@\ndefD@AADfUfZfjBBP\ndefH@LAIVYjjj`@\ndefL@L@YRUeVjjh@@\ndeg@@@FTjWUjjj@@\ndeg@@@JTeV]zjk@@\ndet@@DjYUX^d@@@@@\ndet@@DjYUX^dHbH`@\ndet@@DjYUZ^D`dJ@@\ndet@@DjyZkfyjjj`@\ndetH@@RYVvfFX@Jb@@\ndeu@`Dp@aIeURhYfh@@@@\ndev@@@Re[TjFP@@@@@\ndev@@@rQQJHtpr@@@@@@\ndev@`JADrQQQHyPsPLA@@@\ndevH@JCIEEDceCM@pD@@\ndevH`HX@aJVU|HYjA@`@@\ndg\\B@@Q]rJISQJHr[TDCD@@\ndg\\B@@SSRY[W[FVh@Ih@@\ndg\\B@@pSRf^Y]vzjB``@@\ndg\\D@@QIeUyT{Zjh@@@@\ndg\\D@@QIe]UTjZjjjh@@\ndg\\D@@QIgUYT{Zjh@@@@\ndg\\D@@QIge]hYZjjjh@@\ndg\\D@@QIgfVVSZjjjh@@\ndg\\D@@QIgfWzUZjZih@@\ndg\\D@@SHhdhTdeSMUUUT@@\ndg\\D@@SHhmHcDd]mUT@@@@\ndg\\D@@SHhmHhTmimUUTt@@\ndg\\D@@SHhmHhTmimUUUT@@\ndg\\D@@eIfU_Un`HJj`@@\ndg\\H@@RV^UveVj@Bh@@\ndg\\H@@rQIIQQHg]UPT@@@\ndg\\LPLi`bBBLbRTabRvgUK@A@@@\ndg]D@DdMRYe_[FVjZif@@\ndg]D@DtKrIIQIPiCJt@EJ@@\ndg]D`NV]CDRYWVyih@JZh@@\ndg]H@LlDeYeUuNj``c@@\ndg^D@L@teYee]nj`bH@@\ndg^DPKFPbNBLbbRrdJUM@AKM`Pl@\ndg^D`Ca@BLddLRRVTzm@aG@Pt@\ndg^H@DAIfVU^]Zjjjh@@\ndg^H@DAIf_UTfZjjjh@@\ndg^H@DCHhhhdhb]mTED@@@\ndg^`@@SAgHhihhhb]mTCD@@@\ndgl@P@BH|DfeUUVjjjh@@\ndglB@@RUrJIIIPiMUTuP@@\ndglD@@QIeWUUjjjj@@\ndglD@@QIeuUUjjjfBM`\ndglD@@QImUUUjjjj@@\ndglD@@SHhmHbeDuUUR`@\ndglD@@SHhmHdThuUMM@@\ndglD@@qJ[]UVjjjj@@\ndglDB@QNrJZQEIIMUMUPPt@\ndglDPDlH`xaImUUZZjji@@\ndgm@@@mIUUUVjjjj@@\ndgnB@DBcoHheEiDbuUUU@@\ndg|@P@bFBDfUoeZx{`BH@@@@\ndg|B@@Q[RYvUmgSmj`BJH@@\ndg|D@@SHieEDhcJg[T@@@@@@\ndg|DB@QNrJYQQJHrivuUUUT@@\ndg|D`BTHaIfUm]aNxHHD@@@\ndg|H@@RVYYwySn``@@@@@\ndg|H@@RYe[W[cn@`jjH@@\ndg|H`BBHRYg^U[cnB@`@@@@\ndg}@pBl@c@aPaJYgU|jgZdHHn@@\ndg~D`NTpbDfUvYj[[`@jYa@@\ndg~H@NAIe]YVfNx@J@@@@\ndg\x7F@PBWPbAbLbbbRfaSR]pPTMQ@@\ndiD@@DiUVjj`@\ndiD@`@BDiUVjj`@\ndiD@`@bDee^jj`@\ndiD@`@kDf]Vjj`@\ndiDD@@QIgUZjh@@\ndiDDB@QNR[UfjZ@@\ndiDDPBhLQlQIf_jZh@@\ndiDH@@RVUvjj@@\ndiDH@@RVvjZj@@\ndiDH@@RYuYjjBH`\ndiDH@@rIRHrjj`@\ndiDJ@@PnEInvZjX@@\ndiDJ`HSNDOAJoVZjX@@\ndiDL@@IdfY~jjP@\ndiDL@@kdiVZjj`@\ndiDLPBhPbFbLbbbeiZdHQ@\ndiDLPDp`bH|DfmVZj`@\ndiE@`ND@aJUUjj|@@\ndiEH@BxDeUZjj`@\ndiF@@@RUUZjj@@\ndiF@@@RVUzjj@@\ndiF@@@ReUZjj@@\ndiFD@AADfuUjj`@\ndiFD@BCdf]Zjj`@\ndiFD@FBdiUzjj`@\ndiFD@LCldTatjjh@@\ndiFJAHALXXgNRYvvjj@@\ndiG@@@HTeUWjjp@\ndiG@@@xTiUVjj`@\ndiT@@LdRbQrbYjjZ@@\ndiTH@@RVYz{ajjj`@\ndiTH@@ReUxRfjjf`@\ndiV@@@RVYz{ajjj`@\ndiV@@@RfU|kahDB@@\ndiV@`J@HRfU|kahDB@@\ndiVH@AAJUWaJZjjz@@\ndid@@DjUfaBB`@@\ndid@@DjWffB@h@@\ndid@@DjYUaBHP@@\ndid@@LdbRQk``b@@\ndid@@LdbbQxXF@@@\ndid@P@bNBDfUvf`@h@@\ndid@`@bDf[Waj@@@@\ndid@p@bBbAbDfYun``H@@\ndid@p@bBbFbDfYoa`b@@@\ndidD@@EIe]ih@J@@\ndidD@@EIfU[hBB@@\ndidD@@GHhdhZZ@B`@@\ndidD@@KIDhdZZfjh@@\ndidD@@QIe]Jfjj`@\ndidD@@QInUxV`@@@\ndidD@@iIYgxVB@@@\ndidD@@iJ[gxZB@@@\ndidD`BDHaIf_[hHB@@\ndidDpNDH`ha`cHhhcJZjiX@@\ndidH@@RUe^Fh@@@@\ndidH@@RVUvz`@`@@\ndidH@@RYVZZ@B`@@\ndidH@@RYWZZ@BP@@\ndidH@@RYm^Eh@@@@\ndidH@@RYm^Fh@@@@\ndidH@@RYm^Fjjh@@\ndidH@@RZU~Fjzh@@\ndidHPBBHFHRYgVzB@`@@\ndidHPBBHzHRYgfFBB@@@\ndidH`ACDrJIJFf`@d@@\ndidL@@HTfYun``H@@\ndidL@@IdfYoa`b@@@\ndidL@@KdfYynZej@@\ndidL@@RdfVwaZii@@\ndidL@@rldTTUkjjj`@\ndidlAHJfAAbTyInW[fiY`@\ndie@@@EIYW[n@B@@\ndie@@@EIfW[hBB@@\ndie@@@GHhhdVz@``@@\ndie@PLx@a@cHheDZyjjh@@\ndie@`BDHaIf][hHB@@\ndieD`JXaBPRYgvzejX@@\ndieH@JDDiWTjjjj@@\ndieH@LDDeYWajjj@@\ndif@@@RUe^Gh@@@@\ndif@@@RVUf{`@`@@\ndif@@@RVUv{`@`@@\ndif@@@RfU~Fjjh@@\ndif@@@RfWzXBBP@@\ndif@@@RfWzXBB`@@\ndif@@@rRJEKaj@@@@\ndif@PABHJHRYgVzB@`@@\ndif@`ABHRYWZZ@Bp@@\ndif@`BBHRYgVzB@`@@\ndif@`NBHRYeVFBC@@@\ndifD@D@TfY|fZjj@@\ndifD@J@TiWTjjjj@@\ndifDAHAHeNR[e^Eh@@@@\ndifFAHALkab\\yIgnxVij`@\ndifH@AAJYYhZjj`@\ndifH@DAIfVifjj`@\ndifH@DAInUxVjj`@\ndig@@@xTjU^njjj@@\ndkLB@@SSR[UUYjjjX@@\ndkLB`HSBCprRSEQIYjjjh@@\ndkLD@@QIeUoVjjj`@\ndkLD@@QIeWUVjjj`@\ndkLD@@QIe]Ufijj`@\ndkLD@@QIe^Uvijj`@\ndkLD@@QImUUVjjj`@\ndkLD@@qJ[]UZjzj`@\ndkLFPHSaTpBJBLddJTbNVjVj@@\ndkLH@@RUUUUfjjhHR@\ndkLH@@RVUU^jjjh@@\ndkLH@@rJIJJJrjjjh@@\ndkLH@@rJJQUIIjjjh@@\ndkLLHBi`bL|MbCbDeZuUjZjj@@\ndkLLHLi`bB|MbCbDeVUUjjjj@@\ndkLLXJXPdDdLdEdMdCdLbdLrfdjffj@@\ndkLL`@JdhDiYWUjjjj@@\ndkLL`@idhLdTTJTTjjjj@@\ndkLN@@PiWSR[kVYjjfX@@\ndkM@@@UJeUUZjjj`@\ndkMDHLve@HZPzPFPrQSJQVZfjihHR@\ndkN@@@RUUUVjjjh@@\ndkN@@@ReUUVjjjh@@\ndkND@CALbbRbJRZjjj@@\ndkNH@JAIgeUzjjj`@\ndkNL@BB]RYyVZjjjh@@\ndkO@@@Y\\dbJRRTjjjj@@\ndk\\@`@BDifUWGUN`@@@@@\ndk\\D@@QIenf^WSZjjjh@@\ndk\\D@@SHihheDQgSZjjjh@@\ndk\\H@@RYeg]itvi``D@@\ndk\\H@@RYeg]itxH@@@@@\ndk\\H@@RYm[Watv`@@@@@\ndk\\H@@RfU{WatzB@@@@@\ndk\\H@@rJJZIQDYtv`@@@@@\ndk\\H`ABHRYVvvftx@Jjb@@\ndk]@pLh@a@c`aIegmRkSZj`@P@@\ndk]H@DHDfYyWImMjhHA@@\ndk]LbLjDd@birRPqSIYntuj`PH@@\ndk^@@@RfYU\\]Tzjjjj@@\ndk^D@DBTfYYwimMjdDB@@\ndk^D@DBTfY[W[mNfhHB@@\ndk^HPNHJrXaIf^UvGS``B@@@@\ndk^LpB[a@HpHkprQPjSIILPjjZjj@@\ndkl@`@SLddTlTUgZhBH`@@\ndklB@@Q]RY[WoaZjjZ`@\ndklD@@QIeUY]Mjjjj@@\ndklD@@QIe]URijjjj@@\ndklD@@QIgfU}MjhD@@@\ndklD@@QIgge]Mjj@@@@\ndklD@@QInUvnEh@Jh@@\ndklD@@SHhdhTdjYjjjj@@\ndklD@@eJ[Vvfz`@jh@@\ndklD@@qJ[[m]Njj@@@@\ndklDB@QNR[V[WSZj`@@@\ndklDpLHHbOA@aIkVuVzZjjj@@\ndklH@@RYeZUn`HJj@@\ndklH@@RZWyWSjj`@@@\ndklH@@ReYUfSj``h@@\ndklL@@STf^UvFVh@J`@@\ndklL`AtpbDfUoVih@Ji`@@\ndklLpDppbHBMbLbdLRTtEYX@bP@@\ndkmB@hTDtxeNVS{HihhdUAehBBX@@\ndkmD@DTCRUfWtYV@`e@@\ndkmD@DdCrIJJIPxUV@bE@@\ndknB`NaLt@cIIKEEeiUZB@h@@\ndknD@LCTeYeUTzjBB@@@\ndknD`Ca@BLddLRRVgUhDHpHZ@\ndknD`LIPBDiZuYdzzjfh@@\ndknLPASBBH`HrJPjJKI]ehHB@@\ndknL`LaA@HrRPjIKY]VdDB@@\ndmL@@DjYYVgeBHhb@@\ndmL@@DjYeVdUjjjj@@\ndmLD@@QIe[VfeVj@B@@\ndmLD@@yJY~WJeZYijP@\ndmLH@@RYVuiiV@@@@@@\ndmLH@@RYVuiiVjjjh@@\ndmLH@@RYegXYVhH``@@\ndmLH@@RYiYKnUjjjh@@\ndmLH@@RYiYKnVjjjh@@\ndmLH@@RfYUxYVjifh@@\ndmLH@@RfYUxYVjjjh@@\ndmLH@@rJJIQEneX@@@@@\ndmN@@@RYVuiiVj`@`@@\ndmN@@@ReZ}IiVjjZh@@\ndmN@@@rQQJKFnFP@@H@@\ndmN@@@rQQQQVneP@`@@@\ndmN@`J@HrQQQH{JFZA`J@@\ndmN@`JADrQQQH{JFZA`J@@\ndmNH@BAIfUmiEX@@@@@\ndmN``DkaT@aIefUneVjVjP@\ndmO@@@SdfVUrjUZ`PH@@\ndmT@D@dDdLdJdFdAdIdEdLbdLadjjj`@\ndmT@D@kDdLdJdFdAdIdEdLbdLadjjj`@\ndmT@p@bBBEbDeYyzjjh@@\ndmTB@@RiRYyVZje`@\ndmTB@@SiRYyvZji`@\ndmTB`HSBCpRjuUZjj`@\ndmTB`NFU@HReUfjjY`@\ndmTD@@QIUUUjjj@@\ndmTD@@QIeUyjjj@@\ndmTD@@QIgUUjjj@@\ndmTD@@QImUUjjj@@\ndmTD@@iJU]Vjjj@@\ndmTH@@RUUUYjj`aH\ndmTH@@RUgVYjf`aD\ndmTH@@RUgVijf`aD\ndmTH@@RUgVijf`aH\ndmTH@@RVUWjjj`@\ndmTH@@RZWvjjj`@\ndmTH@@rJQEJJjje@@\ndmUH@JXDiUYjjjh@@\ndmV@@@RUUUjjj`@\ndmV@@@ReUUjjj`@\ndmV@@@RfWYjjj`@\ndmVD@DAdfUufijh@@\ndmVD@LA\\bRbdUjjj`@\ndmVL@IAARYW[Zff`@\ndmWH`DphCpRjyfZjj`@\ndmt@@LdbRdSjP`jH@@\ndmt@H@qAdIdEdDfYVfz@`j@@\ndmt@H@qAdIdEdLbbRadih@Jh@@\ndmtB@@RUR[e^[fjjZ@@\ndmtD@@IIf]yn``J`@@\ndmtD@@QIe]TjZjjh@@\ndmtD@@QIem\\YZ`@`@@\ndmtD@@QIgYVUZh@@@@\ndmtD@@QImYVUZX@@Hr@\ndmtD@@UIfuwaZ@B`@@\ndmtD@@aJye[ahBB`@@\ndmtD@@qJ[[VUjh@@@@\ndmtD@@yJUe^Uj``@@@\ndmtDPAdHc`aIe]jf`@e`@@\ndmtD`JTHaIYe_ihHHP@@\ndmtD`LTHaI[e\\Yi`@`@@\ndmtD`LxLQI[f^Ui``@@@\ndmtD`NTHaIe]Vf`@j`@@\ndmtDpDpHb@aXcHiCDeafV@B@@\ndmtH@@RUfueVZ@@BD`\ndmtH@@RVUv[j@Bh@@\ndmtH@@RYYf[ffjj@@\ndmtH@@RYeZ[hBBh@@\ndmtH@@RYeeZVfjj@@\ndmtH@@RYeeZVjjj@@\ndmtH@@RYeeZZjjj@@\ndmtH@@RYe~[ffjZ@@\ndmtH@@Rfuu[j@BXBAP\ndmtH@@rJJIHhfZZjh@@\ndmtH@@rJJIJEn`HJ`@@\ndmtH@@rJQPiXYjjjh@@\ndmtH@@rJQQHxYjjih@@\ndmtHpEBHZHfHrJIJHrn`BJP@@\ndmtL@@PdfueYUj`@@@\ndmtL@@YTf[gqehHB@@\ndmtL@@hTef^~e``b@@\ndmu@`ATIAIe[^n`BNp@@\ndmuDPHhfBHFHRYf~kjfZj@@\ndmuDPHhfCDFHRYf~kjfZj@@\ndmuH@DTDf^Uqej@B@@\ndmuH@DdDf^UaUj@H@@\ndmuLAHTDZISdfygQehHB@@\ndmv@@@RVUv[n@Bh@@\ndmv@@@RfUWzZjjj@@\ndmv@@@rQQJEJUjh@@@@\ndmvD@BADfueYUjjj`@\ndmvD@DCdf^YyUjB@@@\ndmvD@NADfVyyUjB@@@\ndmvD`La@BLddlTReUhB@@@\ndmvH@BAIVUwaj@B`@@\ndmvH@EAIeYfnZZjx@@\ndmvH@HAIffgejjjh@@\ndmvH@JAJUuTjjjjh@@\ndmvH@LAIUfVUj`H@@@\ndmvL@EAFR[f_FV``H@@\ndmvL@FAIR[ev[fjjZ@@\ndmwLPHkbYPbAbLbbbfezZfZf@@\ndndD`LH@aJZ\\hjijjj`@\ndo\\B`HRUALrRQIIIHjjjej`@\ndo\\D@@QIe]UUZfjjh@@\ndo\\D@@QIe]UVZfjjh@@\ndo\\D@@SHhdddhbfjjjj@@\ndo\\F@@Savtf^UVYjjjY`@\ndo\\J@@PYWHheDhbdfjYjj@@\ndo\\J@@QiuIe[YoZfifT@@\ndo\\L@@ptimyujjjjj`@\ndo\\L@@rTie]mVjjjj`@\ndo^@`L@HReuUUZzjjj@@\ndo^J@JAIW\\bbTRTtqZjifh@@\ndo|B@@RgrJYIQQIGSZZ@`hBL`\ndo|BHDrwBH`HzHNHrJPqQYIRUe``JX@@\ndo|D@@QIeUyemvij@`@`h\ndo|D@@QIe[e~WVjijf`@\ndo|D@@QIgfUYuvj``h@@\ndo|D@@QImYeeVvjjjj`@\ndo|D@@SHhmHbhdk]jjh@@@@\ndo|D@@qJY]YWmzjjB@@@\ndo|DpElH`hbXcHhhidbeFFBBJj`@@\ndo|H@@RV^UviUj`@j`@@\ndo|H@@RYfV{vF@bJj`@@\ndo|H`EBHRYWUnjZ@Bij`@@\ndo}@DDhJs``YBYAY@yByAyCHheEiEYBjfjjjd@@\ndo~@@@RfUYw~V``hj`@@\ndo~@@@RfUYw~Vjjjjh@@\ndo~B@GAEoHheDdhbiMjBBb`@@\ndo~BpIx[\\H``cXcHhheHe]FVBAXi`@@\ndo~D@AADfUvUwSZf`@h@@\ndo~D@LALbbTTJTLWVj`@j@@\ndo~H@DAIVYenEU`Hbj@@\ndo~L@DA]rJJKIJHhfZjjZj@@\neFA@HoBJ@\neFABHiBL@\neFABPiBL@\neFACDlRL@\neFBBHc@@\neFBBlc@@\neFHBJ@\neFHBL@\neFJBhHh@\neFJHbHh@\neF`BL@\neF`BN@\neMB@HRZ@\neMB@Hch@\neMCAD`aBHu@\neMFI@bMP@\neMHAIX@\neMHAIh@\neMLRRWhv@\neM`AId@\neM`AIh@\neM`AIx@\neM`BN`@\neMhDRV@\neOB@Hcfh@\neOHBNZ`@\neO`BNZ`@\nfHbd@@\nfHcp@@\nfHdP@@\nfHd`@@\nfHdp@@\nf`a@`@@FrJIJQQNKTXuUUTu@@@\nf`a@`@@HrJIJQJqIN|uUU@@@@@\nf`a@`D@HQvQSRIIFIIwfjjh@@@@\nf`aA@@@ILjjjj{sUUUUT@@\nf`aA@@@ILsKWRpTADUUP@@\nf`aA@@@YEEDcEEHc^ZZjZj`Pq@\nf`aA@@@YEEDdTddf\\`HJjj@@@\nf`aA@@@YEEEETddfB`Hbjj@@@\nf`aA`@@HqdTRTTbUrMYj`j@d@@\nf`i@A@@AG@cIEEEDedkSdcZ]PA@AE@@@\nf`i@a@FRAD^bDfUm[WirRkN`BJ@BH@@\nf`iA@@@ILrmtvQcAV|uUUUUT@@\nf`iA@@@YEEEEDhdYFTzwf`@``H@@@\nf`iA@@@YHhedhh]HHUhujjjZjj`@@\nf`iA@@@YHhhdXdhhHehujjjZjj`@@\nf`iPB@N^@DYHhhhdldZ\\d[Sj@H@Hx@@@\nf`iQA@B\\|@HpDISLzsnRdcN}TaAQUD@@\nf`q@B@@^BULrjkj\\[uPA@D@@@@\nf`q@`@@^RVYYywbkNXH@BA@@@\nf`qA`@@B|dsLrozkF\\t@P@D@@@\nf`qA`@@H\\drsJkrkF|u@@@P@@@\nf`qA`@@HidTRbRQbREIwfjZh@D@@\nf`qPA@AZ@DPBLdRbbbTVJgV}T@PP@@@@\nf`qPHAHV@cGIF]zDxzt~{HhmEDdThJRmMUUUST@@\nf`qPQ@INxJw`QA@cHhheHeUDqQgPPKELd@@\nf`qPaABLEh@Q@HIKIEvnrQJJJIQXj][uPAA@@@@@\nf`qPb@EZ]xDQdRTRTraTxIVj`@iib@@\nf`qa@@D@rJJJJIZHjl{sUTEDD@@@\nf`qa@@D@rJJJJJFHjl{sULEDD@@@\nf`qp@@@Hpds\\rj~gV|uP@@@@@@\nf`y@@@LdbbbbbfkEBMIuo@@@@@@@@@\nf`~@@@DiUUUUjjjjj@@\nf`~@H@@Hdhug^rJZJIQJFMUTts@@@\nf`~@P@@Ht[HhdhbeLhuUUMT@@\nf`~@P@@HuYIe]UVvijjfP@@\nf`~@`@@BrQQIIQSVUUUUU@@@\nf`~@`@@HRYuUUUjjjjh@@\nf`~@`@@HR[UUUUjjjjh@@\nf`~A@@@YEDeEDcEJjjjj`@@\nf`~Ab@@Dhb`rQJJIEJJUUUUU@@@\nf`~Ab@@TXb`RfYwf^jjjjh@@\nf`~aK@JPQaVcV]xG`KpexYIIUDeHiFjjff`@@\nf`~a`BH@QbDxyIgUUUVjjjj`@@\nf`~ab@IPQ`G`eUnoLsUUUUP@@\nfbe@B@@QBSLjojjheUAAEUUP@@\nfbeAP@@Hu[rSKJkK~lcMPPTUUP@@\nfbea`@D@XgHhhdddcDedFKUUUP@T@@@\nfbm@A@@IdIAJkfYUU|tWICfj@@@@@@@@\nfbmAB@H@SDjnYeUWsQ\\dNZh@@@@@@@@\nfbmAB@H@SDjnYeUWsQ\\dNZhJhJhhh@@\nfbmAB@H@SDjnYeUWsQ\\dNZjjjjjjh@@\nfbmA`@@HsdTRfbTbtjRejshiZjjj`@j@@@\nfbmPa@LN``@QhHrRPqPiSISIJL{r\\kU@DT@A@@@\nfbmp@@@XXeL\x7FLjsobe`iruUPP@@P@@@\nfbu@AP@QAHadPJHeDZbGQ@XaLQfHrJJJ[PrZJ[TyzMLuUULuP@@\nfbu`q@OPQbJ\\`B@aFRRVIJYJYZPTcVV`@@Je`@@\nfbua@@M@rJQEJKJJIKNBgMUUSTBA@HR`\nfby@`@@HRYuUUUUZjjjjj`@@\nfby@`@@HR[UUUUUZjjffj`HTv`\nfby@`@@HR[UUUUUZjjfjj`PT`\nfby@`@@HR[UUUUUZjjjjj`@@\nfbyAb@HHpCpRkUYeU}Zjjjjj`@@\nfbyP`@LR@aIemUUUujjjjjj@@\nfbya@@D@rJJIIIHiIIMUUUUUP@@\nfbya`@D@UkHhhdddbddduUUUUU@@@\nfb}@P@@H]gHheDheEeD\\jugQRgKUUUUUUU@@@\nfb}AA@K@qLxbSLl|}z{ARtYtTp@@QL@P@@@\nfde@@`@QAHa`QhHrJJJZIKPi\\EjpXpQAU@@P@@\nfde@`@@HrIJQQZQXiRljw`mUUUUUT@@\nfde@``ARADDb@qDXaIf]UUUYqVg^``Jh@h@@@\nfdeP`@DNA[HhhhlmMCDRbc^BuPPPuTP@@\nfdea@@H`rQQQIRIZHbmHq`tEAUETT@@\nfdi@`@@BrQQIJHyII\\UPPTUUP@@\nfdi@`@@HRYUUUUPqZjjjjj@@\nfdi@`D@HQrS\\jjrnFKUTuA@@H\\`\nfdiA@@@ILjjjjhXmUUUUT`aE@\nfdiA`@@BCdbbRTTbWTxjjjjjjh@@\nfdiA`@@HLdrjjjjabuUUUUT@@\nfdiA`@@H\\drjjjjabuUUUUT@@\nfdiP`@DVAGHdiEBeCHbkRuP@UR`@@\nfdiXpCHLDXLPHUiwhRcfRJ]]{lbfbbbbbblkMADMCT@@@\nfdi``DF@Pcdfye_UTRf``ajj@PU@\nfdia@@D@RYguUURLZjjjjj@@\nfdiaB@JDAbYEEEIkHhT`qjZjj@@@@\nfdiaPBH@QYwhRclbfbbRLTTpiM@PKLt@@@\nfdq@P@@H]yIeWUUmZjjjfh@@\nfdq@`@@HRYWUUUVijjjj@@\nfdq@`@@HR[UUUUVjjjjj@@\nfdqAa@HRTALhDYIHddiDeTjjjVjj@@\nfdq`@@@ITnjjjmUUUUU@@@\nfduA@@@ILsLjm{AJ\\XOhm@@@@A@@@@\nfdy@`@@XrQQSPzIKJXiYuUUP@UP@@\nfdyA@@@IKLjoSWAF|pDEUUQ@@@\nfdyA@`I@bBQGhbLQdTTVVfbRaYqpXHBh@b@@@\nfdyAB`M@qBXcLPVHkDRYf{}e^RCFBHbh@@@@\nfdyAP@@BUhNQQQQQQDqIdgAZYjZZj`@@\nfdyAP@@BtFNQQQQQJKGISk^Z@HH@h@@@\nfdyA`@@BudbbRTQrbRxhLZBBb`@@@@\nfdyAb@FAb@HrQQJVIKJJDihuRl@ESP@@\nfdyPPBHZ@aShRcdfyVyuwdmFZ@HZfh`@@\nfdyP`@KN@cHihdliEDcpPXm@AUMTP@@\nfdy`R@KPQcp@cIIKDeLeBdJQgKP@@QT@@@\nfdy``@A@Pdrrr\x7FZjLFKTEA@@@@@\nfdya@@D@rJJIIHiQIBTFKUUUP@P@@\nfdya@@D@rJJJIGIQIRTFKUUUP@P@@\nfgA@B@@XbSLrjyFRMjuP@@A@@@\nfgAP@B@BAHirRJIJHkLUpQkAPP@A@@@\nfha@B@@IbURjjjmUUUUT@@\nfha@P@@HiKHhmDeDTdfjjZjh@@\nfha@R@HHehBYddbRRbbSUUULuT@@\nfha@R@HHpPG`eUjjsLuUUUU@@@\nfha@`@@HrJYQZQHiQjjjZZ@@\nfha@p@@HhkQdfV^yuVijVjhDAP\nfhaAP@@HEirSJwkO[TuUTr@@\nfhaAp@@HYqV`iLkWl{mSUSSP@@\nfhaAr@HHpSQk@|DjmUUUVjjjjh@@\nfha`@@@IRlrj~uUUUUP@@\nfha`B@B@bDeUUUUYjjjjhDHP\nfhah@@@\\ExJTjjjjmUUUUT@@\nfhi@R@NJmx@QdbbbbTVadgFCPAAAD@@@\nfhi@a@AFADAbDfUmygZL]z@BhHB@@@\nfhi@a@OAADBBLbTTRbfQtyiUj@B@bX@@@\nfhi@b@HH@DYIHXhhddSiuoKP@P@P@@@\nfhi@c@BFADRb@qFQQQIIPq[LD{tDDPPD@@@\nfhiA@@@ISLwZkf|xMUP@@@@@@\nfhiA@@@YEEDhdcLeKW`m@dP@@@aR@\nfhiAP@@PIqRTmMJ{yoNCUUUUUT@@\nfhiA`@@HBdwLl|zKSoMA@@PP@@@\nfhiIP@DXxHDc^CdTRbfaTVUNZltuUUUL@@\nfhiP@@@ArJJIEJIJYgKN`HJ@B`@@\nfhiPA@B\\@DXBDif]YybTdZiBBbj`@@\nfhiPC@IF@|PBDAFRIJJJIQXz|xMT@PP@@@@\nfhiQ`BHL@`LIQrS\\roZxeN|tDD@A@@@\nfhiQaBOAQbpHb@QzP}dTabbTVeTmF]K@QBrd@@\nfhi`@@@ISLwZkf|xMUP@@@@@@\nfhia@BH@AtIRYeV]WEjwjP`@BH@PP`\nfhiqP@DXxBQoArJIQSPjKJ`mVZZjjjf@@\nfhqA@@@ILkKjj~BuUUUUP@@\nfhqA@@@YEEDhbdddsdAAUUT@@@\nfhqP`@DD@YIeYU_UQfjjjjj@@\nfhqh`@NFmxHJYHihhdkLhIUADQML@@@\nfhyA@@@YHiMEiEbdHtcN}UUUUUU@@@\nfhyA`@@BMdTTTTTTVoMJ|xKPAAPA@@@\nfhyq@@B\\DAIfV]UudeZ|FBBaZeJ@@\nfle@`@@LRfVVW]|QgAhJ@HBF@D@p\nfleA@@@ILwJkttJU`m@PUUUD@@\nfleA`@@HbdssLjoyoNBuT@@@T@@@\nfleP`@DA@eIVUue]^B]yX@J@BT@@@\nfleP`@DA@eIVUue]^B]yZjjjjV`@@\nflea@@D@RYeUmU|RgAZjjBBJ`@@\nfleac@FV|DD`RPEHYEEEMDbeierVcPQDTluD@@\nfleiR@AZlFLDXhdGaddeRRbbRtQk^BuUURmSP@@\nfli@@@LdbVRRbbjTjjjjjj@@\nfli@RBKN|DBX{`|eTjjkLmUUUSMP@@\nfli@`@@HrJIJHiIIIIjjjjjh@@\nfli@`@@HrJKRIEJKIQjjZZihDEH\nfliAb@@Dhb`ReeUU]Vjjjjjh@@\nfli`@@@IJjjjjkUUUUUT@@\nfli`@@@IJjjjjkuUUUUT@@\nflmA@@@YEEEEDeDTbrRgF}EKPLQTKTR@@\nflmAB@K@qDfYYy{uXJVcNbf@@BI`@@@@\nflu@B@@^BULsJjozCZOdu@@@@@@@@@\nflu@B@LTBUJssZktDzHTaAP@@@@@@@\nflu@`@@HrJKRJKRKEISmNCEjjjjjjh@@\nfluP@@@ArQQQQII[ERdmNmxH@@@Bh@@@\nfluh`@DXUipAYEEEELTdedrVcN|uUMSURt@@\nfluqB@AJMDFHdrmrlzmNZ\\EL@PuA@a@@@\nfly@PB@HbeGI{HiiBddddeBTZjjjjf`@@\nfly@`@@XrQQISHqIYQwjjjh@J`@@\nflyAP@@DExNQQQPkIIJYTXtuMUSU@@@\nflyPa@L^bdDQjHrIJJIIHiSQQjjBBJi`@@\nfly`P`LHEIpH`LaFPSHYEMEEDlbeIjtu@pQUP@@\nflya@@D@rJJIIIHiIIVfjjjjjh@@\nflya@@H@RV[WUUTTih@Jjj`DHH\nflyaR@IBlzJPUfQQIYSJJHjTY@AMTuL@@@\nfoA@@@LdbbQbRbDkQfjijjh@@\nfoA@@@LdfbTrbUipufjjjfhDCH\nfoA@P@@Ht[HheMEEEQJeYjjjZi@@\nfoA@R@HHqp@QddebbfRtmV\\mA@@@@@@\nfoA@`@@NrQQQIISDVeVh@@BH@@@\nfoA@`@@VRfUYu^JLz``@B@@@\nfoA@`@@\\RVUoeVBly`BH@@@@\nfoA@b@KN@DISLjohmJMP@@AP@@@\nfoAA@@@IRlrjzkF]U@@@@@@@\nfoAI@BDDXhCIUNQQQSQIHdyYsURE@P@@\nfoAP@@@NRYWVUzLMZ@B`@`@@\nfoAP@@@NRYWVUzLMZjjjjj@@\nfoAPPBHZ@aShRcdfyVywdmFZ@HZd`@@\nfoAQC@DL\\BHiDVbLbbTTJRtYqVhB@@i@@@\nfoA`P@L@IrRRllrkhugLuSMThHHP\nfoA``@D@ydTRbTLRIFFluR`PT@@@\nfoA``@L@QdTVbbbblmV\\u@A@@@@@\nfoA`b@LTX@HRe{VYSVlznZjjZ@@\nfoAa@@A@RfuueV\\lz`@jjb@@\nfoAab@KPQ`@QddebRfR`iF\\m@@A@@@@\nfoAab@NPQ`@QddebbfRtmV\\mA@@@@@@\nfoAqB@JLDPDQdTTTtlRWAV]RuT@B@@@\nfoQ@@@LdbTrdrLYs`cVZjjjjh@@\nfoQ@`@@HRYVvvQRuXsfiZjjY@@\nfoQ@`@@HrJIQJKPjReFlyjjjij`@@\nfoQ`@@@YHhhheDbKTkF]P@PPD@@@\nfoQa@@N@rQQQQJKGbiVLz`BB@D@@@\ngC`D@DPHRfh@\ngC`DADZHRVXRP\ngC`DADZHRVh@\ngC`H@DIKT@@\ngC`HADIKLIH\ngC`HADIKR@@\ngC``@dfZ@@\ngC``AdeY@@\ngC``Adej@@\ngCa@@dkH@\ngCa@@dkP@\ngCa@@dmH@\ngCa@@dmP@\ngCa@@dsP@\ngCa@@duP@\ngCaHH@bNt@@\ngCaHL@aIZ`@\ngCaHLHaIZ`@\ngCah@mJAIj`@\ngCd@Adej@@\ngCe@E`dsP@\ngCh@@dmH@\ngCh@@dmP@\ngCh@@doH@\ngCh@@doP@\ngCh@@duP@\ngChHD@aIU`@\ngChHLHaIZp@\ngCh`LDdsP@\ngCi@DDeV@@\ngCi@DDeZ@@\ngCi@LDej@@\ngFp@DiTt@@@\ngFp@DiTvjh@\ngFpHADILimSP@\ngFpL@DXHHPeJfuU@@\ngFp`@dfTujh@\ngFp`@df_Ejh@\ngFp`AdigVjh@\ngFq@@drfmU@@\ngFq@@drftm@@\ngFq@@drfuM@@\ngFq@@eLzts@@\ngFq@@eLzuU@@\ngFqHJ@aJUMjj@@\ngFt@ATiTvjh@\ngFtHE`DILikUP@\ngFx@@eJf`@@@\ngFx@@eJfuU@@\ngFy@JDiTvjh@\ngFy@LDiWFjh@\ngGP@DiVV`iJ@\ngGP@DiVj`@\ngGPBADJPtaXcHiCUP@\ngGPD@DPHReZj@@\ngGPDADFHRYjY@@\ngGPLADFHlPdrsTID\ngGPLADFHlPdruT@@\ngGPLAbGDlPdruT@@\ngGPP@cTfyi`@\ngGP`@TfYi`@\ngGP`@TfYj`@\ngGP`@dfUjP@\ngGP`@dfUj`@\ngGP`@df]jP@\ngGP`@dfuiaM@\ngGP`@dfuj`@\ngGP`@dfyj`@\ngGP`ADkjj`@\ngGP`ATeVj`@\ngGP`ATeej`@\ngGP`ATf^j`@\ngGP`ATiVj`@\ngGPhH`DYIHUi@@\ngGPhMQDIK]U@@\ngGPlEaDPHlQdTaeh@\ngGQ@@djuT@@\ngGQ@@dkUT@@\ngGQ@@dlmT@@\ngGQ@@drmT@@\ngGQ@@druT@@\ngGQ@@dsML@@\ngGQ@@dsMT@@\ngGQ@@dtuR@@\ngGQ@@eJttHh\ngGQDBHbqBSKTp@\ngGQDJH`qBSKUP@\ngGQDJLPxbSKUP@\ngGQDLHbqBRwSP@\ngGQDLHbqBRwUP@\ngGQHDHaInfh@\ngGQHJ@aJUjh@\ngGQHJHaIUjh@\ngGQLJIARFdLbdMU@@\ngGQ`@bdwMT@@\ngGQ`@jdsmR@@\ngGQdEb@bRFRRVV`@\ngGQh@ZjAJVjh@\ngGT@Ade[j`@\ngGT`EPTfyi`@\ngGX@@dj|tHd\ngGX@@dkUT@@\ngGX@@dk]L@@\ngGX@@dtuV@@\ngGX@@eNuT@@\ngGXHD@aIUVd@\ngGX`hEIWIMkU@@\ngGXhhZ@bS^Rmjf@@\ngGY@BDeVj`@\ngGY@DDeUZP@\ngGY@DDfYj`@\ngGY@HDefZaH`\ngJP@DknX@\ngJPDADFHRYj`@\ngJPDADFHR[f`@\ngJPH@DIJuP@\ngJPH@DIKUP@\ngJPH@DIRuP@\ngJPHADILth@\ngJPXHlPDQzt@@\ngJPXHlPiQzt@@\ngJP`@TeZh@\ngJP`@TfVd@\ngJP`@TfZh@\ngJP`@deVd@\ngJP`@deVh@\ngJP`@dfVh@\ngJP`@dfvd@\ngJP`@dfvh@\ngJP`Adizh@\ngJPdEaDPHRZe`@\ngJQ@@djsBJ@\ngJQ@@dkU@@\ngJQ@@dls@@\ngJQ@@dlu@@\ngJQ@@dmS@@\ngJQ@@drt`@\ngJQ@@dru@@\ngJQ@@dsT`@\ngJQ@@duU@@\ngJQ@@eKS@@\ngJQ@@eKU@@\ngJQHBHaIfj@@\ngJQHBHaInZ@@\ngJQ`@bdvu@@\ngJQhHlOAJmj@@\ngJT@@TeZh@\ngJT@@Te^l@\ngJT`E`TfVh@\ngJU@HPdkU@@\ngJX@@dkU@@\ngJX@@dkt`@\ngJX@@dku@@\ngJX@@dms@@\ngJX@@dmu@@\ngJX@@eKU@@\ngJX`LDdru@@\ngJY@BDeZh@\ngJY@BDeZl@\ngJY@BDizh@\ngJY@DDefh@\ngJY@DDfvd@\ngJY@LDeZh@\ngJYHCabIKTp@\ngJYhCE`DQzt@@\ngKP`Adi\\Zj@@\ngKX@@eKcUP@\ngK\\@ABeKcMH@\ngNpB@DSppPPaJ[Zj`@\ngNpJAbJHLaYArBS]UU@@\ngNpLADXH\\PdjmUP@\ngNpP@jtfvZf@@\ngNp`@TfYZZ@@\ngNp`@dfUZe@@\ngNp`@dfUZf@@\ngNp`@dfUZi@@\ngNp`@dfUZj@@\ngNp`@dfWZfDL@\ngNp`@dfWZj@@\ngNp`@df]Zj@@\ngNp`@df^Zf@@\ngNp`@df^Zj@@\ngNp`@dfvZj@@\ngNp`@tf]jj@@\ngNp`@|bdQjj@@\ngNp`ATf^jj@@\ngNp`ATiUjj@@\ngNp`ATiVjj@@\ngNphDqDYEHcUR@@\ngNphH`DYIHTmJ@@\ngNphJpDIRkUT@@\ngNphJqDIKMTl@@\ngNq@@dlkUP@\ngNq@@dl}MPb`\ngNq@@dr{Vpa`\ngNq@@dr}Vpa`\ngNq@@dsKSP@\ngNq@@eLuUP@\ngNq@AdTbMTpa`\ngNq`@fdskUP@\ngNq`@jdrkUP@\ngNq`@jdssTp@\ngNq`AVeJmUP@\ngNq`AbeMmUP@\ngNqhAbjAJyZj`@\ngNr`ijpiJyImfi`@\ngNt@@TeVzj@@\ngNt@@|dbJjj@@\ngNx@@djmUP@\ngNx@@dsuUP@\ngNx@@eJmTh@\ngNx@@eJmUP@\ngNx@@eRmUP@\ngNx@AddQUUP@\ngNxDHHaQBS]UU@@\ngNx`DFdskUP@\ngNx`LFdjmUP@\ngNy@BDf^jj@@\ngNy@DDfYZi@@\ngNy@FDeYjg@@\ngNyHJPDIJwTt@@\ngOp@DjWkB@@@\ngOp@DjWkjj`@\ngOpH@DILkW@@@@\ngOp`@dfUMZj`@\ngOp`@tiguif`@\ngOp`@tigujj`@\ngOp`ATeekZj`@\ngOq@@drm[ST@@\ngOq@@drm[UT@@\ngOq@@drm\\@@@@\ngOq@@drm]UT@@\ngOq@@eLnmLt@@\ngOq@@eMN]UT@@\ngOq`@ldrikUT@@\ngOqhHl@cIIBjujh@\ngOt@@tjWkjj`@\ngOt@ATiUkjj`@\ngOtLHPDXHhPeLq]UL@@\ngOx@@drm\\@@@@\ngOx@@drm]UT@@\ngOx@@eJqh@P@@\ngOx@@eJyh@P@@\ngOx@@eLmXD@@@\ngOy@DDfYKZj`@\ngOy@JDiWMjj`@\ngOyDEQDDHRYXnZZ@@\ngOyDLpDHHRY\\ujf@@\ngO|HDVHaIeZx@@@',Iq='fHfXAa@\nfJ@FD\neOHBNZ`pge@\neFhHba@\ngFy@DDfXujhB\neFHBJFE@\ngFt@ATiTvjhCAbKF`\ngChHLHaIZ`H\neMhHchLJyH\neFB@JcAaGJ@\ngGPH@EIJmU@P\ngGPHAbIJmU@XTQX|e@\neMBBlRZCAKd`\neMABHXaIhH\ngJPDADXHRUj`LBHmrD\ngJQhHl@bOV`LH^S@\neMHAIdLF^P\ngCa@@dkHFBVyH\ngC`HABQz`H\ngCaDL@b^BTt`P\neMJDbDfP`\ngC`D@DXHRUdB\ngGQhHl@cIHTmPFBMy@\ngJXHD@aIUj@pHVOI@\ngGX@@eMUTA`Uc^P@\ngChDL@aABSM@XHKdp\ngOx@@drm\\@@A`plZp\ngJX`BDdru@P\neMC@thabHzB\ndaDD@@QInXjZjh@`\ndid@@LdbLTifjj`B\ndeTH@@RUYTYY`@@aKChSFxYyE@\ndaD@P@bNbDfUzZ@B@C@`pnxbT\ndeTHPBBHfHrJJIHlLDDP@P\ndaD@`@bDeeVz`@@B\ngFp`@TizJjhCCSGbU@\ndidD@@yIYVXV@H@CBhpf{dB\neFJHbHpP\ngO|@AFeVm]UTAaUcWrP\ngOy@FDiguie`LMj~ID\ndaD@`@bDfUjZ@B@CB`SJ{dL\ndifH@DAInUxV`@@CBdinGdD\ndaEH@HpDeeVz`@@B\ndeVD@AADfVuFVijh@phj[iy@`\ndifH@HAIYexZ@`@B\ndifH@JAJ[gxZB@@CBdJf{dB\ndifH@HAIYW[j@B@B\ngOx@@drm\\@@A`Qj~ID\ngOr@Ajti]qZY`H\ngF|@AbeJfuU@P\ndmN@@@rJJIQEneX@@@@B\ngJX`DBdru@XKGdP\ndid@`@qDeYWaf@@BH\\NABinGdP\ndifD@H@Tee]nh@H@H\ndeVD@LBTeegXV@HP@`\ndaF@@@RYWih@H@LJCBknP`\ndie@@@aJyV[f@B@B\ngGY@JDeUjpLDmrP\ngFt@AdigUjXCAF|Rp\neMBBHRZCAGe@\ngC`HADIKTAa`mrP\nfHcdA@\ngJX`HLdmU@P\ndaz`@@SFyIeYjf@LJAL[nQP\ngNtDLpDHHR[UjhB\ndiTH@@RgeXSaj@B@CAhJfx^P`\ndew@@@pldTTJVTLmP@P@XUaMt|`P\ndiUH@HPDiUWBxZVeh@`\ngNq`AVdlmUPFEfM_DD\ndaxL@@[deiZjh@pXHpf{dT\ngCaHLLQIZ`LLInQ@\ndedBPLinCDKpRUZvjZh@`\ngNplLPDSpPPdukUPFFE_H`\ngJPL@DSpPPdju@XXK\\TH\ndiDDpLH@bOA@aIiUZjh@`\neMIDBKpRYCqH\neMXIDf``\ngG]@DZDfuiPLEkrH\ngGY@JDf]j`LLl^R`\ngJX@@dlu@XZX|PP\ngNxHF@aJZzjPH\ndazD@LADeUffhHr`\ngGQHJLQIUjdCCBKbq@\ndaxD@@QIgUjj@LBrf{dD\ngNy@BDf]jj@phbuxbp\neMXIDfP`\ndiGL@DK`hTfV^ii`B\ndmvD@E@dfYwVzB@j@C@PpjxYT\ndmvL`BaL@HrRRjIJUVjjh@`\ndaFH@NAIeVf`@`@phLInyB@\ndiGD@JxPQIeUZfX@`\ngCi@DDej@ppfyD\ngCa@@dsPFBV@\ngCd@Adej@ppVyH\ndiFD`Li`BDenvjZ`CA`aLxP\neMC@HoABDe`pQyP\ndiDDpLH@bOA@aIkUZjh@`\nfHdpAa@\nfHa@A@\ngFy@JDigVjhCCSGbU@\ndaF@`L@HRe]pjjj@LFBDJfyG@\ngChHL@aJZ`LLHnS@\ndaFD@DCdfVRiji`CCDinQ`\ndaGH@DK`R[e[fii@LLPfyE@\ngOtHLPDYHhckSM@XXI|a@\ngNqhHl@cIICej`H\ngNphJqDIJmSLHf\ngJQ@@dju@XZK\\a@\neMABHYAIhH\ngJP`@TeVhCAQ\\VH\ngCahHlOBNtAaAqK@\ngGQDLHbqBRwUPFEDQoI`\ngJT@ADjzhCAQ\\VH\ndazL@BANR[UZj`CBdrf{dD\ngJX`LHe[U@P\ngNx`LDdssUPFFq_IP\ndaz@@@ReUjj`CC`aLJnyC@\ngNx`LFdjmUPFCDQkyL\ngChHHOAIj`H\ndeUH@JdDin]xZB@`@`\ndieH@BDDfY}n``H@LJpjX^P@\ndaE@@@yIe^f`@`@pHLKnHe@\ngOq@@dlvKUTAaeWqP`\ndid@p@bBbAbDfYun``H@H\ngOy@FDiguie`LMc^IL\ndaDh@DInAIf]nZiX@p`nIe@\ndigH`LH^@HRf_ljfYh@`\ndie@@@EJ[W[j@B@CBlJfGb@`\ndaFH@NAIe^f`@`@pILknHA@\ndieH@BDDfY}n``H@LBpj[b@H\ndeUD`JXiCDrJJJXlMKUL@P\ndaED`LJDCpRjw[fjj@LLinxfD\ndaDH@@RYUih@H@LBCB{bHp\ndaDH@@RVU[f@@@LB`j[bYp\ndmuD@FDDR[f[eV`b@@H\ndaDD@@yIe^fZVX@pSJxbD\ndcN@`LCDRUgeUvjh@@@`\ndk^@@@RYWYVftx@H@@@H\ngO|HEfHaIeZx@@B\ndeT`@@pjrQQIFTpDEP@P\ndif`@@pj[IEEEZxBH`@pPfy@`\ndid@p@bFbAbDfUfn`BH@H\ndaF@@@RYUih@H@LBCB{bHp\ndmuD@JDBRY^YEZjif@H\ndifH@AAJY}rjY[`bf\ndaE@@@yIe^f`@`@pKBknHB@\ndaxD@@QIeejZ@LLJfyG@\ngC`HAbIMLA``nS@\ngJQHH@aJmj@ppQyL\ndaxL`HS@|Dj}fjX@`\nded@@DiUWeiZCHTZ\\IBDpjXYyG@\ngGU@CPdsmJA`e^S@\ndeU@@@aJWeQfZ@@BHlNSJ[agdP\ndmu@pFTH`h`XaIf[oi`bH`@`\neMPBchLF^P\ngGU@DPdrmTAaecrT\ngN}@DZDfVZj@prqxjp\ngGUHEhOAJmZhCAF|pV@\ndmND@DBdfV]RZUZf@@@LBSagdJ\ndev@@@RfV^nFP@@@@LFPj[ayD`\ndcnH@HCHhheEbprl@@A@@P\ngFu@LPeLftu@XYF|f@\ndiT`@@rnRfUjEnBA``B\ndmLD@@IJ]YVDeZj@B@B\ndcnH@DAIYegzUujAHH@H\ndmL@@DjYUVGi@@@`@LNpjxYWdH\ndcn@`N@HRYeevz]`@@`@B\ndcn`@@rfeJY{]kadHHZR@H\ndevH`LX@aJWY\\HYiZZd@`\ngC`HADIKTAaaMrH\ngNpDADXHRV^jhB\ngJPXHlPDQztAxlP\ngNqHLHaIYzj`H\ngGQDJH`qBSKMPdX\ngNq`@jdvsTpFBsyD\ngGQHJHaIYfdRL\ndaFH@NAIe^f`@`@pXDpjx\neMaDBKpRVB\ngOx@@eJqmUTA`xlZ~P@\neFA@H`bLFE\\`\ngKP@Di\\YZ@phbq@\ngJPhIaxIVmPFDODD\ndcLL`NWPbDfUuZf`@jX@LNALJat\ngN|HEb@aIevj`H\ndaE@@@aJyUnx@@@`\ngJQ`@bdru@XS\\a@\ndayL@DhByIeuji@LBSBkdL\ngJY@DDfzhCCSGbB@\ndmVD@JADf^Uvjjh@`\ngJXdIbHaABRuJA@\ndiDDPNDHc@aIUVfjXHb`\ndefHPFDHc@aIUUiji`bJ\ndiEL@DpByIf^ZjX@`\ngNy`LETeUZZDs@\ngNx`BFdskUHFBVxcP\ndmLH@@Rge]aNFh@J`@`\ndeT@@DjU_k``b`@pFDpjXYyG@\ndmNH`Jd@aJYW_JxZA@j@CB`IaWdH\ndmM@@@[IEDhhZXU@@B@@H\nfHdxA@\neM`AIxLI@\ndeVD@DCDenUFVh@@@pZLinGdR\ndeVD@D@dfVuFVh@@@`\ndif@@@RfU~F``@@pYLJf{b@H\ngOx`FDdrikUTA@\ndeVH@IAJYW~F``H@LFSBinyD`\ndig@@@aDiyWaZ@@B@h\ndco@`LK`BLdTTRRITntpTA@Pe@\ngOx@@drm\\@@A`Pl~Ht\ndaE@`LH@aIYZfijd@`\ngOu@LPeLimMTHf\ngJP`@dfvhCCKD\ngKP@Di\\Vi@pLVOH@\ndiVH@BAIfUInFjjZ@H\ngOpH@EILkW@@@P\ngJPDADFHRYj`LBHcrX\ngJPHAbIKUPFBDyeb\ngNthLXaxYIHTuTA@\ngNxdMVLP~BRkUS@P\ngCd@AH}PFBVyH\ngCd@ADkV@pdxe`\ngGU@MPdllrA@\ngJPLADHHLPdwSBJ`\ndeVH@IAIe]ZZ@Bh@LAALJnFP\ndidD@@iJ[gxZB@@CBdJnGdL\ndaFD@FCdig|jfV`bf\ndaE@@@YIeZn`B@@piLiny@@\ndaEH`Dq`BDfUjyjfPC@`SNyE@\ndieH`Dq`BDfUfnZii@H\ndaFH@BAIf]n``@@phLJfyG@\ngFxhMD@cIHXZmTA`c^P@\ndaFH@BAIf]n``@@phLInyE@\ndifD@N@TfUvfZVZ@H\ndeTL`HS@BLddlRPrm@@@FEYS\\OHd\ndmtLPHS@BEBLddlRVFUh@H@H\neMPBchLD^T\ngJPDADFHRYfaIXXH|f@\ndkOL@Dhpg\\bbTTLVZjjj@H\ndcLL@@gTfYe]aZYjf@H\ngNp`@df]Zi@pvMyF\ndmtHPBBHfHRYeUXXHHh@LFCJxYyB`\ndaE@`LX@aIYZfjZd@`\ngGTHE`DILkU@XDKWdH\ndeWD@Li`QImenEjZi@H\ndcOH@DWPRYyWQej@BP@`\ndg}@@@mJYeU|]Tz@@@H@B\ndcND@LADf^YZUZ`HH@H\ndcnH`FDHaIe[gkev@@`@@H\neMbDBDfp|`\ndcMD@DpBRYgmUujh@@@`\ndmuh@DJaePRUeYxZ`@h@H\ndcNd@DJadMrIJIJE\\MP@U@A@\ngChHH@aJz`LDEqS@\ndeVH`BDHaIf]vzB@h@H\ndeU@@@eIYVvG`BL@H\ndcM@pEtIAICICHiCDedLuP@R@D\ndkmH@JTDeVVutvjh@@@`\ndmw@`Ds`BDf^YyUjB@@B\ngGYHMQHIJmT`XJKbq@\neFhHcA`d\nfoA@R@HHqH@QddebRbrPeV\\m@D@@@A@\nfoA@R@HHqp@QddebbrRTeV\\mA@@@@A@\nfoA@R@HHpx@QddebRfRpiFlm@@@P@A`sEFBT{bMF@\ndaE@@@yIe^f`@`@pILJnHG@\ndmuD`LVD@HrREQIXYV`@`@`\nf`i@@@LdbbbTRVHeZ][uHAD@D@@P\nf`i@@@LdbbRRlRSI\\D{qAHT@D@@XJpRgIZMYwnHPH\neMCBthabHzB\ndk\\H@@RfU}WJxZA@fXBIh\nf`iQ@@AR@drnrjmFJt{t@Du@@@@XR@cAJRM[wHN@\ndmM@@@qJY{WJeVj@B@CCdrfzUy@`\nffc@H@@\\eYtRurJJIQPpjYQYLxkR\\m@qAUUUQ@A@\nfmo@H@@B}GKYkrJJJJJIY[EQJKIbmNlGhv`BH@Rjjj``@`\ndk|@@LdbbbRQKauS`@`@`@H\ndcND`La@BLddJTrRzmT@@@XMaMpkqIt\nfc\x7F@H@@\\ddoU{rJJIQQXsIYJIQYLxhpXcftCDP@EUUUA@A@\nfewAP@@Fbt^QQJJJIRI[RKYr`mA}eMVfjjjjjjjj@B\nfheP`@BJ@GHhhdhlhiVRxHuc^ZejjjYh@LNXIsQk^yAp\ndcmH@EtDie_[JxZA@f`@`\nfhy`a@APP@HtDYICDeLeLerVc^Bl@QP@D@A@\nfoQP@@@VRfYe]Z\\d[S`@@@H@@CCaFBUhsnHyX\ndco@@@PTfyVwiWV`BJD@`\ndig@@@aDkYWaZ@@@LJPj[nI@`\nfbmA``BLbCDubNqLxfQQQQJYRqZGdeVBfZfi``bX@B\ndg|L`HS@BLddlTRbLj]kP@@@@FGXpsiwrE@\nflm`@@@ITssJjnsdk^RFHP@@@@@@@D\nf`i@`@@DRYfyU]`mNmyi`@@B@@LN\\dhug^yCP\ndk\\B@@zSrJJIQPenWV`X`h@LLFUyfLX\nedZPK@@@FAEG@dnaec`X\\bbbbbbIbRbTrbJffW@fefUUssPSDQUUUUUUTT@A@\ndmvL`NaL@HrRRqQZUV``@@`\nfa{@a@K]@DNdLdTtbRLRbLfVuAJtFMVuUUUUURsU@AaxAJLxJRcVBgIZ|QA@\nfew@c@OS@DXB@rFRJQZQYFIQFSK[bcVSghzjjjjjjefj@B\nffsA`@@VkdTRTtRabrTRhqQoARUijjijjjj@CArLTxJRcN}GIZ|`d\nfa{@P@@VkwHhdiheCEdheKQbc^BdkSUUSUUTuP@P\ndmL@@DjYeVdU@@`@@LFpjFUxc\\\ndie@@@aJyexV@`@C@binxdB\ndk^D@IAdfYYwz]MjTHB@CBjE]ObYH\nf`ih`@LRt{pAISK[NzRbmV]UUPHEP@FCHDipWnHPX\nfnki`@DTxIPOMrJIQQKYQQIIEgIJSWIVjVZ`hJJb@B\nf`iQP@DX@IQoYEEDeEeBf\\uYwfjffjj`@`\ndetB@@kirQQIRyS]UUUPAaTXPsqDX\ndieH@BDDfYun``H@LLrny`DH\nfoAP@@@XRf^uu\\Tdij`@fXACU@\ndeUH@JDDinUzZBA@@pIBZ^Ip`\ndeUH@AdDimW[j@BP@pIBknId`\ndig@@@pldTRaifij`B\ndmN@@@rRJKQEneP`@@@CCdrfFUyG@\nfb}@B@@AFQQQJQIHqZIgEV\\FNZjjjjjjj@B\nflu`P@@`T[vQQQQQQQIR|EjKdm@DQ@pP@A@\nf`ip@@@V}dbbbbRrQirQmN@@@@a@@B\nfc\x7F``@D@wdTRrbfdRQfRdrvtxLRufmjjjjjjjjjj@CCjLTYpRmFl{pXiKUoQsWodB@\ndmtL@@x\\bbbTUZViVf@H\nfoA@`@@NRYWVUzLMYiY`@hAAcasAJRTYwLHyh\nffcA`@@EjdrnsLjomF^SGKRsP@@AT@@P\nfdy@P@@^CEIe]YVvhpufef@Bf`@H\ndmtL@@QTfyeQehBA@C@jXYxb\\\ndg\\B@@SSRYYYyVvhJA`@H\ndeUD@DpFRUVTYZZZPcJ\ndeVH@IAJYW~F``H@LJPj[nId`\ndk^@@@RfYU\\]Tzjjjj@LKBDpj[ae]L\ngJQHBLQIVi@`\ndieD`Dqa@HRYVZyjfX@`\ndaDH@@RYWih@H@LLrny`HP\ndg^D@LADf^YVyUj@`j@B\ndifH`FDHaIe[kh@b@B\ngGT@Adekj`LD[rX\ndaD`@@pjRfUi`HH@LH{fHS@\ndeTh@DiiAIYmQfj@@@H\ndeVD@NBTfUuifefX@`\ngGYHL`DIMlu@XHwbF@\ndaxD@@QIeYjZ@LBpf{dT\ndcNL@D@mRY^UyUj`@`@`\ndeVh@LKadDimY[j@B`@`\ndmuD@ITDR[e]xV`@d@LBJf{bHh\ndmuH@DpDf[mYUj`@@CAj[ae^Q`\ngOp@DjWkB@@LMc^ZI`\ndcNH@DCHheEDbnmPT@@F@hUMproHt\ndmuL@DHJUIeeTYZijT@`\ndeVH@IAIfu~Eh@H@LJpj[nQH\ndmuD@DpERYyWFVjjj@H\nfgAp@@@XEdbbTTrLiI`mFBBhea`@pYSdmFmrB@\ndiW@@@pdieZfxPHJH@`\ndiW@@@rdiefa[``hH@ppf{bHH\ndev`@@pjYJYejxY@bJH@pPYxc\\\nfgA@@@LdbbbTRoIBMjp@@@@@@FDlTYpTeZmrA`\ndk|H@@RfU{WJx]N`PJjH@`\ndew@@@rdifYj^DHJ``B\nfgAP@@@BrQQJESITXJQk@PT@A@@X\\IPRmFmrB@\ndg|L@@ildTRbrJQTJtEAEL@D\nfoQH@@@XhiJYmg^YpRcVBHJiXX@H\ngOy@LDeYMjj`H\ndet@@DjYUX^d@@@@CAbinF^Hf@\ndie@@@EIe[ih@J@CBdrnGd@\ndmNH@HCHhheDVzU`@@@@LFCJXYxgJ\ndeVD@LBTeegijfjd@`\ndk^HPNHHbXaIf^UvGS``B@@@`\ndclL@@GTf[VwiWV`BJH@pHjy^Q`\ndaEH@DHDfURijZ`CAF{bI`\ndg}H`LHPBDigUm\\h{ifd@H@LAABe]N~Pp\ndaED@LhDRY[Ifjf@H\ndev@PF@HdHrQQJEEPsTlmPA@\ngJPhLQDIKTpD\ngJPhHPDIUmHFDGDl\ngJXhCL@aJyf@`\ndo~L@M@iRYg^ufzB@jj`@pLLJnFUt{\\\neMACDXaIhH\nfjs`@@@ISLjmvj~JBtjp\\iKP@@@@@D@@A@\ndmtB`HSBCprRSFIJUZh@@@pinFUxfD\nfbm@P@@BBeIfYfUWybm^mEh@b@aF`@B\ndclL@@HTfYm]aWViB@`@ppiWfUv`\nf`qAP@@LDivQQQQQRk[`mFZBX`j`@H\nfoQPB@F\\@DYHheBeLdRdeV]Th@@D@@X\\@RfcV]qBI@\nfhypB@F\\T@HrQQJEJIIXeIZ|FiPBB@`@B\ndg\x7FH@LiPRUYYWESnj@@@@@`\nfdeIb@LLDhwh`BLdTrbRJbTOARTxM@PMUUE@AaHBJRu[pX|`x\nfdeP`@E^@{HhhhhhdeBpV`vbt@Q@LD@@XLyhu`qxaW`\ngOx`FJdlm[MLAyib\ndk^D@NBdef^ukmMfi@B@B\nfb}@P@@HM{HihheDeEdbIVk^BgKP@@PT@D@A@\ndidL@@RdfVwaZii@LHKfMp`\nf`i`@@@ISLrj\x7FSdkZ\\@@@DD@@D\nfhyP@@@VrQQQQIYHrgIFtx@@@BJ@@B\ndg}H@DxDfVyVyO[Z`@@@@H\ndcl@@DjYU_egX@@@@@p[J[agaHRz\ndcl@@DjYU_egX@@@@@pYLJngaLJz\nfdyH@@@P]yJVYU_Ud[Uf`@@Bj@D@T\ndeVH@DAIgeQej@@@LJSFx^IT`\nfnsPR@FJtZw`AFRJQQEJIVYQIXYrUVuUUSLt@E@@XN@cAJ\\TXJ\\ejt\ndg}@@@aJVYU^Svv`@@@@`KCj[ae]yaGX\ndazJ`JaLK`BLddNRuR@P\ndmvh@HJfTLbbbdSZZfji@H\ndaDH@@RVU[f@@@LL`fyfUp\nfhyA`@@HBdrrkN~RdcN|uA@@AP@A@\nfdu@@@DjYee]}faRlGtP@@@@@@@@puS`eZMYw`qyC@\nflm@@@LdbbbbbVNVxhUc^bd@@H@@@@@CCUNRUhsoAbgdF@\nfjcA@@@YDhhhhdeElcqVgVcdm@@@@D@@@@P\nflm@@@DjYeeo\x7FyHQm^bd@@@B@@@@CCUNRUjsoAbgdL@\nfduA`@@TCdbbTVbfbLIPQmN}F``@@BH@@H\ndeT`@@aIRge_aV@B@B\ndif@@@RfU~F``@@pxDpj[nPH\ndmND@BA\\bbbTUZ^EjVih@pSayfTp\ngCaHLHaIYaIXHSbV@\ndcLD@@uIfYoXVfZi`CANWbUH\ndmtD@@eIfUTUjBB@@pDHJfxYxfJ\ngNy@LDeVjj@phQkxi`\ndg^F@D@nt{RUefuFVjjjj@LMSJ[ae]N~Q@\ndae@`LH@aJYTj[jZjPCC@RnHw@\nf`ipB@L\\D@HrQQPsISIDXkSoUL`@@P@A@\ndk]D@BXDR[ee]ntvf``H@H\ndk^D@M@dfYmUXUMjP`R@B\nfgA@`@@VRfU{WqS`mFhD@@h@@`\nfoQA@@@ISJkZ~J\\eYu@aA@P@A@\nfhyh@@@XIQVRJJIIZEKLyISo@TEP@P`@P\nfdeI`@LZCDDeYHhhhlhblbTM[pZjjhDBh@B\nf`iH`@LR]x@diemg]IQVkNjjhDBh@CCdBTxHwnHQh\nfHahAa@\nf`q@`@@^RYWUe^cKN`@f@B@DNVB`HpRbmFm{bNB@\ngGP@DiVV`hipJqoDH\ndeTH@@RYWVf`@j@CC`SBhYyG@\ndmu@@@QInUwaZ@BP@pxDJfz^Ph\nda{@`Dp`BDfUvjXHj`\nfoApP@DXxBQgYEDhihUCPVkMMUUTp@X\\hpRcV]qGM@\ndmM@@@{IEEEEZzU@B@@@LBHYWbTp\ndcmH@LXDee[UnWZ@@@@@`\nfoAQ@@NB@drnkJtYYt@Dp@@BGJ\nf`qQ`@DXAKrSLj}klyYsUUH@T@A@\ndclL@@YTfYetngV`hD`@pJag^Ic@\ngNx`LDdvkUHFAqkyD\ndif@@@RYWZZ@BP@piLJnx`B\nfhi`C@I@dDRFICHiCDeMhdaJ|Fj@H`@@@H\nfhia@@I@RfuuYWYqwj@BX@H@PEP\nf`qh@@@XIQfRJJKZJEJgG^ejj`@`@LADYIVcV]rG`\ngC``AdeZ@pdxe`\ndg|d@BE[CDeenUyT{Zj`PH@H\ndk^h@DX]LDf[e]xUvj`PJ@CB`riOdP\nfoAh`@BBUYpLIKW\\kleFluU@`T@A@\ndmu@@@EIYWYnX@J`@pxLJfe^Q@\ndaE@@@aJyUnX@@@pHj[nIF@\ndaDH`NCDRYWih@H@LBALkbEp\nfdyqP@DXxBUoQrJIQSQYQPtDYsSUUUML@D\ndg~H`NX@aJU{YWBNzVZ@B@B\ndaE@@@aJyUnX@@@pPjxTXw@\ndeVD@J@TinUzZBB@@`\ndcMD@JDErQSQIE]MA@e@A@\nfhipb@LBCA@@cIIBhhddmmNmyZ`@@B@@H\ngOx@@eLmXD@AaUcV\nfoQ@@@LdbbbRQrYHRkN@@`@@@@LIXIPRmFl{dN@\nfoI`@@@IRlrj|DTxjqgUP@@AP@A`jBFBdkQkNy@@\nfdu@@@DjYee]}faRlGtP@@@@@@@@pubgEZMYw`qy@P\nfmoA`@@HWdTrTTtbRLrTfVrk^CdhymUUUUUUUUU@AauFBTYrRmFl{pXyKUoQsV\ndifH@AAIfU[hBA@CB`Jf{bXH\nfoAP@@@TRfyWm^RlzP``@@@CBQJLdkQkNxaG@\ndmuD`IVD@HrRFIKKaV@BP@`\ndeVH@DAIgeQej@@@LBSF{fTL`\nfde@`@@DRYg[VU~BLxLViZjB@X@CCfBTyiwnDeXLP\nffsA`@@HjdkJvmk\\oAZMxJRm@AUUUUEP@XNISdeZMYwhirWbHU@\nfa{A`@@HzdrmkZvmrrUYtTeZuP@UUUUQT@D\nfoAQ@@KN@eLrj\x7Fbthu@@@D`@D\nfhy@b@K^@DISLsJntyIq`t@PPID@@XR@rRmF]xOHB@\nf`ia@@A@RYV{UVgEZ]z@Hf`@@@H\nf`yp@@@PudbTRbRbQqRVeNmzZ`@`@`@B\ndknD`La@BLddJTRVtuj``@@`\ndmuL@DpFEIeY~nZifh@pILe^Qp\nf`i@`@@VRYfYU]`eNMyh@`AB@@LNTyKQg^yAP\ndcnD`La@BLddLTJULnkU@A@A@\nfoQ`@@@ISLrkiJRMiu@A@DP@A@\nfoQP@@@LRfV]uZ\\djsjjfh@`@CAPHHIQfcV]r@`\nfoQPB@HH@DYHddheBe``iV]RuJ@E@@P\ndknD@FBTieenGZB`h`@`\ndieD`JXaCDRYgvzejX@pHLi`\ndcND`La@BLddJbRrzmPP@@XMS\\LkqLt\nfoAQ`@EZ@JVQQQQQIWLTZsT`tAP@D\ndcN@@@Rfye~f`Hb`@pzfxYW^HB@\ndg^D@DCdefVUUMjBHb@B\ndg^H@LAJUyfUSjhBH`@`\ndif@pDBHjHFHrJIQEn`HH@H\ndid`@@pjrQQIFf@`h@LDx^XaL\ndcN`@@pjYJYe}k`Hbj@B\ndid@p@qBqAqDfYun``H@H\ndeTL@@jTef_xVB@`@phj[iyD@\nfoAp@@@P\\eKLjorMjsP@@A@B@KBQJLxJRmNyBp\ndeT`@@pjrQQQUMpEAP@XDUCOHX\ndmLh@DkeAIYe]neZhDB@CCBF^HD@\ndk\\d@Dki@|bbTRQR[aejeZf`B\ndeu@PLX@bPaJWY\\HYiZVh@`\ngJ]@EbDfVdCqX`\ndg}@@@aJnYU^Svz`@@@@CARinFUt{yD@\ndcML@DpFEIeY{kfjYj`B\nfoAQ@@KN@eLrj~bthu@@@E@@D\nfoA`@@@ILkjrmFV]@AL@@@ar`\ndaFD@JCdefZyiZ`B\ndclD@@EJY}erevfVfiBJlFpje]xbL\nfoA@`@@VRfum[V\\Uj`@@JP@B\nfoA@@@DjU][VgKNBAJ@@@PiXZ`cAJLdkQk\\`x\ndaFH@HAIYUnh@@@phHJfxf\\\ndeV@@@rRJFICMPD@@XLfTwCrI@\ndk]@`FD@aJY}e\\kSif`@`@`\nf`qQ@@DX@drlrj~gV|uUUUUP@P\ndg}@@@aJVYU^Svv`@@@@`K@Z[ae]obDX\ndmwD@LxPQIeVvUZjZd@`\ndg\x7FH`NW[@HRfum]XYV`BFiH@pXLIewbHX\ndk^H`Lt@aJWYW\\JUiX@J`@pxD[ad~Q@\ndieH@JxLbTTQkfej`CAFGbPP\ndifH@AAJY}rjYY`bj\nfhiA`@@Hddjrm|jIW`mPAD@@@A`kEFBTZsoAyB@\nfgA`@@@ISLrotyHvk@@@@@@@XLYIVcVy`gA@\nf`q@`@@HRUfYU_Sk^Zh@@@@@LELxJRmFl{wDPp\nfoQAB@E@BDiegUteARlz`RB@H@@`\ngGX`JDdjmVA@\ndmtL@@RTeVUaUj@H@C@jXUydYH\ngJX`LDdju@P\nfjs@@@LdbbbbbbfIQxhSiFmGIP@@@@H@@@@LM\\dkQk^bgI^PB\nffkA@@@YDhhhhdeEmE^JuztPiKV`@@@@@@@@@`\nfbc@@@LdbbbbbbQQsEBMKuhir@@@@QE@D@D\nffsA`@@VkdTRTtRabrTRhqQoArUijjjjjfj@B\ndmM@@@iJYewJEYB`R@CBbinWbCH\nfb}H@@@TMYJYewyUWEAJ]yNYB`RZ`B@A@U@\ndiTH@@RfU|kahDB@C@biaxb\\\nf`iQ@@L\\@dmKZzhiJt{sUTu@D@@P\ndg\x7FH@DWPRUe[mzSmh@ajH@`\ndiW@@@rdiefa[``hH@pQnybDH\ndcm@@@sIEDeBeSKhD@@@@XYpkrdiit\ndieH`HHPBLdRbQrjYi`B\ndklF@@pUttiinwZZjjjh@`\nfbu`R@@zM[pHcHhheDeDeNdpTcPDQSTuQ@A`HFJRlxJ^Hcd\nfde@P@@BUYIfYfVwXKWgQZ@H``b@@H\nflu@P@@BBeIfYfWgvJtztV`BH@dh@@pypVcVbgdO@\ndk\\B@@HSRYfu}aWViB@h@LJKiT~QP\nf`iA@@@IKLl~krRmN|p@@@@@@D\ndmM@@@eJYYvfEP@@`@B\ndeUH@LXDeYgXZjYX@`\ndefJ`JaLFP|LddjRcUTpA`rXUMH\nfbc@@@LdbbbbbbejwMB]Hu`ir@@@@@@@@@FFnBehuoAbgOHJ@\nfdm@@@LdbbbbbfjSCBUIuoAb@@@@AP@`A@\ndo|D@@eIfUwU[hBAif@siXz\ndklD@@SHihdhYNEh@Jh@LFPjWSxgB\nfb}A@@@YDhhhhhdhb^ZMyrXyh@@@@@@@@H\nflmA@@@ILroJ}jQ`iV]EKPPA@@D@@D\ndk~H@NAJUmmRkatzVd@J@B\nfgAa@@N@Re[UUIPTcVijeif`@`\ndcNL@M@iRYg^vzB@j`@pdLJnFUt\ndmtL@@SdfVUyUjA@@B\ngFp@DiTvVhC@qX|Ph\ndif@@@rRJEKaj@@@LBSF{fEP`\nfoAp@@@P\\eKLjorMjsP@@A@B@K@QJLyIVg\\PJp\ndk]D@HDMrJJJJHiaTx@HB@@H\ndeTH`NBHRYWVf`@j@CB`SJgbDH\ngGY@DDfUj`LLc^S@\ndknJ@CCNgTieY_hZjjjx@`\nfhiAP@@XyKRTrlovmA^CUUUP@P@FAhDipTmFl{p^HPp\ndeVH@BAIemQfk@@@LFbfxYy@@\nfhyp@@@NBeSMrnjpRkF|D@@@Q@@A`kENBejsoA@\nfleq@@JBBAIYe_e]nR]z``R@Bh@PyP\ndclD@@{HhheEBtkl@D@@@XTGCJ|sFt\ngFx@@eRqUU@P\ngFx@@eSQUU@P\nfoAQ`@DX@pRSJs|kSegMTs@A@@X\\QaVcV]rG@\ndk^L@DBmRYVuvfeVj@Bh@H\ngNs@IcPdmmMHD\nfoAqb@BTYIS`RBSLnmnsakTmTuUP@P\ndaEH@DHDfYVyje`CAJ{bXp\ndieD@DHFRYf^EjiX@`\nf`q@R@HHpP@QddarbbTQtYwej`@@@@B\ndeTL@@J\\bbbRKBuKS@FBUwDq`\ndcLDPBtHbXcHhhddbppPQR@D\ndid@@DjUZnBBH@LBinGfPf@\ngOp`@tiguif`LEksHl\ndaD@P@qBbDfYvzB@@B\ndk\\H@@rIJJIQDYtz`@@@@H\nfgA@@@LdbbbTVKIBMjp@@@@@@FDlTYrRmFmr@`\ndmNH@BAJ]YVDeZj@B@B\ndk~H@LAJYWmRkatzA@@H@B\nfoQh@@@\\UYvRJJJJsITyKQk@A@PtD@D\ndcmD@DpFRYUmRi]Zijj`CAlxYW^P`\ndg}D@DpNrJJIIPiDqvuTuUT@P\nfgAP@@@\\RfYe_irQmV@@@@@@@pxQ`eZM[dI@\ndctDXJXIAICi@YAYCYCHiMEHyUMUPA@\nfjc@`@@VrJIJZIPqYIJcEF|ENVfjjjjjj`@p|cENBdkQkN}EN^HBT\ndieH@BxDfYUa``P@LBCJGbMp\ndg|H`BBHRYe^uXSnBB@@@@pxLhYT~Xvv\ndaE@@@yIeVf`@`@pHLKnHc@\nffk@`@@LrJJJJHsIIPi\\dkWhaRVmPTP@@@B@@A@\nffk@@@DjYee\x7F^uyHQmNCDmP@@@@@@@@@B\nf`q`@@@YIEBedhdnB]zh@J@@@@`\nfhy@`@@\\RYee~uYrVoAZA``@H@@pyPTmFm{fLNB\ndid@@DkYWaz@@@LF`j[ayB@\nfduP@@@^RiggmUvBdjwgAjjfABJh@B\ndcLDpITJsjsZqIfVYXXBHZ@C@`y]xfB\nfhiA@@@ILkkLktYxL@DtD@@H\\lM@Q`eEZMYwnHxD\nfleA@@@ILklrsoQdTp@SUA@@BGJ\ndaF@@@RYWih@H@LBALkbEp\nfdyQ@@@qAdTTRQTRVTYrshBA`@f@DNT\ndaF@@@RYWih@H@LBCBkbIp\nfbmP@@@IrRJJIQJIHdjmhKhu@DPQAT@@X^qS`iJtZsoAbgH\nfhya@@K@rQQQQIYHrgIFtziAhEJh@B\ndmvD@DCDeUeYUj`@@B\ndmL@@DjUgZFUBBbb@LISBinFUyD@\ndmvH`Nd@aIe]Zf`@i`@piLJex`\\\ndeTDPDpHb@cHiCDdLrp@@A@\ndmvD@JBTin_^F``J@C@dIexe\\\ndaE@@@aJmUnjjh@pYLJf{dP\nfl}A@@@IRlrjkoAENJl[tTuP@@@@T@@D\ndknL@LACRYVuvUZh@J@B\nfhyA`@@B|dsLro~pRkF\\t@P@EP@A@\ndk\\B@@z]rJJIQPqnTv`XbH@H\nflu@P@@\\e{HhheEEcDef\\TZJV`Xb@`H@@pebRkN|FJ^PD\ndaGH@LJ`rQQRyULs@D\ndcnD@B@TfYm]aWViB@`@`\ndcn`@@rawIEDhiUTkhPUCD@XDc\\oHt\nfdei`@LJLzHIZrQQQQIJJkDiYpXuUUTBAP@D\nfhypB@M^B@HRfYU[wDYUg^h@@HF``@`\ndmtL`HS@BLddlRRFUh@H@LJqne^QP\nf`qAP@@\\dYvQQQJJDiQgEZZAbBf`@H\nfj}@P@@\\]eIfVWyU^YrV`iZA``HHZ@@H\nfc\x7F@P@@H]oHidhiiDhjdeIemjpXekM[UUUUUUUUUT@F@tXIQgIJtZsoAbdmV}GM^p\nfj}@p@@BMGI\\bbbbbbbtQqxkWgIZ@Hb@dj@@H\nfdeP`@E^@{HhhhhhdeBpV`vbt@Q@LD@@XLxhu`qxcW`\nf`i`@@@ITrlnzpRkN|DA@@@@@F@lTyKQk^yaBg@\ndmtDpFTH`h`XaIf[oi`bH`@phBXYyB`\nfhiP@@@ARYe_YWYqwhBBh@H@@`\ndifH@DAIVUxU`@@C@biny@`\ndifH`Lx@cIEDcBjiYh@phLJayB@\nfoAaB@LJ@DYHhhmMDbQegULpLA@@P\nfhiQP@DX@IW`yEEDeDdbdsakMUMMUU@AaKAJRMYw`|ah\nfhiPP@DX@IWdfYVU_V\\MYjiijZhDKT\ndk\\D@@[HhhiDjbY]ZBhB`@`\ndeu@@@[IEDhcSCH@@@@FB]OBPjD\ndeVD@FADfygFV``@@pjfxYyB@\ndeT@@LdbRTm\\DDT@FCIc\\L|RJ@\ndidD@@iIYmxVfZPB\ndidD@@iIUexVZZPaF\nf`qA@@@ISZzljscoT@Dp@@B@kCRLDhpTmFl{dO@\ndk]H@NTDimYuX]N`B@@@B\nfhip@@@LCdbbTTrLvUIBt@DAEU@@FAdTYpTeFl{p^Pt\ndmuD@DTArJJJFHUZijX@`\ndaE@`LH@aJ[Unh@@@pHD[nQP\ndie@@@aJVuxZ`@@CAlJfx^HD@\ngG]@DZDfUj`LDkqX`\ndet`@@raRfVZi[``jB@H\nfhy@P@@BTGHhhhhhebfBtzwf`BH@J@@H\ndiU@@@aJVupnFZiV@H\nf`ia`@N@LyIfV]UV\\eiwf`XHjh`@`\nfhyaP@E@EIwlbbbbbRjRXhuoAZdF`Jh`@`\ndiDB`HSBCpRjuVji@H\ngGX`BJdsmLIgAZ|b@\nda{D@DJ`yIeUjZBBh\ngGT`EaTf]j`LD[qK@\ndknL`LaM@HrRRqQHjUV`HJ@B\ndcNH@DCHheEBdnmU@@@FGIeMpkqIt\ngChHLLQIZ`LDEqS@\ndmu`@@pjY\\dTTRak`Hbh@H\ndmt@@DjUiZdHJb@CAd[ae^IA@\nfdyhP@DTxIPCAcdTRbbURTRYrRfjVZjjh@LIHEKUg^CD\ndeUH@JXDeVVxVYjdHQ`\nffsA@@@YEDeMDhTihdZLUiwhyZjjjjjjj`@pbcAJ\\EIVcV\\FJ\\uy@h\nfoQP@@@\\rQQQQQKFdeVLx@B@@@@@`\nfj}P@@@LRfVVWV}qF]GIYB`B@@@@ABE@\ngGPP@cTfyj`LInYF`\ndk]D`LFD@HrRPjJIEatujB`H@LJrnE^YNf\ndg|D@@WIEEEDYcPsmAQDT`A@\ndif@@@Rfufz`@`@piLIax`\\\nflepr@M^SAFB\\\\CprRSEQIZJKZR`qZjjjfVh@LIDhpTeVBgdL@\nfjc@`@@LrJJJJHpiQKKdeZ]{IZ`h`@H@P@@`\nfhe@@@LdbbbbbbcxhQiNlD@@@@H@@@pepTeZMYp^P|\nfHbXAy@\nfoAA@@@ILkjrmFV]@AL@@@arpTAFBTUhug\\QpP\nfj}@p@@\\dYI\\bbbTTVaaRaYIUhyZAbJ@BJ@@H\ndcnH@EAIfV][iuhFAH@LJiag^Ib`\nffs`@@@YIDmeEEEDhlZL]DI\\ujjZ@HH@d@@p}AFBTuYw`qSdmP\nfhyH@@@\\EKIEEEEEMUfRU[pYBb`@B@@H\nfjsA`@@TdeJsZvoZpPsfk^CdmA@P@@@@@@D\nfbc@@@DjYU_eUfByKQ`icd@@B`@@@@@H\nfbc@@@LdbbbbbbcJwEB]Hu`ir@@@@DT`D@D\nfdyp`@FRbABS\\rornIW`mA@d@D@HBh\ndmLd@Dqe@TfUeZzUZjUj@LLrUyB@\ndmM@PBx@c@aJYg\\jeZdHB@B\ndaDHPH@HHHRme[fii@H\ndk^@@@RYeg]ntxB@`@@H\ndcmD@DHERYYUZz]Z``R@CBdHYwdL\ndigH@DJ`RYYrYjfh@`\ndigH@LHPRfU~F``@@`\ngC`HAbIMTA@\ngC`DAbZHRVhCCB[dP\ngCahHlNbNlA`a`\ngChHL@aIZ`LBHl\neMFIDRM``',Jq='fHbTA@\nfH`pA@\nfHfXA@\neO`BNZ``\ngKXHL@aJWFj`LBEcrP\ngBP`Adibj`H\ngCh`LDdsPFDWI`\neMBBHRZCAKd`\neFhHcAaWH@\ngC`DADZHRVhB\neMBCDRZCAKd`\ngGP`@dfuiaMX[F|b@\ngCa@@dkHFBbyL\neFABH`bLD\neMDARZCAgd@\ngJQhHl@cIHUhCBGd@\ngJXDB@bABUmTA`Pl^R@\ngJXHD@aIYj@ppqyH\ngGX`DJdsmRA`enP`\neMhDRZCAKd`\neFDBcAaWH@\ngCh`LHe]PD\ngOx@@drm]UTAaqEcV\ngJPH@DISUPFABqyH\ngJX`DBdru@XI[dH\nfHapA@\ngC`DAb[DRVhB\ndaDD@@yIe^fZVX@psBkdH\ndidH@@RUe^Fh@@@pxHPj[nPH\ndig@@@x\\dTRQi`HF@akBdpj[dB\ndaE@@@yIe^f`@`@piLJny@@\ndaEH@HxDeeVz`@@B\ndeUL@DpFGHhdhYWMTuPA@\ngOx@@drm\\@@A`Qb~IT\ngOx@@drm\\@@A`Qc^IL\ngOx@@drm\\@@A`Pm^Hl\ndaD`@@pjRfUi`HH@LDKnHc@\ndaGH@DK`R[e[fiZ@LLQnyE@\ndmN@@@RYVuiiV@@@@@`\ngOq@@eL~mLlAal[qI`\ndid@p@bBbFbDfYoa`b@@H\neMBCDRZCAGe@\ndifH@NAIYVXZ@H@CA`cBX^Qp\ngFt@AdiWEihCAJ|TH\ngCd@ADkj@ppbyL\ngCaHLLQIZ`LDEqS@\ndetH@@RgYVDFZh@H@H\ndiTH@@RgeXSaj@B@B\ngGQ@@eJuRA`Xm^P`\ndiDB@@[aRVeZjj@H\ngGPdMQbGpRUZi@`\neFA@HoBJD\ndaxBPHQn@HhHrRPzKThA@\ngNqhHl@cIHUEj`LLZ~P@\ngGT@ATiVj`LJHl^R`\neMPBcTH\ngJY@BDfVhCAK\\a@\ngGX`JDdsmTA`l^R`\ngGY@BDf^j`LBHmqF`\ngGT`EaTf]jPLDmrD\ngGQHJLQIUjdB\ndefD@FADf]ZZjj@LJJnF^Q`\ngCi@LDeZ@pTwH`\ndaF`@@pjYJYfn@b@@pPfyG@\ngNt`E`tf]Zi@pVODj\ndiFB@BAFEInuZjd@pILJnQp\ngCiHEAxIVtAaMqA@\ngGU@DPdjmTAaekqP`\ngNt`E`tf]Zj@pJM_I`\ngNxhMV@aI[ji`LBHmrL\ngChhMDOBNtA`enP@\ngC``AdeY@pbyH\ndaDH@@RVU[f@@@LJcB[nQP\ndaE@@@aJyUnh@@@pXHpj[d\\\ndaFH@NAIe^f`@`@piLJny@@\ngFxHL@aJYujj@pHboEb\ndaF@@@RYe[hB@@LJ@j[nQ`\ndaFH@HAIYUnh@@@pXHpj[d\\\neFPBca@\ndaE@@@aJyUnX@@@`\ngCe@E`dkPFBbyL\nfH`XA@\ngGPBADZPLaYAIZjhB\ngJPDADFHRUj`H\ngCd@ADiZDE@\ngJY@DDeZhCCSGbB@\ngGX`LDdsmTA`m^P`\neMBBHRYCAKd`\ndkNF@BAIWSR[YVYjjfX@`\neF`BJFE\\`\nda{@`Dq`BLbbbKUR@P\ndaDD@@yIe^fZVX@pSBxbT\ndaD@@DjWZXHB@C@dpnxdL\ndmNH@HCHhheDVzU`@@@@H\ngFy@LDidviXB\ndmv@`LCDRUVUeZj@@@H\ngNxDLLQxbRjuU@P\ndaDH@@RYVih@H@LBrf{b@`\ndaG@`LK`BDimVz`@@B\ndeUD@DhBRY[TYZ`@@B\ndeUD@DdBRY[TYZ`@@B\ndaEH@LXDeYVzje`C@`RnxdD\neMBCDRYCAGe@\ngGPXHlQxIU[U@XR|VH\nfH`TA@\neFJHqHpP\ndaxD@@QIeyjZBBlLpnyC@\ngGP@Di^VaHxTQF|f@\ndctH@@RgUUZYfY@dRVQX\neFBBDcA@\neMBBHRYCAGe@\ngJQHBHaInZ@`\ngGUHEhOAJmZhB\ndmLH@@RgVUaIVj`@`@`\ngOx@@eR}XP@AaXl[rL\ndcnH@JAIYe[yYu`@@@@H\ndev@@@RfV^nFP@@@@LFPj[iy@`\ndmM@@@WIEDhlZxY@@@`@H\ngOx`DBdwQkUTA@\ndclD@@UIfV][iuhFAH@LJiag^Q`\ngJPDADXHRVj`H\ngJPDAbGDRUj`H\ngJQHBLQIVj@`\ngNt`Bpdf{Zj@pJqoHp\neMdDhzB\nfHehA@\nfHfpAa@\nfHchA@\ngOq@@drm[RtA`Uc^Q`\ngFx`LDdrfmU@XHwdp\neMhDRUB\neFbHbHpP\neMJDBDf`pQyP\ndcLL`NWPbDfUuZf`@jX@H\ngNx`JDdskUPFDwLZp\ngJX`DBdru@XS\\RH\ndiDDPNDHc@aIUVjjX@`\ngChHL@aIVPH\ngJXhEBOAJuj@`\neMJDbDf`pQyP\ndaFD@DAdegfyjj`B\ndmNh@DkaTDfVYVzUZiZi@H\ndaDD@@IIf]nfih@`\ndaG@@@kdiVrX@a@C@hPfyG@\ndcnL`LaA@HrRPjIKTrzmPHD@FEYtkh\neOHBNZP`\ndiVH@BAIfUInFjZi@H\ngJPDADFHRUfaIP\ngNqDLHaqBRjuU@P\ngNpDAbODRUVjhB\ngGX`BDdjmTA`m^JD\ndidHPBBHFHRYgvzB@`@`\ngJXHD@aIUZ@`\ndifH`NDIAIe]ih@J@B\ndaE@@@IIf]n``@@pHLJnHw@\ngJPDADFHRYfaIP\ndiFL@JANRY]fjf@LLJayF@\ndie@@@YJYYhP`X@CBdJnGdL\ngGTHE`DIJkSBZppUxhP\neMbDBDfp`\ndeVH@DAIgeQej@@@LJrfx^Hd`\ndcN`@@pjYJYenk`Hbj@B\ndidh@DJaAIVUxZ`@@B\ndcNL@BA]RYe]VEjVih@`\ngNxhIV@aJUji`H\ndcNL@LANRYygiUjB@`@`\ndid`@@pjrQQIFf@`h@LLKaxbL\ndmvH@NAIfUTUZBB@@p[BiaWdR\ndmuH@DpDfWeYUj`@@CAlInF^Hb`\nfH`PAa@\nfHd`Aa@\neFJHaHhP\nfoA@R@HHpx@QddebRfRpiFlm@@@P@A@\ndmLH@@rJJIQEneX@@@@CB`kaWfXt`\ndmLH@@RYVuiiV@@@@@phJxYybXh\nfgA`@@@ISLrotyHvk@@@@@@@P\ndaFH@HAIYUnh@@@pKBinXD\\\nfduH@@@X]GIEEMDdehUbbkF|FKT@@ACPP@P\nf`iH@@@TdyJYW^U|TxJsiB@jjh`@`\ndmM@@@qJY{WJeVj@B@CCdrfxYyB`\ndcmH@AdDfUvVfWX@IjH@`\ndiV@@@RfU|kahDB@CB`JfGbIp\nffsa@@G@RfYU_yUTVeFCDmZ@@@JZ@@@@H\neMJDBDe`pQyP\nflua@`IJAHFdKRCi@hRDQHQDaHaHFJLuzJX@Jfjja`@`\nfhyA`@@B|dsLsKnqVgVBt@Q@BP@A`s`mFlGbIN@\nfbmAP@@BBgJSLsLoOkEZ]ZKPAD@RU@@D\ndaE@@@yJe~fB@`@`\nfoAa@@D@RYYeUuVLyj@@@@@CBRJLxJRmFxbs`\nfdyQ@@DA@drsJkzjlYsT@@@T`@D\nfhiq@@D^BAIefUWuUcNZ`@@B`@B\ndcNL`LaI@HrREQISBjt@Q@AaVTp{q@T\ndg^L`LaM@HrRPiSISJZuPAD@D\nfgA@@@LdbbbTVKIBMjp@@@@@@D\nf`iA`@@F\\eLsJmmARmku@A@QD@@P\nfnk@P@@FbuIefYZWz^`mA|gMVfjjjjjjj@B\ndmLH@@RfY]raV`XD`@phj[axdj\nfoAA@@@ISLjoxmJMP@@AH@AaKAJLxJRk\\PYp\nfnsAH@@LDjwdu[dTTTTTiRbTVWAZMdkPSDEUUUD@D\nfoIQ@@BJ@dsJsllepQkF\\uKUUTt@D\nfoAaB@G\\ADILkkJ}FFm@AP@P@A@\ndk\\H@@RfY]\x7FJEZA`RXBAh\nf`qQ@@G^@eLrj\x7FpmJMP@@AL@@XJpRcNBdiwnPL\nfdyQ@@OA@eLrj~lbthu@@@EL`@D\nfhyi@@DTxIPLbbTTRqbRYrRoAZiYh@J@@`\nf`i@`@@DRYfyU]`mNmyi`@@B@@H\ndaF@@@RYWifef@H\ndaF@@@RYe[hB@@LB@j[bYp\nfhe@@@LdbbbbbqQXXUiJ|DB@H@@B@@`\nfnk@a@CM@DVdLbbdRQbTQdrvIrRoAKVjjjjjiYj`@p]`eF\\EKQoQSdoD@`\nfa{@c@K]@DHBGRFQRZQIFIQFSKZ`eZCFkZjjjjjiYj`@`\ndet@@DjYUX^d@@@@CAlJnF^Hc@\ndeth@DkiAIeeVxYZiZd@pQAxaT\ndcNL`LaE@HrRPzIJSju@D@AaVDwSqEt\ndk]@@@qJY{W\\jUZh@I`@`\ndmN@@@rRJKQEneP`@@@B\nfew@`@@HrJYJJZQIFYJSKJmxNRcfuUUUUUUUUP@XUQ`eF\\dkQkN|FNRu[t_HS@\ndeVH@DAIgeQej@@@LJSJX^It`\ndaFH@NAJY\x7FJifx@`\neMdHTf`|rP\ngJY`DLTefhB\ndiFH`Bp@cIECDjZj@H\ndeTB@@RirJIPqTsUUT@P\ndmtL@@slbTtTLpfjjj@H\ndmNH@BAIfUmiEX@@@@CC`qnFU@\ndg\\D@@eIVUuWaZ@Bj`@peBinWSodP\nfgA@@@DjYU_VByHu`@@@@@@H\ndclH@@rIQQJH|J{P@@@@D\ndk\x7F@@@p\\dTRbfQjVGSBBheF@H\ndk}@@@iJUmUR[atzjZ@H@C@PPje]OdP\ndk}@@@yJUmmRkatp@h@B@B\nfgAp@@@XheLvnjs`iFlDPT@@@A@\nfgAP@@@LrQQJEIITxJQk@QE@@@@P\ndcm@@@EJYYwietB@P@@H\ndif@@@RYWZZ@BP@pXDpj{dB\ndg}B@DpAV|bbbRbrK]imUMMU@FGXU]IwrM@\ndg|H`FBHRYV~Ukcn@H`@@@pdDrne]OdV\ndg|@@DjU_eZx{BAH@@BJlMaLJfe]N~Qp\ngGQhHjOAJmZhCBGfFh\ngGT`EaTf]j`LDkqX`\nfgA`@@@YHheEhTj\\EHu`@@@@@@H\nf`i@`@@VRYfYU]`eNMyh@`AB@@LFRUhso\\QS`\nffc@p@@DEkM\\bbbfbTtRrLKAF|dkTaQ@DDT@@P\nfoQPB@F\\@DYHheBeLdRdeV]Th@@D@@P\ndcoH@NWPRVUmVFUX@ah`B\nf`i@`@@VRYfYU]`eNMyh@`AB@@LF\\Uhso\\QRP\nfhy``@A@|dsLsKnqVgVBt@Q@BP@A@\nfde@P@@J]{HhhhhhdeBpV`vbt@Q@LD@@XLUjw`qxca`\ndaFL@NAFR[e[fff@H\nf`i@`@@VRYfYU]`eNMyh@`AB@@LN\\dkQg^yB`\ndk\\D@@qIY[mZ{SZZ`@P@`\ndk]L@LxDMIe]eRkSZjjjh@`\ndclD@@IIf]U[evjj@@@H\nfgA`@@@ISKJntXKQk@@@@@@@P\nfduP@@@FRe[Y]mRTDhu`q@@BJ@@H@B\ngGPdMQDGpRUYiDe@\ndeT@@DjUghP`h`@pYL[agfPU@\ndg]D@LlMrIJJIIHj]UAAD@D\nf`iA`@@BUdTTTTTRqXKSk^Z@H`@`@B\nfgAP@@@\\RfYe_irQmV@@@@@@@`\ndmw@`Dq`BDeUeYUi`@@cJ\ndiV@`J@HRfU|kahDB@CB`PfGd\\\nfgAPB@LD@DISJ\x7FJdhsakTuT@D`@P\ngGThEXQxIVsS@P\ndk^@@@RfYU\\]Tz@@@@@LECBinFUOdZ\nfgA@`@@\\RfYe_irQmVh@`@H@@`\ndiT`@@rnRfUjEnejfPB\nf`ih@@@XIQVRJJIIZEIgIJ]xB`j@@P@H\ndcn`@@xYUIYeg[ewiBHH@H\ndmtH@@RYWUih@IhBN\\NALJiWb\\H\nf`qAA@A@bOQBSJ{\\ktYYt@EP@P@A@\nfdy@P@@BtGHhhhheEcliuoM@DD@T@@P\ndeu@@@[IEDhcSCH@@@@FAEpsqDh\ndcm@@@{IEEEDcWSh@PD@@P\nfoAP`@DXAIIfU^uYrsfjjP@`@H\ndet@@DjYUX^d@@@@CBlJngfHp`\ndclh@LX]AIVYWxUuj`PH@H\ngJP`@TeVhCCSGdP\ndaxL@@RdfVVjh@ppj[bYp\ndg~@@@RYfUWd}mh@@@@@pdBinFUwbGX\ndifH@HAIfuxZ`@@B\ndmtD`NTLQIe]Vf`@jP@`\ngOqHF@aJW\\ZVXCCA[be@\ndk^@@@RfYU\\]Tz@@@@@LICBinFUxff\ndk}@@@iJUmUR[atpBJ@@@CChPiWSyE@\ndk^@@@RfYU\\]Tz@@@@@LI@j[aeSxfZ\ndeUH@JdDin_xZB@`@pqB[fUt`\ndetH@@RfUWJzZABH@LBInGbIH\nfoQAB@C@BDifYU^gIVtz`B@DH@@`\ngBQHDHaIejhB\ndaG@`D[`bDfUjZ@B@CA@sbUp\ndidHPACDZHrJIJFn`BH@H\ndknDPLa@BABLddJbTVtujAH@@`\ndayDPLZD@HhHrRESKSPAafTwDC`\ngNxHD@aIUVj`LFEcWrP\ndiEDPLZD@HhHrRESIZZ`B\ndayH@DpDf]Vjh@piBinyF@\ndeV@@@RfyWahBB@CBj[agf@a@\ndmtB`HSBCprRSFIJUZh@@@pInF^YaJ\ndif@`D@HRVU^Ejjh@`\ndmL`@@siRfUmhVxHFJH@`\ngO}@EfDfUkZZPH\ngCe@E`dmHD\ndeUL@DpFgHhdhUWMTuHA@\ndclD@@EJY}erevfVfiBJh\ndk\\H@@RYeg]ntxB@P@@H\ndk\\D@@MJ[Vuvy]h@@B@@`\nfoAaB@KN@DISLjoxmJMP@@AP@A@\ndaFH@HAIYUnh@@@phHpfxe\\\ngF|@AbeLzmS@XKWbQ@\ndcoH@DJ`RUeUVy]ZZ`@@CBdhYwbPh\ndmOH@FePRYVukaf@HZH@`\nfhiA`@@Hddjrm|jIW`mPAD@@@A@\ndmL`@@[aRfV[hYVPb``@`\ndeVD@LADfvUFVjjh@p{BinF^P`\nfhi``@L@PdssLjn[s`mU@@@@@A@\ndk|@`@BDie^urnGShD@@`@H\ndeVH@BAIemQfk@@@H\nf`qa@@O@rQQQIISGBtju@@@QX@A@\nfbc@@@LdbbbbbRJvcEBMIpTqrAD@@@@@@@D\ndcl@`@BLbTTRbOBnt@@@@A`LDwCJ{r@@\nfdu@P@@JMGIEEDhhhhgPPkR\\kUUSUUUUP@P\ngFxHH@aJUqiZ@ppMyD\nflu@P@@B\\EIfYfWgvJuztV`BHBDh@@`\nfHbhA@\ndcl@@DjYeeiGP@HH@@`\nf`i@`@@VRYfYU]`eNMyh@`AB@@LF\\ehso\\QQP\nfhy@@@LdbbbRQaREiQoA`@`@@@@@`\nfhyQ@@CA@dsLjm{Jvg^Bt@@DCA@A@\nf`q@Q@HH]x@QXHrRFIJYIZ\\Ehrp@@Dh@A@\ndet@@DjYUX^d@@@@C@binxRXL`\ndg|@P@bCbDfUnuj{[`B@B@@H\ndidHPBBHFHRYgVzB@`@pHLh^HW@\nfoAaB@GDADILroZ{AFmAA@@P@A@\nfoAaB@GLADILk_J}NFm@DP@P@A@\ndc\\@@DjYVYayYt@@H`PB\ndiDB@@SaR[eVfjBLh\nfoAP@@@NRYWVUzLMZ@B`@`@B\ngOq@@drm]SRA@\ndidh@DqaAIe[kffjPB\ndeV@@@rIQIPmw@AH@XUaMp|bP\ndeVB@LAaeJYyzzjjj@LAaBinF^Q`\ndk\\H@@RYm[Watz`@@@@H\nfoAh@@@XIQRTsNj~sdm@P@QT@@P\ndklB@@F]RYe]ynZYif`B\ndaF@@@RYVih@H@LLBnybXp\ndknD@LADf^YWeVhBB`@`\ndcNH`IDLQIe[^ih@Jj@B\ndmO@@@bTieulhUfhB@@H\nfoQh@@@\\eYvRJJJJUKTxkQk@AA@tD@D\nflmAP@@LUyNQQQQEJQJI|Eish~BtuUUUUU@A@\ndkm@@@QInUu^Eh@IdBIh\ndifH@JAJ[gxZB@@C@dJfxgB\nfoA@@@DkfYU]UcNz`@@@@@`\ndk~@@@rRJIJJEi{UNhJHb`@`\ndeT@@DjY]zXFB@@pYLinGbEH\nfoAA@@@ILkjrmFV]@AL@@@arpdAFBTUhunXxIp\ndcMH@EtLbbbRQc]@PLp@P\ndiTH@@RfU|kahDB@CCBknXcB\ndifH@DAInUxZ`@@CAhJfx^HB@\ndeVD@HCDi^UFZh@@@`\ngNx`DJdssTpFBsyD\nf`a``@M@EdbbRTQrdxj``hjh@C@rBTYpTeZl{wHX@\ndeT@@DjUghP`h`@pFDpj[iy@`\ndcl@@DjYU_egX@@@@@pyLJag^XwJ\ndmtL`HS@BLdaRbReUj@@@H\ndaxD@@QIUUjj@LJpj[nQ@\ndcl@@DjYU_egX@@@@@p{J[ag^XaJ\ndaE@@@aJmUnjjh@`\ngOr@ABtiYufjPH\nfgApB@LLx@HRevUUpPTcViYj@H@@`\ndeth@DkiAIeeVxYZiZd@pqJGdD\nfoA@P@@\\e[HhheEBcF\\UihFHIh@B\nfj}@P@@\\ceIfVW{WVYrVoAZA``@JZ@@H\nfbm@p@@BLEN\\bbbbbbbrNwEZ|zKPADPBU@@D\ndg|D@@iJYYytYvzB@@`@B\ndmMH@HxLbbbTQ[iV@@@@@`\ngNt@@\\dbLjj@pJu_I@\ndk\\H@@Rf_YWJtzYX@HBJX\ndg|@@DjU_eZx{BAH@@BJ\\EaLine]ObEX\ndg|@@DjU_eZx{BAH@@BJ\\EaLJfe]Ob]X\nf`q`@@@IKJzljsco\\@Dp@@B@j\ndevH@IAJYW\\kahDB`@`\ndif@@@RYWZZ@B`@pYLJnGd@\ndk_@@@Y\\dTRbJtiad@@Bh@B\ndk\\D`HP@cIHXheDQgSV@@@@@pZnE]ObAH\ndg\x7FL@Ds`f|bbTTtJuCimMUUS@FCTwJ]|P]@\nfoAq`@DXxBSlbbTTtJVhKQffjjjX@LATYpRmFmrC`\ndcnH@IAJ[Vu[ev`@@@@LFSBE]yC@\ndmM@@@[IEDhhZdU@B@@@LJRnF^Hr`\nf`iP`@B\\@aInVUuQRUiwfjij@H@@`\ndc\\@@DjYVYayYt@@H``B\nfgAa`@N@t[HhheDTdsdeFltC@UQ@A@\ndeTD`HP@cIHXhdLk@P@A`UMp|PI@\ndie@@@GHhhdVz@``@pKBinH@`\ndmN@`E@HRfUwrnF`PI`@`\ndk\\L@@R\\bbTTTU[mMjBhB@B\nfbmAR@HHqhH@cIIKDeMEhdlIUgIrm@DM@`D@A@\nfoAP@@@NrJJIHqIYgCV`HH@H@@`\ndcLL@@ETfYUUnZYjZ@H\ndcn@@@rIQQJH|J{p@@@@D\ndaF@@@RZW[jii@H\ndiVL@DBnRYYbZfZjj`B\nfHgPAa@\ngC`DAbZHRVhB\ngCaHLHaIZPLDErP\neMPBcXLJyH\ngCahHlHROTA@',Kq='daD@@DjUZxHD@@\ndaD@@DjUZxHH@@\ndaD@@DjWzXHB@@\ndaD@P@bBbDfYvzB@@@\ndaD@P@bFbDfUjz@H@@\ndaD@P@bNbDfUzZ@B@@\ndaD@P@qFbDfUjz@H@@\ndaD@P@qNbDfUzZ@B@@\ndaD@P@qNqDfUzZ@B@@\ndaD@`@bDfUzZ@B@@\ndaD@`@qDeeVz`@@@\ndaDB`HSJCprRRsTsUU@@\ndaDD@@IIf]nZYX@@\ndaDD@@IIf]nZZh@@\ndaDD@@IIf]n``@@@\ndaDD@@IIf^fZfh@@\ndaDD@@IIfunZfd@@\ndaDD@@QIeUnZjh@@\ndaDD@@YIeZn`B@@@\ndaDD@@qJYfnjjh@@\ndaDD@@yIe^fZVX@@\ndaDD@@yIe^f`@`@@\ndaDH@@RVU[f@@@@\ndaDH@@RVU[j@@@@\ndaDH@@RYVih@H@@\ndaDH@@RYWih@H@@\ndaDH@@RYe[hB@@@\ndaDH@@Rfu[j@@@@\ndaDH`BBHRYg[jfj@@\ndaDL@@PdfyVyjZP@\ndaD``DJnBHrJIPeUMV`@\ndaE@@@IIf]n``@@@\ndaE@@@YIeZn`B@@@\ndaE@@@YJYVf@``@@\ndaE@@@YJY^f@``@@\ndaE@@@aJyUnh@@@@\ndaE@@@yIe^f`@`@@\ndaE@@@yJYfn@b@@@\ndaED@DHNRYWifif@@\ndaED@DXNRYUifej@@\ndaED@DpFRYVkfjY@@\ndaEH@DpDfYbijj`@\ndaEH`Bk`bDfUzZZY`@\ndaF@@@RVU[n@@@@\ndaF@@@RYWih@H@@\ndaF@@@RYe[hB@@@\ndaF@@@RfYk`H`@@\ndaF@@@Rfu[j@@@@\ndaF@@@rQQHtpDD@@\ndaF@`BBHRYg[hH@@@\ndaF@`FBHRYVkh@`@@\ndaF@`JBHRVU[j@@@@\ndaF@`L@HRf^rjjj@@\ndaF@`L@HRf_rjjj@@\ndaF@`N@HRYWih@H@@\ndaF@`NBHRYWih@H@@\ndaF@`NBPRYWih@H@@\ndaF@`NBlRYWih@H@@\ndaF@`NCDRYWih@H@@\ndaFD@BADfyVyjj`@\ndaFD@F@dfYvzB@@@\ndaFD@LADfVZYjj`@\ndaFH@BAIf]n``@@@\ndaFH@BAIfunZfd@@\ndaFH@DAIeUnZjh@@\ndaFH@FAIeZn`B@@@\ndaFH@HAIYUnfjh@@\ndaFH@HAIYUnh@@@@\ndaFH@HAIYUnjjh@@\ndaFH@JAIYVfiih@@\ndaFH@JAJUtjjjh@@\ndaFH@LAIYfnjZd@@\ndaFH@NAJYfnjjh@@\ndaFHPBxHa@aIeTjijX@@\ndaG@@@kdig|jVj`@\ndaG@@@rdifvxH`@@\ndaG@`Bk`bDfUzZ@B@@\ndaG@`Dq`BDfUjyijP@\ndaGD@Dp`yIeVfZiX@@\ndaGH@Dq`RYVkffi@@\ndaTD@@iJUfDifzjjj@@\ndaTH@@ReYaBinjjj`@\ndadH`DBHRZuSFzjjh@@\ndag@@@adie\\infjf@@\ndax@@Dj~fjh@@\ndax@@DkU^jh@@\ndax@@LddUeUT@@\ndaxB@@QnR[VZY`cD\ndaxB@@QnR[VZY`cH\ndaxB@@QnR[VZi`@\ndaxBPLinBHKpRUZjf`@\ndaxB`HSBCpRjuZj`@\ndaxD@@QIeejZ@@\ndaxD@@QImUifALj`\ndaxD@@QImUjj@@\ndaxD@@YIgYjf@@\ndaxD@@iIUmfjBD`\ndaxD@@iJU^jj@@\ndaxD@@yIUUfYADb`\ndaxD@@yIUUfZADb`\ndaxDPFxLPlQIf^ZjBH`\ndaxDPLHHb@aI[Vfj@@\ndaxDPLh@c`aI[njj@@\ndaxD`LhLQIUnjZ@@\ndaxH@@RUUYj`aH\ndaxH@@RUUjj`@\ndaxH@@RV]Zj`@\ndaxJ`HSBxOCIILTuT`@\ndaxL@@cdiZZjh@@\ndaxLRJ[`bBBNyIUjjV@@\ndaxL`HS@BLddJruT@@\ndaxd@DpnAdf]ffXHr@\nday@@@aJnZjj@@\nday@@@yIYwjk@@\nday@PFx@b@aJUVjj@@\nday@`Dp@aIfYjj@@\ndayD@DpFRYVZi`@\ndayD@DxFR[VZy`@\ndayD@HhNRV[ff`RDh\ndayH@DpDfUVjh@@\ndayH@DpDfYfjh@@\ndayH@DpDf]Vjh@@\ndayL@DpFyIeUjj@@\ndazD@BADf{fjh@@\ndazD@LADf]Vjh@@\ndazD@LADf^VihHJ@\ndazD@LADf^fjh@@\ndazD@NADf{Vfl@@\ndazD`B[`BDe[jjX@@\ndazL@BAFR[nZj`@\ndazL@LANRYuZj`@\ndazL`BaL@HrRPzKUP@@\nda{DAJ[`Qi{dfuVjT@@\ndcL@D@dDdJdFdAdIdEdMdLbdtbaa]UUUT@@\ndcL@X@bBbFbAbEbMbDfYn\x7Fijjjj@@\ndcL@X@bBdFdAdEdMdDfYn\x7Fi`bHh@@\ndcLB@@ImrJJJIH}KTsTp@@\ndcLB@@P]R[e[^eh@b`@@\ndcLB@@Q]R[e]nEh@I`@@\ndcLB@@RUR[fVQuhHF@@@\ndcLB@@R]rJYQIPbkT@Q@@@\ndcLB@@RiRYyVQejjjh@@\ndcLB@@SmR[YgYUjB@`@@\ndcLB@@reRf[ujzjjjh@@\ndcLB@@riRfY~Qfjjjh@@\ndcLD@@IIf]z[hHBj@@\ndcLD@@QIgVUWVj`@@@\ndcLD@@QInUvxVfjf`@\ndcLD@@cIMHhTmJuUUU@@\ndcLD@@eJ[W[[j@Bk@@\ndcLD@@iJ[g]xZB@f@bX\ndcLD@@iJ[g]xZB@i@@\ndcLD@@kIEMDdcttDDT@@\ndcLD@@uIfUe[ffZi`@\ndcLD@@uIfUk[ffZi`@\ndcLD@@uIfUk[hBBj@@\ndcLD@@uIfYoXVfZi`@\ndcLD@@uJ[WU[j@BZ@`X\ndcLDHFDH`haXcXaIf[ozYjYjP@\ndcLF@@Rag\\bbTVTILuSUT@@\ndcLF@@RiWTfV[UnZjjj@@\ndcLH@@RYeZvz@`j`@@\ndcLH`ECDRUYWQff@B`HR@\ndcLJ@@PUuInUgzV`BJ@@\ndcLJB@PUuNR[eY~eijjh@@\ndcLL@@GTf[VwiZ@Hh@@\ndcLL@@STfyWWaZ@Bh@@\ndcLL@@yTee[UaX@bh@@\ndcLL`HS@BLddJfRtjmP@P@@\ndcLLpLK`bEbMbDeeYunjjjj@@\ndcL``Dpn@HRYueUujX@@HJ@\ndcMB@DpIWTfeyYaZjeZ@@\ndcMD@DhARYVUuujP`@HJ@\ndcMD`LE]@HRV[]nEjiih@@\ndcMH@DdLbbtRTHju@DP@@\ndcMH@DpLbbbLRRzuT@@@@\ndcMH`BuPBDf[U{aj@BX@@\ndcML@DpByIf]vifjfj`@\ndcM`@@pjY\\dTTRrM\\ADUP@@\ndcN@`IBHRYWVjZ@Bj`@@\ndcND@BADf{YU]Zj@@@@\ndcND@FCTfUiffZjjj@@\ndcND@HC\\bbdTb\x7FKUUUT@@\ndcND@N@dfY{UnZejj@@\ndcNDPBePbMbDfYmoa`bAh@@\ndcND`MIPbLbbbbQcSAAPp@@\ndcNH@FAIe[mkh@bj@@\ndcNH@IAIeYUifejZ`@\ndcNH`Bp@aJY{UWZZjj`@\ndcNH`BpLQI[VUWZZ`@@@\ndcNJ@NCIWTify^ajjjj@@\ndcNL@LAERYyeiuj@b@@@\ndcNL@LAJRY{eUujh@@@@\ndcNL@LAMrJYQIPcKT@E@@@\ndcNL@M@aRYgWVzB@j`@@\ndcNL@NB]RfyUZziZjh@@\ndcO@@@rTie_ZnBBJh@@\ndcOB@DrPY]RYV]jyjfid@@\ndcOD@Ds`wHheELUPmMTra@p\ndc\\D@@iJYeUqavUhJBbD@@\ndc^H@AAJU[VaJZUijZjT@@\ndcl@@DjYU_egX@@@@@@\ndcl@@LdbbbRkRJ`@PD@@@\ndclD@@QIe[WiiUjP@h@@\ndclD@@UIfVW[iv@@B@@@\ndclD@@UIfV][iv@@B@@@\ndclD@@UIfV^XYv@@B@@@\ndclD@@aJyegzUvjBHH@@\ndclD@@iJYW]rnF``IhBI`\ndclD@@iJYW]rnF``JX@@\ndclD@@iJYW]rnF``Jh@@\ndclD@@uIfV_XYV@@@`@@\ndclD@@wHhheEbprl@@A@@@\ndclD`AtLQIe]miaf@Bjb@@\ndclH@@rQSQJH|J{P@@@@@\ndclL@@STfUmVfeVi@B`@@\ndcl`@@`nReeWZY]@``@@@\ndcl`@@i]RfYyvxU@hJhP@\ndcl`@@jURfYwVxY@hFh`@\ndcl`@@j]RfYuVxY@hJhP@\ndcl`@@seRfUeZEnBAbh`@\ndcll@DpYAmRYV]zzUZije`@\ndcm@@@UJfUWyYvhJBH@@\ndcm@@@YJYYwhUtH@@@@@\ndcm@@@iJYeyravPhJH@@\ndcm@@@qJY~Ureujh@H@@\ndcm@@@uIe[UiiV@@@`@@\ndcm@@@uJYU_rnf`Pbh@@\ndcm@@@{IDeCDdUKh@UUB@@\ndcm@@@{IDeCEDUSh@UUD@@\ndcm@`LH@cIEDUDeeKmLp@P@@\ndcmB@DpNE\\bbbTRMBtmULmP@@\ndcmD@hxEChgiWdfYYungX@@H@@@\ndcmH`LJPBLdTTJTvTJtpAE@@@\ndcmh@HqagPRVVY^E]fjVjPbH\ndcn@@@RVUmVy]x@@@@@\ndcn@@@Re]eRi]@B@@@@\ndcn@@@RfV]zE]B@@@@@\ndcn@@@RfYW\\hUhFBJ@@\ndcn@@@rIQQJH|J{p@@@@@\ndcn@@@rQQQIXyPkPLDT@@\ndcn@PABPJPRYgUVy]``@@@@\ndcn@PBBPZPrJJJIP}J{@D@@@@\ndcn@PBBPvPRYe[VFU`@@H@@\ndcn@PIBPJPRYfuvE]``@@@@\ndcn@PNBHvHrJJIQGMtk@@@P@@\ndcn@PNBPVPRYeeVF]`@@`@@\ndcnD@D@TfUe^fWX@B@@@@\ndcnDABCbgnrJJIQPll{@A@@@@\ndcnH@AAIVYWXUvjP`H@@\ndcnH@BAIfUnZQV@@@`@@\ndcnH@HAIYVu[ev`@@@@@\ndcnH@NAIYeUkiviBBH@@\ndcnHAH@NbTfYYungX@@H@@@\ndcnH`DDHaIe[gkev@@`@@@\ndcnH`NtHaIfV]XYV@@@`@@\ndcnH`NtLQIfV]XYV@@@`@@\ndcnL@LAMRYUUrhYZh@h@@\ndcnL`LaA@HrRPjIKTrzmPHD@@\ndcnh@DkaTLbbTTRK\\nmHTA@@@\ndcnh@DkatDfVYengVjVjd@@\ndco@@@JdeYUUnW^fh@@@@\ndco@@@pdif]wJeVhHB`@@\ndco@`FePbDfUmVnFX@ajH@@\ndct@D@dBdFdNdAdIdEdMdDfYkjZjj`bH\ndctB@@I]rJJJIVMRuLDE@\ndctB@@PYRYU{ViijBBP\ndctB@@PYRYU{ViijBB`\ndctB@@PnR[kfVjjj@@\ndctD@@QIeUUZijhHj@\ndctD@@QImUUZjjh@@\ndctF@@IaWTfYn~jjj`@\ndctF@@rngTen{mjjj`@\ndctH@@RUUUZjjj@@\ndctHxJBPRPrPFPfPVPvPrJPqJQUUUT@@\ndctd@BE]ADf{UYjji`@\ndctl@BE]ANR[mUfjjf@@\ndcu@HLd@a@```PcHheDYhuUMP@@\ndcuD@LTIrQPiIQTuUL@@\ndcuDPBwaBPRPrJZSQEUUMT@@\ndcvDBAAAeIeUzZjjh@@\ndcvJ@BAIWTfuiUjjj`@\ndcvJ@LBngTen{mjjj`@\ndcwLaNWXS@|Eg^rRSIIHmUUJ@@\ndc~D@JCldRTTRLHp|J{PA@@P@@\ndeL@@Di[ernDYZjij@@\ndeL@@DjYeIjGijjjj@@\ndeLD@@eJUUXSaFZjjj`@\ndeT@@DjUYk``R`@@\ndeT@@DjWvifjih@@\ndeT@@DjYUXPbDP@@\ndeT@h@bBbFbAbIbDfYokjYjh@@\ndeTB@@IaRYf{aZfj`@\ndeTB@@KiRYg]nZej`@\ndeTD@@EIYWVy`@h@@\ndeTD@@EIfUVz@`h@@\ndeTD@@EIfUvz@`h@@\ndeTD@@EJYU^f```@@\ndeTD@@QIUeQej@@@@\ndeTD@@YIfUqehD@@@\ndeTD@@eIfWVyiff@@\ndeTD@@gHhhhjpmLtp@@\ndeTD@@iIYe^e```@@\ndeTD@@yIYVvE`BH@@\ndeTDPAdH`haIf]VzB@h@@\ndeTDPBDHbXaIf]VzeiZ@@\ndeTD`AdHaIe[jz@HX@@\ndeTD`AdHaIe]jZ@BX@@\ndeTD`NDHaIfVVfBA`@@\ndeTH@@RUYTYY`@@aH\ndeTH@@RUYTYi`@@aH\ndeTH@@RYVZfZZj`@\ndeTH@@rJJIHmtA@pD]@\ndeTH`IBHrJJJJlLADP@@\ndeTHpIBHJHZHrJJJX|LDPP@@\ndeTL@@FTfYW[ffZX@@\ndeTL@@FTf[eXVjjh@@\ndeTL@@KlbbbtK\\uKU@@\ndeTL@@Qdf]eFVjfd@@\ndeTL@@QdfygFV``@@@\ndeTL@@RTfyWxV`@`@@\ndeTLB@QdyInYqehH@@@\ndeTL`BjPkDf[W[jjjh@@\ndeU@@@[IHhXlLuA@@@@\ndeU@@@aJnUqfhH@@@\ndeU@@@aJueQfj@@@@\ndeU@@@gIHhTmpu@A@@@\ndeU@PBdHchaIf^VFBBH@@\ndeU@`Dp@aIgeQej@@@@\ndeUB@DpFFTfUgkfjYX@@\ndeUD@DdBRY[TYZ`@@@\ndeUD@DhBRY[TYZ`@@@\ndeUD@DxBRY[TYZ`@@@\ndeUD@FdDR[VTYZZV`cH\ndeUD@HDIRVUunfef`RKh\ndeUDAHDIFY{HeEDcptAA@@@\ndeUD`DpY@HRYf|fZjj`@\ndeUH@AdDee][j@B`@@\ndeUH@AdLbbRaUM@AT@@\ndeUH@DDDfyVzV`B@@@\ndeUH@HpDenUFYh@@Ha@\ndeUH@JdDie_xZB@`@@\ndeUHBDxNeIgeQej@@@@\ndeUL@DpFgHhdhUWMTuH@@\ndeV@@@RVUenh@J@@\ndeV@@@RYfUa`Hb@@\ndeV@@@RYyTYj`@@@\ndeV@@@rRJEK\\MP@P@@\ndeV@PNBHFHRYeYi``x@@\ndeV@`ABHRYVvf`@j@@\ndeV@`ABHRYWVf`@j@@\ndeV@`D@HRUYTYV`@@@\ndeV@`LADReyTYk`@@@\ndeVD@DCDenUFVh@@@@\ndeVD@FADfygFV``@@@\ndeVD@I@dfYw[hHB`@@\ndeVD@IADfyWxV`@`@@\ndeVD@LB\\bTTQU]TuT`@\ndeVD@N@dfY{[fiZh@@\ndeVD@NADfyUxVfjh@@\ndeVH@AAIfuneh@`@@\ndeVH@AAJYWnf```@@\ndeVH@BAIemQfk@@@@\ndeVH@DAIgeQej@@@@\ndeVH@FAIe[jz@Hh@@\ndeVH@FAIe[zz@Hh@@\ndeVH@HAIYVVz`@h@@\ndeVH@HAJUmQfj@@@@\ndeVH@IAIeYZYiZf@@\ndeVH@IAIe]ZZ@Bl@@\ndeVH@IAJYW~F``H@@\ndeVH@NAIe]jZ@Bh@@\ndeVH@NCIELeBpt@Q@@@\ndeVHPNHHbXaIf^VFBBH@@\ndeV``Dqad@aIeYzyiji@@\ndeVh@DkadDfV[[ffZh@@\ndeW@@@kdigvriZjh@@\ndeW@PH[`bDbDfueFZjZd@@\ndeWH@DJPRY[TYZ`@@@\ndeWH@DZPR[e_aZ@B@@\nded@@DiUUjjj@@\nded@H@RnbAbIbDjUUjjj@@\ndedB@@PYR[UYifXDrj@\ndedB@@PYR[UYjjX@@\ndedD@@QIeVVjjP@\ndedD@@QImUVjj`@\ndedD@@QInUvjj`@\ndedD@@iIUmVZj`aH\ndedH@@RUUUfihDRY@\ndedH@@RUUUfjhHR@\ndedHXDBHjPZPzPFPfPRYZZjjh@@\ndedJ@@IaeIf[zjj`@\ndedJ@@rneI[nvjj`@\nded`@@iiReUvjjh@@\nded`PDqi@HDHRYuYjjX@@\ndeddABxYAAf^R[kujjT@@\ndeddABxYAAf^R[mUjjT@@\ndee@HHxHaH`XbXaImUjjj`@\ndeeD@BdDR[mUjjh@@\ndeeD@BhIRUUUjZdHJ@\ndeeL@BdDEInufji`@\ndefB@BAAeInufji`@\ndefB@BANeInvvjf`@\ndefD`FFPBDiWnjjf@@\ndefD`FFPBDi]nijf@@\ndefJ@BAFFTfy^Zjf@@\ndefL@LAIRYuUjjh@@\ndeg@@@zTiUWjjj@@\ndet@@DjYUX^d@@@@@\ndet@@DjzYKaejjj`@\ndet`@@siRfUna[``ha@@\ndet`@@siRfUnnFP`bb@@\ndet``Dki@HRYYUnFVjVi@@\ndeu@@@[IEDhcSCH@@@@@\ndeu@@@kIEEDceCLaPD@@\ndeu@`Dp@aIeURhYfh@@@@\ndeu@`LH@cIEEBeeCML@D@@\ndeu`@@rfFTigmkadHHhP@\ndev@@@RYyULFZh@H@@\ndev@@@Re[TjFP@@@@@\ndev@@@RfUWJzZABH@@\ndev@@@rQQJHtpr@@@@@@\ndev@PL@HPHRYUTjFVj@@@@\ndev`@@raeJY{ZXYB@iH@@\ndevh@DJndDfVU[af@`hP@\ndevhADIadFf^R[fUnxVijY@@\ndew@@@pldTTJVTLmP@P@@\ndg\\B@@ESrJJIIIFMsSMSM@@\ndg\\B@@Q[R[VUmgVf@HhBL`\ndg\\B@@RSR[YVwEVh@bh@@\ndg\\B@@SSRY[W[FVh@Ih@@\ndg\\B@@i]RVYUuzVBBJh@@\ndg\\B@@kcrIQQQIGBkADQT@@\ndg\\B@@x{RYeewgV`hHh@@\ndg\\D@@cIMHhThYJuT@EP@@\ndg\\D@@eIfU_Un`HJj`@@\ndg\\D`ILIAIme[\\Uf`BJ`@@\ndg\\H@@RYvYiVvfijjBL`\ndg\\HpL@HcpPHRZUYwNvjhH@@@\ndg\\J@@PUMInUeWiZ@Hj`@@\ndg\\L@@H|bbbftbkBtsMUP@@\ndg\\L@@PdfVUZYmjBBh@@\ndg\\L@@k\\dTtRRLwSPPQT`@@\ndg\\L@@yTimYvfF`BJj@@\ndg\\L`HS@BLddJrTRQvmUP@@@@\ndg\\L`NNpbDfUuijZ@Bjf@@\ndg\\``Dpi@HrJJJIRiG[UT`@@@\ndg\\dPHP{YPlmlldabbrV`suKRwP@@\ndg\\d`LF[a@BLddJbbQvfmPLA@@@\ndg]B@HxEv|bTTRRTlgM@aQ@IBt@\ndg]B@NTBNtfY}[fyjVjf`@\ndg]BbInDp`BB{IICiEDUikUUMR@@\ndg]DPCmaBHJHRYf{uzYjYj}@@\ndg]HPAuPbBbDfYw[fzB@ij@@\ndg]H`AuPbDfY_[fz@`ij@@\ndg]L@DDBmIeye]SZjjjX@@\ndg]L`FM]l@cIEEideCeTsUKT@@\ndg]L`LVDL@cIIKDedhLkP@Tp@@\ndg]L`LnDD@cIIBhhd]ikTC@P@@\ndg^@pLCDDHWDRUeYWtzjjZi@@\ndg^B@BAMoHiieDeBimU@DP@@\ndg^D`EppbDenUuaff@Bf@@\ndg^D`NEPbDfUv{jZ@Bij@@\ndg^H@DCHhhhddYimTE@P@@\ndg^L@NAarJJZJUQ]KSUMU@@\ndg^h@BzSlDf^[WUvjiif`cT\ndg_B@DhpK[RY[U{FVh@Ih@@\ndglB@@Q]RYUumZjjZ`@\ndglD@@QImUUUjjjj@@\ndglD@@uIUUUmfYfj@dRVQP\ndglFPHkivpqLqDen{nzjjjh@@\ndglLpDp`BH|MBLbdRRRQSUUUT@@\ndgl`@@`fRefUWejjj`RIH\ndgl`@@z[ReU]Ujjjj`@\ndgm@hNdHR@a@a``PaJeynzjfZj@@\ndgmBPFE]NpQHQDiUf{jjjfX@@\ndgmL@DdAuInUg]jjZj@@\ndgmL`DVCl@cIICDmihmUUS@@\ndgm`PFy]NpBHBDiUUYjjjfX@@\ndgnD@KADfuUUVjjjh@@\ndgnD`H[`BLdTRbJRUUMUT@@\ndgnDpKXpbBBJBDfzUUijjjX@@\ndgnLPFe]CDpHRUVV[jjif`@\ndgnL`Fe]@HRevV[ijif`@\ndg|@@DjU_eZx{BAH@@BJ`\ndg|@@Dke]eVx{x@H@@@@\ndg|@P@bBbDfYweVx{``H@@@@\ndg|@P@bKbDfUvUz[S`@`@`@@\ndg|A@@cJz]mJ{n{nDpjjjjj`@\ndg|A@@rae]oIEEDieFpr]UUUUP@@\ndg|A@@rfE]oIEEhYiCBJ]UUUUP@@\ndg|B@@KSRYfuVXUmjP`X`@@\ndg|D@@MIfV]UivviP`R@@\ndg|D@@MIfV]UivxH@@`@@\ndg|D@@SHhmEDhcJg[T@@@@@@\ndg|D@@WHhheEhbtg\\D@D@@@@\ndg|D@@mIe]e^ftvefef`@\ndg|D@@mIe]e^ftx@H@H@@\ndg|F@@RaW\\bbTbfbi]WMUUMU@@\ndg|H@@RVUvU[cn`@`@@@@\ndg|H@@RVYYwySn``@@@@@\ndg|H@@RYfUWd}mh@@@@@@\ndg|H@@RfYU[GUn`@@H@@@\ndg|J@@RiWHheHmDTWUsUUMUP@@\ndg|L@@HtfYnUvG[ZdJAH@@\ndg|L@@Pdf{UeRh{Zjh@H@@\ndg|L@@t|bbbbbJcBmm@tCD@@\ndg|d@Dq]@\\bbbbfJSSimUSTs@@\ndg|h@DqnAIe[UTYNvjZ`@@@\ndg}@@@OIEEEDTmWSh@PDT@@@\ndg}@@@aJVYU^Svv`@@@@`H\ndg}@@@aJVYU^Svz`@@@@@\ndg}@@@aJnYU^Svz`@@@@@\ndg}@A@LINTjWvUz[SB@`@`@@\ndg}D@DTBrJISQJHrivu@@@@@@\ndg}D@DhBRY[VuFSmj@@@@@@\ndg}D@LHNRfueWXSn`BB@@@@\ndg}DPDp{@H`HRefUWd}njjjjh@@\ndg}DPFES@H`HRenUVrnFijifX@@\ndg}H@DdDfvUmQT{Z`@@@@@\ndg}H@DlDfyVu~GSZ@@@`@@\ndg}H@DxDf^YVyO[Z`@@@@@\ndg}H@IlLdTrTraS]Jt@EMQ@@\ndg}L`JXiTIAIf]VunNzVjZj`@\ndg}L`JXiTIAIf^vunNzVjZj`@\ndg~@`K@HRfYU_GUN`@@B@@@\ndg~D@DCDenUmQd{Z`@@@@@\ndg~D@EADfufUqT{ZZP`HBL`\ndg~H@DAIge[TYNvh@@@@@\ndg~J@DBaW\\bbTbfbi]WMUUMU@@\ndg~L@B@SRYfU]Z]mjPhD`@@\ndg\x7FD@L[`QIeyUTYNvjZ`@@@\ndg\x7FH@DpprJJJJKJmJfuUUUT@@\ndiD@`@RdjeVjj`@\ndiD@x@bLBBdJdFdNdAdDeZjjj`@\ndiDB@@PaR[UVjj@@\ndiDB@@QaR[VvjY@@\ndiDBB@QaSdfumjj`@\ndiDD@@IIUZYjhHR@\ndiDD@@qJYWjjX@@\ndiDH@@RUUVZjBD`\ndiDH@@RUUZjj@@\ndiDH@@RYWZjj@@\ndiDJHDpnDAbHbaahcIIJiIZe@@\ndiDL@@HTeUUffPQJh\ndiDL@@IdfWUij`@\ndiDL@@PTfuUifPSJh\ndiDL@@idiUVjj`@\ndiDL@@idiU^jj`@\ndiDNPHSB[a@XhXrRPzQZe`@\ndiDdAB[aANd^R[gVje@@\ndiDdAB[aANd^R[mVje@@\ndiED@BDDR[mVjj@@\ndiED@BXDR[mfjj@@\ndiED`Dka@HRU[fjf@@\ndiFB@BANEInvZjX@@\ndiFD@BADf{Ujj`@\ndiFD@BADf{Yjj`@\ndiFD@BADf{ejj`@\ndiFD@D@dfWUjj`@\ndiFD@LADf]Yjj`@\ndiFD`JxPBDivzji`@\ndiFD`JxPBLbdTljjX@@\ndiFD`JxPQDivzji`@\ndiFL@DBarJIQEjj`@\ndiFL@LAJRY{fjj@@\ndiF`@@SNGHhhcFjf@@\ndiFd@DpjDNrJJHqji`@\ndiGDAFxPSiGdf]UjiP@\ndiT@@DiUTjxVjih@@\ndiTH@@ReUxRfjef`@\ndiTH@@RfUXQajiV`@\ndiTH@@RfU|kahDB@@\ndiTL@@idfV^nxVfjh@@\ndiU@@@iJYWrnFP`H@@\ndiV@@@RfU|kadHB@@\ndiVH@BAIfUInFjZi@@\ndiVH@FAIfU[fFjjj@@\ndid@@DjUZnBAH@@\ndid@@DjU^nBAHBJ`\ndid@@DjUfaBB`@@\ndid@@DjWffB@h@@\ndid@@DjYUaBHP@@\ndid@@DjYmaBH`@@\ndid@@DjY}nBHH@@\ndid@@LdbRQk``R@@\ndid@@LdbRQk``b@@\ndid@@LdbbQxXF@@@\ndid@H@bBbFbAbDfYnnifj@@\ndid@p@bBbFbDfYoa`b@@@\ndidD@@EIYW[f@B@@\ndidD@@EIe]ih@J@@\ndidD@@EIfU[ffiP@\ndidD@@EIfU[hBB@@\ndidD@@GHhhdVyifX@@\ndidD@@IIf][hHB@@\ndidD@@QInUxV`@@@\ndidD@@iIYgxVB@@@\ndidD@@yIYVXV@H@@\ndidD`BDHaIf_[hHB@@\ndidD`BxHaIf^XZVfP@\ndidD`FDHcHhdhZzZVX@@\ndidH@@RUe^Fh@@@@\ndidH@@RYm^Eh@@@@\ndidH@@RYm^Fh@@@@\ndidH@@Rfu^FfiX@@\ndidHPBBHFHRYgVzB@`@@\ndidHPBBHFHRYgvzB@`@@\ndidH`ABHrJJIEn`HH@@\ndidH`NADRYWZZ@B`@@\ndidL@@HTfYun``H@@\ndidL@@IdfYnfZfj@@\ndidL@@IdfYoaZYf@@\ndidL@@IdfYoa`b@@@\ndidL@@KdfYynZej@@\ndidL@@QdfumnZZjBL`\ndidL@@XTfUfn`BH@@\ndidL@@`TiYeLjjj@@\ndidL@@rldTTqkjjj`@\ndidh@DKaAInV[fiZ`@\ndidh@DqaAIe[kffjP@\ndidh`DJaxHcHheCBjfid@@\ndie@@@EIYW[n@B@@\ndie@@@EIfW[hBB@@\ndie@@@GHhhdVz@``@@\ndie@@@aJvUxZ`@@@\ndie@@@aJyW[j@B@@\ndie@@@iJ[gxZB@@@\ndie@@@yIYVX^@H@@\ndie@`BDHaIf_[hHB@@\ndie@`FDHcHhdhZz@H`@@\ndieD@DpFrJIJFnZii@@\ndieD@JXDR[Y^EjZd@@\ndieDAHxAzQyIYW[j@B@@\ndieD`JXaBPRYgvzejX@@\ndieD`LIN@HRZufFZid@@\ndieH@DDDfyWaZ@@@@\ndieH@DXDfyWaZ@@@@\ndieH@HHDeYWaf@@BH`\ndie`@@pjX\\dTTUk`Hb@@\ndif@@@RUe^Fh@@@@\ndif@@@RUe^Gh@@@@\ndif@@@RYVzZ@B`@@\ndif@@@RYevz@``@@\ndif@@@RYfVF@b@@@\ndif@@@RfU~F``@@@\ndif@@@RfWzXBB`@@\ndif@@@rRJEKaj@@@@\ndif@PABHJHRYgVzB@`@@\ndif@PBBHFHRYgVzB@`@@\ndif@PBBPFPRYgVzB@`@@\ndif@`ABHRYevz@``@@\ndif@`ABPRYVzZ@B`@@\ndif@`ACDRYVzZ@B`@@\ndif@`BBHRYgfzB@`@@\ndif@`D@HRUe^EX@@@@\ndif@`FBHrJIJFn`BH@@\ndif@`L@HRf_\\jjjh@@\ndif@`L@HRf_ljjjh@@\ndifD@BADfyWaZ@@@@\ndifD@BADfyWaZjj@@\ndifD@BADfzwaZjj@@\ndifD@DCdfYTfZiZ@@\ndifD@FADfW[aZjj@@\ndifD@J@TiWTjjjj@@\ndifH@AAIfuxV`@@@\ndifH@AAJ[W[j@B@@\ndifH@BAIfuxV`@@@\ndifH@DAIVUxU`@@@\ndifH@DAIVUxUjj`@\ndifH@DAInUxV`@@@\ndifH@FAIfuxV`@@@\ndifH@HAIYV[j@B@@\ndifH@JAIYekjZf`@\ndifH@JAIYgxZB@@@\ndifH@JAJ[gxZB@@@\ndifH@NAIYVXZ@H@@\ndifH@NAIe]ih@J@@\ndifH@NAJ[VXZ@H@@\ndifL@BANR[efEjjh@@\ndig@@@HTfYun``H@@\ndig@@@h\\dTRQriZj`@\ndig@@@idiggJejj@@\ndig@@@kdigwJejj@@\ndig@@@xTjU^njjj@@\ndig@`D[`bDfUff`@h@@\ndig@`Hi`bDeVWajfi@@\ndigD@DpP[HhdhZyjfd@@\ndigH@DK`R[e^Eh@@@@\ndigH@LxPRVUfz`@`@@\ndigH`BiaBPrJJIEniZf@@\ndkL@`@BDeUUUjjjj@@\ndkLB`HSBCprRSEQIYjjjh@@\ndkLD@@iIUmUVZjj`aH\ndkLF`DkaTpBLbdtVafjjji@@\ndkLH@@ReUUVjjjh@@\ndkLJ@@imMJUUUZjjj`@\ndkLLHLi`bB|MbCbDeVUUjjjj@@\ndkLL`HS@|DjUYfjjjj@@\ndkLP@@anF]OIDhhjeIjYi`RKh\ndkLd@BgSADf{UVZjjf@@\ndkLdABWSAMd~rJZYRHijjjT@@\ndkML@DdAuInUgVjij`@\ndkND@D@dfWUUZjjj@@\ndkNF@BAIWSR[YVYjjfX@@\ndkNH@DAIeuUVjjj`@\ndkNJ`EaLdp|LddrQTrZjij@@\ndkN`@@imMJUUUZjjj`@\ndk\\@@LdbbbRuZQT@B@d@@\ndk\\@`@bDfUvUimN@B@@@@\ndk\\@`@bDfYYwZ]NB@@@@@\ndk\\@`@dDfUnukmN@H@@@@\ndk\\@`@dDfUvUimN@B@@@@\ndk\\@`@dDfYYwZ]NB@@@@@\ndk\\@`@qDfYYwZ]NB@@@@@\ndk\\B@@P]RYVyznTviBBH@@\ndk\\D@@wHhhhhbfESZAhD`@@\ndk\\D`ALHaIe[UzZU`@iiH@@\ndk\\H@@RYeg]itxH@@@@@\ndk\\H@@RZYZ\\SgZifjf@@\ndk\\L@@KTfYmUXUMjP`R@@\ndk\\L@@X|bbbbTa[nEijjf`@\ndk\\L@@htie]urnF``Ij@bX\ndk\\L@@htie]urnF``JZ@aX\ndk\\L@@iTie]{rnF``Jf@@\ndk\\LPDp`BH|LbdRRaRX]Mjj@B@@\ndk\\`@@iSRfYyUnEPJBja@@\ndk\\`@@j]RfYuUnFPJAjb@@\ndk\\`@@kSrQQQSQEngPJBja@@\ndk\\`@@p]rQQJJIFiXP`jJa@@\ndk\\`@@y]RVUmUae^@HZb@@\ndk\\`@@{SRVUnUag^@Hja@@\ndk\\d@Dq]@\\bbbbfJZ]MjjZe`@\ndk\\d@DrUCdfY[uXUujijY`@\ndk\\d@DsmB\\bbbTrQXUujjUj`@\ndk\\d@LHYCdimY^X]N`BJh`@\ndk\\h@DrfAIgeU^GSZfh@@@@\ndk]@@@iJYez\\hSdJBh`@@\ndk]@@@qJYweZ[SB@`@@@@\ndk]@pAt@a@``aIeoeRkSZf`@`@@\ndk]D@DhBRY[UWatvjZ@@@@\ndk]H@ALLbbRfTJiaf@Bfh`@\ndk]H@BtDfYmUXUMjP`R@@\ndk]H@DHDfYyWImMjhHA@@\ndk]H@DTLbTTRaR[mMjj@B@@\ndk]`PLJnTpBIBDiguermNfZijP@\ndk]h@HKiTpRYmfTYtz`HIb@@\ndk]h@HKiTpRgefTYtzZiZiBD`\ndk^@@@RYV{Vntx@`@@@@\ndk^@@@RYeg]itxH@@@@@\ndk^@@@ReUVZ]TziZ@@ABZP\ndk^@@@RfYU\\]Tz@@@@@@\ndk^@@@Ri_YVftzjX@H@@\ndk^@pFBPJPFPrJJJZIGitxHB@@@@\ndk^@pIBPJPvPRYfuUaTxH@B@@@\ndk^D@BADfyVux]Mh@@@@@\ndk^D@BCTfYmUXUMjP`R@@\ndk^D@E@dfYuY[eNB@jh`@\ndk^D@IADfvYWz]MjdHB@@\ndk^D@MBdie]vrnF``Jj@@\ndk^H@DAIVUm^GSV@@@@@@\ndk^HPJHJqjqIf[u^GS`b@@@@@\ndk^HPNHJrXaIf^UvGS``B@@@@\ndk^HPNHLR\\QIf^UvGS``B@@@@\ndk^H`BPHaInUm^GSj@@@@@@\ndk^d@DkaTCrJIQQHynWVjVZZ@@\ndk^d@LzULDRZUnTYtvjiZi@@\ndk^h@DkatLbbTTTU[iujejjP@\ndk_@@@Y\\daRRRTaMNjjjj`@\ndk_@@@pdif]uriUjB@j@@\ndk_@@@pldTTtRJrmMjfjU`@\ndk_@@@x|dTTTRJkit@HBH@@\ndk_H@NWPRfumUaeZ@HZb@@\ndklB@@F]RYe]ynZYif`@\ndklB@@P]R[e[_iZ@Hj@@\ndklB@@QSrJYJIJF]ZX@b@cH\ndklB@@QmR[fUxUZBBF@@\ndklB@@RURYYV|jZfij`@\ndklB@@SSR[VUVUZZfjPcH\ndklD@@QIe]URijjjj@@\ndklD@@QIe]e]Mjz@@@@\ndklD@@QIgfU}MjhD@@@\ndklD@@[Hhdhbdjz@Hjh@@\ndklD@@eJ[Vvfz`@jh@@\ndklD@@kIEMDdbnf``bd@@\ndklD`AtIAIVUm~fh@bh@@\ndklF@@PittfzugEVjjjh@@\ndklL@@PtfyWUxV`@j`@@\ndkl`@@kaRe[vTf@HZj@ah\ndkm@`MLHcHhhdcDfz@`Z|@@\ndkmD@DHCRYvUvUZh@J@@\ndkmD@DTCRUfWtYV@`e@@\ndkmD@DdCrIJJIPxUV@bE@@\ndkmH`FVPbDfUonkh@bZ`@@\ndkmH`NVPbDfUunih@JZ`@@\ndkmL@DpAuIfYeVejjZZ@@\ndkmL@DpIOHhhdddVyjiVfBN`\ndkn@@@RYfVya`Hbj@@\ndkn@pLCDDHWDRUeYWSjjijP@\ndknD@BADfyWUxV`@j`@@\ndknD@D@tf^UeFVh@J`@@\ndknD@D@|bbtTTUgVhBH`@@\ndknD@LALbbRbaRtvjh@@@@\ndknDPAWPbFbLbbRbJfkh@bf`@@\ndknD`ItpdDfUunih@Ji`@@\ndknD`ItpdLbbRaTtih@Ji`@@\ndknD`ItpdLbbbRLt[hBBi`@@\ndknH`It@aJYWunF``JX@@\ndknL@BACR[me]]Zj@B@@\ndknL@D@SRUVYu]ZXHB@cH\ndknL@D@]R[e[_iZ@Hj@@\ndknL@FASrJJIIJT]ZBBb@@\ndknL@M@iRYg^un``Jj@@\ndkoD@DJPwHheLeCAej@BX@@\ndkoH@DHpRYYUfSZ``h@@\ndk|@`@BLdTRbtQJagShD@@`@@\ndk|D@@QIeeUrYngVjijjP@\ndk|H@@RYegYngUN@@BHH@@\ndk|H@@RfYU_JGUN`@@B@@@\ndk}@@@iJUmUR[atpBJ@@@@\ndk~@@@RYkUdUgYNjjjjh@@\ndk~@@@RfYU_JGUN`@@B@@@\ndk\x7F@@@p\\dTRbfQjVGSBBheF@@\ndmL@@DjYUVGi@@@`@@\ndmL@@DjYeVdU@@`@@@\ndmL@`@bLbbbTQ[iV@@@@@@\ndmLB@@IerJJJZEnxVijf`@\ndmLB@@RURYUVJaejVjh@@\ndmLD@@QIe[VfeVi@B@@\ndmLD@@UJYWwJxZA@j@@\ndmLH@@RYVuiiV@@@@@@\ndmLH@@RYiYKnUjjjh@@\ndmLH@@RYiYKnVjjjh@@\ndmLH@@RffeZVVjjjh@@\ndmLH@@Rge[aNFh@JP@@\ndmLH@@rJJIQEneX@@@@@\ndmLL@@yTee[VFUX@bb@@\ndmL`@@pURfV[jVDHJbD@@\ndmL`@@siRfUmhVxHFJH@@\ndmL``Dke@HRYYeXYUjfVh@@\ndmLd@Dqe@TfUeZzUZZeZ@@\ndmLh@DkeAIefUaeVjYZ`@\ndmM@@@UJYY~nFP@@H@@\ndmM@@@eJYY^fEP@@`@@\ndmM@@@iJYewJEYB`b@@\ndmM@@@qJYm^neP`@@@@\ndmM@@@{IDeCDbjU@Bjb@@\ndmM@@@{IEEEEZzU@B@@@@\ndmMB@DpNETfYYZEiZjYZ@@\ndmMH@DHDfVUzzUZ`PH@@\ndmMHAHx@zIrJJIQEneX@@@@@\ndmMh@DkaePRYYe[iUifjd@@\ndmN@@@RYVuiiV@@@@@@\ndmN@@@RfVWiaT@@H@@@\ndmN@@@RfV]iaT@@H@@@\ndmN@@@RfV^kad@@B@@@\ndmN@@@RfV_kad@@B@@@\ndmN@@@RfWUiiTB@@@@@\ndmN@@@RfYWraV`XH`@@\ndmN@@@rQQJJFfEP@@`@@\ndmN@@@rQQQQVneP@`@@@\ndmN@`BBHRYe[ZQV@@@@@@\ndmN@`BBPRYUUiiV@@@@@@\ndmN@`BBPRYe[ZQV@@@@@@\ndmN@`DBHRYVuiiV@@@@@@\ndmN@`FBPRYVukiV@@@@@@\ndmN@`NBPRYegXYV@@@@@@\ndmN@`NBPrJJIQEneX@@@@@\ndmN@pN@H`HPHrRPqIZneUhDB@@\ndmND@BA\\bbbReInFjZjd@@\ndmNH@BAIfUmiEX@@@@@\ndmNH@DAIe[VfeX@@@@@\ndmNH@FAIe[VneX@@@@@\ndmNH@NAIYe^neZdHB@@\ndmNH@NAIfV]aeX@@@@@\ndmNH@NCHhheDVzU`@@@@@\ndmNHAH@Nb\\bbbTQ[iV@@@@@@\ndmN``DkaT@aIefUneVjVjP@\ndmNd@DqiTARYVUkiUijUh@@\ndmNh@DkaTDfVYVzUZiZi@@\ndmO@@@pdif]|jUZ``H@@\ndmT@p@AHbDbLddJRRjjj`@\ndmT@p@bBBEbDeYyzjjh@@\ndmTB`HSBCpRjuUZjj`@\ndmTD@@QIeU]jjZ@@\ndmTD@@iIUmUfjjBD`\ndmTDPDpOBICIHUHiIjjh@@\ndmTD`HdDiJiUZjjf@@\ndmTH@@RUgVYjf`aD\ndmTH@@rJQHiRjZZ@@\ndmU@HLT@a@c`bPaIggYjjf@@\ndmU@pLD@a@c`cHheEKFjfh@@\ndmUL`LZiT@aI[mZjZf@@\ndmVD@D@dfWUVjjh@@\ndmVDbNp`qNeI[Ujfjj@@\ndmVL@LAERYuUZjj`@\ndmWLaNeXS@|Ie^Rju_ZjiP@\ndmlH@@RfVvaFGajjjj`@\ndmn@@@RZUYX^gejjjj`@\ndmnH@JAJUe^D[iV`@@@@@\ndmt@H@bAdIdEdDfUvjZ@Bj@@\ndmt@h@bJbNbIbEbDfV[vFZjj`@\ndmtB@@PUrJZIJGiZ@H`@@\ndmtB@@QeR[f]FV``H@@\ndmtD@@QIe]TjZjjh@@\ndmtD@@QIgYVUZh@@@@\ndmtD@@UIfuwaZ@B`@@\ndmtD@@WHhidh^Eh@J@@\ndmtD@@aJye[ahBB`@@\ndmtD@@iJ[g_ahHBP@@\ndmtDPAdH`haIf]Un``J`@@\ndmtDPBTHbXaIf^un``JP@@\ndmtDPITICiAIfVYa`Ha`@@\ndmtDPNDHaXcHhheEVfBAb@@\ndmtH@@RVUv[fZjZ@@\ndmtH@@RVUv[jZjZ@@\ndmtH@@RYWZih@Jh@@\ndmtH@@RYeY[hBBh@@\ndmtH@@RYe[[hBBh@@\ndmtH@@RYeeZVjjj@@\ndmtH@@RYeeZZjjj@@\ndmtH@@RYe~[ffjZ@@\ndmtH@@RYvUeVf@@BL`\ndmtH@@RZWv[jjjZ@@\ndmtH@@Rfuu[j@BXBA`\ndmtH@@Rfuu[j@Bd@@\ndmtH@@Rfuu[j@Bh@@\ndmtH@@rIQQQWiXBH`@@\ndmtHPNBHHHRUUeeZjii@@\ndmtH`IBHRUe]xY`@hBH`\ndmtH`JBHRYVueYi@@BHP\ndmtL@@ETfYUVyifi`@\ndmtL@@F\\bbfbQXVjjj@@\ndmtL@@QTfyW^Eh@I@@\ndmtL@@QTfyW^Eh@J@@\ndmtL@@QdfygQehHB@@\ndmtL@@SdfVyyUjB@@@\ndmtL@@yTee[VE`BJ@@\ndmu@@@GHhhhhvf@bb@@\ndmu@@@aJyyVUjh@@@@\ndmuB@DpAeTfUejyjiV`@\ndmuB@DpFFTfUgjyjfj`@\ndmuB@DpFe\\bbRaRkfjZi@@\ndmuD@ITDR[e]xV`@h@@\ndmuDAHDEFU{HeEDh^F`HJ@@\ndmuDPHhfBHFHRYf~kjfZj@@\ndmuD`LVD@HrRRqIXYV`@`@@\ndmuD`NEeCDRUf_FYiZfBH`\ndmuH@DHDf]eYUj`@@@\ndmuH@DTDf^Uqej@B@@\ndmuH@DpDfUmYUj`@@@\ndmuH`BhPbDfUmzZZZ[`@\ndmuH`BhPdDfUmzZZZ[`@\ndmuH`BhPkDfUmzZZZ[`@\ndmuH`BhPqDfUmzZZZ[`@\ndmuH`NVPbDfUujZ@Bf@@\ndmuL@DpByIf]zfZjZh@@\ndmuL`DxiTHaIevxjfjih@@\ndmv@@@RYfVXXBHh@@\ndmv@@@RYfWXXBHh@@\ndmv@@@rJQFIJUjh@@@@\ndmv@@@rRJIIFUjB`@@@\ndmv@HDBHFPfPVPRYWZih@Jh@@\ndmv@`D@HRUVUeUj@@@@\ndmv@pEBPRPrPrJPqIXYj`@`@@\ndmvD@DATf^Uqej@B@@\ndmvD@DATfvUqej@B@@\ndmvD@DCdf^YyUjB@@@\ndmvD@EADfvUqej@B@@\ndmvDPBS@BHBLddjbReUjjj@@\ndmvD`La@BLddlRVFUh@H@@\ndmvD`La@BLddlTReUhB@@@\ndmvD`NePBDig{ljfZf`@\ndmvH@DAIgfVUZ`H@@@\ndmvH@FAIe[Vn`BJ`@@\ndmvH@FCHhhhdYUhJ@@@\ndmvH@IAIeYVfZVih@@\ndmvH`NdHaIe]Zf`@i`@@\ndmvL`IaL@HrRRqIPUV`B@@@\ndmw@`Dp`BDf]eYUj`@@@\ndmwD@ByPQInvVUZjeh@@\ndmwH@Dp`RYvUeVj@@@@\ndnDH@@ReVDijiZ@@\ndnDH@@ReVDijij@@\ndnDH@@ReVDijji@@\ndnDH@@ReVDijjj@@\ndndD@@iJUhPfijjj`@\ndo\\BPHP{@HtHrRPzIKJQZjje`@\ndo\\F@@Savtf^UVYjjjY`@\ndo\\NPHQevw@DHDrRYJKJHqjfii`@\ndo]B@L\\Dv|bbtTTUTZjjjX@@\ndo]LpN_JXH`o@xaIUkUWjfjjh@@\ndo^B@NAK]ImgV^ZjjjX@@\ndo^D@D@dfWUUUjjjj`@\ndo^J`CaLKP|LddqbbTlZijfhHi@\ndo^LpGZ[BHHHhHrJZQIIHrZjji`@\ndo^LpNWSCDpHJPRV]efyjjfZBH`\ndo|D@@QIgfUYuvj``h@@\ndo|D@@QIn[WVeVj@Bj@@\ndo|J@@S[_HheDdeDYMjBBb`@@\ndo|L@@RtfUVuwSZjp@h@@\ndo|L@@RtfvYWwSZjA@h@@\ndo|L@@iTinU]_ihHHjXBCP\ndo|L@@iTinU]_ihHHjXBC`\ndo|L@@kldTtTRQrEZBHbi@@\ndo|LPDp`BH|LbdRRbRQmvjhJ@@@\ndo}D@DTIrJZQQQIVwZjAh@@@\ndo}DPAuSBHJHRYg]nVzB@ij`@@\ndo}D`ILmBHrJJIIESQa``bfh@@\ndo}H@DHLbbTbjbRmvjj`@@@\ndo}H@DpLbbbLRVTWVj`@j@@\ndo}L`LvDl@cIIBidleIUZ`@i`@@\ndo~B@BAC_HiiMhdh]mjj@H`@@\ndo~B@BAK]InVuf[fjjif`@\ndo~B@I@moHhhdddcFEjVjYh@@\ndo~B@JAC_HheBhdh]mjj@H`@@\ndo~D@DAlbbbRQbRmvjj`@@@\ndo~D`MNpqLbbbRLRa[hBAkZ@cd\ndo~J`IaLNpaLddlRTrTEVh@bf@@\ndo~J`LXUMpBLdTVafffrjjjjj`@\ndo~L@LANRYye[iujBBJ`@@\ndo\x7FD@Lk`QIgff{yVjZjjP@\ndo\x7FH@DxPRYYe}ZyiZjjh@@\neFA@lhbL@\neFAAD`bJ@\neFAADhbL@\neFABDhbL@\neFABHhbL@\neFACDlRL@\neFB@Hc@@\neFB@Jc@@\neFBBHc@@\neFBBlc@@\neFBCDc@@\neFDBb`@\neFDBc@@\neFHBL@\neFJHaHh@\neFJHbHh@\neFJhUKhh@\neFPBc@@\neF`BL@\neFbHbHp@\neM@H~@\neM@ghz@\neMA@HPaIT@\neMABD`bM`@\neMABlZqIh@\neMACD\\QIh@\neMB@Jch@\neMBBHRZ@\neMBBPRY@\neMBBlRZ@\neMDARV@\neMFI@bMP@\neMHAIX@\neMHAIh@\neMJHaHv@\neMLRRWhv@\neMPBch@\neM`BN`@\neMhDRV@\neMhHRY@\neMhHRZ@\neMhHch@\neMpBXv@\neO@Hyj@\neOBBHcfh@\neOBCDcfh@\neOHBNZ`@\neOPBcfX@\neO`BNZ`@\nfH`T@@\nfH`X@@\nfH`p@@\nfHap@@\nfHbT@@\nfHbh@@\nfHc`@@\nfHcp@@\nfHdH@@\nfHdP@@\nfHdT@@\nfHd`@@\nfHe@@@\nfHf@@@\nfHfp@@\nfHg`@@\nfHgp@@\nfHpXT@\nfHppT@\nfHpp\\@\nfHppt@\nfHrpT@\nfHtHh@\nfHt`x@\nfI@@\nf`a@`@@FrJJIJQJrLxtsUSU@@@\nf`a@`@@ZrJJJQIPfJltEPTT@@@\nf`aA@@@YEEEEETdfB`Hbjj@@@\nf`aAb@NFlBHrJIKIRJUTY@AUST@@@\nf`aAc@DXH@HaxNBLbdRRaRfvcMUP@U@@@\nf`aIB@DXHYp@aIgVUueZZf@Bf@PI@\nf`aPP@B\\@bwlbfdRLRTIgMUT@Q@@@\nf`aQB@NJdBHRYWVyncH@JZjp@@\nf`aQC@IVLBPQHXdLbdLRTvf`eUPADu@@@\nf`a`R@NBT{pHaIe]nujL`@jZf@@@\nf`a`b@KRTCDRVUvuYgJ@BZj`@@\nf`a`b@KRTCDrIQISPiJLEPACUT@@@\nf`a`r@MPQ`ag@|LddqTTRbNkMU@QC@@@\nf`aaB@NRADYEDediiDZL`@ijj@@@\nf`aaP@I@QPQlbbTabaaeNMUSUUP@@\nf`aa`@D@IqInvuURTZjZff`HUV`\nf`axB@DXIhso@BDf]e]VdijZVZZ@`dj@\nf`i@@@DjYeW_YHRg^@@`B@@@@\nf`i@@@DjYeeoyHQm^@@@@h@@@\nf`i@@@LdbbbRVQSEBd{p@DA@@@@@\nf`i@@@LdbbbTRISIBt{p@DD@@@@@\nf`i@@@LdbbbbbRkEB\\[p@@PD@@@@\nf`i@@@LdbbbbbSwEBlYp@@P@P@@@\nf`i@@@LdbbbbbT_EBl[p@@PA@@@@\nf`i@@@LdbbbbbrKEBl[p@@PA@@@@\nf`i@P@@HD[HhdihdhUSbkN|uHPTUB@@\nf`i@P@@HTYIe[VUZLeiwfi@HjhP@@\nf`i@`@@BrJJJJJIXlUiuoP@D@@@@@@\nf`i@`@@DRYfyU]`mNmyi`@@B@@@\nf`i@`@@DRYfyU]`mNmzB@@@@@@@\nf`i@`@@LRYfYuw`eV]zhJ`@@@@@\nf`i@`@@VRYfYU]`eNMyh@`AB@@@\nf`iA@@@YHhhheLUfBdYwjBb@@H@@@\nf`iH@@@RM[IEDhmBdeS`eZL@@DMTD@@\nf`iH@@@RM[IEDhmBdmS`eZL@@DMTD@@\nf`iH@@@RlyJYY{WZ\\Dhu`@@bZ``@@\nf`iH@@@Rl{IEDhmCDcS`eZl@@DUTB@@\nf`iH@@@TUkIEEEhdd\\s`eV\\APUTlL@@\nf`iH@@@XDiJYYe_jRXHs``jJYa`@@\nf`iH@@@XDkIEDhhdmcRSAF\\DEQSLL@@\nf`iH@@@\\myJYfUuZ\\d[S`@@@H`@@@\nf`iH@@@\\myJYfUuz\\d[S`@@@H`@@@\nf`iH@@@\\m{IEEEDedkSdcZ\\@@@AD@@@\nf`iI@@DJlxBSJvmjtYKQkMR@PuQ@@@\nf`iI@@HXIXBRsKkZtxkSoSUP@@@@aA@\nf`iI@@H\\EhFQJJIPiIKdeZ]z`HajTh@@\nf`iP@@@NRifyWU`eVMzjPZBbP@@\nf`iP@@@VRifyWm`eNMzjPZBbP@@\nf`iP@@@^RifyW}`eVLzjPZBJP@@\nf`iPB@LD@DYHhbhheD\\TtYwjY`@@H@@@\nf`iPBDK^@D@eJYfUuZ\\d[S`@@@H`@@@\nf`iP`@BJ@aInUUm^BUiwfjjh@F@@@\nf`iP`@BR@IIf]UUV\\UiwfiZZ@H@@@\nf`iP`@B\\@aInVUuQRUiwfjij@H@@@\nf`iP`@FB@QIf[eU~Rtzwff@@@H@@@\nf`iQ@DK^@BTifYWUirQmN@@@@b@@@\nf`i`@@@IKLlj{pRkN}p@@@@@@@\nf`i``@M@YdTTTTQaRyIQg^jBb@@P@@\nf`ih@@@LDYvRJIQPiIJcIJtX@@`ZhH@@\nf`ih@@@\\UYvRJJJJsIJgIZMX@HBJhH@@\nf`ih@@@\\eYvRJJJJUKJgEZMX@HHFhH@@\nf`ip@@@N}dbbbbRQrhJQmV@@@@J@@@\nf`ip@@@V|eLsJzmNRMip@@@DP@@@\nf`ip@@@XTeLwOvmNJt{pPQT@A@@@\nf`ip@@@XUdbbRRRaRiqVg^BBJ`@H@@\nf`ip@@@XdeLv{ZmNRUYpPQP@Q@@@\nf`ip@@@XdeLv{Z}NRUYpPQP@Q@@@\nf`ip@@@XheLwKkkNRT{pQAU@@@@@\nf`ip@@@XmdbbbTTqVXJQg^@H@@@@@@\nf`ipB@BBdBHRYeY_Y`nRlzBAJiXh@@\nf`iq@@BBdAIf]UUV\\UiwhH@@@@@@@\nf`iqB@HLEpBHeJkOrnBBd{uKSL@E@@@\nf`q@@`@QBhbtQzHRYegeudcN``@`b@@@\nf`q@A@@QGhaIe]fuzLlz@Bh@H@@@\nf`q@P@@HMyIgUe[]V]yjZ@Hh`PI@\nf`q@P@@HMyIgfU{]V]yjhDHh`@@\nf`q@P@@HTyIgUfUc^lyjZBB``PI@\nf`q@P@@HU{HhlhheEbhvkMULttr@@\nf`q@`@@LRYfWg^Qg^ZB`@@@@@\nf`q@`@@LrJIJHiQITx{t@QL@@@`J@\nf`q@`@@^rJIJEJIKTYYt@Dp@P@ar@\nf`q@`@@^rJJJJHpkTXkSPTPUT@@@\nf`q@a@JZADFbLbTTTtRROIN}TmPLA@@@\nf`qA@@@ILrj}jISoM@a@@@BEH\nf`qA@@@YDhhdmECJ\\UYfiVjZP@@\nf`qA@@@YDhhdmECJ\\UZfiVjZP@@\nf`qAB@C@bLbRbbbebpmV]TAECTP@@\nf`qAP@@De[vQQQSJJFTVcNZdJHbP@@\nf`qA`@@FmdbbdTUrRiIVjBh`jd@@\nf`qA`@@N|dsLjolcZlt@@@T@@@\nf`qAp@@HiIUgYEDhdeDUBLEIjfjjZ`@@\nf`qHB@DXIp@PdskKJzcN|uL@@@@HD`\nf`qHP@DXxHDg^rJIQSPiITyhsSUUUS@@@\nf`qH`@DLeXELbbTVRTQfcN|uMUUUP@@\nf`qH`@DXEXFdfYfWVYHufjjYZj@@\nf`qH`@DXUXAlbbbJRRfaJN|uUUSUP@@\nf`qH`@DXtX@lbbbbbQRSIF\\uUMULp@@\nf`qIB@DXIHP@aIgYfuyJ]yjZXHB@PI@\nf`qI`@DXEXBmYEEEEDldVRMYjjfVj`@@\nf`qP@@@PRkfUWvQmVj@@@J@@@\nf`qP@@@^RYWUe^cKN`@j@B@@@\nf`qP@@@^Rfuue]gKNh@J@B@@@\nf`qPA@NBAD^bDfUuYWhrsh@I`@`ACd@\nf`qPB@DX@DILwLjoiuoMUUUUT@@\nf`qP`@DBAKHheHdhbmPSoMUTa@P@@\nf`qP`@DXAqIfWnUrL]yjihHB@@@\nf`qQ@@FRAdTRTbdTLMYwfjjjjj@@\nf`qQB@G^`@HRkfUWvQmVZ@@@J@@@\nf`qQP@DX@ISoILrl{jsfcMUMMUT@@\nf`qQb@LD``P@cIICDhhlUrSoKUST@@@@\nf`q`P@L@PHrSJwJwJmhsUAAQQ@@@\nf`q`P@L@PZvQQISQIPiVg^ZhBBjH@@\nf`q``@D@EdTRdTTqVYIwfji``H@@\nf`q``@D@HdsOZvjmN|uT@@@@@@\nf`q`b@ILD@HrQQJJMQITX{uMTt@D@@@\nf`qa@@D@RYyeg^Qg^Z``@@@@@\nf`qa@@E@rQQIZYQF\\EjuA@TuQ@@@\nf`qaA@ODADRbLbbbRRaRsAV]AAD@D@@@\nf`qaP@D@hIVdfVfy{IrgFjjfjj@@\nf`qa`@L@PKHhdheEhchsoMTD@@@@@\nf`qa`@L@QqIgfUnyZTYjBBJh`@@\nf`qi`@DTxIPC^rJIQQJiILyISUKMUU@@@\nf`qi`@FRt{pDDR[mfUuVkNZjhDB`@@\nf`qp@@@Hpds\\rj~gV|uP@@@@@@\nf`qpP@DXxBQoYEDhihUDZBuYijjji`@@\nf`qpP@DXxBSoYEDhihTdjBtYijjji`@@\nf`qpa@LR}A@@`PQddaTTRRUtxuejBPha@@\nf`qq@@BBdAIf]ueV\\]zB@f@@@PE@\nf`qqa@FRUAE`AE@cIICDedeDJQgKTuMMR@@\nf`y@@@DjYeU\x7FYHRmF]x@B@D@H@@\nf`y@@@DjYegvYHQmNmx@@@BHA@@\nf`y@@@LdbbbbbfkEBMIuo@@@@@@@@@\nf`y@@@LdbbbbbfkEB]Jqo@@A@PA@@@\nf`y`@@@ISKJmmQ`cZ][p@@@@@@@@\nf`y`@@@ISKONwS`iZl{p@@@@@@@@\nf`y`@@@YHdhddhTbTeiQg^@@@@@@@@@\nf`~@`@@HR[UUUUjjjjh@@\nf`~@`@@TReWUUVjjjjh@@\nf`~@`@@XRfUv}ZjijihDFP\nf`~A@@@IJskNmSUMTtABGb@\nf`~A@@@ISjjjkSLsMTBbDbqDjP\nf`~AJ@HHpSdmV|CprRSESSSSMUUUT`@@\nf`~P@`LZ@DHBBAB`cHheDdiXduUTuT@@\nf`~PQ@FBUhu`AD@aJUYn{zjjfZ`@@\nf`~X@@@TYhwdiUUU]jjjjj@@\nf`~`J@NPQafcN|CprRRjYQQFMTuSS@@@\nf`~`b`KLLBHHDTBNA@`aIneUYYjjjZ`@@\nfbc@@@DjYeUU]VJEKQ`ybd@B@D@@B@@@\nfbc@@@DjYeUoUvRDkQoIcd@B@@@@B@@@\nfbc@@@DjYeU\x7FUVRDkQg^Sd@B@D@@B@@@\nfbc@@@DjYeeU_VZEJt_ICd@@@@`D@`@@\nfbc@@@DjYeeou~RD[W`qSd@@@B@P@`@@\nfbc@@@LdbbbTRrTroIBMks`ir@@@@Q@@H@@\nfbc@@@LdbbbTVraVSIBMkugIr@@@AP@@@@@\nfbc@@@LdbbbTVtRqSEBMjsoIr@@@@T@@H@@\nfbc@@@LdbbbbbVVrOEBlXKoIr@@D@B@@P@@\nfbc@@@LdbbbbbbQQsEBMKuhir@@@@@@@@@@\nfbc@@@LdbbbbbbQQsEBMKuhir@@@@QC@D@@\nfbc@@@LdbbbbbbcJwEB]Hu`ir@@@@@@@@@@\nfbc@@@LdbbbbbbcJwEB]Hu`ir@@@@DS@D@@\nfbc@@@LdbbbbbbcJwEB]Jq`ir@@DA@@A@@@\nfbcA@@@ILk\\l|{txhug^BgKRtBPBAQ@@@\nfbcA@@@ILsLjmwLDkqhnRGKP@@@@@@@@@\nfbcA@@@IRlrrk{pQSbgVBgM@@A@@@@@@@\nfbcP@@@JrQQJJZYIEITxKQkN|gH@@TUP@D@@@\nfbc`@@@IRlrrk{pQSbgVBgMT@A@@AP@@@\nfbc`@@@YHheDdihhTZLD[P[^Sd@@@@@@@@@@\nfbc`@@@YHheEMDXlijREhuoQSd@@@@@@@@@@\nfbeA`@@^rdrjljmkXKSUAAQT@ab@\nfbe`P@J@QjNQQJZZJJJYWrVjjAbJH@@\nfbe`a@KDtBlEVYEEDeDh]ddfR``hfzj`@@\nfbeaS@Lfl{rPPd@a@QddRTRRTTWTtFj`hFXf@@@\nfbepP@LZbAC`iLjmmnvjtuUMUTuP@@\nfbepb@BZ}IPHaIf]V_VUgHHBfjjh@@\nfbm@P@@BMGHhhhhhhm\\esRoNBt@PT@UP@@@\nfbm@P@@HD[HhdihdiddYSbk^BuHPTUMP`@@\nfbm@Q@NV\\DDPfHrJJIQPqIJFLUKT\\pPM@LQA@@@\nfbm@`@@HRYvYfUWTYrPyZh@@@H@@@@\nfbm@`@@YRYVvW[WhpThi`@ij`@b@@@\nfbm@a`FZADDbAQ@XaLQfHrJJJZIFYJK\\ehwhpQDT@QQ@@@\nfbm@a`FZADDbAQChbLQfHrJJJZIFZKP|ejpTpQDTDDQ@@@\nfbm@p@@JMxJ\\bbbbbbbRegAZb{KPADPLQ@@@\nfbm@p@@JmxJ\\bbbbbbTRacAZb{KPADDLQ@@@\nfbmAP@@BLENQQQQQQQYG[bm^]Eh@bHAJ`@@\nfbmAP@@BLGNQQQQQQQYGGbm^]Eh@bHAJ`@@\nfbmAP@@BUGNQQQQQQPeIG`mNBehIb@@b`@@\nfbmAP@@HCGNQQJIYJJKSDeEZBejjjjjjj@@\nfbmAP@@HeENQQJJJIKZIUgMNmyj@H@Bj`@@\nfbmAP@@JMxNQQQQQQQIRw`mQ\\eh@bHFB`@@\nfbmAP@@\\e{NQQQJJKFIKYgEFbehFH`HB`@@\nfbmAP@@\\e{NQQQJJKFIKYgEFbehJH`HB`@@\nfbmH@@@\\myJYfUu^zgIFtx@@@BJj@@@\nfbmH`@EVBdGlbbbbbbTRacAZb{KPADDLQ@@@\nfbmH`@JRUDBTin]UUU~Bt[pZB@@@Bh@@@\nfbmHp@J\\DhBC^SgHheEDjheEdsdeNBuRsUUUU@@@\nfbmIR@LLDiAoIr@HrJSJIHjIQY\\EIS`uUJuUUUP@@\nfbmPB@NA@DYHhhhddmcEJ\\e[Sj@H@`jj@@@\nfbmPB@NQ@DYHhhhddmbeJ\\e[Sj@H@`jj@@@\nfbmPBDK^@D@eJYfUu^zgIFtx@@@BJj@@@\nfbmPP@KA@kr\\bbbbbbTRacAZb{KPADDLQ@@@\nfbmPQ@ARMyNPQAHaIf]UyYugIJ\\FB@ijjjJ@@\nfbmPa@NVtx@PHHrQQQJIPjXiTD[r\\uK]S[@D@@@\nfbmPb@AJ|dDPdrmrljoSfgQs@DMPTBD@@@\nfbmQB@AJRBHRYVyV[WisShy`BFhJAB@@@\nfbmQB@EZBBHrJIJZIQQYHtxhr\\pACUAHA@@@\nfbmQ`@DX@sJSJmsJkyJNbEKUM@@@E@@@\nfbm`P@A@BgJSLsLoOsEZ}ZKPADABM@@@\nfbm``@A@BdsLsJk|qVoVbt@Q@PeP@@@\nfbmh@@@XIQRTsNj~rkNRtgHD@DU@@@@@\nfbmhR@DXyJroIrBPrJIQQJIYPyTDZs`sUMUUURp@@\nfbmi@@DTxIPLbbTTRqbRfSNRUxKUKM@AUP@@\nfbmiP@DTxIPC^SgHheEDjheEdsdeNBuRsUUUU@@@\nfbmiP@DTxIPMVSgHheEDbmDddsde^BuRsUUUU@@@\nfbmiP@DXXKPC^SgHhdheDYeEeSdmNBtuUUUUK@@@\nfbmiP@LF]EHIJugIEEEDdhjllRegAcUUUPHET@@@\nfbu@AP@QAHadPJHeDZbGQ@XaLQfHrJJJ[PrZJ[TyzMLuUULuP@@\nfbu@S@FZCDDPRH\\DILs\\}|r|dgLDQDuKTH@@\nfbu@`@@IRYV{eg^irNX@aj```@Py@\nfbu@`@@YRYWYeg_hrJVeffjef`@@\nfbu@b@AVADYEDeCMHiXdjLSf@Bji`@@@@\nfbu@b@AVADYEDeeMHiXdjLSf@Bji`@@@@\nfbu@b@AVADYEEDeKHiXdf\\Sf@`ji`@@@@\nfbuAP@@BtFNQQQQQJKGIISk^Z@HH@j`@@\nfbuAb`BJ\\BHxDRbHqFXcHhhmDiXhmDpPTpPQLDDP@@@\nfbuHP@DTtdCHyrJIJHqSSIQTxEKUMUUTsP@@\nfbuI`@DXtX@dyEEEEDbddUFRLyjjZjZfX@@\nfbuIp@DXyHBcV|gNQQJJIQYJGJ`gAZZjjjjj@@\nfbuPJ@CA`cANJ}dGaddfJbTQbbvejJVjjjjYf`@@\nfbuPb@GQxDDPdrnlkk]QakP@T@EMP@@@\nfbuPr@FVQSaoIAHIM\\w[KntDxMUUULuMP@@\nfbuPr@NBYhwdyAbIJsNk[NbdxLtmTuUUPHP`\nfbu`P@D`EGNQQQQQUJYJIb`iZjjjjfZ@@\nfbu`a@HLx@JpDYHheHiXddmgQBfje`@@B@@@\nfbuaP@L`pIUdifV^}}ZR|FjjjB@k`@@\nfbua`@D@uYIge]nyTTdyZ`@jh@@@@\nfbupb@IZCrpHaIe]YfWjcOA`@jjYZb@@\nfbupr@DDhJRgIrBHRYYymY_IRwiiijjjj`@@\nfby@`@@HR[UUUUUZjjfjj`PT`\nfby@`@@HR[UUUUUZjjjjj`@@\nfby@`@@YRUUUUUUjknknj`@@\nfbyA`@@LdeJl{jjkSUMLsU@haxjSKd\nfbyPH@BZ@bR`qSgHieDeDXdhhuUSUTs@@@\nfb}@@@DjYeeou~RD[Wdy@@@@`H@@@@\nfb}@@@DjYee\x7F]^RD[S`q@@@@@B`@@@\nfb}@@@LdbbbTRrJROEBMipT`@@@@DP@@@\nfb}@@@LdbbbbbRJvcCBTzHT`@A@@@P@@@\nfb}@@@LdbbbbbRVNSCBTzJ\\`@QDTB@P@@\nfb}@@@LdbbbbbRVNSCBTzJ\\`@QDTD@P@@\nfb}@@@LdbbbbbT^VgEBl[tT`@A@DTpD@@\nfb}@@@LdbbbbbT^VgEBl[tT`@A@DUPD@@\nfb}@@@LdbbbbbTqrOEBlXLT`@A@EEPB@@\nfb}@@@LdbbbbbqfvSMBlYr\\`@A@A@@@@@\nfb}@B@@RBSJwKON}NJMYr\\p@@@@P@@@@\nfb}@B@@RBSLl|{vkART{r\\p@@A@@@@@@\nfb}@B@@RBSLl}kvkART{r\\p@@A@@@@@@\nfb}@B@@RBSLsKoN{AJlxJ\\p@P@@@@@@@\nfb}@B@@RFQQISQJJZUJgEFmyNX@@@B@@@@@\nfb}@B@@RFQQQJIIQIKEgIF|ENX@@@H@@@@@\nfb}@B@@RFQQQJIIZJYEgIFmyNX@@@@H@@@@\nfb}@B@@RFQQQJJIJIIEgIVCfJX@@@B@@@@@\nfb}@B@@RFQQQJKPiIII`iJlzJX@@@@@H@@@\nfb}@B@@RFQQQJKPjJEI`iJmyNX@@@B@@@@@\nfb}@B@@XfQJJJJIIPy[bmNlENZ@@@@@@@@@\nfb}@B@@XfQJJJJJIJHgfc^\\fNZ@@@@@@@@@\nfb}@P@@HuYIeYUWUVcEF]yNViBBJ`B@@@\nfb}@P@@HuYIeYUWUVcEF]yNViBBJjjF@@\nfb}@`@@DRYfu{VuXIVg^Sei`P@@@`@@@\nfb}@`@@DrJJJZIQIIHlEksdqsA@@@@@@@@@\nfb}@`@@IRYe^um]HpTcVCEhH@@@J`@@@\nfb}@`@@LrJJJJHpiII\\dkSoIsTED@@@H@@@\nfb}@`@@LrJJJJHpiSI\\dkSoIsTED@@@H@@@\nfb}@`@@LrJJJJHpqJH|dhw`isTED@@@H@@@\nfb}@`@@LrJJJJHyISI\\dkSoIsTED@@@H@@@\nfb}@`@@LrJJJJHyJIK\\dkS`isTED@@@H@@@\nfb}@`@@LrJJJJJFHpi\\UjsoIrtCA@@@H@@@\nfb}@`@@LrJJJJJFHpi\\UjsoIs@D@@@@@@@@\nfb}@`@@LrJJJJKIIIH|Djw`isTEP@@@@@@@\nfb}@`@@QrJJJJIIQ[Rldiq`~Rt@@@@QP@@@\nfb}@`@@YRYVumUehrVkNBeijZfjjj`@@\nfb}@`@@YRYVum[ehrRkNBf@BjjjjJ`@@\nfb}@`@@YrJJQIFYIIKDyKUgQRuUTuPAT@@@\nfb}@`@@YrJJQIFYJYKDyIUgQRuUUU@AU@@@\nfb}A@@@IKLlrjk|Dip\\qSP@@PEAA@@@\nfb}A@@@ILrrkszlyHuoIs@@@@P@@@@@\nfb}A@@@ILrsvoZlEISoIrtDPD@@P@@@\nfb}A@@@ISZvmkZlxkSoQSP@@@AU@H@@\nfb}A@@@ISZvmkZlxkSoQSP@@@AU@P@@\nfb}A@@@IS\\lrr}|UjwlacP@@DAL@P@@\nfb}A@@@YDhhhhdeCenJtzpTyX@@@@@@@@@\nfb}A@@@YDhhhhdeCenJtzpTyh@@@@@@@@@\nfb}AP@@HeGNQQJJIQKIIVgIJlGrVi``@BJ`@@\nfb}A`@@^sdTTTTRRRTqXIScQBeh@`@HB@@@@\nfb}Ab@HHp@HrRRqQJKQYSRtZwhirm@@@@@@@@@\nfb}P`@FJ@yIfY^]oTcIFmyNV`hJI`B@@@\nfb}QB@I^DBHrJIKJSQQIHtYjw`is@AUTC@E@@@\nfb}Q`@EV@`rSJvlkjmFRt{r\\mR@A@@D@@@\nfb}`@@@ISLsJr{tdjqh~r@@@@AD@@@@\nfb}`@@@YHhhhhdddbj\\dxtPy@@@@`B@@@@\nfb}`@@@YHhhhhdhebZRUXt_I@@@@@``@@@\nfb}`@@@YHhhhhdhehZRUXt_Y@@@@@b@@@@\nfb}`@@@YHhhhheEDmZRTxrPy@@@@`H@@@@\nfb}`@@@YHhhhhecLej\\eXsdy@@@`@`@@@@\nfb}`B@I@qDefYe^}~RtzpTyh@@@@@@@@@\nfb}a@@A@RYVvWVuhrVg^Sf@@B@@@@@@@\nfb}a@@B@rJJIJZKQZHldKW`is@@@D@@@@@@\nfb}a@@I@RVYfU{wyKSkASe`@@@@@@@@@\nfb}a@@I@RVYfU{wyKSkASf`@@@@@@@@@\nfb}i@@DTxIPLbbTTRqbRfSNRUxJ\\mTlt@EUD@@\nfb}pB@ARlBHRYWYee}hpsk^Sf@Bjih@B@@@\nfde@@@DjYeew]daZCD@BBBH@@@@\nfde@P@@BLGHhhhhhhlcqVoNbt@QD@d@@@\nfde@P@@B\\EIfYfWgXkWkQZ@H`HR@@@\nfde@P@@B|EIfYfW^XkSkQZ@H`Ab@@@\nfde@P@@B}EIfYfW_XkSkAZ@H`AJ@@@\nfde@P@@HmEIeYVuuhqQkNZdHHfjH@@\nfde@P@@HuyIe[UWZhrRcVZd@JjjH@@\nfde@P@@J]{HhhhhhdeBpV`vbt@Q@LD@@@\nfde@P@@J]{HhhhhhdeBpV`vc@@A@T@@@@\nfde@P@@V}GHhhhhdddfpRgFBt@P@dT@@@\nfde@``ARADDb@qDXaIf]UUUYqVg^``Jh@h@@@\nfde@b@AZAHILkmmrmFBdFL@ETp@D@@@\nfdeAB`A@bBQCH`LQFHrJJJZEIII\\dkSoPQDT@T@@@\nfdeAB`E@bBQCH`LQFHrJJJZESII\\DkSoPQAP@TP@@\nfdeAC@A@bKQDXaIe]ke_hpRoA`@j``H`@@\nfdeA`@@BBdsLsJklUkuhm@DPDI@@@\nfdeA`@@DMdTTTtRRVLXKSk^ZXFBbjH@@\nfdeA`@@ZcdbbTRVTJVipRoAhHbH@J@@@\nfdeA`@@\\udTTRbVfaRYrQ`qZA`h@B@@@\nfdeH@@@\\mEJYfUuWirQmN@@@@bh@@@\nfdeH`@JZ|DBDfV]VUugEF]yjeifZj`@@\nfdeI@@DRUXFQQZJIPiIJVcVCEj@bFiR`@@\nfdeI@@DTYxBSOKKropQgAbuMT@@@@@@\nfdeI@@D\\tXBSOLkK{KRcAbuA@eTiP@@\nfdeI@@HXHHBRrr{lkQbmAcSU@D@@@@@\nfdeI@@HXIHBRsKj|kSbmAcMU@D@@@HPP\nfdeI@@JRUDBTwNjjkpVc^CPP@@@P@@@\nfdeI`@LZCDDeYHhhidhblbTM[pZjjhDBh@@\nfdeI`BLTEhDGZewIWNlkJ|DZLFMUMUUUU@aR@\nfdeP`@AN@EIfYfWgXkWkQZ@H`HR@@@\nfdeP`@E^@{HhhhhhdeBpV`vbt@Q@LD@@@\nfdeP`@E^@{HhhhhhdeBpV`vbuHSPME@@@\nfdeQ@@DXAdTTtTTRbqUXp_QZh@``A@@@\nfde``@A@BdsLsKslUiuhm@DPAI@@@\nfde`a@AZlBHYDYEDeDXhiij\\Lywh@bifjJ@@\nfde`b@H\\d@HRfYfWwQJMxLZjfZ``H@@\nfdea@@D@RYyVuUQRMXLVh@@@@@@@\nfdea`@D@DGHhhhieeLbrVkN|uPL@AP@@@\nfdea`@O@T{HhhhhhdeBpV`vbt@Q@LD@@@\nfdei@@JJMYpLbTTRbbQVwMBUxMAATDAB@@@\nfdei`@LZlFHIJrQQQQYQEYDhZw`uUUPHEP@@\nfdei`@LZlFHIJrQQQSIQEYDhZw`uUUPHEP@@\nfdep@@@NCdbbbbRVTqhJQm^@@@@hh@@@\nfdep@@@XMdbbTRRReRhHu`qBBbj@@`@@\nfdepA@BLBBlmVYEEEDhXXmNJd[tX@`@HH@@@\nfdepB@DXb@HrJIQIIHqIDhkPXsUAA@@@@@\nfdi@P@@Hl{HhehiheDT[sTuUSMPHD`\nfdi@S@HRMX@PkptDIUvnmZoAM@ASUT@@@\nfdiAP@@HtxJSOLkkMISTAATs@@@\nfdiPa@LR``@PHHrRPjJJKHyF\\mPLQEP@@\nfdiQ@@DT@drsMsJvbuTEDU@@@\nfdiQQ@LJ``ug@BM^FRRVIKIRIHTeh@Ijj@@@\nfdiQb@LR``P@cIIBhheedmjru@pQU@@@\nfdiQb@LR``P@cIIBhhhlcdYru@qDU@@@\nfdiY`@D\\Eiw`TbdloJskjuUT@EUP@@\nfdi`a@OJtCDMbIKJrwzsAT@QMUT@@@\nfdiqb@JLtRR`RFQQQYIPiIIgJVjVjj`@@\nfdiqc@DXipu`QBoCPQdTRbtQdrRIRZjjZjd@@\nfdm@@@LdbbRbaRRsA\\UKU`~b@@@@@@P@@@\nfdm@@@LdbbbbbfjSCBUIuoAb@@@@APA@@@\nfdq@P@@HMYIeUU{UZjjijh@@\nfdqPP@JV@bT\\bbTrTQRRSUUMUL``z@\nfdqPQ@FBUhu`AD@aJUYn{ujjjYjh@@\nfdu@@@DjYUUUU`nR\\GtP@@@@@@@@@\nfdu@@@DjYee]}daVuxLP`I@H@H@@@\nfdu@@@DjYee]}daZlGtP@@@@@@@@@\nfdu@@@DjYee]}faRlGtP@@@@@@@@@\nfdu@@@DjYee\x7F_daFtxLP@@@@@@@@@\nfdu@@@LdbbbTRerwIBMiw``@@@QT@`@@\nfdu@@@LdbbbTRerwIBMiw``@@@QTA@@@\nfdu@@@LdbbbTRrJwEBMipXp@PP@A@@@@\nfdu@@@LdbbbTVtjSEBMjpXa@RA@@P@@@\nfdu@@@LdbbbbbTUHeRlXOh`@@@@@@@@@\nfdu@@@LdbbbbbVVOEBlXOh`@A@EBA@@@\nfdu@A@@Q@XaIfV^]nXJRc^]F@@@@B@@@@\nfdu@A@@Q@XaIfV^^^XJVcN}F@@@@B@@@@\nfdu@A@@QBhcHhheDeDULpTmV|FL@@D@@@@@@\nfdu@B@@ABTsJk^wb`mJmzM@@@@DP@@@\nfdu@B@@ABTsJk^\x7Fb`mJmxM@@@@AP@@@\nfdu@B@@AFRJJJIIEYHeIFuXLZ@B@@@`@@@\nfdu@B@@AFRJJJIIQIDeIFtEtZ@@@@``@@@\nfdu@B@@AFRJJJIIQSDeIFtEtZ@@@@``@@@\nfdu@P@@HLGHheDdhidYJqk^bFKUUUUUUP@@\nfdu@`@@ARYegg[fBdhwgQZA@@hZA@@\nfdu@`@@ARYegggfBehsoQZA@@bZA@@\nfduA@@@IKJrkN{NJMipXt@A@@@P@@@\nfduA@@@IKLl{[oAJmhpXt@@E@@@@@@\nfduA@@@ILroZviFBdZpXmA@@@A@@@@\nfduA@@@YDhdhddhTf\\uXs`qh@@@@@@@@\nfduA@@@YEEDhmEhTjLehu`qZB@@@B@@@\nfduP@@@ARige[uVBeip_Qjjj`Pbh@@\nfduQ@@FRAdbbTJbRrRIQVcVCFiPHJ@H@@@\nfdu`@@@IRmjjviJBdZpX`@@@@@@@@@\nfdu`@@@IRmkZviJBdZpX`@@@@@@@@@\nfdu`@@@IRmrsziFRtZpX`@@@@@@@@@\nfdu`@@@YHheEMDXTjREhu`q@@@@@@@@@\nfdu`@@@YHheEMD\\djREhu`q@@@@@@@@@\nfdu`@@@YHheEhTddj\\EHu`q@@@@@@@@@\nfdu`@@@YHheEhTidj\\EHu`q@@@@@@@@@\nfdu`@@@YHhhhhdhbZRUXp_Q@@@@@@@@@\nfdu`@@@YHhhhhdhcjRTxp_QB@@`H@`@@\nfdu``@C@PdrmljmtYKUoAbuHAEUTT@@\nfdua@@E@RfVU]uZLEhu`qhJ@@@@@@@\nfdua`@E@QyIefUW]iqUcNCEj@@@JhD@@\nfduh@@@XUkrTs\\zsoSbcZ\\FHPQTa@E@@@\nfdy@P@@BtGHhhhheEcliuoM@DD@T@@@\nfdy@P@@HHkHheDdhcLeHpXmPPL@@@@@\nfdy@a@CVADQbLbbRrbbbQuF^C@ATC@P@@@\nfdy@c@OAAD\\BAABSJzlrsQchp@UA@p@@@\nfdyA`@@FcdTTTbbveVhsUf`jHbj`@@\nfdyAa@OAbBPUHYEEDhddYeFRMZB@B@i`@@\nfdyH@@@PL{IDhhddih]JvoM@@ADT@H@`\nfdyHQ@DDiIU`qADBBLbbTVbaTVQJF]MMMUUU@@@\nfdyIP@DXTx@cAbdsJslolyjsUSSUUTBDh\nfdyIP@DXxHDkAbdrs]k^txYsSUUUUL@@\nfdyIR@DDiHRkAbBHrJIQZJEQYDhYtttuUUT@@\nfdyP@@@HrJSQQUQIQ\\d{uTuT@EH@@\nfdyP@`ARADDbOQDXcHhhlleDeBsc`pPDpAD@HB`\nfdyPA@FN@DRBDiee\x7Ff\\R`qhJBHH@@@@\nfdyPP@LQ@bSdfUYiUmSkQZjBhJb@@@\nfdyPQ@LN``u`AF`cIICDmdiCEQg`kUUTmUP@@\nfdyP`@AR@EJ[WUe]Yqwj@Bh@J@@@\nfdyP`@NJ@aIeeU~UehLVhHD`@@DIP\nfdyPp@DX@IU`qrJJIJIIPyLxYsUSSUUT@@\nfdyPr@JLIIU`qAHILsjvskNF]RuSUUU@@@\nfdy``@F@PdwLk{JbTFKPPQP@@@@\nfdya@@D@RYyWUeQRCEj@Bh@@@@\nfdya@@D@rJJJFJIIKVCzKUT@@@@@@\nfdya`@I@UyIfUVyUEhLV``hjjH@@\nfdyaa@AZMXDPrHrJIJHqQRqTxYt@QTsTd@@\nfdyhP@AZlFHDDXdw]lrnfmxKUUTBAP@@\nfdyhR@AZlFLDXhCprRRiIQQIYFmxKUUTBAP@@\nfdypB@FR\\@HrQQHqQIIZXi[uLuAAEP@@\nfdypP@DT|@chyEDidhdeDQRlEj@@`X`@@\nfdypQ@FVcSao@dDRFQRZJZYQKF`gAjjjifZ`@@\nfdyp`@BRb@RSLzkrnsg`mRtt@E@HB`\nfdyp`@LLxABSKrjszJQgMTu@QT@@@\nfdypp@DXT@QgAbdsJslolyjsUSSUUTBDh\nfdypp@DXxBVkAbdrs]k^tDYsSUUUUL@@\nfdyppBDDDCbkAbR^]dTtRaTRbrXHsfiZjjjh@@\nfdyq`@LLxA@TfWe^UtTg^Zij`@i@@@\nfdyq`@LLxACdfWeUWTTcVZij@Jh@@@\nfdyqb@IFmcpPABT{JvnvIVcUPACTtP@@\nfdyqb@LFcA@`AFRREQQIJH{WcVVhIBhhP@@\nfd}@`@@AReYU[WhHirQmV}DLZ@`P@HBB@@\nfgA@@@DjYU_VByHu`@@@@@@@\nfgA@@@LdbbbTVKIBMjp@@@@@@@\nfgA@`@@\\RfYe_irQmVh@`@H@@@\nfgAAB@F@bDfUmUZ\\EHuh@b`@@@@\nfgAA`@@HLdrmlkQdmFluHADq@@@\nfgAA`@@HtdrmjkQdeFluH@UP`@@\nfgAH@@@RM[IEDhmBeS`eZL@@DMA@@@\nfgAH@@@XDYJYYmvdfBuXHJbYF@@\nfgAH@@@XDkIEDhhdkRSAFlDEQRc@@@\nfgAH@@@XhiJYmg]gAJMXH`jeF@@\nfgAP@@@LrQQJEIITxJQk@QE@@@@@\nfgAP@@@LrQQJJYPtdKQk@@@@@@@@\nfgAP@@@\\RfYe_irQmV@@@@@@@@\nfgAP@B@BAHirRJIJHkLUpQkAPP@A@@@\nfgAPB@LD@DISJ\x7FJdhsakTuT@E@@@\nfgAQ@@D\\@drjjkQ`iFm@AT@@@@@\nfgAQ@@F\\AdbbTJRRIPTcVjUj@H@@@\nfgA`@@@ISKJntXKQk@@@@@@@@\nfgA`@@@ISLrotyHvk@@@@@@@@\nfgA`@@@ITsZkdYHvkU@@@D@@@\nfgA`@@@YHheEMEZRDhu`@@@@@@@\nfgA`@@@YHheEhTj\\EHu`@@@@@@@\nfgA`B@N@BDifYWz\\d[Uj@H@B@@@\nfgAa@@K@RfYU_qPVeFh@@@h@@@\nfgAp@@@RMdbbTVaVipRmF@@BF``@@\nfgAp@@@XDeLzmkQ`iFlDAP@D@@@\nfgAp@@@XIdbbaRRRqPTcVZ`@@H@@@\nfgAp@@@XheLvnjs`iFlDPT@@@@@\nfgAp@@@XxeLzjkQ`iFlDAT@@@@@\nfgApB@LLx@HRevUUpPTcViYj@H@@@\nfha@P`BLxHDQSp}DAbDeUkUUZjZjjh@@\nfha@P`BLxHFISp}bAqDeUkUUZjZjjh@@\nfha@R@HHpPG`eUjjjjuUUUU@@@\nfha@s@DXHJR`ADOBgadTbRTbVQSUTuUT@@\nfhaA`@@HhdvvjjjuUUUU@@@\nfhaP`@DFAkHhdddeDTfjjijh@@\nfhaPq@BA`cAg^@DVBLddNRRRRdVjjje`@@\nfhaYA@DXHKPY@BOABSNrrzsTtuULAAKd@\nfhe@@@DjYfYuvJEXpW^@@B@I@H@@\nfhe@@@DjYfYuvJEXpW^@@B@J@H@@\nfhe@@@LdbbRbbvJYKfaJ|D@@bH@@@@@\nfhe@@@LdbbbTReryHQmN|D@@@BI@`@@\nfhe@@@LdbbbVTJrYHQeZ|D@`B`@@@@@\nfhe@@@LdbbbbbbcxhQiNlD@@@@H@@@@\nfhe@B@@ABTrnjjnJ\\Ehw`tB@@@D@@@\nfhe@B@@ABTrnkJ~J\\Ejw`tB@@@D@@@\nfhe@B@@ABTroZjnJ\\Ehw`tB@@@D@@@\nfhe@B@@ABTroZvnJ\\Ehw`tB@@@D@@@\nfhe@B@@ABTrssziJBUhw`tB@@@D@@@\nfhe@B@@AFRJIKSQJHxirQk^CPH@@@P@@@\nfhe@B@@AFRJIQZKPiDhIVc^CPH@@@P@@@\nfhe@B@@QBSLl|{[ARTYw`p@@@@@@@@\nfheA@@@IRlk\\~pQSdm^CUAAD@T@@@\nfheA`@@B|eJkjrmAEFNCxMUUUUUU@@@\nfheA`@@TdeJsZvoACNZmxMA@P@@@@@@\nfhe`@@@ILkJl{tYKRk^C@@P@@D@@@\nfhe`@@@ILrsZklxJQoNC@P@@@D@@@\nfhe`@@@YEDeDdeBeQbcZmxL@@@@@@@@@\nfhea@@D@RYYU[WirVcV|F@@@@@@@@@\nfhea@@D@rJIJIIJEJcEFu[pX@@@@@@@@\nfhea`@E@iIJUfv]nBFBL{pZB@`@@@@@\nfhep`@BLT@NQQQQQQRq\\tJQkN|uLuUTuP@@\nfhi@B`@QAHadPzHCDYEEELUDeCpUoPQA@DP@@@\nfhi@R@HHpPG`eUoLjn[s`mU@@@@@@@\nfhi@`@@ArJIKJQPiZcG^`@j`@`@@\nfhi@a@AFADAbDfUmygZL]z@BhHB@@@\nfhi@c@BFADRb@qFQQQIIPqYLD{tDDPPD@@@\nfhi@c@BFADRb@qFQQQIIPq[LD{tDDPPD@@@\nfhiAA@C@b@qBSJ{Lk}FN}@ATDA@@@\nfhiAP@@Lt[vQQQRJKYE\\d{sSUTruP@@\nfhiA`@@B|dsLro~jqgM@D@AP@@@\nfhiA`@@HBdrszroScoMSLp@T@`J@\nfhiA`@@H|dwMkZ{IUgMT@@AP@@@\nfhiA`@@LCdTRTQRbRuNN}@DS@A@B@h\nfhiE@BBTYIQoA@rY\\drkKN{hsoUTt``T@@@\nfhiHA@MFlDDPRHrJJKKQQQUgK^``HBH`@@\nfhiHP@DXyhDoARYYn]uhHuffjjji`@@\nfhiHR@DDhJQoAADILl{[NdijttuuUUP@@\nfhiH`@DXtX@lbbbbbQRRYHsfjijii`@@\nfhiH`@LTEhBDfWVUm]F|Ejf@@@@DJP\nfhiHa@LTdFB@A@`cIIBhhddmmNMYZ`hJJD@@\nfhiIB@DXHI`@aIg[ee_UoAZi`@@@@@\nfhiIP@DDhpDc^BdwJ|lzsfkMMMUUU@@@\nfhiIP@DXxHDc^CdTRbbeTVSNZltuUUUT@@\nfhiIP@DXxHDc^CdTRbfaTVUNZltuUUUL@@\nfhiI`@DXYH@`yEDeDbhdmScoMUJp@T@`J@\nfhiP@@@^Rfuue]Yrsj@B`@h@@@\nfhiPA@FV@DTBDiee~yqJ|F`hHB@@@@\nfhiPP@DB@SblbfdtQfbRHspVjff@B@@@\nfhiPR@AFIsPPRBSL{sNiFV]UUULtt@@\nfhiPR@JLIIpPRBSLzmj{NFmRuSUUT@@\nfhiP`@DZAyIgVYW^VkNZjdHBh@@\nfhiPb@OA``@QddabbRRvRkF\\m@@@E@@@\nfhiPp@DX@IQoARYeYv}YsUfjffjj`@@\nfhiQ@@DX@drlrloks`mTEATQ@@@\nfhiQ@@INAdTTTbRarrk^BtEP@@@@@\nfhiQP@DN@Qg@iLoJjnbdZsUUPAR@@@\nfhiQP@DX@IW`yEEDeDdbdsakMUMMUU@@@\nfhiQR@IV`Sa``dLbRbTQbbvXIwjjjjYf`@@\nfhiQR@JLIIW``BDfYu[UV\\MZejfjjh@@\nfhiQR@JLIIW``kDfYu[UV\\MZejfjjh@@\nfhiQR@LVXhw``BLdTVbRbbqqRVjYjjZj`@@\nfhiQR@LVXhw``QLdTVbRbbqqRVjYjjZj`@@\nfhiX@@@PHJpTiZUWYz\\MyjjbAb@DCH\nfhi`P@L@yHrTrlm{}AN}UUTDAP@@\nfhi`P@O@YivQQQRIV[SDXhsPUADu@@@\nfhi``@A@|eMkkJzsegT@E@AP@@@\nfhi`p@D@Dhw`iLrl{^lyjsUKUUUP@@\nfhia@@E@rQQQHyJIHToAhJ@h@@@@\nfhia@@F@rJJIIJZIDRoAjBBB`@@@\nfhiaB@ARADILk[OZtYxL@EM@@@@@\nfhiaB@K^@DISLjo{Bthu@@@EL@@@\nfhih@@@XIQRTsNj\x7FKNRtA@AEU@@@\nfhii`@DTxIPGArJIQQJiJYgIJZiYjjj@@\nfhii`@DTxIPK^rJIQQHpsQgEZZPhBjh@@\nfhipP@DXx@w`yEDhihUEePVkMMUSUS@@@\nfhipP@DXxBW`iLlwZztDZsSUUUTp@@\nfhipR@BTYIW``qDfY][UV\\MZejfjjh@@\nfhipS@IZCpSo@bBAA`cHhheHeTiFJlFBAXiVH@@\nfhipS@IZCpSo@bBAA`cHhheHeTiFJlFBAXjVH@@\nfhipb@AZL{pLQIe]e]ncEZ`@iZfb@@\nfhipp@DXH@Rc^BdsJsmzsfkMUMMUU@@@\nfhiq@@DZBCHheEDeDceNmyj@@`B@@@\nfhiqB@IF]hDadTTTbfLVWE^CMUTmRt@@\nfhiqP@DXH@RoArJJIJIIEIgCVZjZZjj@@\nfhiqP@DXxBQoArJIQSPjKJ`mVZZjjjf@@\nfhiqP@DXxBVoARYYnuuhHuffjjji`@@\nfhiq`@DXX@SdfUgyVj\\lyjiZjjh@@\nfhq@`@@AJDaHRDRHrHPVB`Hbjj`@@\nfhq@`@@JRVUeZZYq`@jjj@@@\nfhq@p@@HdkQlbbTLTTvRHqjZfZj`@@\nfhq@p@@LdkUlbbbbVRLrxIijjif`P]@\nfhqA``BZLBHyxRbEAFQQQIIRUILEAADuU@@@\nfhqH`@L\\MXBLbbTTJbbLTYjh@bj@@@\nfhqP@@@XRe]UUUpQijjjj`@@\nfhqQ@@DT@drlsLjXKUAQEP@@@\nfhqQa@EVEHpHc`QdTRVTaTrUFP@RtuP@@\nfhqXB@J\\dZpPAFRIJIJJqIBeUSMMUpHJ`\nfhuA@@@IRlrjwPQSbgFC[u@@@@@@@@@\nfhy@@@LdbbbTRrNYHVoA@@```@@@@\nfhy@B@@XbRsLrj\x7FIZ}XM@@@A@@@@\nfhy@`@@TrQQRJKSE[deV|F`jH@@`@@\nfhy@`@@\\RYee~uYrVoAZA``@H@@@\nfhyA@`A@bBQGh`LPdsNjjjsdeV]A@U@E@@@@\nfhyA@`F@bMQCh`LQdTTTTQbVSN\\U[tDPAAE@@@\nfhyA@`G@bEQEh`LPdsKOKrrQoN}A@qSUE@@@\nfhyA`@@BMdTTTTTTVoEZ|xL@A@P@@@@\nfhyA`@@BMdTTTTTTVoMJ|xKPAAPA@@@\nfhyA`@@B|dsLsKnqVgVBt@Q@BP@@@\nfhyA`@@DtdsMksjpRc^CA@@P@@@@@\nfhyA`@@Hldrlk[oQbcN|uHPQUP`@@\nfhyH@@@RM[IEDhmBdeJ\\DkQ`@@ajhH@@\nfhyH@@@RlyJYY{W^gAJMX@@HfjB@@\nfhyH@@@Rl{IEDhmCDcJ\\DkU`@@bjhD@@\nfhyH@@@XEkIEDhheYdjRXKpP`jJ@@P@@\nfhyH@@@XxkIEDeDehTjBYspP`Xj@@`@@\nfhyHB@FR|DDPdrmkZkSdeV]@@@@E@@@@\nfhyI@@DTeXBSOJwJxeZmxKT@QUJL@@\nfhyI@@H\\EhFQJJIPiIQ\\dkS`tADMSIP@@\nfhyI@@LDmxFRJIHqQJKTXYuoTuPDPtP@@\nfhyI`@LJMxD`yHhhhUEedRfcN}UUTBAT@@@\nfhyP@@@VRe]VYWISUoAjjjP`J@@@\nfhyP`@BZ@aInUe{W`mF|Ejjj`@X@@\nfhy`@@@ISLjm{btjw`t@@@@@@@@\nfhya@@O@RYVyY\x7FiqQkN`@@@B`@@@\nfhyh@@@\\eYvRJJJJUKKTxkQk@AA@uPP@@\nfhyh@@@\\e[vRJJJJUIITxkQk@AA@uPP@@\nfhyh`@LJMxHIYHhhhUEedRfcN}UUTBAT@@@\nfhyi@@HXIipDefWVuj\\TYwijh@@H@DHH\nfhyp@@@XheLwKjjsde^BDPUPP@@@@\nfhyp@@@^CdbbbbRRqUNRmip@@@PT@@@\nfhyq@@DTTAIgeUUTTmF|Ej@Bh@@@@\nfkAA`@@TTeJwsLDDXkVcUUUUUP@@\nfle@P@@HlgHhdeECDhl[^\\EjjifffP@@\nfle@`@@IrJJIERJEJYgOA`HJh@J@@@\nfleA@@@YEEMDleHiJJQoM@PUUUD@@\nfleA`@@HhdskLl{~cNBuT@@DR@@@\nfleA`@@QSdTTTRRTqfdcZ\\t@@AEL@@@\nfleHP@DTuDC@irJIJHqZIJJgOQZijjjYh@@\nfleH`@G^BdELdTtRRJbbQyHLZBBJ@``@@\nfleIP@DXxHDgQSdTRbffbTNUNV|tuUUUTp@@\nfleIb@LB\\fBK@BLddJbbRTVVuypUj@@HB`@@\nfleP@@@DrRJPqQPiIYtPijjh@@@@@@\nflePR@MA`cA@^FRRXqQJJMJJMjsU@QCUQ@@@\nflePq@JLIIShiAVVdLbbbrRfTbJYpwjVjZjjj@@\nflePr@FVQSaoIAHIM\\w[KoPS`uUUTsTt@@\nfleQP@EQ@RsoYEEEleDhiXRcNMUUULuT@@\nfleQa@AVQ``@bGadTbRTRRLrUgQRuUAP@@@@@\nfleQp@DXArVcQSdTTRbbtRjSAN|uUUMUUP@@\nfle`Q@LPQxH@bpQddebbTVVdleF\\m@@ABt@@@\nfle``@G@edbbTRbbNv`mV}PPE@QT@@@\nfle``@O@bdsJl{zpmV\\tD@@EL@@@\nflea@@E@RfVYUV\\tWIjh@`B`@@@\nflea`@I@U{HhhddihddQZbehHJJjh`@@\nflehPBDX}EHAJCFT}dTTRTRRJRsNFluTtuUKP@@\nflep@@@XleLl{rkMA^RDEEAAP@@@\nflepR@F\\EhtT`qDeYeumUEIwiiZijjjABD@\nflep`@IJR@bSLj~rnlExKTmS@AT@aJ@\nflepr@DDhJRgQRBHRYYym[teK^fffjjjh@@\nflepr@M^SAFB\\\\CprRSEQIZJKZR`qZjjjfVh@@\nfleqQ@BTYIShiAVVdLbbbVRfTbJYpwjVjZjjj@@\nfli@`@@HRYUoV]YjZfZjhFBPiT@\nfliAP`BJt{p@bG``HpDYEHdeCEeEluUTuMT@@\nflm@@@DjYee]{YhTkA|d@@@@@`@@@\nflm@@@DjYee]\x7FYhTkA}D@@@@@H@@@\nflm@@@DjYeeo\x7FyHQm^bd@@@B@@@@@\nflm@@@DjYee\x7F[yHQmNBd@@@@@`@@@\nflm@@@DjYee\x7F]yHQmNCD@@@@@H@@@\nflm@@@DjYee\x7F^yHQmNCD@@@@@H@@@\nflm@@@DjYee\x7F_yHQmNCD@@@@@H@@@\nflm@@@LdbbbTRRJrYHVkAbd@@B@@@@@@\nflm@@@LdbbbTRfvQyHQm^CD@@@BbhA@@\nflm@@@LdbbbTRrJRxhQmNBd@@@@@`@@@\nflm@@@LdbbbTRrVQyHQm^bd@@@@`@@@@\nflm@@@LdbbbTRraQxhQmNbd@@@@H@@@@\nflm@@@LdbbbTRvRQyHQmNbd@@@@H@@@@\nflm@@@LdbbbTVRcRXhQmNbd@@@B@@@@@\nflm@@@LdbbbbbT^VxhUc^CD@@H@bhA@@\nflm@@@LdbbbbbT^VxhUc^CD@@H@bhB@@\nflm@@@LdbbbbbVNVxhUc^bd@@H@@@@@@\nflm@@@LdbbbbbqeRYhUcNbd@@H@@@@@@\nflm@@@LdbbbbbqvRYhUcNbd@@H@@@@@@\nflm@B@@AFRJJJIISKHdiHvoQSP@P@PA@@@@\nflm@B@@QFQQQQQIJGKLdjqoQS@@@@@@@@@@\nflm@`@@QRYfYUgydeNLGrV`@@@BH@@@\nflmA@@@ILkZvmkQdeV]EL@@@@@@@@@\nflmA@@@ILroZrnQ`iFmEKPP@@PD@@@\nflmA@@@ILroZvjQ`iFlFKPP@@@T@@@\nflmA@@@ILroZvlQ`iFlFKPP@@@T@@@\nflmA@@@ILroZvnQ`iFlFKPP@@@T@@@\nflmA@@@ILrrknjsdcV]EL@@@@@@@@@\nflmA@@@ILrsssoQdmFmEKPP@@PD@@@\nflmA@@@ILrszmkQdeV]EKPQ@@@D@@@\nflmA@@@YEDeMDhmBeSbcV]EL@@@@@@@@@\nflmA@@@YEEDeieEBdQ`iF]EKPP@D@D@@@\nflmA@@@YEEDhihcBdqTcV]EL@@@@@@@@@\nflmA@@@YEEDhmBeLdpTeV]EKPQ@@@D@@@\nflmA@@@YEEDhmBeLdpTeV]EL@@@@@@@@@\nflmA@@@YEELeDdeBepQgV}EMP@@@@@@@@\nflmPB@JQAbIKJr|jlsfk^RFMUUTDCED@@\nflmQb@LA`bp@cIICMDeMDUARd[S`qVjZjZfi`@@\nflm`@@@ILkZvmkQdeV]EL@@@@@@@@@\nflm`@@@YHhhhhdddkSdgFbEH@@DDTHD@@\nflm`@@@YHhhhhdhb]RRkFCzH@@@@@P@@@\nflm`@@@YHhhhhdhbmRRkFCyH@@@@A@@@@\nflm`@@@YHhhhhdhecRRkFcyH@@@@D@@@@\nflm`@@@YHhhhhdhecRRkFcyH@A@ALPD@@\nflm`@@@YHhhhhdhecRRkFcyH@PAHE@P@@\nflm`@@@YHhhhhdhecRRkFcyH@PAPD`P@@\nflm`@@@YHhhhheED[RRgFbEH@@@D@@@@@\nflm``@J@SdbbbbLTtQqJRl[pXtC@@@AP@@@\nflm``@K@PdrlkZk]FJLxLTmRDDUUQP@@\nflma@@@`rJIJZIQPiZgEZmxLX@@@BI`H@@\nflma@@E@RfVU^_VcAZLzJZB`@H@@@@@\nflmaA@KQ@DPdLdtTRbtJfRyIVkNbfifZiZej@@\nflma`@E@PGHheEDdhlbj\\TxwhiZ`@@`j`P@@\nflmp@@@VbeJvsOviJBUYtT`@AEP@A@@@\nflmq@@DVbAIee^WmZ\\djshiZ`PJ`@b@@@\nflu@B@@YbULsJjozCZOdu@@@@@@@@@\nflu@P@@BBgHhhhhhehlVJtztV`BH@dh@@@\nflu@P@@BLGHhhhhhhlcnJuytV`BH`Dh@@@\nflu@P@@B|gHheIhcEDiRLu[pVfjjjjZ`@@\nflu@P@@HL{HheEhhdhh~BlGtVijjijj`@@\nflu@P@@\\\\eIfVWyWv\\ehLV`hHB@h@@@\nflu@P@@\\uyIfVWnuv\\dZJV`XJ@`H@@@\nflu@`@@IRYfUV}~QmNCEh@@@@H@@@\nfluAA@A@bDqFQQIYRZJIGTYYpXp@UPTAD@@@\nfluAA`A@bBQCHatQFHSDYEEEMBdeEerRm^CADQPPT@@@\nfluAP@@B\\ENQQQQQQKQXlUkuhm@DPDIP@@@\nfluAP@@B\\FJSLsLoJlUkudm@DPDLP@@@\nfluAP@@HdkrSKLvzjlxYtTmT`pQUP`@@\nfluA`@@FmdbfbTTRbawEZ|EM@@@TA@@@@\nfluA`@@LbdsLn}kkNRuxKPTD@AT@@@\nfluH@D@\\|DDiSLrno}SdcZ\\@@@AET@@@\nfluHP@CNbdDeZrQQQQIJJkHeKNCFjjj`PJ`@@\nfluHP@DLtdGGArJJIQRzJJY`iJbejijjYjh@@\nfluHP@DTuDC@irJIJHqZIJJgIJ}Ejfjjifh@@\nfluH`@DXYh@lbbRbTVaRRisQoAZZfd@Jh@@\nfluH`@DXlDGlbbbLTRvTtISQhiZj@`ijb@@\nfluHb@BJtFW`QFQQQIZKQFZLDZtTpPPuTuE@@@\nfluHb@LLxxLPABTnrjkmxHJQkTlu@EMP@@@\nfluI@@DRUXFQQZJIPiIIRtZpXmPDPuKE@@@\nfluI@@DTXDFQQZIQKQQH|DYtTmSU@A@@@@@\nfluI@@D\\tXBSOLkKzrthpXmPPIUKE@@@\nfluIP@DTxHFgQSdTRbRJbTVSNBe[sSMUUUTp@@\nfluIP@LF]DDeZSdbbbbRTUVQJV\\FMUUU@`U@@@\nfluI`@FNbdAMILs[[N}SagAcUUUPHEP@@\nfluPB@AJADYEDeMDeELcSfgQS@DMPT@P@@\nfluPP@B\\@bslbfbTTRTrJDkQhiZjfjX@`@@\nfluPP@DTAsUlbbTTTlTvRXHshiZjZfX@`@@\nfluPR@JLHkU`RBSL|kMjlyxLTtmUMUUU@@@\nfluP`@DIAEIeUYU}jLehufiZYfij`@@\nfluP`@DX@qIeVyeUrT]DJVjZ@@@H@@@\nfluQB@AJbBHRYVyV]VgMNBf@HZ``P`@@\nfluXP@L\\MZJPRUkIEEEDdeFdbTLxLZjjjP@j@@@\nflu``@C@PdwLl|{XmNbEKPP@@@@@@@\nflu`a@BJ|CDebYEEEheDXdlpQkQSAADtB@P@@\nflua`@A@\\EIfYfWgvJuztV`BHBDh@@@\nflua`@O@TyIfYfUYvBtFrV`BHAbH@@@\nfluab@FF|D@QdbbTVbtjfQJBUjuRsUMMTh@@\nfluhP@DTxIPCQSdTRbbUTRRSNRT{sUKMUUUP@@\nfluhP@LF]EHIJudbbbbRTUVQJV\\FMUUU@`U@@@\nfluhP@NFmEHIJudbbbbRRcRQJF\\FMUSUP@U@@@\nfluib@DTxIQhi@DYEDhheUDddsdeN|uRsUUUT@@\nfluqB@AJMDFHdrmrlzmNZ\\EL@PuA@a@@@\nfluq`@LLxAGdfWeUWnEHuhiZij@IjB@@\nflux@@@XIRRmYHhhldiUDdsdcQRA@AQ@@@@@\nflyAP@@HLxJSJkZrkzluT`DUT@@@\nflyP`@EN@cHheDdehbdiZZ``bjj@@@\nflyPa@LR``@PHHrRPjJJKHiXsejAbHj`@@\nflyQ`@D^@jrSOJzmkBdu@AUUT@@@\nflyQ`@D^@jrS[JzmkBdu@AUUT@@@\nflyQa@ANQ``@bGadTbRTRRLVegMUPTAU@@@\nflyQa@LR`cp@`PQddaTTRrrRVkKTCAEU@@@\nflyQb@LR``P@cIIBhhhlbecNVhFHbj@@@\nfly`R@IB]zH@aJYYY]gcAjhJ`fF@@@\nfly`R@JRU[p@cIDeDeEECEXijji`Z@@@@\nfly`R@LPQEH@cIICEDddeUC^VjB`hf@@@\nflyqa@AFlcA@ADOCHiDdhThmLmFZj`@jj@@@\nfl}@@@LdbbbTVtQVXhQmV]zJZPZHPPJH@@\nfl}@B@@AFRJJIIJISXxhKRoF]EM@@@@@@@@@@\nfl}A@@@ILrjvjvQ`mJLxLTmA@@@@@@@@\nfl}A@@@IRlkjmkAENRUYtTuPH@@@T@@@\nfl}A@@@IRlkvjkAENRUYtTuPP@@@T@@@\nfl}A@@@IRlkvmkAENRUYtTuPP@@@T@@@\nfl}A@@@IRlrjkoAENJl[tTt@@@@@@@@@\nfl}A@@@YHdhdmMDhcAENRt{tTt@@@@@@@@@\nfl}A@@@YHdhheDddcAENJ]ZHTt@@@@@@@@@\nfl}A@@@YHhhhddjdeAGERuYtTuUAL@@U@@@\nfl}`@@@IRlrjkoAENJl[tTuP@@@@T@@@\nfl}a@@A@rQIQJKJZEJBF\\UXshij``@@@h@@@\nfl}a@@E@ReYWm[VBJ\\ejshij``@@@h@@@\nfl}a@@M@ReYYy}VBF\\TZshij``@@@h@@@\nfl}aA@HF@DIdLdRbbRRtJR``cIFuXLZjjZjeif@@\nfoA@@@DjYeUlVeN`@`B`@@@\nfoA@C@@QAHadQdTTTqTRWAV]ADD@@@@@\nfoA@C@@QBhatQdTTRbRQsIFmA@A@P@@@\nfoA@P@@HTYIe[VWhrVfi@Hj`@@\nfoAA@@@ILkjrmFV]@AL@@@ar@\nfoAAB@A@bDfUmyVcKN`@j@@@@@\nfoAAP@@HiIRSKSrvSeNMUTuUP@@\nfoAA`@@HXdrl{ZlyYsUUP@P@@\nfoAAa@JJtBHxDYDhhmEDcpUgUKPJA@@@\nfoAHB@MFlx@QdbdRjbbcNZmUUTmU@@@\nfoAH`@D\\DXDdfVYe_ScVZ`ZBP`@@\nfoAH`BDRLx@cR{dfYeW\\cMVZhFAF@@@\nfoAHb@DXx[U`QFQQJIZJHdiHtuUTsL@@\nfoAI@@MFlxFRJQJjJJLyjuUURuT@@\nfoAP@@@NRYWVUzLMZ@B`@`@@\nfoAPB@KN@DISLjohmJMP@@AP@@@\nfoAPB`ARADDbCREiAt`dsNkmLyITDASUD@@\nfoAPR@JLIIU`RBSLzmlsagTmTuUP@@\nfoAP`@BZ@aInvYWejsfjiB@`@@\nfoAP`@CNAQJYW]YqS`jB@jj`@@\nfoAPa@LDlx@P`HRYYeUuVLyjjjjj@@\nfoAQA@BLD@HpDISOmrnJV]TuL@D@@@\nfoAX@`LTXKU`AB@aPPxHRYWUeRTMYjjjZV@@\nfoA`@@@ILkjrmFV]@AT@@@@@\nfoA`@@@IML|{wEFmUSA@R@@@\nfoA`R@EPQap@cIIKEDhliZMYZ@@@`@@\nfoA``@M@PdvkLkruYsUTa@P@@\nfoAaB@GDADILroZ{AFmAA@@P@@@\nfoAaB@GLADILk_J}NFm@DP@P@@@\nfoAaB@G\\ADILkkJ}FFm@AP@P@@@\nfoAaB@KN@DISLjoxmJMP@@AP@@@\nfoAaP@D@hIVlbbTbfbjIrgFjjfjh@@\nfoAa`@M@PQInvYWejsffiB@`@@\nfoAh`@JBUYpBILoLkleFmUU@`T@@@\nfoAp@@@P\\eKLjorMjsP@@A@B@H\nfoAp`@BZLABS]lroKUgMURDA@@@\nfoApb@HDihp@cIDhideDRLUZi``Xh@@\nfoAqB@JLtPFHdsM|klEYuKTa@P@@\nfoAq`@DXxBSlbbTTtJVhKQffjjjX@@\nfoI@@@LdbbRbaQXKbiZlx@@@`@`@@\nfoIA@@@IRlkZ|DTyKUgUSUUJs@@@\nfoIA@@@IRlkjlDTyIUgUPP@AP@@@\nfoIA@@@IRlkvlDTyIUgT@@@@@@@@\nfoIA@@@IRlrj|DTxjqgUP@@AP@@@\nfoI`@@@IRlrj|DTxjqgUP@@AP@@@\nfoQ@@@DjYUU}`nRlx@@@`@@@@\nfoQ@@@LdbbbRQrYHRkN@@`@@@@@\nfoQ@@@LdbbbTRQYHUmN@@@B@@@@\nfoQ@@@LdbbbTReyHQmN@@@@`@@@\nfoQ@B@@QBSJvjkQdeV]@@@@@@@@\nfoQ@B@@QBSLljzsdcV]@@@@@@@@\nfoQ@B@@QBSLsJ~pRkF]@@@@@@@@\nfoQ@B@@RBSLljzsdcV]@@@@@@@@\nfoQ@B@@XbSLljzsdcV]@@@@@@@@\nfoQ@`@@VRfYeUzBd[Sj@B@H`@@\nfoQA@@@ILrrkkNRMYt@@@@@@@@\nfoQA@@@YEEEEDdbrRkF\\t@D@D@@@\nfoQA`@@HldrlkZtXhugMRDDSD@@\nfoQA`@@HldrmlrtYKUgMR@QSD@@\nfoQH@@@BMYJUme]IRRmF@@HFhH@@\nfoQH@@@BMYJYYo]irRmF@@HFhH@@\nfoQH@@@RMYJUm[UIPRmF@@BFhH@@\nfoQH@@@RM[IEDhmBdj\\DkQ`@@ajB@@\nfoQH@@@Rl{IEDhmCDj\\DkU`@@bjA@@\nfoQH@@@TUkIEEEhddV\\Djs`JBjeF@@\nfoQH@@@TuYJYfU}YrQeZ@H@XhH@@\nfoQH@@@XD[IEDhiedZRXKU``jIfF@@\nfoQH@@@XDiJYYe_iI`cNBBhidX@@\nfoQH@@@XhiJYmg]YpRcVBHJiXX@@\nfoQH@@@XhiJYmg_YpRcVBHJiXX@@\nfoQH@@@XhiJYmggYrRkNBHJfhT@@\nfoQH@@@XhkIEELhddV\\Djs`bBjeF@@\nfoQH@@@XikIEEMEDUf\\ejs`bBjeF@@\nfoQI@@HTEHBRrjjjsdeV]P@TtiP@@\nfoQP@@@FRfYeUz\\e[SjdF`hd@@\nfoQP@@@NrQQJIIYFcAZMX@@@@`@@@\nfoQP@@@NrQQQQIHz`iFuX@@@@`@@@\nfoQP@@@RrQQJIJUJgAJlx@@H@@@@@\nfoQP@@@TRfY^]|TEjsiB`hjb@@\nfoQP@@@TrQQQQPqUdeZlx@`@@@@@@\nfoQP@@@VRfYe]Z\\d[S`@@@H@@@\nfoQP@@@XRf[WmZ\\djs``@@@@@@\nfoQP@@@XRf[eUz\\UXs``@@@@@@\nfoQQ@BHB@ChRdsLrolDjqgP@@@@@@@\nfoQ`@@@ILkZjmFRUYt@@@@@@@@\nfoQ`@@@IRlrj|DJlYuT@@@P@@@\nfoQ`@@@ISKON}NBeYp@@A@@@@@\nfoQ`@@@ISLro}NRMjp@@@A@@@@\nfoQ`@@@ITsJkyFJlYuPA@@P@@@\nfoQ`@@@ITs\\k{AJlYuT`tDR@@@\nfoQ`@@@YHheEMDcPTeF\\@@@A@@@@\nfoQ`A@F@bIABSJvmkSdeV]@DT@A@@@\nfoQ``@B@mdTTTTRQSEBdYsTmUSU@@@\nfoQ``@M@HdsMjoLDhugMLBDUD@@\nfoQa@@A@RYfYWvBUXsf`B@@`@@\nfoQh@@@LDZvRJIQPiYTYIVc@@DCTD@@\nfoQh@@@TejrTsLkzsdcJt@P@qPP@@\nfoQh@@@TtZvRJJJIIELEHrm@D@LTD@@\nfoQh@@@TyjvRJJJIIELEHrm@D@LTD@@\nfoQh@@@XijvRJJIQIELEHrm@D@LTD@@\nfoQp@@@XDeLvzjtyIUgAAEP@@@@\nfoQp@@@XDeLzmjtXJQkA@T@DP@@\nfoQp@@@XDeLzmktXJQkA@T@DP@@\nfoQp@@@XHeLsnjxiIUgMPP@A@@@\nfoQp@@@XdeL|oZtYIUgA@U@@P@@\nfoQp@@@XidbbfbQRSNRUYpQAT@@@@@\ngBQ@@eJuT@@\ngBX@@eLUT@@\ngBXHHHaIejh@\ngC`@Die@@\ngC`@Dkz@@\ngC`DADHHRVh@\ngC`DADZHRVXRH\ngC`DADZHRVXRP\ngC`DADZHRVh@\ngC`DADZHRVx@\ngC`DAHJPRZd@\ngC`DAHJPRZh@\ngC`DAHZPRVh@\ngC`DAVRDRZh@\ngC`DAbKDRZh@\ngC`DAbZHRVh@\ngC`DAb[DRVx@\ngC`HADIKLIH\ngC`HADIKR@@\ngC`HAbIKJ@@\ngC`HAbIKLIH\ngC`HAbIML@@\ngC`LAHJPt`duP@\ngC`LAVJluXduP@\ngC`LAbKDvHduP@\ngC``@dfZ@@\ngC``AdeY@@\ngC``AdeZ@@\ngC``Adej@@\ngC``Adij@@\ngCa@@dkH@\ngCa@@dkP@\ngCa@@dkX@\ngCa@@dmH@\ngCa@@dmP@\ngCa@@dmX@\ngCa@@duP@\ngCaHDHaIZ`@\ngCaHH@bNt@@\ngCaHL@aIZ`@\ngCaHLHaIZP@\ngCaHLHaIZ`@\ngCaHLLQIZP@\ngCaHLLQIZ`@\ngCah@mJAIj`@\ngCahHl@bNt@@\ngCahHlGBNt@@\ngCahHlHRNj@@\ngCahhlAa]ncm@@\ngCd@Adej@@\ngCdDI`BHDRZh@\ngCe@E`dmH@\ngCe@E`dsP@\ngCe@I`dlpd`\ngCh@@dmH@\ngCh@@doH@\ngCh@@doP@\ngChHD@aIU`@\ngChHLHaIZp@\ngChHLLQIZ`@\ngChHLLQIZp@\ngCh`LDdsP@\ngChhMDHbNl@@\ngCi@DDeV@@\ngCi@DDeZ@@\ngCi@DDfZ@@\ngCi@HDefDd@\ngCi@LDej@@\ngCiHEAxIVt@@\ngFp@DiTt@@@\ngFp@DiTvjh@\ngFpD@DTHReSZj`@\ngFpD@DXHReSZj`@\ngFp`@dfTujX@\ngFp`@dfTujh@\ngFp`@df_Ejh@\ngFp`ATiTvjh@\ngFp`AdiTvjh@\ngFq@@drfmM@@\ngFq@@eLqUU@@\ngFq@@eLzts@@\ngFqDDHbqBSZmUT@@\ngFqDLLRqBRwcUT@@\ngFq`AbeJfuU@@\ngFr@ACTiTt@@@\ngFr@ACTi[FZd@\ngFt@ATiTt@@@\ngFt@AdiTt@@@\ngFtHCQbILimKP@\ngFtHE`DILikUP@\ngFtHLPDISNmLp@\ngFu@E`drfmT`@\ngFx@@eJf`@@@\ngFx@@eJfuU@@\ngFx@@eRfuK@@\ngFxHJ@aJ]bjn@@\ngFxHL@aJYujj@@\ngFx`HJeJxtm@@\ngFy@DDfTujh@\ngFy@JDiTvjh@\ngFy@LDeXvjh@\ngF|@AbeJf`@@@\ngF|HLZ@aJYuif@@\ngGP@DiVj`@\ngGP@Di]jQ@`\ngGPBABRHTQXcIHUTpf`\ngGPDADVHRUZZDT@\ngGPDAbWDRUZj@@\ngGPP@cTfyi`@\ngGPPACTiVj`@\ngGPXHlQxIU[U@@\ngGP`@TfYi`@\ngGP`@df]j`@\ngGP`@dfuiaM@\ngGP`@dfuj`@\ngGP`@dfyi`@\ngGP`ADkjj`@\ngGP`ATeUfQD`\ngGP`ATeUfaE@\ngGP`ATeVj`@\ngGP`ATfVi`@\ngGP`ATiVj`@\ngGPhCQDIJuS@@\ngGPhH`DYIHUi@@\ngGPhMPDIK]U@@\ngGPhMQDIK]U@@\ngGQ@@djltHh\ngGQ@@djuT@@\ngGQ@@dlmT@@\ngGQ@@drmT@@\ngGQ@@dsMT@@\ngGQ@@dsmT@@\ngGQ@@eMUT@@\ngGQDB@jQBUkUP@\ngGQDL@aAFQRFj`@\ngGQDLHbqBRwSP@\ngGQDLHbqBRwUP@\ngGQHBHaIeiXQH\ngGQHBLQIeiXQH\ngGQHJ@aJ]jd@\ngGQHJHaIUjh@\ngGQLJIARFdLbdMU@@\ngGQXHlZHROjj@@\ngGQ`@ZdrmJ@@\ngGQ`@ZdrmT@@\ngGQ`@ZdruT@@\ngGQ`@bdvmT@@\ngGQ`@bdwMT@@\ngGQ`AjdmmR@@\ngGQ`AjdmmT@@\ngGQdEb@bRFRRVV`@\ngGQhHlLSIHTmP@\ngGR@aEQJyJ]Vh@\ngGT@ADiVj`@\ngGT@ATeVj`@\ngGT@ATeWfqE@\ngGT`E`TfUj`@\ngGT`E`TfYj`@\ngGT`HQTeVY`iZ@\ngGU@E`drmT@@\ngGU@E`dsmT@@\ngGU`E`ZdrmT@@\ngGXHD@aIUVd@\ngGXXKEb@cIIBmP@\ngGX`BDdvmT@@\ngGX`BDdwMT@@\ngGX`LDdsmT@@\ngGXhKB@aJvZX@\ngGY@DDfYj`@\ngGY@HDeVZaI@\ngGY@HDefZaH`\ngGY@HDfVj`@\ngGYDIQDJHR[Zj@@\ngGYHEAxIVkU@@\ngGYHEAxIVsU@@\ngGY`BETfuj`@\ngG\\DEj@`aBSJuP@\ngG]HHjPDIJuS@@\ngGp`AldTLRuUT@@\ngGq@@eJqFuUP@\ngJP@DjYd@\ngJPBAVJPt`YAIjj@@\ngJPDADFHRUfaH`\ngJPDADFHRUfaI@\ngJPDADFHRUj`@\ngJPDADFHR[f`@\ngJPDADQpRZj`@\ngJPDAbGDRUj`@\ngJPD`DHHGKdfzd@\ngJPH@EIRuP@\ngJPHADIJsHd`\ngJPHADIJsPd`\ngJPHADIJtpb`\ngJPHADIJuP@\ngJPHADIKSP@\ngJPHADILth@\ngJPHAVILuP@\ngJPHAVIMUP@\ngJPHAbIJuP@\ngJPHAbIKUP@\ngJPLADZPL`dmM@@\ngJPLAHJPt`duU@@\ngJPXHlPDQzt@@\ngJPXHlQBQ{T@@\ngJPXHlQxQ{T@@\ngJP`@TeVd@\ngJP`@TeVh@\ngJP`@TeZh@\ngJP`@TfVd@\ngJP`@TfZh@\ngJP`@dfVh@\ngJPdEaDPHRZe`@\ngJQ@@djsBJ@\ngJQ@@dju@@\ngJQ@@dkSBJ@\ngJQ@@dkU@@\ngJQ@@dls@@\ngJQ@@dmU@@\ngJQ@@dru@@\ngJQ@@dsT`@\ngJQ@@duU@@\ngJQ@@eKU@@\ngJQDDH`qBRmT@@\ngJQDHG@nBUMT@@\ngJQHBHaIVYDd@\ngJQHBHaIVj@@\ngJQHBHaIfe@@\ngJQHBHaInZ@@\ngJQHBLQIVj@@\ngJQHBOAJfj@@\ngJQ`@bdvt`@\ngJQhHl@bOV`@\ngJQhHlHROZ`@\ngJT@@Te^h@\ngJT`E`TfVh@\ngJU@DPdru@@\ngJU@E`dru@@\ngJU@LPdkU@@\ngJUHLT@aIjZ@@\ngJX@@dkU@@\ngJX@@eKU@@\ngJXDBLQXbS]V@@\ngJXHD@aIUZ@@\ngJX`BDdvu@@\ngJX`DBdjt`@\ngJX`DBdru@@\ngJX`LDdru@@\ngJY@BDfZl@\ngJY@DDeUh@\ngJY@DDeVh@\ngJY@DDefh@\ngJY@DDfVh@\ngJY@LDijX@\ngJYHC`DIKTp@\ngJYHCaHIKTp@\ngJYHEAxIVmP@\ngJ\\@@ldru@@\ngJ]@EbDfVh@\ngJ]HEcAxQ{T@@\ngKP`@df\\Vj@@\ngKP`Adi\\Zj@@\ngKQ@@eKEUP@\ngKQ@@eKcRp@\ngKQ@@eKcUP@\ngKQHDHaImZj`@\ngKQHH@aJoEj`@\ngKT@Adi\\Vf@@\ngKX@@eKcRp@\ngKX@@eKcUP@\ngKX@@eMEUP@\ngKXHL@aJWFe`@\ngKXHL@aJWFj`@\ngNpH@DIRkUJ@@\ngNpJAHJPtaYArFQRFUU@@\ngNpN@iRHTQh`qEbGDQ\x7FUT@@\ngNpOAYRHTQh`qEbGDAh\x7Fjj@@\ngNpP@jtfvZf@@\ngNp`@TfWZZ@@\ngNp`@dfVZf@@\ngNp`ATeVYjDT@\ngNp`ATf^Zf@@\ngNphJqDIKMTl@@\ngNphMQLISMSL@@\ngNpiJqbOdeVjf@@\ngNplJqHJPtadTaeTp@\ngNplLQDSpPPduluP@\ngNq@@djkLpRb`\ngNq@@djkLpRc@\ngNq@@djkMPb`\ngNq@@dl{UP@\ngNq@@dskUP@\ngNq@@eJmTh@\ngNqDLHaqBRjuU@@\ngNqHFHaIfYjP@\ngNqhHl@cIIBej`@\ngNqhHlOAJkVj`@\ngNqhMV@aI[jY`@\ngNt@@teUjj@@\ngNtHAePIKMuT@@\ngNtHE`DILruT@@\ngNtHE`DILzuR@@\ngNt`E`TfWZj@@\ngNt`E`TfYZj@@\ngNu@E`dl{Tpa`\ngNu@E`dskUH@\ngNu@LpdjmUP@\ngNx@@eLuUP@\ngNx@@eRmUP@\ngNxHD@aIUUj`@\ngNxHD@aIUej`@\ngNx`LDdskUH@\ngNx`LDdskUP@\ngNx`LDdssUP@\ngNxhMV@aI[ji`@\ngNyHEAxIVjuT@@\ngNy`BDtf{Zj@@\ngNy`LDtf]Zj@@\ngNyhGE`DYIITmT@@\ngN|@ADeJkUP``\ngOp@DiUkZV`@\ngOp@DjWkB@@@\ngOpHADILkW@@@@\ngOpHADILkWTu@@\ngOpHAbILkW@@@@\ngOp`@dfUMZj`@\ngOp`@dfVqZj`@\ngOp`@tiguif`@\ngOp`ATiekjj`@\ngOp`A\\dTQbjj`@\ngOphH`DYIHUVmT`@\ngOq@@drm[ST@@\ngOq@@drm[UT@@\ngOq@@drm\\@@@@\ngOq@@drm]UT@@\ngOq@@eLvmLt@@\ngOq`@VdsNkTl@@\ngOq`@fdrikTl@@\ngOq`@fdrikUL@@\ngOt@@TeVkzj`@\ngOt@@tiek@`@@\ngOt@@|daRbjj`@\ngOtHE`DILl[MT`@\ngOt`EPTfVKZj`@\ngOu@Dpdrm[Rt@@\ngOx@@drm]UT@@\ngOx@@eSRuUT@@\ngOx@AddJTUUT@@\ngOxHBHaIeZx@@@\ngOxHDHaIeZx@@@\ngOxHDIAIeZx@@@\ngOxHHHaIeZzjh@\ngOxHL@aJY}Zjh@\ngOy@DDfYKZj`@\ngOy@LDiisjj`@\ngOy@hAteIeZx@@@\ngOz@ACVeKNLuR@@',Lq='fHe`A@\neO`BNZ``\ngF|HLZ@aJYuif@pPWHP\nfHgdAy@\neMBBHRZCAKd`\ngNxDLHaqBRjuU@P\neMABHXaIhH\ngCaHLHaIZ`H\ngJQ@@dls@XKGd`\ngJP`AdejhC@qX|`@\ngGP`ATiVj`LCEcrT\ngGP`ATeVj`LJEc^R@\ngGQ`@jdjmTAal[rP\ngGQhHl@cIIJmPFBMyL\neFBCPcA@\ndid@p@bBbFbDfYoa`b@@LJ@fx^QP\ndaDD@@IIf]n``@@pkBiny@@\ngChDL@aABSM@XHKdp\ndazD`N``|DjwVjh@`\nfHbTAa@\ndaxB@@InRYgZj`CBdpf{dT\ngNxLL@aAABDfVZj@`\ngGPB@DHHpQPaIUZdB\neMC@HoABDe`pccH\ngC`L@DSpPPdkPD\ngJU@LPdjr`XIGd`\ndaG@@@kdig|jVj`CAdpj[`\ndaD@@DjUZxHH@CAdpj[nQ@\ngKP`@Tixjj@ptVyB\neMhDRZCAKd`\ngJQ@@dkU@XDQGdp\ngC`D@DXHRVhCAC\\TP\ngGQ@@djmTA`Xl^R`\ngGXXKEb@cIIBmPD\neF`BNFE@\ndidHPBBHFHRYgVzB@`@`\ngCe@E`dsPFBV@\ndaFL@HABR[e[fii@H\ndaE@@@yIe^f`@`@piLJny@@\ngO|HEfHaIeZx@@B\ndcN`@@pjYJYe}k`Hbj@B\ngN|HEb@cHhhVj`H\ngNp`@deUZj@pvkxhP\ndax`@@bnRn[ff`RHj\ngJT@@deVhCCKGbD@\ndigDPLXXP@b`cIIKDnEZfd@`\neFBBDcAaWH@\ndeVL@HAIR[e_aZ@B@B\ndaDH@@RVU[f@@@LJ`j[nQ`\ndcoH@DJ`RUeUVy]ZZ`@@B\ndmO@@@pldTTJRraej@B`@`\ngNp`@dfUZf@pvMyF\neMBCDRZCAKd`\ngJPDAbGDRUj`H\neMBBHRYCAKd`\ngCe@H`dtpD\ngCi@LDek@`\neMaDBKpRVB\nfHbXAa@\nfHfpAa@\neMIDBKpRYB\ngJ\\HEb@aIej@pPWEb\ndcNL@M@aRYgWVzB@j`@`\ngNpXHmQxYIIXuTA`Uj~P`\ndeT@P@bNBDfUuih@J`@pxDpjx^QH\ndidHPABHJHRYgVzB@`@`\ndcl`@@`nReeWZY]@``@@B\ngOr@Ajti]qZY`H\ndeV@@@ReWTj@Bj@CCdrfxYyD@\ndeV@`IBHRYWVf`@j@CC`SBkiy@`\ngChHLIAIZ`LEMrH\ngNxhMVIAI[ji`H\ndeWH@DJPRY[TYZ`@@B\ngGPhCPDILuS@XXR|f@\ngGP`ATivj`LMc^JH\ngNp`ATiUjj@pLQoEf\ngJQ@@djsBJptQyL\ngNu@H`dlkUPFBqxjp\ndeV@@@RYyTYj`@@CBPaLinFP\ngCah@mJAIj`LLE`\neMFIHRM`pgd@\ndaD@P@bNbDfUzZ@B@B\ndifL@LAARYYRijjX@`\ndaxD@@QIUYjj@LJpj[nQ@\ndcN@@@rQQHqIKmUP@@A`hPfEprn|b`\ndmtH@@RYWUih@Jh@LAALJnF^Ph\ndmv@@@rRJIIFUjB`@@pyLJnF^Hr`\ngF|@AbeJfuU@P\ndiFB@LANEIn^ZjX@`\ndeTL@@Z\\bbRbM]@DT@FC@fMt|`P\ndcNL@MAirJIJHsUt@QU@A@\ndeT`@@pjrQQIFTpDEP@P\nfHd`Aa@\neFJHbHpP\nfdu@@@LdbbTRrLtSI\\TZwhaDTUUTU@A`YFBTYpRmF\\GfI^b\ndiT@@DjYVnXPbD`@pIBinYgB\nf`iA@@@YHheDdYdj\\DiwjBH`@H@@`\nf`iA@@@YHheEdXdj\\EIwjBHH@H@@`\ndmLH`I@HRfUwrnF`PJ`@phBXUxbT\ndedF@@PfFTf{nZjf@LFSBii@\nfoQPP@DLAhulbbbTRTLYpTmNZjjjZh@H\nf`i@@@DjUg^uhIVg^BBb`@H@B\nflmA@@@ILkZvmsQdeV\\EL@EUUUQT@FEdXipTeZM[tT|P@H\ndet@@DjYUX^dHbH`CAdJfF^Ig@\ngGY@JDeejpH\nfb}`@@@YEDeMDhTihjLUiwhyZjjjjjjh@LOHpRgARUhugASgb@Q@\ndclH@@rQQRJItJ{PUPB@FEDwSoLhU@\ndg|D@@OIEEHhfmPkmAU@T`A@\ndcMB@HDDWTfyV{iZ@HX@H\nfak`A@H`dDrBSJ{LjnmkQc`m{@ATDDUUPP@P\neghP@@@LdbbdRdRQbLf|VRcEdbgfPRSUUUUUUUUMUP@D\ndidD@@iJ[gxZB@@C@dkaxaL\ndeTD@@QImYQejjj@LNpj[agdH\nfdu@@@LdbTRTjTtSA\\dZwhmUUUUUU@A@\nfoQ@B@@AFRJJJIQX|UJqgT@PP@`@D\nflmA@@@ILk[J}sQdmV\\EL@EUUUQT@D\nfgA`B@N@BDifYWz\\d[Uj@H@B@@H\ngGUHEj@aIUZdB\ngGX`LDdjmLIkBWba@\nfnk@`@@MrJJZQIFYJSKYdmF|DmVjjjjjjjf@B\nfoQA@@@YHhhheEcqTkF]PAA@B@@P\ndiFD@HCDemVZj`bJ\ndcMH@LDDeYWWajjjj@LIBDJfzUyF`\neFDBcAaWH@\nfewAP@@Vkv^QQISQJFKIQJVZcEF|EIVejjjjjijj@CBJLTxJRmFl{t\\ejOdI`\ndmTJ@@PYUIeUYjjV@LBpjWdL\ndmLH@@RYVuiiV@BiH@`\ndk\\D@@OIEEHhkje]hJbB`@`\nffsA`@@VkdTRTtRabrTRhqQoARUiZjijjjj@B\ndg|@@DjYV}~T{`d@@@@H\ndeTH@@Rfuunh@I@C@lJfybCH\nfoA@P@@HuYIe[UvhrRfi@Bf`@H\ndeTD@@QIgeQej@@@LJrfx^QH\ndif@`D@HRUe^Eh@@@pXDinGbPP\ndg~@pIBPRPrPrJPqISIBiwU@@@@@P\ndmL@@DjYe^TYHI`h@H\nfb}A@@@YEEHdcLeIeb\\djsdyZjjjjjjh@LOXIQgARUjsoAbgdL`\nfc\x7FAR@CQ[OXICHhhdiDimCEDYmBrPk^SfmSUUMUUUSTt@D\ndmLH@@Re^UpiVe`@`@`\ngNp`@deYZj@pNM_I@\ndmtH`NBHRYWVih@Jh@LNALJiWdB\ngJQHBLQIVi@`\ndmtD`NTHaIe]Vf`@jP@`\ndmtD`NTHaIe]Vf`@j`@`\ndklLPDp`BH|LbdRRaRtvjh@@@`\ndeU@@@aJWeQej@@B@h\ndg|`@@pjrQQQGQIMqvAPP@@@P\ndeT`@@pjrQQQUMpEAP@P\ndaED@DXNRYUifjf@H\ndmND@DCdfVUrjUZjZi@H\ndg}@@@aJVYU^Svv`@@@@`J\ndg\\@@LdbbbTicRDQIS@Qs@\ndk\\d@DsmB\\bbbTrQXUujjUj`C@liS@\ndazH@DAImeiiBL\\Lrfxb\\\ndkNH@DAImgVVffZPsEhV\ndefH@HAIUVYfZPrKJj\ndaF`@@pjYJYfn@b@@pPfyG@\nfhi@@@LdbbQfTbdxKShF@jjh`@`\ndmt@`@bDfUuZZ@Bj@CC`SBkaxdj\ngOt`DPtfWMZi`OLlP\ndaED@DHNRY[Jfjf@H\nfluP`@DD@iIfUeeurL]DJVjhD@@H@@`\ngNy@BDeUZj@pJu_EB\ndcNDpBWPbNBIdDfYUYa``bT@H\ndmvDPIxPbBBDfYZVVBAX@CB`rfgdJ\nfoA@@@DjYW{WdkN`XH@@@Ppp\nfluA`@@HddrsJoN}NF}EKTpT@@D@@XJEKQkN|FJ^PX\nfhi@P@@HeyIefU]zgCNZfB`ah@B\ndcML`EvDp@cIIKDedLkP@R@D\ndieD@DpARYevyjiX@pRnxfD\nfhyA`@@HBdrrkN~RdcN|uA@@AP@A@\nfbmQA@AJRBH}HILk\\kOKtyit\\pACTE@a@@P\ndg^B@HAEmInUvwaZ@Bj`@`\ndg^L@DAaRUf^uNvjj@@@H\ndcNH@DAIgfYgVhBH@CBnFUwbI`\nfle`a@BAbBHeDYEEEDhXXeDqWdpPU@PL@@D\ndefL`BpYBHRVvZfjX@`\nfoAa@@F@rJJJIJINUeFZBPjHP@pDRgARtYwHJ@\ngOtHDpDILrWST`XHwbF@\ndigD@HK`QInW[fiZ`CCDX^QP\ndeWD`D[ad@cHhhhXWTuSHA@\ndg}@@@mJYeU|]Tz@@@H@CAPpj[aeSyE`\ndiU@@@aJUUpnFZ@@@LBinGbP`\nfoQP@@@XRfYV}|LEjsfjYA@h@B\nfhe@@@LdbbTbTurXkfaF|DHbhh@B@@`\ndaFD@DAdfUjyjf`CCDsnQP\ndax@@DjWzjh@pxHPj[nQ`\ndmO@@@`diWe\\JUfh@D@LFanFUyB@\nfdep@@@PHeKOLjoxH[s`sU@@@D`@D\ngGP@DiVj`LBl[qB@\nfdy@`@@QRYWYVuzL|F@Bj`@`@B\ndeTH@@RYWZf`@f@CA`SBj^PH\ndmtH@@RYWYih@IhBNh\ndg\\HpMBPRPrPrJPqIPsCMT@ET@D\ndmv@`EBHRYWYih@Jl@H\nfoApB@EZ\\BHRYWYyZLlz@BjiR@B\ndaG@@@rdifzxHH@CAB[bYp\ndev@`NBHRYUVfFX@Ja@LJCBz^H``\ndcnH`IHHaIfU]XUvBBJa@H\nfoAQB@FF\\BHRYVwUZ\\EJ@HjjQ@B\ndg]L`LfDD@cIIBhhd]ikTC@P@XMc\\nwrI@\ndaGH@DJ`rJIPdsTt`D\ngNq@@dsUUPFEDM_I`\ndklHPABHvHrJJIJQEn`HJj@B\ndeTH@@RfUWihHH@CCdpfxYyE@\nfduA@@@ILroZviFBdZpXmUUU@AT@A@\nfhya@@O@rJIJYRIKFcIF]z@BjjjJ@CCPDXhJRmFl{d@`\ndk\x7F@@@bdie^urnGSfj`@X@H\ndg|@@DjU^yZx{BAJjbBJX\ndid@@DjU^nBAHBJ\\FaLJaxc\\\ndeUD`LjD@HrREQICJt@@@XUc\\L|`@\nf`qPb@LB``@QddaTTRRVtzwej@@@@@B\nfdu@@@DjUe_[e`nRM[tYfjjjjj`@`\nfb}@@@LdtTbRLrTrUAJt{r\\}UUUTuUT@D\nflmA@@@ILkZvmsQdeV\\EMUMUUUUT@D\nfnk@`@@MrJJQJJYERKH{`iJ|dmZjjjjjjjj@C@hHXIQgARUYw`qSdmV|ah\ngJPhLQbIKTpD\nf`i@@@DjYU[wxjQg^`@@@@@@CBRLDirVkNy`Qo@\ngGUHEj@aIgZdCAA\\VH\ndcn@@@RieU~V]jjjj`CATpj[ae]yD@\ndmO@@@`diYuZZU@`@@@LJsae^PP\ndieD`JXaCDRYgvzejX@`\ndeVH@ACHhhdbwPDE@AaPXUMsNFP\ndknH@FAIfUWaVhHHh@LMBLJfxYt~Ph\ndg]H@FlDfYU^QVhHHj@B\ndknH`Bp@aJYvU}NfhH@@H\ngNqdEb@b^FQRHmU@XXU|VH\ndg|@P@bEbLbbbTRRKR]pP@P@@D\ndmtL@@QTfyeQehBB@B\neMdHTf`pQyP\ngN|@@ldssUPFEBM_I`\ndef`@@SFyIeYfjZ`B\ndidHPBCDGDRYgVzB@`@`\ndigD@HHPQInYXVijPCCDknPp\ndmvD@EADfzUQUjjj`B\ndifH`JD@cIEDd\\jfZh@`\ngOq@@dlvKUTAaeWqP`\ndaDD@@QInXjZjh@`\ndmMH@LhDeemZZUYh@H@H\ndeVH@ICHhdhUSP@U`A`RYS\\PN@\ndclL@@STfUmVfeVi@B`@pKaW^Ie@\nfheA`@@JLeLrsLsPPkR\\kUUSUUUU@A@\nflui@@DXYYHDfYewvUHrVhiZjifjjZ@CCUNBdkQg^CGdD`\nfoQa@@N@rQQQQJKGbiVLz`BB@D@@pxQ`iZM[bEN@\ndifH@BAIfuxZ`@@CB`inGbQ`\nfoAP@@@NrQQQUQIZgCVifh@J@@`\nfdy``@N@PdrwL}kjtFKTuUA@H@FDbLxHug^CqBh`\ndig@@@`Tke]nX@H@H\ngOq@AdbbLUUTA`Uc^Q`\ndk\\L@@jTie]urnF``Jf@B\nf`i`B@F@bDfUmYWiqQg^`BJB@H@B\ndeVD`NFPBDfUvih@I`@pKBj^HB@\nfde@@@LdbbTRtQRSIJuzHQEDDDD@A@\ndidHPABHZHRYVZz@H`@`\ndidH@@RfU~F``@@pYLJnGdL\ndmMH@NTDee[VFUX@ab@LBpfWbQH\ndcmH@DhDfUe]aWVjA@`@piLzUx`z\nfluQB@NNbBHRYeg[Y]dkNBfBAbB`P`@H\ndid@@Ldbbq[`bB@CBdJfGdX\ndefJ`JaLFP|LddjRcUTpA@\ndifD@DCdfVZaZjZ@LLRaxfL\ndeT@@DjU_k``b`@pzDJfz^I``\ndeWH@DJ`RYYTfZij`CCNF^He@\ndmwH@DJ`RYYYIfjZj@H\nfb}@`@@YRYUUm[ehrRkNBf@BjjjjJ`@p|CANBdkQk^CEN^PP\ndcnH@BAIf^u[evB@@@@H\ndcLB`HSFCpRj}mUujh@@@`\ndidH`HApRkm^Fh@@@`\ngJPhH`xQ}TA@\ndmL@@DjYUVgaHIBh@H\ndidD@@GIEEE\\jfZh@`\ngFp@DiTvVhC@qX|Ph\nfj}P@@@IrRJJIQJIHdqUmA]Fh@bBHJh@@p}bgARUhug^CGIP\ndknD@BADfyvUtvjh@@@pELXYWSxaT\ndeTH@@RYWVf`@j@CC`SBkiy@`\ndg|@@DjU_eZx{BAH@@BJlAaBjUt~Ec]X\nfj}a@@O@rJIJYRIKEIRcIF]z@Bjjjfb`@`\ngNy@LDegjj@plVKxaP\ndo~D@AA|dTRRbqbyZBBbj@B\ndidD@@yJ[VXZ@H@CBdpf{dB\ndmLD@@iJYW_JxZB@f@CBdIagdJ\ndaF@@@RYe[hB@@LJCBinQp\ndieH@FxDigwJiejBJlLpjyb\\H\ndg^L@LAER[eYWNvjhD@@LAPjF]N~Qp\ngNr@Ab|dTQeiDe@\ndaF@`H@HRVU[jjj@LNaLJf{d@\nfnk@`@@ErQSISISISQJUgEZ]zNZ@Bjjjjlj@B\nfb}A@@@YHiiDdYdidjBUiwdyjjjjh@J`@H\ndmLD@@eIfUTfEV``R@B\ngJQ@@dru@XS\\RH\ndeVD`NFPbDfU^ih@I`@phLIiyB@\ndmuD`BFUBHRYg^[hHBX@H\ndid@`@qDf[Wai@@@H\ndcmH@DhDfUe]aWVjA@`@pKNe^X`z\ndmvD@HADfueYUi`@@cF\ndeT@@LdbRQU\\DBT@FCPaTL|SF@\nf`q@@@DjU]YWiqw``R``@ABe@\nflmA@@@YEDeMDhTieQbmN}EKUUUUUUT@D\nfduA@@@ILkZvmmFRUYpXmMUUUUU@A@\nfoAA@@@ISZzrkNV]P@T@@@AaIFBTXkQk\\PIp\nfhi@`@@^Rfuue]Yrsj@B`@h@@`\nfdyA@@@ISZvnsKNJuP@@EUP@A@\ndo\x7FD@DgP]IefU\x7FTvhJBI@B\ndkm@@@kIHbdhdmNjj@@@LIaLIae]xeF\nf`qa`@D@EyIfYUwrLuYj`XHf@@pYrRmFmqEG@\ndk]H@FDLbbbbbey]N@H`@@B\ndeT@@DjU_k``b`@pFDpj[iy@`\ndg^L`LaC@HrRPzIKJRju@AL@FCXS\\L|bp\nfoAa@@I@RYeU{QJlyhHHjb@B\nfoA@B@@RFQQIPiQITYYt@Dp@@BGF\ndid@P@qAqDfY]n`HH@H\ndidHPFBHJHRYf~FBH@@`\nfdei`@LJLzHIZrQQQQIJJkDiYpXuUUTBAP@D\ndeTH@@RYWZf`@f@CB`SBkbCH\ndeTH@@RYWVf`@f@ckA`SBknHL`\ndmMH@HxLbbbTQ[iV@@@@@`\ndk\\@`@BDifUWGUN`@@@@CBPpj[ae^Ii`\ndeVD`LFPBDieWrjfiXHr`\ndg|L@@Z|dTRbfQuCNtEDQR@D\ndg]H@IlDfYU^QUhHHj@B\ndcLD`BtHaIfUuXXHHj@B\ndk_H@DhpRYVU}aWVjA@h@H\nf`a@@@DjU_gVipHDffXAbbItp\ndo}L@HdDuInUunxV`@ij@B\ndcvB@LANGHieEDXuTuPA@\ndkn@@@rQIJJIGSjhJ@@CATrfxYWSx`P\ndaxD@@QIeUjZBB\\Bpf{dT\ngCh`DLdlpd\\InHP\nfoAA@@@ILkjrmFV]@AL@@@ar`\ndk\\D@@WIEEEDYjF]hJHbP@`\ndmtH@@Rfuu[j@BXBAlJSBinXLj\ndklH@@Rfuuenh@Ij@`V\ndidB@@[aRfU\\jjjh@p[BinGbP`\nfhi@H@@XyIQgArQQIQJ[JV`k^jjjjjj@CCTBTxJVcN|GfDLl\nfhiA@@@ILrsZkf|xKUP@@@@@XFps`iJtZsoAyA@\nfdyA@@@YEEHeEBdekA}Ejj@@@@@B\ndmVH@NAIe^Yjjj@LNpjXYWdH\nfhqa@@C@RYegV^tyjhHJJ@@`\nfoAQ@@DJ@drrkNrMYsTDEED@FDdDYpVcV]sEIJ@\nfjc@`@@YRYVum[ezLdjs`i`@jjjjhj@B\ndeTH@@RfUWihHH@CCdpjxYyC@\nfbmA@@@YEEMHdcLeHVRt[pVjjjjjjh@H\nfbmA@@@YEEeHhmEhXmV]zNVjjjjjjh@H\ndk^L@IANRY[f~]tvjjjj@LApjXUt~It`\ndo|H@@rJJIPjHqaZZijZBAh\nfgAP@@@LrQQJJYPtdKQk@@@@@@@XBQPTeZMX\nfH`pA@\neMJH|Df`pgd@\neOB@HcfhH',Mq='gFq`@ldrfmU@XR|a@\ngKXHL@aJWFj`LBHmrD\ngFyHL`DILikUPD\ngKTHLPDIRxtlA``nR@\nfHgdA@\ngNp`@dfzZj@pJMX\neFA@H`bLFE\\`\ngCh`LDdsPFDWI`\neMHAIdLF^P\ngCa@@dkHFBVyH\ngJXLL@aAABDfVhB\ngJXHD@aIUZ@`\ngJXHD@aIYj@ppqyH\neMXIDe``\ngOx@@drm\\@@A`plZp\ngOx@@drm]UTAaqEcV\nfHapA@\ndid@@DjUfaBB`@LNaLinGdD\ndidD@@IIf]xViZPCBlJfGd\\\ndidH`DBHR[e^FX@@@`\ndid@p@bBbFbDfYoa`b@@H\ndidH@@RUe^Fh@@@pxHPj[nPH\ngCd@ADkj@ppbyL\ndifD@BADfyWaZ@@@LBQnGdT\ndeV@@@RgYTYj`@@C@PaLinGdR\ndiDD@@QIn]Zjh@pYLInGdT\ngNp`@dfuZi@pJM_I`\ngC`hIaxIVtA`enP@\ngGYHKAxIVkU@XYZ|PH\ndiFL`JcaCpRm[fjf@H\nfHbTA@\ngNtDLpDHHR[UjhB\ngGYHEAxIVsU@XEXwd@\ngNq`@jdvkSPf\\Ll~P`\ndaxB`HSN@HrRPyKUPAaaTwHx\ngCaHDIAIZ`LINS@\ngJXDBIARBS]TA`h^S@\ndiE@PBDIAIAInujjh@`\ndiED@BDDR[mVjj@LFPj[ayF@\nday@`Lx@aIUVjj@LNaLJf{d@\ndmv@@@Rf~UeZj@@@LEBDpfxYT\ndmvD@E@dfYwVzB@j@C@PpjxYT\ndmvL`BaL@HrRRjIJUVjjh@`\ndmvL@EAFR[f_FV``H@LFRfxU@\neFA@HoBJD\ndiGD@JxPQIeUZfX@`\ndiFB@BAFEInuZjd@pILJnQp\ngCd@ADij@ppfyD\ndaFD@DCdfVRiji`CCDinQ`\ndaGH`LJn@HRf_rjYj@H\ngC`@Die@ptVy@\ngNq`AVeJmUPFEbq_IP\ngC`DAb[DRVhB\ngCaHLLQIZ`LDEqS@\ngJP`@dfzhCCA[ba@\ngNp`@dfuZZDu`TZ~S@\ndaxD@@QInuij@LBRf{dD\ngGPBADZPLaYAIZjhB\neMABHYAIhH\ngNqBLIAREdGHIMmUTAaCrP\ngNt`BpdfyZZDu@\ngJX@@dmu@XKGdP\ngGX`LDdsmTA`m^P`\ngNx@@eZoUHD\ngJY@HDeVhCAK\\a@\ngJXDB@bABUmTA`Pl^R@\nfHbxA@\ngGU@DPdrmRAaecrT\ndedd@DHYCdfW[ZfZBBX\ndcLL@@STfVyVUZ`HD@H\ngNp@DkVzZ@`\ndaD@`@bDfUjZ@B@CB`SJ{dL\ndcLL@@QTfVUV]Z```@LNJfFUwd\\\ndaxh@LInAIUYffBDh\ndaDh@DInAIf]nZiX@`\ndmu@`LTHaIUe\\Yj`@`@`\ndaFD@N@dfYvyjV`C@lJfyG@\ndeVDPL[`bB|DeYgFZjZh@`\ndk^@@@Ri_YVftzjX@H@LEaLJae]L\ngC`HADIMTAa`mrP\ndeT@`@bLbdTJPsU@@@D\ndcND`La@BLddJbRrzmPP@@X]S\\LkoHp\ngNphJpDIKMULA`QD^Q`\ndiDL`JXPBLbTTJjfd@`\ndaG@`LK`BDimVz`@@B\ndeU@PAdH`haIf]vzB@h@H\ndmuD@FdBRgmeeZjjV@LNaLJaWd\\\ngNp`ATeejj@pXfM_H`\ngNqhHjOAJmVj`LLZ~HX\nfH`TA@\ndaxL`HS@BLddNRuT@XXUMrN@\ngGR@ACTkVY`iZ`\ngN|@ABeJ|mPTdxJofVH\ndax`@@`fRegej`RIJ\ndayD@LXJRVUYk`aJ\ndigH`LJn@HRf_ljfYh@ppBGdH\ndaFH@BAIeZf`@`@piLiny@@\ngNuHDx@cIHTej`H\ndknH@EAJYYfQN`bJ`@LIRfxYWS@\ndifH`BDHaIe]ih@I@B\ndmOH`LJQ@HRf^yriVfZZh@`\ndmNH@HCHhheDVzU`@@@@LNCJXYWd\\\ndcnH@AAIVYWXUvjP`H@LABDkag^Ph\ndew@@@pldTTJVTLmP@P@XUaMt|`P\ndaEH@DhDfVVyje`C@dJfyG@\ngNp`@deVZf@pvkyH\ngGQ@@dluTAaQEkrD\ngNq@@djuUPFCDqkyD\ngNq`@fdskUPFFq_IP\ngJU@HPdkU`P\neMhHRVCBP\ngCi@HDfZ@`\ngCh@@dmHFFBwI@\ngGYHLQDIJuU@P\nfH`XA@\nfHbXA@\ngNpHAxITkUTAaqEcV\ngKP@Di\\YZ@phbq@\ngJQHH@aJmj@ppQyL\neMJDBDf`pQyP\ngJT@@TfVdCAK\\PH\ndiG@`Dp`BDfWejjPB\nfHcdAa@\ngCh`LDdkPFDwI@\ndcNL@M@iRYg^vzB@j`@pdLJnFUt\ndigD@Dq`EIeYkfje`CCDqnQp\ndiT`@@rnRfUjEnBA``B\ndmNH@DAIe[VfeVj@B@CAdkae^Q`\ndedHXJBPRPrPzPFPfPrJPsRUSU@D\nfHcDA@\nfHdXA@\nfHdxA@\nfHgHA@\ndeTH@@RYeYn`HJ@CC`JfxYyF@\ngGQHJ@aIUjhCBqXwd@\ndaG@@@kdiVrX@a@C@hPfyG@\ndcll@DpYAmRYV]zzUZije`B\ngNtPHPjtfvZf@`\ndeVD@HADfvUFVh@@@`\ngGT`EaTf]jPLDmrD\ndaGH@DK`R[e[fiZ@LLQnyE@\ndmvD@BADfyW^Eh@J@CBdingdJ\ngGY@LDeVj`LJEc^R@\ndcnL`LaA@HrRPjIKTrzmPHD@FEYtkh\ndaE@@@aJyUnX@@@pkBinyD@\ndie@`NDHaIe]ih@J@CA`SBknPH\ndeU@@@aJueQfj@@@LAaLinF^P@\ngN}@EbDf]Zj@pJqoHp\ndeU@@@aJueQfj@@@LAALinF^Q@\ngJT@@Te^hCCKD\ndcwD@BWPQInuVZjjX@`\ndiVH@BAIfUInFjZi@H\ndaz@pL@HPHxHRYyZj`B\ndaz@pL@HPHXHRYvZj`B\ngGYHE@DYIHUj@prqyJ\ndayH`Dr@|Djyfjh@`\ngNp`@dfuZi@pJq_IP\nfHdpAa@\ngNv@ALZtkefj@`\ndieH`Dq`BDfUfnZii@H\ndaF@PDBHzHRYWih@H@H\ndaF@`N@HRf_rjYj@LFALJnyC@\ndaFH@BAIf]n``@@phLInyE@\ngGQDB@baBSKTpD\ngJPhLPDIKTpFDGEB\ndeV@@@RiU\\Yjjj`C@XSBkagdL\ndg\\`@@SFRYueUNvjjjj@H\ndg]@@@IIUYyT{zjh@@@`\ndaxH@@RV^jj`CAhPj[nHF@\ngGX`JDdjmTA`c^JX\ndg^D@IADfvYW}MjhHB@B\ngNq@@dl{MPblLZ~JX\ndctJ`HSBtOCIILXhduSUHQS@\ndkl@@LdbRdSRjP`jJ`@pfDJfxUOdZ\ndidH@@RYeVz@``@pXLJf{dB\ndigH`LJa@HRf^|jfZX@`\ngJQhHb@cIHUhCB[bA@\ndaD@P@qBbDfYvzB@@B\ndkm@@@IIUfU}OjhH@@H\ndid`@@pjRfUjXBB`@pSaxbL\ndcND@E@dfYwYn``Jh@LACBkaW^QH\ndaF@@@RevpjjZ@LJAB[nIF@\ndmvD@BADfy]fyjjf`B\ndifD@NADfvWaZjj@LJSBx^IS@\ndcND@IADfvYXYZjjZ@H\ndcND@BADfuYU]Zj@@@LASJxYW^PP\ndeTH@@RYe\\YjB@@C@Ppj[agdP\ndmvH@FCHhhhdYVhJ@@C@PcJ[aWbEH\ndcNL@LANRYygiUjB@`@`\ndifH@HCIDhYLJejh@`\ndaED@DHFRUe[fff@LDRnIG@\nfH`PAa@\nfHdHA@\nfHchA@\neFHBJFE@\nffs@`@@HRYyVWm[W`cVCDmVjjjjh@J`@H\ndcLD@@UIUe]FVX@J@aKCdrfx]yD@\nfewI@@DXldBSLjj\x7FJzriFJMEHymUDDUUT@@P@D\ndeVD@HADfyeFV`H@@piJ[ayD`\ndet@@DjYUZ^D`dJ@CBlJngfHp`\ndaFH@LAIVUnjjh@pZDJf{d@\ndg~H@AAJYU^wiNzBBJh`B\ndet@@DjYUZ^D`dJ@CCdpjxYyC@\nf`qa@@D@RYeU{TRg^Zij`@`DJVBbTxIVcV]{bLL@\ndmLH@@Rfum[iV`@jD@`\ndclD@@WIEDhbmSCMAPEH@XTa]N|Rf@\ndclD@@UJYY\x7Fkaf`bBh@LJPfg^IW@\nfgAA@@@ISLrotyHvk@@@@@@@P\ndeUD@DhFrJIJHusUMT@XYgCqDd\ndaE@@@iJUdf@H`@pRnydHp\nflu@P@@\\u{HhheEEMGEfJdzJV`hJ``H@@pecdeZMzJ^PB\ndet@@DjYUX^dHbH`CAdJfx^Id`\nfbmAP@@BU{NQQQQQQFVKYbmNbehJbHHB`@B\ndeTB@@pYRf[^njjj`CBXSBinFP\nfc\x7FAR@AU[OXICHiEHhhhblXdicBkvkNbdmUUULuUMSUT@FD\\DYpTeZMYp\\eFo_HM@\ne`\\PM@@@AHdbjnikgolbbbbbbTtRQTNTtrSBfgefRvsPADADqUUUUA@@P\ngNpP@btf{Zi@prqyJ\nffsA@@@YEHeHlcDeidfRtXLRujnjjj`@j@@`\nfb}`@@@YEDeMDhTihjLUiwhyZjjjjjjh@H\nfjcP@@@FrQQQQII[FIJgIVt{IP@@@``@@@@`\ndknL`LaE@HrRPzIJZ]Vh@b@CAlInd~Qh\ndeUH@JdDin]xZB@`@pIBX^IT`\ndg]L`LnDT@cIIChdieNkT@QP@P\ndclH@@rQQRJJuJ{PUDD@FEDwCoLjU@\nfbmq@@DX|CHhhdhddlTdQfgV|uTTEEUA@Aa[ENBdkQkNCEOHY@\nffs@@@LdbbdRdRQbLcKRU{HTmUUUUUUUP@XNqQgARtZwhiKWb@y@\nfewAa@AMfBPSHYEHiEEEDUeheme^uYt\\ujjjifjjZZh@L@xHs`iJtZs`yKT_H]@\neghPJD@AFimk`IAJPrJQRJJJHkKQK[PiJ^ZVNQYUUUUTsUUMMSU@@P\ndedL@@PTfUUZjf@LJpjz^Pp\nfnkA`@@N[dTRTtTTlVRbUFJlFNZmKUUUUUUT@F@TXipTeZMYw`iKWbLm@\nfig@P@@NZOHhdihhiXleDbjLUXL\\uZVjjjjjij`@`\nechX@@@@LdbbbRRQbvfRdrqC@bfadRUSRvuUUUUUUUUUTp@D\ndg|L@@v|dTTbbluRitEQDR@FCTwSm|Pc@\ndeTL@@zTfY[XXBH`@pHJz^I``\nf`q@P@@Ty{HeEEDhYDIVcLDQEUB@D\nfi{@p@@Zldh|dTTbblvrbbfaiIQjmhJbHjBHJP@B\ndk]D@BXCrJJJIPkiWV`d`X@H\nfhiP`@CAAQJ[gnyodgAhHHBB@@B\ndidHPBBlFlRYgVzB@`@phLKayE@\ndmvH`NTHaIeU^f`@jP@phLKaxbR\ndeV@`ABHRYeun`HJ@C@`pixeB\ndmv@`EBHRYe][hBBl@H\nf`q@A@@QGhaIfUvuv\\lz@`h@H@@pdCAIJtZwnHhp\ndeV@`IBHRYWVf`@k@CB`SBkd\\\nfbu`c@MDTBleVYbLbbbtReTTRsAQSAADtAA@@D\nfhi`A@I@b@qFQQQIIJEKLx{tAAT@D@@P\nfns`P@H`EYvQQQQQVZEIQILDhuVm@tDMUUTP@XZiqQoASdmV|`p\nflmA@@@ILkZjmsSdeV\\ELuMUUUUT@D\ndmtD@@QImYVUZh@@@p{B[ae^P`\nf`i`B@A@bDfUmvUhpTg^`@i`@H@B\ngOx@@eJimUTA`plZ~R@\ndeVH@LAIVYqejfj@H\nfb}A@@@ISZvmk\\lxkSoQsP@UUUUEP@P\nfig@p@@HdFJ\\bbTTRRtVbarRipwlezMjXJBQj`B@@B\nfoA@P@@HeYIefU]ipsfi`hD`@H\ndcLL@@G\\dTRRbOKPPTP@XRfES\\L{rE@\nfdu@@@LdbbTRrLtSI\\TZwhaDTUUTU@A`[AJLxIVcNCGfDR|\nflmA@@@YEEDhhTmMDsdmV\\ELAAUUUQT@D\nfdu@@@LdbbTRrLtSI\\TZwhaDTUUTU@A`YFBTxIVcNCGfFR|\nfnkAb@AUZBPrJQRJJJHkFIIWmV]EIZjjjYjjZj@CCvBLxJRmFlxNRmyAH\ndid@`@qDeYWaf@@BHlJ@j[nYBB\ndid@`@qDeYWaf@@BH\\F@j[axdH\ndeTH`ICDRUe_af@B@bJ\ndmu@PNTHbPaIe]Zf`@i`@pxDpjz^Ph\ndmN@@@RfV_Jaf`RB`@phFzUxbT\ngOy@HDfUkjj`LJlZ~P@\nfhy@@@LdbbTRtjVYHRg^BHh`@b@@`\ngNy`LETeVZZDs@\ndazD@BADfUvjXHj`\nfb}A@@@ILk[Kk\\tyKSoQsUSUUUUUP@X^AAFJLEIVcV\\EODQr\ndmMH`LVPBDiiu\\XUjjij@H\nfc\x7F@R@AUFlDadTbdTTTQVLRVtJU{UgQRVjjjfZjffjj@CBNBLxJRmFlxNR}godD`\nfnk@`@@UrJIJYQKISQHrgIZ]zNZjZjjjjjj@B\ndmtD@@QImYVUZjjh@pGBinFUyB@\ndk\\L@@STfUmWiiUjP@j@C@nESxfT\nfhy@@@DjUg^UzBu[pP`hhH@`@H\nfb}a@@L`rJJQIFYJYKDyIUgQSUUUUUTuP@P\ngNqhHf@cIICejPLH^KT\nf`aab@APQp@QddafaRbR]yZij@@@B\ndid@`@dDfYUn`HD@LJpj[nPH\nfoA@B@@RFQQQIEQILyYtA@p@@BGJ\ndmvH`Dd@aIVUwaZ@B`@pxDJfy^Q@\nfj}P`@GEAKIEDheED]idaZm{IZB@hBJjB@CCrBTxJRmFmxLTyy@h\nfigp`@DZBAJSJvnrnwJtYIW`}FuP@UUS@AD@A@\nfoQ`@@@ISLk^~FBuYuAH@@P@AaHANBdkQk\\Pqp\nfc\x7F@Q@AMzlD`fPrJQRJJJHkKQK[ZJ}jshykUUUSMUTuMU@A`GAF\\EIVcV\\DmQ{yCH\nfc\x7F@Q@AMglD`fPrJQRJJJHkKQK[Pj}jshykUUUSMUTtuU@A`GAF\\EIVcV\\DmQkxcNp\ngF|LHjOC^A|DiTt@@B\nf`qPa@LV``@QzPrRPqQIIXjlYruPTDT@A@\nfluQ@@A^AdRbRrTTRbKNN}EK@AP@@@@@P\nfa{A`@@E[dTTbRLrTfVrTIrRkNSejjjjjjijh@H\nfhiP@@@\\RiU{WZisQjjh@Jj@@`\nfj}PB@AVADILkkJzrmFFmdl@EUUP@D@A@\nfbmP@@@\\RiUyU^uipsdyjj```@H@@p]AAJBdkQg^CEODXr\nfdyAp@@XyHwhiSJrwoMANCUUTDDT@AajAJ\\EKQoAcqBE@\nfbu@H@@TyjwdyrIQQQHqISYBTFKADQUUTP@P\ndg~HPItHciCHhheEJlpz\\ADMLPAaPXpwrE@\nf`q`a@IZ|BHyHYEEDhiULfBtZ@bFfd`@`\ndmLH@@Rfum[iV`@jH@piLIaxdj\nffsA@@@ILklsJwoQbk^RVl@EPT@@A@@D\nfnkA@@@ISZvjsJklyIRZe[P@UT@@@D@@P\nfdeA@@@ISJ~sJ\x7FAF\\FMA@UA@D@A@\ndetH@@rQQJHtpsPT@`AaQMpsqDh\nfdeIB@LZmzH@cIEEEEeEDnRfkF}TmULsMP@XLAQbcNxbPP\nflm@@@LddTTVbvebyjUcNBfjjjjjjj@C@LDXIQgARUhuoQSr@P\ngJU@H`dkSBI`\ndmtD@@QIee^UZ``@@phfFUxe\\\ndmtD@@QIee^UZ``@@pinFUxaT\ndmtL@@QTfyeQehBA@C@jXYxb\\\ndcLH@@RYeZvz@`j`@pxBXYW^Q`\neMhDRUB\ndmuD@HTDrJZIJGaZ@B`@piBkix`j\ndmv@`LBHRUVUeZj@@@H\ndmv@`LCDRUVUeZj@@@H\ndaDH@@RYU[fei@LJpj[nP`\ndiFB@LAFEIg[Zjh@pkBiny@`\ndknD@MADfym]eVj@B`@pELinF]yA`\ndmvD@D@dfWeYUj`@@CCdpfxYyB`\ndifD@LADfyWaZjj@LBinGbHP\ndifH@BAIVUxZjj`C@PaLinGdD\neMPBcXLIyP\ndcNH@DCHheEDbnmPT@@F@hUMproHt\ndcMH@ITDee]UnX@Jh@H\ndaG@`LK`BLdTRIUSUpPea`IWHH\ndaDh@DqnAIeZfZZd@`\ndaF@`BBHRYg[hH@@LJCBknPp\ndaEL@LhDYIe\\jZjd@pcnQ`\ndk^@@@RfYU\\]Tz@@@@@LECBinFUOdZ\ndidHPBBHzHRYgfFBB@@phL[ayA@\ngOt@AdiguZj`LBmWrD\ndk]H@HdDeVYWz]NjdHB@B\ndcnH`Jd@aJYW]rnF`PJh@H\ndcnD`HI`BDfYoVnWZfX@@@`\ndk^@PIBHDHRYWY^fWX@JZa@H\ndg~H@ACHhdieDYSSl@EMTPA@\ndg~D@DClbbbTlRUL]mUL@D@FGYSRngrK@\ndmM@PAdI@iAIf]UneXHBd`B\ndieDPLZD@HhHrREQKaVii@LBrnGd@\ngNq@@dsmUPFEbq_DJ\ndctB`HSE@HrRPzKJKUUT@P\ngCi@DDfj@pPwED\ndmvH`IDHaIfUin`HJ`@`\ndo~L@M@iRYg^ufzB@jj`@pLLJnFUt{\\\ndo~L@CB]RfyV~^F`HJj`@`\ndg\x7F@`LIPBLdTTJRrVUJtsUUM@FD`S\\Lz]`\ndaxD@@QIgUjfBJlBpj{dL\ndmt@@DjYefdHbR@CBlinWbCH\ndklH@@rQSQJYGahBBi@B\ndkNH@DAIeV}VijY`pdYZ\ngJY@HDeYdRVFBwHP\nf`i@P@@HD[HhdihdhUSbkN|uHPTUB@FALT]{bEZ@\ndaxL@@rdenjjh@pZDrf{d@\nf`qpB@MF|BHRYWUVVcG^`@jjeH@H\ndclD@@[HhheBeSKkTp@P@XEpkoLri@\ndclD@@kHheCDdUKkSP@P@P\ndaFH@FAIe^fZVh@ppj{bI`\nfdyaC@OVAHHdLRFQRFIJZZKHRoAj`BIZX`@`\nfbma@@E@rJJJJISIIH|TyL\\mSHPA@D@@P\ngFp@DiTujdCBqXwd@\ndeVH@LAI[eQejZi@H\ndk]H`FTpbDfUmYkiV@Hjd`B\nfj}qB@BVBxDadTTRbRRTVasIFCzLDEETuUD`A`kENBtZJ\\eyCp\ndid@@LdbbQxXF@@CBdrf{bDH\ndklD@@QIge]aej@BX@H\ndg\\D`EtICHhdhXhUSP@UMP@P\ndmtH@@rJJJJia`HbP@phBinxdj\ndeeL`LZDh@cIHULhmMT@XUc\\L|`@\nfbuqb@LRuA@`AFRREQJKKSISUhiVhFB@@`@B\ndmtHPBBHfPRYeUXXHHh@H\ndk^@@@rQQHjIKJtzjX@HBJX\ndg^L@FAMrJZJIKPbkPPPr@FGIS\\JgrK@\nff}`P@F@QkNQSQQIZEJIJXRdeZBBFX@b@DCT\nff}`P@F@QhnQSQQIZEJIKHRdyZBBFX@J@DCT\ngK\\@ABeKcMHFDgDl\ndmO@@@`diWe\\JUfh@D@H\ndmO@@@`ldRabRpiVZ`@P@`\ndg~DPFvpbEBLbbRfbRM\\JpAESI@D\ngC`@Dio@pfyD\neMBBHR[B\nfjc@`@@ERYYefUez\\UyvPiZZjjjjjj@CCrLTXJVcV]zJ\\ey@H\nfb}@`@@YRYVum[ehrRkNBf@BjjjjJ`@p\\AFJ\\EIQk^bgODZB\ngGU@EPdvmVAaErX\nfb}A@@@YEDeMDhTihjLUiwhyZVjjjjjh@LGHqS`iJtZs`isrDP\ndk}@@@aJVyY|hUtyjBHJ@B\ngGY@DDfvZPLLC^IX\ngNx`HLdmluHdX\ndkl@@LdbbTbsyYfjihDaV`\ngNp@DjWZeDCaXl[rL\ndet@@DjYWkafefiPCBdKagfYS@\nfhiAb@HHp@HrRRqQIYPsdcVVfiVji@B\nfoA@@@LdbRbbJvDkQfijZjT@H\ndknDpNWPdDdLdLbdLTRbdzhHF@@`\nfhq`c@IF\\BPQHXdLbdLRTvfdDjj@Hfj@@`\ndif@@@Rfy^F`H@@pYBinGd@\ndmvL@HAER[e]xV`@h@LFPj[eyD`\ndkmL@HLDUInUv~Eh@Jh@H\ndknH@LAIYyfUMj`H`@H\ndiFL@JAARY[vji@LLpfxe\\\ndiFL@BAAR[fvji@LLSJxc\\\ndeVD@JADfvYzVjZd@piLXYyG@\nflea`@F@bgHhhhdhebmIVTYhIBhjXP@`\ndeU@`Dp@aIUeQej@@@H\ndg|@@DjU_eZx{BAH@@BJ\\MaBine]N~Q`\nfhy@`@@BrJJJJIQ[DQeNmyh@`HB`@B\nfdeIb@LLDhwh`BLdTrbRJbTOARTxM@PMUUE@AaHBJRu[pX|`x\ndg}@@@yJeWe^nNzjhHB@CAXPjFUt{yF@\ndk]@@@yJeWmZ{Sjj`@`@pVDJnFUt~Q`\ndcmH@DxLbbTTRM\\nmPLA@AaSCJ{rD@\nfhyPB@ARADILkkKZtYhw`p@UR@A@@P\nfhyPB@ARADYEEDdXeLdsfc^C@PUH@D@A@\ndk]H`AdpqDfUmUiev@Bfd`B\nf`y@@@LdbbTRfNSI\\TKSoADTT@A@@XZPRcARtZso\\a`\ndg|H@@RfuvU[cn`@`@@@peBiiWSodX\ndg|H`ABHRYW[ficn@BjjH@`\ndk_H@FdprJJISPkatzjjiZ@H\ndid@p@qBqFqDfYoa`b@@H\ngGP@DiUja@xEXwbD@\ngGP@DiVV`iJpJqoDH\ndeTH@@RYWVf`@j@CC`SBhYyG@\ndmtH@@RYWUih@Jh@LNALJiWb\\H\ndcLD@@SHihheJLmUTu@D\ndg\x7F@`IVpbDfUuYZX{`@jiR@LAALJat~Qp\ndigDPLXXP@b`cIHUDnEZfd@`\ndk]@`AdHaIe[VZYS`@ijH@pxDret~Qp\nfoAaB@G\\ADILkkJ}FFm@AP@P@A@\ngOxhMDOAJmZvjhCCXobE@\ndk^D`La@|Dju^UimMjj@B@B\nfdyQR@NBYhwh`qLbRbbrRfTHiNCMKUMUU@aB`\nfluQ@@DJAdTTTTTRTQsEVbEKU@C@PD@@P\nfoA`P@D@DZrSLrkrQfgMTC@q@@P\nf`q@`@@HR[fYU_Sk^Zh@@@@@LMHYpTeZMYwnPp\ndg]LPElYLHc`cHhdliBiSP@Rtp@XLBXUKrK@\ndcNL@JAMRYYYyUjijX@`\ndcnH`LD@aJY{grevj``H@LABDxYW^Ie@\ndg~D`LePBLdTRTrbNT]uU@pD@D\nfoAab@LRMX@PeL|lmxiiuUAPLP@D\ndeVHPLX@bPaJU{\\Jjfj@LJANF^IA@\ndaFH`Dx@aIYUnZjh@`\ndeTH`N@HRVUmaX@b@B\ndknH@CAIe]eZZ@Bjl@H\ndk^@@@RfYU\\]Tz@@@@@LI@j[aWSxfR\nfhi`C@I@dDRFICHiCDeMhdaJ|Fj@H`@@@H\ndif@@@RYevz@``@pkBkay@@\ndew@@@pldTTJVTLmP@P@XEc]ODHP\ndaGH@LK`RUe[fjV@H\ndaGH@LK`RUe[jjV@LLBfxdH\ndeWH@BZPRYVtYjje`CB`SNGbQP\ndmtH`EBHRYWUih@Jh@LNALJnWbCH\ndo|H`EBHRYWUnjZ@Bjj`@`\ndcMH@AtDieV_ehHJH@H\ndk}@@@qJYVu|kitvi`PJ@B\nfgAP@@@XRf^~UqQgCVZjh@J@@`\ndk^H@LAJYeUquSh@@@@@pzfxUt~HT`\nfjcY`@LJCENPBedbbTRTQRRTqJZLxLTuUU@aDuQ@AaJAJ\\EJvR|``\nfmoYr@CA[t\\xHpSbk@|LddqTRbrfbRRRfejJRmzMjjjjhDHfjH@H\ngGPhMQHIJmT`XR\\RH\nfewh`@DRMdXHyEDhhhddmiUDenRuzp^cZ`B@HFjjhH@H\nf`qPC@MNADBBEQFQQQIQGQIgG^`HJjeH@H\nfdy@@`@Q@`bTQFHRYeWvUv\\|F@`j`@`@B\nf`qPA@MNADBBDfUuUehqwh@JjiR@B\nfdu@@@LdbbbRRltSIBTZwhmUUUUUU@A@\nfigA`@@U{dTRTtTTjTtRJeFJMxJRmKUUUUUSUP@P\nffcP@@@\\RiUyU_YVgCNKVjjBBF@B@DGVDXHIPTeZMYwhirVoHX@\ndg}@@@yJeWmZntzjh@J@CAXSBhYt{yG@\nfdea@@D@RYZymUr\\lGtVjh@@@`@B\nfoAQ@@KN@eLrj\x7Fbthu@@@D`@D\ndg}@@@cIDmEDhcJg\\l@@@@QE@\ndcn`@@rnuJYVuhUtHJIb@LBSIWbMh\nf`ipB@DX|@HRYeYU\\cMNmyjhhJHH@H\ndmvH`DH@aIVyVUZX@@@pxD[ae^Q@\ndeTD@@QInUnEh@H@LBRnGbEH\ndeT@@DjU^k``b`@pyLinF^IA@\nfoQH@@@XHYJYg]UirRkN@bJjhh@H\nfgApB@LLx@HRevUUpPTcViYj@H@@`\ngOy@DDfUkZj`LBcWrX\nfoQa@@H@RVUmUV\\djsj@Bh@@@B\ndclH`ABHRYVyzy]`BFh`B\ndcLDPND@c@aJY~vrjefj`B\ndidHPACDZHrJIJFn`BH@H\ngO|@AfeJykSTA@\ngNt`MPdfUZj@pvKyJ\ngNt`IPdfuZj@pJu_HP\ndg_@`DGPbDfUueZZ@Bjj@B\ndmWD@Hi`QImfUjjj@H\nfhiP`@C^@aImYWVVRoAZX@JiRACEaiFBLxJVkNB\ndmUL`LZDh@cIHULdeijh@p[FxYWd@\ngOp@DiYKZj`LFHlWrT\ndcLh@DxYCHdeEDcJmPDD@FEIcCoH`\ndcLL@@QTfVUV]Z```@LFJfFUxgZ\ndg^H@LAJUyfUSjhBH`@`\ndeTH`L@HReyTYi`@@CC`RfxYyD@\ngGYHL`DIMlu@XR|pq@\ngGYHL`DIMlu@XHwbF@\ndiFH`Bp@cIECDjZj@H\ndcLH`BBHRYg^fzB@j`@`\ndeTDPBdHcoAIfUvFBBD@H\ndg\x7F@PBWPbAbLbbbRfaSR]pPTMQ@F@acCNg_HH\nfoApA@BZLBHeDYEEDhiiTqUgPPUCTP@P\ndcN@PBBHFHrJJIJXmLDED`A`HXwCJ{rB@\ndig@@@pdigwJZjf@LFSJ[ayA@\ngOu@E`drm[SRAalWrT\ndkmL`LNDD@cIIBhdmeuZ`PH@LFqneSyE@\ndg_L`LfxPPBLddJbbQvfmPLA@A@\ndaxB`HRnCpRk]Zj`B\ndeV@PIBHzHRYeea`Ha@CA`rfx^QH\ndmu@PITHchcHhheEVF@bF@`Z\ngChhMDOBNtA`enP@\ndid@p@qBqAqDfYun``H@H\ndidD@@qJYWkjjj`CAhPj[nPH\ndeTL@@rTie_kjjjh@pZDJnGbMH\ndcLD@@QIeVuWVj`@@CClInFUyF`\nfoA@@@DjYUgWfiF`XJbh`@`\ndidH@@RYWZZ@B`@pXDpjGd\\\ndaD@P@qNqDfUZZ@B@CB`pf{dH\ndklL@@Ptfym]eVj@BP@pyLinF^Qh\ngOpH@DISOkSM@XDZobQ@\ndk\\H`F@HrQQJVIHjtzeX@H@LA@jFUt~P`\ndcLL`NWPBDif~{JifZZ@H\ndid`@@pjrQQQUn@hH@LDJnHp`\ndaG@@@pdifZxBH@CCJ[nPP\ndg~H@IAJ[VvuneZ@Bfh`CCdpfgSyE`\nf`qQB@FF\\BPRYVwUZgAR`BJjfD@H\ngNu@Dpdr{UPFAVMyF\ndg~D@EBldTtVRTOBntDASQ@D\ndig@@@pldTTQk`Ha@C@biny@`\nfoAh@@@XHkVRJJHyQJTxZpDPsUBBGJ\ndaF@@@Re]Jjjj@LFBDpnxeL\nfoA`@@@YHhhhdeCHvgT@D@P@@XJ@Q`iJtZsnPh\ndeT@@DjUijP`j@@pzDinF^PP\nfhy``@E@xeLsO}kSdm^CPUAEUQ@A`YAJLEKQkN|D\ndiT@@DjYVnXPbD`@pYLinGdD\ndmLH@@RfyYxYV`HJD@`\nfoA`B@N@bLbbbTNbRXJshBH`@@@H\ndcLH`L@HrQJSPiKmMP@@A@\ndaFL@LANRYUJfjf@LH[bXP\ndie@`NDHaIe]ih@I@B\ndmO@@@pldTTLRqiUj@``@pKNe^HP`\ndmM@@@iJYYtjFYBHJ@B\ndaF@`J@HRfWrjYj@LLAFxa\\\ndmOH@FePRYVukaf@HZH@`\nfdyAb@HHpCpRkVYU_]Nmyj`@@B`@B\nfj}A`@@HTdssKMzsRuYp\\mPDQTuTT@D\nfhiA`@@Hddjrm|jIW`mPAD@@@A@\nfhiQ@@JJ@djlklkIW`mUUT@D@AaiFBLxJUg^CqDC@\nfakQB@CQDBPRYW[Ue]eZLCEoX@Ijjj`@H@B\nfde@`@@JrJJJJJIIGLUXOhmR@T@A@@D\nfdya`@O@x[HhheDlTdhsdmM@pQUU@@P\ndmLD@@qJY{WJeZj@B@CAj[ae^IB@\ndg^L@HACR[e]VxV`@jX@H\nfdea`@D@dkHheEDheDUSfoVBuLDpTT`@P\ndmtL@@RTeVUaUj@H@C@jXUydYH\nf`qa@@O@rQQQIISGBtju@@@QX@A@\ndcNH@FAIfYWgZ`hH@CAPcBinF]yB`\nfnsqa@BLSFkPRBiCHhhhhhddTmeNZlGrXHbhjejXh@LABehsoE[rC@\nffsA@@@ILvuoZunsfc^BVkSUUUUUUL@D\ndcnD`I[PdDfUm^neX@biH@`\ndk^H`MdICHhheCEfGS`HaiH@prfOfXn`\nfbmqB@BNcxD`dsLrzlvqUhasAATUSUD`A@\ndid@@DjU^nBAHBJ\\FaBiaxf\\\ndmtDPNDHaYAIfVUi``X`@`\ndknH@JAIUYfUNjjjj@H\nfhipR@LBBhug@BLdTRbrRfeIrcJjjfYfP@pxBBBTE[bOA@\nfbmha@LBm{Lm@BDaBTrs\x7F]riFTyEMUULsRuL@FD`PPR``isqGK@\ndaE@@@aJyUnX@@@pKB[nIE@\ngGP@Djuj`LLm^JD\ndmtH@@rQQEJZDjjfh@`\nfhy@`@@\\RYeeo]xjQg^ZB`@@h@@`\nf`iA@@@ISZvmjsbmN}P@@ETD@D\ndclH@@RfumVy]h@Jh`CCdpfFUyF`\nflu@`@@IRfum]fUgEZ}F`@jjjhh@H\nfbm@`@@FrJJQIFIKSIDx{r\\uUUUT@E@@P\ndmtDPBTLSlSHhhmDVFBBI@CB`pnGdR\ndmvH@DAIe[VUZh@@@pyJ[ae^Q`\ndk^@PABHtHRYWYZfWX@JVb@LFALJmxbf\ndifD`JxPBDig{JifZ@H\ndmNH`NTHaIfVUazXBHX`B\neMJDBDe`pQyP\ndmOH@NFPRYegXYUhD``@`\ndaDD@@iJY\x7FJjjh@pjLInyB@\nf`q@a@ARADJbDfYV]}YsQhBBjjb@B\ndg\x7FD@NfpqJYoYWJNzjjjV`B\ndg}B@HTDf|bfbbTThfnmA@tI@D\nfluPP@DTAsUlbbTTTlTvRXHshiZjZfX@`@CCbLLDFJ^P|\nfhyQ@@CA@dsLjm{Jvg^Bt@@DCA@A@\ndmMH`NV@BDke[VFUX@ab@LBpfWbAH\ndiU@@@aJmUJnFjjj@H\ndo}D@DTNrJIQSEQJaZjZjj@H\nfoAP`@DNAsHheEEJefBMYji`@h@B\ndiF@@@rQRJVjj`CC`SBiny@`\nf`qa@@O@rJJJIIKERMjsPTEEQ@A@\nfb}A@@@YEEHdcLeLeB\\djshijjjj`@j`@H\nfoQHA@FJuXFH`HRYUYYbTT[UjZjefZ@B\ndcND@E@TfUvvf`@jh@LEALJnFUt\ndeTD`AdHaIe]jZ@BX@LJALj^PP\nfhi`b@FN|BHRYV~YYirpX@bBA`@@`\ndmtD@@qIUe\\YZ`@`@piFxUxdT\nf`a@P@@HhyIeWfUsVZj`HF@@`\ndklH@@rJJJIEQa`Hbj@C@PIne]OdT\ndeTH@@rJJIEMtAAP@XTFES\\PY@\nfdy@`@@QrJJIEIQISLx{tAAT@E@@D\ndg\\H`BBPRYg^U[hHBjh@H\ndcLH`IBPRYeWvz@`j`@`\ndaD@`@\\DimVz`@@CAdpj[nP@\ngGQDHG@nFRQFj`H\ndet@@DjYUZ^D`dJ@CAdpjz^Hp`\ndaxD@@QInejj@LLRfxa\\\ngOq@@eLvmLtA`l~Ht\ndidH`ACDRYWZZ@B`@pXDpjGd\\\ndaFH@NAIe^fiZh@`\nfoAPB@AVADILkkJ}FFm@AUUQ@A`hBLDiIVcNyC`\ndmvD@LADfvUaej@B@CAliae^P`\ndmvD@LADfWeYUj`@@CAlIne^PH\ndknL`LaM@HrRRqIYPYV`@f@B\ndeVH@BAIf_VzB@h@LJCBinHL`\nf`q@a@AZAD^bDfUmn]hrsh@JB@`@B\ndaD@P@bNbDfUZZ@B@C@`pnxbL\ndcND@NADfUyU]Zj@@@H\ndeV@@@RV[TYzP@@C@j[axaR\ndeT@@DjU_k``R`HipZDJagfXw@\ngGX@@dj|tHkAF|Rp\ndeT@@DjY]zXFB@@pYLinGbEH\ndcN@@@RVYUYwiBB@Hr`\ndif@@@rRQHqajjZ@H\ndazH@LAIYzjj@LFALJnxdL\nfbmPb@B^dDDadTTRbbRRTQXjrXy``jJjZbP@pyc`kNSgfC^b\ndazL`BaL@HrRPyKUPA`a]qJX\ndmVL`BaL@HrRPzJIZjj@H\ndcvL`BaL@HrRPyQZKUUT@P\ngNs@EjpdssTpFFU_I`\ndmtD@@QIn[VUZh@@@pyLkae^Pp\nfoA@`@@BRYfYWuVLyh@`@@@CBUF\\EIVcNyBp\ndmt@@DjYnvDHbJ@CBlkiWdB\ndg}H@ATLbbbbbJpj]mMULm@D\ndg}H@JxDefY^Zy[YjZZZ@H\ndmvH@ACHhhhdYUhJ@@CAbine^Q`\ndaF@`B@HRYg[hH@@LJCJ[nPP\ndg}B@DpAV|bbbRbrK]imUMMU@FGXU]IwrM@\ndcnL`LZ]@HrQQZYPyPkSMUU@D\nfoAQ`@EZ@JVQQQQQIWLTZsT`tAP@D\ndk\\D`HP@cIICDhl^e]V`bB`@`\ndg~D@D@dfVumQd{Z`@@@@H\ndg~D@FCldTRQTRUL]uUT@D@FBpfESKiwrI@\nfoA`b@LRL@HRf[g^|TuZjhDHh@B\ndif@@@RfU~F``@@pxDpj[nPH\nfoAqb@BTYIS`RBSLnmnsakTmTuUP@P\nfoAqP@DXH@RcNrJJIJKJUgMVZjZZj`@pybgIFl{dM@\ndaEH@DHDfYVyje`CAJ{bXp\ndieD@DHFRYf^EjiX@`\ngG]@EbDf]jPLInIT\ndmuD@FxERf_UrjYZjBJh\nfoQP@@@JRi^yWZ\\DhujjhHBh@B\ndo}B@HtDUtfyW[WaZ@BZh@LAPj[mM~QX\nf`qQA@CVHBHEDYEEDeMCLfR]zBBbfh`@`\ndmtH@@rJJIJEn`HJ`@pxBiae^Ig@\ndmvH@EAJYYtYhH``@pYBinGbRh\ndg^L@LAER[eYWNvjhD@@H\ndcND@BADfuYU]Zj@@@LNrnFUwbHP\ndazB`LaFx@cIHThmS@FBXWDC`\ndcND@EADfvUtYZ`@h@LFSJ[axbz\ndcN@pBBlFlVlRYfyVfBBb`@`\ndiFH@BAIUgjjl@pZDrfGd@\ndiFL@JANRY[vjj@H\ndcLL@@G\\bbbRfK]@PUP@P\ndo~D@E@\\bbbReRd[hBBjj@B\ngJ\\@AbeMK@XKGdP\ngJ\\@ABe[LbQpfyhb\ndaz`@@bfyJUfZYBHX\ndmM@PBx@c@aJYg\\jeZdHB@B\nfoQh@@@XdhrTsMlkSdeFluMTsUP@P\ndk}@@@YJUUURkatp@j@@@CCdqiWSyD@\nfnkA`@@X[dRdtbRLrTfVsIZMxIZmUUUUUUUT@D\nflm@@@DiUem^yIRVkNBejjjjjjj@B\nfhyA@@@YEHhheDecKQoNCU@DP@P@A@\nffs@`@@VrQSRIQSJJpqBt[pXykUUUUUUUU@A@\nfhyI@@DDlDBSJkrmmFF]xKUM@DMD@D\ndk]H@LdDie[Wz]N``@@@B\nfbmAP@@\\}GNRJJJJIEJIN`cA\\dBB@@Hb@@CBSNBdkpXyxbs`\nfakAp@@\\}EIVyHhhhhdTheETjBLEvPHH@@bbj@@B\ndk]L@LxDMIe]eRkSZjjjh@`\ndcNDPBePbNBDfYWia``aX@LFCBimy@`\nfnsHJ@CASVBLDxjsP^FRRXjIQYSQIYRuEIZmUUUU@`U@@XFpRcARUjwlm@\ngGX`DJdvsTA`a^YX`\ngNy`HLtemfjDc@\nfdy`a@BRlBHcDYEEDehXhlfB|FBBFB@`@B\ndg}D@DXIrJJIQIH|ivuTB`P@P\ndk]@PAdH`paIe[UZY]`@iiH@`\ndaxB@@QnRUeZj`CBlJf{dP\ndk]H`LiPBDigm[rivfPBJ@B\nfoQAB@C@BDifYU^gIVtz`B@DH@@`\ndmtH@@RVYWeVhH@@LNRfxYWbPP\ndknL@HALRYyyWSZj`@@B\nfoAPb@ARLxDPdrmjktXkT@EMLd@D\ndg|@@DjU_eZx{BAH@@BJlIaBjUt{yfMp\ndifL@NAJR[Y^Ejfd@`\ndk^D@JADfV[W[mMjdHB@B\ndmOH@FyPRf_eriVjjeh@pkBhYxdH\nflmA@@@ILkZvmsQdeV\\EL@EUUUQR@FC`HqS`iJtZwhiy@H\ndkn@`ECDRYWUZf`@fj@cj\ndmvH`BTHaIf]un``J`@`\nflea`@A@bgHhhhhdmL]EVMyh@`@bX@@puQgAJtZsoI@\ndcLHPBCD{DrJJKQFLLDDU@A@\ngJPdLQDHHRViPH\ndg|@@DjU_eZx{BAH@@BJ\\IaBiiWSyf]X\nfleA`@@HcdTtbLRTrVVc^BuU@@@T@@P\ndaEH@LHDimVz`@@CAD{bTp\ndeWD`LjXP@cIHUDdLkP@@A`fMs@xP\ndigH@LxPRfuVz`@`@pQBxeB\nfhiHb@LB|FB@AFRREQQII[Z][ru@@@D@@P\nfhiHb@LBdFB@AFRREQQII[Z][ru@@@D@@X\\psbm^CqAk@\ndigH@LhPRfuVz`@`@`\ndeUD@LhARfuunh@J@B\ndkoH@DgPRYYe]SZ`hH@C@TJfxUOdZ\ndkoH@DEPRYYUfSZ``h@C@TJfz]OdJ\nf`qq@@IZLAJ[WUeV\\]z`@j@@@@`\nfbmP@@@IrRJJIQJIHdjmhKhu@DPQAT@@X^qS`iJtZsoAbgH\ngOt@ATiVKjj`LFHlWrT\ndid@P@bNbDfYYa`H`@LJBfx^Q`\nfhiAP@@XyKRTrlovmA^CUUUP@P@FAhDipTmFl{p^HPp\ndcM@pItIAICICHiCDeLJuPAB@D\ndg_H@FzpRYeUyEZ``bh@H\ndmtH`EBHRYWUih@Jh@LAALJnE^QH\ndcOH@NGPrJJQVIKmUURpAapXUMJ|b@\ndg_H@AVprJQRIVYG]UUUK@D\ndeWL@DpPzTfY[ifjfX@`\ndaF@@@RYe[hB@@LJCBknPp\nfhyq@@DZ|AIe[yVvcCV|EjVffej@B\ndmuH@DTDfYUQUjjj`CClJfxYyB`\nfbuQ@@DJAdTTRRQRRbRDiNVjjjf@B@PMP\ndknD@HALbfbRadxV`@j`@`\nfoAAB@A@bLbbbRebRYrshBB`@@@LM@pQgAJtZsnPh\ngJXHD@aIUj@pHVOI@\ngOq@@dsI[UTA`Ub~R`\nfoAA@@@ILs|kiFZLuU@AT@A@\ndeTD`FdICHhdhUWPAD`A`RYS\\`P\ngOp@DjukZj`LJlZ~R@\nf`iA@@@YEEHdcLeB\\djsjjjjjjP@`\nf`iAB@K@BDifYYWyjSc^h@`HP`@CBPBLEIVcNyCp\nfdeAB@K@BLdTTTRbRsWMR\\[u@DABEP@A@\ndcoH@DGPRUe[^e]Z@HX`CBdJnwbPh\nfb}P@@@NrRJJIQJKEKB\\hvoQsTEQ@PDT@@XNQbgARUjw`qSd\ndeU@@@{IEEEKP`QP`A@\ndevH`Bp@cIEEBdUCMR@D@D\ndmMH`Ad`bDfYuVzU``JR@H\ndcm@PNDHaPaIeU]iaf@Bfb@LFCBE]yB@\nfluQA@N^bBH\\DILrsmlnrUoNCA@qMTlT@D\ndcn`@@rawIEDhiUTkhPUCD@XDc\\oHt\ndeVH@JAJ[e^f``P@LJPfz^PH\nf`iQB@C^X@HrQQJVIJYDijsoUJpACQ@A@\ndeTL@@zTif[rjjjh@`\ndeTL@@SdfUUifiZh@`\ndcLD@@QIe[]ifiZY``Z\ndaE@`BhHaIfUn`H@@ppL{bTp\ndmv@HBBHFPfPVPRYUzih@Jh@H\ndmuD`FFUBHrJIJERn`BI`@`\ndknD@D@tfUe_WViB@`HJpEL[aeSyE@\ndo~D@DA|bbRbbRqvviB`R@`f\ndid@@LdfbQ[fji`B\nfakA@@@YEDeLlhheHmXjL}DNZjefjjjij@aeJprBBLTxJVcN|FJ\\ejwdN`\nf`qI`@DXUh@oYEEEELTlfRtYjjZfi`@piSakNyCp\nfdyHP@DXUh@oQrJJJJXiYILehsUTuMTl@D\ndiDL@@RdfVejj`CCB[nIP`\ndmtHPIBHVHRYfUXXBHX@H\ndieH@JxLbTTQkfej`CAFGbPP\ndeUD@FxJRVYmnYjZ`B\nf`qa@@M@RZYU_tVeFj`hHj`@LM@pRcNBTZsnHex\nflma@@@`rJIJYRIIIFcIF\\GtX@Jjjjbh@LG@Qb`iJtZsoQR\ndcl@@LdbbRRkCtaDQTPA`JXU\\LkrM@\ndid@p@bFbAbDfUfn`BH@H\ndeTL@@PdfzuFVjfd@`\nf`qA@@@ILroljJSoMA@P@@@FBlDipVcV]{dF@\nfhi@p@@XyHpTieY[wzB]zjjhHB`@LCPIS`iJtZs`|``\ndig@@@pldTTqkiZjPB\ndaE@@@qJYfn@b@@ppj[d\\\ndigH@LhPRfuvz`@`@`\nfhiHb@LBlFB@AFRREQQII[Z][ru@@@D@@P\ndeWD`LjXP@cIHUDdLkP@@AafMwD@P\ngJT@ADiZhCCBGbV@\ndaFH`Lx@aJY_JjZX@ppD[dD\nfhiPb@KN|D@PeLrj~lKRcT@@@Tp@A@\ndifD`Na@BLddJT[ejj`B\ndaEH@DpDeYVyjj`B\ndiD`@@iarQQHijVPPJj\ndk}@@@qJYWm\\katvid@J@B\ndknL`JaLCprRRiQIKSZjjj`B\ndet@`@bDfUUiaf@Bh`CB`pagbIp\ndclH@@RYfmj[iZZjj`B\nfoQA`@@\\ldsKKnlyKUgM@pQSD@FALJL{bJB@\nfgAa`@N@t[HhheDTdsdeFltC@UQ@A@\ndmMD@DhNRYYWIiUjA@`@`\ndieD@HXDR[eVyiih@`\ndeU@@@qJYejxBHh@LJJfF^Qp\ndmuH@DXDfUgjZ@Bj@B\ndg}H@NTDefV]qT{Zi``H@H\ndcNH@DCHheEBdnmU@@@FGIeMpkqIt\nf`a``@L@PdwJrkIoMUPPT@@P\ndazD@DAdeYvjh@pKBknIB@\ndefD`BpPbDemgijj@H\ndo~D@D@|bTTTTRqvvhJHb@B\ndmvD`La@|Dj~Uaej@B@B\ndaE@@@yJeVnBB@@`\ndmwH@DHPRYe~[fjVj@H\ndazL@BANR[UZj`CBdrf{dD\nfHgPAa@\nfle@`@@^rJJJIRYqQY`hi`HbjB@@@H\ndaF@@@Ri][jjj@LNaLJfx\ndif@@@rRQHkajjj@H\ndet@P@bBBDfYVX^fBBHPCA`pfx^QH\ndmL@P@bBBLbbbRaZ^FBBbD@`\nfoQa@@A@RYfYWvBUXsf`B@@`@CCeNRUhsnHEX\nfhy``@A@|dsLro~pRkF\\t@P@EP@A@\ngJY@HDeYhRVFFOH`\nfoAAA@F@bGQBSJwroSakPAD@D@@P\nfle@P@@HBeIeUUYVsN}EjjZZjj`HxjX\ndeT@`@dDfYWFZP`@HR`\ndid@`@bDf[Wai@@@LJ@j[nI``\ndmtH`IBHRYm]xZP@h@H\nf`qPB@E^ADILkjlkQcoP@SUUDBGJ\nflua@@E@rQSISIQJEIgEFbf`@jj@@`@H\ndeu@`DhHaIeURhYiZeh@`\nfb}a@@E@RfuUe]UYqQg^Sf`@jjh@H@@pT`eARtZpT|qbNP\nfa{a`@C@hwIEMEDdbdeLeeqQkAbVmAAQUT@EE@@XFPpTeZlxJ\\}y@h\nffs``@C@heMrklkjoIF\\FIZtDDTD@@D@@P\ngFx@@eLzuU@XFHwbF@\nfdyqB@BZ|hDQdTTRRLrbbXKtXHHfjjH@H\nfnsQB@IVDBHrJIJZEQQHjIJgK^jv@HZjj`@H@B\ndcm@@@uJfUWyYvjjjj@LISBinFUxdZ\nfoAq@@DDlAInfUWUXsfjjjjh@H\ndg]D@DDKR[e]UxVjjjf@H\nfoAQB@ARlBHRYVueZLTZ@Bfib@B\ndg}H@DHLbbbQTRS\\]mTl@D@FERng_Lrf@\nf`qP`@DD@iIfUumV\\]yjei@B@@`\ndmuL@DHIWHhhdcFyjef`B\nfhi@P@@HD[Hhdihdibj\\UYjP`hf`DGL\ndiD@@Dj}Yjj`CCdpj[ayD@\nfdy@`@@FrQQRJKZDjTdjuATQQMH@D\nfhi@@@LdbbRRNbRxkpXFBb@@@Ppp\ndmt@@DjU_ZxHHZ@cgChSBjUyc\\H\nfle@@@LdbbRbVQfRyjJXFBhh@@ACE@\ndaE@@@qJYVnBB@@pHj[nHf@\ndidH`D@HRUe^Eh@@@phBkaxdL\ndmuH`LY@BDfymYUj`@@B\ndcMH`DjPBDefulUZ`BH@H\ngOt@AdiWqZY`LCDWrT\ndg}H@NLLbbbTVaSR]pPRTi@D\nfhiAB@A@bLbbRaRbbRhqwh@J`HD@@`\nfhy`@@@YIEEDdXllQbk^CTCP@@P@A`kENRUhsoAxae`\nflu`@@@ILksLjoQchaSUUP@@A@@D\nfbu`c@@drBlEVQkLbbbRbTJRTkI^SAAQAAQ@@FF`XEIVc^CGOHV@\ndcLL@@WTfZYWijjjj@H\nfle`b@LPP@HrRPjIIKQIYQhiVhHHB@@@LCXYpTeV]xLT|aP\ngNq`@fdwKUPFAQkyL\nfhiA`@@HBdrsJrziuoMP@A@H@A`JcARUjwnHLx\nfdya`@H`PEIefUewTzwfh@@`F@@H\nflep`@OQRABSKLjo{JlYsT@@@Tp@AajcARUhsoQSrG@\ndifH@DAIfU[fjV`CBlJf{dB\nfjc@`@@ERYeu[VyV\\Uiwhy`HJjjjhj@B\nfb}A@@@ILk[Kk\\tyKSoQruSUUUUUP@P\nfb}@`@@YRfum[VwYqVg^bf`@jjjjJ`@`\nfhyh@@@XisVRJIJYIPiTDkW`aAQL@A@@XRQabcV]xOHZ@\ngCi@DDfj@ppfyD\ndeT@P@qIdDeYWxY`@`Ha`\nfdeIB@BTuEpIAIfUeUeYrt_A`HJjfhd@LNXEHu`qxcW`\ndg}H`BFPdLbbRRTTmLnp@TuP`FAYiwqDT\nfbmIB@BTtepICHhhdhddhhfsed~C@PUUMUD`A@\ndk}@@@YJUUURkatp@j@@@C@dZUxTYi`\nfhyP`@LRAyIYUVuugIJly`@h@JH@B\nfhyQ@@M^@eMk\\lzpRmN}PA@QAD@A@\nfhyQ@@NVAdbfRfRRTKARUYu@DT@T@@D\nfhyP`@EAAQJ[eWUWdmF|F``bh@@@B\nf`yP@@@FRifUV|cAVTxwj`YBhjH@LITxJRmV]{bLF@\ndk\\H`D@HRUe[Watv`@jb@LA@jzUt~Q@\ndaDHPNBHPHrJIPeUMT`FE@aTwHx\ndeTHPNBHPHrJIPiEUMUR@P\ndclD`HP@cIHXdiepkjp@UD@XLU]J{r@@\ndk^D`La@BLddlRTrF]MZ@Bh`B\ndcl``Ae]BHRYWYZY]`@fhPB\ndg~@PNBHTHRYUueicn@BjZD@pdLKaWSobHp\nfdyPc@IZxHDPPH\\DYEEDiDfeDfJCFBAXjfa@B\ndaF@@@RiUkjjj@LFaBinyF@\ndg|H@@rJJIHqIMqw@PP@@@XKaTwRng_H@\ndidH@@RVe~Fjjh@pzLJfx^P@\ngJT`LPdfvhCAX|b@\ndcLH@@RVYVnffjZh@`\ndeTH@@rJJIFMsSUT@XUeLL|``\ngOt@ADiYqfjPLMB~R`\nfhia@@I@RYeU{UDkpV``bjh`@`\ndcnH@DAIfYwXUuj`PH@LFJne]y@`\ndk^D@DCTfYgUXUMj`Pb@B\ngCi@HDff@pPwH`\ndknDpItpdDdLdLbdLRTtEZh@a`@`\ndmuD@HXDR[fUEV```@LLRexS\\L`\ndk\\L@@STfUm]iiUjP@j@CBlkaOdR\ndieD@DpFRYUrijfh@pILx^QP\ndcnL@LAFRYV]Zy]Zi`@@B\ndmL@@DjYUVGi@@@`@LNpjxYWdH\ndet`@@iiRfYunFVifi@LLpayG@\ndidh@HJfAIfYxZie`B\neMbDBDfp`\ngJPXHlQLQzlA@\ngC`DADZHRVhB\ngC`DAbZHRVhB\ngCahHlNbNlA@\ngJQhHl@cIHUhCBGd@\ngKP`@Ti\\Zj@pRwI@',Nq='daD@@DiYZYji`@\ndaD@@DjUZxHD@@\ndaD@@DjUZxHH@@\ndaD@@DjWjXHB@@\ndaD@@DjWzXHB@@\ndaD@@DjYvxH`@@\ndaD@@DkeVyjj`@\ndaD@P@bBbDfYvzB@@@\ndaDD@@IIf]nZZh@@\ndaDD@@IIf]n``@@@\ndaDD@@YJZUnjjh@@\ndaDD@@aJVdnjjh@@\ndaDD@@yIe^f`@`@@\ndaDH@@RVU[j@@@@\ndaDH@@RYVih@H@@\ndaDH@@RYe[hB@@@\ndaDH`BCDRYg[hH@@@\ndaDH`NBPRYWih@H@@\ndaDH`NCDRYWih@H@@\ndaE@@@YIeZn`B@@@\ndaE@@@aJyUnh@@@@\ndaED@DpFRYVkfjY@@\ndaF@@@RYe[hB@@@\ndaFD`L[`BDi]lJjf`@\ndaFH@BAIf]n``@@@\ndaFH@DAIYUnZjh@@\ndaFH@DAIeUnZjh@@\ndaFH@FAIeZn`B@@@\ndaFH@LAIVUnjjh@@\ndaFH@NAIe^fZVh@@\ndaFH@NAIe^f`@`@@\ndaFH`LP@aIe\\jZjX@@\ndax@H@SHbDbLbLddjUUT@@\ndax@P@SHqDjmZjh@@\ndaxB@@rnRV{Zj`@\ndaxD@@QIgUjj@@\ndaxD@@YIgYjf@@\ndaxD@@iJU^jj@@\ndaxD@@yIUVjj@@\ndaxHpJBPRPrPrJPiUUH@@\ndaxL@@SDfUVjh@@\ndaxL@@idiUZjh@@\ndaxL`HS@BLddNRuT@@\ndaxL`HS@BLddNbuT@@\ndax``HJn@JRgmjZp@\nday@@@yIUVjj@@\nday@@@yIUWjk@@\nday@PDh@c`aIUUZe@@\nday@`Dp@aIfYjj@@\ndayD@HPNRVUZj`@\ndayD@LHNRUUZZPcH\ndayH@DpDfYfjh@@\ndazB@BAJyImmji@@\ndazD@LADf]Vjh@@\ndazD@LADf^fjh@@\ndazD@LCdeYzjh@@\ndazD@NADf{Vfl@@\ndazH@DAIeYjZ@@\nda{D@HQ`iIYejZ@@\nda{H`Dpn@HRYeZjP@\ndcL@@LdbRbceBDEEP@@\ndcL@X@bBbFbAbEbMbDfYn\x7Fijjjj@@\ndcL@X@bBdFdAdEdMdDfYn\x7Fi`bHh@@\ndcLB@@imReUUUvjjjh@@\ndcLB@@riRfY~Qfjjjh@@\ndcLD@@IIf]^xViZjP@\ndcLD@@uIfUk[hBBj@@\ndcLDPDtHaXaInUvxY`@i@@\ndcLDPEt@b@cIDhTdiaUUUM@@\ndcLF@@Rag\\bbTVTILuSUT@@\ndcLH@@RYeUqVhHH`@@\ndcLJ@@PUuInUgzV`BJ@@\ndcLL@@STfue^UZX@HBL`\ndcLL@@STfyWWaZ@Bh@@\ndcL`@@qaRfU~ZxHDj`Hj@\ndcM@@@aJufYgZhBH@@\ndcM@PEtH`haIf]u[hHBj@@\ndcM@PEtHchcHhheETppDQT@@\ndcMD`LE]@HRV[]nEjiih@@\ndcMH`BuPBDf[U{aj@BX@@\ndcN@@@Ri]mUvjl@@@@\ndcNB@BAEuInVVFV`HF@@\ndcNBAHAEvISdfyV{iZ@HX@@\ndcNBAHAEvISdfyW[aZ@BX@@\ndcND@BADf{YU]Zj@@@@\ndcND@E@\\bbRafUM@AUP@@\ndcND@E@dfYwYn``Jh@@\ndcND@MBldTtRROSPPQP@@\ndcNH@EAJ[WU[j@Bj@@\ndcNH@FAIfYWgZ`hH@@\ndcNH`Bp@aJY{UWZZjj`@\ndcNL@FA]rJJJJJlLADU@@@\ndcNL@LB]RVV]vE`BFP@@\ndcNL@M@iRYg^vzB@j`@@\ndcNLaIQ]BHf]yInUvxY`@f@@\ndcOD@DgPGHheDhkrmUSU@@\ndcOH@NFPrQISQHrmTAA@@@\ndcOLAHHPSXeNrJZJIHlKTuUH@@\ndclD@@GIEEHhfUKmAU@P@@\ndclD@@UIfV][iuhFAH@@\ndclD@@iJYW]rnF``JX@@\ndcll@DpYAmRYV]zzUZije`@\ndcll@Dsm@iRYgeVE]ZjeZ`@\ndcm@@@YJYYwhUtH@@@@@\ndcm@@@{IDeCDdUKh@UUD@@\ndcm@@@{IDeCEDUSh@UUD@@\ndcmD@DHERYYUZz]Z``R@@\ndcmD@DpFRYV]Zy]Zi`@@@\ndcnD@D@TfYg\\fWVjYjh@@\ndcnDPN[PBABLdTRbauCtuTuSP@@\ndcnL@LAFRYV]Zy]Zi`@@@\ndcnL@LAMRYUUrhYZh@h@@\ndcnL`LaA@HrRPjIKTrzmPHD@@\ndcndADkaTMz]yIefU[iUjeji@@\ndco@@@bdigyWJWYj`@P@@\ndctB@@RURY]VvjjZ@@\ndctB@@ieReUWZjjj@@\ndctF@@rngTen{mjjj`@\ndctFPHSNePBABLddNRTbuUK@@\ndctH@@RYuUVjjj@@\ndctHxJBPRPrPFPfPVPvPrJPqJQUUUT@@\ndcu@`Dp@aIeUUZjjh@@\ndcuFPBvDsiT@`PcIICdeHmURp@@\ndcwD@HQ`iIYgmZfjh@@\ndeL@@Di[ernDYZjij@@\ndeL@@DjYeIjGijjjj@@\ndeLD@@EJUWhRfFZjjj`@\ndeT@@DjWvifjih@@\ndeT@@LdbTRoBuUM@@\ndeTD@@EIe]jZ@Bh@@\ndeTD@@QImeQej@@@@\ndeTD@@eIfu~Eh@H@@\ndeTD@@gHhhhjppDQ@@@\ndeTD`AdHaIe]jZ@BX@@\ndeTD`NDHaIfVVfBA`@@\ndeTH@@RUYTYY`@@aH\ndeTH@@RUYTYi`@@aH\ndeTH@@RV[TYjP@@@\ndeTH@@RYVZfZZj`@\ndeTH`IBHrJJJJlLADP@@\ndeTL@@PTfyVzV`B@@@\ndeTL@@QdfygFV``@@@\ndeTd@HRi@TefUzVjih@@\ndeU@@@aJWeQfj@@@@\ndeU@@@aJueQfj@@@@\ndeU@`Dp@aIgeQej@@@@\ndeUD@DdARYgxfZjf`@\ndeUH@HPDefuFVh@@@@\ndeUH@JdDin_xZB@`@@\ndeV@@@RgYTYj`@@@\ndeVD@BADf{YxVjjX@@\ndeVD@D@TfVgJfjih@@\ndeVD@DAdfygFV``@@@\ndeVD@FADfygFV``@@@\ndeVD@IADfyWxV`@`@@\ndeVDAHAHeNR[e[aZ@B@@\ndeVDAHAHeNR[e_aZ@B@@\ndeVDAHAHeNR[fTYZ@`@@\ndeVDPJJPqFqDfYokjVjX@@\ndeVDaBxPbBd^RYeei``X@@\ndeVDaNFPbNf^RYWZf`@f@@\ndeVH@IAIe]ZZ@Bh@@\ndeVH@LAIUeQfjjj@@\ndeVHpDxHbhbXcHheEJptuMH@@\ndeW@@@RTfyWxV`@`@@\ndeWH@DZPR[e_aZ@B@@\nded@@DiUUjjj@@\nded@@Dj_VfZZ@@\ndedB@@RiR[UUijhHr@\ndedB@@iiReUvjjh@@\ndedD@@QIkWZjj`@\ndedD@@QInUVfj`cH\ndedD@@QInUVjj`@\ndedD@@aJVfjjj`@\ndedH@@RUUVjjh@@\ndedHXDBHjPZPzPFPfPRYZZjjh@@\ndedHXJBPRPrPzPFPfPrJPsRUSU@@\ndedH`HALRkUUjjh@@\ndedJ@@rneI[nvjj`@\nded`@@iiReUvjjh@@\ndedd@BXYADfyvZZfBL`\ndee@@@gHeHiCuUV@@\ndeeD@BdDR[mUjjh@@\ndeeL@HdBEIVUfji`@\ndefD@LADf^]Zjj@@\ndefD`FFPBDiWnjjf@@\ndefH@DAIfWvjj`@\ndefH@DAIfuVjj`@\ndefL@IAAR[UYjjX@@\ndefL`F`Y@HRkvyjjX@@\ndefL`JHY@HrJJEJUMS@@\ndeflAHrfxDfISdfvyZfi@@\ndegD@HZPEIUUfji`@\ndet@@DjYUX^d@@@@@\ndet``Dki@HRYYUnFVjVi@@\ndet``DkiBHRYYUnFXBBa@@\ndeth@DkiAIeeVxYZiZd@@\ndeth@DxYAIeeZxYZfYh@@\ndeu@@@kIEEDceCLaPD@@\ndeu@`Dp@aIeURhYfh@@@@\ndev@@@RYyULFZh@H@@\ndev`@@rfeJY{ZxYBBJD@@\ndew@@@pldTTJVTLmP@P@@\ndg\\B@@Q[R[VUmgVf@HhBL`\ndg\\D@@eIfU_Un`HJj`@@\ndg\\H@@RYWY^ih@Jmh@@\ndg\\H@@RYvYiVvfZjjBL`\ndg\\d`LF[a@BLddJbbQvfmPLA@@@\ndg\\h@Fd{AIfUWxUZABI`Hj@\ndg\\h`LJfd@aJ[Ywkaijejh@@\ndg]B@NTBNtfY}[fyjVjf`@\ndg]D`NV]BHRYWVyih@JZh@@\ndg]HPAuPbBbDfYw[fzB@ij@@\ndg]H`AuPbDfY_[fz@`ij@@\ndg]L`LnDD@cIIBhhd]ikTC@P@@\ndg^@@@RigvuNzjk@@@@\ndg^B@BAEoHiihdhUNmT@QP@@\ndg^B@BAMoHiieDeBimU@DP@@\ndg^DPNI`BDBDfumeSmZf`@@@\ndg^D`INpbLbdTJTv`kU@DS@@@\ndg^D`La@BLddLTRRIvmTE@@@@\ndg^H@DCHhhhddYimTE@P@@\ndg^L@EACR[YWVFVh@JX@@\ndg^L`Iia@HrIQPjZIG]TtuU@@\ndg^LpDx{BHjHVHrJIQPiZLLADUL@@\ndg_@pNfpqDBCBDfUyWbjfijfPad\ndg_BAHJ`QSbT{HiheBidwMMUUR@@\ndg_L`LfxPPBLddJbbQvfmPLA@@@\ndgl@P@SHbDjmUUZjjjh@@\ndglB@@Q]RYUumZjjZ`@\ndglFPHSItpBEBLddNRRTbuUTl@@\ndglP@@anESoIDhhd\\iMSMS@dWP\ndgll@DrSAmRYV^uZiik``X\ndgmFPBnDr]L@aPcIICddeHmUUK@@\ndgmL`DVCl@cIICDmihmUUS@@\ndgnNPKaLXYTOBoCIIKHiMTuUMU@@\ndgoL`LD~gPBLdTtTjTuMSTt@@\ndg|@@DjU_eZx{BAH@@BJ`\ndg|@P@bEbLbbbTVbKR]pP@P@@@\ndg|D@@MIfV]UivviP`R@@\ndg|H@@RVUvU[cm`@`@@@@\ndg|H@@RVUvU[cn`@`@@@@\ndg|L@@t|bbbbbJcBmm@tCD@@\ndg|d@Dq]@\\bbbbfJSSimUSTs@@\ndg|laLJce]@HjUyJYn[gFNzZZYj`@\ndg}DPEtYBHXHrJIJZQGUvg@DKUB@@\ndg}H@IlLdTrTraS]Jt@EMQ@@\ndg}L@hdDNISid~R[fV]F]mh@@H@@@\ndg}L`JXiTIAIf]VunNzVjZj`@\ndg~D@D@\\bbTbTjQL]mUR@D@@\ndg~D@IB|dTrTraS]Jt@EMQ@@\ndg~D`LS`BDfUueRh{Zjh@H@@\ndg~L`EF[BHRYVyukev@HZiH@@\ndg\x7F@PBWPbAbLbbbRfaSR]pPTMQ@@\ndg\x7FH@DjprJIQEQIMqvuRp@P@@\ndiD@X@dDdLdJdNdAdLbdLtjjh@@\ndiD@`@RdjeVjj`@\ndiDB`HSB@HrRPiIZj`@\ndiDD@@QImVZjh@@\ndiDD@@QImZZjh@@\ndiDD@@QIn]Zjh@@\ndiDH@@RVUjjj@@\ndiDH@@RYWZjj@@\ndiDH@@RYuVjj@@\ndiDL@@Pdf]ijj`@\ndiDL@@RdfUUjj`@\ndiDL@@idiUVjj`@\ndiDL@@xTiUVjj`@\ndiDLB@PTyInuZZh@@\ndiE@@@aJuUjjh@@\ndiEL@HPJyIYvZjX@@\ndiEL@HhNEIgVZjX@@\ndiFD@AADfuUjj`@\ndiFD@BADf{Yjj`@\ndiFD@FADfumjjP@\ndiFD`JxPBDivzji`@\ndiFD`JxPBLbdTljjX@@\ndiFD`JxPQDivzji`@\ndiFHPF`@`PaJoUZjT@@\ndiFJAHALXXgNRYvvjj@@\ndiF`@@SNGHhhcFjf@@\ndiG@@@HTeUWjjp@\ndiGD@HQ`kHeECFij@@\ndiGD@HhPyIUVZjX@@\ndiGL@HQ`h\\bTTQZfh@@\ndiT@@DjV]nxVjihHB@\ndiTH@@ReUxRfjjf`@\ndiTH@@RgeXSaj@B@@\ndiV@@@RfU|kadHB@@\ndid@@DjYUaBHP@@\ndid@@LdbRQk``b@@\ndid@`@bDfUff`@h@@\ndidD@@EIfU[ffiP@\ndidD@@EIfU[hBB@@\ndidD@@GHhdhZZ@B`@@\ndidD`BD@aIf][hHB@@\ndidD`BDHaIf][hHB@@\ndidD`BDHaIf_[hHB@@\ndidD`BDLQIf][hHB@@\ndidH@@RUe^Fh@@@@\ndidH@@RYVZZ@B`@@\ndidH@@RYWZZ@BP@@\ndidH@@RYm^Eh@@@@\ndidHPNBHPHRYWbjfid@@\ndidL@@IdfYoa`b@@@\ndidL@@XTfUfn`BH@@\ndie@@@EIfW[hBB@@\ndie@`BDHaIf][hHB@@\ndie@`Dh@aIefkfif`@\ndieD@DpFRYVZyjfd@@\ndieD`JXaBPRYgvzejX@@\ndieD`LIN@HRZufFZid@@\ndieH`DpPBDfYeaZjj@@\ndif@@@rQQIFfZjj@@\ndif@@@rQQIFfjjj@@\ndif@`D@HRUe^EX@@@@\ndifDAHAHeNR[e^Eh@@@@\ndifDPJHPqFqLbbbekjVi`@\ndig@@@xTjU^njjj@@\ndigDaLHVx@cI{dieljjUj@@\ndigH@DK`R[e^Eh@@@@\ndkL@`@SDjUUUjjjj@@\ndkLD@@QIe[WVijj`@\ndkLFPHSAWPBIBLddNRRdVjjV@@\ndkLF`DkaTpBLbdtVafjjji@@\ndkLH@@RUUUVjnzl@@\ndkLJ@@imMJUUUZjjj`@\ndkLP@@anF]OIDhhjeIjYi`RKh\ndkLP@@bfF]OIDeDTeIiii`RKh\ndkMFPBNDpUt@bPcIICddiEjje`@\ndkNLPFi]@HLHrJQRjJJjZih@@\ndkNLpFjUBPRPrPrJPsISJjZfd@@\ndk\\@`@bDfYYwZ]NB@@@@@\ndk\\D@@wHhhhhbfESZAhD`@@\ndk\\D`ELLSHhdidhZZ]`@iYH@@\ndk\\H@@RYeg]itxH@@@@@\ndk\\``LyS@HRfUeTUtzBHba@@\ndk\\d@Dq]@\\bbbbfJZ]MjjZe`@\ndk\\d@DrUCdfY[uXUujijY`@\ndk]@@@uJYU_lkihDHj`@@\ndk]D@DxJRYY{Tjtvfh@H@@\ndk]DAHDMF]yIYeg^ESh@@H@@@\ndk]H@ALLbbRfTJiaf@Bfh`@\ndk]H@DHDfYyWImMjhHA@@\ndk^@@@RfYU\\]Tz@@@@@@\ndk^@@@Ri_YVftzjX@H@@\ndk^@PDCDNHRVf]}aWZjB@h@@\ndk^D@A@|bbRfTJiaf@Bfh`@\ndk^D@IADfvYWz]MjdHB@@\ndk^FAHALzSb\\{HhmEDdQe]Zjjjh@@\ndk^HPNHJrXaIf^UvGS``B@@@@\ndk^L@LANrJIQQEJatvjf@@@@\ndklB@@PcR[me]]Zj@B@@\ndklB@@QSrJYJIJF]ZX@b@cH\ndklB@@QmR[fUxUZBBF@@\ndklD@@QInUvnEh@Jh@@\ndklD@@eJ[Vvfz`@jh@@\ndklDpBhJqJs@cHiiCDeeviX@H@@\ndklH@@RYfVya`Hbj@@\ndklH@@rJJJQJy]ZBhB@@\ndkl`@@SFRYvUWSZjjj`@\ndkl`@@kaRe[vTf@HZj@ah\ndkm@`ATHaIe[ujZ@BfhBA`\ndkm@`MLHcHhhdcDfz@`Z|@@\ndkmB@hHD[heNJ^{HihhhTmMjf@@@@\ndkmD@DpJRY{YWSZf`@@@\ndkmDBLdDf\\bbTTJTWVj`@`@@\ndkmDpDgSBHjHVHRYYmYn`HJf@@\ndkmH`NVPbDfUunih@JZ`@@\ndknD@DATf^UuFVh@J`@@\ndknD@E@\\bbRafRih@Jj`@@\ndknD@LALbbRbaRtvjh@@@@\ndknD`J``BDjy{Utvfjjh@@\ndknL@BACR[me]]Zj@B@@\ndknL@BAMR[meYSZj@H@@\ndknL@M@iRYg^un``Jj@@\ndknLBLANSlbfbTaRtvjzjh@@\ndknh@HF]LDee]V[iiYihDbz@\ndko@@@aLdRbbRQtvjB`@HB@\ndkoH@Dp`RYuYWSZj`@@@\ndk}@@@iJUmUR[atpBJ@@@@\ndmLB@@RURYUVJaejVjh@@\ndmLD@@QIe[VfeVi@B@@\ndmLH@@RYiYKnUjjjh@@\ndmLH@@RffeZVVjjjh@@\ndmLH@@rJJIQEneX@@@@@\ndmLL@@SdfVUrjUZ`PH@@\ndmM@@@[IEDhhZFU@B@@@@\ndmM@@@aJUe^neYhDB@@\ndmMD@DdNRYYWJiUjA@`@@\ndmMD@DhNRYYWIiUjA@`@@\ndmMh@DkaePRYYe[iUifjd@@\ndmN@@@RYVuiiVj`@`@@\ndmN@pN@H`HPHrRPqIZneUhDB@@\ndmND@BA\\bbbReInFjZjd@@\ndmND@DBdfV]RZUZf@@@@\ndmND@DCdfVUrjUZ`PH@@\ndmND@DCdfVUrjUZfZj@@\ndmND@DCdfVUrjUZjZj@@\ndmNDPNYPBABLdTRbch^fjfih@@\ndmNH@NAIYe^neZdHB@@\ndmN`@@raeJY{VfFP`JX`@\ndmN``DkaT@aIefUneVjVjP@\ndmNh@HqaTDeee^FUhBBa@@\ndmOD@DiPyIee\\feVhDB@@\ndmOD@DjPyIee\\feVhDB@@\ndmOH@DK`RYYWkiUiB@`@@\ndmT@p@SHbDbLddJRRjjj`@\ndmTB@@ieReUWjjj`@\ndmTB`HSBCpRjuUZjj`@\ndmTBpHin@JspDJrQIRGJjZj@@\ndmTFPHSFFPBNBLddNRdVje`@\ndmTH@@RUUUjjj`@\ndmTJ@@IaUIf[~jjj@@\ndmTL@@xTiUUZjjh@@\ndmTLPDp`|HdLddJTTfjj`@\ndmTl@HQeBiRVY~Zfi`@\ndmU@HLT@a@c`bPaIggYjjf@@\ndmU@pLD@a@c`cHheEKFjfh@@\ndmUD@DXARYW[ZffP@\ndmUH@JXDiUYjjjh@@\ndmULAFDDIiGdf{evfYh@@\ndmUL`DzIT@cIICeMEjjX@@\ndmVB`NiiT@cIDeEKJijX@@\ndmVH@DAIYUUjjj@@\ndmVL`JFU@HReYyjjf`@\ndmWH`DphCpRjyfZjj`@\ndmtB@@QeR[f]FV``H@@\ndmtB@@QeR[f_FV``H@@\ndmtD@@UIfU|YZB@`@@\ndmtD@@WHhhhiVF@bJ@@\ndmtD@@iJ[g_ahHBP@@\ndmtD@@kIEEHeYV`j@@@\ndmtH@@RYeY[hBBh@@\ndmtH@@Rfuu[j@BXBA`\ndmtL@@SdfVyyUjB@@@\ndmtL@@YTf[gqehHB@@\ndmtL@@yTee[VE`BJ@@\ndmt`@@kaRe[vI`BFhBF`\ndmu@@@aJueXYj`@`@@\ndmu@@@aJue\\Yj`@`@@\ndmu@@@aJufVUj`H@@@\ndmu@`ITICHhhdbfz@`j@@\ndmu@`ITLSHhhdbfz@`j@@\ndmuD`LVD@HrRRqIXYV`@`@@\ndmuDaBWaBHJQ{HhheEVfBAb@@\ndmuH@HPDefUyUjB@@@\ndmuL@DTAeIf]XfZjih@@\ndmuL@DpFUIeY~nZijd@@\ndmv@@@Rf~UeZj@@@@\ndmvD@DATf^Uqej@B@@\ndmvD@E@dfYwVzB@j@@\ndmvD@EADfyW^Ejjj`@\ndmvD@NADfVyyUjB@@@\ndmvDPDePbAbDfUmjZ@Bf@@\ndmvD`EJPqDfY}VzB@j@@\ndmvD`La@BLddlTReUhB@@@\ndmvD`NePBDig{ljfZf`@\ndmvH@FCHhhhdYUhJ@@@\ndmvL@FAIR[f^EV```@@\ndmvLaFFUCDZYyIe[fn`BF`@@\ndmw@@@aDiUmYUj`@@`H\ndmw@@@aLdRbaReVj@@B@`\ndmwDAHEPRISlbfbRazV`BH@@\ndmwDAHYPRISdfygQehHB@@\ndo\\@T@dDdLdJdFdAdEdCdKdGdLbdLaffdjjjjh@@\ndo\\J@@imMJUUUUjjjjh@@\ndo\\L@@idiUUUVjjjj`@\ndo\\NPHQevw@DHDrRYJKJHqjfii`@\ndo]B@L\\Dv|bbtTTUTZjjjX@@\ndo^J`CaLKP|LddqbbTlZijfhHi@\ndo^bAHs[\\DxXeNrJYQQEIQjiji`@\ndo|D@@QIn[WVeVj@Bj@@\ndo|D@@uJYU]]zZBBJj@@\ndo|H@@RV^UviUj`@j`@@\ndo|H@@rJIKJQFRf`@jjh@@\ndo|H@@rJJIHiIIn`HJjh@@\ndo|J@@S[_HheDdeDYMjBBb`@@\ndo|L@@iTinU]_ihHHjXBC`\ndo|L`BmpbDfYU]Ya``biX@@\ndo|hHLJkL@bYAYCYCIEEMMHi|jfZjfh@@\ndo|h`LInT@aJYV}zFZZjVj`@\ndo}B`K^DxpBLddLRRVbgUhHHZ@@\ndo}DPAuSBHJHRYg]nVzB@ij`@@\ndo}DPAuSCDJHRYgunVzB@ij`@@\ndo}H@DpLbbbLRVTWVj`@j@@\ndo}`@@pjX|dTTRfb\\k`Hbjj@@\ndo~B@BAC_HiiMhdh]mjj@H`@@\ndo~B@BAK]InVuf[fjjif`@\ndo~D@DAdfuke^wZjhH@@@\ndo~J@BAFM|bfffbRavvjh@b@@\ndo~L@CB]RfyV~^F`HJj`@@\ndo~L@GBnrQSQQIGHUhHbJh@@\ndo~L@LANRYye[iujBBJ`@@\ndo~L@M@iRYg^ufzB@jj`@@\ndo\x7FBAHXpSAb\\yIgUnUmvji`@@@\ndo\x7FD@Ds`gHhhcDliIUjiijh@@\ndo\x7FD@HPpEIYunUmvji`@@@\neFA@lbrJ@\neFA@lhbL@\neFAADdRJ@\neFABHhbL@\neFABljrL@\neFACDlRL@\neFB@Jc@@\neFBBHc@@\neFBBPc@@\neFBCDc@@\neFDBb`@\neFHBJ@\neFHBL@\neFJHbHh@\neFPBc@@\neFhXNic@@\neMA@HPaIT@\neMA@JXaIh@\neMABPYAIh@\neMBBHRZ@\neMCALhabHz@\neMDARU@\neMDARZ@\neMFI@bMP@\neMFIBrMP@\neMHAIX@\neMHAIh@\neMHBN`@\neMIHdFPRZ@\neMPAR^@\neMPBch@\neM`AIx@\neMdXZWde`@\neMhDRV@\neOBBHcfh@\neOHBNZ`@\nfH@@\nfH`P@@\nfH`T@@\nfH`X@@\nfHap@@\nfHbd@@\nfHbh@@\nfHcP@@\nfHcd@@\nfHdH@@\nfHdP@@\nfHe@@@\nfHep@@\nfHf@@@\nfI@@\nf`aA@@@ILsKWRpTADUUP@@\nf`aA@@@YEEDdTddf\\`HJjj@@@\nf`aA`D@HdclbfdQTtRSoMUTp@@@@\nf`aI`@BJ]xBAIMrlvkYsUTuRt@@\nf`aPP@B\\@bwlbfdRLRTIgMUT@Q@@@\nf`aPP`M^djs`QAHa`PXHrJJJYRZr\\EADBtl@@@\nf`aPRDL^`cQ`AAnrRPjIIQUJlmPPLT@@@\nf`aP`@DJA{HdiEDXhkKRuSUTl@@\nf`aP`@LV@cHhmEDidYISTAAUP@@\nf`aQB@NJdBHRYWVyncH@JZjp@@\nf`a`C@C@dDRFICHiCDeBhlJUT@EUP@@\nf`a`P@B@TZvQQQIIRsHTfiZfZd@@\nf`a``@M@PdwMkjteMT@ES@@@\nf`aa`@D@IiIeyWUYJZh@Jj@@@\nf`aa`@M@D[IEDdhbeiqUAAQUP@@\nf`apRBNJMACo@BES{lddJVbbbfcKUURsP@@\nf`aq@@DV\\CHheEDcddkSTE@UP@@\nf`i@P@@HD[HhdihdhUSbkN|uHPTUB@@\nf`i@P@@HTYIe[VUZLeiwfi@HjhP@@\nf`i@a@FRAD^bDfUm[WirRkN`BJ@BH@@\nf`iP@@@TRfY]{WEAJ]zPhH@B@@@\nf`iQA@B\\|@HpDISLzsnRdcN}TaAQUD@@\nf`ia`@D@M{HhhhddYdJVkN|uPTDUD@@\nf`q@@@DjYU[oxjQh@@BJ`@@\nf`q@@`@QAHbtQzHrJJKKQIPlxYtDA@DP@@@\nf`q@A`@XaLQfHKDMb^qDfYn}gyJshHb@``@@\nf`q@`@@HR[YWYTTg^Z`@`@@@@\nf`qA`@@HddrsLk~gV|u@D@@@@@\nf`qHA@DBMxDQhHRYWVY~cMV`@jikHDMH\nf`qHP@DXxHDg^rJIQSPiITyhsSUUUS@@@\nf`qI`@DXUh@oYEEEELTlfRtYjjZfi`@@\nf`qI`BHDDhBOPeGIMsOljpSoMSUUUL@@\nf`qPA@B^ADVbLbbbrtRTKNF]A@PAD@@@\nf`qPHAHV@cGIF]zDxzt~{HhmEDdThJRmMUUUST@@\nf`qPP@LT@cfdfV^V]IRwfjhJBD@@\nf`qP`@DV@YIfYmYZLuyj`hFb@@@\nf`qP`@DXAqIfWnUrL]yjihHB@@@\nf`qPaABLEh@Q@HIKIEvnrQJJJIQXj][uPAA@@@@@\nf`qQA@CVHBHEDYEEDeMCLfR]zBBbfh`@@\nf`qQB@G^`@HRkfUWvQmVZ@@@J@@@\nf`qQ`@BF@abS\\oLkjL{sUUHPD@@@\nf`qQ`@BZ@crS]lrnruYsUTa@T@@@\nf`qQa@LL`bp@b`QddebRRRrxKQejZ@J`@@\nf`qQb@FBt{p@aJYUW]lTEJfZjZZ`@@\nf`qQb@LD``P@cIICDhhlUrSoKUST@@@@\nf`q`P@L@PZvQQISQIPiVg^ZhBBjH@@\nf`q`a@APX@H}HIUwOJ~JUgMRtDAP@@\nf`q`a@IN|BHDDILklnmQbmP@TuSB@@\nf`q`b@ILD@HrQQJJMQITX{uMTt@D@@@\nf`qa@@E@rQQIZYQF\\EjuA@TuQ@@@\nf`qaB@AFADYEEDeHUDV\\]z@`jfh`@@\nf`qaB@BXAHYDihiheD]F]zfh@@@@@@\nf`qa`@D@x{HhheKDeeQakMUL@EP@@\nf`qa`@L@QqIgfUnyZTYjBBJh`@@\nf`qi`@DTxIPC^rJIQQJiILyISUKMUU@@@\nf`qp@@@Hpds\\rj~gV|uP@@@@@@\nf`qpP@DXxBSoYEDhihTdjBtYijjji`@@\nf`qpa@LR}A@@`PQddaTTRRUtxuejBPha@@\nf`qqaBLDhJR`AA`bdr{IEELiiXeacoTtsMMR@@\nf`~@P@@HtYIeU]mVjjij`@@\nf`~@`@@^rIRHiIJQUUUUU@@@\nf`~@p@@HYHsdf]UUVZjjjf@@\nf`~@p@@TYhwdiUUU]jjjjj@@\nf`~@q@HHpjs`AA`cIIChlhmHmUUTl@@\nf`~A@@@ILjjjjuUUUT@@\nf`~PA@F\\@DPBDiUVYUjjjjj@@\nf`~`J@NPQafcN|CprRRjYQQFMTuSS@@@\nf`~a@@D@rJIIJHiIMSUMU@PRe@\nf`~hPBHXhsPDV]BT\\dvvjw[TuUSP@@\nfbc@@@LdbbbbbbQQsEBMKuhir@@@@@@@@@@\nfbc@@@LdbbbbbbcJwEB]Hu`ir@@@@@@@@@@\nfbca`@D@BeIfYeUuu`iJMytTyZh@P@bjB@@\nfbca`DK@iftiSJzljlxipVk^bEMA@@@DU@P@@\nfbe@@@LdbbTbsbrTWESMUSSML@dESY|`\nfbe@P@@LU[HhhhehcBdeHsPTEEUT@@@\nfbe`p@L@PkPTiLjlr{{V|uU@pQS@@@\nfbei`@LTDjpDIRYUyWYVtYjj@Bjj@@@\nfbepb@BZ}IPHaIf]V_VUgHHBfjjh@@\nfbepb@BZ}IPLQIf]V_VUgHHBfjjh@@\nfbexc@LDXXJ\\pP@HdD^BLddNfRRbRrcoKTuUKTt@@\nfbm@A@@IdHaJkfYUU|tWICfj@@@@@@@@\nfbm@P@@BSgHhhhheEl\\hHrgV|t@PDASP@@@\nfbmAB@H@SDjnYeUWsQ\\dNZh@@@@@@@@\nfbmAB@H@SDjnYeUWsQ\\dNZjjjjjjh@@\nfbmAP@@BLENQQQQQQQYG[bm^]Eh@bHAJ`@@\nfbmAP@@HeENQQJJJIKZIUgMNmyj@H@Bj`@@\nfbmA`@@ZMdRbbTVTJRRxIwdyXBBjZ@@@@@\nfbmHR@LLDhwdy@DYHieDdUDhlnBdipZ@`Zjjhh@@\nfbmIR@AZ|fBKNb@HrRPrJI[QQETEHr\\kUUUMRmP@@\nfbmIR@LLDiAoIr@HrJSJIHjIQY\\EIS`uUJuUUUP@@\nfbmP@@@FRieYgUWxYrPyjj`B@HH@@@\nfbmP@@@IrRJJIQJIHdjmhKhu@DPQAT@@@\nfbmPB@AIADILklsJktYyL\\p@UTEAA@@@\nfbmP`@DJ@{HhhhhhdhbRqUhaRuP@pPAP@@@\nfbmP`@DJAgHhhhhhdhcjqUhaRuP@pDAP@@@\nfbmPb@AJ|dDPdrmrljoSfgQs@DMPTBD@@@\nfbmQ`@DX@sJSJmsJkyJNbEKUM@@@E@@@\nfbm`B@I@bDfYYu^yV\\uYNXB@P`B@@DAP\nfbma@@I@rQSQIYQQIH|ExJ\\tAAUAPA@@@\nfbmaB@KDADILsjlmklxZr\\pPEUPPA@@@\nfbmhR@DXXJQoIrBHrJIJKJFIQYTyIS`pAEUUUPp@@\nfbmiP@DTxIPC^SgHheEDjheEdsdeNBuRsUUUU@@@\nfbmiP@LF]EHIJugIEEEDdhjllRegAcUUUPHET@@@\nfbmpC@KQRBHDDZBLbbbRaRtTRSNJtgLAAUMUjTBFh\nfbmp`@DXHBRSL|lvj}F^SGKUUH@@A@@@\nfbmp`@LFbAFQQISISJQIFRcVbej`@jjjJ@@\nfbmq@@ARLAJ[eg_fWbmFSf``Xj@`D@@\nfbmq@@DBBAIeeU{e^VgASejBBIB@H@@\nfbmq@@IJLAJ[V{WeY`mFbf`BFj@BD@@\nfbmq@@IJLAJ[V{WfU`mFSf`BFj@`D@@\nfbu@P`HHdhCaapXxZ\\LdTbdQLblRRmIsUUUUP@P@@\nfbuA@@@ILsOsJk~SGKUUP@@@@@@\nfbuAb@KIrBHrJJIJKIQGJLx[tAAUSULP@@\nfbuHa@LRlFB@A@`cIIBhhddmeIiqkKTDaQUP`@@\nfbuI`@DXtX@dyEEEEDbddUFRLyjjZjZfX@@\nfbuPP@DD@xL\\bbTRTRVReUNRLu@aQTuP@@\nfbuPb@GQxDDPdrnlkk]QakP@T@EMP@@@\nfbuPc@IFSdD`bPqHYEHXdimDheAJCFj@H`HF@@@\nfbuPr@NBYhwdyAbIJsNk[NbdxLtmTuUUPHP`\nfbuQP@DA@JR\\yEEEDhdUmMBLLyj`VBJf`@@\nfbuQP@DTAqTTyEDhhkMEDhfBCejii``X`@@\nfbuQP@DV@JR\\yEEEDhdXceBLMyj`VBbf`@@\nfbuQ`@DX@ZvQQZQQQQIH{rXyZj`A`@@@@\nfbuQa@EVEHpHc`QdTRVTaTrbRhqNX@IZZjh`@@\nfbu`Q@LJMxH@cHQdbbRbTTRfuqQrZj`f`ZJ@@@\nfbu`R@DFmzH@aIfYU^]^QmNZhJBbjh@@\nfbu`c@IQRBPQHXdLbdLRTvbbTPeAsU@DPPL@@@\nfbuaP@D@M[tTfYeUyuyFtyj`hJJj`@@\nfbuaPBH@PrQhRclbfbbtQfbRpeQRtDDS@A@@@\nfby@A@@IdHaJkUUUUVjjjjjj@@\nfby@`@@HR[UUUUUZjjffj`HTvP\nfby@`@@YRUUUUUUjknknj`@@\nfbyPH@BZ@bR`qSgHieDeDXdhhuUSUTs@@@\nfbyPq@BY`cGhi@DABLddNRRRRRTbuUUUTl@@\nfbyQ@@DX@dsLjjjjuUUUUU@@@\nfb}@@@DjYee\x7F]^RD[S`q@@@@@B`@@@\nfb}@@@LdbbbbbTqvwEBlXJ\\`@A@D@@@@@\nfb}@`@@YRYVumUWhrRkNbf@@@@@B@@@@\nfb}@`@@YRYVum[ehrRkNBf@BjjjjJ`@@\nfb}@`D@HQvQSRJIJUJYHRiZ]yNVjjjjjjj@@\nfb}AB@L`BLbbdRQfRfRqNRUYtTuUUUP@UP@@\nfb}A`@@JsdbbTTRrLrVqPVkNbfjjjj@Bj@@@\nfb}A`@@\\sdTTbRLrTrVIrRkNbejjjj@Bj@@@\nfb}P@@@NrRJJIQJKEKB\\hvoQsTEQ@PDT@@@\nfb}P@@@PRfYVVwWppQoIcfZj``@@X@@@\nfb}P@@@XrQQQJJYIKXxXHshaRuS@D@AT@@@\nfb}PB@NZADILkljlktXkS`is@ATu@`E@@@\nfb}Q@@AJ@drmvjroSbmNBgL@PttD@T@@@\nfb}Q@@EI@drmrsvoSdeV\\gL@PDPD@P@@@\nfb}Q@@IZ@eMk]jmjpVcN|gM@DMM@@T@@@\nfb}QB@I^DBHrJIKJSQQIHtYjw`is@AUTC@E@@@\nfb}``@H`YdTTTTTLQbJxkUg^bfBH@bBjH`@@\nfb}a@@H`rQSISJQIIHlxjwhaSP@UUUUEP@@\nfb}pB@ARTBHRYWYUVuhqVg^Sf@BZZP@J@@@\nfb}pB@ARTBPRYWYUVuhqVg^Sf@BZZP@J@@@\nfb}q@@DDtAIf[eUo\\cMNlENVjdJ@@A`@@\nfb}q@@IJtAJ[V{UUU`mF]yNZ@HZZ`@`@@\nfde@``ARAbDq@xd\\QIf]UUUYqVg^``Jh@h@@@\nfdeH@@@XH[IEEEeddhUSdeV\\ADUUUQP@@\nfdeIB@LLxxH@aJWYUWVpPTcViYj@JZ@@@\nfdeP@@@PReenUe|TLGtYj`@B@P@@\nfdeP@@@RRi^VV}z\\txLZjB`@@`@@\nfdeQ@@DXAdTTtTTRbqUXp_QZh@``A@@@\nfdeQ`@DT@hvQQJIPqSYIgMVCEjjif@H@@@\nfdeQ`@DTAsVQQJJJVKII`cVCEjijZ@H@@@\nfde`a@AZlBHYDYEDeDXhiij\\Lywh@bifjJ@@\nfde`b@H\\d@HRfYfWwQJMxLZjfZ``H@@\nfdeaB@OBADILk\\}vmNJLFL@PuP@D@@@\nfdehB@L\\tYp@cIEDeEEDdfKUgAcPQ@EP@`@@\nfdei`@LZlFHIJrQQQQYQEYDhZw`uUUPHEP@@\nfdi@Q@ARUh@Q`HRf^yoYlTjVfZjj@@\nfdi@`@@FrJJIJQIFPly@PUUUP@@\nfdiA`@@HedTtbJRbbV|Ejjj@``@@\nfdiHbBDBuxLPQGicvQQIYJJIHrcH@Jjjf@@@\nfdi``@D@TdrrvokFluUA@UP@@\nfdiqb@JLtRR`RFQQQYIPiIIgJVjVjj`@@\nfdiqb@LR}A@`AFRREQQKKI[UejA`bj@@@\nfdq@P@@HMYIeUU{UZjjijh@@\nfdq@P@@HM[HhddhbidhuUUSUP@@\nfdqPP@LQ@`p\\bbRTTaVRcUUUUS@@@\nfduP@@@PRfYVuU\\LxIt_Afjh@@A@@@\nfduP@@@PRfYemurLyJp_QfjB@@A`@@\nfduP@@@XrQQQJJEIKXYpU`~buS@@@E@@@\nfduQ@@DVAdTTRbLRrRYrVcVCEjejjjjh@@\nfdup@@@PheLr|jjxipSh~CMUP@@B@@@\nfdy@A`@QAHadPzHCDQbLbbbfJbRbwAV}ADD@QP@@@\nfdy@`@@HRYUm]eUZCEjh@H@@@@\nfdy@`@@HrJZJIPiSIBdFKUUUP@H@@\nfdy@`@@QRYWYnuzL|F@Bj`@`@@\nfdy@a@CVADQbLbbRrbbbQuF^C@ATC@P@@@\nfdy@b@@qAHILrrkssIF}A@AAPp@@@\nfdyA@@@IMLrroioNcU@DPA@@@@\nfdyAP@@BUhNQQQQQQDqIdgAZYjZZj`@@\nfdyA`@@DMdTTTtRRQtXKSfiB`XZ`@@\nfdyA`@@HedTtbabRrRMXLVjjfh@@@@\nfdyAa@OAbBPUHYEEDhddYeFRMZB@B@i`@@\nfdyAb@FAb@HrQQJVIKJJDihuRl@ESP@@\nfdyAb@NV\\@HrQQQPsQSHdhFMRsTC@P@@\nfdyH`@DXxHDlbbTTRbTrmNAbtuSUUU@@@\nfdyHb@LBeFB@AFRREQQIIZ[Sk^Vh@@@h@@@\nfdyP`@DR@IIf_mfYHstVjihBH`@@\nfdyP`@DZ@GHhleEEDbmKSoMURE@U@@@\nfdyP`@MNAKIEDheEDZhHvkTDQTMTP@@\nfdyPb@INHDD`dsJnnkLDkTDDUUSB@@\nfdyPb@NJCDDQdTTRatRbTXIwhBH`HF@@@\nfdyQ@@DQAdTRbbRbbuUzpVh@@bB@@@\nfdyQ`@BN@`VQSQIQPqIXs`qZjfdHB@@@\nfdyQ`@DLAHrSKjrmjZpXmUUU@A@@@\nfdyQ`@DZAxJSNlrnkKUgMURDAU@@@\nfdy``@F@PdwLkwKBUxKPPPAA@@@\nfdy```LHx@H}HAdHrFQQIYJIJFPeCVZjh@Jj@@@\nfdy`b@LHx@HrJIKIQIPrDhZsUU@AUP@@\nfdy`b@N^BBHrJJIPzJIQLDzLADPPL@@@\nfdyppBDDDCbkAbR^]dTtRaTRbrXHsfiZjjjh@@\nfdyqb@LFcA@`AFRREQQIJH{WcVVhIBhhP@@\nfgAAB@A@bDfUmmZLEHuh@J@@`@@\nfgAQ@@BT@dsJ}kQ`iFmTuP@T@@@\nfgApB@LLx@HRevUUpPTcViYj@H@@@\nfha@B@@IbURjjjmUUUUT@@\nfha@R@HHehBYddbRRbbSUUULuT@@\nfha@R@HHpPG`eUjjjjuUUUU@@@\nfha@b@HR@fYIHddhdbeUUSUU@@@\nfhaPq@BA`cAg^@DVBLddNbbbjdVjjje`@@\nfha`@@@IJjjjjuUUUUP@@\nfheA@@@YEEDhmEhbpTmFmxKPI@TAB@@@\nfhe`@@@YHhhhdeEbRdcZ|xM@@@@A@@@@\nfhi@B@@QFQQQJGQSILyxLAAU@@@@@\nfhi@B`@QAHadPzHCDYEEELUDeCpUoPQA@DP@@@\nfhi@`@@ArJIKJQPiZcG^`@j`@`@@\nfhi@`@@HR[fyYwU[pVjB@@@@@@\nfhi@`@@^RYvUm]djsfj@@@h@@@\nfhiA`@@H|dwMkZ{IUgMT@@AP@@@\nfhiAc@HHDApHxT\\DjUgf{zBuZj`bHj`@@\nfhiIP@DDhpDc^BdwJ|lzsfkMMMUUU@@@\nfhiIP@DXxHDc^CdTRbfaTVUNZltuUUUL@@\nfhiP`@DZAyIgVYW^VkNZjdHBh@@\nfhiPb@OA``@QddabbRRvRkF\\m@@@E@@@\nfhiPp@DX@IQoARYeYv}YsUfjffjj`@@\nfhiQ@@AJAdTTTTTJtPeVBtuTttt@@\nfhiQ@@INAdTTTbRarrk^BtEP@@@@@\nfhiQB@BJd@HrJJIKPjII`oA``ah@@@@\nfhiQP@DX@IW`yEEDeDdbdsakMUMMUU@@@\nfhiQR@IV`Sa``dLbRbTQbbvXIwjjjjYf`@@\nfhiQR@JLIIW``kDfYu[UV\\MZejfjjh@@\nfhiQR@LVXhw``BLdTVbRbbqqRVjYjjZj`@@\nfhiQ`@EZ@JVQQQQQIEYLTYsT`tDU@@@\nfhiQ`@HHAIvQJJQJIHqBdYsUUUUUP@@\nfhi`P@D@ExJSLrkoqFZluPLDSP@@\nfhi`Q@@pP{p@`pQddabTVrTgIJtmADEKP@@\nfhi``@B@PdwWLjn[s`mU@@@@@@@\nfhi``@N@PdrkLjn[s`mUTEAD@@@\nfhia`@D@HkHheDdhiTiF|EjBA`@@@@\nfhia`@F@ykIEDdUEDejLmzjj``b`@@\nfhiab@LRMx@PeLwN}obfkUU@aET@@@\nfhihB@NJtXHHcHhdlhhiDzLMz@BijjD@@\nfhihbBDBtYw``bGRGlbbRrTTRLhpuh@JijX`PT`\nfhipR@BTYIW``dDfY][UV\\MZejfjjh@@\nfhipS@IZCpSo@bBAA`cHhheHeTiFJlFBAXiVH@@\nfhipS@IZCpSo@bBAA`cHhheHeTiFJlFBAXjVH@@\nfhip`@DDh@VQQJIIFIJBd{sTuL@E@@@\nfhipp@DXH@Rc^BdsJsmzsfkMUMMUU@@@\nfhiqP@DXxBQoArJIQSPjKJ`mVZZjjjf@@\nfhiq`@DDTA`TfYg]m^B|EjiZjjh@@\nfhq@R@HHehBYddbRRbbSRcUUTsUR@@\nfhqAb@AV\\BHRYV}Umhr@Bjij@@@\nfhqQ@@DT@drlsLjXKUAQEP@@@\nfhqQa@EVEHpHc`QdTRVTaTrUFP@RtuP@@\nfhqQa@EVEHpHc`QdTRVTaTreFP@RtuP@@\nfhqXB@J\\dZpPAFRIJIJJqIBeUSMMUpHJ`\nfhq`P@L@PkRSJkLn}oMUPLDH@@\nfhq`R@LPQhp@cIIBhdeDTeVVhHFJP@@\nfhq`c@BN|BlEVJkLbbbfbUbbYJBBbjZ@@@\nfhqa`@J@QqIenvyepVjj`H`@@\nfhqa`@M@DYIfUoWigHBBjjh@@\nfhyA@`A@qBXglPNHdsNjjjsdeV]A@U@E@@@@\nfhyA`@@H}dTRbTTTUQFJ]xKTuUMST@@\nfhyH@@@XxkIEDeDehTjBYspP`Xj@@`@@\nfhyHB@L\\uX@QdbbRbbRRhmFmxMAD@T@H@@\nfhy``@D@TdsLrj\x7FQaoNBuP@`@P@@@\nfhyh@@@XdhvRJJJYQIWTyIQkMSULuUP@@\nfhyh@@@XhrRTsvzrmNBexHPPu@@P@@\nfhyq@@DTtAIemU[TTmF|Ej@B`@H@@\nfle@`@@HrJZQEJIJYJshiZjj@BjH@@\nfle@`@@VRfYU_uTVeFh@@@jk@@@\nfleA@@@ILr{LjnzHTmUT@@@@@@\nfleA@@@YEEMDleHiJJQoM@PUUUD@@\nfleAC@N@bERBYCHhheDlUDefRCFBAbX@H@PU@\nfleAP@@HMEJSKLjt\x7FJ|xKTCAUBq@@@\nfleA`@@HXdvmrnrncQRuS@A@@@@@\nfleA`@@QSdTTTRRTqfdcZ\\t@@AEL@@@\nfleAp@@XDkSdyHhheDeLUmRPXuUUTDAP@@\nfleH@@@XIsIEEDbeMMDjBbfVjfhBA@HPbP\nfleHP@DTuDC@irJIJHqZIJJgOQZijjjYh@@\nfleHb@LBdfB@AFRREQQIIZIZ][ru@@@EP@@@\nfleI`BHDihBEPeGYEMDhULeDdsdTmMUMP@H@@\nfleIb@LBLFBH`BLddJbbRRvVtzwej@@@J`@@\nfleIb@LB\\fBK@BLddJbbRTVVuypUj@@HB`@@\nfleIb@LB}FB@`BLddJbbRRrVtzwej@@@J`@@\nfleP@@@PReevYU]tPijjh@@@@@@\nflePB@K^AbILklrj}Qc`p@TaAE@B@h\nfleP`@DA@eIVUue]^B]yX@J@BT@@@\nfleQ@@DLAdTTTTtRbTI`~RuUATQQ@@@\nfleQ@@DR@dsK[JsioARuUTEDQ@@@\nfleQA@CQDBHUDILrl{ljsdTpDEUUJd@@\nfle``@G@edbbTRbbNv`mV}PPE@QT@@@\nfle`a@BQRBHEDYEEDeMCDdXrSoPPT@PS@@@\nflea@@E@RfVYUV\\tWIjh@`B`@@@\nfleaR@LPQYw`AFRRFJJYIQW\\dXJuPLDqT@@@\nflehPBDX}EHAJCFT}dTTRTRRJRsNFluTtuUKP@@\nflepa@LRmA@@`PQddaTTRRVLfgFlmPREETH@@\nflepb@BZSIPHaIf]V_eV\\bfB@ijjh`@@\nfleqa@EV\\JQ`QG@cHhdliBieDjLbf@BVfjh`@@\nfli@q@HHpxLPAG`cIICdddddiEjjjje`@@\nfli`@@@IJjjjjkUUTuUTBDd\nflmA@@@ILkZvmsQdeV\\EL@EUUUQT@@\nflmP@@@PRfYV{U_CAZ]DJYji@@@F@@@\nflmP@@@PRfYV{U_CAZ]DJYjj@@@F@@@\nflmP@@@PRfYen]|cIF]DJYj```@E@@@\nflmP@@@PrQQQQJZKIDcIV]DJYj``H@F@@@\nflmP@@@TRfYmUo\\eEVMzJYf`h@@J@@@\nflmP@@@XRfYYmU\x7FCAF]DJVjX@`@J@@@\nflma@@@`rJIJYRIIIFcIF\\GtX@Jjjjbh@@\nflma@@I@rJJIQZEJYI`iJlzJV`b@@@H@@@\nflmp@@@PheLr|vj~JBtzHTsUR@@@L@@@\nflmpB@ARTBHRYWYUUVcEZlzJX@Iij@B@@@\nflu@b@HiADILrljm|sfoVC@PUUULT@@\nfluA@@@IRlkjtrRdeFmUA@UUUD@@\nfluAP@@HeEJSKLrorlyiuoMPA@@UP@@@\nfluAP@@\\dZJSLlsr{lxjtTm@qEETq@@@\nfluHR@DTyKShi@DILlk|lzs`iV}@PUUUTL@@\nfluHR@IV|fANJ@DIJrsmrnpSoQSUUUTBAP@@\nfluH`@DXTDAdfYfYu_Yit_AZj`eAbb@@\nfluHb@LLxxLPABTnrjkmxHJQkTlu@EMP@@@\nfluPB@AIADYEEDdUEDhlsc`qS@PL@@@@BGH\nfluPB@AJADYEDeMDeELcSfgQS@DMPT@P@@\nfluPP@DTAsUlbbTTTlTvRXHshiZjZfZjf@@\nfluPP@DVApt\\bbTTtabbJhJs`iZfZiifi@@\nfluPQ@N^YirPA@`cIEDhh\\ihmjByHLZjZiijZ`@@\nfluP`@DX@qIeVyeUrT]DJVjZ@@@H@@@\nfluP`@NQ@EIfVVWV~Je[rV`h@`FH@@@\nfluQ@@DJAdTTTTTRTQsEVbEKU@C@PD@@@\nfluQ@@LT@eLw\\jmvJF|zMLp@@DT@@@\nfluQ@@LT@eLw\\jonJF|xMLp@@AT@@@\nfluQ@@LTAdbbffbbTRNJVbEMLp@PPD@@@\nfluQP@BBAIWhyEEEDhdbddpQgQRtpSAEUD@@\nfluQ`@DQ@KrSLrkssDYiwhmTC@qSQ@@@\nfluQ`@DTAjvQQJQJVJ[ILDYtTmUUSL@P@@@\nflu``@O@HdsNkKN{NF}ELDAUPPA@@@\nflu`a@NZlBHHDILro\\wjpQgQSUUULpA@@@\nflu`b@AQRBHRYVvuvZcARMZ@B`@jX`@@\nflu`b@LQR@HRf[eU\x7F[EEVLzf@@@JZ@@@\nflu`p@E@EJwhyEEEEDdTddqQgQRuHMAEUD@@\nflua`@@`PYIn[VuUYJlzJVj@Bj@@@@@\nfluhP@LF]EHIJudbbbbRTUVQJV\\FMUUU@`U@@@\nfluib@DXXJQhiADYEDeEeCDdeSdeN}@DUUUTL@@\nflup@@@JLeLlj||}AZ]EHP@PTUPP@@\nflup@@@PIdbTVbbbRQNBVCzLuP@P@R@@@\nflup@@@PIdbTVbbbRQfBVCzLuP@P@R@@@\nflup@@@PIdbTVbbbTRNBVbELuP@PPB@@@\nflupB@ARtBHRYV{WfVgEZbf@HZhB@P@@\nflupB@ARtBHRYWYWfVcEZbf@BZhB@P@@\nflupB@ARtBPRYWYWfVcEZbf@BZhB@P@@\nflupB@NBtBHRYeg_fUdeZbfBAbhB@P@@\nfly@R@HHpPG`eUjrjkm`mUPTEE@@@\nflyIPAHDTXBMV]BT\\jW]dTtTVRbfaSAMSUSMUH@@\nflyQP@BI@`uhyEMDeEEDkhkSUUUMLt@@\nflyQ`@D^@jrSOJzmkBdu@AUUT@@@\nflyQ`@D^@jrS[JzmkBdu@AUUT@@@\nflyaP@B@QjJTfvUfW^\\Ejj`XHX@@\nflya`@M@EGIEDdhbdeD^JhHJJjj@@@\nflyaa@KRtX@PkprQSQHiRYIS`j@`ijj`@@\nflypqBE^CAFNBdCpixRyTyIIXiEXdhX{sUST@EP@@\nflyqP@DDTA`XiRYf]vug`fjejjfh@@\nflyqP@DDTAchiRYf]wU[`fjejjih@@\nflyyB@J\\dZwhi@DYHdhdhkDhdJUUMUMUjBFd\nfoA@C@@IdHaDQddabbRRrkF]U@@@@@@@\nfoA@P@@HtYIe[UUhrRfi@Bj`@@\nfoA@`@@HR[e]e^Blyh@H@@@@\nfoA@`@@HR[fV]qZlyhH@@@@@\nfoA@`@@\\RVUoeVBly`BH@@@@\nfoAHB@L\\TX@QdbbRbRRheZMPQDUP`@@\nfoAIB@HXHK`@aJV~vUqRsjYf`@`@@\nfoAPP@DTAsUlbbTTTlTXHsfjfifh@@\nfoAP`@BZ@aInvYWejsfjiB@`@@\nfoAP`@DX@qIeVyWIRsfjZB@`@@\nfoAP`@DXAsHhhddXdbLlyjih@H@@\nfoAPa@LDlx@P`HRYYeUuVLyjjjjj@@\nfoAQPBHF@aVgPeGIMrms\x7FIFlt@PuI@@@\nfoA`P@D@DZrSLrkrQfgMTC@q@@@\nfoA``@D@edTTTbtjQNV\\uUR@D@@@\nfoA``@L@QdTVbbbblmV\\u@A@@@@@\nfoA`a@NRTBH]DYEEDhid\\pQkPD@pD@@@\nfoAh@@@TYIVRJJKIFILxJSTsUUR@@\nfoAh`@JBUYpBILoLkleFmUU@`T@@@\nfoAi@@HDhrPLbbbbjbRxJsjiejjh@@\nfoAq@@DD\\AIeUyWhpufjZ@B`@@\nfoAq@DLTtCDiS\\kN\x7FAV]PD@P@@@@\nfoAqB@JLDPDQdTTTtlRWAV]RuT@B@@@\nfoAqB@JLtPFHdsM|klEYuKTa@P@@\nfoAq`@DXxBSlbbTTtJVhKQffjjjX@@\nfoAqb@NZ]ab`ABTsoL\x7FbakUUULmP@@\nfoQ`@@@ILkZjmFRUYt@@@@@@@@\nfoQa@@N@rQQQQJKGbiVLz`BB@D@@@\ngBQ@@eLUT@@\ngC`@Die@@\ngC`HAbIMT@@\ngC`LADJPt`duP@\ngC`LAbJPt`duP@\ngC`LAbKDvHduP@\ngC``@dfZ@@\ngC``Adej@@\ngCa@@dkH@\ngCa@@dmP@\ngCa@@dsP@\ngCa@@duP@\ngCa@@eMH@\ngCaHH@bNt@@\ngCahHlHRNj@@\ngCb@AEdij@@\ngCd@ADij@@\ngCd@ADkZ@@\ngCd@Adej@@\ngCd@Adez@@\ngCdDI`BHDRZh@\ngCe@h`NdiKLIH\ngCh@@doH@\ngCh@@doP@\ngCh@AGj@@\ngChHD@aIU`@\ngChHL@aIZ`@\ngCh`DLdkP@\ngCi@DDeV@@\ngCi@DDeZ@@\ngCi@DDfZ@@\ngCi@LDej@@\ngFp@DiTt@@@\ngFp`@dfTujX@\ngFp`@dfTujh@\ngFp`ATiTvjh@\ngFq`AbeJfuU@@\ngFr@ACTi[FZd@\ngFt@ADiTt@@@\ngFt@ADi[FZd@\ngFtHE`DILikUP@\ngFu@E`drfmU@@\ngFx@@eJfuU@@\ngFx`DBdwFmU@@\ngFy@DDfTujh@\ngF|HLZ@aJYuif@@\ngGP@DiVj`@\ngGPBAHJPLaYAInih@\ngGPBAHJPtaYCHiCUP@\ngGPHAbIL}U@@\ngGPP@cTfyi`@\ngGP`@TfYj`@\ngGP`@dfUjP@\ngGP`@df]jP@\ngGP`ATeVn`@\ngGP`ATefj`@\ngGP`ATiVj`@\ngGP`Adinj`@\ngGPhMPDIK]U@@\ngGQ@@eMUT@@\ngGQDB@jQBUkUP@\ngGQDL@aAFQRFj`@\ngGQHHGCIHcUP@\ngGQLJIARFdLbdMU@@\ngGQXHlZHROjj@@\ngGQ`@ZdruT@@\ngGQ`@bdvmT@@\ngGQ`@bdwMT@@\ngGQdLZ@j^BTuSP@\ngGQhHl@cIIBmP@\ngGT@ADiVj`@\ngGT@ATeVj`@\ngGT@ATeWjp@\ngGTHE`DILluBR@\ngGTHE`DIL{U@@\ngGU@CPdrm^@@\ngGU@E`drmT@@\ngGU@E`dsmT@@\ngGXHD@aIUVd@\ngGX`DBdsMT@@\ngGX`LDdsmT@@\ngGX`hEIWIMsU@@\ngGY@DDfUjP@\ngGY@DDfYj`@\ngGY@HDfVj`@\ngGYHI`xISkTbB@\ngGp`@TizHvjj@@\ngJPB@fRHTQhbOj`@\ngJPBAbJHt`YAIjj@@\ngJPD@eRHczh@\ngJPD@fRHczh@\ngJPH@DIJuP@\ngJPH@EIRuP@\ngJPH@eQ}T@@\ngJPHADILth@\ngJPLADJPt`duU@@\ngJPLAHJPt`duU@@\ngJPXHlPDQzt@@\ngJP`@TeVh@\ngJP`@TeZh@\ngJP`@TfZh@\ngJP`@dfvd@\ngJPhLQVIKTp@\ngJPlLPDRPTaGth@\ngJQ@@djsBJ@\ngJQ@@dkU@@\ngJQ@@dlu@@\ngJQ@@dru@@\ngJQ@@dsU@@\ngJQ@@duU@@\ngJQ`@bdvu@@\ngJT@@TeZh@\ngJT@@Te^h@\ngJTHE`DILmP@\ngJT`E`TfVh@\ngJU@H`dlu@@\ngJX@@dkU@@\ngJX@@dkt`@\ngJX@@dku@@\ngJX@@dms@@\ngJXDBLQQBS]V@@\ngJXDBLQXbS]V@@\ngJX`DBdru@@\ngJY@BDeZh@\ngJY@BDfZh@\ngJY@DDefh@\ngJY@DDfVh@\ngJY@DDfvd@\ngJYHC`DIKTp@\ngKP`@Tixjj@@\ngKP`Adi\\Zj@@\ngKQ@@eKcRp@\ngKQ@@eKcUP@\ngKQ@@eOEUP@\ngKT@ADi\\Yi@@\ngKT@Adi\\Vf@@\ngKX@@eKcUP@\ngKY@HDi\\ZV@@\ngNpJAbJHLaYArBS]UU@@\ngNpP`jusyImfi`@\ngNpXHlPDYIHTmT@@\ngNp`@dfUZf@@\ngNp`@dfuZj@@\ngNp`ATf^Zf@@\ngNp`ATiUjj@@\ngNplJqHJPtadTaeTp@\ngNq@@dskUP@\ngNq@@dssUP@\ngNq@@eLuUP@\ngNq@@eM]UP@\ngNqTHmVDPhbU_MS@@\ngNq`@bdvkUP@\ngNq`AVeJmUP@\ngNqhHl@cIICej`@\ngNt@@Te[zj@@\ngNtDLpDDHRevnl@\ngNu@HPdjkUH@\ngNu@H`dlkUP@\ngNu@H`dl{UP@\ngNxHD@aIUUfaM@\ngNx`LDdskUP@\ngNx`LDdssUP@\ngNx`LFdjmUP@\ngNxhMV@aI[ji`@\ngNy@FDfWj[@@\ngNy@LDeVjj@@\ngNyhGE`DYIITmT@@\ngN|@ADeJkUP``\ngN|`HfJdlsTp@\ngOp@DjWkB@@@\ngOpH@DILkW@@@@\ngOpHADILkW@@@@\ngOp`@tigujj`@\ngOp`AdeekZZP@\ngOq@@drm[UT@@\ngOq@@drm\\@@@@\ngOq@@drm]UT@@\ngOqHL@aIYZzfd@\ngOq`@fdrikTl@@\ngOq`@fdrikUL@@\ngOq`AVeL~mUT@@\ngOqhHl@cIIBjujh@\ngOrHEcPDILl[MT`@\ngOt@AdigkB@@@\ngOtHE`DILjZuU@@\ngOtHE`DILl[MT`@\ngOt`E`tfUMZi`@\ngOx@@drm\\@@@@\ngOx@@drm]UT@@\ngOxHHHaIeZzjh@\ngOx`BDdwM[UT@@\ngOx`DFdrikTl@@\ngOx`FDdwM[UT@@\ngOy@DDfYKZj`@\ngOy@JDiWMjj`@',Oq,Pq=false,Qq,Rq,Sq='eOHBNZ`pge@\neFDBcA@\neFhHbayP\ngFy@DDfXujhCAF|f@\nfHcpAa@\neMPBchLDnR\neMBBHRZCAKd`\neFJHbHpXI@\ngCa@@dkHFBbyL\neMJDbDf``\neMHAIdLF^P\ngJXHD@aIUj@pHVOI@\ngGY@BDeVj`LJHm^P`\neMPBcTH\ngCd@ADkj@pTWDL\ndeTH`ABHRVUunh@J@B\ndidL@@SdfU{aZjj@H\ndiDD@@QIeuZjh@pkJ[ayA@\ngOy@HDfUk`@@H\ngJ\\@ADe\\u@XKGd`\ndeV@@@Rffeijjj`CClJfxYy@@\ndeVH`NdLQIe]ZZ@Bd@LFALJny@`\ngOx@@drm\\@@A`Qb~IT\ngF|@AbeJfuU@P\ndmN@@@rJJIQEneX@@@@B\ndifH@BAIfuxV`@@C@j[ayF@\ngOx@@drm\\@@A`PZ~IX\ndaF@@@Rfu[j@@@LJSBinHG@\ndeU@@@yJUuRh@Jh@LNaLIagd\\\ngNxLL@aAABDfVZj@`\ngGPB@DHHpQPaIUZdB\ngGY@BDf^j`LBHmqF`\ndazH@LAIUUjj@LJPj[nIF@\neMBCDRZCAKd`\ngGQHJLQIUjdB\ndefD@B@TfUvZjf@H\neFA@HoBJFD\\h\ndayD@DpNRYWZj`C@lInxbT\ngFt@AdiTvjhC@qF|Pp\ngOtHLPDYHhckSM@XXI|a@\ngNqBLIAREdGHIMmUTA@\neMXIDfP`\ngJT@ADiVhPT\ngGY@JDeVj`LJHl^R`\ngCi@DDeZ@pTwH`\ngNx@@eJmUPFCDVKyJ\ngNxDLHaqBRjuU@P\ngGP@DjVj`LKEc^P@\ngNx`LFdjmUPFCDQkyL\nfHbpAa@\ndidHPBBHFHRYgVzB@`@`\ndaD@P@bNbDfUzZ@B@B\ndaD@`@bDfUjZ@B@CB`SJ{dL\ngCe@E`dsPFBV@\ngOy@FDiguie`LMc^IL\ndaDh@DInAIf]nZiX@`\ndigH`LH^@HRf_ljfYh@ppBkbHH\ndaD@`@QDeeVz`@@B\ngCahHlHRNtA@\ngOxhMDOAJmZvjhC@xu|`@\ndaED`LJDCpRjw[fjj@LLinxfD\ndcMDPLzfBHKpRUZvYvjZjh@`\ndeT@`@bLbdTJPsU@@@D\ndcNH@DAIee^eVhHB@CCbine]yB@\ngJX`DBdru@XKGdP\ndif``DqnDHaIeYih@J@B\ngOx@@drm\\@@A`Pm^Hl\ngCa`@ldsPFBV@\ngNqhHl@cIICej`LD[rL\ngO|HLV@aJY}ZYhB\ngF}@LZDigVfXB\ndaE@@@aJyUnh@@@pXLJf{dP\ndmN@@@RYVuiiV@@@@@`\ngFx@@eJfuT`XEXwbD@\ngFt@ADiWFZdCCPwdP\ngGPhMQDIJmU@P\ngNqHLHaIYzj`H\ngNx`JBdr}UPFCDVKyJ\ndaF@`LBHRVU[j@@@H\ngGX`HJeZuTA`Uc^R@\nfHfpA@\ngOy@FDfUkZj`LBl[rL\ngGQ@@dluTA`Xm^HD\ndaz@@@Revjj`CAdpj[nHD@\ndiTD@@yJ]VaNFjjj@H\ndmLL@@iTie]|kahHBh@LJPfE^QH\ngOu@DPdrykURA`l~Q@\ndiT`@@rnRfUjEnBA``B\nfHgHA@\ndid@@DjU^nBBD@LFaLiaxa\\\neFACr`bLD\ndcl`@@`nReeWZY]@``@@B\ndid@`@bDfYUn`HD@LFCBiny@`\ndmv@`EBHrJJIHin`HFpHy`\ngOq@@eL~mLlAal[qI`\ndeUH@DpDf^UFVh@@@pYLinGdR\ndeVH@IAJYW~F``H@LFSBinyD`\ndeV`@@rfyJWUlIifjBJh\ndig@@@aDiyWaZ@@B@h\ndeWH@DJPRY[TYZ`@@B\neMBBHRZCAGe@\ngNy@FDfUZk@pRwIP\ndiFL@LANRYufjf@H\ngGU@E`dkMTAad[rP\ndayH`Dr@|Djyfjh@`\ngNtDLpDHHR[UjhB\ngJXDL@aABS[TA``nKD\nfHe`Aa@\ngJQ@@drt`XI[dH\ndid@@DjU^nBBD@LNaLJfGd\\\ndifH@FAIYU[fZj`CBlJnGdP\ngJ]@HbDfZdB\ndie@@@EJf[XZjj`CChSBiny@`\ndk^@@@RYWYVftx@H@@@H\neFA@H`bLFE\\`\ndeTH@@RYWZf`@j@CC`SBhYyG@\ngOu@JpeJymUTA@\ndcN`@@pjYJYenk`Hbj@B\ndmu@@@QIUYVUZh@@@pdDpfxYWdP\ndmvH@NAIfUTUZBB@@p[BiaWdR\ndkmH@NtDfYU_EZ``b`@`\ndkmH@JTDeVVutvjh@@@`\ndieH@DHLbbbaIjZjPB\nfHaXA@\nfHbHAa@\nfHdXA@\nfoA@R@HHpx@QddebRfRpiFlm@@@P@A@\nf`i@@@LdbbRRlRSI\\D{qAHT@D@@P\ngNq`@VdsMUPFCDQkyL\nflmA@@@YEEHdcLdddSdmV]EKUUSU@E@@D\nffsA`@@YkdTTbRLrTrVRIrRkNbejjjj@Bfh@B\ndaE@@@yIe^f`@`@pILJnHG@\nfeg@`@@SrJIJYYQQJQZqYTYzH\\uUKMUUUUUPDLiVEPPQbgARtYw`qSdmV}gr@p\ndcMH@HdDf[W[ai@Bh@H\nfc\x7F@H@@JbdoU{RYfYU_uyYubc^\\D\\vi@H@BjjjhH@H\nfb}A@@@YEDeMDhTihjLUiwhyZZjjjjjh@LKHqPTeZMYpTyxcdP\nfnk@P@@FbuIefYZWz^`mA|gMVfjjjjjjj@B\ndmLH@@rJJIQEneX@@@@CA`iae^Ig@\nfleA@@@YDddhiiDdYipTtsUUUSTABHe@\ndk|@@LdbbbRQKauS`@`@`@H\nfjc@`@@ERYYefUez\\UyvPiZZjjjjjj@CCrLTXJVcV]zJ\\ey@H\nfjc@`@@JrJIJZJJUJZJcMF|EIVfjjjjjj`@`\nffsA`@@VkdTRTtRabrTRhqQoARUijjijjjj@B\ndie@@@aJyW[f@B@CBlJfGdP\ndeT`@@`YRnUunX@I@B\ndmL@@DjYeVdU@@`@@H\nfgA`@@@YIEDeDUfJxHujjhHB`@LNPQbgIZxQbcV@\nf`ia@@E@Rfume]gEF]z`@jjhh@H\nfoQ`@@@YIEDeDTdpWEF]UUAAE@@P\ndcnH@EAJ[UU[ev`@jb@LFPje]xfB\ndcnH@MAJ[UU[ev`@jb@LFPjE]xfR\nfjc@`@@ErJIJZIPiSQJcEZ]zNVfjjjjjj`@p\\cENBdkQkNBgODQJ\nfb}A@@@ILk[Kk\\tyKSoQsUSUMUUUP@XN`cEFBdhugASgb@Q@\nfb}@@@LdbbTbrNTtWEBt{t\\sMUUUUUT@D\ndg\\L@@{\\dTrTRtKCPAET`A@\ngOtHLPDISKWSU@XDJodH\ndknH@ECHhheLcFz@`jh@H\ngJP@DknhCCSGdP\ndcLh@LKaCIELeBhwT@ET@D\ndmv@@@RfUYyZBB`@LAALInF^IR`\nfgA@@@DjYU_VByHu`@@@@@@LIXISdeZM[dD@\ndidD@@yIfVXXBH@C@`[axfT\ngCa@@dmHFBVxa@\ndaxL`HS@\\DjeZjh@`\ngNphH`xITkUTA@\ndaGH@Dp`RYeifjf@H\ndg\x7FH`LKU@HrQQSIJHyTgSA@qP@P\ndev@@@RfU\x7FJxZA@h@LFAB[iyE@\ngOqHJ@aJUcZZhB\ndaF@@@Rfu[j@@@LFALJfyG@\ndieH@NDDfUvf`@d@H\ngFx`LDdrfmU@XHwdp\ngNuhJxl@cIIJej`H\ndifD@JAdiWTjjjj@H\ngGT`DQTfuiaMXIWdp\nffcAp@@DEiNZyEEEMDihbdiVBMzIViBb@@hh@@pTQVcV]ENZ`\ndo^H@DAIeV}YZfifdLIFU`\ndg\x7FH@LiPRUYYWESnj@@@@@`\nfoAQ`@DD@irSLwN~QakMUPHE@@XLpJVcNx`g@\ndg}D@DHERYfywIcmjjA@`@p[AgSobEp\ndieD@DHNRY[Rijih@pSFybDH\ndg}D`Ad{BHRYVuUiev@BffH@pXDrmob\\H\nfoAPb@ARlxDPdrmkKtXht@EMSD@D\ndclD@@[HhheBeSKkTp@P@P\ndcNH@MAJ[WY[j@Bg@B\nf`qA@@@ISZzljscoT@Dp@@B@kARLDhsQkN}q@i@\nfbmpB@NV|@HRfYe]WnirQmN@@@@bj`@@pt@cAJLxLTy@\ndg}L@LhDWHheDbhdSG[UTtuHQu@\ndmLH@@rJJIQEneX@@@@CC`rfFUyG@\nfbmPb@FRTX@QdbdTjRtrbRps`dyjjjYh@J@@pLAAJ\\e[pXisrF@\ngJ]@EbDeVdCBKba@\ndcl@@DjYU_egX@@@@@pELJae]xc\\\ngJP@Di[xCBGdp\nffsA`@@IkdbfbTTTTRtwxkW`ecfjjBjHjJb@B\ngJYhMDPDYIBm@P\ndmM@@@iJYWWJxYB@j@CBbiiWdH\ndaD@@DjYvxH`@CCLJnXc\\\ndg^@@@rQQQIIRcmAPTT@D\ndiV@`J@HRfU|kahDB@CB`PfGd\\\ndmL``LXU@HRfvUxYV`HFH@pHDj^Ph\ngC`DADJHRZhCCBKdp\ndmtH@@RYWYih@Jh@LAALJnE^QH\ndmtH@@RYWUih@Jh@LNALJiWb\\H\ndeUD@DhFrJIJHusUMT@XDf]OHD\ndg}@@@aJVYU^Svv`@@@@`K@Z[ae]ObEX\ndmwH@DePRUe]xV`@d@H\ndcO@@@b\\dtTRJgBpDEP@P\ndaD@@DjYvxH`@C@linxbD\ngJX`BDdru@XS\\a@\ndie@@@aJyW[f@B@C@lJfxdB\ndid@`@bDfYUn`HH@LFCBiayG@\nfoAPB@NBADILkjrmFV]@AL@@@ar`\ndeUH@BxLbbbRKCA@d@D\ndk}@@@iJUmUR[atpBJ@@@CChPiWSyE@\ndeUD`LjD@HrREQICJt@@@XUc\\L|`@\ndmuD`IVD@HrRFIKKaV@BP@`\ndmtD`ATIAIe]Vf`@jP@pILJnIt`\ngGP`@deUjpLI^JX\ndmvD@LADfvUYUjjj`B\ndcND`MEPdDfUuZf`@jX@H\nfjc@`@@HRUe[Vun^Bt[pTeZ@Bjjjhj@CBqJ\\dkQkN}GI^H`d\nfnk@`@@HRYVum[VydjshiJuj`@jjjjJ`@`\nfa{@P@@HkwHdhdididihdUpVc^BdkP@UUUUTeP@P\nfjcAB@C@bLbbbRdRbLJbXJUoQRVejjfjjjh@LG@pQ`iJuYpXyKrBP\nfb}A@@@ILk[Kk\\tyKSoQruSUUUUUP@P\nfjcP@@@FrQQQQII[FIJgIVt{IP@@@``@@@@`\ndaD@@DjUZxHH@CAL[ePaUp\ndcNH`BpHaI[fuWZZ`@@B\ndcND@DCdefV]]Z`b@@H\ngJQDD@b^BRkTA@\ndg|@P@bEbLbbbTRRKR]pP@P@@D\nfoA@A@@QChaIe]YWhpuh@J@B@@LE@Q`eIJtYwH\\@\nfoA@@@DjYUgWfiF`XJbh`@`\ndeTL`HS@BLddlRPrm@@@FETwCOHp\ngOx@@eLmXD@A`Uc^Q`\ndclD@@EJY}erevfVfiBJX\ndcNH@DAIge]FVh@I@B\nfoA`@@@ILkjrmFV]@AL@@@arptcAJ\\DkQkNy@@\nf`qqa@FF}Sf`RBICHiMEMCEmPUgUUUSKT@D\ndg}@@@aJVYU^Svv`@@@@`K@Rine]obBX\ndg}@@@aJyeg_eNvB@@@@B\ndmND@LAdige\\jUifZi@H\ngKP@H~Jj`LEcqQ@\ndidD@@QIenifjj`B\ndklD@@QIVVU}MjlH@@H\nf`q@`@@HRUfYU_Sk^Zh@@@@@H\ndeVD@BADfVuFVh@@@`\ndmwH@HePRYm]xZ`@d@H\nf`qP@@@^RYWUe^cKN`@f@B@DNT\nffsA`@@VkdTRTtRabrTRhqQoARUiZjijjjj@CArLTxJRcN}GIZ|`d\ngNy`DBtf]Zi@`\ndeV@pJBPRPrPrJPqICMT@@@^DwCH\nfaw@@@LdbbbbbbeeRRwMB]Jq`y[Wh@@@@@@@@@@A@\ndmN@@@rQQJJFfEP@@`@B\ndmN@@@RfV_kad@@B@@`\ndco@@@g\\dTRbtM\\J`@AB@A@\nfdep@@@PheL}sJkxiXOhsUP@@@`@D\nfhiP@@@ARVUue]Yrsn@B`@l@@`\ndcn@@@RfV]zXU@@BD@CCdJnFUyF`\nfbmP@@@IRfV]}e]hpR`q@@BF@B`@PUP\ndeWH@DJ`RYYTfZij`B\nfb}@`@@YRYVum[ehrRkNBejfjjjjj`@`\ndk~H@NAJUmmRkatzVd@J@B\ndklD@@IIf]vVzB@jh@LICBkagSyB`\nfbuA@@@ISZ{Ll{lxgM@AMTD@@HBh\ndaF@@@RYe[hB@@LBpj[bAp\ndco@@@bdigyWJWYj`@P@pYLKawbRh\ndaGD@Bi`QInUnZZX@`\ndkmD@DpJRY{e]]Zf@B@B\ndaDH@@RZW[jii@H\ndazL`FaL@HrRRiKUPA`VEW\\a`\ndaG@@@[diWRh@I@C@hSByG@\ndigH`LJn@HRf_ljfYh@ppBGdH\ndmL@`@bLbbbTQ[iV@@@@@`\ndeVD@DCdfV[hVjfh@piJXYyF@\ndaDH@@RYe[hB@@LJCBinQp\nfoAPR@JLIIS`UbSLzmnsakTmTuUP@X\\Ab`mV]rA`\ndcmH@NDLbbbTVKCJm@dE@A@\ndg~H@FAIfU{TYNz``@@@B\ndklF@@sittieY_hZjjjh@pfBinF]ObHh\ndcLH@@RfUWnf``R`@pyBinW^Id`\ndo|H@@rQQIIZHkihHDjd@H\ndo~L`MaL@HrRPjIIKISVhHHh@H\ndg]LPElYLHc`cHhdliBiSP@Rtp@XLBXUKrK@\ndcNH`LdHaIUe[EZh@b@B\ndk|H@@RfYU_JGUN`@@B@@pxfzUt~PH\ndklH`I@HRYeuin`HJj@B\nfdyQR@FJtZwh`BLdTbbJTRlvF\\eUUULsU@AaHBLDiqQ`|aD\ndmu@pNTIAICICHiCDeafj@B@CAJ[bBh\ndkmHpNTpdDdLdLbdLRRrFZh@JP@`\ndkoL`LhXPpBLdaTRbQeUhBBP@`\nfdy`a@BRlBHcDYEEDehXhlfB|FBBFB@`@B\nfb}A@@@YEDeMDhTihjLUiwhyZZjjjjjd@H\ndco@@@rldTRTJU\\naA@@@A`Prn|qLX\nf`qA@@@IMKMjoiuoUV@@@@@D\nflm`@@@YHhhhhdhb]RRkFCzH@@@@@P@@D\ndg\\L`AWPbDfUv{ZZ@Bij@B\ndg|H`KBHRYWYWimN@B@B@@`\ndaF@@@RYVih@H@LBSJ{b@p\ndeTD@@QIVYQehB@@LBinGbQH\ndeTD`NDHaIfVVfBA`@LJCJX^QH\neMhDRZCAKd`\ndeU@@@cIMDeBwL@E@A@\ndaG@`D[`bDfUZZ@B@CA@sbEp\ndcMDPDy]BHXHRYV}jz@Hi`@`\ndmt@@DjYZVTHbh@CBdxYWfES@\ndev@@@Re[TjFP@@@@LFaFxYyA@\ngOpHADILkW@@@XHWaVK@\ndmv@HBBHFPfPVPRYUzih@Jh@H\nfhy@`@@\\RYee~uYrVoAZA``@H@@pI`uo\\YPTmNB\ndk\\L@@{\\bbbTQr[iUhF@j@CALOaTpr\nfoA@P@@\\ekHhheECEf\\TYhFHbh@B\nfoA@P@@\\d[HhheELiVRUYhFHhh@CCFMNyabk@\nfoA@P@@\\T[HhheDXef\\eihFBJh@B\nfj}@P@@\\tEIfVWnU]YrQdeZA`hB`b@@H\nfj}@P@@\\ueIfVWn]}YrQhiZA`hB@j@@H\nfj}@P@@\\MeIfVWVu}YrVhiZA`bB@j@@H\ndmL@@DjYUVGi@@@`@LFrnFUxbD\ndclH@@RfumVy]h@@@@CCdpfE]yD`\ndmMH@HxLbbbTQ[iV@@@@@`\ndeVH@HCHiEBdLuP@@A@\nfoA@b@HH@DYIHXeEhbkUgKPP@@@@P\ngJU@LPdjt`XQ\\TH\ndeUH@JDDinUzZBA@@`\ndclD@@IJ]YWaIVj`@h@H\ngNpXHlQxIUZuTAaGqF`\ndiDB`HSBCpRjuVji@H\nf`a`R@LPP{p@cIIBhddbhhruAADp@F@jLxJRm^yC@\ndk]@`FD@aJY}e\\kSif`@`@`\nfoA@B@@QBRsLjnjqgSP@@@@@P\ndmtH`IBHRUe]xY`@hBH\\F@j[exdH\nfhiAB@O@bDfYUg_QZlzdH@@J@DDL\nfhi@`@@HRYyWVUEKpVh@J@@@@pebcNBMYwnHTD\nfjm@P@@HceIge]YWVEKpVh@J@Bf@@H\ndmt@@DkYV~Gh@J@B\nfjm@h@@XDkSoQRWIEEDhdicDXZRcfjjjjjjj@CBLDDhs`iJtYw`irWdK@\nfleab@OPQD@QddabbRRvrbkF\\m@@@EL@@XJhs`iJt[r^P\\\nfbuQb@OA`aH@cIICEDdelehjqgKP@@ATp@A@\ndkld`LWSa@BLddlRVRFUh@JP@`\ndg}H`AfpbDfUmYZYS`@ijR@LJALh~Iv`\neFJHSHp^I@\ngOpH@fILkW@@@P\ngJPhI`fIKTpD\ndaE@@@aJyUnX@@@ppj{fQc@\ndmvDPLa@BNdLddlTVeUhH@@LBqagdJ\ndcnH@LCIEDUDeeKmLuUT@P\ndaFD@LADfyVyjj`C@binyB@\ndcnD@BCTfYVuiEX@@B@@`\ngC`LADJPt`duPFDGI@\neMACDXaIhH\neMACD\\QIhH\neMBBlRZCAKd`\neFJHAHhP\neM`AIdLF^P\nfHbXAa@\neMBCDRZCAGe@\ngJQdEbOBRD_M@P\ngChHHDQIj`H',Tq='fHbTA@\nfH`pA@\ngFp`@dfTujXCAZ|a@\ngFx@@eJftu@XVKF|`@\neO`BNZ``\nfH`XA@\nfHdpAa@\ngNxHLHaIVjj`H\neFJHbHpP\neMABHXaIhH\ngJXHD@aIYj@ppqyH\ngCi@HDej@pRwDH\ngCi@LDeZ@pTWI`\ngCd@ADiZDE@\ngOx@@drm\\@@A`plZp\ngGX`LDdsmTA`m^P`\ngCiHLaDIMLA@\nfHapA@\ndeTH@@RY[TYjp@@B\ngCa@@dkHFBVyH\ndeTD@@eIYWVy`@h@LFpjXYyD@\ndeT@`@qDeVUFZX@@HR`\ngJQhHl@cIHUhCBGd@\ndifH@DAInUxV`@@CBdinGdD\ndeV@@@Rge[aj@B@CChPjxYy@@\ndeV@@@RgfTYj@`@CChRfxYy@@\ndeVD@AADfVuFVijh@phj[iy@`\ndeVD@ICDiieZZjjh@`\ndaGH@DK`R[e[fiZ@LLQnyE@\ngOr@Ajti]qZY`H\ngGPhIPDIU{T`XXK\\a@\ndid@`@qDeYWaf@@BH\\NABinGdP\nfHa@A@\ngNq`@jdvkSPf\\Ll~P`\ndedB`LkiCDRV{njjh@`\ndiGH@Dr`RY{fjj@H\ngJY@BDfZhC@bK\\a@\ngJY@BDfVhCAK\\a@\ngGY@JDf]j`LLl^R`\ngGY@BDfUj`LLm^P`\ngJPH@DISUPFABqyH\ngJX@@dlu@XZX|PP\ngNxHF@aJZzjPH\ndazD@LADeUffhHr`\ngGT`EaTf]jPLDmrD\ngCh`LDdsPFDWI`\ngGX`JDdsmTA`l^R`\ndmv@@@Rf~UeZj@@@LEBDpfxYT\ngOx@@drm]UTAaqEcV\ngOx@@drm\\@@A`Qc^IL\ndmvL`BaL@HrRRjIJUVjjh@`\neFA@HoBJD\ndiFB@BAFEInuZjd@pILJnQp\ndayH@DpDf]Vjh@pKBinHg@\ngNuHLzHaIUji`H\ngNt`E`tf]Zj@pJM_I`\neMJD|Df`pYy@\ngJPhLQxIRuPD\ndaDL@@KdfYvyjV`CCLJnPp\neMBBlRZCAGe@\ngOq`AVeL~mUTA`Yb~Q@\neMPBchLD^T\ndaF@@@ReYJjjj@LNaLJf{d@\ndaE@@@aJyUnX@@@`\ngCe@E`dkPFBbyL\ngCahHlGBNtA@\ngC`@Die@ptVy@\ngC`DAb[DRVhB\ngCaHLLQIZ`LDEqS@\ngGPBADZPLaYAIZjhB\neMABHYAIhH\ngJX@@dkT`XFKGd@\ngJY@DDeZhCCSGbB@\ngGT`CPdfuj`LLl^R`\ngGX`DJdsmRA`enP`\ngFq@@eMqUW@P\ndkNF@BAIWSR[YVYjjfX@`\ndeVDPL[`bB|DeYgFZjZh@`\ndeVL`LxY@HRf][JjZV`cJ\ngO|HEfHaIeZx@@B\ndaxD@@QInuij@LBRf{dD\ndaxD@@QImUijBLlBRf{bXP\ndedB@@PYR[fyijXHqpQIxe\\\nfH`TA@\ndaxL`HS@BLddNRuR@P\neFJHqHpP\ndaxL`Lk`qDenzjh@pXDpj{bPp\neFPBca@\ngG]@EjDfUj`LEcqJ`\ndedd@DpaCdfU{ZjZ@H\ndmOH`LJQ@HRf^yriVfZZh@`\ndaE@@@yIe^f`@`@piLJny@@\ndevH`LX@aJWY\\HYiZZd@`\ndaEH@DHDfURijZ`CCL[nP`\ndaFH@HAIYUnh@@@pXHpj[d\\\ngFt`CQdidviXB\ngJPXHlPDQztAxlP\ngJPDAbGDRUj`H\ngNx@@eJmThFCbqky@\neMA@HPaIXLD^T\ngGYHLQDIJuU@P\ngGP`ATf]j`LLl^R`\ngFp`AdeoEjhCCHwd`\ngOp`AdeekZZPLMB~R`\ngF|@@ZeJxru@XYF|d@\ngOy@FDiekjj`LKEc^Q`\ngOx@@eJqmUTA`xlZ~P@\ngOx@@eLvmUTA`xlZ~P@\ngOxHBHaIeZzjhB\ngKP@Di\\YZ@phbq@\ndiDB`HSB@HrRPyIZj`CCBknHp`\ngNq@@djmUPFEfM_DD\ndcLL`HS@BLddjbRtjmP@P@P\ngJPJAHR`Tai@rBSUTA@\neMCArhabHzCCI@\ngNy@BDf[Zj@pruxbp\ngJY@DDfzhCCSGbB@\ngNx`JDdskUPFDwLZp\ndmVD@JADf^Uvjjh@`\ngChHL@aIVPH\ngNy`LETeUZZDs@\ngNt`LPdfUZi@pexlp\ndiEH@DpDfYUjj`C@bkaxfL\ndidD@@EJ[W[j@B@CBdJfGdX\ndmtH@@rJIJFRf`@j`@pxDrae^Pp\ndaG@@@kdiVrX@a@B\neMhDRUB\ngOx`FDdrikUTA@\ngJXHD@aIUZ@`\ndcnL`LaA@HrRPjIKTrzmPHD@FEYtkh\ngCi@LDeZ@pTwH`\ngFq@@eNqUU@XZX|Rh\ngKP@Di\\Vi@pLVOH@\ndiVH@BAIfUInFjZi@H\ngNqhHl@cIHUEj`LLZ~P@\ndaxD@@QIe]ji@LBpf{dT\ndiFL@J@aRY]Zjj@LABDpfx^QP\ngNx@@eLmUPFEfM_DD\ngNy@LDeUji@phQkxi`\ngGX`BDdjmTA`m^JD\ndazD@FCdfUVjx@`\ngCd@AH}PFBVyH\ngChHHGBOTA@\ngC`DADZHRVhB\ndeVB`BaLd@cIIBeDwKULpA@\ngJT@@deVhCCCGbb@\ngNu@E`drkUPFFM_Dl\ngGXHJGAJijhC@qX|e@\ndidH@@RYeVz@``@pXLJf{dB\ndaEH`Dq`BDfUjyjfPC@`SNyE@\ndieH`Dq`BDfUfnZii@H\ngJPDADFHRYfaHp\ndmuL@DpFEIeY~nZifh@pILe^Qp\ndklL@@STfue]eVj@B`@pyL[ad~PP\ndifH@HAIVUxY`@@bGA`Pjx^Pp\ndid@p@bBbFbDfYoa`b@@H\ndid`@@pjRfUjXBB`@pSaxbL\ndcNH@DCHheEDbnmPT@@F@hUMproHt\ndaz`@@SFyIeYjf@LJAL[nQP\ndaG@`LK`BDimVz`@@B\ndeT`@@pjrQQIFTpDEP@P\ndid`@@pjrQQIFf@`h@LLKaxbL\ngNy@JDeUjj@phVKxiP\ndigD@DpP[HhdhZyjfd@pqLkdB\ndmuH@DpDfWeYUj`@@CAlInF^Hb`\nfH`PAa@\nfHdHA@\nfHchA@\neFHBJFE@\ndmLH@@RYVuiiV@@@@@phJxYybXh\ndid@@DjY}nBHH@LJSBh^Qp\nf`i@@@LdbbbTRVHeZ][uHAD@D@@XRhs`iJtZwLX]x\nfhyA@@@ILklrstXYw`p@TaDA@@P\nfdeA@@@ISZvmvkNJLFM@AUH@D@A@\ngOu@HpeK^MKTAaKqY`\nffsA`@@LudTTTeRdVTtLIps`ySeijjjZjjj@B\ndklF@@XUttief_kjjjjh@`\nfi{@h@@LDipTmzOHhhhhhbXdiidnBu[IV`fHbjjjh`@pmIVc^CENRmzObNU@\ndayL@DpFyIgYjf@LLPnyF@\nfb}A@@@ISZvmk\\lxkSoQsP@UUUUEP@XNQ`cARUhugASgdH`\nffs``@L@QdTVbdRQfRfRxhu`qKUjjjjj@Bh@B\nfnkA`@@HkdRTTbRLrTrVSIZMxNRmUUUU@AUP@D\nfa{A@@@YEDeMDhTihdXjLUiwhyZZjjjjjjj@B\ndaxH@@RYWjZPcKA`SBinIG@\ndet@@DjYUX^d@@@@CAlJnF^Hc@\ngF|LHjOC^A|DiTt@@B\nfnk@`@@UrJIJZIPiSQHrcEZ]zNVfjjjjjZj@CCrLTxJRmFlxJ\\myBH\nfc\x7F@P@@E{OHhdiheBeMDhUEKQbmN}GKSUUUUUTsUT@D\nfnk@@@LdbTTRbTRQNfiKak^cejjfjjjjj`@`\ndmtL@@QTfyeQehBA@C@jXYxb\\\ndklL@@PtfVV]WVhH`P@`\ndeTD@@eJ[WVz`@h@LBRfgfXP`\ndcLL@@{TimY]ah@bh@LJPaW^Id`\ndcLH@@RYeZvz@`j`@pxBXYW^Q`\nfoAPB@LD@DYHhdcEEEQagTuPDQ@@P\ndeVH@DAIgeQej@@@LJSJX^It`\ndeVH@IAJYW~F``H@LJPj[nId`\ndk^@@@RfYU\\]Tzjjjj@LKBDpj[ae]L\ndid@@DjYmaBH`@LBrnGbDp\ngJX@@dmu@XKGdP\ndg^L@D@[rJIJJIGZ[UAPD@D\ndeVh@LKadDimY[j@B`@paNHL`\neMPBcXLIyP\ndcMH@ITDee]UnX@Jh@H\ndg]D@DpCRYVuveVj@BX@H\nfgA@@@DjYU_VByHu`@@@@@@H\ndclL@@{TivY~DeZhHB`@`\ndid@@LdbRQk``b@CAhPjX^Q`\neMJH\\Df`pgd@\ndeTD@@YIfUqehH@@LFJfxYyF@\ndclH@@rQQRJJuJ{PUDB@FCES\\L|Q]@\ndif@@@RYWZZ@BP@pXDpj{dB\ndaF`@@pjYJYfn@b@@pPfyG@\ndeV@@@RYV~f`@i@B\ndeV@`BBHRYg]n``I@B\ndmNH@NCHhheDVzU`@@@@LJCJe^Hp`\ndcnD`HI`BDfYoVnWZfX@@@`\ndmND@DCdfVUrjUZjZi@LFrfFUyB@\ndeVH@HAIYWVz`@d@H\ndaG@@@rdifvxH`@C@linyA@\ndieDPLZD@HhHrREQKaVii@LBrnGd@\ndk]L@LhDeIeoYR[SZjZZdHZ`\ndmTJ`HSNd@cIICeMEjjh@`\ndaxB@@rnRV{jj`CAhSBinHG@\ndeeD@DHFR[eyiihHr`\ndaxD@@QIgUjfBJlBpj{dL\nfoQPB@F\\@DYHheBeLdRdeV]Th@@D@@P\ndk^@@@RfYU\\]Tz@@@@@LECBinFUOdZ\ndclD@@kHheCDdUKkSP@P@P\ndeV@@@RV[TYzP@@C@j[axaR\ndklHPBBPzPrJJKQEIa``bZ@bV\ngGP@DiUjaAXEXwbD@\ndmtH@@RYWUih@IhBN\\NALJiWb\\H\ndk]L@LxDMIe]eRkSZjjjh@`\ngN}@DVDfUZi@`\nfb}@`@@YRYVum[ehrRkNBf@BjjjjJ`@p\\AFJ\\EIQk^bgODZB\nffs@`@@URYVumfv^cIJmyIX@Jjjjjbh@H\nfjc@`@@ErJJQIFYJSKXgIJlyNVjjjjjjj`@pC`eF\\EIVkN|FJ\\|`T\ndg|@@DjU_eZx{BAH@@BJ\\MaBine]N~Q`\neFJHSHpP\ndiV@`J@HRfU|kahDB@CB`QnGdD\nfhy`B@J@BDifW_e\\TDkpZA`R@B@DAVF`DYpTeFl{wHA@\nfhy`B@N@BLdTTTRRVqirUmNh@`BBh@@pt@cARUhugAyCp\ngNx`DJdssTpFBsyD\nfgA`@@@ISLrotyHvk@@@@@@@P\nfdy@@@LdbbRbVbJwMAc@pU@P@@ab`\nf`q@`@@^RYWUe^cKN`@f@B@DNVB`HpRbmFm{bNB@\ngGP@DiVV`iJpJqoDH\ndeTH@@RYWVf`@j@CC`SBhYyG@\nf`qAA@A@bOQBSJ{\\ktYYt@EP@P@A@\ndif@@@RfU~F``@@pYLJf{b@H\ngJQ@@dsU`XKGbD@\ndeVD@LADfvUFVjjh@p{BinF^P`\ndmtD@@QIV[VUZh@@@pZfxYWbQ`\ndeT@`@bDfUuih@Jp@`\ndid@@DkYWaz@@@LF`j[ayB@\ngNx`LDdskUXD\nfoA@@@DkfYU]UcNz`@@@@@`\ngJX`DBdru@XS\\RH\ndeVH@DAIgeQej@@@LJrfx^Hd`\ndg~@@@RYfUWd}mh@@@@@pdBinFUwbGX\ndid@@DjYmaBH`@LJpj[nPH\ndmvD@EBdin]~F``I@CAdJfF^Ph\ndetD@@eIYe~DYZjjh@`\nfmoA`@@HWdTrTTtbRLrTfVrk^CdhymUUUUUUUUU@AauFBTYrRmFl{pXyKUoQsV\nf`qAA@A@dORBSJ{\\ktYYt@EP@P@A`IFBTUhunHxH\nfb}@`@@YRYVum[ehrRkNBf@BjjjjJ`@p\\cENBdkQk^CEOD@r\nfb}A@@@YEDeMDhTihjLUiwhyjjjjjjih@L@``cAJ\\EIQkNBgOHQ@\ndknL`LaE@HrRPzIJZ]Vh@b@CAlInd~Qh\ndeUH@JdDin]xZB@`@pIBX^IT`\ndg]L`LnDT@cIIChdieNkT@QP@P\nfgAAB@H@BDjyeUrLd[Uf`@@B@@H\nfgA@@@LdbbbTVKIBMjp@@@@@@D\ndieD`JXaCDRYgvzejX@pHLi`\ndmUL`LZDh@cIHULdeijh@p[FxYWd@\ndeV@@@RfyWahBB@CBj[agf@a@\ndcNH@DAIVYeEZ`Hb@B\ndaF@@@Rfu[j@@@LJABinIg@\ndo}D@LlMrIJJIIHiSjhHHh@H\ndid@p@qBqAqDfYun``H@H\ndcn@@@RieU~V]jB`b@C@\\JfxYWbIh\nfoA``@I@PdwJ{J|EYsP@UUD@D\ndeUH@AdDim][j@B`@pYLJfGdR\ndid@@DjUZnBAH@LFaLkayC@\nfoA@@@DjU][VgKNBAJ@@@PiXZ`cAJLdkQk\\`x\ndid@@DjYUaBHP@LFSBiayG@\nfoA@B@@QBSJzlkQegP@S@@@H\\lM@Q`eNBUhunP\\\ndif@@@RYWZZ@BP@piLJnx`B\nfoA`@@@ILkjrmFV]@AL@@@ar`\ndeUD@DXIRYvTYZjV`B\nfdyAb@HHpCpRkVYU_]Nmyj`@@B`@B\nfhiA`@@Hddjrm|jIW`mPAD@@@A@\nfhi``@L@PdrkLjn[s`mUTEAD@A@\ndeVD@LADeVUFVh@@@`\ngJX`LDdju@XP|Tp\ndifH@HAIfuxZ`@@B\ndif@PBBPRPR[e^Fh@@@pbNES\\H\ndcNHpJtIAICICHiCDedLuP@R@D\ndkn@hJBPRPrPvPNPrJPqIPpYj`@j@B\ndeV@@@Re]Xj@Bj@C@fxYxTIa@\ndmvHPBTIAIAInV_ij@HP@`\ndcNH@NCHeEDdYplAAT@D\ndkn@PBBPRPR[eW[aj@Bj@B\nfHbhA@\nf`q`@@@YIEBedhdnB]zh@J@@@@`\ndeVD@D@dfVuFVh@@@pXj[agdH\nfoAP@@@NRYWVUzLMZ@B`@`@B\nfde`@@@ISKN~rmFBTFH@@Pp@@@`j`\nfoQHA@FJuXFH`HRYUYYbTT[UjZjefZ@B\ndid@P@BJdDef_ahH@@LJ`fx^P@\ndeT@@DjY]zXFB@@pYLinGbEH\ndaF@@@RYWih@H@LBCB{bIP\ndcl@@DjYU_egX@@@@@pYLJngaLJz\ndknL`LaM@HrRRqIYPYV`@f@B\ndeTh@DiiAIgeQej@@@H\ndg^@@@rIRJEJFRoU@AT@D\nfoAp@@@PidbbvbRafJJuMT@ER@A@\ndk_@PBxpbEbDfYYUZ]NBBbT`B\ndcvB@FAEuIm[VZijX@`\ndknL@JABR[mfWSZZ@`@B\ngOq@@drm]SRA@\nfoA@`@@BRYfYWuVLyh@`@@@CBUF\\EIVcNyBp\ndaG@@@[diWRh@I@B\ndeT`@@pjrQQQUMpEAP@XDUCOHX\ndaDh@DqnAIeZfZZd@`\nfoAP`@DNAsHheEEJefBMYji`@h@B\nfoAQ`@DX@pRSJs|kSegMTs@A@@X\\QaVcV]rG@\ndg|@@DjU_eZx{BAH@@BJ\\MaLJfe]N~Qp\ndk\\H@@RYm[Watz`@@@@H\ndmvH@JAJ[g_ahHBP@`\ndmvH`BdIAIfUya``a`@pkBhUyC@\ngOy@HDfUkjj`LJlZ~P@\ndifD@LADeYWaZjj@H\ndif@@@RfU~F``@@pxDpj[nPH\ndmtH@@rJJIJEn`HJ`@pxBiae^Ig@\ndeVH@NAJ[VvF`BH@LBPfGfYt`\nfoQH`@LBUXCldTRabRrIRRcVjj`PJ`@LFPISfg\\``\ndmM@PBx@c@aJYg\\jeZdHB@B\ndaF@`H@HRVU[jjj@LNaLJf{d@\nfgA`B@K@BDifUW|TEiQj@@@J@@H\nflmAP@@LUyNQQQQEJQJI|Eish~BtuUUUUU@A@\nfnkA`@@U[dTRTrTrTtRJeFRUYpTmMUUUUTuT@D\ndaDH@@RYUih@H@LBCB{bHp\nfoAQ@@DZ@drsJkjlYsUPTDP@D\nfoQAB@C@BDifYU^gIVtz`B@DH@@`\nf`ip@@@F}dbbTRRRQhpVkN@@@BB@@B\ndklD@@QIe]e]MjZ@@@LApfFUt~QP\ndg|H@@rIQQJQZ}NfuUUTt@P\ndmtD`ATHaIe]nf`@jP@pXDpjGb]H\nfbmPB@NA@DYHhhhddmcEJ\\e[Sj@H@`jj@@H\ndg~H@HAIYeg_eNzB@@@@B\ndmvD`La@BLdabRbxUjjj@H\ngGQ@@djutA`c^HP\ndaD@@DjUZxHH@CBdpj[bQp\nf`qpB@DLxBHRYUvYZcKN`@hBB@@H\ndigDPLXXP@b`cIHUDnEZfd@`\ndaE@@@aJmUnjjh@`\nf`iA@@@YHhhhheEZBdxwj@BHHH@@`\nfluA`@@HRdrvmkZxiFlFKUUUP@U@@P\ndev@@@Re[TjFP@@@@LFaB[iy@`\ndif`@@pjGIEEEVxB``@ppjGdL\ndaE@`BhHaIfUn`H@@ppL{bTp\ndg~L@BAER[e[gzSmh@B@@@`\nfj}@P@@\\teIfVWn]{YrQhyZA`hBBJ@@H\nfj}@P@@\\LeIfVWVu{YrVhyZA`bBBJ@@H\ndmvH@DAIVUVUZjjh@peLJfxYWdP\ngNp`@dfzZj@pJMX\ndk\\D`HP@cIHXheDQgSV@@@@@pZnE]ObAH\ndeVD@FADfygFV``@@pjfxYyB@\nf`qa@@H@RVUuYUgG^h@J`@@@H\ndaD@`@qDfYVz@`@CB`pj[d\\\ndif@@@rRJEKaj@@@LJSJ[nHP`\ndifD`Na@BLddJT[ejj`B\nfhiQ`@DX@pRSJswJ}N^BuSMA@P@D\nfgAa`@N@t[HhheDTdsdeFltC@UQ@A@\neMBBHRYCAKd`\ngGPdMQDGpRUYiDe@\ngGPXHlQxIU[U@XR|VH\ndeU@@@qJYejxBHh@LJJfF^Qp\ndmuH@DXDfUgjZ@Bj@B\ndcNH@DCHheEBdnmU@@@FGIeMpkqIt\nfnk@`@@UrJIJZIPiSQHrcEZ]zNVfjjjjjZj@CCrLTxJRmFlxLTyyAh\neMJDBDe`pQyP\ndeTH`ICDRUe_af@B@bGB`Jf{fPd`\nf`qAB@O@qDfYUg_EjsjP`@@`AAC@\ndcLH`ICDRUe]^FX@J`Ha`\ndid@`@bDf[Wai@@@LJ@j[nI``\nfhi``@C@PdsrnljJW`mP@T@@@A@\nfiwpP@DVz@wliLsLj{[klrPeFCEKQmTtuAPHDQ@@FDj\\d{qV}GrF`\nfhiPb@OA``@QddabbRRvRkF\\m@@@E@@F@jLxJRmFxcpP\ndeVH`IDIAIe[ZZ@Bd@LBSJGb]H\ndmvH`IDIAIe[^f`@i`HF`\ndmv@`ABHRYWUih@IhBN\\FALJixgJ\ndaE@`FxLQIfVfifx@`\ndaFH`BxLQIe\\jffh@`\ndig@@@aDkYWaZ@@@LJPj[nI@`\neMBBHR[B\ndmtD@@QIn[VUZd@@@pYFxYWdT\nfhi@`@@HR[YfUWMypVf`@@@@@ptc`iJtZsoAyC@\ndk^H`MDIAIe[mZy]`BIjD@`\ndaF@@@RYWifef@H\ngOx@@drm[RtA`Uc^HL\ndidL@@QdfU\\jZff@LBpf{dB\ndcmH@HDLbTTRbOBnt@@@@A@\ndaF@@@RZW[jii@H\ndaFH@NAIe^f`@`@piLJny@@\nfHgPAa@\ngC`DAbZHRVhB\neMBBHRYCAGe@\neFJH\\HpXQr`\neMPBchLF^P\neMbDBDfp`',Uq='daD@@DiYZYji`@\ndaD@@DjUZxHD@@\ndaD@@DjUZxHH@@\ndaD@@DjWjXHB@@\ndaD@@DjWzXHB@@\ndaD@@DjYvxH`@@\ndaD@P@bBbDfYvzB@@@\ndaD@P@bFbDfUjz@H@@\ndaD@P@bNBDfUzZ@B@@\ndaD@`@BDeeVz`@@@\ndaDB@@InRYgrfiZ@@\ndaDD@@IIf]n``@@@\ndaDD@@QIeZfZfh@@\ndaDD@@QIe\\jZehHj@\ndaDD@@QIe\\jZfh@@\ndaDD@@QIe\\jZihHj@\ndaDD@@YIeZn`B@@@\ndaDD@@iIeenjZd@@\ndaDD@@yIe^fZVX@@\ndaDD@@yIe^f`@`@@\ndaDH@@RVU[f@@@@\ndaDH@@RVU[j@@@@\ndaDH@@RYVih@H@@\ndaDH@@RYWih@H@@\ndaDH@@RYe[hB@@@\ndaDH`NBPRYWih@H@@\ndaDH`NBlRYWih@H@@\ndaDH`NCDRYWih@H@@\ndaDL@@SdfURijZ`@\ndaDL`HS@BLddJS\\mUP@@\ndaE@@@YIeZn`B@@@\ndaE@@@yIe^f`@`@@\ndaED@DHNRYWifif@@\ndaED@DpFRYVkfjY@@\ndaEH@DXDf[Vyje`@\ndaF@@@RYe[hB@@@\ndaF@`BBHRYg[hH@@@\ndaF@`FBHRYVkh@`@@\ndaF@`NBHRYWih@H@@\ndaFD`HI`bDfYjzif`@\ndaFD`JK`BLbbbMMTtp@@\ndaFH@DAIeUnZjh@@\ndaFH@HAIYUnfjh@@\ndaGD@Dp`yIeVfZiX@@\ndadL`HS`BLddJULwKUU@@\ndax@X@bDbLbJbFbNbLbdLeUT@@\ndax@X@bDbLbJdFdNdLbdLeUT@@\ndax@X@bDbLdJbFbNdLbdLeUT@@\ndaxB@@QnR[VZY`cH\ndaxB`HSBCpRjuZj`@\ndaxD@@IIeujj@@\ndaxD@@QImUifALj`\ndaxD@@iJU^jj@@\ndaxL@@SDfUVjh@@\ndaxL`HS@BLddNbuT@@\nday@`Dp@aIfYjj@@\ndazD@LADf]Vjh@@\ndazD@LADf^Vjh@@\ndazD@LADf^fjh@@\ndazH@DAIfujj@@\ndazH@DAImUjj@@\ndazH`LPHaInVZj@@\ndcL@@DjYn}aBHbh@@\ndcL@X@bBbFbAbEbMbDfYn\x7Fijjjj@@\ndcL@X@bDbJbAbEbMbDfn^_ijjjj@@\ndcLB@@RURYYejyjieh@@\ndcLB@@RiRYyVQejjjh@@\ndcLD@@IIf]z[hHBj@@\ndcLD@@eJ[W[[j@Bk@@\ndcLD@@iJ[g]xZB@f@bX\ndcLD@@uIfUk[hBBj@@\ndcLD@@uIfU}FV`PJ@@\ndcLDHFDH`haXcXaIf[ozYjYjP@\ndcLF@@IaWTfYn\x7Fijjjj@@\ndcLF@@Rag\\bbTVTILuSUT@@\ndcLH@@RYWUZZ@Bj`@@\ndcLH@@RYWYzZ@Bj`@@\ndcLH@@RYeZvz@`j`@@\ndcLH@@rJJIJGMtAAU@@@\ndcLJB@PUuNR[eY~eijjh@@\ndcLL@@QTfvUtYZ`@h@@\ndcLL`HS@BLddJfRtjmP@P@@\ndcM@@@WIDeBddU@AMTACP\ndcMD@DTIR[fVQuhHF@@@\ndcMH`BuPBDf[U{aj@BX@@\ndcMh@DKaePR[eoVEjVfhHF@\ndcNB`BaLtOCIILeBdnmUUU@@\ndcND@DATfyeXYZ@`h@@\ndcND@LADfU[U]Zj@@@@\ndcNH@ICHhdhdYSP@UT@@\ndcO@@@aDiUm^UZh@HB@`\ndcl@@DjYU_egX@@@@@@\ndcl@@DjYn}BXVjjjd@@\ndclD@@UIfV][iuhFAH@@\ndclD@@iJYW]rnF``IhBI`\ndclD@@iJYW]rnF``Jh@@\ndcll@Dsm@iRYgeVE]ZjeZ`@\ndcllADqe@]R]{HhdhdcWRkURmT@@\ndcm@@@YJYYwhUtH@@@@@\ndcmH@DpLbbbLRQTnmU@A@@@\ndcmH@DpLbbbLRQTnmUUUP@@\ndcnL`LaA@HrRPjIKTrzmPHD@@\ndcndADkatIjYyIefY[eujeji@@\ndcnl@DsetBeIf^UXUujjUj@@\ndct@@DiUUVjjj`@\ndctB@@I]rJJJIVMRuLDE@\ndctB@@PYRYU{ViijBBP\ndctB@@RURY]VvjjZ@@\ndctBHFxYBHRHrHkprJPqREUUMT@@\ndctF@@IaWTfYn~jjj`@\ndctF@@rngTen{mjjj`@\ndctd@DrmATf^VYjji`@\ndcuD`FWi@HrQRXiSUTttDU@\ndcvB`JFUt@aJUgfjjfX@@\ndcvD@LADf^eujjj`@\ndcvD@LADf^fYjjj`@\ndcvHPF`G@WCIIEXmIUSUPQR@\ndeL@@DjYeIjGijjjj@@\ndeT@@DjWvifjih@@\ndeT@@DjYUXPbDP@@\ndeT@p@bDbLbLbdLRPsU@@@@\ndeTB@@KiRYg]nZej`@\ndeTD@@SHheDYaMUMP@@\ndeTD@@eIff\\Ijjf@@\ndeTDB@YnRYe\\YZB@@@\ndeTD`AdHaIe[jz@HX@@\ndeTD`NDHaIfVVfBA`@@\ndeTH@@RYVZfZZj`@\ndeTH@@RYe\\YZA@@@\ndeTH@@rJJIHmtA@pD]@\ndeTH@@rJJIHmtAAH@@\ndeTL@@JTfYoXXHH`@@\ndeTL`BjPkDf[W[jjjh@@\ndeU@@@EIYe^g``p@@\ndeU@@@aJWeQfj@@@@\ndeU@PBdHchaIf^VFBBH@@\ndeU@`Dp@aIgeQej@@@@\ndeUD@HDIRVUunfef`RKh\ndeV@@@RVUenh@J@@\ndeV@@@RYeun`HJ@@\ndeV@PNBHFHRYeYi``x@@\ndeV@pBBHzHfHRYgea``b@@\ndeVD@FADfygFV``@@@\ndeVD@IADfyWxV`@`@@\ndeVH@BAIf]VzB@h@@\ndeVH@FAIfUqehH@@@\ndeVH@HAIYf^f`H`@@\ndeVH@IAJ[Vvz`@h@@\ndeVH`Ax@aIfVVfBA`@@\nded@@DiUUjjj@@\nded@@Dj_VfZZ@@\nded@X@bDbBbFbAbIbDf{nijZ@@\ndedD@@QIeVVjjP@\ndeeD@DdAR[UYjjX@@\ndefD@LADf^]Zjj@@\ndefD`FFPBDiWnjjf@@\ndefD`FFPBDi]nijf@@\ndefJ`JaLFP|LddjRcUTp@@\ndet@@DjYUX^d@@@@@\ndet@@DjYUX^dHbH`@\ndetL@@jTie]rnF``J@@\ndet``Dki@HRYYUnFVjVi@@\ndev@@@rQQJHtpr@@@@@@\ndev@PL@HPHRYUTjFVj@@@@\ndevh@DJndDfVU[af@`hP@\ndevhADIadFf^R[fUnxVijY@@\ndg\\B@@Q[R[VUmgVf@HhBL`\ndg\\B@@SSRY[W[FVh@Ih@@\ndg\\D@@eIfU_Un`HJj`@@\ndg\\H`ABHRYVwUih@Jjh@@\ndg]HPAuPbBbDfYw[fzB@ij@@\ndg^B@BAMoHiieDeBimU@DP@@\ndg^L`LxY@HrQQYJEIYUSRuUAFP\ndgl@@DiUUUZjjjh@@\ndglBA@RUSe{HihheDbtuSUAFH\ndglBPHRU@HhHrRPrIRIkUTmP@@\ndglD@@QIgV]YjfjZBJ`\ndglD@@QImUUUjjjj@@\ndglFPHkivpqLqDen{nzjjjh@@\ndgmB@LxDWTfU[{Vjjfh@@\ndgnD@KADfuUUVjjjh@@\ndgnD`H[`BLdTRbJRUUMUT@@\ndg|D@@OIEEHhfmPkmAU@T`@@\ndg|DPFDH`haIf[oWiNyjY`@@@\ndg|H@@RVUvU[cn`@`@@@@\ndg|H@@RYfUWd}mh@@@@@@\ndg|L@@ildTRbrJQTJtEAEL@@\ndg|d@Dq]@\\bbbbfJSSimUSTs@@\ndg|l@Dq]@[rJJJJXiMNfuUMSL@@\ndg}@@@aJVYU^Svv`@@@@`H\ndg}@@@mJYeU|]Tz@@@H@@\ndg}D@AlBRYgU][iVB@jjD@@\ndg}D@DpCrJJHqIYERzuT@EP@@\ndg}D`LHU@HrQQYJJEYQwSMMMR@@\ndg}H@DHDfYV]rX{Zi``H@@\ndg}HPAuPbBbDfYw[fx{``JZb@@\ndg}L`FWSl@cIEHhbeEc]MUUTsP@@\ndg}L`JXiTIAIf]VunNzVjZj`@\ndg~D@EADfufUqT{ZZP`HBL`\ndg~L@BAER[e[gzSmh@B@@@@\ndg~L@IAKR[Ye]z]MjdHB`@@\ndiD@`@RdjeVjj`@\ndiDB`HSB@HrRPiIZj`@\ndiDD@@GIEHjjjj@@\ndiDJHDpnDAbHbaahcIIJiIZe@@\ndiDL@@xTiUVjj`@\ndiDLPBhPbFbLbbbeiZdHQ@\ndiDNPHSB[a@XhXrRPzQZe`@\ndiE@@@sIDhcFZj@@\ndiF@PHApiprRQVRjj`@\ndiFD@AADfuUjj`@\ndiFD@LADf^Yjj`@\ndiFD@LADf^]jj`@\ndiFD`JxPBLbdTljjX@@\ndiFDpAk`bDbLbLbdJTjjX@@\ndiTH@@RfU|kahDB@@\ndiTL@@X\\dRRaaJzjZj@@\ndiV@@@RfU|kahDB@@\ndid@@DjUZnBAH@@\ndid@@DjUfaBB`@@\ndid@@DjYUaBHP@@\ndid@@LdbRQk``R@@\ndid@@LdbbQxXF@@@\ndidD@@IIf][hHB@@\ndidH@@RYVZZ@B`@@\ndidH@@RYVzZ@B`@@\ndidH@@RYeVz@``@@\ndidH@@RYevz@``@@\ndidH@@RYfVF@b@@@\ndidH@@RYm^Fh@@@@\ndidHHFBHJHzHFHRYgljZfh@@\ndidH`DBHR[e^FX@@@@\ndidL@@IdfYoa`b@@@\ndidL@@RdfV^fZjj@@\ndidL@@SdfVTjZfZ@@\ndidL@@pTee^fZZi@@\ndidh@DKaAInV[fiZ`@\ndie@`HxGCIIEJnFjiX@@\ndieD@DpFRYVZyjfd@@\ndieD`JXaBPRYgvzejX@@\ndieD`LIN@HRZufFZid@@\ndieH`Dq`BDfUfnZii@@\ndie`@@pjX\\dTTUk`Hb@@\ndif@@@rJJIEn`HH@@\ndif@@@rRJEKaj@@@@\ndifH@AAJ[W[j@B@@\ndifH@BAIfuxV`@@@\ndifH@JAJ[gxZB@@@\ndigH@DK`R[e^Eh@@@@\ndigL@Ds`XTfUfn`BH@@\ndkLB`HSB@HrRPiIIIZjjh@@\ndkLH@@RUUUVjjjh@@\ndkMB@LxDeTfU]mZjij@@\ndkNF@BAIWSR[YVYjjfX@@\ndk\\@`@bDfYYwZ]NB@@@@@\ndk\\B@@SSRYVuVfeVi@Bh@@\ndk\\D@@QIee}RkaZfjjh@@\ndk\\D@@wHhhhhbfESZAhD`@@\ndk\\D@@wHhhhhbfESZBhD`@@\ndk\\H@@RYWYVftx@H@@@@\ndk\\H@@RYeg]itxH@@@@@\ndk\\L@@x|bbbTTJZUuhFHF@@\ndk\\b@Dsm@iMIf^UvE]Zjeih@@\ndk\\d@Dq]@\\bbbbfJZ]MjjZe`@\ndk\\d@DsmB\\bbbTrQXUujjUj`@\ndk]D@JxCRe]YTjtzjjjj@@\ndk]D`LHY@HRf]VwJtzYifi@@\ndk]H@DpLbbbLRVJeujh@J@@\ndk]H`FVPbDfUonkmN@Hfh`@\ndk^@@@RfYU\\]Tzjjjj@@\ndk^@@@rQQRJJjaTzBjBD@@\ndk^D@IADfvYWz]MjdHB@@\ndk^d@DXYtCRYf[WaWVjfji@@\ndk^d@DkaTMRYYe]neVjZji@@\ndk^d@DkatCRYYfUngVjZji@@\ndk_D@DHPuIeeevySZ`hD`@@\ndklB@@PcR[me]]ZZ@B@@\ndklB@@PcR[me]]Zj@B@@\ndklB@@QSrJYJIJF]ZX@b@cH\ndklB@@QmR[fUxUZBBF@@\ndklD@@MJ[eZ~F`HJh@@\ndklD@@eJ[Vvfz`@jh@@\ndklH@@RYfWua`Hbe@@\ndklL`HS@BLddJbRvWUjB@`@@\ndkl`@@kaRe[vTf@HZj@ah\ndkm@@@GHhhhhdvf@bbh@@\ndkmD@DHCRYvUvUZh@J@@\ndkmD@DTCRUfWtYV@`e@@\ndkmD@DdCrIJJIPxUV@bE@@\ndkmDpDgSBHjHVHRYYmYn`HJf@@\ndkmH`NVPbDfUunih@JZ`@@\ndkn@`D@HRUUYWSVj`@@@\ndknD@LALbbRbaRtvjh@@@@\ndknH@DCHhmEEEYuj@bH@@\ndknL`IaMADRge][aj@Bf@@\ndk~@@@RfYU_JGUN`@@B@@@\ndmL@@DjYUVGi@@@`@@\ndmL@`@VDifU^FUifje@@\ndmL@`@ZDifU^FUifje@@\ndmLB@@RURYUVJaejVjh@@\ndmLD@@QIe[VfeVi@B@@\ndmLD@@QIe[VfeVj@B@@\ndmLH@@RYVuiiV@BjH@@\ndmLH@@RYe~Ifyjjjh@@\ndmLH@@RYiYKnUjjjh@@\ndmLL@@SdfVUrjUZ`PH@@\ndmLd@DpYBdfV]VzUZZjV@@\ndmLd@Dqe@TfUeZzUZZeZ@@\ndmM@@@yJUntfePBIhP@\ndmN@pN@H`HPHrRPqIZneUhDB@@\ndmND@BA\\bbbReInFjZjd@@\ndmNH@BAIfUmiEX@@@@@\ndmNH@NAIYe^neZdHB@@\ndmNh@DkaTDfVYVzUZiZi@@\ndmO@@@SdfVUrjUZ`PH@@\ndmTB`HSB@HrRPiQQZjj@@\ndmTH@@RUUUjjj`@\ndmU@pLD@a@c`cHheEKFjfh@@\ndmV@HLBHQpYpVHrIRHrJjjj@@\ndmt@H@bAdIdEdDfUvjZ@Bj@@\ndmt@X@bBbFbAbIbEbDfYojzfZj`@\ndmtD@@QIgYVUZh@@@@\ndmtH@@RYWUih@Jh@@\ndmtH@@RYe[[hBBh@@\ndmtH@@RYeeZVfjj@@\ndmtH@@RYeeZZjjj@@\ndmtH@@RYe~[ffjZ@@\ndmtH@@Rfuu[j@BXBA`\ndmtH@@rJJIHin`HJ`@@\ndmuD`LVD@HrRRqIXYV`@`@@\ndmuH`Dq`BLbbRbJkfjZi@@\ndmvD@DATf^Uqej@B@@\ndmvD`La@BLddlTReUhB@@@\ndmvHPEHJsjsHhhmDVFBBK@@\ndmvHPEHLSlSHhhmDVFBBK@@\ndmvH`ITICHhdhdZZ@Bj@@\ndmvL`BaL@HrRRjIJUVjjh@@\ndnDH@@ReVDijiZ@@\ndnDH@@ReVDijjj@@\ndo\\H@@RUUUUZjjjj@@\ndo^HpAxH`hb`aIevyffjZjhHa@\ndo|J@@S[_HheDdeDYMjBBb`@@\ndo|L@@RtfUVuwSZjp@h@@\ndo|L@@RtfvYWwSZjA@h@@\ndo|L@@iTinU]_ihHHjXBCP\ndo|L@@iTinU]_ihHHjXBC`\ndo~L@MAER[e]mnEh@Ij`@@\ndo~``LKad@aJ[V]Y[j@Bjj@@\neF@Hh@\neFAAD`bJ@\neFAADdRJ@\neFABD`bJ@\neFABHhbL@\neFACDlRL@\neFBBHc@@\neFBBlc@@\neFBCDc@@\neFHBJ@\neFPBc@@\neFbHbHp@\neFhXNic@@\neM@Hv@\neMA@JXaIh@\neMABHXaIh@\neMABHYAIh@\neMB@Jch@\neMBBHRZ@\neMBBPRY@\neMCALhabHz@\neMDARV@\neMFI@bMP@\neMFIGBMP@\neMFiDqzN`@\neMHAIX@\neMHAId@\neMHAIh@\neMIdEJHIcd@\neMJDBDeP@\neMPBch@\neM`AIx@\neMbDbDfp@\neMhDRV@\neO@Hyj@\neOB@Hcfh@\neOBBHcfh@\neOBCDcfh@\neOHBNZ`@\neOPBcfX@\neO`BNZ`@\nfH`D@@\nfH`T@@\nfH`X@@\nfH`p@@\nfHa@@@\nfHap@@\nfHbT@@\nfHcD@@\nfHcP@@\nfHcT@@\nfHcd@@\nfHdH@@\nfHdP@@\nfHdd@@\nfHdp@@\nfHep@@\nfHfX@@\nfHgP@@\nfHgd@@\nfHgh@@\nfHpXT@\nfHpp\\@\nfI@@\nfJ@@\nf`a@P@@Ht[HheDhmD\\jsTE@qP@@\nf`a@`@@FrJJIJQJrLy@PUUT@@@\nf`aA@@@ILsKWRpTADUUP@@\nf`aAb@NFlBHrJIKIRJUTY@AUST@@@\nf`aQC@IVLBPQHXdLbdLRTvf`eUPADu@@@\nf`ahB@LDxJP@aJ[V]Yf\\h@Jjj@@@\nf`ahB@LDxJP@aJ[V]ZV\\h@Jjj@@@\nf`aq@@DV\\CHheEDcddkSTE@UP@@\nf`i@P@@HD[HhdihdhUSbkN|uHPTUB@@\nf`i@P@@HTYIe[VUZLeiwfi@HjhP@@\nf`i@`@@DRYfyU]`mNmyi`@@B@@@\nf`i@`@@HrJSQQQIH|UiuoMP@P@P@@@\nf`i@`@@VRYfYU]`eNMyh@`AB@@@\nf`i@a@ARADDbDfYuUUYqVg^``Jh@@@@\nf`i@a@BBADNbLbbbRfaRcIBU[tDE@AA@@@\nf`i@a@FRAD^bDfUm[UirRkN`BJ@BH@@\nf`i@a@FRAD^bDfUm[WirRkN`BJ@BH@@\nf`iA@@@YHhhheLUfBdYwjBb@@H@@@\nf`iA`@@HldrlkZuFJMYsTaADtP@@\nf`iA`@@HmdTRTRTrQhqQkNZdHHfb@@\nf`iQA@B\\|@HpDISLzsnRdcN}TaAQUD@@\nf`ih@@@\\eYvRJJJJUKJgEZMX@HHFhH@@\nf`ip@@@XTeLwOvmNJt{pPQT@A@@@\nf`q@@@DjUgm_hJs``hB@`@@\nf`q@`@@LRYfWg^Qg^ZB`@@@@@\nf`qA`@@FmdbbdTUrRiIVjBh`jd@@\nf`qA`@@Hpdrlrj~gV|uP@@@@@@\nf`qHC@DXxHDPrHMDYEDeDUEEj\\]z@HhH@@@@\nf`qP@@@PrQJJJIIFJlYsP@P@P@`B@\nf`qPB@DX@DILwLjoiuoMUUUUT@@\nf`qP`@DBAKHheHdhbmPSoMUTa@P@@\nf`q`@@@YEEDbdhdf\\]z@`f@@@PE@\nf`q`B@B@bDfYwVUYqwhHBX@@A@T@\nf`q`B@O@dDfUuYWhrsh@I`@`ACd@\nf`q```JBTBHQDXbBAFQREQQSYJVg^jjfX@@@@\nf`q```JBTBHQDXbBAFQREQQSYJVg^jjfZjd@@\nf`q`b@LPP@HrRPjJIIKZ][ru@@@@@@@\nf`qa@@D@RYyV{TRg^Z`B@@@@@\nf`qa@@D@RYyeg^Qg^Z``@@@@@\nf`qh@@@XirVRJJZGJII`g^BHa`@@DAP\nf`qi`@DTxIPC^rJIQQJiILyISUKMUU@@@\nf`qp@@@Hpds\\rj~gV|uP@@@@@@\nf`qpa@LR}A@@`PQddaTTRRUtxuejBPha@@\nf`y@@@LdbbbbbfkEBMIuo@@@@@@@@@\nf`~`b`KLLBHHDTBNA@`aIneUYYjjjZ`@@\nfbc@@@LdbbbTRLRqWEBMjsoIs@AA@UUPT@@\nfbc@@@LdbbbbbbQQsEBMKuhir@@@@@@@@@@\nfbc@@@LdbbbbbbcJwEB]Hu`ir@@@@@@@@@@\nfbc`@@@ISLrj}{dirVgVBgLmUTrB@m@@@\nfbc`@@@YHheEMDXlijREhuoQSd@@@@@@@@@@\nfbe@P@@HM[HheDhdhmbdisTEALUT@@@\nfbe@`@@HRYWUUUUIQjjjjjj`@@\nfbe`P@N@P[vQQJJKQRzZIUfjijjZi@@\nfbm@@@DjYVWV}~ZlENXI@H@@@@@@\nfbmAP@@BLENQQQQQQQYG[bm^]Eh@bHAJ`@@\nfbmAP@@BUGNQQQQQQPeIG`mNBehIb@@b`@@\nfbmH`@EVBdGlbbbbbbTRacAZb{KPADDLQ@@@\nfbmI@@DTdhFQQIQIQHqIY`iJSejfYjB@h@@\nfbmPB@NA@DYHhhhddmcEJ\\e[Sj@H@`jj@@@\nfbmPb@AJ|dDPdrmrljoSfgQs@DMPTBD@@@\nfbmQB@AJRBHRYVyV[WisShy`BFhJAB@@@\nfbmp@@@V|eLsJzo]SdcZ\\@@@AEU@@@@\nfbu@@@LdbbRbVbrQwMQS@pU@TA@BFH\nfbu@@@LdbbRbbtRJOCIs@pUPE@@BFH\nfbu@`@@YRYWYeg_hrJX@Ij``H@Py@\nfbu@`@@YRYWYeg_hrJX@JY``H@HDj`\nfbua@@D@rJJJPjJYIK^SGKUUX@@@@@@\nfbupb@LVcA@@cIIBmDeLhThkUgKU@ASUQ@@@\nfby@`@@HR[UUUUUZjjfjj`PT`\nfby@`@@HR[UUUUUZjjjjj`@@\nfbya@@D@R[UUUUUZjjjjj`@@\nfb}@@@DjYee\x7F]^RD[S`q@@@@@B`@@@\nfb}@P@@H]gHheDheEeD\\jugQRgKUUUUUUU@@@\nfb}@`@@LrJJJJHyISI\\dkSoIsTED@@@H@@@\nfb}@`@@LrJJJJHyJIK\\dkS`isTED@@@H@@@\nfb}@`@@LrJJJJKIIIH|Djw`isTEP@@@@@@@\nfb}@`@@YRYVum[ehrRkNBf@BjjjjJ`@@\nfb}@`@@YrJJQIFYJYKDyIUgQRuUUU@AU@@@\nfb}@`D@HQvQSRJIJUJYHRiZ]yNVjjjjjjj@@\nfb}A@@@IS\\lj{j|DjsoIsP@@@@@@@@@\nfb}A@@@IS\\l~kZ|dkSoIsP@@@@@@@@@\nfb}A@@@YHihhhdeCenJtzpTyh@@@@@@@@@\nfb}P@@@RRVYfU{wyKSkASg`@@@@@@@@@\nfb}`@@@YHhhhhdhecjRUXt_I@@@@@``@@@\nfb}`B@A@dDfUmeumZLeiwdy`@@`@@@@@@\nfb}`B@A@dDfUmevUzLeipTy`@@`@@@@@@\nfb}`B@B@dLbbbRfbtVbKIBuxJ\\p@@A@@@@@@\nfb}a@@I@RVYfU{wyKSkASe`@@@@@@@@@\nfde@@P@QAHadQJHuDFb@qDXcHhhiMMEciS`eN}MUUUUUP@@\nfde@@P@QAHadQJHuDFbOQDXcHhhiMMEcCS`eNCMUUUUUP@@\nfde@P@@BLGHhhhhhhlcqVoNbt@QD@d@@@\nfde@``ARADDb@qDXaIf]UUUYqVg^``Jh@h@@@\nfde@``BRADLb@qDXcHhhheCBdeqTmN}ADT@AQ@@@\nfdeQ@@DFAdTRbRRTJbIQVg^ZiBBZhDDJP\nfde``@A@BdsLsKslUiuhm@DPAI@@@\nfde`b@H\\d@HRfYfWwQJMxLZjfZ``H@@\nfdi@P@@HM[HheDhdhkdisTE@lT`@@\nfdiA`@@HedTtbJRbbV|Ejjj@``@@\nfdiQb@LR``P@cIIBhheedmjru@pQU@@@\nfdq@`@@HR[UUUUVjjjjj@@\nfdqA`@@LdeJl{jjtuSSLtBbGbiLnP\nfdqPP@LQ@`p\\bbRTTaVRcUUUUS@@@\nfdqaQ@JDT{pPQEoCHhhihldYMRuUUMPHHP\nfdu@@@DjYee]}daZlGtP@@@@@@@@@\nfdu@@@DjYee]}faRlGtP@@@@@@@@@\nfdu@@@DjYee\x7F_daFtxLP@@@@@@@@@\nfdu@@@LdbbbbbTUHeRlXOh`@@@@@@@@@\nfdu@A`@QAHadQJHCDQbLbbbdtqTVeA\\TYuoSUUUUUT`@@\nfduAC`H`bBQCHbTQjHMD^bLbbbdttVLUNBTxHXsUUUUUU@@@\nfdu`@@@YHheEhTddj\\EHu`q@@@@@@@@@\nfdu`@@@YHhhhdeCejBd[S`q@@@@@@@@@\nfdu`@@@YHhhhhdhbZRUXp_Q@@@@@@@@@\nfdy@A`@XaLPVH[DCbQqDfYuygUgG^``JB@h@@@\nfdy@A`@XaLQfHKDUbNqDfYn\x7Ff_d`q`bHhH@@@@\nfdyAP@@BUhNQQQQQQDqIdgAZYjZZj`@@\nfdyAa@MAbBHQDIM_JztvBRUUP@UTt@@@\nfdyP`@AR@EJ[WUe]Yqwj@Bh@J@@@\nfdyhP@DTxIPCAcdTRbbURTRYrRfjVZjjh@@\nfdyi`@DTxIP@qrJIQQJiJILyISUKMUUT@@\nfdyqb@LFcA@`AFRREQQIJH{WcVVhIBhhP@@\nfgA@@@DjYU_VByHu`@@@@@@@\nfgA@@@LdbbbTVKIBMjp@@@@@@@\nfgA@P@@HEkHheHeEBRdmFluUMUU@@@\nfgA@`@@\\RfYe_irQmVh@`@H@@@\nfgAA`@@HLdrmlkQdmFluHADq@@@\nfgAH@@@XhiJYmg]gAJMXH`jeF@@\nfgAP@@@\\RfYe_irQmV@@@@@@@@\nfgAp@@@XheLvnjs`iFlDPT@@@@@\nfgApB@LLx@HRevUUpPTcViYj@H@@@\nfha@R@HHpPG`eUjjjjuUUUU@@@\nfhaH@@@\\DyJUUUWVjjjjj@@\nfhep`@BLT@NQQQQQQRq\\tJQkN|uLuUTuP@@\nfhi@B`@QAHadPzHCDYEEELUDeCpUoPQA@DP@@@\nfhi@`@@RrJJIIQFYHVoAZBA`@@@@\nfhi@c@ARAD\\B@qFQQIYRFIKTX{t@Dl@D@@@\nfhiA`@@B|dsLro~jqgM@D@AP@@@\nfhiHB@EZLDDQdTRVTTTQUFF]@A@QE@@@\nfhiHa@LTdFB@A@`cIIBhhddmmNMYZ`hJJD@@\nfhiIP@DXxHDc^CdTRbfaTVUNZltuUUUL@@\nfhiP@@@ArJJIEJIJYgKN`HJ@B`@@\nfhiPA@BAADNbLbbbrrbRaYrwhHB`BH@@@\nfhiP`@DZAyIgVYW^VkNZjdHBh@@\nfhiPb@OA``@QddabbRRvRkF\\m@@@E@@@\nfhiQA@BADBH]DYEEEeeDeBseoPPE@DP@@@\nfhiQ`@DX@pRSJswJ}N^BuSMA@P@@\nfhia@@E@rQQQHyJIHToAhJ@h@@@@\nfhia@@J@RfywVUxKpZB@f@@@Pe@\nfhipS@IZCpSo@bBAA`cHhheHeTiFJlFBAXjVH@@\nfhiq@@DZBCHheEDeDceNmyj@@`B@@@\nfhiqB@IF]hDadTTTbfLVWE^CMUTmRt@@\nfhiqP@DXxBQoArJIQSPjKJ`mVZZjjjf@@\nfhiqb@LJMAC`AFRRVIKIQKBd{rt@EMLP@@\nfhq@`@@NrQQJIJYHxRjBHbjh@@\nfhqXB@J\\dZpPAFRIJIJJqIBeUSMMUpHJ`\nfhy@C`@QAHadQJHuDFbOQFQQQRZZKFTxIS`sUUUUMP@@\nfhy@`@@HRY[fUWpPwgAZj@@@H@@@\nfhy@c@ARADDb@qBSLzjjkNJt{tDAU@D@@@\nfhyA@`A@bBQGh`LPdsNjjjsdeV]A@U@E@@@@\nfhyA`@@BMdTTTTTTVoEZ|xL@A@P@@@@\nfhyA`@@B|dsLsKnqVgVBt@Q@BP@@@\nfhyA`@@Hldrlk[oQbcN|uHPQUP`@@\nfhyH@@@XxkIEDeDehTjBYspP`Xj@@`@@\nfhyI`@LJMxD`yHhhhUEedRfcN}UUTBAT@@@\nfhyPA@B\\@DXBDif]WmRTekpZdHHjjH@@\nfhy``@A@|dsLsKnqVgVBt@Q@BP@@@\nfhyaB@K^@DISLjo{XhKRcT@@@Tt@@@\nfhyh@@@\\e[vRJJJJUIITxkQk@AA@uPP@@\nfkAA`@@TTeJwsLDDXkVcUUUUUP@@\nfle@B`@QAhbtPzHSDYEEDeCEHe\\s`XpDEATA@@@\nfle@Q`OAbdDPRHYDRbMQAhcHhhiMMEcEJ\\DjZjjjfZ`@@\nfleHb@LBdfB@AFRREQQIIZIZ][ru@@@EP@@@\nfleI@BHTDh@NaNQJJIIKQPi\\eEMA@eT@@@ar@\nfleP`@DA@eIVUue]^B]yX@J@BT@@@\nfle`@@@YIEDeDThll[tTuSPP@@@BAH\nfle``@A@Pdrrr\x7FrjcQRuAPI@@@HJ`\nfle``@D`TeMrnkZ\x7FAAcPDET@D@@@\nflehPBDX}EHAJCFT}dTTRTRRJRsNFluTtuUKP@@\nfli@@@LdbVRRbbjTjjjjjj@@\nfliA@@@IJjjjjkUUUUUT@@\nfliAq@LDhkSo@BH^FRJQSIQUQ[TuUMSU@@@\nflm@@@DjYee]\x7FYhTkA}D@@@@@H@@@\nflm@@@DjYee\x7F]yHQmNCD@@@@@H@@@\nflm@@@DjYee\x7F_yHQmNCD@@@@@H@@@\nflm@@@LdbbbTRrJRxhQmNBd@@@@@`@@@\nflm@@@LdbbbTRvRQyHQmNbd@@@@H@@@@\nflm@@@LdbbbTVRcRXhQmNbd@@@B@@@@@\nflm@@@LdbbbbbTQnEjUcA}D@@@@@H@@@\nflm@@@LdbbbbbTUVDjUcA|d@@@@@`@@@\nflm@B@@AFQQQQQIJGKLdjqoQS@@@@@@@@@@\nflm@B@@RFQQQJKPiIILEIUgQS@@@@@@@@@@\nflm@B@@XfQQQQQIJGKLdjqoQS@@@@@@@@@@\nflmA@@@ILrrknjsdcV]EL@@@@@@@@@\nflm`@@@YHhhhdeElcPTcZ}EH@@@A@@@@@\nflm`@@@YHhhhhdhecRRkFcyH@@@@D@@@@\nflmaA@KQ@DPdLdtTRbtJfRyIVkNbfifZiZej@@\nfluAP@@B\\FJSLsLoJlUkudm@DPDLP@@@\nfluAP@@HdkrSKLvzjlxYtTmT`pQUP`@@\nfluAc`OARBHIDLbIQFh`tQDHrJJJSSQXiZgAJCFZjjjjfh@@\nfluPP@DTAsUlbbTTTlTvRXHshiZjZfX@`@@\nfluP`@DX@qIeVyeUrT]DJVjZ@@@H@@@\nfluQB@DXX@HRYUnYU\\eGQBejf`@@B@@@\nflu`a@BJ|CDebYEEEheDXdlpQkQSAADtB@P@@\nfluaP@E@EKt\\bbbbbRJRRXhshiZdF`bjb@@\nflua`@O@QGHheEDdebeJ\\UXsfh@@@ih@@@\nfluib@DTxIQhi@DYEDhheUDddsdeN|uRsUUUT@@\nflux@@@XIRRmYHhhldiUDdsdcQRA@AQ@@@@@\nflyAP@@HLxJSJkZrkzluT`DUT@@@\nfl}A@@@IRlrjkoAENJl[tTuP@@@@T@@@\nfl}A@@@YEEDhih]DdpTmZlxLTmSSTuJsT@@\nfl}A@@@YHdhheDddcAENJ]ZHTt@@@@@@@@@\nfl}AA`D`bBQCHbTPFHcDYEEEIibheEPWEFl{wdsUUUUUUT@@\nfl}aA@HF@DIdLdRbbRRtJR``cIFuXLZjjZjeif@@\nfoA@`@@VRfUYu^JLz``@B@@@\nfoAAB@A@bDfUmyVcKN`@j@@@@@\nfoAA`@@HXdrkkJdiYsTtp@PBGH\nfoAP@@@NRYWVUzLMZ@B`@`@@\nfoAP@@@XReeV]qZlyfjjjj@@\nfoAPB@KN@DISLjohmJMP@@AP@@@\nfoAPQ@LF`aV`AD`cIICDTiCJLlyZjiZj@@\nfoAP`@BZ@aInvYWejsfjiB@`@@\nfoA`@@@IKLrjzkF]u@@@@@@@\nfoA`@@@ILkjrmFV]@AL@@@ar@\nfoA`@@@IML|{wEFmUSA@R@@@\nfoA``@L@QdTVbbbblmV\\u@A@@@@@\nfoA`a@AZlBHYDYEDeDXhiSagPAESLP@@\nfoAaB@G\\ADILkkJ}FFm@AP@P@@@\nfoAa`@J@PIImeeWyJsfjZb@@@@\nfoAp@@@P\\eKLjorMjsP@@A@B@H\nfoAq`@DXxBSlbbTTtJVhKQffjjjX@@\nfoIA@@@IRlkZ|DTyKUgUSUUJs@@@\nfoIA@@@IRlrj|DTxjqgUP@@AP@@@\nfoQ@b@BBADYEEDeMBdrPeV]AAP@A@@@\nfoQ@b@FRADILkZvmNRUYt@QP@D@@@\nfoQA`@@HldrmlrtYKUgMR@QSD@@\nfoQH@@@RM[IEDhmBdj\\DkQ`@@ajB@@\nfoQH@@@XhiJYmg^YpRcVBHJiXX@@\nfoQP@@@FRfYeUz\\e[S`@@@`@@@\nfoQP@@@FRfYeUz\\e[S`@`JbB@@\nfoQa@@N@rQQQQJKGbiVLz`BB@D@@@\nfoQp@@@XdeLv{ZtyIUgAAE@@P@@\nfoQp@@@XidbbfbQRSNRUYpQAT@@@@@\ngBQ@@eJuT@@\ngBX@@eLUT@@\ngC`@Die@@\ngC`@H{P@\ngC`DADJHRZd@\ngC`DADZHRVXRP\ngC`DADZHRVx@\ngC`DAHJPRZd@\ngC`HADIKLIH\ngC`HAVIMT@@\ngC`HAbIKJ@@\ngC`LADJHtPduP@\ngC`LADJPt`duP@\ngC`LAVJluXduP@\ngC``Adej@@\ngCa@@dkH@\ngCa@@dmH@\ngCa@@dmP@\ngCa@@dmX@\ngCa@@dsP@\ngCa@@duP@\ngCaA@NRVd@\ngCaHH@bNt@@\ngCah@mJAIj`@\ngCahHl@bNj@@\ngCahHlGBNt@@\ngCahHlHRNj@@\ngCahhlAa]ncm@@\ngCaihlLr\\nwQz`@\ngCd@ADij@@\ngCd@ADkZ@@\ngCdDI`BHDRZh@\ngCh@@doH@\ngCh@@duP@\ngCi@DDfZ@@\ngFp@DiTt@@@\ngFp`@dfTujh@\ngFp`@df_Ejh@\ngFp`ATiTvjh@\ngFq@@eOKUU@@\ngFq`@ldrfmT`@\ngFr@ACTi[FZd@\ngFr@ACTi_FVh@\ngFt@ATigVVh@\ngFt@AdigUjX@\ngFtHE`DILikUP@\ngFu@E`drfmU@@\ngFx@@eJf`@@@\ngFx`LDdrfmU@@\ngFy@DDfXujh@\ngFy@JDiTvjh@\ngFy@LDeXvjh@\ngFy@LDi^Jnh@\ngGP@DiVj`@\ngGPBADJHLQXaInih@\ngGPBADJHtQXcHiCUp@\ngGPBAHJPLaYAInih@\ngGPBAHJPtaYCHiCUP@\ngGP`@TfYi`@\ngGP`ATeVj`@\ngGP`ATeVn`@\ngGP`ATiVj`@\ngGQ@@dkUT@@\ngGQ@@eMUT@@\ngGQLJHaQFbLbdMT`@\ngGQXHlZHROjj@@\ngGQ`@bdwMT@@\ngGQhHl@cIIBmP@\ngGQhHlLSIHTmP@\ngGQhHlOAJmZh@\ngGT@ADiVj`@\ngGU@E`dmmT@@\ngGXHD@aIUVd@\ngGXLJHaQFbLbdMU`@\ngGX`LDdsmT@@\ngGY@HDeVZaI@\ngGY@HDefZaH`\ngGYHLaDIMtu@@\ngG]HHjPDIJuS@@\ngJP@DjYd@\ngJPBADJHtPXaIjj@@\ngJPBADJHtPYAIjj@@\ngJPDADQpRZj`@\ngJPHADIKSP@\ngJPHADILth@\ngJPHAVILuP@\ngJPLADJHLPdwS@@\ngJPXHlPDQzt@@\ngJPXHlPLQzt@@\ngJP`@TeVd@\ngJP`@TeZh@\ngJP`@TfVd@\ngJP`@deVh@\ngJP`@dfvd@\ngJPlLQDPHTPduR`@\ngJQ@@dlu@@\ngJQ@@dmU@@\ngJQ@@duU@@\ngJQ@@eKU@@\ngJQDHG@nBUMT@@\ngJQHBHaIfe@@\ngJQHBLQIfe@@\ngJQ`@bdvu@@\ngJT@ADiYhRP\ngJT`E`TfVh@\ngJX@@dkU@@\ngJX@@dms@@\ngJX@@eKU@@\ngJX`LDdru@@\ngJY@DDfvd@\ngJYHC`DIKTp@\ngKP`@df\\Vj@@\ngKQ@@eKcRp@\ngKQ@@eKcUP@\ngKX@@eKcUP@\ngK\\@ADeKbuH@\ngNpXHlPDYIHTmT@@\ngNp`@dfVZf@@\ngNp`@df]Zi@@\ngNp`@dfzZj@@\ngNp`@teUZi@@\ngNplJqDJHtQdTaeUP@\ngNpmJqDJHtP~rJPrjX@\ngNq@@dssUP@\ngNq`AVeJmUP@\ngNqhHl@cIICej`@\ngNqhHlOAJkVj`@\ngNtDLpDDHRevnl@\ngNtHE`DILruT@@\ngNx@@eRmUP@\ngNx`LDdskUH@\ngNx`LDdskUP@\ngNx`LDdssUP@\ngNx`LFdjmUP@\ngNy`LDtf]Zj@@\ngN|@ADeJkUP``\ngN}HEbpDILzuR@@\ngOp@DjWkB@@@\ngOpHADILkW@@@@\ngOpXHlPDYIHUVmU@@\ngOp`@dfUMZf`@\ngOp`@dfVqZj`@\ngOp`@tiguif`@\ngOp`AdeekZZP@\ngOphH`DYIHUVmT`@\ngOq@@drm[UT@@\ngOq@@drm\\@@@@\ngOq`@fdrikTl@@\ngOq`@fdrikUL@@\ngOqhHl@cIIBjujh@\ngOtHE`DILl[MT`@\ngOx@@drm\\@@@@\ngOx@@drm]UT@@\ngOx@@eJqh@P@@\ngOxHDHaIeZx@@@\ngOy@DDfYKZj`@\ngOz@ACVeKNLuR@@\ngO|HDVHaIeZx@@@';var ZD=fI(134);DG(169,1,{},_q);var $D=fI(169);DG(170,1,{},er);_.assessDruglikeness=function fr(a){var b;return fq(this.a,(b=a.a,bA(cr),b));};_.getDetail=function gr(){return eA(this.a.a);};_.getDruglikenessString=function hr(a){return gq(a.a);};var br=-999,cr;var _D=fI(170);DG(39,1,{39:1},qu,ru);_.addAtom=function su(a){return Ih(this.a,a);};_.addBond=function tu(a,b,c){return Jh(this.a,a,b,c);};_.addFragment=function uu(a,b,c){Fk(this.a,a.a,b,c);};_.ib=function vu(){Jp(this.a);};_.jb=function wu(a){Kp(this.a,a);};_.addImplicitHydrogens=function xu(a){a===undefined?this.ib():this.jb(a);};_.addMissingChirality=function yu(){rp(this.a);};_.addMolecule=function zu(a){return Kh(this.a,a.a);};_.addOrChangeAtom=function Au(a,b,c,d,e,f){return Lh(this.a,a,b,c,d,e,f);};_.addOrChangeBond=function Bu(a,b,c){return Mh(this.a,a,b,c);};_.addRing=function Cu(a,b,c,d){return Nh(this.a,a,b,c,d);};_.addRingToAtom=function Du(a,b,c){return Oh(this.a,a,b,c);};_.addRingToBond=function Eu(a,b,c){return Ph(this.a,a,b,c);};_.addSubstituent=function Fu(a,b){return Qh(this.a,a.a,b);};_.canonizeCharge=function Gu(a){return Hk(this.a,a);};_.changeAtom=function Hu(a,b,c,d,e){return Rh(this.a,a,b,c,d,e);};_.changeAtomCharge=function Iu(a,b){return Sh(this.a,a,b);};_.changeBond=function Ju(a,b){return Th(this.a,a,b);};_.convertStereoBondsToSingleBonds=function Ku(a){Ik(this.a,a);};_.copyAtom=function Lu(a,b,c,d){return Vh(this.a,a.a,b,c,d);};_.copyBond=function Mu(a,b,c,d,e,f){return Wh(this.a,a.a,b,c,d,e,f);};_.copyMolecule=function Nu(a){Xh(this.a,a.a);};_.copyMoleculeByAtoms=function Ou(a,b,c,d){Jk(this.a,a.a,b,c,d);};_.copyMoleculeByBonds=function Pu(a,b,c,d){return Kk(this.a,a.a,b,c,d);};_.copyMoleculeProperties=function Qu(a){Wo(this.a,a.a);};_.deleteAtom=function Ru(a){Zh(this.a,a);};_.deleteAtomOrBond=function Su(a,b){return $h(this.a,a,b);};_.deleteAtoms=function Tu(a){return _h(this.a,a);};_.deleteBond=function Uu(a){ai(this.a,a);};_.deleteBondAndSurrounding=function Vu(a){bi(this.a,a);};_.deleteMarkedAtomsAndBonds=function Wu(){return ci(this.a);};_.deleteMolecule=function Xu(){di(this.a);};_.deleteSelectedAtoms=function Yu(){return ei(this.a);};_.ensureHelperArrays=function Zu(a){Xo(this.a,a);};_.findAlleneCenterAtom=function $u(a){return Mk(this.a,a);};_.findAtom=function _u(a,b){return fi(this.a,a,b);};_.findBINAPChiralityBond=function av(a){return Nk(this.a,a);};_.findBond=function bv(a,b){return gi(this.a,a,b);};_.findRingSystem=function cv(a,b,c,d){Ok(this.a,a,b,c,d);};_.getAbsoluteAtomParity=function gv(a){return Yo(this.a,a);};_.getAbsoluteBondParity=function hv(a){return Zo(this.a,a);};_.getAllAtoms=function iv(){return this.a.o;};_.getAllBonds=function jv(){return this.a.p;};_.getAllConnAtoms=function kv(a){return Qk(this.a,a);};_.getAllHydrogens=function lv(a){return Rk(this.a,a);};_.getAromaticRingCount=function ov(){return Sk(this.a);};_.getAtomAbnormalValence=function pv(a){return hi(this.a,a);};_.getAtomCIPParity=function qv(a){return ii(this.a,a);};_.getAtomCharge=function rv(a){return ji(this.a,a);};_.getAtomColor=function sv(a){return ki(this.a,a);};_.getAtomCustomLabel=function tv(a){return li(this.a,a);};_.getAtomESRGroup=function uv(a){return ni(this.a,a);};_.getAtomESRType=function vv(a){return oi(this.a,a);};_.getAtomLabel=function wv(a){return pi(this.a,a);};_.getAtomList=function xv(a){return qi(this.a,a);};_.getAtomListString=function yv(a){return ri(this.a,a);};_.getAtomMapNo=function zv(a){return si(this.a,a);};_.getAtomMass=function Av(a){return ti(this.a,a);};_.getAtomParity=function Bv(a){return ui(this.a,a);};_.getAtomPi=function Cv(a){return Tk(this.a,a);};_.getAtomPreferredStereoBond=function Dv(a){return Uk(this.a,a);};_.getAtomQueryFeatures=function Ev(a){return vi(this.a,a);};_.getAtomRadical=function Fv(a){return wi(this.a,a);};_.getAtomRingBondCount=function Gv(a){return Vk(this.a,a);};_.getAtomRingCount=function Hv(a,b){return Wk(this.a,a,b);};_.getAtomRingSize=function Iv(a){return Xk(this.a,a);};_.getAtomX=function Jv(a){return xi(this.a,a);};_.getAtomY=function Kv(a){return yi(this.a,a);};_.getAtomZ=function Lv(a){return zi(this.a,a);};_.getAtomicNo=function Mv(a){return Ai(this.a,a);};_.getAtoms=function Ov(){return this.a.d;};_.getAverageBondLength=function Pv(a){return Yk(this.a,a);};_.getAverageTopologicalAtomDistance=function Qv(){return Zk(this.a);};_.getBond=function Rv(a,b){return $k(this.a,a,b);};_.getBondAngle=function Sv(a,b){return Di(this.a,a,b);};_.getBondAtom=function Tv(a,b){return Ei(this.a,a,b);};_.getBondBridgeMaxSize=function Uv(a){return Fi(this.a,a);};_.getBondBridgeMinSize=function Vv(a){return Gi(this.a,a);};_.getBondCIPParity=function Wv(a){return Hi(this.a,a);};_.getBondESRGroup=function Xv(a){return Ii(this.a,a);};_.getBondESRType=function Yv(a){return Ji(this.a,a);};_.getBondLength=function Zv(a){return Ki(this.a,a);};_.getBondOrder=function $v(a){return Mi(this.a,a);};_.getBondParity=function _v(a){return Ni(this.a,a);};_.getBondPreferredStereoBond=function aw(a){return Sl(this.a,a);};_.getBondQueryFeatures=function bw(a){return Oi(this.a,a);};_.getBondRingSize=function cw(a){return _k(this.a,a);};_.getBondType=function dw(a){return Pi(this.a,a);};_.getBondTypeSimple=function ew(a){return Qi(this.a,a);};_.getBonds=function fw(){return this.a.e;};_.getChiralText=function gw(){return $o(this.a);};_.getChirality=function hw(){return this.a.G;};_.getCompactCopy=function iw(){return new ru(_o(this.a));};_.getConnAtom=function jw(a,b){return al(this.a,a,b);};_.getConnAtoms=function kw(a){return bl(this.a,a);};_.getConnBond=function lw(a,b){return cl(this.a,a,b);};_.getConnBondOrder=function mw(a,b){return dl(this.a,a,b);};_.getDefaultMaxValenceUncharged=function nw(a){return Ri(this.a,a);};_.getDiastereotopicAtomIDs=function ow(){return sp(this.a);};_.getElectronValenceCorrection=function pw(a){return Si(this.a,a);};_.getFisherProjectionParity=function qw(a,b,c,d){return gl(this.a,a,b,c,d);};_.getFragmentAtoms=function rw(a){return hl(this.a,a);};_.getFragmentNumbers=function sw(a,b){return il(this.a,a,b);};_.getFragments=function tw(){var a,b,c;a=ap(this.a);c=ZB(aE,fR,39,a.length,0,1);for(b=0;b<a.length;b++){c[b]=new ru(a[b]);}return c;};_.getFreeValence=function uw(a){return jl(this.a,a);};_.getHandleHydrogenMap=function vw(){return kl(this.a);};_.getHoseCodes=function ww(a){a=a||{};var b=(typeof a.maxSphereSize===TT?5:a.maxSphereSize)|0;var c=(typeof a.type===TT?0:a.type)|0;return up(this.a,b,c);};_.getIDCode=function xw(){return cp(this.a);};_.getIDCodeAndCoordinates=function yw(){return{idCode:this.getIDCode(),coordinates:this.getIDCoordinates()};};_.getIDCoordinates=function zw(){return dp(this.a);};_.getImplicitHigherValence=function Aw(a,b){return ll(this.a,a,b);};_.getImplicitHydrogens=function Bw(a){return ml(this.a,a);};_.getIndex=function Cw(){return _n($z(nu),this.a);};_.getMaxAtoms=function Dw(){return this.a.K;};_.getMaxBonds=function Ew(){return this.a.L;};_.getMaxValence=function Fw(a){return Ti(this.a,a);};_.getMaxValenceUncharged=function Gw(a){return Ui(this.a,a);};_.getMolecularFormula=function Hw(){return new Hz(this.a);};_.getMolweight=function Iw(){return nl(this.a);};_.getName=function Jw(){return this.a.M;};_.getNumberOfHydrogens=function Kw(){return Lp(this.a);};_.getOccupiedValence=function Lw(a){return ol(this.a,a);};_.getPath=function Mw(a,b,c,d,e){return pl(this.a,a,b,c,d,e);};_.getPathBonds=function Nw(a,b,c){ql(this.a,a,b,c);};_.getPathLength=function Ow(a,b){return rl(this.a,a,b);};_.getProperties=function Pw(){return new Lz(this.a);};_.getRotatableBondCount=function Qw(){return ul(this.a);};_.kb=function Rw(a,b,c){var d;d=new zo(this.a,c);Dd(d,new tH(0,0,a,b));xd(d);return wo(d);};_.getSortedConnMap=function Sw(a){return vl(this.a,a);};_.getStereoBond=function Tw(a){return wl(this.a,a);};_.getStereoCenterCount=function Uw(){return ep(this.a);};_.getStereoProblem=function Vw(a){return Wi(this.a,a);};_.getSubstituent=function Ww(a,b,c,d,e){return yl(this.a,a,b,c,d.a,e);};_.getSubstituentSize=function Xw(a,b){return zl(this.a,a,b);};_.getSymmetryRank=function Yw(a){return fp(this.a,a);};_.invalidateHelperArrays=function Zw(a){Zi(this.a,a);};_.inventCoordinates=function $w(){pu(this);};_.isAllylicAtom=function _w(a){return Dl(this.a,a);};_.isAromaticAtom=function ax(a){return El(this.a,a);};_.isAromaticBond=function bx(a){return Fl(this.a,a);};_.isAtomConfigurationUnknown=function cx(a){return $i(this.a,a);};_.isAtomMarkedForDeletion=function dx(a){return _i(this.a,a);};_.isAtomParityPseudo=function ex(a){return aj(this.a,a);};_.isAtomStereoCenter=function fx(a){return bj(this.a,a);};_.isAutoMappedAtom=function ix(a){return cj(this.a,a);};_.isBINAPChiralityBond=function jx(a){return Gl(this.a,a);};_.isBondBackgroundHilited=function kx(a){return dj(this.a,a);};_.isBondBridge=function lx(a){return ej(this.a,a);};_.isBondForegroundHilited=function mx(a){return fj(this.a,a);};_.isBondMarkedForDeletion=function nx(a){return gj(this.a,a);};_.isBondParityPseudo=function ox(a){return hj(this.a,a);};_.isBondParityUnknownOrNone=function px(a){return ij(this.a,a);};_.isDelocalizedBond=function qx(a){return Hl(this.a,a);};_.isElectronegative=function rx(a){return jj(this.a,a);};_.isElectropositive=function sx(a){return kj(this.a,a);};_.isFlatNitrogen=function tx(a){return Il(this.a,a);};_.isFragment=function ux(){return this.a.I;};_.isMarkedAtom=function vx(a){return lj(this.a,a);};_.isMetalAtom=function wx(a){return mj(this.a,a);};_.isNaturalAbundance=function xx(a){return nj(this.a,a);};_.isOrganicAtom=function yx(a){return oj(this.a,a);};_.isPurelyOrganic=function zx(){return pj(this.a);};_.isRingAtom=function Ax(a){return Kl(this.a,a);};_.isRingBond=function Bx(a){return Ll(this.a,a);};_.isSelectedAtom=function Cx(a){return qj(this.a,a);};_.isSelectedBond=function Dx(a){return rj(this.a,a);};_.isSimpleHydrogen=function Ex(a){return Ml(this.a,a);};_.isSmallRingAtom=function Fx(a){return Nl(this.a,a);};_.isSmallRingBond=function Gx(a){return Ol(this.a,a);};_.isStabilizedAtom=function Hx(a){return Pl(this.a,a);};_.isStereoBond=function Ix(a){return sj(this.a,a);};_.markAtomForDeletion=function Jx(a){uj(this.a,a);};_.markBondForDeletion=function Kx(a){vj(this.a,a);};_.normalizeAmbiguousBonds=function Lx(){return Ql(this.a);};_.removeAtomColors=function Mx(){yj(this.a);};_.removeAtomCustomLabels=function Nx(){this.a.r=null;};_.removeAtomMarkers=function Ox(){zj(this.a);};_.removeAtomSelection=function Px(){Aj(this.a);};_.removeBondHiliting=function Qx(){Bj(this.a);};_.removeExplicitHydrogens=function Rx(){Ul(this.a);};_.renumberESRGroups=function Sx(a){return Dj(this.a,a);};_.scaleCoords=function Tx(a){Ej(this.a,a);};_.setAllAtoms=function Ux(a){Fj(this.a,a);};_.setAllBonds=function Vx(a){Gj(this.a,a);};_.setAssignParitiesToNitrogen=function Wx(a){gp(this.a,a);};_.setAtomAbnormalValence=function Xx(a,b){Hj(this.a,a,b);};_.setAtomCIPParity=function Yx(a,b){Ij(this.a,a,b);};_.setAtomCharge=function Zx(a,b){Jj(this.a,a,b);};_.setAtomColor=function $x(a,b){Kj(this.a,a,b);};_.setAtomConfigurationUnknown=function _x(a,b){Lj(this.a,a,b);};_.setAtomCustomLabel=function ay(a,b){Mj(this.a,a,b);};_.setAtomESR=function by(a,b,c){Oj(this.a,a,b,c);};_.setAtomList=function cy(a,b,c){Qj(this.a,a,b,c);};_.setAtomMapNo=function dy(a,b,c){Rj(this.a,a,b,c);};_.setAtomMarker=function ey(a,b){Sj(this.a,a,b);};_.setAtomMass=function fy(a,b){Tj(this.a,a,b);};_.setAtomParity=function gy(a,b,c){Uj(this.a,a,b,c);};_.setAtomQueryFeature=function hy(a,b,c){Vj(this.a,a,b,c);};_.setAtomRadical=function iy(a,b){Wj(this.a,a,b);};_.setAtomSelection=function jy(a,b){Xj(this.a,a,b);};_.setAtomX=function ky(a,b){Zj(this.a,a,b);};_.setAtomY=function ly(a,b){$j(this.a,a,b);};_.setAtomZ=function my(a,b){_j(this.a,a,b);};_.setAtomicNo=function ny(a,b){ak(this.a,a,b);};_.setBondAtom=function oy(a,b,c){bk(this.a,a,b,c);};_.setBondBackgroundHiliting=function py(a,b){ck(this.a,a,b);};_.setBondCIPParity=function qy(a,b){dk(this.a,a,b);};_.setBondESR=function ry(a,b,c){ek(this.a,a,b,c);};_.setBondForegroundHiliting=function sy(a,b){fk(this.a,a,b);};_.setBondParity=function ty(a,b,c){gk(this.a,a,b,c);};_.setBondParityUnknownOrNone=function uy(a){hk(this.a,a);};_.setBondQueryFeature=function vy(a,b,c){ik(this.a,a,b,c);};_.setBondType=function wy(a,b){jk(this.a,a,b);};_.setChirality=function xy(a){kk(this.a,a);};_.setFragment=function yy(a){lk(this.a,a);};_.setHydrogenProtection=function zy(a){mk(this.a,a);};_.setMaxAtoms=function Ay(a){nk(this.a,a);};_.setMaxBonds=function By(a){ok(this.a,a);};_.setName=function Cy(a){pk(this.a,a);};_.setParitiesValid=function Dy(a){Xl(this.a,a);};_.setStereoBondFromAtomParity=function Ey(a){Yl(this.a,a);};_.setStereoBondFromBondParity=function Fy(a){Zl(this.a,a);};_.setStereoBondsFromParity=function Gy(){$l(this.a);};_.setToRacemate=function Hy(){this.a.J=true;};_.setUnknownParitiesToExplicitlyUnknown=function Iy(){hp(this.a);};_.shareSameFragment=function Jy(a,b){return rl(this.a,a,b)!=-1;};_.stripIsotopInfo=function Ky(){return sk(this.a);};_.stripSmallFragments=function Ly(){return _l(this.a);};_.stripStereoInformation=function My(){tk(this.a);};_.supportsImplicitHydrogen=function Ny(a){return am(this.a,a);};_.toMolfile=function Oy(){var a;a=new Hm(this.a);return a.b.a;};_.toSVG=function Py(a,b,c){if(!$doc.createElement){throw new Error("Molecule#toSVG cannot be used outside of a browser's Window environment");}return this.kb(a,b,c);};_.toSmiles=function Qy(){return Bo(_z(nu),this.a);};_.translateCoords=function Ry(a,b){vk(this.a,a,b);};_.validate=function Sy(){ip(this.a);};_.zoomAndRotate=function Ty(a,b,c){xk(this.a,a,b,c);};_.zoomAndRotateInit=function Uy(a,b){yk(this.a,a,b);};var ir=UT,jr=UT,kr,lr=TR,mr=SR,nr=ST,or=0,pr=3,qr=1,rr=2,sr=0,tr=64,ur=384,vr=448,wr=192,xr=256,yr=320,zr=128,Ar,Br=1,Cr=2,Dr=4,Er=0,Fr=3,Gr=1,Hr=6,Ir=2,Jr=1,Kr=2,Lr=PQ,Mr=3,Nr=25,Or=XQ,Pr=1920,Qr=4,Rr=7,Sr=iR,Tr=xQ,Ur=PR,Vr=5,Wr=17,Xr=TQ,Yr=YQ,Zr=29,$r=128,_r=qR,as=yQ,bs=256,cs=cR,ds=32768,es=512,fs=tQ,gs=zQ,hs=16,is=OQ,js=1048576,ks=32,ls=UQ,ms=64,ns=4,os=8,ps=ER,qs=33554432,rs=FR,ss=3,ts=14,us=SQ,vs=WQ,ws=3,xs=22,ys=120,zs=4,As=3,Bs=PR,Cs=48,Ds=32,Es=0,Fs=16,Gs=4,Hs=48,Is=1,Js=0,Ks=3,Ls=2,Ms=1,Ns=0,Os=3,Ps=2,Qs=VT,Rs=bR,Ss=2,Ts=18,Us=cR,Vs=15,Ws=4,Xs=0,Ys=$Q,Zs=8,$s=960,_s=4,at=6,bt=6,ct=15360,dt=4,et=10,ft=8,gt=2,ht=qR,it=786480,jt=20,kt=tQ,lt=16,mt=32,nt=SQ,ot=3,pt=14,qt=48,rt=2,st=4,tt=786495,ut=1,vt=4,wt=26,xt=128,yt=64,zt=2,At=9,Bt=127,Ct=32,Dt=1,Et=4,Ft=17,Gt=458752,Ht=VQ,It=AQ,Jt=cR,Kt=qR,Lt=zQ,Mt=196608,Nt=0,Ot=327680,Pt=6,Qt=5,Rt=32,St=0,Tt=1,Ut=2,Vt=8,Wt=128,Xt=1,Yt=4,Zt=2,$t=32,_t=64,au=16,bu=252,cu=15,du=1,eu=0,fu=7,gu=3,hu=47,iu=79,ju=31,ku=190,lu=16,mu,nu;var aE=fI(39);DG(171,1,{},Vy);_.getField=function Wy(a){var b,c;c=Sp(this.a);for(b=0;b<c.length;b++){if(DJ(c[b],a)){return Rp(this.a,b);}}return null;};_.getFieldData=function Xy(a){return Rp(this.a,a);};_.getFieldNames=function Yy(a){return Tp(this.a,a);};_.getMolecule=function Zy(){return new ru(Up(this.a));};_.getNextFieldData=function $y(){var a;return a=this.a.a.a,a;};_.getNextMolFile=function _y(){var a;return a=this.a.f.a,a;};_.next=function az(){return Op(this.a);};var bE=fI(171);DG(173,1,{},bz);_.isFragmentInMolecule=function cz(){return Jn(this.a);};_.setFragment=function dz(a){On(this.a,a.a);};_.setMol=function ez(a,b){Pn(this.a,b.a);On(this.a,a.a);};_.setMolecule=function fz(a){Pn(this.a,a.a);};var dE=fI(173);DG(174,1,{},gz);_.createIndex=function iz(a){return _n(this.a,a.a);};_.isFragmentInMolecule=function oz(){return bo(this.a);};_.setFragment=function pz(a,b){eo(this.a,a.a,b);};_.setMolecule=function qz(a,b){fo(this.a,a.a,b);};var cE=fI(174);DG(175,1,{},Cz);_.assessRisk=function Dz(a,b){return Wq(a.a,(bA(Az),b));};_.getDetail=function Ez(a,b){return eA(Xq(a.a,b));};var rz=3,sz=2,tz,uz=1,vz=0,wz=2,xz=0,yz=3,zz=1,Az;var eE=fI(175);DG(176,1,{},Fz);var fE=fI(176);DG(124,126,{},Hz);EG(_,{absoluteWeight:{'get':function Iz(){return Bm(this);}}});EG(_,{formula:{'get':function Jz(){return Cm(this);}}});EG(_,{relativeWeight:{'get':function Kz(){return Dm(this);}}});var gE=fI(124);DG(125,127,{},Lz);EG(_,{acceptorCount:{'get':function Mz(){return an(this);}}});EG(_,{donorCount:{'get':function Nz(){return bn(this);}}});EG(_,{logP:{'get':function Oz(){return cn(this);}}});EG(_,{logPString:{'get':function Pz(){return eA(_p((new aq(),this.a)));}}});EG(_,{logS:{'get':function Qz(){return Bq((Aq(),this.a));}}});EG(_,{logSString:{'get':function Rz(){return eA(Cq((Aq(),this.a)));}}});EG(_,{polarSurfaceArea:{'get':function Sz(){return tq((rq(),this.a));}}});EG(_,{polarSurfaceAreaString:{'get':function Tz(){return eA(uq((rq(),this.a)));}}});EG(_,{rotatableBondCount:{'get':function Uz(){return ul(this.a);}}});EG(_,{stereoCenterCount:{'get':function Vz(){return ep(this.a);}}});var hE=fI(125);DG(81,1,{},cA);_.a=null;_.b=null;_.c=null;_.d=null;_.e=null;_.f=null;_.g=null;_.i=null;var Wz=null;var jE=fI(81);DG(130,1,{},dA);var iE=fI(130);var hA;DG(72,1,{4:1},nA);_.eb=function oA(a,b){return mA(a,b);};_.ab=function pA(a){return this===a;};var kE=fI(72);DG(119,1,{},sA);var lE=fI(119);DG(141,61,wQ);var pE=fI(141);DG(75,141,{75:1,4:1,12:1,14:1},MA);_.mb=function NA(){LA(this);return this.c;};_.ob=function OA(){return VC(this.b)===VC(JA)?null:this.b;};var JA;var mE=fI(75);DG(165,1,{});var oE=fI(165);var QA=0,RA=0,SA=-1;DG(158,165,{},eB);var aB;var qE=fI(158);DG(161,1,{},AB);var xB;var vE=fI(161);DG(120,1,{},RB);_.b=0;_.c=false;_.d=0;_.e=0;_.f=3;_.g=false;_.i=3;_.j=40;_.k=0;_.n=0;_.o=1;_.p=1;_.q='-';_.r='';_.t='';_.u='';_.v=false;var wE=fI(120);DG(163,1,{},TB);var xE=fI(163);var cC;var FC,GC,HC,IC;DG(79,14,vQ);var UE=fI(79);DG(25,79,vQ);var PE=fI(25);DG(121,25,vQ,NG);var yE=fI(121);DG(16,1,{16:1},VG,WG,XG);_.ab=function $G(a){return PC(a,16)&&a.c==this.c;};_.cb=function _G(){return this.c;};_.db=function bH(){return aI(zE),zE.k+'[r='+(this.c>>16&255)+',g='+(this.c>>8&255)+',b='+(this.c&255)+']';};_.a=0;_.b=null;_.c=0;var PG,QG,RG,SG;var zE=fI(16);DG(118,1,{},eH);_.b=0;var cH=null;var AE=fI(118);DG(107,1,{107:1});_.ab=function fH(a){var b;if(PC(a,40)){b=a;return this.a==b.a&&this.b==b.b;}return this===a;};_.cb=function gH(){var a;a=AI(this.a);a=vG(a,mG(AI(this.b),31));return tG(a)^tG(qG(a,32));};var CE=fI(107);DG(40,107,{107:1,40:1,4:1},hH,iH);_.db=function jH(){return'Point2D.Double['+this.a+', '+this.b+']';};_.a=0;_.b=0;var BE=fI(40);DG(188,1,{});var FE=fI(188);DG(108,188,{108:1});_.ab=function oH(a){var b;if(a===this){return true;}if(PC(a,19)){b=a;return this.c==b.c&&this.d==b.d&&this.b==b.b&&this.a==b.a;}return false;};_.cb=function pH(){var a;a=AI(this.c);a=eG(a,mG(AI(this.d),37));a=eG(a,mG(AI(this.b),43));a=eG(a,mG(AI(this.a),47));return tG(a)^tG(qG(a,32));};var EE=fI(108);DG(19,108,{108:1,19:1},sH,tH);_.db=function uH(){return aI(DE),DE.k+'[x='+this.c+',y='+this.d+',w='+this.b+',h='+this.a+']';};_.a=0;_.b=0;_.c=0;_.d=0;var DE=fI(19);DG(189,1,{});var KE=fI(189);DG(55,189,{},xH);_.a=0;var GE=fI(55);DG(168,1,{});var IE=fI(168);DG(167,168,{});var HE=fI(167);DG(122,167,{},yH);var JE=fI(122);DG(54,189,{},AH);_.a=0;var LE=fI(54);DG(71,1,{94:1});_.db=function GH(){return this.a;};var ME=fI(71);DG(152,26,wQ,HH);var NE=fI(152);DG(48,26,wQ,IH,JH);var YE=fI(48);DG(160,48,wQ,KH);var OE=fI(160);var BI,CI;DG(53,1,{4:1,30:1,53:1});_.fb=function FI(a){return this.b-a.b;};_.compareTo=function EI(a){return this.b-a.b;};_.equals=function GI(a){return this===a;};_.ab=function(a){return this.equals(a);};_.hashCode=function HI(){return YP(this);};_.cb=function(){return this.hashCode();};_.name=function II(){return this.a!=null?this.a:''+this.b;};_.ordinal=function JI(){return this.b;};_.toString=function KI(){return this.a!=null?this.a:''+this.b;};_.db=function(){return this.toString();};_.b=0;var TE=fI(53);DG(21,26,wQ,MI,NI);var WE=fI(21);DG(142,26,wQ,OI);var XE=fI(142);var ZI;DG(45,70,{4:1,30:1,45:1,70:1},aJ);_.fb=function cJ(a){return bJ(this.a,a.a);};_.ab=function dJ(a){return PC(a,45)&&iG(a.a,this.a);};_.cb=function eJ(){return tG(this.a);};_.db=function gJ(){return''+uG(this.a);};_.a=0;var _E=fI(45);var iJ;DG(159,26,wQ,oJ);var aF=fI(159);DG(63,21,wQ,tJ);var cF=fI(63);DG(43,1,{4:1,43:1},uJ);_.ab=function vJ(a){var b;if(PC(a,43)){b=a;return this.c==b.c&&IN(this.d,b.d)&&IN(this.a,b.a)&&IN(this.b,b.b);}return false;};_.cb=function wJ(){return rN(aC(VB(eF,1),fR,1,5,[YI(this.c),this.a,this.d,this.b]));};_.db=function xJ(){return this.a+'.'+this.d+'('+(this.b!=null?this.b:'Unknown Source')+(this.c>=0?':'+this.c:'')+')';};_.c=0;var gF=fI(43);DG(97,71,{94:1},SJ);var hF=fI(97);DG(38,71,{94:1},ZJ,$J,_J);var iF=fI(38);DG(143,48,wQ,aK);var jF=fI(143);var bK;DG(50,26,wQ,eK,fK);var mF=fI(50);DG(51,1,fU);_.fb=function gK(a){return AJ(this.a,a.a);};_.ab=function hK(a){var b;if(a===this){return true;}if(!PC(a,51)){return false;}b=a;return DJ(this.a,b.a);};_.cb=function iK(){return cQ(this.a);};_.db=function jK(){return this.a;};var nF=fI(51);DG(192,1,{});var pF=fI(192);DG(69,192,{},mK,nK);var oF=fI(69);DG(33,1,gU);_.add=function sK(a){throw dG(new fK('Add not supported on this collection'));};_.addAll=function tK(a){var b,c,d;OP(a);b=false;for(d=a.yb();d.Bb();){c=d.Cb();b=b|this.add(c);}return b;};_.clear=function uK(){var a;for(a=this.yb();a.Bb();){a.Cb();a.Db();}};_.contains=function vK(a){return oK(this,a,false);};_.containsAll=function wK(a){return pK(this,a);};_.isEmpty=function xK(){return this.size()==0;};_.remove=function yK(a){return oK(this,a,true);};_.removeAll=function zK(a){return qK(this,a);};_.retainAll=function AK(a){var b,c,d;OP(a);b=false;for(c=this.yb();c.Bb();){d=c.Cb();if(!a.contains(d)){c.Db();b=true;}}return b;};_.toArray=function BK(){return this.zb(ZB(eF,fR,1,this.size(),5,1));};_.zb=function CK(a){var b,c,d,e;e=this.size();a.length<e&&(a=(d=new Array(e),bC(d,a)));c=this.yb();for(b=0;b<e;++b){a[b]=c.Cb();}a.length>e&&(a[e]=null);return a;};_.db=function DK(){return rK(this);};var qF=fI(33);DG(74,33,hU);_.addAtIndex=function EK(a,b){throw dG(new fK('Add not supported on this list'));};_.add=function FK(a){this.addAtIndex(this.size(),a);return true;};_.addAllAtIndex=function GK(a,b){var c,d,e;OP(b);c=false;for(e=b.yb();e.Bb();){d=e.Cb();this.addAtIndex(a++,d);c=true;}return c;};_.clear=function HK(){this.Ab(0,this.size());};_.ab=function IK(a){var b,c,d,e,f;if(a===this){return true;}if(!PC(a,93)){return false;}f=a;if(this.size()!=f.size()){return false;}e=f.yb();for(c=this.yb();c.Bb();){b=c.Cb();d=e.Cb();if(!(VC(b)===VC(d)||b!=null&&pc(b,d))){return false;}}return true;};_.cb=function JK(){return BN(this);};_.indexOf=function KK(a){var b,c;for(b=0,c=this.size();b<c;++b){if(IN(a,this.getAtIndex(b))){return b;}}return-1;};_.yb=function LK(){return new UK(this);};_.lastIndexOf=function MK(a){var b;for(b=this.size()-1;b>-1;--b){if(IN(a,this.getAtIndex(b))){return b;}}return-1;};_.removeAtIndex=function NK(a){throw dG(new fK('Remove not supported on this list'));};_.Ab=function OK(a,b){var c,d;d=new YK(this,a);for(c=a;c<b;++c){MP(d.a<d.c.size());d.c.getAtIndex(d.b=d.a++);TK(d);}};_.setAtIndex=function PK(a,b){throw dG(new fK('Set not supported on this list'));};_.subList=function QK(a,b){return new $K(this,a,b);};var uF=fI(74);DG(111,1,{},UK);_.Bb=function VK(){return RK(this);};_.Cb=function WK(){return SK(this);};_.Db=function XK(){TK(this);};_.a=0;_.b=-1;var rF=fI(111);DG(112,111,{},YK);_.Db=function ZK(){TK(this);};var sF=fI(112);DG(113,74,hU,$K);_.addAtIndex=function _K(a,b){QP(a,this.b);this.c.addAtIndex(this.a+a,b);++this.b;};_.getAtIndex=function aL(a){NP(a,this.b);return this.c.getAtIndex(this.a+a);};_.removeAtIndex=function bL(a){var b;NP(a,this.b);b=this.c.removeAtIndex(this.a+a);--this.b;return b;};_.setAtIndex=function cL(a,b){NP(a,this.b);return this.c.setAtIndex(this.a+a,b);};_.size=function dL(){return this.b;};_.a=0;_.b=0;var tF=fI(113);DG(190,1,{164:1});_.getOrDefault=function nL(a,b){var c;return c=this.get(a),c==null&&!this.containsKey(a)?b:c;};_.putIfAbsent=function tL(a,b){var c;return c=this.get(a),c!=null?c:this.put(a,b);};_.replace=function vL(a,b){return this.containsKey(a)?this.put(a,b):null;};_.clear=function hL(){this.Eb().clear();};_.containsKey=function iL(a){return!!fL(this,a,false);};_.containsValue=function jL(a){return eL(this,a);};_.ab=function kL(a){var b,c,d;if(a===this){return true;}if(!PC(a,46)){return false;}d=a;if(this.c!=d.c){return false;}for(c=new qO(new vO(d).b);RK(c.a);){b=c.b=SK(c.a);if(!aM(this,b)){return false;}}return true;};_.get=function lL(a){return mL(fL(this,a,false));};_.cb=function oL(){return AN(this.Eb());};_.isEmpty=function pL(){return this.c==0;};_.keySet=function qL(){return new CL(this);};_.put=function rL(a,b){throw dG(new fK('Put not supported on this map'));};_.putAll=function sL(a){var b,c;OP(a);for(c=new qO(a.Eb().b);RK(c.a);){b=c.b=SK(c.a);aO(this,b.Fb(),b.Gb());}};_.remove=function uL(a){return mL(fL(this,a,true));};_.size=function wL(){return this.Eb().b.c;};_.db=function xL(){var a,b,c;c=new UN('{','}');for(b=new qO(this.Eb().b);RK(b.a);){a=b.b=SK(b.a);TN(c,gL(this,a.Fb())+'='+gL(this,a.Gb()));}return!c.a?c.c:c.e.length==0?c.a.a:c.a.a+(''+c.e);};_.values=function yL(){return new ML(this);};var BF=fI(190);DG(166,33,iU);_.ab=function zL(a){var b;if(a===this){return true;}if(!PC(a,62)){return false;}b=a;if(b.size()!=this.size()){return false;}return pK(this,b);};_.cb=function AL(){return AN(this);};_.removeAll=function BL(a){var b,c,d,e;OP(a);e=this.size();if(e<a.size()){for(b=this.yb();b.Bb();){c=b.Cb();a.contains(c)&&b.Db();}}else{for(d=a.yb();d.Bb();){c=d.Cb();this.remove(c);}}return e!=this.size();};var GF=fI(166);DG(149,166,iU,CL);_.clear=function DL(){WN(this.a);};_.contains=function EL(a){return bM(this.a,a);};_.yb=function FL(){var a;a=new qO(new vO(this.a).b);return new IL(a);};_.remove=function GL(a){if(bM(this.a,a)){bO(this.a,a);return true;}return false;};_.size=function HL(){return this.a.c;};var wF=fI(149);DG(150,1,{},IL);_.Bb=function JL(){return RK(this.a.a);};_.Cb=function KL(){var a;a=oO(this.a);return a.Fb();};_.Db=function LL(){pO(this.a);};var vF=fI(150);DG(116,33,gU,ML);_.clear=function NL(){WN(this.a);};_.contains=function OL(a){return eL(this.a,a);};_.yb=function PL(){var a;return a=new qO(new vO(this.a).b),new RL(a);};_.size=function QL(){return this.a.c;};var yF=fI(116);DG(117,1,{},RL);_.Bb=function SL(){return RK(this.a.a);};_.Cb=function TL(){var a;return a=oO(this.a),a.Gb();};_.Db=function UL(){pO(this.a);};var xF=fI(117);DG(76,1,{76:1,78:1});_.ab=function WL(a){var b;if(!PC(a,78)){return false;}b=a;return IN(this.c,b.Fb())&&IN(this.d,b.Gb());};_.Fb=function XL(){return this.c;};_.Gb=function YL(){return this.d;};_.cb=function ZL(){return JN(this.c)^JN(this.d);};_.db=function $L(){return this.c+'='+this.d;};var zF=fI(76);DG(77,76,{76:1,77:1,78:1},_L);var AF=fI(77);DG(191,190,{164:1});_.containsKey=function dM(a){return bM(this,a);};_.Eb=function eM(){return new hM(this);};_.get=function fM(a){return cM(this,a);};_.keySet=function gM(){return new mM(this);};var FF=fI(191);DG(115,166,iU,hM);_.contains=function iM(a){return PC(a,78)&&aM(this.b,a);};_.yb=function jM(){return new qO(this.b);};_.remove=function kM(a){var b;if(PC(a,78)){b=a;return cO(this.b,b);}return false;};_.size=function lM(){return this.b.c;};var CF=fI(115);DG(91,166,iU,mM);_.clear=function nM(){WN(this.a);};_.contains=function oM(a){return bM(this.a,a);};_.yb=function pM(){var a;return a=new qO(new vO(this.a).b),new sM(a);};_.remove=function qM(a){if(bM(this.a,a)){bO(this.a,a);return true;}return false;};_.size=function rM(){return this.a.c;};var EF=fI(91);DG(92,1,{},sM);_.Bb=function uM(){return RK(this.a.a);};_.Cb=function vM(){var a;return a=oO(this.a),a.Fb();};_.Db=function wM(){pO(this.a);};var DF=fI(92);DG(17,74,{4:1,5:1,37:1,33:1,74:1,17:1,36:1,93:1,193:1},LM);_.addAtIndex=function MM(a,b){xM(this,a,b);};_.add=function NM(a){return yM(this,a);};_.addAllAtIndex=function OM(a,b){return zM(this,a,b);};_.addAll=function PM(a){return AM(this,a);};_.clear=function QM(){this.a=ZB(eF,fR,1,0,5,1);};_.contains=function RM(a){return CM(this,a,0)!=-1;};_.getAtIndex=function SM(a){return BM(this,a);};_.indexOf=function TM(a){return CM(this,a,0);};_.isEmpty=function UM(){return this.a.length==0;};_.yb=function VM(){return new dN(this);};_.lastIndexOf=function WM(a){return DM(this,a);};_.removeAtIndex=function XM(a){return FM(this,a);};_.remove=function YM(a){return GM(this,a);};_.Ab=function ZM(a,b){HM(this,a,b);};_.setAtIndex=function $M(a,b){return IM(this,a,b);};_.size=function _M(){return this.a.length;};_.toArray=function aN(){return JM(this);};_.zb=function bN(a){return KM(this,a);};var IF=fI(17);DG(44,1,{},dN);_.Bb=function eN(){return this.a<this.c.a.length;};_.Cb=function fN(){return cN(this);};_.Db=function gN(){SP(this.b!=-1);FM(this.c,this.a=this.b);this.b=-1;};_.a=0;_.b=-1;var HF=fI(44);var CN;DG(157,1,{4:1},EN);_.eb=function FN(a,b){return OP(a),PH(a,(OP(b),b));};_.ab=function GN(a){return this===a;};var JF=fI(157);DG(162,26,wQ,HN);var KF=fI(162);DG(83,1,{},RN,SN);_.a=0;_.b=0;var KN,LN,MN=0;var LF=fI(83);DG(95,1,{},UN);_.db=function VN(){return!this.a?this.c:this.e.length==0?this.a.a:this.a.a+(''+this.e);};var MF=fI(95);DG(46,191,{4:1,164:1,46:1},hO,iO);_.clear=function jO(){WN(this);};_.Eb=function kO(){return new vO(this);};_.put=function lO(a,b){return aO(this,a,b);};_.remove=function mO(a){return bO(this,a);};_.size=function nO(){return this.c;};_.c=0;var VF=fI(46);DG(32,1,{},qO);_.Cb=function tO(){return oO(this);};_.Bb=function sO(){return RK(this.a);};_.Db=function uO(){pO(this);};var NF=fI(32);DG(41,115,iU,vO);_.clear=function wO(){WN(this.a);};var OF=fI(41);DG(60,77,{76:1,77:1,78:1,60:1},xO);_.b=false;var PF=fI(60);DG(90,1,{},yO);_.db=function zO(){return'State: mv='+this.c+' value='+this.d+' done='+this.a+' found='+this.b;};_.a=false;_.b=false;_.c=false;var QF=fI(90);DG(42,53,lU,FO);_.Hb=function GO(){return false;};_.Ib=function HO(){return false;};var AO,BO,CO,DO;var UF=gI(42,IO);DG(146,42,lU,JO);_.Ib=function KO(){return true;};var RF=gI(146,null);DG(147,42,lU,LO);_.Hb=function MO(){return true;};_.Ib=function NO(){return true;};var SF=gI(147,null);DG(148,42,lU,OO);_.Hb=function PO(){return true;};var TF=gI(148,null);DG(73,166,{4:1,37:1,33:1,36:1,62:1},SO);_.add=function TO(a){return QO(this,a);};_.clear=function UO(){WN(this.a);};_.contains=function VO(a){return RO(this,a);};_.yb=function WO(){var a;return a=new qO(new vO(new mM(this.a).a).b),new sM(a);};_.remove=function XO(a){return bO(this.a,a)!=null;};_.size=function YO(){return this.a.c;};var WF=fI(73);DG(151,74,{4:1,5:1,37:1,33:1,74:1,36:1,93:1,193:1},_O);_.addAtIndex=function aP(a,b){eP(a,this.a.a.length+1);xM(this.a,a,b);};_.add=function bP(a){return ZO(this,a);};_.addAllAtIndex=function cP(a,b){eP(a,this.a.a.length+1);return zM(this.a,a,b);};_.addAll=function dP(a){return AM(this.a,a);};_.clear=function fP(){this.a.a=ZB(eF,fR,1,0,5,1);};_.contains=function gP(a){return CM(this.a,a,0)!=-1;};_.containsAll=function hP(a){return pK(this.a,a);};_.getAtIndex=function iP(a){eP(a,this.a.a.length);return BM(this.a,a);};_.indexOf=function jP(a){return CM(this.a,a,0);};_.isEmpty=function kP(){return this.a.a.length==0;};_.yb=function lP(){return new dN(this.a);};_.lastIndexOf=function mP(a){return DM(this.a,a);};_.removeAtIndex=function nP(a){eP(a,this.a.a.length);return FM(this.a,a);};_.removeAll=function oP(a){return qK(this.a,a);};_.Ab=function pP(a,b){HM(this.a,a,b);};_.setAtIndex=function qP(a,b){eP(a,this.a.a.length);return IM(this.a,a,b);};_.size=function rP(){return this.a.a.length;};_.subList=function sP(a,b){return new $K(this.a,a,b);};_.toArray=function tP(){return JM(this.a);};_.zb=function uP(a){return $O(this,a);};_.db=function vP(){return rK(this.a);};var XF=fI(151);DG(104,51,fU);var $F=fI(104);DG(105,104,fU,FP);var YF=fI(105);DG(137,104,fU,JP);var ZF=fI(137);var YC=hI('C');var _C=hI('I');var aD=hI('J');var _F=hI('S');var aG=hI('Z');var ZC=hI('D');var $C=hI('F');var XC=hI('B');_=HG('OCL.DrugScoreCalculator',_q);_.calculate=ar;dr();_=HG('OCL.DruglikenessPredictor',er);_.DRUGLIKENESS_UNKNOWN=br;ou();_=HG('OCL.Molecule',ru);_.FISCHER_PROJECTION_LIMIT=ir;_.STEREO_ANGLE_LIMIT=jr;_.VALIDATION_ERRORS_STEREO=kr;_.VALIDATION_ERROR_AMBIGUOUS_CONFIGURATION=lr;_.VALIDATION_ERROR_ESR_CENTER_UNKNOWN=mr;_.VALIDATION_ERROR_OVER_UNDER_SPECIFIED=nr;_.cAtomCIPParityNone=or;_.cAtomCIPParityProblem=pr;_.cAtomCIPParityRorM=qr;_.cAtomCIPParitySorP=rr;_.cAtomColorBlack=sr;_.cAtomColorBlue=tr;_.cAtomColorDarkGreen=ur;_.cAtomColorDarkRed=vr;_.cAtomColorGreen=wr;_.cAtomColorMagenta=xr;_.cAtomColorOrange=yr;_.cAtomColorRed=zr;_.cAtomLabel=Ar;_.cAtomParity1=Br;_.cAtomParity2=Cr;_.cAtomParityIsPseudo=Dr;_.cAtomParityNone=Er;_.cAtomParityUnknown=Fr;_.cAtomQFAny=Gr;_.cAtomQFAromState=Hr;_.cAtomQFAromStateBits=Ir;_.cAtomQFAromStateShift=Jr;_.cAtomQFAromatic=Kr;_.cAtomQFCharge=Lr;_.cAtomQFChargeBits=Mr;_.cAtomQFChargeShift=Nr;_.cAtomQFFlatNitrogen=Or;_.cAtomQFHydrogen=Pr;_.cAtomQFHydrogenBits=Qr;_.cAtomQFHydrogenShift=Rr;_.cAtomQFMatchStereo=Sr;_.cAtomQFMoreNeighbours=Tr;_.cAtomQFNarrowing=Ur;_.cAtomQFNeighbourBits=Vr;_.cAtomQFNeighbourShift=Wr;_.cAtomQFNeighbours=Xr;_.cAtomQFNoMoreNeighbours=Yr;_.cAtomQFNoOfBits=Zr;_.cAtomQFNot0Hydrogen=$r;_.cAtomQFNot0Neighbours=_r;_.cAtomQFNot0PiElectrons=as;_.cAtomQFNot1Hydrogen=bs;_.cAtomQFNot1Neighbour=cs;_.cAtomQFNot1PiElectron=ds;_.cAtomQFNot2Hydrogen=es;_.cAtomQFNot2Neighbours=fs;_.cAtomQFNot2PiElectrons=gs;_.cAtomQFNot2RingBonds=hs;_.cAtomQFNot3Hydrogen=is;_.cAtomQFNot3Neighbours=js;_.cAtomQFNot3RingBonds=ks;_.cAtomQFNot4Neighbours=ls;_.cAtomQFNot4RingBonds=ms;_.cAtomQFNotAromatic=ns;_.cAtomQFNotChain=os;_.cAtomQFNotCharge0=ps;_.cAtomQFNotChargeNeg=qs;_.cAtomQFNotChargePos=rs;_.cAtomQFPiElectronBits=ss;_.cAtomQFPiElectronShift=ts;_.cAtomQFPiElectrons=us;_.cAtomQFRingSize=vs;_.cAtomQFRingSizeBits=ws;_.cAtomQFRingSizeShift=xs;_.cAtomQFRingState=ys;_.cAtomQFRingStateBits=zs;_.cAtomQFRingStateShift=As;_.cAtomQFSimpleFeatures=Bs;_.cAtomRadicalState=Cs;_.cAtomRadicalStateD=Ds;_.cAtomRadicalStateNone=Es;_.cAtomRadicalStateS=Fs;_.cAtomRadicalStateShift=Gs;_.cAtomRadicalStateT=Hs;_.cBondCIPParityEorP=Is;_.cBondCIPParityNone=Js;_.cBondCIPParityProblem=Ks;_.cBondCIPParityZorM=Ls;_.cBondParityEor1=Ms;_.cBondParityNone=Ns;_.cBondParityUnknown=Os;_.cBondParityZor2=Ps;_.cBondQFAllFeatures=Qs;_.cBondQFAromState=Rs;_.cBondQFAromStateBits=Ss;_.cBondQFAromStateShift=Ts;_.cBondQFAromatic=Us;_.cBondQFBondTypes=Vs;_.cBondQFBondTypesBits=Ws;_.cBondQFBondTypesShift=Xs;_.cBondQFBridge=Ys;_.cBondQFBridgeBits=Zs;_.cBondQFBridgeMin=$s;_.cBondQFBridgeMinBits=_s;_.cBondQFBridgeMinShift=at;_.cBondQFBridgeShift=bt;_.cBondQFBridgeSpan=ct;_.cBondQFBridgeSpanBits=dt;_.cBondQFBridgeSpanShift=et;_.cBondQFDelocalized=ft;_.cBondQFDouble=gt;_.cBondQFMatchStereo=ht;_.cBondQFNarrowing=it;_.cBondQFNoOfBits=jt;_.cBondQFNotAromatic=kt;_.cBondQFNotRing=lt;_.cBondQFRing=mt;_.cBondQFRingSize=nt;_.cBondQFRingSizeBits=ot;_.cBondQFRingSizeShift=pt;_.cBondQFRingState=qt;_.cBondQFRingStateBits=rt;_.cBondQFRingStateShift=st;_.cBondQFSimpleFeatures=tt;_.cBondQFSingle=ut;_.cBondQFTriple=vt;_.cBondTypeCross=wt;_.cBondTypeDeleted=xt;_.cBondTypeDelocalized=yt;_.cBondTypeDouble=zt;_.cBondTypeDown=At;_.cBondTypeIncreaseOrder=Bt;_.cBondTypeMetalLigand=Ct;_.cBondTypeSingle=Dt;_.cBondTypeTriple=Et;_.cBondTypeUp=Ft;_.cChiralityDiastereomers=Gt;_.cChiralityEpimers=Ht;_.cChiralityIsomerCountMask=It;_.cChiralityKnownEnantiomer=Jt;_.cChiralityMeso=Kt;_.cChiralityNotChiral=Lt;_.cChiralityRacemic=Mt;_.cChiralityUnknown=Nt;_.cChiralityUnknownEnantiomer=Ot;_.cDefaultAtomValence=Pt;_.cESRGroupBits=Qt;_.cESRMaxGroups=Rt;_.cESRTypeAbs=St;_.cESRTypeAnd=Tt;_.cESRTypeOr=Ut;_.cHelperBitCIP=Vt;_.cHelperBitIncludeNitrogenParities=Wt;_.cHelperBitNeighbours=Xt;_.cHelperBitParities=Yt;_.cHelperBitRings=Zt;_.cHelperBitSymmetryDiastereotopic=$t;_.cHelperBitSymmetryEnantiotopic=_t;_.cHelperBitSymmetrySimple=au;_.cHelperBitsStereo=bu;_.cHelperCIP=cu;_.cHelperNeighbours=du;_.cHelperNone=eu;_.cHelperParities=fu;_.cHelperRings=gu;_.cHelperSymmetryDiastereotopic=hu;_.cHelperSymmetryEnantiotopic=iu;_.cHelperSymmetrySimple=ju;_.cMaxAtomicNo=ku;_.cMaxConnAtoms=lu;_.cRoundedMass=mu;_.fromIDCode=dv;_.fromMolfile=ev;_.fromSmiles=fv;_.getAngle=mv;_.getAngleDif=nv;_.getAtomicNoFromLabel=Nv;_.isAtomicNoElectronegative=gx;_.isAtomicNoElectropositive=hx;_=HG('OCL.SDFileParser',Vy);_=HG('OCL.SSSearcher',bz);_=HG('OCL.SSSearcherWithIndex',gz);_.bitCount=hz;_.getHexStringFromIndex=jz;_.getIndexFromHexString=kz;_.getKeyIDCode=lz;_.getSimilarityAngleCosine=mz;_.getSimilarityTanimoto=nz;Bz();_=HG('OCL.ToxicityPredictor',Cz);_.RISK_HIGH=rz;_.RISK_LOW=sz;_.RISK_NAMES=tz;_.RISK_NO=uz;_.RISK_UNKNOWN=vz;_.TYPE_IRRITANT=wz;_.TYPE_MUTAGENIC=xz;_.TYPE_REPRODUCTIVE_EFFECTIVE=yz;_.TYPE_TUMORIGENIC=zz;_=HG('OCL.Util',Fz);_.getHoseCodesFromDiastereotopicID=Gz;MH();_=HG('java.lang.Boolean');_.$isInstance=NH;_=HG('java.lang.CharSequence');_.$isInstance=QH;_=HG('java.lang.Comparable');_.$isInstance=tI;_=HG('java.lang.Double');_.$isInstance=yI;_=HG('java.lang.Number');_.$isInstance=vI;_=HG('java.lang.String');_.$isInstance=HJ;_=HG('java.lang.Throwable');_.of=DA;var eQ=(TA(),WA);var gwtOnLoad=gwtOnLoad=zG;xG(JG);AG('permProps',[[['locale','default'],['user.agent',_T]]]);$sendStats('moduleStartup','moduleEvalEnd');gwtOnLoad(__gwtModuleFunction.__errFn,__gwtModuleFunction.__moduleName,__gwtModuleFunction.__moduleBase,__gwtModuleFunction.__softPermutationId,__gwtModuleFunction.__computePropValue);$sendStats('moduleStartup','end');$gwt&&$gwt.permProps&&__gwtModuleFunction.__moduleStartupDone($gwt.permProps);// End GWT code
var toReturn=$wnd["OCL"];toReturn.version='4.4.0';return toReturn;}var isBrowser,globalEnv,document;if(typeof window!=='undefined'){// usual browser window
isBrowser=true;globalEnv=window;document=window.document;}else if(typeof self!=='undefined'){// Web Worker
isBrowser=true;globalEnv=self;document={};}else if(typeof global!=='undefined'){// Node.js
isBrowser=false;globalEnv=global;document={};}else{// Other environment (example: CouchDB)
isBrowser=false;globalEnv=root;document={};}var fakeWindow;if(isBrowser&&!true){fakeWindow=globalEnv;}else{fakeWindow={};fakeWindow.setTimeout=globalEnv.setTimeout?globalEnv.setTimeout.bind(globalEnv):noop;fakeWindow.clearTimeout=globalEnv.clearTimeout?globalEnv.clearTimeout.bind(globalEnv):noop;fakeWindow.setInterval=globalEnv.setInterval?globalEnv.setInterval.bind(globalEnv):noop;fakeWindow.clearInterval=globalEnv.clearInterval?globalEnv.clearInterval.bind(globalEnv):noop;// required since GWT 2.8.0
fakeWindow.Error=globalEnv.Error;fakeWindow.Math=globalEnv.Math;fakeWindow.RegExp=globalEnv.RegExp;fakeWindow.TypeError=globalEnv.TypeError;}if(!fakeWindow.document){fakeWindow.document=document;}var exportedApi=getExports(fakeWindow);if(true){// NodeJS
fillExports(exportedApi,exports);}else { var i, obj, l, path; }function fillExports(obj,exports){var keys=Object.keys(obj);for(var i=0;i<keys.length;i++){exports[keys[i]]=obj[keys[i]];}}function noop(){}})(this);
/* WEBPACK VAR INJECTION */}.call(this, __webpack_require__(4)))

/***/ }),
/* 3 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


var Matrix = __webpack_require__(0);
/**
 * @private
 * Check that a row index is not out of bounds
 * @param {Matrix} matrix
 * @param {number} index
 * @param {boolean} [outer]
 */


exports.checkRowIndex = function checkRowIndex(matrix, index, outer) {
  var max = outer ? matrix.rows : matrix.rows - 1;

  if (index < 0 || index > max) {
    throw new RangeError('Row index out of range');
  }
};
/**
 * @private
 * Check that a column index is not out of bounds
 * @param {Matrix} matrix
 * @param {number} index
 * @param {boolean} [outer]
 */


exports.checkColumnIndex = function checkColumnIndex(matrix, index, outer) {
  var max = outer ? matrix.columns : matrix.columns - 1;

  if (index < 0 || index > max) {
    throw new RangeError('Column index out of range');
  }
};
/**
 * @private
 * Check that the provided vector is an array with the right length
 * @param {Matrix} matrix
 * @param {Array|Matrix} vector
 * @return {Array}
 * @throws {RangeError}
 */


exports.checkRowVector = function checkRowVector(matrix, vector) {
  if (vector.to1DArray) {
    vector = vector.to1DArray();
  }

  if (vector.length !== matrix.columns) {
    throw new RangeError('vector size must be the same as the number of columns');
  }

  return vector;
};
/**
 * @private
 * Check that the provided vector is an array with the right length
 * @param {Matrix} matrix
 * @param {Array|Matrix} vector
 * @return {Array}
 * @throws {RangeError}
 */


exports.checkColumnVector = function checkColumnVector(matrix, vector) {
  if (vector.to1DArray) {
    vector = vector.to1DArray();
  }

  if (vector.length !== matrix.rows) {
    throw new RangeError('vector size must be the same as the number of rows');
  }

  return vector;
};

exports.checkIndices = function checkIndices(matrix, rowIndices, columnIndices) {
  var rowOut = rowIndices.some(r => {
    return r < 0 || r >= matrix.rows;
  });
  var columnOut = columnIndices.some(c => {
    return c < 0 || c >= matrix.columns;
  });

  if (rowOut || columnOut) {
    throw new RangeError('Indices are out of range');
  }

  if (typeof rowIndices !== 'object' || typeof columnIndices !== 'object') {
    throw new TypeError('Unexpected type for row/column indices');
  }

  if (!Array.isArray(rowIndices)) rowIndices = Array.from(rowIndices);
  if (!Array.isArray(columnIndices)) rowIndices = Array.from(columnIndices);
  return {
    row: rowIndices,
    column: columnIndices
  };
};

exports.checkRange = function checkRange(matrix, startRow, endRow, startColumn, endColumn) {
  if (arguments.length !== 5) throw new TypeError('Invalid argument type');
  var notAllNumbers = Array.from(arguments).slice(1).some(function (arg) {
    return typeof arg !== 'number';
  });
  if (notAllNumbers) throw new TypeError('Invalid argument type');

  if (startRow > endRow || startColumn > endColumn || startRow < 0 || startRow >= matrix.rows || endRow < 0 || endRow >= matrix.rows || startColumn < 0 || startColumn >= matrix.columns || endColumn < 0 || endColumn >= matrix.columns) {
    throw new RangeError('Submatrix indices are out of range');
  }
};

exports.getRange = function getRange(from, to) {
  var arr = new Array(to - from + 1);

  for (var i = 0; i < arr.length; i++) {
    arr[i] = from + i;
  }

  return arr;
};

exports.sumByRow = function sumByRow(matrix) {
  var sum = Matrix.Matrix.zeros(matrix.rows, 1);

  for (var i = 0; i < matrix.rows; ++i) {
    for (var j = 0; j < matrix.columns; ++j) {
      sum.set(i, 0, sum.get(i, 0) + matrix.get(i, j));
    }
  }

  return sum;
};

exports.sumByColumn = function sumByColumn(matrix) {
  var sum = Matrix.Matrix.zeros(1, matrix.columns);

  for (var i = 0; i < matrix.rows; ++i) {
    for (var j = 0; j < matrix.columns; ++j) {
      sum.set(0, j, sum.get(0, j) + matrix.get(i, j));
    }
  }

  return sum;
};

exports.sumAll = function sumAll(matrix) {
  var v = 0;

  for (var i = 0; i < matrix.rows; i++) {
    for (var j = 0; j < matrix.columns; j++) {
      v += matrix.get(i, j);
    }
  }

  return v;
};

/***/ }),
/* 4 */
/***/ (function(module, exports) {

var g; // This works in non-strict mode

g = function () {
  return this;
}();

try {
  // This works if eval is allowed (see CSP)
  g = g || new Function("return this")();
} catch (e) {
  // This works if the window reference is available
  if (typeof window === "object") g = window;
} // g can still be undefined, but nothing to do about it...
// We return undefined, instead of nothing here, so it's
// easier to handle this case. if(!global) { ...}


module.exports = g;

/***/ }),
/* 5 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


exports.hypotenuse = function hypotenuse(a, b) {
  var r;

  if (Math.abs(a) > Math.abs(b)) {
    r = b / a;
    return Math.abs(a) * Math.sqrt(1 + r * r);
  }

  if (b !== 0) {
    r = a / b;
    return Math.abs(b) * Math.sqrt(1 + r * r);
  }

  return 0;
}; // For use in the decomposition algorithms. With big matrices, access time is
// too long on elements from array subclass
// todo check when it is fixed in v8
// http://jsperf.com/access-and-write-array-subclass


exports.getEmpty2DArray = function (rows, columns) {
  var array = new Array(rows);

  for (var i = 0; i < rows; i++) {
    array[i] = new Array(columns);
  }

  return array;
};

exports.getFilled2DArray = function (rows, columns, value) {
  var array = new Array(rows);

  for (var i = 0; i < rows; i++) {
    array[i] = new Array(columns);

    for (var j = 0; j < columns; j++) {
      array[i][j] = value;
    }
  }

  return array;
};

/***/ }),
/* 6 */
/***/ (function(module, exports, __webpack_require__) {

/* WEBPACK VAR INJECTION */(function(global, process) {(function (global, undefined) {
  "use strict";

  if (global.setImmediate) {
    return;
  }

  var nextHandle = 1; // Spec says greater than zero

  var tasksByHandle = {};
  var currentlyRunningATask = false;
  var doc = global.document;
  var setImmediate;

  function addFromSetImmediateArguments(args) {
    tasksByHandle[nextHandle] = partiallyApplied.apply(undefined, args);
    return nextHandle++;
  } // This function accepts the same arguments as setImmediate, but
  // returns a function that requires no arguments.


  function partiallyApplied(handler) {
    var args = [].slice.call(arguments, 1);
    return function () {
      if (typeof handler === "function") {
        handler.apply(undefined, args);
      } else {
        new Function("" + handler)();
      }
    };
  }

  function runIfPresent(handle) {
    // From the spec: "Wait until any invocations of this algorithm started before this one have completed."
    // So if we're currently running a task, we'll need to delay this invocation.
    if (currentlyRunningATask) {
      // Delay by doing a setTimeout. setImmediate was tried instead, but in Firefox 7 it generated a
      // "too much recursion" error.
      setTimeout(partiallyApplied(runIfPresent, handle), 0);
    } else {
      var task = tasksByHandle[handle];

      if (task) {
        currentlyRunningATask = true;

        try {
          task();
        } finally {
          clearImmediate(handle);
          currentlyRunningATask = false;
        }
      }
    }
  }

  function clearImmediate(handle) {
    delete tasksByHandle[handle];
  }

  function installNextTickImplementation() {
    setImmediate = function setImmediate() {
      var handle = addFromSetImmediateArguments(arguments);
      process.nextTick(partiallyApplied(runIfPresent, handle));
      return handle;
    };
  }

  function canUsePostMessage() {
    // The test against `importScripts` prevents this implementation from being installed inside a web worker,
    // where `global.postMessage` means something completely different and can't be used for this purpose.
    if (global.postMessage && !global.importScripts) {
      var postMessageIsAsynchronous = true;
      var oldOnMessage = global.onmessage;

      global.onmessage = function () {
        postMessageIsAsynchronous = false;
      };

      global.postMessage("", "*");
      global.onmessage = oldOnMessage;
      return postMessageIsAsynchronous;
    }
  }

  function installPostMessageImplementation() {
    // Installs an event handler on `global` for the `message` event: see
    // * https://developer.mozilla.org/en/DOM/window.postMessage
    // * http://www.whatwg.org/specs/web-apps/current-work/multipage/comms.html#crossDocumentMessages
    var messagePrefix = "setImmediate$" + Math.random() + "$";

    var onGlobalMessage = function onGlobalMessage(event) {
      if (event.source === global && typeof event.data === "string" && event.data.indexOf(messagePrefix) === 0) {
        runIfPresent(+event.data.slice(messagePrefix.length));
      }
    };

    if (global.addEventListener) {
      global.addEventListener("message", onGlobalMessage, false);
    } else {
      global.attachEvent("onmessage", onGlobalMessage);
    }

    setImmediate = function setImmediate() {
      var handle = addFromSetImmediateArguments(arguments);
      global.postMessage(messagePrefix + handle, "*");
      return handle;
    };
  }

  function installMessageChannelImplementation() {
    var channel = new MessageChannel();

    channel.port1.onmessage = function (event) {
      var handle = event.data;
      runIfPresent(handle);
    };

    setImmediate = function setImmediate() {
      var handle = addFromSetImmediateArguments(arguments);
      channel.port2.postMessage(handle);
      return handle;
    };
  }

  function installReadyStateChangeImplementation() {
    var html = doc.documentElement;

    setImmediate = function setImmediate() {
      var handle = addFromSetImmediateArguments(arguments); // Create a <script> element; its readystatechange event will be fired asynchronously once it is inserted
      // into the document. Do so, thus queuing up the task. Remember to clean up once it's been called.

      var script = doc.createElement("script");

      script.onreadystatechange = function () {
        runIfPresent(handle);
        script.onreadystatechange = null;
        html.removeChild(script);
        script = null;
      };

      html.appendChild(script);
      return handle;
    };
  }

  function installSetTimeoutImplementation() {
    setImmediate = function setImmediate() {
      var handle = addFromSetImmediateArguments(arguments);
      setTimeout(partiallyApplied(runIfPresent, handle), 0);
      return handle;
    };
  } // If supported, we should attach to the prototype of global, since that is where setTimeout et al. live.


  var attachTo = Object.getPrototypeOf && Object.getPrototypeOf(global);
  attachTo = attachTo && attachTo.setTimeout ? attachTo : global; // Don't get fooled by e.g. browserify environments.

  if ({}.toString.call(global.process) === "[object process]") {
    // For Node.js before 0.9
    installNextTickImplementation();
  } else if (canUsePostMessage()) {
    // For non-IE10 modern browsers
    installPostMessageImplementation();
  } else if (global.MessageChannel) {
    // For web workers, where supported
    installMessageChannelImplementation();
  } else if (doc && "onreadystatechange" in doc.createElement("script")) {
    // For IE 6–8
    installReadyStateChangeImplementation();
  } else {
    // For older browsers
    installSetTimeoutImplementation();
  }

  attachTo.setImmediate = setImmediate;
  attachTo.clearImmediate = clearImmediate;
})(typeof self === "undefined" ? typeof global === "undefined" ? this : global : self);
/* WEBPACK VAR INJECTION */}.call(this, __webpack_require__(4), __webpack_require__(18)))

/***/ }),
/* 7 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


module.exports = __webpack_require__(0).Matrix;
module.exports.Decompositions = module.exports.DC = __webpack_require__(49);

/***/ }),
/* 8 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


module.exports = abstractMatrix;

var LuDecomposition = __webpack_require__(9);

var SvDecomposition = __webpack_require__(10);

var arrayUtils = __webpack_require__(37);

var util = __webpack_require__(3);

var MatrixTransposeView = __webpack_require__(42);

var MatrixRowView = __webpack_require__(43);

var MatrixSubView = __webpack_require__(44);

var MatrixSelectionView = __webpack_require__(45);

var MatrixColumnView = __webpack_require__(46);

var MatrixFlipRowView = __webpack_require__(47);

var MatrixFlipColumnView = __webpack_require__(48);

function abstractMatrix(superCtor) {
  if (superCtor === undefined) superCtor = Object;
  /**
   * Real matrix
   * @class Matrix
   * @param {number|Array|Matrix} nRows - Number of rows of the new matrix,
   * 2D array containing the data or Matrix instance to clone
   * @param {number} [nColumns] - Number of columns of the new matrix
   */

  class Matrix extends superCtor {
    static get [Symbol.species]() {
      return this;
    }
    /**
     * Constructs a Matrix with the chosen dimensions from a 1D array
     * @param {number} newRows - Number of rows
     * @param {number} newColumns - Number of columns
     * @param {Array} newData - A 1D array containing data for the matrix
     * @return {Matrix} - The new matrix
     */


    static from1DArray(newRows, newColumns, newData) {
      var length = newRows * newColumns;

      if (length !== newData.length) {
        throw new RangeError('Data length does not match given dimensions');
      }

      var newMatrix = new this(newRows, newColumns);

      for (var row = 0; row < newRows; row++) {
        for (var column = 0; column < newColumns; column++) {
          newMatrix.set(row, column, newData[row * newColumns + column]);
        }
      }

      return newMatrix;
    }
    /**
     * Creates a row vector, a matrix with only one row.
     * @param {Array} newData - A 1D array containing data for the vector
     * @return {Matrix} - The new matrix
     */


    static rowVector(newData) {
      var vector = new this(1, newData.length);

      for (var i = 0; i < newData.length; i++) {
        vector.set(0, i, newData[i]);
      }

      return vector;
    }
    /**
     * Creates a column vector, a matrix with only one column.
     * @param {Array} newData - A 1D array containing data for the vector
     * @return {Matrix} - The new matrix
     */


    static columnVector(newData) {
      var vector = new this(newData.length, 1);

      for (var i = 0; i < newData.length; i++) {
        vector.set(i, 0, newData[i]);
      }

      return vector;
    }
    /**
     * Creates an empty matrix with the given dimensions. Values will be undefined. Same as using new Matrix(rows, columns).
     * @param {number} rows - Number of rows
     * @param {number} columns - Number of columns
     * @return {Matrix} - The new matrix
     */


    static empty(rows, columns) {
      return new this(rows, columns);
    }
    /**
     * Creates a matrix with the given dimensions. Values will be set to zero.
     * @param {number} rows - Number of rows
     * @param {number} columns - Number of columns
     * @return {Matrix} - The new matrix
     */


    static zeros(rows, columns) {
      return this.empty(rows, columns).fill(0);
    }
    /**
     * Creates a matrix with the given dimensions. Values will be set to one.
     * @param {number} rows - Number of rows
     * @param {number} columns - Number of columns
     * @return {Matrix} - The new matrix
     */


    static ones(rows, columns) {
      return this.empty(rows, columns).fill(1);
    }
    /**
     * Creates a matrix with the given dimensions. Values will be randomly set.
     * @param {number} rows - Number of rows
     * @param {number} columns - Number of columns
     * @param {function} [rng=Math.random] - Random number generator
     * @return {Matrix} The new matrix
     */


    static rand(rows, columns, rng) {
      if (rng === undefined) rng = Math.random;
      var matrix = this.empty(rows, columns);

      for (var i = 0; i < rows; i++) {
        for (var j = 0; j < columns; j++) {
          matrix.set(i, j, rng());
        }
      }

      return matrix;
    }
    /**
     * Creates a matrix with the given dimensions. Values will be random integers.
     * @param {number} rows - Number of rows
     * @param {number} columns - Number of columns
     * @param {number} [maxValue=1000] - Maximum value
     * @param {function} [rng=Math.random] - Random number generator
     * @return {Matrix} The new matrix
     */


    static randInt(rows, columns, maxValue, rng) {
      if (maxValue === undefined) maxValue = 1000;
      if (rng === undefined) rng = Math.random;
      var matrix = this.empty(rows, columns);

      for (var i = 0; i < rows; i++) {
        for (var j = 0; j < columns; j++) {
          var value = Math.floor(rng() * maxValue);
          matrix.set(i, j, value);
        }
      }

      return matrix;
    }
    /**
     * Creates an identity matrix with the given dimension. Values of the diagonal will be 1 and others will be 0.
     * @param {number} rows - Number of rows
     * @param {number} [columns=rows] - Number of columns
     * @param {number} [value=1] - Value to fill the diagonal with
     * @return {Matrix} - The new identity matrix
     */


    static eye(rows, columns, value) {
      if (columns === undefined) columns = rows;
      if (value === undefined) value = 1;
      var min = Math.min(rows, columns);
      var matrix = this.zeros(rows, columns);

      for (var i = 0; i < min; i++) {
        matrix.set(i, i, value);
      }

      return matrix;
    }
    /**
     * Creates a diagonal matrix based on the given array.
     * @param {Array} data - Array containing the data for the diagonal
     * @param {number} [rows] - Number of rows (Default: data.length)
     * @param {number} [columns] - Number of columns (Default: rows)
     * @return {Matrix} - The new diagonal matrix
     */


    static diag(data, rows, columns) {
      var l = data.length;
      if (rows === undefined) rows = l;
      if (columns === undefined) columns = rows;
      var min = Math.min(l, rows, columns);
      var matrix = this.zeros(rows, columns);

      for (var i = 0; i < min; i++) {
        matrix.set(i, i, data[i]);
      }

      return matrix;
    }
    /**
     * Returns a matrix whose elements are the minimum between matrix1 and matrix2
     * @param {Matrix} matrix1
     * @param {Matrix} matrix2
     * @return {Matrix}
     */


    static min(matrix1, matrix2) {
      matrix1 = this.checkMatrix(matrix1);
      matrix2 = this.checkMatrix(matrix2);
      var rows = matrix1.rows;
      var columns = matrix1.columns;
      var result = new this(rows, columns);

      for (var i = 0; i < rows; i++) {
        for (var j = 0; j < columns; j++) {
          result.set(i, j, Math.min(matrix1.get(i, j), matrix2.get(i, j)));
        }
      }

      return result;
    }
    /**
     * Returns a matrix whose elements are the maximum between matrix1 and matrix2
     * @param {Matrix} matrix1
     * @param {Matrix} matrix2
     * @return {Matrix}
     */


    static max(matrix1, matrix2) {
      matrix1 = this.checkMatrix(matrix1);
      matrix2 = this.checkMatrix(matrix2);
      var rows = matrix1.rows;
      var columns = matrix1.columns;
      var result = new this(rows, columns);

      for (var i = 0; i < rows; i++) {
        for (var j = 0; j < columns; j++) {
          result.set(i, j, Math.max(matrix1.get(i, j), matrix2.get(i, j)));
        }
      }

      return result;
    }
    /**
     * Check that the provided value is a Matrix and tries to instantiate one if not
     * @param {*} value - The value to check
     * @return {Matrix}
     */


    static checkMatrix(value) {
      return Matrix.isMatrix(value) ? value : new this(value);
    }
    /**
     * Returns true if the argument is a Matrix, false otherwise
     * @param {*} value - The value to check
     * @return {boolean}
     */


    static isMatrix(value) {
      return value != null && value.klass === 'Matrix';
    }
    /**
     * @prop {number} size - The number of elements in the matrix.
     */


    get size() {
      return this.rows * this.columns;
    }
    /**
     * Applies a callback for each element of the matrix. The function is called in the matrix (this) context.
     * @param {function} callback - Function that will be called with two parameters : i (row) and j (column)
     * @return {Matrix} this
     */


    apply(callback) {
      if (typeof callback !== 'function') {
        throw new TypeError('callback must be a function');
      }

      var ii = this.rows;
      var jj = this.columns;

      for (var i = 0; i < ii; i++) {
        for (var j = 0; j < jj; j++) {
          callback.call(this, i, j);
        }
      }

      return this;
    }
    /**
     * Returns a new 1D array filled row by row with the matrix values
     * @return {Array}
     */


    to1DArray() {
      var array = new Array(this.size);

      for (var i = 0; i < this.rows; i++) {
        for (var j = 0; j < this.columns; j++) {
          array[i * this.columns + j] = this.get(i, j);
        }
      }

      return array;
    }
    /**
     * Returns a 2D array containing a copy of the data
     * @return {Array}
     */


    to2DArray() {
      var copy = new Array(this.rows);

      for (var i = 0; i < this.rows; i++) {
        copy[i] = new Array(this.columns);

        for (var j = 0; j < this.columns; j++) {
          copy[i][j] = this.get(i, j);
        }
      }

      return copy;
    }
    /**
     * @return {boolean} true if the matrix has one row
     */


    isRowVector() {
      return this.rows === 1;
    }
    /**
     * @return {boolean} true if the matrix has one column
     */


    isColumnVector() {
      return this.columns === 1;
    }
    /**
     * @return {boolean} true if the matrix has one row or one column
     */


    isVector() {
      return this.rows === 1 || this.columns === 1;
    }
    /**
     * @return {boolean} true if the matrix has the same number of rows and columns
     */


    isSquare() {
      return this.rows === this.columns;
    }
    /**
     * @return {boolean} true if the matrix is square and has the same values on both sides of the diagonal
     */


    isSymmetric() {
      if (this.isSquare()) {
        for (var i = 0; i < this.rows; i++) {
          for (var j = 0; j <= i; j++) {
            if (this.get(i, j) !== this.get(j, i)) {
              return false;
            }
          }
        }

        return true;
      }

      return false;
    }
    /**
     * Sets a given element of the matrix. mat.set(3,4,1) is equivalent to mat[3][4]=1
     * @abstract
     * @param {number} rowIndex - Index of the row
     * @param {number} columnIndex - Index of the column
     * @param {number} value - The new value for the element
     * @return {Matrix} this
     */


    set(rowIndex, columnIndex, value) {
      // eslint-disable-line no-unused-vars
      throw new Error('set method is unimplemented');
    }
    /**
     * Returns the given element of the matrix. mat.get(3,4) is equivalent to matrix[3][4]
     * @abstract
     * @param {number} rowIndex - Index of the row
     * @param {number} columnIndex - Index of the column
     * @return {number}
     */


    get(rowIndex, columnIndex) {
      // eslint-disable-line no-unused-vars
      throw new Error('get method is unimplemented');
    }
    /**
     * Creates a new matrix that is a repetition of the current matrix. New matrix has rowRep times the number of
     * rows of the matrix, and colRep times the number of columns of the matrix
     * @param {number} rowRep - Number of times the rows should be repeated
     * @param {number} colRep - Number of times the columns should be re
     * @return {Matrix}
     * @example
     * var matrix = new Matrix([[1,2]]);
     * matrix.repeat(2); // [[1,2],[1,2]]
     */


    repeat(rowRep, colRep) {
      rowRep = rowRep || 1;
      colRep = colRep || 1;
      var matrix = new this.constructor[Symbol.species](this.rows * rowRep, this.columns * colRep);

      for (var i = 0; i < rowRep; i++) {
        for (var j = 0; j < colRep; j++) {
          matrix.setSubMatrix(this, this.rows * i, this.columns * j);
        }
      }

      return matrix;
    }
    /**
     * Fills the matrix with a given value. All elements will be set to this value.
     * @param {number} value - New value
     * @return {Matrix} this
     */


    fill(value) {
      for (var i = 0; i < this.rows; i++) {
        for (var j = 0; j < this.columns; j++) {
          this.set(i, j, value);
        }
      }

      return this;
    }
    /**
     * Negates the matrix. All elements will be multiplied by (-1)
     * @return {Matrix} this
     */


    neg() {
      return this.mulS(-1);
    }
    /**
     * Returns a new array from the given row index
     * @param {number} index - Row index
     * @return {Array}
     */


    getRow(index) {
      util.checkRowIndex(this, index);
      var row = new Array(this.columns);

      for (var i = 0; i < this.columns; i++) {
        row[i] = this.get(index, i);
      }

      return row;
    }
    /**
     * Returns a new row vector from the given row index
     * @param {number} index - Row index
     * @return {Matrix}
     */


    getRowVector(index) {
      return this.constructor.rowVector(this.getRow(index));
    }
    /**
     * Sets a row at the given index
     * @param {number} index - Row index
     * @param {Array|Matrix} array - Array or vector
     * @return {Matrix} this
     */


    setRow(index, array) {
      util.checkRowIndex(this, index);
      array = util.checkRowVector(this, array);

      for (var i = 0; i < this.columns; i++) {
        this.set(index, i, array[i]);
      }

      return this;
    }
    /**
     * Swaps two rows
     * @param {number} row1 - First row index
     * @param {number} row2 - Second row index
     * @return {Matrix} this
     */


    swapRows(row1, row2) {
      util.checkRowIndex(this, row1);
      util.checkRowIndex(this, row2);

      for (var i = 0; i < this.columns; i++) {
        var temp = this.get(row1, i);
        this.set(row1, i, this.get(row2, i));
        this.set(row2, i, temp);
      }

      return this;
    }
    /**
     * Returns a new array from the given column index
     * @param {number} index - Column index
     * @return {Array}
     */


    getColumn(index) {
      util.checkColumnIndex(this, index);
      var column = new Array(this.rows);

      for (var i = 0; i < this.rows; i++) {
        column[i] = this.get(i, index);
      }

      return column;
    }
    /**
     * Returns a new column vector from the given column index
     * @param {number} index - Column index
     * @return {Matrix}
     */


    getColumnVector(index) {
      return this.constructor.columnVector(this.getColumn(index));
    }
    /**
     * Sets a column at the given index
     * @param {number} index - Column index
     * @param {Array|Matrix} array - Array or vector
     * @return {Matrix} this
     */


    setColumn(index, array) {
      util.checkColumnIndex(this, index);
      array = util.checkColumnVector(this, array);

      for (var i = 0; i < this.rows; i++) {
        this.set(i, index, array[i]);
      }

      return this;
    }
    /**
     * Swaps two columns
     * @param {number} column1 - First column index
     * @param {number} column2 - Second column index
     * @return {Matrix} this
     */


    swapColumns(column1, column2) {
      util.checkColumnIndex(this, column1);
      util.checkColumnIndex(this, column2);

      for (var i = 0; i < this.rows; i++) {
        var temp = this.get(i, column1);
        this.set(i, column1, this.get(i, column2));
        this.set(i, column2, temp);
      }

      return this;
    }
    /**
     * Adds the values of a vector to each row
     * @param {Array|Matrix} vector - Array or vector
     * @return {Matrix} this
     */


    addRowVector(vector) {
      vector = util.checkRowVector(this, vector);

      for (var i = 0; i < this.rows; i++) {
        for (var j = 0; j < this.columns; j++) {
          this.set(i, j, this.get(i, j) + vector[j]);
        }
      }

      return this;
    }
    /**
     * Subtracts the values of a vector from each row
     * @param {Array|Matrix} vector - Array or vector
     * @return {Matrix} this
     */


    subRowVector(vector) {
      vector = util.checkRowVector(this, vector);

      for (var i = 0; i < this.rows; i++) {
        for (var j = 0; j < this.columns; j++) {
          this.set(i, j, this.get(i, j) - vector[j]);
        }
      }

      return this;
    }
    /**
     * Multiplies the values of a vector with each row
     * @param {Array|Matrix} vector - Array or vector
     * @return {Matrix} this
     */


    mulRowVector(vector) {
      vector = util.checkRowVector(this, vector);

      for (var i = 0; i < this.rows; i++) {
        for (var j = 0; j < this.columns; j++) {
          this.set(i, j, this.get(i, j) * vector[j]);
        }
      }

      return this;
    }
    /**
     * Divides the values of each row by those of a vector
     * @param {Array|Matrix} vector - Array or vector
     * @return {Matrix} this
     */


    divRowVector(vector) {
      vector = util.checkRowVector(this, vector);

      for (var i = 0; i < this.rows; i++) {
        for (var j = 0; j < this.columns; j++) {
          this.set(i, j, this.get(i, j) / vector[j]);
        }
      }

      return this;
    }
    /**
     * Adds the values of a vector to each column
     * @param {Array|Matrix} vector - Array or vector
     * @return {Matrix} this
     */


    addColumnVector(vector) {
      vector = util.checkColumnVector(this, vector);

      for (var i = 0; i < this.rows; i++) {
        for (var j = 0; j < this.columns; j++) {
          this.set(i, j, this.get(i, j) + vector[i]);
        }
      }

      return this;
    }
    /**
     * Subtracts the values of a vector from each column
     * @param {Array|Matrix} vector - Array or vector
     * @return {Matrix} this
     */


    subColumnVector(vector) {
      vector = util.checkColumnVector(this, vector);

      for (var i = 0; i < this.rows; i++) {
        for (var j = 0; j < this.columns; j++) {
          this.set(i, j, this.get(i, j) - vector[i]);
        }
      }

      return this;
    }
    /**
     * Multiplies the values of a vector with each column
     * @param {Array|Matrix} vector - Array or vector
     * @return {Matrix} this
     */


    mulColumnVector(vector) {
      vector = util.checkColumnVector(this, vector);

      for (var i = 0; i < this.rows; i++) {
        for (var j = 0; j < this.columns; j++) {
          this.set(i, j, this.get(i, j) * vector[i]);
        }
      }

      return this;
    }
    /**
     * Divides the values of each column by those of a vector
     * @param {Array|Matrix} vector - Array or vector
     * @return {Matrix} this
     */


    divColumnVector(vector) {
      vector = util.checkColumnVector(this, vector);

      for (var i = 0; i < this.rows; i++) {
        for (var j = 0; j < this.columns; j++) {
          this.set(i, j, this.get(i, j) / vector[i]);
        }
      }

      return this;
    }
    /**
     * Multiplies the values of a row with a scalar
     * @param {number} index - Row index
     * @param {number} value
     * @return {Matrix} this
     */


    mulRow(index, value) {
      util.checkRowIndex(this, index);

      for (var i = 0; i < this.columns; i++) {
        this.set(index, i, this.get(index, i) * value);
      }

      return this;
    }
    /**
     * Multiplies the values of a column with a scalar
     * @param {number} index - Column index
     * @param {number} value
     * @return {Matrix} this
     */


    mulColumn(index, value) {
      util.checkColumnIndex(this, index);

      for (var i = 0; i < this.rows; i++) {
        this.set(i, index, this.get(i, index) * value);
      }

      return this;
    }
    /**
     * Returns the maximum value of the matrix
     * @return {number}
     */


    max() {
      var v = this.get(0, 0);

      for (var i = 0; i < this.rows; i++) {
        for (var j = 0; j < this.columns; j++) {
          if (this.get(i, j) > v) {
            v = this.get(i, j);
          }
        }
      }

      return v;
    }
    /**
     * Returns the index of the maximum value
     * @return {Array}
     */


    maxIndex() {
      var v = this.get(0, 0);
      var idx = [0, 0];

      for (var i = 0; i < this.rows; i++) {
        for (var j = 0; j < this.columns; j++) {
          if (this.get(i, j) > v) {
            v = this.get(i, j);
            idx[0] = i;
            idx[1] = j;
          }
        }
      }

      return idx;
    }
    /**
     * Returns the minimum value of the matrix
     * @return {number}
     */


    min() {
      var v = this.get(0, 0);

      for (var i = 0; i < this.rows; i++) {
        for (var j = 0; j < this.columns; j++) {
          if (this.get(i, j) < v) {
            v = this.get(i, j);
          }
        }
      }

      return v;
    }
    /**
     * Returns the index of the minimum value
     * @return {Array}
     */


    minIndex() {
      var v = this.get(0, 0);
      var idx = [0, 0];

      for (var i = 0; i < this.rows; i++) {
        for (var j = 0; j < this.columns; j++) {
          if (this.get(i, j) < v) {
            v = this.get(i, j);
            idx[0] = i;
            idx[1] = j;
          }
        }
      }

      return idx;
    }
    /**
     * Returns the maximum value of one row
     * @param {number} row - Row index
     * @return {number}
     */


    maxRow(row) {
      util.checkRowIndex(this, row);
      var v = this.get(row, 0);

      for (var i = 1; i < this.columns; i++) {
        if (this.get(row, i) > v) {
          v = this.get(row, i);
        }
      }

      return v;
    }
    /**
     * Returns the index of the maximum value of one row
     * @param {number} row - Row index
     * @return {Array}
     */


    maxRowIndex(row) {
      util.checkRowIndex(this, row);
      var v = this.get(row, 0);
      var idx = [row, 0];

      for (var i = 1; i < this.columns; i++) {
        if (this.get(row, i) > v) {
          v = this.get(row, i);
          idx[1] = i;
        }
      }

      return idx;
    }
    /**
     * Returns the minimum value of one row
     * @param {number} row - Row index
     * @return {number}
     */


    minRow(row) {
      util.checkRowIndex(this, row);
      var v = this.get(row, 0);

      for (var i = 1; i < this.columns; i++) {
        if (this.get(row, i) < v) {
          v = this.get(row, i);
        }
      }

      return v;
    }
    /**
     * Returns the index of the maximum value of one row
     * @param {number} row - Row index
     * @return {Array}
     */


    minRowIndex(row) {
      util.checkRowIndex(this, row);
      var v = this.get(row, 0);
      var idx = [row, 0];

      for (var i = 1; i < this.columns; i++) {
        if (this.get(row, i) < v) {
          v = this.get(row, i);
          idx[1] = i;
        }
      }

      return idx;
    }
    /**
     * Returns the maximum value of one column
     * @param {number} column - Column index
     * @return {number}
     */


    maxColumn(column) {
      util.checkColumnIndex(this, column);
      var v = this.get(0, column);

      for (var i = 1; i < this.rows; i++) {
        if (this.get(i, column) > v) {
          v = this.get(i, column);
        }
      }

      return v;
    }
    /**
     * Returns the index of the maximum value of one column
     * @param {number} column - Column index
     * @return {Array}
     */


    maxColumnIndex(column) {
      util.checkColumnIndex(this, column);
      var v = this.get(0, column);
      var idx = [0, column];

      for (var i = 1; i < this.rows; i++) {
        if (this.get(i, column) > v) {
          v = this.get(i, column);
          idx[0] = i;
        }
      }

      return idx;
    }
    /**
     * Returns the minimum value of one column
     * @param {number} column - Column index
     * @return {number}
     */


    minColumn(column) {
      util.checkColumnIndex(this, column);
      var v = this.get(0, column);

      for (var i = 1; i < this.rows; i++) {
        if (this.get(i, column) < v) {
          v = this.get(i, column);
        }
      }

      return v;
    }
    /**
     * Returns the index of the minimum value of one column
     * @param {number} column - Column index
     * @return {Array}
     */


    minColumnIndex(column) {
      util.checkColumnIndex(this, column);
      var v = this.get(0, column);
      var idx = [0, column];

      for (var i = 1; i < this.rows; i++) {
        if (this.get(i, column) < v) {
          v = this.get(i, column);
          idx[0] = i;
        }
      }

      return idx;
    }
    /**
     * Returns an array containing the diagonal values of the matrix
     * @return {Array}
     */


    diag() {
      var min = Math.min(this.rows, this.columns);
      var diag = new Array(min);

      for (var i = 0; i < min; i++) {
        diag[i] = this.get(i, i);
      }

      return diag;
    }
    /**
     * Returns the sum by the argument given, if no argument given,
     * it returns the sum of all elements of the matrix.
     * @param {string} by - sum by 'row' or 'column'.
     * @return {Matrix|number}
     */


    sum(by) {
      switch (by) {
        case 'row':
          return util.sumByRow(this);

        case 'column':
          return util.sumByColumn(this);

        default:
          return util.sumAll(this);
      }
    }
    /**
     * Returns the mean of all elements of the matrix
     * @return {number}
     */


    mean() {
      return this.sum() / this.size;
    }
    /**
     * Returns the product of all elements of the matrix
     * @return {number}
     */


    prod() {
      var prod = 1;

      for (var i = 0; i < this.rows; i++) {
        for (var j = 0; j < this.columns; j++) {
          prod *= this.get(i, j);
        }
      }

      return prod;
    }
    /**
     * Computes the cumulative sum of the matrix elements (in place, row by row)
     * @return {Matrix} this
     */


    cumulativeSum() {
      var sum = 0;

      for (var i = 0; i < this.rows; i++) {
        for (var j = 0; j < this.columns; j++) {
          sum += this.get(i, j);
          this.set(i, j, sum);
        }
      }

      return this;
    }
    /**
     * Computes the dot (scalar) product between the matrix and another
     * @param {Matrix} vector2 vector
     * @return {number}
     */


    dot(vector2) {
      if (Matrix.isMatrix(vector2)) vector2 = vector2.to1DArray();
      var vector1 = this.to1DArray();

      if (vector1.length !== vector2.length) {
        throw new RangeError('vectors do not have the same size');
      }

      var dot = 0;

      for (var i = 0; i < vector1.length; i++) {
        dot += vector1[i] * vector2[i];
      }

      return dot;
    }
    /**
     * Returns the matrix product between this and other
     * @param {Matrix} other
     * @return {Matrix}
     */


    mmul(other) {
      other = this.constructor.checkMatrix(other);

      if (this.columns !== other.rows) {
        // eslint-disable-next-line no-console
        console.warn('Number of columns of left matrix are not equal to number of rows of right matrix.');
      }

      var m = this.rows;
      var n = this.columns;
      var p = other.columns;
      var result = new this.constructor[Symbol.species](m, p);
      var Bcolj = new Array(n);

      for (var j = 0; j < p; j++) {
        for (var k = 0; k < n; k++) {
          Bcolj[k] = other.get(k, j);
        }

        for (var i = 0; i < m; i++) {
          var s = 0;

          for (k = 0; k < n; k++) {
            s += this.get(i, k) * Bcolj[k];
          }

          result.set(i, j, s);
        }
      }

      return result;
    }

    strassen2x2(other) {
      var result = new this.constructor[Symbol.species](2, 2);
      const a11 = this.get(0, 0);
      const b11 = other.get(0, 0);
      const a12 = this.get(0, 1);
      const b12 = other.get(0, 1);
      const a21 = this.get(1, 0);
      const b21 = other.get(1, 0);
      const a22 = this.get(1, 1);
      const b22 = other.get(1, 1); // Compute intermediate values.

      const m1 = (a11 + a22) * (b11 + b22);
      const m2 = (a21 + a22) * b11;
      const m3 = a11 * (b12 - b22);
      const m4 = a22 * (b21 - b11);
      const m5 = (a11 + a12) * b22;
      const m6 = (a21 - a11) * (b11 + b12);
      const m7 = (a12 - a22) * (b21 + b22); // Combine intermediate values into the output.

      const c00 = m1 + m4 - m5 + m7;
      const c01 = m3 + m5;
      const c10 = m2 + m4;
      const c11 = m1 - m2 + m3 + m6;
      result.set(0, 0, c00);
      result.set(0, 1, c01);
      result.set(1, 0, c10);
      result.set(1, 1, c11);
      return result;
    }

    strassen3x3(other) {
      var result = new this.constructor[Symbol.species](3, 3);
      const a00 = this.get(0, 0);
      const a01 = this.get(0, 1);
      const a02 = this.get(0, 2);
      const a10 = this.get(1, 0);
      const a11 = this.get(1, 1);
      const a12 = this.get(1, 2);
      const a20 = this.get(2, 0);
      const a21 = this.get(2, 1);
      const a22 = this.get(2, 2);
      const b00 = other.get(0, 0);
      const b01 = other.get(0, 1);
      const b02 = other.get(0, 2);
      const b10 = other.get(1, 0);
      const b11 = other.get(1, 1);
      const b12 = other.get(1, 2);
      const b20 = other.get(2, 0);
      const b21 = other.get(2, 1);
      const b22 = other.get(2, 2);
      const m1 = (a00 + a01 + a02 - a10 - a11 - a21 - a22) * b11;
      const m2 = (a00 - a10) * (-b01 + b11);
      const m3 = a11 * (-b00 + b01 + b10 - b11 - b12 - b20 + b22);
      const m4 = (-a00 + a10 + a11) * (b00 - b01 + b11);
      const m5 = (a10 + a11) * (-b00 + b01);
      const m6 = a00 * b00;
      const m7 = (-a00 + a20 + a21) * (b00 - b02 + b12);
      const m8 = (-a00 + a20) * (b02 - b12);
      const m9 = (a20 + a21) * (-b00 + b02);
      const m10 = (a00 + a01 + a02 - a11 - a12 - a20 - a21) * b12;
      const m11 = a21 * (-b00 + b02 + b10 - b11 - b12 - b20 + b21);
      const m12 = (-a02 + a21 + a22) * (b11 + b20 - b21);
      const m13 = (a02 - a22) * (b11 - b21);
      const m14 = a02 * b20;
      const m15 = (a21 + a22) * (-b20 + b21);
      const m16 = (-a02 + a11 + a12) * (b12 + b20 - b22);
      const m17 = (a02 - a12) * (b12 - b22);
      const m18 = (a11 + a12) * (-b20 + b22);
      const m19 = a01 * b10;
      const m20 = a12 * b21;
      const m21 = a10 * b02;
      const m22 = a20 * b01;
      const m23 = a22 * b22;
      const c00 = m6 + m14 + m19;
      const c01 = m1 + m4 + m5 + m6 + m12 + m14 + m15;
      const c02 = m6 + m7 + m9 + m10 + m14 + m16 + m18;
      const c10 = m2 + m3 + m4 + m6 + m14 + m16 + m17;
      const c11 = m2 + m4 + m5 + m6 + m20;
      const c12 = m14 + m16 + m17 + m18 + m21;
      const c20 = m6 + m7 + m8 + m11 + m12 + m13 + m14;
      const c21 = m12 + m13 + m14 + m15 + m22;
      const c22 = m6 + m7 + m8 + m9 + m23;
      result.set(0, 0, c00);
      result.set(0, 1, c01);
      result.set(0, 2, c02);
      result.set(1, 0, c10);
      result.set(1, 1, c11);
      result.set(1, 2, c12);
      result.set(2, 0, c20);
      result.set(2, 1, c21);
      result.set(2, 2, c22);
      return result;
    }
    /**
     * Returns the matrix product between x and y. More efficient than mmul(other) only when we multiply squared matrix and when the size of the matrix is > 1000.
     * @param {Matrix} y
     * @return {Matrix}
     */


    mmulStrassen(y) {
      var x = this.clone();
      var r1 = x.rows;
      var c1 = x.columns;
      var r2 = y.rows;
      var c2 = y.columns;

      if (c1 !== r2) {
        // eslint-disable-next-line no-console
        console.warn(`Multiplying ${r1} x ${c1} and ${r2} x ${c2} matrix: dimensions do not match.`);
      } // Put a matrix into the top left of a matrix of zeros.
      // `rows` and `cols` are the dimensions of the output matrix.


      function embed(mat, rows, cols) {
        var r = mat.rows;
        var c = mat.columns;

        if (r === rows && c === cols) {
          return mat;
        } else {
          var resultat = Matrix.zeros(rows, cols);
          resultat = resultat.setSubMatrix(mat, 0, 0);
          return resultat;
        }
      } // Make sure both matrices are the same size.
      // This is exclusively for simplicity:
      // this algorithm can be implemented with matrices of different sizes.


      var r = Math.max(r1, r2);
      var c = Math.max(c1, c2);
      x = embed(x, r, c);
      y = embed(y, r, c); // Our recursive multiplication function.

      function blockMult(a, b, rows, cols) {
        // For small matrices, resort to naive multiplication.
        if (rows <= 512 || cols <= 512) {
          return a.mmul(b); // a is equivalent to this
        } // Apply dynamic padding.


        if (rows % 2 === 1 && cols % 2 === 1) {
          a = embed(a, rows + 1, cols + 1);
          b = embed(b, rows + 1, cols + 1);
        } else if (rows % 2 === 1) {
          a = embed(a, rows + 1, cols);
          b = embed(b, rows + 1, cols);
        } else if (cols % 2 === 1) {
          a = embed(a, rows, cols + 1);
          b = embed(b, rows, cols + 1);
        }

        var halfRows = parseInt(a.rows / 2);
        var halfCols = parseInt(a.columns / 2); // Subdivide input matrices.

        var a11 = a.subMatrix(0, halfRows - 1, 0, halfCols - 1);
        var b11 = b.subMatrix(0, halfRows - 1, 0, halfCols - 1);
        var a12 = a.subMatrix(0, halfRows - 1, halfCols, a.columns - 1);
        var b12 = b.subMatrix(0, halfRows - 1, halfCols, b.columns - 1);
        var a21 = a.subMatrix(halfRows, a.rows - 1, 0, halfCols - 1);
        var b21 = b.subMatrix(halfRows, b.rows - 1, 0, halfCols - 1);
        var a22 = a.subMatrix(halfRows, a.rows - 1, halfCols, a.columns - 1);
        var b22 = b.subMatrix(halfRows, b.rows - 1, halfCols, b.columns - 1); // Compute intermediate values.

        var m1 = blockMult(Matrix.add(a11, a22), Matrix.add(b11, b22), halfRows, halfCols);
        var m2 = blockMult(Matrix.add(a21, a22), b11, halfRows, halfCols);
        var m3 = blockMult(a11, Matrix.sub(b12, b22), halfRows, halfCols);
        var m4 = blockMult(a22, Matrix.sub(b21, b11), halfRows, halfCols);
        var m5 = blockMult(Matrix.add(a11, a12), b22, halfRows, halfCols);
        var m6 = blockMult(Matrix.sub(a21, a11), Matrix.add(b11, b12), halfRows, halfCols);
        var m7 = blockMult(Matrix.sub(a12, a22), Matrix.add(b21, b22), halfRows, halfCols); // Combine intermediate values into the output.

        var c11 = Matrix.add(m1, m4);
        c11.sub(m5);
        c11.add(m7);
        var c12 = Matrix.add(m3, m5);
        var c21 = Matrix.add(m2, m4);
        var c22 = Matrix.sub(m1, m2);
        c22.add(m3);
        c22.add(m6); //Crop output to the desired size (undo dynamic padding).

        var resultat = Matrix.zeros(2 * c11.rows, 2 * c11.columns);
        resultat = resultat.setSubMatrix(c11, 0, 0);
        resultat = resultat.setSubMatrix(c12, c11.rows, 0);
        resultat = resultat.setSubMatrix(c21, 0, c11.columns);
        resultat = resultat.setSubMatrix(c22, c11.rows, c11.columns);
        return resultat.subMatrix(0, rows - 1, 0, cols - 1);
      }

      return blockMult(x, y, r, c);
    }
    /**
     * Returns a row-by-row scaled matrix
     * @param {number} [min=0] - Minimum scaled value
     * @param {number} [max=1] - Maximum scaled value
     * @return {Matrix} - The scaled matrix
     */


    scaleRows(min, max) {
      min = min === undefined ? 0 : min;
      max = max === undefined ? 1 : max;

      if (min >= max) {
        throw new RangeError('min should be strictly smaller than max');
      }

      var newMatrix = this.constructor.empty(this.rows, this.columns);

      for (var i = 0; i < this.rows; i++) {
        var scaled = arrayUtils.scale(this.getRow(i), {
          min,
          max
        });
        newMatrix.setRow(i, scaled);
      }

      return newMatrix;
    }
    /**
     * Returns a new column-by-column scaled matrix
     * @param {number} [min=0] - Minimum scaled value
     * @param {number} [max=1] - Maximum scaled value
     * @return {Matrix} - The new scaled matrix
     * @example
     * var matrix = new Matrix([[1,2],[-1,0]]);
     * var scaledMatrix = matrix.scaleColumns(); // [[1,1],[0,0]]
     */


    scaleColumns(min, max) {
      min = min === undefined ? 0 : min;
      max = max === undefined ? 1 : max;

      if (min >= max) {
        throw new RangeError('min should be strictly smaller than max');
      }

      var newMatrix = this.constructor.empty(this.rows, this.columns);

      for (var i = 0; i < this.columns; i++) {
        var scaled = arrayUtils.scale(this.getColumn(i), {
          min: min,
          max: max
        });
        newMatrix.setColumn(i, scaled);
      }

      return newMatrix;
    }
    /**
     * Returns the Kronecker product (also known as tensor product) between this and other
     * See https://en.wikipedia.org/wiki/Kronecker_product
     * @param {Matrix} other
     * @return {Matrix}
     */


    kroneckerProduct(other) {
      other = this.constructor.checkMatrix(other);
      var m = this.rows;
      var n = this.columns;
      var p = other.rows;
      var q = other.columns;
      var result = new this.constructor[Symbol.species](m * p, n * q);

      for (var i = 0; i < m; i++) {
        for (var j = 0; j < n; j++) {
          for (var k = 0; k < p; k++) {
            for (var l = 0; l < q; l++) {
              result[p * i + k][q * j + l] = this.get(i, j) * other.get(k, l);
            }
          }
        }
      }

      return result;
    }
    /**
     * Transposes the matrix and returns a new one containing the result
     * @return {Matrix}
     */


    transpose() {
      var result = new this.constructor[Symbol.species](this.columns, this.rows);

      for (var i = 0; i < this.rows; i++) {
        for (var j = 0; j < this.columns; j++) {
          result.set(j, i, this.get(i, j));
        }
      }

      return result;
    }
    /**
     * Sorts the rows (in place)
     * @param {function} compareFunction - usual Array.prototype.sort comparison function
     * @return {Matrix} this
     */


    sortRows(compareFunction) {
      if (compareFunction === undefined) compareFunction = compareNumbers;

      for (var i = 0; i < this.rows; i++) {
        this.setRow(i, this.getRow(i).sort(compareFunction));
      }

      return this;
    }
    /**
     * Sorts the columns (in place)
     * @param {function} compareFunction - usual Array.prototype.sort comparison function
     * @return {Matrix} this
     */


    sortColumns(compareFunction) {
      if (compareFunction === undefined) compareFunction = compareNumbers;

      for (var i = 0; i < this.columns; i++) {
        this.setColumn(i, this.getColumn(i).sort(compareFunction));
      }

      return this;
    }
    /**
     * Returns a subset of the matrix
     * @param {number} startRow - First row index
     * @param {number} endRow - Last row index
     * @param {number} startColumn - First column index
     * @param {number} endColumn - Last column index
     * @return {Matrix}
     */


    subMatrix(startRow, endRow, startColumn, endColumn) {
      util.checkRange(this, startRow, endRow, startColumn, endColumn);
      var newMatrix = new this.constructor[Symbol.species](endRow - startRow + 1, endColumn - startColumn + 1);

      for (var i = startRow; i <= endRow; i++) {
        for (var j = startColumn; j <= endColumn; j++) {
          newMatrix[i - startRow][j - startColumn] = this.get(i, j);
        }
      }

      return newMatrix;
    }
    /**
     * Returns a subset of the matrix based on an array of row indices
     * @param {Array} indices - Array containing the row indices
     * @param {number} [startColumn = 0] - First column index
     * @param {number} [endColumn = this.columns-1] - Last column index
     * @return {Matrix}
     */


    subMatrixRow(indices, startColumn, endColumn) {
      if (startColumn === undefined) startColumn = 0;
      if (endColumn === undefined) endColumn = this.columns - 1;

      if (startColumn > endColumn || startColumn < 0 || startColumn >= this.columns || endColumn < 0 || endColumn >= this.columns) {
        throw new RangeError('Argument out of range');
      }

      var newMatrix = new this.constructor[Symbol.species](indices.length, endColumn - startColumn + 1);

      for (var i = 0; i < indices.length; i++) {
        for (var j = startColumn; j <= endColumn; j++) {
          if (indices[i] < 0 || indices[i] >= this.rows) {
            throw new RangeError('Row index out of range: ' + indices[i]);
          }

          newMatrix.set(i, j - startColumn, this.get(indices[i], j));
        }
      }

      return newMatrix;
    }
    /**
     * Returns a subset of the matrix based on an array of column indices
     * @param {Array} indices - Array containing the column indices
     * @param {number} [startRow = 0] - First row index
     * @param {number} [endRow = this.rows-1] - Last row index
     * @return {Matrix}
     */


    subMatrixColumn(indices, startRow, endRow) {
      if (startRow === undefined) startRow = 0;
      if (endRow === undefined) endRow = this.rows - 1;

      if (startRow > endRow || startRow < 0 || startRow >= this.rows || endRow < 0 || endRow >= this.rows) {
        throw new RangeError('Argument out of range');
      }

      var newMatrix = new this.constructor[Symbol.species](endRow - startRow + 1, indices.length);

      for (var i = 0; i < indices.length; i++) {
        for (var j = startRow; j <= endRow; j++) {
          if (indices[i] < 0 || indices[i] >= this.columns) {
            throw new RangeError('Column index out of range: ' + indices[i]);
          }

          newMatrix.set(j - startRow, i, this.get(j, indices[i]));
        }
      }

      return newMatrix;
    }
    /**
     * Set a part of the matrix to the given sub-matrix
     * @param {Matrix|Array< Array >} matrix - The source matrix from which to extract values.
     * @param {number} startRow - The index of the first row to set
     * @param {number} startColumn - The index of the first column to set
     * @return {Matrix}
     */


    setSubMatrix(matrix, startRow, startColumn) {
      matrix = this.constructor.checkMatrix(matrix);
      var endRow = startRow + matrix.rows - 1;
      var endColumn = startColumn + matrix.columns - 1;
      util.checkRange(this, startRow, endRow, startColumn, endColumn);

      for (var i = 0; i < matrix.rows; i++) {
        for (var j = 0; j < matrix.columns; j++) {
          this[startRow + i][startColumn + j] = matrix.get(i, j);
        }
      }

      return this;
    }
    /**
     * Return a new matrix based on a selection of rows and columns
     * @param {Array<number>} rowIndices - The row indices to select. Order matters and an index can be more than once.
     * @param {Array<number>} columnIndices - The column indices to select. Order matters and an index can be use more than once.
     * @return {Matrix} The new matrix
     */


    selection(rowIndices, columnIndices) {
      var indices = util.checkIndices(this, rowIndices, columnIndices);
      var newMatrix = new this.constructor[Symbol.species](rowIndices.length, columnIndices.length);

      for (var i = 0; i < indices.row.length; i++) {
        var rowIndex = indices.row[i];

        for (var j = 0; j < indices.column.length; j++) {
          var columnIndex = indices.column[j];
          newMatrix[i][j] = this.get(rowIndex, columnIndex);
        }
      }

      return newMatrix;
    }
    /**
     * Returns the trace of the matrix (sum of the diagonal elements)
     * @return {number}
     */


    trace() {
      var min = Math.min(this.rows, this.columns);
      var trace = 0;

      for (var i = 0; i < min; i++) {
        trace += this.get(i, i);
      }

      return trace;
    }
    /*
     Matrix views
     */

    /**
     * Returns a view of the transposition of the matrix
     * @return {MatrixTransposeView}
     */


    transposeView() {
      return new MatrixTransposeView(this);
    }
    /**
     * Returns a view of the row vector with the given index
     * @param {number} row - row index of the vector
     * @return {MatrixRowView}
     */


    rowView(row) {
      util.checkRowIndex(this, row);
      return new MatrixRowView(this, row);
    }
    /**
     * Returns a view of the column vector with the given index
     * @param {number} column - column index of the vector
     * @return {MatrixColumnView}
     */


    columnView(column) {
      util.checkColumnIndex(this, column);
      return new MatrixColumnView(this, column);
    }
    /**
     * Returns a view of the matrix flipped in the row axis
     * @return {MatrixFlipRowView}
     */


    flipRowView() {
      return new MatrixFlipRowView(this);
    }
    /**
     * Returns a view of the matrix flipped in the column axis
     * @return {MatrixFlipColumnView}
     */


    flipColumnView() {
      return new MatrixFlipColumnView(this);
    }
    /**
     * Returns a view of a submatrix giving the index boundaries
     * @param {number} startRow - first row index of the submatrix
     * @param {number} endRow - last row index of the submatrix
     * @param {number} startColumn - first column index of the submatrix
     * @param {number} endColumn - last column index of the submatrix
     * @return {MatrixSubView}
     */


    subMatrixView(startRow, endRow, startColumn, endColumn) {
      return new MatrixSubView(this, startRow, endRow, startColumn, endColumn);
    }
    /**
     * Returns a view of the cross of the row indices and the column indices
     * @example
     * // resulting vector is [[2], [2]]
     * var matrix = new Matrix([[1,2,3], [4,5,6]]).selectionView([0, 0], [1])
     * @param {Array<number>} rowIndices
     * @param {Array<number>} columnIndices
     * @return {MatrixSelectionView}
     */


    selectionView(rowIndices, columnIndices) {
      return new MatrixSelectionView(this, rowIndices, columnIndices);
    }
    /**
    * Calculates and returns the determinant of a matrix as a Number
    * @example
    *   new Matrix([[1,2,3], [4,5,6]]).det()
    * @return {number}
    */


    det() {
      if (this.isSquare()) {
        var a, b, c, d;

        if (this.columns === 2) {
          // 2 x 2 matrix
          a = this.get(0, 0);
          b = this.get(0, 1);
          c = this.get(1, 0);
          d = this.get(1, 1);
          return a * d - b * c;
        } else if (this.columns === 3) {
          // 3 x 3 matrix
          var subMatrix0, subMatrix1, subMatrix2;
          subMatrix0 = this.selectionView([1, 2], [1, 2]);
          subMatrix1 = this.selectionView([1, 2], [0, 2]);
          subMatrix2 = this.selectionView([1, 2], [0, 1]);
          a = this.get(0, 0);
          b = this.get(0, 1);
          c = this.get(0, 2);
          return a * subMatrix0.det() - b * subMatrix1.det() + c * subMatrix2.det();
        } else {
          // general purpose determinant using the LU decomposition
          return new LuDecomposition(this).determinant;
        }
      } else {
        throw Error('Determinant can only be calculated for a square matrix.');
      }
    }
    /**
     * Returns inverse of a matrix if it exists or the pseudoinverse
     * @param {number} threshold - threshold for taking inverse of singular values (default = 1e-15)
     * @return {Matrix} the (pseudo)inverted matrix.
     */


    pseudoInverse(threshold) {
      if (threshold === undefined) threshold = Number.EPSILON;
      var svdSolution = new SvDecomposition(this, {
        autoTranspose: true
      });
      var U = svdSolution.leftSingularVectors;
      var V = svdSolution.rightSingularVectors;
      var s = svdSolution.diagonal;

      for (var i = 0; i < s.length; i++) {
        if (Math.abs(s[i]) > threshold) {
          s[i] = 1.0 / s[i];
        } else {
          s[i] = 0.0;
        }
      } // convert list to diagonal


      s = this.constructor[Symbol.species].diag(s);
      return V.mmul(s.mmul(U.transposeView()));
    }

  }

  Matrix.prototype.klass = 'Matrix';
  /**
   * @private
   * Check that two matrices have the same dimensions
   * @param {Matrix} matrix
   * @param {Matrix} otherMatrix
   */

  function checkDimensions(matrix, otherMatrix) {
    // eslint-disable-line no-unused-vars
    if (matrix.rows !== otherMatrix.rows || matrix.columns !== otherMatrix.columns) {
      throw new RangeError('Matrices dimensions must be equal');
    }
  }

  function compareNumbers(a, b) {
    return a - b;
  }
  /*
   Synonyms
   */


  Matrix.random = Matrix.rand;
  Matrix.diagonal = Matrix.diag;
  Matrix.prototype.diagonal = Matrix.prototype.diag;
  Matrix.identity = Matrix.eye;
  Matrix.prototype.negate = Matrix.prototype.neg;
  Matrix.prototype.tensorProduct = Matrix.prototype.kroneckerProduct;
  Matrix.prototype.determinant = Matrix.prototype.det;
  /*
   Add dynamically instance and static methods for mathematical operations
   */

  var inplaceOperator = `
(function %name%(value) {
    if (typeof value === 'number') return this.%name%S(value);
    return this.%name%M(value);
})
`;
  var inplaceOperatorScalar = `
(function %name%S(value) {
    for (var i = 0; i < this.rows; i++) {
        for (var j = 0; j < this.columns; j++) {
            this.set(i, j, this.get(i, j) %op% value);
        }
    }
    return this;
})
`;
  var inplaceOperatorMatrix = `
(function %name%M(matrix) {
    matrix = this.constructor.checkMatrix(matrix);
    checkDimensions(this, matrix);
    for (var i = 0; i < this.rows; i++) {
        for (var j = 0; j < this.columns; j++) {
            this.set(i, j, this.get(i, j) %op% matrix.get(i, j));
        }
    }
    return this;
})
`;
  var staticOperator = `
(function %name%(matrix, value) {
    var newMatrix = new this[Symbol.species](matrix);
    return newMatrix.%name%(value);
})
`;
  var inplaceMethod = `
(function %name%() {
    for (var i = 0; i < this.rows; i++) {
        for (var j = 0; j < this.columns; j++) {
            this.set(i, j, %method%(this.get(i, j)));
        }
    }
    return this;
})
`;
  var staticMethod = `
(function %name%(matrix) {
    var newMatrix = new this[Symbol.species](matrix);
    return newMatrix.%name%();
})
`;
  var inplaceMethodWithArgs = `
(function %name%(%args%) {
    for (var i = 0; i < this.rows; i++) {
        for (var j = 0; j < this.columns; j++) {
            this.set(i, j, %method%(this.get(i, j), %args%));
        }
    }
    return this;
})
`;
  var staticMethodWithArgs = `
(function %name%(matrix, %args%) {
    var newMatrix = new this[Symbol.species](matrix);
    return newMatrix.%name%(%args%);
})
`;
  var inplaceMethodWithOneArgScalar = `
(function %name%S(value) {
    for (var i = 0; i < this.rows; i++) {
        for (var j = 0; j < this.columns; j++) {
            this.set(i, j, %method%(this.get(i, j), value));
        }
    }
    return this;
})
`;
  var inplaceMethodWithOneArgMatrix = `
(function %name%M(matrix) {
    matrix = this.constructor.checkMatrix(matrix);
    checkDimensions(this, matrix);
    for (var i = 0; i < this.rows; i++) {
        for (var j = 0; j < this.columns; j++) {
            this.set(i, j, %method%(this.get(i, j), matrix.get(i, j)));
        }
    }
    return this;
})
`;
  var inplaceMethodWithOneArg = `
(function %name%(value) {
    if (typeof value === 'number') return this.%name%S(value);
    return this.%name%M(value);
})
`;
  var staticMethodWithOneArg = staticMethodWithArgs;
  var operators = [// Arithmetic operators
  ['+', 'add'], ['-', 'sub', 'subtract'], ['*', 'mul', 'multiply'], ['/', 'div', 'divide'], ['%', 'mod', 'modulus'], // Bitwise operators
  ['&', 'and'], ['|', 'or'], ['^', 'xor'], ['<<', 'leftShift'], ['>>', 'signPropagatingRightShift'], ['>>>', 'rightShift', 'zeroFillRightShift']];
  var i;

  for (var operator of operators) {
    var inplaceOp = eval(fillTemplateFunction(inplaceOperator, {
      name: operator[1],
      op: operator[0]
    }));
    var inplaceOpS = eval(fillTemplateFunction(inplaceOperatorScalar, {
      name: operator[1] + 'S',
      op: operator[0]
    }));
    var inplaceOpM = eval(fillTemplateFunction(inplaceOperatorMatrix, {
      name: operator[1] + 'M',
      op: operator[0]
    }));
    var staticOp = eval(fillTemplateFunction(staticOperator, {
      name: operator[1]
    }));

    for (i = 1; i < operator.length; i++) {
      Matrix.prototype[operator[i]] = inplaceOp;
      Matrix.prototype[operator[i] + 'S'] = inplaceOpS;
      Matrix.prototype[operator[i] + 'M'] = inplaceOpM;
      Matrix[operator[i]] = staticOp;
    }
  }

  var methods = [['~', 'not']];
  ['abs', 'acos', 'acosh', 'asin', 'asinh', 'atan', 'atanh', 'cbrt', 'ceil', 'clz32', 'cos', 'cosh', 'exp', 'expm1', 'floor', 'fround', 'log', 'log1p', 'log10', 'log2', 'round', 'sign', 'sin', 'sinh', 'sqrt', 'tan', 'tanh', 'trunc'].forEach(function (mathMethod) {
    methods.push(['Math.' + mathMethod, mathMethod]);
  });

  for (var method of methods) {
    var inplaceMeth = eval(fillTemplateFunction(inplaceMethod, {
      name: method[1],
      method: method[0]
    }));
    var staticMeth = eval(fillTemplateFunction(staticMethod, {
      name: method[1]
    }));

    for (i = 1; i < method.length; i++) {
      Matrix.prototype[method[i]] = inplaceMeth;
      Matrix[method[i]] = staticMeth;
    }
  }

  var methodsWithArgs = [['Math.pow', 1, 'pow']];

  for (var methodWithArg of methodsWithArgs) {
    var args = 'arg0';

    for (i = 1; i < methodWithArg[1]; i++) {
      args += `, arg${i}`;
    }

    if (methodWithArg[1] !== 1) {
      var inplaceMethWithArgs = eval(fillTemplateFunction(inplaceMethodWithArgs, {
        name: methodWithArg[2],
        method: methodWithArg[0],
        args: args
      }));
      var staticMethWithArgs = eval(fillTemplateFunction(staticMethodWithArgs, {
        name: methodWithArg[2],
        args: args
      }));

      for (i = 2; i < methodWithArg.length; i++) {
        Matrix.prototype[methodWithArg[i]] = inplaceMethWithArgs;
        Matrix[methodWithArg[i]] = staticMethWithArgs;
      }
    } else {
      var tmplVar = {
        name: methodWithArg[2],
        args: args,
        method: methodWithArg[0]
      };
      var inplaceMethod2 = eval(fillTemplateFunction(inplaceMethodWithOneArg, tmplVar));
      var inplaceMethodS = eval(fillTemplateFunction(inplaceMethodWithOneArgScalar, tmplVar));
      var inplaceMethodM = eval(fillTemplateFunction(inplaceMethodWithOneArgMatrix, tmplVar));
      var staticMethod2 = eval(fillTemplateFunction(staticMethodWithOneArg, tmplVar));

      for (i = 2; i < methodWithArg.length; i++) {
        Matrix.prototype[methodWithArg[i]] = inplaceMethod2;
        Matrix.prototype[methodWithArg[i] + 'M'] = inplaceMethodM;
        Matrix.prototype[methodWithArg[i] + 'S'] = inplaceMethodS;
        Matrix[methodWithArg[i]] = staticMethod2;
      }
    }
  }

  function fillTemplateFunction(template, values) {
    for (var value in values) {
      template = template.replace(new RegExp('%' + value + '%', 'g'), values[value]);
    }

    return template;
  }

  return Matrix;
}

/***/ }),
/* 9 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


var Matrix = __webpack_require__(0); // https://github.com/lutzroeder/Mapack/blob/master/Source/LuDecomposition.cs


function LuDecomposition(matrix) {
  if (!(this instanceof LuDecomposition)) {
    return new LuDecomposition(matrix);
  }

  matrix = Matrix.Matrix.checkMatrix(matrix);
  var lu = matrix.clone(),
      rows = lu.rows,
      columns = lu.columns,
      pivotVector = new Array(rows),
      pivotSign = 1,
      i,
      j,
      k,
      p,
      s,
      t,
      v,
      LUrowi,
      LUcolj,
      kmax;

  for (i = 0; i < rows; i++) {
    pivotVector[i] = i;
  }

  LUcolj = new Array(rows);

  for (j = 0; j < columns; j++) {
    for (i = 0; i < rows; i++) {
      LUcolj[i] = lu[i][j];
    }

    for (i = 0; i < rows; i++) {
      LUrowi = lu[i];
      kmax = Math.min(i, j);
      s = 0;

      for (k = 0; k < kmax; k++) {
        s += LUrowi[k] * LUcolj[k];
      }

      LUrowi[j] = LUcolj[i] -= s;
    }

    p = j;

    for (i = j + 1; i < rows; i++) {
      if (Math.abs(LUcolj[i]) > Math.abs(LUcolj[p])) {
        p = i;
      }
    }

    if (p !== j) {
      for (k = 0; k < columns; k++) {
        t = lu[p][k];
        lu[p][k] = lu[j][k];
        lu[j][k] = t;
      }

      v = pivotVector[p];
      pivotVector[p] = pivotVector[j];
      pivotVector[j] = v;
      pivotSign = -pivotSign;
    }

    if (j < rows && lu[j][j] !== 0) {
      for (i = j + 1; i < rows; i++) {
        lu[i][j] /= lu[j][j];
      }
    }
  }

  this.LU = lu;
  this.pivotVector = pivotVector;
  this.pivotSign = pivotSign;
}

LuDecomposition.prototype = {
  isSingular: function isSingular() {
    var data = this.LU,
        col = data.columns;

    for (var j = 0; j < col; j++) {
      if (data[j][j] === 0) {
        return true;
      }
    }

    return false;
  },

  get determinant() {
    var data = this.LU;

    if (!data.isSquare()) {
      throw new Error('Matrix must be square');
    }

    var determinant = this.pivotSign,
        col = data.columns;

    for (var j = 0; j < col; j++) {
      determinant *= data[j][j];
    }

    return determinant;
  },

  get lowerTriangularMatrix() {
    var data = this.LU,
        rows = data.rows,
        columns = data.columns,
        X = new Matrix.Matrix(rows, columns);

    for (var i = 0; i < rows; i++) {
      for (var j = 0; j < columns; j++) {
        if (i > j) {
          X[i][j] = data[i][j];
        } else if (i === j) {
          X[i][j] = 1;
        } else {
          X[i][j] = 0;
        }
      }
    }

    return X;
  },

  get upperTriangularMatrix() {
    var data = this.LU,
        rows = data.rows,
        columns = data.columns,
        X = new Matrix.Matrix(rows, columns);

    for (var i = 0; i < rows; i++) {
      for (var j = 0; j < columns; j++) {
        if (i <= j) {
          X[i][j] = data[i][j];
        } else {
          X[i][j] = 0;
        }
      }
    }

    return X;
  },

  get pivotPermutationVector() {
    return this.pivotVector.slice();
  },

  solve: function solve(value) {
    value = Matrix.Matrix.checkMatrix(value);
    var lu = this.LU,
        rows = lu.rows;

    if (rows !== value.rows) {
      throw new Error('Invalid matrix dimensions');
    }

    if (this.isSingular()) {
      throw new Error('LU matrix is singular');
    }

    var count = value.columns;
    var X = value.subMatrixRow(this.pivotVector, 0, count - 1);
    var columns = lu.columns;
    var i, j, k;

    for (k = 0; k < columns; k++) {
      for (i = k + 1; i < columns; i++) {
        for (j = 0; j < count; j++) {
          X[i][j] -= X[k][j] * lu[i][k];
        }
      }
    }

    for (k = columns - 1; k >= 0; k--) {
      for (j = 0; j < count; j++) {
        X[k][j] /= lu[k][k];
      }

      for (i = 0; i < k; i++) {
        for (j = 0; j < count; j++) {
          X[i][j] -= X[k][j] * lu[i][k];
        }
      }
    }

    return X;
  }
};
module.exports = LuDecomposition;

/***/ }),
/* 10 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


var Matrix = __webpack_require__(0);

var util = __webpack_require__(5);

var hypotenuse = util.hypotenuse;
var getFilled2DArray = util.getFilled2DArray; // https://github.com/lutzroeder/Mapack/blob/master/Source/SingularValueDecomposition.cs

function SingularValueDecomposition(value, options) {
  if (!(this instanceof SingularValueDecomposition)) {
    return new SingularValueDecomposition(value, options);
  }

  value = Matrix.Matrix.checkMatrix(value);
  options = options || {};
  var m = value.rows,
      n = value.columns,
      nu = Math.min(m, n);
  var wantu = true,
      wantv = true;
  if (options.computeLeftSingularVectors === false) wantu = false;
  if (options.computeRightSingularVectors === false) wantv = false;
  var autoTranspose = options.autoTranspose === true;
  var swapped = false;
  var a;

  if (m < n) {
    if (!autoTranspose) {
      a = value.clone(); // eslint-disable-next-line no-console

      console.warn('Computing SVD on a matrix with more columns than rows. Consider enabling autoTranspose');
    } else {
      a = value.transpose();
      m = a.rows;
      n = a.columns;
      swapped = true;
      var aux = wantu;
      wantu = wantv;
      wantv = aux;
    }
  } else {
    a = value.clone();
  }

  var s = new Array(Math.min(m + 1, n)),
      U = getFilled2DArray(m, nu, 0),
      V = getFilled2DArray(n, n, 0),
      e = new Array(n),
      work = new Array(m);
  var nct = Math.min(m - 1, n);
  var nrt = Math.max(0, Math.min(n - 2, m));
  var i, j, k, p, t, ks, f, cs, sn, max, kase, scale, sp, spm1, epm1, sk, ek, b, c, shift, g;

  for (k = 0, max = Math.max(nct, nrt); k < max; k++) {
    if (k < nct) {
      s[k] = 0;

      for (i = k; i < m; i++) {
        s[k] = hypotenuse(s[k], a[i][k]);
      }

      if (s[k] !== 0) {
        if (a[k][k] < 0) {
          s[k] = -s[k];
        }

        for (i = k; i < m; i++) {
          a[i][k] /= s[k];
        }

        a[k][k] += 1;
      }

      s[k] = -s[k];
    }

    for (j = k + 1; j < n; j++) {
      if (k < nct && s[k] !== 0) {
        t = 0;

        for (i = k; i < m; i++) {
          t += a[i][k] * a[i][j];
        }

        t = -t / a[k][k];

        for (i = k; i < m; i++) {
          a[i][j] += t * a[i][k];
        }
      }

      e[j] = a[k][j];
    }

    if (wantu && k < nct) {
      for (i = k; i < m; i++) {
        U[i][k] = a[i][k];
      }
    }

    if (k < nrt) {
      e[k] = 0;

      for (i = k + 1; i < n; i++) {
        e[k] = hypotenuse(e[k], e[i]);
      }

      if (e[k] !== 0) {
        if (e[k + 1] < 0) {
          e[k] = 0 - e[k];
        }

        for (i = k + 1; i < n; i++) {
          e[i] /= e[k];
        }

        e[k + 1] += 1;
      }

      e[k] = -e[k];

      if (k + 1 < m && e[k] !== 0) {
        for (i = k + 1; i < m; i++) {
          work[i] = 0;
        }

        for (j = k + 1; j < n; j++) {
          for (i = k + 1; i < m; i++) {
            work[i] += e[j] * a[i][j];
          }
        }

        for (j = k + 1; j < n; j++) {
          t = -e[j] / e[k + 1];

          for (i = k + 1; i < m; i++) {
            a[i][j] += t * work[i];
          }
        }
      }

      if (wantv) {
        for (i = k + 1; i < n; i++) {
          V[i][k] = e[i];
        }
      }
    }
  }

  p = Math.min(n, m + 1);

  if (nct < n) {
    s[nct] = a[nct][nct];
  }

  if (m < p) {
    s[p - 1] = 0;
  }

  if (nrt + 1 < p) {
    e[nrt] = a[nrt][p - 1];
  }

  e[p - 1] = 0;

  if (wantu) {
    for (j = nct; j < nu; j++) {
      for (i = 0; i < m; i++) {
        U[i][j] = 0;
      }

      U[j][j] = 1;
    }

    for (k = nct - 1; k >= 0; k--) {
      if (s[k] !== 0) {
        for (j = k + 1; j < nu; j++) {
          t = 0;

          for (i = k; i < m; i++) {
            t += U[i][k] * U[i][j];
          }

          t = -t / U[k][k];

          for (i = k; i < m; i++) {
            U[i][j] += t * U[i][k];
          }
        }

        for (i = k; i < m; i++) {
          U[i][k] = -U[i][k];
        }

        U[k][k] = 1 + U[k][k];

        for (i = 0; i < k - 1; i++) {
          U[i][k] = 0;
        }
      } else {
        for (i = 0; i < m; i++) {
          U[i][k] = 0;
        }

        U[k][k] = 1;
      }
    }
  }

  if (wantv) {
    for (k = n - 1; k >= 0; k--) {
      if (k < nrt && e[k] !== 0) {
        for (j = k + 1; j < n; j++) {
          t = 0;

          for (i = k + 1; i < n; i++) {
            t += V[i][k] * V[i][j];
          }

          t = -t / V[k + 1][k];

          for (i = k + 1; i < n; i++) {
            V[i][j] += t * V[i][k];
          }
        }
      }

      for (i = 0; i < n; i++) {
        V[i][k] = 0;
      }

      V[k][k] = 1;
    }
  }

  var pp = p - 1,
      iter = 0,
      eps = Math.pow(2, -52);

  while (p > 0) {
    for (k = p - 2; k >= -1; k--) {
      if (k === -1) {
        break;
      }

      if (Math.abs(e[k]) <= eps * (Math.abs(s[k]) + Math.abs(s[k + 1]))) {
        e[k] = 0;
        break;
      }
    }

    if (k === p - 2) {
      kase = 4;
    } else {
      for (ks = p - 1; ks >= k; ks--) {
        if (ks === k) {
          break;
        }

        t = (ks !== p ? Math.abs(e[ks]) : 0) + (ks !== k + 1 ? Math.abs(e[ks - 1]) : 0);

        if (Math.abs(s[ks]) <= eps * t) {
          s[ks] = 0;
          break;
        }
      }

      if (ks === k) {
        kase = 3;
      } else if (ks === p - 1) {
        kase = 1;
      } else {
        kase = 2;
        k = ks;
      }
    }

    k++;

    switch (kase) {
      case 1:
        {
          f = e[p - 2];
          e[p - 2] = 0;

          for (j = p - 2; j >= k; j--) {
            t = hypotenuse(s[j], f);
            cs = s[j] / t;
            sn = f / t;
            s[j] = t;

            if (j !== k) {
              f = -sn * e[j - 1];
              e[j - 1] = cs * e[j - 1];
            }

            if (wantv) {
              for (i = 0; i < n; i++) {
                t = cs * V[i][j] + sn * V[i][p - 1];
                V[i][p - 1] = -sn * V[i][j] + cs * V[i][p - 1];
                V[i][j] = t;
              }
            }
          }

          break;
        }

      case 2:
        {
          f = e[k - 1];
          e[k - 1] = 0;

          for (j = k; j < p; j++) {
            t = hypotenuse(s[j], f);
            cs = s[j] / t;
            sn = f / t;
            s[j] = t;
            f = -sn * e[j];
            e[j] = cs * e[j];

            if (wantu) {
              for (i = 0; i < m; i++) {
                t = cs * U[i][j] + sn * U[i][k - 1];
                U[i][k - 1] = -sn * U[i][j] + cs * U[i][k - 1];
                U[i][j] = t;
              }
            }
          }

          break;
        }

      case 3:
        {
          scale = Math.max(Math.max(Math.max(Math.max(Math.abs(s[p - 1]), Math.abs(s[p - 2])), Math.abs(e[p - 2])), Math.abs(s[k])), Math.abs(e[k]));
          sp = s[p - 1] / scale;
          spm1 = s[p - 2] / scale;
          epm1 = e[p - 2] / scale;
          sk = s[k] / scale;
          ek = e[k] / scale;
          b = ((spm1 + sp) * (spm1 - sp) + epm1 * epm1) / 2;
          c = sp * epm1 * (sp * epm1);
          shift = 0;

          if (b !== 0 || c !== 0) {
            shift = Math.sqrt(b * b + c);

            if (b < 0) {
              shift = -shift;
            }

            shift = c / (b + shift);
          }

          f = (sk + sp) * (sk - sp) + shift;
          g = sk * ek;

          for (j = k; j < p - 1; j++) {
            t = hypotenuse(f, g);
            cs = f / t;
            sn = g / t;

            if (j !== k) {
              e[j - 1] = t;
            }

            f = cs * s[j] + sn * e[j];
            e[j] = cs * e[j] - sn * s[j];
            g = sn * s[j + 1];
            s[j + 1] = cs * s[j + 1];

            if (wantv) {
              for (i = 0; i < n; i++) {
                t = cs * V[i][j] + sn * V[i][j + 1];
                V[i][j + 1] = -sn * V[i][j] + cs * V[i][j + 1];
                V[i][j] = t;
              }
            }

            t = hypotenuse(f, g);
            cs = f / t;
            sn = g / t;
            s[j] = t;
            f = cs * e[j] + sn * s[j + 1];
            s[j + 1] = -sn * e[j] + cs * s[j + 1];
            g = sn * e[j + 1];
            e[j + 1] = cs * e[j + 1];

            if (wantu && j < m - 1) {
              for (i = 0; i < m; i++) {
                t = cs * U[i][j] + sn * U[i][j + 1];
                U[i][j + 1] = -sn * U[i][j] + cs * U[i][j + 1];
                U[i][j] = t;
              }
            }
          }

          e[p - 2] = f;
          iter = iter + 1;
          break;
        }

      case 4:
        {
          if (s[k] <= 0) {
            s[k] = s[k] < 0 ? -s[k] : 0;

            if (wantv) {
              for (i = 0; i <= pp; i++) {
                V[i][k] = -V[i][k];
              }
            }
          }

          while (k < pp) {
            if (s[k] >= s[k + 1]) {
              break;
            }

            t = s[k];
            s[k] = s[k + 1];
            s[k + 1] = t;

            if (wantv && k < n - 1) {
              for (i = 0; i < n; i++) {
                t = V[i][k + 1];
                V[i][k + 1] = V[i][k];
                V[i][k] = t;
              }
            }

            if (wantu && k < m - 1) {
              for (i = 0; i < m; i++) {
                t = U[i][k + 1];
                U[i][k + 1] = U[i][k];
                U[i][k] = t;
              }
            }

            k++;
          }

          iter = 0;
          p--;
          break;
        }
      // no default
    }
  }

  if (swapped) {
    var tmp = V;
    V = U;
    U = tmp;
  }

  this.m = m;
  this.n = n;
  this.s = s;
  this.U = U;
  this.V = V;
}

SingularValueDecomposition.prototype = {
  get condition() {
    return this.s[0] / this.s[Math.min(this.m, this.n) - 1];
  },

  get norm2() {
    return this.s[0];
  },

  get rank() {
    var eps = Math.pow(2, -52),
        tol = Math.max(this.m, this.n) * this.s[0] * eps,
        r = 0,
        s = this.s;

    for (var i = 0, ii = s.length; i < ii; i++) {
      if (s[i] > tol) {
        r++;
      }
    }

    return r;
  },

  get diagonal() {
    return this.s;
  },

  // https://github.com/accord-net/framework/blob/development/Sources/Accord.Math/Decompositions/SingularValueDecomposition.cs
  get threshold() {
    return Math.pow(2, -52) / 2 * Math.max(this.m, this.n) * this.s[0];
  },

  get leftSingularVectors() {
    if (!Matrix.Matrix.isMatrix(this.U)) {
      this.U = new Matrix.Matrix(this.U);
    }

    return this.U;
  },

  get rightSingularVectors() {
    if (!Matrix.Matrix.isMatrix(this.V)) {
      this.V = new Matrix.Matrix(this.V);
    }

    return this.V;
  },

  get diagonalMatrix() {
    return Matrix.Matrix.diag(this.s);
  },

  solve: function solve(value) {
    var Y = value,
        e = this.threshold,
        scols = this.s.length,
        Ls = Matrix.Matrix.zeros(scols, scols),
        i;

    for (i = 0; i < scols; i++) {
      if (Math.abs(this.s[i]) <= e) {
        Ls[i][i] = 0;
      } else {
        Ls[i][i] = 1 / this.s[i];
      }
    }

    var U = this.U;
    var V = this.rightSingularVectors;
    var VL = V.mmul(Ls),
        vrows = V.rows,
        urows = U.length,
        VLU = Matrix.Matrix.zeros(vrows, urows),
        j,
        k,
        sum;

    for (i = 0; i < vrows; i++) {
      for (j = 0; j < urows; j++) {
        sum = 0;

        for (k = 0; k < scols; k++) {
          sum += VL[i][k] * U[j][k];
        }

        VLU[i][j] = sum;
      }
    }

    return VLU.mmul(Y);
  },
  solveForDiagonal: function solveForDiagonal(value) {
    return this.solve(Matrix.Matrix.diag(value));
  },
  inverse: function inverse() {
    var V = this.V;
    var e = this.threshold,
        vrows = V.length,
        vcols = V[0].length,
        X = new Matrix.Matrix(vrows, this.s.length),
        i,
        j;

    for (i = 0; i < vrows; i++) {
      for (j = 0; j < vcols; j++) {
        if (Math.abs(this.s[j]) > e) {
          X[i][j] = V[i][j] / this.s[j];
        } else {
          X[i][j] = 0;
        }
      }
    }

    var U = this.U;
    var urows = U.length,
        ucols = U[0].length,
        Y = new Matrix.Matrix(vrows, urows),
        k,
        sum;

    for (i = 0; i < vrows; i++) {
      for (j = 0; j < urows; j++) {
        sum = 0;

        for (k = 0; k < ucols; k++) {
          sum += X[i][k] * U[j][k];
        }

        Y[i][j] = sum;
      }
    }

    return Y;
  }
};
module.exports = SingularValueDecomposition;

/***/ }),
/* 11 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


exports.array = __webpack_require__(12);
exports.matrix = __webpack_require__(39);

/***/ }),
/* 12 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


function compareNumbers(a, b) {
  return a - b;
}
/**
 * Computes the sum of the given values
 * @param {Array} values
 * @returns {number}
 */


exports.sum = function sum(values) {
  var sum = 0;

  for (var i = 0; i < values.length; i++) {
    sum += values[i];
  }

  return sum;
};
/**
 * Computes the maximum of the given values
 * @param {Array} values
 * @returns {number}
 */


exports.max = function max(values) {
  var max = values[0];
  var l = values.length;

  for (var i = 1; i < l; i++) {
    if (values[i] > max) max = values[i];
  }

  return max;
};
/**
 * Computes the minimum of the given values
 * @param {Array} values
 * @returns {number}
 */


exports.min = function min(values) {
  var min = values[0];
  var l = values.length;

  for (var i = 1; i < l; i++) {
    if (values[i] < min) min = values[i];
  }

  return min;
};
/**
 * Computes the min and max of the given values
 * @param {Array} values
 * @returns {{min: number, max: number}}
 */


exports.minMax = function minMax(values) {
  var min = values[0];
  var max = values[0];
  var l = values.length;

  for (var i = 1; i < l; i++) {
    if (values[i] < min) min = values[i];
    if (values[i] > max) max = values[i];
  }

  return {
    min: min,
    max: max
  };
};
/**
 * Computes the arithmetic mean of the given values
 * @param {Array} values
 * @returns {number}
 */


exports.arithmeticMean = function arithmeticMean(values) {
  var sum = 0;
  var l = values.length;

  for (var i = 0; i < l; i++) {
    sum += values[i];
  }

  return sum / l;
};
/**
 * {@link arithmeticMean}
 */


exports.mean = exports.arithmeticMean;
/**
 * Computes the geometric mean of the given values
 * @param {Array} values
 * @returns {number}
 */

exports.geometricMean = function geometricMean(values) {
  var mul = 1;
  var l = values.length;

  for (var i = 0; i < l; i++) {
    mul *= values[i];
  }

  return Math.pow(mul, 1 / l);
};
/**
 * Computes the mean of the log of the given values
 * If the return value is exponentiated, it gives the same result as the
 * geometric mean.
 * @param {Array} values
 * @returns {number}
 */


exports.logMean = function logMean(values) {
  var lnsum = 0;
  var l = values.length;

  for (var i = 0; i < l; i++) {
    lnsum += Math.log(values[i]);
  }

  return lnsum / l;
};
/**
 * Computes the weighted grand mean for a list of means and sample sizes
 * @param {Array} means - Mean values for each set of samples
 * @param {Array} samples - Number of original values for each set of samples
 * @returns {number}
 */


exports.grandMean = function grandMean(means, samples) {
  var sum = 0;
  var n = 0;
  var l = means.length;

  for (var i = 0; i < l; i++) {
    sum += samples[i] * means[i];
    n += samples[i];
  }

  return sum / n;
};
/**
 * Computes the truncated mean of the given values using a given percentage
 * @param {Array} values
 * @param {number} percent - The percentage of values to keep (range: [0,1])
 * @param {boolean} [alreadySorted=false]
 * @returns {number}
 */


exports.truncatedMean = function truncatedMean(values, percent, alreadySorted) {
  if (alreadySorted === undefined) alreadySorted = false;

  if (!alreadySorted) {
    values = [].concat(values).sort(compareNumbers);
  }

  var l = values.length;
  var k = Math.floor(l * percent);
  var sum = 0;

  for (var i = k; i < l - k; i++) {
    sum += values[i];
  }

  return sum / (l - 2 * k);
};
/**
 * Computes the harmonic mean of the given values
 * @param {Array} values
 * @returns {number}
 */


exports.harmonicMean = function harmonicMean(values) {
  var sum = 0;
  var l = values.length;

  for (var i = 0; i < l; i++) {
    if (values[i] === 0) {
      throw new RangeError('value at index ' + i + 'is zero');
    }

    sum += 1 / values[i];
  }

  return l / sum;
};
/**
 * Computes the contraharmonic mean of the given values
 * @param {Array} values
 * @returns {number}
 */


exports.contraHarmonicMean = function contraHarmonicMean(values) {
  var r1 = 0;
  var r2 = 0;
  var l = values.length;

  for (var i = 0; i < l; i++) {
    r1 += values[i] * values[i];
    r2 += values[i];
  }

  if (r2 < 0) {
    throw new RangeError('sum of values is negative');
  }

  return r1 / r2;
};
/**
 * Computes the median of the given values
 * @param {Array} values
 * @param {boolean} [alreadySorted=false]
 * @returns {number}
 */


exports.median = function median(values, alreadySorted) {
  if (alreadySorted === undefined) alreadySorted = false;

  if (!alreadySorted) {
    values = [].concat(values).sort(compareNumbers);
  }

  var l = values.length;
  var half = Math.floor(l / 2);

  if (l % 2 === 0) {
    return (values[half - 1] + values[half]) * 0.5;
  } else {
    return values[half];
  }
};
/**
 * Computes the variance of the given values
 * @param {Array} values
 * @param {boolean} [unbiased=true] - if true, divide by (n-1); if false, divide by n.
 * @returns {number}
 */


exports.variance = function variance(values, unbiased) {
  if (unbiased === undefined) unbiased = true;
  var theMean = exports.mean(values);
  var theVariance = 0;
  var l = values.length;

  for (var i = 0; i < l; i++) {
    var x = values[i] - theMean;
    theVariance += x * x;
  }

  if (unbiased) {
    return theVariance / (l - 1);
  } else {
    return theVariance / l;
  }
};
/**
 * Computes the standard deviation of the given values
 * @param {Array} values
 * @param {boolean} [unbiased=true] - if true, divide by (n-1); if false, divide by n.
 * @returns {number}
 */


exports.standardDeviation = function standardDeviation(values, unbiased) {
  return Math.sqrt(exports.variance(values, unbiased));
};

exports.standardError = function standardError(values) {
  return exports.standardDeviation(values) / Math.sqrt(values.length);
};
/**
 * IEEE Transactions on biomedical engineering, vol. 52, no. 1, january 2005, p. 76-
 * Calculate the standard deviation via the Median of the absolute deviation
 *  The formula for the standard deviation only holds for Gaussian random variables.
 * @returns {{mean: number, stdev: number}}
 */


exports.robustMeanAndStdev = function robustMeanAndStdev(y) {
  var mean = 0,
      stdev = 0;
  var length = y.length,
      i = 0;

  for (i = 0; i < length; i++) {
    mean += y[i];
  }

  mean /= length;
  var averageDeviations = new Array(length);

  for (i = 0; i < length; i++) averageDeviations[i] = Math.abs(y[i] - mean);

  averageDeviations.sort(compareNumbers);

  if (length % 2 === 1) {
    stdev = averageDeviations[(length - 1) / 2] / 0.6745;
  } else {
    stdev = 0.5 * (averageDeviations[length / 2] + averageDeviations[length / 2 - 1]) / 0.6745;
  }

  return {
    mean: mean,
    stdev: stdev
  };
};

exports.quartiles = function quartiles(values, alreadySorted) {
  if (typeof alreadySorted === 'undefined') alreadySorted = false;

  if (!alreadySorted) {
    values = [].concat(values).sort(compareNumbers);
  }

  var quart = values.length / 4;
  var q1 = values[Math.ceil(quart) - 1];
  var q2 = exports.median(values, true);
  var q3 = values[Math.ceil(quart * 3) - 1];
  return {
    q1: q1,
    q2: q2,
    q3: q3
  };
};

exports.pooledStandardDeviation = function pooledStandardDeviation(samples, unbiased) {
  return Math.sqrt(exports.pooledVariance(samples, unbiased));
};

exports.pooledVariance = function pooledVariance(samples, unbiased) {
  if (typeof unbiased === 'undefined') unbiased = true;
  var sum = 0;
  var length = 0,
      l = samples.length;

  for (var i = 0; i < l; i++) {
    var values = samples[i];
    var vari = exports.variance(values);
    sum += (values.length - 1) * vari;
    if (unbiased) length += values.length - 1;else length += values.length;
  }

  return sum / length;
};

exports.mode = function mode(values) {
  var l = values.length,
      itemCount = new Array(l),
      i;

  for (i = 0; i < l; i++) {
    itemCount[i] = 0;
  }

  var itemArray = new Array(l);
  var count = 0;

  for (i = 0; i < l; i++) {
    var index = itemArray.indexOf(values[i]);
    if (index >= 0) itemCount[index]++;else {
      itemArray[count] = values[i];
      itemCount[count] = 1;
      count++;
    }
  }

  var maxValue = 0,
      maxIndex = 0;

  for (i = 0; i < count; i++) {
    if (itemCount[i] > maxValue) {
      maxValue = itemCount[i];
      maxIndex = i;
    }
  }

  return itemArray[maxIndex];
};

exports.covariance = function covariance(vector1, vector2, unbiased) {
  if (typeof unbiased === 'undefined') unbiased = true;
  var mean1 = exports.mean(vector1);
  var mean2 = exports.mean(vector2);
  if (vector1.length !== vector2.length) throw 'Vectors do not have the same dimensions';
  var cov = 0,
      l = vector1.length;

  for (var i = 0; i < l; i++) {
    var x = vector1[i] - mean1;
    var y = vector2[i] - mean2;
    cov += x * y;
  }

  if (unbiased) return cov / (l - 1);else return cov / l;
};

exports.skewness = function skewness(values, unbiased) {
  if (typeof unbiased === 'undefined') unbiased = true;
  var theMean = exports.mean(values);
  var s2 = 0,
      s3 = 0,
      l = values.length;

  for (var i = 0; i < l; i++) {
    var dev = values[i] - theMean;
    s2 += dev * dev;
    s3 += dev * dev * dev;
  }

  var m2 = s2 / l;
  var m3 = s3 / l;
  var g = m3 / Math.pow(m2, 3 / 2.0);

  if (unbiased) {
    var a = Math.sqrt(l * (l - 1));
    var b = l - 2;
    return a / b * g;
  } else {
    return g;
  }
};

exports.kurtosis = function kurtosis(values, unbiased) {
  if (typeof unbiased === 'undefined') unbiased = true;
  var theMean = exports.mean(values);
  var n = values.length,
      s2 = 0,
      s4 = 0;

  for (var i = 0; i < n; i++) {
    var dev = values[i] - theMean;
    s2 += dev * dev;
    s4 += dev * dev * dev * dev;
  }

  var m2 = s2 / n;
  var m4 = s4 / n;

  if (unbiased) {
    var v = s2 / (n - 1);
    var a = n * (n + 1) / ((n - 1) * (n - 2) * (n - 3));
    var b = s4 / (v * v);
    var c = (n - 1) * (n - 1) / ((n - 2) * (n - 3));
    return a * b - 3 * c;
  } else {
    return m4 / (m2 * m2) - 3;
  }
};

exports.entropy = function entropy(values, eps) {
  if (typeof eps === 'undefined') eps = 0;
  var sum = 0,
      l = values.length;

  for (var i = 0; i < l; i++) sum += values[i] * Math.log(values[i] + eps);

  return -sum;
};

exports.weightedMean = function weightedMean(values, weights) {
  var sum = 0,
      l = values.length;

  for (var i = 0; i < l; i++) sum += values[i] * weights[i];

  return sum;
};

exports.weightedStandardDeviation = function weightedStandardDeviation(values, weights) {
  return Math.sqrt(exports.weightedVariance(values, weights));
};

exports.weightedVariance = function weightedVariance(values, weights) {
  var theMean = exports.weightedMean(values, weights);
  var vari = 0,
      l = values.length;
  var a = 0,
      b = 0;

  for (var i = 0; i < l; i++) {
    var z = values[i] - theMean;
    var w = weights[i];
    vari += w * (z * z);
    b += w;
    a += w * w;
  }

  return vari * (b / (b * b - a));
};

exports.center = function center(values, inPlace) {
  if (typeof inPlace === 'undefined') inPlace = false;
  var result = values;
  if (!inPlace) result = [].concat(values);
  var theMean = exports.mean(result),
      l = result.length;

  for (var i = 0; i < l; i++) result[i] -= theMean;
};

exports.standardize = function standardize(values, standardDev, inPlace) {
  if (typeof standardDev === 'undefined') standardDev = exports.standardDeviation(values);
  if (typeof inPlace === 'undefined') inPlace = false;
  var l = values.length;
  var result = inPlace ? values : new Array(l);

  for (var i = 0; i < l; i++) result[i] = values[i] / standardDev;

  return result;
};

exports.cumulativeSum = function cumulativeSum(array) {
  var l = array.length;
  var result = new Array(l);
  result[0] = array[0];

  for (var i = 1; i < l; i++) result[i] = result[i - 1] + array[i];

  return result;
};

/***/ }),
/* 13 */
/***/ (function(module, exports, __webpack_require__) {

/**
 * Created by acastillo on 7/5/16.
 */
const SpinSystem = __webpack_require__(14);

const AutoAssigner = __webpack_require__(15);

const OCLE = __webpack_require__(17);

const DEBUG = false;

function autoAssign(entry, options) {
  if (entry.spectra.h1PeakList) {
    return assignmentFromPeakPicking(entry, options);
  } else {
    return assignmentFromRaw(entry, options);
  }
}

function assignmentFromRaw(entry, options) {//TODO Implement this method

  /*var molfile = entry.molfile;
  var spectra = entry.spectra;
   var molecule =  OCLE.Molecule.fromMolfile(molfile);
   molecule.addImplicitHydrogens();
   entry.molecule = molecule;
  entry.diaIDs = molecule.getGroupedDiastereotopicAtomIDs();
   //Simulate and process the 1H-NMR spectrum at 400MHz
  var jcampFile = molFiles[i].replace("mol_","h1_").replace(".mol",".jdx");
  var spectraData1H = SD.load(spectra.h1);//
    var signals = spectraData1H.nmrPeakDetection({nStddev:3, baselineRejoin:5, compute:false});
  spectra.solvent = spectraData1H.getParamString(".SOLVENT NAME", "unknown");
  entry.diaID = molecule.toIDCode();
   signals = integration(signals, molecule.countAtom("H"));
   for(var j=0;j< signals.length;j++){
      signals[j]._highlight=[-(j+1)];
  }
   spectra.h1PeakList = signals;
   return assignmentFromPeakPicking(entry,options);
  */
}

function assignmentFromPeakPicking(entry, options) {
  const predictor = options.predictor;
  var molecule, diaIDs, molfile;
  var spectra = entry.spectra;

  if (!entry.molecule) {
    molecule = OCLE.Molecule.fromMolfile(entry.molfile);
    molecule.addImplicitHydrogens();
    diaIDs = molecule.getGroupedDiastereotopicAtomIDs();

    for (var j = 0; j < diaIDs.length; j++) {
      diaIDs[j].nbEquivalent = diaIDs[j].atoms.length;
    }

    diaIDs.sort(function (a, b) {
      if (a.atomLabel == b.atomLabel) {
        return b.nbEquivalent - a.nbEquivalent;
      }

      return a.atomLabel < b.atomLabel ? 1 : -1;
    });
    entry.molecule = molecule;
    entry.diaIDs = diaIDs;
    entry.diaID = molecule.getIDCode();
  } else {
    molecule = entry.molecule;
    diaIDs = entry.diaIDs;
  } //H1 prediction


  var h1pred = predictor.predict(molecule, {
    group: true,
    ignoreLabile: false
  });
  if (!h1pred || h1pred.length === 0) return null; //console.log(h1pred);

  var optionsError = {
    iteration: options.iteration || 1,
    learningRatio: options.learningRatio || 1
  };

  for (var j = 0; j < h1pred.length; j++) {
    h1pred[j].error = getError(h1pred[j], optionsError);
  }

  h1pred.sort(function (a, b) {
    if (a.atomLabel == b.atomLabel) {
      return b.integral - a.integral;
    }

    return a.atomLabel < b.atomLabel ? 1 : -1;
  });

  try {
    spectra.h1PeakList.sort(function (a, b) {
      return b.integral - a.integral;
    });
    const spinSystem = new SpinSystem(h1pred, spectra.h1PeakList); //console.log(spinSystem);

    const autoAssigner = new AutoAssigner(spinSystem, options);
    return autoAssigner.getAssignments();
  } catch (e) {
    console.log("Could not assign this molecule.");
    return null;
  }
}

function getError(prediction, param) {
  //console.log(prediction)
  //Never use predictions with less than 3 votes
  if (prediction.std == 0 || prediction.ncs < 3) {
    return 20;
  } else {
    //factor is between 1 and +inf
    //console.log(prediction.ncs+" "+(param.iteration+1)+" "+param.learningRatio);
    var factor = 3 * prediction.std / Math.pow(prediction.ncs, (param.iteration + 1) * param.learningRatio); //(param.iteration+1)*param.learningRatio*h1pred[indexSignal].ncs;

    return 3 * prediction.std + factor;
  }

  return 20;
}

module.exports = autoAssign;

/***/ }),
/* 14 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";
/**
 * Created by acastillo on 9/2/16.
 */


const DEBUG = false;

class SpinSystem {
  constructor(diaIDsArray, signalsArray, opt) {
    var options = Object.assign({}, opt);
    this.diaIDsArray = diaIDsArray;
    this.signalsArray = signalsArray;
    this.cosy = options.cosySignals || null;
    this.connCosy = options.cosyPaths || null;
    this.connHmbc = options.hmbcPaths || null;
    this.hmbc = options.hmbcSignals || null;
    this.init();
  }

  init() {
    const nDiaIds = this.diaIDsArray.length;
    const nSignals = this.signalsArray.length;
    const diaIDByAtomLabel = {};
    const indexByAtomLabel = {};
    var shiftsH = [];
    var shiftsC = [];
    var windowH = [];
    var windowC = [];
    var signals1D = this.signalsArray;
    var nH = 0,
        nC = 0,
        i = 0;

    try {
      this.chemicalShiftsT = new Array(nDiaIds);
      this.chemicalShiftsTError = new Array(nDiaIds);
      this.diaList = new Array(nDiaIds);
      var dia = null;

      for (i = 0; i < nDiaIds; i++) {
        dia = this.diaIDsArray[i];

        if (diaIDByAtomLabel[dia.atomLabel]) {
          diaIDByAtomLabel[dia.atomLabel].push(dia.diaIDs[0]);
          indexByAtomLabel[dia.atomLabel].push(i);
        } else {
          diaIDByAtomLabel[dia.atomLabel] = [dia.diaIDs[0]];
          indexByAtomLabel[dia.atomLabel] = [i];
        }

        this.diaList[i] = dia.integral;
        this.chemicalShiftsT[i] = dia.delta;
        this.chemicalShiftsTError[i] = dia.error || 0;
      } // We can't have more signals than different protons in the molecule if the integral
      // matches the nH


      this.signals = new Array(nSignals);
      this.chemicalShiftsE = new Array(nSignals);
      this.signalsWidth = new Array(nSignals);

      for (i = 0; i < nSignals; i++) {
        var from = signals1D[i].from;
        var to = signals1D[i].to;
        this.chemicalShiftsE[i] = (from + to) / 2.0;
        this.signals[i] = Math.round(signals1D[i].integral);
        shiftsH.push((from + to) / 2.0);
        windowH.push(Math.abs(from - to));
        this.signalsWidth[i] = Math.abs(from - to);
      } //System.out.println(diaIDsH.size());
      //var row,col;

      /*if(cosy.length()>0&&connCosy.length()>0){
          cosyT = new byte[nH][nH];
          //To parse the theoretical COSY
          for(int i=connCosy.length()-1;i>=0;i--){
              JSONObject pair = (JSONObject) connCosy.get(i);
              row=diaIDsH.indexOf(pair.getString("diaID1"));
              col=diaIDsH.indexOf(pair.getString("diaID2"));
              if(row>=0&&col>=0){
                  cosyT[row][col]=(byte)4;
                  if(pair.getInt("distance")==4)
                      cosyT[row][col]=(byte)2;
              }
          }
          cosyE = new byte[nSignals][nSignals];
          double resolutionX = Double.MAX_VALUE/6;//((JSONObject)hmbc.get(0)).getDouble("resolutionX");
          if(((JSONObject)cosy.get(0)).has("resolutionX"))
          resolutionX = ((JSONObject)cosy.get(0)).getDouble("resolutionX");
           for(int i=cosy.length()-1;i>=0;i--){
              double x = 0,y=0;
              //JSONObject crossPeak = null;
              if(cosy.get(i) instanceof JSONObject){
                  x=((JSONObject)cosy.get(i)).getDouble("shiftX");
                  y=((JSONObject)cosy.get(i)).getDouble("shiftY");
              }
              else{
                  if(cosy.get(i) instanceof NMRSignal2D){
                      x =((NMRSignal2D)cosy.get(i)).toJSON().getDouble("shiftX");
                      y =((NMRSignal2D)cosy.get(i)).toJSON().getDouble("shiftY");
                  }
              }
              row = getIndex(x,shiftsH,windowH, resolutionX);
              col = getIndex(y,shiftsH,windowH, resolutionX);
              if(row>=0&&col>=0)
                  cosyE[row][col]=1;
          }
          //To complete the missing diagonal signals
          for(int i=0;i<cosyE.length;i++)
          cosyE[i][i]=1;
      }
       if(hmbc.length()>0&&connHmbc.length()>0){
          hmbcT = new byte[nH][nC];
          //To parse the theoretical HMBC
          for(int i=connHmbc.length()-1;i>=0;i--){
              JSONObject pair = (JSONObject) connHmbc.get(i);
              row=diaIDsH.indexOf(pair.getString("diaID1"));
              col=diaIDsC.indexOf(pair.getString("diaID2"));
              if(row>=0&&col>=0)
                  hmbcT[row][col]=1;
          }
            for(int i=hmbc.length()-1;i>=0;i--){
              //JSONObject crossPeak = null;
              double y=0;//x=0;
              if(hmbc.get(i) instanceof JSONObject){
                  //x=((JSONObject)hmbc.get(i)).getDouble("shiftX");
                  y=((JSONObject)hmbc.get(i)).getDouble("shiftY");
              }
              else{
                  if(hmbc.get(i) instanceof NMRSignal2D){
                      //crossPeak = (JSONObject) ((NMRSignal2D)hmbc.get(i)).toJSON().get("peaks");
                      //x =((NMRSignal2D)hmbc.get(i)).toJSON().getDouble("shiftX");
                      y =((NMRSignal2D)hmbc.get(i)).toJSON().getDouble("shiftY");
                  }
              }
              if(shiftsC.size()==0)
                  shiftsC.add(y);
              else{
                  int index = shiftsC.binarySearch(y);
                  if(index<0)
                      shiftsC.beforeInsert(-(index+1), y);
              }
          }
          if(DEBUG) System.out.println("shifts C : "+shiftsC);
          if(DEBUG)  System.out.println("diaIDs C : "+diaIDsC);
          shiftsC=fitCount(diaIDsC,shiftsC);
          if(DEBUG)  System.out.println("shifts C : "+shiftsC);
          hmbcE = new byte[nSignals][shiftsC.size()];
          double resolutionX = Double.MAX_VALUE/6;//((JSONObject)hmbc.get(0)).getDouble("resolutionX");
          if(((JSONObject)hmbc.get(0)).has("resolutionX"))
          resolutionX = ((JSONObject)hmbc.get(0)).getDouble("resolutionX");
          double resolutionY = Double.MAX_VALUE/6;//((JSONObject)hmbc.get(0)).getDouble("resolutionY");
          if(((JSONObject)hmbc.get(0)).has("resolutionY"))
          resolutionY = ((JSONObject)hmbc.get(0)).getDouble("resolutionY");
          for(int i=hmbc.length()-1;i>=0;i--){
              //JSONObject crossPeak = null;
              double x=0,y=0;
              if(hmbc.get(i) instanceof JSONObject){
                  x=((JSONObject)hmbc.get(i)).getDouble("shiftX");
                  y=((JSONObject)hmbc.get(i)).getDouble("shiftY");
              }
              else{
                  if(hmbc.get(i) instanceof NMRSignal2D){
                      //crossPeak = (JSONObject) ((NMRSignal2D)hmbc.get(i)).toJSON().get("peaks");
                      x =((NMRSignal2D)hmbc.get(i)).toJSON().getDouble("shiftX");
                      y =((NMRSignal2D)hmbc.get(i)).toJSON().getDouble("shiftY");
                  }
              }
              //System.out.println(resolutionX+" "+resolutionY+" "+windowH);
              row = getIndex(x,shiftsH, windowH, resolutionX);
              col = getIndex(y,shiftsC, null, resolutionY);
              //System.out.println(row+" , "+col);
              if(row>=0&&col>=0)
                  hmbcE[row][col]=1;
          }
      }*/

    } catch (e) {
      // TODO Auto-generated catch block
      console.log("Exception in SpinSystem " + e);
    }
  }

  _getIndex(value, shifts, windows, resolution) {
    var minDiff = Number.MAX_VALUE;
    var index = shifts.length - 1;

    for (var i = shifts.length - 1; i >= 0; i--) {
      if (Math.abs(value - shifts[i]) < minDiff) {
        minDiff = Math.abs(value - shifts[i]);
        index = i;
      }
    }

    if (windows) {
      if (minDiff <= windows[index] / 2) return index;
    }

    if (minDiff <= Math.abs(resolution * 4)) return index;
    return -1;
  }

}

module.exports = SpinSystem;

/***/ }),
/* 15 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";

/**
 * Created by acastillo on 9/2/16.
 */

const TreeSet = __webpack_require__(16);

const defaultOptions = {
  minScore: 1,
  maxSolutions: 100,
  errorCS: -1,
  onlyCount: false,
  timeout: 20000,
  condensed: true
};
const DEBUG = false;

class Assignment {
  constructor(spinSystem, opt) {
    var options = Object.assign({}, defaultOptions, opt);
    this.spinSystem = spinSystem;
    this.minScore = options.minScore;
    this.maxSolutions = options.maxSolutions;
    this.errorCS = options.errorCS;
    this.onlyCount = options.onlyCount;
    this.timeout = options.timeout;
    this.MAXERRORSHMBC = 1;
    this.condensed = options.condensed;
    this.timeoutTerminated = 0;
    this.score = 0;
    this.nSolutions = 0;
    this.nSteps = 0;
    this.lowerBound = 0;
    this.solutions = null;

    this.comparator = function (a, b) {
      return b.score - a.score;
    };
  }

  getAssignments() {
    var date = new Date();
    this.timeStart = date.getTime();
    var i, j, k, nSignals, nDiaIDs;
    if (DEBUG) console.log(this.spinSystem);
    this.lowerBound = this.minScore;

    do {
      this.nSolutions = 0;
      this.nSteps = 0;
      this.solutions = new TreeSet(this.comparator);
      if (this.spinSystem.hmbcE != null) this.spinSystem.hmbcLines = {};
      if (this.spinSystem.cosyE != null) this.spinSystem.cosyLines = {};
      nSignals = this.spinSystem.signals.length;
      nDiaIDs = this.spinSystem.diaList.length;
      this.scores = new Array(nSignals);
      let partial = new Array(nSignals);

      for (i = 0; i < nSignals; i++) {
        this.scores[i] = 1;
        partial[i] = [];
      }

      var diaMask = new Array(nDiaIDs);

      for (i = diaMask.length - 1; i >= 0; i--) diaMask[i] = true;

      try {
        this.exploreTreeRec(this.spinSystem.signals, this.spinSystem.diaList, nSignals - 1, nDiaIDs - 1, diaMask, partial);
      } catch (e) {
        console.log("Exception in assignment: " + e);
      }

      this.lowerBound -= 0.1;
      if (DEBUG) console.log("Decreasing lowerBound: " + this.lowerBound);
    } while (this.solutions.isEmpty() && this.lowerBound >= 0.4); //Format the result


    this._formatAssignmentOutput();

    return this.solutions.elements;
  }

  _formatAssignmentOutput(format) {
    var nSignals = this.spinSystem.signalsArray.length;
    var i, j, k;
    var assignment = this.solutions.elements;
    var nSolutions = this.solutions.length;

    for (i = 0; i < nSolutions; i++) {
      var assignment = this.solutions.elements[i].assignment;
      this.solutions.elements[i].index = i + "";
      var assignmentNew = {};

      for (j = 0; j < nSignals; j++) {
        var diaIDs = assignment[j];
        var tmp = new Array(diaIDs.length);

        for (k = 0; k < diaIDs.length; k++) {
          tmp[k] = this.spinSystem.diaIDsArray[diaIDs[k]].diaIDs[0];
        }

        if (this.condensed) assignmentNew[this.spinSystem.signalsArray[j].signalID] = tmp;else {
          assignment[j] = {
            signalID: this.spinSystem.signalsArray[j].signalID,
            delta: Math.round(this.spinSystem.signalsArray[j].signal[0].delta * 100) / 100,
            diaID: tmp
          };
        }
      }

      if (this.condensed) this.solutions.elements[i].assignment = assignmentNew;
    }
  }

  exploreTreeRec(signals, diaList, indexSignal, indexDia, diaMask, partial) {
    //If this happens, we can assign this atom group to this signal
    while (indexDia >= 0 && signals[indexSignal] >= diaList[indexDia]) {
      //Force a return if the loop time is longer than the given timeout
      const d = new Date();

      if (d.getTime() - this.timeStart > this.timeout) {
        this.timeoutTerminated = true;
        return;
      } //We can speed up it by checking the chemical shift first


      if (diaMask[indexDia] && this._isWithinCSRange(indexSignal, indexDia)) {
        this.nSteps++;
        const sizePartial = partial[indexSignal].length; //Assign the atom group to the signal

        diaMask[indexDia] = false; //Mark this atom group as used

        partial[indexSignal][sizePartial] = indexDia; //Add the atom group index to the assignment list

        signals[indexSignal] -= diaList[indexDia]; //Subtract the group from signal integral
        //If this signal is completely assigned, we have to verify all the restrictions

        if (signals[indexSignal] == 0) {
          let keySum = this._accomplishCounts(indexSignal, partial);

          if (DEBUG) console.log("Accomplish count: " + keySum);

          if (keySum != 0) {
            //Verify the restrictions. A good solution should give a high score
            this.score = this._solutionScore(partial, indexSignal, keySum);
            if (DEBUG) console.log(this.score + " " + partial); //This is a solution

            if (this.score > 0) {
              if (indexSignal == 0) {
                //We found a new solution
                this.nSolutions++;
                var solution = {
                  assignment: this._cloneArray(partial),
                  score: this.score
                };

                if (this.solutions.length >= this.maxSolutions) {
                  if (this.score > this.solutions.last().score) {
                    this.solutions.pollLast();
                    this.solutions.add(solution);
                  }
                } else {
                  this.solutions.add(solution);
                }
              } else {
                //Each new signal that we assign will produce a new level on the tree.
                indexSignal--; //Lets go forward with the next signal

                indexDia = diaList.length;

                while (!diaMask[--indexDia]); //SolutionTree newLevel = new SolutionTree(solution);


                this.exploreTreeRec(signals, diaList, indexSignal, indexDia, diaMask, partial);
                indexSignal++;
              }
            }
          }
        } else {
          //It says that the signal should be assigned by combining 2 or more signals
          const previousIndexDia = indexDia;

          while (indexDia > 0 && !diaMask[--indexDia]);

          if (indexDia >= 0) this.exploreTreeRec(signals, diaList, indexSignal, indexDia, diaMask, partial);
          indexDia = previousIndexDia;
        } //Deallocate this atom group to try the next one.


        indexDia = partial[indexSignal].splice(sizePartial, 1); //Add the atom group index to the assignment list

        diaMask[indexDia] = true; //Mark this atom group as available

        signals[indexSignal] += diaList[indexDia]; //Subtract the group from signal integral

        this.scores[indexSignal] = 1;
      }

      indexDia--;
    }
  }

  _cloneArray(data) {
    return JSON.parse(JSON.stringify(data));
  }

  _isWithinCSRange(indexSignal, indexDia) {
    if (this.spinSystem.chemicalShiftsE != null && this.spinSystem.chemicalShiftsT != null) {
      if (this.errorCS == 0) return true;
      var cfAtoms = this.spinSystem.chemicalShiftsT[indexDia];
      if (cfAtoms == -9999999) return true;
      var cfSignal = this.spinSystem.chemicalShiftsE[indexSignal];
      var error = this.spinSystem.chemicalShiftsTError[indexDia];
      if (error < Math.abs(this.errorCS)) error = this.errorCS;
      var csError = Math.abs(this.spinSystem.signalsWidth[indexSignal] / 2.0 + Math.abs(error));
      if (Math.abs(cfSignal - cfAtoms) <= csError) return true;else return false;
    }

    return true;
  }

  _accomplishCounts(indexSignal, partial) {
    //Check the chemical shift
    var keySum = -1;
    var keySumCOSY = 1;
    var keySumHMBC = 1;

    if (this.spinSystem.cosyE != null) {
      keySumCOSY = this._accomplishCount(indexSignal, partial[indexSigna + l], this.spinSystem.cosyT, this.spinSystem.cosyE, this.spinSystem.cosyLines, true, MAXERRORSCOSY);
      keySum = keySumCOSY;
    }

    if (keySum != 0) {
      if (this.spinSystem.hmbcE != null) {
        keySumHMBC = this._accomplishCount(indexSignal, partial[indexSignal], this.spinSystem.hmbcT, this.spinSystem.hmbcE, this.spinSystem.hmbcLines, false, MAXERRORSHMBC);
        keySum = keySumHMBC;
      }
    }

    return keySum;
  }
  /**
   * This function calculates the expected connectivity pattern each time that a experimental
   * signal is completely assigned to a set of theoretical signals. Once the new patter is calculate,
   * it verifies if this patter is potentially a candidate or not. The function stores the
   * calculated pattern in a hash-map so that each of those expensive operations is performed just once.
   * If the pattern accomplish with the count of expected and observed signals, it will return the key
   * of the patter in the hash map, or 0 in other case.
   * @param index
   * @param signals
   * @param theoretical
   * @param experimental
   * @param hashMap
   * @param isSymmetryc
   * @return
   */


  _accomplishCount(index, signals, theoretical, experimental, hashMap, isSymmetryc, maxErrors) {
    return 1;
  }
  /**
   * This function calculates the score of a given partial assignment.
   * @param partial
   * @param current
   * @param keySingalAsg
   * @return The score of the given function.
   * @throws JSONException
   */


  _solutionScore(partial, current, keySingalAsg) {
    var score = 0;
    var expLH = 0;

    if (this.spinSystem.cosyLines != null) {
      expLH++;
      score = this._cosyScore(partial, current, keySingalAsg);
    }

    if (this.spinSystem.hmbcLines != null) {
      expLH++;
      score += this._hmbcScore(partial, current, keySingalAsg);
    }

    if (this.spinSystem.chemicalShiftsT != null && this.errorCS > 0) {
      expLH++;
      score += this._chemicalShiftScore(partial, current, keySingalAsg);
    }

    if (expLH == 0) {
      expLH = 3;
      score = 3;
    }

    this.scores[current] = score / expLH;
    var sumLh = 0;
    var count = 0;

    for (var i = this.scores.length - 1; i >= 0; i--) {
      if (this.scores[i] != -1) {
        sumLh += this.scores[i];
        count++;
      }
    }

    if (sumLh < this.scores.length * this.lowerBound) return -sumLh / count;
    return sumLh / count;
  }
  /**
   * This function calculates the assignment score for the chemical shift restrictions.
   * @param partial
   * @param current
   * @param keySingalAsg
   * @return
   */


  _chemicalShiftScore(partial, current, keySingalAsg) {
    if (this.errorCS <= 0) return 1;
    var csSignal = this.spinSystem.chemicalShiftsE[current];
    var widthSignal = this.spinSystem.signalsWidth[current] / 2.0;
    var score = 0;
    var csGroup = 0;
    var diff = 0;
    var nbGroups = 0;

    try {
      var assignedGroups = partial[current];

      for (var i = assignedGroups.length - 1; i >= 0; i--) {
        csGroup = this.spinSystem.chemicalShiftsT[assignedGroups[i]];

        if (csGroup != -9999999) {
          nbGroups++;
          diff = Math.abs(csSignal - csGroup);
          if (diff <= widthSignal) score += 1;else {
            diff = Math.abs(diff - widthSignal);
            score += -0.25 / this.errorCS * diff + 1;
          }
        }
      }

      if (nbGroups == 0) return 1.0;
      return score / nbGroups;
    } catch (e) {
      console.log("Exception in chemical shift score function " + e);
    }

    return 1;
  }
  /**
   * This function calculates the assignment score for the COSY restrictions.
   * @param partial
   * @param current
   * @param keySingalAsg
   * @return
   */


  _cosyScore(partial, current, keySingalAsg) {
    var goodness = 0;
    var cosyLine = this.spinSystem.cosyLines[keySingalAsg];
    var count1 = 0;
    var count0 = 0;
    var size = cosyLine.length - 1;

    for (var i = partial.length - 1; i >= 0; i--) {
      try {
        var signal2 = partial[i];

        if (i != current && signal2.length > 0) {
          var key = 0; //The unique key for the union of those signals

          try {
            for (var j = signal2.length - 1; j >= 0; j--) {
              key |= 1 << signal2.getInt(j);
            }
          } catch (ex) {
            console.log("Exception in cosy score function " + ex);
          }

          var cosyLine2 = this.spinSystem.cosyLines[key];
          var crossPeak = false;

          for (var j = size; j >= 0; j--) {
            if (cosyLine[j] == 6 && cosyLine2[j] != 0) crossPeak = true;
          }

          if (crossPeak) count1++;else count0++;
          if (this.spinSystem.cosyE[current][i] == 0 && crossPeak) goodness -= 0.5;
          if (this.spinSystem.cosyE[current][i] == 1 && !crossPeak) goodness -= 0.5;
          if (this.spinSystem.cosyE[current][i] == 1 && crossPeak) goodness += 1;
          if (this.spinSystem.cosyE[current][i] == 0 && !crossPeak) goodness += 0.5;
        }
      } catch (e1) {
        console.log("Exception in cosy score function " + e);
      }
    }

    return Math.exp(-Math.abs(count1 + count0 / 2.0 - goodness) / 2.0);
  }
  /**
   * This function calculates the assignment score for the HMBC restrictions.
   * @param partial
   * @param current
   * @param keySingalAsg
   * @return
   */


  _hmbcScore(partial, current, keySingalAsg) {
    var hmbcLine = this.spinSystem.hmbcLines[keySingalAsg];
    var sizeT = hmbcLine.length - 1;
    var sizeE = this.spinSystem.hmbcE[0].length - 1;
    var freedom = sizeT - sizeE + this.MAXERRORSHMBC;
    var crossPeaks = 0;

    for (var j = sizeT; j >= 0; j--) {
      if (hmbcLine[j] == 1) crossPeaks++;
    }

    for (var j = sizeE; j >= 0; j--) {
      if (this.spinSystem.hmbcE[current][j] == 1) crossPeaks--;
    }

    if (crossPeaks < freedom) crossPeaks = freedom;
    return Math.exp(-Math.abs(crossPeaks - freedom) / (sizeT + 1));
  }

}

module.exports = Assignment;

/***/ }),
/* 16 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";

/**
 * Created by acastillo on 9/3/16.
 */

class TreeSet {
  constructor(compatator) {
    this.length = 0;
    this.elements = [];
    if (compatator) this.compatator = compatator;else this.compatator = function (a, b) {
      return a - b;
    };
  }

  size() {
    return this.elements.length;
  }

  last() {
    return this.elements[this.length - 1];
  }

  first() {
    return this.elements[0];
  }

  isEmpty() {
    return this.size() === 0;
  }

  pollLast() {
    if (this.length > 0) {
      this.length--;
      return this.elements.splice(this.length, 1);
    }

    return null;
  }

  pollFirst() {
    if (this.length > 0) {
      this.length--;
      return this.elements.splice(0, 1);
    }

    return null;
  }

  add(element) {
    let index = this.binarySearch(element);

    if (index < 0) {
      index = -index - 1;
    }

    this.elements.splice(index, 0, element);
    this.length++;
  }
  /**
   * Performs a binary search of value in array
   * @param {number[]} array - Array in which value will be searched. It must be sorted.
   * @param {number} value - Value to search in array
   * @return {number} If value is found, returns its index in array. Otherwise, returns a negative number indicating where the value should be inserted: -(index + 1)
   */


  binarySearch(value) {
    var low = 0;
    var high = this.elements.length - 1;

    while (low <= high) {
      var mid = low + high >>> 1;
      var midValue = this.elements[mid];
      var cmp = this.compatator(midValue, value);

      if (cmp < 0) {
        low = mid + 1;
      } else if (cmp > 0) {
        high = mid - 1;
      } else {
        return mid;
      }
    }

    return -(low + 1);
  }

}

module.exports = TreeSet;

/***/ }),
/* 17 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


__webpack_require__(6);

var OCL = __webpack_require__(2);

module.exports = exports = OCL;
exports.DB = __webpack_require__(19);
exports.RXN = __webpack_require__(25);
OCL.Molecule.prototype.getGroupedDiastereotopicAtomIDs = __webpack_require__(27);
OCL.Molecule.prototype.getExtendedDiastereotopicAtomIDs = __webpack_require__(28);
OCL.Molecule.prototype.toVisualizerMolfile = __webpack_require__(29);
OCL.Molecule.prototype.getGroupedHOSECodes = __webpack_require__(30);
OCL.Molecule.prototype.getNumberOfAtoms = __webpack_require__(31);
OCL.Molecule.prototype.toDiastereotopicSVG = __webpack_require__(32);
OCL.Molecule.prototype.getAtomsInfo = __webpack_require__(33);
OCL.Molecule.prototype.getAllPaths = __webpack_require__(34);
OCL.Molecule.prototype.getConnectivityMatrix = __webpack_require__(53);

/***/ }),
/* 18 */
/***/ (function(module, exports) {

// shim for using process in browser
var process = module.exports = {}; // cached from whatever global is present so that test runners that stub it
// don't break things.  But we need to wrap it in a try catch in case it is
// wrapped in strict mode code which doesn't define any globals.  It's inside a
// function because try/catches deoptimize in certain engines.

var cachedSetTimeout;
var cachedClearTimeout;

function defaultSetTimout() {
  throw new Error('setTimeout has not been defined');
}

function defaultClearTimeout() {
  throw new Error('clearTimeout has not been defined');
}

(function () {
  try {
    if (typeof setTimeout === 'function') {
      cachedSetTimeout = setTimeout;
    } else {
      cachedSetTimeout = defaultSetTimout;
    }
  } catch (e) {
    cachedSetTimeout = defaultSetTimout;
  }

  try {
    if (typeof clearTimeout === 'function') {
      cachedClearTimeout = clearTimeout;
    } else {
      cachedClearTimeout = defaultClearTimeout;
    }
  } catch (e) {
    cachedClearTimeout = defaultClearTimeout;
  }
})();

function runTimeout(fun) {
  if (cachedSetTimeout === setTimeout) {
    //normal enviroments in sane situations
    return setTimeout(fun, 0);
  } // if setTimeout wasn't available but was latter defined


  if ((cachedSetTimeout === defaultSetTimout || !cachedSetTimeout) && setTimeout) {
    cachedSetTimeout = setTimeout;
    return setTimeout(fun, 0);
  }

  try {
    // when when somebody has screwed with setTimeout but no I.E. maddness
    return cachedSetTimeout(fun, 0);
  } catch (e) {
    try {
      // When we are in I.E. but the script has been evaled so I.E. doesn't trust the global object when called normally
      return cachedSetTimeout.call(null, fun, 0);
    } catch (e) {
      // same as above but when it's a version of I.E. that must have the global object for 'this', hopfully our context correct otherwise it will throw a global error
      return cachedSetTimeout.call(this, fun, 0);
    }
  }
}

function runClearTimeout(marker) {
  if (cachedClearTimeout === clearTimeout) {
    //normal enviroments in sane situations
    return clearTimeout(marker);
  } // if clearTimeout wasn't available but was latter defined


  if ((cachedClearTimeout === defaultClearTimeout || !cachedClearTimeout) && clearTimeout) {
    cachedClearTimeout = clearTimeout;
    return clearTimeout(marker);
  }

  try {
    // when when somebody has screwed with setTimeout but no I.E. maddness
    return cachedClearTimeout(marker);
  } catch (e) {
    try {
      // When we are in I.E. but the script has been evaled so I.E. doesn't  trust the global object when called normally
      return cachedClearTimeout.call(null, marker);
    } catch (e) {
      // same as above but when it's a version of I.E. that must have the global object for 'this', hopfully our context correct otherwise it will throw a global error.
      // Some versions of I.E. have different rules for clearTimeout vs setTimeout
      return cachedClearTimeout.call(this, marker);
    }
  }
}

var queue = [];
var draining = false;
var currentQueue;
var queueIndex = -1;

function cleanUpNextTick() {
  if (!draining || !currentQueue) {
    return;
  }

  draining = false;

  if (currentQueue.length) {
    queue = currentQueue.concat(queue);
  } else {
    queueIndex = -1;
  }

  if (queue.length) {
    drainQueue();
  }
}

function drainQueue() {
  if (draining) {
    return;
  }

  var timeout = runTimeout(cleanUpNextTick);
  draining = true;
  var len = queue.length;

  while (len) {
    currentQueue = queue;
    queue = [];

    while (++queueIndex < len) {
      if (currentQueue) {
        currentQueue[queueIndex].run();
      }
    }

    queueIndex = -1;
    len = queue.length;
  }

  currentQueue = null;
  draining = false;
  runClearTimeout(timeout);
}

process.nextTick = function (fun) {
  var args = new Array(arguments.length - 1);

  if (arguments.length > 1) {
    for (var i = 1; i < arguments.length; i++) {
      args[i - 1] = arguments[i];
    }
  }

  queue.push(new Item(fun, args));

  if (queue.length === 1 && !draining) {
    runTimeout(drainQueue);
  }
}; // v8 likes predictible objects


function Item(fun, array) {
  this.fun = fun;
  this.array = array;
}

Item.prototype.run = function () {
  this.fun.apply(null, this.array);
};

process.title = 'browser';
process.browser = true;
process.env = {};
process.argv = [];
process.version = ''; // empty string to avoid regexp issues

process.versions = {};

function noop() {}

process.on = noop;
process.addListener = noop;
process.once = noop;
process.off = noop;
process.removeListener = noop;
process.removeAllListeners = noop;
process.emit = noop;
process.prependListener = noop;
process.prependOnceListener = noop;

process.listeners = function (name) {
  return [];
};

process.binding = function (name) {
  throw new Error('process.binding is not supported');
};

process.cwd = function () {
  return '/';
};

process.chdir = function (dir) {
  throw new Error('process.chdir is not supported');
};

process.umask = function () {
  return 0;
};

/***/ }),
/* 19 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";
/* WEBPACK VAR INJECTION */(function(setImmediate) {

var OCL = __webpack_require__(2);

var Molecule = OCL.Molecule;

var parseSDF = __webpack_require__(21);

var Papa = __webpack_require__(22);

var extend = __webpack_require__(23);

var moleculeCreator = __webpack_require__(24);

var defaultDBOptions = {
  length: 0,
  computeProperties: false
};

function DB(options) {
  options = extend({}, defaultDBOptions, options);
  this.data = new Array(options.length);
  this.molecules = new Array(options.length);
  this.statistics = null;
  this.length = 0;
  this.computeProperties = !!options.computeProperties;
  this.searcher = null;
}

var defaultSDFOptions = {
  onStep: function onStep(current, total) {}
};

DB.parseSDF = function (sdf, options) {
  if (typeof sdf !== 'string') {
    throw new TypeError('sdf must be a string');
  }

  options = extend({}, defaultSDFOptions, options);
  return new Promise(function (resolve, reject) {
    var parsed = parseSDF(sdf);
    var molecules = parsed.molecules;
    var db = new DB(options);
    db.statistics = parsed.statistics;
    var i = 0,
        l = molecules.length;
    parseNext();

    function parseNext() {
      if (i === l) {
        return resolve(db);
      }

      try {
        db.push(Molecule.fromMolfile(molecules[i].molfile.value), molecules[i]);
      } catch (e) {
        return reject(e);
      }

      options.onStep(++i, l);
      setImmediate(parseNext);
    }
  });
};

var defaultCSVOptions = {
  header: true,
  dynamicTyping: true,
  skipEmptyLines: true,
  onStep: function onStep(current, total) {}
};

DB.parseCSV = function (csv, options) {
  if (typeof csv !== 'string') {
    throw new TypeError('csv must be a string');
  }

  options = extend({}, defaultCSVOptions, options);
  return new Promise(function (resolve, reject) {
    var parsed = Papa.parse(csv, options);
    var fields = parsed.meta.fields;
    var stats = new Array(fields.length);
    var firstElement = parsed.data[0];
    var datatype, datafield;

    for (var i = 0; i < fields.length; i++) {
      stats[i] = {
        label: fields[i],
        isNumeric: typeof firstElement[fields[i]] === 'number'
      };
      var lowerField = fields[i].toLowerCase();

      if (moleculeCreator.has(lowerField)) {
        datatype = moleculeCreator.get(lowerField);
        datafield = fields[i];
      }
    }

    if (!datatype) {
      throw new Error('this document does not contain any molecule field');
    }

    var db = new DB(options);
    db.statistics = stats;
    var i = 0,
        l = parsed.data.length;
    parseNext();

    function parseNext() {
      if (i === l) {
        return resolve(db);
      }

      try {
        db.push(datatype(parsed.data[i][datafield]), parsed.data[i]);
      } catch (e) {
        return reject(e);
      }

      options.onStep(++i, l);
      setImmediate(parseNext);
    }
  });
};

DB.prototype.push = function (molecule, data) {
  if (data === undefined) data = {};
  this.molecules[this.length] = molecule;
  var molecularFormula = molecule.getMolecularFormula();

  if (!molecule.index) {
    molecule.index = molecule.getIndex();
    molecule.idcode = molecule.getIDCode();
    molecule.mw = molecularFormula.relativeWeight;
  }

  this.data[this.length++] = data;

  if (this.computeProperties) {
    var properties = molecule.getProperties();
    data.properties = {
      absoluteWeight: molecularFormula.absoluteWeight,
      relativeWeight: molecule.mw,
      formula: molecularFormula.formula,
      acceptorCount: properties.acceptorCount,
      donorCount: properties.donorCount,
      logP: properties.logP,
      logS: properties.logS,
      polarSurfaceArea: properties.polarSurfaceArea,
      rotatableBondCount: properties.rotatableBondCount,
      stereoCenterCount: properties.stereoCenterCount
    };
  }
};

var defaultSearchOptions = {
  format: 'oclid',
  mode: 'substructure',
  limit: 0
};

DB.prototype.search = function (query, options) {
  options = extend({}, defaultSearchOptions, options);

  if (typeof query === 'string') {
    query = moleculeCreator.get(options.format.toLowerCase())(query);
  } else if (!(query instanceof Molecule)) {
    throw new TypeError('toSearch must be a Molecule or string');
  }

  var result;

  switch (options.mode.toLowerCase()) {
    case 'exact':
      result = this.exactSearch(query, options.limit);
      break;

    case 'substructure':
      result = this.subStructureSearch(query, options.limit);
      break;

    case 'similarity':
      result = this.similaritySearch(query, options.limit);
      break;

    default:
      throw new Error('unknown search mode: ' + options.mode);
  }

  return result;
};

DB.prototype.exactSearch = function (query, limit) {
  var queryIdcode = query.getIDCode();
  var result = new DB();
  limit = limit || Number.MAX_SAFE_INTEGER;

  for (var i = 0; i < this.length; i++) {
    if (this.molecules[i].idcode === queryIdcode) {
      result.push(this.molecules[i], this.data[i]);
      if (result.length >= limit) break;
    }
  }

  return result;
};

DB.prototype.subStructureSearch = function (query, limit) {
  var needReset = false;

  if (!query.isFragment()) {
    needReset = true;
    query.setFragment(true);
  }

  var queryIndex = query.getIndex();
  var queryMW = query.getMolecularFormula().relativeWeight;
  var searcher = this.getSearcher();
  searcher.setFragment(query, queryIndex);
  var searchResult = [];

  for (var i = 0; i < this.length; i++) {
    searcher.setMolecule(this.molecules[i], this.molecules[i].index);

    if (searcher.isFragmentInMolecule()) {
      searchResult.push([this.molecules[i], i]);
    }
  }

  searchResult.sort(function (a, b) {
    return Math.abs(queryMW - a[0].mw) - Math.abs(queryMW - b[0].mw);
  });
  var length = limit || searchResult.length;
  var result = new DB({
    length: length
  });

  for (var i = 0; i < length; i++) {
    result.push(this.molecules[searchResult[i][1]], this.data[searchResult[i][1]]);
  }

  if (needReset) {
    query.setFragment(false);
  }

  return result;
};

DB.prototype.similaritySearch = function (query, limit) {
  var queryIndex = query.getIndex();
  var queryMW = query.getMolecularFormula().relativeWeight;
  var queryIDCode = query.getIDCode();
  var searchResult = new Array(this.length);
  var similarity;

  for (var i = 0; i < this.length; i++) {
    if (this.molecules[i].idcode === queryIDCode) {
      similarity = 1e10;
    } else {
      similarity = OCL.SSSearcherWithIndex.getSimilarityTanimoto(queryIndex, this.molecules[i].index) * 100000 - Math.abs(queryMW - this.molecules[i].mw) / 1000;
    }

    searchResult[i] = [similarity, i];
  }

  searchResult.sort(function (a, b) {
    return b[0] - a[0];
  });
  var length = limit || searchResult.length;
  var result = new DB({
    length: length
  });

  for (var i = 0; i < length; i++) {
    result.push(this.molecules[searchResult[i][1]], this.data[searchResult[i][1]]);
  }

  return result;
};

DB.prototype.getSearcher = function () {
  return this.searcher || (this.searcher = new OCL.SSSearcherWithIndex());
};

module.exports = DB;
/* WEBPACK VAR INJECTION */}.call(this, __webpack_require__(20).setImmediate))

/***/ }),
/* 20 */
/***/ (function(module, exports, __webpack_require__) {

/* WEBPACK VAR INJECTION */(function(global) {var scope = typeof global !== "undefined" && global || typeof self !== "undefined" && self || window;
var apply = Function.prototype.apply; // DOM APIs, for completeness

exports.setTimeout = function () {
  return new Timeout(apply.call(setTimeout, scope, arguments), clearTimeout);
};

exports.setInterval = function () {
  return new Timeout(apply.call(setInterval, scope, arguments), clearInterval);
};

exports.clearTimeout = exports.clearInterval = function (timeout) {
  if (timeout) {
    timeout.close();
  }
};

function Timeout(id, clearFn) {
  this._id = id;
  this._clearFn = clearFn;
}

Timeout.prototype.unref = Timeout.prototype.ref = function () {};

Timeout.prototype.close = function () {
  this._clearFn.call(scope, this._id);
}; // Does not start the time, just sets up the members needed.


exports.enroll = function (item, msecs) {
  clearTimeout(item._idleTimeoutId);
  item._idleTimeout = msecs;
};

exports.unenroll = function (item) {
  clearTimeout(item._idleTimeoutId);
  item._idleTimeout = -1;
};

exports._unrefActive = exports.active = function (item) {
  clearTimeout(item._idleTimeoutId);
  var msecs = item._idleTimeout;

  if (msecs >= 0) {
    item._idleTimeoutId = setTimeout(function onTimeout() {
      if (item._onTimeout) item._onTimeout();
    }, msecs);
  }
}; // setimmediate attaches itself to the global object


__webpack_require__(6); // On some exotic environments, it's not clear which object `setimmediate` was
// able to install onto.  Search each possibility in the same order as the
// `setimmediate` library.


exports.setImmediate = typeof self !== "undefined" && self.setImmediate || typeof global !== "undefined" && global.setImmediate || this && this.setImmediate;
exports.clearImmediate = typeof self !== "undefined" && self.clearImmediate || typeof global !== "undefined" && global.clearImmediate || this && this.clearImmediate;
/* WEBPACK VAR INJECTION */}.call(this, __webpack_require__(4)))

/***/ }),
/* 21 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";
 // options: an object

function parse(sdf, options) {
  // we will find the delimiter in order to be much faster and not use regular expression
  var header = sdf.substr(0, 1000);
  var crlf = '\n';

  if (header.indexOf('\r\n') > -1) {
    crlf = '\r\n';
  } else if (header.indexOf('\r') > -1) {
    crlf = '\r';
  }

  var sdfParts = sdf.split(crlf + '$$$$' + crlf);
  var molecules = [];
  var labels = {};
  var start = Date.now();
  var i = 0,
      ii = sdfParts.length,
      sdfPart,
      parts,
      molecule,
      j,
      jj,
      lines,
      from,
      to,
      label,
      k,
      kk;

  for (; i < ii; i++) {
    sdfPart = sdfParts[i];
    parts = sdfPart.split(crlf + '>');

    if (parts.length > 0 && parts[0].length > 5) {
      molecule = {};
      molecules.push(molecule);
      molecule.molfile = {
        type: 'mol2d',
        value: parts[0] + crlf
      };
      jj = parts.length;

      for (j = 1; j < jj; j++) {
        lines = parts[j].split(crlf);
        from = lines[0].indexOf('<');
        to = lines[0].indexOf('>');
        label = lines[0].substring(from + 1, to);

        if (labels[label]) {
          labels[label].counter++;
        } else {
          labels[label] = {
            counter: 1,
            isNumeric: true
          };
        }

        kk = lines.length - 1;

        for (k = 1; k < kk; k++) {
          if (molecule[label]) {
            molecule[label] += crlf + lines[k];
          } else {
            molecule[label] = lines[k];
          }
        }

        if (labels[label].isNumeric) {
          if (!isFinite(molecule[label])) {
            labels[label].isNumeric = false;
          }
        }
      }
    }
  } // all numeric fields should be converted to numbers


  var numericFields = [];

  for (var label in labels) {
    var currentLabel = labels[label];

    if (currentLabel.isNumeric) {
      currentLabel.minValue = Number.MAX_VALUE;
      currentLabel.maxValue = Number.MIN_VALUE;

      for (var j = 0; j < molecules.length; j++) {
        if (molecules[j][label]) {
          var value = parseFloat(molecules[j][label]);
          molecules[j][label] = value;
          if (value > currentLabel.maxValue) currentLabel.maxValue = value;
          if (value < currentLabel.minValue) currentLabel.minValue = value;
        }
      }
    }
  } // we check that a label is in all the records


  for (var key in labels) {
    if (labels[key].counter == molecules.length) {
      labels[key].always = true;
    } else {
      labels[key].always = false;
    }
  }

  var statistics = [];

  for (var key in labels) {
    var statistic = labels[key];
    statistic.label = key;
    statistics.push(statistic);
  }

  return {
    time: Date.now() - start,
    molecules: molecules,
    labels: Object.keys(labels),
    statistics: statistics
  };
}

module.exports = parse;

/***/ }),
/* 22 */
/***/ (function(module, exports, __webpack_require__) {

var __WEBPACK_AMD_DEFINE_RESULT__;/*!
	Papa Parse
	v4.1.2
	https://github.com/mholt/PapaParse
*/
(function (global) {
  "use strict";

  var IS_WORKER = !global.document && !!global.postMessage,
      IS_PAPA_WORKER = IS_WORKER && /(\?|&)papaworker(=|&|$)/.test(global.location.search),
      LOADED_SYNC = false,
      AUTO_SCRIPT_PATH;
  var workers = {},
      workerIdCounter = 0;
  var Papa = {};
  Papa.parse = CsvToJson;
  Papa.unparse = JsonToCsv;
  Papa.RECORD_SEP = String.fromCharCode(30);
  Papa.UNIT_SEP = String.fromCharCode(31);
  Papa.BYTE_ORDER_MARK = "\ufeff";
  Papa.BAD_DELIMITERS = ["\r", "\n", "\"", Papa.BYTE_ORDER_MARK];
  Papa.WORKERS_SUPPORTED = !IS_WORKER && !!global.Worker;
  Papa.SCRIPT_PATH = null; // Must be set by your code if you use workers and this lib is loaded asynchronously
  // Configurable chunk sizes for local and remote files, respectively

  Papa.LocalChunkSize = 1024 * 1024 * 10; // 10 MB

  Papa.RemoteChunkSize = 1024 * 1024 * 5; // 5 MB

  Papa.DefaultDelimiter = ","; // Used if not specified and detection fails
  // Exposed for testing and development only

  Papa.Parser = Parser;
  Papa.ParserHandle = ParserHandle;
  Papa.NetworkStreamer = NetworkStreamer;
  Papa.FileStreamer = FileStreamer;
  Papa.StringStreamer = StringStreamer;

  if ( true && module.exports) {
    // Export to Node...
    module.exports = Papa;
  } else if (isFunction(global.define) && global.define.amd) {
    // Wireup with RequireJS
    !(__WEBPACK_AMD_DEFINE_RESULT__ = (function () {
      return Papa;
    }).call(exports, __webpack_require__, exports, module),
				__WEBPACK_AMD_DEFINE_RESULT__ !== undefined && (module.exports = __WEBPACK_AMD_DEFINE_RESULT__));
  } else {
    // ...or as browser global
    global.Papa = Papa;
  }

  if (global.jQuery) {
    var $ = global.jQuery;

    $.fn.parse = function (options) {
      var config = options.config || {};
      var queue = [];
      this.each(function (idx) {
        var supported = $(this).prop('tagName').toUpperCase() == "INPUT" && $(this).attr('type').toLowerCase() == "file" && global.FileReader;
        if (!supported || !this.files || this.files.length == 0) return true; // continue to next input element

        for (var i = 0; i < this.files.length; i++) {
          queue.push({
            file: this.files[i],
            inputElem: this,
            instanceConfig: $.extend({}, config)
          });
        }
      });
      parseNextFile(); // begin parsing

      return this; // maintains chainability

      function parseNextFile() {
        if (queue.length == 0) {
          if (isFunction(options.complete)) options.complete();
          return;
        }

        var f = queue[0];

        if (isFunction(options.before)) {
          var returned = options.before(f.file, f.inputElem);

          if (typeof returned === 'object') {
            if (returned.action == "abort") {
              error("AbortError", f.file, f.inputElem, returned.reason);
              return; // Aborts all queued files immediately
            } else if (returned.action == "skip") {
              fileComplete(); // parse the next file in the queue, if any

              return;
            } else if (typeof returned.config === 'object') f.instanceConfig = $.extend(f.instanceConfig, returned.config);
          } else if (returned == "skip") {
            fileComplete(); // parse the next file in the queue, if any

            return;
          }
        } // Wrap up the user's complete callback, if any, so that ours also gets executed


        var userCompleteFunc = f.instanceConfig.complete;

        f.instanceConfig.complete = function (results) {
          if (isFunction(userCompleteFunc)) userCompleteFunc(results, f.file, f.inputElem);
          fileComplete();
        };

        Papa.parse(f.file, f.instanceConfig);
      }

      function error(name, file, elem, reason) {
        if (isFunction(options.error)) options.error({
          name: name
        }, file, elem, reason);
      }

      function fileComplete() {
        queue.splice(0, 1);
        parseNextFile();
      }
    };
  }

  if (IS_PAPA_WORKER) {
    global.onmessage = workerThreadReceivedMessage;
  } else if (Papa.WORKERS_SUPPORTED) {
    AUTO_SCRIPT_PATH = getScriptPath(); // Check if the script was loaded synchronously

    if (!document.body) {
      // Body doesn't exist yet, must be synchronous
      LOADED_SYNC = true;
    } else {
      document.addEventListener('DOMContentLoaded', function () {
        LOADED_SYNC = true;
      }, true);
    }
  }

  function CsvToJson(_input, _config) {
    _config = _config || {};

    if (_config.worker && Papa.WORKERS_SUPPORTED) {
      var w = newWorker();
      w.userStep = _config.step;
      w.userChunk = _config.chunk;
      w.userComplete = _config.complete;
      w.userError = _config.error;
      _config.step = isFunction(_config.step);
      _config.chunk = isFunction(_config.chunk);
      _config.complete = isFunction(_config.complete);
      _config.error = isFunction(_config.error);
      delete _config.worker; // prevent infinite loop

      w.postMessage({
        input: _input,
        config: _config,
        workerId: w.id
      });
      return;
    }

    var streamer = null;

    if (typeof _input === 'string') {
      if (_config.download) streamer = new NetworkStreamer(_config);else streamer = new StringStreamer(_config);
    } else if (global.File && _input instanceof File || _input instanceof Object) // ...Safari. (see issue #106)
      streamer = new FileStreamer(_config);

    return streamer.stream(_input);
  }

  function JsonToCsv(_input, _config) {
    var _output = "";
    var _fields = []; // Default configuration

    /** whether to surround every datum with quotes */

    var _quotes = false;
    /** delimiting character */

    var _delimiter = ",";
    /** newline character(s) */

    var _newline = "\r\n";
    unpackConfig();
    if (typeof _input === 'string') _input = JSON.parse(_input);

    if (_input instanceof Array) {
      if (!_input.length || _input[0] instanceof Array) return serialize(null, _input);else if (typeof _input[0] === 'object') return serialize(objectKeys(_input[0]), _input);
    } else if (typeof _input === 'object') {
      if (typeof _input.data === 'string') _input.data = JSON.parse(_input.data);

      if (_input.data instanceof Array) {
        if (!_input.fields) _input.fields = _input.data[0] instanceof Array ? _input.fields : objectKeys(_input.data[0]);
        if (!(_input.data[0] instanceof Array) && typeof _input.data[0] !== 'object') _input.data = [_input.data]; // handles input like [1,2,3] or ["asdf"]
      }

      return serialize(_input.fields || [], _input.data || []);
    } // Default (any valid paths should return before this)


    throw "exception: Unable to serialize unrecognized input";

    function unpackConfig() {
      if (typeof _config !== 'object') return;

      if (typeof _config.delimiter === 'string' && _config.delimiter.length == 1 && Papa.BAD_DELIMITERS.indexOf(_config.delimiter) == -1) {
        _delimiter = _config.delimiter;
      }

      if (typeof _config.quotes === 'boolean' || _config.quotes instanceof Array) _quotes = _config.quotes;
      if (typeof _config.newline === 'string') _newline = _config.newline;
    }
    /** Turns an object's keys into an array */


    function objectKeys(obj) {
      if (typeof obj !== 'object') return [];
      var keys = [];

      for (var key in obj) keys.push(key);

      return keys;
    }
    /** The double for loop that iterates the data and writes out a CSV string including header row */


    function serialize(fields, data) {
      var csv = "";
      if (typeof fields === 'string') fields = JSON.parse(fields);
      if (typeof data === 'string') data = JSON.parse(data);
      var hasHeader = fields instanceof Array && fields.length > 0;
      var dataKeyedByField = !(data[0] instanceof Array); // If there a header row, write it first

      if (hasHeader) {
        for (var i = 0; i < fields.length; i++) {
          if (i > 0) csv += _delimiter;
          csv += safe(fields[i], i);
        }

        if (data.length > 0) csv += _newline;
      } // Then write out the data


      for (var row = 0; row < data.length; row++) {
        var maxCol = hasHeader ? fields.length : data[row].length;

        for (var col = 0; col < maxCol; col++) {
          if (col > 0) csv += _delimiter;
          var colIdx = hasHeader && dataKeyedByField ? fields[col] : col;
          csv += safe(data[row][colIdx], col);
        }

        if (row < data.length - 1) csv += _newline;
      }

      return csv;
    }
    /** Encloses a value around quotes if needed (makes a value safe for CSV insertion) */


    function safe(str, col) {
      if (typeof str === "undefined" || str === null) return "";
      str = str.toString().replace(/"/g, '""');
      var needsQuotes = typeof _quotes === 'boolean' && _quotes || _quotes instanceof Array && _quotes[col] || hasAny(str, Papa.BAD_DELIMITERS) || str.indexOf(_delimiter) > -1 || str.charAt(0) == ' ' || str.charAt(str.length - 1) == ' ';
      return needsQuotes ? '"' + str + '"' : str;
    }

    function hasAny(str, substrings) {
      for (var i = 0; i < substrings.length; i++) if (str.indexOf(substrings[i]) > -1) return true;

      return false;
    }
  }
  /** ChunkStreamer is the base prototype for various streamer implementations. */


  function ChunkStreamer(config) {
    this._handle = null;
    this._paused = false;
    this._finished = false;
    this._input = null;
    this._baseIndex = 0;
    this._partialLine = "";
    this._rowCount = 0;
    this._start = 0;
    this._nextChunk = null;
    this.isFirstChunk = true;
    this._completeResults = {
      data: [],
      errors: [],
      meta: {}
    };
    replaceConfig.call(this, config);

    this.parseChunk = function (chunk) {
      // First chunk pre-processing
      if (this.isFirstChunk && isFunction(this._config.beforeFirstChunk)) {
        var modifiedChunk = this._config.beforeFirstChunk(chunk);

        if (modifiedChunk !== undefined) chunk = modifiedChunk;
      }

      this.isFirstChunk = false; // Rejoin the line we likely just split in two by chunking the file

      var aggregate = this._partialLine + chunk;
      this._partialLine = "";

      var results = this._handle.parse(aggregate, this._baseIndex, !this._finished);

      if (this._handle.paused() || this._handle.aborted()) return;
      var lastIndex = results.meta.cursor;

      if (!this._finished) {
        this._partialLine = aggregate.substring(lastIndex - this._baseIndex);
        this._baseIndex = lastIndex;
      }

      if (results && results.data) this._rowCount += results.data.length;
      var finishedIncludingPreview = this._finished || this._config.preview && this._rowCount >= this._config.preview;

      if (IS_PAPA_WORKER) {
        global.postMessage({
          results: results,
          workerId: Papa.WORKER_ID,
          finished: finishedIncludingPreview
        });
      } else if (isFunction(this._config.chunk)) {
        this._config.chunk(results, this._handle);

        if (this._paused) return;
        results = undefined;
        this._completeResults = undefined;
      }

      if (!this._config.step && !this._config.chunk) {
        this._completeResults.data = this._completeResults.data.concat(results.data);
        this._completeResults.errors = this._completeResults.errors.concat(results.errors);
        this._completeResults.meta = results.meta;
      }

      if (finishedIncludingPreview && isFunction(this._config.complete) && (!results || !results.meta.aborted)) this._config.complete(this._completeResults);
      if (!finishedIncludingPreview && (!results || !results.meta.paused)) this._nextChunk();
      return results;
    };

    this._sendError = function (error) {
      if (isFunction(this._config.error)) this._config.error(error);else if (IS_PAPA_WORKER && this._config.error) {
        global.postMessage({
          workerId: Papa.WORKER_ID,
          error: error,
          finished: false
        });
      }
    };

    function replaceConfig(config) {
      // Deep-copy the config so we can edit it
      var configCopy = copy(config);
      configCopy.chunkSize = parseInt(configCopy.chunkSize); // parseInt VERY important so we don't concatenate strings!

      if (!config.step && !config.chunk) configCopy.chunkSize = null; // disable Range header if not streaming; bad values break IIS - see issue #196

      this._handle = new ParserHandle(configCopy);
      this._handle.streamer = this;
      this._config = configCopy; // persist the copy to the caller
    }
  }

  function NetworkStreamer(config) {
    config = config || {};
    if (!config.chunkSize) config.chunkSize = Papa.RemoteChunkSize;
    ChunkStreamer.call(this, config);
    var xhr;

    if (IS_WORKER) {
      this._nextChunk = function () {
        this._readChunk();

        this._chunkLoaded();
      };
    } else {
      this._nextChunk = function () {
        this._readChunk();
      };
    }

    this.stream = function (url) {
      this._input = url;

      this._nextChunk(); // Starts streaming

    };

    this._readChunk = function () {
      if (this._finished) {
        this._chunkLoaded();

        return;
      }

      xhr = new XMLHttpRequest();

      if (!IS_WORKER) {
        xhr.onload = bindFunction(this._chunkLoaded, this);
        xhr.onerror = bindFunction(this._chunkError, this);
      }

      xhr.open("GET", this._input, !IS_WORKER);

      if (this._config.chunkSize) {
        var end = this._start + this._config.chunkSize - 1; // minus one because byte range is inclusive

        xhr.setRequestHeader("Range", "bytes=" + this._start + "-" + end);
        xhr.setRequestHeader("If-None-Match", "webkit-no-cache"); // https://bugs.webkit.org/show_bug.cgi?id=82672
      }

      try {
        xhr.send();
      } catch (err) {
        this._chunkError(err.message);
      }

      if (IS_WORKER && xhr.status == 0) this._chunkError();else this._start += this._config.chunkSize;
    };

    this._chunkLoaded = function () {
      if (xhr.readyState != 4) return;

      if (xhr.status < 200 || xhr.status >= 400) {
        this._chunkError();

        return;
      }

      this._finished = !this._config.chunkSize || this._start > getFileSize(xhr);
      this.parseChunk(xhr.responseText);
    };

    this._chunkError = function (errorMessage) {
      var errorText = xhr.statusText || errorMessage;

      this._sendError(errorText);
    };

    function getFileSize(xhr) {
      var contentRange = xhr.getResponseHeader("Content-Range");
      return parseInt(contentRange.substr(contentRange.lastIndexOf("/") + 1));
    }
  }

  NetworkStreamer.prototype = Object.create(ChunkStreamer.prototype);
  NetworkStreamer.prototype.constructor = NetworkStreamer;

  function FileStreamer(config) {
    config = config || {};
    if (!config.chunkSize) config.chunkSize = Papa.LocalChunkSize;
    ChunkStreamer.call(this, config);
    var reader, slice; // FileReader is better than FileReaderSync (even in worker) - see http://stackoverflow.com/q/24708649/1048862
    // But Firefox is a pill, too - see issue #76: https://github.com/mholt/PapaParse/issues/76

    var usingAsyncReader = typeof FileReader !== 'undefined'; // Safari doesn't consider it a function - see issue #105

    this.stream = function (file) {
      this._input = file;
      slice = file.slice || file.webkitSlice || file.mozSlice;

      if (usingAsyncReader) {
        reader = new FileReader(); // Preferred method of reading files, even in workers

        reader.onload = bindFunction(this._chunkLoaded, this);
        reader.onerror = bindFunction(this._chunkError, this);
      } else reader = new FileReaderSync(); // Hack for running in a web worker in Firefox


      this._nextChunk(); // Starts streaming

    };

    this._nextChunk = function () {
      if (!this._finished && (!this._config.preview || this._rowCount < this._config.preview)) this._readChunk();
    };

    this._readChunk = function () {
      var input = this._input;

      if (this._config.chunkSize) {
        var end = Math.min(this._start + this._config.chunkSize, this._input.size);
        input = slice.call(input, this._start, end);
      }

      var txt = reader.readAsText(input, this._config.encoding);
      if (!usingAsyncReader) this._chunkLoaded({
        target: {
          result: txt
        }
      }); // mimic the async signature
    };

    this._chunkLoaded = function (event) {
      // Very important to increment start each time before handling results
      this._start += this._config.chunkSize;
      this._finished = !this._config.chunkSize || this._start >= this._input.size;
      this.parseChunk(event.target.result);
    };

    this._chunkError = function () {
      this._sendError(reader.error);
    };
  }

  FileStreamer.prototype = Object.create(ChunkStreamer.prototype);
  FileStreamer.prototype.constructor = FileStreamer;

  function StringStreamer(config) {
    config = config || {};
    ChunkStreamer.call(this, config);
    var string;
    var remaining;

    this.stream = function (s) {
      string = s;
      remaining = s;
      return this._nextChunk();
    };

    this._nextChunk = function () {
      if (this._finished) return;
      var size = this._config.chunkSize;
      var chunk = size ? remaining.substr(0, size) : remaining;
      remaining = size ? remaining.substr(size) : '';
      this._finished = !remaining;
      return this.parseChunk(chunk);
    };
  }

  StringStreamer.prototype = Object.create(StringStreamer.prototype);
  StringStreamer.prototype.constructor = StringStreamer; // Use one ParserHandle per entire CSV file or string

  function ParserHandle(_config) {
    // One goal is to minimize the use of regular expressions...
    var FLOAT = /^\s*-?(\d*\.?\d+|\d+\.?\d*)(e[-+]?\d+)?\s*$/i;
    var self = this;
    var _stepCounter = 0; // Number of times step was called (number of rows parsed)

    var _input; // The input being parsed


    var _parser; // The core parser being used


    var _paused = false; // Whether we are paused or not

    var _aborted = false; // Whether the parser has aborted or not

    var _delimiterError; // Temporary state between delimiter detection and processing results


    var _fields = []; // Fields are from the header row of the input, if there is one

    var _results = {
      // The last results returned from the parser
      data: [],
      errors: [],
      meta: {}
    };

    if (isFunction(_config.step)) {
      var userStep = _config.step;

      _config.step = function (results) {
        _results = results;
        if (needsHeaderRow()) processResults();else // only call user's step function after header row
          {
            processResults(); // It's possbile that this line was empty and there's no row here after all

            if (_results.data.length == 0) return;
            _stepCounter += results.data.length;
            if (_config.preview && _stepCounter > _config.preview) _parser.abort();else userStep(_results, self);
          }
      };
    }
    /**
     * Parses input. Most users won't need, and shouldn't mess with, the baseIndex
     * and ignoreLastRow parameters. They are used by streamers (wrapper functions)
     * when an input comes in multiple chunks, like from a file.
     */


    this.parse = function (input, baseIndex, ignoreLastRow) {
      if (!_config.newline) _config.newline = guessLineEndings(input);
      _delimiterError = false;

      if (!_config.delimiter) {
        var delimGuess = guessDelimiter(input);
        if (delimGuess.successful) _config.delimiter = delimGuess.bestDelimiter;else {
          _delimiterError = true; // add error after parsing (otherwise it would be overwritten)

          _config.delimiter = Papa.DefaultDelimiter;
        }
        _results.meta.delimiter = _config.delimiter;
      }

      var parserConfig = copy(_config);
      if (_config.preview && _config.header) parserConfig.preview++; // to compensate for header row

      _input = input;
      _parser = new Parser(parserConfig);
      _results = _parser.parse(_input, baseIndex, ignoreLastRow);
      processResults();
      return _paused ? {
        meta: {
          paused: true
        }
      } : _results || {
        meta: {
          paused: false
        }
      };
    };

    this.paused = function () {
      return _paused;
    };

    this.pause = function () {
      _paused = true;

      _parser.abort();

      _input = _input.substr(_parser.getCharIndex());
    };

    this.resume = function () {
      _paused = false;
      self.streamer.parseChunk(_input);
    };

    this.aborted = function () {
      return _aborted;
    };

    this.abort = function () {
      _aborted = true;

      _parser.abort();

      _results.meta.aborted = true;
      if (isFunction(_config.complete)) _config.complete(_results);
      _input = "";
    };

    function processResults() {
      if (_results && _delimiterError) {
        addError("Delimiter", "UndetectableDelimiter", "Unable to auto-detect delimiting character; defaulted to '" + Papa.DefaultDelimiter + "'");
        _delimiterError = false;
      }

      if (_config.skipEmptyLines) {
        for (var i = 0; i < _results.data.length; i++) if (_results.data[i].length == 1 && _results.data[i][0] == "") _results.data.splice(i--, 1);
      }

      if (needsHeaderRow()) fillHeaderFields();
      return applyHeaderAndDynamicTyping();
    }

    function needsHeaderRow() {
      return _config.header && _fields.length == 0;
    }

    function fillHeaderFields() {
      if (!_results) return;

      for (var i = 0; needsHeaderRow() && i < _results.data.length; i++) for (var j = 0; j < _results.data[i].length; j++) _fields.push(_results.data[i][j]);

      _results.data.splice(0, 1);
    }

    function applyHeaderAndDynamicTyping() {
      if (!_results || !_config.header && !_config.dynamicTyping) return _results;

      for (var i = 0; i < _results.data.length; i++) {
        var row = {};

        for (var j = 0; j < _results.data[i].length; j++) {
          if (_config.dynamicTyping) {
            var value = _results.data[i][j];
            if (value == "true" || value == "TRUE") _results.data[i][j] = true;else if (value == "false" || value == "FALSE") _results.data[i][j] = false;else _results.data[i][j] = tryParseFloat(value);
          }

          if (_config.header) {
            if (j >= _fields.length) {
              if (!row["__parsed_extra"]) row["__parsed_extra"] = [];
              row["__parsed_extra"].push(_results.data[i][j]);
            } else row[_fields[j]] = _results.data[i][j];
          }
        }

        if (_config.header) {
          _results.data[i] = row;
          if (j > _fields.length) addError("FieldMismatch", "TooManyFields", "Too many fields: expected " + _fields.length + " fields but parsed " + j, i);else if (j < _fields.length) addError("FieldMismatch", "TooFewFields", "Too few fields: expected " + _fields.length + " fields but parsed " + j, i);
        }
      }

      if (_config.header && _results.meta) _results.meta.fields = _fields;
      return _results;
    }

    function guessDelimiter(input) {
      var delimChoices = [",", "\t", "|", ";", Papa.RECORD_SEP, Papa.UNIT_SEP];
      var bestDelim, bestDelta, fieldCountPrevRow;

      for (var i = 0; i < delimChoices.length; i++) {
        var delim = delimChoices[i];
        var delta = 0,
            avgFieldCount = 0;
        fieldCountPrevRow = undefined;
        var preview = new Parser({
          delimiter: delim,
          preview: 10
        }).parse(input);

        for (var j = 0; j < preview.data.length; j++) {
          var fieldCount = preview.data[j].length;
          avgFieldCount += fieldCount;

          if (typeof fieldCountPrevRow === 'undefined') {
            fieldCountPrevRow = fieldCount;
            continue;
          } else if (fieldCount > 1) {
            delta += Math.abs(fieldCount - fieldCountPrevRow);
            fieldCountPrevRow = fieldCount;
          }
        }

        if (preview.data.length > 0) avgFieldCount /= preview.data.length;

        if ((typeof bestDelta === 'undefined' || delta < bestDelta) && avgFieldCount > 1.99) {
          bestDelta = delta;
          bestDelim = delim;
        }
      }

      _config.delimiter = bestDelim;
      return {
        successful: !!bestDelim,
        bestDelimiter: bestDelim
      };
    }

    function guessLineEndings(input) {
      input = input.substr(0, 1024 * 1024); // max length 1 MB

      var r = input.split('\r');
      if (r.length == 1) return '\n';
      var numWithN = 0;

      for (var i = 0; i < r.length; i++) {
        if (r[i][0] == '\n') numWithN++;
      }

      return numWithN >= r.length / 2 ? '\r\n' : '\r';
    }

    function tryParseFloat(val) {
      var isNumber = FLOAT.test(val);
      return isNumber ? parseFloat(val) : val;
    }

    function addError(type, code, msg, row) {
      _results.errors.push({
        type: type,
        code: code,
        message: msg,
        row: row
      });
    }
  }
  /** The core parser implements speedy and correct CSV parsing */


  function Parser(config) {
    // Unpack the config object
    config = config || {};
    var delim = config.delimiter;
    var newline = config.newline;
    var comments = config.comments;
    var step = config.step;
    var preview = config.preview;
    var fastMode = config.fastMode; // Delimiter must be valid

    if (typeof delim !== 'string' || Papa.BAD_DELIMITERS.indexOf(delim) > -1) delim = ","; // Comment character must be valid

    if (comments === delim) throw "Comment character same as delimiter";else if (comments === true) comments = "#";else if (typeof comments !== 'string' || Papa.BAD_DELIMITERS.indexOf(comments) > -1) comments = false; // Newline must be valid: \r, \n, or \r\n

    if (newline != '\n' && newline != '\r' && newline != '\r\n') newline = '\n'; // We're gonna need these at the Parser scope

    var cursor = 0;
    var aborted = false;

    this.parse = function (input, baseIndex, ignoreLastRow) {
      // For some reason, in Chrome, this speeds things up (!?)
      if (typeof input !== 'string') throw "Input must be a string"; // We don't need to compute some of these every time parse() is called,
      // but having them in a more local scope seems to perform better

      var inputLen = input.length,
          delimLen = delim.length,
          newlineLen = newline.length,
          commentsLen = comments.length;
      var stepIsFunction = typeof step === 'function'; // Establish starting state

      cursor = 0;
      var data = [],
          errors = [],
          row = [],
          lastCursor = 0;
      if (!input) return returnable();

      if (fastMode || fastMode !== false && input.indexOf('"') === -1) {
        var rows = input.split(newline);

        for (var i = 0; i < rows.length; i++) {
          var row = rows[i];
          cursor += row.length;
          if (i !== rows.length - 1) cursor += newline.length;else if (ignoreLastRow) return returnable();
          if (comments && row.substr(0, commentsLen) == comments) continue;

          if (stepIsFunction) {
            data = [];
            pushRow(row.split(delim));
            doStep();
            if (aborted) return returnable();
          } else pushRow(row.split(delim));

          if (preview && i >= preview) {
            data = data.slice(0, preview);
            return returnable(true);
          }
        }

        return returnable();
      }

      var nextDelim = input.indexOf(delim, cursor);
      var nextNewline = input.indexOf(newline, cursor); // Parser loop

      for (;;) {
        // Field has opening quote
        if (input[cursor] == '"') {
          // Start our search for the closing quote where the cursor is
          var quoteSearch = cursor; // Skip the opening quote

          cursor++;

          for (;;) {
            // Find closing quote
            var quoteSearch = input.indexOf('"', quoteSearch + 1);

            if (quoteSearch === -1) {
              if (!ignoreLastRow) {
                // No closing quote... what a pity
                errors.push({
                  type: "Quotes",
                  code: "MissingQuotes",
                  message: "Quoted field unterminated",
                  row: data.length,
                  // row has yet to be inserted
                  index: cursor
                });
              }

              return finish();
            }

            if (quoteSearch === inputLen - 1) {
              // Closing quote at EOF
              var value = input.substring(cursor, quoteSearch).replace(/""/g, '"');
              return finish(value);
            } // If this quote is escaped, it's part of the data; skip it


            if (input[quoteSearch + 1] == '"') {
              quoteSearch++;
              continue;
            }

            if (input[quoteSearch + 1] == delim) {
              // Closing quote followed by delimiter
              row.push(input.substring(cursor, quoteSearch).replace(/""/g, '"'));
              cursor = quoteSearch + 1 + delimLen;
              nextDelim = input.indexOf(delim, cursor);
              nextNewline = input.indexOf(newline, cursor);
              break;
            }

            if (input.substr(quoteSearch + 1, newlineLen) === newline) {
              // Closing quote followed by newline
              row.push(input.substring(cursor, quoteSearch).replace(/""/g, '"'));
              saveRow(quoteSearch + 1 + newlineLen);
              nextDelim = input.indexOf(delim, cursor); // because we may have skipped the nextDelim in the quoted field

              if (stepIsFunction) {
                doStep();
                if (aborted) return returnable();
              }

              if (preview && data.length >= preview) return returnable(true);
              break;
            }
          }

          continue;
        } // Comment found at start of new line


        if (comments && row.length === 0 && input.substr(cursor, commentsLen) === comments) {
          if (nextNewline == -1) // Comment ends at EOF
            return returnable();
          cursor = nextNewline + newlineLen;
          nextNewline = input.indexOf(newline, cursor);
          nextDelim = input.indexOf(delim, cursor);
          continue;
        } // Next delimiter comes before next newline, so we've reached end of field


        if (nextDelim !== -1 && (nextDelim < nextNewline || nextNewline === -1)) {
          row.push(input.substring(cursor, nextDelim));
          cursor = nextDelim + delimLen;
          nextDelim = input.indexOf(delim, cursor);
          continue;
        } // End of row


        if (nextNewline !== -1) {
          row.push(input.substring(cursor, nextNewline));
          saveRow(nextNewline + newlineLen);

          if (stepIsFunction) {
            doStep();
            if (aborted) return returnable();
          }

          if (preview && data.length >= preview) return returnable(true);
          continue;
        }

        break;
      }

      return finish();

      function pushRow(row) {
        data.push(row);
        lastCursor = cursor;
      }
      /**
       * Appends the remaining input from cursor to the end into
       * row, saves the row, calls step, and returns the results.
       */


      function finish(value) {
        if (ignoreLastRow) return returnable();
        if (typeof value === 'undefined') value = input.substr(cursor);
        row.push(value);
        cursor = inputLen; // important in case parsing is paused

        pushRow(row);
        if (stepIsFunction) doStep();
        return returnable();
      }
      /**
       * Appends the current row to the results. It sets the cursor
       * to newCursor and finds the nextNewline. The caller should
       * take care to execute user's step function and check for
       * preview and end parsing if necessary.
       */


      function saveRow(newCursor) {
        cursor = newCursor;
        pushRow(row);
        row = [];
        nextNewline = input.indexOf(newline, cursor);
      }
      /** Returns an object with the results, errors, and meta. */


      function returnable(stopped) {
        return {
          data: data,
          errors: errors,
          meta: {
            delimiter: delim,
            linebreak: newline,
            aborted: aborted,
            truncated: !!stopped,
            cursor: lastCursor + (baseIndex || 0)
          }
        };
      }
      /** Executes the user's step function and resets data & errors. */


      function doStep() {
        step(returnable());
        data = [], errors = [];
      }
    };
    /** Sets the abort flag */


    this.abort = function () {
      aborted = true;
    };
    /** Gets the cursor position */


    this.getCharIndex = function () {
      return cursor;
    };
  } // If you need to load Papa Parse asynchronously and you also need worker threads, hard-code
  // the script path here. See: https://github.com/mholt/PapaParse/issues/87#issuecomment-57885358


  function getScriptPath() {
    var scripts = document.getElementsByTagName('script');
    return scripts.length ? scripts[scripts.length - 1].src : '';
  }

  function newWorker() {
    if (!Papa.WORKERS_SUPPORTED) return false;
    if (!LOADED_SYNC && Papa.SCRIPT_PATH === null) throw new Error('Script path cannot be determined automatically when Papa Parse is loaded asynchronously. ' + 'You need to set Papa.SCRIPT_PATH manually.');
    var workerUrl = Papa.SCRIPT_PATH || AUTO_SCRIPT_PATH; // Append "papaworker" to the search string to tell papaparse that this is our worker.

    workerUrl += (workerUrl.indexOf('?') !== -1 ? '&' : '?') + 'papaworker';
    var w = new global.Worker(workerUrl);
    w.onmessage = mainThreadReceivedMessage;
    w.id = workerIdCounter++;
    workers[w.id] = w;
    return w;
  }
  /** Callback when main thread receives a message */


  function mainThreadReceivedMessage(e) {
    var msg = e.data;
    var worker = workers[msg.workerId];
    var aborted = false;
    if (msg.error) worker.userError(msg.error, msg.file);else if (msg.results && msg.results.data) {
      var abort = function abort() {
        aborted = true;
        completeWorker(msg.workerId, {
          data: [],
          errors: [],
          meta: {
            aborted: true
          }
        });
      };

      var handle = {
        abort: abort,
        pause: notImplemented,
        resume: notImplemented
      };

      if (isFunction(worker.userStep)) {
        for (var i = 0; i < msg.results.data.length; i++) {
          worker.userStep({
            data: [msg.results.data[i]],
            errors: msg.results.errors,
            meta: msg.results.meta
          }, handle);
          if (aborted) break;
        }

        delete msg.results; // free memory ASAP
      } else if (isFunction(worker.userChunk)) {
        worker.userChunk(msg.results, handle, msg.file);
        delete msg.results;
      }
    }
    if (msg.finished && !aborted) completeWorker(msg.workerId, msg.results);
  }

  function completeWorker(workerId, results) {
    var worker = workers[workerId];
    if (isFunction(worker.userComplete)) worker.userComplete(results);
    worker.terminate();
    delete workers[workerId];
  }

  function notImplemented() {
    throw "Not implemented.";
  }
  /** Callback when worker thread receives a message */


  function workerThreadReceivedMessage(e) {
    var msg = e.data;
    if (typeof Papa.WORKER_ID === 'undefined' && msg) Papa.WORKER_ID = msg.workerId;

    if (typeof msg.input === 'string') {
      global.postMessage({
        workerId: Papa.WORKER_ID,
        results: Papa.parse(msg.input, msg.config),
        finished: true
      });
    } else if (global.File && msg.input instanceof File || msg.input instanceof Object) // thank you, Safari (see issue #106)
      {
        var results = Papa.parse(msg.input, msg.config);
        if (results) global.postMessage({
          workerId: Papa.WORKER_ID,
          results: results,
          finished: true
        });
      }
  }
  /** Makes a deep copy of an array or object (mostly) */


  function copy(obj) {
    if (typeof obj !== 'object') return obj;
    var cpy = obj instanceof Array ? [] : {};

    for (var key in obj) cpy[key] = copy(obj[key]);

    return cpy;
  }

  function bindFunction(f, self) {
    return function () {
      f.apply(self, arguments);
    };
  }

  function isFunction(func) {
    return typeof func === 'function';
  }
})(typeof window !== 'undefined' ? window : this);

/***/ }),
/* 23 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


var hasOwn = Object.prototype.hasOwnProperty;
var toStr = Object.prototype.toString;

var isArray = function isArray(arr) {
  if (typeof Array.isArray === 'function') {
    return Array.isArray(arr);
  }

  return toStr.call(arr) === '[object Array]';
};

var isPlainObject = function isPlainObject(obj) {
  if (!obj || toStr.call(obj) !== '[object Object]') {
    return false;
  }

  var hasOwnConstructor = hasOwn.call(obj, 'constructor');
  var hasIsPrototypeOf = obj.constructor && obj.constructor.prototype && hasOwn.call(obj.constructor.prototype, 'isPrototypeOf'); // Not own constructor property must be Object

  if (obj.constructor && !hasOwnConstructor && !hasIsPrototypeOf) {
    return false;
  } // Own properties are enumerated firstly, so to speed up,
  // if last one is own, then all properties are own.


  var key;

  for (key in obj) {
    /**/
  }

  return typeof key === 'undefined' || hasOwn.call(obj, key);
};

module.exports = function extend() {
  var options,
      name,
      src,
      copy,
      copyIsArray,
      clone,
      target = arguments[0],
      i = 1,
      length = arguments.length,
      deep = false; // Handle a deep copy situation

  if (typeof target === 'boolean') {
    deep = target;
    target = arguments[1] || {}; // skip the boolean and the target

    i = 2;
  } else if (typeof target !== 'object' && typeof target !== 'function' || target == null) {
    target = {};
  }

  for (; i < length; ++i) {
    options = arguments[i]; // Only deal with non-null/undefined values

    if (options != null) {
      // Extend the base object
      for (name in options) {
        src = target[name];
        copy = options[name]; // Prevent never-ending loop

        if (target !== copy) {
          // Recurse if we're merging plain objects or arrays
          if (deep && copy && (isPlainObject(copy) || (copyIsArray = isArray(copy)))) {
            if (copyIsArray) {
              copyIsArray = false;
              clone = src && isArray(src) ? src : [];
            } else {
              clone = src && isPlainObject(src) ? src : {};
            } // Never move original objects, clone them


            target[name] = extend(deep, clone, copy); // Don't bring in undefined values
          } else if (typeof copy !== 'undefined') {
            target[name] = copy;
          }
        }
      }
    }
  } // Return the modified object


  return target;
};

/***/ }),
/* 24 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


var Molecule = __webpack_require__(2).Molecule;

var fields = new Map();
module.exports = fields;
fields.set('oclid', Molecule.fromIDCode);
fields.set('idcode', Molecule.fromIDCode);
fields.set('smiles', Molecule.fromSmiles);
fields.set('molfile', Molecule.fromMolfile);

/***/ }),
/* 25 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


var OCL = __webpack_require__(2);

var parseRXN = __webpack_require__(26);

function RXN(rxn, options) {
  if (!rxn) {
    this.reagents = [];
    this.products = [];
  } else {
    var parsed = parseRXN(rxn);
    this.reagents = generateInfo(parsed.reagents);
    this.products = generateInfo(parsed.products);
  }
}

RXN.prototype.addReagent = function (molfile) {
  this.reagents.push(getMolfileInfo(molfile));
};

RXN.prototype.addProduct = function (molfile) {
  this.products.push(getMolfileInfo(molfile));
};

RXN.prototype.toRXN = function () {
  var result = [];
  result.push("$RXN");
  result.push("");
  result.push("");
  result.push("Openchemlib");
  result.push(format3(this.reagents.length) + format3(this.products.length));

  for (var i = 0; i < this.reagents.length; i++) {
    result.push("$MOL");
    result.push(getMolfile(this.reagents[i].molfile));
  }

  for (var i = 0; i < this.products.length; i++) {
    result.push("$MOL");
    result.push(getMolfile(this.products[i].molfile));
  }

  return result.join("\n");
};

function getMolfile(molfile) {
  var lines = ~molfile.indexOf("\r\n") ? molfile.split("\r\n") : molfile.split(/[\r\n]/);
  return lines.join("\n");
}

function format3(number) {
  var length = (number + "").length;
  return "   ".substring(0, 3 - length) + number;
}

function generateInfo(molecules) {
  for (var i = 0; i < molecules.length; i++) {
    molecules[i] = getMolfileInfo(molecules[i]);
  }

  return molecules;
}

function getMolfileInfo(molfile) {
  var ocl = OCL.Molecule.fromMolfile(molfile);
  return {
    molfile: molfile,
    smiles: ocl.toSmiles(),
    mf: ocl.getMolecularFormula().formula,
    mw: ocl.getMolecularFormula().relativeWeight,
    idCode: ocl.getIDCode
  };
}

module.exports = RXN;

/***/ }),
/* 26 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


function parse(rxn) {
  if (typeof rxn !== 'string') {
    throw new TypeError('Parameter "rxn" must be a string');
  } // we will find the delimiter in order to be much faster and not use regular expression


  var header = rxn.substr(0, 1000);
  var crlf = '\n';

  if (header.indexOf('\r\n') > -1) {
    crlf = '\r\n';
  } else if (header.indexOf('\r') > -1) {
    crlf = '\r';
  }

  var rxnParts = rxn.split(crlf + '$MOL' + crlf);
  var reagents = [];
  var products = [];
  var result = {};
  result.reagents = reagents;
  result.products = products; // the first part is expected to contain the number of reagents and products
  // First part should start with $RXN
  // and the fifth line should contain the number of reagents and products

  if (rxnParts.length === 0) throw new Error('file looks empty');
  var header = rxnParts[0];
  if (header.indexOf("$RXN") != 0) throw new Error('file does not start with $RXN');
  var lines = header.split(crlf);
  if (lines.length < 5) throw new Error('incorrect number of lines in header');
  var numberReagents = lines[4].substring(0, 3) >> 0;
  var numberProducts = lines[4].substring(3, 6) >> 0;
  if (numberReagents + numberProducts != rxnParts.length - 1) throw new Error('not the correct number of molecules');

  for (var i = 1; i < rxnParts.length; i++) {
    if (i <= numberReagents) {
      reagents.push(rxnParts[i]);
    } else {
      products.push(rxnParts[i]);
    }
  }

  return result;
}

module.exports = parse;

/***/ }),
/* 27 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


module.exports = function getGroupedDiastereotopicAtomIDs(options) {
  var options = options || {};
  var label = options.atomLabel;
  var diaIDs = this.getDiastereotopicAtomIDs(options);
  var diaIDsObject = {};

  for (var i = 0; i < diaIDs.length; i++) {
    if (!label || this.getAtomLabel(i) === label) {
      var diaID = diaIDs[i];

      if (!diaIDsObject[diaID]) {
        diaIDsObject[diaID] = {
          counter: 1,
          atoms: [i],
          oclID: diaID,
          atomLabel: this.getAtomLabel(i),
          _highlight: [diaID]
        };
      } else {
        diaIDsObject[diaID].counter++;
        diaIDsObject[diaID].atoms.push(i);
      }
    }
  }

  var diaIDsTable = [];

  for (var key of Object.keys(diaIDsObject)) {
    diaIDsTable.push(diaIDsObject[key]);
  }

  return diaIDsTable;
};

/***/ }),
/* 28 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


module.exports = function getExtendedDiastereotopicAtomIDs() {
  var molecule = this.getCompactCopy();
  molecule.addImplicitHydrogens();
  var diaIDs = molecule.getDiastereotopicAtomIDs();
  var newDiaIDs = [];

  for (var i = 0; i < diaIDs.length; i++) {
    var diaID = diaIDs[i];
    var newDiaID = {
      oclID: diaID,
      hydrogenOCLIDs: [],
      nbHydrogens: 0
    };

    for (var j = 0; j < molecule.getAllConnAtoms(i); j++) {
      var atom = molecule.getConnAtom(i, j);

      if (molecule.getAtomicNo(atom) === 1) {
        newDiaID.nbHydrogens++;

        if (newDiaID.hydrogenOCLIDs.indexOf(diaIDs[atom]) === -1) {
          newDiaID.hydrogenOCLIDs.push(diaIDs[atom]);
        }
      }
    }

    newDiaIDs.push(newDiaID);
  }

  return newDiaIDs;
};

/***/ }),
/* 29 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


module.exports = function toVisualizerMolfile() {
  var diaIDs = this.getGroupedDiastereotopicAtomIDs();
  var highlight = [];
  var atoms = {};
  diaIDs.forEach(function (diaID) {
    atoms[diaID.oclID] = diaID.atoms;
    highlight.push(diaID.oclID);
  });
  var molfile = {
    type: 'mol2d',
    value: this.toMolfile(),
    _highlight: highlight,
    _atoms: atoms
  };
  return molfile;
};

/***/ }),
/* 30 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


var Util = __webpack_require__(2).Util;

module.exports = function getGroupedHOSECodes(options) {
  var options = options || {};
  var diaIDs = this.getGroupedDiastereotopicAtomIDs(options);
  diaIDs.forEach(function (diaID) {
    var hoses = Util.getHoseCodesFromDiastereotopicID(diaID.oclID, options);
    diaID.hoses = [];
    var level = 1;

    for (var hose of hoses) {
      diaID.hoses.push({
        level: level++,
        oclID: hose
      });
    }
  });
  return diaIDs;
};

/***/ }),
/* 31 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


module.exports = function getNumberOfAtoms(options) {
  var options = options || {};
  var label = options.atomLabel;
  var mf = this.getMolecularFormula().formula;
  var parts = mf.split(/(?=[A-Z])/);

  for (var part of parts) {
    var atom = part.replace(/[0-9]/g, '');

    if (atom == label) {
      return part.replace(/[^0-9]/g, '') * 1 || 1;
    }
  }

  return 0;
};

/***/ }),
/* 32 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


module.exports = function toDiastereotopicSVG(options) {
  var options = options || {};
  var width = options.width || 300;
  var height = options.width || 200;
  var prefix = options.width || "ocl";
  var svg = options.svg;
  var diaIDs = this.getDiastereotopicAtomIDs();
  if (!svg) svg = this.toSVG(width, height, prefix);
  svg = svg.replace(/Atom:[0-9]+\"/g, function (value) {
    var atom = value.replace(/[^0-9]/g, '');
    return value + ' data-atomid="' + diaIDs[atom] + '"';
  });
  return svg;
};

/***/ }),
/* 33 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


module.exports = function getAtomsInfo() {
  var diaIDs = this.getDiastereotopicAtomIDs();
  var results = [];

  for (var i = 0; i < diaIDs.length; i++) {
    var result = {
      oclID: diaIDs[i],
      extra: {
        singleBonds: 0,
        doubleBonds: 0,
        tripleBonds: 0,
        aromaticBonds: 0,
        cnoHybridation: 0 // should be 1 (sp), 2 (sp2) or 3 (sp3)

      }
    };
    var extra = result.extra;
    results.push(result);
    result.abnormalValence = this.getAtomAbnormalValence(i); // -1 is normal otherwise specified

    result.charge = this.getAtomCharge(i);
    result.cipParity = this.getAtomCIPParity(i);
    result.color = this.getAtomColor(i);
    result.customLabel = this.getAtomCustomLabel(i); //        result.esrGroup=this.getAtomESRGroup(i);
    //        result.esrType=this.getAtomESRType(i);

    result.atomicNo = this.getAtomicNo(i);
    result.label = this.getAtomLabel(i); //        result.list=this.getAtomList(i);
    //        result.listString=this.getAtomListString(i);
    //        result.mapNo=this.getAtomMapNo(i);

    result.mass = this.getAtomMass(i); //        result.parity=this.getAtomParity(i);
    //        result.pi=this.getAtomPi(i);
    //        result.preferredStereoBond=this.getAtomPreferredStereoBond(i);
    //        result.queryFeatures=this.getAtomQueryFeatures(i);

    result.radical = this.getAtomRadical(i);
    result.ringBondCount = this.getAtomRingBondCount(i); //        result.ringCount=this.getAtomRingCount(i);

    result.ringSize = this.getAtomRingSize(i);
    result.x = this.getAtomX(i);
    result.y = this.getAtomY(i);
    result.z = this.getAtomZ(i);
    result.allHydrogens = this.getAllHydrogens(i);
    result.connAtoms = this.getConnAtoms(i);
    result.allConnAtoms = this.getAllConnAtoms(i);
    result.isAromatic = this.isAromaticAtom(i);
    result.isAllylic = this.isAllylicAtom(i);
    result.isStereoCenter = this.isAtomStereoCenter(i);
    result.isRing = this.isRingAtom(i);
    result.isSmallRing = this.isSmallRingAtom(i);
    result.isStabilized = this.isStabilizedAtom(i);
    result.extra.singleBonds = result.allHydrogens;

    for (var j = 0; j < this.getAllConnAtoms(i); j++) {
      var bond = this.getConnBond(i, j);
      var bondOrder = this.getBondOrder(bond);

      if (this.isAromaticBond(bond)) {
        extra.aromaticBonds++;
      } else if (bondOrder === 1) {
        extra.singleBonds++;
      } else if (bondOrder === 2) {
        extra.doubleBonds++;
      } else if (bondOrder === 3) {
        extra.tripleBonds++;
      }
    }

    result.extra.totalBonds = result.extra.singleBonds + result.extra.doubleBonds + result.extra.tripleBonds + result.extra.aromaticBonds;

    if (result.atomicNo === 6) {
      result.extra.cnoHybridation = result.extra.totalBonds - 1;
    } else if (result.atomicNo === 7) {
      result.extra.cnoHybridation = result.extra.totalBonds;
    } else if (result.atomicNo === 8) {
      result.extra.cnoHybridation = result.extra.totalBonds + 1;
    }
  }

  return results;
};

/***/ }),
/* 34 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


var floydWarshall = __webpack_require__(35);

var Matrix = __webpack_require__(7);

module.exports = function getAllPaths(options) {
  var options = options || [];
  var fromLabel = options.fromLabel || '';
  var toLabel = options.toLabel || '';
  var minLength = options.minLength === undefined ? 1 : options.minLength;
  var maxLength = options.maxLength === undefined ? 4 : options.maxLength; // we need to find all the atoms 'fromLabel' and 'toLabel'

  var results = {};
  var diaIDs = this.getDiastereotopicAtomIDs();
  var connectivityMatrix = this.getConnectivityMatrix(); // TODO have a package that allows to convert the connectivityMatrix to a distanceMatrix

  var pathLengthMatrix = floydWarshall(new Matrix(connectivityMatrix));

  for (var from = 0; from < this.getAllAtoms(); from++) {
    for (var to = 0; to < this.getAllAtoms(); to++) {
      if (!fromLabel || this.getAtomLabel(from) === fromLabel) {
        if (!toLabel || this.getAtomLabel(to) === toLabel) {
          var key = diaIDs[from] + "_" + diaIDs[to];
          var pathLength = pathLengthMatrix[from][to];

          if (pathLength >= minLength && pathLength <= maxLength) {
            if (!results[key]) {
              results[key] = {
                fromDiaID: diaIDs[from],
                toDiaID: diaIDs[to],
                fromAtoms: [],
                toAtoms: [],
                fromLabel: this.getAtomLabel(from),
                toLabel: this.getAtomLabel(to),
                pathLength: pathLength
              };
            }

            if (results[key].fromAtoms.indexOf(from) === -1) results[key].fromAtoms.push(from);
            if (results[key].toAtoms.indexOf(to) === -1) results[key].toAtoms.push(to);
          }
        }
      }
    }
  }

  var finalResults = [];

  for (var key in results) {
    finalResults.push(results[key]);
  }

  return finalResults;
};

/***/ }),
/* 35 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


const Matrix = __webpack_require__(7);
/**
 * Algorithm that finds the shortest distance from one node to the other
 * @param {Matrix} adjMatrix - A squared adjacency matrix
 * @return {Matrix} - Distance from a node to the other, -1 if the node is unreachable
 */


function floydWarshall(adjMatrix) {
  if (Matrix.isMatrix(adjMatrix) && adjMatrix.columns !== adjMatrix.rows) throw new TypeError('The adjacency matrix should be squared');
  const numVertices = adjMatrix.columns;
  let distMatrix = new Matrix(numVertices, numVertices);
  distMatrix.apply((row, column) => {
    // principal diagonal is 0
    if (row === column) distMatrix.set(row, column, 0);else {
      let val = adjMatrix.get(row, column); // edges values remain the same

      if (val) distMatrix.set(row, column, val); // 0 values become infinity
      else distMatrix.set(row, column, Number.POSITIVE_INFINITY);
    }
  });

  for (let k = 0; k < numVertices; ++k) for (let i = 0; i < numVertices; ++i) for (let j = 0; j < numVertices; ++j) {
    let dist = distMatrix.get(i, k) + distMatrix.get(k, j);
    if (distMatrix.get(i, j) > dist) distMatrix.set(i, j, dist);
  } // When there's no connection the value is -1


  distMatrix.apply((row, column) => {
    if (distMatrix.get(row, column) === Number.POSITIVE_INFINITY) distMatrix.set(row, column, -1);
  });
  return distMatrix;
}

module.exports = floydWarshall;

/***/ }),
/* 36 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


if (!Symbol.species) {
  Symbol.species = Symbol.for('@@species');
}

/***/ }),
/* 37 */
/***/ (function(module, exports, __webpack_require__) {

module.exports = exports = __webpack_require__(38);
exports.getEquallySpacedData = __webpack_require__(40).getEquallySpacedData;
exports.SNV = __webpack_require__(41).SNV;

/***/ }),
/* 38 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


const Stat = __webpack_require__(11).array;
/**
 * Function that returns an array of points given 1D array as follows:
 *
 * [x1, y1, .. , x2, y2, ..]
 *
 * And receive the number of dimensions of each point.
 * @param array
 * @param dimensions
 * @returns {Array} - Array of points.
 */


function coordArrayToPoints(array, dimensions) {
  if (array.length % dimensions !== 0) {
    throw new RangeError('Dimensions number must be accordance with the size of the array.');
  }

  var length = array.length / dimensions;
  var pointsArr = new Array(length);
  var k = 0;

  for (var i = 0; i < array.length; i += dimensions) {
    var point = new Array(dimensions);

    for (var j = 0; j < dimensions; ++j) {
      point[j] = array[i + j];
    }

    pointsArr[k] = point;
    k++;
  }

  return pointsArr;
}
/**
 * Function that given an array as follows:
 * [x1, y1, .. , x2, y2, ..]
 *
 * Returns an array as follows:
 * [[x1, x2, ..], [y1, y2, ..], [ .. ]]
 *
 * And receives the number of dimensions of each coordinate.
 * @param array
 * @param dimensions
 * @returns {Array} - Matrix of coordinates
 */


function coordArrayToCoordMatrix(array, dimensions) {
  if (array.length % dimensions !== 0) {
    throw new RangeError('Dimensions number must be accordance with the size of the array.');
  }

  var coordinatesArray = new Array(dimensions);
  var points = array.length / dimensions;

  for (var i = 0; i < coordinatesArray.length; i++) {
    coordinatesArray[i] = new Array(points);
  }

  for (i = 0; i < array.length; i += dimensions) {
    for (var j = 0; j < dimensions; ++j) {
      var currentPoint = Math.floor(i / dimensions);
      coordinatesArray[j][currentPoint] = array[i + j];
    }
  }

  return coordinatesArray;
}
/**
 * Function that receives a coordinate matrix as follows:
 * [[x1, x2, ..], [y1, y2, ..], [ .. ]]
 *
 * Returns an array of coordinates as follows:
 * [x1, y1, .. , x2, y2, ..]
 *
 * @param coordMatrix
 * @returns {Array}
 */


function coordMatrixToCoordArray(coordMatrix) {
  var coodinatesArray = new Array(coordMatrix.length * coordMatrix[0].length);
  var k = 0;

  for (var i = 0; i < coordMatrix[0].length; ++i) {
    for (var j = 0; j < coordMatrix.length; ++j) {
      coodinatesArray[k] = coordMatrix[j][i];
      ++k;
    }
  }

  return coodinatesArray;
}
/**
 * Tranpose a matrix, this method is for coordMatrixToPoints and
 * pointsToCoordMatrix, that because only transposing the matrix
 * you can change your representation.
 *
 * @param matrix
 * @returns {Array}
 */


function transpose(matrix) {
  var resultMatrix = new Array(matrix[0].length);

  for (var i = 0; i < resultMatrix.length; ++i) {
    resultMatrix[i] = new Array(matrix.length);
  }

  for (i = 0; i < matrix.length; ++i) {
    for (var j = 0; j < matrix[0].length; ++j) {
      resultMatrix[j][i] = matrix[i][j];
    }
  }

  return resultMatrix;
}
/**
 * Function that transform an array of points into a coordinates array
 * as follows:
 * [x1, y1, .. , x2, y2, ..]
 *
 * @param points
 * @returns {Array}
 */


function pointsToCoordArray(points) {
  var coodinatesArray = new Array(points.length * points[0].length);
  var k = 0;

  for (var i = 0; i < points.length; ++i) {
    for (var j = 0; j < points[0].length; ++j) {
      coodinatesArray[k] = points[i][j];
      ++k;
    }
  }

  return coodinatesArray;
}
/**
 * Apply the dot product between the smaller vector and a subsets of the
 * largest one.
 *
 * @param firstVector
 * @param secondVector
 * @returns {Array} each dot product of size of the difference between the
 *                  larger and the smallest one.
 */


function applyDotProduct(firstVector, secondVector) {
  var largestVector, smallestVector;

  if (firstVector.length <= secondVector.length) {
    smallestVector = firstVector;
    largestVector = secondVector;
  } else {
    smallestVector = secondVector;
    largestVector = firstVector;
  }

  var difference = largestVector.length - smallestVector.length + 1;
  var dotProductApplied = new Array(difference);

  for (var i = 0; i < difference; ++i) {
    var sum = 0;

    for (var j = 0; j < smallestVector.length; ++j) {
      sum += smallestVector[j] * largestVector[i + j];
    }

    dotProductApplied[i] = sum;
  }

  return dotProductApplied;
}
/**
 * To scale the input array between the specified min and max values. The operation is performed inplace
 * if the options.inplace is specified. If only one of the min or max parameters is specified, then the scaling
 * will multiply the input array by min/min(input) or max/max(input)
 * @param input
 * @param options
 * @returns {*}
 */


function scale(input, options) {
  var y;

  if (options.inPlace) {
    y = input;
  } else {
    y = new Array(input.length);
  }

  const max = options.max;
  const min = options.min;

  if (typeof max === "number") {
    if (typeof min === "number") {
      var minMax = Stat.minMax(input);
      var factor = (max - min) / (minMax.max - minMax.min);

      for (var i = 0; i < y.length; i++) {
        y[i] = (input[i] - minMax.min) * factor + min;
      }
    } else {
      var currentMin = Stat.max(input);
      var factor = max / currentMin;

      for (var i = 0; i < y.length; i++) {
        y[i] = input[i] * factor;
      }
    }
  } else {
    if (typeof min === "number") {
      var currentMin = Stat.min(input);
      var factor = min / currentMin;

      for (var i = 0; i < y.length; i++) {
        y[i] = input[i] * factor;
      }
    }
  }

  return y;
}

module.exports = {
  coordArrayToPoints: coordArrayToPoints,
  coordArrayToCoordMatrix: coordArrayToCoordMatrix,
  coordMatrixToCoordArray: coordMatrixToCoordArray,
  coordMatrixToPoints: transpose,
  pointsToCoordArray: pointsToCoordArray,
  pointsToCoordMatrix: transpose,
  applyDotProduct: applyDotProduct,
  scale: scale
};

/***/ }),
/* 39 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


var arrayStat = __webpack_require__(12);

function compareNumbers(a, b) {
  return a - b;
}

exports.max = function max(matrix) {
  var max = -Infinity;

  for (var i = 0; i < matrix.length; i++) {
    for (var j = 0; j < matrix[i].length; j++) {
      if (matrix[i][j] > max) max = matrix[i][j];
    }
  }

  return max;
};

exports.min = function min(matrix) {
  var min = Infinity;

  for (var i = 0; i < matrix.length; i++) {
    for (var j = 0; j < matrix[i].length; j++) {
      if (matrix[i][j] < min) min = matrix[i][j];
    }
  }

  return min;
};

exports.minMax = function minMax(matrix) {
  var min = Infinity;
  var max = -Infinity;

  for (var i = 0; i < matrix.length; i++) {
    for (var j = 0; j < matrix[i].length; j++) {
      if (matrix[i][j] < min) min = matrix[i][j];
      if (matrix[i][j] > max) max = matrix[i][j];
    }
  }

  return {
    min: min,
    max: max
  };
};

exports.entropy = function entropy(matrix, eps) {
  if (typeof eps === 'undefined') {
    eps = 0;
  }

  var sum = 0,
      l1 = matrix.length,
      l2 = matrix[0].length;

  for (var i = 0; i < l1; i++) {
    for (var j = 0; j < l2; j++) {
      sum += matrix[i][j] * Math.log(matrix[i][j] + eps);
    }
  }

  return -sum;
};

exports.mean = function mean(matrix, dimension) {
  if (typeof dimension === 'undefined') {
    dimension = 0;
  }

  var rows = matrix.length,
      cols = matrix[0].length,
      theMean,
      N,
      i,
      j;

  if (dimension === -1) {
    theMean = [0];
    N = rows * cols;

    for (i = 0; i < rows; i++) {
      for (j = 0; j < cols; j++) {
        theMean[0] += matrix[i][j];
      }
    }

    theMean[0] /= N;
  } else if (dimension === 0) {
    theMean = new Array(cols);
    N = rows;

    for (j = 0; j < cols; j++) {
      theMean[j] = 0;

      for (i = 0; i < rows; i++) {
        theMean[j] += matrix[i][j];
      }

      theMean[j] /= N;
    }
  } else if (dimension === 1) {
    theMean = new Array(rows);
    N = cols;

    for (j = 0; j < rows; j++) {
      theMean[j] = 0;

      for (i = 0; i < cols; i++) {
        theMean[j] += matrix[j][i];
      }

      theMean[j] /= N;
    }
  } else {
    throw new Error('Invalid dimension');
  }

  return theMean;
};

exports.sum = function sum(matrix, dimension) {
  if (typeof dimension === 'undefined') {
    dimension = 0;
  }

  var rows = matrix.length,
      cols = matrix[0].length,
      theSum,
      i,
      j;

  if (dimension === -1) {
    theSum = [0];

    for (i = 0; i < rows; i++) {
      for (j = 0; j < cols; j++) {
        theSum[0] += matrix[i][j];
      }
    }
  } else if (dimension === 0) {
    theSum = new Array(cols);

    for (j = 0; j < cols; j++) {
      theSum[j] = 0;

      for (i = 0; i < rows; i++) {
        theSum[j] += matrix[i][j];
      }
    }
  } else if (dimension === 1) {
    theSum = new Array(rows);

    for (j = 0; j < rows; j++) {
      theSum[j] = 0;

      for (i = 0; i < cols; i++) {
        theSum[j] += matrix[j][i];
      }
    }
  } else {
    throw new Error('Invalid dimension');
  }

  return theSum;
};

exports.product = function product(matrix, dimension) {
  if (typeof dimension === 'undefined') {
    dimension = 0;
  }

  var rows = matrix.length,
      cols = matrix[0].length,
      theProduct,
      i,
      j;

  if (dimension === -1) {
    theProduct = [1];

    for (i = 0; i < rows; i++) {
      for (j = 0; j < cols; j++) {
        theProduct[0] *= matrix[i][j];
      }
    }
  } else if (dimension === 0) {
    theProduct = new Array(cols);

    for (j = 0; j < cols; j++) {
      theProduct[j] = 1;

      for (i = 0; i < rows; i++) {
        theProduct[j] *= matrix[i][j];
      }
    }
  } else if (dimension === 1) {
    theProduct = new Array(rows);

    for (j = 0; j < rows; j++) {
      theProduct[j] = 1;

      for (i = 0; i < cols; i++) {
        theProduct[j] *= matrix[j][i];
      }
    }
  } else {
    throw new Error('Invalid dimension');
  }

  return theProduct;
};

exports.standardDeviation = function standardDeviation(matrix, means, unbiased) {
  var vari = exports.variance(matrix, means, unbiased),
      l = vari.length;

  for (var i = 0; i < l; i++) {
    vari[i] = Math.sqrt(vari[i]);
  }

  return vari;
};

exports.variance = function variance(matrix, means, unbiased) {
  if (typeof unbiased === 'undefined') {
    unbiased = true;
  }

  means = means || exports.mean(matrix);
  var rows = matrix.length;
  if (rows === 0) return [];
  var cols = matrix[0].length;
  var vari = new Array(cols);

  for (var j = 0; j < cols; j++) {
    var sum1 = 0,
        sum2 = 0,
        x = 0;

    for (var i = 0; i < rows; i++) {
      x = matrix[i][j] - means[j];
      sum1 += x;
      sum2 += x * x;
    }

    if (unbiased) {
      vari[j] = (sum2 - sum1 * sum1 / rows) / (rows - 1);
    } else {
      vari[j] = (sum2 - sum1 * sum1 / rows) / rows;
    }
  }

  return vari;
};

exports.median = function median(matrix) {
  var rows = matrix.length,
      cols = matrix[0].length;
  var medians = new Array(cols);

  for (var i = 0; i < cols; i++) {
    var data = new Array(rows);

    for (var j = 0; j < rows; j++) {
      data[j] = matrix[j][i];
    }

    data.sort(compareNumbers);
    var N = data.length;

    if (N % 2 === 0) {
      medians[i] = (data[N / 2] + data[N / 2 - 1]) * 0.5;
    } else {
      medians[i] = data[Math.floor(N / 2)];
    }
  }

  return medians;
};

exports.mode = function mode(matrix) {
  var rows = matrix.length,
      cols = matrix[0].length,
      modes = new Array(cols),
      i,
      j;

  for (i = 0; i < cols; i++) {
    var itemCount = new Array(rows);

    for (var k = 0; k < rows; k++) {
      itemCount[k] = 0;
    }

    var itemArray = new Array(rows);
    var count = 0;

    for (j = 0; j < rows; j++) {
      var index = itemArray.indexOf(matrix[j][i]);

      if (index >= 0) {
        itemCount[index]++;
      } else {
        itemArray[count] = matrix[j][i];
        itemCount[count] = 1;
        count++;
      }
    }

    var maxValue = 0,
        maxIndex = 0;

    for (j = 0; j < count; j++) {
      if (itemCount[j] > maxValue) {
        maxValue = itemCount[j];
        maxIndex = j;
      }
    }

    modes[i] = itemArray[maxIndex];
  }

  return modes;
};

exports.skewness = function skewness(matrix, unbiased) {
  if (typeof unbiased === 'undefined') unbiased = true;
  var means = exports.mean(matrix);
  var n = matrix.length,
      l = means.length;
  var skew = new Array(l);

  for (var j = 0; j < l; j++) {
    var s2 = 0,
        s3 = 0;

    for (var i = 0; i < n; i++) {
      var dev = matrix[i][j] - means[j];
      s2 += dev * dev;
      s3 += dev * dev * dev;
    }

    var m2 = s2 / n;
    var m3 = s3 / n;
    var g = m3 / Math.pow(m2, 3 / 2);

    if (unbiased) {
      var a = Math.sqrt(n * (n - 1));
      var b = n - 2;
      skew[j] = a / b * g;
    } else {
      skew[j] = g;
    }
  }

  return skew;
};

exports.kurtosis = function kurtosis(matrix, unbiased) {
  if (typeof unbiased === 'undefined') unbiased = true;
  var means = exports.mean(matrix);
  var n = matrix.length,
      m = matrix[0].length;
  var kurt = new Array(m);

  for (var j = 0; j < m; j++) {
    var s2 = 0,
        s4 = 0;

    for (var i = 0; i < n; i++) {
      var dev = matrix[i][j] - means[j];
      s2 += dev * dev;
      s4 += dev * dev * dev * dev;
    }

    var m2 = s2 / n;
    var m4 = s4 / n;

    if (unbiased) {
      var v = s2 / (n - 1);
      var a = n * (n + 1) / ((n - 1) * (n - 2) * (n - 3));
      var b = s4 / (v * v);
      var c = (n - 1) * (n - 1) / ((n - 2) * (n - 3));
      kurt[j] = a * b - 3 * c;
    } else {
      kurt[j] = m4 / (m2 * m2) - 3;
    }
  }

  return kurt;
};

exports.standardError = function standardError(matrix) {
  var samples = matrix.length;
  var standardDeviations = exports.standardDeviation(matrix);
  var l = standardDeviations.length;
  var standardErrors = new Array(l);
  var sqrtN = Math.sqrt(samples);

  for (var i = 0; i < l; i++) {
    standardErrors[i] = standardDeviations[i] / sqrtN;
  }

  return standardErrors;
};

exports.covariance = function covariance(matrix, dimension) {
  return exports.scatter(matrix, undefined, dimension);
};

exports.scatter = function scatter(matrix, divisor, dimension) {
  if (typeof dimension === 'undefined') {
    dimension = 0;
  }

  if (typeof divisor === 'undefined') {
    if (dimension === 0) {
      divisor = matrix.length - 1;
    } else if (dimension === 1) {
      divisor = matrix[0].length - 1;
    }
  }

  var means = exports.mean(matrix, dimension);
  var rows = matrix.length;

  if (rows === 0) {
    return [[]];
  }

  var cols = matrix[0].length,
      cov,
      i,
      j,
      s,
      k;

  if (dimension === 0) {
    cov = new Array(cols);

    for (i = 0; i < cols; i++) {
      cov[i] = new Array(cols);
    }

    for (i = 0; i < cols; i++) {
      for (j = i; j < cols; j++) {
        s = 0;

        for (k = 0; k < rows; k++) {
          s += (matrix[k][j] - means[j]) * (matrix[k][i] - means[i]);
        }

        s /= divisor;
        cov[i][j] = s;
        cov[j][i] = s;
      }
    }
  } else if (dimension === 1) {
    cov = new Array(rows);

    for (i = 0; i < rows; i++) {
      cov[i] = new Array(rows);
    }

    for (i = 0; i < rows; i++) {
      for (j = i; j < rows; j++) {
        s = 0;

        for (k = 0; k < cols; k++) {
          s += (matrix[j][k] - means[j]) * (matrix[i][k] - means[i]);
        }

        s /= divisor;
        cov[i][j] = s;
        cov[j][i] = s;
      }
    }
  } else {
    throw new Error('Invalid dimension');
  }

  return cov;
};

exports.correlation = function correlation(matrix) {
  var means = exports.mean(matrix),
      standardDeviations = exports.standardDeviation(matrix, true, means),
      scores = exports.zScores(matrix, means, standardDeviations),
      rows = matrix.length,
      cols = matrix[0].length,
      i,
      j;
  var cor = new Array(cols);

  for (i = 0; i < cols; i++) {
    cor[i] = new Array(cols);
  }

  for (i = 0; i < cols; i++) {
    for (j = i; j < cols; j++) {
      var c = 0;

      for (var k = 0, l = scores.length; k < l; k++) {
        c += scores[k][j] * scores[k][i];
      }

      c /= rows - 1;
      cor[i][j] = c;
      cor[j][i] = c;
    }
  }

  return cor;
};

exports.zScores = function zScores(matrix, means, standardDeviations) {
  means = means || exports.mean(matrix);
  if (typeof standardDeviations === 'undefined') standardDeviations = exports.standardDeviation(matrix, true, means);
  return exports.standardize(exports.center(matrix, means, false), standardDeviations, true);
};

exports.center = function center(matrix, means, inPlace) {
  means = means || exports.mean(matrix);
  var result = matrix,
      l = matrix.length,
      i,
      j,
      jj;

  if (!inPlace) {
    result = new Array(l);

    for (i = 0; i < l; i++) {
      result[i] = new Array(matrix[i].length);
    }
  }

  for (i = 0; i < l; i++) {
    var row = result[i];

    for (j = 0, jj = row.length; j < jj; j++) {
      row[j] = matrix[i][j] - means[j];
    }
  }

  return result;
};

exports.standardize = function standardize(matrix, standardDeviations, inPlace) {
  if (typeof standardDeviations === 'undefined') standardDeviations = exports.standardDeviation(matrix);
  var result = matrix,
      l = matrix.length,
      i,
      j,
      jj;

  if (!inPlace) {
    result = new Array(l);

    for (i = 0; i < l; i++) {
      result[i] = new Array(matrix[i].length);
    }
  }

  for (i = 0; i < l; i++) {
    var resultRow = result[i];
    var sourceRow = matrix[i];

    for (j = 0, jj = resultRow.length; j < jj; j++) {
      if (standardDeviations[j] !== 0 && !isNaN(standardDeviations[j])) {
        resultRow[j] = sourceRow[j] / standardDeviations[j];
      }
    }
  }

  return result;
};

exports.weightedVariance = function weightedVariance(matrix, weights) {
  var means = exports.mean(matrix);
  var rows = matrix.length;
  if (rows === 0) return [];
  var cols = matrix[0].length;
  var vari = new Array(cols);

  for (var j = 0; j < cols; j++) {
    var sum = 0;
    var a = 0,
        b = 0;

    for (var i = 0; i < rows; i++) {
      var z = matrix[i][j] - means[j];
      var w = weights[i];
      sum += w * (z * z);
      b += w;
      a += w * w;
    }

    vari[j] = sum * (b / (b * b - a));
  }

  return vari;
};

exports.weightedMean = function weightedMean(matrix, weights, dimension) {
  if (typeof dimension === 'undefined') {
    dimension = 0;
  }

  var rows = matrix.length;
  if (rows === 0) return [];
  var cols = matrix[0].length,
      means,
      i,
      ii,
      j,
      w,
      row;

  if (dimension === 0) {
    means = new Array(cols);

    for (i = 0; i < cols; i++) {
      means[i] = 0;
    }

    for (i = 0; i < rows; i++) {
      row = matrix[i];
      w = weights[i];

      for (j = 0; j < cols; j++) {
        means[j] += row[j] * w;
      }
    }
  } else if (dimension === 1) {
    means = new Array(rows);

    for (i = 0; i < rows; i++) {
      means[i] = 0;
    }

    for (j = 0; j < rows; j++) {
      row = matrix[j];
      w = weights[j];

      for (i = 0; i < cols; i++) {
        means[j] += row[i] * w;
      }
    }
  } else {
    throw new Error('Invalid dimension');
  }

  var weightSum = arrayStat.sum(weights);

  if (weightSum !== 0) {
    for (i = 0, ii = means.length; i < ii; i++) {
      means[i] /= weightSum;
    }
  }

  return means;
};

exports.weightedCovariance = function weightedCovariance(matrix, weights, means, dimension) {
  dimension = dimension || 0;
  means = means || exports.weightedMean(matrix, weights, dimension);
  var s1 = 0,
      s2 = 0;

  for (var i = 0, ii = weights.length; i < ii; i++) {
    s1 += weights[i];
    s2 += weights[i] * weights[i];
  }

  var factor = s1 / (s1 * s1 - s2);
  return exports.weightedScatter(matrix, weights, means, factor, dimension);
};

exports.weightedScatter = function weightedScatter(matrix, weights, means, factor, dimension) {
  dimension = dimension || 0;
  means = means || exports.weightedMean(matrix, weights, dimension);

  if (typeof factor === 'undefined') {
    factor = 1;
  }

  var rows = matrix.length;

  if (rows === 0) {
    return [[]];
  }

  var cols = matrix[0].length,
      cov,
      i,
      j,
      k,
      s;

  if (dimension === 0) {
    cov = new Array(cols);

    for (i = 0; i < cols; i++) {
      cov[i] = new Array(cols);
    }

    for (i = 0; i < cols; i++) {
      for (j = i; j < cols; j++) {
        s = 0;

        for (k = 0; k < rows; k++) {
          s += weights[k] * (matrix[k][j] - means[j]) * (matrix[k][i] - means[i]);
        }

        cov[i][j] = s * factor;
        cov[j][i] = s * factor;
      }
    }
  } else if (dimension === 1) {
    cov = new Array(rows);

    for (i = 0; i < rows; i++) {
      cov[i] = new Array(rows);
    }

    for (i = 0; i < rows; i++) {
      for (j = i; j < rows; j++) {
        s = 0;

        for (k = 0; k < cols; k++) {
          s += weights[k] * (matrix[j][k] - means[j]) * (matrix[i][k] - means[i]);
        }

        cov[i][j] = s * factor;
        cov[j][i] = s * factor;
      }
    }
  } else {
    throw new Error('Invalid dimension');
  }

  return cov;
};

/***/ }),
/* 40 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";

/**
 *
 * Function that returns a Number array of equally spaced numberOfPoints
 * containing a representation of intensities of the spectra arguments x
 * and y.
 *
 * The options parameter contains an object in the following form:
 * from: starting point
 * to: last point
 * numberOfPoints: number of points between from and to
 * variant: "slot" or "smooth" - smooth is the default option
 *
 * The slot variant consist that each point in the new array is calculated
 * averaging the existing points between the slot that belongs to the current
 * value. The smooth variant is the same but takes the integral of the range
 * of the slot and divide by the step size between two points in the new array.
 *
 * @param x - sorted increasing x values
 * @param y
 * @param options
 * @returns {Array} new array with the equally spaced data.
 *
 */

function getEquallySpacedData(x, y, options) {
  if (x.length > 1 && x[0] > x[1]) {
    x = x.slice().reverse();
    y = y.slice().reverse();
  }

  var xLength = x.length;
  if (xLength !== y.length) throw new RangeError("the x and y vector doesn't have the same size.");
  if (options === undefined) options = {};
  var from = options.from === undefined ? x[0] : options.from;

  if (isNaN(from) || !isFinite(from)) {
    throw new RangeError("'From' value must be a number");
  }

  var to = options.to === undefined ? x[x.length - 1] : options.to;

  if (isNaN(to) || !isFinite(to)) {
    throw new RangeError("'To' value must be a number");
  }

  var reverse = from > to;

  if (reverse) {
    var temp = from;
    from = to;
    to = temp;
  }

  var numberOfPoints = options.numberOfPoints === undefined ? 100 : options.numberOfPoints;

  if (isNaN(numberOfPoints) || !isFinite(numberOfPoints)) {
    throw new RangeError("'Number of points' value must be a number");
  }

  if (numberOfPoints < 1) throw new RangeError("the number of point must be higher than 1");
  var algorithm = options.variant === "slot" ? "slot" : "smooth"; // default value: smooth

  var output = algorithm === "slot" ? getEquallySpacedSlot(x, y, from, to, numberOfPoints) : getEquallySpacedSmooth(x, y, from, to, numberOfPoints);
  return reverse ? output.reverse() : output;
}
/**
 * function that retrieves the getEquallySpacedData with the variant "smooth"
 *
 * @param x
 * @param y
 * @param from - Initial point
 * @param to - Final point
 * @param numberOfPoints
 * @returns {Array} - Array of y's equally spaced with the variant "smooth"
 */


function getEquallySpacedSmooth(x, y, from, to, numberOfPoints) {
  var xLength = x.length;
  var step = (to - from) / (numberOfPoints - 1);
  var halfStep = step / 2;
  var start = from - halfStep;
  var output = new Array(numberOfPoints);
  var initialOriginalStep = x[1] - x[0];
  var lastOriginalStep = x[x.length - 1] - x[x.length - 2]; // Init main variables

  var min = start;
  var max = start + step;
  var previousX = Number.MIN_VALUE;
  var previousY = 0;
  var nextX = x[0] - initialOriginalStep;
  var nextY = 0;
  var currentValue = 0;
  var slope = 0;
  var intercept = 0;
  var sumAtMin = 0;
  var sumAtMax = 0;
  var i = 0; // index of input

  var j = 0; // index of output

  function getSlope(x0, y0, x1, y1) {
    return (y1 - y0) / (x1 - x0);
  }

  main: while (true) {
    while (nextX - max >= 0) {
      // no overlap with original point, just consume current value
      var add = integral(0, max - previousX, slope, previousY);
      sumAtMax = currentValue + add;
      output[j] = (sumAtMax - sumAtMin) / step;
      j++;
      if (j === numberOfPoints) break main;
      min = max;
      max += step;
      sumAtMin = sumAtMax;
    }

    if (previousX <= min && min <= nextX) {
      add = integral(0, min - previousX, slope, previousY);
      sumAtMin = currentValue + add;
    }

    currentValue += integral(previousX, nextX, slope, intercept);
    previousX = nextX;
    previousY = nextY;

    if (i < xLength) {
      nextX = x[i];
      nextY = y[i];
      i++;
    } else if (i === xLength) {
      nextX += lastOriginalStep;
      nextY = 0;
    } // updating parameters


    slope = getSlope(previousX, previousY, nextX, nextY);
    intercept = -slope * previousX + previousY;
  }

  return output;
}
/**
 * function that retrieves the getEquallySpacedData with the variant "slot"
 *
 * @param x
 * @param y
 * @param from - Initial point
 * @param to - Final point
 * @param numberOfPoints
 * @returns {Array} - Array of y's equally spaced with the variant "slot"
 */


function getEquallySpacedSlot(x, y, from, to, numberOfPoints) {
  var xLength = x.length;
  var step = (to - from) / (numberOfPoints - 1);
  var halfStep = step / 2;
  var lastStep = x[x.length - 1] - x[x.length - 2];
  var start = from - halfStep;
  var output = new Array(numberOfPoints); // Init main variables

  var min = start;
  var max = start + step;
  var previousX = -Number.MAX_VALUE;
  var previousY = 0;
  var nextX = x[0];
  var nextY = y[0];
  var frontOutsideSpectra = 0;
  var backOutsideSpectra = true;
  var currentValue = 0; // for slot algorithm

  var currentPoints = 0;
  var i = 1; // index of input

  var j = 0; // index of output

  main: while (true) {
    if (previousX >= nextX) throw new Error('x must be an increasing serie');

    while (previousX - max > 0) {
      // no overlap with original point, just consume current value
      if (backOutsideSpectra) {
        currentPoints++;
        backOutsideSpectra = false;
      }

      output[j] = currentPoints <= 0 ? 0 : currentValue / currentPoints;
      j++;
      if (j === numberOfPoints) break main;
      min = max;
      max += step;
      currentValue = 0;
      currentPoints = 0;
    }

    if (previousX > min) {
      currentValue += previousY;
      currentPoints++;
    }

    if (previousX === -Number.MAX_VALUE || frontOutsideSpectra > 1) currentPoints--;
    previousX = nextX;
    previousY = nextY;

    if (i < xLength) {
      nextX = x[i];
      nextY = y[i];
      i++;
    } else {
      nextX += lastStep;
      nextY = 0;
      frontOutsideSpectra++;
    }
  }

  return output;
}
/**
 * Function that calculates the integral of the line between two
 * x-coordinates, given the slope and intercept of the line.
 *
 * @param x0
 * @param x1
 * @param slope
 * @param intercept
 * @returns {number} integral value.
 */


function integral(x0, x1, slope, intercept) {
  return 0.5 * slope * x1 * x1 + intercept * x1 - (0.5 * slope * x0 * x0 + intercept * x0);
}

exports.getEquallySpacedData = getEquallySpacedData;
exports.integral = integral;

/***/ }),
/* 41 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


exports.SNV = SNV;

var Stat = __webpack_require__(11).array;
/**
 * Function that applies the standard normal variate (SNV) to an array of values.
 *
 * @param data - Array of values.
 * @returns {Array} - applied the SNV.
 */


function SNV(data) {
  var mean = Stat.mean(data);
  var std = Stat.standardDeviation(data);
  var result = data.slice();

  for (var i = 0; i < data.length; i++) {
    result[i] = (result[i] - mean) / std;
  }

  return result;
}

/***/ }),
/* 42 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


var BaseView = __webpack_require__(1);

class MatrixTransposeView extends BaseView {
  constructor(matrix) {
    super(matrix, matrix.columns, matrix.rows);
  }

  set(rowIndex, columnIndex, value) {
    this.matrix.set(columnIndex, rowIndex, value);
    return this;
  }

  get(rowIndex, columnIndex) {
    return this.matrix.get(columnIndex, rowIndex);
  }

}

module.exports = MatrixTransposeView;

/***/ }),
/* 43 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


var BaseView = __webpack_require__(1);

class MatrixRowView extends BaseView {
  constructor(matrix, row) {
    super(matrix, 1, matrix.columns);
    this.row = row;
  }

  set(rowIndex, columnIndex, value) {
    this.matrix.set(this.row, columnIndex, value);
    return this;
  }

  get(rowIndex, columnIndex) {
    return this.matrix.get(this.row, columnIndex);
  }

}

module.exports = MatrixRowView;

/***/ }),
/* 44 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


var BaseView = __webpack_require__(1);

var util = __webpack_require__(3);

class MatrixSubView extends BaseView {
  constructor(matrix, startRow, endRow, startColumn, endColumn) {
    util.checkRange(matrix, startRow, endRow, startColumn, endColumn);
    super(matrix, endRow - startRow + 1, endColumn - startColumn + 1);
    this.startRow = startRow;
    this.startColumn = startColumn;
  }

  set(rowIndex, columnIndex, value) {
    this.matrix.set(this.startRow + rowIndex, this.startColumn + columnIndex, value);
    return this;
  }

  get(rowIndex, columnIndex) {
    return this.matrix.get(this.startRow + rowIndex, this.startColumn + columnIndex);
  }

}

module.exports = MatrixSubView;

/***/ }),
/* 45 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


var BaseView = __webpack_require__(1);

var util = __webpack_require__(3);

class MatrixSelectionView extends BaseView {
  constructor(matrix, rowIndices, columnIndices) {
    var indices = util.checkIndices(matrix, rowIndices, columnIndices);
    super(matrix, indices.row.length, indices.column.length);
    this.rowIndices = indices.row;
    this.columnIndices = indices.column;
  }

  set(rowIndex, columnIndex, value) {
    this.matrix.set(this.rowIndices[rowIndex], this.columnIndices[columnIndex], value);
    return this;
  }

  get(rowIndex, columnIndex) {
    return this.matrix.get(this.rowIndices[rowIndex], this.columnIndices[columnIndex]);
  }

}

module.exports = MatrixSelectionView;

/***/ }),
/* 46 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


var BaseView = __webpack_require__(1);

class MatrixColumnView extends BaseView {
  constructor(matrix, column) {
    super(matrix, matrix.rows, 1);
    this.column = column;
  }

  set(rowIndex, columnIndex, value) {
    this.matrix.set(rowIndex, this.column, value);
    return this;
  }

  get(rowIndex) {
    return this.matrix.get(rowIndex, this.column);
  }

}

module.exports = MatrixColumnView;

/***/ }),
/* 47 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


var BaseView = __webpack_require__(1);

class MatrixFlipRowView extends BaseView {
  constructor(matrix) {
    super(matrix, matrix.rows, matrix.columns);
  }

  set(rowIndex, columnIndex, value) {
    this.matrix.set(this.rows - rowIndex - 1, columnIndex, value);
    return this;
  }

  get(rowIndex, columnIndex) {
    return this.matrix.get(this.rows - rowIndex - 1, columnIndex);
  }

}

module.exports = MatrixFlipRowView;

/***/ }),
/* 48 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


var BaseView = __webpack_require__(1);

class MatrixFlipColumnView extends BaseView {
  constructor(matrix) {
    super(matrix, matrix.rows, matrix.columns);
  }

  set(rowIndex, columnIndex, value) {
    this.matrix.set(rowIndex, this.columns - columnIndex - 1, value);
    return this;
  }

  get(rowIndex, columnIndex) {
    return this.matrix.get(rowIndex, this.columns - columnIndex - 1);
  }

}

module.exports = MatrixFlipColumnView;

/***/ }),
/* 49 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


var Matrix = __webpack_require__(0).Matrix;

var SingularValueDecomposition = __webpack_require__(10);

var EigenvalueDecomposition = __webpack_require__(50);

var LuDecomposition = __webpack_require__(9);

var QrDecomposition = __webpack_require__(51);

var CholeskyDecomposition = __webpack_require__(52);

function inverse(matrix) {
  matrix = Matrix.checkMatrix(matrix);
  return solve(matrix, Matrix.eye(matrix.rows));
}
/**
 * Returns the inverse
 * @memberOf Matrix
 * @static
 * @param {Matrix} matrix
 * @return {Matrix} matrix
 * @alias inv
 */


Matrix.inverse = Matrix.inv = inverse;
/**
 * Returns the inverse
 * @memberOf Matrix
 * @static
 * @param {Matrix} matrix
 * @return {Matrix} matrix
 * @alias inv
 */

Matrix.prototype.inverse = Matrix.prototype.inv = function () {
  return inverse(this);
};

function solve(leftHandSide, rightHandSide) {
  leftHandSide = Matrix.checkMatrix(leftHandSide);
  rightHandSide = Matrix.checkMatrix(rightHandSide);
  return leftHandSide.isSquare() ? new LuDecomposition(leftHandSide).solve(rightHandSide) : new QrDecomposition(leftHandSide).solve(rightHandSide);
}

Matrix.solve = solve;

Matrix.prototype.solve = function (other) {
  return solve(this, other);
};

module.exports = {
  SingularValueDecomposition: SingularValueDecomposition,
  SVD: SingularValueDecomposition,
  EigenvalueDecomposition: EigenvalueDecomposition,
  EVD: EigenvalueDecomposition,
  LuDecomposition: LuDecomposition,
  LU: LuDecomposition,
  QrDecomposition: QrDecomposition,
  QR: QrDecomposition,
  CholeskyDecomposition: CholeskyDecomposition,
  CHO: CholeskyDecomposition,
  inverse: inverse,
  solve: solve
};

/***/ }),
/* 50 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


const Matrix = __webpack_require__(0).Matrix;

const util = __webpack_require__(5);

const hypotenuse = util.hypotenuse;
const getFilled2DArray = util.getFilled2DArray;
const defaultOptions = {
  assumeSymmetric: false
}; // https://github.com/lutzroeder/Mapack/blob/master/Source/EigenvalueDecomposition.cs

function EigenvalueDecomposition(matrix, options) {
  options = Object.assign({}, defaultOptions, options);

  if (!(this instanceof EigenvalueDecomposition)) {
    return new EigenvalueDecomposition(matrix, options);
  }

  matrix = Matrix.checkMatrix(matrix);

  if (!matrix.isSquare()) {
    throw new Error('Matrix is not a square matrix');
  }

  var n = matrix.columns,
      V = getFilled2DArray(n, n, 0),
      d = new Array(n),
      e = new Array(n),
      value = matrix,
      i,
      j;
  var isSymmetric = false;

  if (options.assumeSymmetric) {
    isSymmetric = true;
  } else {
    isSymmetric = matrix.isSymmetric();
  }

  if (isSymmetric) {
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        V[i][j] = value.get(i, j);
      }
    }

    tred2(n, e, d, V);
    tql2(n, e, d, V);
  } else {
    var H = getFilled2DArray(n, n, 0),
        ort = new Array(n);

    for (j = 0; j < n; j++) {
      for (i = 0; i < n; i++) {
        H[i][j] = value.get(i, j);
      }
    }

    orthes(n, H, ort, V);
    hqr2(n, e, d, V, H);
  }

  this.n = n;
  this.e = e;
  this.d = d;
  this.V = V;
}

EigenvalueDecomposition.prototype = {
  get realEigenvalues() {
    return this.d;
  },

  get imaginaryEigenvalues() {
    return this.e;
  },

  get eigenvectorMatrix() {
    if (!Matrix.isMatrix(this.V)) {
      this.V = new Matrix(this.V);
    }

    return this.V;
  },

  get diagonalMatrix() {
    var n = this.n,
        e = this.e,
        d = this.d,
        X = new Matrix(n, n),
        i,
        j;

    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        X[i][j] = 0;
      }

      X[i][i] = d[i];

      if (e[i] > 0) {
        X[i][i + 1] = e[i];
      } else if (e[i] < 0) {
        X[i][i - 1] = e[i];
      }
    }

    return X;
  }

};

function tred2(n, e, d, V) {
  var f, g, h, i, j, k, hh, scale;

  for (j = 0; j < n; j++) {
    d[j] = V[n - 1][j];
  }

  for (i = n - 1; i > 0; i--) {
    scale = 0;
    h = 0;

    for (k = 0; k < i; k++) {
      scale = scale + Math.abs(d[k]);
    }

    if (scale === 0) {
      e[i] = d[i - 1];

      for (j = 0; j < i; j++) {
        d[j] = V[i - 1][j];
        V[i][j] = 0;
        V[j][i] = 0;
      }
    } else {
      for (k = 0; k < i; k++) {
        d[k] /= scale;
        h += d[k] * d[k];
      }

      f = d[i - 1];
      g = Math.sqrt(h);

      if (f > 0) {
        g = -g;
      }

      e[i] = scale * g;
      h = h - f * g;
      d[i - 1] = f - g;

      for (j = 0; j < i; j++) {
        e[j] = 0;
      }

      for (j = 0; j < i; j++) {
        f = d[j];
        V[j][i] = f;
        g = e[j] + V[j][j] * f;

        for (k = j + 1; k <= i - 1; k++) {
          g += V[k][j] * d[k];
          e[k] += V[k][j] * f;
        }

        e[j] = g;
      }

      f = 0;

      for (j = 0; j < i; j++) {
        e[j] /= h;
        f += e[j] * d[j];
      }

      hh = f / (h + h);

      for (j = 0; j < i; j++) {
        e[j] -= hh * d[j];
      }

      for (j = 0; j < i; j++) {
        f = d[j];
        g = e[j];

        for (k = j; k <= i - 1; k++) {
          V[k][j] -= f * e[k] + g * d[k];
        }

        d[j] = V[i - 1][j];
        V[i][j] = 0;
      }
    }

    d[i] = h;
  }

  for (i = 0; i < n - 1; i++) {
    V[n - 1][i] = V[i][i];
    V[i][i] = 1;
    h = d[i + 1];

    if (h !== 0) {
      for (k = 0; k <= i; k++) {
        d[k] = V[k][i + 1] / h;
      }

      for (j = 0; j <= i; j++) {
        g = 0;

        for (k = 0; k <= i; k++) {
          g += V[k][i + 1] * V[k][j];
        }

        for (k = 0; k <= i; k++) {
          V[k][j] -= g * d[k];
        }
      }
    }

    for (k = 0; k <= i; k++) {
      V[k][i + 1] = 0;
    }
  }

  for (j = 0; j < n; j++) {
    d[j] = V[n - 1][j];
    V[n - 1][j] = 0;
  }

  V[n - 1][n - 1] = 1;
  e[0] = 0;
}

function tql2(n, e, d, V) {
  var g, h, i, j, k, l, m, p, r, dl1, c, c2, c3, el1, s, s2, iter;

  for (i = 1; i < n; i++) {
    e[i - 1] = e[i];
  }

  e[n - 1] = 0;
  var f = 0,
      tst1 = 0,
      eps = Math.pow(2, -52);

  for (l = 0; l < n; l++) {
    tst1 = Math.max(tst1, Math.abs(d[l]) + Math.abs(e[l]));
    m = l;

    while (m < n) {
      if (Math.abs(e[m]) <= eps * tst1) {
        break;
      }

      m++;
    }

    if (m > l) {
      iter = 0;

      do {
        iter = iter + 1;
        g = d[l];
        p = (d[l + 1] - g) / (2 * e[l]);
        r = hypotenuse(p, 1);

        if (p < 0) {
          r = -r;
        }

        d[l] = e[l] / (p + r);
        d[l + 1] = e[l] * (p + r);
        dl1 = d[l + 1];
        h = g - d[l];

        for (i = l + 2; i < n; i++) {
          d[i] -= h;
        }

        f = f + h;
        p = d[m];
        c = 1;
        c2 = c;
        c3 = c;
        el1 = e[l + 1];
        s = 0;
        s2 = 0;

        for (i = m - 1; i >= l; i--) {
          c3 = c2;
          c2 = c;
          s2 = s;
          g = c * e[i];
          h = c * p;
          r = hypotenuse(p, e[i]);
          e[i + 1] = s * r;
          s = e[i] / r;
          c = p / r;
          p = c * d[i] - s * g;
          d[i + 1] = h + s * (c * g + s * d[i]);

          for (k = 0; k < n; k++) {
            h = V[k][i + 1];
            V[k][i + 1] = s * V[k][i] + c * h;
            V[k][i] = c * V[k][i] - s * h;
          }
        }

        p = -s * s2 * c3 * el1 * e[l] / dl1;
        e[l] = s * p;
        d[l] = c * p;
      } while (Math.abs(e[l]) > eps * tst1);
    }

    d[l] = d[l] + f;
    e[l] = 0;
  }

  for (i = 0; i < n - 1; i++) {
    k = i;
    p = d[i];

    for (j = i + 1; j < n; j++) {
      if (d[j] < p) {
        k = j;
        p = d[j];
      }
    }

    if (k !== i) {
      d[k] = d[i];
      d[i] = p;

      for (j = 0; j < n; j++) {
        p = V[j][i];
        V[j][i] = V[j][k];
        V[j][k] = p;
      }
    }
  }
}

function orthes(n, H, ort, V) {
  var low = 0,
      high = n - 1,
      f,
      g,
      h,
      i,
      j,
      m,
      scale;

  for (m = low + 1; m <= high - 1; m++) {
    scale = 0;

    for (i = m; i <= high; i++) {
      scale = scale + Math.abs(H[i][m - 1]);
    }

    if (scale !== 0) {
      h = 0;

      for (i = high; i >= m; i--) {
        ort[i] = H[i][m - 1] / scale;
        h += ort[i] * ort[i];
      }

      g = Math.sqrt(h);

      if (ort[m] > 0) {
        g = -g;
      }

      h = h - ort[m] * g;
      ort[m] = ort[m] - g;

      for (j = m; j < n; j++) {
        f = 0;

        for (i = high; i >= m; i--) {
          f += ort[i] * H[i][j];
        }

        f = f / h;

        for (i = m; i <= high; i++) {
          H[i][j] -= f * ort[i];
        }
      }

      for (i = 0; i <= high; i++) {
        f = 0;

        for (j = high; j >= m; j--) {
          f += ort[j] * H[i][j];
        }

        f = f / h;

        for (j = m; j <= high; j++) {
          H[i][j] -= f * ort[j];
        }
      }

      ort[m] = scale * ort[m];
      H[m][m - 1] = scale * g;
    }
  }

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      V[i][j] = i === j ? 1 : 0;
    }
  }

  for (m = high - 1; m >= low + 1; m--) {
    if (H[m][m - 1] !== 0) {
      for (i = m + 1; i <= high; i++) {
        ort[i] = H[i][m - 1];
      }

      for (j = m; j <= high; j++) {
        g = 0;

        for (i = m; i <= high; i++) {
          g += ort[i] * V[i][j];
        }

        g = g / ort[m] / H[m][m - 1];

        for (i = m; i <= high; i++) {
          V[i][j] += g * ort[i];
        }
      }
    }
  }
}

function hqr2(nn, e, d, V, H) {
  var n = nn - 1,
      low = 0,
      high = nn - 1,
      eps = Math.pow(2, -52),
      exshift = 0,
      norm = 0,
      p = 0,
      q = 0,
      r = 0,
      s = 0,
      z = 0,
      iter = 0,
      i,
      j,
      k,
      l,
      m,
      t,
      w,
      x,
      y,
      ra,
      sa,
      vr,
      vi,
      notlast,
      cdivres;

  for (i = 0; i < nn; i++) {
    if (i < low || i > high) {
      d[i] = H[i][i];
      e[i] = 0;
    }

    for (j = Math.max(i - 1, 0); j < nn; j++) {
      norm = norm + Math.abs(H[i][j]);
    }
  }

  while (n >= low) {
    l = n;

    while (l > low) {
      s = Math.abs(H[l - 1][l - 1]) + Math.abs(H[l][l]);

      if (s === 0) {
        s = norm;
      }

      if (Math.abs(H[l][l - 1]) < eps * s) {
        break;
      }

      l--;
    }

    if (l === n) {
      H[n][n] = H[n][n] + exshift;
      d[n] = H[n][n];
      e[n] = 0;
      n--;
      iter = 0;
    } else if (l === n - 1) {
      w = H[n][n - 1] * H[n - 1][n];
      p = (H[n - 1][n - 1] - H[n][n]) / 2;
      q = p * p + w;
      z = Math.sqrt(Math.abs(q));
      H[n][n] = H[n][n] + exshift;
      H[n - 1][n - 1] = H[n - 1][n - 1] + exshift;
      x = H[n][n];

      if (q >= 0) {
        z = p >= 0 ? p + z : p - z;
        d[n - 1] = x + z;
        d[n] = d[n - 1];

        if (z !== 0) {
          d[n] = x - w / z;
        }

        e[n - 1] = 0;
        e[n] = 0;
        x = H[n][n - 1];
        s = Math.abs(x) + Math.abs(z);
        p = x / s;
        q = z / s;
        r = Math.sqrt(p * p + q * q);
        p = p / r;
        q = q / r;

        for (j = n - 1; j < nn; j++) {
          z = H[n - 1][j];
          H[n - 1][j] = q * z + p * H[n][j];
          H[n][j] = q * H[n][j] - p * z;
        }

        for (i = 0; i <= n; i++) {
          z = H[i][n - 1];
          H[i][n - 1] = q * z + p * H[i][n];
          H[i][n] = q * H[i][n] - p * z;
        }

        for (i = low; i <= high; i++) {
          z = V[i][n - 1];
          V[i][n - 1] = q * z + p * V[i][n];
          V[i][n] = q * V[i][n] - p * z;
        }
      } else {
        d[n - 1] = x + p;
        d[n] = x + p;
        e[n - 1] = z;
        e[n] = -z;
      }

      n = n - 2;
      iter = 0;
    } else {
      x = H[n][n];
      y = 0;
      w = 0;

      if (l < n) {
        y = H[n - 1][n - 1];
        w = H[n][n - 1] * H[n - 1][n];
      }

      if (iter === 10) {
        exshift += x;

        for (i = low; i <= n; i++) {
          H[i][i] -= x;
        }

        s = Math.abs(H[n][n - 1]) + Math.abs(H[n - 1][n - 2]);
        x = y = 0.75 * s;
        w = -0.4375 * s * s;
      }

      if (iter === 30) {
        s = (y - x) / 2;
        s = s * s + w;

        if (s > 0) {
          s = Math.sqrt(s);

          if (y < x) {
            s = -s;
          }

          s = x - w / ((y - x) / 2 + s);

          for (i = low; i <= n; i++) {
            H[i][i] -= s;
          }

          exshift += s;
          x = y = w = 0.964;
        }
      }

      iter = iter + 1;
      m = n - 2;

      while (m >= l) {
        z = H[m][m];
        r = x - z;
        s = y - z;
        p = (r * s - w) / H[m + 1][m] + H[m][m + 1];
        q = H[m + 1][m + 1] - z - r - s;
        r = H[m + 2][m + 1];
        s = Math.abs(p) + Math.abs(q) + Math.abs(r);
        p = p / s;
        q = q / s;
        r = r / s;

        if (m === l) {
          break;
        }

        if (Math.abs(H[m][m - 1]) * (Math.abs(q) + Math.abs(r)) < eps * (Math.abs(p) * (Math.abs(H[m - 1][m - 1]) + Math.abs(z) + Math.abs(H[m + 1][m + 1])))) {
          break;
        }

        m--;
      }

      for (i = m + 2; i <= n; i++) {
        H[i][i - 2] = 0;

        if (i > m + 2) {
          H[i][i - 3] = 0;
        }
      }

      for (k = m; k <= n - 1; k++) {
        notlast = k !== n - 1;

        if (k !== m) {
          p = H[k][k - 1];
          q = H[k + 1][k - 1];
          r = notlast ? H[k + 2][k - 1] : 0;
          x = Math.abs(p) + Math.abs(q) + Math.abs(r);

          if (x !== 0) {
            p = p / x;
            q = q / x;
            r = r / x;
          }
        }

        if (x === 0) {
          break;
        }

        s = Math.sqrt(p * p + q * q + r * r);

        if (p < 0) {
          s = -s;
        }

        if (s !== 0) {
          if (k !== m) {
            H[k][k - 1] = -s * x;
          } else if (l !== m) {
            H[k][k - 1] = -H[k][k - 1];
          }

          p = p + s;
          x = p / s;
          y = q / s;
          z = r / s;
          q = q / p;
          r = r / p;

          for (j = k; j < nn; j++) {
            p = H[k][j] + q * H[k + 1][j];

            if (notlast) {
              p = p + r * H[k + 2][j];
              H[k + 2][j] = H[k + 2][j] - p * z;
            }

            H[k][j] = H[k][j] - p * x;
            H[k + 1][j] = H[k + 1][j] - p * y;
          }

          for (i = 0; i <= Math.min(n, k + 3); i++) {
            p = x * H[i][k] + y * H[i][k + 1];

            if (notlast) {
              p = p + z * H[i][k + 2];
              H[i][k + 2] = H[i][k + 2] - p * r;
            }

            H[i][k] = H[i][k] - p;
            H[i][k + 1] = H[i][k + 1] - p * q;
          }

          for (i = low; i <= high; i++) {
            p = x * V[i][k] + y * V[i][k + 1];

            if (notlast) {
              p = p + z * V[i][k + 2];
              V[i][k + 2] = V[i][k + 2] - p * r;
            }

            V[i][k] = V[i][k] - p;
            V[i][k + 1] = V[i][k + 1] - p * q;
          }
        }
      }
    }
  }

  if (norm === 0) {
    return;
  }

  for (n = nn - 1; n >= 0; n--) {
    p = d[n];
    q = e[n];

    if (q === 0) {
      l = n;
      H[n][n] = 1;

      for (i = n - 1; i >= 0; i--) {
        w = H[i][i] - p;
        r = 0;

        for (j = l; j <= n; j++) {
          r = r + H[i][j] * H[j][n];
        }

        if (e[i] < 0) {
          z = w;
          s = r;
        } else {
          l = i;

          if (e[i] === 0) {
            H[i][n] = w !== 0 ? -r / w : -r / (eps * norm);
          } else {
            x = H[i][i + 1];
            y = H[i + 1][i];
            q = (d[i] - p) * (d[i] - p) + e[i] * e[i];
            t = (x * s - z * r) / q;
            H[i][n] = t;
            H[i + 1][n] = Math.abs(x) > Math.abs(z) ? (-r - w * t) / x : (-s - y * t) / z;
          }

          t = Math.abs(H[i][n]);

          if (eps * t * t > 1) {
            for (j = i; j <= n; j++) {
              H[j][n] = H[j][n] / t;
            }
          }
        }
      }
    } else if (q < 0) {
      l = n - 1;

      if (Math.abs(H[n][n - 1]) > Math.abs(H[n - 1][n])) {
        H[n - 1][n - 1] = q / H[n][n - 1];
        H[n - 1][n] = -(H[n][n] - p) / H[n][n - 1];
      } else {
        cdivres = cdiv(0, -H[n - 1][n], H[n - 1][n - 1] - p, q);
        H[n - 1][n - 1] = cdivres[0];
        H[n - 1][n] = cdivres[1];
      }

      H[n][n - 1] = 0;
      H[n][n] = 1;

      for (i = n - 2; i >= 0; i--) {
        ra = 0;
        sa = 0;

        for (j = l; j <= n; j++) {
          ra = ra + H[i][j] * H[j][n - 1];
          sa = sa + H[i][j] * H[j][n];
        }

        w = H[i][i] - p;

        if (e[i] < 0) {
          z = w;
          r = ra;
          s = sa;
        } else {
          l = i;

          if (e[i] === 0) {
            cdivres = cdiv(-ra, -sa, w, q);
            H[i][n - 1] = cdivres[0];
            H[i][n] = cdivres[1];
          } else {
            x = H[i][i + 1];
            y = H[i + 1][i];
            vr = (d[i] - p) * (d[i] - p) + e[i] * e[i] - q * q;
            vi = (d[i] - p) * 2 * q;

            if (vr === 0 && vi === 0) {
              vr = eps * norm * (Math.abs(w) + Math.abs(q) + Math.abs(x) + Math.abs(y) + Math.abs(z));
            }

            cdivres = cdiv(x * r - z * ra + q * sa, x * s - z * sa - q * ra, vr, vi);
            H[i][n - 1] = cdivres[0];
            H[i][n] = cdivres[1];

            if (Math.abs(x) > Math.abs(z) + Math.abs(q)) {
              H[i + 1][n - 1] = (-ra - w * H[i][n - 1] + q * H[i][n]) / x;
              H[i + 1][n] = (-sa - w * H[i][n] - q * H[i][n - 1]) / x;
            } else {
              cdivres = cdiv(-r - y * H[i][n - 1], -s - y * H[i][n], z, q);
              H[i + 1][n - 1] = cdivres[0];
              H[i + 1][n] = cdivres[1];
            }
          }

          t = Math.max(Math.abs(H[i][n - 1]), Math.abs(H[i][n]));

          if (eps * t * t > 1) {
            for (j = i; j <= n; j++) {
              H[j][n - 1] = H[j][n - 1] / t;
              H[j][n] = H[j][n] / t;
            }
          }
        }
      }
    }
  }

  for (i = 0; i < nn; i++) {
    if (i < low || i > high) {
      for (j = i; j < nn; j++) {
        V[i][j] = H[i][j];
      }
    }
  }

  for (j = nn - 1; j >= low; j--) {
    for (i = low; i <= high; i++) {
      z = 0;

      for (k = low; k <= Math.min(j, high); k++) {
        z = z + V[i][k] * H[k][j];
      }

      V[i][j] = z;
    }
  }
}

function cdiv(xr, xi, yr, yi) {
  var r, d;

  if (Math.abs(yr) > Math.abs(yi)) {
    r = yi / yr;
    d = yr + r * yi;
    return [(xr + r * xi) / d, (xi - r * xr) / d];
  } else {
    r = yr / yi;
    d = yi + r * yr;
    return [(r * xr + xi) / d, (r * xi - xr) / d];
  }
}

module.exports = EigenvalueDecomposition;

/***/ }),
/* 51 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


var Matrix = __webpack_require__(0).Matrix;

var hypotenuse = __webpack_require__(5).hypotenuse; //https://github.com/lutzroeder/Mapack/blob/master/Source/QrDecomposition.cs


function QrDecomposition(value) {
  if (!(this instanceof QrDecomposition)) {
    return new QrDecomposition(value);
  }

  value = Matrix.checkMatrix(value);
  var qr = value.clone(),
      m = value.rows,
      n = value.columns,
      rdiag = new Array(n),
      i,
      j,
      k,
      s;

  for (k = 0; k < n; k++) {
    var nrm = 0;

    for (i = k; i < m; i++) {
      nrm = hypotenuse(nrm, qr[i][k]);
    }

    if (nrm !== 0) {
      if (qr[k][k] < 0) {
        nrm = -nrm;
      }

      for (i = k; i < m; i++) {
        qr[i][k] /= nrm;
      }

      qr[k][k] += 1;

      for (j = k + 1; j < n; j++) {
        s = 0;

        for (i = k; i < m; i++) {
          s += qr[i][k] * qr[i][j];
        }

        s = -s / qr[k][k];

        for (i = k; i < m; i++) {
          qr[i][j] += s * qr[i][k];
        }
      }
    }

    rdiag[k] = -nrm;
  }

  this.QR = qr;
  this.Rdiag = rdiag;
}

QrDecomposition.prototype = {
  solve: function solve(value) {
    value = Matrix.checkMatrix(value);
    var qr = this.QR,
        m = qr.rows;

    if (value.rows !== m) {
      throw new Error('Matrix row dimensions must agree');
    }

    if (!this.isFullRank()) {
      throw new Error('Matrix is rank deficient');
    }

    var count = value.columns;
    var X = value.clone();
    var n = qr.columns;
    var i, j, k, s;

    for (k = 0; k < n; k++) {
      for (j = 0; j < count; j++) {
        s = 0;

        for (i = k; i < m; i++) {
          s += qr[i][k] * X[i][j];
        }

        s = -s / qr[k][k];

        for (i = k; i < m; i++) {
          X[i][j] += s * qr[i][k];
        }
      }
    }

    for (k = n - 1; k >= 0; k--) {
      for (j = 0; j < count; j++) {
        X[k][j] /= this.Rdiag[k];
      }

      for (i = 0; i < k; i++) {
        for (j = 0; j < count; j++) {
          X[i][j] -= X[k][j] * qr[i][k];
        }
      }
    }

    return X.subMatrix(0, n - 1, 0, count - 1);
  },
  isFullRank: function isFullRank() {
    var columns = this.QR.columns;

    for (var i = 0; i < columns; i++) {
      if (this.Rdiag[i] === 0) {
        return false;
      }
    }

    return true;
  },

  get upperTriangularMatrix() {
    var qr = this.QR,
        n = qr.columns,
        X = new Matrix(n, n),
        i,
        j;

    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        if (i < j) {
          X[i][j] = qr[i][j];
        } else if (i === j) {
          X[i][j] = this.Rdiag[i];
        } else {
          X[i][j] = 0;
        }
      }
    }

    return X;
  },

  get orthogonalMatrix() {
    var qr = this.QR,
        rows = qr.rows,
        columns = qr.columns,
        X = new Matrix(rows, columns),
        i,
        j,
        k,
        s;

    for (k = columns - 1; k >= 0; k--) {
      for (i = 0; i < rows; i++) {
        X[i][k] = 0;
      }

      X[k][k] = 1;

      for (j = k; j < columns; j++) {
        if (qr[k][k] !== 0) {
          s = 0;

          for (i = k; i < rows; i++) {
            s += qr[i][k] * X[i][j];
          }

          s = -s / qr[k][k];

          for (i = k; i < rows; i++) {
            X[i][j] += s * qr[i][k];
          }
        }
      }
    }

    return X;
  }

};
module.exports = QrDecomposition;

/***/ }),
/* 52 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


var Matrix = __webpack_require__(0).Matrix; // https://github.com/lutzroeder/Mapack/blob/master/Source/CholeskyDecomposition.cs


function CholeskyDecomposition(value) {
  if (!(this instanceof CholeskyDecomposition)) {
    return new CholeskyDecomposition(value);
  }

  value = Matrix.checkMatrix(value);

  if (!value.isSymmetric()) {
    throw new Error('Matrix is not symmetric');
  }

  var a = value,
      dimension = a.rows,
      l = new Matrix(dimension, dimension),
      positiveDefinite = true,
      i,
      j,
      k;

  for (j = 0; j < dimension; j++) {
    var Lrowj = l[j];
    var d = 0;

    for (k = 0; k < j; k++) {
      var Lrowk = l[k];
      var s = 0;

      for (i = 0; i < k; i++) {
        s += Lrowk[i] * Lrowj[i];
      }

      Lrowj[k] = s = (a[j][k] - s) / l[k][k];
      d = d + s * s;
    }

    d = a[j][j] - d;
    positiveDefinite &= d > 0;
    l[j][j] = Math.sqrt(Math.max(d, 0));

    for (k = j + 1; k < dimension; k++) {
      l[j][k] = 0;
    }
  }

  if (!positiveDefinite) {
    throw new Error('Matrix is not positive definite');
  }

  this.L = l;
}

CholeskyDecomposition.prototype = {
  get lowerTriangularMatrix() {
    return this.L;
  },

  solve: function solve(value) {
    value = Matrix.checkMatrix(value);
    var l = this.L,
        dimension = l.rows;

    if (value.rows !== dimension) {
      throw new Error('Matrix dimensions do not match');
    }

    var count = value.columns,
        B = value.clone(),
        i,
        j,
        k;

    for (k = 0; k < dimension; k++) {
      for (j = 0; j < count; j++) {
        for (i = 0; i < k; i++) {
          B[k][j] -= B[i][j] * l[k][i];
        }

        B[k][j] /= l[k][k];
      }
    }

    for (k = dimension - 1; k >= 0; k--) {
      for (j = 0; j < count; j++) {
        for (i = k + 1; i < dimension; i++) {
          B[k][j] -= B[i][j] * l[i][k];
        }

        B[k][j] /= l[k][k];
      }
    }

    return B;
  }
};
module.exports = CholeskyDecomposition;

/***/ }),
/* 53 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


module.exports = function getConnectivityMatrix(options) {
  // TODO remove this line when the bug is fixed ... (in OCL addImplicitHydrogens)
  // this.ensureHelperArrays(this.cHelperNeighbours)
  this.toMolfile();
  var nbAtoms = this.getAllAtoms();
  var result = new Array(nbAtoms);

  for (var i = 0; i < nbAtoms; i++) {
    result[i] = new Array(nbAtoms).fill(0);
    result[i][i] = 1;
  }

  for (var i = 0; i < nbAtoms; i++) {
    for (var j = 0; j < this.getAllConnAtoms(i); j++) {
      result[i][this.getConnAtom(i, j)] = 1;
    }
  }

  return result;
};

/***/ })
/******/ ]);
});
//# sourceMappingURL=nmr-auto-assignment.js.map