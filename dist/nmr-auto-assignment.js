(function webpackUniversalModuleDefinition(root, factory) {
	if(typeof exports === 'object' && typeof module === 'object')
		module.exports = factory();
	else if(typeof define === 'function' && define.amd)
		define([], factory);
	else if(typeof exports === 'object')
		exports["nmrAutoAssignment"] = factory();
	else
		root["nmrAutoAssignment"] = factory();
})(this, function() {
return /******/ (function(modules) { // webpackBootstrap
/******/ 	// The module cache
/******/ 	var installedModules = {};

/******/ 	// The require function
/******/ 	function __webpack_require__(moduleId) {

/******/ 		// Check if module is in cache
/******/ 		if(installedModules[moduleId])
/******/ 			return installedModules[moduleId].exports;

/******/ 		// Create a new module (and put it into the cache)
/******/ 		var module = installedModules[moduleId] = {
/******/ 			exports: {},
/******/ 			id: moduleId,
/******/ 			loaded: false
/******/ 		};

/******/ 		// Execute the module function
/******/ 		modules[moduleId].call(module.exports, module, module.exports, __webpack_require__);

/******/ 		// Flag the module as loaded
/******/ 		module.loaded = true;

/******/ 		// Return the exports of the module
/******/ 		return module.exports;
/******/ 	}


/******/ 	// expose the modules object (__webpack_modules__)
/******/ 	__webpack_require__.m = modules;

/******/ 	// expose the module cache
/******/ 	__webpack_require__.c = installedModules;

/******/ 	// __webpack_public_path__
/******/ 	__webpack_require__.p = "";

/******/ 	// Load entry module and return exports
/******/ 	return __webpack_require__(0);
/******/ })
/************************************************************************/
/******/ ([
/* 0 */
/***/ function(module, exports, __webpack_require__) {

	/**
	 * Created by acastillo on 7/5/16.
	 */
	const SpinSystem = __webpack_require__(1);
	const AutoAssigner = __webpack_require__(2);

	function autoAssign(entry, options){
	    if(entry.spectra.h1PeakList){
	        return assignmentFromPeakPicking(entry, options);
	    }
	    else{
	        return assignmentFromRaw(entry, options);
	    }
	}

	function assignmentFromRaw(entry, options){
	    var molfile = entry.molfile;
	    var spectra = entry.spectra;

	    var molecule=ACT.load(molfile);

	    molecule.expandHydrogens();

	    entry.molecule = molecule;
	    entry.diaIDs = molecule.getDiastereotopicAtomIDs();

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
	}

	function assignmentFromPeakPicking(entry, options) {
	    const predictor = options.predictor;
	    var molecule, diaIDs, molfile;

	    var spectra = entry.spectra;
	    if(!entry.molecule) {
	        molecule=ACT.load(entry.molfile);
	        molecule.expandHydrogens();
	        diaIDs=molecule.getDiastereotopicAtomIDs();

	        for (var j = 0; j < diaIDs.length; j++) {
	            diaIDs[j].nbEquivalent=diaIDs[j].atoms.length;
	        }

	        diaIDs.sort(function(a,b) {
	            if (a.element == b.element) {
	                return b.nbEquivalent-a.nbEquivalent;
	            }
	            return a.element<b.element?1:-1;
	        });
	        entry.molecule = molecule;
	        entry.diaIDs = diaIDs;
	        entry.diaID = molecule.toIDCode();
	    }
	    else {
	        molecule = entry.molecule;
	        diaIDs = entry.diaIDs;
	    }

	    //H1 prediction
	    var h1pred = predictor.predict(molecule, {group:true});
	    if(!h1pred || h1pred.length === 0)
	        return null;

	    var optionsError = {iteration:options.iteration || 1, learningRatio:options.learningRatio || 1};

	    for (var j=0; j<h1pred.length; j++) {
	        h1pred[j].error = getError(h1pred[j], optionsError);
	    }

	    h1pred.sort(function(a,b) {
	        if (a.atomLabel==b.atomLabel) {
	            return b.integral-a.integral;
	        }
	        return a.atomLabel<b.atomLabel?1:-1;
	    });

	    try{
	        const spinSystem = new SpinSystem(h1pred, spectra.h1PeakList);
	        const autoAssigner = new AutoAssigner(spinSystem, {minScore:1 ,maxSolutions:3000, errorCS:-1});
	        return autoAssigner.getAssignments();
	    }
	    catch(e){
	        console.log("Could not assign this molecule.");
	        return null;
	    }
	}

	function  getError(prediction, param){
	    //console.log(prediction)
	    //Never use predictions with less than 3 votes
	    if(prediction.std==0||prediction.ncs<3){
	        return 20;
	    }
	    else{
	        //factor is between 1 and +inf
	        //console.log(prediction.ncs+" "+(param.iteration+1)+" "+param.learningRatio);
	        var factor = 3*prediction.std/
	            (Math.pow(prediction.ncs,(param.iteration+1)*param.learningRatio));//(param.iteration+1)*param.learningRatio*h1pred[indexSignal].ncs;
	        return 3*prediction.std+factor;
	    }
	    return 20;
	}

	module.exports = autoAssign;


/***/ },
/* 1 */
/***/ function(module, exports) {

	/**
	 * Created by acastillo on 9/2/16.
	 */
	'use strict'
	const DEBUG  = true;
	class SpinSystem {
	    constructor(diaIDsArray, signalsArray, opt){
	        var options = Object.assign({}, opt);
	        this.diaIDsArray = diaIDsArray;
	        this.signalsArray = signalsArray;
	        this.cosy = options.cosySignals || null;
	        this.connCosy = options.cosyPaths || null;
	        this.connHmbc = options.hmbcPaths || null;
	        this.hmbc = options.hmbcSignals || null;
	        this.init();
	    }

	    init(){
	        const nDiaIds = this.diaIDsArray.length;
	        const nSignals = this.signalsArray.length;
	        const diaIDByAtomLabel = {};
	        const indexByAtomLabel = {};
	        var shiftsH = [];
	        var shiftsC = [];
	        var windowH = [];
	        var windowC = [];

	        var signals1D = this.signalsArray;
	        var nH = 0, nC = 0, i = 0;

	        try {
	            this.chemicalShiftsT = new Array(nDiaIds);
	            this.chemicalShiftsTError = new Array(nDiaIds);
	            this.diaList = new Array(nDiaIds);
	            var dia = null;
	            for(i = 0; i < nDiaIds; i++) {
	                dia = this.diaIDsArray[i];
	                if(diaIDByAtomLabel[dia.atomLabel]) {
	                    diaIDByAtomLabel[dia.atomLabel].push(dia.diaIDs[0]);
	                    indexByAtomLabel[dia.atomLabel].push(i);
	                }
	                else {
	                    diaIDByAtomLabel[dia.atomLabel] = [dia.diaIDs[0]];
	                    indexByAtomLabel[dia.atomLabel] = [i];
	                }
	                this.diaList[i] = dia.integral;
	                this.chemicalShiftsT[i] = dia.delta;
	                this.chemicalShiftsTError[i] = dia.error || 0;
	            }
	            // We can't have more signals than different protons in the molecule if the integral
	            // matches the nH
	            this.signals = new Array(nSignals);
	            this.chemicalShiftsE = new Array(nSignals);
	            this.signalsWidth = new Array(nSignals);

	            for(i = 0; i < nSignals; i++) {
	                var from = signals1D[i].from;
	                var to = signals1D[i].to;
	                this.chemicalShiftsE[i] = (from+to)/2.0;
	                this.signals[i] = Math.round(signals1D[i].integral);
	                shiftsH.push((from+to)/2.0);
	                windowH.push(Math.abs(from-to));
	                this.signalsWidth[i] = Math.abs(from-to);
	            }
	            //System.out.println(diaIDsH.size());
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
	        for(var i = shifts.length - 1; i >= 0; i--){
	            if(Math.abs(value-shifts[i]) < minDiff) {
	                minDiff = Math.abs(value - shifts[i]);
	                index = i;
	            }
	        }
	        if(windows) {
	            if(minDiff <= windows[index]/2)
	                return index;
	        }
	        if(minDiff <= Math.abs(resolution*4))
	            return index;

	        return -1;
	    }

	}

	module.exports = SpinSystem;


/***/ },
/* 2 */
/***/ function(module, exports, __webpack_require__) {

	'use strict'
	/**
	 * Created by acastillo on 9/2/16.
	 */
	const TreeSet = __webpack_require__(3);

	const defaultOptions = {minScore:1, maxSolutions: 100, errorCS:-1, onlyCount: false, timeout:20000};

	const DEBUG = false;

	class Assignment{
	    constructor(spinSystem, opt){
	        var options = Object.assign({}, defaultOptions, opt);
	        this.spinSystem = spinSystem;
	        this.minScore = options.minScore;
	        this.maxSolutions = options.maxSolutions;
	        this.errorCS = options.errorCS;
	        this.onlyCount = options.onlyCount;
	        this.timeout = options.timeout;
	        this.MAXERRORSHMBC = 1;

	        this.timeoutTerminated = 0;
	        this.score = 0;
	        this.nSolutions = 0;
	        this.nSteps = 0;
	        this.lowerBound = 0;

	        this.solutions = null;

	        this.comparator = function(a, b){
	            return a.score - b.score;
	        }
	    }

	    getAssignments(){
	        var date = new Date();
	        this.timeStart = date.getTime();
	        var i, j, k, nSignals, nDiaIDs;

	        if(DEBUG) console.log(this.spinSystem);

	        this.lowerBound = this.minScore;
	        do{
	            this.nSolutions = 0;
	            this.nSteps = 0;
	            this.solutions = new TreeSet(this.comparator);

	            if(this.spinSystem.hmbcE != null)
	                this.spinSystem.hmbcLines = {};
	            if(this.spinSystem.cosyE != null)
	                this.spinSystem.cosyLines = {};

	            nSignals = this.spinSystem.signals.length;
	            nDiaIDs = this.spinSystem.diaList.length;
	            this.scores = new Array(nSignals);
	            let partial = new Array(nSignals);

	            for (i = 0; i < nSignals; i++) {
	                this.scores[i] = 1;
	                partial[i] = [];
	            }

	            var diaMask = new Array(nDiaIDs);

	            for(i = diaMask.length - 1; i >= 0; i--)
	                diaMask[i] = true;

	            try {
	                this.exploreTreeRec(this.spinSystem.signals,
	                    this.spinSystem.diaList,nSignals - 1, nDiaIDs - 1, diaMask, partial);

	            } catch (e) {
	                console.log("Exception in assignment: " + e);
	            }
	            this.lowerBound -= 0.1;
	            if(DEBUG)	console.log("Decreasing lowerBound: " + this.lowerBound);
	        }while(this.solutions.isEmpty() && this.lowerBound >= 0.4);

	        //Format the result
	        var assignment = this.solutions.elements;
	        var nSolutions = this.solutions.length;
	        for(i = 0; i < nSolutions; i++){
	            var assignment = this.solutions.elements[i].assignment;
	            for(j = 0; j < nSignals; j++){
	                var diaIDs = assignment[j];
	                var tmp = new Array(diaIDs.length);
	                for(k = 0; k < diaIDs.length; k++){
	                    tmp[k] = this.spinSystem.diaIDsArray[diaIDs[k]].diaIDs[0];
	                }
	                assignment[j] = {signalID : this.spinSystem.signalsArray[j].signalID,
	                                 delta : Math.round(this.spinSystem.signalsArray[j].signal[0].delta*100)/100,
	                                 //integral: Math.round(this.spinSystem.signalsArray[j].integral*100)/100,
	                                 diaID : tmp}
	            }
	        }

	        return this.solutions.elements;
	    }

	    exploreTreeRec(signals, diaList, indexSignal, indexDia, diaMask, partial) {
	        //If this happens, we can assign this atom group to this signal
	        while(indexDia >= 0 && signals[indexSignal] >= diaList[indexDia]) {
	            //Force a return if the loop time is longer than the given timeout
	            const d = new Date();
	            if((d.getTime() - this.timeStart) > this.timeout){
	                this.timeoutTerminated=true;
	                return;
	            }
	            //We can speed up it by checking the chemical shift first
	            if(diaMask[indexDia] && this._isWithinCSRange(indexSignal, indexDia)) {
	                this.nSteps++;
	                const sizePartial = partial[indexSignal].length;
	                //Assign the atom group to the signal
	                diaMask[indexDia] = false;//Mark this atom group as used
	                partial[indexSignal][sizePartial] = indexDia;//Add the atom group index to the assignment list
	                signals[indexSignal] -= diaList[indexDia];//Subtract the group from signal integral
	                //If this signal is completely assigned, we should verify all the restrictions
	                if(signals[indexSignal] == 0){
	                    let keySum  = this._accomplishCounts(indexSignal, partial);
	                    //System.out.println("Accomplish count: "+keySum);
	                    if(keySum != 0){
	                        //Verify the restrictions. A good solution should give a high score
	                        this.score = this._solutionScore(partial, indexSignal, keySum);
	                        if(DEBUG) console.log(this.score+" "+partial);
	                        if(DEBUG) console.log("score: "+this.score);
	                        //This is a solution
	                        if(this.score>0){
	                            if(indexSignal == 0){//We found a new solution
	                                this.nSolutions++;
	                                var solution = {assignment: this._cloneArray(partial), score: this.score};
	                                if (this.solutions.length >= this.maxSolutions) {
	                                    if (this.score > this.solutions.last().score) {
	                                        this.solutions.pollLast();
	                                        this.solutions.add(solution);
	                                    }
	                                } else {
	                                    this.solutions.add(solution);
	                                }
	                            }
	                            else{
	                                //Each new signal that we assign will produce a new level on the tree.
	                                indexSignal--;//Lets go forward with the next signal
	                                indexDia = diaList.length;
	                                while(!diaMask[--indexDia]);
	                                //SolutionTree newLevel = new SolutionTree(solution);
	                                this.exploreTreeRec(signals, diaList, indexSignal, indexDia, diaMask, partial);
	                                indexSignal++;
	                            }
	                        }
	                    }

	                }
	                else{
	                    //It says that the signal should be assigned by combining 2 or more signals
	                    const previousIndexDia = indexDia;
	                    while(indexDia > 0 && !diaMask[--indexDia]);
	                    if(indexDia>=0)
	                        this.exploreTreeRec(signals, diaList, indexSignal, indexDia, diaMask, partial);
	                    indexDia = previousIndexDia;
	                }
	                //Deallocate this atom group to try the next one.
	                indexDia = partial[indexSignal].splice(sizePartial,1);//Add the atom group index to the assignment list
	                diaMask[indexDia] = true;//Mark this atom group as available
	                signals[indexSignal] += diaList[indexDia];//Subtract the group from signal integral
	                this.scores[indexSignal] = 1;
	            }
	            indexDia--;
	        }
	    }

	    _cloneArray(data){
	        return JSON.parse(JSON.stringify(data));
	    }

	    _isWithinCSRange(indexSignal, indexDia) {
	        if(this.spinSystem.chemicalShiftsE != null && this.spinSystem.chemicalShiftsT != null){
	            if(this.errorCS == 0)
	                return true;
	            var cfAtoms = this.spinSystem.chemicalShiftsT[indexDia];

	            if(cfAtoms==-9999999)
	                return true;
	            var cfSignal = this.spinSystem.chemicalShiftsE[indexSignal];
	            var error = this.spinSystem.chemicalShiftsTError[indexDia];
	            if(error === 0)
	                error = this.errorCS;

	            var csError = Math.abs(this.spinSystem.signalsWidth[indexSignal]/2.0+Math.abs(error));
	            if(Math.abs(cfSignal - cfAtoms) <= csError)
	                return true;
	            else
	                return false;
	        }
	        return true;
	    }

	    _accomplishCounts(indexSignal, partial){
	        //Check the chemical shift
	        var keySum = -1;
	        var keySumCOSY = 1;
	        var keySumHMBC = 1;
	        if(this.spinSystem.cosyE != null){
	            keySumCOSY = this._accomplishCount(indexSignal, partial[indexSigna+l],
	                this.spinSystem.cosyT, this.spinSystem.cosyE, this.spinSystem.cosyLines, true, MAXERRORSCOSY);
	            keySum = keySumCOSY;
	        }
	        if(keySum!=0){
	            if(this.spinSystem.hmbcE != null){
	                keySumHMBC = this._accomplishCount(indexSignal, partial[indexSignal], this.spinSystem.hmbcT,
	                    this.spinSystem.hmbcE, this.spinSystem.hmbcLines, false, MAXERRORSHMBC);
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
	    _accomplishCount(index, signals,theoretical, experimental, hashMap, isSymmetryc, maxErrors) {
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
	        var expLH=0;
	        if(this.spinSystem.cosyLines != null) {
	            expLH++;
	            score=this._cosyScore(partial, current, keySingalAsg);
	        }
	        if(this.spinSystem.hmbcLines != null) {
	            expLH++;
	            score+=this._hmbcScore(partial, current, keySingalAsg);
	        }
	        if(this.spinSystem.chemicalShiftsT != null && this.errorCS > 0) {
	            expLH++;
	            score+=this._chemicalShiftScore(partial, current, keySingalAsg);
	        }

	        if(expLH==0){
	            expLH=3;
	            score=3;
	        }

	        this.scores[current] = score/expLH;
	        var sumLh=0;
	        var count=0;
	        for(var i = this.scores.length-1; i >= 0; i--) {
	            if(this.scores[i] != -1) {
	                sumLh += this.scores[i];
	                count++;
	            }
	        }

	        if(sumLh<this.scores.length*this.lowerBound)
	            return -sumLh/count;
	        return sumLh/count;
	    }

	    /**
	     * This function calculates the assignment score for the chemical shift restrictions.
	     * @param partial
	     * @param current
	     * @param keySingalAsg
	     * @return
	     */
	    _chemicalShiftScore(partial, current, keySingalAsg) {

	        if(this.errorCS <= 0)
	            return 1;

	        var csSignal = this.spinSystem.chemicalShiftsE[current];
	        var widthSignal = this.spinSystem.signalsWidth[current]/2.0;

	        var score = 0;
	        var csGroup = 0;
	        var diff=0;
	        var nbGroups = 0;
	        try {
	            var assignedGroups = partial[current];
	            for(var i = assignedGroups.length-1; i >= 0; i--) {
	                csGroup = this.spinSystem.chemicalShiftsT[assignedGroups[i]];
	                if(csGroup != -9999999) {
	                    nbGroups++;
	                    diff = Math.abs(csSignal-csGroup);
	                    if(diff <= widthSignal)
	                        score+=1;
	                    else{
	                        diff = Math.abs(diff-widthSignal);
	                        score += (-0.25/this.errorCS)*diff+1;
	                    }
	                }
	            }
	            if(nbGroups==0)
	                return 1.0;
	            return score/nbGroups;
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
	        var goodness=0;
	        var cosyLine = this.spinSystem.cosyLines[keySingalAsg];
	        var count1 =0;
	        var count0 =0;
	        var size = cosyLine.length-1;

	        for(var i = partial.length-1;i >= 0; i--) {
	            try {
	                var signal2 = partial[i];
	                if(i != current && signal2.length > 0) {
	                    var key=0;
	                    //The unique key for the union of those signals
	                    try {
	                        for(var j = signal2.length-1; j >= 0; j--) {
	                            key|=1<<signal2.getInt(j);
	                        }
	                    } catch (ex) {
	                        console.log("Exception in cosy score function " + ex);
	                    }

	                    var cosyLine2 = this.spinSystem.cosyLines[key];
	                    var crossPeak = false;
	                    for(var j = size; j >= 0; j--) {
	                        if(cosyLine[j] == 6 && cosyLine2[j] != 0)
	                            crossPeak = true;
	                    }
	                    if(crossPeak)
	                        count1++;
	                    else
	                        count0++;

	                    if(this.spinSystem.cosyE[current][i] == 0 && crossPeak)
	                        goodness-=0.5;
	                    if(this.spinSystem.cosyE[current][i] == 1 && !crossPeak)
	                        goodness-=0.5;
	                    if(this.spinSystem.cosyE[current][i] == 1 && crossPeak)
	                        goodness+=1;
	                    if(this.spinSystem.cosyE[current][i] == 0 && !crossPeak)
	                        goodness+=0.5;
	                }
	            } catch (e1) {
	                console.log("Exception in cosy score function " + e);
	            }
	        }
	        return Math.exp(-Math.abs((count1+count0/2.0)-goodness)/2.0);
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
	        var sizeT = hmbcLine.length-1;
	        var sizeE = this.spinSystem.hmbcE[0].length-1;
	        var freedom = sizeT-sizeE + this.MAXERRORSHMBC;
	        var crossPeaks = 0;
	        for(var j = sizeT; j >= 0; j--) {
	            if(hmbcLine[j] == 1)
	                crossPeaks++;
	        }
	        for(var j = sizeE; j >= 0; j--) {
	            if(this.spinSystem.hmbcE[current][j] == 1)
	                crossPeaks--;
	        }

	        if(crossPeaks < freedom)
	            crossPeaks=freedom;

	        return Math.exp(-Math.abs(crossPeaks-freedom)/(sizeT+1));
	    }
	}
	module.exports = Assignment;

/***/ },
/* 3 */
/***/ function(module, exports) {

	'use strict';
	/**
	 * Created by acastillo on 9/3/16.
	 */

	class TreeSet{

	    constructor(compatator){
	        this.length = 0;
	        this.elements = [];
	        if(compatator)
	            this.compatator = compatator;
	        else
	            this.compatator = function(a, b){ return a - b };
	    }

	    size(){
	        return this.elements.length;
	    }

	    last(){
	        return this.elements[this.length-1];
	    }

	    first(){
	        return this.elements[0];
	    }

	    isEmpty(){
	        return this.size()===0;
	    }

	    pollLast(){
	        if(this.length>0){
	            this.length--;
	            return this.elements.splice(this.length, 1);
	        }
	        return null;
	    }

	    pollFirst(){
	        if(this.length>0) {
	            this.length--;
	            return this.elements.splice(0, 1);
	        }
	        return null;
	    }

	    add(element){
	        let index = this.binarySearch(element);
	        if(index < 0){
	            index = -index-1;
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
	            var mid = (low + high) >>> 1;
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

/***/ }
/******/ ])
});
;