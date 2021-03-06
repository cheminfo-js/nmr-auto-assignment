/**
 * Created by acastillo on 5/7/16.
 */
'use strict';

const SD = require('spectra-data');
const FS = require('fs');
const OCLE = require("openchemlib-extended");
const autoassigner = require('../src/index');
const Predictor = require("nmr-predictor");


function createSpectraData(filename, label, data) {
    var spectrum = SD.NMR.fromJcamp(
        FS.readFileSync(__dirname + filename).toString()
    );
    return spectrum;
};

function loadFile(filename){
    return FS.readFileSync(__dirname + filename).toString();
}


describe('Auto-assignment ethylvinylether', function () {
    var molfile = loadFile("/data/ethylvinylether/structure.mol");
    var molecule = OCLE.Molecule.fromMolfile(molfile);
    molecule.addImplicitHydrogens();
    var nH = molecule.getMolecularFormula().formula.replace(/.*H([0-9]+).*/,"$1")*1;
    var diaIDs = molecule.getGroupedDiastereotopicAtomIDs();
    for (var j=0; j<diaIDs.length; j++) {
        diaIDs[j].nbEquivalent=diaIDs[j].atoms.length;
    }
    diaIDs.sort(function(a,b) {
        if (a.atomLabel==b.atomLabel) {
            return b.nbEquivalent-a.nbEquivalent;
        }
        return a.atomLabel<b.atomLabel?1:-1;
    });
    //var diaID = molecule.getIDCode();
    const predictor = new Predictor(JSON.parse(loadFile("/../src/h1_database.json")));

    var spectrum = createSpectraData("/data/ethylvinylether/1h.jdx");
    var peakPicking = spectrum.nmrPeakDetection({
        "nH": nH,
        realTop: true,
        thresholdFactor: 1,
        clean: true,
        compile: true,
        idPrefix: "1H",
        format:"new"
    });
    //console.log(peakPicking);
    //console.log("PP", spectrum.sd.spectra[0].data);
    it('Known patterns for ethylvinylether', function () {
        var result = autoassigner({molecule:molecule, diaIDs:diaIDs,
                spectra:{h1PeakList:peakPicking, solvent:spectrum.getParamString(".SOLVENT NAME", "unknown")}},
            {predictor:predictor}
        );
        //console.log(JSON.stringify(result));
    });
});