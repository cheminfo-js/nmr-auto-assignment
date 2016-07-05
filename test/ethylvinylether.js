/**
 * Created by acastillo on 5/7/16.
 */
'use strict';

var SD = require('spectra-data');
var FS = require('fs');
var OCLE = require("openchemlib-extended");

function createSpectraData(filename, label, data) {
    var spectrum = SD.NMR.fromJcamp(
        FS.readFileSync(__dirname + filename).toString()
    );
    return spectrum;
};

function loadMolfile(filename){
    return FS.readFileSync(__dirname + filename).toString();

}

describe('spectra-data test peak picking', function () {
        var molfile = loadMolfile("/data/ethylvinylether/structure.mol");
    var molecule = OCLE.Molecule.fromMolfile(molfile);
    molecule.addImplicitHydrogens();
    var nH = molecule.getNumberOfAtoms('H');
    var diaIDs = molecule.getDiastereotopicAtomIDs();
    for (var j=0; j<diaIDs.length; j++) {
        diaIDs[j].nbEquivalent=diaIDs[j].atoms.length;
    }
    diaIDs.sort(function(a,b) {
        if (a.element==b.element) {
            return b.nbEquivalent-a.nbEquivalent;
        }
        return a.element<b.element?1:-1;
    });
    var diaID = molecule.toIDCode();
    console.log(diaIDs);

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
    it('Known patterns for ethylvinylether', function () {

    });
});