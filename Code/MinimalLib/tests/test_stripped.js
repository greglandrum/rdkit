//
//
//  Copyright (C) 2019 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

const assert = require('assert');
const {
    performance
  } = require('perf_hooks');
var initRDKitModule = require("../demo/RDKit_minimal.js");
var RDKitModule;
const fs       = require('fs');
const readline = require('readline');

// the goal here isn't to be comprehensive (the RDKit has tests for that),
// just to make sure that the wrappers are working as expected
function test_basics() {
    var bmol = RDKitModule.get_mol("c1ccccc","");
    assert.equal(bmol.is_valid(),0);
    
    var mol = RDKitModule.get_mol("c1ccccc1O","");
    assert.equal(mol.is_valid(),1);
    var qmol = RDKitModule.get_mol("c1ccccc1","");
    assert.equal(qmol.is_valid(),1);

    var match = mol.get_substruct_match(qmol);

    var pmatch = JSON.parse(match);
    assert.equal(pmatch.atoms.length,6);

    var qmol2 = RDKitModule.get_qmol("ccc");
    assert.equal(qmol2.is_valid(),1);

    var match2 = mol.get_substruct_match(qmol2);

    var pmatch2 = JSON.parse(match2);
    assert.equal(pmatch2.atoms.length,3);

    // var fp = mol.get_morgan_fp_as_uint8array(2,1024);
    // console.log(fp);

    // assert.equal(mol.get_to msmiles(),"Oc1ccccc1");
    // assert.equal(mol.get_inchi(),"InChI=1S/C6H6O/c7-6-4-2-1-3-5-6/h1-5,7H");
    // assert.equal(RDKitModule.get_inchikey_for_inchi(mol.get_inchi()),"ISWSIDIOOBJBQZ-UHFFFAOYSA-N");

    // var mb = mol.get_molblock();
    // assert(mb.search("M  END")>0);
    // var mol2 = RDKitModule.get_mol(mb);
    // assert.equal(mol2.is_valid(),1);
    // assert.equal(mol2.get_smiles(),"Oc1ccccc1");
    
    // var mjson = mol.get_json();
    // assert(mjson.search("commonchem")>0);
    // var mol3 = RDKitModule.get_mol(mjson);
    // assert.equal(mol3.is_valid(),1);
    // assert.equal(mol3.get_smiles(),"Oc1ccccc1");
    
    // var descrs = JSON.parse(mol.get_descriptors());
    // assert.equal(descrs.NumAromaticRings,1);
    // assert.equal(descrs.NumRings,1);
    // assert.equal(descrs.amw,94.11299);

    // var checkStringBinaryFpIdentity = (stringFp, binaryFp) => {
    //     assert.equal(binaryFp.length, stringFp.length / 8);
    //     for (var i = 0, c = 0; i < binaryFp.length; ++i) {
    //         var byte = 0;
    //         for (var j = 0; j < 8; ++j, ++c) {
    //             if (stringFp[c] === "1") {
    //                 byte |= (1 << j);
    //             }
    //         }
    //         assert.equal(byte, binaryFp[i]);
    //     }
    // };

    // var fp1 = mol.get_morgan_fp();
    // assert.equal(fp1.length,2048);
    // assert.equal((fp1.match(/1/g)||[]).length,11);
    // var fp1Uint8Array = mol.get_morgan_fp_as_uint8array();
    // checkStringBinaryFpIdentity(fp1, fp1Uint8Array);
    // var fp2 = mol.get_morgan_fp(0,512);
    // assert.equal(fp2.length,512);
    // assert.equal((fp2.match(/1/g)||[]).length,3);
    // var fp2Uint8Array = mol.get_morgan_fp_as_uint8array(0, 512);
    // checkStringBinaryFpIdentity(fp2, fp2Uint8Array);
    
    // var svg = mol.get_svg();
    // assert(svg.search("svg")>0);

    // var qmol = RDKitModule.get_qmol("Oc(c)c");
    // assert.equal(qmol.is_valid(),1);
    // var match = mol.get_substruct_match(qmol);
    // var pmatch = JSON.parse(match);
    // assert.equal(pmatch.atoms.length,4);
    // assert.equal(pmatch.atoms[0],6);
    // var svg2 = mol.get_svg_with_highlights(match);
    // assert(svg2.search("svg")>0);
    // assert(svg.search("#FF7F7F")<0);
    // assert(svg2.search("#FF7F7F")>0);
}

initRDKitModule().then(function(instance) {
    var done = {};
    const waitAllTestsFinished = () => {
        const poll = resolve => {
            if (Object.values(done).every(v => v)) {
                resolve();
            } else {
                setTimeout(() => poll(resolve), 100);
            }
        }
        return new Promise(poll);
    }
    RDKitModule = instance;
    console.log(RDKitModule.version());
    test_basics();
    waitAllTestsFinished().then(() =>
        console.log("Tests finished successfully")
    );
});
