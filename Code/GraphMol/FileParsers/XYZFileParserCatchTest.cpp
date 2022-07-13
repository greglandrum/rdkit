#include "RDGeneral/test.h"
#include "catch.hpp"
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <RDGeneral/FileParseException.h>

using namespace RDKit;

TEST_CASE("xyz file parser") {
    SECTION("basics") {
        std::string rdbase = getenv("RDBASE");
        std::string fName = rdbase +
                            "/Code/GraphMol/FileParsers/test_data/"
                            "acetate.xyz";
        std::unique_ptr<RWMol> mol(XYZFileToMol(fName, 0));
        REQUIRE(mol);
        
        REQUIRE(mol->getNumAtoms() == 7);
        
        size_t ind = 0;
        for (auto a : {6u, 6u, 1u, 1u, 1u, 8u, 8u}) {
            CHECK(mol->getAtomWithIdx(ind)->getAtomicNum() == a);
            ind++;
        }
        
        RDGeom::POINT3D_VECT positions = {
            RDGeom::Point3D{-4.71686, 0.89919, 0.05714},
            RDGeom::Point3D{-3.24898, 0.98400, -0.22830},
            RDGeom::Point3D{-5.04167, 1.74384, 0.67862},
            RDGeom::Point3D{-5.01710, -0.02205, 0.56344},
            RDGeom::Point3D{-5.21076, 0.96874, -0.91208},
            RDGeom::Point3D{-2.65909, 2.05702, -0.34025},
            RDGeom::Point3D{-2.63413, -0.18702, -0.48679}
        };
        
        CHECK(mol->hasProp("_FileComments"));
        CHECK(mol->getProp<std::string>("_FileComments") == "charge=-1=");
        
        auto conf = &mol->getConformer();
        REQUIRE(conf);
        REQUIRE(conf->getNumAtoms() == 7);
        
        ind = 0;
        for (auto p : positions) {
            CHECK(conf->getAtomPos(ind).x == p.x);
            CHECK(conf->getAtomPos(ind).y == p.y);
            CHECK(conf->getAtomPos(ind).z == p.z);
            ind++;
        }
    }
    SECTION("basics w/o comment and w/ nonstandard spacing") {
        std::string rdbase = getenv("RDBASE");
        std::string fName = rdbase +
                            "/Code/GraphMol/FileParsers/test_data/"
                            "ethane.xyz";
        std::unique_ptr<RWMol> mol(XYZFileToMol(fName, 0));
        REQUIRE(mol);
        
        REQUIRE(mol->getNumAtoms() == 8);
        
        size_t ind = 0;
        for (auto a : {6u, 6u, 1u, 1u, 1u, 1u, 1u, 1u}) {
            CHECK(mol->getAtomWithIdx(ind)->getAtomicNum() == a);
            ind++;
        }
        
        RDGeom::POINT3D_VECT positions = {
            RDGeom::Point3D{-4.58735, 0.92696, 0.00000},
            RDGeom::Point3D{-3.11050, 0.92696, 0.00000},
            RDGeom::Point3D{-4.93786, 1.78883, 0.58064},
            RDGeom::Point3D{-4.93786, -0.00682, 0.45608},
            RDGeom::Point3D{-4.93786, 0.99888, -1.03672},
            RDGeom::Point3D{-2.75999, 0.85505, 1.03672},
            RDGeom::Point3D{-2.75998, 1.86075, -0.45608},
            RDGeom::Point3D{-2.75998, 0.06509, -0.58064}
        };
        
        CHECK(!mol->hasProp("_FileComments"));
        
        auto conf = &mol->getConformer();
        REQUIRE(conf);
        REQUIRE(conf->getNumAtoms() == 8);
        
        ind = 0;
        for (auto p : positions) {
            CHECK(conf->getAtomPos(ind).x == p.x);
            CHECK(conf->getAtomPos(ind).y == p.y);
            CHECK(conf->getAtomPos(ind).z == p.z);
            ind++;
        }
    }
}
