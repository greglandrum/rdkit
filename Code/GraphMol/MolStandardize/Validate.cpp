#include "Validate.h"
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/MolStandardize/FragmentCatalog/FragmentCatalogParams.h>
#include <GraphMol/MolStandardize/FragmentCatalog/FragmentRemover.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <iostream>
#include <vector>
#include <string>

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

using namespace std;
using namespace RDKit;

namespace RDKit {
class RWMol;
class ROMol;

namespace MolStandardize {

std::vector<ValidationErrorInfo> RDKitValidation::validate(
    const ROMol &mol, bool reportAllFailures) const {
  ROMol molCopy = mol;
  std::vector<ValidationErrorInfo> errors;

  unsigned int na = mol.getNumAtoms();

  if (!na) {
    errors.push_back(
        ValidationErrorInfo("ERROR: [NoAtomValidation] Molecule has no atoms"));
  }

  // loop over atoms
  for (size_t i = 0; i < na; ++i) {
    if (!reportAllFailures) {
      if (errors.size() >= 1) {
        break;
      }
    }
    Atom *atom = molCopy.getAtomWithIdx(i);
    try {
      int explicitValence = atom->calcExplicitValence();
    } catch (const MolSanitizeException &e) {
      // std::cout << e.message() << std::endl;
      errors.push_back(ValidationErrorInfo("INFO: [ValenceValidation] " +
                                           std::string(e.message())));
    }
  }
  return errors;
}

void NoAtomValidation::run(
    const ROMol &mol, bool reportAllFailures,
    std::vector<ValidationErrorInfo> &errors) const {
  unsigned int na = mol.getNumAtoms();

  if (!na) {
    errors.push_back(
        ValidationErrorInfo("ERROR: [NoAtomValidation] Molecule has no atoms"));
  }
}

void FragmentValidation::run(
    const ROMol &mol, bool reportAllFailures,
    std::vector<ValidationErrorInfo> &errors) const {
  std::string rdbase = getenv("RDBASE");
  std::string fgrpFile = rdbase +
                         "/Code/GraphMol/MolStandardize/FragmentCatalog/"
                         "data/fragmentPatterns.txt";
  std::shared_ptr<FragmentCatalogParams> fparams(
      new FragmentCatalogParams(fgrpFile));
  FragmentCatalog fcat(fparams.get());

  const std::vector<std::shared_ptr<ROMol>> &fgrps = fparams->getFuncGroups();
  INT_VECT mapping;
  VECT_INT_VECT atom_mapping;
  std::vector<ROMOL_SPTR> frags =
      MolOps::getMolFrags(mol, true, &mapping, &atom_mapping);

  for (auto &fgrp : fgrps) {
    std::string fname;
    fgrp->getProp(common_properties::_Name, fname);
    std::vector<RDKit::MatchVectType> res;
    unsigned int matches = SubstructMatch(mol, *fgrp, res);
    //		std::cout << fname << " matches " << matches << std::endl;
    if (matches != 0 && frags.size() != 0) {
      VECT_INT_VECT substructmap;  // store idxs of frag from substructmatch
      for (const auto &match : res) {
        std::vector<int> vec;
        for (const auto &pair : match) {
          vec.push_back(pair.second);
        }
        substructmap.push_back(vec);
      }

      // to stop the same fragment being reported many times if present
      // multiple times in molecule
      bool fpresent = false;

      for (auto &molfragidx : atom_mapping) {
        std::sort(molfragidx.begin(), molfragidx.end());
        for (auto &substructidx : substructmap) {
          std::sort(substructidx.begin(), substructidx.end());
          //					// help to debug...
          //					std::cout << "molfragidx: "  <<
          // std::endl; 					for (const auto &i : molfragidx) {
          // std::cout << i; }; 					std::cout << std::endl; 					std::cout <<
          //"substructidx: "  << std::endl;
          // for (const auto &i : substructidx) { std::cout << i; }; 					std::cout <<
          //std::endl;
          //					//
          if (molfragidx == substructidx & !fpresent) {
            std::string msg = fname + " is present";
            errors.push_back(
                ValidationErrorInfo("INFO: [FragmentValidation] " + msg));
            fpresent = true;
          }
        }
      }
    }
  }
}

void NeutralValidation::run(
    const ROMol &mol, bool reportAllFailures,
    std::vector<ValidationErrorInfo> &errors) const {
  int charge = RDKit::MolOps::getFormalCharge(mol);
  if (charge != 0) {
    std::string charge_str;
    if (charge > 0) {
      charge_str = "+" + std::to_string(charge);
    } else {
      charge_str = std::to_string(charge);
    }
    std::string msg = "Not an overall neutral system (" + charge_str + ')';
    //		std::cout << msg << std::endl;
    errors.push_back(ValidationErrorInfo("INFO: [NeutralValidation] " + msg));
  }
}

void IsotopeValidation::run(
    const ROMol &mol, bool reportAllFailures,
    std::vector<ValidationErrorInfo> &errors) const {
  unsigned int na = mol.getNumAtoms();
  std::set<string> isotopes;

  // loop over atoms
  for (size_t i = 0; i < na; ++i) {
    if (!reportAllFailures) {
      if (errors.size() >= 1) {
        break;
      }
    }
    const Atom *atom = mol.getAtomWithIdx(i);
    unsigned int isotope = atom->getIsotope();
    if (isotope != 0) {
      std::string symbol = atom->getSymbol();
      isotopes.insert(std::to_string(isotope) + symbol);
    }
  }

  for (auto &isotope : isotopes) {
    errors.push_back(ValidationErrorInfo(
        "INFO: [IsotopeValidation] Molecule contains isotope " + isotope));
  }
}

// constructor
MolVSValidation::MolVSValidation() {
	std::vector<MolVSValidations*> validations = 
	{new NoAtomValidation(), 
	new FragmentValidation(), 
	new NeutralValidation(), 
	new IsotopeValidation()};
	this->d_validations = validations;
}

// overloaded constructor
MolVSValidation::MolVSValidation(const std::vector<MolVSValidations*> validations) {
	this->d_validations = validations;
}

// copy constructor
MolVSValidation::MolVSValidation(const MolVSValidation &other) {
	d_validations = other.d_validations;
}

MolVSValidation::~MolVSValidation() {
	for (const auto p_method : this->d_validations) {
		delete p_method;
	}
};

std::vector<ValidationErrorInfo> MolVSValidation::validate(
    const ROMol &mol, bool reportAllFailures) const {
  std::vector<ValidationErrorInfo> errors;

	for (const auto method : this->d_validations) {
		method->run(mol, reportAllFailures, errors); 
	}

  return errors;
}

std::vector<ValidationErrorInfo> AllowedAtomsValidation::validate(
		const ROMol &mol, bool reportAllFailures) const {
  std::vector<ValidationErrorInfo> errors;
  unsigned int na = mol.getNumAtoms();

  for (size_t i = 0; i < na; ++i) {
    if (!reportAllFailures) {
      if (errors.size() >= 1) {
        break;
      }
    }
    const Atom *qatom = mol.getAtomWithIdx(i);
    bool match = false;
    // checks to see qatom matches one of list of allowedAtoms
    for (const auto &allowedAtom : this->d_allowedList) {
      if (allowedAtom->Match(qatom)) {
        match = true;
        break;
      }
    }
    // if no match, append to list of errors.
    if (!match) {
      std::string symbol = qatom->getSymbol();
      errors.push_back(
          ValidationErrorInfo("INFO: [AllowedAtomsValidation] Atom " + symbol +
                              " is not in allowedAtoms list"));
    }
  }
  return errors;
}

std::vector<ValidationErrorInfo> DisallowedAtomsValidation::validate(
    const ROMol &mol, bool reportAllFailures) const {
  std::vector<ValidationErrorInfo> errors;
  unsigned int na = mol.getNumAtoms();

  for (size_t i = 0; i < na; ++i) {
    if (!reportAllFailures) {
      if (errors.size() >= 1) {
        break;
      }
    }
    const Atom *qatom = mol.getAtomWithIdx(i);
    bool match = false;
    // checks to see qatom matches one of list of allowedAtoms
    for (const auto &disallowedAtom : this->d_disallowedList) {
      if (disallowedAtom->Match(qatom)) {
        match = true;
      }
    }
    // if no match, append to list of errors.
    if (match) {
      std::string symbol = qatom->getSymbol();
      errors.push_back(
          ValidationErrorInfo("INFO: [DisallowedAtomsValidation] Atom " +
                              symbol + " is in disallowedAtoms list"));
    }
  }
  return errors;
}

}  // namespace MolStandardize
}  // namespace RDKit