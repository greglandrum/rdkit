#include "TransformCatalogParams.h"
#include "TransformCatalogUtils.h"
#include <GraphMol/RDKitBase.h>
#include <sstream>

namespace RDKit {
namespace MolStandardize {

TransformCatalogParams::TransformCatalogParams(const std::string &transformFile) {
	d_transformations.clear();
	d_transformations = readTransformations(transformFile);
}

TransformCatalogParams::TransformCatalogParams(const TransformCatalogParams &other) {
	d_typeStr = other.d_typeStr;
	d_transformations.clear();

	const std::vector<std::shared_ptr<ChemicalReaction>> &otransforms = 
					other.getTransformations();
	for (auto &transi : otransforms) {
		std::shared_ptr<ChemicalReaction> transform( new ChemicalReaction(*transi) );
		d_transformations.push_back(transform);
	}
}

TransformCatalogParams::~TransformCatalogParams() {}

const std::vector<std::shared_ptr<ChemicalReaction>> &TransformCatalogParams::getTransformations() const {
  return d_transformations;
}

const ChemicalReaction *TransformCatalogParams::getTransformation(unsigned int fid) const {
  URANGE_CHECK(fid, d_transformations.size());
  // return d_transformations[fid];
  return d_transformations[fid].get();
}

void TransformCatalogParams::toStream(std::ostream &ss) const {
	ss << d_transformations.size() << "\n";
}

std::string TransformCatalogParams::Serialize() const {
	std::stringstream ss;
	toStream(ss);
	return ss.str();
}

void TransformCatalogParams::initFromStream(std::istream &ss) {

}

void TransformCatalogParams::initFromString(const std::string &text) {

}

} // namespace MolStandardize
} // namespace RDKit
