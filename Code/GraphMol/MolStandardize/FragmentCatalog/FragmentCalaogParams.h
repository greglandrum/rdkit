#ifndef __RD_FRAGMENT_CATALOG_PARAMS_H__
#define __RD_FRAGMENT_CATALOG_PARAMS_H__

#include <Catalogs/CatalogParams.h>
#include "FragmentCalaogUtils.h"
#include <GraphMol/RDKitBase.h>
#include <string>
#include <vector>
#include <iostream>

namespace RDKit {
class ROMol;

namespace MolStandardize {
class FragmentCatalogParams : public RDCatalog::CatalogParams {

	public:
		FragmentCatalogParams() {
			d_typeStr = "Fragment Catalog Parameters";
			d_funcGroups.clear();
		}

		FragmentCatalogParams(const std::string &fgroupFile);
		// copy constructor
		FragmentCatalogParams(const FragmentCatalogParams &other);

		~FragmentCatalogParams();

		unsigned int getNumFuncGroups() const {
		       return static_cast<unsigned int>(d_funcGroups.size()); }

		const std::vector<std::shared_ptr<ROMol>> &getFuncGroups() const;

		const ROMol *getFuncGroup(unsigned int fid) const;
		
		void toStream(std::ostream &) const;
		std::string Serialize() const;
		void initFromStream(std::istream &ss);
		void initFromString(const std::string &text);

	private: 
		std::vector<std::shared_ptr<ROMol>> d_funcGroups;

}; // class FragmentCatalogParams

} // namespace MolStandardize
} // namespace RDKit

#endif