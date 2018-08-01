#include <RDBoost/Wrap.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolStandardize/FragmentCatalog/FragmentRemover.h>

namespace python = boost::python;
using namespace RDKit;

namespace {

ROMol* removeHelper(MolStandardize::FragmentRemover &self, const ROMol &mol){
	return self.remove(mol);
}

ROMol* chooseHelper(MolStandardize::LargestFragmentChooser &self, const ROMol &mol){
	return self.choose(mol);
}

} // namespace

BOOST_PYTHON_MODULE(Fragment) {
	python::scope().attr("__doc__") = 
					"Module containing tools for dealing with molecules with more than \
					covalently bonded unit";

	std::string docString = "";

	python::class_<MolStandardize::FragmentRemover, boost::noncopyable>(
					"FragmentRemover", python::init<>() )
					.def(python::init<std::string, bool>())
					.def("remove", &removeHelper,
							(python::arg("self"), python::arg("mol")),
							"",
							python::return_value_policy<python::manage_new_object>())
					;					

	python::class_<MolStandardize::LargestFragmentChooser, boost::noncopyable>(
					"LargestFragmentChooser", python::init<bool>(
									(python::arg("preferOrganic") = false) ))
					.def("choose", &chooseHelper,
							(python::arg("self"), python::arg("mol")),
							"",
							python::return_value_policy<python::manage_new_object>())
					;

//					.def(python 



}
