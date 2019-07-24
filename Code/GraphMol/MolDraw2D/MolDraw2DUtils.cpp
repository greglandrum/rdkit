//
//  Copyright (C) 2016-2019 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>

#include <GraphMol/RWMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <RDGeneral/BoostEndInclude.h>
#include <limits>
#include <cmath>
#include <Numerics/Conrec.h>

namespace RDKit {
namespace MolDraw2DUtils {

namespace {
bool isAtomCandForChiralH(const RWMol &mol, const Atom *atom) {
  // conditions for needing a chiral H:
  //   - stereochem specified
  //   - in at least two rings
  if ((!mol.getRingInfo()->isInitialized() ||
       mol.getRingInfo()->numAtomRings(atom->getIdx()) > 1) &&
      (atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW ||
       atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW)) {
    return true;
  }
  return false;
}
}  // end of anonymous namespace

void prepareMolForDrawing(RWMol &mol, bool kekulize, bool addChiralHs,
                          bool wedgeBonds, bool forceCoords) {
  if (kekulize) {
    MolOps::Kekulize(mol, false);  // kekulize, but keep the aromatic flags!
  }
  if (addChiralHs) {
    std::vector<unsigned int> chiralAts;
    for (RWMol::AtomIterator atIt = mol.beginAtoms(); atIt != mol.endAtoms();
         ++atIt) {
      if (isAtomCandForChiralH(mol, *atIt)) {
        chiralAts.push_back((*atIt)->getIdx());
      }
    }
    if (chiralAts.size()) {
      bool addCoords = false;
      if (!forceCoords && mol.getNumConformers()) addCoords = true;
      MolOps::addHs(mol, false, addCoords, &chiralAts);
    }
  }
  if (forceCoords || !mol.getNumConformers()) {
    // compute 2D coordinates in a standard orientation:
    const bool canonOrient = true;
    RDDepict::compute2DCoords(mol, nullptr, canonOrient);
  }
  if (wedgeBonds) {
    WedgeMolBonds(mol, &mol.getConformer());
  }
}

void prepareAndDrawMolecule(MolDraw2D &drawer, const ROMol &mol,
                            const std::string &legend,
                            const std::vector<int> *highlight_atoms,
                            const std::vector<int> *highlight_bonds,
                            const std::map<int, DrawColour> *highlight_atom_map,
                            const std::map<int, DrawColour> *highlight_bond_map,
                            const std::map<int, double> *highlight_radii,
                            int confId) {
  RWMol cpy(mol);
  prepareMolForDrawing(cpy);
  drawer.drawMolecule(cpy, legend, highlight_atoms, highlight_bonds,
                      highlight_atom_map, highlight_bond_map, highlight_radii,
                      confId);
}

void updateDrawerParamsFromJSON(MolDraw2D &drawer, const char *json) {
  PRECONDITION(json, "no parameter string");
  updateDrawerParamsFromJSON(drawer, std::string(json));
};
#define PT_OPT_GET(opt) opts.opt = pt.get(#opt, opts.opt)

void get_colour_option(boost::property_tree::ptree *pt, const char *pnm,
                       DrawColour &colour) {
  PRECONDITION(pnm && strlen(pnm), "bad property name");
  if (pt->find(pnm) == pt->not_found()) return;

  boost::property_tree::ptree::const_iterator itm = pt->get_child(pnm).begin();
  colour.r = itm->second.get_value<float>();
  ++itm;
  colour.g = itm->second.get_value<float>();
  ++itm;
  colour.b = itm->second.get_value<float>();
  ++itm;
}

void updateDrawerParamsFromJSON(MolDraw2D &drawer, const std::string &json) {
  if (json == "") return;
  std::istringstream ss;
  ss.str(json);
  MolDrawOptions &opts = drawer.drawOptions();
  boost::property_tree::ptree pt;
  boost::property_tree::read_json(ss, pt);
  PT_OPT_GET(atomLabelDeuteriumTritium);
  PT_OPT_GET(dummiesAreAttachments);
  PT_OPT_GET(circleAtoms);
  PT_OPT_GET(continuousHighlight);
  PT_OPT_GET(flagCloseContactsDist);
  PT_OPT_GET(includeAtomTags);
  PT_OPT_GET(clearBackground);
  PT_OPT_GET(legendFontSize);
  PT_OPT_GET(multipleBondOffset);
  PT_OPT_GET(padding);
  PT_OPT_GET(additionalAtomLabelPadding);
  get_colour_option(&pt, "highlightColour", opts.highlightColour);
  get_colour_option(&pt, "backgroundColour", opts.backgroundColour);
  get_colour_option(&pt, "legendColour", opts.legendColour);
  if (pt.find("atomLabels") != pt.not_found()) {
    BOOST_FOREACH (boost::property_tree::ptree::value_type const &item,
                   pt.get_child("atomLabels")) {
      opts.atomLabels[boost::lexical_cast<int>(item.first)] =
          item.second.get_value<std::string>();
    }
  }
}

void contourAndDrawGrid(MolDraw2D &drawer, const double *grid,
                        const std::vector<double> &xcoords,
                        const std::vector<double> &ycoords, size_t nContours,
                        std::vector<double> &levels, bool dashNegative) {
  PRECONDITION(grid, "no data");
  PRECONDITION(nContours > 0, "no contours");

  size_t nX = xcoords.size();
  size_t nY = ycoords.size();
  double dX = (xcoords.back() - xcoords.front()) / nX;
  double dY = (ycoords.back() - ycoords.front()) / nY;
  if (!levels.size()) {
    levels.resize(nContours);
    double minV = std::numeric_limits<double>::max();
    double maxV = std::numeric_limits<double>::min();
    for (size_t i = 0; i < nX; ++i) {
      for (size_t j = 0; j < nY; ++j) {
        minV = std::min(minV, grid[i * nY + j]);
        maxV = std::max(maxV, grid[i * nY + j]);
      }
    }
    for (size_t i = 0; i < nContours; ++i) {
      levels[i] = minV + i * (maxV - minV) / (nContours - 1);
    }
  }
  std::vector<conrec::ConrecSegment> segs;
  conrec::Contour(grid, 0, nX - 1, 0, nY - 1, xcoords.data(), ycoords.data(),
                  nContours, levels.data(), segs);
  const auto olw = drawer.lineWidth();
  const auto odash = drawer.dash();
  const auto ocolor = drawer.colour();
  drawer.setLineWidth(1.0);
  drawer.setColour(DrawColour(0.5, 0.5, 0.5));
  static DashPattern negDash = {2, 6};
  static DashPattern posDash;
  std::cerr << "  DRAWING " << segs.size() << " SEGMENTS" << std::endl;
  for (const auto &seg : segs) {
    if (dashNegative && seg.isoVal < 0) {
      drawer.setDash(negDash);
    } else {
      drawer.setDash(posDash);
    }
    drawer.drawLine(seg.p1, seg.p2);
  }
  drawer.setDash(odash);
  drawer.setLineWidth(olw);
  drawer.setColour(ocolor);
};

void contourAndDrawGaussians(MolDraw2D &drawer,
                             const std::vector<Point2D> &locs,
                             const std::vector<double> &weights,
                             const std::vector<double> &widths,
                             size_t nContours, std::vector<double> &levels,
                             double resolution, bool dashNegative) {
  PRECONDITION(locs.size() == weights.size(), "size mismatch");
  PRECONDITION(locs.size() == widths.size(), "size mismatch");

  // start by setting up the grid
  size_t nx = (size_t)ceil(drawer.range().x / resolution) + 1;
  size_t ny = (size_t)ceil(drawer.range().y / resolution) + 1;
  std::vector<double> xcoords(nx);
  for (size_t i = 0; i < nx; ++i) {
    xcoords[i] = drawer.minPt().x + i * resolution;
  }
  std::vector<double> ycoords(ny);
  for (size_t i = 0; i < ny; ++i) {
    ycoords[i] = drawer.minPt().y + i * resolution;
  }
  std::unique_ptr<double[]> grid(new double[nx * ny]);

  // populate the grid from the gaussians:
  for (size_t ix = 0; ix < nx; ++ix) {
    auto px = drawer.minPt().x + ix * resolution;
    for (size_t iy = 0; iy < ny; ++iy) {
      auto py = drawer.minPt().y + iy * resolution;
      Point2D pt(px, py);
      double accum = 0.0;
      for (size_t ig = 0; ig < locs.size(); ++ig) {
        auto d2 = (pt - locs[ig]).lengthSq();
        auto contrib = weights[ig] / widths[ig] *
                       exp(-0.5 * d2 / (widths[ig] * widths[ig]));
        accum += contrib;
      }
      grid[ix * ny + iy] = accum / (2 * M_PI);
    }
  }

  // and render it:
  contourAndDrawGrid(drawer, grid.get(), xcoords, ycoords, nContours, levels,
                     dashNegative);
};
}  // namespace MolDraw2DUtils
}  // namespace RDKit
