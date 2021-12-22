//
//  Copyright (C) 2014-2021 David Cosgrove and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Original author: David Cosgrove (AstraZeneca)
// 27th May 2014
//
// This class makes a 2D drawing of an RDKit molecule.
// It draws heavily on $RDBASE/GraphMol/MolDrawing/MolDrawing.h.
// One purpose of this is to make it easier to overlay annotations on top of
// the molecule drawing, which is difficult to do from the output of
// MolDrawing.h
// The class design philosophy echoes a standard one:
// a virtual base class defines the interface and does all
// the heavy lifting and concrete derived classes implement
// library-specific drawing code such as drawing lines, writing strings
// etc.

#include <RDGeneral/export.h>
#ifndef RDKITMOLDRAW2D_H
#define RDKITMOLDRAW2D_H

#include <vector>

#include <Geometry/point.h>
#include <Geometry/Transform2D.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/MolDraw2D/MolDraw2DHelpers.h>

// ****************************************************************************
using RDGeom::Point2D;

namespace RDKit {

class DrawMol;
class DrawText;

//! MolDraw2D is the base class for doing 2D renderings of molecules
class RDKIT_MOLDRAW2D_EXPORT MolDraw2D {
 public:
  //! constructor for a particular size
  /*!
    \param width       : width (in pixels) of the rendering
    \param height      : height (in pixels) of the rendering
    \param panelWidth  : (optional) width (in pixels) of a single panel
    \param panelHeight : (optional) height (in pixels) of a single panel

    The \c panelWidth and \c panelHeight arguments are used to provide the
    sizes of the panels individual molecules are drawn in when
    \c drawMolecules() is called.
  */
  MolDraw2D(int width, int height, int panelWidth, int panelHeight);
  virtual ~MolDraw2D();

  //! \name Methods that must be provided by child classes
  //@{
 private:
  virtual void initDrawing() = 0;
  virtual void initTextDrawer(bool noFreetype) = 0;

 public:
  //! clears the contents of the drawing
  virtual void clearDrawing() = 0;
  //! draws a line from \c cds1 to \c cds2 using the current drawing style
  /// in atom coords.
  virtual void drawLine(const Point2D &cds1, const Point2D &cds2) = 0;
  //! draw a polygon.  Note that if fillPolys() returns false, it
  //! doesn't close the path.  If you want it to in that case, you
  //! do it explicitly yourself.
  virtual void drawPolygon(const std::vector<Point2D> &cds) = 0;
  //@}

  //! draw a single molecule
  /*!
    \param mol             : the molecule to draw
    \param legend          : the legend (to be drawn under the molecule)
    \param highlight_atoms : (optional) vector of atom ids to highlight
    \param highlight_atoms : (optional) vector of bond ids to highlight
    \param highlight_atom_map   : (optional) map from atomId -> DrawColour
    providing the highlight colors. If not provided the default highlight colour
    from \c drawOptions() will be used.
    \param highlight_bond_map   : (optional) map from bondId -> DrawColour
    providing the highlight colors. If not provided the default highlight colour
    from \c drawOptions() will be used.
    \param highlight_radii : (optional) map from atomId -> radius (in molecule
    coordinates) for the radii of atomic highlights. If not provided the default
    value from \c drawOptions() will be used.
    \param confId          : (optional) conformer ID to be used for atomic
    coordinates

  */
  virtual void drawMolecule(
      const ROMol &mol, const std::string &legend,
      const std::vector<int> *highlight_atoms,
      const std::vector<int> *highlight_bonds,
      const std::map<int, DrawColour> *highlight_atom_map = nullptr,
      const std::map<int, DrawColour> *highlight_bond_map = nullptr,
      const std::map<int, double> *highlight_radii = nullptr, int confId = -1);

  //! \overload
  virtual void drawMolecule(
      const ROMol &mol, const std::vector<int> *highlight_atoms = nullptr,
      const std::map<int, DrawColour> *highlight_map = nullptr,
      const std::map<int, double> *highlight_radii = nullptr, int confId = -1);

  //! \overload
  virtual void drawMolecule(
      const ROMol &mol, const std::string &legend,
      const std::vector<int> *highlight_atoms = nullptr,
      const std::map<int, DrawColour> *highlight_map = nullptr,
      const std::map<int, double> *highlight_radii = nullptr, int confId = -1);

  //! \overload
  virtual void drawMolecule(
      const ROMol &mol, const std::vector<int> *highlight_atoms,
      const std::vector<int> *highlight_bonds,
      const std::map<int, DrawColour> *highlight_atom_map = nullptr,
      const std::map<int, DrawColour> *highlight_bond_map = nullptr,
      const std::map<int, double> *highlight_radii = nullptr, int confId = -1);

  //! draw molecule with multiple colours allowed per atom.
  /*!
    \param mol             : the molecule to draw
    \param legend          : the legend (to be drawn under the molecule)
    \param highlight_atom_map   : map from atomId -> DrawColours
    providing the highlight colours.
    \param highlight_bond_map   : map from bondId -> DrawColours
    providing the highlight colours.
    \param highlight_radii : map from atomId -> radius (in molecule
    coordinates) for the radii of atomic highlights. If not provided for an
    index, the default value from \c drawOptions() will be used.
    \param confId          : (optional) conformer ID to be used for atomic
    coordinates
  */
  virtual void drawMoleculeWithHighlights(
      const ROMol &mol, const std::string &legend,
      const std::map<int, std::vector<DrawColour>> &highlight_atom_map,
      const std::map<int, std::vector<DrawColour>> &highlight_bond_map,
      const std::map<int, double> &highlight_radii,
      const std::map<int, int> &highlight_linewidth_multipliers,
      int confId = -1);

  //! draw multiple molecules in a grid
  /*!
    \param mols             : the molecules to draw
    \param legends          : (optional) the legends (to be drawn under the
    molecules)
    \param highlight_atoms  : (optional) vectors of atom ids to highlight
    \param highlight_atoms  : (optional) vectors of bond ids to highlight
    \param highlight_atom_map   : (optional) maps from atomId -> DrawColour
    providing the highlight colors. If not provided the default highlight colour
    from \c drawOptions() will be used.
    \param highlight_bond_map   : (optional) maps from bondId -> DrawColour
    providing the highlight colors. If not provided the default highlight colour
    from \c drawOptions() will be used.
    \param highlight_radii  : (optional) maps from atomId -> radius (in molecule
    coordinates) for the radii of atomic highlights. If not provided the default
    value from \c drawOptions() will be used.
    \param confId           : (optional) conformer IDs to be used for atomic
    coordinates

    The \c panelWidth and \c panelHeight values will be used to determine the
    number of rows and columns to be drawn. Theres not a lot of error checking
    here, so if you provide too many molecules for the number of panes things
    are likely to get screwed up.
    If the number of rows or columns ends up being <= 1, molecules will be
    being drawn in a single row/column.
  */
  virtual void drawMolecules(
      const std::vector<ROMol *> &mols,
      const std::vector<std::string> *legends = nullptr,
      const std::vector<std::vector<int>> *highlight_atoms = nullptr,
      const std::vector<std::vector<int>> *highlight_bonds = nullptr,
      const std::vector<std::map<int, DrawColour>> *highlight_atom_maps =
          nullptr,
      const std::vector<std::map<int, DrawColour>> *highlight_bond_maps =
          nullptr,
      const std::vector<std::map<int, double>> *highlight_radii = nullptr,
      const std::vector<int> *confIds = nullptr);

  //! draw a ChemicalReaction
  /*!
    \param rxn                 : the reaction to draw
    \param highlightByReactant : (optional) if this is set, atoms and bonds will
    be highlighted based on which reactant they come from. Atom map numbers
    will not be shown.
    \param highlightColorsReactants : (optional) provide a vector of colors for
    the
    reactant highlighting.
    \param confIds   : (optional) vector of confIds to use for rendering. These
    are numbered by reactants, then agents, then products.
  */
  virtual void drawReaction(
      const ChemicalReaction &rxn, bool highlightByReactant = false,
      const std::vector<DrawColour> *highlightColorsReactants = nullptr,
      const std::vector<int> *confIds = nullptr);

  //! \name Transformations
  //@{
  // transform a set of coords in the molecule's coordinate system
  // to drawing system coordinates and vice versa. Note that the coordinates
  // have
  // the origin in the top left corner, which is how Qt and Cairo have it, no
  // doubt a holdover from X Windows. This means that a higher y value will be
  // nearer the bottom of the screen. This doesn't really matter except when
  // doing text superscripts and subscripts.

  //! transform a point from the molecule coordinate system into the drawing
  //! coordinate system
  virtual Point2D getDrawCoords(const Point2D &mol_cds) const;
  //! returns the drawing coordinates of a particular atom
  virtual Point2D getDrawCoords(int at_num) const;
  virtual Point2D getAtomCoords(const std::pair<int, int> &screen_cds) const;
  //! transform a point from drawing coordinates to the molecule coordinate
  //! system
  virtual Point2D getAtomCoords(
      const std::pair<double, double> &screen_cds) const;
  //! returns the molecular coordinates of a particular atom
  virtual Point2D getAtomCoords(int at_num) const;
  //@}
  //! return the width of the drawing area.
  virtual int width() const { return width_; }
  //! return the height of the drawing area.
  virtual int height() const { return height_; }
  //! return the width of the drawing panels.
  virtual int panelWidth() const { return panel_width_; }
  //! return the height of the drawing panels.
  virtual int panelHeight() const { return panel_height_; }
  virtual int drawHeight() const { return panel_height_ - legend_height_; }

  //! returns the drawing scale (conversion from molecular coords -> drawing
  /// coords)
  double scale() const { return scale_; }
  //! calculates the drawing scale (conversion from molecular coords -> drawing
  /// coords)
  void calculateScale(int width, int height, const ROMol &mol,
                      const std::vector<int> *highlight_atoms = nullptr,
                      const std::map<int, double> *highlight_radii = nullptr,
                      int confId = -1);
  //! overload
  /// calculate a single scale that will suit all molecules.  For use by
  /// drawMolecules primarily.
  void calculateScale(int width, int height, const std::vector<ROMol *> &mols,
                      const std::vector<std::vector<int>> *highlight_atoms,
                      const std::vector<std::map<int, double>> *highlight_radii,
                      const std::vector<int> *confIds,
                      std::vector<std::unique_ptr<RWMol>> &tmols);
  // set [xy]_trans_ to the middle of the draw area in molecule coords
  void centrePicture(int width, int height);

  //! explicitly sets the scaling factors for the drawing
  void setScale(int width, int height, const Point2D &minv, const Point2D &maxv,
                const ROMol *mol = nullptr);
  //! sets the drawing offset (in drawing coords)
  void setOffset(int x, int y) {
    x_offset_ = x;
    y_offset_ = y;
  }
  //! returns the drawing offset (in drawing coords)
  Point2D offset() const { return Point2D(x_offset_, y_offset_); }

  //! returns the minimum point of the drawing (in molecular coords)
  Point2D minPt() const { return Point2D(x_min_, y_min_); }
  //! returns the width and height of the grid (in molecular coords)
  Point2D range() const { return Point2D(x_range_, y_range_); }

  //! font size in drawing coordinate units. That's probably pixels.
  virtual double fontSize() const;
  virtual void setFontSize(double new_size);

  //! sets the current draw color
  virtual void setColour(const DrawColour &col) { curr_colour_ = col; }
  //! returns the current draw color
  virtual DrawColour colour() const { return curr_colour_; }
  //! sets the current dash pattern
  virtual void setDash(const DashPattern &patt) { curr_dash_ = patt; }
  //! returns the current dash pattern
  virtual const DashPattern &dash() const { return curr_dash_; }

  //! sets the current line width
  virtual void setLineWidth(int width) { drawOptions().bondLineWidth = width; }
  //! returns the current line width
  virtual int lineWidth() const { return drawOptions().bondLineWidth; }

  //! using the current scale, work out the size of the label in molecule
  //! coordinates.
  /*!
     Bear in mind when implementing this, that, for example, NH2 will appear as
     NH<sub>2</sub> to convey that the 2 is a subscript, and this needs to
     accounted for in the width and height.
   */
  virtual void getStringSize(const std::string &label, double &label_width,
                             double &label_height) const;
  // get the overall size of the label, allowing for it being split
  // into pieces according to orientation.
  void getLabelSize(const std::string &label, OrientType orient,
                    double &label_width, double &label_height) const;
  // return extremes for string in molecule coords.
  void getStringExtremes(const std::string &label, OrientType orient,
                         const Point2D &cds, double &x_min, double &y_min,
                         double &x_max, double &y_max) const;

  //! drawString centres the string on cds.
  virtual void drawString(const std::string &str, const Point2D &cds);
  // unless the specific drawer over-rides this overload, it will just call
  // the first one.  SVG for one needs the alignment flag.
  virtual void drawString(const std::string &str, const Point2D &cds,
                          TextAlignType align);
  //! draw a triangle
  virtual void drawTriangle(const Point2D &cds1, const Point2D &cds2,
                            const Point2D &cds3);
  //! draw an ellipse
  virtual void drawEllipse(const Point2D &cds1, const Point2D &cds2);
  // draw the arc of a circle between ang1 and ang2.  Note that 0 is
  // at 3 o-clock and 90 at 12 o'clock as you'd expect from your maths.
  // ang2 must be > ang1 - it won't draw backwards.  This is not enforced.
  // Angles in degrees.
  virtual void drawArc(const Point2D &centre, double radius, double ang1,
                       double ang2);
  // and a general ellipse form
  virtual void drawArc(const Point2D &centre, double xradius, double yradius,
                       double ang1, double ang2);
  //! draw a rectangle
  virtual void drawRect(const Point2D &cds1, const Point2D &cds2);
  //! draw a line indicating the presence of an attachment point (normally a
  //! squiggle line perpendicular to a bond)
  virtual void drawAttachmentLine(const Point2D &cds1, const Point2D &cds2,
                                  const DrawColour &col, double len = 1.0,
                                  unsigned int nSegments = 16);
  //! draw a wavy line like that used to indicate unknown stereochemistry
  virtual void drawWavyLine(const Point2D &cds1, const Point2D &cds2,
                            const DrawColour &col1, const DrawColour &col2,
                            unsigned int nSegments = 16,
                            double vertOffset = 0.05);
  //! draw a line where the ends are different colours
  virtual void drawLine(const Point2D &cds1, const Point2D &cds2,
                        const DrawColour &col1, const DrawColour &col2);
  //! adds additional information about the atoms to the output. Does not make
  //! sense for all renderers.
  virtual void tagAtoms(const ROMol &mol) { RDUNUSED_PARAM(mol); }
  //! set whether or not polygons are being filled
  virtual bool fillPolys() const { return fill_polys_; }
  //! returns either or not polygons should be filled
  virtual void setFillPolys(bool val) { fill_polys_ = val; }

  //! returns our current drawing options
  MolDrawOptions &drawOptions() { return options_; }
  //! \overload
  const MolDrawOptions &drawOptions() const { return options_; }

  //! returns the coordinates of the atoms of the current molecule in molecular
  //! coordinates
  const std::vector<Point2D> &atomCoords() const {
    PRECONDITION(activeMolIdx_ >= 0, "no index");
    return at_cds_[activeMolIdx_];
  }
  //! returns the atomic symbols of the current molecule
  const std::vector<std::pair<std::string, OrientType>> &atomSyms() const {
    PRECONDITION(activeMolIdx_ >= 0, "no index");
    return atom_syms_[activeMolIdx_];
  }
  //! Draw an arrow with either lines or a filled head (when asPolygon is true)
  virtual void drawArrow(const Point2D &cds1, const Point2D &cds2,
                         bool asPolygon = false, double frac = 0.05,
                         double angle = M_PI / 6);

  // reset to default values all the things the c'tor sets
  void tabulaRasa();

  virtual bool supportsAnnotations() { return true; }
  virtual void drawAnnotation(const AnnotationType &annotation);

  bool hasActiveAtmIdx() const { return activeAtmIdx1_ >= 0; }
  int getActiveAtmIdx1() const { return activeAtmIdx1_; }
  int getActiveAtmIdx2() const { return activeAtmIdx2_; }
  void setActiveAtmIdx(int at_idx1 = -1, int at_idx2 = -1) {
    at_idx1 = (at_idx1 < 0 ? -1 : at_idx1);
    at_idx2 = (at_idx2 < 0 ? -1 : at_idx2);
    if (at_idx2 >= 0 && at_idx1 < 0) {
      std::swap(at_idx1, at_idx2);
    }
    activeAtmIdx1_ = at_idx1;
    activeAtmIdx2_ = at_idx2;
  }
  bool hasActiveBndIdx() const { return activeBndIdx_ >= 0; }
  int getActiveBndIdx() const { return activeBndIdx_; }
  void setActiveBndIdx(int bnd_idx = -1) {
    activeBndIdx_ = (bnd_idx < 0 ? -1 : bnd_idx);
  }
  void setActiveClass(std::string actClass = std::string("")) {
    d_activeClass = actClass;
  }
  std::string getActiveClass() const { return d_activeClass; }

 protected:
  std::unique_ptr<DrawText> text_drawer_;
  std::string d_activeClass;

 private:
  bool needs_scale_;
  int width_, height_, panel_width_, panel_height_, legend_height_;
  double scale_;
  double x_min_, y_min_, x_range_, y_range_;
  double x_trans_, y_trans_;
  int x_offset_, y_offset_;  // translation in screen coordinates
  bool fill_polys_;
  int activeMolIdx_;
  int activeAtmIdx1_;
  int activeAtmIdx2_;
  int activeBndIdx_;
  std::vector<std::unique_ptr<DrawMol>> draw_mols_;

  DrawColour curr_colour_;
  DashPattern curr_dash_;
  MolDrawOptions options_;

  std::vector<std::vector<Point2D>> at_cds_;  // from mol
  std::vector<std::vector<int>> atomic_nums_;
  std::vector<std::vector<std::pair<std::string, OrientType>>> atom_syms_;
  // by the time annotations_ are drawn, we're only ever using the trans_ member
  // of the StringRect, but it is convenient to keep the whole thing rather than
  // just a StringPos for the position for calculating the scale of the drawing.
  // Went a long way down the rabbit hole before realising this, hence this
  // note.
  std::vector<std::vector<AnnotationType>> annotations_;
  std::vector<std::vector<std::pair<std::shared_ptr<StringRect>, OrientType>>>
      radicals_;
  Point2D bbox_[2];
  std::vector<std::vector<MolDrawShape>> pre_shapes_;
  std::vector<std::vector<MolDrawShape>> post_shapes_;
  bool needs_init_ = true;

  // return a DrawColour based on the contents of highlight_atoms or
  // highlight_map, falling back to atomic number by default
  DrawColour getColour(
      int atom_idx, const std::vector<int> *highlight_atoms = nullptr,
      const std::map<int, DrawColour> *highlight_map = nullptr);
  DrawColour getColourByAtomicNum(int atomic_num);

  // set the system up to draw the molecule including calculating the scale.
  std::unique_ptr<RWMol> setupDrawMolecule(
      const ROMol &mol, const std::vector<int> *highlight_atoms,
      const std::map<int, double> *highlight_radii, int confId, int width,
      int height);
  // copies of atom coords, atomic symbols etc. are stashed for convenience.
  // these put empty collections onto the stack and pop the off when done.
  void pushDrawDetails();
  void popDrawDetails();

  // Do the drawing, the new way
  void startDrawing();
  void drawAllMolecules();

  // do the initial setup bits for drawing a molecule.
  std::unique_ptr<RWMol> initMoleculeDraw(
      const ROMol &mol, const std::vector<int> *highlight_atoms,
      const std::map<int, double> *highlight_radii, int confId = -1);
  void setupTextDrawer();

  // if bond_colours is given, it must have an entry for every bond, and it
  // trumps everything else.  First in pair is bonds begin atom, second is
  // end atom.
  void drawBonds(const ROMol &draw_mol,
                 const std::vector<int> *highlight_atoms = nullptr,
                 const std::map<int, DrawColour> *highlight_atom_map = nullptr,
                 const std::vector<int> *highlight_bonds = nullptr,
                 const std::map<int, DrawColour> *highlight_bond_map = nullptr,
                 const std::vector<std::pair<DrawColour, DrawColour>>
                     *bond_colours = nullptr);
  // do the finishing touches to the drawing
  void finishMoleculeDraw(const ROMol &draw_mol,
                          const std::vector<DrawColour> &atom_colours);
  void drawLegend(const std::string &legend);
  // draw a circle in the requested colour(s) around the atom.
  void drawHighlightedAtom(int atom_idx, const std::vector<DrawColour> &colours,
                           const std::map<int, double> *highlight_radii);
  // calculate the rectangle that goes round the string, taking its
  // orientation into account.  Centre of StringRect
  // won't be the same as label_coords, necessarily, as the string might
  // be offset according to orient.
  StringRect calcLabelRect(const std::string &label, OrientType orient,
                           const Point2D &label_coords) const;
  // calculate parameters for an ellipse that roughly goes round the label
  // of the given atom.
  void calcLabelEllipse(int atom_idx,
                        const std::map<int, double> *highlight_radii,
                        Point2D &centre, double &xradius,
                        double &yradius) const;
  // annot.rect_ will have a width of -1.0 if there's a problem.
  void calcAnnotationPosition(const ROMol &mol, const Atom *atom,
                              AnnotationType &annot) const;
  void calcAnnotationPosition(const ROMol &mol, const Bond *bond,
                              AnnotationType &annot) const;
  void calcAnnotationPosition(const ROMol &mol, AnnotationType &annot) const;
  // find where to put the given annotation around an atom.  Starting
  // search at angle start_ang, in degrees.
  void calcAtomAnnotationPosition(const ROMol &mol, const Atom *atom,
                                  double start_ang,
                                  AnnotationType &annot) const;

  // draw 1 or more coloured line along bonds
  void drawHighlightedBonds(
      const ROMol &mol,
      const std::map<int, std::vector<DrawColour>> &highlight_bond_map,
      const std::map<int, int> &highlight_linewidth_multipliers,
      const std::map<int, double> *highlight_radii);
  int getHighlightBondWidth(
      int bond_idx,
      const std::map<int, int> *highlight_linewidth_multipliers) const;
  // move p2 so that the line defined by p1 to p2 touches the ellipse for the
  // atom highlighted.
  void adjustLineEndForHighlight(int at_idx,
                                 const std::map<int, double> *highlight_radii,
                                 Point2D p1, Point2D &p2) const;

  void extractAtomCoords(const ROMol &mol, int confId, bool updateBBox);
  void extractAtomSymbols(const ROMol &mol);
  void extractMolNotes(const ROMol &mol);
  void extractAtomNotes(const ROMol &mol);
  void extractBondNotes(const ROMol &mol);
  void extractRadicals(const ROMol &mol);
  void extractSGroupData(const ROMol &mol);
  void extractVariableBonds(const ROMol &mol);
  void extractBrackets(const ROMol &mol);
  void extractLinkNodes(const ROMol &mol);

  void drawAtomLabel(int atom_num,
                     const std::vector<int> *highlight_atoms = nullptr,
                     const std::map<int, DrawColour> *highlight_map = nullptr);
  OrientType calcRadicalRect(const ROMol &mol, const Atom *atom,
                             StringRect &rad_rect);
  void drawRadicals(const ROMol &mol);
  // find a good starting point for scanning round the annotation
  // atom.  If we choose well, the first angle should be the one.
  // Returns angle in radians.
  double getNoteStartAngle(const ROMol &mol, const Atom *atom) const;
  // see if the note will clash with anything else drawn on the molecule.
  // Returns 0 if no clash, 1-3 if there is a clash, denoting what clashed.
  int doesAtomNoteClash(const Point2D &note_pos,
                        const std::vector<std::shared_ptr<StringRect>> &rects,
                        const ROMol &mol, unsigned int atom_idx) const;
  int doesBondNoteClash(const Point2D &note_pos,
                        const std::vector<std::shared_ptr<StringRect>> &rects,
                        const ROMol &mol, const Bond *bond) const;
  bool doesNoteClashNbourBonds(
      const Point2D &note_pos,
      const std::vector<std::shared_ptr<StringRect>> &rects, const ROMol &mol,
      const Atom *atom) const;
  // does the note intersect with atsym, and if not, any other atom symbol.
  bool doesNoteClashAtomLabels(
      const Point2D &note_pos,
      const std::vector<std::shared_ptr<StringRect>> &rects, const ROMol &mol,
      unsigned int atom_idx) const;
  bool doesNoteClashOtherNotes(
      const Point2D &note_pos,
      const std::vector<std::shared_ptr<StringRect>> &rects) const;

  // take the coords for atnum, with neighbour nbr_cds, and move cds out to
  // accommodate
  // the label associated with it.
  void adjustBondEndForLabel(const std::pair<std::string, OrientType> &lbl,
                             const Point2D &nbr_cds, Point2D &cds) const;

  // adds LaTeX-like annotation for super- and sub-script.
  std::pair<std::string, OrientType> getAtomSymbolAndOrientation(
      const Atom &atom) const;
  std::string getAtomSymbol(const Atom &atom, OrientType orientation) const;
  OrientType getAtomOrientation(const Atom &atom) const;

  // things used by calculateScale.
  void adjustScaleForAtomLabels(const std::vector<int> *highlight_atoms,
                                const std::map<int, double> *highlight_radii);
  void adjustScaleForRadicals(const ROMol &mol);
  void adjustScaleForAnnotation(const std::vector<AnnotationType> &notes);

 private:
  virtual void updateMetadata(const ROMol &mol, int confId) {
    RDUNUSED_PARAM(mol);
    RDUNUSED_PARAM(confId);
  }
  virtual void updateMetadata(const ChemicalReaction &rxn) {
    RDUNUSED_PARAM(rxn);
  }

 protected:
  std::vector<std::pair<std::string, std::string>> d_metadata;
  unsigned int d_numMetadataEntries = 0;

  virtual void doContinuousHighlighting(
      const ROMol &mol, const std::vector<int> *highlight_atoms,
      const std::vector<int> *highlight_bonds,
      const std::map<int, DrawColour> *highlight_atom_map,
      const std::map<int, DrawColour> *highlight_bond_map,
      const std::map<int, double> *highlight_radii);

  virtual void highlightCloseContacts();
  // if bond_colours is given, it must have an entry for every bond, and it
  // trumps everything else.  First in pair is bonds begin atom, second is
  // end atom.
  virtual void drawBond(
      const ROMol &mol, const Bond *bond, int at1_idx, int at2_idx,
      const std::vector<int> *highlight_atoms = nullptr,
      const std::map<int, DrawColour> *highlight_atom_map = nullptr,
      const std::vector<int> *highlight_bonds = nullptr,
      const std::map<int, DrawColour> *highlight_bond_map = nullptr,
      const std::vector<std::pair<DrawColour, DrawColour>> *bond_colours =
          nullptr);
  virtual void drawAtomLabel(int atom_num, const DrawColour &draw_colour);
  //! DEPRECATED
  virtual void drawAnnotation(const std::string &note,
                              const StringRect &note_rect) {
    AnnotationType annot;
    annot.text_ = note;
    annot.rect_ = note_rect;
    drawAnnotation(annot);
  }

  // calculate the width to draw a line in draw coords.
  virtual double getDrawLineWidth() const;

  // sort out coords and scale for drawing reactions.
  void get2DCoordsForReaction(ChemicalReaction &rxn, Point2D &arrowBegin,
                              Point2D &arrowEnd, std::vector<double> &plusLocs,
                              double spacing, const std::vector<int> *confIds);
  // despite the name, this is only ever used for molecules in a reaction.
  void get2DCoordsMol(RWMol &mol, double &offset, double spacing, double &maxY,
                      double &minY, int confId, bool shiftAgents,
                      double coordScale);
};

// return true if the line l1s->l1f intersects line l2s->l2f.  If ip is not
// nullptr, the intersection point is stored in it.
RDKIT_MOLDRAW2D_EXPORT bool doLinesIntersect(const Point2D &l1s,
                                             const Point2D &l1f,
                                             const Point2D &l2s,
                                             const Point2D &l2f,
                                             Point2D *ip = nullptr);
// return true if line ls->lf intersects (or is fully inside) the
// rectangle of the string.
RDKIT_MOLDRAW2D_EXPORT bool doesLineIntersectLabel(const Point2D &ls,
                                                   const Point2D &lf,
                                                   const StringRect &lab_rect,
                                                   double padding = 0.0);
inline void setDarkMode(MolDrawOptions &opts) {
  assignDarkModePalette(opts.atomColourPalette);
  opts.backgroundColour = DrawColour{0, 0, 0, 1};
  opts.annotationColour = DrawColour{0.9, 0.9, 0.9, 1};
  opts.legendColour = DrawColour{0.9, 0.9, 0.9, 1};
  opts.symbolColour = DrawColour{0.9, 0.9, 0.9, 1};
  opts.variableAttachmentColour = DrawColour{0.3, 0.3, 0.3, 1};
}
inline void setDarkMode(MolDraw2D &d2d) { setDarkMode(d2d.drawOptions()); }
inline void setMonochromeMode(MolDrawOptions &opts, const DrawColour &fgColour,
                              const DrawColour &bgColour) {
  opts.atomColourPalette.clear();
  opts.atomColourPalette[-1] = fgColour;
  opts.backgroundColour = bgColour;
  opts.annotationColour = fgColour;
  opts.legendColour = fgColour;
  opts.symbolColour = fgColour;
  opts.variableAttachmentColour = fgColour;
}
inline void setMonochromeMode(MolDraw2D &opts, const DrawColour &fgColour,
                              const DrawColour &bgColour) {
  setMonochromeMode(opts.drawOptions(), fgColour, bgColour);
}
}  // namespace RDKit

#endif  // RDKITMOLDRAW2D_H
