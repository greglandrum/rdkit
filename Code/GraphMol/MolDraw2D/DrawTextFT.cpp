//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//
// Original author: David Cosgrove (CozChemIx) on 06/05/2020.
//

#include <cstdio>
#include <iostream>

#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/DrawTextFT.h>

using namespace std;

namespace RDKit {

// ****************************************************************************
  DrawTextFT::DrawTextFT(double max_fnt_sz, double min_fnt_sz,
                         const string &font_file) :
    DrawText(max_fnt_sz, min_fnt_sz), library_(nullptr), face_(nullptr),
    x_trans_(0), y_trans_(0), string_y_max_(0) {

  string fontfile = getFontFile();

  int err_code = FT_Init_FreeType(&library_);
  if(err_code != FT_Err_Ok) {
    throw runtime_error(string("Couldn't initialise Freetype."));
  }
//  // take the first face
//  err_code = FT_New_Face(library_, fontfile.c_str(), 0, &face_);
//  if(err_code != FT_Err_Ok) {
//    throw runtime_error(string("Font file ") + fontfile + string(" not found."));
//  }
//  em_scale_ = 1.0 / face_->units_per_EM;
  setFontFile(font_file);

}

// ****************************************************************************
DrawTextFT::~DrawTextFT() {

  FT_Done_Face(face_);
  FT_Done_FreeType(library_);

}

// ****************************************************************************
void DrawTextFT::drawChar(char c, const Point2D &cds) {

  FT_Load_Char(face_, c, FT_LOAD_NO_SCALE | FT_LOAD_NO_BITMAP);
  x_trans_ = cds.x;
  y_trans_ = cds.y;
  extractOutline();

}

// ****************************************************************************
double DrawTextFT::fontCoordToDrawCoord(FT_Pos fc) const {

  double pc = fontSize() * fc * em_scale_;

  return pc;

}

// ****************************************************************************
void DrawTextFT::fontPosToDrawPos(FT_Pos fx, FT_Pos fy,
                                  double &dx, double &dy) const {

  dx = x_trans_ + fontCoordToDrawCoord(fx);
  // freetype has the origin at bottom left, we want it in the top left.
  FT_Pos rev_y = string_y_max_ - fy;
  dy = y_trans_ + fontCoordToDrawCoord(rev_y);

}

// ****************************************************************************
double DrawTextFT::extractOutline() {

  FT_Outline_Funcs callbacks;

  callbacks.move_to = moveToFunction;
  callbacks.line_to = lineToFunction;
  callbacks.conic_to = conicToFunction;
  callbacks.cubic_to = cubicToFunction;

  callbacks.shift = 0;
  callbacks.delta = 0;

  FT_GlyphSlot slot = face_->glyph;
  FT_Outline &outline = slot->outline;

  FT_Error error = FT_Outline_Decompose(&outline, &callbacks, this);
  if(error != FT_Err_Ok) {
    /* not sure what to do in this case */ ;
  }
  return fontCoordToDrawCoord(slot->advance.x);

}

// ****************************************************************************
string DrawTextFT::getFontFile() const {

  if(!font_file_.empty()) {
    return font_file_;
  }
  string ff_name = getenv("RDBASE") ? getenv("RDBASE") : "";
  if(ff_name.empty()) {
    return "NO_FONT_FILE"; // the c'tor will throw when this file isn't found.
  }
  ff_name += "/Code/GraphMol/MolDraw2D/Telex-Regular.ttf";
  return ff_name;

}

// ****************************************************************************
void DrawTextFT::setFontFile(const string &font_file) {

  if(face_ && font_file == font_file_) {
    return;
  }

  font_file_ = font_file;
  string ff = getFontFile();
  FT_Done_Face(face_);
  // take the first face
  int err_code = FT_New_Face(library_, ff.c_str(), 0, &face_);
  if(err_code != FT_Err_Ok) {
    throw runtime_error(string("Font file ") + ff + string(" not found."));
  }
  em_scale_ = 1.0 / face_->units_per_EM;

}

// ****************************************************************************
void DrawTextFT::getStringRects(const string &text,
                                vector<shared_ptr<StringRect>> &rects,
                                vector<TextDrawType> &draw_modes,
                                vector<char> &draw_chars) const {

  TextDrawType draw_mode = TextDrawType::TextDrawNormal;
  double running_x = 0.0, max_y = 0.0;
  for(size_t i = 0; i < text.length(); ++i) {
    // setStringDrawMode moves i along to the end of any <sub> or <sup>
    // markup
    if ('<' == text[i] && setStringDrawMode(text, draw_mode, i)) {
      continue;
    }
    draw_chars.push_back(text[i]);
    FT_Pos this_x_min, this_y_min, this_x_max, this_y_max, advance;
    calcGlyphBBox(text[i], this_x_min, this_y_min, this_x_max, this_y_max,
                  advance);
    double oscale = selectScaleFactor(text[i], draw_mode);
    double p_x_min = oscale * fontCoordToDrawCoord(this_x_min);
    double p_y_min = oscale * fontCoordToDrawCoord(this_y_min);
    double p_x_max = oscale * fontCoordToDrawCoord(this_x_max);
    double p_y_max = oscale * fontCoordToDrawCoord(this_y_max);
    double width = p_x_max - p_x_min;
    double height = p_y_max - p_y_min;
    Point2D offset(p_x_min + width / 2.0, p_y_max / 2.0);
    Point2D g_centre(offset.x, p_y_max - height / 2.0);
    rects.push_back(shared_ptr<StringRect>(new StringRect(offset, g_centre, width, height)));
    rects.back()->trans_.x = running_x;
    draw_modes.push_back(draw_mode);
    running_x += p_x_max;
    max_y = max(max_y, p_y_max);
  }
  for(auto r: rects) {
    r->g_centre_.y = max_y - r->g_centre_.y;
    r->offset_.y = max_y / 2.0;
  }

  adjustStringRectsForSuperSubScript(draw_modes, rects);

}

// ****************************************************************************
void DrawTextFT::calcGlyphBBox(char c, FT_Pos &x_min, FT_Pos &y_min,
                               FT_Pos &x_max , FT_Pos &y_max,
                               FT_Pos &advance) const {

  FT_Load_Char(face_, c, FT_LOAD_NO_SCALE | FT_LOAD_NO_BITMAP);
  FT_GlyphSlot slot = face_->glyph;
  FT_Outline &outline = slot->outline;
  FT_BBox bbox;

  FT_Outline_Get_BBox(&outline, &bbox);

  x_min = bbox.xMin;
  y_min = bbox.yMin;
  x_max = bbox.xMax;
  y_max = bbox.yMax;

  advance = slot->advance.x;

}

// ****************************************************************************
int moveToFunction(const FT_Vector *to, void *user) {

  auto *rdft = static_cast<DrawTextFT *>(user);
  return rdft->MoveToFunctionImpl(to);

}

// ****************************************************************************
int lineToFunction(const FT_Vector *to, void *user) {

  auto *rdft = static_cast<DrawTextFT *>(user);
  return rdft->LineToFunctionImpl(to);

}

// ****************************************************************************
int conicToFunction(const FT_Vector *control, const FT_Vector *to,
                    void *user) {

  auto *rdft = static_cast<DrawTextFT *>(user);
  return rdft->ConicToFunctionImpl(control, to);

}

// ****************************************************************************
int cubicToFunction(const FT_Vector *controlOne,
                    const FT_Vector *controlTwo,
                    const FT_Vector *to, void *user) {

  auto *rdft = static_cast<DrawTextFT *>(user);
  return rdft->CubicToFunctionImpl(controlOne, controlTwo, to);

}

} // namespace RDKit
