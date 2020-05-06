//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Original author: David Cosgrove (CozChemIx) on 29/04/2020.
//
// A concrete class derived from DrawText that uses SVG
// to draw text onto a picture.

#include <sstream>
#include <boost/algorithm/string.hpp>

#include <GraphMol/MolDraw2D/DrawTextSVG.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <GraphMol/MolDraw2D/MolDraw2DDetails.h>

using namespace std;

namespace RDKit {

string DrawColourToSVG(const RDKit::DrawColour &col);

// ****************************************************************************
DrawTextSVG::DrawTextSVG(ostream &oss, string &d_act_class) :
    oss_(oss), d_active_class_(d_act_class) {

}

// ****************************************************************************
void DrawTextSVG::getStringSize(const string &label, double &label_width,
                                double &label_height) const {

  label_width = 0.0;
  label_height = 0.0;

  MolDraw2D::TextDrawType draw_mode = MolDraw2D::TextDrawNormal;

  bool had_a_super = false;
  bool had_a_sub = false;

  double act_font_size = fontSize();
  for (int i = 0, is = label.length(); i < is; ++i) {
    // setStringDrawMode moves i along to the end of any <sub> or <sup>
    // markup
    if ('<' == label[i] && setStringDrawMode(label, draw_mode, i)) {
      continue;
    }

    label_height = act_font_size;
    double char_width =
        act_font_size *
        static_cast<double>(MolDraw2D_detail::char_widths[(int)label[i]]) /
        MolDraw2D_detail::char_widths[(int)'M'];
    if (MolDraw2D::TextDrawSubscript == draw_mode) {
      char_width *= 0.5;
      had_a_sub = true;
    } else if (MolDraw2D::TextDrawSuperscript == draw_mode) {
      char_width *= 0.5;
      had_a_super = true;
    }
    label_width += char_width;
  }

  // subscript keeps its bottom in line with the bottom of the bit chars,
  // superscript goes above the original char top by a bit (empirical)
  if (had_a_super) {
    label_height *= 1.1;
  }
  if (had_a_sub) {
    label_height *= 1.1;
  }

}

namespace {
void escape_xhtml(std::string &data) {
  boost::algorithm::replace_all(data, "&", "&amp;");
  boost::algorithm::replace_all(data, "\"", "&quot;");
  boost::algorithm::replace_all(data, "\'", "&apos;");
  boost::algorithm::replace_all(data, "<", "&lt;");
  boost::algorithm::replace_all(data, ">", "&gt;");
}
}  // namespace

// ****************************************************************************
void DrawTextSVG::drawString(const std::string &str, const Point2D &cds,
                             MolDraw2D::AlignType align) {

  std::cout << "DrawTextSVG::drawString " << str << " at "
      << cds.x << ", " << cds.y << std::endl;

  unsigned int fontSz = fontSize();
  double string_width, string_height;
  getStringSize(str, string_width, string_height);
  std::cout << "string size : " << string_width << " by " << string_height << std::endl;
  std::string col = DrawColourToSVG(colour());

  std::string text_anchor = "middle";
  double tmult = 0.0;
  if(align == MolDraw2D::END) {
    text_anchor = "end";
    tmult = -1.0;
  } else if(align == MolDraw2D::START) {
    text_anchor = "start";
    tmult = 1.0;
  }
  Point2D draw_coords = Point2D(cds.x, cds.y);
  // fonts are laid out with room for wider letters like W and hanging bits like g.
  // Very few atomic symbols need to care about this, and common ones look a bit
  // out of line.  For example O sits to the left of a double bond.  This is an
  // empirical tweak to push it back a bit.  Use the string_height in x and y
  // as it's just a correction factor based on the scaled font size.
  draw_coords.x += string_height * tmult * 0.1;
  draw_coords.y += string_height * 0.15;

  oss_ << "<text dominant-baseline=\"central\" text-anchor=\""
       << text_anchor << "\"";
  oss_ << " x='" << draw_coords.x;
  oss_ << "' y='" << draw_coords.y << "'";

  if (!d_active_class_.empty()) {
    oss_ << " class='" << d_active_class_ << "'";
  }
  oss_ << " style='font-size:" << fontSz
       << "px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;"
          "font-family:sans-serif;"
       << "fill:" << col << "'";
  oss_ << " >";

  MolDraw2D::TextDrawType draw_mode = MolDraw2D::TextDrawNormal;
  std::string span;
  bool first_span = true;
  auto write_span = [&]() {
    if (!first_span) {
      escape_xhtml(span);
      oss_ << span << "</tspan>";
      span = "";
    }
    first_span = false;
  };

  for (int i = 0, is = str.length(); i < is; ++i) {
    // setStringDrawMode moves i along to the end of any <sub> or <sup>
    // markup
    if ('<' == str[i] && setStringDrawMode(str, draw_mode, i)) {
      write_span();
      oss_ << "<tspan";
      switch (draw_mode) {
        // To save people time later - on macOS Catalina, at least, Firefox
        // renders the superscript as a subscript.  It's fine on Safari.
        case MolDraw2D::TextDrawSuperscript:
          oss_ << " style='baseline-shift:super;font-size:" << fontSz * 0.75
               << "px;"
               << "'";
          break;
        case MolDraw2D::TextDrawSubscript:
          oss_ << " style='baseline-shift:sub;font-size:" << fontSz * 0.75
               << "px;"
               << "'";
          break;
        default:
          break;
      }
      oss_ << ">";
      continue;
    }
    if (first_span) {
      first_span = false;
      oss_ << "<tspan>";
      span = "";
    }
    span += str[i];
  }
  escape_xhtml(span);
  oss_ << span << "</tspan>";
  oss_ << "</text>\n";

}

// ****************************************************************************
// draw the char, with the bottom left hand corner at cds
void DrawTextSVG::drawChar(char c, const Point2D &cds) {
  unsigned int fontSz = fontSize();
  std::string col = DrawColourToSVG(colour());

  oss_ << "<text";
  oss_ << " x='" << cds.x;
  oss_ << "' y='" << cds.y << "'";
  oss_ << " style='font-size:" << fontSz
       << "px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;"
          "font-family:sans-serif;text-anchor:start;"
       << "fill:" << col << "'";
  oss_ << " >";
  oss_ << c;
  oss_ << "</text>";

}

} // namespace RDKit