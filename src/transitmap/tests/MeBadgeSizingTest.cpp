#include <algorithm>
#include <cmath>
#include <sstream>
#include <string>

#include "shared/rendergraph/RenderGraph.h"
#include "transitmap/output/SvgRenderer.h"
#include "transitmap/tests/MeBadgeSizingTest.h"
#include "util/Misc.h"
#include "util/String.h"
#include "util/geo/Geo.h"

using shared::rendergraph::RenderGraph;
using transitmapper::config::Config;
using transitmapper::output::SvgRenderer;

namespace {

double extractFontSize(const std::string &svg, const std::string &label) {
  size_t labelPos = svg.find(">" + label + "<");
  TEST(labelPos != std::string::npos);
  size_t textStart = svg.rfind("<text", labelPos);
  TEST(textStart != std::string::npos);
  size_t fontPos = svg.find("font-size=\"", textStart);
  TEST(fontPos != std::string::npos);
  fontPos += std::string("font-size=\"").size();
  size_t fontEnd = svg.find("\"", fontPos);
  TEST(fontEnd != std::string::npos);
  return std::stod(svg.substr(fontPos, fontEnd - fontPos));
}

double extractRectAttribute(const std::string &svg, const std::string &fillValue,
                            const std::string &attribute) {
  std::string fillMarker = "fill=\"" + fillValue + "\"";
  size_t fillPos = svg.find(fillMarker);
  TEST(fillPos != std::string::npos);
  size_t rectStart = svg.rfind("<rect", fillPos);
  TEST(rectStart != std::string::npos);
  std::string attrMarker = attribute + "=\"";
  size_t attrPos = svg.find(attrMarker, rectStart);
  TEST(attrPos != std::string::npos);
  attrPos += attrMarker.size();
  size_t attrEnd = svg.find("\"", attrPos);
  TEST(attrEnd != std::string::npos);
  return std::stod(svg.substr(attrPos, attrEnd - attrPos));
}

double computeBadgeHeight(double starPx, double labelHeightPx) {
  double textHeightForPadding = labelHeightPx > 0.0 ? labelHeightPx : starPx;
  double padTop = textHeightForPadding * 0.28;
  double padBottom = textHeightForPadding * 0.12;
  double contentHeight = std::max(starPx, textHeightForPadding);
  return padTop + padBottom + contentHeight;
}

double computeBadgeWidth(double starPx, double labelHeightPx,
                         size_t labelCpCount) {
  double textHeightForPadding = labelHeightPx > 0.0 ? labelHeightPx : starPx;
  double padX = textHeightForPadding * 0.6;
  double starGapPx = starPx * 0.2;
  double labelWidthPx = labelCpCount * (labelHeightPx * 0.6);
  return padX * 2.0 + starPx + starGapPx + labelWidthPx;
}

}  // namespace

void MeBadgeSizingTest::run() {
  Config baseCfg;
  baseCfg.outputResolution = 1.0;
  baseCfg.renderMe = true;
  baseCfg.renderMeLabel = true;
  baseCfg.meStationWithBg = true;
  baseCfg.meStationBgFill = "#abc123";
  baseCfg.meStationBgStroke = "#000000";
  baseCfg.meStationTextColor = "#123456";
  baseCfg.meStationFill = "#654321";
  baseCfg.meLandmark.coord = util::geo::DPoint(0, 0);
  baseCfg.meLandmark.label = "Here";
  baseCfg.meLandmark.fontSize = baseCfg.meLabelSize;
  baseCfg.meLandmark.color = baseCfg.meStationFill;

  RenderGraph g;

  Config autoCfg = baseCfg;
  autoCfg.meStarSize = 75.0;
  autoCfg.meStarSizeExplicit = true;
  autoCfg.meLabelSizeExplicit = false;
  autoCfg.meLandmark.fontSize = autoCfg.meLabelSize;

  std::ostringstream autoSvgOut;
  SvgRenderer autoRenderer(&autoSvgOut, &autoCfg);
  autoRenderer.print(g);
  std::string autoSvg = autoSvgOut.str();
  double autoFontSize = extractFontSize(autoSvg, autoCfg.meLandmark.label);
  TEST(baseCfg.meStarSize > 0.0);
  double expectedAutoFont = baseCfg.meLabelSize *
                            (autoCfg.meStarSize / baseCfg.meStarSize);
  TEST(std::abs(autoFontSize - expectedAutoFont) < 1e-6);
  size_t labelCpCount = util::toWStr(autoCfg.meLandmark.label).size();
  double starPx = autoCfg.meStarSize * autoCfg.outputResolution;
  double expectedAutoHeight = computeBadgeHeight(starPx, autoFontSize);
  double expectedAutoWidth =
      computeBadgeWidth(starPx, autoFontSize, labelCpCount);
  double rectAutoHeight =
      extractRectAttribute(autoSvg, autoCfg.meStationBgFill, "height");
  double rectAutoWidth =
      extractRectAttribute(autoSvg, autoCfg.meStationBgFill, "width");
  TEST(std::abs(rectAutoHeight - expectedAutoHeight) < 1e-6);
  TEST(std::abs(rectAutoWidth - expectedAutoWidth) < 1e-6);

  Config explicitCfg = autoCfg;
  explicitCfg.meLabelSizeExplicit = true;
  explicitCfg.meLabelSize = 60.0;
  explicitCfg.meLandmark.fontSize = explicitCfg.meLabelSize;

  std::ostringstream explicitSvgOut;
  SvgRenderer explicitRenderer(&explicitSvgOut, &explicitCfg);
  explicitRenderer.print(g);
  std::string explicitSvg = explicitSvgOut.str();
  double explicitFontSize =
      extractFontSize(explicitSvg, explicitCfg.meLandmark.label);
  TEST(std::abs(explicitFontSize - explicitCfg.meLabelSize) < 1e-6);
  double explicitStarPx = explicitCfg.meStarSize * explicitCfg.outputResolution;
  double expectedExplicitHeight =
      computeBadgeHeight(explicitStarPx, explicitFontSize);
  double expectedExplicitWidth =
      computeBadgeWidth(explicitStarPx, explicitFontSize, labelCpCount);
  double rectExplicitHeight =
      extractRectAttribute(explicitSvg, explicitCfg.meStationBgFill, "height");
  double rectExplicitWidth =
      extractRectAttribute(explicitSvg, explicitCfg.meStationBgFill, "width");
  TEST(std::abs(rectExplicitHeight - expectedExplicitHeight) < 1e-6);
  TEST(std::abs(rectExplicitWidth - expectedExplicitWidth) < 1e-6);
}
