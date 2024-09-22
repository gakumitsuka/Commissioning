#include <iostream>
#include <fstream>
#include <cassert>

#include <boost/array.hpp>
#include <boost/format.hpp>
#include <boost/date_time.hpp>
#include <boost/filesystem.hpp>
#include <boost/multi_array.hpp>
#include <boost/program_options.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include <TAxis.h>
#include <TCanvas.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TH1.h>
#include <TRint.h>
#include <TStyle.h>

using boost::program_options::options_description;
using boost::program_options::value;
using boost::program_options::variables_map;
using boost::program_options::store;
using boost::program_options::parse_command_line;
using boost::program_options::notify;


std::vector < std::vector<std::string> > ParseCsv(const std::string &filepath)
{
  std::string line;
  std::ifstream ifs(filepath);

  std::vector< std::vector<std::string> > cells;
  while (std::getline(ifs, line))
    {
      std::vector<std::string> data;

      boost::tokenizer< boost::escaped_list_separator< char > > tokens(line);
      for (const std::string& it : tokens)
        data.push_back(it);

      cells.push_back(data);
    }
  cells.erase(cells.begin()); // remove header line

  return cells;
}


int main(int argc,char *argv[])
{
  /* ROOT initialization */
  int tmp=1;
  TApplication theApp("App",&tmp,argv);

  int fontid = 42;
  gStyle->SetStatFont(fontid);
  gStyle->SetLabelFont(fontid,"XYZ");
  gStyle->SetLabelFont(fontid,"");
  gStyle->SetTitleFont(fontid,"XYZ");
  gStyle->SetTitleFont(fontid,"");
  gStyle->SetTextFont(fontid);
  gStyle->SetLegendFont(fontid);
  gStyle->SetOptStat(0);

  /*
   * Command line option
   */
  options_description opt("Option");
  // define format
  opt.add_options()
  ("help,h", "display help")
  ("input,i",  value<std::string>(), "input csv file");

  // analyze command line
  variables_map argmap;
  store(parse_command_line(argc, argv, opt), argmap);
  notify(argmap);

  // if no matching option, show help
  if (argmap.count("help") || !argmap.count("input"))
    {
      std::cerr << opt << std::endl;
      return 1;
    }

  const std::string sinput = argmap["input"].as<std::string>();

  static boost::posix_time::ptime epoch(boost::gregorian::date(1970,  1,  1));
  static boost::posix_time::ptime day0( boost::gregorian::date(2024, 10,  9));
  static boost::posix_time::ptime day1( boost::gregorian::date(2024, 10, 22));
  static boost::posix_time::ptime day2( boost::gregorian::date(2024, 11, 29));
  static boost::posix_time::ptime day3( boost::gregorian::date(2024, 12, 27));

  const auto cells = ParseCsv(sinput);

  std::vector<double> vsler, viler, vbler;
  std::vector<double> vsher, viher, vbher;
  std::vector<double> vscol, vbpro, vlumi;

  const auto her2ler = 1.83/2.58;
  const auto nb2346  = 2346.;
  const auto nb393   =  393.;

  auto Lsph = [](const double &mA2)
              {
                return (4.77727 * exp(-2.4878 * mA2) + 4.33312) * 1.0e+31;
              };

  auto Lspl = [](const double &mA2)
              {
                return (4.77727 * exp(-2.4878 * mA2) + 4.33312) * 1.0e+31 / (1.0e+34);
              };

  for (const auto& data : cells)
    {
      assert(data.size()==2 && "data's column size should be 2.");

      auto *formatd = new boost::posix_time::time_input_facet("%Y-%m-%d");
      boost::posix_time::ptime ptimed;

      const double iler = boost::lexical_cast<double>(data[0]);
      const double iher = iler * her2ler;

      std::istringstream iss(data[1]);
      iss.imbue(std::locale(iss.getloc(), formatd));
      iss >> ptimed;
      delete formatd;

      auto to_time_t = [](const boost::posix_time::ptime &date)
                       {
                         return time_t((date - epoch).total_seconds());
                       };
      const auto sec = boost::lexical_cast<double>(to_time_t(ptimed) - 788918400);

      vsler.push_back(sec);  // ref to 1 Jan, 1995
      viler.push_back(iler);
      vbler.push_back(iler/nb2346 * 1000.);

      vsher.push_back(sec);
      viher.push_back(iher);
      vbher.push_back(iher/nb2346 * 1000.);

      // Collision starts on Day 1
      vscol.push_back(sec);

      if (ptimed > day1)
        {
          // Filling scheme
          //          if (iler < 1.4)
          //            {
          //              const auto mA2 = iler/nb393  * iher/nb393  * 1.0e+6;
          //              vlumi.push_back(mA2  * Lsph(mA2) * nb393);
          //              vluml.push_back(mA2  * Lspl(mA2) * nb393);
          //            }
          //          else
          //            {
          vbpro.push_back(iler/nb2346 * 1000. * iher/nb2346 * 1000.);

          const auto mA2 = iler/nb2346 * iher/nb2346 * 1.0e+6;
          vlumi.push_back(mA2  * Lsph(mA2) * nb2346 / 1.0e+35); // normalized to 1e+35
          //            }
        }
      else
        {
          vbpro.push_back(-100);
          vlumi.push_back(-100);
        }
    }

  //--- Beam current ---
  TCanvas *c0 = new TCanvas("c0", "c0", 700, 500);
  c0->cd();

  // LER beam current
  auto giler = boost::shared_ptr<TGraph>(new TGraph(vsler.size(), &vsler[0], &viler[0]));
  giler->SetName("");
  giler->SetTitle("2024C");
  giler->SetMarkerStyle(24);
  giler->SetMarkerColor(2);
  giler->SetLineColor(2);
  giler->Draw("AL");
  giler->GetHistogram()->SetMinimum(0.);
  giler->GetHistogram()->SetMaximum(3.);
  auto *ax0 = giler->GetXaxis();
  ax0->SetLimits(vsler.front()-86400, vsler.back()+86400);
  giler->GetXaxis()->SetTimeDisplay(1);
  giler->GetXaxis()->SetNdivisions(507);
  giler->GetXaxis()->SetTimeFormat("%b %d");
  giler->GetXaxis()->SetLabelSize(0.05);
  giler->GetXaxis()->SetLabelOffset(0.02);
  giler->GetYaxis()->SetTitle("Beam current (mA)");
  giler->GetYaxis()->SetTitleSize(0.05);
  giler->GetYaxis()->SetTitleOffset(0.9);
  giler->GetYaxis()->CenterTitle(true);
  giler->GetYaxis()->SetLabelSize(0.05);

  // HER beam current
  auto giher = boost::shared_ptr<TGraph>(new TGraph(vsher.size(), &vsher[0], &viher[0]));
  giher->SetMarkerStyle(24);
  giher->SetMarkerColor(4);
  giher->SetLineColor(4);
  giher->Draw("L,same");

  auto* leg0 = new TLegend(0.1446991,0.6659664,0.4212034,0.8823529,NULL,"brNDC");
  leg0->AddEntry(giler.get(), "LER", "l");
  leg0->AddEntry(giher.get(), "HER", "l");
  leg0->SetBorderSize(0);
  leg0->SetFillColor(0);
  leg0->SetTextFont(42);
  leg0->Draw();

  c0->Update();

  //--- Bunch current ---
  TCanvas *c1 = new TCanvas("c1", "c1", 700, 500);
  c1->cd();

  // LER bunch current
  auto gbler = boost::shared_ptr<TGraph>(new TGraph(vsler.size(), &vsler[0], &vbler[0]));
  gbler->SetName("");
  gbler->SetTitle("2024C");
  gbler->SetMarkerStyle(24);
  gbler->SetMarkerColor(2);
  gbler->SetLineColor(2);
  gbler->Draw("AL");
  gbler->GetHistogram()->SetMinimum(0.);
  gbler->GetHistogram()->SetMaximum(1.2);
  auto *ax1 = gbler->GetXaxis();
  ax1->SetLimits(vsler.front()-86400, vsler.back()+86400);
  gbler->GetXaxis()->SetTimeDisplay(1);
  gbler->GetXaxis()->SetNdivisions(507);
  gbler->GetXaxis()->SetTimeFormat("%b %d");
  gbler->GetXaxis()->SetLabelSize(0.05);
  gbler->GetXaxis()->SetLabelOffset(0.02);
  gbler->GetYaxis()->SetTitle("Bunch current (mA)");
  gbler->GetYaxis()->SetTitleSize(0.05);
  gbler->GetYaxis()->SetTitleOffset(0.9);
  gbler->GetYaxis()->CenterTitle(true);
  gbler->GetYaxis()->SetLabelSize(0.05);

  // HER bunch current
  auto gbher = boost::shared_ptr<TGraph>(new TGraph(vsher.size(), &vsher[0], &vbher[0]));
  gbher->SetMarkerStyle(24);
  gbher->SetMarkerColor(4);
  gbher->SetLineColor(4);
  gbher->Draw("L,same");

  // Bunch current product
  auto gbpro = boost::shared_ptr<TGraph>(new TGraph(vscol.size(), &vscol[0], &vbpro[0]));
  gbpro->SetLineStyle(7);
  gbpro->SetLineColor(1);
  gbpro->Draw("L,same");

  // Bunch current product axis
  TGaxis *axisr = new TGaxis(vscol.back()+86400, 0.0, vscol.back()+86400, 1.2, 0, 1.2, 510, "+L");
  axisr->SetTitle("Bunch current product (mA^{2})");
  axisr->SetTitleSize(0.05);
  axisr->SetTitleOffset(0.9);
  axisr->SetTitleFont(42);
  axisr->CenterTitle(true);
  axisr->SetLabelSize(0.05);
  axisr->SetLabelFont(42);
  axisr->Draw();

  auto* leg1 = new TLegend(0.1446991,0.6659664,0.6346705,0.8823529,NULL,"brNDC");
  leg1->AddEntry(gbler.get(), "LER", "l");
  leg1->AddEntry(gbher.get(), "HER", "l");
  leg1->AddEntry(gbpro.get(), "Bunch current product", "l");
  leg1->SetBorderSize(0);
  leg1->SetFillColor(0);
  leg1->SetTextFont(42);
  leg1->Draw();

  c1->Update();

  //--- Luminosity ---
  TCanvas *c2 = new TCanvas("c2", "c2", 700, 500);
  c2->cd();

  // Luminosity
  auto glumi = boost::shared_ptr<TGraph>(new TGraph(vscol.size(), &vscol[0], &vlumi[0]));
  glumi->SetName("");
  glumi->SetTitle("2024C");
  glumi->SetMarkerStyle(24);
  glumi->SetMarkerColor(2);
  glumi->SetLineColor(44);
  glumi->Draw("AL");
  glumi->GetHistogram()->SetMinimum(0.);
  glumi->GetHistogram()->SetMaximum(1.2);
  auto *ax2 = glumi->GetXaxis();
  ax2->SetLimits(vscol.front()-86400, vscol.back()+86400);
  glumi->GetXaxis()->SetTimeDisplay(1);
  glumi->GetXaxis()->SetNdivisions(507);
  glumi->GetXaxis()->SetTimeFormat("%b %d");
  glumi->GetXaxis()->SetLabelSize(0.05);
  glumi->GetXaxis()->SetLabelOffset(0.02);
  glumi->GetYaxis()->SetTitle("Luminosity (10^{35} cm^{-2} s^{-1})");
  glumi->GetYaxis()->SetTitleSize(0.05);
  glumi->GetYaxis()->SetTitleOffset(0.9);
  glumi->GetYaxis()->CenterTitle(true);
  glumi->GetYaxis()->SetLabelSize(0.05);

  // Bunch current product
  gbpro->Draw("L,same");

  // Bunch current product axis
  axisr->Draw();

  auto* leg2 = new TLegend(0.1446991,0.6659664,0.6848138,0.8823529,NULL,"brNDC");
  leg2->AddEntry(glumi.get(), "Luminosity", "l");
  leg2->AddEntry(gbpro.get(), "Bunch current product", "l");
  leg2->SetBorderSize(0);
  leg2->SetFillColor(0);
  leg2->SetTextFont(42);
  leg2->Draw();

  c2->Update();

  theApp.Run();
  return 0;
}
