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
#include <TGraph.h>
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
  static boost::posix_time::ptime day1( boost::gregorian::date(2024, 10, 21));
  static boost::posix_time::ptime day2( boost::gregorian::date(2024, 11, 29));
  static boost::posix_time::ptime day3( boost::gregorian::date(2024, 12, 27));

  const auto cells = ParseCsv(sinput);

  std::vector<double> vsler, viler;
  std::vector<double> vsher, viher;
  std::vector<double> vlumh, vluml;

  const auto her2ler = 1.83/2.58;
  const auto Lsph    = 5.0e+31 * 1.0e+6 / (1.0e+34);
  const auto Lspl    = 4.0e+31 * 1.0e+6 / (1.0e+34);
  const auto nb2346  = 2346.;
  const auto nb393   =  393.;

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

      // Collision starts on Day 1
      if (ptimed > day1)
        {
          vsher.push_back(sec);
          viher.push_back(iher);

          // Filling scheme
          if (iler < 1.4)
            {
        	  vlumh.push_back(iler/nb393  * iher/nb393  * Lsph * nb393);
              vluml.push_back(iler/nb393  * iher/nb393  * Lspl * nb393);
            }
          else
            {
              vlumh.push_back(iler/nb2346 * iher/nb2346 * Lsph * nb2346);
              vluml.push_back(iler/nb2346 * iher/nb2346 * Lspl * nb2346);
            }
        }
    }

  TCanvas *c0 = new TCanvas("c0", "c0", 1000, 500);
  c0->cd();

  auto giler = boost::shared_ptr<TGraph>(new TGraph(vsler.size(), &vsler[0], &viler[0]));
  giler->SetName("");
  giler->SetTitle("2024C");
  giler->SetMarkerStyle(24);
  giler->SetMarkerColor(2);
  giler->SetLineColor(2);
  giler->Draw("AP");
  giler->GetHistogram()->SetMinimum( 0.);
  giler->GetHistogram()->SetMaximum(11.);
  giler->GetXaxis()->SetTimeDisplay(1);
  giler->GetXaxis()->SetTimeFormat("%b %d");

  auto giher = boost::shared_ptr<TGraph>(new TGraph(vsher.size(), &vsher[0], &viher[0]));
  giher->SetMarkerStyle(24);
  giher->SetMarkerColor(4);
  giher->SetLineColor(4);
  giher->Draw("P,same");

  auto glumh = boost::shared_ptr<TGraph>(new TGraph(vsher.size(), &vsher[0], &vlumh[0]));
  glumh->SetMarkerStyle(24);
  glumh->SetMarkerColor(4);
  glumh->SetLineColor(4);
  glumh->Draw("L,same");

  c0->Update();

  theApp.Run();
  return 0;
}
