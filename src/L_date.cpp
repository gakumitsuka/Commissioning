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

#include <TFile.h>
#include <TGraph.h>

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
  /*
   * Command line option
   */
  options_description opt("Option");
  // define format
  opt.add_options()
  ("help,h", "display help")
  ("ini,i",  value<std::string>(), "ini file");

  // analyze command line
  variables_map argmap;
  store(parse_command_line(argc, argv, opt), argmap);
  notify(argmap);

  // if no matching option, show help
  if (argmap.count("help") || !argmap.count("ini"))
    {
      std::cerr << opt << std::endl;
      return 1;
    }

  const std::string ini = argmap["ini"].as<std::string>();

  boost::property_tree::ptree pt;
  try
    {
      boost::property_tree::read_ini(ini.c_str(), pt);
    }
  catch (boost::property_tree::ptree_error& e)
    {
      std::cout << "ptree_error " << e.what() << std::endl;
      exit(-1);
    }

  std::vector<std::string> sINI = {pt.get<std::string>("PV.RATE"), pt.get<std::string>("PV.CURRENT"), pt.get<std::string>("PV.CHARGE1"), pt.get<std::string>("PV.EFFICIENCY")};

  auto fout = boost::shared_ptr<TFile>(new TFile("analysis01.root", "recreate"));

  static boost::posix_time::ptime epoch(boost::gregorian::date(1970, 1, 1));

  /* Search for EPICS PV csv files */
  for (auto &INI: sINI)
    {
      const int it = &INI - &sINI[0];

      boost::filesystem::path pd(INI);

      /* Process good files only */
      if (boost::filesystem::exists(pd) && boost::filesystem::file_size(pd) > 1000 && pd.extension() == ".csv")
        {
          const auto cells = ParseCsv(pd.string());

          std::vector<double> vsec;
          std::vector<double> vval;

          for (const auto& data : cells)
            {
              assert(data.size()==2 && "data's column size should be 2.");
              if (it > 1 && data[1] == "0")
                continue;

              auto *formatd = new boost::posix_time::time_input_facet("%Y/%m/%d %H:%M:%s.%f");
              boost::posix_time::ptime ptimed;

              std::istringstream iss(data[0]);
              iss.imbue(std::locale(iss.getloc(), formatd));
              iss >> ptimed;
              delete formatd;

              auto to_time_t = [](const boost::posix_time::ptime &date)
                               {
                                 return time_t((date - epoch).total_seconds());
                               };

              vsec.push_back(boost::lexical_cast<double>(to_time_t(ptimed) - 788918400)); // ref to 1 Jan, 1995
              vval.push_back(boost::lexical_cast<double>(data[1]));
            }

          std::cout << "Loaded " << pd.stem().string() << std::endl;
          auto gr = boost::shared_ptr<TGraph>(new TGraph(vsec.size(), &vsec[0], &vval[0]));
          gr->SetName(pd.stem().string().c_str());
          gr->SetTitle(pd.stem().string().c_str());
          gr->SetMarkerStyle(20+it);
          gr->SetMarkerColor(1+it);
          gr->SetLineColor(1+it);
          gr->Write();
        }
    }

  fout->Write();
  fout->Close();

  return 0;
}
