{
  std::ifstream ifs("Lsp5e+31.dat");

  const unsigned int np = 11;
  auto *g = new TGraph(np);
  g->SetMarkerStyle(24);
  g->SetMarkerSize(2);

  unsigned int i = 0;
  while(!ifs.eof())
    {
      double current, lsp;
      ifs >> current >> lsp;
      g->SetPoint(i, current, lsp);
      i++;
      if (i==np)
	break;
    }

  auto *func4 = new TF1("func", "[0]*exp([1]*x)+[2]", 0.05, 0.9);
  func4->SetLineColor(4);
  func4->SetLineStyle(7);
  g->Fit(func4, "NOME", "", 0, 0.4);

  auto *func5 = new TF1("func", "[0]*exp([1]*x)+[2]", 0.05, 0.9);
  func5->SetLineColor(2);
  func5->SetLineStyle(7);
  g->Fit(func5, "NOME", "", 0, 0.8);
  
  g->SetTitle("");
  g->GetXaxis()->SetTitle("Bunch current product (mA^{2})");
  g->GetXaxis()->CenterTitle(true);
  g->GetYaxis()->SetTitle("Specific luminosity (10^{31} cm^{-2} s^{-1} mA^{-2})");
  g->GetYaxis()->CenterTitle(true);
  g->Draw("ALP");
  //func4->Draw("C,same");
  func5->Draw("C,same");
  auto *axis = g->GetXaxis();
  axis->SetLimits(0., 0.9);
  g->GetHistogram()->SetMinimum(0.);
  g->GetHistogram()->SetMaximum(10.);
}
