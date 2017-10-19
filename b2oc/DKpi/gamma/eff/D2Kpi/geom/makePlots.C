{
	gROOT->ProcessLine(".L lhcbStyle.C");
	lhcbStyle();
	gStyle->SetPalette(1,0);
	const Int_t NRGBs = 5;
	const Int_t NCont = 255;
	Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00};
	Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51};
	Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00};
	Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00};
	TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
	gStyle->SetNumberContours(NCont);
	gStyle->SetOptStat(0000);

	TFile* f = TFile::Open("hists/geom_sq15_all.root");
	TH2D* eff = (TH2D*)f->Get("efficiency");
	TH2D* errp = (TH2D*)f->Get("errorHi");
	TH2D* errm = (TH2D*)f->Get("errorLo");
	TH2D* errppc = errp->Clone("errorHiPc");
	TH2D* errmpc = errm->Clone("errorLoPc");

	errppc->Scale(100.);
	errppc->Divide(eff);
	errmpc->Scale(100.);
	errmpc->Divide(eff);

	TCanvas * canvas = new TCanvas("c2","c2");
	eff->SetLabelFont(62,"x");
	eff->SetLabelFont(62,"y");
	eff->SetTitleFont(62,"x");
	eff->SetTitleFont(62,"y");
	eff->SetTitleSize(0.06,"x");
	eff->SetTitleSize(0.06,"y");
	eff->SetLabelSize(0.05,"x");
	eff->SetLabelSize(0.05,"y");
	eff->SetXTitle("m'");
	eff->SetYTitle("#theta'");
	eff->GetZaxis()->SetRangeUser(0.31,0.37);
	canvas->SetRightMargin(0.17);
	eff->Draw("colztext30");
	canvas->SaveAs("../latex/figs/geom_eff.pdf");
	eff->Draw("colz");
	canvas->SaveAs("../latex/figs/geom_eff_nolabel.pdf");

	gStyle->SetPaintTextFormat("0.2f %%");

	errp->SetLabelFont(62,"x");
	errp->SetLabelFont(62,"y");
	errp->SetTitleFont(62,"x");
	errp->SetTitleFont(62,"y");
	errp->SetTitleSize(0.06,"x");
	errp->SetTitleSize(0.06,"y");
	errp->SetLabelSize(0.05,"x");
	errp->SetLabelSize(0.05,"y");
	errp->SetXTitle("m'");
	errp->SetYTitle("#theta'");
	canvas->SetRightMargin(0.17);
	errp->Draw("colz");
	errppc->Draw("text30same");
	canvas->SaveAs("../latex/figs/geom_errp.pdf");

	errm->SetLabelFont(62,"x");
	errm->SetLabelFont(62,"y");
	errm->SetTitleFont(62,"x");
	errm->SetTitleFont(62,"y");
	errm->SetTitleSize(0.06,"x");
	errm->SetTitleSize(0.06,"y");
	errm->SetLabelSize(0.05,"x");
	errm->SetLabelSize(0.05,"y");
	errm->SetXTitle("m'");
	errm->SetYTitle("#theta'");
	canvas->SetRightMargin(0.17);
	errm->Draw("colz");
	errmpc->Draw("text30same");
	canvas->SaveAs("../latex/figs/geom_errm.pdf");

	errppc->Scale(0.01);
	errmpc->Scale(0.01);

	errppc->SetLabelFont(62,"x");
	errppc->SetLabelFont(62,"y");
	errppc->SetTitleFont(62,"x");
	errppc->SetTitleFont(62,"y");
	errppc->SetTitleSize(0.06,"x");
	errppc->SetTitleSize(0.06,"y");
	errppc->SetLabelSize(0.05,"x");
	errppc->SetLabelSize(0.05,"y");
	errppc->SetXTitle("m'");
	errppc->SetYTitle("#theta'");
	errppc->GetZaxis()->SetRangeUser(0.008,0.031);
	canvas->SetRightMargin(0.17);
	errppc->Draw("colztext30");
	canvas->SaveAs("../latex/figs/geom_errp_frac.pdf");
	errppc->Draw("colz");
	canvas->SaveAs("../latex/figs/geom_errp_frac_nolabel.pdf");

	errmpc->SetLabelFont(62,"x");
	errmpc->SetLabelFont(62,"y");
	errmpc->SetTitleFont(62,"x");
	errmpc->SetTitleFont(62,"y");
	errmpc->SetTitleSize(0.06,"x");
	errmpc->SetTitleSize(0.06,"y");
	errmpc->SetLabelSize(0.05,"x");
	errmpc->SetLabelSize(0.05,"y");
	errmpc->SetXTitle("m'");
	errmpc->SetYTitle("#theta'");
	errmpc->GetZaxis()->SetRangeUser(0.008,0.031);
	canvas->SetRightMargin(0.17);
	errmpc->Draw("colztext30");
	canvas->SaveAs("../latex/figs/geom_errm_frac.pdf");
	errmpc->Draw("colz");
	canvas->SaveAs("../latex/figs/geom_errm_frac_nolabel.pdf");

   std::cout << errppc->GetBinContent(errppc->GetMinimumBin()) << "-" << errppc->GetBinContent(errppc->GetMaximumBin()) << std::endl;
   std::cout << errmpc->GetBinContent(errmpc->GetMinimumBin()) << "-" << errmpc->GetBinContent(errmpc->GetMaximumBin()) << std::endl;
}
