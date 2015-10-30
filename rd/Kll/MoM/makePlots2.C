Double_t getMax(Int_t n, Double_t* vals) {
	Double_t max(-999.);

	for(Int_t i=0; i<n; ++i) {
		if(vals[i] > max) max=vals[i];
	}

	return max;
}

Double_t getMin(Int_t n, Double_t* vals) {
	Double_t min(999.);

	for(Int_t i=0; i<n; ++i) {
		if(vals[i] < min) min=vals[i];
	}

	return min;
}

void makePlots2() {
	gROOT->ProcessLine(".x lhcbStyle.C");
	gStyle->SetOptStat(0);

	Double_t ignore, eff, massfit1, massfit2, bias1, bias2;

	Double_t value[6][19], stat[6][19], syst[6][19];
	Int_t ng(4);

	TH1D pull("pull","",15,-3.0,3.0);

	std::ifstream fin("errs.dat");
	fin.ignore(256,'\n');
	for(Int_t q=0; q<19; ++q) {
		for(Int_t g=1; g<ng+1; ++g) {	
			fin >> ignore >> ignore >> value[g+1][q] >> stat[g+1][q] >> eff >> massfit1 >> massfit2 >> bias1 >> bias2;
			syst[g+1][q] = TMath::Sqrt(stat[g+1][q]*stat[g+1][q] + eff*eff + massfit1*massfit1 + massfit2*massfit2 + bias1*bias1 + bias2*bias2);
			//plot syst errors on top of stats

			if(q<17) pull.Fill(value[g+1][q]/syst[g+1][q]);
		}
		value[0][q] = value[2][q];
		stat[0][q]  = stat[2][q];
		syst[0][q]  = syst[2][q];

		value[1][q] = 2*value[3][q]+1;
		stat[1][q]  = 2*stat[3][q];
		syst[1][q]  = 2*syst[3][q];
	}

	fin.close();

	Double_t  xs[19] = { 0.54, 1.55, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 11.375, 12.125, 15.5, 16.5, 17.5, 18.5, 19.5, 20.5, 21.5, 3.55, 18.5};
	Double_t exs[19] = { 0.44, 0.45, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,  0.375,  0.375,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5, 2.45,  3.5};

	TCanvas c;

	TF1 gaus("gaus","gaus(0)",-3.,3.);
	gaus.SetParameter(0,68);
	gaus.SetParameter(1,pull.GetBinCenter(pull.GetMaximumBin()));
	gaus.SetParameter(2,pull.GetStdDev());
	gaus.SetParLimits(1,-2.0,2.0);
	gaus.SetParLimits(2,0.3,3.0);
	pull.Fit(&gaus,"S");
	pull.Draw();
	c.SaveAs("gPulls.pdf");

	for(Int_t i=0; i<ng+2; ++i) {
            TString varName("");
            TString varName2("");

	    switch(i) {
		    case 0:
			    varName="A_{FB}";
			    varName2="AFB";
			    break;
		    case 1:
			    varName="F_{H}";
			    varName2="FH";
			    break;
		    default:
			    varName="G^{("; varName+=(i-1); varName+=")}";
			    varName2="G"; varName2+=(i-1);
	    }
	    
	    Double_t yMin = getMin(19,value[i])-getMax(19,syst[i]);
	    Double_t yMax = getMax(19,value[i])+getMax(19,syst[i]);

	    if(yMax-yMin == 0) continue;

	    TH2D h("h","",1,0.0,22.,1,yMin,yMax);

	    TLine l(0.0,0.0,22.0,0.0);
	    l.SetLineStyle(kDashed);
	    l.SetLineColor(kRed);
	    
	    h.GetYaxis()->SetTitle(varName);
	    h.GetXaxis()->SetTitle("q^{2} [GeV^{2}/c^{4}]");

	    TGraphErrors g1(17,xs,   value[i],   exs,   stat[i]   );
	    TGraphErrors g2( 2,xs+17,value[i]+17,exs+17,stat[i]+17);
	    TGraphErrors g3(17,xs,   value[i],   exs,   syst[i]   );
	    TGraphErrors g4( 2,xs+17,value[i]+17,exs+17,syst[i]+17);

	    g2.SetLineColor(kBlue);
	    g2.SetMarkerColor(kBlue);
	    g3.SetLineColor(kRed);
	    g4.SetLineColor(kRed);
	    g2.SetFillColor(kBlue);
	    g4.SetFillColor(kRed);
	    g2.SetFillStyle(3001);
	    g4.SetFillStyle(3002);

	    h.Draw();
	    l.Draw();
	    g4.Draw("E2same");
	    g2.Draw("E2same");
	    g3.Draw("PEsame");
	    g1.Draw("PEsame");
	    c.SaveAs("plots/fromPatrickNew/"+varName2+".pdf");
	    c.SaveAs("plots/fromPatrickNew/"+varName2+".png");
	}

}
