#include "TGraph.h"
#include <vector>
#include "Riostream.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TCanvas.h"

void plot(int n=10)
{
	std::cout << std::string(80,'*') << std::endl;
	std::cout << "* Results for " << n << " runs" << std::endl;
	std::cout << std::string(80,'*') << std::endl;

	TGraph* gv[n];

	for ( int i = 0; i < n; ++i )
	{
		gv[i] = new TGraph(Form("result.%d",i));
	}

	TGraphErrors* g = new TGraphErrors(gv[0]->GetN());

	for ( int i = 0; i < g->GetN(); ++i )
	{
		vector<double> vy;
		
		for ( int j = 0; j < n; ++j ) 
		{
			vy.push_back(gv[j]->GetY()[i]);
		}

		Double_t x = gv[0]->GetX()[i];

		Double_t y = TMath::Mean(vy.size(),&vy[0]);
		Double_t ey = TMath::RMS(vy.size(),&vy[0]);

		std::cout << Form("%2d parallel processes : mean RootMarks %7.2f",TMath::Nint(x),y);

		if ( n > 1 )
		{
			std::cout << Form("rms %7.2f ",ey);
		}
		std::cout << endl;
		g->SetPoint(i,x,y);
		g->SetPointError(i,0,ey);
	}
	
	std::cout << std::string(80,'*') << std::endl;

	g->SetLineColor(1);
	g->SetMarkerStyle(20);

	g->SetTitle("SubRootMarks vs # of proc");

	g->GetXaxis()->SetTitle("Number of stress processes");

	g->GetYaxis()->SetTitle("SubRootMarks");
	g->GetYaxis()->SetTitleOffset(1.2);

	TCanvas* c =  new TCanvas("SubRootMarksVsNProc","SubRootMarksVsNProc");

	g->Draw("ALEP");

	c->SaveAs("SubRootMarksVsNProc.pdf");


}