#include <TMathBase.h>
#include <TString.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TLine.h>
#include <TEllipse.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TFile.h>
#include <TMultiGraph.h>
#include <cassert>
#include <cmath>
#include <ctime>
#include <math.h>
#include <iostream>
#include <TFile.h>

using namespace std;

const int n1 = 9;
const int n2 = 113;
const int resolution = 1000;

void SaveAndGenerateProjections(TFile*file,TH2*th2);
void MovingSum(TH1*th1);
	
void LinearEventGenerator(double tracker_cells[][n2], double a, double b);
void CircularEventGenerator(double tracker_cells[][n2], double x0, double y0, double r0);
void CircularEventGeneratorFromFoil(double tracker_cells[][n2], double y0, double sinus, double r0);

int legendre()
{
	double tracker_cells[n1][n2] = {0.0};
	// LinearEventGenerator( tracker_cells, -1.0, 95.2);
	// LinearEventGenerator( tracker_cells, -1.18, 43.0);
	// LinearEventGenerator( tracker_cells, -16328, 602000);
	// LinearEventGenerator( tracker_cells, 6.43, 18.5);
	//LinearEventGenerator( tracker_cells, -0.97, 40.0);
	//LinearEventGenerator( tracker_cells, -0.97, 25.0);
	// CircularEventGenerator( tracker_cells, 40.46, -13.45, 36.73);
	CircularEventGeneratorFromFoil(tracker_cells, 56.,-0.8, +20.);
		
	TH2F *tracker_hits = new TH2F("tracker", "tracker; x; y", n1, -0.5, n1-0.5, n2, -0.5, n2-0.5);
	for(int i = 0; i < n1; i++)
		for(int j = 0; j < n2; j++)
		{
			tracker_hits->Fill(i, j, tracker_cells[i][j]);
		}

	clock_t start = clock();
	vector<double> hits;	
	for(int i = 0; i < n1; i++)
		for(int j = 0; j < n2; j++)
			if(tracker_cells[i][j] != 0.0)
			{
				hits.push_back(i);
				hits.push_back(j);
				hits.push_back(tracker_cells[i][j]);
			}
			
	
	TCanvas *can1 = new TCanvas("can1", "canvas", 1000, 1000);
	can1->Range(-2,-2,n1+2,n2+2);
	TEllipse *el[hits.size()];
	for(int i = 0; i < hits.size()/3; i++)
	{
		el[i] = new TEllipse(hits.at(3*i), hits.at((3*i) + 1), hits.at((3*i)+2), hits.at((3*i)+2));
		el[i]->SetLineWidth(2);
   		el[i]->Draw("SAME");
	}
	// can1->SaveAs("lines_without_rec.png");
	
	cout << "number of hits: " << hits.size()/3 << endl;		

	vector<double> reconstructed_lines;
	int hit_count_last = 0;


	TFile*file = new TFile("data.root","RECREATE");

	
	//while( hits.size()/2 > 5 && hits.size() != hit_count_last)
	
	{
		
		double r, theta;
		double phi1 = 0.0;
		double phi2 = 2*M_PI;
		double R1 = -80;
		double R2 = 80;
		double gaussian_theta = 2;
		double gaussian_r = 2;
		
		double peak_Theta;
		double peak_R;
		double delta_phi, delta_R;
		for(int q = 0; q < 4; q++)
		{	
			double offset = (phi2 - phi1)/(2.0*resolution);
			TH2F *sinograms = new TH2F(Form("legender_track_%lu_iter_%d.png", reconstructed_lines.size()/2, q), "sinograms; theta; r", resolution, phi1+offset, phi2+offset, resolution, R1, R2);		
			for(int i = 0; i < hits.size()/3; i++)
			{
				for(int k = 0; k < resolution; k++)
				{
					theta = phi1 + (k * (phi2 - phi1) / resolution);
					r = ( -(hits.at(3*i))*cos(theta) ) - ( -(hits.at(3*i+1))*sin(theta) );
					for(size_t j = 0; j<400; j++) {
						sinograms->Fill( gRandom->Gaus(theta,gaussian_theta*(phi2- phi1)/resolution),gRandom->Gaus(r - hits.at(3*i+2),gaussian_r*(R2-R1)/resolution));
						sinograms->Fill( gRandom->Gaus(theta,gaussian_theta*(phi2- phi1)/resolution),gRandom->Gaus(r + hits.at(3*i+2),gaussian_r*(R2-R1)/resolution));
					}
				}	
			}
				
			//if(hit_count_last == 0 && q == 0) sinograms->Draw("COLZ");
			
			sinograms->SetStats(0);
			sinograms->Draw("COLZ");
			// can1->SetLogz();
			// can1->SaveAs(Form("legender_track-%lu_iter-%d.png", reconstructed_lines.size()/2, q));
			SaveAndGenerateProjections(file,sinograms);
			
			double maximum = 0.0;
			for(int i = 1; i < resolution; i++)
			{
				for(int j = 1; j < resolution; j++)
				{
					if(maximum < sinograms->GetBinContent(i,j))
					{
						maximum = sinograms->GetBinContent(i,j);
						peak_Theta = i;
						peak_R = j;
					}
				}
			}
			delete sinograms;
			
			delta_phi = phi2 - phi1;
			delta_R = R2 - R1;
			
			peak_Theta = phi1 + delta_phi * peak_Theta / (double)resolution;
			peak_R = R1 + delta_R * (peak_R-0.5) / (double)resolution;
			
			//cout << "iteration " << q  << ": " << peak_Theta << "	" << peak_R << endl;
			
			phi1 = peak_Theta - 0.05*delta_phi;
			phi2 = peak_Theta + 0.05*delta_phi;
			R1 = peak_R - 0.05*delta_R;
			R2 = peak_R + 0.05*delta_R;
		
		}
		cout << "peak found:" << peak_Theta << "	" << peak_R << endl;
		
		double a = 1.0 / tan(peak_Theta);
		double b = peak_R / sin(peak_Theta);
		
		
		double distance;
		double denominator = sqrt((a*a) + 1);	
		
		hit_count_last = hits.size();
		int iter = 0;
		int hit_counter = 0;
		while(iter < hits.size()/3)
		{
			distance = abs(hits.at(3*iter+1) - (a*hits.at(3*iter)) - b) / denominator;
			if( abs(distance - hits.at(3*iter+2)) < 0.05 )
			{
				hits.erase(hits.begin() + 3*iter, hits.begin() + 3*iter+3);
				hit_counter++;
			}	
			else
				iter++;
		}
		if(hit_counter > 3)
		{
			reconstructed_lines.push_back(a);
			reconstructed_lines.push_back(b);
		}	
		//cout << hits.size()/3 << " hits remaining" << endl << endl;
	}
	
	std::cout<<"Finished after " << (clock() - start) / (double) CLOCKS_PER_SEC << " seconds.";
	
	for(int i = 0; i < reconstructed_lines.size()/2; i++)
	{
		std::cout << "line: " << reconstructed_lines.at(2*i) << " " << reconstructed_lines.at(2*i+1) << endl; 
		
		TLine *line = new TLine(-1.0, -1.0*reconstructed_lines.at(2*i) + reconstructed_lines.at(2*i+1), 51.0, reconstructed_lines.at(2*i)*51.0 + reconstructed_lines.at(2*i+1));
		line->SetLineColor(kRed);
		line->Draw("Same");
		line->SetLineWidth(2);
		//can1->SaveAs(Form("lines_with_rec_%d.png", i));
	}
	
	return EXIT_SUCCESS;
}
void SaveAndGenerateProjections(TFile*file,TH2*th2)
{
  size_t nBinsX = th2->GetNbinsX();
  size_t nBinsY = th2->GetNbinsY();
  double xMin = th2->GetXaxis()->GetXmin();
  double xMax = th2->GetXaxis()->GetXmax();
  double yMin = th2->GetYaxis()->GetXmin();
  double yMax = th2->GetYaxis()->GetXmax();

  TH1D* projX = new TH1D(Form("%s_x",th2->GetName()), "Projection on X axis (Max)", nBinsX, xMin, xMax);
  TH1D* projY = new TH1D(Form("%s_y",th2->GetName()), "Projection on Y axis (Max)", nBinsY, yMin, yMax);

  for (size_t ix = 1; ix <= nBinsX; ++ix) {
    double maxValX = 0;
    for (int iy = 1; iy <= nBinsY; ++iy) {
      double binContent = th2->GetBinContent(ix, iy);
      if (binContent > maxValX) {
          maxValX = binContent;
      }
    }
    projX->SetBinContent(ix, maxValX);
  }

  for (size_t iy = 1; iy <= nBinsY; ++iy) {
    double maxValY = 0;
    for (int ix = 1; ix <= nBinsX; ++ix) {
      double binContent = th2->GetBinContent(ix, iy);
      if (binContent > maxValY) {
        maxValY = binContent;
      }
    }
    projY->SetBinContent(iy, maxValY);
  }

  file->cd();
  th2->Write();
  projX->Write();
  projY->Write();
}
void LinearEventGenerator( double tracker_cells[][n2], double a, double b)
{
	double distance;
	double denominator = sqrt((a*a) + 1);
	
	for(int i = 0; i < n1; i++)
		for(int j = 0; j < n2; j++)
		{
			distance = abs(j - (a*i) - b) / denominator;
			if( tracker_cells[i][j] == 0.0 && distance < 0.5 )
				tracker_cells[i][j] = distance;
		}
}
void CircularEventGenerator( double tracker_cells[][n2], double x0, double y0, double r0)
{
	double distance;
	double deltaxy;
	
	for(int i = 0; i < n1; i++)
		for(int j = 0; j < n2; j++)
		{
			deltaxy = (i-x0)*(i-x0) + (j-y0)*(j-y0);
			distance = abs(r0 - sqrt(deltaxy));
			if( tracker_cells[i][j] == 0.0 && distance < 0.5 )
				tracker_cells[i][j] = distance;
		}
}
void CircularEventGeneratorFromFoil(double tracker_cells[][n2], double y0, double sinus, double r0)
{
	assert(abs(sinus)<=1.);
	CircularEventGenerator(tracker_cells, sinus*r0, y0+sqrt(1-sinus*sinus)*r0,TMath::Abs(r0));
}
/// @brief returns number of bin in the middle of window of len
int MovingSum(TH1*th1,size_t len)
{
  int nBins = th1->GetNbinsX();

	double max = -1000.;
	int max_i = 1;
  for (int i = 1+len; i <= nBins; ++i) {
    int count = 0;
    double sum = 0.0;

    for (int j = i - len; j <= i; ++j) {
      if (j >= 1 && j <= nBins) {
        sum += th1->GetBinContent(j);
        ++count;
      }
  	}

  	if (sum>max) {
  		max = sum;
  		max_i = i;
  	}
  }

  return max_i-len/2;
}
