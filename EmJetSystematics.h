#include <string>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TStopwatch.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TParameter.h>

#ifndef OUTPUT
#define OUTPUT(x) std::cout<<#x << ": " << x << std::endl
#endif

using std::string;

class EmJetSystematics
{
 public:
  void SetFilenames(string filename_data, string filename_qcd);
  double CalculateIpXYShift();
 private:
  string m_filename_data;
  string m_filename_qcd;
};

void EmJetSystematics::SetFilenames(string filename_data, string filename_qcd)
{
  m_filename_data = filename_data;
  m_filename_qcd = filename_qcd;
}

double EmJetSystematics::CalculateIpXYShift()
{
  TDirectory* currentfile = gDirectory;
  double x_data = TMath::Log10(0.25);
  TFile* file_data = new TFile(m_filename_data.c_str());
  TH1F* hist_data = (TH1F*)file_data->Get("sys_log_track_ipXY");
  int nbins_data = hist_data->GetNbinsX();
  int binx_data = hist_data->FindBin(x_data);
  // Include overflow and underflow bins in integral caculation
  double fraction = hist_data->Integral(0, binx_data) / hist_data->Integral(0, nbins_data+1);
  OUTPUT(fraction);

  TFile* file_qcd = new TFile(m_filename_qcd.c_str());
  TH1F* hist_qcd = (TH1F*)file_qcd->Get("sys_log_track_ipXY");
  double total = hist_qcd->Integral();
  OUTPUT(total);
  double target = total * fraction; // Target integral
  OUTPUT(target);
  int nbins_qcd = hist_qcd->GetNbinsX();
  OUTPUT(nbins_qcd);
  int binx_qcd;
  for (int n=0; n < hist_qcd->GetNbinsX()+1; n++) {
    if (hist_qcd->Integral(0, n) > target) {
      OUTPUT(hist_qcd->Integral(0, n));
      OUTPUT("Breaking");
      binx_qcd = n;
      break;
    }
  }
  double x_qcd = hist_qcd->GetBinLowEdge(binx_qcd);
  OUTPUT(x_data);
  OUTPUT(x_qcd);
  double shift = TMath::Power(10, x_data-x_qcd);

  file_qcd->Close();
  file_data->Close();

  currentfile->cd();

  return shift;
}

// double EmJetSystematics::CalculateIpXYSigShift()
// {
//   TDirectory* currentfile = gDirectory;
//   double x_data = TMath::Log10(0.25);
//   TFile* file_data = new TFile(m_filename_data.c_str());
//   TH1F* hist_data = (TH1F*)file_data->Get("sys_log_track_ipXYSig");
//   int nbins_data = hist_data->GetNbinsX();
//   int binx_data = hist_data->FindBin(x_data);
//   // Include overflow and underflow bins in integral caculation
//   double fraction = hist_data->Integral(0, binx_data) / hist_data->Integral(0, nbins_data+1);
//   OUTPUT(fraction);

//   TFile* file_qcd = new TFile(m_filename_qcd.c_str());
//   TH1F* hist_qcd = (TH1F*)file_qcd->Get("sys_log_track_ipXYSig");
//   double total = hist_qcd->Integral();
//   OUTPUT(total);
//   double target = total * fraction; // Target integral
//   OUTPUT(target);
//   int nbins_qcd = hist_qcd->GetNbinsX();
//   OUTPUT(nbins_qcd);
//   int binx_qcd;
//   for (int n=0; n < hist_qcd->GetNbinsX()+1; n++) {
//     if (hist_qcd->Integral(0, n) > target) {
//       OUTPUT(hist_qcd->Integral(0, n));
//       OUTPUT("Breaking");
//       binx_qcd = n;
//       break;
//     }
//   }
//   double x_qcd = hist_qcd->GetBinLowEdge(binx_qcd);
//   OUTPUT(x_data);
//   OUTPUT(x_qcd);
//   double shift = TMath::Power(10, x_data-x_qcd);

//   file_qcd->Close();
//   file_data->Close();

//   currentfile->cd();

//   return shift;
// }
