#include <string>
#include <cassert>
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
const double cut_3dSig = 8.0;

using std::string;


class EmJetSystematics
{
 public:
  EmJetSystematics();
  void SetFilenames(string filename_data, string filename_qcd);
  double CalculateIpXYShift();
  double GetShiftedIpXY(double ipXY_in);
  double Calculate3dSigShift();
  double GetShifted3dSig(double in_3dSig);
  void SetTesting(int testing) {m_testing = testing;};
  void SetDirectionModeling(int direction) {m_direction_modeling = direction;};
  void SetDirectionPdf(int direction) {m_directions_pdf = direction;};
 private:
  int m_testing = 0; // Set to non-zero to enable testing

  // Directions for systematic shifts
  //  0 : no shift
  //  1 : positive shift
  // -1 : negative shift
  int m_directions_pdf;
  int m_direction_modeling;

  int m_pdf_ready = 0; // Only calculate pdf shifts if not zero
  string m_filename_data;
  string m_filename_qcd;
  double m_ipXY_shift;
  double m_3dSig_shift;
};

EmJetSystematics::EmJetSystematics()
{
  m_directions_pdf      = 0;
  m_direction_modeling = 0;
}

void EmJetSystematics::SetFilenames(string filename_data, string filename_qcd)
{
  m_filename_data = filename_data;
  m_filename_qcd = filename_qcd;
  m_ipXY_shift = CalculateIpXYShift();
  m_3dSig_shift = Calculate3dSigShift();
  m_pdf_ready = 1;
}

double EmJetSystematics::CalculateIpXYShift()
{
  OUTPUT("================================\n");
  OUTPUT("CalculateIpXYShift\n");
  OUTPUT("================================\n");
  TDirectory* currentfile = gDirectory;
  double x_data = TMath::Log10(0.25);
  TFile* file_data = new TFile(m_filename_data.c_str());
  TH1F* hist_data = (TH1F*)file_data->Get("sys_log_track_ipXY");
  int nbins_data = hist_data->GetNbinsX();
  int binx_data = hist_data->FindBin(x_data);
  // Include overflow and underflow bins in integral caculation
  // WARNING: Integral calculation does not work well if overflow/underflow bins have high number of entries
  // Bin histograms so that number of entries in overflow/underflow bins are low
  double fraction = hist_data->Integral(0, binx_data) / hist_data->Integral(0, nbins_data+1);
  OUTPUT(fraction);

  TFile* file_qcd = new TFile(m_filename_qcd.c_str());
  TH1F* hist_qcd;
  if (m_testing==0) {
    hist_qcd = (TH1F*)file_qcd->Get("sys_log_track_ipXY");
  }
  else {
    hist_qcd = (TH1F*)file_qcd->Get("systest_log_track_ipXY");
  }
  int nbins_qcd = hist_qcd->GetNbinsX();
  OUTPUT(nbins_qcd);
  double total = hist_qcd->Integral(0, nbins_qcd+1);
  OUTPUT(total);
  double target = total * fraction; // Target integral
  OUTPUT(target);
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

double EmJetSystematics::GetShiftedIpXY(double ipXY_in)
{
  int direction = m_direction_modeling;
  assert(direction==0 || direction ==1 || direction ==-1);
  assert(m_pdf_ready);
  double ipXY_out = ipXY_in * TMath::Power(m_ipXY_shift, m_direction_modeling);
  return ipXY_out;
}

double EmJetSystematics::Calculate3dSigShift()
{
  OUTPUT("================================\n");
  OUTPUT("Calculate3dSigShift\n");
  OUTPUT("================================\n");
  TDirectory* currentfile = gDirectory;
  double x_data = cut_3dSig;
  TFile* file_data = new TFile(m_filename_data.c_str());
  TH1F* hist_data = (TH1F*)file_data->Get("sys_track_3dSig");
  int nbins_data = hist_data->GetNbinsX();
  int binx_data = hist_data->FindBin(x_data);
  // Include overflow and underflow bins in integral caculation
  // WARNING: Integral calculation does not work well if overflow/underflow bins have high number of entries
  // Bin histograms so that number of entries in overflow/underflow bins are low
  double fraction = hist_data->Integral(0, binx_data) / hist_data->Integral(0, nbins_data+1);
  OUTPUT(fraction);

  TFile* file_qcd = new TFile(m_filename_qcd.c_str());
  TH1F* hist_qcd;
  if (m_testing==0) {
    hist_qcd = (TH1F*)file_qcd->Get("sys_track_3dSig");
  }
  else {
    hist_qcd = (TH1F*)file_qcd->Get("systest_track_3dSig");
  }
  int nbins_qcd = hist_qcd->GetNbinsX();
  OUTPUT(nbins_qcd);
  double total = hist_qcd->Integral(0, nbins_qcd+1);
  OUTPUT(total);
  double target = total * fraction; // Target integral
  OUTPUT(target);
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
  double shift = x_data / x_qcd;

  file_qcd->Close();
  file_data->Close();
  currentfile->cd();
  return shift;
}

double EmJetSystematics::GetShifted3dSig(double in_3dSig)
{
  int direction = m_direction_modeling;
  assert(direction==0 || direction ==1 || direction ==-1);
  assert(m_pdf_ready);
  double out_3dSig = in_3dSig * TMath::Power(m_3dSig_shift, m_direction_modeling);
  return out_3dSig;
}
