#include "TH1F.h"
#include "TH3.h"
#include "TFile.h"
#include "TRandom1.h"
#include "TRandom2.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

class  TriggerReWeighting {
 public:

   TriggerReWeighting ( ) { } ;

 TriggerReWeighting( std::string generatedFile,
                     // std::string dataFile,
                     std::string GenHistName,
                     std::string DataHistName) :
   generatedFileName_( generatedFile),
     // dataFileName_     ( dataFile ),
     GenHistName_      ( GenHistName ),
     DataHistName_     ( DataHistName )
     {
       generatedFile_ = new TFile( generatedFileName_.c_str() ) ; //MC distribution
       // dataFile_      = new TFile( dataFileName_.c_str() );       //Data distribution

       // Data_distr_ = new TH1F(  *(static_cast<TH1F*>(dataFile_->Get( DataHistName_.c_str() )->Clone() )) );
       Data_distr_ = new TH1F(  *(static_cast<TH1F*>(generatedFile_->Get( DataHistName_.c_str() )->Clone() )) );
       auto MC_distr_temp = new TH1F(  *(static_cast<TH1F*>(generatedFile_->Get( GenHistName_.c_str() )->Clone() )) );

       MC_distr_ = MC_distr_temp;

       // normalize both histograms first
       // Data_distr_->Scale( 1.0/ Data_distr_->Integral() );
       // MC_distr_->Scale( 1.0/ MC_distr_->Integral() );

       // Data_distr_->Print();
       // std::cout << "Dumping data dist" << std::endl;
       // for(int ibin = 1; ibin<Data_distr_->GetNbinsX()+1; ++ibin){
       //   std::cout << "   " << ibin-1 << ", " << Data_distr_->GetBinContent(ibin) << std::endl;
       // }
       // MC_distr_->Print();
       // std::cout << "Dumping MC dist" << std::endl;
       // for(int ibin = 1; ibin<MC_distr_->GetNbinsX()+1; ++ibin){
       //   std::cout << "   " << ibin-1 << ", " << MC_distr_->GetBinContent(ibin) << std::endl;
       // }

       weights_ = new TH1F( *(Data_distr_)) ;

       // MC * data/MC = data, so the weights are data/MC:

       weights_->SetName("triggerWeights");

       TH1F* den = new TH1F(*(MC_distr_));

       // weights_->Divide( den );  // so now the average weight should be 1.0

       int NBins = weights_->GetNbinsX();
       for(int ibin = 1; ibin<NBins+1; ++ibin){
         double ratio = Data_distr_->GetBinContent(ibin) / MC_distr_->GetBinContent(ibin);
         weights_->SetBinContent(ibin, ratio);
       }

       // std::cout << "Data_distr_->Integral(): " << Data_distr_->Integral() << std::endl;
       // std::cout << "MC_distr_->Integral(): " << MC_distr_->Integral() << std::endl;
       // std::cout << " Lumi/Pileup Reweighting: Computed Weights per In-Time Nint " << std::endl;


       for(int ibin = 1; ibin<NBins+1; ++ibin){
         // std::cout << "   " << ibin-1 << "\t" << weights_->GetBinContent(ibin) << "\t" << weights_->GetBinContent(ibin) * MC_distr_->GetBinContent(ibin)<< std::endl;
       }

     }

   double weight( float n_int ){
     int bin = weights_->GetXaxis()->FindBin( n_int );
     return weights_->GetBinContent( bin );
   }

 protected:

   std::string generatedFileName_;
   std::string dataFileName_;
   std::string GenHistName_;
   std::string DataHistName_;
   TFile *generatedFile_;
   TFile *dataFile_;
   TH1F  *weights_;

   //keep copies of normalized distributions:
   TH1F*      MC_distr_;
   TH1F*      Data_distr_;

};
