#include "MINERvA_CC0pi_XSec_1DQ2_nue.h"

//******************************************************************** 
MINERvA_CC0pi_XSec_1DQ2_nue::MINERvA_CC0pi_XSec_1DQ2_nue(std::string inputfile, FitWeight *rw, std::string  type, std::string fakeDataFile){
//******************************************************************** 

  // Define Measurement
  fName = "MINERvA_CC0pi_XSec_1DQ2_nue";
  fPlotTitles = "; E_{e} (GeV); d#sigma/dE_{e} (cm^{2}/GeV)";
  EnuMin = 0.0;
  EnuMax = 15.0;
  fNormError = 0.101;
  fDefaultTypes = "FIX/FULL";
  fAllowedTypes = "FIX,FREE,SHAPE/DIAG,FULL/NORM/MASK";
  Measurement1D::SetupMeasurement(inputfile, type, rw, fakeDataFile);

  // Setup Data File
  std::string datafile = FitPar::GetDataBase()+"/MINERvA/CC0pi/MINERvA_CC0pi_nue_Data_ARX1509_05729.root";
  std::string dist_name = "";
  
  dist_name = "1DQ2";
  fPlotTitles = "; Q_{QE}^{2} (GeV^{2}); d#sigma/dQ_{QE}^{2} (cm^{2}/GeV^{2})";
  
  SetDataFromFile(datafile, "Data_" + dist_name);
  SetCovarFromDataFile(datafile, "Covar_" + dist_name);

  // Setup Default MC Hists
  SetupDefaultHist();

  // Different generators require slightly different rescaling factors.
  fScaleFactor = (this->fEventHist->Integral("width")*100.0*1E-38/(fNEvents+0.))/this->TotalIntegratedFlux(); 

};

//********************************************************************
void MINERvA_CC0pi_XSec_1DQ2_nue::FillEventVariables(FitEvent *event){
//********************************************************************

  Enu_rec = 0.0;
  
  // Get the relevant signal information
  for (UInt_t j = 2; j < event->Npart(); ++j){

    int PID = (event->PartInfo(j))->fPID;
    if (!event->PartInfo(j)->fIsAlive) continue;

    if (abs(PID) == 11){

      Thetae = ((event->PartInfo(0))->fP.Vect().Angle((event->PartInfo(j))->fP.Vect()));
      Ee = (event->PartInfo(j))->fP.E()/1000.0;

      Enu_rec = FitUtils::EnuQErec((event->PartInfo(j))->fP, cos(Thetae), 34.,true);
      Q2QEe   = FitUtils::Q2QErec((event->PartInfo(j))->fP, cos(Thetae), 34.,true);

      break;
    }
  }
  
  fXVar = Q2QEe;
  LOG(EVT) << "fXVar = "<<fXVar<<std::endl;
  return;
}



//********************************************************************
bool MINERvA_CC0pi_XSec_1DQ2_nue::isSignal(FitEvent *event){
//*******************************************************************

  // Get Nue/NueBar events
  if (fabs(event->PDGnu()) != 12)  return false;
  
  // CC0Pi
  if (!event->IsFS0Pi() || !event->IsCC()) return false;
      
  // restrict energy range
  if (event->Enu()/1000.0 < this->EnuMin || event->Enu()/1000.0 > this->EnuMax) return false;

  // Restrct Ee
  if (Ee < 0.5) return false;
  
  return true;
};

