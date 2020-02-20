// Copyright 2016 L. Pickering, P Stowell, R. Terri, C. Wilkinson, C. Wret

/*******************************************************************************
*    This file is part of NUISANCE.
*
*    NUISANCE is free software: you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    NUISANCE is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with NUISANCE.  If not, see <http://www.gnu.org/licenses/>.
*******************************************************************************/

#include "MicroBooNE_CC1MuNp_XSec_1Dthetamup_nu.h"
#include "MicroBooNE_SignalDef.h"

//********************************************************************
MicroBooNE_CC1MuNp_XSec_1Dthetamup_nu::MicroBooNE_CC1MuNp_XSec_1Dthetamup_nu(nuiskey samplekey) {
//********************************************************************

  // Sample overview ---------------------------------------------------
  std::string descrip = "MicroBooNE_CC1MuNp_XSec_1Dthetamup_nu sample. \n" \
                        "Target: Ar \n" \
                        "Flux: BNB FHC numu \n" \
                        "Signal: CC 1 mu + Np \n";

  // Setup common settings
  fSettings = LoadSampleSettings(samplekey);
  fSettings.SetDescription(descrip);
  fSettings.SetXTitle("#theta_{p#mu}^{reco} (radians)");
  fSettings.SetYTitle("d#sigma/d#theta_{p#mu}^{reco} (cm^{2}/radian/nucleon)");
  fSettings.SetAllowedTypes("FULL,DIAG/FREE,SHAPE,FIX/SYSTCOV/STATCOV","FIX/FULL");
  fSettings.SetEnuRange(0.0, 6.8);
  fSettings.DefineAllowedTargets("Ar");

  // Plot information
  fSettings.SetTitle("MicroBooNE_CC1MuNp_XSec_1Dthetamup_nu");
  fSettings.DefineAllowedSpecies("numu");

  FinaliseSampleSettings();

  // Scaling Setup ---------------------------------------------------
  // ScaleFactor automatically setup for DiffXSec/cm2/Nucleon
  fScaleFactor = GetEventHistogram()->Integral("width") / fNEvents * 1E-38 / TotalIntegratedFlux();

  // Plot Setup -------------------------------------------------------
  std::string inputFile = FitPar::GetDataBase() + "/MicroBooNE/CC1MuNp/CCNp_data_MC_cov_dataRelease.root";
  SetDataFromRootFile(inputFile, "DataXsec_thetamup");
  ScaleData(1E-38);
  SetCovarFromRootFile(inputFile, "CovarianceMatrix_thetamup");

  // Load smearing matrix ---------------------------------------------
  TFile* inputRootFile = TFile::Open(inputFile.c_str());
  assert(inputRootFile && inputRootFile->IsOpen());
  TH2D* tempsmear = (TH2D*) inputRootFile->Get("SmearingMatrix_thetamup");
  assert(tempsmear);

  // Normalize columns
  TH1D* hpx = tempsmear->ProjectionX("_smearing_px");
  for (int i=1; i<tempsmear->GetNbinsX()+1; i++) {
    for (int j=1; j<tempsmear->GetNbinsY()+1; j++) {
      double v = tempsmear->GetBinContent(i, j) / hpx->GetBinContent(i);
      tempsmear->SetBinContent(i, j, v);
    }
  }
  delete hpx;

  fSmearingMatrix = new TMatrixDSym(fDataHist->GetNbinsX());
  assert(fSmearingMatrix);
  for (int i=0; i<fDataHist->GetNbinsX(); i++) {
    for (int j=0; j<fDataHist->GetNbinsX(); j++) {
      (*fSmearingMatrix)(i,j) = tempsmear->GetBinContent(i+1, j+1);
    }
  }
  inputRootFile->Close();

  // Final setup  -----------------------------------------------------
  FinaliseMeasurement();
};


bool MicroBooNE_CC1MuNp_XSec_1Dthetamup_nu::isSignal(FitEvent* event) {
  return SignalDef::isCC1MuNp(event, EnuMin, EnuMax);
};


void MicroBooNE_CC1MuNp_XSec_1Dthetamup_nu::FillEventVariables(FitEvent* event) {
  if (event->NumFSParticle(13) == 0) return;
  if (event->NumFSParticle(2212) == 0) return;
  TVector3 vpmu = event->GetHMFSParticle(13)->fP.Vect();
  TVector3 vppl = event->GetHMFSParticle(2212)->fP.Vect();
  fXVar = vpmu.Angle(vppl);
};


void MicroBooNE_CC1MuNp_XSec_1Dthetamup_nu::ConvertEventRates() {
  // Do standard conversion
  Measurement1D::ConvertEventRates();

  // Apply MC truth -> reco smearing
  TH1D* truth = (TH1D*) fMCHist->Clone(TString(fMCHist->GetName()) + "_truth");

  for (int ireco=1; ireco<fMCHist->GetNbinsX()+1; ireco++) {
    double total = 0;
    for (int itrue=1; itrue<fMCHist->GetNbinsX()+1; itrue++) {
      total += truth->GetBinContent(itrue) * truth->GetBinWidth(itrue) * fSmearingMatrix->operator()(ireco-1, itrue-1);
    }
    fMCHist->SetBinContent(ireco, total / fMCHist->GetBinWidth(ireco));
  }
}
