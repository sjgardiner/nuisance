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

#include "MicroBooNE_CC1MuNp_XSec_1Dcosp_nu.h"
#include "MicroBooNE_SignalDef.h"

//********************************************************************
MicroBooNE_CC1MuNp_XSec_1Dcosp_nu::MicroBooNE_CC1MuNp_XSec_1Dcosp_nu(nuiskey samplekey) {
//********************************************************************

  // Sample overview ---------------------------------------------------
  std::string descrip = "MicroBooNE_CC1MuNp_XSec_1Dcosp_nu sample. \n" \
                        "Target: Ar \n" \
                        "Flux: BNB FHC numu \n" \
                        "Signal: CC 1 mu + Np \n";

  // Setup common settings
  fSettings = LoadSampleSettings(samplekey);
  fSettings.SetDescription(descrip);
  fSettings.SetXTitle("cos#theta_{p}^{reco}");
  fSettings.SetYTitle("d#sigma/dcos#theta_{p}^{reco} (cm^{2}/nucleon)");
  fSettings.SetAllowedTypes("FULL,DIAG/FREE,SHAPE,FIX/SYSTCOV/STATCOV","FIX/FULL");
  fSettings.SetEnuRange(0.0, 6.8);
  fSettings.DefineAllowedTargets("Ar");

  // Plot information
  fSettings.SetTitle("MicroBooNE_CC1MuNp_XSec_1Dcosp_nu");
  fSettings.DefineAllowedSpecies("numu");

  FinaliseSampleSettings();

  // Scaling Setup ---------------------------------------------------
  // ScaleFactor automatically setup for DiffXSec/cm2/Nucleon
  fScaleFactor = GetEventHistogram()->Integral("width") / fNEvents * 1E-38 / TotalIntegratedFlux();

  // Plot Setup -------------------------------------------------------
  std::string inputFile = FitPar::GetDataBase() + "/MicroBooNE/CC1MuNp/CCNp_data_MC_cov_dataRelease.root";
  SetDataFromRootFile(inputFile, "DataXsec_pangle");
  ScaleData(1E-38);
  SetCovarFromRootFile(inputFile, "CovarianceMatrix_pangle");

  // Load smearing matrix ---------------------------------------------
  TFile* inputRootFile = TFile::Open(inputFile.c_str());
  assert(inputRootFile && inputRootFile->IsOpen());
  TH2D* hsmear = (TH2D*) inputRootFile->Get("SmearingMatrix_pangle");
  assert(hsmear);
  fSmearingMatrix = (TH2D*) hsmear->Clone("_ub_ccnp_smearing");
  fSmearingMatrix->SetDirectory(0);
  inputRootFile->Close();
  assert(fSmearingMatrix);

  // Normalize columns
  TH1D* hpx = fSmearingMatrix->ProjectionX("_smearing_py");
  for (int i=1; i<fSmearingMatrix->GetNbinsX()+1; i++) {
    for (int j=1; j<fSmearingMatrix->GetNbinsY()+1; j++) {
      double v = fSmearingMatrix->GetBinContent(i, j) / hpx->GetBinContent(i);
      fSmearingMatrix->SetBinContent(i, j, v);
    }
  }
  delete hpx;

  // Final setup  -----------------------------------------------------
  FinaliseMeasurement();
};


bool MicroBooNE_CC1MuNp_XSec_1Dcosp_nu::isSignal(FitEvent* event) {
  return SignalDef::isCC1MuNp(event, EnuMin, EnuMax);
};


void MicroBooNE_CC1MuNp_XSec_1Dcosp_nu::FillEventVariables(FitEvent* event) {
  if (event->NumFSParticle(2212) == 0) return;
  fXVar = event->GetHMFSParticle(2212)->fP.Vect().CosTheta();
};


void MicroBooNE_CC1MuNp_XSec_1Dcosp_nu::ConvertEventRates() {
  // Do standard conversion
  Measurement1D::ConvertEventRates();

  // Apply MC truth -> reco smearing
  TH1D* truth = (TH1D*) fMCHist->Clone(TString(fMCHist->GetName()) + "_truth");

  for (int ireco=1; ireco<fMCHist->GetNbinsX()+1; ireco++) {
    double total = 0;
    for (int itrue=1; itrue<fMCHist->GetNbinsX()+1; itrue++) {
      total += truth->GetBinContent(itrue) * truth->GetBinWidth(itrue) * fSmearingMatrix->GetBinContent(ireco, itrue);
    }
    fMCHist->SetBinContent(ireco, total / fMCHist->GetBinWidth(ireco));
  }
}

