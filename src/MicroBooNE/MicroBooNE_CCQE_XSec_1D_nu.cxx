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

#include <cmath>
#include "TSpline.h"
#include "MicroBooNE_CCQE_XSec_1D_nu.h"
#include "MicroBooNE_SignalDef.h"

//********************************************************************
MicroBooNE_CCQE_XSec_1D_nu::MicroBooNE_CCQE_XSec_1D_nu(nuiskey samplekey) {
//********************************************************************
  fSettings = LoadSampleSettings(samplekey);
  std::string name = fSettings.GetS("name");
  std::string objSuffix;

  if (name.find("MicroBooNE_CCQE_XSec_1DPmu_nu") != std::string::npos) {
    fDist = kPmu;
    objSuffix = "pmu";
    fSettings.SetXTitle("P_{#mu} (GeV/c)");
    fSettings.SetYTitle("d#sigma/dP_{#mu} (cm^{2}/(GeV/c)/^{40}Ar)");
  }
  else if (name.find("MicroBooNE_CCQE_XSec_1Dcosmu_nu") != std::string::npos) {
    fDist = kCosMu;
    objSuffix = "ctmu";
    fSettings.SetXTitle("cos#theta_{#mu}");
    fSettings.SetYTitle("d#sigma/dcos#theta_{#mu} (cm^{2}/^{40}Ar)");
  }
  else if (name.find("MicroBooNE_CCQE_XSec_1DPp_nu") != std::string::npos) {
    fDist = kPp;
    objSuffix = "pp";
    fSettings.SetXTitle("P_{p} (GeV/c)");
    fSettings.SetYTitle("d#sigma/dP_{p} (cm^{2}/(GeV/c)/^{40}Ar)");
  }
  else if (name.find("MicroBooNE_CCQE_XSec_1Dcosp_nu") != std::string::npos) {
    fDist = kCosP;
    objSuffix = "ctp";
    fSettings.SetXTitle("cos#theta_{p}");
    fSettings.SetYTitle("d#sigma/dcos#theta_{p} (cm^{2}/^{40}Ar)");
  }
  else {
    assert(false);
  }

  // Full or partial phase space
  if (name.find("FullPS") != std::string::npos) {
    fFullPS = true;
  }
  else if (name.find("PartPS") != std::string::npos) {
    fFullPS = false;
  }
  else {
    assert(false);
  }

  // Sample overview ---------------------------------------------------
  std::string descrip = name + " sample, arxiv:2006.00108.\n" \
                        "Target: 40Ar\n" \
                        "Flux: BNB FHC numu\n" \
                        "Signal: CCQE\n";

  fSettings.SetDescription(descrip);
  fSettings.SetTitle(name);
  fSettings.SetAllowedTypes("FULL,DIAG/FREE,SHAPE,FIX/SYSTCOV/STATCOV", "FIX/FULL");
  fSettings.SetEnuRange(0.0, 6.8);
  fSettings.DefineAllowedTargets("Ar");
  fSettings.DefineAllowedSpecies("numu");
  FinaliseSampleSettings();

  // Load data ---------------------------------------------------------
  std::string inputFile = FitPar::GetDataBase() + "/MicroBooNE/CCQE/ccqe_data_release.root";

  std::string dirname = std::string(fFullPS ? "fullps" : "partps");
  SetDataFromRootFile(inputFile, dirname + "/data_" + objSuffix);
  ScaleData(1E-38);

  // ScaleFactor for DiffXSec/cm2/Nucleus (with factor of 40 for Ar40)
  fScaleFactor = GetEventHistogram()->Integral("width") / fNEvents * 1E-38 * 40 / TotalIntegratedFlux();

  SetCovarFromRootFile(inputFile, dirname + "/cov_" + objSuffix);

  // Set up proton range approximation
  size_t npbins = 31;
  float ke_p[npbins] = {
    10, 15, 20, 30, 40, 80, 100, 150, 200, 250, 300, 350, 400,
    450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000,
    1500, 2000, 2500, 3000, 4000, 5000
  };
  float range_p[npbins] = {
    1.887E-1, 3.823E-1, 6.335E-1, 1.296, 2.159, 7.375, 1.092E1,
    2.215E1, 3.627E1, 5.282E1, 7.144E1, 9.184E1, 1.138E2,
    1.370E2, 1.614E2, 1.869E2, 2.132E2, 2.403E2, 2.681E2,
    2.965E2, 3.254E2, 3.548E2, 3.846E2, 4.148E2, 4.454E2,
    7.626E2, 1.090E3, 1.418E3, 1.745E3, 2.391E3, 3.022E3
  };

  TGraph* range_proton = new TGraph(npbins);
  for (size_t i=0; i<npbins; i++) {
      float E = 938 + ke_p[i];
      range_proton->SetPoint(i, std::sqrt(E*E-938*938), range_p[i]/1.396);
  }
  fProtonRange = new TSpline3("range_vs_p_proton", range_proton);

  // Final setup ------------------------------------------------------
  FinaliseMeasurement();
};


bool MicroBooNE_CCQE_XSec_1D_nu::isSignal(FitEvent* event) {
  return SignalDef::MicroBooNE::isCCQE(event, EnuMin, EnuMax, fFullPS, fProtonRange);
};


void MicroBooNE_CCQE_XSec_1D_nu::FillEventVariables(FitEvent* event) {
  if (event->NumFSParticle(13) == 0 || event->NumFSParticle(2212) == 0) return;

  if (fDist == kPmu) {
    fXVar = event->GetHMFSParticle(13)->p() / 1000;
  }
  else if (fDist == kCosMu) {
    fXVar = event->GetHMFSParticle(13)->P3().CosTheta();
  }
  else if (fDist == kPp) {
    fXVar = event->GetHMFSParticle(2212)->p() / 1000;
  }
  else if (fDist == kCosP) {
    fXVar = event->GetHMFSParticle(2212)->P3().CosTheta();
  }
}

