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

#include "TRandom.h"
#include "TSpline.h"
#include "MicroBooNE_SignalDef.h"
#include "FitUtils.h"

namespace SignalDef {
  namespace MicroBooNE {

bool isCC1MuNp(FitEvent* event, double EnuMin, double EnuMax) {
  // Check CC inclusive
  if (!SignalDef::isCCINC(event, 14, EnuMin, EnuMax)) return false;

  // Veto events which don't have exactly 1 FS muon
  if (event->NumFSMuon() != 1) return false;

  // Veto events with FS mesons
  if (event->NumFSPions() != 0) return false;

  // Veto events with FS electrons
  if (event->NumFSElectron() != 0) return false;

  // Veto events with FS photons
  if (event->NumFSPhoton() != 0) return false;

  // Muon momentum above threshold
  if (event->GetHMFSParticle(13)->fP.Vect().Mag() < 100) return false;

  // Leading proton within momentum range
  if (event->NumFSParticle(2212) == 0) return false;
  double plead = event->GetHMFSParticle(2212)->fP.Vect().Mag();
  if (plead > 300 && plead < 1200) return true;

  return false;
}


bool isCCQE(FitEvent* event, double EnuMin, double EnuMax, bool fullPS, TSpline3* prange) {
  // Check CC inclusive
  if (!SignalDef::isCCINC(event, 14, EnuMin, EnuMax)) return false;

  // Veto events which don't have exactly 1 FS muon with p > 100 MeV
  if (event->NumFSMuon() != 1) return false;
  if (event->GetHMFSParticle(13)->fP.Vect().Mag() < 100) return false;

  // Veto events with FS pi0s, electrons, or gammas
  if (event->NumFSPiZero() != 0) return false;
  if (event->NumFSElectron() != 0) return false;
  if (event->NumFSPhoton() != 0) return false;

  // No FS charged pions with p > 70
  if (event->NumFSChargePions() > 0) {
    if (event->GetHMFSParticle(PhysConst::pdg_charged_pions)->p() > 70) return false;
  }

  // One proton over momentum threshold
  size_t np = 0;
  for (size_t i=0; i<event->NParticles(); i++) {
    if (event->GetParticleState(i) == kFinalState &&
        event->GetParticlePDG(i) == 2212 &&
        event->GetParticleP3(i).Mag() >= 300) {
      np++;
    }
  }
  if (np != 1) return false;

  // Muon/proton angles
  FitParticle& muon = *event->GetHMFSMuon();
  FitParticle& proton = *event->GetHMFSProton();

  float dtmp = proton.P3().Angle(muon.P3()) * 180 / TMath::Pi();
  dtmp = dtmp <   0 ? dtmp + 180 : dtmp;
  dtmp = dtmp > 180 ? dtmp - 180 : dtmp;

  float dpmp = (muon.P3().Phi() - proton.P3().Phi()) * 180 / TMath::Pi();
  dpmp = dpmp <   0 ? dpmp + 360 : dpmp;
  dpmp = dpmp > 360 ? dpmp - 360 : dpmp;

  // Phase space cuts
  if (muon.p() < 100 || muon.p() > 1500) return false;
  float ctmax = fullPS ? 0.95 : 0.8;
  if (muon.P3().CosTheta() < -0.65 || muon.P3().CosTheta() > ctmax) return false;

  if (proton.p() < 300 || proton.p() > 1000) return false;
  if (proton.P3().CosTheta() < 0.15 || proton.P3().CosTheta() > 0.95) return false;

  if ((proton.P4() + muon.P4()).Perp() > 350) return false;
  if (std::abs(dtmp - 90) > 55) return false;
  if (std::abs(dpmp - 180) > 35) return false;

  // Fiducial volume for approximate proton containment correction
  TRandom3* rand = new TRandom3(42);
  float rvx = rand->Uniform(   3,  256);
  float rvy = rand->Uniform(-115,  115);
  float rvz = rand->Uniform(   5, 1037);
  TVector3 rv(rvx, rvy, rvz);
  float range = prange->Eval(proton.p());
  TVector3 re = rv + range * proton.P3().Unit();

  if (!(re.X() >=    3 && re.X() <=  256) &&
       (re.Y() >= -115 && re.Y() <=  115) &&
       (re.Z() >=    5 && re.Z() <= 1037)) {
    return false;
  }

  return true;
}

  }  // namespace MicroBooNE
}  // namespace SignalDef

