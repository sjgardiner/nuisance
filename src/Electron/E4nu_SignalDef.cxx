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

#include "E4nu_SignalDef.h"
#include "FitUtils.h"

bool SignalDef::e4nu::isEM1p0pi( FitEvent* event, double EnuMin,
  double EnuMax )
{
  const double PROTON_MASS = 0.938; // GeV / c^2
  const int ELECTRON = 11;
  const int PHOTON = 22;
  const int PI_PLUS = 211;
  const int PI_MINUS = -211;
  const int PROTON = 2212;

  // First check that this is an electron scattering event
  if ( !event->HasISParticle(ELECTRON) ) return false;
  if ( !event->HasFSParticle(ELECTRON) ) return false;

  // Get the initial and final electron 4-momenta
  const TLorentzVector& initial_electron
    = event->GetHMISParticle( ELECTRON )->fP;

  const TLorentzVector& final_electron
    = event->GetHMFSParticle( ELECTRON )->fP;

  // Store some components of these 4-momenta for easy retrieval
  double my_Ev = initial_electron.E() / 1e3; // GeV

  double my_El = final_electron.E() / 1e3; // GeV
  double my_pxl = final_electron.X() / 1e3; // GeV
  double my_pyl = final_electron.Y() / 1e3; // GeV
  double my_pzl = final_electron.Z() / 1e3; // GeV

  // Calculates the value of Q^2 for the given event
  TLorentzVector q = initial_electron - final_electron;
  double Q2 = -1. * ( q[3]*q[3] - q[2]*q[2]
    - q[1]*q[1] - q[0]*q[0] ) / 1e6; // GeV^2

  // Calculates W for the event
  double W = std::sqrt( -Q2 + PROTON_MASS*PROTON_MASS
    + 2.*PROTON_MASS*(my_Ev - my_El) ); // GeV

  // Calculates total momentum (GeV) and angle (degrees) WRT z-axis
  double p_tot = std::sqrt( my_pxl*my_pxl + my_pyl*my_pyl + my_pzl*my_pzl );
  double theta_z = std::acos( my_pzl / p_tot ) * 180. / M_PI;

  // Sets counters for each type of hadron
  int proton_counter = 0;

  int pion_counter = 0;

  int photon_counter = 0;

  // Get the total number of particles in the event
  int num_particles = event->NParticles();

  // Goes through events with 1.161 GeV to check electron requirements
  if ( my_Ev == 1.161 && Q2 > 0.1 && p_tot > 0.40 && theta_z > (17. + 7./p_tot)
    && W < 2. )
  {
    // Checks individual requirements for final-state hadrons
    for ( int n = 0; n < num_particles; ++n ) {

      // Skip particles that aren't part of the final state
      int state = event->GetParticleState( n );
      if ( state != kFinalState ) continue;

      // Retrieve the particle 4-momentum and PDG code
      const TLorentzVector& h_v = event->GetParticleP4( n );

      int my_pdgf = event->GetParticlePDG( n );

      //Determines type of hadron and checks related requirements
      //Proton requirements
      if ( my_pdgf == PROTON ) {
        double proton_momentum = std::sqrt( h_v[0]*h_v[0] + h_v[1]*h_v[1]
          +h_v[2]*h_v[2] ) / 1e3; // GeV
        double theta_proton = std::acos( h_v[2] / 1e3
          / proton_momentum ) * 180. / M_PI;
        if ( proton_momentum > 0.3 && theta_proton > 12. ) {
          proton_counter += 1;
        }
      }

      //Photon requirements
      else if ( my_pdgf == PHOTON ) {
        double photon_momentum = std::sqrt( h_v[0]*h_v[0] + h_v[1]*h_v[1]
          + h_v[2]*h_v[2] ) / 1e3; // GeV

        if ( photon_momentum > 0.3 ){
          photon_counter += 1;
        }
      }

      //Positive pion requirements
      else if ( my_pdgf == PI_PLUS ) {
        double pion_momentum = std::sqrt( h_v[0]*h_v[0] + h_v[1]*h_v[1]
          + h_v[2]*h_v[2] ) / 1e3; // GeV

        double theta_pion = std::acos( h_v[2] / 1e3 / pion_momentum ) * 180. / M_PI;
        if ( pion_momentum > 0.15 && theta_pion > 12. ) {
          pion_counter += 1;
        }
      }

      //Negative pion requirements
      else if ( my_pdgf == PI_MINUS ) {
        double pion_momentum = std::sqrt( h_v[0]*h_v[0] + h_v[1]*h_v[1]
          + h_v[2]*h_v[2] ) / 1e3;
        double theta_pion = std::acos( h_v[2] / 1e3 / pion_momentum ) * 180. / M_PI;
        if ( pion_momentum > 0.15 && theta_pion > (17. + 4./pion_momentum) ) {
          pion_counter += 1;
        }
      }
    } // Loop over particles
  } // Electron cuts

  // Counts the event if it meets the requirements
  if ( proton_counter == 1 && pion_counter == 0 && photon_counter == 0 ) {
    return true;
  }

  return false;
}
