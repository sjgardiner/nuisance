// Copyright 2016-2021 L. Pickering, P Stowell, R. Terri, C. Wilkinson, C. Wret

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

#include "ElectronScattering_e4nu.h"
#include "E4nu_SignalDef.h"

ElectronScattering_e4nu::ElectronScattering_e4nu( nuiskey samplekey ) {

  // Sample overview
  std::string descrip =
      "Electron Scattering e4nu sample. \n"
      "Target: Multiple \n"
      "Flux: Energy should match data being handled \n"
      "Signal: 1-proton 0-pion electron scattering events"
      " (with phase space cuts) \n";
  fSettings = LoadSampleSettings( samplekey );
  fSettings.SetDescription( descrip );
  fSettings.DefineAllowedSpecies( "electron" );
  fSettings.SetTitle( "Electron" );
  fSettings.SetAllowedTypes( "FIX/DIAG", "FIX,FREE,SHAPE/DIAG/NORM/MASK" );
  fSettings.SetXTitle( "E_{QE}" );
  fSettings.SetYTitle( "d#sigma/dE_{QE} (#mub/GeV)" );
  //fIsNoWidth = true;

  FinaliseSampleSettings();

  // Plot Setup -------------------------------------------------------
  SetDataFromName( fSettings.GetS("originalname") );
  SetCovarFromDiagonal();

  // Scaling Setup ---------------------------------------------------
  // ScaleFactor automatically setup for DiffXSec/cm2/Nucleon
  fScaleFactor = GetEventHistogram()->Integral() * 1E-38 / double(fNEvents)
    / GetFluxHistogram()->Integral();
    // / TotalIntegratedFlux();

  // Convert from cm^2 / nucleon to microbarn / nucleus
  // TODO: revisit the handling of the mass number here
  fScaleFactor *= 2. / 1e-30;

  std::cout << "Event Integral = " << GetEventHistogram()->Integral() << '\n';
  std::cout << "Flux Integral = " << GetFluxHistogram()->Integral() << '\n';
  std::cout << "FNEvents = " << fNEvents << '\n';
  std::cout << "ScaleFactor = " << fScaleFactor << std::endl;

  // Finish up
  FinaliseMeasurement();
};

void ElectronScattering_e4nu::SetDataFromName( const std::string& name ) {

  if ( name != "e4nu_EQE" ) NUIS_ABORT( "Unrecognized e4nu data set" );

  std::string in_data_file_name = FitPar::GetDataBase()
    + "/Electron/41586_2021_4046_MOESM2_ESM.csv";

  std::ifstream in_data_file( in_data_file_name );

  if ( !in_data_file.good() ) {
    NUIS_ABORT("Failed to open e-scattering database file: "
      + in_data_file_name );
  }

  // Pull from file and define vectors
  std::vector<double> low_bin_edge;
  std::vector<double> high_bin_edge;
  std::vector<double> cross_section;
  std::vector<double> error;
  std::string line;

  //Accounts for the lines prior to actually starting the dataset
  std::getline( in_data_file, line );
  std::getline( in_data_file, line );
  std::getline( in_data_file, line );

  // Goes through dataset and takes values
  while ( std::getline(in_data_file, line) ) {

    // Defines variables for loop
    std::stringstream line_ss( line );
    char dummy;
    int bin_index;
    double value_1, value_2, value_3, value_4;

    // Goes through file and assigns values to variables
    line_ss >> bin_index >> dummy >> value_1 >> dummy >> value_2
      >> dummy >> value_3 >> dummy >> value_4;

    //Puts variables in their respective vectors
    low_bin_edge.push_back( value_1 );
    high_bin_edge.push_back( value_2 );
    cross_section.push_back( value_3 );
    error.push_back( value_4 );

  }

  std::vector<double> bin_edges = low_bin_edge;

  int n_bins = low_bin_edge.size();

  bin_edges.push_back( high_bin_edge[n_bins - 1] );

  double* bin_edges_array = &bin_edges[0];

  // Form the histograms for data and MC
  fDataHist = new TH1D( (fName + "_data").c_str(), (fName + "_data").c_str(),
    n_bins, bin_edges_array );

  fMCHist = dynamic_cast< TH1D* >( fDataHist->Clone("MC") );

  // Set the bin contents for the data histogram
  for( int i = 0; i < n_bins; i++ ) {
    fDataHist->SetBinContent( i + 1, cross_section[i] );
    fDataHist->SetBinError( i + 1, error[i] );
  }

}

void ElectronScattering_e4nu::FillEventVariables(FitEvent* event) {

  const int ELECTRON = 11;
  const double PROTON_MASS = 0.938; // GeV / c^2
  const double ELECTRON_MASS = 0.000511; // GeV / c^2
  // Note that there is a misprint in the paper. Value confirmed with Afro
  // Papadopoulou
  const double AVG_REMOVAL_ENERGY_12C = 0.021; // GeV

  if ( event->NumFSParticle(ELECTRON) == 0 ) return;

  FitParticle* eout = event->GetHMFSParticle( ELECTRON );

  const TLorentzVector& p4eout = eout->fP;
  double El = p4eout.E() / 1e3; // GeV
  double kl = p4eout.P() / 1e3; // GeV
  double cosl = p4eout.Vect().CosTheta();
  double ml = ELECTRON_MASS;

  double E_qe = ( 2.*PROTON_MASS*AVG_REMOVAL_ENERGY_12C
    + 2.*PROTON_MASS*El - ml*ml ) / 2. / ( PROTON_MASS - El + kl*cosl );

  fXVar = E_qe;
}

bool ElectronScattering_e4nu::isSignal( FitEvent* event ) {
  double dummy = 0.;
  return SignalDef::e4nu::isEM1p0pi( event, dummy, dummy );
};

int ElectronScattering_e4nu::GetNDOF() {
  return fDataHist->GetNbinsX();
}

double ElectronScattering_e4nu::GetLikelihood() { return 0.0; }

void ElectronScattering_e4nu::SetFitOptions(std::string opt) { return; }

TH1D* ElectronScattering_e4nu::GetMCHistogram(void) { return fMCHist; }

TH1D* ElectronScattering_e4nu::GetDataHistogram(void) {
  return fDataHist;
}
