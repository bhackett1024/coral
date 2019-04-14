/* Copyright 2019 Brian Hackett. Released under the MIT license. */

"use strict";

const { Terms } = require("./units");

const UniversalGasConstant = Terms.JoulesPerKelvinMole(8.3145);

// The charge on a mole of electrons.
const Faraday = Terms.CoulombsPerMole(96485);

// Compute the disassociation constant for H2O <-> H+ + OH-
function waterDisassociation() {
  // FIXME approximating here with the value for 25C. This also has a different
  // value in seawater due to different ionic activities.
  return Terms.MolesSquaredPerLiterSquared(1e-14);
}

// For species of interest, track the following information:
//
// conductivity: limiting equivalent ionic conductivities (S cm^2/mol)
// diffusion: diffusion coefficient (m^2/s)
// charge: charge on the species
//
// Except where noted below, conductivity and diffusion coefficients are from:
// http://www.aqion.de/site/194
const Species = {
  "Na+": { conductivity: 50, diffusion: 1.33 * 1e-9, charge: 1 },
  "Cl-": { conductivity: 76.2, diffusion: 2.03 * 1e-9, charge: -1 },
  "Mg^{2+}": { conductivity: 53, diffusion: 0.705 * 1e-9, charge: 2 },
  "Ca^{2+}": { conductivity: 59.6, diffusion: 0.793 * 1e-9, charge: 2 },
  "K+": { conductivity: 73.6, diffusion: 1.96 * 1e-9, charge: 1 },
  "SO4^{2-}": { conductivity: 80.4, diffusion: 1.07 * 1e-9, charge: -2 },
  "OH-": { conductivity: 197.9, diffusion: 5.27 * 1e-9, charge: -1 },
  "H+": { conductivity: 349.6, diffusion: 9.31 * 1e-9, charge: 1 },
  "HCO3-": { conductivity: 44.3, diffusion: 1.18 * 1e-9, charge: -1 },
  "CO3^{2-}": { conductivity: 71.7, diffusion: 0.955 * 1e-9, charge: -2 },

  // The CO2 diffusion coefficient is from:
  //
  // Shane P. Cadogan, Geoffrey C. Maitland, J. P. Martin Trusler
  // Diffusion Coefficients of CO2 and N2 in Water at Temperatures between
  //   298.15 K and 423.15 K at Pressures up to 45 MPa
  // J. Chem. Eng. Data 2014, 59, 519−525
  //
  // This is at a high pressure (14 MPa = 138 atm), which will affect things. Oh well.
  "CO2": { diffusion: 2.233 * 1e-9 },

  // Estimate the diffusion coefficient for OCl- using that of CL-.
  "OCl-": { diffusion: 2.03 * 1e-9 }
};
for (const entry of Object.values(Species)) {
  if (entry.conductivity) {
    entry.conductivity = Terms.SiemensSquareCentimersPerMole(entry.conductivity);
  }
  if (entry.diffusion) {
    entry.diffusion = Terms.SquareMetersPerSecond(entry.diffusion);
  }
  if (entry.charge) {
    entry.charge = Terms.Number(entry.charge);
  }
}

// Concentrations of the major ionic species in seawater (99.9% of ions,
// S=35, T=25C), in mol/kg-H2O, are from:
//
// Frank J. Millero
// Chemical Oceanography (CRC Press, 2013) page 67
const SeawaterIons = {
  "Na+": Terms.MolesPerSeawaterKg(0.486),
  "Cl-": Terms.MolesPerSeawaterKg(0.567),
  "Mg^{2+}": Terms.MolesPerSeawaterKg(0.055),
  "Ca^{2+}": Terms.MolesPerSeawaterKg(0.011),
  "K+": Terms.MolesPerSeawaterKg(0.011),
  "SO4^{2-}": Terms.MolesPerSeawaterKg(0.029),
  "B(OH)3": Terms.MolesPerSeawaterKg(3.3e-4),
  "B(OH)4-": Terms.MolesPerSeawaterKg(1.05e-4)
};

function seawaterConcentration(name, S) {
  return SeawaterIons[name].mul(S).div(Terms.Salinity(35));
}

module.exports = {
  UniversalGasConstant, Faraday, waterDisassociation, seawaterConcentration, SeawaterIons, Species
}
