/* Copyright 2019 Brian Hackett. Released under the MIT license. */

"use strict";

const { carbonateConcentrations } = require("./carbonate");
const { Units, Terms } = require("./units");

function aragoniteSaturation(T, S, DIC, pH) {
  // Compute CO3 concentration based on the environment.
  const CO3 = carbonateConcentrations(T, S, DIC, pH).CO3;

  // Ca concentration varies little relative to other salts in seawater, and can
  // be estimated using the water's salinity.
  //
  // Frank J. Millero
  // Chemical Oceanography (CRC Press, 2013), pp 296
  const Ca = Terms.MolesPerSeawaterKg(S.normalize(Units.Salinity) * 0.0002934);

  // Compute the solubility product of aragonite based on the environment.
  const Ksp = associationKspAragonite(T, S);

  // In solution, CO3 and Ca ions will dissolve or precipitate in order to
  // reach the following equilibrium:
  //
  // Ksp = [CO3][Ca]
  //
  // Concentrations with this Ksp are based on mol/kg-sw instead of mol/l.
  // When this equality holds, the solution is saturated.
  //
  // If the concentration of CO3 or Ca were to decrease then the solution is
  // undersaturated and any available CaCO3 will eventually dissolve to restore
  // the equilibrium.
  //
  // If the concentration of CO3 or Ca were to increase then the solution is
  // supersaturated and CaCO3 will eventually precipitate to restore the
  // equilibrium.
  //
  // The saturation state is the ratio [CO3][Ca]/Ksp. If this is 1 then the
  // solution is saturated, if below 1 the solution is undersaturated, and
  // if above 1 the solution is supersaturated.
  return CO3.mul(Ca).div(Ksp).number();
}

// The method below computes the solubility product of Aragonite for a given
// temperature and salinity. These calculations are from:
//
// Frank J. Millero
// Chemical Oceanography (CRC Press, 2013) pp 272
//
// This uses ln instead of log in some places, which doesn't produce correct
// results. The original paper (my copy is illegible in places, alas) cited
// below, uses log(), so that's what is used here.
//
// Alfonso Mucci
// The solubility of calcite and aragonite in seawater at various salinities, temperatures, and one atmosphere total pressure
// American Journal of Science (September 1983)
function associationKspAragonite(T, S) {
  T = T.normalize(Units.Kelvin);
  S = S.normalize(Units.Salinity);
  const logKsp_thermodynamic = -171.945 - 0.077993 * T + 2903.293 / T + 71.595 * Math.log10(T);
  const A = -0.068393 + 0.0017276 * T + 88.135 / T;
  const B = -0.10018;
  const C = 0.0059415;
  const logKsp = logKsp_thermodynamic + A * Math.sqrt(S) + B * S + C * Math.pow(S, 1.5);
  return Terms.MolesSquaredPerSeawaterKgSquared(Math.pow(10, logKsp));
}

module.exports = { aragoniteSaturation };
