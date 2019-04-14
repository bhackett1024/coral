/* Copyright 2019 Brian Hackett. Released under the MIT license. */

const { assert } = require("./utils");
const { Units, Terms } = require("./units");

// Compute the concentration of carbonate species.
function carbonateConcentrations(T, S, DIC, pH) {
  // FIXME per the paper cited below defining K1/K2, instead of [H+] this should be using
  // [H+] + [HSO4-] + [HF], but trying to make this correction does not seem to
  // give correct results.
  const H = pH.concentrationH();

  // At equilibrium, the concentrations of the compounds in each reaction are
  // related according to the following equations:
  //
  // K1 = [H+][HCO3-]/[CO2]
  // K2 = [H+][CO3^{2-}]/[HCO3-]
  //
  // K1 and K2 are constants determined by the seawater's temperature and
  // salinity. [H2O] is omitted from the first equation as it is essentially
  // constant and will be folded into the computation of K1. Given K1, K2, and
  // the total DIC, the concentrations [CO2], [HCO3^-], and [CO3^{2-}] can be
  // determined from [H+] according to the derivations below.
  const K1 = associationK1(T, S);
  const K2 = associationK2(T, S);

  // DIC = [CO2] + [HCO3-] + [CO3^{2-}]
  // DIC = [CO2] + K1*[CO2]/[H+] + K2*K1*[CO2]/[H+]^2 (substituting equations below)
  // DIC = [CO2](1 + K1/[H+] + K2*K1/[H+]^2)
  // [CO2] = DIC / (1 + K1/[H+] + K2*K1/[H+]^2)
  const CO2 = DIC.div(Terms.Number(1).add(K1.div(H)).add(K1.mul(K2).div(H.mul(H))));

  // [H+][HCO3-]/[CO2] = K1
  // [HCO3-] = K1*[CO2]/[H+]
  const HCO3 = K1.mul(CO2).div(H);

  // [H+][CO3^{2-}]/[HCO3-] = K2
  // [CO3^{2-}] = K2*[HCO3-]/[H+]
  // [CO3^{2-}] = K2*(K1*[CO2]/[H+])/[H+]
  // [CO3^{2-}] = K2*K1*[CO2]/[H+]^2
  const CO3 = K2.mul(K1).mul(CO2).div(H.mul(H));

  return { CO2, HCO3, CO3 };
}

// Generate data for a Bjerrum plot of carbonate concentrations against pH.
function generateBjerrumData(T, S, startPH, endPH, numPoints) {
  startPH = startPH.normalize(Units.pH);
  endPH = endPH.normalize(Units.pH);
  const pHs = [];
  const concentrations = [];
  for (let i = 0; i <= numPoints; i++) {
    const pH = Terms.pH(startPH + (i / numPoints) * (endPH - startPH));
    pHs.push(pH);
    concentrations.push(carbonateConcentrations(T, S, Terms.Molarity(1), pH));
  }
  return {
    pH: pHs,
    CO2: concentrations.map(v => v.CO2),
    HCO3: concentrations.map(v => v.HCO3),
    CO3: concentrations.map(v => v.CO3)
  };
}

// The methods below compute values of K1 and K2 for a given temperature and
// salinity. These calculations are from:
//
// Frank J. Millero, Taylor B. Graham, Fen Huang, Hector Bustos-Serrano, Denis Pierrot
// Dissociation constants of carbonic acid in seawater as a function of salinity and temperature
// Marine Chemistry 100 (2006) 80–94

// Compute the stoichiometric association constant for the first ionization of
// carbonic acid in seawater, as a function of temperature and salinity.
function associationK1(T, S) {
  T = T.normalize(Units.Kelvin);
  S = S.normalize(Units.Salinity);
  const pK1_0 = -126.34048 + 6320.813 / T + 19.568224 * Math.log(T);
  const A = 13.4191 * Math.sqrt(S) + 0.0331 * S - 0.0000533 * Math.pow(S, 2);
  const B = -530.123 * Math.sqrt(S) - 6.103 * S;
  const C = -2.06950 * Math.sqrt(S);
  const pK1 = pK1_0 + A + B / T + C * Math.log(T);
  return Terms.MolesPerSeawaterKg(Math.pow(10, -pK1));
}

// Compute the stoichiometric association constant for the second ionization of
// carbonic acid in seawater, as a function of temperature and salinity.
function associationK2(T, S) {
  T = T.normalize(Units.Kelvin);
  S = S.normalize(Units.Salinity);
  const pK2_0 = -90.18333 + 5143.692 / T + 14.613358 * Math.log(T);
  const A = 21.0894 * Math.sqrt(S) + 0.1248 * S - 0.0003687 * Math.pow(S, 2);
  const B = -772.483 * Math.sqrt(S) - 20.051 * S;
  const C = -3.3336 * Math.sqrt(S);
  const pK2 = pK2_0 + A + B / T + C * Math.log(T);
  return Terms.MolesPerSeawaterKg(Math.pow(10, -pK2));
}

module.exports = { generateBjerrumData, carbonateConcentrations };
