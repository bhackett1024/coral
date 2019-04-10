/* Copyright 2019 Brian Hackett. Released under the MIT license. */

// Return the density of seawater in kg/l.
function density(TC, S) {
  // FIXME approximating here, but this could be calculated more precisely.
  return 1.025;
}

// Molarity is mol/l, while molality is mol/kg-H2O
// These have similar values in solutions but different sources might use one or
// the other so we have to convert between them.
function molarityToMolality(TC, S, M) {
  // Subtract the weight of the salts from the weight of a liter of saltwater
  // to get the weight of H2O per liter. Using S = g/kg of sea water.
  const D = density(TC, S); // kg/l
  const S2 = S * D; // g/kg * kg/l => g/l
  const D_H2O = D - S2 / 1000; // kg-H2O/l

  // mol/l / (kg-H2O/l) => mol/kg-H2O
  return M / D_H2O;
}

// Compute the concentration of carbonate species in mol/kg.
// TC = temperature in degrees C
// S = salinity (g/kg)
// DIC = dissolved inorganic carbon, [CO2] + [HCO3-] + [CO3^{2-}], mol/kg
function carbonateConcentrations(TC, S, DIC, pH) {
  const Hmolarity = Math.pow(10, -pH);
  const H = Hmolarity / density(TC, S);

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
  const K1 = associationK1(TC, S);
  const K2 = associationK2(TC, S);

  // DIC = [CO2] + [HCO3-] + [CO3^{2-}]
  // DIC = [CO2] + K1*[CO2]/[H+] + K2*K1*[CO2]/[H+]^2 (substituting equations below)
  // DIC = [CO2](1 + K1/[H+] + K2*K1/[H+]^2)
  // [CO2] = DIC / (1 + K1/[H+] + K2*K1/[H+]^2)
  const CO2 = DIC / (1 + K1 / H + K1 * K2 / (H * H));

  // [H+][HCO3-]/[CO2] = K1
  // [HCO3-] = K1*[CO2]/[H+]
  const HCO3 = K1 * CO2 / H;

  // [H+][CO3^{2-}]/[HCO3-] = K2
  // [CO3^{2-}] = K2*[HCO3-]/[H+]
  // [CO3^{2-}] = K2*(K1*[CO2]/[H+])/[H+]
  // [CO3^{2-}] = K2*K1*[CO2]/[H+]^2
  const CO3 = K2 * K1 * CO2 / (H * H);

  return { CO2, HCO3, CO3 };
}

// Generate the data used for comparing CO3^{2-} concentrations according to pH.
// print(carbonateConcentrations(28, 35, 1, 7).CO3 /
//       carbonateConcentrations(28, 35, 1, 8.1).CO3);
// 0.08573205565286607
// print(carbonateConcentrations(28, 35, 1, 9).CO3 /
//       carbonateConcentrations(28, 35, 1, 8.1).CO3);
// 4.114304857123481

//print(carbonateConcentrations(28, 36, 2.252298, 7.76).CO3);

function densityH2O(TC, S) {
  // Subtract the weight of the salts from the weight of a liter of saltwater
  // to get the weight of H2O per liter. Using S = g/kg of sea water.
  const D = density(TC, S); // kg/l
  const S2 = S * D; // g/kg * kg/l => g/l
  return D - S2 / 1000; // kg-H2O/l
}

// Generate data for a Bjerrum plot of carbonate concentrations against pH.
function generateBjerrumData(TC, S, start_pH, end_pH, numPoints) {
  const pHs = [];
  const concentrations = [];
  for (let i = 0; i <= numPoints; i++) {
    const pH = start_pH + (i / numPoints) * (end_pH - start_pH);
    pHs.push(pH);
    concentrations.push(carbonateConcentrations(TC, S, 1, pH));
  }
  return {
    pH: pHs,
    CO2: concentrations.map(v => v.CO2),
    HCO3: concentrations.map(v => v.HCO3),
    CO3: concentrations.map(v => v.CO3)
  };
}

// Generate the data used in the carbonate species chart.
// generateBjerrumData(28, 35, 4, 10, 20);

// The methods below compute values of K1 and K2 for a given temperature and
// salinity. These calculations are from:
//
// Frank J. Millero, Taylor B. Graham, Fen Huang, Hector Bustos-Serrano, Denis Pierrot
// Dissociation constants of carbonic acid in seawater as a function of salinity and temperature
// Marine Chemistry 100 (2006) 80–94

// Compute the stoichiometric association constant for the first ionization of
// carbonic acid in seawater, as a function of temperature and salinity.
function associationK1(TC, S) {
  const T = TC + 273.15;
  const pK1_0 = -126.34048 + 6320.813 / T + 19.568224 * Math.log(T);
  const A = 13.4191 * Math.sqrt(S) + 0.0331 * S - 0.0000533 * Math.pow(S, 2);
  const B = -530.123 * Math.sqrt(S) - 6.103 * S;
  const C = -2.06950 * Math.sqrt(S);
  const pK1 = pK1_0 + A + B / T + C * Math.log(T);
  return Math.pow(10, -pK1);
}

// Compute the stoichiometric association constant for the second ionization of
// carbonic acid in seawater, as a function of temperature and salinity.
function associationK2(TC, S) {
  const T = TC + 273.15;
  const pK2_0 = -90.18333 + 5143.692 / T + 14.613358 * Math.log(T);
  const A = 21.0894 * Math.sqrt(S) + 0.1248 * S - 0.0003687 * Math.pow(S, 2);
  const B = -772.483 * Math.sqrt(S) - 20.051 * S;
  const C = -3.3336 * Math.sqrt(S);
  const pK2 = pK2_0 + A + B / T + C * Math.log(T);
  return Math.pow(10, -pK2);
}

module.exports = { generateBjerrumData, carbonateConcentrations, density };
