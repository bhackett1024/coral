/* Copyright 2019 Brian Hackett. Released under the MIT license. */

function assert(v) {
  if (!v) {
    throw new Error("Assertion failed!");
  }
}

// Return the density of the H2O in seawater in kg-H2O/l.
function densityH2O(TC, S) {
  // Subtract the weight of the salts from the weight of a liter of saltwater
  // to get the weight of H2O per liter. Using S = g/kg of sea water.
  const D = density(TC, S); // kg/l
  const S2 = S * D; // g/kg * kg/l => g/l
  return D - S2 / 1000; // kg-H2O/l
}

// Universal gas constant.
const R = 8.3145; // J/(K*mol)

// Faraday, the charge on a mole of electrons.
const F = 96485; // C/mol

// Avogadro's number, the number of things in a mole.
const Avogadro = 6.022e23;

// The number of electrons transferred by a coulomb (1 ampere-second).
const C = Avogadro / F;

// Compute the disassociation constant for H2O <-> H+ + OH-
function waterDisassociation(TC) {
  // FIXME approximating here with the value for 25C.
  return 1e-14;
}

// Compute the disassociation constant for B(OH)3 + H2O <-> H+ + B(OH)4-
function boricAcidDisassociation(TC) {
  // FIXMe approximating here with the value for 25C.
  // S. Zumdahl, S. Zumdahl, D. DeCoste
  // Chemistry, Tenth Edition (Cengage Learning, 2014), pp A24
  return 5.8e-10;
}

// Compute the potential in volts required to generate hydrogen and chlorine
// gas via electrolysis of seawater.
function electrolysisPotential(TC, S, pH) {
  const T = TC + 273.15;

  // Seawater electrolysis is based on the following half-reactions, with
  // standard reduction potentials:
  //
  // Oxidation: 2 Cl- -> Cl2 + 2e-  (E_SRP 1.36V)
  // Reduction: 2 H2O + 2e- -> H2 + 2 OH-  (E_SRP 0.83V)
  // Overall: 2 Cl- + 2H2O -> Cl2 + H2 + 2 OH-
  //
  // The E_SRP of the entire reaction is -1.36V - 0.83V = -2.19V
  //
  // The actual potential of the reaction depends on the concentrations of the
  // reaction components.
  //
  // E = E_SRP - RT/nF ln(Q)
  // units: (J/C) = (J/C) - (J/(K*mol)) * K / (C/mol)

  // F. Millero, Chemical Oceanography (CRC Press, 2013) pp 67
  const Cl_molality = 0.566; // for S=35, mol/kg-H2O
  const Cl = Cl_molality * densityH2O(TC, S); // mol/kg-H2O * kg-H2O/l = mol/l

  // pH is -log([H+]). [OH-] can be calculated from this using the dissociation
  // constant for water, Kw.
  const Kw = waterDisassociation(TC);
  const H = Math.pow(10, -pH);
  const OH = Kw / H;

  // Compute the actual potential using the Nernst equation.
  const Q = Math.pow(OH, 2) / Math.pow(Cl, 2);
  return -2.19 - R * T / (2 * F) * Math.log(Q);
}

// Compute the number of hydroxide ions (mol OH-) needed to raise the pH of a
// liter of water from startPH to endPH.
function hydroxideRequirement(TC, S, DIC, startPH, endPH) {
  // Compute the concentrations (mol/l) of different species at the initial and
  // final pH. For each of these species, determine the number of moles of OH-
  // ions needed to effect the change.
  let newOH = 0;

  // [H+] decreases as pH increases. Each removed H+ ion must have been
  // neutralized by an introduced OH- ion.
  const startH = Math.pow(10, -startPH);
  const endH = Math.pow(10, -endPH);
  assert(endH < startH);
  newOH += startH - endH;

  // [OH-] increases as pH increases. As H+ ions are removed, additional OH-
  // ions are needed so that [H+] and [OH-] are in equilibrium according to the
  // disassociation constant Kw of water.
  const Kw = waterDisassociation(TC);
  const startOH = Kw / startH;
  const endOH = Kw / endH;
  assert(endOH > startOH);
  newOH += endOH - startOH;

  // As pH increases, the concentrations of CO2 and HCO3 will decrease, adding
  // H+ ions to the water which need to be neutralized by OH-. Compute the total
  // number of hydrogen atoms in the carbonate species before and after the pH
  // change.
  const startCarbonate = carbonateConcentrations(TC, S, DIC, startPH);
  const endCarbonate = carbonateConcentrations(TC, S, DIC, endPH);

  // Computed concentrations are in mol/kg-sw, multiply by the density kg-sw/l
  // to get concentrations in mol/l. Ignore CO3, which has no hydrogen.
  const startCO2 = startCarbonate.CO2 * density(TC, S);
  const endCO2 = endCarbonate.CO2 * density(TC, S);
  const startHCO3 = startCarbonate.HCO3 * density(TC, S);
  const endHCO3 = endCarbonate.HCO3 * density(TC, S);

  // CO2 is effectively H2CO3 (carbonic acid), which it must form before it can
  // disassociate to HCO3-. Carbonic acid is diprotic and two H+ ions are
  // generated when it completely disassociates.
  newOH += 2 * (startCO2 - endCO2);

  // HCO3- is the conjugate base of H2CO3 but is an acid itself, releasing an H+
  // ion when it disassociates which must be neutralized.
  newOH += startHCO3 - endHCO3;

  // Boric acid is also present in seawater, in low concentrations compared to
  // carbonate species.
  const Kb = boricAcidDisassociation(TC);

  // The total amount of boric acid species at S=35, in mol/kg-H2O
  // F. Millero, Chemical Oceanography (CRC Press, 2013) pp 67
  const totalBoricKg = 0.0001045 + 0.0003259;
  const totalBoric = totalBoricKg * densityH2O(TC, S);

  // As pH increases, [B(OH)4-] increases. Each such ion needs an OH- ion to
  // form from B(OH)3.
  //
  // Kb = [H+][B(OH)4-]/[B(OH)3]
  // Kb = [H+][B(OH)4-]/([TotalBoric] - [B(OH)4-])
  // Kb[TotalBoric] - Kb[B(OH)4-] = [H+][B(OH)4-]
  // Kb[TotalBoric] = [H+][B(OH)4-] + Kb[B(OH)4-]
  // Kb[TotalBoric] = ([H+] + Kb)[B(OH)4-]
  // [B(OH)4-] = Kb[TotalBoric] / ([H+] + Kb)
  const startBorate = Kb * totalBoric / (startH + Kb);
  const endBorate = Kb * totalBoric / (endH + Kb);
  assert(startBorate < endBorate);
  newOH += endBorate - startBorate;

  return newOH;
}

// Compute the theoretical limit for the amount of seawater that can have its pH
// raised from startPH to endPH using electrolysis, using a power source of
// W watts. Result is in liters per second.
function electrolysisLimit(W, TC, S, DIC, startPH, endPH) {
  // Compute the maximum amperage possible with W watts. In practice the voltage
  // required will be greater than that produced by electrolysisPotential().
  // The actual voltage depends on the environment and will need to be
  // determined experimentally.
  const Amps = W / -electrolysisPotential(TC, S, startPH);

  // Each amp will move one coulomb of electrons per second. Each electron will
  // reduce one H2O molecule to OH-.
  const OH = Amps * C;

  return OH / Avogadro / hydroxideRequirement(TC, S, DIC, startPH, endPH);
}

// Determine the amount of water that can be alkalized by a 100W solar panel in
// conditions expected for the end of the century. Use the average values of
// temperature, salinity, DIC, and pH for the grid coordinate (66,130) closest
// to Opunohu Bay in 2096-2100.
//
// Data source:
//
// Dunne, John; John, Jasmin; Shevliakova, Elena; Stouffer, Ronald; Griffies, Stephen; Malyshev, Sergey; Milly, P.; Sentman, Lori; Adcroft, Alistair; Cooke, William; Dunne, Krista; Hallberg, Robert; Harrison, Matthew; Krasting, John; Levy, Hiram; Phillips, Peter; Samuels, Bonita; Spelman, Michael; Winton, Michael; Wittenberg, Andrew; Zadeh, Niki
// NOAA GFDL GFDL-ESM2M, rcp85 experiment output for CMIP5 AR5, served by ESGF
//
// Tracking IDs:
//
// Temperature: 9cfab26d-0591-4c50-83fc-76f4501d4771
// Salinity: d8ccebdd-c8dd-46ed-a876-382a7ec16c61
// DIC: 2307f360-a117-4a6a-8541-c6099b5f5976
// pH: 9c5ddc62-cdc4-42d8-bfea-5c7800031bc1
//
// Looking at the contributions of the various ions in hydroxideRequirement()
// gives a good sense of the buffering that happens in seawater and makes it
// resistant to changes in pH. Of the hydroxide ions required to change the
// pH to 8.2:
//
// 0.004% (1 in 24107) neutralizes an H+ ion.
// 0.37% (1 in 264) is needed to maintain the equilibrium between [H+] and [OH-].
// 91.3% is needed to neutralize H+ ions added as CO2 and HCO3- disassociate.
// 8.3% is needed to bring boric acid species into equilibrium.
//
// > print(electrolysisLimit(100, 29.5, 35.4, 0.002229376 / density(29.5, 35.4),
//                           7.76, 8.2));
// 2.1220809101172438

// Calculate the theoretical potential needed for seawater electrolysis.
// > print(electrolysisPotential(29.5, 35.4, 7.76));
// -1.8304110483470388

// Calculate the number of moles of hydroxide ions needed to alkalize one liter.
// > print(hydroxideRequirement(29.5, 35.4, 0.002229376 / density(29.5, 35.4),
//                              7.76, 8.2));
// 0.00026682691299919205
