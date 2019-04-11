/* Copyright 2019 Brian Hackett. Released under the MIT license. */

const { assert } = require("./utils");
const { Units, Terms } = require("./units");
const { carbonateConcentrations, density } = require("./carbonate");

// Universal gas constant.
const R = 8.3145; // J/(K*mol)

// Faraday, the charge on a mole of electrons.
const F = 96485; // C/mol

// Avogadro's number, the number of things in a mole.
const Avogadro = 6.022e23;

// The number of electrons transferred by a coulomb (1 ampere-second).
const C = Avogadro / F;

// Compute the disassociation constant for H2O <-> H+ + OH-
function waterDisassociation() {
  // FIXME approximating here with the value for 25C.
  return 1e-14;
}

// Compute the disassociation constant for B(OH)3 + H2O <-> H+ + B(OH)4-
function boricAcidDisassociation() {
  // FIXMe approximating here with the value for 25C.
  // S. Zumdahl, S. Zumdahl, D. DeCoste
  // Chemistry, Tenth Edition (Cengage Learning, 2014), pp A24
  return 5.8e-10;
}

// Compute the potential in volts required to generate hydrogen and chlorine
// gas via electrolysis of seawater.
function electrolysisPotential(T, S, pH) {
  T = T.normalize(Units.Kelvin);
  S = S.normalize(Units.Salinity);

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
  const ClStandard = Terms.MolesPerSeawaterKg(0.566);

  const Cl = ClStandard.normalize(Units.Molarity) * S / 35;

  // pH is -log([H+]). [OH-] can be calculated from this using the dissociation
  // constant for water, Kw.
  const Kw = waterDisassociation();
  const H = pH.concentrationH().normalize(Units.Molarity);
  const OH = Kw / H;

  // Compute the actual potential using the Nernst equation.
  const Q = Math.pow(OH, 2) / Math.pow(Cl, 2);
  return Terms.Volts(-2.19 - R * T / (2 * F) * Math.log(Q));
}

// Compute the number of hydroxide ions (mol OH-) needed to raise the pH of a
// liter of water from startPH to endPH.
function hydroxideRequirement(T, S, DIC, startPH, endPH) {
  // Compute the concentrations (mol/l) of different species at the initial and
  // final pH. For each of these species, determine the number of moles of OH-
  // ions needed to effect the change.
  let newOH = 0;

  // [H+] decreases as pH increases. Each removed H+ ion must have been
  // neutralized by an introduced OH- ion.
  const startH = startPH.concentrationH().normalize(Units.Molarity);
  const endH = endPH.concentrationH().normalize(Units.Molarity);
  assert(endH < startH);
  newOH += startH - endH;

  // [OH-] increases as pH increases. As H+ ions are removed, additional OH-
  // ions are needed so that [H+] and [OH-] are in equilibrium according to the
  // disassociation constant Kw of water.
  const Kw = waterDisassociation();
  const startOH = Kw / startH;
  const endOH = Kw / endH;
  assert(endOH > startOH);
  newOH += endOH - startOH;

  // As pH increases, the concentrations of CO2 and HCO3 will decrease, adding
  // H+ ions to the water which need to be neutralized by OH-. Compute the total
  // number of hydrogen atoms in the carbonate species before and after the pH
  // change.
  const startCarbonate = carbonateConcentrations(T, S, DIC, startPH);
  const endCarbonate = carbonateConcentrations(T, S, DIC, endPH);

  const startCO2 = startCarbonate.CO2.normalize(Units.Molarity);
  const endCO2 = endCarbonate.CO2.normalize(Units.Molarity);
  const startHCO3 = startCarbonate.HCO3.normalize(Units.Molarity);
  const endHCO3 = endCarbonate.HCO3.normalize(Units.Molarity);

  // CO2 is effectively H2CO3 (carbonic acid), which it must form before it can
  // disassociate to HCO3-. Carbonic acid is diprotic and two H+ ions are
  // generated when it completely disassociates.
  newOH += 2 * (startCO2 - endCO2);

  // HCO3- is the conjugate base of H2CO3 but is an acid itself, releasing an H+
  // ion when it disassociates which must be neutralized.
  newOH += startHCO3 - endHCO3;

  // Boric acid is also present in seawater, in low concentrations compared to
  // carbonate species.
  const Kb = boricAcidDisassociation();

  // The total amount of boric acid species at S=35, in mol/kg-H2O
  // F. Millero, Chemical Oceanography (CRC Press, 2013) pp 67
  const totalBoricStandard = Terms.MolesPerSeawaterKg(0.0001045 + 0.0003259);

  const totalBoric =
    totalBoricStandard.normalize(Units.Molarity) * S.normalize(Units.Salinity) / 35;

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

  return Terms.Molarity(newOH);
}

// Compute the theoretical limit for the amount of seawater that can have its pH
// raised from startPH to endPH using electrolysis, using a power source W.
function electrolysisLimit(W, T, S, DIC, startPH, endPH) {
  // Compute the maximum amperage possible with W watts. In practice the voltage
  // required will be greater than that produced by electrolysisPotential().
  // The actual voltage depends on the environment and will need to be
  // determined experimentally.
  const potential = -electrolysisPotential(T, S, startPH).normalize(Units.Volts);
  const Amps = W.normalize(Units.Watts) / potential;

  // Each amp will move one coulomb of electrons per second. Each electron will
  // reduce one H2O molecule to OH-.
  const OH = Amps * C;

  const requirement = hydroxideRequirement(T, S, DIC, startPH, endPH).normalize(Units.Molarity);

  return Terms.LitersPerSecond(OH / Avogadro / requirement);
}

module.exports = { electrolysisLimit, electrolysisPotential, hydroxideRequirement, C, Avogadro };
