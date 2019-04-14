/* Copyright 2019 Brian Hackett. Released under the MIT license. */

"use strict";

const { assert } = require("./utils");
const { Units, Terms } = require("./units");
const { carbonateConcentrations, density } = require("./carbonate");
const { waterDisassociation, Faraday, UniversalGasConstant, seawaterConcentration } = require("./constants");

// Compute the disassociation constant for B(OH)3 + H2O <-> H+ + B(OH)4-
function boricAcidDisassociation() {
  // FIXME approximating here with the value for 25C.
  // S. Zumdahl, S. Zumdahl, D. DeCoste
  // Chemistry, Tenth Edition (Cengage Learning, 2014), pp A24
  return Terms.Molarity(5.8e-10);
}

// Compute the disassociation constant for HSO4- <-> H+ + SO4^{2-}
function bisulfateDisassociation() {
  // FIXME approximating here with the value for 25C.
  // S. Zumdahl, S. Zumdahl, D. DeCoste
  // Chemistry, Tenth Edition (Cengage Learning, 2014), pp A24
  return Terms.Molarity(1.2e-2);
}

// Compute the potential in volts required to generate hydrogen and chlorine
// gas via electrolysis of seawater.
function electrolysisPotential(T, S, pH) {
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

  const Cl = seawaterConcentration("Cl-", S);

  // pH is -log([H+]). [OH-] can be calculated from this using the dissociation
  // constant for water, Kw.
  const Kw = waterDisassociation();
  const OH = Kw.div(pH.concentrationH());

  // Compute the actual potential using the Nernst equation.
  const E_SRP = Terms.Volts(-2.19);
  const logQ = OH.mul(OH).div(Cl.mul(Cl)).naturalLogarithm();
  return E_SRP.sub(logQ.mul(UniversalGasConstant).mul(T).div(Faraday.mul(Terms.Number(2))));
}

// Compute the number of hydroxide ions (mol OH-) needed to raise the pH of a
// liter of water from startPH to endPH, with the contribution sorted by species.
function hydroxideContributors(T, S, DIC, startPH, endPH) {
  // Compute the concentrations (mol/l) of different species at the initial and
  // final pH. For each of these species, determine the number of moles of OH-
  // ions needed to effect the change.
  const contributors = {};

  // [H+] decreases as pH increases. Each removed H+ ion must have been
  // neutralized by an introduced OH- ion.
  const startH = startPH.concentrationH();
  const endH = endPH.concentrationH();
  assert(endH.lessThan(startH));
  contributors["H+"] = startH.sub(endH);

  // [OH-] increases as pH increases. As H+ ions are removed, additional OH-
  // ions are needed so that [H+] and [OH-] are in equilibrium according to the
  // disassociation constant Kw of water.
  const Kw = waterDisassociation();
  const startOH = Kw.div(startH);
  const endOH = Kw.div(endH);
  assert(startOH.lessThan(endOH));
  contributors["OH-"] = endOH.sub(startOH);

  // As pH increases, the concentrations of CO2 and HCO3 will decrease, adding
  // H+ ions to the water which need to be neutralized by OH-. Compute the total
  // number of hydrogen atoms in the carbonate species before and after the pH
  // change.
  const startCarbonate = carbonateConcentrations(T, S, DIC, startPH);
  const endCarbonate = carbonateConcentrations(T, S, DIC, endPH);

  // CO2 is effectively H2CO3 (carbonic acid), which it must form before it can
  // disassociate to HCO3-. Carbonic acid is diprotic and two H+ ions are
  // generated when it completely disassociates.
  contributors["CO2"] = startCarbonate.CO2.sub(endCarbonate.CO2).mul(Terms.Number(2));

  // HCO3- is the conjugate base of H2CO3 but is an acid itself, releasing an H+
  // ion when it disassociates which must be neutralized.
  contributors["HCO3-"] = startCarbonate.HCO3.sub(endCarbonate.HCO3);

  // Boric acid is also present in seawater, in low concentrations compared to
  // carbonate species.
  const Kb = boricAcidDisassociation();
  const totalBoric = seawaterConcentration("B(OH)3", S).add(seawaterConcentration("B(OH)4-", S));

  // As pH increases, [B(OH)4-] increases. Each such ion needs an OH- ion to
  // form from B(OH)3.
  //
  // Kb = [H+][B(OH)4-]/[B(OH)3]
  // Kb = [H+][B(OH)4-]/([TotalBoric] - [B(OH)4-])
  // Kb[TotalBoric] - Kb[B(OH)4-] = [H+][B(OH)4-]
  // Kb[TotalBoric] = [H+][B(OH)4-] + Kb[B(OH)4-]
  // Kb[TotalBoric] = ([H+] + Kb)[B(OH)4-]
  // [B(OH)4-] = Kb[TotalBoric] / ([H+] + Kb)
  const startBorate = Kb.mul(totalBoric).div(startH.add(Kb));
  const endBorate = Kb.mul(totalBoric).div(endH.add(Kb));
  assert(startBorate.lessThan(endBorate));
  contributors["B(OH)4-"] = endBorate.sub(startBorate);

  // Sulfuric acid is also present in seawater, with sulfate concentrations about
  // 10 times greater than DIC. This should be accounted for as well, though it
  // does not have as large a buffering effect because H2SO4 completely disassociates
  // and HSO4- is also a strong acid and mostly disassociates.
  //
  // The calculations here are analogous to those for boric acid.
  const Ks = Terms.Molarity(1.2e-2);
  const totalSulfate = seawaterConcentration("SO4^{2-}", S);
  const startSulfate = Ks.mul(totalSulfate).div(startH.add(Ks));
  const endSulfate = Ks.mul(totalSulfate).div(endH.add(Ks));
  assert(startSulfate.lessThan(endSulfate));
  contributors["SO4^{2-}"] = endSulfate.sub(startSulfate);

  return contributors;
}

// Compute the total amount of hydroxide ions needed between the various contributors.
function hydroxideRequirement(T, S, DIC, startPH, endPH) {
  const contributors = hydroxideContributors(T, S, DIC, startPH, endPH);
  let rv = Terms.Molarity(0);
  for (const v of Object.values(contributors)) {
    rv = rv.add(v);
  }
  return rv;
}

// Return how many amp-seconds are required to raise a volume V from startPH to endPH.
function electrolysisRequirement(T, S, DIC, V, startPH, endPH) {
  const requirement = hydroxideRequirement(T, S, DIC, startPH, endPH);
  const OH = requirement.mul(V);

  // Each amp will move one coulomb of electrons per second. Each electron will
  // reduce one H2O molecule to OH-.
  return OH.mul(Faraday);
}

// Compute the theoretical limit for the amount of seawater that can have its pH
// raised from startPH to endPH using electrolysis, using a power source W.
function electrolysisLimit(W, T, S, DIC, startPH, endPH) {
  // Compute the maximum amperage possible with W watts. In practice the voltage
  // required will be greater than that produced by electrolysisPotential().
  // The actual voltage depends on the environment and will need to be
  // determined experimentally.
  const potential = electrolysisPotential(T, S, startPH).negate();
  const amps = W.div(potential);

  const V = Terms.Liters(1);
  const requirement = electrolysisRequirement(T, S, DIC, V, startPH, endPH);
  return V.mul(amps).div(requirement);
}

module.exports = {
  electrolysisLimit,
  electrolysisPotential,
  electrolysisRequirement,
  hydroxideContributors,
  hydroxideRequirement
};
