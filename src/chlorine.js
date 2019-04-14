/* Copyright 2019 Brian Hackett. Released under the MIT license. */

"use strict";

const { Units, Terms } = require("./units");

// Compute the steady state amount of a substance N which exists due to a
// single source emitting the substance at releaseRate (N/s), with the specified
// half life.
function steadyStateAmount(releaseRate, halfLife) {
  // The average lifetime of a species is its half life divided by ln(2), which
  // works out to 1.44 times its half life.
  const lifetime = halfLife.div(Terms.Number(Math.log(2)));

  // Multiply by the release rate to get the steady state amount of the substance
  // due to this source.
  return releaseRate.mul(lifetime);
}

// Compute the steady state amount of chlorine (g) in the environment with design #1:
// the anode is placed directly into the environment.
function design1(anodeRate, halfLife) {
  return steadyStateAmount(anodeRate, halfLife);
}

// Compute the steady state amount of chlorine in the environment with design #2:
// the anode is placed into a compartment with a fixed volume and outflow rate.
// We assume the compartment is well mixed.
function design2(anodeRate, halfLife, volume, outflow) {
  // The steady state amount of chlorine in the compartment is the same as would
  // be produced by design #1, ignoring any chlorine lost to outflow. Ignoring
  // the outflow here is conservative to do, and since designs will have a very
  // small outflow compared to the compartment volume, this is a good approximation.
  const compartmentChlorine = design1(anodeRate, halfLife);
  const compartmentConcentration = compartmentChlorine.div(volume);

  // The amount of chlorine released to the environment depends on the outflow
  // rate and the concentration of chlorine in the compartment.
  const releaseRate = compartmentConcentration.mul(outflow);

  return steadyStateAmount(releaseRate, halfLife);
}

// Compute the steady state amount of chlorine in the environment with design #3:
// the anode is placed into a compartment with a fixed primaryVolume and outflow rate,
// which then goes through a passage of size secondaryVolume before reaching the
// environment.
//
// We assume the primary compartment is well mixed, and that water passes through
// the secondary compartment in an orderly, laminar fashion. The effects of
// diffusion (chlorine species migrating from higher to lower concentration regions)
// are ignored.
function design3(anodeRate, halfLife, primaryVolume, outflow, secondaryVolume) {
  // The steady state amount of chlorine produced by design #2.
  const chlorine = design2(anodeRate, halfLife, primaryVolume, outflow);

  // How many half lifes has the chlorine emitted from the primary compartment
  // gone through before reaching the environment?
  const lifetimeCount = secondaryVolume.div(outflow).div(halfLife);

  // The fraction of chlorine which still remains by the time water reaches the
  // environment.
  const fraction = Terms.Number(Math.pow(0.5, lifetimeCount.number()));

  return chlorine.mul(fraction);
}

module.exports = { design1, design2, design3, steadyStateAmount };
