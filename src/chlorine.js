/* Copyright 2019 Brian Hackett. Released under the MIT license. */

// Compute the steady state amount of chlorine which exists due to a
// single source emitting releaseRate chlorine (g/s), with the specified
// halfLife (s) for chlorine.
function steadyStateChlorine(releaseRate, halfLife) {
  // The average lifetime of a species is its half life divided by ln(2), which
  // works out to 1.44 times its half life.
  const lifetime = halfLife / Math.log(2);

  // Multiply by the release rate to get the steady state amount of chlorine due
  // to this source.
  return releaseRate * lifetime;
}

// Compute the steady state amount of chlorine (g) in the environment with design #1:
// the anode is placed directly into the environment.
function design1(anodeRate, halfLife) {
  return steadyStateChlorine(anodeRate, halfLife);
}

// Compute the steady state amount of chlorine (g) in the environment with design #2:
// the anode is placed into a compartment with a fixed volume (l) and outflow
// rate (l/s). We assume the compartment is well mixed.
function design2(anodeRate, halfLife, volume, outflow) {
  // The steady state amount of chlorine in the compartment is the same as would
  // be produced by design #1, ignoring any chlorine lost to outflow. Ignoring
  // the outflow here is conservative to do, and since designs will have a very
  // small outflow compared to the compartment volume, this is a good approximation.
  const compartmentConcentration = design1(anodeRate, halfLife) / volume;

  // The amount of chlorine released to the environment depends on the outflow
  // rate and the concentration of chlorine in the compartment.
  const releaseRate = compartmentConcentration * outflow;

  return steadyStateChlorine(releaseRate, halfLife);
}

// Compute the steady state amount of chlorine (g) in the environment with design #3:
// the anode is placed into a compartment with a fixed primaryVolume (l) and outflow
// rate (l/s), which then goes through a passage of size secondaryVolume (l)
// before reaching the environment.
//
// We assume the primary compartment is well mixed, and that water passes through
// the secondary compartment in an orderly, laminar fashion. The effects of
// diffusion (chlorine species migrating from higher to lower concentration regions)
// are ignored.
function design3(anodeRate, halfLife, primaryVolume, outflow, secondaryVolume) {
  // How many half lifes has the chlorine emitted from the primary compartment
  // gone through before reaching the environment?
  const lifetimeCount = secondaryVolume / outflow / halfLife;

  // The fraction of chlorine which still remains by the time water reaches the
  // environment.
  const fraction = Math.pow(0.5, lifetimeCount);

  return fraction * design2(anodeRate, halfLife, primaryVolume, outflow);
}

module.exports = { design1, design2, design3 };
