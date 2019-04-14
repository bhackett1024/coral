/* Copyright 2019 Brian Hackett. Released under the MIT license. */

"use strict";

const { carbonateConcentrations, density, densityH2O } = require("./carbonate");
const { hydroxideRequirement } = require("./electrolysis");
const { steadyStateAmount } = require("./chlorine");
const { Units, Terms } = require("./units");
const { Species, SeawaterIons, waterDisassociation, Faraday } = require("./constants");

function solveForParameter(start, end, fn) {
  for (var i = 0; i < 30; i++) {
    const middle = start.add(end).div(Terms.Number(2));
    const value = fn(middle);
    if (value.isNegative()) {
      end = middle;
    } else {
      start = middle;
    }
  }
  return start;
}

function integrate(start, end, fn) {
  // How many points to include when integrating.
  const count = 10000;

  const diff = end - start;
  let total = 0;
  for (let i = 0; i < count; i++) {
    const x = start + diff * i / count;
    total += fn(x) / count;
  }
  return total;
}

// Compute the amount of a species that has diffused across an interface by the
// given time, in mol/m^2 of interface area.
//
// concentration (mol/m^3): fixed concentration at the interface
// coefficient (m^2/s): diffusion coefficient for the species
// time (s): elapsed time

// This implementation uses separation of variables. While not hard to determine,
// the resulting integrals and sums have to be approximated, making this very
// computationally expensive to determine accurately and vulnerable to roundoff
// errors.
function diffusion1(interfaceConcentration, coefficient, time) {
  interfaceConcentration = interfaceConcentration.normalize(Units.MolesPerCubicMeter);
  coefficient = coefficient.normalize(Units.SquareMetersPerSecond);
  time = time.normalize(Units.Seconds);

  // Diffusion in one dimension proceeds according to the PDE:
  //
  // dC(x,t)/dt = D * d^2C(x,t)/dx^2
  //
  // Instead of an infinite solution, consider things along a 1m gradient.
  // The interface concentration is Ci, and the concentration at the end of
  // the gradient is fixed at zero.
  //
  // Boundary conditions:
  //
  // C(0,t) = Ci  for all t
  // C(1,t) = 0   for all t
  //
  // Initial conditions:
  //
  // C(x,0) = 0   for x > 0
  //
  // The steady state of this will ultimately be the function
  // C(x,infinity) = Ci(1-x)
  //
  // To solve this using separation of variables, we need to transform the
  // conditions so that the boundary conditions are zero:
  //
  // C'(x) = C(x,t) - Ci(1-x)  or  C(x,t) = C'(x) + Ci(1-x)
  //
  // Boundary conditions:
  //
  // C'(0,t) = C'(1,t) = 0  for all t
  //
  // Initial conditions:
  //
  // C'(x,0) = -Ci(1-x)
  //
  // The solution for a zero-boundary PDE like this is below.
  //
  // Stanley J. Farlow
  // Partial Differential Equations for Scientists and Engineers (Dover Publications Inc, 1993)
  // page 39
  //
  // C'(x,t) = Sum_n=1..inf a_n exp(-t(n pi sqrt(D))^2) sin(n pi x)
  // where a_n = 2 * Integral_0..1 C'(x,0) sin(n pi x) dx

  // Number of sums to perform.
  const numsums = 10000;

  // Compute a_n coefficients.
  const coefficients = [];
  for (let n = 1; n <= numsums; n++) {
    coefficients[n] = 2 * integrate(0, 1, x => {
      return -interfaceConcentration * (1 - x) * Math.sin(n * Math.PI * x);
    });
  }

  function Cprime(x, time) {
    let total = 0;
    for (let n = 1; n <= numsums; n++) {
      const powerTerm = -time * Math.pow(n * Math.PI * Math.sqrt(coefficient), 2);
      total += coefficients[n] * Math.exp(powerTerm) * Math.sin(n * Math.PI * x);
    }
    return total;
  }

  // Compute Integral_0..1 C(x,t) dx
  const rv = integrate(0, 1, x => {
    return Cprime(x, time) + interfaceConcentration * (1 - x);
  });
  return Terms.MolesPerSquareMeter(rv);
}

// This implementation uses Laplace transforms. The solution is more complicated
// to determine but it is much simpler than diffusion1() and can be computed exactly.
function diffusion2(interfaceConcentration, coefficient, time) {
  // As for diffusion1(), diffusion in one dimension proceeds according to the PDE:
  //
  // dC(x,t)/dt = D * d^2C(x,t)/dx^2
  //
  // Boundary conditions:
  //
  // C(0,t) = Ci  for all t
  // lim(x->infinity) C(x,t) = 0   for all t
  //
  // Initial conditions:
  //
  // C(x,0) = 0   for x > 0
  //
  // Compute the Laplace transform of the PDE wrt t:
  //
  // L{dC(x,t)/dt} = L{D * d^2C(x,t)/dx^2}
  // s * LC(x,s) - C(x,0) = D * d^2LC(x,s)/dx^2
  // d^2LC(x,s)/dx^2 - s/D * LC(x,s) = 0
  //
  // This is an ODE, with the below solution.
  // ref: http://eqworld.ipmnet.ru/en/solutions/ode/ode0201.pdf
  //
  // LC(x,s) = c1 * sinh(x*sqrt(s/D)) + c2 * cosh(x*sqrt(s/D))
  //
  // Do some transformations to eliminate sinh/cosh, using identities for
  // sinh/cosh.
  // ref: http://math2.org/math/trig/hyperbolics.htm
  //
  // with A = x*sqrt(s/D)
  // LC(x,s) = c1 * sinh(A) + c2 * cosh(A)
  // LC(x,s) = c1 * (exp(A) - exp(-A)) / 2 + c2 * (exp(A) + exp(-A)) / 2
  // LC(x,s) = c3 * (exp(A) - exp(-A)) + c4 * (exp(A) + exp(-A))
  // LC(x,s) = c3 exp(A) - c3 exp(-A) + c4 exp(A) + c4 exp(-A)
  // LC(x,s) = (c3 + c4) exp(A) + (c4 - c3)exp(-A)
  // LC(x,s) = c5 exp(A) + c6 exp(-A)
  //
  // Using a fresh pool of constants:
  //
  // LC(x,s) = c1 * exp(x*sqrt(s/D)) + c2 * exp(-x*sqrt(s/D))
  //
  // Compute the Laplace transform of the boundary conditions:
  //
  // LC(0,s) = L{Ci}
  // LC(0,s) = Ci/s
  //
  // LC(inf,s) = L{0}
  // LC(inf,s) = 0
  //
  // Substituting to find c1/c2:
  //
  // c1 * exp(0*sqrt(s/D)) + c2 * exp(-0*sqrt(s/D)) = Ci/s
  // c1 + c2 = Ci/s
  //
  // c1 * exp(inf*sqrt(s/D)) + c2 * exp(-inf*sqrt(s/D)) = 0
  //
  // c1 = 0
  // c2 = Ci/s
  // LC(x,s) = Ci/s * exp(-x*sqrt(s/D))
  //
  // Compute the inverse transform to get back to C(x,t). First, transform the
  // equation so that it fits into a transformation table.
  //
  // L^-1{LC(x,s)} = L^-1{Ci/s * exp(-x*sqrt(s/D))}
  // C(x,t) = Ci * L^-1{1/s * exp(-sqrt(s*(x^2/D)))}
  //
  // Using an entry in the following table, we get the equation below.
  // ref: http://eqworld.ipmnet.ru/en/auxiliary/inttrans/LapInv5.pdf
  //
  // C(x,t) = Ci * erfc(sqrt(x^2/D)/(2*sqrt(t)))
  // C(x,t) = Ci * erfc(sqrt(x^2)*sqrt(1/D)/(2*sqrt(t)))
  // C(x,t) = Ci * erfc(x/(2*sqrt(D*t)))
  //
  // We can compute the exact integral of this function:
  //
  // Integral_0..infinity erfc(ax) dx = 1/(a*sqrt(pi))
  //
  // reference:
  //
  // Edward W. Ng and Murray Geller
  // A Table of Integrals of the Error Functions
  // Journal of Research of the National Bureau of Standards --
  //   B. Mathematical Sciences Vol. 73B, No. 1, January-March 1969
  //
  // Therefore:
  //
  // Integral_0..infinity C(x,t) dx = Ci / (sqrt(pi) * 1/(2*sqrt(D*t)))
  // Integral_0..infinity C(x,t) dx = Ci * 2 * sqrt(D*t) / sqrt(pi)

  const factor = Terms.Number(2 / Math.sqrt(Math.PI));
  return interfaceConcentration.mul(factor).mul(coefficient.mul(time).sqrt());
}

// Compute the approximate contribution of an ion to the conductivity of a
// solution, based on its molarity and limiting equivalent conductivity.
// The constants for limiting conductivity are based on extrapolation to a
// solution where the ion has infinite dilution. At higher concentrations
// this method will overapproximate the conductivity, because ions in
// concentrated solutions interact with each other and reduce the overall
// conductivity of the solution. There are methods to calculate conductivity
// which account for this but these depend on additional constants determined
// experimentally which are hard to find values for.
function ionConductivity(name, concentration) {
  // Ignore boric acid species (unknown conductivity, low concentration).
  if (/B\(OH\)/.test(name)) {
    return Terms.SiemensPerMeter(0);
  }
  const { conductivity, charge } = Species[name];
  return concentration.mul(conductivity).mul(charge.abs());
}

let seawaterConductivity = Terms.SiemensPerMeter(0);
for (const [name, concentration] of Object.entries(SeawaterIons)) {
  seawaterConductivity = seawaterConductivity.add(ionConductivity(name, concentration));
}
// > console.log(seawaterConductivity.normalize(SiemensPerMeter));
// 8.212238499999998
//
// Seawater is actually 5.3 S/m. The calculated value is 55% higher due to the
// overapproximation made in ionConductivity().

// Calculate movement due to migration and diffusion of ions and other species
// when an anode and cathode compartment are connected together.
//
// amps (A): How much current is moving between the compartments.
// halfLife (s): Half-life of chlorine compounds.
// volume (l): Size of the anode compartment.
// outflow (l/s): Exchange between the anode compartment and environment.
// interfaceSize (cm^2): Area of the interface between the compartments, across which
//   species can flow.
// current (cm/s): Speed of convection in both compartments along the interface.
// startPH (pH): Environmental pH.
// endPH (pH): pH which the cathode compartment is being raised to.
function electrolysisIonMovement({ T, S, DIC,
                                   amps, halfLife, volume, outflow, interfaceSize, current,
                                   startPH, endPH }) {
  // Here are the species that can have different concentrations in the two
  // compartments, and whose movements we need to account for:
  //
  // H+/OH-
  // CO2/HCO3-/CO3^{2-}
  // Cl2/HOCl/OCl-

  // How much charge moves between compartments.
  const chargeRate = amps.div(Faraday); // mol/s

  // Compute the steady state pH of the anode compartment: enough H+ is produced
  // each second to lower the inflow of water (equal to the outflow) from the
  // environmental pH to that of the steady state.
  //
  // The amount of H+ produced at the anode should equal the amount of OH-
  // produced at the cathode, which also equals the rate of charge transfer
  // between the compartments. If the anode produces O2 then this happens
  // immediately, while if the anode produces Cl2 then the H+ forms as the
  // Cl2 disassociates.
  const anodePH = solveForParameter(Terms.pH(1), startPH, (pH) => {
    const requirement = hydroxideRequirement(T, S, DIC, pH, startPH);
    return requirement.mul(outflow).sub(chargeRate);
  });

  // Concentrations in the cathode compartment.
  const cathode = carbonateConcentrations(T, S, DIC, endPH);
  const cathodeH = endPH.concentrationH();
  const cathodeOH = waterDisassociation().div(cathodeH);

  // Concentrations in the anode compartment.
  const anode = carbonateConcentrations(T, S, DIC, anodePH);
  const anodeH = anodePH.concentrationH();
  const anodeOH = waterDisassociation().div(anodeH);

  // Compute how many ions migrate across the interface, in mol/s.
  function migration(name, concentration) {
    // The transport number for an ion is the fraction of current which it
    // carries. In general these should be determined experimentally, but they can
    // be approximated using the fraction of their equivalent ionic conductivity
    // out of the total across all ions in solution:
    //
    // Allen J. Bard, Larry R. Faulkner
    // Electrochemical methods: fundamentals and applications, 2nd ed. (John Wiley & Sons Inc, 2001), page 67
    const { charge } = Species[name];
    const transport = ionConductivity(name, concentration).div(seawaterConductivity);
    return transport.mul(chargeRate).div(charge.abs());
  }

  // The anode consumes Cl- and/or produces H+, while the cathode produces OH-.
  // This creates a charge imbalance, and ions will migrate between the
  // compartments to balance out the charge. Negatively charged ions will move
  // to the anode, and positively charged ions will move to the cathode.
  // Migration that can occur:
  //
  // Cathode -> Anode: OH-, HCO3-, CO3^{2-}
  // Anode -> Cathode: H+
  //
  // Note that we assume the concentration of chlorine compounds in the cathode
  // compartment is zero, so none will migrate to the anode (any such migration
  // would be beneficial, in any case).

  const migrationOH = migration("OH-", cathodeOH);
  const migrationHCO3 = migration("HCO3-", cathode.HCO3);
  const migrationCO3 = migration("CO3^{2-}", cathode.CO3);
  const migrationH = migration("H+", anodeH);

  // All of these migrating ions have the effect of lowering the pH in the
  // cathode compartment, undoing the work done at the cathode to alkalize the
  // water.
  //
  // Since every electron transferred generates OH- at the cathode, the work
  // lost due to migration is the same as the fraction of current which
  // transports these ions.
  const migrationTotal =
    migrationOH.add(migrationHCO3).add(migrationCO3.cmul(2)).add(migrationH);
  const migrationLoss = migrationTotal.div(chargeRate).number();

  // A species will diffuse from a compartment where it is at higher concentration
  // to the compartment where it is at a lower concentration. The rate of
  // diffusion depends on several factors, and will speed up under these
  // conditions:
  //
  // - When the concentration difference is greater.
  // - When the species has a higher diffusion coefficient, and is more mobile
  //   in the water.
  // - When convection around the interface between the compartments is stronger.
  //   If the water in both compartments is still, the species will build up
  //   around the interface, flattening out the concentration gradient and
  //   reducing the rate of diffusion between the compartments. As a current
  //   flows, lower concentration water will be replenished at the interface and
  //   increase the overall rate of diffusion.
  //
  // Simultaneously modeling both diffusion and convection around the interface
  // is tricky, so we simplify things. Compute the time taken for water to sweep
  // past the interface, refreshing the water at the interface with that of the
  // well mixed compartment as a whole. Calculate the amount of diffusion that
  // will happen in this time, and multiply this by the refresh frequency with
  // which the water will refresh to get the rate of diffusion.

  // The time taken for water to sweep past the interface. We assume the interface
  // region is square.
  const refreshTime = interfaceSize.sqrt().div(current); // s

  // Compute the rate of diffusion across the interface, in mol/s.
  //
  // concentrationDifference (mol/l)
  // coefficient (m^2/s)
  function diffusion(name, concentrationDifference) {
    if (concentrationDifference.isNegative()) {
      return Terms.MolesPerSecond(0);
    }
    const { diffusion: coefficient } = Species[name];

    // Change this to diffusion1 to compare implementations.
    const rate = diffusion2(concentrationDifference, coefficient, refreshTime); // mol/m^2
    return rate.mul(interfaceSize).div(refreshTime); // mol/s
  }

  // Diffusing species we have to worry about:
  //
  // Cathode -> Anode: OH-, HCO3-, CO3^{2-}
  // Anode -> Cathode: H+, CO2, Cl2/HOCl/OCl-

  const diffusionH = diffusion("H+", anodeH.sub(cathodeH));
  const diffusionOH = diffusion("OH-", cathodeOH.sub(anodeOH));
  const diffusionCO2 = diffusion("CO2", anode.CO2.sub(cathode.CO2));
  const diffusionHCO3 = diffusion("HCO3-", cathode.HCO3.sub(anode.HCO3));
  const diffusionCO3 = diffusion("CO3^{2-}", cathode.CO3.sub(anode.CO3));

  // As with migration, any diffusion of these ions represents work lost while
  // alkalizing water in the compartment.
  const diffusionTotal = diffusionH.add(diffusionOH).add(diffusionCO2.cmul(2)).add(diffusionHCO3).add(diffusionCO3.cmul(2));
  const diffusionLoss = diffusionTotal.div(chargeRate).number();

  // Chlorine compounds do not affect pH, but we need to keep track of them to
  // manage diffusion of these toxic compounds into the cathode compartment.
  // While chlorine could be either Cl2, HOCl, or OCl-, we simplify things by
  // treating them all as OCl-. This is equivalent to the more general problem
  // except for effects from different diffusion coefficients for the species
  // (none of which we know precisely anyways).

  // For every two electrons conducted, at most one molecule of Cl2 will be
  // produced. The rate of active chlorine atom generation is then at most the
  // rate of conducted electrons.
  const anodeCl = steadyStateAmount(chargeRate, halfLife).div(volume);
  const diffusionCl = diffusion("OCl-", anodeCl); // mol/s

  // The diffusion of chlorine into the cathode compartment is a point source.
  // Compute how much seawater is needed to dilute this chlorine to PNEC concentration.
  const chlorineWeight = diffusionCl.mul(Terms.GramsPerMole(35.45)); // g/s
  const chlorineRelease = steadyStateAmount(chlorineWeight, halfLife);

  return { migrationLoss, diffusionLoss, chlorineRelease };
}

module.exports = { electrolysisIonMovement };
