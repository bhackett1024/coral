/* Copyright 2019 Brian Hackett. Released under the MIT license. */

const { carbonateConcentrations } = require("./carbonate");

// Return how much OH- is needed to neutralize any acid introduced due to flux
// of CO2 from the atmosphere into the water after raising the water pH from
// startPH to endPH. Result is in mol m^-2 h^-1.
function neutralizeRequirement(TC, S, DIC, startPH, endPH, windSpeed) {
  const startCarbonate = carbonateConcentrations(TC, S, DIC, startPH);
  const endCarbonate = carbonateConcentrations(TC, S, DIC, endPH);

  // Atmospheric flux depends on the CO2 concentration in the atmosphere and its
  // solubility in water, which can be derived to find an imbalance in CO2
  // concentration vs. the equilibrium.
  //
  // We can simplify this by assuming that the atmosphere is in equilibrium with
  // the amount of CO2 in the water at startPH. The CO2 imbalance is then the
  // difference in CO2 concentrations between the two different values of pH.
  const imbalance = (startCarbonate.CO2 - endCarbonate.CO2) / 1000; // mol cm^-3

  // The rate at which atmospheric CO2 will enter the water to correct this
  // imbalance depends on the sea state. The windier is, the more turbulent the
  // water surface will be, and the faster the rate. The calculations below are
  // based on the following reference to account for wind effects:
  //
  // Peter S. Liss, Liliane Merlivat
  // Air-Sea Gas Exchange Rates: Introduction and Synthesis
  // The Role of Air-Sea Exchange in Geochemical Cycling 113-127, 1986 (D Reidel, Dordrecht)
  //
  // The rate is the product of the concentration imbalance and a transfer
  // velocity term kw for how fast CO2 will transfer in the water at a given
  // wind speed (there is also a transfer velocity term for air, but for CO2
  // the term for water is dominant, and it is conservative to ignore the one
  // for air).

  // Use the formulae for calculating kw as given in the above reference,
  // normalized to a Schmidt number of 600 (CO2 at 20C).
  let kw; // cm/h
  if (windSpeed < 3.6) {
    kw = 0.17 * windSpeed;
  } else if (windSpeed < 13) {
    kw = 2.85 * windSpeed - 9.65;
  } else {
    kw = 5.9 * windSpeed - 49.3;
  }

  // Correct to the Schmidt number of CO2 at 30C.
  if (windSpeed < 3.6) {
    kw = kw * Math.pow(360, -2/3) / Math.pow(600, -2/3);
  } else {
    kw = kw * Math.pow(360, -0.5) / Math.pow(600, -.5);
  }

  const CO2flux = imbalance * kw * 100 * 100; // mol m^-2 h^-1

  // The CO2 entering the water will disassociate according to the carbonate
  // concentrations at the pH the water was raised to. In order to keep the pH
  // stable, we need OH to neutralize the H+ ion added when H2CO3 disassociates
  // to HCO3, and the two H+ ions added when H2CO3 fully disassociates to CO3.
  return CO2flux * (endCarbonate.HCO3 / DIC + endCarbonate.CO3 * 2 / DIC);
}

module.exports = { neutralizeRequirement };
