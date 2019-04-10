/* Copyright 2019 Brian Hackett. Released under the MIT license. */

const { density } = require("../carbonate");

// Make sure that a computed value matches what the data in the wiki.
function expect(actual, expected) {
  const sa = JSON.stringify(actual);
  const se = JSON.stringify(expected);
  if (sa != se) {
    throw new Error(`Mismatch! Expected ${se}, got ${sa}`);
  }
}

// Get some constants for environmental conditions to consider while modeling.

// Get the average values of temperature, salinity, DIC, and pH for the grid
// coordinate (66,130) closest to Opunohu Bay on Moorea in 2006-2010 and 2096-2100.
//
// Data source:
//
// Dunne, John; John, Jasmin; Shevliakova, Elena; Stouffer, Ronald; Griffies, Stephen; Malyshev, Sergey; Milly, P.; Sentman, Lori; Adcroft, Alistair; Cooke, William; Dunne, Krista; Hallberg, Robert; Harrison, Matthew; Krasting, John; Levy, Hiram; Phillips, Peter; Samuels, Bonita; Spelman, Michael; Winton, Michael; Wittenberg, Andrew; Zadeh, Niki
// NOAA GFDL GFDL-ESM2M, rcp85 experiment output for CMIP5 AR5, served by ESGF
//
// Tracking IDs:
//
// Temperature: 23fc9a5a-1991-4f34-aae7-ccb51fb4a7cd, 9cfab26d-0591-4c50-83fc-76f4501d4771
// Salinity: 1d57bf2b-d5f7-4539-9b77-d6878ad161be, d8ccebdd-c8dd-46ed-a876-382a7ec16c61
// DIC: 2f917c29-567b-495d-812a-a70d8b758484, 2307f360-a117-4a6a-8541-c6099b5f5976
// pH: 1b93c61d-3a60-44a8-9677-902a1bb53944, 9c5ddc62-cdc4-42d8-bfea-5c7800031bc1
const Temp_2010 = 27.6, Temp_2100 = 29.5;
const Salinity_2010 = 35.4, Salinity_2100 = 35.4;
const DIC_2010 = 0.002037928 / density(27.6, 35.4), DIC_2100 = 0.002229376 / density(29.5, 35.4);
const pH_2010 = 8.08, pH_2100 = 7.76;

// Target pH we want to alkalize water to.
const pH_Target = 8.2;

module.exports = {
  expect,
  Temp_2010, Temp_2100,
  Salinity_2010, Salinity_2100,
  DIC_2010, DIC_2100,
  pH_2010, pH_2100,
  pH_Target
};
