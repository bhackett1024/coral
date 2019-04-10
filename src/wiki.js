/* Copyright 2019 Brian Hackett. Released under the MIT license. */

const { generateBjerrumData, density } = require("./carbonate");
const { aragoniteSaturation } = require("./saturation");

// This file checks that all data included on http://coral.wiki is in sync with
// the data produced by this repository.

function expect(actual, expected) {
  const sa = JSON.stringify(actual);
  const se = JSON.stringify(expected);
  if (sa != se) {
    throw new Error(`Mismatch! Expected ${se}, got ${sa}`);
  }
}

// https://coral.wiki/wiki/index.php?title=Carbonate-Script
//
// Computer Bjerrum plot data for carbonate concentrations as a function of pH.
expect(generateBjerrumData(28, 35, 4, 10, 20), {
  pH: [
    4,4.3,4.6,4.9,5.2,5.5,5.8,6.1,6.4,6.7,7,7.300000000000001,7.6,7.9,8.2,8.5,
    8.8,9.1,9.4,9.7,10
  ],
  CO2: [
    0.984521075526212,0.969583696419006,0.9410933224090265,0.8889700081106523,
    0.8004947579655667,0.6678369658227451,0.5018130225484696,0.3353017888106607,
    0.20155084443942506,0.11199784775270286,0.059111131571207697,0.03016884139968606,
    0.014996011673134002,0.0072358410873917985,0.0033421465853907114,
    0.0014422309477453872,0.0005649467721451165,0.00019708664968350493,
    0.00006147748027210752,0.000017571226527866345,0.000004741199864460002
  ],
  HCO3: [
    0.015478732301200735,0.030415550135701457,0.05890376620951796,0.1110190434169304,
    0.19946599336739915,0.3320326764218774,0.49779702797022,0.663660916176523,
    0.7959668739977576,0.8825108405413283,0.9293507313926659,0.9463875181866063,
    0.9386120851525606,0.9036482091746433,0.8327906985033544,0.7170429702231562,
    0.5604255200611997,0.39009248044595485,0.24278756027814646,0.13845621325410828,
    0.07454158708511233
  ],
  CO3: [
    1.9217258716646177e-7,7.534452924028192e-7,0.000002911381455470726,
    0.000010948472417309967,0.00003924866703411113,0.00013035775537761956,
    0.00038994948131036125,0.0010372950128162079,0.0024822815628172544,
    0.00549131170596871,0.011538137036126519,0.023443640413707827,0.04639190317430532,
    0.08911594973796505,0.1638671549112548,0.2815147988290984,0.4390095331666552,
    0.6097104329043617,0.7571509622415813,0.8615262155193639,0.9254536717150231
  ]
});

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

// https://coral.wiki/wiki/index.php?title=Climate
//
// Compute aragonite saturation states around Moorea in 2006-2010 and 2096-2100.
expect(aragoniteSaturation(Temp_2010, Salinity_2010, DIC_2010, pH_2010),
       4.094835416162626);
expect(aragoniteSaturation(Temp_2100, Salinity_2100, DIC_2100, pH_2100),
       2.461015911331112);
