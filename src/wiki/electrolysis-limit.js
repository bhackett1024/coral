/* Copyright 2019 Brian Hackett. Released under the MIT license. */

const { electrolysisLimit } = require("../electrolysis");
const {
  expect,
  Temp_2100, Salinity_2100, DIC_2100, pH_2100,
  pH_Target
} = require("./utils");

// https://coral.wiki/wiki/index.php?title=Electrolysis
//
// Determine the amount of water that can be alkalized by a 100W solar panel in
// conditions expected for the end of the century.
expect(electrolysisLimit(100, Temp_2100, Salinity_2100, DIC_2100, pH_2100, pH_Target),
       2.1220809101172438);
