/* Copyright 2019 Brian Hackett. Released under the MIT license. */

const { electrolysisLimit } = require("../electrolysis");
const { Units, Terms } = require("../units");
const {
  expect,
  Temp_2100, Salinity_2100, DIC_2100, pH_2100,
  pH_Target
} = require("./utils");

// https://coral.wiki/wiki/index.php?title=Electrolysis
//
// Determine the amount of water that can be alkalized by a 100W solar panel in
// conditions expected for the end of the century.
expect(electrolysisLimit(Terms.Watts(100),
                         Temp_2100, Salinity_2100, DIC_2100, pH_2100,
                         pH_Target).normalize(Units.LitersPerSecond),
       2.146058637396675);
