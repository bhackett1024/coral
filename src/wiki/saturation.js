/* Copyright 2019 Brian Hackett. Released under the MIT license. */

"use strict";

const { aragoniteSaturation } = require("../saturation");
const {
  expect,
  Temp_2010, Temp_2100,
  Salinity_2010, Salinity_2100,
  DIC_2010, DIC_2100,
  pH_2010, pH_2100
} = require("./utils");

// https://coral.wiki/wiki/index.php?title=Climate
//
// Compute aragonite saturation states around Moorea in 2006-2010 and 2096-2100.
expect(aragoniteSaturation(Temp_2010, Salinity_2010, DIC_2010, pH_2010),
       4.094835416162627);
expect(aragoniteSaturation(Temp_2100, Salinity_2100, DIC_2100, pH_2100),
       2.4610159113311125);
