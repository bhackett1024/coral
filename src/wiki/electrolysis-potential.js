/* Copyright 2019 Brian Hackett. Released under the MIT license. */

"use strict";

const { electrolysisPotential } = require("../electrolysis");
const { Units } = require("../units");
const { expect, Temp_2100, Salinity_2100, DIC_2100, pH_2100 } = require("./utils");

// https://coral.wiki/wiki/index.php?title=Electrolysis
//
// Calculate the theoretical potential needed for seawater electrolysis.
expect(electrolysisPotential(Temp_2100, Salinity_2100, pH_2100).normalize(Units.Volts),
       -1.829128647204989);
