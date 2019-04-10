/* Copyright 2019 Brian Hackett. Released under the MIT license. */

const { electrolysisPotential } = require("../electrolysis");
const { expect, Temp_2100, Salinity_2100, DIC_2100, pH_2100 } = require("./utils");

// https://coral.wiki/wiki/index.php?title=Electrolysis
//
// Calculate the theoretical potential needed for seawater electrolysis.
expect(electrolysisPotential(Temp_2100, Salinity_2100, pH_2100),
       -1.8304110483470388);
