/* Copyright 2019 Brian Hackett. Released under the MIT license. */

const { neutralizeRequirement } = require("../exchange");
const { Units, Terms } = require("../units");
const {
  expect,
  Temp_2100, Salinity_2100, DIC_2100, pH_2100,
  pH_Target
} = require("./utils");

// https://coral.wiki/wiki/index.php?title=Electrolysis
//
// Compute how much hydroxide is required to neutralize the acid introduced
// due to flux of CO2 from the atmosphere into the water after raising the
// water's pH with an electrolysis cell (mol m^-2 h^-1). Compute this for a
// wind speed of 6 m/s.
expect(neutralizeRequirement(Temp_2100, Salinity_2100, DIC_2100, pH_2100, pH_Target,
                             Terms.MetersPerSecond(6)).normalize(Units.MolesPerSquareMeterHour),
       0.001666281067108141);
