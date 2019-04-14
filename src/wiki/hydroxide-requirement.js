/* Copyright 2019 Brian Hackett. Released under the MIT license. */

"use strict";

const { hydroxideContributors } = require("../electrolysis");
const { Units, Terms } = require("../units");
const {
  expect,
  Temp_2100, Salinity_2100, DIC_2100, pH_2100,
  pH_Target
} = require("./utils");

// https://coral.wiki/wiki/index.php?title=Electrolysis
//
// Calculate the number of moles of hydroxide ions needed to alkalize one liter.
const contributors =
  hydroxideContributors(Temp_2100, Salinity_2100, DIC_2100, pH_2100, pH_Target);

let requirement = Terms.Molarity(0);
for (const v of Object.values(contributors)) {
  requirement = requirement.add(v);
}
expect(requirement.normalize(Units.Molarity), 0.0002681740304786918);

// Looking at the contributions of the various ions in hydroxideContributors()
// gives a good sense of the buffering that happens in seawater and makes it
// resistant to changes in pH.
const contributorFractions = {};
for (const [name,v] of Object.entries(contributors)) {
  contributorFractions[name] = v.div(requirement).number();
}
expect(contributorFractions, {
  "H+": 0.00004127332845367099,
  "OH-": 0.003764172292604455,
  "CO2": 0.11048097382654368,
  "HCO3-": 0.7983541697869926,
  "B(OH)4-": 0.08725600506722427,
  "SO4^{2-}": 0.00010340569818143602
});
