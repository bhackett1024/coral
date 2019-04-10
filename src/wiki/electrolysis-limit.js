/* Copyright 2019 Brian Hackett. Released under the MIT license. */

const { electrolysisLimit } = require("../electrolysis");
const { expect, Temp_2100, Salinity_2100, DIC_2100, pH_2100 } = require("./utils");

// https://coral.wiki/wiki/index.php?title=Electrolysis
//
// Determine the amount of water that can be alkalized by a 100W solar panel in
// conditions expected for the end of the century.
//
// Looking at the contributions of the various ions in hydroxideRequirement()
// gives a good sense of the buffering that happens in seawater and makes it
// resistant to changes in pH. Of the hydroxide ions required to change the
// pH to 8.2:
//
// 0.004% (1 in 24107) neutralizes an H+ ion.
// 0.37% (1 in 264) is needed to maintain the equilibrium between [H+] and [OH-].
// 91.3% is needed to neutralize H+ ions added as CO2 and HCO3- disassociate.
// 8.3% is needed to bring boric acid species into equilibrium.
expect(electrolysisLimit(100, Temp_2100, Salinity_2100, DIC_2100, pH_2100, 8.2),
       2.1220809101172438);
