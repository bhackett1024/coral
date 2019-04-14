/* Copyright 2019 Brian Hackett. Released under the MIT license. */

const { design1, design2, design3 } = require("../chlorine");
const { Units } = require("../units");
const {
  expect,
  Chlorine_HalfLife, Chlorine_PNEC,
  Cell_ChlorineRate, Cell_PrimaryCompartment, Cell_SecondaryCompartment, Cell_Outflow
} = require("./utils");

// https://coral.wiki/wiki/index.php?title=Electrolysis
//
// Compute the amounts of chlorine which escapes to the environment with each
// cell design (g).
const Chlorine1 = design1(Cell_ChlorineRate, Chlorine_HalfLife).normalize(Units.Grams);
const Chlorine2 = design2(Cell_ChlorineRate, Chlorine_HalfLife,
                          Cell_PrimaryCompartment, Cell_Outflow).normalize(Units.Grams);
const Chlorine3 = design3(Cell_ChlorineRate, Chlorine_HalfLife,
                          Cell_PrimaryCompartment, Cell_Outflow,
                          Cell_SecondaryCompartment).normalize(Units.Grams);
expect(Chlorine1, 63.36316619584327);
expect(Chlorine2, 4.875398701107535);
expect(Chlorine3, 0.06405789516709288);

// Compare the amounts of chlorine produced by the designs:
expect(Chlorine2 / Chlorine1, 0.07694373551407804);
expect(Chlorine3 / Chlorine2, 0.013139006488339296);
expect(Chlorine3 / Chlorine1, 0.0010109642401565342);

function chlorineToBleach(Cl) {
  return Cl / 0.48 / 0.0525;
}

// Compute the approximate amount of bleach (ml) with an equivalent amount of chlorine.
expect(chlorineToBleach(Chlorine1), 2514.411356977908);
expect(chlorineToBleach(Chlorine2), 193.46820242490222);
expect(chlorineToBleach(Chlorine3), 2.541979966948131);

// Compute the volume needed to dilute 0.064 g of chlorine to the PNEC concentration.
// 23m * 23m * 3m = 1587 m^3
expect(Chlorine3 / Chlorine_PNEC.normalize(Units.GramsPerLiter) / 1000, 1525.1879801688783); // m^3
