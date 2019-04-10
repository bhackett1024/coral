/* Copyright 2019 Brian Hackett. Released under the MIT license. */

const { electrolysisIonMovement } = require("../migration");
const {
  expect,
  Temp_2100, Salinity_2100, DIC_2100, pH_2100,
  pH_Target, Chlorine_HalfLife,
  Cell_Amps, Cell_PrimaryCompartment, Cell_Outflow,
  Cell_InterfaceSize, Cell_InterfaceCurrent
} = require("./utils");

// https://coral.wiki/wiki/index.php?title=Electrolysis
//
// Compute alkalization losses from ion migration and diffusion in the example
// electrolysis cell, and the effect of chlorine diffusing into the cathode.
const movement = electrolysisIonMovement({
  TC: Temp_2100, S: Salinity_2100, DIC: DIC_2100,
  amps: Cell_Amps,
  halfLife: Chlorine_HalfLife,
  volume: Cell_PrimaryCompartment,
  outflow: Cell_Outflow,
  interface: Cell_InterfaceSize,
  current: Cell_InterfaceCurrent,
  startPH: pH_2100, endPH: pH_Target
});
expect(movement.migrationLoss + movement.diffusionLoss, 0.004190706572475512);
expect(movement.chlorineDilution, 94.47230599266902);
