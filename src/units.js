/* Copyright 2019 Brian Hackett. Released under the MIT license. */

"use strict";

const { assert } = require("./utils");

// Framework for keeping track of the unit associated with each term, allowing
// easy interconversion between different units, and making sure that the
// units used in a calculation match up with each other.

function Term(value, unit) {
  assert(typeof value == "number" && !Number.isNaN(value), "Bad value");
  if (unit instanceof Array) {
    // Internal constructor which supplies the primitive components directly.
    this.value = value;
    this.components = unit;
  } else if (unit.offset) {
    // Units offset from a base are only supported for terms with a single dimension.
    assert(!unit.scale);
    assert(unit.components.length == 1 &&
           !unit.components[0].power &&
           !unit.components[0].unit.components);

    this.value = value + unit.offset;
    this.components = [ { unit: unit.components[0].unit, power: 1 } ];
  } else {
    const { scale, components } = reduceUnit(unit);
    this.value = value * scale;
    this.components = components;
  }
}

Term.prototype = {
  normalize(unit) {
    const { scale, components } = reduceUnit(unit);
    assert(componentsMatch(this.components, components),
           `Normalize failed: ${JSON.stringify(this.components)} => ${unit.name}`);
    return this.value / scale;
  },

  div(v) {
    assert(v instanceof Term);
    const newComponents = cloneComponents(v.components);
    for (const entry of newComponents) {
      entry.power *= -1;
    }
    return this.mul(new Term(1 / v.value, newComponents));
  },

  mul(v) {
    assert(v instanceof Term);
    const newComponents = cloneComponents(this.components);
    mergeComponents(newComponents, v.components, 1);
    return new Term(this.value * v.value, newComponents);
  },

  cmul(n) {
    return this.mul(Terms.Number(n));
  },

  sqrt() {
    const newComponents = cloneComponents(this.components);
    for (const entry of newComponents) {
      assert(entry.power % 2 == 0);
      entry.power /= 2;
    }
    return new Term(Math.sqrt(this.value), newComponents);
  },

  add(v) {
    assert(v instanceof Term);
    assert(componentsMatch(this.components, v.components),
           `Components must match for add: ${JSON.stringify(this.components)} ${JSON.stringify(v.components)}`);
    return new Term(this.value + v.value, this.components);
  },

  sub(v) {
    return this.add(v.negate());
  },

  negate(v) {
    return new Term(-this.value, this.components);
  },

  lessThan(v) {
    assert(v instanceof Term);
    assert(componentsMatch(this.components, v.components),
           `Components must match for lessThan: ${JSON.stringify(this.components)} ${JSON.stringify(v.components)}`);
    return this.value < v.value;
  },

  isNegative() {
    return this.value < 0;
  },

  // Convert a dimensionless quantity to its number.
  number() {
    return this.normalize(Units.Number);
  },

  // Get the natural logarithm of a dimensionless quantity.
  naturalLogarithm() {
    return Terms.Number(Math.log(this.number()));
  },

  // Get the absolute value of a dimensionless quantity.
  abs() {
    return Terms.Number(Math.abs(this.number()));
  },

  // Get the concentration of H+ for a given pH.
  concentrationH() {
    assert(this.components.length == 1 &&
           this.components[0].power == 1 &&
           this.components[0].unit == Units.pH);
    return Terms.Molarity(Math.pow(10, -this.value));
  }
};

// Convert a unit to a factor of base unit components.
function reduceUnit(unit) {
  if (unit == Base.Number) {
    return { scale: 1, components: [] };
  }
  if (!unit.components) {
    return { scale: 1, components: [{ unit, power: 1 }] };
  }

  assert(!unit.offset);
  const rv = { scale: unit.scale || 1, components: [] };
  for (const { unit: componentUnit, power } of unit.components) {
    const { scale, components } = reduceUnit(componentUnit);
    rv.scale *= Math.pow(scale, power || 1);
    mergeComponents(rv.components, components, power || 1);
  }

  return rv;
}

function mergeComponents(components, newComponents, power) {
  assert(power);
  for (const { unit, power: componentPower } of newComponents) {
    assert(!unit.components && componentPower);
    const newPower = power * componentPower;
    mergeSingleComponent(components, unit, newPower);
  }
}

function mergeSingleComponent(components, unit, power) {
  assert(!unit.components && power);
  for (let i = 0; i < components.length; i++) {
    if (components[i].unit == unit) {
      components[i].power += power;
      if (components[i].power == 0) {
        components[i] = components[components.length - 1];
        components.pop();
      }
      return;
    }
  }
  components.push({ unit, power });
}

function cloneComponents(components) {
  const rv = [];
  for (const entry of components) {
    rv.push({ ...entry });
  }
  return rv;
}

function componentsMatch(acomponents, bcomponents) {
  if (acomponents.length != bcomponents.length) {
    return false;
  }
  for (const { unit, power } of acomponents) {
    if (!bcomponents.some(({ unit: bunit, power: bpower }) => {
      return unit == bunit && power == bpower;
    })) {
      return false;
    }
  }
  return true;
}

// All other units are expressed in terms of the base units.
const Base = {
  // Special dimensionless unit.
  Number: {
    name: ""
  },

  Kelvin: {
    name: "K"
  },

  Salinity: {
    name: "S"
  },

  pH: {
    name: "pH"
  },

  Moles: {
    name: "mol"
  },

  Grams: {
    name: "g"
  },

  Meters: {
    name: "m"
  },

  Seconds: {
    name: "s"
  },

  Coulombs: {
    name: "C"
  },

  Joules: {
    name: "J"
  }
};

// Approximate weight of a liter of seawater, in kilograms.
const SeawaterDensity = 1.025;

const Derived1 = {
  Celsius: {
    name: "°C",
    components: [ { unit: Base.Kelvin } ],
    offset: 273.15
  },

  Decimeters: {
    name: "dm",
    components: [ { unit: Base.Meters } ],
    scale: 0.1
  },

  Centimeters: {
    name: "cm",
    components: [ { unit: Base.Meters } ],
    scale: 0.01
  },

  GramsPerSecond: {
    name: "g/s",
    components: [
      { unit: Base.Grams },
      { unit: Base.Seconds, power: -1 }
    ]
  },

  GramsPerMole: {
    name: "g/mol",
    components: [
      { unit: Base.Grams },
      { unit: Base.Moles, power: -1 }
    ]
  },

  CubicMeters: {
    name: "m^3",
    components: [ { unit: Base.Meters, power: 3 } ]
  },

  MetersPerSecond: {
    name: "m/s",
    components: [
      { unit: Base.Meters },
      { unit: Base.Seconds, power: -1 }
    ]
  },

  SquareMetersPerSecond: {
    name: "m^2/s",
    components: [
      { unit: Base.Meters, power: 2 },
      { unit: Base.Seconds, power: -1 }
    ]
  },

  MolesPerSecond: {
    name: "mol/s",
    components: [
      { unit: Base.Moles },
      { unit: Base.Seconds, power: -1 }
    ]
  },

  MolesPerSquareMeter: {
    name: "mol/m^2",
    components: [
      { unit: Base.Moles },
      { unit: Base.Meters, power: -2 }
    ]
  },

  MolesPerCubicMeter: {
    name: "mol/m^3",
    components: [
      { unit: Base.Moles },
      { unit: Base.Meters, power: -3 }
    ]
  },

  Hours: {
    name: "h",
    components: [ { unit: Base.Seconds } ],
    scale: 1 / 3600
  },

  Amperes: {
    name: "A",
    components: [
      { unit: Base.Coulombs },
      { unit: Base.Seconds, power: -1 }
    ]
  },

  Volts: {
    name: "V",
    components: [
      { unit: Base.Joules },
      { unit: Base.Coulombs, power: -1 }
    ]
  },

  JoulesPerKelvinMole: {
    name: "J K^-1 mol^-1",
    components: [
      { unit: Base.Joules },
      { unit: Base.Kelvin, power: -1 },
      { unit: Base.Moles, power: -1 }
    ]
  },

  CoulombsPerMole: {
    name: "C/mol",
    components: [
      { unit: Base.Coulombs },
      { unit: Base.Moles, power: -1 }
    ]
  },
};

const Derived2 = {
  Liters: {
    name: "l",
    components: [ { unit: Derived1.Decimeters, power: 3 } ]
  },

  SquareCentimeters: {
    name: "cm^2",
    components: [ { unit: Derived1.Centimeters, power: 2 } ]
  },

  CentimetersPerSecond: {
    name: "cm/s",
    components: [
      { unit: Derived1.Centimeters },
      { unit: Base.Seconds, power: -1 }
    ]
  },

  CentimetersPerHour: {
    name: "cm/h",
    components: [
      { unit: Derived1.Centimeters },
      { unit: Derived1.Hours, power: -1 }
    ]
  },

  MolesPerSquareMeterHour: {
    name: "mol m^-2 h^-1",
    components: [
      { unit: Base.Moles },
      { unit: Base.Meters, power: -2 },
      { unit: Derived1.Hours, power: -1 }
    ]
  },

  Watts: {
    name: "W",
    components: [
      { unit: Derived1.Amperes },
      { unit: Derived1.Volts }
    ]
  },

  AmpSeconds: {
    name: "A s",
    components: [
      { unit: Derived1.Amperes },
      { unit: Base.Seconds }
    ]
  },

  Ohms: {
    name: "Ω",
    components: [
      { unit: Derived1.Volts },
      { unit: Derived1.Amperes, power: -1 }
    ]
  },
};

const Derived3 = {
  Molarity: {
    name: "mol/l",
    components: [
      { unit: Base.Moles },
      { unit: Derived2.Liters, power: -1 }
    ]
  },

  GramsPerLiter: {
    name: "g/l",
    components: [
      { unit: Base.Grams },
      { unit: Derived2.Liters, power: -1 }
    ]
  },

  LitersPerSecond: {
    name: "l/s",
    components: [
      { unit: Derived2.Liters },
      { unit: Base.Seconds, power: -1 }
    ]
  },

  MolesSquaredPerLiterSquared: {
    name: "mol^2 l^-2",
    components: [
      { unit: Base.Moles, power: 2 },
      { unit: Derived2.Liters, power: -2 }
    ]
  },

  Siemens: {
    name: "S",
    components: [
      { unit: Derived2.Ohms, power: -1 }
    ]
  },

  // A kilogram of seawater as a unit of volume.
  SeawaterKg: {
    name: "kg-sw",
    components: [ { unit: Derived2.Liters } ],
    scale: 1 / SeawaterDensity
  },
};

const Derived4 = {
  SiemensPerMeter: {
    name: "S/m",
    components: [
      { unit: Derived3.Siemens },
      { unit: Base.Meters, power: -1 }
    ]
  },

  SiemensSquareCentimersPerMole: {
    name: "S cm^2 mol^-1",
    components: [
      { unit: Derived3.Siemens },
      { unit: Derived1.Centimeters, power: 2 },
      { unit: Base.Moles, power: -1 }
    ]
  },

  MolesPerSeawaterKg: {
    name: "mol/kg-sw",
    components: [
      { unit: Base.Moles },
      { unit: Derived3.SeawaterKg, power: -1 }
    ]
  },
};

const Units = {
  ...Base,
  ...Derived1,
  ...Derived2,
  ...Derived3,
  ...Derived4
};

const Terms = {};
for (const [id,unit] of Object.entries(Units)) {
  Terms[id] = value => new Term(value, unit);
}

module.exports = { Units, Terms };
