/* Copyright 2019 Brian Hackett. Released under the MIT license. */

const { assert } = require("./utils");

// Framework for keeping track of the unit associated with each term, allowing
// easy interconversion between different units, and making sure that the
// units used in a calculation match up with each other.

function Term(value, unit) {
  this.value = value;
  this.unit = unit;
}

Term.prototype = {
  normalize(unit) {
    assert(this.unit && unit);
    if (this.unit == unit) {
      return this.value;
    }
    if (this.unit.base == unit) {
      if (this.unit.offset) {
        return this.value + this.unit.offset;
      }
    }
    const factor = dimensionFactor(this.unit, unit);
    if (factor) {
      return this.value * factor;
    }
    console.error(`Normalize failed: ${this.unit.name} => ${unit.name}`);
    assert(false);
  },

  // Get the concentration of H+ for a given pH.
  concentrationH() {
    assert(this.unit == Units.pH);
    return Terms.Molarity(Math.pow(10, -this.value));
  }
};

function dimensionFactor(unit1, unit2) {
  assert(unit1 && unit2);

  if (unit1 == unit2) {
    return 1;
  }

  if (unit1.base == unit2) {
    return unit1.scale;
  }

  if (unit2.base == unit1) {
    return 1 / unit2.scale;
  }

  const { components: c1 } = unit1;
  const { components: c2 } = unit2;
  if (!c1 || !c2 || c1.length != c2.length) {
    return;
  }
  let factor = 1;
  for (let i = 0; i < c1.length; i++) {
    if (c1[i].power != c2[i].power) {
      return;
    }
    const dfactor = dimensionFactor(c1[i].unit, c2[i].unit);
    if (!dfactor) {
      return;
    }
    factor *= Math.pow(dfactor, c1[i].power || 1);
  }
  return factor;
}

const Base = {
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

  Liters: {
    name: "l"
  },

  Seconds: {
    name: "s"
  },

  Coulombs: {
    name: "C"
  },

  Volts: {
    name: "V"
  }
};

// Approximate weight of a liter of seawater, in kilograms.
const SeawaterDensity = 1.025;

const Scaled = {
  Celsius: {
    name: "°C",
    base: Base.Kelvin,
    offset: 273.15
  },

  Centimeter: {
    name: "cm",
    base: Base.Meters,
    scale: 0.01
  },

  // A kilogram of seawater as a unit of volume.
  SeawaterKg: {
    name: "kg-sw",
    base: Base.Liters,
    scale: 1 / SeawaterDensity
  }
};

const Derived = {
  SquareCentimeters: {
    name: "cm^2",
    components: [ { unit: Scaled.Centimeter, power: 2 } ]
  },

  Molarity: {
    name: "mol/l",
    components: [
      { unit: Base.Moles },
      { unit: Base.Liters, power: -1 }
    ]
  },

  MolesPerSeawaterKg: {
    name: "mol/kg-sw",
    components: [
      { unit: Base.Moles },
      { unit: Scaled.SeawaterKg, power: -1 }
    ]
  },

  GramsPerLiter: {
    name: "g/l",
    components: [
      { unit: Base.Grams },
      { unit: Base.Liters, power: -1 }
    ]
  },

  GramsPerSecond: {
    name: "g/s",
    components: [
      { unit: Base.Grams },
      { unit: Base.Seconds, power: -1 }
    ]
  },

  LitersPerSecond: {
    name: "l/s",
    components: [
      { unit: Base.Liters },
      { unit: Base.Seconds, power: -1 }
    ]
  },

  MetersPerSecond: {
    name: "m/s",
    components: [
      { unit: Base.Meters },
      { unit: Base.Seconds, power: -1 }
    ]
  },

  CentimetersPerSecond: {
    name: "cm/s",
    components: [
      { unit: Scaled.Centimeters },
      { unit: Base.Seconds, power: -1 }
    ]
  },

  MolesPerSquareMeterHour: {
    name: "mol m^-2 h^-1",
    components: [
      { unit: Base.Moles },
      { unit: Base.Meters, power: -2 },
      { unit: Scaled.Hours, power: -1 }
    ]
  },

  Amperes: {
    name: "A",
    components: [
      { unit: Base.Coulombs },
      { unit: Base.Seconds, power: -1 }
    ]
  },
};

const Derived2 = {
  Watts: {
    name: "W",
    components: [
      { unit: Derived.Amperes },
      { unit: Base.Volts }
    ]
  }
};

const Units = {
  ...Base,
  ...Scaled,
  ...Derived,
  ...Derived2
};

const Terms = {};
for (const [id,unit] of Object.entries(Units)) {
  Terms[id] = value => new Term(value, unit);
}

module.exports = { Units, Terms };
