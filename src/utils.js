/* Copyright 2019 Brian Hackett. Released under the MIT license. */

function assert(v) {
  if (!v) {
    throw new Error("Assertion failed!");
  }
}

module.exports = { assert };
