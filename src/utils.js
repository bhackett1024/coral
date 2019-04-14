/* Copyright 2019 Brian Hackett. Released under the MIT license. */

"use strict";

function assert(v, why) {
  if (!v) {
    throw new Error(`Assertion failed: ${why}`);
  }
}

module.exports = { assert };
