// `clippy` is a code linting tool for improving code quality by catching
// common mistakes or strange code patterns. If the `clippy` feature is
// provided, it is enabled and all compiler warnings are prohibited.
#![cfg_attr(feature = "clippy", deny(warnings))]
#![cfg_attr(feature = "clippy", feature(plugin))]
#![cfg_attr(feature = "clippy", plugin(clippy))]
#![cfg_attr(feature = "clippy", allow(inline_always))]
#![cfg_attr(feature = "clippy", allow(too_many_arguments))]
#![cfg_attr(feature = "clippy", allow(unreadable_literal))]

// The compiler provides `test` (on nightly) for benchmarking tools, but
// it's hidden behind a feature flag. Enable it if we're testing.
#![cfg_attr(test, feature(test))]
#[cfg(test)]
extern crate test;

extern crate rand;

extern crate pairing;

mod bls12_381;
