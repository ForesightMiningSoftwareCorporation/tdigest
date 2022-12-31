//! T-Digest algorithm in Rust
//!
//! A data structure for approximating the [quantile
//! function](https://en.wikipedia.org/wiki/Quantile_function) of a sample
//! distribution. T-digests have bounded memory requirements, making them useful
//! for streaming data analysis.
//!
//! [t-digest
//! paper](https://github.com/tdunning/t-digest/blob/main/docs/t-digest-paper/histo.pdf)
//!
//! ## Example
//!
//! ```rust
//! use tdigest::TDigest;
//!
//! let t = TDigest::default();
//! let values: Vec<f64> = (1..=1_000_000).map(f64::from).collect();
//!
//! let t = t.merge_sorted(100, &values);
//!
//! let ans = t.estimate_quantile(0.99).unwrap();
//! let expected: f64 = 990_000.0;
//!
//! let percentage: f64 = (expected - ans).abs() / expected;
//! assert!(percentage < 0.01);
//! ```

mod cluster;
mod t_digest;

pub use cluster::*;
pub use t_digest::*;
