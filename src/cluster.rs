use float_ord::FloatOrd;
use std::{
    cmp::Ordering,
    ops::{Add, AddAssign},
};

/// A cluster of samples represented as the mean of those samples.
#[derive(Debug, Clone, Default)]
pub struct Cluster {
    pub(crate) count: usize,
    pub(crate) sum: f64,
}

impl PartialEq for Cluster {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.sum == other.sum && self.count == other.count
    }
}
impl Eq for Cluster {}

impl PartialOrd for Cluster {
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Cluster {
    #[inline]
    fn cmp(&self, other: &Self) -> Ordering {
        FloatOrd(self.mean()).cmp(&FloatOrd(other.mean()))
    }
}

impl Add for Cluster {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        Self {
            count: self.count + rhs.count,
            sum: self.sum + rhs.sum,
        }
    }
}

impl AddAssign for Cluster {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        self.count += rhs.count;
        self.sum += rhs.sum;
    }
}

impl Cluster {
    #[inline]
    pub fn new(sum: f64, count: usize) -> Self {
        Self { sum, count }
    }

    #[inline]
    pub fn count(&self) -> usize {
        self.count
    }

    #[inline]
    pub fn sum(&self) -> f64 {
        self.sum
    }

    #[inline]
    pub fn mean(&self) -> f64 {
        self.sum / self.count as f64
    }
}
