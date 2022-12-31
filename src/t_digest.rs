use crate::cluster::*;
use float_ord::FloatOrd;

/// Approximation of a sample distribution's [quantile
/// function](https://en.wikipedia.org/wiki/Quantile_function).
///
/// As more samples are merged into the digest, it uses a scaling function to
/// bound the number of [Cluster]s.
#[derive(Debug, Clone)]
pub struct TDigest {
    clusters: Vec<Cluster>,
    count: usize,
    sum: f64,
    min: f64,
    max: f64,
}

impl Default for TDigest {
    #[inline]
    fn default() -> Self {
        Self {
            clusters: Vec::new(),
            count: 0,
            sum: 0.0,
            min: std::f64::INFINITY,
            max: std::f64::NEG_INFINITY,
        }
    }
}

impl TDigest {
    /// Total number of samples.
    #[inline]
    pub fn count(&self) -> usize {
        self.count
    }

    /// Sum of all samples.
    #[inline]
    pub fn sum(&self) -> f64 {
        self.sum
    }

    /// Mean of all samples.
    #[inline]
    pub fn mean(&self) -> Option<f64> {
        (self.count > 0).then_some(self.sum / self.count as f64)
    }

    /// Minimum of all samples.
    #[inline]
    pub fn min(&self) -> f64 {
        self.min
    }

    /// Maximum of all samples.
    #[inline]
    pub fn max(&self) -> f64 {
        self.max
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.clusters.is_empty()
    }

    /// Convenience method that sorts `values` before calling `merge_sorted`.
    pub fn merge_unsorted(self, max_clusters: usize, values: &mut [f64]) -> Self {
        values.sort_by_key(|v| FloatOrd(*v));
        self.merge_sorted(max_clusters, values)
    }

    /// Creates a cluster for each of `sorted_values` and merges them with the
    /// clusters of `self`, using the scaling function to keep the total
    /// resulting cluster count below `max_clusters`.
    pub fn merge_sorted(self, max_clusters: usize, sorted_values: &[f64]) -> Self {
        if sorted_values.is_empty() {
            return self;
        }

        let total_count = self.count + sorted_values.len();
        let mut sum = self.sum;
        let mut min = self.min;
        let mut max = self.max;
        for value in sorted_values.iter().cloned() {
            sum += value;
            min = min.min(value);
            max = max.max(value);
        }

        let mut clusters = self.clusters;
        clusters.extend(
            sorted_values
                .iter()
                .cloned()
                .map(|val| Cluster::new(val, 1)),
        );
        clusters.sort();
        let merged_clusters = Self::merge_clusters(max_clusters, total_count, &clusters);

        Self {
            clusters: merged_clusters,
            count: total_count,
            sum,
            min,
            max,
        }
    }

    /// Merges all of the clusters of `digests`, using the scaling function to
    /// keep the total resulting cluster count below `max_clusters`.
    pub fn merge_digests(max_clusters: usize, digests: Vec<Self>) -> Self {
        let n_clusters = digests.iter().map(|d| d.clusters.len()).sum();
        if n_clusters == 0 {
            return Self::default();
        }

        let mut clusters = Vec::with_capacity(n_clusters);
        let mut total_count = 0;
        let mut sum = 0.0;
        let mut min = std::f64::INFINITY;
        let mut max = std::f64::NEG_INFINITY;
        for mut digest in digests.into_iter() {
            total_count += digest.count;
            sum += digest.sum;
            min = min.min(digest.min);
            max = max.max(digest.max);
            clusters.append(&mut digest.clusters);
        }
        clusters.sort();
        let mut merged = Self {
            clusters,
            count: total_count,
            sum,
            min,
            max,
        };
        merged.clusters = Self::merge_clusters(max_clusters, merged.count, &merged.clusters);
        merged
    }

    /// Merge adjacent clusters while keeping clusters smaller than the upper
    /// bound determined by the scaling function.
    fn merge_clusters(
        max_clusters: usize,
        total_count: usize,
        sorted_clusters: &[Cluster],
    ) -> Vec<Cluster> {
        if total_count == 0 {
            return Vec::new();
        }

        let mut merged_clusters: Vec<Cluster> = Vec::with_capacity(max_clusters);
        let mut current_merge = Cluster::default();
        {
            let mut count_so_far = 0;
            let mut k_limit = 1.0;
            let mut upper_bound = k_to_q(k_limit, max_clusters as f64) * total_count as f64;
            for cluster in sorted_clusters {
                count_so_far += cluster.count;
                if count_so_far as f64 > upper_bound {
                    // Finish the current merge.
                    merged_clusters.push(current_merge);
                    // Start a new merge.
                    current_merge = Cluster::default();
                    k_limit += 1.0;
                    upper_bound = k_to_q(k_limit, max_clusters as f64) * total_count as f64;
                }
                current_merge += cluster.clone();
            }
        }
        merged_clusters.push(current_merge);
        merged_clusters.sort(); // Only necessary to fix float imprecision.
        merged_clusters
    }

    /// Returns an estimate for
    /// [quantile](https://en.wikipedia.org/wiki/Quantile) `q` where `0.0 <= q
    /// <= 1.0`.
    ///
    /// For example:
    ///   - `q=0.0` returns the _minimum_
    ///   - `q=0.5` returns the _median_
    ///   - `q=1.0` returns the _maximum_
    pub fn estimate_quantile(&self, q: f64) -> Option<f64> {
        if self.clusters.is_empty() {
            return None;
        }

        let rank = q * self.count as f64;

        // Find cluster containing desired rank.
        let mut pos;
        let mut t;
        if q > 0.5 {
            if q >= 1.0 {
                return Some(self.max());
            }

            // Scan back to front until finding the cluster containing the
            // desired rank.
            pos = 0;
            t = self.count;
            for (k, cluster) in self.clusters.iter().enumerate().rev() {
                t -= cluster.count;
                if rank >= t as f64 {
                    pos = k;
                    break;
                }
            }
        } else {
            if q <= 0.0 {
                return Some(self.min());
            }

            // Scan front to back until finding the cluster containing the
            // desired rank.
            pos = self.clusters.len() - 1;
            t = 0;
            for (k, cluster) in self.clusters.iter().enumerate() {
                if rank < (t + cluster.count) as f64 {
                    pos = k;
                    break;
                }
                t += cluster.count;
            }
        }

        // Interpolate between clusters.
        let center = self.clusters[pos].mean();
        let mut delta = 0.0;
        let mut min = self.min;
        let mut max = self.max;
        if self.clusters.len() > 1 {
            if pos == 0 {
                max = self.clusters[pos + 1].mean();
                delta = max - center;
            } else if pos == (self.clusters.len() - 1) {
                min = self.clusters[pos - 1].mean();
                delta = center - min;
            } else {
                min = self.clusters[pos - 1].mean();
                max = self.clusters[pos + 1].mean();
                delta = (max - min) / 2.0;
            }
        }
        let value = center + ((rank - t as f64) / self.clusters[pos].count as f64 - 0.5) * delta;
        Some(value.clamp(min, max))
    }
}

/*
 * A good biased scaling function has the following properties:
 *   - The value of the function k(0, delta) = 0, and k(1, delta) = delta.
 *     This is a requirement for any t-digest function.
 *   - The limit of the derivative of the function dk/dq at 0 is inf, and at
 *     1 is inf. This provides bias to improve accuracy at the tails.
 *   - For any q <= 0.5, dk/dq(q) = dk/dq(1-q). This ensures that the accuracy
 *     of upper and lower quantiles are equivalent.
 *
 * The scaling function used here is...
 *   k(q, d) = (IF q >= 0.5, d - d * sqrt(2 - 2q) / 2, d * sqrt(2q) / 2)
 *
 *   k(0, d) = 0
 *   k(1, d) = d
 *
 *   dk/dq = (IF q >= 0.5, d / sqrt(2-2q), d / sqrt(2q))
 *   limit q->1 dk/dq = inf
 *   limit q->0 dk/dq = inf
 *
 *   When plotted, the derivative function is symmetric, centered at q=0.5.
 *
 * Note that FMA has been tested here, but benchmarks have not shown it to be a
 * performance improvement.
 */

// Scale function
//
// fn q_to_k(q: f64, d: f64) -> f64 {
//     if q >= 0.5 {
//         d - d * (0.5 - 0.5 * q).sqrt()
//     } else {
//         d * (0.5 * q).sqrt()
//     }
// }

/// Inverse of the scale function
fn k_to_q(k: f64, d: f64) -> f64 {
    let k_div_d = k / d;
    if k_div_d >= 0.5 {
        let base = 1.0 - k_div_d;
        1.0 - 2.0 * base * base
    } else {
        2.0 * k_div_d * k_div_d
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const ERR: f64 = 0.01;

    #[test]
    fn test_merge_against_uniform_distribution() {
        let t = TDigest::default();
        let values: Vec<_> = (1..=1_000_000).map(f64::from).collect();

        let t = t.merge_sorted(100, &values);

        let ans = t.estimate_quantile(1.0).unwrap();
        let expected = 1_000_000.0;

        let percentage = (expected - ans).abs() / expected;
        assert!(percentage < ERR);

        let ans = t.estimate_quantile(0.99).unwrap();
        let expected = 990_000.0;

        let percentage = (expected - ans).abs() / expected;
        assert!(percentage < ERR);

        let ans = t.estimate_quantile(0.01).unwrap();
        let expected = 10_000.0;

        let percentage = (expected - ans).abs() / expected;
        assert!(percentage < ERR);

        let ans = t.estimate_quantile(0.0).unwrap();
        let expected = 1.0;

        let percentage = (expected - ans).abs() / expected;
        assert!(percentage < ERR);

        let ans = t.estimate_quantile(0.5).unwrap();
        let expected = 500_000.0;

        let percentage = (expected - ans).abs() / expected;
        assert!(percentage < ERR);
    }

    #[test]
    fn test_merge_against_skewed_distribution() {
        let t = TDigest::default();
        let mut values: Vec<_> = (1..=600_000).map(f64::from).collect();
        for _ in 0..400_000 {
            values.push(1_000_000.0);
        }

        let t = t.merge_sorted(100, &values);

        let ans = t.estimate_quantile(0.99).unwrap();
        let expected = 1_000_000.0;
        let percentage = (expected - ans).abs() / expected;
        assert!(percentage < ERR);

        let ans = t.estimate_quantile(0.01).unwrap();
        let expected = 10_000.0;

        let percentage = (expected - ans).abs() / expected;
        assert!(percentage < ERR);

        let ans = t.estimate_quantile(0.5).unwrap();
        let expected = 500_000.0;

        let percentage = (expected - ans).abs() / expected;
        assert!(percentage < ERR);
    }

    #[test]
    fn test_merge_digests() {
        let mut digests: Vec<TDigest> = Vec::new();

        for _ in 1..=100 {
            let t = TDigest::default();
            let values: Vec<_> = (1..=1_000).map(f64::from).collect();
            let t = t.merge_sorted(100, &values);
            digests.push(t)
        }

        let t = TDigest::merge_digests(100, digests);

        let ans = t.estimate_quantile(1.0).unwrap();
        let expected = 1000.0;

        let percentage = (expected - ans).abs() / expected;
        assert!(percentage < ERR);

        let ans = t.estimate_quantile(0.99).unwrap();
        let expected = 990.0;

        let percentage = (expected - ans).abs() / expected;
        assert!(percentage < ERR);

        let ans = t.estimate_quantile(0.01).unwrap();
        let expected = 10.0;

        let percentage = (expected - ans).abs() / expected;
        assert!(percentage < 0.2);

        let ans = t.estimate_quantile(0.0).unwrap();
        let expected = 1.0;

        let percentage = (expected - ans).abs() / expected;
        assert!(percentage < ERR);

        let ans = t.estimate_quantile(0.5).unwrap();
        let expected = 500.0;

        let percentage = (expected - ans).abs() / expected;
        assert!(percentage < ERR);
    }
}
