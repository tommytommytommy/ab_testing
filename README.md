# A/B Testing Simulation
This simulation explores the interactions between conversion rates, minimum effect sizes, maximum probability for false positives/false negatives, and sample sizes with respect to A/B testing and categorical metrics (metrics that can only be answered "yes" or "no"; for example, flipping coins, examining user clicks, and other behaviors that can be modeled as a binomial distribution).

The simulation takes in the following parameters:
1) `alpha`, the maximum probability for false positives
2) `beta`, where the maximum probability for false negatives is 1-beta
3) `p`, the underlying odds for the control group; this is the parameter for the binomial distribution
4) `minimum_effect_size`, the minimum desired change to detect, if a change exists
5) `test_odds_delta`, the difference between the control and the test's binomial distribution parameter

This simulation always flips two coins, one with odds `p` and one with odds `p + test_odds_delta`. Flips are aggregated into runs, and each run compares the number of heads and tails in the control versus the test using a Chi-square test for independence. 

If `test_odds_delta` is 0, then the objective truth is that the control and the test are the same, and any statistically significant differences at the p=0.05 level are false positives.

If `test_odds_delta` is not 0, then the objective truth is that the control and the test are different, and any not statistically significant differences are false negatives.

The simulation builds a histogram of false positives or false negatives based on the analysis mode. To help visualize the variance in the percentage of runs that are false positives/false negatives, coin flips are aggregated in runs, runs are aggregated in trials, and the number of false positives/negatives are aggregated at the trial level. By default, the simulation generates 100 runs per trial (each run containing the dynamically calculated sample size for the test and the control), and 100 total trials.
