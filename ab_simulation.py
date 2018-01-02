import matplotlib.pyplot as plt
import numpy as np
import random
from scipy.stats import chi2, norm


def histogram(parameters):

    hist, bins = np.histogram(parameters['data'], parameters['bins'])
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2

    plt.figure(parameters['figure'])
    plt.xlabel(parameters['xLabel'])
    plt.ylabel(parameters['yLabel'])
    plt.suptitle(parameters['title'])
    plt.bar(center, hist, align='center', width=width)

# run_simulation()
# simulate an A/B experiment by modeling a series of coin flips (draw samples from a binomial distribution)
#
# inputs
#   baseline_odds: float between 0 and 1; the control threshold for heads
#   test_odds_delta: float between 0 and 1; baseline_odds + test_odds_delta indicates the test threshold for heads
#   alpha: float between 0 and 1; false positive rate
#   beta: float between 0 and 1; false negative rate
#   minimum_effect_size: float between 0 and 1; min. desired effect size to detect if the control and treatment are different
#   total_trials: integer; total number of trials to run (each trial draws a dynamically calculated number of samples)
#
# outputs
#   significant_differences: the number of trials that resulted in a statistically different number of heads
def run_simulation(baseline_odds=0.01, test_odds_delta=0, alpha=0.05, beta=0.2, minimum_effect_size=0.005, total_trials=100):

    # f(alpha, beta, baseline_odds, minimum_effect_size) = samples_per_trial
    # http://sphweb.bumc.bu.edu/otlt/MPH-Modules/BS/BS704_Power/BS704_Power_print.html
    samples_per_trial = int(2 * baseline_odds * (1 - baseline_odds) * ((norm.ppf(1-alpha/2) + norm.ppf(1-beta)) / minimum_effect_size)**2)

    # chi-squared inputs
    # expected counts, assuming independence
    # -----------------------------------
    # |         | test | control | TOTAL
    # -----------------------------------
    # | heads   | 50    | 50     | 100
    # -----------------------------------
    # | tails   | 50    | 50     | 100
    # -----------------------------------
    # | TOTAL   | 100   | 100    | 200
    #
    # degrees of freedom: (r - 1)(c - 1) = (2 - 1)(2 - 1) = 1
    degrees_of_freedom = 1

    significant_differences = 0

    for _ in range(total_trials):

        control_heads_per_trial = 0
        for _ in range(samples_per_trial):
            if random.random() < baseline_odds:
                control_heads_per_trial += 1

        control_tails_per_trial = samples_per_trial - control_heads_per_trial

        test_heads_per_trial = 0
        for _ in range(samples_per_trial):
            if random.random() < baseline_odds + test_odds_delta:
                test_heads_per_trial += 1

        test_tails_per_trial = samples_per_trial - test_heads_per_trial

        total_number_of_samples = samples_per_trial * 2

        expected_heads_percentage = (control_heads_per_trial + test_heads_per_trial) / (total_number_of_samples)
        expected_tails_percentage = (control_tails_per_trial + test_tails_per_trial) / (total_number_of_samples)

        expected_heads = expected_heads_percentage * samples_per_trial
        expected_tails = expected_tails_percentage * samples_per_trial
        X_2_control_heads = ((control_heads_per_trial - expected_heads) ** 2) / expected_heads
        X_2_control_tails = ((control_tails_per_trial - expected_tails) ** 2) / expected_tails
        X_2_test_heads = ((test_heads_per_trial - expected_heads) ** 2) / expected_heads
        X_2_test_tails = ((test_tails_per_trial - expected_tails) ** 2) / expected_tails

        X_2 = X_2_control_heads + X_2_control_tails + X_2_test_heads + X_2_test_tails

        p = chi2.sf(X_2, degrees_of_freedom)
        significant = p < alpha

        if significant:
            significant_differences += 1

    return significant_differences


def main():

    # number of trials to perform for each simulation; each
    # trial draws SAMPLE_SIZE samples (dynamically calculated based on A/B experiment parameters)
    total_trials = 100

    # number of times to run the simulation; each run contains total_trials number of trials
    # multiple runs allow us to observe the variance in false positives/negatives across runs
    # on average, the number of false positives/negatives should correspond with the selected alpha/beta values
    total_runs = 100

    false_positives_simulation = {'data': [], 'bins': 100,
                  'figure': 0,
                  'xLabel': 'False positives',
                  'yLabel': 'Counts',
                  'title': 'A/B Experiments'}

    # false positives exploration
    for runs in range(total_runs):
        significant = run_simulation(
            baseline_odds=0.01,
            test_odds_delta=0,
            alpha=0.05,
            beta=0.2,
            minimum_effect_size=0.005,
            total_trials=total_trials
        )
        false_positives_simulation['data'].append(significant)

    average = sum(false_positives_simulation['data']) / len(false_positives_simulation['data'])
    print("Mean false positives per run: %s" % average)
    histogram(false_positives_simulation)

    # false negatives exploration
    false_negatives_simulation = {'data': [], 'bins': 100,
                  'figure': 1,
                  'xLabel': 'False negatives',
                  'yLabel': 'Counts',
                  'title': 'A/B Experiments'}

    for runs in range(total_runs):
        significant = run_simulation(
            baseline_odds=0.01,
            test_odds_delta=0.005,
            alpha=0.05,
            beta=0.2,
            minimum_effect_size=0.005,
            total_trials=total_trials
        )
        false_negatives_simulation['data'].append(total_trials - significant)

    average = sum(false_negatives_simulation['data']) / len(false_negatives_simulation['data'])
    print("Mean false negatives per run: %s" % average)
    histogram(false_negatives_simulation)

    plt.show()

if __name__ == "__main__":
    main()
