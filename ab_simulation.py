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
    plt.show()


def run_simulation():

    # a/b experiment parameters
    baseline_odds = 0.01
    test_odds_delta = 0
    alpha = 0.05
    beta = 0.2
    minimum_effect_size = 0.005

    # f(alpha, beta, baseline_odds, minimum_effect_size) = samples_per_trial
    # http://sphweb.bumc.bu.edu/otlt/MPH-Modules/BS/BS704_Power/BS704_Power_print.html
    samples_per_trial = int(2 * baseline_odds * (1 - baseline_odds) * ((norm.ppf(1-alpha/2) + norm.ppf(1-beta)) / minimum_effect_size)**2)

    # number of experiments to run
    total_trials = 100

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

    false_positives = 0

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

        p = chi2.isf(X_2, degrees_of_freedom)
        significant = p < 0.05

        if significant:
            false_positives += 1

    # number of times that the experiment produced a false positive--there was no difference between the test and the
    # control, but the experiment suggests that there is
    print("False positives: %s / %s = %s" % (false_positives, total_trials, false_positives / total_trials))
    return false_positives


def main():

    simulation = {'data': [], 'bins': 100,
                  'figure': 0,
                  'xLabel': 'False positives',
                  'yLabel': 'Counts',
                  'title': 'A/B Experiments'}

    for runs in range(100):
        false_positives = run_simulation()
        simulation['data'].append(false_positives)

    average = sum(simulation['data']) / len(simulation['data'])
    print("Mean false positives per run: %s" % average)
    histogram(simulation)


if __name__ == "__main__":
    main()
