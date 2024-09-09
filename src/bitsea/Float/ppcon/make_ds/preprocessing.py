def missing_extreme_test(pres_variable_df, min_extreme=10, max_extreme=180):
    """
    remove float where there are no measurements at the top/bottom of the interval
    :param pres_variable_df: df of the pressure of the variable
    :param min_extreme: minimum value where there must be at least one measurement
    :param max_extreme: maximum value where there must be at least one measurement
    :return: flag 1 if the measurement is accepted, 0 of it need to be removed
    """
    if pres_variable_df[0] > min_extreme or pres_variable_df[-1] < max_extreme:
        return 0
    else:
        return 1


def counting_measurement_test(variable_df, tradeoff=50):
    """
    removing float which do not have enough measurements (tradeoff)
    :param variable_df: df of the variable
    :param tradeoff: minimum number of measurements required
    :return: flag 1 if the measurement is accepted, 0 of it need to be removed
    """
    if len(variable_df) >= tradeoff:
        return 1
    else:
        return 0


def noisy_profile_test(discrete_variable, percentage_tradeoff=0.1):
    """
    removing float with too noise
    :param discrete_variable: df of the variable
    :param percentage_tradeoff: tradeoff percentage of noisy points
    :return: flag 1 if the measurement is accepted, 0 of it need to be removed
    """
    return


