import numpy as np
import pandas as pd
import pytest

from eNMRpy.Measurement.eNMR_Methods import _eNMR_Methods

X_AXIS = "U / [V]"
SLOPE = 0.5
INTERCEPT = 3.0
OUTLIER_IDX = [3, 14]


def make_measurement():
    """Bypasses __init__, which needs raw Bruker data, and sets only what lin_huber uses."""
    rng = np.random.default_rng(0)
    u = np.linspace(-100, 100, 21)
    ph0 = SLOPE * u + INTERCEPT + rng.normal(0, 0.2, u.size)
    ph0[OUTLIER_IDX] += [40, -35]

    # shuffled rows, since lin_huber sorts by the x-axis itself
    order = rng.permutation(u.size)
    meas = object.__new__(_eNMR_Methods)
    meas.eNMRraw = pd.DataFrame({X_AXIS: u[order], "ph0": ph0[order]})
    meas._x_axis = X_AXIS
    meas.lin_res_dic = {}
    return meas, u


def test_lin_huber_recovers_line_despite_outliers():
    meas, u = make_measurement()

    meas.lin_huber()
    res = meas.lin_res_dic["ph0"]

    assert np.ravel(res["m"])[0] == pytest.approx(SLOPE, abs=0.01)
    assert np.ravel(res["b"])[0] == pytest.approx(INTERCEPT, abs=0.3)
    assert 0 < res["sig_m"] < 0.01
    assert 0 < res["r_square"] <= 1
    np.testing.assert_allclose(res["x"], u)
    assert res["y_fitted"].shape == u.shape


def test_lin_huber_marks_outliers():
    meas, u = make_measurement()

    meas.lin_huber()

    flagged = meas.eNMRraw.loc[meas.eNMRraw["outlier"], X_AXIS]
    np.testing.assert_allclose(sorted(flagged), u[OUTLIER_IDX])


def test_lin_huber_respects_ulim():
    meas, _ = make_measurement()

    meas.lin_huber(ulim=(-40, 40))
    res = meas.lin_res_dic["ph0"]

    np.testing.assert_allclose(res["x"], np.arange(-40, 41, 10))
    # m is used as mu[0] downstream in lin_results_df, so it must stay indexable
    assert np.ndim(res["m"]) == 1
