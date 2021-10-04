"""test_production.py

This contains validation scripts  comparing circuit-based Z2 simulation to
Erik's numerics code. To run the pre-arranged test suite just call

    python3 -m pytest test_production.py

Alternatively, you can arrange a custom test suite using the __name__==__main__
logic at the end of this module.

Notes:

    - There's currently a discrepancy between the numerics code and cirq
        outputs. However it can be completely fixed in postprocessing.
    - This will not run for 5x5 grid; the sparse matrix operations become
        to expensive to actually validate against the qsimcirq outputs.
"""
import numpy as np
import pytest

from z2_sim.src.QuantumCircuits.Cirq_Code import production
from z2_sim.src.NumericsPython.Z2DualCorrelation import z2_dual_correlater_numeric

# Visually verify agreement between the two simulation methods.
PRODUCE_PLOTS = False
if PRODUCE_PLOTS:
    import matplotlib.pyplot as plt


# Tolerances for comparison to numerics. These are intentionally well above
# machine precision to allow for the error accumulated due to imperfect
# gate decompositions.
ATOL = 1e-4
RTOL = 1e-2


def postprocess_numerics(truth):
    """
    FIXME: There is a complex conjugate discrepancy between the
    numerics and circuit simulation. Less importantly, there's a different
    grid indexing scheme. Both of these can be corrected in postprocessing.
    """
    truth = truth[:,:,1:].transpose(2, 1, 0)
    truth = truth[:,:,::-1] # reflect the grid
    truth = truth.conj() # ????
    return truth


@pytest.mark.parametrize('n', [3, 4])
@pytest.mark.parametrize('trotter_steps', [10])
@pytest.mark.parametrize('obc', [True, False])
@pytest.mark.parametrize('decompbasis', [None, "ISP"])
def test_compute_obs_with_intermediate_state_vector_all_observables(n, trotter_steps, obc, decompbasis):
    """Test that the main qsim routine is consistent with numerics."""

    jcoup = 10
    dt = 0.05

    out = production.compute_obs_with_intermediate_state_vector(
        n, trotter_steps, jcoup, dt, all_observables=True, obc=obc, decompbasis=decompbasis
    )
    assert out.shape == (trotter_steps, n, n)

    truth = z2_dual_correlater_numeric(n, jcoup, 1 / jcoup, dt, trotter_steps, obc=obc, localgroundstate=False)
    truth = postprocess_numerics(truth)

    np.testing.assert_allclose(truth, out, atol=ATOL, rtol=RTOL)
    if PRODUCE_PLOTS:
        fig_re, axes_re = plt.subplots(n, n, sharex=True, constrained_layout=True)
        fig_im, axes_im = plt.subplots(n, n, sharex=True, constrained_layout=True)
        for i in range(n):
            for j in range(n):
                axes_re[i,j].plot(truth[:,i,j].real, c='r', label='Re[numeric]')
                axes_re[i,j].plot(out[:,i,j].real, c='pink',  label='Re[circuit]')
                axes_im[i,j].plot(truth[:,i,j].imag, c='b', label='Im[numeric]')
                axes_im[i,j].plot(out[:,i,j].imag, c='cyan', label='Im[circuit]')
        axes_re[-1, -1].legend()
        axes_im[-1, -1].legend()
        for axes in [axes_re, axes_im]:
            for ax in axes[-1,:]:
                ax.set_xlabel("trotter step")
        for ax in axes_re[:,0]:
            ax.set_ylabel(r"$Re\langle X_i (X+iY)_a\rangle$")
        for ax in axes_im[:,0]:
            ax.set_ylabel(r"$Im\langle X_i (X+iY)_a\rangle$")
        fig_re.savefig("plot_test_production_n{}_re".format(n))
        fig_im.savefig("plot_test_production_n{}_im".format(n))


@pytest.mark.parametrize('n', [3])
@pytest.mark.parametrize('trotter_steps', [10])
@pytest.mark.parametrize('func', [
    production.compute_obs_with_intermediate_state_vector,
    production.compute_ancillafree_obs_with_intermediate_state_vector
])
@pytest.mark.parametrize('obc', [True])
def test_run_parameter_sweep_all_observables(n, trotter_steps, func, obc):
    """This just verifies the sweep outputs a reasonable object."""

    jcoup_arr = [10, 8, 6]
    dt_arr = [0.05, 0.01]

    res = production.run_noiseless_parameter_sweep(
        n,
        trotter_steps,
        jcoup_arr=jcoup_arr,
        dt_arr=dt_arr,
        all_observables=True,
        func=func,
        obc=obc,
        decompbasis=None,
    )
    assert res.shape == (len(jcoup_arr), len(dt_arr), trotter_steps, n, n)

    for j in range(len(jcoup_arr)):
        for k in range(len(dt_arr)):
            truth = z2_dual_correlater_numeric(n, jcoup_arr[j], 1 / jcoup_arr[j], dt_arr[k], trotter_steps, obc=obc, localgroundstate=False)
            truth = postprocess_numerics(truth)
            np.testing.assert_allclose(truth, res[j,k], atol=ATOL, rtol=RTOL)


@pytest.mark.parametrize('n', [3, 4])
@pytest.mark.parametrize('trotter_steps', [10])
@pytest.mark.parametrize('obc', [True, False])
@pytest.mark.parametrize('decompbasis', [None, "ISP"])
def test_compute_ancillafree_obs_with_intermediate_state_vector_all_observables(n, trotter_steps, obc, decompbasis):
    """Test that the ancilla-free qsim routine is consistent with ancilla."""

    jcoup = 10
    dt = 0.05

    out = production.compute_ancillafree_obs_with_intermediate_state_vector(
        n, trotter_steps, jcoup, dt, all_observables=True, obc=obc, decompbasis=decompbasis
    )
    truth = production.compute_obs_with_intermediate_state_vector(
        n, trotter_steps, jcoup, dt, all_observables=True, obc=obc, decompbasis=decompbasis
    )

    assert out.shape == (trotter_steps, n, n)
    np.testing.assert_allclose(truth, out, atol=ATOL, rtol=RTOL)


@pytest.mark.parametrize('n', [3])
@pytest.mark.parametrize('decompbasis', [None, "SIS"])
def test_compute_noisy_obs(n, decompbasis):
    """Test that the no-checkpoint qsim is compatible with checkpoint."""

    jcoup = 10
    dt = 0.05
    obc = True

    trotter_steps = 3
    trotter_start, trotter_stop = (1, trotter_steps + 1)  # EXCLUSIVE indexing
    out = production.compute_noisy_obs(
        n=n,
        trotter_start=trotter_start,
        trotter_stop=trotter_stop,
        jcoup=jcoup,
        dt=dt,
        all_observables=True,
        qsim_options=dict(t=32, f=4, r=1),
        noise_models=None,
        obc=obc,
        decompbasis=decompbasis
    )
    assert out.shape == (trotter_stop - trotter_start, n, n)

    truth = production.compute_obs_with_intermediate_state_vector(
        n, trotter_steps, jcoup, dt, all_observables=True, obc=obc, decompbasis=decompbasis
    )

    np.testing.assert_allclose(truth, out, atol=ATOL, rtol=RTOL)


@pytest.mark.parametrize('n', [3])
@pytest.mark.parametrize('decompbasis', [None, "SIS"])
def test_compute_noisy_obs_composition(n, decompbasis):
    """Test that the no-checkpoint qsim concatenates correctly"""
    jcoup = 10
    dt = 0.05
    obc = True

    start_stops = [
        (1, 4),
        (4, 6),
        (6, 7)
    ]
    all_res = []
    for (trotter_start, trotter_stop) in start_stops:
        out = production.compute_noisy_obs(
            n=n,

            trotter_start=trotter_start,
            trotter_stop=trotter_stop,
            jcoup=jcoup,
            dt=dt,
            all_observables=True,
            qsim_options=dict(t=32, f=4, r=1),
            noise_models=None,
            obc=obc,
            decompbasis=decompbasis
        )
        all_res.append(out)
    out = np.vstack(all_res)

    truth = production.compute_obs_with_intermediate_state_vector(
        n, 6, jcoup, dt, all_observables=True, obc=obc, decompbasis=decompbasis
    )

    np.testing.assert_allclose(truth, out, atol=ATOL, rtol=RTOL)


if __name__ == "__main__":


    # test_compute_ancillafree_obs_with_intermediate_state_vector_all_observables(3, 10, True, "MS")
    test_compute_obs_with_intermediate_state_vector_all_observables(3, 10, True, "MS")
    # test_compute_noisy_obs_composition(3, None)
    # test_compute_noisy_obs(3, None)
     # test_run_parameter_sweep_all_observables(3, 10)

     # FIXME: numerics cant actually support this computation, or at least it
     # hangs for a looong time.
     # test_compute_obs_with_intermediate_state_vector_all_observables(5, 20)
