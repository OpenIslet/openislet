"""Round-3 numerical oracles and adversarial controls; all fixtures are synthetic."""
from __future__ import annotations

import math
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from scripts import run_layer4_tissue_and_metabolic_twin as l4


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def test_rk4_reused_rhs_buffer_retains_independent_stage_values():
    work = np.empty(2)
    def reused_rhs(t, y):
        work[:] = (-y[0], -2.0*y[1])
        return work
    times = [0., .5, 1.]
    reused = l4.rk4_integrate(reused_rhs, (0., 1.), np.array([1., 2.]), times, dt_max=.05)
    fresh = l4.rk4_integrate(lambda t, y: np.array([-y[0], -2.*y[1]]), (0., 1.), np.array([1., 2.]), times, dt_max=.05)
    expected = np.array([[1., math.exp(-.5), math.exp(-1.)], [2., 2.*math.exp(-1.), 2.*math.exp(-2.)]])
    require(np.array_equal(reused, fresh), "RHS memory reuse changed RK4 mathematics")
    require(np.max(np.abs(reused-expected)) < 1e-6, "RK4 disagrees with the independent exponential solution")


def test_rk4_has_fourth_order_step_convergence_for_smooth_decay():
    errors = [abs(l4.rk4_integrate(lambda t, y: -y, (0., 1.), np.array([1.]), [1.], dt_max=h)[0, 0]-math.exp(-1.))
              for h in (.2, .1, .05)]
    require(all(14. < a/b < 18. for a, b in zip(errors, errors[1:])), "smooth RK4 refinement lost fourth-order convergence")


@pytest.mark.parametrize("law", ["first_order", "shrinking_core", "ginstling_brounshtein"])
@pytest.mark.parametrize("half_life", [.01, 10., 68.])
def test_inventory_exact_half_life_nonnegative_and_integrated_source_closure(law, half_life):
    a = 1.-math.pow(.5, 1./3.)
    f_gb_half = 1. - (2./3.)*.5 - (.5**(2./3.))
    depletion = (half_life / (3. * f_gb_half)) if law == "ginstling_brounshtein" else (half_life/a)
    times = np.array([0., half_life, 2.*half_life, depletion, 2.*depletion])
    actual = l4._cao2_inventory_fractions(times, half_life, law)
    require(actual[0] == 1. and abs(actual[1]-.5) < 2e-15, "declared half-life or initial stock changed")
    require(np.all(actual >= 0.) and np.all(np.diff(actual) <= 0.), "inventory left its monotone nonnegative domain")
    # Independent Gauss-Legendre integration of -dn/dt: the SC source is a
    # quadratic until finite-time depletion, first-order is exponential, and
    # Ginstling-Brounshtein uses u = sqrt(t) substitution for the t^(-1/2) start.
    nodes, weights = np.polynomial.legendre.leggauss(24)
    for stop, remaining in zip(times[1:], actual[1:]):
        upper = min(stop, depletion) if law in {"shrinking_core", "ginstling_brounshtein"} else stop
        if law == "ginstling_brounshtein":
            u_max = math.sqrt(upper)
            u_grid = .5 * u_max * (nodes + 1.)
            t_grid = u_grid ** 2
            rem_grid = l4._cao2_inventory_fractions(t_grid, half_life, law)
            k_gb = f_gb_half / half_life
            rate_grid = np.where(rem_grid > 1e-15, 1.5 * k_gb / ((np.maximum(1e-15, rem_grid) ** (-1./3.)) - 1.), 0.)
            released = float(np.dot(weights, .5 * u_max * 2. * u_grid * rate_grid))
            tol = 1e-5 if stop >= depletion else 1e-12
            require(abs(released+remaining-1.) < tol, "integrated release does not close the initial-stock ledger")
        else:
            grid = .5*upper*(nodes+1.)
            rate = (3.*a/half_life*(1.-a*grid/half_life)**2 if law == "shrinking_core"
                    else math.log(2.)/half_life*np.exp(-math.log(2.)*grid/half_life))
            released = .5*upper*float(np.dot(weights, rate))
            require(abs(released+remaining-1.) < 2e-14, "integrated release does not close the initial-stock ledger")
    if law in {"shrinking_core", "ginstling_brounshtein"}:
        require(np.max(actual[-2:]) < 1e-43, "stock survives finite-time depletion")
    if law == "ginstling_brounshtein" and half_life == 68.:
        # AUDIT-REF: verify exact_piecewise_conserving vs early_shell_clamped regularization integral at 336 h
        f_reg_336 = l4._ginstling_brounshtein_piecewise_regularized_hydrolyzed_fraction(336., 68., 0.05)
        require(0.907 < f_reg_336 < 0.908, "piecewise-conserving 336-h GB hydrolyzed fraction drifted")


@pytest.mark.parametrize("times,half_life,law", [([True], 68., "first_order"), ([float("nan")], 68., "first_order"),
    ([-1.], 68., "first_order"), ([0.], True, "first_order"), ([0.], 0., "first_order"), ([0.], 68., "unknown")])
def test_inventory_rejects_invalid_inputs(times, half_life, law):
    with pytest.raises(ValueError):
        l4._cao2_inventory_fractions(np.asarray(times), half_life, law)


@pytest.mark.parametrize("n", [3, 11, 31])
def test_thomas_matches_independent_dense_solve(n):
    rng = np.random.default_rng(729+n)
    lower, upper = rng.uniform(-.7, .7, (2, n-1))
    diagonal = 2.+rng.random(n)
    rhs = rng.normal(size=n)
    matrix = np.diag(diagonal)+np.diag(lower, -1)+np.diag(upper, 1)
    actual = l4._tridiag_thomas_solve(lower, diagonal, upper, rhs)
    require(np.max(np.abs(actual-np.linalg.solve(matrix, rhs))) < 2e-15, "Thomas elimination changed the linear system")


@pytest.mark.parametrize("reaction_derivative", [0., .3])
def test_mean_constrained_newton_matches_full_bordered_system_including_singular_block(reaction_derivative):
    n = 11
    volumes = np.linspace(.2, 2., n)
    weights = volumes/np.sum(volumes)
    conductance = np.linspace(.5, 1.5, n-1)
    lower, upper = conductance/volumes[1:], conductance/volumes[:-1]
    diagonal = -np.r_[upper, 0.]-np.r_[0., lower]
    derivative = np.r_[np.full(5, reaction_derivative), np.zeros(n-5)]
    diagonal -= derivative
    supply = np.r_[np.zeros(5), np.ones(n-5)]
    residual, mean_error = np.linspace(-1., 2., n), .17
    jacobian = np.diag(diagonal)+np.diag(lower, -1)+np.diag(upper, 1)
    bordered = np.block([[jacobian, supply[:, None]], [weights[None, :], np.zeros((1, 1))]])
    expected = np.linalg.solve(bordered, np.r_[-residual, -mean_error])
    dp, ds = l4._mean_constrained_newton_step(lower, diagonal, upper, derivative, weights, supply, residual, mean_error)
    require(np.max(np.abs(np.r_[dp, ds]-expected)) < 3e-13, "reduced constrained solve differs from independent full bordering")
    require(np.max(np.abs(bordered@np.r_[dp, ds]+np.r_[residual, mean_error])) < 1e-12, "bordered Newton residual failed")


@pytest.mark.parametrize("law", ["first_order", "shrinking_core", "ginstling_brounshtein"])
def test_pouch_inventory_remains_defined_for_short_half_life(law):
    result = l4._run_single_pouch_pde_regime("ROUND3-SHORT-LIVED", .45, "central", .01, True, .97,
                                            release_kinetics_law=law, nz_override=5, _assess_numerical_gate=False)
    stock = [row["cao2_molar_inventory_remaining_mmol"] for row in result["trajectory_0h_to_336h"]]
    require(stock[0] > 0. and stock[1] == 0. and all(value >= 0. for value in stock), "production pouch did not use exact nonnegative stock")


def test_spherical_solver_and_uptake_refine_toward_first_order_analytic_oracle(monkeypatch):
    radius, surface, x, km = .006, 20., 2.5, 1e12
    diffusion, alpha = l4.D_O2_ISLET_TISSUE_CM2_S, l4.ALPHA_BUNSEN_MOL_CM3_MMHG
    monkeypatch.setattr(l4, "KM_O2_MMHG", km)
    monkeypatch.setattr(l4, "VMAX_ISLET_PURE_MOL_CM3_S", alpha*diffusion*(x/radius)**2*km)
    core_exact = surface*x/math.sinh(x)
    uptake_exact = 4.*math.pi*diffusion*alpha*surface*radius*(x/math.tanh(x)-1.)
    errors = []
    for nodes in (31, 61, 121):
        core, uptake, diagnostics = l4._solve_fixed_surface_spherical_bvp_with_diagnostics(surface, radius, nodes)
        require(diagnostics["converged"] is True, "linear-limit spherical BVP did not converge")
        errors.append(abs(core-core_exact))
        require(abs(uptake/uptake_exact-1.) < .003, "radial uptake disagrees with the independent boundary-flux integral")
    require(all(3.8 < a/b < 4.2 for a, b in zip(errors, errors[1:])), "spherical field does not refine at second order")


def saturated_two_zone_oracle(mean, radius, phi, closure):
    """Continuum integration of constant reaction in each spherical region."""
    outer = radius/phi**(1./3.)
    q = l4.VMAX_ISLET_PURE_MOL_CM3_S/l4.ALPHA_BUNSEN_MOL_CM3_MMHG
    supply = q*phi/(1. if closure == "uniform" else 1.-phi)
    inner_reaction = q-supply if closure == "uniform" else q
    a = inner_reaction*radius**2/(6.*l4.D_O2_ISLET_TISSUE_CM2_S)
    b = (q if closure == "uniform" else q+supply)*radius**3/(3.*l4.D_O2_HYDROGEL_CM2_S)
    moment_inner = inner_reaction*radius**5/(30.*l4.D_O2_ISLET_TISSUE_CM2_S)
    moment_shell = (a*(outer**3-radius**3)/3.
        + b*((outer**3-radius**3)/(3.*radius)-(outer**2-radius**2)/2.)
        - supply/(6.*l4.D_O2_HYDROGEL_CM2_S)*((outer**5-radius**5)/5.-radius**2*(outer**3-radius**3)/3.))
    core = mean-3.*(moment_inner+moment_shell)/outer**3
    return core, core+a


@pytest.mark.parametrize("closure", ["uniform", "gel_only"])
def test_two_zone_refines_to_piecewise_analytic_solution_and_conserves_flux(monkeypatch, closure):
    monkeypatch.setattr(l4, "KM_O2_MMHG", 1e-12)
    core_exact, surface_exact = saturated_two_zone_oracle(40., .006, .1262247, closure)
    errors = []
    for nodes in (12, 24, 48):
        actual = l4._solve_two_zone_fv_wigner_seitz_bvp(40., .006, n_per_zone=nodes, supply_closure=closure)
        require(actual["converged"] is True, "saturated two-zone solve failed")
        errors.append(abs(actual["core_po2_unrounded_mmhg"]-core_exact))
        require(abs(actual["surface_po2_mmhg"]-surface_exact) < .05, "interface pressure disagrees with the continuum oracle")
        require(abs(actual["mean_constraint_error_mmhg"]) < 1e-11, "unit-cell mean oxygen constraint lost")
        require(abs(actual["integrated_mass_balance_error_mmhg_s"]) < 1e-9, "closed-cell source and consumption do not balance")
        require(abs(actual["interface_flux_balance_error_mmhg_cm3_s"]) < 1e-15, "interface flux fails the inner-volume reaction balance")
    require(all(3.8 < a/b < 4.2 for a, b in zip(errors, errors[1:])), "two-zone field does not refine at second order")


@pytest.mark.parametrize("km", [None, 25000.])
def test_peroxide_nodewise_mass_and_oxygen_stoichiometry(km):
    source = np.linspace(.1, 1e5, 21)
    spacing, catalytic, perfusion = .035/22., 18., .017
    result = l4._solve_peroxide_slab_bvp(source, spacing, catalytic, perfusion, km_uM=km)
    h, decomp = result["concentration_uM"], result["decomposition_uM_s"]
    boundary_outflow = l4.D_H2O2_HYDROGEL_CM2_S*(h[0]+h[-1])/spacing
    input_rate = float(np.sum(source)*spacing)
    removed_rate = float(np.sum(decomp+perfusion*h)*spacing)+boundary_outflow
    require(abs(input_rate-removed_rate) < 1e-7, "peroxide source is not equal to reaction, washout and host-boundary flux")
    expected_decomp = catalytic*h if km is None else catalytic*h/(1.+h/km)
    require(np.max(np.abs(decomp-expected_decomp)) < 1e-9, "oxygen source used a different catalytic reaction")
    require(np.all(.5*decomp >= 0.), "stoichiometric O2 source became negative")


def test_stalled_spherical_newton_cannot_claim_convergence(monkeypatch):
    monkeypatch.setattr(l4, "_tridiag_thomas_solve", lambda lo, diag, up, rhs: np.zeros_like(rhs))
    _, _, result = l4._solve_fixed_surface_spherical_bvp_with_diagnostics(25., .006)
    require(result["converged"] is False and result["residual_inf_mmhg_s"] > 1., "zero Newton step masquerades as oxygen convergence")


def test_stalled_two_zone_newton_cannot_claim_convergence(monkeypatch):
    monkeypatch.setattr(l4, "_mean_constrained_newton_step", lambda lo, diag, up, derivative, weights, supply, residual, mean: (np.zeros_like(residual), 0.))
    result = l4._solve_two_zone_fv_wigner_seitz_bvp(20., .006, n_per_zone=8)
    require(result["converged"] is False and result["residual_inf_mmhg_s"] > 1., "stalled constrained solve passed on its mean constraint alone")


def synthetic_metabolic_integrator(first_glucose_mutation=None, capture=None, low_lead=4., low_oat=3.8):
    calls = [0]
    def integrate(rhs, span, y0, times, dt_max=.25):
        closure = dict(zip(rhs.__code__.co_freevars, (cell.cell_contents for cell in rhs.__closure__)))
        viable = closure["viable_endocrine_fraction"]
        if capture is not None:
            capture.append((viable, float(y0[0])))
        states = np.tile(np.asarray(y0)[:, None], (1, len(times)))
        glucose = 5.
        if closure["ex_start_min"] == 448.:
            glucose = low_lead if closure["enable_low_g_sst_clamp"] else low_oat
        states[0] = glucose
        if calls[0] in {0, 3} and first_glucose_mutation is not None:
            first_glucose_mutation(states[0])
        calls[0] += 1
        return states
    return integrate


@pytest.mark.parametrize("failure", ["severe", "nadir", "time_in_range"])
def test_metabolic_gate_checks_all_solved_points_not_minute_downsample(monkeypatch, failure):
    def mutate(glucose):
        if failure == "time_in_range":
            indices = np.flatnonzero(np.arange(len(glucose)) % 4 != 0)
            glucose[indices[:1000]] = 3.8
        else:
            glucose[1] = 2.99 if failure == "severe" else 3.49
    monkeypatch.setattr(l4, "rk4_integrate", synthetic_metabolic_integrator(mutate))
    row = l4.solve_24h_paracrine_and_whole_body_metabolic_ode()["graft_architectures_compared"][0]
    require(row["min_glucose_mM"] == 5. and row["time_in_tight_range_3p9_to_10p0_mM_pct"] == 100., "fixture did not isolate the downsampling defect")
    require(row["passes_clinical_glycemic_safety_targets"] is False, "quarter-minute safety failure disappeared during minute downsampling")


def coupling_row(**updates):
    row = dict(regime_id="ROUND3-SYNTHETIC", min_120um_cluster_core_po2_mmhg=4., min_slab_po2_across_14_days_mmhg=20.,
        max_interstitial_h2o2_uM=1., max_interstitial_ph=7.4, coupled_dynamic_14d_viable_endocrine_fraction_v80_066=1.,
        passes_all_oxygen_peroxide_and_ph_gates=True, passes_strict_gel_only_two_zone_fv_2p0mmHg_core_gate=True)
    row.update(updates)
    return row


def test_pde_viability_transfer_does_not_round_across_initial_state_branch(monkeypatch):
    seen = []
    monkeypatch.setattr(l4, "rk4_integrate", synthetic_metabolic_integrator(capture=seen))
    source = coupling_row(coupled_dynamic_14d_viable_endocrine_fraction_v80_066=.6,
                          coupled_dynamic_14d_viable_endocrine_fraction_unrounded=.59999)
    row = l4.solve_24h_paracrine_and_whole_body_metabolic_ode([source])["matched_pde_ode_qualification_scenarios"][0]
    actual = [(v, g) for v, g in seen if .59 < v < .61]
    require(actual == [(.59999, 17.2-12.1*.59999)], "display rounding crossed the existing v=0.6 initialization branch")
    require(row["coupled_viable_endocrine_fraction"] == .6 and row["coupled_viable_endocrine_fraction_unrounded"] == .59999,
            "displayed and transferred viability are not distinguished")


def test_pde_chemical_transfer_uses_raw_threshold_inputs(monkeypatch):
    monkeypatch.setattr(l4, "rk4_integrate", synthetic_metabolic_integrator())
    source = coupling_row(max_interstitial_h2o2_uM=5., max_interstitial_h2o2_unrounded_uM=5.00049,
                          max_interstitial_ph=7.6, max_interstitial_ph_unrounded=7.60049)
    row = l4.solve_24h_paracrine_and_whole_body_metabolic_ode([source])["matched_pde_ode_qualification_scenarios"][0]
    expected = math.exp(-.00049/12.)*math.exp(-.00049/.35)
    require(abs(row["coupled_viable_endocrine_fraction_unrounded"]-expected) < 1e-15, "rounding discarded chemical deficits before state transfer")


@pytest.mark.parametrize("field,value", [
    ("coupled_dynamic_14d_viable_endocrine_fraction_unrounded", True),
    ("coupled_dynamic_14d_viable_endocrine_fraction_unrounded", 1.01),
    ("coupled_dynamic_14d_viable_endocrine_fraction_unrounded", float("nan")),
    ("max_interstitial_h2o2_unrounded_uM", -1.),
    ("max_interstitial_h2o2_unrounded_uM", float("inf")),
    ("max_interstitial_ph_unrounded", 14.1),
])
def test_raw_pde_state_is_validated_before_ode_execution(monkeypatch, field, value):
    def forbidden(*args, **kwargs):
        raise AssertionError("invalid PDE state reached integration")
    monkeypatch.setattr(l4, "rk4_integrate", forbidden)
    with pytest.raises(ValueError):
        l4.solve_24h_paracrine_and_whole_body_metabolic_ode([coupling_row(**{field: value})])


def test_low_glucose_runtime_invariant_uses_unrounded_minimum(monkeypatch):
    monkeypatch.setattr(l4, "rk4_integrate", synthetic_metabolic_integrator(low_lead=3.899999))
    with pytest.raises(RuntimeError, match="Low-glucose challenge invariant"):
        l4.solve_24h_paracrine_and_whole_body_metabolic_ode()


@pytest.mark.parametrize("comparison,call_index,message", [("OAT", 3, "OAT_1 invariant"), ("zero_SST", 13, "Zero-SST collapse")])
def test_exact_comparison_rejects_unrounded_state_difference_hidden_by_all_reported_metrics(monkeypatch, comparison, call_index, message):
    original = synthetic_metabolic_integrator()
    calls = [0]
    def altered(*args, **kwargs):
        states = original(*args, **kwargs)
        if calls[0] == call_index:
            # UCN3 at t=0.25 is absent from every reported glucose/shutoff metric.
            states[3, 1] += 1e-12
        calls[0] += 1
        return states
    monkeypatch.setattr(l4, "rk4_integrate", altered)
    with pytest.raises(RuntimeError, match=message):
        l4.solve_24h_paracrine_and_whole_body_metabolic_ode()


def test_fixed_step_metabolic_output_does_not_claim_error_certification(monkeypatch):
    monkeypatch.setattr(l4, "rk4_integrate", synthetic_metabolic_integrator())
    diagnostics = l4.solve_24h_paracrine_and_whole_body_metabolic_ode()["numerical_time_integration_diagnostics"]
    require(diagnostics["status"] == "FIXED_STEP_NOT_ERROR_CERTIFIED" and diagnostics["continuous_time_safety_certified"] is False,
            "finite reporting grid masquerades as continuous-time or error certification")


@pytest.mark.parametrize("step", [.25, .125, .0625])
def test_rk4_inplace_rhs_cannot_mutate_the_base_or_any_saved_stage(step):
    def inplace(t, y):
        np.negative(y, out=y)
        return y
    actual = l4.rk4_integrate(inplace, (0., 1.), np.array([1., 2.]), [0., .5, 1.], dt_max=step)
    fresh = l4.rk4_integrate(lambda t, y: -y, (0., 1.), np.array([1., 2.]), [0., .5, 1.], dt_max=step)
    require(np.array_equal(actual, fresh), "an in-place callback changed the mathematical RK4 state")
    require(abs(actual[0, -1]-math.exp(-1.)) < 2e-5, "in-place decay disagrees with the exponential oracle")


@pytest.mark.parametrize("value", [-1., float("nan"), float("inf"), True, np.bool_(False)])
def test_top_level_raw_core_is_domain_validated_before_coupling(monkeypatch, value):
    def forbidden(*args, **kwargs):
        raise AssertionError("invalid raw legacy core reached integration")
    monkeypatch.setattr(l4, "rk4_integrate", forbidden)
    with pytest.raises(ValueError):
        l4.solve_24h_paracrine_and_whole_body_metabolic_ode([coupling_row(min_120um_cluster_core_po2_unrounded_mmhg=value)])


def test_legacy_core_fallback_prefers_its_exported_raw_companion(monkeypatch):
    seen = []
    monkeypatch.setattr(l4, "rk4_integrate", synthetic_metabolic_integrator(capture=seen))
    source = coupling_row(min_120um_cluster_core_po2_mmhg=1.2, min_120um_cluster_core_po2_unrounded_mmhg=1.19999)
    result = l4.solve_24h_paracrine_and_whole_body_metabolic_ode([source])["matched_pde_ode_qualification_scenarios"][0]
    require([(v, g) for v, g in seen if .59 < v < .61] == [(.599995, 17.2-12.1*.599995)],
            "raw legacy core was replaced by its rounded displayed value")
    require(result["core_coupling_source"] == "legacy_dirichlet_upper_bound_not_qualified", "raw transport promoted legacy scientific qualification")


def test_conservative_core_retains_precedence_over_valid_legacy_raw_core(monkeypatch):
    monkeypatch.setattr(l4, "rk4_integrate", synthetic_metabolic_integrator())
    source = coupling_row(min_120um_cluster_core_po2_mmhg=1.2, min_120um_cluster_core_po2_unrounded_mmhg=1.19999,
        wigner_seitz_120um_cluster_pericellular_shell_audit_v80_001={"two_zone_fv_gel_only_supply_120um_core_po2_unrounded_mmhg": 2.})
    row = l4.solve_24h_paracrine_and_whole_body_metabolic_ode([source])["matched_pde_ode_qualification_scenarios"][0]
    require(row["core_coupling_source"] == "gel_only_unrounded" and row["coupled_viable_endocrine_fraction_unrounded"] == 1.,
            "legacy core overrode the conservative gel-only source")


@pytest.mark.parametrize("oat_nadir", [3.899999, 3.894])
def test_low_glucose_raw_success_is_not_vetoed_by_display_rounding(monkeypatch, oat_nadir):
    monkeypatch.setattr(l4, "rk4_integrate", synthetic_metabolic_integrator(low_lead=4., low_oat=oat_nadir))
    row = l4.solve_24h_paracrine_and_whole_body_metabolic_ode()["low_glucose_challenge_below_4p2mM_oat1_verification_ultracode_v80_15"]
    require(row["low_glucose_challenge_oat1_without_lt_4p2mM_clamp_min_glucose_0p25min_mM"] == oat_nadir,
            "valid raw threshold result was replaced by display rounding")


@pytest.mark.parametrize("contradiction", ["legacy_core", "low_core", "low_surface", "low_slab", "high_peroxide", "high_pH"])
def test_forged_pde_booleans_cannot_override_missing_or_failing_raw_complete_product_evidence(monkeypatch, contradiction):
    monkeypatch.setattr(l4, "rk4_integrate", synthetic_metabolic_integrator())
    source = coupling_row(
        passes_cluster_surface_15mmHg_full_first_phase_gsis_gate_gel_only=True,
        wigner_seitz_120um_cluster_pericellular_shell_audit_v80_001={
            "two_zone_fv_gel_only_supply_120um_core_po2_unrounded_mmhg": 2.,
            "two_zone_fv_gel_only_supply_surface_po2_unrounded_mmhg": 15.,
        },
    )
    if contradiction == "legacy_core":
        del source["wigner_seitz_120um_cluster_pericellular_shell_audit_v80_001"]
    elif contradiction == "low_core":
        source["wigner_seitz_120um_cluster_pericellular_shell_audit_v80_001"]["two_zone_fv_gel_only_supply_120um_core_po2_unrounded_mmhg"] = 1.99999
    elif contradiction == "low_surface":
        source["wigner_seitz_120um_cluster_pericellular_shell_audit_v80_001"]["two_zone_fv_gel_only_supply_surface_po2_unrounded_mmhg"] = 14.99999
    elif contradiction == "low_slab":
        source.update(min_slab_po2_across_14_days_mmhg=15., min_slab_po2_unrounded_mmhg=14.99999)
    elif contradiction == "high_peroxide":
        source.update(max_interstitial_h2o2_uM=5., max_interstitial_h2o2_unrounded_uM=5.00049)
    else:
        source.update(max_interstitial_ph=7.6, max_interstitial_ph_unrounded=7.60049)
    row = l4.solve_24h_paracrine_and_whole_body_metabolic_ode([source])["matched_pde_ode_qualification_scenarios"][0]
    require(row["pde_passes_all_oxygen_peroxide_and_ph_gates"] and row["pde_passes_strict_gel_only_qualification"]
            and row["coupled_metabolic_glycemic_gate_passed"], "fixture did not isolate forged imported flags from raw evidence")
    require(row["coupled_complete_product_pass"] is False and row["qualified_exact_product"] is False,
            "conditional complete-product verdict contradicted its represented numerical evidence")


def test_complete_product_numeric_thresholds_are_closed_without_claiming_exact_product_qualification(monkeypatch):
    monkeypatch.setattr(l4, "rk4_integrate", synthetic_metabolic_integrator())
    source = coupling_row(
        min_slab_po2_across_14_days_mmhg=15.,
        min_slab_po2_unrounded_mmhg=15.,
        max_interstitial_h2o2_uM=5.,
        max_interstitial_h2o2_unrounded_uM=5.,
        max_interstitial_ph=7.6,
        max_interstitial_ph_unrounded=7.6,
        passes_cluster_surface_15mmHg_full_first_phase_gsis_gate_gel_only=True,
        wigner_seitz_120um_cluster_pericellular_shell_audit_v80_001={
            "two_zone_fv_gel_only_supply_120um_core_po2_unrounded_mmhg": 2.,
            "two_zone_fv_gel_only_supply_surface_po2_unrounded_mmhg": 15.,
        },
    )
    row = l4.solve_24h_paracrine_and_whole_body_metabolic_ode([source])["matched_pde_ode_qualification_scenarios"][0]
    require(row["coupled_complete_product_pass"] is True and row["qualified_exact_product"] is False,
            "closed numerical thresholds were confused with empirical exact-product qualification")
