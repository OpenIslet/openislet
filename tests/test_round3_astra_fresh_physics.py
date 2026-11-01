"""Independent analytical and manufactured-solution checks for the fresh physics audit."""
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


@pytest.mark.parametrize("step", [.25, .125, .0625, .17])
@pytest.mark.parametrize("closed_endpoints", [False, True])
def test_declared_pulse_has_exact_duration_independent_of_rk4_endpoint_values(step, closed_endpoints):
    start, stop = .83, 2.11
    def rhs(t, y):
        on = start <= t <= stop if closed_endpoints else start < t < stop
        return np.array([float(on)])
    rhs.integration_breakpoints = (start, stop)
    result = l4.rk4_integrate(rhs, (0., 3.), np.zeros(1), [0., 1.2, 3.], dt_max=step)
    expected = np.array([[0., 1.2-start, stop-start]])
    require(np.max(np.abs(result-expected)) < 2e-15,
            "RK4 added or discarded exposure at an endpoint or stepped across a declared switch")


def test_piecewise_forced_decay_matches_analytic_convolution():
    start, stop, rate = .83, 2.11, .7
    def rhs(t, y):
        return -rate*y + float(start <= t <= stop)
    rhs.integration_breakpoints = (start, stop)
    expected = math.exp(-rate*3.) + (1.-math.exp(-rate*(stop-start)))/rate*math.exp(-rate*(3.-stop))
    errors = [abs(l4.rk4_integrate(rhs, (0., 3.), np.ones(1), [3.], dt_max=h)[0, 0]-expected)
              for h in (.1, .05, .025)]
    require(errors[-1] < 2e-9, "piecewise forcing disagrees with its exact exponential convolution")
    require(all(12. < a/b < 20. for a, b in zip(errors, errors[1:])),
            "known forcing changes reduced smooth-segment RK4 to first order")


@pytest.mark.parametrize("bad", [[float("nan")], [True], [-1.], [3.1], ["1.0"]])
def test_rk4_rejects_invalid_declared_breakpoints(bad):
    def rhs(t, y):
        return -y
    rhs.integration_breakpoints = bad
    with pytest.raises(ValueError, match="integration_breakpoints"):
        l4.rk4_integrate(rhs, (0., 3.), np.ones(1), [3.])


@pytest.mark.parametrize("nodes", [51, 101, 201])
@pytest.mark.parametrize("width", [.34, 100./350.])
@pytest.mark.parametrize("scale", [.1754**1.5, .42**1.5, 1.])
def test_sharp_lamina_flux_matches_independent_piecewise_linear_resistance(monkeypatch, nodes, width, scale):
    monkeypatch.setattr(l4, "CENTRAL_LAMINA_HALF_WIDTH_FRACTION", width)
    half = .0175
    z = np.linspace(-half, half, nodes)
    # Integral from -L to z of 1/D, independently constructed as three segments.
    inside_length = np.clip(z+width*half, 0., 2.*width*half)
    pressure = (z+half-inside_length)/l4.D_O2_HYDROGEL_CM2_S
    pressure += inside_length/(scale*l4.D_O2_HYDROGEL_CM2_S)
    flux = l4._slab_oxygen_face_diffusivities(z, scale)*np.diff(pressure)/np.diff(z)
    require(np.max(np.abs(flux-1.)) < 2e-12,
            "face operator does not preserve constant flux through the actual material interface")


@pytest.mark.parametrize("surface", [.00011, .001, .01, .1, 1., 40.])
def test_ewald_diagnostic_solves_its_shell_equation_without_creating_oxygen(surface):
    radius = .006
    result = l4.solve_spherical_cluster_po2_bvp(surface, radius)
    diagnostic = result["ewald_analytic_shell_solver_diagnostics"]
    p_surface = diagnostic["surface_po2_unrounded_mmhg"]
    require(0. <= p_surface <= surface, "a consuming shell created oxygen above the prescribed mean")
    core, uptake, spherical = l4._solve_fixed_surface_spherical_bvp_with_diagnostics(p_surface, radius)
    x = .1262247**(1./3.)
    factor = 1.-1.8*x+x**3-.2*x**6
    resistance = factor/(4.*math.pi*radius*l4.D_O2_HYDROGEL_CM2_S*l4.ALPHA_BUNSEN_MOL_CM3_MMHG)
    residual = p_surface + resistance*uptake - surface
    require(abs(residual) <= 1.01e-10*max(1., surface), "Ewald shell did not close its own supply-drop equation")
    require(diagnostic["converged"] and spherical["converged"], "reported diagnostic did not solve its coupled equations")
    require(core <= p_surface, "consuming spherical tissue created oxygen")


@pytest.mark.parametrize("phi", [.8, .99])
def test_ewald_shell_resistance_vanishes_as_extracellular_shell_disappears(phi):
    result = l4.solve_spherical_cluster_po2_bvp(40., .006, phi_islet_in_slab=phi)
    x = phi**(1./3.)
    exact = 1.-1.8*x+phi-.2*phi**2
    require(abs(result["ewald_shell_form_factor_exact"]-exact) <= 5.1e-8,
            "an artificial positive Ewald floor prevents the no-shell limit")
    if phi == .99:
        require(40.-result["ewald_analytic_shell_solver_diagnostics"]["surface_po2_unrounded_mmhg"] < 1e-4,
                "vanishing shell retained a finite oxygen drop")


def test_slab_minimum_is_reduced_from_solved_states_without_a_100mmhg_cap():
    result = l4._run_single_pouch_pde_regime(
        "FRESH-HIGH-OXYGEN-CONTROL", 0., "none", 68., True, 0.,
        boundary_po2_mmhg=250., demand_multiplier_extra=.01,
        k_perf_max_s_override=0., nz_override=5, _assess_numerical_gate=False)
    reported = [row["min_slab_po2_mmhg"] for row in result["trajectory_0h_to_336h"]]
    require(min(reported) > 200., "control did not isolate a fully oxygenated steady slab")
    require(abs(result["min_slab_po2_unrounded_mmhg"]-min(reported)) < .0051,
            "aggregate minimum is an initializer rather than any solved slab state")


@pytest.mark.parametrize("peak", [1., 25000., 250000.])
def test_nonlinear_peroxide_recovers_independent_manufactured_parabola(peak):
    n, length, catalytic, perfusion, km = 39, .035, 18., .017, 25000.
    spacing = length/(n+1)
    z = np.arange(1, n+1)*spacing
    expected = 4.*peak*z*(length-z)/length**2
    reaction = catalytic*expected/(1.+expected/km)
    source = 8.*l4.D_H2O2_HYDROGEL_CM2_S*peak/length**2 + reaction + perfusion*expected
    actual = l4._solve_peroxide_slab_bvp(source, spacing, catalytic, perfusion, km_uM=km)
    require(np.max(np.abs(actual["concentration_uM"]-expected)) < 1e-8*max(1., peak),
            "nonlinear catalase solution differs from a manufactured field with analytically exact curvature")
    boundary_flux = l4.D_H2O2_HYDROGEL_CM2_S*(expected[0]+expected[-1])/spacing
    removal = spacing*sum(actual["decomposition_uM_s"]+perfusion*actual["concentration_uM"])+boundary_flux
    require(abs(removal-spacing*sum(source)) < 1e-7*max(1., peak),
            "manufactured peroxide field loses mass through clearance or the boundaries")


def test_two_zone_linear_limit_matches_independent_continuum_moments(monkeypatch):
    radius, phi, mean, x, km = .006, .18, 20., 1.3, 1e8
    di, dg = l4.D_O2_ISLET_TISSUE_CM2_S, l4.D_O2_HYDROGEL_CM2_S
    outer = radius*phi**(-1./3.)
    # p_inner(r)=C*sinh(lambda*r)/(lambda*r); normalize C=1,
    # use flux continuity to determine shell supply and then integrate <p>.
    p_surface = math.sinh(x)/x
    supply = 3.*di*radius*(math.cosh(x)-math.sinh(x)/x)/(outer**3-radius**3)
    inner_moment = radius**3*(math.cosh(x)/x**2-math.sinh(x)/x**3)
    shell_moment = p_surface*(outer**3-radius**3)/3. + supply/(3.*dg)*(
        outer**3*((outer**3-radius**3)/(3.*radius)-(outer**2-radius**2)/2.)
        -.5*((outer**5-radius**5)/5.-radius**2*(outer**3-radius**3)/3.))
    exact_core = mean/(3.*(inner_moment+shell_moment)/outer**3)
    monkeypatch.setattr(l4, "KM_O2_MMHG", km)
    monkeypatch.setattr(l4, "VMAX_ISLET_PURE_MOL_CM3_S", l4.ALPHA_BUNSEN_MOL_CM3_MMHG*di*(x/radius)**2*km)
    errors = []
    for nodes in (12, 24, 48):
        actual = l4._solve_two_zone_fv_wigner_seitz_bvp(mean, radius, phi, nodes, "gel_only")
        require(actual["converged"], "first-order two-zone BVP failed to converge")
        errors.append(abs(actual["core_po2_unrounded_mmhg"]-exact_core))
        require(abs(actual["integrated_mass_balance_error_mmhg_s"]) < 1e-8,
                "no-flux unit cell violated supply-consumption balance")
    require(errors[-1] < .002 and all(3.8 < a/b < 4.2 for a, b in zip(errors, errors[1:])),
            "two-zone solution failed second-order convergence to independent continuum moments")


@pytest.mark.parametrize("mean", [.01, 2., 40.])
def test_homogeneous_unit_cell_has_constant_exact_solution(mean):
    actual = l4._solve_two_zone_fv_wigner_seitz_bvp(mean, .006, 1., 8, "uniform")
    require(actual["converged"], "homogeneous unit-cell control did not converge")
    require(abs(actual["core_po2_unrounded_mmhg"]-mean) < 1e-10,
            "homogeneous matched supply generated a spurious spatial oxygen gradient")
    exact = l4.VMAX_ISLET_PURE_MOL_CM3_S/l4.ALPHA_BUNSEN_MOL_CM3_MMHG*mean/(l4.KM_O2_MMHG+mean)
    require(abs(actual["oxygen_supply_volume_mean_mmhg_s"]-exact) < 1e-9,
            "homogeneous supply differs from exact Michaelis-Menten consumption")
