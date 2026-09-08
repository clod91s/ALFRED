"""
═══════════════════════════════════════════════════════════════════════════════
 Positron induced dose and cell survival in liver tissue
═══════════════════════════════════════════════════════════════════════════════

 1.  THE GENERAL G FACTOR — INTEGRAL DEFINITION  [Lea & Catcheside 1942; Dale
     1985; Brenner & Sachs 1998]

     The standard LQ survival model writes:
         SF = exp[ −α·D  −  β · h(Ḋ_peak) · G · D² ]

     where G is the generalized Lea-Catcheside factor:

         G = (2/D²) ∫₀ᵀ Ḋ(t₂) [∫₀^{t₂} Ḋ(t₁) e^{−μ(t₂−t₁)} dt₁] dt₂   (1)

     μ = ln 2 / T_half  [s⁻¹]  is the first-order DSB repair rate.
     This double integral sums, for every pair of time points (t₁ < t₂),
     the product of their dose contributions weighted by the probability
     e^{−μ(t₂−t₁)} that the first DSB has NOT yet been repaired by the time
     the second DSB arrives.

     For CONSTANT dose rate Ḋ over a total time T = D/Ḋ, the integral (1)
     evaluates to the closed form [Dale 1985]:

         G_analytic(T) = 2[μT − 1 + e^{−μT}] / (μT)²              (2)



 ─────────────────────────────────────────────────────────────────────────────
 2.  DISCRETE TIME-STEP FORMULATION

     For an arbitrary time-varying Ḋ(t) discretised into N steps of width Δt:

         ΔD_n  = Ḋ_n · Δt                           [dose deposited in step n]
         w_n   = ∫₀^{t_n} Ḋ(s) e^{−μ(t_n−s)} ds     [sub-lethal-damage memory]

     The memory variable satisfies the exact recursion:

         w_{n+1} = e^{−μΔt} · (w_n + ΔD_n)                         (3)

     (This is the exact solution of dw/dt = Ḋ(t) − μw for a piecewise-
     constant dose rate; it introduces NO approximation beyond discretisation.)

     The quadratic damage accumulator is updated at each step as:

         E₂  +=  2 · ΔD_n · w_n            [interaction: new × past damage]

     At the end of irradiation the G factor is recovered:

         G_discrete = E₂ / D_total²

     and the survival fraction is:

         SF_discrete = exp(−α·D − β·h·E₂)

     This reproduces G_analytic to machine precision for constant Ḋ

 ─────────────────────────────────────────────────────────────────────────────
 3.  HP "Δt = (period)/3" — PHYSICAL AND COMPUTATIONAL CONSIDERATIONS

     Time step has been imposed equal to 1/3 of the shortest beam period.

     SOURCE       Period           Δt = period/3    Steps for 1 Gy @ Ḋ_avg
     ─────────────────────────────────────────────────────────────────────
     Chopped      T_chop = 20 ns   6.67 ns           ~2 × 10¹³   → no time variation required
     Bunched      T_rep  = 120 s   40 s              ~10³–10⁵    → time variation required

     The key Nyquist-Shannon requirement is that Δt should resolve the
     physically relevant timescale, i.e. the one that affects the G factor.
     Two scales matter:

       (a)  The repair half-time T_half = 1800 s  →  μ⁻¹ = 2598 s.
            Any Δt ≪ T_half is sufficient to track repair kinetics.

       (b)  The pulse period (T_chop or T_rep):
            Only relevant if μ·T_period is NOT negligible.

     For the CHOPPED source:
         μ·T_chop = (ln2/1800) × 20 × 10⁻⁹ = 7.7 × 10⁻¹²
         → Negligible repair within one chopper period by 12 orders of
           magnitude.  Brenner & Hall (1991) Int J Radiat Oncol Biol Phys prove analytically that pulsed ≡ continuous when
           τ_ON / T_half → 0.  The discrete stepper with Δt = T_chop/3
           is computationally infeasible AND physically unnecessary.

     For the BUNCHED source:
         μ·T_rep = (ln2/1800) × 120 = 0.0462
         → ~4.5 % of sub-lethal damage repairs between consecutive bunches.
           This is small but non-zero; it accumulates coherently over many
           repetitions and leads to a quantifiable correction
         → v7 uses Δt = T_rep/3 = 33 s for the bunched source.
           Each bunch (τ_b = 10 ns) is modelled as instantaneous (τ_b/Δt
           = 2.5 × 10⁻¹⁰; no physics is lost).

 ─────────────────────────────────────────────────────────────────────────────
 4.  ANALYTIC FORMULA FOR n INSTANTANEOUS BUNCHES

     For n equal instantaneous doses ΔD_b = Ḋ_b_peak · τ_b delivered every
     T_rep with repair factor r = e^{−μ·T_rep} between each pair:

         E₂(n) = 2·ΔD_b² · r/(1−r) · [n − (1 − r^n)/(1−r)]         (4)

         G_bunch(n) = E₂(n) / (n·ΔD_b)²

     Asymptotic behaviours (total dose D = n·ΔD_b, time T = n·T_rep):

       n→1:  G_bunch → 0   (single bunch; no past damage to interact with)
       n→∞:  G_bunch → 2r/[(1−r)·n] = 2r·T_rep/[(1−r)·T]

     Compare the analytic (continuous) formula for the same T:
       G_analytic → 2/(μ·T) = 2·T_rep/(μ·T_rep·T)    [for μT ≫ 1]

     Ratio  G_bunch / G_analytic  →  r·μ·T_rep/(1−r) = r·μ·T_rep/(1−e^{−μ·T_rep})

     Numerically (μ·T_rep = 0.0462, r = 0.9549):
         Ratio → 0.9549 × 0.0462 / (1 − 0.9549) = 0.04412/0.04510 = 0.978

     → a G_discrete is ~ 2.2 % LOWER than the G_analytic for large n.
     → β · G_discrete · D² < β · G_analytic · D²

     Authors consideration: the bunched source, with long inter-bunch gaps,
     allows more sub-lethal repair than a continuous source at the same
     average dose rate.

 ─────────────────────────────────────────────────────────────────────────────
 5.  CONTINUOUS AND CHOPPED (PULSED) SOURCES

     • Continuous: Ḋ = const → G_discrete ≡ G_analytic
     • Chopped:    μ·T_chop ≈ 10⁻¹¹ → identical to continuous at Ḋ_avg.
       The discrete stepper with any Δt ≫ T_chop gives the same result.

 ─────────────────────────────────────────────────────────────────────────────
 REFERENCES
 [1]  ICRU Report 37 (1984) — Stopping powers for positrons
 [2]  ICRU Report 44 (1989) — Tissue substitutes; liver composition
 [3]  Berger & Seltzer (1964) NASA SP-3012
 [4]  Son et al. (2013) Radiat Oncol 8:61 — α/β = 2 Gy liver
 [5]  Roberts & Holt (1988) Radiat Res 113:51
 [6]  Dale (1985) Br J Radiol 58:515 — G factor, T_half
 [7]  Lea & Catcheside (1942) J Genet 44:216
 [8]  Thames (1985) Int J Radiat Biol 47:319
 [9]  Pratx & Fung (2019) Clin Cancer Res 25:3414 — FLASH O₂ model
 [10] Bourhis et al. (2019) Radiother Oncol 139:18
 [11] Vozenin et al. (2019) Clin Cancer Res 25:35
 [12] Montay-Gruel et al. (2021) PNAS 118:e2015521118
 [13] Jansen et al. (2021) Int J Radiat Oncol Biol Phys 111:S90
 [14] Limoli & Vozenin (2023) Annu Rev Cancer Biol 7:1-21
 [15] Browne & Firestone (1986) Table of Radioactive Isotopes — ²²Na decay
 [16] Bé et al. (2004) Monographies BIPM-5 Vol.1 — mean β⁺ energy ²²Na
 [17] Liu et al. (2004) Nucl Instrum Methods B — W mesh efficiency 1.2×10⁻³
 [18] Schultz & Lynn (1988) Rev Mod Phys 60:701 — slow positron beam review
 [19] Schultz et al. (1991) Rare-gas moderator — Ar/Kr/Xe efficiencies
 [20] NIST ESTAR — https://physics.nist.gov/Star
 [21] NNDC NuDat 3 — https://www.nndc.bnl.gov/nudat3/
 [22] Brenner & Hall (1991) Int J Radiat Oncol Biol Phys 20:181 — PDR ≡ LDR
      conditions:  τ_pulse ≪ T_half  and  Δ_rep ≲ T_half
 [23] Fowler & Mount (1992) Int J Radiat Oncol Biol Phys 23:661 — PDR design
 [24] Sachs et al. (1998) Radiat Res 150:83 — G factor universality in LQ
 [25] Dale (1985) Br J Radiol 58:515 — G factor, repair kinetics [same as [6]]
 [26] Thames (1985) Int J Radiat Biol 47:319 — G, fractionation [same as [8]]
═══════════════════════════════════════════════════════════════════════════════
"""

import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import warnings
import pandas as pd
warnings.filterwarnings("ignore")

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  USER CONFIGURATION                                                      ║
# ╚══════════════════════════════════════════════════════════════════════════╝

ENERGIES_KEV    = [5, 10, 15, 20, 50, 100]
E_MAX_NA22  = 545.7    # [keV]  β⁺ endpoint energy
E_MEAN_NA22 = 215.5    # [keV]  mean β⁺ kinetic energy [16]
# ─────────────SOURCE ACTIVITY─────────
A_NA22          = 1.85e9        # [Bq]
BR_POS          = 0.9034        # [-]
ETA_MOD         = 1.2e-3        # [-]
ETA_TRANS       = 0.50          # [-]
# ────────────CHOPPER───────────────────
F_CHOP_HZ       = 50.0e6        # [Hz]
TAU_ON_S        = 2.0e-9        # [s]
CHOPPER_TRANS   = TAU_ON_S * F_CHOP_HZ   # = 0.10
# ────────────BUNCHER────────────────────
N_BUNCH         = 1.0e7        # [e⁺/bunch]
TAU_BUNCH       = 10.0e-9       # [s] # BUNCHER PULSE DURATION
#F_REP           = 8.3e-3        # [Hz]   → T_rep = 120.5 s
ETA_BUNCH       = 0.1           # [-]
# ─────────────SPOT SIZE─────────────────
FWHM_MM         = 2.0           # [mm]
# ────────────LQ PARAMETERS LIVER────────
ALPHA           = 0.10          # [Gy⁻¹]
BETA            = 0.05          # [Gy⁻²]
T_HALF_S        = 0.5 * 3600.0  # [s]
# ───────────FLASH REGIME────────────────
FLASH_ONSET     = 40.0          # [Gy/s]
FLASH_SAT       = 120.0         # [Gy/s]
FLASH_H_MIN     = 0.10          # [-]
# ───────────DOSE THRESHOLDS─────────────
D_MAX           = 10.0          # [Gy]
D_THRESHOLD     = 1.0           # [Gy]

# ╔══════════════════════════════════════════════════════════════════════════╗
# ── Derived source quantities ───────────────────────────────────────────────
# ╚══════════════════════════════════════════════════════════════════════════╝

ETA_EFF         = ETA_MOD * ETA_TRANS #1E-4
R_POSITRON      = BR_POS * ETA_EFF * A_NA22 #EFFECTIVE ACTIVITY, AFTER LOSSES ~1E6
AREA_MM2        = np.pi * (FWHM_MM / 2.0)**2 #BEAM WAIST
N_PER_PERIOD    = R_POSITRON / F_CHOP_HZ #~2E-2
R_CHOPPED       = CHOPPER_TRANS * R_POSITRON #CHOPPER LOSSES, EFFECTIVE ACTIVITY ~1E5
A_EQUIV_BUNCH   = R_POSITRON * ETA_BUNCH #BUNCHER LOSSES, EFFECTIVE ACTIVITY ~1E5
T_REP       = N_BUNCH / A_EQUIV_BUNCH
F_REP           = 1.0 / T_REP    # [s]  bunch repetition period
DC_EFF_BUNCH    = TAU_BUNCH * F_REP

# ── Physical constants ────────────────────────────────────────────────────
r_e    = 2.8179e-15 # classical electron radius [m]
m_e    = 9.1093826e-31 #electron mass [kg]
c      = const.c #speed of light [km/s]
N_A    = const.Avogadro #Na → Avogadro
e_J    = 1.602176634e-19 #eV/J
E_rest = m_e * c**2 / e_J # kg*(km/s)^2/(eV/J)
CONV   = 1.602176634e-8 #conversion factor for u.m. consistency

# ── Liver tissue (ICRU 44) ────────────────────────────────────────────────
_w  = np.array([0.102, 0.139, 0.030, 0.716,
                0.002, 0.003, 0.003, 0.002, 0.003]) #relative weights of elements → [H, C, N, O, Na, P, S, Cl, K]
_Z  = np.array([1, 6, 7, 8, 11, 15, 16, 17, 19], dtype=float) #atomic numbers
_A  = np.array([1.0079, 12.011, 14.0067, 15.994,
                22.98977, 30.97376, 32.06, 35.453, 39.0983]) #atomic mass
_I  = np.array([19.2, 81.0, 82.0, 106.0,
                149.0, 195.5, 180.0, 180.0, 190.0]) #mean excitation energy https://physics.nist.gov/cgi-bin/Star/compos.pl?matno=007
_ZA_i         = _Z / _A
ZAratio_liver = float(np.sum(_w * _ZA_i))
Imean_liver   = float(np.exp(
    np.sum(_w * _ZA_i * np.log(_I)) / ZAratio_liver)) #Bragg's additivity rule

# ── Repair rate ───────────────────────────────────────────────────────────
MU = np.log(2.0) / T_HALF_S     # [s⁻¹] #sublethal radiation damage repair rate

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  GLOBAL PLOT STYLE DEFINITION                                            ║
# ╚══════════════════════════════════════════════════════════════════════════╝

plt.rcParams.update({
    "font.family"      : "DejaVu Sans",
    "font.size"        : 10.0,
    "axes.grid"        : True,
    "grid.linestyle"   : ":",
    "grid.alpha"       : 0.38,
    "figure.dpi"       : 140,
    "axes.spines.top"  : False,
    "axes.spines.right": False,
    "xtick.direction"  : "in",
    "ytick.direction"  : "in",
    "legend.framealpha": 0.95,
    "legend.edgecolor" : "lightgray",
})

def _make_visual_style(energies):
    n   = len(energies)
    pts = np.linspace(0.10, 0.90, max(n, 2))
    cmap = plt.cm.plasma
    LS  = ["-", "--", "-.", ":", (0,(3,1,1,1)), (0,(5,2))]
    LW  = [2.2, 2.0, 1.8, 1.6, 1.8, 1.6]
    colors     = {E: cmap(pts[i])       for i, E in enumerate(energies)}
    linestyles = {E: LS[i % len(LS)]    for i, E in enumerate(energies)}
    linewidths = {E: LW[i % len(LW)]    for i, E in enumerate(energies)}
    return colors, linestyles, linewidths

PALETTE, LINESTYLE, LW_map = _make_visual_style(ENERGIES_KEV)
VALIDITY_MIN = 2.55

def _energy_handles():
    return [Line2D([0],[0], color=PALETTE[E], ls=LINESTYLE[E], lw=2.2,
                   label=f"{E} keV" + ("  ⚠" if E == 3 else ""))
            for E in ENERGIES_KEV]

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  STOPPING POWER AND LOREENTZ FACTOR                                      ║
# ╚══════════════════════════════════════════════════════════════════════════╝

# ══════════════════════════════════════════════════════════════════════════════
# PHYSICS FUNCTIONS: positron stopping power and Lorentz factor
# ══════════════════════════════════════════════════════════════════════════════

def stopping_power(E_keV):
    E_keV = np.asarray(E_keV, dtype=float)
    tau   = E_keV * 1e3 / E_rest
    beta  = np.sqrt(tau * (tau + 2.0)) / (tau + 1.0)
    Fp    = (2.0 * np.log(2.0)
             - (beta**2 / 12.0) * (23.0
                                    + 14.0 / (tau + 2.0)
                                    + 10.0 / (tau + 2.0)**2
                                    + 4.0  / (tau + 2.0)**3))
    bracket   = (2.0 * np.log(E_keV * 1e3 / Imean_liver)
                 + np.log(1.0 + tau / 2.0) + Fp)
    prefactor = (2.0 * np.pi * r_e**2 * m_e * c**2 * N_A * ZAratio_liver) / beta**2
    return prefactor * bracket * 6.24150907e16

def beta_gamma(E_keV):
    tau  = np.asarray(E_keV, dtype=float) * 1e3 / E_rest
    beta = np.sqrt(tau * (tau + 2.0)) / (tau + 1.0)
    return beta / np.sqrt(1.0 - beta**2)

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  DOSE RATES COMPUTATIONS                                                 ║
# ╚══════════════════════════════════════════════════════════════════════════╝

def dose_rates(E_keV):
    S = float(stopping_power(E_keV))
    Phi_cont  = R_POSITRON / AREA_MM2
    DR_cont   = S * Phi_cont * CONV
    DR_p_avg  = CHOPPER_TRANS * DR_cont
    Phi_b_avg  = A_EQUIV_BUNCH / AREA_MM2
    Phi_b_inst = N_BUNCH / (TAU_BUNCH * AREA_MM2)
    DR_b_avg   = S * Phi_b_avg  * CONV
    DR_b_peak  = S * Phi_b_inst * CONV
    return {
        "cont_avg": DR_cont,
        "p_avg"   : DR_p_avg,
        "p_peak"  : DR_cont,
        "b_avg"   : DR_b_avg,
        "b_peak"  : DR_b_peak,
    }

# Create a list to store data rows
dose_rate_data = []

for E in ENERGIES_KEV:
    dr = dose_rates(E)
    S = float(stopping_power(E))
    dose_rate_data.append({
        "Energy [keV]": E,
        "Stopping Power [MeV*cm2/g]": f"{S:.4f}",
        "Continuous [Gy/s]": f"{dr['cont_avg']:.4e}",
        "Pulsed Avg [Gy/s]": f"{dr['p_avg']:.4e}",
        "Pulsed Peak [Gy/s]": f"{dr['p_peak']:.4e}",
        "Bunched Avg [Gy/s]": f"{dr['b_avg']:.4e}",
        "Bunched Peak [Gy/s]": f"{dr['b_peak']:.4e}"
    })

# Convert to DataFrame for a clean table display
dr_table = pd.DataFrame(dose_rate_data)

# Export to CSV
dr_table.to_csv('dose_rates_summary.csv', index=False)
print("Data exported to 'dose_rates_summary.csv'.")

print("SUMMARY TABLE: DOSE RATES AND STOPPING POWER PER ENERGY")
print(dr_table)

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  RADIOBIOLOGICAL MODEL                                                   ║
# ╚══════════════════════════════════════════════════════════════════════════╝

def G(T_s):
    """Analytic Lea-Catcheside G [Dale 1985] for constant dose rate."""
    muT = MU * np.asarray(T_s, dtype=float)
    return np.where(muT < 1e-8,
                    1.0 - muT / 3.0,
                    2.0 * (muT - 1.0 + np.exp(-muT)) / muT**2)

def flash_h(ddot_peak):
    dp = np.asarray(ddot_peak, dtype=float).ravel()
    h  = np.ones_like(dp)
    tr = (dp > FLASH_ONSET) & (dp <= FLASH_SAT)
    fl = dp > FLASH_SAT
    h[tr] = 1.0 - (1.0 - FLASH_H_MIN) * (dp[tr] - FLASH_ONSET) \
                                       / (FLASH_SAT - FLASH_ONSET)
    h[fl] = FLASH_H_MIN
    return float(h[0]) if h.size == 1 else h

def SF_curve(D_arr, ddot_avg, ddot_peak):
    """Analytic SF"""
    D   = np.asarray(D_arr, dtype=float)
    T_s = np.where(D > 0, D / ddot_avg, 0.0)
    Gv  = G(T_s)
    h   = flash_h(ddot_peak)
    return np.exp(-ALPHA * D - BETA * h * Gv * D**2)

# ── Pre-compute dose rates and analytic SF  ────────────────
DR   = {E: dose_rates(E) for E in ENERGIES_KEV}
N_PTS = 1000 #changed from 800
D_arr = np.linspace(0.0, D_MAX, N_PTS)

SF_v6 = {E: {
    "cont"   : SF_curve(D_arr, DR[E]["cont_avg"], DR[E]["cont_avg"]),
    "pulsed" : SF_curve(D_arr, DR[E]["p_avg"],    DR[E]["p_peak"]),
    "bunched": SF_curve(D_arr, DR[E]["b_avg"],     DR[E]["b_peak"]),
} for E in ENERGIES_KEV}

# ══════════════════════════════════════════════════════════════════════════════
# DISCRETE TIME-STEP ENGINE FOR G FACTOR AND SURVIVAL FRACTION
# ══════════════════════════════════════════════════════════════════════════════


# w -> quantity of sub-lethal damage present in a cell at a given moment
def _memory_step(w, dD, decay):
    """
    Exact single-step update of the sub-lethal-damage memory variable w:
        w_{n+1} = exp(−μΔt)·(w_n + ΔD_n)
    Returns new w.
    """
    return decay * (w + dD)

def discrete_G_continuous(D_target, ddot_avg, dt_s, N_MIN=200):
    """
    Compute G_discrete for CONTINUOUS irradiation at constant dose rate.

    Should reproduce G_analytic(T) = 2[μT−1+e^{−μT}]/(μT)² to machine
    precision when N is large (validated in console output).
    N_MIN=200 ensures convergence even for short irradiations (T ≪ T_half).

    Parameters
    ----------
    D_target : total dose [Gy]
    ddot_avg : constant dose rate [Gy/s]
    dt_s     : preferred time step [s]
    N_MIN    : minimum number of steps (default 200)
    """
    T_total  = D_target / ddot_avg
    N        = max(int(round(T_total / dt_s)), N_MIN)
    dt_exact = T_total / N          # adjust Δt so steps divide T exactly
    dD       = ddot_avg * dt_exact  # dose per step [Gy]
    decay    = np.exp(-MU * dt_exact)

    w  = 0.0
    E2 = 0.0
    for _ in range(N):
        E2 += 2.0 * dD * w
        w   = _memory_step(w, dD, decay)

    G_disc = E2 / D_target**2 if D_target > 0 else 1.0
    return G_disc

def discrete_SF_bunched_array(D_arr_target, E_keV):
    """
    Compute SF_discrete for the BUNCHED source over an array of dose values.

    Algorithm
    ---------
    Each bunch is modelled as an INSTANTANEOUS dose increment
        ΔD_b = Ḋ_b_peak · τ_b   [Gy]
    delivered at t = 0, T_rep, 2·T_rep, ...

    The user's "Δt = T_rep / 3" is implemented by splitting each bunch
    period into 3 sub-steps of length T_rep/3:
        sub-step 0: dose ΔD_b  (bunch fires at the START of the period)
        sub-steps 1,2: dose 0  (gap; only repair)

    This is mathematically equivalent to using Δt = T_rep with a single
    instantaneous dose per period — the extra zero-dose sub-steps only
    apply the repair decay (verified numerically in the console).

    Parameters
    ----------
    D_arr_target : 1-D array of dose values [Gy] (must be monotone)
    E_keV        : positron beam energy [keV]

    Returns
    -------
    G_arr_disc : G_discrete at each D value
    SF_arr_disc: SF_discrete at each D value
    """
    D_arr_target = np.asarray(D_arr_target, dtype=float)
    dr    = dose_rates(E_keV)
    dD_b  = dr["b_peak"] * TAU_BUNCH        # dose per bunch [Gy]
    h_val = flash_h(dr["b_peak"])

    if dD_b <= 0:
        # degenerate: no bunched dose
        sf = np.exp(-ALPHA * D_arr_target)
        return np.ones_like(D_arr_target), sf

    # ── Sub-step decay factors (T_rep / 3 sub-step) ───────────────────────
    dt_sub  = T_REP / 3.0
    decay_s = np.exp(-MU * dt_sub)   # decay per sub-step (40 s)

    # ── Build output arrays ───────────────────────────────────────────────
    G_out  = np.zeros_like(D_arr_target)
    SF_out = np.zeros_like(D_arr_target)

    # ── Single-pass accumulation ──────────────────────────────────────────
    w    = 0.0    # sub-lethal-damage memory
    E2   = 0.0    # quadratic damage accumulator
    D_n  = 0.0    # accumulated dose

    idx  = 0      # current index in D_arr_target
    n_total = int(D_arr_target[-1] / dD_b) + 2   # max bunches needed

    for _ in range(n_total):
        # ── Sub-step 0: BUNCH fires ───────────────────────────────────────
        E2  += 2.0 * dD_b * w         # interaction of new dose with past
        D_n += dD_b
        w    = _memory_step(w, dD_b, decay_s)   # decay T_rep/3

        # ── Sub-steps 1, 2: zero dose, only repair ────────────────────────
        w = _memory_step(w, 0.0, decay_s)
        w = _memory_step(w, 0.0, decay_s)

        # ── Record SF at all target doses that have been reached ──────────
        while idx < len(D_arr_target) and D_n >= D_arr_target[idx]:
            # Optionally interpolate for the fractional last bunch
            frac = D_arr_target[idx] / D_n if D_n > 0 else 1.0
            # For the fractional dose use linear interpolation of E2
            E2_interp = E2 * (D_arr_target[idx] / D_n)**2  if D_n > 0 else E2
            G_out[idx]  = (E2_interp / D_arr_target[idx]**2
                           if D_arr_target[idx] > 0 else 1.0)
            SF_out[idx] = np.exp(-ALPHA * D_arr_target[idx]
                                 - BETA * h_val * E2_interp)
            idx += 1
        if idx >= len(D_arr_target):
            break

        # Fill any remaining points (doses above last bunch)
    if idx < len(D_arr_target):
        G_out[idx:]  = G_out[idx - 1] if idx > 0 else 0.0
        SF_out[idx:] = SF_out[idx - 1] if idx > 0 else 1.0

    return G_out, SF_out

    def analytic_G_bunched_formula(n, r):
        """
        Closed-form G for n equal instantaneous bunches separated by T_rep,
        with inter-bunch repair factor r = e^{-μ·T_rep}.  Equation (★★★★).

        G_bunch(n) = 2r(1-r) · [n(1-r) - (1-r^n)] / [n²(1-r)²]
               = 2r / [n(1-r)] · [1 - (1-r^n)/(n(1-r))]
        """
        if n == 0:
           return 1.0
        if abs(1.0 - r) < 1e-14:          # r ≈ 1 (no repair): G_bunch → (n-1)/n
           return (n - 1.0) / n
        rn  = r**n
        num = 2.0 * r * (1.0 - r) * (n * (1.0 - r) - (1.0 - rn))
        den = (n * (1.0 - r))**2
        return num / den

# ══════════════════════════════════════════════════════════════════════════════
# PRE-COMPUTE v7 DISCRETE SF AND BUNCHED DATA
# ══════════════════════════════════════════════════════════════════════════════
T_HI_BUNCH = 43200    # [s]

r_rep = np.exp(-MU * T_REP) # inter-bunch repair factor

# 1. Closed-form vectorised SF for Bunched Source
def SF_bunch_fast(n_arr, dD_b, r, h_val):
    n   = np.asarray(n_arr, dtype=float)
    D   = n * dD_b
    if abs(1.0 - r) < 1e-14:
        E2 = dD_b**2 * n * (n - 1.0)
    else:
        log_rn = n * np.log(r)
        rn     = np.where(log_rn > -700.0, np.exp(log_rn), 0.0)
        E2     = (2.0 * dD_b**2 * r * (n * (1.0 - r) - (1.0 - rn)) / (1.0 - r)**2)
    E2 = np.maximum(E2, 0.0)
    return np.exp(-ALPHA * D - BETA * h_val * E2)

# 2. Pre-compute discrete SF for all energies (D_MAX grid)
D_plot  = np.linspace(0.0, D_MAX, N_PTS)
D_plot_nz = D_plot.copy(); D_plot_nz[0] = 1e-6

G_disc_plot  = {}
SF_disc_plot = {}

for E in ENERGIES_KEV:
    G_d, SF_d = discrete_SF_bunched_array(D_plot_nz, E)
    G_disc_plot[E]  = G_d
    SF_disc_plot[E] = SF_d

# 3. Pre-compute bunched data up to T_HI_BUNCH using log-sampling to save RAM
print("  Pre-computing bunched SF (sparse log-sampling) … ", end="", flush=True)
_BUNCH_DATA = {}
N_SAMPLES = 5000

for E in ENERGIES_KEV:
    dr    = DR[E]
    dD_b  = dr["b_avg"] / F_REP
    h_val = float(flash_h(dr["b_peak"]))
    n_total_max = int(T_HI_BUNCH / T_REP) + 2

    n_indices = np.unique(np.logspace(0, np.log10(n_total_max), N_SAMPLES).astype(int))
    n_arr_b = np.concatenate(([0], n_indices)).astype(float)
    T_arr_b = n_arr_b * T_REP
    D_arr_b = n_arr_b * dD_b

    sf_arr_b    = SF_bunch_fast(n_arr_b, dD_b, r_rep, h_val)
    sf_arr_b[0] = 1.0

    _BUNCH_DATA[E] = {
        "n"  : n_arr_b,
        "T"  : T_arr_b,
        "D"  : D_arr_b,
        "SF" : sf_arr_b,
        "dD" : dD_b,
    }
print("done.")

# ══════════════════════════════════════════════════════════════════════════════
# SHARED TIME-AXIS CONFIG  (Figs A, B, C)
# t_lo = D_START/Ḋ_max: fastest curve starts within 0.01% of D=0 / t=0
# ══════════════════════════════════════════════════════════════════════════════
DIV_MIN  = 60.0
D_START  = 1e-3
_dr3_cont = DR[ENERGIES_KEV[0]]["cont_avg"]
_dr3_puls = DR[ENERGIES_KEV[0]]["p_avg"]
_T_LO_S = {"cont": D_START/_dr3_cont, "pulsed": D_START/_dr3_puls, "bunched": T_REP/10.0}
_T_HI_S = {"cont": 6e3, "pulsed": 6e4, "bunched": T_HI_BUNCH}
_SRC_CFG = [
    ("cont",    "Continuous", "cont_avg", "cont_avg"),
    ("pulsed",  "Pulsed",     "p_avg",    "p_peak"),
    ("bunched", "Bunched",    "b_avg",    "b_peak"),
]

# Reference SF levels
_SF_LEVELS = [0.50, 0.10, 0.01]
_SF_COLORS = ["#888888", "#888888", "#aaaaaa"]
_SF_LS     = ["--", ":", "-."]
_SF_LABELS = ["SF = 50%", "SF = 10%", "SF = 1%"]


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE A — EXPOSURE TIME (log) vs ABSORBED DOSE  [minutes, all panels]
# ══════════════════════════════════════════════════════════════════════════════
# global y-range: union of all three panels (smallest t_lo → largest t_hi)
_A_Y_LO = min(_T_LO_S.values()) / DIV_MIN          # ≈ 5e-5 min (cont)
_A_Y_HI = (_BUNCH_DATA[ENERGIES_KEV[-1]]["T"][      # slowest energy at D_MAX
    np.searchsorted(_BUNCH_DATA[ENERGIES_KEV[-1]]["D"], D_MAX)] * 3) / DIV_MIN

fig_A, axes_A = plt.subplots(1, 3, figsize=(17, 6.2),
                              sharey=True,
                              gridspec_kw={"wspace": 0.20})
fig_A.subplots_adjust(left=0.07, right=0.975, top=0.84, bottom=0.27)

for ax, (src, label, da_k, dp_k) in zip(axes_A, _SRC_CFG):
    t_lo_s = _T_LO_S[src]; t_hi_s = _T_HI_S[src]
    ax.axvline(D_THRESHOLD, color="dimgray", ls="-.", lw=1.0, zorder=2)

    for E in ENERGIES_KEV:
        col  = PALETTE[E]; ls = LINESTYLE[E]; lw = LW_map[E]
        davg = DR[E][da_k]

        if src == "bunched":
            bd  = _BUNCH_DATA[E]
            T_s = bd["T"]; D_s = bd["D"]
            mask = (D_s <= D_MAX * 1.02)
            T_pl = T_s.copy(); T_pl[0] = t_lo_s
            ax.step(D_s[mask], T_pl[mask] / DIV_MIN,
                    color=col, lw=lw, where="post", label=f"{E} keV",
                    solid_joinstyle="miter")
        else:
            D_plt = np.linspace(0, D_MAX, 800)
            T_plt = D_plt / davg
            mask  = (T_plt >= t_lo_s) & (T_plt <= t_hi_s * 1.05)
            ax.semilogy(D_plt[mask], T_plt[mask] / DIV_MIN,
                        color=col, ls=ls, lw=lw, label=f"{E} keV")
            t_thr = D_THRESHOLD / davg
            if t_lo_s <= t_thr <= t_hi_s:
                ax.plot(D_THRESHOLD, t_thr / DIV_MIN,
                        "x", color=col, ms=8, mew=1.8, zorder=6)

    ax.set_xlim(0, D_MAX)
    ax.set_xlabel("Absorbed dose  $D$  [Gy]", fontsize=10.5)
    ax.set_title(label, fontsize=11, fontweight="bold", pad=4)

axes_A[0].set_yscale("log")             # sharey propagates log scale to all
axes_A[0].set_ylim(_A_Y_LO, _A_Y_HI)   # single call sets all three via sharey
axes_A[0].set_ylabel("Exposure time  $t$  [min]  (log)", fontsize=10.5)

fig_A.legend(handles=_energy_handles() + [
    Line2D([0],[0], color="dimgray", ls="-.", lw=1.0,
           label=f"$D$ = {D_THRESHOLD} Gy  threshold"),
    Line2D([0],[0], color="k", marker="x", ls="none",
           ms=8, mew=1.8, label="× time to threshold"),
    Line2D([0],[0], color="k", ls="-", lw=2.0,
           label=rf"Bunched staircase  ($T_{{rep}}$ = {T_REP:.0f} s)"),
], loc="lower center", ncol=4, fontsize=8.0, bbox_to_anchor=(0.5, 0.00))
fig_A.suptitle(
    "Exposure time  $t$  [min, log]  vs  Absorbed dose  $D$  — all three regimes\n"
    r"$t_\mathrm{lo} = D_\mathrm{start}/\dot{D}_\mathrm{max}$: "
    r"fastest curve (3 keV) visible from $D \approx 0$  ·  all panels in minutes",
    fontsize=9.5, y=0.995)
# plt.savefig("/mnt/user-data/outputs/figA_time_vs_dose.png",
#             dpi=155, bbox_inches="tight")
# print("  Saved: figA_time_vs_dose.png")



# ══════════════════════════════════════════════════════════════════════════════
# FIGURE B — SURVIVING FRACTION  vs  EXPOSURE TIME (log, minutes)
# ──────────────────────────────────────────────────────────────────────────────
#  ALIGNMENT: sharey=True + dynamic t_lo
#  ─────────────────────────────────────────────────────────────────────────────
#  (a) sharey=True  →  single shared y-axis; identical ticks/limits locked.
#  (b) t_lo = D_start / Ḋ_avg_max   with D_start = 1 mGy
#      ensures the fastest energy (3 keV) has SF > 0.9999 at the left edge.
#          Continuous:  t_lo = 0.001/0.324 = 3.1 ms  → 5.1×10⁻⁵ min
#          Pulsed:      t_lo = 0.001/0.0324 = 31 ms  → 5.1×10⁻⁴ min
#          Bunched:     t_lo = T_rep/10 = 12 s         → 0.20     min
#  (c) wspace = 0 (panels touch; sharey makes the shared axis seamless).
# ══════════════════════════════════════════════════════════════════════════════

DIV_MIN  = 60.0
D_START  = 1e-3   # Gy → SF = exp(−0.1×0.001) = 0.9999

_dr3 = DR[ENERGIES_KEV[0]]   # E = 3 keV  →  highest dose rate
_T_LO_S = {
    "cont"   : D_START / _dr3["cont_avg"],
    "pulsed" : D_START / _dr3["p_avg"],
    "bunched": T_REP / 10.0,
}
_T_HI_S = {
    "cont"   : 6e3,
    "pulsed" : 6e4,
    "bunched": T_HI_BUNCH,
}
_FIGB_CFG = [
    ("cont",    "Continuous",  "cont_avg",  "cont_avg"),
    ("pulsed",  "Pulsed",      "p_avg",     "p_peak"),
    ("bunched", "Bunched",     "b_avg",     "b_peak"),
]

fig_B, axes_B = plt.subplots(
    1, 3, figsize=(17, 7.0),
    sharey=True,
    gridspec_kw={"wspace": 0.20})
fig_B.subplots_adjust(left=0.07, right=0.975, top=0.83, bottom=0.29)

for ax, (src, label, da_k, dp_k) in zip(axes_B, _FIGB_CFG):
    t_lo_s = _T_LO_S[src]
    t_hi_s = _T_HI_S[src]

    for sfv, col_r, ls_r in zip(_SF_LEVELS, _SF_COLORS, _SF_LS):
        ax.axhline(sfv, color=col_r, lw=0.85, ls=ls_r, alpha=0.60, zorder=1)

    for E in ENERGIES_KEV:
        col  = PALETTE[E]; ls = LINESTYLE[E]; lw = LW_map[E]
        davg = DR[E][da_k]; dpeak = DR[E][dp_k]

        if src == "bunched":
            bd = _BUNCH_DATA[E]
            T_s = bd["T"]; sf_s = bd["SF"]
            T_plot = T_s.copy()
            T_plot[0] = t_lo_s          # replace T=0 with t_lo_s (inside axis)
            mask = (T_plot <= t_hi_s * 1.02) & (sf_s >= 1e-8)
            ax.step(T_plot[mask] / DIV_MIN, sf_s[mask],
                    color=col, lw=lw, where="post", label=f"{E} keV",
                    solid_joinstyle="miter")
            sf_m = sf_s[mask]; t_m = T_plot[mask]
            if sf_m[0] > 0.50 and sf_m[-1] < 0.50:
                i50 = int(np.searchsorted(-sf_m, -0.50))
                ax.plot(t_m[i50] / DIV_MIN, 0.50, "o", color=col,
                        ms=7, zorder=7, markeredgecolor="k", markeredgewidth=0.7)
        else:
            t_arr = np.logspace(np.log10(t_lo_s), np.log10(t_hi_s), 1600)
            sf_t  = SF_curve(davg * t_arr, davg, dpeak)
            ax.semilogy(t_arr / DIV_MIN, sf_t,
                        color=col, ls=ls, lw=lw, zorder=3, label=f"{E} keV")
            if sf_t[0] > 0.50 > sf_t[-1]:
                t50 = float(np.interp(0.50, sf_t[::-1], t_arr[::-1]))
                ax.plot(t50 / DIV_MIN, 0.50, "o", color=col,
                        ms=7, zorder=7, markeredgecolor="k", markeredgewidth=0.7)
            t_1Gy = D_THRESHOLD / davg if davg > 0 else 1e99
            if t_lo_s <= t_1Gy <= t_hi_s:
                sf_1 = float(np.interp(t_1Gy, t_arr, sf_t))
                ax.plot(t_1Gy / DIV_MIN, sf_1, "x", color=col,
                        ms=8, mew=1.8, zorder=7)

    ax.set_xscale("log")
    ax.set_xlim(t_lo_s / DIV_MIN, t_hi_s / DIV_MIN)
    ax.set_xlabel("Exposure time  $t$  [min]  (log)", fontsize=10.5)
    ax.set_title(label, fontsize=11, fontweight="bold", pad=4)

    for sfv, lbl_r, col_r in zip(_SF_LEVELS, _SF_LABELS, _SF_COLORS):
        ax.text(ax.get_xlim()[1] * 0.82, sfv * 1.15,
                lbl_r, ha="right", fontsize=7, color=col_r)

    _ann = {
        "bunched": (r"v7 step fn  (SF const between bunches)" "\n"
                    r"FLASH $h=0.10$,  $G\approx 0$",
                    "#6a3d9a", "#f2e6ff"),
        "pulsed" : (r"Analytic $G$ — exact" "\n"
                    r"($\mu T_\mathrm{chop}=7.7\!\times\!10^{-12}$)",
                    "#1f77b4", "#e6f2ff"),
        "cont"   : (r"Analytic $G(T)$ — exact" "\n"
                    r"(const $\dot{D}$, Dale 1985)",
                    "#2ca02c", "#e6f9ee"),
    }
    txt, ec, fc = _ann[src]
    ax.text(0.03, 0.10, txt, transform=ax.transAxes,
            fontsize=7.5, color=ec,
            bbox=dict(boxstyle="round,pad=0.3", fc=fc, ec=ec, lw=0.8, alpha=0.9))

axes_B[0].set_ylabel("Cell surviving fraction  SF  (log)", fontsize=10.5)
axes_B[0].set_ylim(3e-4, 1.6)      # applied to all via sharey

fig_B.legend(
    handles=_energy_handles() + [
        Line2D([0],[0], color="k", marker="o", ls="none", ms=7,
               markeredgecolor="k", label="● LD₅₀  (SF = 50 %)"),
        Line2D([0],[0], color="k", marker="x", ls="none",
               ms=8, mew=1.8, label=f"× $D$ = {D_THRESHOLD:.0f} Gy  threshold"),
        Line2D([0],[0], color="k", ls="-", lw=2.0,
               label=rf"Bunched step fn  ($T_\mathrm{{rep}}$ = {T_REP:.0f} s)"),
    ],
    loc="lower center", ncol=5, fontsize=8.2, bbox_to_anchor=(0.5, 0.00))

fig_B.suptitle(
    r"Surviving fraction  SF  vs  Exposure time  $t$  [min, log]  — v7"
    "\n"
    r"sharey · all panels in minutes · "
    r"$t_\mathrm{lo}$ = D$_\mathrm{start}$/\dot{D}$_\mathrm{max}$ ensures SF = 1 at left edge  ·  "
    rf"$\alpha={ALPHA}$ Gy⁻¹,  $\beta={BETA}$ Gy⁻²,  $T_{{1/2}}={T_HALF_S/3600:.1f}$ h",
    fontsize=9.2, y=0.995)

# plt.savefig("/mnt/user-data/outputs/figB_SF_vs_time.png",
#             dpi=155, bbox_inches="tight")
# print("  Saved: figB_SF_vs_time.png")

# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════
# FIGURE C — ACCUMULATED DOSE  vs  EXPOSURE TIME (log) [minutes, sharey]
# sharey=True: dose axis identical across all panels (0 → D_MAX)
# t_lo dynamic: fastest curve starts at D ≈ 0 at left edge
# ══════════════════════════════════════════════════════════════════════════════
fig_C, axes_C = plt.subplots(1, 3, figsize=(17, 6.2),
                              sharey=True,
                              gridspec_kw={"wspace": 0.20})
fig_C.subplots_adjust(left=0.07, right=0.975, top=0.84, bottom=0.27)

for ax, (src, label, da_k, dp_k) in zip(axes_C, _SRC_CFG):
    t_lo_s = _T_LO_S[src]; t_hi_s = _T_HI_S[src]
    ax.axhline(D_THRESHOLD, color="dimgray", ls="-.", lw=1.1, zorder=2)

    for E in ENERGIES_KEV:
        col  = PALETTE[E]; ls = LINESTYLE[E]; lw = LW_map[E]
        davg = DR[E][da_k]

        if src == "bunched":
            bd  = _BUNCH_DATA[E]
            T_s = bd["T"]; D_s = bd["D"]
            # cap at D_MAX so shared y-axis works cleanly
            mask = (T_s <= t_hi_s * 1.02) & (D_s <= D_MAX * 1.02)
            T_pl = T_s.copy(); T_pl[0] = t_lo_s
            ax.step(T_pl[mask] / DIV_MIN, D_s[mask],
                    color=col, lw=lw, where="post", label=f"{E} keV",
                    solid_joinstyle="miter")
        else:
            t_lin = np.logspace(np.log10(t_lo_s), np.log10(t_hi_s), 800)
            ax.plot(t_lin / DIV_MIN, davg * t_lin,
                    color=col, ls=ls, lw=lw, label=f"{E} keV")
            t_thr = D_THRESHOLD / davg
            if t_lo_s <= t_thr <= t_hi_s:
                ax.axvline(t_thr / DIV_MIN, color=col, ls=":", lw=0.9, alpha=0.65)

    ax.set_xscale("log")
    ax.set_xlim(t_lo_s / DIV_MIN, t_hi_s / DIV_MIN)
    ax.set_xlabel("Exposure time  $t$  [min]  (log)", fontsize=10.5)
    ax.set_title(label, fontsize=11, fontweight="bold", pad=4)
    for i, E_a in enumerate([ENERGIES_KEV[0], ENERGIES_KEV[-1]]):
        ax.text(0.03, 0.96 - i*0.09,
                rf"$\dot{{D}}$ ({E_a} keV) = {DR[E_a][da_k]:.2e} Gy/s",
                transform=ax.transAxes, fontsize=7.2,
                color=PALETTE[E_a], va="top")

axes_C[0].set_ylabel("Accumulated dose  $D$  [Gy]", fontsize=10.5)
axes_C[0].set_ylim(0, D_MAX * 1.05)    # propagates to all via sharey

fig_C.legend(handles=_energy_handles() + [
    Line2D([0],[0], color="dimgray", ls="-.", lw=1.1,
           label=f"$D$ = {D_THRESHOLD} Gy  threshold"),
    Line2D([0],[0], color="gray",    ls=":",  lw=0.9,
           label="Vertical dotted = time to threshold (cont./pulsed)"),
    Line2D([0],[0], color="k",       ls="-",  lw=2.0,
           label=rf"Bunched staircase  ($T_{{rep}}$ = {T_REP:.0f} s)"),
], loc="lower center", ncol=4, fontsize=8.0, bbox_to_anchor=(0.5, 0.00))
fig_C.suptitle(
    r"Accumulated dose  $D$  vs  Exposure time  $t$  [min, log]  — all three regimes"
    "\n"
    r"sharey · all panels minutes · $t_\mathrm{lo}$: fastest curve starts at $D \approx 0$"
    r"  ·  bunched capped at $D_\mathrm{MAX}$ = 10 Gy for shared dose axis",
    fontsize=9.5, y=0.995)
# plt.savefig("/mnt/user-data/outputs/figC_dose_vs_time.png",
#             dpi=155, bbox_inches="tight")
# print("  Saved: figC_dose_vs_time.png")

# CONSOLE SUMMARY
# ══════════════════════════════════════════════════════════════════════════════
LINE = "═" * 88
print()
print(LINE)
print("  v7 BUG FIX — BUNCHED SF WAS TRUNCATED AT D_MAX")
print(LINE)
print(f"""
  ROOT CAUSE:  n_max was set by D_MAX / ΔD_b.  At D_MAX = 10 Gy the
               bunched SF reaches only exp(−α·D_MAX) = exp(−1) ≈ 0.37.
               FLASH effect (h = 0.10) suppresses the β term,
               G ≈ 0 suppresses it further → SF_min ≈ 0.37 within D_MAX.
               Lower 60% of log y-axis was empty.

  FIX:         n_max now covers the full time window T_HI = {T_HI_BUNCH:.2e} s
               (≈ {T_HI_BUNCH/86400:.0f} days).  Closed-form vectorised formula
               used for efficiency (no iteration).

  PHYSICS:     Bunched source has SF_min determined by time (not D_MAX).
               At T_HI ≈ 400 d, slowest energy (100 keV):
                 D ≈ {DR[ENERGIES_KEV[-1]]['b_avg'] * T_HI_BUNCH:.1f} Gy
                 → SF ≈ exp(−{ALPHA * DR[ENERGIES_KEV[-1]]['b_avg'] * T_HI_BUNCH:.1f})
                       = {np.exp(-ALPHA * DR[ENERGIES_KEV[-1]]['b_avg'] * T_HI_BUNCH):.3f}
               Fastest energy (3 keV):
                 D ≈ {DR[ENERGIES_KEV[0]]['b_avg'] * T_HI_BUNCH:.0f} Gy → SF → 0
""")
print("  FIGURES SAVED:")
for f in ["figA_time_vs_dose.png", "figB_SF_vs_time.png", "figC_dose_vs_time.png"]:
    print(f"    {f}")
print(LINE)

# ═══════════════════════════════════════════════════════════════════════════
# PRE-COMPUTE ALL DOSE RATES AND SURVIVAL CURVES
# ═══════════════════════════════════════════════════════════════════════════
DR   = {E: dose_rates(E) for E in ENERGIES_KEV}
N_PTS = 800
D_arr = np.linspace(0.0, D_MAX, N_PTS)

SF = {E: {
    "cont"   : SF_curve(D_arr, DR[E]["cont_avg"], DR[E]["cont_avg"]),
    "pulsed" : SF_curve(D_arr, DR[E]["p_avg"],    DR[E]["p_peak"]),
    "bunched": SF_curve(D_arr, DR[E]["b_avg"],     DR[E]["b_peak"]),
} for E in ENERGIES_KEV}

# ═══════════════════════════════════════════════════════════════════════════
# CONSOLE OUTPUT
# ═══════════════════════════════════════════════════════════════════════════
def print_console_summary():
    """Prints the physical parameters and radiobiological summary to the console."""
    LINE = "═" * 90
    SEP  = "─" * 90

    print(LINE)
    print("  ²²Na SOURCE MODEL")
    print(LINE)
    print(f"  Total activity             A_Na22  = {A_NA22:.2e}  Bq")
    print(f"  β⁺ branching ratio         BR_POS  = {BR_POS:.4f}  (NNDC NuDat 3 [21])")
    print(f"  β⁺ endpoint / mean energy           = {E_MAX_NA22} / {E_MEAN_NA22} keV  [16]")
    print(f"  Moderator efficiency       ETA_MOD  = {ETA_MOD:.2e}  (W mesh, Liu 2004 [17])")
    print(f"  Transport efficiency       ETA_TRANS= {ETA_TRANS:.2f}")
    print(f"  Combined efficiency        ETA_EFF  = {ETA_EFF:.4e}")
    print(f"  Positrons at sample        R_sample = {R_POSITRON:.3e}  e⁺/s  (= BR×η_mod×η_trans×A = {BR_POS}×{ETA_MOD:.2e}×{ETA_TRANS}×{A_NA22:.2e})")
    print(f"  Beam area (FWHM = {FWHM_MM} mm)             = {AREA_MM2:.4f}  mm²")
    print()
    print("  PULSED BEAM (²²Na + electrical chopper)")
    print(f"  Frequency / duty cycle              = {F_CHOP_HZ/1e6:.0f} MHz  CHOPPER_TRANS={CHOPPER_TRANS:.2f}")
    print(f"  ON-duration per period    τ_ON       = {TAU_ON_S*1e9:.2f} ns")
    print(f"  ⟨positrons per period⟩   N_per_period = {N_PER_PERIOD:.3e}  (λ«1 → 0 or 1 per period)")
    print(f"  D_chopped = {CHOPPER_TRANS:.2f} × D_continuous  (first principles, Schultz & Lynn 1988 [18])")
    print()
    print("   BUNCHED BEAM")
    print(f"  N_b = {N_BUNCH:.0e} e⁺/bunch   τ_b = {TAU_BUNCH*1e9:.0f} ns   f_rep = {F_REP:.0e} Hz")
    print(f"  A_equiv = {A_EQUIV_BUNCH:.2e} Bq   DC_eff = {DC_EFF_BUNCH:.2e}")
    print()

    print(LINE)
    print("  DOSE RATE TABLE  (all Gy/s)")
    print(LINE)
    print(f"  {'E[keV]':>7}  {'βγ':>7}  {'S/ρ':>12}  {'Ḋ_cont':>12}  "
          f"{'Ḋ_p_peak':>12}  {'Ḋ_p_avg':>12}  {'Ḋ_b_avg':>12}  {'Ḋ_b_peak':>14}")
    print(f"  {'':>7}  {'':>7}  {'[MeV·cm²/g]':>12}  {'[Gy/s]':>12}  "
          f"{'[Gy/s]':>12}  {'[Gy/s]':>12}  {'[Gy/s]':>12}  {'[Gy/s]':>14}")
    print("  " + SEP[:88])
    for E in ENERGIES_KEV:
        dr  = DR[E];  bg = beta_gamma(E);  S = float(stopping_power(E))
        tag = "  ⚠" if E < 5 else ""
        print(f"  {E:>7}  {bg:>7.4f}  {S:>12.4f}  {dr['cont_avg']:>12.4e}  "
              f"{dr['p_peak']:>12.4e}  {dr['p_avg']:>12.4e}  "
              f"{dr['b_avg']:>12.4e}  {dr['b_peak']:>14.4e}{tag}")

    print()
    print("  NOTE: Ḋ_p_peak = Ḋ_cont  (chopper does NOT compress positrons)")
    print("        Ḋ_p_avg  = DC × Ḋ_cont  (positrons LOST during OFF period)")
    print()
    print(LINE)
    print("  SURVIVAL FRACTION  AT KEY DOSES")
    print(LINE)
    print(f"  {'E[keV]':>7}  {'Source':>12}  {'h':>6}  {'G@5Gy':>8}  "
          f"{'t_safe[s]':>10}  {'SF@1Gy':>9}  {'SF@5Gy':>9}  {'SF@10Gy':>10}")
    print("  " + SEP[:80])
    for E in ENERGIES_KEV:
        dr = DR[E]
        for lbl, sk, da_k, dp_k in [
                ("Continuous", "cont",    "cont_avg", "cont_avg"),
                ("Pulsed",     "pulsed",  "p_avg",    "p_peak"),
                ("Bunched",    "bunched", "b_avg",    "b_peak")]:
            davg  = dr[da_k]; dpeak = dr[dp_k]
            G5    = float(G(5.0 / davg))
            h_val = flash_h(dpeak)
            t_s   = D_THRESHOLD / davg
            sf1   = float(SF[E][sk][np.argmin(np.abs(D_arr - 1.0))])
            sf5   = float(SF[E][sk][np.argmin(np.abs(D_arr - 5.0))])
            sf10  = float(SF[E][sk][np.argmin(np.abs(D_arr - 10.0))])
            flt   = " !" if dpeak > FLASH_ONSET else "  "
            print(f"  {E:>7}  {lbl:>12}  {h_val:>6.2f}{flt}  {G5:>8.5f}  "
                  f"{t_s:>10.1f}  {sf1:>9.5f}  {sf5:>9.5f}  {sf10:>10.6f}")
        print()
    print(LINE)

# Call the function to print the summary to console
print_console_summary()

# ═══════════════════════════════════════════════════════════════════════════

plt.show()

# export1: dose history and time required to reach 1 Gy for all three source types 

# Set export resolution
N_EXPORT_PTS = 10000

# 1. Determine a Global Common Time Axis
# We take the union of all relevant time bounds to ensure consistency
t_min = min(_T_LO_S['cont'], _T_LO_S['pulsed'], _T_LO_S['bunched'])
t_max = max(_T_HI_S['cont'], _T_HI_S['pulsed'], _T_HI_S['bunched'])

t_common = np.geomspace(t_min, t_max, N_EXPORT_PTS)

# 2. Continuous Source (using common axis)
data_cont = {'Time [s]': t_common}
for E in ENERGIES_KEV:
    data_cont[f'Dose_{E}keV [Gy]'] = DR[E]['cont_avg'] * t_common
pd.DataFrame(data_cont).to_csv('dose_history_continuous.csv', index=False)

# 3. Pulsed Source (using common axis)
data_pulsed = {'Time [s]': t_common}
for E in ENERGIES_KEV:
    data_pulsed[f'Dose_{E}keV [Gy]'] = DR[E]['p_avg'] * t_common
pd.DataFrame(data_pulsed).to_csv('dose_history_pulsed.csv', index=False)

# 4. Bunched Source (Synchronized via resampling/forward-fill)
# Since bunched data is discrete, we map it onto the common timeline using forward-fill
bunched_frames = []
for E in ENERGIES_KEV:
    df_orig = pd.DataFrame({
        'Time [s]': _BUNCH_DATA[E]['T'],
        f'Dose_{E}keV [Gy]': _BUNCH_DATA[E]['D']
    })
    bunched_frames.append(df_orig.set_index('Time [s]'))

# Combine original discrete points
dose_bunch_raw = pd.concat(bunched_frames, axis=1).sort_index()

# Reindex to the common time axis to make the files structurally identical
dose_history_bunched = dose_bunch_raw.reindex(dose_bunch_raw.index.union(t_common)).ffill().loc[t_common]
dose_history_bunched.index.name = 'Time [s]'
dose_history_bunched.to_csv('dose_history_bunched.csv')

# 5. Time to 1 Gy Summary Table
time_1gy_data = []
for E in ENERGIES_KEV:
    time_1gy_data.append({
        "Energy [keV]": E,
        "Continuous [s]": round(1.0 / DR[E]['cont_avg'], 4),
        "Pulsed [s]": round(1.0 / DR[E]['p_avg'], 4),
        "Bunched [s]": round(1.0 / DR[E]['b_avg'], 4)
    })

time_1gy_df = pd.DataFrame(time_1gy_data)
time_1gy_df.to_csv('time_to_1Gy_summary.csv', index=False)

print(f"Files generated with {N_EXPORT_PTS} synchronized time points:")
print(" - dose_history_continuous.csv (Shared Time Axis)")
print(" - dose_history_pulsed.csv     (Shared Time Axis)")
print(" - dose_history_bunched.csv    (Shared Time Axis + Forward Fill)")
print(" - time_to_1Gy_summary.csv")

print("\nTIME REQUIRED TO DEPOSIT 1 Gy [seconds]:")
print(time_1gy_df.head(len(ENERGIES_KEV)))

#export 2: time to specific SF thresholds for all three source types

# 1. Reuse the common time axis from the previous step
# (t_common and N_EXPORT_PTS are already defined in the kernel)

# 2. Continuous SF History
sf_data_cont = {'Time [s]': t_common}
for E in ENERGIES_KEV:
    dr = DR[E]
    sf_data_cont[f'SF_{E}keV'] = SF_curve(dr['cont_avg'] * t_common, dr['cont_avg'], dr['cont_avg'])
pd.DataFrame(sf_data_cont).to_csv('sf_history_continuous.csv', index=False)

# 3. Pulsed SF History
sf_data_pulsed = {'Time [s]': t_common}
for E in ENERGIES_KEV:
    dr = DR[E]
    sf_data_pulsed[f'SF_{E}keV'] = SF_curve(dr['p_avg'] * t_common, dr['p_avg'], dr['p_peak'])
pd.DataFrame(sf_data_pulsed).to_csv('sf_history_pulsed.csv', index=False)

# 4. Bunched SF History (using forward-fill to maintain the step-function nature)
bunched_sf_frames = []
for E in ENERGIES_KEV:
    df_orig = pd.DataFrame({
        'Time [s]': _BUNCH_DATA[E]['T'],
        f'SF_{E}keV': _BUNCH_DATA[E]['SF']
    })
    bunched_sf_frames.append(df_orig.set_index('Time [s]'))

sf_bunch_raw = pd.concat(bunched_sf_frames, axis=1).sort_index()
sf_history_bunched = sf_bunch_raw.reindex(sf_bunch_raw.index.union(t_common)).ffill().loc[t_common]
sf_history_bunched.index.name = 'Time [s]'
sf_history_bunched.to_csv('sf_history_bunched.csv')

print("SF history files generated (synchronized with dose history time axis):")
print(" - sf_history_continuous.csv")
print(" - sf_history_pulsed.csv")
print(" - sf_history_bunched.csv")

# --- Original Threshold Table Logic ---
sf_targets = [0.99, 0.95, 0.90, 0.50, 0.10]
sf_summary_data = []

for E in ENERGIES_KEV:
    dr = DR[E]
    for mode in ['cont', 'pulsed']:
        da_k = 'cont_avg' if mode == 'cont' else 'p_avg'
        dp_k = 'cont_avg' if mode == 'cont' else 'p_peak'
        t_search = np.geomspace(1e-3, 1e6, 50000)
        sf_vals = SF_curve(dr[da_k] * t_search, dr[da_k], dr[dp_k])
        row = {"Energy [keV]": E, "Mode": mode.capitalize()}
        for target in sf_targets:
            if sf_vals[-1] <= target:
                t_target = np.interp(target, sf_vals[::-1], t_search[::-1])
                row[f"Time to SF={target}"] = f"{t_target:.2f}s"
            else:
                row[f"Time to SF={target}"] = ">1e6s"
        sf_summary_data.append(row)

    bd = _BUNCH_DATA[E]
    row_b = {"Energy [keV]": E, "Mode": "Bunched"}
    for target in sf_targets:
        idx = np.where(bd['SF'] <= target)[0]
        if len(idx) > 0:
            row_b[f"Time to SF={target}"] = f"{bd['T'][idx[0]]:.2f}s"
        else:
            row_b[f"Time to SF={target}"] = ">T_max"
    sf_summary_data.append(row_b)

sf_summary_df = pd.DataFrame(sf_summary_data)
sf_summary_df.to_csv('sf_threshold_summary.csv', index=False)
print(sf_summary_df)

# export 3: dose delivered in 60 seconds for each energy and mode

# Calculate dose delivered in 60 seconds for each energy and mode
dose_60s_data = []
EXPOSURE_TIME = 60.0 # seconds

for E in ENERGIES_KEV:
    dr = DR[E]
    dose_60s_data.append({
        'Energy [keV]': E,
        'Continuous Dose [Gy]': round(dr['cont_avg'] * EXPOSURE_TIME, 4),
        'Pulsed Dose [Gy]': round(dr['p_avg'] * EXPOSURE_TIME, 4),
        'Bunched Dose [Gy]': round(dr['b_avg'] * EXPOSURE_TIME, 4)
    })

dose_60s_df = pd.DataFrame(dose_60s_data)
print(f"TOTAL DOSE DELIVERED IN {EXPOSURE_TIME} SECONDS")
print(dose_60s_df)

# export 4: bunched G-factor and repair per cycle

def analytic_G_bunched_val(n, r):
    """Closed-form G for n instantaneous bunches separated by T_rep."""
    if n <= 1: return 0.0
    term1 = (2.0 * r) / (n * (1.0 - r))
    term2 = 1.0 - (1.0 - r**n) / (n * (1.0 - r))
    return term1 * term2

# Cycle-specific repair parameters
r_val = np.exp(-MU * T_REP)
repair_per_cycle_pct = (1.0 - r_val) * 100

print(f"--- INTER-BUNCH REPAIR PER CYCLE ---")
print(f"Repetition Period (T_rep): {T_REP:.2f} s")
print(f"Repair Factor (r):         {r_val:.4f}")
print(f"Repair per Cycle (%):      {repair_per_cycle_pct:.2f} %\n")

# Parameters for the detailed table
dose_points = [1.0, 5.0, 10.0]
bunched_g_data = []

for E in ENERGIES_KEV:
    dr = DR[E]
    d_bunch = dr['b_peak'] * TAU_BUNCH

    for D in dose_points:
        n_bunches = int(np.ceil(D / d_bunch))
        g_fact = analytic_G_bunched_val(n_bunches, r_val)

        bunched_g_data.append({
            'Energy [keV]': E,
            'Target Dose [Gy]': D,
            'Dose/Bunch [Gy]': round(d_bunch, 4),
            'n (Bunches)': n_bunches,
            'r (Repair Factor)': round(r_val, 4),
            'Repair/Cycle (%)': round(repair_per_cycle_pct, 2),
            'G factor (Bunched)': round(g_fact, 6)
        })

# Create and display DataFrame
bunched_g_df = pd.DataFrame(bunched_g_data)
print("BUNCHED REGIME: G-FACTOR AND CALCULATION PARAMETERS")
print(bunched_g_df)

# Export to CSV
bunched_g_df.to_csv('bunched_g_factor_parameters.csv', index=False)

#export 5: G-factor vs accumulated dose at 20 keV

# 1. Configuration for 20 keV
E_TARGET = 20
dr = DR[E_TARGET]
time_axis = t_common  # Using the shared axis from previous cells

# 2. Calculate Dose, G and SF histories
# Dose histories
dose_cont = dr['cont_avg'] * time_axis
dose_pulsed = dr['p_avg'] * time_axis
# For bunched, dose is a staircase based on integer number of pulses
n_arr = (time_axis / T_REP).astype(int)
dose_bunched = n_arr * (dr['b_avg'] / F_REP)

# G-factors
g_cont = G(time_axis)
g_pulsed = G(time_axis)

def analytic_G_bunched(n_arr, r):
    n = np.asarray(n_arr, dtype=float)
    res = np.ones_like(n)
    mask = n > 1
    term1 = (2.0 * r) / (n[mask] * (1.0 - r))
    term2 = 1.0 - (1.0 - r**n[mask]) / (n[mask] * (1.0 - r))
    res[mask] = term1 * term2
    res[n <= 1] = 0.0
    return res

g_bunched = analytic_G_bunched(n_arr, r_rep)

# SF histories (optional reference)
sf_cont = SF_curve(dose_cont, dr['cont_avg'], dr['cont_avg'])
sf_pulsed = SF_curve(dose_pulsed, dr['p_avg'], dr['p_peak'])

def SF_bunch_v(n_arr, dD_b, r, h_val):
    n = np.asarray(n_arr, dtype=float)
    D = n * dD_b
    log_rn = n * np.log(r)
    rn = np.where(log_rn > -700.0, np.exp(log_rn), 0.0)
    E2 = (2.0 * dD_b**2 * r * (n * (1.0 - r) - (1.0 - rn)) / (1.0 - r)**2)
    return np.exp(-ALPHA * D - BETA * h_val * np.maximum(E2, 0.0))

sf_bunched = SF_bunch_v(n_arr, dr['b_avg'] / F_REP, r_rep, flash_h(dr['b_peak']))

# 3. Create DataFrame and Export
g_dose_history_df = pd.DataFrame({
    'Time [s]': time_axis,
    'Dose_Continuous [Gy]': dose_cont,
    'G_Continuous': g_cont,
    'Dose_Pulsed [Gy]': dose_pulsed,
    'G_Pulsed': g_pulsed,
    'Dose_Bunched [Gy]': dose_bunched,
    'G_Bunched': g_bunched
})
g_dose_history_df.to_csv('g_vs_dose_history_20keV.csv', index=False)

# 4. Plotting
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

# Plot 1: G vs Time (Reference)
ax1.semilogx(time_axis, g_cont, label='Continuous / Pulsed', color=PALETTE[E_TARGET], lw=2.5)
ax1.semilogx(time_axis, g_bunched, label='Bunched (Step-wise)', color='black', ls='--', lw=2)
ax1.set_title(f'G-factor vs Time ({E_TARGET} keV)')
ax1.set_xlabel('Time [s] (log)')
ax1.set_ylabel('G-factor')
ax1.grid(True, which='both', linestyle=':', alpha=0.5)
ax1.legend()
ax1.set_ylim(-0.05, 1.05)

# Plot 2: G vs Accumulated Dose
ax2.plot(dose_cont, g_cont, label='Continuous', color=PALETTE[E_TARGET], lw=2.5)
ax2.plot(dose_pulsed, g_pulsed, label='Pulsed', color=PALETTE[E_TARGET], ls=':', lw=2)
ax2.step(dose_bunched, g_bunched, label='Bunched', color='black', ls='--', lw=2, where='post')
ax2.set_title(f'G-factor vs Accumulated Dose ({E_TARGET} keV)')
ax2.set_xlabel('Accumulated Dose [Gy]')
ax2.set_ylabel('G-factor')
ax2.set_xlim(0, D_MAX)
ax2.grid(True, linestyle=':', alpha=0.5)
ax2.legend()

plt.tight_layout()
plt.show()

print("G-factor vs Dose relationship data exported to 'g_vs_dose_history_20keV.csv'.")

# export 6: FLASH sparing factor summary for 20 keV

# 1. Define the range for the theoretical curve
dot_peak_axis = np.logspace(-1, 9, 1000)
h_curve = flash_h(dot_peak_axis)

# 2. Collect dots for the three beam configurations, filtered for 20 keV
E_FILTER = 20
points = []
dr = DR[E_FILTER]
# Continuous
points.append({'Type': 'DataPoint', 'Energy': E_FILTER, 'Regime': 'Continuous', 'Ḋ_peak': dr['cont_avg'], 'h': flash_h(dr['cont_avg'])})
# Pulsed
points.append({'Type': 'DataPoint', 'Energy': E_FILTER, 'Regime': 'Pulsed', 'Ḋ_peak': dr['p_peak'], 'h': flash_h(dr['p_peak'])})
# Bunched
points.append({'Type': 'DataPoint', 'Energy': E_FILTER, 'Regime': 'Bunched', 'Ḋ_peak': dr['b_peak'], 'h': flash_h(dr['b_peak'])})

# Create a DataFrame for the theoretical curve to include in the CSV
curve_df = pd.DataFrame({
    'Type': 'TheoreticalCurve',
    'Energy': np.nan,
    'Regime': 'FLASH_Function',
    'Ḋ_peak': dot_peak_axis,
    'h': h_curve
})

pts_df = pd.DataFrame(points)

# Combine data points and the curve for export
export_df = pd.concat([pts_df, curve_df], ignore_index=True)
export_df.to_csv('flash_sparing_factor_summary_20keV.csv', index=False)

# 3. Plotting
plt.figure(figsize=(12, 6))
plt.semilogx(dot_peak_axis, h_curve, color='gray', lw=2, alpha=0.5, label='FLASH Function $h(\dot{D}_{peak})$')

# Add markers for each regime (only for the selected energy)
markers = {'Continuous': 'o', 'Pulsed': 's', 'Bunched': '^'}
for regime, m in markers.items():
    subset = pts_df[pts_df['Regime'] == regime]
    plt.scatter(subset['Ḋ_peak'], subset['h'], marker=m, s=120, label=f"{regime} ({E_FILTER} keV)", edgecolors='k', zorder=5)

# Annotate the onset and saturation
plt.axvline(FLASH_ONSET, color='r', ls='--', alpha=0.3)
plt.axvline(FLASH_SAT, color='r', ls='--', alpha=0.3)
plt.text(FLASH_ONSET, 0.5, f' Onset\n({FLASH_ONSET} Gy/s)', color='r', rotation=90, va='center', ha='right', fontsize=9)
plt.text(FLASH_SAT, 0.5, f' Saturation\n({FLASH_SAT} Gy/s)', color='r', rotation=90, va='center', ha='left', fontsize=9)

plt.title(f'FLASH Sparing Factor ($h$) vs. Peak Dose Rate - {E_FILTER} keV Only')
plt.xlabel('Peak Dose Rate $\dot{D}_{peak}$ [Gy/s] (log scale)')
plt.ylabel('Sparing Factor $h$')
plt.grid(True, which='both', linestyle=':', alpha=0.6)
plt.legend()
plt.ylim(0, 1.1)
plt.tight_layout()
plt.show()

print(f"FLASH sparing factor data (points and curve) for {E_FILTER} keV exported to 'flash_sparing_factor_summary_20keV.csv'.")
print(pts_df.sort_values(by=['Regime']))

# export 7: FLASH sparing factor summary for all energies

# 1. Define the range for the theoretical curve
dot_peak_axis = np.logspace(-1, 9, 1000)
h_curve = flash_h(dot_peak_axis)

# 2. Collect dots for the three beam configurations, filtered for 20 keV
E_FILTER = 20
points = []
dr = DR[E_FILTER]
# Continuous
points.append({'Type': 'DataPoint', 'Energy': E_FILTER, 'Regime': 'Continuous', 'Ḋ_peak': dr['cont_avg'], 'h': flash_h(dr['cont_avg'])})
# Pulsed
points.append({'Type': 'DataPoint', 'Energy': E_FILTER, 'Regime': 'Pulsed', 'Ḋ_peak': dr['p_peak'], 'h': flash_h(dr['p_peak'])})
# Bunched
points.append({'Type': 'DataPoint', 'Energy': E_FILTER, 'Regime': 'Bunched', 'Ḋ_peak': dr['b_peak'], 'h': flash_h(dr['b_peak'])})

# Create a DataFrame for the theoretical curve to include in the CSV
curve_df = pd.DataFrame({
    'Type': 'TheoreticalCurve',
    'Energy': np.nan,
    'Regime': 'FLASH_Function',
    'Ḋ_peak': dot_peak_axis,
    'h': h_curve
})

pts_df = pd.DataFrame(points)

# Combine data points and the curve for export
export_df = pd.concat([pts_df, curve_df], ignore_index=True)
export_df.to_csv('flash_sparing_factor_summary_20keV.csv', index=False)

# 3. Plotting
plt.figure(figsize=(12, 6))
plt.semilogx(dot_peak_axis, h_curve, color='gray', lw=2, alpha=0.5, label='FLASH Function $h(\dot{D}_{peak})$')

# Add markers for each regime (only for the selected energy)
markers = {'Continuous': 'o', 'Pulsed': 's', 'Bunched': '^'}
for regime, m in markers.items():
    subset = pts_df[pts_df['Regime'] == regime]
    plt.scatter(subset['Ḋ_peak'], subset['h'], marker=m, s=120, label=f"{regime} ({E_FILTER} keV)", edgecolors='k', zorder=5)

# Annotate the onset and saturation
plt.axvline(FLASH_ONSET, color='r', ls='--', alpha=0.3)
plt.axvline(FLASH_SAT, color='r', ls='--', alpha=0.3)
plt.text(FLASH_ONSET, 0.5, f' Onset\n({FLASH_ONSET} Gy/s)', color='r', rotation=90, va='center', ha='right', fontsize=9)
plt.text(FLASH_SAT, 0.5, f' Saturation\n({FLASH_SAT} Gy/s)', color='r', rotation=90, va='center', ha='left', fontsize=9)

plt.title(f'FLASH Sparing Factor ($h$) vs. Peak Dose Rate - {E_FILTER} keV Only')
plt.xlabel('Peak Dose Rate $\dot{D}_{peak}$ [Gy/s] (log scale)')
plt.ylabel('Sparing Factor $h$')
plt.grid(True, which='both', linestyle=':', alpha=0.6)
plt.legend()
plt.ylim(0, 1.1)
plt.tight_layout()
plt.show()

print(f"FLASH sparing factor data (points and curve) for {E_FILTER} keV exported to 'flash_sparing_factor_summary_20keV.csv'.")
print(pts_df.sort_values(by=['Regime']))

#export 8: Dose required to reach SF=0.90 for all three source types

# We will use the previously calculated survival data to find the dose at SF=0.90
sf_target = 0.90
comparison_results = []

for E in ENERGIES_KEV:
    dr = DR[E]

    # 1. Continuous Dose at SF=0.90
    t_search = np.geomspace(1e-3, 1e5, 100000)
    sf_cont = SF_curve(dr['cont_avg'] * t_search, dr['cont_avg'], dr['cont_avg'])
    dose_cont_90 = np.interp(sf_target, sf_cont[::-1], (dr['cont_avg'] * t_search)[::-1])

    # 2. Pulsed Dose at SF=0.90
    sf_pulsed = SF_curve(dr['p_avg'] * t_search, dr['p_avg'], dr['p_peak'])
    dose_pulsed_90 = np.interp(sf_target, sf_pulsed[::-1], (dr['p_avg'] * t_search)[::-1])

    # 3. Bunched Dose at SF=0.90
    bd = _BUNCH_DATA[E]
    # Find first index where SF drops below 0.90
    idx_b = np.where(bd['SF'] <= sf_target)[0][0]
    dose_bunched_90 = bd['D'][idx_b]

    comparison_results.append({
        'Energy [keV]': E,
        'Dose Cont [Gy]': round(dose_cont_90, 4),
        'Dose Pulsed [Gy]': round(dose_pulsed_90, 4),
        'Dose Bunched [Gy]': round(dose_bunched_90, 4),
        'Add. Dose vs Cont [Gy]': round(dose_bunched_90 - dose_cont_90, 4),
        'Add. Dose vs Pulsed [Gy]': round(dose_bunched_90 - dose_pulsed_90, 4),
        '% Increase vs Cont': round((dose_bunched_90 / dose_cont_90 - 1) * 100, 2)
    })

viability_comparison_df = pd.DataFrame(comparison_results)
print(viability_comparison_df)

# Brief summary for the user
print(f"\nOn average, the bunched source requires {viability_comparison_df['% Increase vs Cont'].mean():.1f}% more dose than the continuous beam to reach 90% viability.")

# export 9: Survival fraction at 30 minutes for bunched source

t_target_s = 30 * 60  # 30 minutes in seconds
sf_30min_results = []

for E in ENERGIES_KEV:
    bd = _BUNCH_DATA[E]

    # Find the SF at 30 minutes using interpolation or finding the nearest index
    # Since the bunched SF is a step function, we find the last bunch delivered before or at 30 min
    idx = np.searchsorted(bd['T'], t_target_s, side='right') - 1
    sf_val = bd['SF'][idx]
    dose_val = bd['D'][idx]

    sf_30min_results.append({
        'Energy [keV]': E,
        'Exposure Time [min]': 30,
        'Total Bunches': int(bd['n'][idx]),
        'Accumulated Dose [Gy]': round(dose_val, 4),
        'Survival Fraction (SF)': f"{sf_val:.6e}"
    })

sf_30min_df = pd.DataFrame(sf_30min_results)

# Export for reference
sf_30min_df.to_csv('sf_at_30min_bunched.csv', index=False)

# export 10: Effective Quadratic Coefficient B = beta * h * G for all energies and modes

analysis_results = []

# We analyze the parameters over the standard 0-10 Gy range
for E in ENERGIES_KEV:
    dr = DR[E]

    for mode in ['Continuous', 'Pulsed', 'Bunched']:
        if mode == 'Continuous':
            d_avg, d_peak = dr['cont_avg'], dr['cont_avg']
        elif mode == 'Pulsed':
            d_avg, d_peak = dr['p_avg'], dr['p_peak']
        else:
            d_avg, d_peak = dr['b_avg'], dr['b_peak']

        # 1. FLASH sparing factor h (constant for a given peak dose rate)
        h_val = float(flash_h(d_peak))

        # 2. Mean G factor over the irradiation time (0 to 10 Gy)
        # G varies with time/dose, so we take the average over the 10 Gy delivery
        d_vals = np.linspace(0.01, 10, 100)
        if mode == 'Bunched':
            # Use the discrete engine results for bunched G
            g_vals, _ = discrete_SF_bunched_array(d_vals, E)
        else:
            t_vals = d_vals / d_avg
            g_vals = G(t_vals)

        g_mean = np.mean(g_vals)

        # 3. Effective Quadratic Coefficient B = beta * h * G
        # This determines the 'curvature' of the survival curve
        b_eff = BETA * h_val * g_mean

        analysis_results.append({
            'Energy [keV]': E,
            'Regime': mode,
            'Peak DR [Gy/s]': f"{d_peak:.2e}",
            'Mean G': round(g_mean, 4),
            'Sparing h': round(h_val, 3),
            'Effective Beta [Gy^-2]': f"{b_eff:.4e}",
            'Quadratic Reduction [%]': round((1 - (b_eff / BETA)) * 100, 1)
        })

analysis_df = pd.DataFrame(analysis_results)
print(analysis_df)

# Export the analysis
analysis_df.to_csv('lq_quadratic_factor_analysis.csv', index=False)

#export 11: Dose required to reach SF=0.50 for all three source types

analysis_results = []

# We analyze the parameters over the standard 0-10 Gy range
for E in ENERGIES_KEV:
    dr = DR[E]

    for mode in ['Continuous', 'Pulsed', 'Bunched']:
        if mode == 'Continuous':
            d_avg, d_peak = dr['cont_avg'], dr['cont_avg']
        elif mode == 'Pulsed':
            d_avg, d_peak = dr['p_avg'], dr['p_peak']
        else:
            d_avg, d_peak = dr['b_avg'], dr['b_peak']

        # 1. FLASH sparing factor h (constant for a given peak dose rate)
        h_val = float(flash_h(d_peak))

        # 2. Mean G factor over the irradiation time (0 to 10 Gy)
        # G varies with time/dose, so we take the average over the 10 Gy delivery
        d_vals = np.linspace(0.01, 10, 100)
        if mode == 'Bunched':
            # Use the discrete engine results for bunched G
            g_vals, _ = discrete_SF_bunched_array(d_vals, E)
        else:
            t_vals = d_vals / d_avg
            g_vals = G(t_vals)

        g_mean = np.mean(g_vals)

        # 3. Effective Quadratic Coefficient B = beta * h * G
        # This determines the 'curvature' of the survival curve
        b_eff = BETA * h_val * g_mean

        analysis_results.append({
            'Energy [keV]': E,
            'Regime': mode,
            'Peak DR [Gy/s]': f"{d_peak:.2e}",
            'Mean G': round(g_mean, 4),
            'Sparing h': round(h_val, 3),
            'Effective Beta [Gy^-2]': f"{b_eff:.4e}",
            'Quadratic Reduction [%]': round((1 - (b_eff / BETA)) * 100, 1)
        })

analysis_df = pd.DataFrame(analysis_results)
print(analysis_df)

# Export the analysis
analysis_df.to_csv('lq_quadratic_factor_analysis.csv', index=False)

#export 12: Time to reach specific survival fraction thresholds for all three source types

sf_targets = [0.99, 0.95, 0.90]
results = []

print("--- TIME TO REACH VIABILITY THRESHOLDS (ALL REGIMES) ---")

for E in ENERGIES_KEV:
    dr = DR[E]

    # 1. Continuous & Pulsed (Interpolation search)
    for mode in ['Continuous', 'Pulsed']:
        davg = dr['cont_avg'] if mode == 'Continuous' else dr['p_avg']
        dpeak = dr['cont_avg'] if mode == 'Continuous' else dr['p_peak']

        t_search = np.linspace(0, 500, 200000)
        sf_vals = SF_curve(davg * t_search, davg, dpeak)

        for target in sf_targets:
            t_target = np.interp(target, sf_vals[::-1], t_search[::-1])
            results.append({
                'Energy [keV]': E,
                'Regime': mode,
                'Target Viability': f"{target*100:.0f}%",
                'Time [s]': round(t_target, 4),
                'Dose [Gy]': round(davg * t_target, 4)
            })

    # 2. Bunched (Discrete search in pre-computed data)
    bd = _BUNCH_DATA[E]
    for target in sf_targets:
        idx = np.where(bd['SF'] <= target)[0]
        if len(idx) > 0:
            t_bunch = bd['T'][idx[0]]
            d_bunch = bd['D'][idx[0]]
        else:
            t_bunch = np.nan
            d_bunch = np.nan

        results.append({
            'Energy [keV]': E,
            'Regime': 'Bunched',
            'Target Viability': f"{target*100:.0f}%",
            'Time [s]': round(t_bunch, 2) if not np.isnan(t_bunch) else ">T_max",
            'Dose [Gy]': round(d_bunch, 4) if not np.isnan(d_bunch) else np.nan
        })

viability_summary_df = pd.DataFrame(results)
pivot_df = viability_summary_df.pivot_table(
    index=['Energy [keV]', 'Regime'],
    columns='Target Viability',
    values='Time [s]',
    aggfunc='first'
)

print(pivot_df)

#export 13: Survival fraction at specific time points for all three source types

# Updated target times to include the newly requested points
target_times = [5, 10, 30, 100, 200, 300, 400, 500, 600, 700, 800]  # seconds
sf_comparison_rows = []

for E in ENERGIES_KEV:
    dr = DR[E]
    bd = _BUNCH_DATA[E]

    for t in target_times:
        # Continuous
        sf_c = SF_curve(dr['cont_avg'] * t, dr['cont_avg'], dr['cont_avg'])
        # Pulsed
        sf_p = SF_curve(dr['p_avg'] * t, dr['p_avg'], dr['p_peak'])
        # Bunched (Step-function: find status at time t)
        idx = np.searchsorted(bd['T'], t, side='right') - 1
        # Safety check for index out of bounds if t is larger than max precomputed T
        idx = max(0, min(idx, len(bd['SF']) - 1))
        sf_b = bd['SF'][idx]

        sf_comparison_rows.append({
            'Energy [keV]': E,
            'Time [s]': t,
            'SF Continuous': round(float(sf_c), 6),
            'SF Pulsed': round(float(sf_p), 6),
            'SF Bunched': round(float(sf_b), 6)
        })

sf_time_table = pd.DataFrame(sf_comparison_rows)
sf_pivot = sf_time_table.pivot(index=['Energy [keV]'], columns='Time [s]')

# Exporting the data to an Excel file
sf_pivot.to_excel('liver_sf_comparison.xlsx')
print("Data exported to 'liver_sf_comparison.xlsx'.")

# Displaying with a pivot for a clean Energy vs Time matrix
print(sf_pivot)