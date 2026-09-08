# Numerical dosimetric analysis of low-energy positrons in living biological samples

Python code supporting the paper:

> M. Vicini, P. Folegati, G. Maero, L. Resciniti, R. Ferragut, C. Conci,
> **Numerical dosimetric analysis of low-energy positrons in living biological samples**,
> *Radiation Physics and Chemistry* **251** (2027) 114352.
> [https://doi.org/10.1016/j.radphyschem.2026.114352](https://doi.org/10.1016/j.radphyschem.2026.114352)

Part of the special issue *IFARP-7*. Open access under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

---

## What this code does

Positron Annihilation Spectroscopy (PAS) is a promising non-destructive probe of the
nanoscale free volume and electronic environment of cells and tissues. Moving PAS from
bulk, drop-cast sources to a depth-resolved *intravital* setting requires moderated,
monoenergetic positron beams in the 5–20 keV range — but the absorbed dose such beams
deliver to living matter had never been quantified.

This script implements the numerical framework of the paper. It couples positron
transport to a time-resolved radiobiological model in order to estimate the limit of
cellular viability during a PAS measurement. Specifically it:

1. **Computes the mass collision stopping power** *S(E)/ρ* of a liver sinusoidal
   endothelial cell monolayer over 5–100 keV, using the Berger–Seltzer formulation
   (ICRU Report 37) with the elemental composition of human liver tissue from ICRU
   Report 44, combined through Bragg's additivity rule.
2. **Derives absorbed dose rates** for a nominal 1.85 GBq (50 mCi) <sup>22</sup>Na
   source in three beam-delivery regimes — **continuous**, **pulsed** (electrostatic
   chopper) and **bunched** (positron trap) — tracking both time-averaged
   (*Ḋ*<sub>avg</sub>) and instantaneous peak (*Ḋ*<sub>peak</sub>) dose rates.
3. **Propagates dose through a Linear-Quadratic (LQ) survival model** extended in time
   to include sub-lethal damage repair (the Lea–Catcheside protraction factor
   *G*(μ,*t*)) and ultra-high dose rate tissue sparing (a phenomenological FLASH
   sparing factor *h*(*Ḋ*<sub>peak</sub>)):

   *SF*(*D*) = exp[ −α·*D* − β·*h*(*Ḋ*<sub>peak</sub>)·*G*(μ,*t*)·*D*² ]

4. **Exports the resulting dose and survival histories** as CSV/XLSX tables and produces
   Figs. 2–5 of the paper. (Fig. 1, the comparison of this stopping-power model against
   Gumus et al. 2018, is not generated here.)

The central question the code answers is operational: *for how long can an intravital PAS
measurement run before the surviving fraction of the inspected sample drops below the
99% / 95% / 90% viability thresholds?*

## Repository contents

| File | Description |
| --- | --- |
| `positron_dosimetry.py` | The complete model: stopping power, source dosimetry, LQ survival model, discrete time-step engine, figures and data export. Single self-contained script. |
| `requirements.txt` | Python dependencies. |
| `CITATION.cff` | Machine-readable citation metadata (GitHub renders a "Cite this repository" button from it). |
| `LICENSE` | CC BY 4.0, matching the licence of the article. |

## Requirements

Python 3.11 or later (developed and published with Python 3.12). Dependencies:

```
numpy
scipy
matplotlib
pandas
openpyxl
```

Install them with:

```bash
pip install -r requirements.txt
```

## Running the model

```bash
python positron_dosimetry.py
```

The script is non-interactive and takes no command-line arguments: all inputs are
constants in the `USER CONFIGURATION` block at the top of the file. It prints a set of
summary tables to the console, writes the data files listed below **into the current
working directory**, and opens the figures in interactive matplotlib windows.

Runtime is a few seconds: every regime is evaluated either in closed form or with a
vectorised recursion, so no step is computationally heavy.

> **Note on figures.** All figures are shown on screen with `plt.show()` and **none are
> written to disk**. Figures A, B and C (paper Figs. 3 and 4) carry `plt.savefig(...)`
> calls that are commented out — uncomment them and adjust the output path to save those
> three. The *G*-factor and FLASH figures (paper Figs. 5 and 2) have no `savefig` call at
> all; add one if you need the image files. Be aware that the console summary prints a
> `FIGURES SAVED:` block listing three PNG filenames even though nothing is saved.
> To run the script on a headless machine, set the non-interactive backend first:
> `MPLBACKEND=Agg python positron_dosimetry.py` (the data exports are unaffected).

## Configuration

The physical and radiobiological inputs are grouped at the top of the script. The values
committed here are exactly those used for the published results:

| Parameter | Value | Meaning |
| --- | --- | --- |
| `ENERGIES_KEV` | `[5, 10, 15, 20, 50, 100]` | Positron implantation energies [keV] |
| `A_NA22` | `1.85e9` | <sup>22</sup>Na source activity [Bq] (50 mCi) |
| `E_MAX_NA22`, `E_MEAN_NA22` | `545.7`, `215.5` | <sup>22</sup>Na β<sup>+</sup> endpoint and mean kinetic energy [keV] (reported in the console summary) |
| `BR_POS` | `0.9034` | β<sup>+</sup> branching ratio |
| `ETA_MOD` | `1.2e-3` | Tungsten-mesh moderation efficiency |
| `ETA_TRANS` | `0.50` | Electrostatic transport efficiency |
| `F_CHOP_HZ`, `TAU_ON_S` | `50e6`, `2e-9` | Chopper frequency [Hz] and ON-window [s] → geometric transmission 0.10 |
| `N_BUNCH`, `TAU_BUNCH`, `ETA_BUNCH` | `1e7`, `10e-9`, `0.1` | Positrons per bunch, compressed bunch duration [s], trap efficiency → *T*<sub>rep</sub> ≈ 100 s |
| `FWHM_MM` | `2.0` | Gaussian beam FWHM [mm] → exposed area π mm² |
| `ALPHA`, `BETA` | `0.10`, `0.05` | LQ coefficients [Gy⁻¹], [Gy⁻²] → α/β = 2 Gy (late-responding normal tissue, hepatic) |
| `T_HALF_S` | `1800` | Molecular repair half-time of liver cells [s] → μ = 3.85 × 10⁻⁴ s⁻¹ |
| `FLASH_ONSET`, `FLASH_SAT`, `FLASH_H_MIN` | `40`, `120`, `0.10` | FLASH sparing: onset and saturation peak dose rate [Gy/s], and maximum sparing |
| `D_MAX`, `D_THRESHOLD` | `10.0`, `1.0` | Upper dose bound and reference dose of the plots [Gy] |

Adding energies to `ENERGIES_KEV` propagates automatically through the tables, the
exports and the figures (colours and line styles are generated from the list). **Do not
remove 20 keV from the list**, however: the single-energy analyses (`E_TARGET` /
`E_FILTER`, the two exports named `*_20keV.csv`) look that value up directly and the
script would stop with a `KeyError` partway through the exports.

Note also that the worked numerical example in the module docstring dates from an earlier
configuration with *T*<sub>rep</sub> = 120 s (*r* = 0.955, μ*T*<sub>rep</sub> = 0.046).
The parameters committed here give *T*<sub>rep</sub> = 99.7 s, *r* = 0.962 and
μ*T*<sub>rep</sub> = 0.038 — these are the values used in the published results, and the
ones the script prints at run time. The derivations in that docstring are unaffected.

## Generated data files

Running the script produces the following files in the working directory. They are the
numerical content behind the paper's figures and tables:

| Output file | Content | Related paper item |
| --- | --- | --- |
| `dose_rates_summary.csv` | Stopping power and average/peak dose rates per energy, for the three regimes | **Table 1** |
| `time_to_1Gy_summary.csv` | Wall-clock time to deposit 1 Gy, per energy and regime | Sec. 3.1 |
| `dose_history_continuous.csv`<br>`dose_history_pulsed.csv`<br>`dose_history_bunched.csv` | Accumulated dose vs. time on a shared logarithmic time axis (10⁴ points); the bunched file is a forward-filled staircase | **Fig. 4** |
| `sf_history_continuous.csv`<br>`sf_history_pulsed.csv`<br>`sf_history_bunched.csv` | Surviving fraction vs. time on the same shared axis | **Fig. 3** |
| `sf_threshold_summary.csv` | Time required to cross *SF* = 99/95/90/50/10% | Sec. 3.2–3.3 |
| `g_vs_dose_history_20keV.csv` | Lea–Catcheside factor *G* and accumulated dose vs. time at 20 keV, three regimes | **Fig. 5** |
| `bunched_g_factor_parameters.csv` | Dose per bunch, number of bunches, inter-bunch repair factor *r* and resulting *G* at 1/5/10 Gy | Sec. 2.3.1, 3.3 |
| `flash_sparing_factor_summary_20keV.csv` | *h*(*Ḋ*<sub>peak</sub>) curve plus the operating point of each regime at 20 keV | **Fig. 2** |
| `lq_quadratic_factor_analysis.csv` | Effective quadratic coefficient β·*h*·⟨*G*⟩ and its reduction relative to β | Sec. 3.3 |
| `sf_at_30min_bunched.csv` | Bunches delivered, dose and *SF* after 30 min of bunched irradiation | Sec. 3.3 |
| `liver_sf_comparison.xlsx` | *SF* matrix (energy × elapsed time) for the three regimes at 5–800 s | **Table 2** |

## Model outline

**Stopping power.** `stopping_power(E_keV)` evaluates the Berger–Seltzer mass collision
stopping power (ICRU 37) with the Rohrlich–Carlson positron-specific correction
*F*<sup>+</sup>(τ), for a composite liver medium whose ⟨*Z/A*⟩ and mean excitation energy
⟨*I*⟩ are obtained from ICRU 44 mass fractions via Bragg's additivity rule.

**Source dosimetry.** `dose_rates(E_keV)` converts the useful positron rate at the sample
*R*<sub>sample</sub> = *A* · BR · η<sub>mod</sub> · η<sub>trans</sub> into a fluence rate
over the beam area, and then into dose rates for each regime. Note that an active chopper
acts purely as a geometric shutter and *cannot* compress the beam, so
*Ḋ*<sub>peak</sub><sup>pulsed</sup> = *Ḋ*<sub>cont</sub> while the average is scaled by the
duty cycle. The trap-based bunched regime, in contrast, compresses 10⁷ positrons into a
10 ns burst and reaches *Ḋ*<sub>peak</sub> ~ 10⁷–10⁸ Gy/s.

**Protraction factor.** `G(T_s)` is the closed-form Dale (1985) expression
*G*(*t*) = 2[μ*t* − 1 + exp(−μ*t*)]/(μ*t*)², with a second-order Taylor expansion below
μ*t* < 10⁻⁸ where the analytic form becomes numerically unstable.

**FLASH sparing.** `flash_h(ddot_peak)` implements the piecewise-linear sparing factor:
*h* = 1 below 40 Gy/s, decreasing linearly to *h*<sub>min</sub> = 0.10 at 120 Gy/s, and
constant thereafter.

**Discrete time-step engine.** For the bunched regime the dose arrives in discrete
instantaneous increments, so *G* cannot be taken from the continuous formula. The script
integrates the sub-lethal-damage memory variable *w* with the exact recursion
*w*<sub>*n*+1</sub> = e<sup>−μΔ*t*</sup>(*w*<sub>*n*</sub> + Δ*D*<sub>*n*</sub>) and
accumulates *E*₂ += 2·Δ*D*<sub>*n*</sub>·*w*<sub>*n*</sub>, from which
*G* = *E*₂/*D*². This reproduces the analytic *G* to machine precision for a constant
dose rate, and is also available as a closed form for *n* equal instantaneous bunches
(`SF_bunch_fast`). The extended docstring at the top of the script derives these
expressions and discusses the choice of time step.

## Main result

Across all three regimes the achievable measurement duration scales as 1/*S*(*E*), and the
absorbed dose is governed by the **time-averaged** flux reaching the sample — essentially
decoupled from the temporal microstructure of the beam. The continuous and pulsed
configurations differ by at most one order of magnitude in delivered dose, and the
99% viability boundary (*D* ≈ 100 mGy) is crossed within seconds: 0.44 s at 5 keV to
4.4 s at 100 keV for the continuous beam, roughly one decade later for the pulsed one.

The bunched regime is qualitatively different. It activates two protective mechanisms at
once — FLASH sparing at *h* = *h*<sub>min</sub> = 0.10, since the peak dose rate reaches
10⁷–10⁸ Gy/s, and significant sub-lethal repair across the ≈100 s inter-bunch interval —
which together roughly double the effective *LD*<sub>50</sub>, from ≈2.86 Gy to 5.7–5.98 Gy.
That protection comes at the price of dose *quantization*: each bunch deposits a minimum
increment of 216 mGy (100 keV) to 2176 mGy (5 keV), so *SF* = 99% is unreachable at any
energy and the 95% threshold is met only marginally, at 50 and 100 keV. In the 5–20 keV
window relevant to depth-resolved PAS, no regime satisfies the strict *SF* ≥ 95%
non-destructive criterion.

The practical conclusion is that sample viability is set by the time-averaged flux, not by
the choice of time structure. Collecting the ~10⁶ annihilation events needed for a
statistically robust lifetime spectrum before the surviving fraction drops below 90%
therefore calls for a 30–50 mCi <sup>22</sup>Na source together with a digital detection
system reaching at least 10⁵ counts/s — which frames depth-resolved PAS of living cells as
an instrumental challenge for next-generation positron beam systems.

## Citation

If you use this code, please cite the article:

```bibtex
@article{Vicini2027,
  title   = {Numerical dosimetric analysis of low-energy positrons in living biological samples},
  author  = {Vicini, Matteo and Folegati, Paola and Maero, Giancarlo and
             Resciniti, Leonardo and Ferragut, Rafael and Conci, Claudio},
  journal = {Radiation Physics and Chemistry},
  volume  = {251},
  pages   = {114352},
  year    = {2027},
  doi     = {10.1016/j.radphyschem.2026.114352},
  issn    = {0969-806X}
}
```

## License

The code in this repository is released under the
[Creative Commons Attribution 4.0 International licence (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/),
the same licence as the article. See [`LICENSE`](LICENSE).

## Acknowledgements

This research is funded by the European Research Council (ERC), project **ALFRED**,
G.A. 101221323. Views and opinions expressed are those of the authors only and do not
necessarily reflect those of the European Union or the European Research Council.

## Contact

Matteo Vicini — matteo.vicini@polimi.it
Department of Chemistry, Materials and Chemical Engineering "G. Natta", Politecnico di Milano

Rafael Ferragut — rafael.ferragut@polimi.it · Claudio Conci — claudio.conci@polimi.it
