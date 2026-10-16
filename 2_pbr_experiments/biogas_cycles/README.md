# Biogas cycling experiment

Data and analysis for three countercurrent synthetic-biogas cycles on 31 May 2023, using a ceria packed bed at 930 °C (1203 K). Reduction uses CH4/CO2 diluted in argon; oxidation uses CO2 diluted in argon.

| File | Contents / purpose |
| --- | --- |
| `rx_test_027.lvm` | Original instrument data used by `plot_data.py`: tab-separated measurements with decimal commas and no header. |
| `plot_data.py` | Extracts the required instrument channels, corrects the H2 signal for background and CO/CO2 cross-sensitivity, writes the processed CSV, and plots the experiment with shaded reduction/oxidation intervals. |
| `31052023_biogas_930C_dataset.csv` | Processed time series: time (s), bed temperature (°C), pressure relative to ambient (bar), inlet stream flows and mass-flow-controller set points (standard L/min), and measured gas concentrations (%). CH4 and CO2 stream flows include their argon diluent. |
| `analysis_biogas.py` | Reads the processed CSV, identifies the three reduction/oxidation cycles, integrates molar flows, and prints conversions, yields, elemental balances, carbon-balance deficits, and mean per-cycle feed amounts for model inputs. Also saves an overview plot. |
| `analysis_biogas_output.txt` | Saved console output from the analysis, including per-cycle values and three-cycle totals. |
| `analysis_biogas_summary.csv` | Compact metric/value summary of the analysis. Inlet gas amounts are totals over three cycles; `n_CeO2_mol` is the bed inventory. This file is not automatically written by the current analysis script. |
| `reactor_bed` | LaTeX table fragment listing bed dimensions, ceria loading (292.3 g), operating temperature and pressure, and void fractions. |
| `plots/31052023_Biogas_930C_dataset.png` | Raster overview of gas concentrations, inlet flows, and bed temperature. |
| `plots/31052023_Biogas_930C_dataset.pdf` | PDF version of the overview plot. |

The `.idea/` directory contains editor settings: `.name` (project name), `biogas_cycles.iml` (module definition), `modules.xml` (module list), `misc.xml` (project configuration), `workspace.xml` (local workspace state), `inspectionProfiles/profiles_settings.xml` (inspection settings), and `.gitignore` (editor-file exclusions). These are not experimental inputs.

Run the scripts from this directory with Python and NumPy, pandas, Matplotlib, and SciPy installed. Use `plot_data.py` to regenerate the processed data, then `analysis_biogas.py` for the integrated results. Both scripts write the same plot filenames, so the last script run determines the saved plot style. The analysis prints its results; saving the text output requires console redirection or capture.
